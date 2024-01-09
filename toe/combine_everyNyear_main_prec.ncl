;###########################################################################
;---Author: Yi-Xian Li, modified by Hien Bui
;---Goal: Calculate time of emergence of MJO amplitude
;---Input: MJO precipitation or wind amplitude timeseries
;---For details, please check our paper: https://doi.org/10.1038/s41612-023-00344-z
;---Please send me an email if you have any questions: hien.bui@monash.edu
;###########################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/kf_filter.ncl"
begin

;============================================================================
; Step 1: Settings of variables
;============================================================================

; path and variable settings
  parent_dir     = "./data/" 				; parent directory of all model members' folders
  sub_dir        = "/day/" 						; the path between the .nc file and the corresponding model folder
  read_var       = "U850_stddev" 							; name of the output variable for the output .nc file
  diro           = "./" 								; output directory
  output_var     = "U850_stddev"

; domain as the input of kf_filtet (fed grids)
  lat_n          =  25.
  lat_s          = -25.

;============================================================================
; Step 2: Find names and paths of the .nc files for each model
;============================================================================

  models_hist_dirs  = systemfunc("ls -d "+parent_dir+"b.e21.BHIST*.*.LE2-????.???") 	; model folders' paths ; ONLY for Historical

  nmodels           = dimsizes(models_hist_dirs)
  models_hist_names = new(nmodels,string) 						; model folders' names
  models_ssp_names  = new(nmodels,string) 						; model folders' names
  models_names 	    = new(nmodels,string) 						; model names for output ; diro/filo

  do i = 0,nmodels-1
     str                  = str_split(models_hist_dirs(i),"/")
     models_hist_names(i) = str(dimsizes(str)-1) 					; derive model folders' names only
     str_ssp              = str_split(models_hist_names(i),"Bsc")
     cmip_or_smbb_detect  = str_split(str_ssp(dimsizes(str_ssp)-1),".")

     if (str_ssp(1).eq."HIST")
        str_ssp(1) = "SSP370"
        if (cmip_or_smbb_detect(0).eq."mbb")
           models_ssp_names(i) = str_ssp(0) + "B" + str_ssp(1) + "s" + str_ssp(2)
           models_names(i)     = "s" + str_ssp(2)
           else if (cmip_or_smbb_detect(0).eq."mip6")
                   models_ssp_names(i) = str_ssp(0) + "B" + str_ssp(1) + "c" + str_ssp(2)
                   models_names(i)     = "c" + str_ssp(2)
           end if
        end if
     end if
  end do

  models_ssp_dirs = parent_dir + models_ssp_names
  
;============================================================================
; Step 3: Combining
;============================================================================

  print("")
  print("This thread is processing the model: " + models_names(current_nmodel))

  dir  = diro+"everyNyear/"+N+"/" + models_names(current_nmodel)
  diro = diro+"everyNyear/"+N+"/"

  print("Directory: "+diro)
  print("")
  print("=============================================")

  files  = systemfunc("ls "+dir+"/*.nc")
  print(files)
  nfiles = dimsizes(files)
  f      = addfiles(files,"r")
  ListSetType (f, "join")
  var    = f[:]->$read_var$

  ending_year = new(nfiles,integer) 							; model names for output ; diro/filo

  do i = 0,nfiles-1
     str := str_split(files(i),".-")
     if (i.eq.0)
        start_year =  str(dimsizes(str)-3)
     end if
     if (i.eq.nfiles-1)
        end_year =  str(dimsizes(str)-2)
     end if
     ending_year(i) = toint(str(dimsizes(str)-2))
  end do
    
  print("Years: "+start_year+"-"+end_year)

  ending_year!0   = "ending_year"
  var!0           = "ending_year"
  var&ending_year = ending_year
  printVarSummary(var)

;===================================================================
; Step 4: Output .nc files
;===================================================================

  print("Now is dealing with creating a new nc. file.")

  nclpath = systemfunc("pwd")
  nclpath = nclpath + "/" + systemfunc("ls *.ncl -tr |tail -n 1")
  
;===================================================================

  lat   = var&lat
  lon   = var&lon

  nlat  = dimsizes(lat)
  nlon  = dimsizes(lon)
  nyear = dimsizes(ending_year)

  filo  = models_names(current_nmodel)+"."+start_year+"-"+end_year+".nc"

  system("/bin/rm -f " + diro + filo)  							; remove if exists
  fout  = addfile (diro + filo, "c")  							; open output file

  ; ------------------
      setfileoption(fout,"DefineMode",True)
  ; ------------------    
      fAtt               = True            						; assign file attributes
      fAtt@title         = "NCL Efficient Approach to netCDF Creation"
      fAtt@source_file   =  "original-file.nc"
      fAtt@Conventions   = "None"
      fAtt@nclpath       = nclpath
      fAtt@creation_date = systemfunc ("date")
      fileattdef( fout, fAtt )           						; copy file attributes
  ; ------------------    
      dimNames = (/"ending_year", "lat",  "lon"/)
      dimSizes = (/ nyear,         nlat,   nlon/)
      dimUnlim = (/ False,        False,  False/)
      filedimdef(fout,dimNames,dimSizes,dimUnlim)
  ; ------------------
     filevardef(fout, "ending_year" ,typeof(ending_year),getvardims(ending_year))
     filevardef(fout, "lat"         ,typeof(lat)        ,getvardims(lat))
     filevardef(fout, "lon"         ,typeof(lon)        ,getvardims(lon))
     filevardef(fout, output_var    ,typeof(var)        ,getvardims(var))
  ; ------------------
     filevarattdef(fout,output_var   ,var)
     filevarattdef(fout,"ending_year",ending_year)
     filevarattdef(fout,"lat"        ,lat)
     filevarattdef(fout,"lon"        ,lon)
  ; ------------------
     setfileoption(fout,"DefineMode",False)
  ; ------------------
     fout->ending_year  = (/ending_year/)
     fout->lat          = (/lat/)
     fout->lon          = (/lon/)
     fout->$output_var$ = (/var/)

  print("NCL file path: "+nclpath)
  print(".nc file path: "+diro+filo)

end
