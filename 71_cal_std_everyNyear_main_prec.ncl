load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/kf_filter.ncl"
begin

;============================================================================
; Step 1: Settings of variables
;============================================================================

; path and variable settings
  parent_dir     = "./data/" 					; parent directory of all model members' folders
  sub_dir        = "/day/" 							; the path between the .nc file and the corresponding model folder
  read_var       = "U850" 									; name of the desired variable in CESM2 .nc files
  output_var     = "U850_stddev" 								; name of the output variable for the output .nc file
  diro           = "./" 									; output directory

; domain as the input of kf_filtet (fed grids)
  lat_n     =  25.
  lat_s     = -25.

; kf_filter settings
  mis            = 0
  mis@_FillValue = 0 										; just for searching all equivalent depth wave signals

  tMin      = 30  										; 20-100 days
  tMax      = 80
  kMin      = 1   										; wave number from 1 - 5
  kMax      = 5
  hMin      = mis
  hMax      = mis
  waveName  = "none"
  obsPerDay = 1   										; 1 time step for a day as daily data

; Make sure the folder structure is established
  detect_subdir1 = systemfunc("ls "+diro+"everyNyear")
  if (all(ismissing(detect_subdir1)))
    system("mkdir "+diro+"everyNyear")
  end if
  detect_subdir2 = systemfunc("ls "+diro+"everyNyear/"+N)
  if (all(ismissing(detect_subdir2)))
    system("mkdir "+diro+"everyNyear/"+N)
  end if

  diro = diro+"everyNyear/"+N+"/"

;============================================================================
; Step 2: Find names and paths of the .nc files for each model
;============================================================================

  models_hist_dirs  = systemfunc("ls -d "+parent_dir+"b.e21.BHIST*.*.LE2-????.???") 		; model folders' paths ---- ONLY for Historical

  nmodels           = dimsizes(models_hist_dirs)
  models_hist_names = new(nmodels,string) 							; model folders' names
  models_ssp_names  = new(nmodels,string) 							; model folders' names
  models_names      = new(nmodels,string) 							; model names for output ---- diro/filo
  
  do i = 0,nmodels-1
    str                  = str_split(models_hist_dirs(i),"/")
    models_hist_names(i) = str(dimsizes(str)-1) 						; derive model folders' names only
    str_ssp              = str_split(models_hist_names(i),"Bsc")
    cmip_or_smbb_detect  = str_split(str_ssp(dimsizes(str_ssp)-1),".")
  
    if(str_ssp(1).eq."HIST")
       str_ssp(1) = "SSP370"
       if(cmip_or_smbb_detect(0).eq."mbb")
          models_ssp_names(i) = str_ssp(0) + "B" + str_ssp(1) + "s" + str_ssp(2)
          models_names(i)     = "s" + str_ssp(2)
          else if(cmip_or_smbb_detect(0).eq."mip6")
                  models_ssp_names(i) = str_ssp(0) + "B" + str_ssp(1) + "c" + str_ssp(2)
                  models_names(i)     = "c" + str_ssp(2)
          end if
       end if
    end if
  end do

  models_ssp_dirs = parent_dir + models_ssp_names
  print("")
  print("This thread is processing the model: " + models_names(current_nmodel))
  print("")
  print("=============================================")

  days_per_year = 365 										; 365 days per year (without existence of leap years in the model outputs)
  years_per_nc  = 10 										; each nc file contains 10 years
  nc_hist_dir   = models_hist_dirs(current_nmodel) + sub_dir 					; folder containing all nc files for currently-picked  model
  nc_ssp_dir    = models_ssp_dirs(current_nmodel)  + sub_dir 					; folder containing all nc files for currently-picked  model
  nc_hist_files = systemfunc("ls "+nc_hist_dir+"*U850.*.nc")
  nc_ssp_files  = systemfunc("ls "+nc_ssp_dir+"*U850.*.nc")
  nc_files      = array_append_record(nc_hist_files,nc_ssp_files,0)
  
;============================================================================
; Step 3: Judge HIST or SSP to choose the corresponding start/end year
        ; and to skip nc files not covering the desired years
;============================================================================

  nncfiles       = dimsizes(nc_files)
  nskipped_files = 0

  do i = 0,nncfiles-1
     str := str_split(nc_files(i),"/")
     str := str_split(str(dimsizes(str)-1),".")
     str := str(dimsizes(str)-2) 									; year duration, e.g., 18500101-18591231
     year_str      = str_split_by_length(str,(/4,4,1,4,4/))
     nc_start_year = year_str(0) 									; e.g., 1850 in the same example
     nc_end_year   = year_str(3) 									; e.g., 1859 in the same example

     if nc_end_year.lt.start_year .or. nc_start_year.gt.end_year
        print("nc file skipped (not in the desired years: "+start_year+"-"+end_year+").")
        nskipped_files = nskipped_files + 1
        continue
     end if

;============================================================================
; Step 4: Generate time points based on the stary/end years of the .nc files
        ; will be used later in Step 5
;============================================================================

     print("=============================================")
     print("Reading years: "+nc_start_year+"-"+nc_end_year)

     years_skipped_head = 0 									; for most of the nc files whose years are fully covered within the desired years ; default
     years_skipped_tail = 0 									; same as above but for detecting ending year
 
     if (nc_start_year.lt.start_year)
        years_skipped_head = start_year - toint(nc_start_year) 					; for nc files with years only partly covered
        if (years_skipped_head.ne.0)
            max_years_skipped_head = years_skipped_head 					; record the only non-zero years_skipped_head
        end if
     end if
  
     if (nc_end_year.gt.end_year)
        years_skipped_tail = toint(nc_end_year) - end_year 					; for nc files with years only partly covered
        if (years_skipped_tail.ne.0)
            max_years_skipped_tail = years_skipped_tail 					; record the ony non-zero years_skipped_tail
        end if
     end if

     if (toint(nc_end_year)-toint(nc_start_year).ne.years_per_nc-1)
        if (toint(nc_end_year).ne.2100)
            print("Success: HIST/SSP switched. (HIST nc-files reading finished)")
        end if
        years_skipped_tail = years_per_nc-1 - (toint(nc_end_year)-toint(nc_start_year))		; default years_per_nc = 10 but there's ONE nc files that does not have these many years, 
												; so I manually adjust its years_per_nc to be 10 to calculate the "year_skipped_tail" 
												; (won't affect reading results but just easier to process it here.).
     end if
  
     print("head skipped years: "+years_skipped_head+" ; tail skipped years: "+years_skipped_tail)
  
     doy1 = day_of_year((/2021,2021/),(/ 1,11/),(/ 1, 1/)) 					; Starting Julian days (boreal winter: 01/01 & 11/01)
     doy2 = day_of_year((/2021,2021/),(/ 4,12/),(/30,31/)) 					; Ending Julian days (boreal winter: 04/30 & 12/31)
  
     time_points_start := years_skipped_head * days_per_year + doy1
     time_points_end   := years_skipped_head * days_per_year + doy2

     j = 1
     do while(years_skipped_head+years_skipped_tail+j.le.years_per_nc-1)
        time_points_start := array_append_record(time_points_start, (years_skipped_head+j)* days_per_year + doy1,0)
        time_points_end   := array_append_record(time_points_end  , (years_skipped_head+j)* days_per_year + doy2,0)
        j = j+1
     end do

     time_points_start = time_points_start -1 							; NCL starts from 0 (but Julian days start from 1)
     time_points_end   = time_points_end -1 							; NCL starts from 0 

     if (.not.isdefined("time_points"))
        time_points_nc_shift = 0
        time_points         := ispan(time_points_nc_shift+time_points_start(0),time_points_nc_shift+time_points_end(0),1)
     else
        total_years_read     = toint(nc_start_year) - previous_nc_start_year
        time_points_nc_shift = time_points_nc_shift + total_years_read*days_per_year 
        time_points         := array_append_record(time_points,ispan(time_points_nc_shift+time_points_start(0),time_points_nc_shift+time_points_end(0),1),0)
     end if

     previous_nc_start_year  = toint(nc_start_year)

     j = 1
     do while(j.lt.dimsizes(time_points_start))
        time_points := array_append_record(time_points,ispan(time_points_nc_shift+time_points_start(j),time_points_nc_shift+time_points_end(j),1),0)
        j = j+1
     end do

;============================================================================
; Step 5: read the variable according for the time points and specified lat range first
;============================================================================

     f = addfile(nc_files(i),"r")
     
     if (.not.isdefined("var"))
        var := f->$read_var$(time_points_start(0):time_points_end(dimsizes(time_points_end)-1),{lat_s:lat_n},:)
     else
        var := array_append_record(var,f->$read_var$(time_points_start(0):time_points_end(dimsizes(time_points_end)-1),{lat_s:lat_n},:),0) 			; attach the var array along the first dimension (time)
     end if
  end do

  print("=============================================")
  printVarSummary(var)
  time_points = time_points - time_points(0) 							; when reading var, it starts from here

;============================================================================
; Step 6: Then apply the kf_filter and refine into output
;============================================================================
 
  lat       = var&lat
  lon       = var&lon
  filtData  = var
  do y = 0, dimsizes(lat)-1
     filtData(:,{lat(y)},:) = kf_filter(var(:,{lat(y)},:), obsPerDay, tMin, tMax, kMin, kMax, hMin, hMax, waveName )
  end do

  printVarSummary(filtData)

  varo     = filtData(time_points,:,:)
  delete (filtData)
  delete (var)
  day      = time_points
  day!0    = "day"
  varo!0   = "day"
  varo&day = day

;============================================================================
; Step 7: Calculate the STDDEV and then output
;============================================================================

  varo             := dim_stddev_n_Wrap(varo,0)
;  varo              = varo * 86400 * 1000 ; Convert units meter/second to mm/day
  varo@units        = "mm/day"
  varo@long_name    = "Standard deviation of "+varo@long_name+" at each grid filtered by kf_filter"
  delete(varo@lat)
  varo@cell_methods = "daily data from CESM2L, period:"+start_year+"-"+end_year
  printVarSummary(varo)

  nlat = dimsizes(lat)
  nlon = dimsizes(lon)
  
  ;-------------------

  print("Now is dealing with creating a new nc. file.")

  nclpath = systemfunc("pwd")
  nclpath = nclpath + "/" + systemfunc("ls *.ncl -tr |tail -n 1")

  ;-------------------

  filo= output_var+ "." + start_year+"-"+end_year + ".nc"
  if (all(ismissing(systemfunc("ls "+diro+models_names(current_nmodel)))))
    system("mkdir "+diro+models_names(current_nmodel))    					; create subdir to accomondate many nc files
  end if
  
  system("/bin/rm -f " + diro +models_names(current_nmodel) + "/" + filo)    			; remove if exists
  fout  = addfile (diro + models_names(current_nmodel) + "/" + filo, "c")  			; open output file

  ;-------------------    
  setfileoption(fout,"DefineMode",True)
  ;-------------------
      fAtt               = True   							         ; assign file attributes
      fAtt@title         = "NCL Efficient Approach to netCDF Creation"
      fAtt@source_file   =  "original-file.nc"
      fAtt@Conventions   = "None"
      fAtt@nclpath       = nclpath
      fAtt@creation_date = systemfunc ("date")
      fileattdef( fout, fAtt )            							; copy file attributes
  ;-------------------    
      dimNames           = (/"lat",  "lon"/)
      dimSizes           = (/ nlat,  nlon /)
      dimUnlim           = (/False,  False/)
      filedimdef(fout,dimNames,dimSizes,dimUnlim)
  ;-------------------
     filevardef(fout, "lat"  ,typeof(lat),getvardims(lat))
     filevardef(fout, "lon"  ,typeof(lon),getvardims(lon))
     filevardef(fout, output_var ,typeof(varo)  ,getvardims(varo))
  ;-------------------   
     filevarattdef(fout,output_var,varo)
     filevarattdef(fout,"lat"     ,lat)
     filevarattdef(fout,"lon"     ,lon)
  ;-------------------
      setfileoption(fout,"DefineMode",False)
  ;-------------------
     fout->lat           = (/lat/)
     fout->lon           = (/lon/)
     fout->$output_var$  = (/varo/)

  print("NCL file path: "+nclpath)
  print(".nc file path: "+diro+filo)
 
end
