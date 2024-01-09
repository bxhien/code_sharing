;###########################################################################
;---Author: Yi-Xian Li, modified by Hien Bui
;---Goal: Calculate time of emergence of MJO amplitude
;---Input: MJO precipitation or wind amplitude timeseries
;---For details, please check our paper: https://doi.org/10.1038/s41612-023-00344-z
;---Please send me an email if you have any questions: hien.bui@monash.edu
;###########################################################################

begin

;============================================================================
; Step 0: Settings of variables
;============================================================================

; path and variable settings
  dir       = "./"
  dir       = dir + N + "/"
  read_var  = "U850_stddev" 							;stddev = Standard Deviation ; name of the output variable for the output .nc file
  diro      = "./output/" 	; output directory ; make sure you create this folder in advance


; year settings
  start_year = 1861
  end_year   = 2100
  pre_year   = 2005 ; the last year of definition of current climate (i.e., start_year to pre_year; e.g., 1861-2005)

  consecutive_year = (/1,3,10,1/); criteria of consecutive years exceeding the upper bound ; last "1" means for requiring all subsequent years to pass

;============================================================================
; Step 1: Find names of the .nc files for each model
;============================================================================

  files  = systemfunc("ls "+dir+"*"+start_year+"-"+end_year+".nc") ; model folders' paths ; ONLY for Historical
  print(files)

; Make sure the folder structure is established
  detect_subdir1 = systemfunc("ls "+diro+"runningNyear")
  if (all(ismissing(detect_subdir1)))
    system("mkdir "+diro+"runningNyear")
  end if
  detect_subdir2 = systemfunc("ls "+diro+"runningNyear/"+N)
  if (all(ismissing(detect_subdir2)))
    system("mkdir "+diro+"runningNyear/"+N)
  end if
  diro = diro+"runningNyear/"+N+"/"

  models_nc_names = files
  nmodels = dimsizes(models_nc_names)
  do i = 0,nmodels-1
    str = str_split(files(i),"/")
    models_nc_names(i) = str(dimsizes(str)-1)
  end do
  print("This thread is processing the model: " + models_nc_names(current_nmodel))
  f = addfile(files(current_nmodel),"r")

  var    = f->$read_var$
  year   = f->ending_year - N/2 ; middle year of the running years
  year!0 = "year"

;============================================================================
; Step 2: Define historical bounds (natural variability)
;============================================================================

  ts      = var
  ts!0    = "year"
  ts&year = year
  printVarSummary(ts)

  lat  = ts&lat
  lon  = ts&lon
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)
  nconsecutive_year = dimsizes(consecutive_year)
  sz  = dimsizes(ts)

  var := var({start_year:pre_year},:,:)
  do i = 0, nlat-1
    do j = 0, nlon-1
      qsort(var(:,i,j))
    end do
  end do

  p = 95 ; percentile
  ind_p = dimsizes(var(:,0,0)) * (p / 100.)     ; the location of the percentile wanted
  fl = ind_p - toint(ind_p)			; float part
  if fl.gt.0.5
    ind_p := toint(ind_p) 
    ub = var(ind_p,:,:)
  else if fl.lt.0.5
    ind_p := toint(ind_p) - 1
    ub = var(ind_p,:,:)
  else ; if.fl.eq.0.5
    ind_p := toint(ind_p) - 1
    ub = ( var(ind_p,:,:) + var(ind_p+1,:,:) ) /2.
  end if
  end if

/;
  mean = dim_avg_n_Wrap(var({start_year:pre_year},:,:),0) ; mean of input_var within current climatology
  sd = dim_stddev_n_Wrap(var({start_year:pre_year},:,:),0) ; standard deviation of input_var within current climatology
  ub = mean + 1*sd ; upper bound of current climatology ; critical value (tipping point) ;<----- currently use 1SD
  copy_VarCoords(mean,ub)
;/


;============================================================================
; Step 3: Calculate (4) and (5)
;============================================================================

  ToE = new((/nconsecutive_year,sz(1),sz(2)/),float)
  min_ts = new(sz,float)
  do i = 0 ,nconsecutive_year-1
    do y = 0 , dimsizes(year)-1 - (consecutive_year(i)-1)
      if (i.ne.nconsecutive_year-1)
          min_ts(y,:,:) = dim_min_n_Wrap (ts(y:y+consecutive_year(i)-1,:,:), 0) ; min of [y,y+consecutive_year(i)-1] years ; test for if this > ub
      else ; for consecutive_year.eq.-99 ; last case (all the years after the year must exceed ub)
          min_ts(y,:,:) = dim_min_n_Wrap (ts(y:,:,:), 0) ; min of [y,end_year] years ; test for if this > ub
      end if
    end do

    do j = 0,nlat-1
      do k = 0,nlon-1
        gt_yr := N*ind(min_ts(:,j,k).gt.ub(j,k)) + year(0) ; all the years > ub
        gt_yr_future_ind := ind(gt_yr.gt.pre_year) ; ind of years after pre_year (future climatology, 2005 as default)
        if(.not.all(ismissing(gt_yr_future_ind)))
          ToE(i,j,k) = gt_yr(gt_yr_future_ind(0))
        end if
      end do
    end do
  end do
  consecutive_year!0 = "consecutive_year"
  ToE!0 = "consecutive_year"
  ToE&consecutive_year = consecutive_year
  copy_VarCoords(ub,ToE(0,:,:))
  printVarSummary(ToE)

  ;===================================================================

  print("Now is dealing with creating a new nc. file.")

  nclpath = systemfunc("pwd")
  nclpath = nclpath + "/" + systemfunc("ls *.ncl -tr |tail -n 1")
  ;...
  ;===================================================================

  nyear = dimsizes(year)
  filo = models_nc_names(current_nmodel)

  system("/bin/rm -f " + diro + filo)    ; remove if exists
  fout  = addfile (diro + filo, "c")  ; open output file

  ;===================================================================
  ; explicitly declare file definition mode. Improve efficienconsecutive_year.
  ;===================================================================
      setfileoption(fout,"DefineMode",True)
  ;===================================================================
  ; create global attributes of the file
  ;===================================================================
      fAtt               = True            ; assign file attributes
      fAtt@title         = "NCL Efficient Approach to netCDF Creation"
      fAtt@source_file   =  "original-file.nc"
      fAtt@Conventions   = "None"
      fAtt@nclpath       = nclpath
      fAtt@creation_date = systemfunc ("date")
      fileattdef( fout, fAtt )            ; copy file attributes
  ;===================================================================
      dimNames = (/"consecutive_year","year", "lat",  "lon"/)
      dimSizes = (/ nconsecutive_year,  nyear, nlat ,  nlon /)
      dimUnlim = (/False, False,  False, False/)
      filedimdef(fout,dimNames,dimSizes,dimUnlim)
  ;===================================================================
     filevardef(fout, "ub" ,typeof(ub),getvardims(ub))
     filevardef(fout, "ts" ,typeof(ts),getvardims(ts))
     filevardef(fout, "ToE" ,typeof(ToE),getvardims(ToE))
     filevardef(fout, "consecutive_year" ,typeof(consecutive_year),getvardims(consecutive_year))
     filevardef(fout, "year" ,typeof(year),getvardims(year))
     filevardef(fout, "lat"  ,typeof(lat),getvardims(lat))
     filevardef(fout, "lon"  ,typeof(lon),getvardims(lon))
  ;===================================================================
     filevarattdef(fout,"ub" ,ub)
     filevarattdef(fout,"ts" ,ts)
     filevarattdef(fout,"ToE" ,ToE)
     filevarattdef(fout,"consecutive_year",consecutive_year)
     filevarattdef(fout,"year" ,year)
     filevarattdef(fout,"lat"  ,lat)
     filevarattdef(fout,"lon"  ,lon)
  ;===================================================================
  ; explicitly exit file definition mode. **NOT REQUIRED**
  ;===================================================================
      setfileoption(fout,"DefineMode",False)
  ;====================================================================
     fout->consecutive_year   = (/consecutive_year/)
     fout->year   = (/year/)
     fout->lat    = (/lat/)
     fout->lon    = (/lon/)
     fout->ub     = (/ub/) ; (3)
     fout->ts     = (/ts/) ; (4)
     fout->ToE    = (/ToE/); (5)

  print("NCL file path: "+nclpath)
  print(".nc file path: "+diro+filo)

end
