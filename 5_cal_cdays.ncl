;###########################################################################
;---Goal: Calculate the consecutive high-pressue days
;---Input: omega at 500 hPa and KBDI from CESM2-LE
;---For details, please check our paper: https://doi.org/10.1029/2022GL099017
;---Please send me an email if you have any questions: hien.bui@monash.edu
;###########################################################################


begin

; --- data path and variables
dir_w = "./"
dir_f = "./"
v_w   = "w500" ; variable name in the nc file
v_f   = "kbdi" ; variable name in the nc file

; --- reading data

  file1 = systemfunc("ls "+dir_w+"*.nc")
  f_w   = addfile(file1(jj),"r")

  file2 = systemfunc("ls "+dir_f+"*.nc")
  f_f   = addfile(file2(jj),"r")

  read_dim =  getfilevardimsizes(f_w,v_w) ; suppose that variables w and f have the same dimensions
  time_dim = read_dim(0) -1

  v_lat = "lat"
  v_lon = "lon"

; --- domain
  time_step= 1
  lat_n    = 90.
  lat_s    = 0.
  lon_w    = 0
  lon_e    = 360

; --- convert w500 to [0,1] array -- positive=1, zero and negative =0
  wpos = where(f_w->$v_w$(::time_step,{lat_s:lat_n},{lon_w:lon_e}).gt.0,1,0)
  fire = f_f->$v_f$(::time_step,{lat_s:lat_n},{lon_w:lon_e})
  clim = dim_avg_n_Wrap(fire,0)
  clim = where(clim.ne.0,clim,clim@_FillValue)

  lat = f_w->$v_lat$({lat_s:lat_n})
  lon = f_w->$v_lon$({lon_w:lon_e})
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)
  
; --- convert wpos from [0,1] to [0,consecutive-days] array

  ;wtemp = wpos ; turn on if you want to compare [0,1] and [0,c-days] arrays
  print("Converting wpos starts")

  do i = 0, nlat-1
    do j = 0, nlon-1
      x = ind(wpos(:,i,j).eq.0) ; locations where wpos=0 ; i.e., negative w500

;==============================================================================
      if(.not.any(ismissing(x))) ; prevent Missing (i.e., all(wpos)=1) occurs
        k = 0 ; in case dimsizes(x)=0 and thus not replaced

;     --head
        if(x(0).ne.0)
          wpos(0:x(0)-1,i,j) = x(0)
        end if

;     --between
        if(dimsizes(x).ge.2)
          do k = 0,dimsizes(x)-2 ; -2 for NCL's 0-start and two adjacent elements required at a time for one loop
            wpos(x(k)+1:x(k+1)-1,i,j) = x(k+1)-x(k)-1
          end do
        end if

;     --tail
        if(x(dimsizes(x)-1).ne.time_dim)
          wpos(x(dimsizes(x)-1)+1:,i,j) = time_dim-x(dimsizes(x)-1)
        end if

      end if
;==============================================================================

      delete(x) ; dimensions of x can vary
    end do
   ; print("lat: "+lat(i))
  end do

  print("Converting wpos finished")

;---------
; put the comment mark below at the beginning if you want to compare [0,1] and [0,c-days] arrays by printing them out
/;
  do i = 0, nlat-1
    do j = 0, nlon-1
      print(wpos(:,i,j)+", "+wtemp(:,i,j))
      print("---------------------------")
    end do
  end do
;/

;----------
; calculate n-day averages for each n

  n_max = 20;max(wpos) ; the longest consecutive days
  n_day = ispan(0,n_max,1)
  n_day!0 = "n_day"
  nday = dimsizes(n_day)

  print("Creating nd_avg array")
  nd_avg = new((/n_max+1,nlat,nlon/),float)
  print("nd_avg array created")

  do n = 0, n_max ; 0-day means the collection of negative-500w events
    fire_temp = where(wpos.eq.n,fire,fire@_FillValue)
    nd_avg(n,:,:) = dim_avg_n(fire_temp,0)/clim  ; n-day average
    nd_avg(n,:,:) = where(ismissing(nd_avg(n,:,:)).and..not.ismissing(fire(0,:,:)),0,nd_avg(n,:,:)) ; add 0 to Missing area on lands

  ;  print("n-day: "+n)
  end do


  print("calculation of nd_avg array finished")

;----------
; feed Metadata (n-day,lat,lon)

  nd_avg!0 = "n_day"
  nd_avg&n_day = n_day
  nd_avg!1 = "lat"
  nd_avg&lat = lat
  nd_avg!2 = "lon"
  nd_avg&lon = lon
  printVarSummary(nd_avg)

;===============================================
 print("=================================================================")
 print("Writing out NetCDF files ....")
;---- Creating filename and output variable

  str = str_split(file2(jj),"/")
  str:= str_split(str(dimsizes(str)-1),"._")
  output_filename = str_join(str(5:6),".")  ; extract the ensemble name
  print("Now is dealing with file: "+output_filename)

;-----
  output_dir = "./"
  filo       = "cdays_kbdi_his_" + output_filename + ".nc" ; entire output .nc filename

  system("/bin/rm -f " + output_dir + filo)    ; remove if exists
  fout  = addfile (output_dir + filo, "c")  ; open output file

    ;===================================================================
    ; explicitly declare file definition mode. Improve efficiency.
    ;===================================================================
        setfileoption(fout,"DefineMode",True)
    ;===================================================================
    ; create global attributes of the file
    ;===================================================================
        fAtt               = True            ; assign file attributes
        fAtt@title         = "NCL Efficient Approach to netCDF Creation"
        fAtt@source_file   =  "original-file.nc"
        fAtt@Conventions   = "None"
        fAtt@creation_date = systemfunc ("date")
        fileattdef( fout, fAtt )            ; copy file attributes
    ;===================================================================
        dimNames = (/"n_day", "lat" , "lon"/)
        dimSizes = (/ nday  ,  nlat , nlon /)
        dimUnlim = (/False  , False , False/)
        filedimdef(fout,dimNames,dimSizes,dimUnlim)
    ;===================================================================
       filevardef(fout, "n_day" ,typeof(n_day),getvardims(n_day))
       filevardef(fout, "lat"  ,typeof(lat),getvardims(lat))
       filevardef(fout, "lon"  ,typeof(lon),getvardims(lon))
       filevardef(fout, "nd_avg"    ,typeof(nd_avg)  ,getvardims(nd_avg))
    ;===================================================================
       filevarattdef(fout,"nd_avg",nd_avg)
       filevarattdef(fout,"n_day" ,n_day)
       filevarattdef(fout,"lat"  ,lat)
       filevarattdef(fout,"lon"  ,lon)
    ;===================================================================
    ; explicitly exit file definition mode. **NOT REQUIRED**
    ;===================================================================
        setfileoption(fout,"DefineMode",False)
    ;====================================================================
       fout->n_day   = (/n_day/)
       fout->lat     = (/lat/)
       fout->lon     = (/lon/)
       fout->nd_avg  = (/nd_avg/)

    print(".nc file path: "+output_dir+filo)

end
