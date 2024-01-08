;###########################################################################
;---Goal: Calculate the eastward shift of the MJO
;---Input: MJO precipitation from 100 ensemble members from CESM2-LE
;---For details, please check our paper: https://doi.org/10.1175/JCLI-D-23-0213.1
;---Please send me an email if you have any questions: hien.bui@monash.edu
;###########################################################################


load "./centroid_lon.ncl"
load "./centroid_lat.ncl"

begin
  dir  = "./"
  output_dir = "./output/"
  input_var_name   = "PRECT"
  output_var_name1 = "c_lat"
  output_var_name2 = "c_lon"
  
;************************************************************
; READ DATA
;************************************************************
  latN = 8.
  latS = -12.
  lonW = 60
  lonE = 210

  file1  = systemfunc("ls "+dir+"*.nc")
  f      = addfile(file1(i),"r")

  var    = f->$input_var_name$({latS:latN},{lonW:lonE})
  lat    = tofloat(f->lat({latS:latN}))
  lon    = tofloat(f->lon({lonW:lonE}))
  nlat   = dimsizes(lat)
  nlon   = dimsizes(lon)
  printVarSummary(var)
  print("=====================================")
  print("Starting to calculate stddev ....")

;**********************************************
; climatological
;**********************************************
  c_lon = centroid_lon(lat,lon,var({latS:latN},{lonW:lonE}))
  c_lat = centroid_lat(lat,lon,var({latS:latN},{lonW:lonE}))
delete(var)
  print(c_lon)
  print(c_lat)

;**********************************************
; Output .nc files
;**********************************************
;---- Creating filename and output variable

  str = str_split(file1(i),"/")
  str:= str_split(str(dimsizes(str)-1),"._")
  output_filename = str_join(str(1:7),".")  ; extract the ensemble name
  print("Now is dealing with file: "+output_filename)

;-----
  filo       = "shift_" + output_filename + ".nc" ; entire output .nc filename

  system("/bin/rm -f " + output_dir + filo)    ; remove if exists
  fout  = addfile (output_dir + filo, "c")  ; open output file

;-----
  setfileoption(fout,"DefineMode",True)

  fAtt               = True            ; assign file attributes
  fAtt@title         = "NCL Efficient Approach to netCDF Creation"
  fAtt@source_file   =  "original-file.nc"
  fAtt@Conventions   = "None"
  fAtt@creation_date = systemfunc ("date")
  fileattdef( fout, fAtt )            ; copy file attributes
 ;===================================================================
       setfileoption(fout,"DefineMode",False)
       fout->$output_var_name1$      = (/c_lat/)
       fout->$output_var_name2$      = (/c_lon/)

  print(".nc file path: "+output_dir+filo)
end
