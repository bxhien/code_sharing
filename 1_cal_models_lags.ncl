;###########################################################################
;---Goal: Calculate linear regression models to estimate feedback between
;soil moisture and next-day precipitation in Australia
;---Input: CYGNSS soil moisture and GPM precipitation
;---For details, please check our paper: https://doi.org/10.1088/1748-9326/ad15b7
;---Please send me an email if you have any questions: hien.bui@monash.edu
;###########################################################################


begin

  dir   = "./"
  fil_p = "GPM_prec_day_2018-2022.nc4"          ;--- precipitation
  fil_s = "CYGNSS_sm_day_2018-2022.nc"          ;--- soil moisture
  fx    = addfile(dir+fil_s,"r")                
  fy    = addfile(dir+fil_p,"r")                
  filo  = "model"+model+"_multivar_Yest_4-lag.nc"

  varx_name = "SM_daily"
  vary_name = "precipitationCal"
  
;--- find out overlapped time
  time_x = fx->time - 0.5 			;--- shift phase
  time_y = dble2flt(fy->time)

  same_time = venn2_intersection(time_x,time_y)
  timeind   = get1Dindex(time_y,same_time)

;--- create xlat/lon as not existent in the nc file
  sz   = getfilevardimsizes(fx,varx_name)
  latN = 38
  latS = -38
  lonW = -135
  lonE = 164
  xlat = fspan(latS,latN,sz(1))
  xlon = fspan(lonW,lonE,sz(2))
  
;---  subset domain
  latN = -10
  latS = -38
  lonW = 110
  lonE = 160

  totalind = ind(xlat.ge.latS.and.xlat.le.latN)
  xlat_sta = totalind(0)
  xlat_end = totalind(dimsizes(totalind)-1) 
  totalind:= ind(xlon.ge.lonW.and.xlon.le.lonE)
  xlon_sta = totalind(0)
  xlon_end = totalind(dimsizes(totalind)-1)
  xlat    := xlat(xlat_sta:xlat_end)
  xlon    := xlon(xlon_sta:xlon_end)

  x        = fx->$varx_name$(:,xlat_sta:xlat_end,xlon_sta:xlon_end)
  y        = fy->$vary_name$(timeind,{lonW:lonE},{latS:latN})
  y       := y (time|:,lat|:,lon|:) ; re-arrange dimensions

  ylat     = fy->lat({latS:latN})
  ylon     = fy->lon({lonW:lonE})

;--- regrid
  y:= linint2_Wrap(ylon,ylat,y,False,xlon,xlat,0)

;--- lead-lag shift
  y!1  = "lat"
  y!2  = "lon"
  lat  = y&lat
  lon  = y&lon

  avgx = dim_avg_n_Wrap(x,0)
  avgx:= conform(y,avgx,(/1,2/))
  y    = (/where(ismissing(avgx),y@_FillValue,y)/) 		;--- mask out ocean

  day_lags = 4 							;--- include -1 to -4 days
  ntime    = dimsizes(y&time) - day_lags 			;--- length used for every time series

  nlat = dimsizes(lat)
  nlon = dimsizes(lon)

  if (model.eq.1) 						;--- without soil
     xx = new((/ntime,day_lags,nlat,nlon/),float) 		;--- depedant variables
  else if (model.eq.2) 						;--- with soil (yesterday)
     xx = new((/ntime,day_lags+1,nlat,nlon/),float)
  end if
  end if

  ismiss = new((/ntime,nlat,nlon/),float)
  ismiss = 0. 							;--- later to detect where is missing

  i = 0
  do while (i.lt.day_lags) 					;--- from day -4 to day -1 (i=0 -> day=-4; i=3 -> day=-1)
    xx(:,i,:,:) = y(i:ntime-1+i,:,:)
    ismiss = ismiss + xx(:,i,:,:)
    i = i + 1
  end do
  iszero = ismiss 

;--- pp = previous precipitation
  pp1   = xx(:,0,:,:)  ;--- ntime,nlat,nlon
  ppavg = dim_avg_n_Wrap(xx,1)

  printVarSummary(xx)

  y:= y(day_lags:,:,:) ;--- today (day=0)
  time = y&time

  if (model.eq.2) 						;--- with soil (yesterday)
     i                  = day_lags -1
     x                 := x(i:ntime-1+i,:,:) 
     x&time            := y&time
     x@long_name        = "Yesterday soil moisture"
     copy_VarCoords(y,x)
     xx(:,day_lags,:,:) = (/x/) 				;--- i=4 ->day=-1 soil
     ismiss = ismiss + xx(:,day_lags,:,:)
     x_std  = conform(x,dim_stddev_n_Wrap(x,0),(/1,2/))
     ismiss = where(x_std.eq.0,ismiss@_FillValue,ismiss) 
  end if

;--- multi-variable regression
  sz          = dimsizes(y)
  Yest        = new(sz,float)
  Yest_nzprec = new(sz,float)
  Yest_zprec  = new(sz,float)
  RSE         = new((/nlat,nlon/),float) 			;--- residual standard error, similar to RMSE but considering degree of freedom
  RSE_nzprec  = new((/nlat,nlon/),float)
  RSE_zprec   = new((/nlat,nlon/),float)

  opt = True

  nvars = day_lags 						;--- how many days = how many dependant vars
  if (model.eq.2) nvars = day_lags + 1 end if 			;--- plus yesterday soil 

  nsample        = new((/nlat,nlon/),float) 			;--- number of time points used
  nsample_nzprec = new((/nlat,nlon/),float) 
  nsample_zprec  = new((/nlat,nlon/),float) 

  cor            = new((/nlat,nlon/),float)
  cor_nzprec     = new((/nlat,nlon/),float)
  cor_zprec      = new((/nlat,nlon/),float)

  reg            = new((/nvars+1,nlat,nlon/),float)
  reg_nzprec     = new((/nvars+1,nlat,nlon/),float)
  reg_zprec      = new((/nvars+1,nlat,nlon/),float)

  pval           = new((/nvars+1,nlat,nlon/),float) 
  pval_nzprec    = new((/nvars+1,nlat,nlon/),float) 
  pval_zprec     = new((/nvars+1,nlat,nlon/),float) 

  do i = 0,nlat-1
     do j = 0,nlon-1
        vvind:= ind(.not.ismissing(ismiss(:,i,j)))
        if all(ismissing(vvind)) continue end if
        if dimsizes(vvind).le.nvars+1 continue end if 		;--- causes error (and meaningless) when too small
        b = reg_multlin_stats(y(vvind,i,j),xx(vvind,:,i,j),opt) ;--- Yest = b(0) + b(1)*X1 + b(2)*X2 + ...+b(M)*XM
        cor(i,j)        = b@r
        reg(:,i,j)      = (/b/)
        nsample(i,j)    = b@N
        Yest(vvind,i,j) = b@Yest
        RSE(i,j)        = b@RSE
        pval(:,i,j)     = b@pval

;--- conditional regression when any of the days' (-4 to -1) prec is not zero (i.e., acummulated rainfall >0 in the past 4 days) for model_1
;--- conditional regression when any of the days' (-4 to -1) prec and yesterday's soil is not zero for model_2
       nzind:= ind(ismiss(:,i,j).gt.0) 				;--- non-zero-precs locations
       if all(ismissing(nzind)) continue end if
       if dimsizes(nzind).le.nvars+1 continue end if 		;--- causes error (and meaningless) when too small
       b = reg_multlin_stats(y(nzind,i,j),xx(nzind,:,i,j),opt)
       cor_nzprec(i,j)        = b@r
       reg_nzprec(:,i,j)      = (/b/)
       nsample_nzprec(i,j)    = b@N
       Yest_nzprec(nzind,i,j) = b@Yest
       RSE_nzprec(i,j)        = b@RSE
       pval_nzprec(:,i,j)     = b@pval

;--- conditional regression when all of the days' (-4 to -1) prec are zero (i.e., accumulated rainfall = 0 in the past 4 days) and yesterday soil has value (non-missing)
       zind:= ind(.not.ismissing(ismiss(:,i,j)).and.iszero(:,i,j).eq.0) 	;--- zero-precs locations
       if all(ismissing(zind)) continue end if
       if dimsizes(zind).le.nvars+1 continue end if 		;--- causes error (and meaningless) when too small
       b = reg_multlin_stats(y(zind,i,j),xx(zind,:,i,j),opt)
       cor_zprec(i,j)       = b@r
       reg_zprec(:,i,j)     = (/b/)
       nsample_zprec(i,j)   = b@N
       Yest_zprec(zind,i,j) = b@Yest
       RSE_zprec(i,j)       = b@RSE
       pval_zprec(:,i,j)    = b@pval
    end do
  end do

  copy_VarCoords(y(0,:,:),RSE)
  copy_VarCoords(y(0,:,:),RSE_nzprec)
  copy_VarCoords(y(0,:,:),RSE_zprec)
  copy_VarCoords(y(0,:,:),cor)
  copy_VarCoords(y(0,:,:),cor_nzprec)
  copy_VarCoords(y(0,:,:),cor_zprec)
  copy_VarCoords(y(0,:,:),nsample)
  copy_VarCoords(y(0,:,:),nsample_nzprec)
  copy_VarCoords(y(0,:,:),nsample_zprec)
  copy_VarCoords(y,Yest)
  copy_VarCoords(y,Yest_nzprec)
  copy_VarCoords(y,Yest_zprec)
  
  R = Yest - y
  copy_VarCoords(y,R)

  lag_day = ispan(-(day_lags+1),-1,1)
  lag_day(0) = 0 						;--- lag_day = (/0,-4,-3,-2,-1/) where lag_day= 0 means y-intercept
  if (model.eq.2) 
    lag_day:= array_append_record(lag_day,-1,0) 		;--- lag_day = (/0,-4,-3,-2,-1,-1/) where the last -1 means yestady soil
  end if
  lag_day!0 = "lag_day"
  reg!0 = "lag_day"
  reg&lag_day = lag_day

  copy_VarCoords(y(0,:,:),reg(0,:,:))
  copy_VarCoords(reg,reg_nzprec)
  copy_VarCoords(reg,reg_zprec)
  copy_VarCoords(reg,pval)
  copy_VarCoords(reg,pval_nzprec)
  copy_VarCoords(reg,pval_zprec)
  printVarSummary(cor)
  printVarSummary(reg)
  printVarSummary(reg_nzprec)

  nlag_day       = dimsizes(lag_day)
  R@long_name    = "Residual (Yest - y)"
  y@long_name    = "Today precipitation"
  pval@long_name = "p-value"

;--- output to nc file
  diro = "./"
  system("/bin/rm -f " + diro + filo) 		   ;--- remove if exists
  fout  = addfile (diro + filo, "c")  		   ;--- open output file

  ;===================================================================
      setfileoption(fout,"DefineMode",True)
  ;===================================================================
      fAtt               = True            ; assign file attributes
      fAtt@title         = "NCL Efficient Approach to netCDF Creation"
      fAtt@source_file   =  "original-file.nc"
      fAtt@Conventions   = "None"
      fAtt@creation_date = systemfunc ("date")
      fileattdef( fout, fAtt )            ; copy file attributes
  ;===================================================================
      dimNames = (/"time", "lag_day",  "lat", "lon"/)
      dimSizes = (/ntime, nlag_day , nlat ,  nlon /)
      dimUnlim = (/True, False, False, False/)
      filedimdef(fout,dimNames,dimSizes,dimUnlim)
  ;===================================================================
      filevardef(fout, "time"  ,typeof(time),getvardims(time))
      filevardef(fout, "lag_day"  ,typeof(lag_day),getvardims(lag_day))
      filevardef(fout, "lat",typeof(lat),getvardims(lat))
      filevardef(fout, "lon",typeof(lon),getvardims(lon))
      filevardef(fout, "cor"  ,typeof(cor),getvardims(cor))
      filevardef(fout, "cor_nzprec"  ,typeof(cor_nzprec),getvardims(cor_nzprec))
      filevardef(fout, "cor_zprec"  ,typeof(cor_zprec),getvardims(cor_zprec))
      filevardef(fout, "reg",typeof(reg),getvardims(reg))
      filevardef(fout, "reg_nzprec"  ,typeof(reg_nzprec),getvardims(reg_nzprec))
      filevardef(fout, "reg_zprec"  ,typeof(reg_zprec),getvardims(reg_zprec))
      filevardef(fout, "pval",typeof(pval),getvardims(pval))
      filevardef(fout, "pval_nzprec"  ,typeof(pval_nzprec),getvardims(pval_nzprec))
      filevardef(fout, "pval_zprec"  ,typeof(pval_zprec),getvardims(pval_zprec))
      filevardef(fout, "nsample"  ,typeof(nsample),getvardims(nsample))
      filevardef(fout, "nsample_nzprec"  ,typeof(nsample_nzprec),getvardims(nsample_nzprec))
      filevardef(fout, "nsample_zprec"  ,typeof(nsample_zprec),getvardims(nsample_zprec))
      filevardef(fout, "RSE",typeof(RSE),getvardims(RSE))
      filevardef(fout, "RSE_nzprec"  ,typeof(RSE_nzprec),getvardims(RSE_nzprec))
      filevardef(fout, "RSE_zprec"  ,typeof(RSE_zprec),getvardims(RSE_zprec))
      filevardef(fout, "Yest",typeof(Yest),getvardims(Yest))
      filevardef(fout, "Yest_nzprec"  ,typeof(Yest_nzprec),getvardims(Yest_nzprec))
      filevardef(fout, "Yest_zprec"  ,typeof(Yest_zprec),getvardims(Yest_zprec))
      filevardef(fout, "y",typeof(y),getvardims(y))
      filevardef(fout, "R",typeof(R),getvardims(R))
      if (model.eq.2)
         filevardef(fout, "x",typeof(x),getvardims(x))
      end if
      filevardef(fout, "pp1"  ,typeof(pp1),getvardims(pp1))
      filevardef(fout, "ppavg"  ,typeof(ppavg),getvardims(ppavg))
  ;===================================================================
      filevarattdef(fout,"nsample"       ,nsample)
      filevarattdef(fout,"nsample_nzprec",nsample_nzprec)
      filevarattdef(fout,"nsample_zprec" ,nsample_zprec)
      filevarattdef(fout,"reg",reg)
      filevarattdef(fout,"reg_nzprec" ,reg_nzprec)
      filevarattdef(fout,"reg_zprec"  ,reg_zprec)
      filevarattdef(fout,"pval",pval)
      filevarattdef(fout,"pval_nzprec",pval_nzprec)
      filevarattdef(fout,"pval_zprec" ,pval_zprec)
      filevarattdef(fout,"Yest",Yest)
      filevarattdef(fout,"Yest_nzprec",Yest_nzprec)
      filevarattdef(fout,"Yest_zprec" ,Yest_zprec)
      filevarattdef(fout,"RSE",RSE)
      filevarattdef(fout,"RSE_nzprec" ,RSE_nzprec)
      filevarattdef(fout,"RSE_zprec"  ,RSE_zprec)
      filevarattdef(fout,"cor",cor)
      filevarattdef(fout,"cor_nzprec" ,cor_nzprec)
      filevarattdef(fout,"cor_zprec"  ,cor_zprec)
      filevarattdef(fout,"lat",lat)
      filevarattdef(fout,"lon",lon)
      filevarattdef(fout,"time",time)
      filevarattdef(fout,"lag_day",lag_day)
      filevarattdef(fout,"R",R)
      if (model.eq.2)
         filevarattdef(fout,"x",x)
      end if
      filevarattdef(fout,"pp1",pp1)
      filevarattdef(fout,"ppavg",ppavg)
  ;===================================================================
      setfileoption(fout,"DefineMode",False)
  ;====================================================================
      fout->lat            = (/lat/)
      fout->lon            = (/lon/)
      fout->time           = (/time/)
      fout->lag_day        = (/lag_day/)
      fout->y              = (/y/)
      fout->R              = (/R/)
      fout->reg            = (/reg/)
      fout->reg_nzprec     = (/reg_nzprec/)
      fout->reg_zprec      = (/reg_zprec/)
      fout->pval           = (/pval/)
      fout->pval_nzprec    = (/pval_nzprec/)
      fout->pval_zprec     = (/pval_zprec/)
      fout->Yest           = (/Yest/)
      fout->Yest_nzprec    = (/Yest_nzprec/)
      fout->Yest_zprec     = (/Yest_zprec/)
      fout->RSE            = (/RSE/)
      fout->RSE_nzprec     = (/RSE_nzprec/)
      fout->RSE_zprec      = (/RSE_zprec/)
      fout->nsample        = (/nsample/)
      fout->nsample_nzprec = (/nsample_nzprec/)
      fout->nsample_zprec  = (/nsample_zprec/)
      fout->cor            = (/cor/)
      fout->cor_nzprec     = (/cor_nzprec/)
      fout->cor_zprec      = (/cor_zprec/)
      if (model.eq.2)
         fout->x = (/x/)
      end if
      fout->pp1            = (/pp1/)
      fout->ppavg          = (/ppavg/)

  print(".nc file path: "+diro+filo)

end
