;###########################################################################
;---Author: Yi-Xian Li, modified by Hien Bui
;---Goal: Make animation for stationary waves and wildfire
;---For details, please check our paper: https://doi.org/10.1029/2022GL099017
;---Please send me an email if you have any questions: hien.bui@monash.edu
;###########################################################################

begin

  ntimes   = 81
  dir = "./data/"
;--
  
  inf_fire1  = addfile(dir+"LE2-1231.011_NFIRE_1850-2100_ANN.nc","r")
  inf_fire2  = addfile(dir+"LE2-1231.011_burn_1850-2100_ANN.nc","r")
  inf_wave1  = addfile(dir+"LE2-1231.011_burn_1850-2100_ANN.nc","r")
  inf_wave2  = addfile(dir+"LE2-1231.011_burn_1850-2100_ANN.nc","r")

  fire1  = inf_fire1->NFIRE(170:251,:,:)
  fire1  = fire1*86400*365*1000 
  printVarSummary(fire1)
  mf1    = inf_fire1->NFIRE(150:169,:,:)
  mf1    = mf1*86400*365*1000
  mean1  = dim_avg_n_Wrap(mf1,0)
  mean11 = conform_dims(dimsizes(fire1),mean1,(/1,2/))
  printVarSummary(mean11)
  fire2  = inf_fire2->FAREA_BURNED(170:251,:,:)
  fire2  = fire2*86400*365*1000
  mf2    = inf_fire2->FAREA_BURNED(150:169,:,:)
  mf2    = mf2*86400*365*1000
  mean2  = dim_avg_n(mf2,0)
  mean22 = conform_dims(dimsizes(fire2),mean2,(/1,2/))
  fire1 = where(fire1.le.0,fire1@_FillValue,fire1)
  fire2 = where(fire2.le.0,fire2@_FillValue,fire2)
  wave1  = inf_wave1->FAREA_BURNED(170:251,:,:)
  wave2  = inf_wave1->FAREA_BURNED(170:251,:,:)
  ;wave1 = wave1*86400*365*1000
  ;wave2 = wave2*86400*365*1000
  fire1_ano = fire1
  fire2_ano = fire2
  fire1_ano = (fire1-mean11)
  fire2_ano = (fire2-mean22)

  fire1_ano = fire1_ano>0
  fire2_ano = fire2_ano>0


  time   = inf_fire1->time 

;-- set Max/Min/Int
  fireMin = 0.
  fireMax = 60
  fireInt = 10
  npts    = ((fireMax-fireMin)/fireInt)+1

  bMin = 0.
  bMax = 300
  bInt = 50

  waveMin = 0.0
  waveMax = 50.0
  waveInt = 10.0
  nconts  = ((waveMax-waveMin)/waveInt)+1

;-- set plot output file name
  plotout = "./burn_ano"

;-- set workstation resources
  wks_type                   = "x11"                  ;-- plot output type
  wks_type@wkBackgroundColor = "black"                ;-- set workstation background to black (or grey18)
  wks_type@wkWidth           =  2500                  ;-- for presentations 16:9 ratio
  wks_type@wkHeight          =  2500                
  wks  = gsn_open_wks(wks_type,plotout)
  cmap = "MPL_OrRd" 
  gsn_define_colormap(wks,cmap)                     ;-- define colormap
  colors = gsn_retrieve_colormap(wks)               ;-- retrieve color map for common labelbar

;-- set resources
  res                      =  True                  ;-- set resources for plot
  res@cnInfoLabelOn        =  False                 ; Yi-Xian helps you turn it off :)
  res@gsnDraw              =  False                 ;-- don't draw plot yet
  res@gsnFrame             =  False                 ;-- don't advance frame
  res@gsnAddCyclic         =  True                  ;-- add cyclic point
  res@gsnRightString       =  ""                    ;-- don't draw right string
  res@gsnLeftString        =  ""                    ;-- don't draw left string


  res@lbLabelBarOn         =  False                 ;-- add common colorbar
  


  res@tiMainFontHeightF    =  0.02                  ;-- main title font size
  res@tiMainOffsetYF       =  0.06                  ;-- move title upward
  res@tiMainFontColor      = "white"                ;-- set font color to white
  
  res@vpYF                 =  0.78                  ;-- set viewport Y position
  res@vpHeightF            =  0.45                  ;-- set viewport width
  res@vpWidthF             =  0.8                   ;-- set viewport height


  res3                     = res
  res4                     = res

  res@cnFillOn             =  True                  ;-- enable color fill
  res@cnLinesOn            =  False                 ;-- disable contour lines
  res@cnLineLabelsOn       =  False                 ;-- disable contour lines
  res@cnFillMode           =  "RasterFill"           ;-- set fill mode
  res@cnRasterSmoothingOn  =  True                  ;-- smooth contouring
  res@cnLevelSelectionMode = "ManualLevels"         ;-- use manual contour line levels
  res@cnMinLevelValF       =  fireMin                ;-- contour min. value
  res@cnMaxLevelValF       =  fireMax                ;-- contour max. value
  res@cnLevelSpacingF      =  fireInt                ;-- contour interval

  res@mpOceanFillColor     = "White"
  res@mpInlandWaterFillColor = "White"
  res@mpProjection         = "Orthographic"         ;-- use projection
  res@mpCenterLatF         =  35                    ;-- center latitude at 53 deg north
  res@mpPerimOn            =  False                 ;-- don't draw the box around the plot

  res1                     =  res                   ;-- resource setting for left plot
  res1@vpXF                =  -0.1                  ;-- set viewport X position 1. panel

  res2                     =  res1                   ;-- resource setting for right plot
  res2@vpXF                =  0.365                 ;-- set viewport X position 2. panel
  res2@cnMinLevelValF      =  bMin                ;-- contour min. value
  res2@cnMaxLevelValF      =  bMax                ;-- contour max. value
  res2@cnLevelSpacingF     =  bInt
;-- resources for the contours
  res3@cnFillOn            = False
  res3@cnLinesOn           = True
  res3@cnMonoLineColor     = False
  res3@cnMinLevelValF      = waveMin             
  res3@cnMaxLevelValF	   = waveMax           
  res3@cnLevelSpacingF     = waveInt               
  res3@cnLineThicknessF    = 5.0
  res3@cnLevelSpacingF     = 1
  res3@cnLineLabelsOn      = False

  res4@cnFillOn            = False
  res4@cnLinesOn           = True
  res4@cnMinLevelValF      = waveMin
  res4@cnMaxLevelValF      = waveMax
  res4@cnLevelSpacingF     = waveInt
  res4@cnLineThicknessF    = 5.0
  res4@cnLevelSpacingF     = 1
  res4@cnLineLabelsOn      = False 

;-- text resources: labelbar units
  tires1                   =  True                  ;-- text resources title string
  tires1@txFontColor       = "white"                ;-- change to white
  tires1@txFontHeightF     =  0.018                 ;-- text font size
  tires1@txJust            = "CenterCenter"         ;-- text justification

  tires2                   =  True                  ;-- text resources year string
  tires2@txFontColor       = "white"                ;-- change to white
  tires2@txFontHeightF     =  0.015                 ;-- text font size
  tires2@txJust            = "CenterCenter"         ;-- text justification
  tires2@txFont            =  30                    ;-- text font "courier-bold"

  tires3                   =  True                  ;-- text resources copyright string
  tires3@txFontColor       = "white"                ;-- change to white
  tires3@txFontHeightF     =  0.013
  tires3@txJust            = "BottomRight"          ;-- text justification
  
;-- create plots
  incr     = 10;0.5 ;-- turn ?? degree to the east for each time step
  startlon = 15.0                                 
  tstart   = 0                                 

;-- assign plot array for the data sets displayed on two globes
  plotf = new(2,graphic)
  plotw = new(2,graphic)

;-- loop over time steps and rotate the globe 1 degrees per frame
  ;j=0
  do i=tstart,ntimes-1
     ;irot = where((j*incr).le.-180., (j*incr)+360., (j*incr)) ;-- rotation longitude step
     irot = where((i*incr).le.-180., (i*incr)+360., (i*incr)) ;-- rotation longitude step
     rlon = startlon - irot                                   ;-- define mpCenterLonF value
     res1@mpCenterLonF = rlon                                 ;-- center lon position
     res2@mpCenterLonF = rlon   

     plotf(0) = gsn_csm_contour_map(wks,fire1_ano(i,:,:),res1)
;     plotw(0) = gsn_csm_contour(wks,wave1(i,:,:),res3)
;     plotw(0) = ColorNegDashZeroPosContour(plotw(0),"blue","blue","blue")
;     overlay(plotf(0),plotw(0))
     plotf(1) = gsn_csm_contour_map(wks,fire2_ano(i,:,:),res2)
;     plotw(1) = gsn_csm_contour(wks,wave2(i,:,:),res4)
;     plotw(1) = ColorNegDashZeroPosContour(plotw(1),"blue","blue","blue")
;     overlay(plotf(1),plotw(1))
 
     draw(plotf)

;####################
; Add label bar
;####################
undef("add_labelbar1")
procedure add_labelbar1(wks,map)
begin
;-- retrive values from plot
  getvalues map@contour
     "cnFillColors"              :   cols
     "cnLevels"                  :   levs
     "vpWidthF"                  :   vpw
     "vpHeightF"                 :   vph
     "vpXF"                      :   vpx
     "vpYF"                      :   vpy
  end getvalues
         
;-- set the labels and the labelbar width, height and position
  labs     = levs+""
  nlevs    = dimsizes(labs)
  lbwidth  = 0.5
  lbheight = 0.06
  lbx      = ((1.0-(lbwidth+vpx))/2) + vpx/2
  lby      = ((1.0-(vpy-vph))/2) - 0.01

;-- set labelbar resources
  lbres                       =  True
  lbres@gsnFrame              =  False                 ;-- don't advance frame
  lbres@vpWidthF              =  lbwidth               ;-- width of labelbar
  lbres@vpHeightF             =  lbheight              ;-- height of labelbar
  lbres@vpXF                  =  lbx                   ;-- labelbar x-position
  lbres@vpYF                  =  lby                   ;-- labelbar y-position
  lbres@lbPerimOn             =  False                 ;-- no label bar box
  lbres@lbOrientation         =  "Horizontal"          ;-- orientation
  lbres@lbLabelFontHeightF    =  0.015                 ;-- label font height
  lbres@lbLabelAlignment      =  "InteriorEdges"       ;-- where to label
  lbres@lbLabelPosition       = "Top"                  ;-- write labels on top of the labelbar
  lbres@lbMonoFillPattern     =  True                  ;-- fill solid
  lbres@lbFillColors          =  cols                  ;-- use colors
  lbres@lbFillPattern         =  0                     ;-- fill solid
  lbres@lbLabelFontColor      = "white"                ;-- label color
  lbres@lbMonoFillColor       =  False                 ;-- no mono fill color
  
  gsn_labelbar_ndc(wks,nlevs+1,labs,lbx,lby,lbres)     ;-- draw transparent labelbar
  
end

    add_labelbar1(wks,plotf)
   
     setvalues plotf@data
       "sfDataArray" : (/inf_fire1->NFIRE(i,:,:)*86400*365*1000/)
       "sfDataArray" : (/inf_fire2->FAREA_BURNED(i,:,:)*86400*365*1000/)
     end setvalues
     setvalues plotw@data
       "sfDataArray" : (/inf_fire1->NFIRE(i,:,:)*86400*365*1000/)
       "sfDataArray" : (/inf_fire2->FAREA_BURNED(i,:,:)*86400*365*1000/)
     end setvalues

     gsn_text_ndc(wks,"           LE2-1231.011~C~                ",0.5,0.82,tires1) 
     gsn_text_ndc(wks,"Fire count",0.1, 0.8, tires1)           
     gsn_text_ndc(wks,"Burned area", 0.9, 0.8, tires1)        
     ;gsn_text_ndc(wks,"[fraction/year]", 0.9, 0.301, tires3)      
     gsn_text_ndc(wks,"Year:"+(i+2020),0.525, 0.348, tires2) 

     frame(wks)
     
;     j = j+1
  end do

end
