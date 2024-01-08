undef("centroid_lon")
function centroid_lon(lat:numeric,lon:numeric,var:numeric)
local dims,rank,dphi,dtheta,dx,dy,i,lon2d
begin

;--- check if "var" is valid
  dims = dimsizes(var)
  rank = dimsizes(dims)
  if(rank.ne.2) 
    print("Error: The input data must be 2-dimensional")
  else if (dims(0).ne.dimsizes(lat))
    print("Error: The number of elements of 1st dimension of input data isn't identical to nlat.")
  else if (dims(1).ne.dimsizes(lon))
    print("Error: The number of elements of 1st dimension of input data isn't identical to nlon.")
  end if
  end if
  end if
;--- check complete

;--- calculate 1-d dphi and dtheta
;  nlat = dims(0)
;  nlon = dims(1)

  dphi   = new(dims(0),typeof(var))
  dtheta = new(dims(1),typeof(var))
  
  dphi(1:dims(0)-2) = ( lat(2:dims(0)-1) - lat(0:dims(0)-3) ) / 2.
  dphi(0) = lat(1) - lat(0)
  dphi(dims(0)-1) = lat(dims(0)-1) - lat(dims(0)-2)

  dtheta(1:dims(1)-2) = ( lon(2:dims(1)-1) - lon(0:dims(1)-3) ) / 2.
  dtheta(0) = lon(1) - lon(0)
  dtheta(dims(1)-1) = lon(dims(1)-1) - lon(dims(1)-2)

;--- calculation complete

;--- compute 2-d dx and dy
  dx = new(dims,typeof(var))
  dy = new(dims,typeof(var))

  dy = conform(var,dphi,0) ; consider that earth radius is a constant; thus dy = constant
  do i = 0,dims(0)-1
    dx(i,:) = dtheta * cos(lat(i)*4.*atan(1.0)/180.) ; dx is larger towards the equator
  end do
;--- computation complete

 
;--- calculate centroid of longitude (c_lon)
  lon2d = conform(var,lon,1)
  c_lon = sum (var*dx*dy*lon2d) / sum (var*dx*dy) ; centroid weighted by mass*area (var*dx*dy)
;--- calculation complete
  

  return(c_lon)
 end
