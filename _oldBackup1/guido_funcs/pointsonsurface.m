function [newpoints,center,radius,coeffs]=pointsonsurface(vc,points,order);


[center,radius]=sphfit(vc)
[coeffs,err]=harmfit(vc,center,order);
newpoints=mk_vcharm(points,center,coeffs);
  
  
  return;