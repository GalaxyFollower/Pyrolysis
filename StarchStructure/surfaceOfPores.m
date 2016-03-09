function [Sp_calc, np]=surfaceOfPores(eps,Sp,Dp,Lp,R)

%Volume of sphere
Vs=4*pi*R^3/3;
%mean r1
r=(R-Lp)/2;
vp=poreVolume(Dp,Lp,r);
np=eps*Vs./vp;
Vp=np.*poreVolume(Dp,Lp,r);
Sp_calc=np.*poreSurface(Dp,Lp,r);







