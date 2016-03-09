function Vp=poreVolume(Dp,Lp,r)
%volume of the cone
Vc=pi*Dp^2*Lp/12;
%Volumr of the sphere cap
tanphi=Dp./(Lp+r)/2;
Vsc=pi*(Lp+r).^3/3.*(2*(tanphi.^2+1).^1.5-3*tanphi.^2-2);
Vp=Vc+Vsc;
%Vp=Vc;

