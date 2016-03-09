function Sp=poreSurface(Dp,Lp,r)
%Surface of the cone
tanalpha=Dp./(2*Lp);
Sc=pi*Dp*Lp.*(tanalpha+sqrt(1+tanalpha.^2))/2;

%surface of the sphere cap
tanphi=Dp./(Lp+r)/2;
Ssc=2*pi*(Lp+r).^2.*(tanphi.^2-sqrt(tanphi.^2+1)+1);

Sp=Sc+Ssc;
