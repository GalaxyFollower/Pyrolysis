function F=shapeFactor(Dp,Lp,r1)
alpha=atan(Dp/(2*Lp));
phi=atan(Dp./(2*(r1+Lp)));

F=((Lp+r1)./cos(phi)+r1*sin(alpha)./(1-cos(phi)).*(log(sin(alpha-phi)/sin(alpha))*sin(alpha)...
    +phi.*cos(alpha)))./(sqrt((Lp).^2 +Dp^2/4));


end