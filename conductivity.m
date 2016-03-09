
function lambda=conductivity(rho,T,porosity,cord,param)
Lp=param.Lp;
Dp=param.Dp;
dc=param.dc;
emis=param.emis;
stefBoltz=param.stefBoltz;
gamma=shapeFactor(Dp,Lp,cord);
porosity0=param.porosity0;
phi=atan(Dp./(2*(cord+Lp)));
lambda=(rho(:,1)./sum(rho(:,1:2),2)*param.lambda(1)+ ...
       rho(:,2)./sum(rho(:,1:2),2)*param.lambda(2)).*(1-porosity)+ ...
       param.lambda(3)*porosity;
   %4*emis*stefBoltz.*T.^3.*(gamma.*(cord.*(1-cos(phi))+Lp).*porosity0...
   %   +dc*(porosity-porosity0));
end
