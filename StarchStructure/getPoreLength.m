function Lp=getPoreLength(param)
unit=1E9;
rho=param.rhob0*1000/unit^3;
Sp=param.Sp*unit^2/1000;
Dp=param.Dp*unit;
R=45E-6*unit;% mesh size of the particules analised by Juszczak (2002)
eps=param.porosity0;
Sp=Sp*rho*4*pi*(R)^3/3;
options = optimoptions('fsolve','tolfun',1E-10);
Lp=fsolve(@(Lp)abs(Sp- ...
    surfaceOfPores(parameps,Sp,Dp,Lp,R))/Sp,104,options)/unit;