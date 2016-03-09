function [Tend, rhogasesend]=limitTemperature(CM,T,rho, ygas,r,param)


Tend0=mean([T(1); param.Tinf]);
options = optimoptions(@fsolve,'Display','off');
warning('off','all')
Tend=fsolve(@(Tend)energyBalance(Tend,rho,ygas,CM,T,r,param),Tend0,options);


    
warning('on','all')
[~, rho]=energyBalance(Tend,rho,ygas,CM,T,r,param);
rhogasesend=rho(end,3:4);

end

function [fo, rho]=energyBalance(Tend,rho,ygas,CM,T,r,param)
T=[T;Tend];
nv=param.nv;
W=param.W;
porosity=1-(1-param.porosity0)*sum(rho(:,1:2),2)/param.rhob0;
rhotot=param.Pinf*porosity(end)/param.Rgas/Tend/(ygas/W(4) +(1-ygas)/W(3));
rho(nv+1,3)=(1-ygas)*rhotot;
rho(nv+1,4)=ygas*rhotot;

%...................................
%polynomials for rho_i
%...................................
a_rho=zeros([param.np param.order+1 param.nc]);
for i=1:param.nc
    [a_rho(:,:,i), intervals] =...
        getPolynomials(CM,rho(:,i),r,param.unit,param.order);
end

%mean density of the whole particle
param.rho_mean=3*sum(volumeIntegral(sum(a_rho,3),intervals,intervals))/param.r^3;
%convection coefficient
param.h_conv=hconv(param.Tinf,param);


[a_T, intervalsT]=getPolynomials(CM,T,r,param.unit,param.order);
a_T=a_T(end,:);
a_dTdx=differentiatePolynomials(a_T);
lambda=conductivity(rho,T,porosity,CM,param);
h=param.h_conv;
Tinf=param.Tinf;
sigma=param.stefBoltz;
eps=param.emis;
fo=h*(Tinf-Tend)+sigma*eps*(Tinf^4-Tend^4)-lambda(end)*polyval(a_dTdx,param.r);
end

