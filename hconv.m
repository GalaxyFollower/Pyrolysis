
function [h,Re]=hconv(T,param,r)

P=param.Pinf;
Rgas=param.Rgas;
W=param.air.W;%Kg/mol
e_k=param.air.e_k;
sigma=param.air.sigma;
if nargin<=2
    r=param.r;
else
    %r=r(:);
end
k=1.3806488E-23;%Boltzman constant
Na=6.02214129E23;%Avogadro's number
g=9.8;%m/s^2 gravity
n=param.vter_mult;% Multiplier for the terminal velocity of the particle;

%density
rhog=W*P/Rgas./T;

%viscosity
kT_e=T/e_k;
A=1.16145;
B=0.14874;
C=0.52487;
D=0.7732;
E=2.16178;
F=2.43787;
Omega=A*kT_e.^-B+C*exp(-D*kT_e)+E*exp(-F*kT_e);
C_mu=5/16*sqrt(k/Na/pi)*1E20;%1E20 because sigma is in amstrongs

mu=C_mu*sqrt(W*T)/sigma^2./Omega;%Pa*s;

%Conductivity[1]
%   [1] Handbook of Thermal Conductivity of Liquids and Gases;
%       Natan B. Vargaftik
%       url: goo.gl/uXWIVe
lambda=(-0.9474+11.961*(T/100)-2.3632*(T/100).^2+0.8406*(T/100).^3 ...
       -0.1747*(T/100).^4+1.904E-2*(T/100).^5-1.035E-3*(T/100).^6 ...
       +2.228E-5*(T/100).^7)*1E-3;%W/(m*K)

%terminal velocity
rhos=param.rho_mean;
v=n.*sqrt(4*g*r.^2.*(rhos-rhog)./(18.*mu));

% %Heat Capacity
Cp=polyval(param.air.Cp,T);

%Prandlt Number
Pr=Cp.*mu./lambda;
%Reynolds Number
Re=2*rhog.*v.*r./mu;
%Nusselt Number
Nu=2+0.03.*Pr.^0.33.*Re.^0.54+0.35.*Pr.^0.36.*Re.^0.58;
h=lambda.*Nu./(2*r);%W/(m^2*K)

