clc
clear
addpath('polynomialOperations')
addpath('AuxiliarScripts')
addpath('StarchStructure')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global param
%Initial temperature[K]
T0=298.15;
%Temperature at the border [K]
Tinf=2500;
%Initial pressure [Pa]
P0=1E5;
%Initial particle radious[m]ans
r=25E-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INTEGRATION AND DISCRETIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%number of discretization volumes
nv=12;
%integration time [s]
param.tf=0.5;
%Order of polynomials for the regression of density and temperature
order=3;
%Fit until the second derivative
deriv2=true;
%Changing order of magnitude of the radius 
param.unit=10^-((-2*(log10(r)<0)+1).*ceil(abs(log10(r))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DECLARATION OF PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Particule radius
param.r=r;
%Ideal gas constant [J/(mol*K)]:
param.Rgas=8.314;
%Molar mass [kg/mol]:
param.W=[0 ;0;128e-3;44e-3]; %[Biomass ???; tar; gas ]
%Initial Density of biomass [Kg/m^3]
param.rhob0=1492; %[Karathanos, 1993]
%Initial Pressure
param.P0=P0;
%Permeability of solids [m^2]
param.K=[5e-8; 1e-5]; %Biomass ; charcoal
%Viscosity of both gases [kg /(m*s)]
param.mu=3e-5;
%pressure at the infinity [Pa]
param.Pinf=1E5;
%number of reactions
param.nrxn=4;
%number of compounds
param.nc=4;
%number of solid compounds
param.nsc=2;
%number of gas compounds
param.ngc=2;
%Heat Capacity [J/(Kg*K)]
param.Cp=1000*[0 1e-3 1.5
          -6.7e-7 2e-3 0.44
          -2e-6 4.6e-3 -0.162
          -2e-7 7e-4 0.761 ];
param.Cv=param.Cp;
param.Cv(3:4,3)=param.Cp(3:4,3)-param.Rgas./param.W(3:4);
%Reference temperature
param.Tref=293.15;
%Thermal conductivity [W/(m*K)]
param.lambda=[0.35; 0.1; 0.26];
%Emisivity of biomass
param.eps=0.15;
%Temperature at the surface
param.Tinf=Tinf;
%Initial Temperature 
param.T0=T0;

%--------Biomass Structure---------------
%Volume of pores
param.Vp=1.1E-6;%m3/kg [Juszczak,2002]
%Surface of pores
param.Sp=475; %m2/kg[Juszczak,2002]
param.Dp=sqrt(3)*9.28*1E-9;%m
%Porosity of biomass
param.porosity0=param.rhob0*param.Vp;
%param.Lp=getPoreLength(param);
%Length of conical part of starch pore
param.Lp=104.5E-9;
%pore diameter in char
param.dc=0.34E-9;%m
%--------Reaction Kinetics---------------
%[1]:Biomass->gas
%[2]:Biomass->tar
%[3]:Biomass->char
%[4]:tar->gas

%Frequency factor [1/s]:

param.kinetics.A=[6.8e8 8.23e8 0 1e5];%torrado
param.kinetics.A=[1.3e8 2e8 1.08e7 0];%Chan et al
param.kinetics.A=[1.3e8 2e8 1.08e7 1.48e6];%Chan et al2
param.kinetics.name='chan et al2';
param.kinetics.A=[ 0 0 0 0];
%Activation Energy [J/mol]:
param.kinetics.E=[139e3 119e3 0 85e3];%torrado
param.kinetics.E=[140e3 133e3 121e3 144e3];%chan et al
%Enthalpies of reaction [J/kg]
param.kinetics.Hrxn=[25e3 25e3 25e3 0 ] ;%[Haseli,2011]
param.kinetics.Hrxn=[538e3 0 -2E6 0 538e3] ;%[Milosavljevic; 1996]
%param.kinetics.Hrxn=[418e3 418e3 418e3 0 ];%[Authier,2009]


%--------Discretization parameters---------------
%Polynomial order
param.order=order;
%Fit until the second derivative
param.deriv2=deriv2;
%Number of volumes
param.nv=nv;
%number of polynomials
if(deriv2)
    param.np=(nv-1)/(order-2)+2;
else
    param.np=(nv)/(order);
end


%-----------External gas properties------------
%                (Those of air)
param.air.W=28.97E-3;%Kg/mol 
param.air.Cp=[-1.291418E-7 3.230959E-4 -4.703434E-2 9.917485E2];
%Lennard-Jones parameters for viscosity;
param.air.sigma=3.617;%A°
param.air.e_k=97;%K
% Multiplier for the terminal velocity of the particle
param.vter_mult=1;
%Stefan boltzman constant
param.stefBoltz=5.670373E-8; %W m^2 K^4
%emisivity
param.emis=1;


ri=(0:r/(nv):r)'*1E6;
% ri=[(0:4*r/(3*nv):r/3)'
%     (r/3+2*r/(3*nv):2*r/(3*nv):2*r/3)'
%     (2*r/3+4*r/(3*nv):4*r/(3*nv):r)']*1E6;

%  ri=[(0:3*r/(2*nv):r/2)'
%       (r/2+3*r/(4*nv):3*r/(4*nv):r)']*1E6;
n=length(ri);
CM0=0.75*(ri(2:n).^4-ri(1:n-1).^4)./( ri(2:n).^3-ri(1:n-1).^3);
%ri=[ri(1:n-1);ri(n-1)+0.75*(ri(n)-ri(n-1));ri(end)];
%CM0=mean([ri(2:nv+1) ri(1:nv)],2);
T= ones([nv+1,1])*T0;
%T=fixLast(CM0,T, Tinf, ri);

rhog0=P0*param.porosity0*param.W(4)./(param.Rgas*T);

%%
X0=[ones([nv+1 1])*param.rhob0
    zeros([nv+1 1]) 
%    zeros([nv+1 1]) 
%    rhog0(1:nv+1)
    ones([nv+1 1])
    ones([nv,1])*P0
    T(1:nv)
    CM0*1E-6
    ];


%matlabpool open local 4
figure(3)
for i=1:4
    hdl(i)=subplot(2,2,i);
end
figure(4)
for i=1:4
    hdl(i+4)=subplot(2,2,i);
end
figure(5)
for i=1:4
    hdl(i+8)=subplot(2,2,i);
end

tic
plotem=true;

odeoptions=odeset('Stats','on','OutputFCn', ...
    @(t,X,init)plotDataAtStep(t,ri*1E-6,X,param,init,hdl,plotem),'NonNegative',1);
[t X]=ode15s(@(t,X)massAndEnergyBalances(t,X,ri*1e-6,param),[ 0 param.tf],X0,odeoptions);
param.duration=toc;


%matlabpool close
%[X t]=rungekutta(@(t,X)massAndEnergyBalances(t,X,ri*1e-6,param),[0 tf],X0,0.00001);
yield=evaluateData(t,ri*1E-6,X,param,true);

