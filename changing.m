function changing(Tinf, kinetics,r,vter_mult)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global param
%Initial temperature[K]
T0=298.15; 
%Temperature at the border [K]
%Tinf=1000;
%Initial pressure [Pa]
P0=1E5;
%Initial particle radious[m]
%r=2.5E-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INTEGRATION AND DISCRETIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%number of discretization volumes

nv=ceil(8*r/35e-6+60/7);
%nv=100;
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
param.rhob0=1440; %[Karathanos, 1993]
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
param.lambda=[0.35; 0.1; 0.026];

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
param.kinetics=kinetics;


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
param.vter_mult=vter_mult;
%Stefan boltzman constant
param.stefBoltz=5.670373E-8; %W m^2 K^4
%emisivity
param.emis=0.95;

ri=(0:r/(nv):r)';
% ri=[(0:4*r/(3*nv):r/3)'
%     (r/3+2*r/(3*nv):2*r/(3*nv):2*r/3)'
%     (2*r/3+4*r/(3*nv):4*r/(3*nv):r)']*1E6;

%  ri=[(0:3*r/(2*nv):r/2)'
%       (r/2+3*r/(4*nv):3*r/(4*nv):r)']*1E6;
n=length(ri);
CM0=0.75*(ri(2:n).^4-ri(1:n-1).^4)./( ri(2:n).^3-ri(1:n-1).^3);
%ri=[ri(1:n-1);ri(n-1)+0.75*(ri(n)-ri(n-1));ri(end)];
%CM0=mean([ri(2:nv+1) ri(1:nv)],2);
T= ones([nv,1])*T0;
%T=fixLast(CM0,T, Tinf, ri);

rhog0=param.porosity0*P0*param.W(4)./(param.Rgas*T);

%
X0=[ones([nv+1 1])*param.rhob0
    zeros([nv+1 1]) 
    %zeros([nv 1]) 
    %rhog0(1:nv)
    ones(nv+1,1)
    ones([nv,1])*(P0-param.Pinf)
    T
    ];



figure(getFigureHdl('Densities'));
for i=1:4
    hdl(i)=subplot(2,2,i);
end
figure(getFigureHdl('Densities Gradients'));
for i=1:4
    hdl(i+4)=subplot(2,2,i);
end
figure(getFigureHdl('T & P'));
for i=1:4
    hdl(i+8)=subplot(2,2,i);
end
figure(getFigureHdl('CM'));
for i=1:2
    hdl(i+12)=subplot(2,1,i);
end

tic
plotem=true;

odeoptions=odeset('Stats','on','OutputFCn', ...
    @(t,X,init)plotDataAtStep(t,ri,X,param,init,hdl,plotem),'NonNegative',1);


time=0;
XX=X0';
t0=0;
tf=param.tf;
succes=true;
%while mean(X0(1:nv+1))>0.01*param.rhob0
    %try
     
        [t X]=ode15s(@(t,X)massAndEnergyBalances(t,X,ri,param),[ t0 tf],X0,odeoptions);

        time=[time;t(2:end)];
        XX=[XX;X(2:end,:)];
        if(t(end)<tf)
            succes=false;
        end
        t0=tf;
        tf=1.5*tf;
        X0=X(end,:)';
        
%     catch err
%         ERROR=['ERROR at: ' kinetics.name '_Tinf=' num2str(param.Tinf) '_r=' num2str(param.r) ...
%     '_nvter=' num2str(param.vter_mult) '\n'];
%         fprintf(ERROR)
%         fileID = fopen('Simulations\Errors.txt','w');
%         fid = fopen('changing.txt','w');
%         fprintf(fid, ERROR);
%         fclose(fid);
%         succes=false;
%     end
    if ~succes
        X0(:)=0;
    end
%end

i=find(mean(XX(:,1:nv+1),2)<.01*param.rhob0);
if isempty(i)
    i=length(time);
end

time=time(1:i(1));
XX=XX(1:i(1),:);

param.tf=time(end);
param.duration=toc;
save(['Simulations\' kinetics.name '_Tinf=' num2str(param.Tinf) '_r=' num2str(param.r) ...
    '_nvter=' num2str(param.vter_mult) '.mat']);





