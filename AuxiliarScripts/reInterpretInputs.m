



n=length(r_ev);
nv=param.nv;%number of volumes
nc=param.nc;%number of components
nsc=param.nsc;%number of solid components
ngc=nc-nsc;%number of gaseous components
nrxn=param.nrxn;%number of reacitons
order=param.order;%order of polynomials
np=param.np;%number of polynomials
Rgas=param.Rgas;%Ideal gas constant
W=param.W;%molas masses
unit=param.unit;

%Component indices
bio=1; %biomass
cha=2;
gou=3;%tar=goudron (tar is a matlab function)
gas=4;
gases=param.nsc+1:param.nc;
solids=1:param.nsc;
%nodes indices
bionodes=1:nv+1;
charnodes=max(bionodes)+1:max(bionodes)+nv+1;
solidnodes=[bionodes charnodes];
%gounodes=max(charnodes)+1:max(charnodes)+nv;
%gasnodes=max(gounodes)+1:max(gounodes)+nv;
ygasnode=max(charnodes)+1:max(charnodes)+nv+1;
%gasesnodes=[gounodes gasnodes ];
Pnodes=max(ygasnode)+1:max(ygasnode)+nv;
Tnodes=max(Pnodes)+1:max(Pnodes)+nv;
CMnodes=max(Tnodes)+1:max(Tnodes)+nv;
CMTnodes=max(CMnodes)+1:max(CMnodes)+nv;

%% ------------------------------------------------------------
% Reinterpretation of inputs
%------------------------------------------------------------

rho=zeros(nv+1,nc);
rho(:,bio)=X(bionodes);
rho(:,cha)=X(charnodes);
porosity=1-(1-param.porosity0)*sum(rho(:,solids),2)/param.rhob0;
%rho(1:nv,gou)=X(gounodes);
%rho(1:nv,gas)=X(gasnodes);
ygas=X(ygasnode);
%ygas=rho(end,gas)/sum(rho(end,gases));

P=X(Pnodes);

T=X(Tnodes);

rhotot=(P+param.Pinf).*porosity(1:nv)/Rgas./T./(ygas(1:nv)/W(gas)+(1-ygas(1:nv))/W(gou));
rho(1:nv,gou)=(1-ygas(1:nv)).*rhotot;
rho(1:nv,gas)=ygas(1:nv).*rhotot;


CM=0.75*(r_ev(2:nv+1).^4-r_ev(1:nv).^4)./( r_ev(2:nv+1).^3-r_ev(1:nv).^3);
%CM=X(CMnodes);
r=[r_ev(1:nv);CM(end);r_ev(nv+1)];


%CM=mean([r(1:n-1) r(2:n)],2);
try
    cord=[CM; r(end)];
    [Tend rhogasesend]=limitTemperature(CM,T,rho, ygas(end),r,param);
catch err
    CM=0.75*(r_ev(2:nv+1).^4-r_ev(1:nv).^4)./( r_ev(2:nv+1).^3-r_ev(1:nv).^3);
    r=[r_ev(1:nv);CM(end);r_ev(nv+1)];
    [Tend rhogasesend]=limitTemperature(CM,T,rho, ygas(end),r,param);
end
    
T=[T; Tend];
P=[P;0]; 
%rhotot=P.*porosity/Rgas./T./(ygas/W(gas)+(1-ygas)/W(gou));
rho(end,gases)=rhogasesend;


%% Fixing external pressure and temperature
% if(t<0.1)
%     Tend=(param.Tinf-param.T0)/0.1*t+param.T0;
% else
%     Tend=param.Tinf;
% end
% Tend=param.Tinf;
% T=[T; Tend];
% rhotot=(1/(ygas/W(gas) +(1-ygas)/W(gou))*param.Pinf/Rgas/T(end));
% rho(nv+1,gou)=(1-ygas)*rhotot;
% rho(nv+1,gas)=ygas*rhotot;
%P=(Rgas*(rho(:,gases)*(1./W(gases))).*T);

pressure=[rho(1:nv,gases)*(1./W(gases))*Rgas.*T(1:nv)./porosity(1:nv);param.Pinf];
