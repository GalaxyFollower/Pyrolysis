% Returns the time derivatives for the component densities (?i),
% temperature (T), and the position of the centers of mass (CM) 
% for each component.
% INPUT:
%   t: time
%   X: vector containing, component densities and temperature vectors for 
%       each node.
%   r: vector containing the radial positions of the  
%   param: structure variable, containing the thermochemical, 
%       thermophysical and kinectic parameters of the system among other 
%       constants.
% OUTPUT:
%   dXdt: Column vector containing the time 
function [dXdt]=massAndEnergyBalancesOld(t,X,r,param)
t
n=length(r);

X=reshape(X,[length(r) param.nc+1]);
rho=X(:,1:param.nc);
T=round(X(:,param.nc+1)*10000)/10000;

gases=param.nsc+1:param.nc;
solids=1:param.nsc;


%------------------------------------------------------------
%Calculation of the pressure; its gradient and the velocity at the nodes
%------------------------------------------------------------
P=round((rho(:,gases)*(1./param.M(gases))).*T*param.R*100)/100;
% P(n+1)=param.Pinf;
% aP=polynomiame([r ;2*r(end)-r(n-1)],P,2);
% daPdr=derivatethem(aP);
% dPdr=evaluatethemsimply(daPdr,[r ;2*r(n)-r(n-1)]);

P(n)=param.Pinf;
aP=polynomiame(round(r*1e10)/1e5,P,2);
aP(:,1)=aP(:,1)/1e10;
aP(:,2)=aP(:,2)/1e5;
daPdr=derivatethem(aP);
dPdr=evaluatethemsimply(daPdr,r);


dPdr(1)=0;
P=P(1:n);
dPdr=dPdr(1:n);
daPdr=daPdr(1:n-2,:);
aP=aP(1:n-2,:);

upwind=dPdr<0;
if(t>0)
    [P T];
end
% P(n)=param.Pinf;
% aP=polynomiame(r,P,1);
% daPdr=derivatethem(aP);
% dPdr=evaluatethem(daPdr,r,upwind);

porosity=1-(1-param.porosity0)*sum(rho(:,solids),2)/param.rhob0;
K=porosity/param.porosity0.*(rho(:,solids)*param.K)./sum(rho(:,solids),2);
u=-K/param.mu.*dPdr;

%---------------------------------------------------------------
%calculation of the center of mass
%---------------------------------------------------------------

rhototal=sum(rho,2);
arhototal=polynomiame(r,rhototal,2);
cm=centerOfMass(arhototal,r,upwind);

%---------------------------------------------------------------
% Mass Balances
%---------------------------------------------------------------

uu=repmat(u,[1 2]);
rr=repmat(r,[1 2]);

massconvection(:,gases)=rr(2:n,:).^2.*uu(2:n,:).*rho(2:n,gases)...
    -rr(1:n-1,:).^2.*uu(1:n-1,:).*rho(1:n-1,gases);
[RxnRate,Hrxn]=reactionRate(r,rho,T,param,upwind); %[n,3]=size(RxnRate]
rr=repmat(r,[1 param.nc]);%auxiliar matrix for next operation
dRhodt=3*(-massconvection+RxnRate)./(rr(2:n,:).^3-rr(1:n-1,:).^3);

%-----------------------------------------------------------------
% Energy Balance
%-----------------------------------------------------------------
aT=polynomiame(round(r*1e10)/1e5,T,2);
aT(:,1)=aT(:,1)/1e10;
aT(:,2)=aT(:,2)/1e5;

daTdr=derivatethem(aT);
rhoCv=rho.*([T.^2 T ones([n 1])]*param.Cv');
asumrho34=polynomiame(r, (K/param.mu).*sum(rhoCv(:,3:4),2),2);
athermalconvection=[multipleConvolution(multipleConvolution(daPdr,daTdr),asumrho34)  zeros([n-2 2]) ];
thermalconvection=integratethem(athermalconvection,r,upwind);
aPdPdr=multipleConvolution([daPdr zeros([n-2 2])],aP);
akmu=polynomiame(r,K/param.mu,2);
PdPdr= volumemean(akmu,r,upwind).*integratethem(aPdPdr,r,upwind);
dTdr=evaluatethemsimply(daTdr,r);
dTdr(1)=0;
conduction=param.lambda*(r(2:n).^2.*dTdr(2:n) - r(1:n-1).^2.*dTdr(1:n-1));
asumrho=polynomiame(r, (K/param.mu).*sum(rhoCv,2),2);
sumrho=volumemean(asumrho,r,upwind);
dTdt=(thermalconvection+PdPdr+conduction -Hrxn)./sumrho;




%-----------------------------------------------------------------
% Extrapolate to nodes
%-----------------------------------------------------------------

%calculate dRhodt at the nodes
adRhodt1=polynomiame([- cm(1); cm],[dRhodt(1,1) ;dRhodt(:,1)],2);% simmetry around r=0
adRhodt2=polynomiame([- cm(1); cm],[dRhodt(1,2) ;dRhodt(:,2)],2);% simmetry around r=0
adRhodt3=polynomiame([- cm(1); cm],[dRhodt(1,3) ;dRhodt(:,3)],2);% simmetry around r=0
adRhodt4=polynomiame([- cm(1); cm],[dRhodt(1,4) ;dRhodt(:,4)],2);% simmetry around r=0
adTdt=polynomiame([- cm(1); cm],[dTdt(1) ;dTdt(:)],2);% simmetry around r=0

dRho1dt=evaluatethem([adRhodt1;adRhodt1(end,:)],r,upwind);
dRho2dt=evaluatethem([adRhodt2;adRhodt2(end,:)],r,upwind);
dRho3dt=evaluatethem([adRhodt3;adRhodt3(end,:)],r,upwind);
dRho4dt=evaluatethem([adRhodt4;adRhodt4(end,:)],r,upwind);
dTdt=evaluatethem([adTdt;adTdt(end,:)],r,upwind);
dTdt(end)=0;
dXdt=[dRho1dt;dRho2dt;dRho3dt;dRho4dt;dTdt];
end

function [RxnRate,Hrxn]=reactionRate(r,rho,T,param,upwind)
    %[1]:Biomass->gas
    %[2]:Biomass->tar
    %[3]:Biomass->char
    %[4]:tar->gas
    
    A=param.kinetics.A;
    E=param.kinetics.E;
    R=param.R;
    n=length(r);
    k=repmat(A,[n 1]).*exp(-(1./(R*T))*E);
    
    rxn(:,1)=k(:,1).*rho(:,1);
    rxn(:,2)=k(:,2).*rho(:,1);
    rxn(:,3)=k(:,3).*rho(:,1);
    rxn(:,4)=k(:,4).*rho(:,3);
    
    TT=[T.^3 T.^2 T];
    Cvrxn1=param.Cv(4,:)-param.Cv(1,:);
    Cvrxn2=param.Cv(3,:)-param.Cv(1,:);
    Cvrxn3=param.Cv(2,:)-param.Cv(1,:);
    Cvrxn4=param.Cv(4,:)-param.Cv(3,:);
   
    
    Hrxn=(TT*Cvrxn1').*rxn(:,1)+(TT*Cvrxn2').*rxn(:,2) ...
        +(TT*Cvrxn3').*rxn(:,3)+(TT*Cvrxn4').*rxn(:,4);
    
    Rxn(:,1)=volumeintegral(r,rxn(:,1),upwind);
    Rxn(:,2)=volumeintegral(r,rxn(:,2),upwind);
    Rxn(:,3)=volumeintegral(r,rxn(:,3),upwind);
    Rxn(:,4)=volumeintegral(r,rxn(:,1),upwind);
    Hrxn=volumeintegral(r,Hrxn,upwind);
    
    RxnRate(:,1)=-Rxn(:,1)-Rxn(:,2)-Rxn(:,3);
    RxnRate(:,2)=Rxn(:,3);   
    RxnRate(:,3)=Rxn(:,2)-Rxn(:,4);
    RxnRate(:,4)=Rxn(:,1)+Rxn(:,4);
end
