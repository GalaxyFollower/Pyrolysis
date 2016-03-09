function [dRhodt]=massbalance(t,r,rho,T,param)
t

r=r(:);
T=T(:);


n=length(r);
rho=reshape(rho,[length(r) param.nc]);
gases=param.nsc+1:param.nc;
solids=1:param.nsc;


%------------------------------------------------------------
%Calculation of the pressure and its gradient and the velocity at the nodes
%------------------------------------------------------------
P=(rho(:,gases)*(1./param.M(gases))).*T*param.R;
P(n+1)=param.Pinf;
a=polynomiame([r ;r(end)+1e-6],P);
daPdr=derivatethem(a);
dPdr=evaluatethemsimply(daPdr,[r ;2*r(n)-r(n-1)]);
% P(n)=param.Pinf;
% a=polynomiame(r,P);
% daPdr=derivatethem(a);
% dPdr=evaluatethemsimply(daPdr,r);

dPdr(1)=0;
dPdr=dPdr(1:n);
upwind=dPdr<0;
% for i=1:param.nsc
%     arhosolid=polynomiame(r,rho(:,i));
%     [rhoSolid(:,1)]=volumemean(arhosolid,r,upwind);
% end
porosity=1-(1-param.porosity0)*sum(rho(:,solids),2)/param.rhob0;
K=porosity/param.porosity0.*(rho(:,solids)*param.K)./sum(rho(:,solids),2);
%K=rho(:,1)/param.rhob0*param.K;
u=-K/param.mu.*dPdr;

%---------------------------------------------------------------
%calculation of the center of mass
%---------------------------------------------------------------

rhototal=sum(rho,2);
arhototal=polynomiame(r,rhototal);
cm=centerOfMass(arhototal,r,upwind);

uu=repmat(u,[1 2]);
rr=repmat(r,[1 2]);

massconvection(:,gases)=rr(2:n,:).^2.*uu(2:n,:).*rho(2:n,gases)...
    -rr(1:n-1,:).^2.*uu(1:n-1,:).*rho(1:n-1,gases);
[RxnRate]=reactionRate(r,rho,T,param.kinetics,upwind); %[n,3]=size(RxnRate]
rr=repmat(r,[1 param.nc]);%auxiliar matrix for next operation
dRhodt=3*(-massconvection+RxnRate)./(rr(2:n,:).^3-rr(1:n-1,:).^3);

%calculate dRhodt at the nodes
adRhodt1=polynomiame([- cm(1); cm],[dRhodt(1,1) ;dRhodt(:,1)]);% simmetry around r=0
adRhodt2=polynomiame([- cm(1); cm],[dRhodt(1,2) ;dRhodt(:,2)]);% simmetry around r=0
adRhodt3=polynomiame([- cm(1); cm],[dRhodt(1,3) ;dRhodt(:,3)]);% simmetry around r=0
adRhodt4=polynomiame([- cm(1); cm],[dRhodt(1,4) ;dRhodt(:,4)]);% simmetry around r=0


dRho1dt=evaluatethem([adRhodt1;adRhodt1(end,:)],r,upwind);
dRho2dt=evaluatethem([adRhodt2;adRhodt2(end,:)],r,upwind);
dRho3dt=evaluatethem([adRhodt3;adRhodt3(end,:)],r,upwind);
dRho4dt=evaluatethem([adRhodt4;adRhodt4(end,:)],r,upwind);


dRhodt=[dRho1dt;dRho2dt;dRho3dt;dRho4dt];
end

function [RxnRate]=reactionRate(r,rho,T,kinetics,upwind)
    %[1]:Biomass->gas
    %[2]:Biomass->tar
    %[3]:Biomass->char
    %[4]:tar->gas
    
    A=kinetics.A;
    E=kinetics.E;
    R=kinetics.R;
    n=length(r);
    k=repmat(A,[n 1]).*exp(-repmat(E,[n 1])./(R*T));
    %k(:,1)=A(1).*exp(-E(1)./(R*T));
    %k(:,2)=A(2).*exp(-E(2)./(R*T));
    %k(:,3)=A(3).*exp(-E(3)./(R*T));
      
    rxn(:,1)=volumeintegral(r,k(:,1).*rho(:,1),upwind);
    rxn(:,2)=volumeintegral(r,k(:,2).*rho(:,1),upwind);
    rxn(:,3)=volumeintegral(r,k(:,3).*rho(:,1),upwind);
    rxn(:,4)=volumeintegral(r,k(:,4).*rho(:,3),upwind);
    
    RxnRate(:,1)=-rxn(:,1)-rxn(:,2)-rxn(:,3);
    RxnRate(:,2)=rxn(:,3);   
    RxnRate(:,3)=rxn(:,2)-rxn(:,4);
    RxnRate(:,4)=rxn(:,1)+rxn(:,4);
end

