%...................................
%polynomials for rho_i
%...................................
a_rho=zeros([param.np order+1 nc]);
a_drhodt=zeros([param.np order nc]);

for i=1:nc
    [a_rho(:,:,i), intervals] =getPolynomials(CM,rho(:,i),r,unit,order);
    a_drhodr(:,:,i)=differentiatePolynomials(a_rho(:,:,i));
    drhodr_cm=evaluatePolynomials(a_drhodr(:,:,i),intervals,cord);
    [a_drhodr_cm(:,:,i), intervals] =getPolynomials(CM,drhodr_cm,r,unit,order);  
end

a_drhodr=a_drhodr_cm;


rho_r=fun(a_rho,intervals,r_ev,1:nc);
drhodr_r=fun(a_drhodr,intervals,r_ev,1:nc);


%...................................
%Polynomial for T
%...................................


[a_T, intervals]=getPolynomials(CM,T,r,unit,order);
a_dTdr=differentiatePolynomials(a_T);
a_d2Tdr2=differentiatePolynomials(a_dTdr);
dTdr_r=evaluatePolynomials(a_dTdr,intervals,r_ev);


%...................................
%Polynomial for P
%...................................
[a_P, intervals] =getPolynomials(CM,P,r,unit,order);
P_r=evaluatePolynomials(a_P,intervals,r_ev)+param.Pinf;
a_dPdr=differentiatePolynomials(a_P);
a_d2Pdr2=differentiatePolynomials(a_dPdr);
dPdr_r=evaluatePolynomials(a_dPdr,intervals,r_ev);
d2Pdr2_r=evaluatePolynomials(a_d2Pdr2,intervals,r_ev);
%...................................
%Permeability
%...................................

a_K=(param.K(bio)-param.K(cha))*a_rho(:,:,bio)/param.rhob0+ ...
    [zeros([np order]) ones([np 1])*param.K(cha)];
a_dKdr=(param.K(bio)-param.K(cha))*a_drhodr(:,:,bio)/param.rhob0;
K_r=evaluatePolynomials(a_K,intervals,r_ev);

dKdr_r=evaluatePolynomials(a_dKdr,intervals,r_ev);

%...................................
% Heat capacity
%...................................
a_rhoCv=zeros([param.np 2*order+1 nc]);
for i=1:nc
    Cv(:,i)=polyval(param.Cv(i,:),T);
    [a_Cv(:,:,i) ] =getPolynomials(CM,...
        Cv(:,i),r,unit,order);
    [a_rhoCv(:,:,i) ] =multipleConvolution(a_rho(:,:,i),a_Cv(:,:,i));
end

%int_rw^re[rhoCv*r^2]dr
RHOCVR2=volumeIntegral(sum(a_rhoCv,3),intervals,r);
RHOCVR2=[RHOCVR2(1:nv-1);sum(RHOCVR2(nv:nv+1))];
    
%int_rw^re[rhoCv*r^3]dr
a_rhocvr=[sum(a_rhoCv,3) zeros(np,1)];
RHOCVR3=volumeIntegral(a_rhocvr,intervals,r);
RHOCVR3=[RHOCVR3(1:nv-1);sum(RHOCVR3(nv:nv+1))];

CMT=RHOCVR3./RHOCVR2;
%...................................
% velocity
%...................................
u_r=-K_r.*dPdr_r/param.mu;

%...................................
%Rates of reaction
%...................................
[a_rxn, a_Rrxn] = rxnRate(a_rho,T,r,CM,param);


%...................................
%  Conduction
%...................................
lambda=conductivity(rho,T,porosity,cord,param);
[a_lambda, intervals] =getPolynomials(CM,lambda,r,unit,order);

lambda_r=evaluatePolynomials(a_lambda,intervals,r_ev);
%...................................
%  Energy of formation
%...................................
deltaUrxn=Urxn(T,param);
a_Urxn=zeros([np order+1 nrxn]);
for i=1:nrxn
   a_Urxn(:,:,i)=getPolynomials(CM,deltaUrxn(:,i),r,unit,order);
end
