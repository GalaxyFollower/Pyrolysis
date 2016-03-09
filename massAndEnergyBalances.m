% Returns the time derivatives for the component densities (?i),
% temperature (T), and the position of the centers of mass (CM) 
% for each component.
% INPUT:
%   t: time
%   X: vector containing, component densities and temperature vectors for 
%       each node.
%   r_ev: vector containing the radial positions of the  
%   param: structure variable, containing the thermochemical, 
%       thermophysical and kinectic parameters of the system among other 
%       constants.
% OUTPUT:
%   dXdt: Column vector containing the time 
function [dXdt]=massAndEnergyBalances(logt,X,r_ev,param)

% ------------------------------------------------------------
% Reinterpretation of inputs
%------------------------------------------------------------
reInterpretInputs

% ------------------------------------------------------------
% Principal Polynomials
%------------------------------------------------------------
calculateAllPolynomials



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               MASS BALANCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rho*u*r^2(r)
%zero for all
rhour2_r=zeros(size(rho_r));
%Not zero for gases
rhour2_r(:,gases)=rho_r(:,gases).*repmat(u_r.*r_ev.^2,[1 2]); 



for i=1:nc
	RXNaux(:,:,i)=volumeIntegral(a_Rrxn(:,:,i),intervals,r);    
end




RXN_RATES=[RXNaux(1:nv-1,:);sum(RXNaux(nv:nv+1,:),1)]; 
% RXN_RATES(:,bio)=-sum(RXN(:,1:3),2);
% RXN_RATES(:,cha)=RXN(:,3);
% RXN_RATES(:,gou)=RXN(:,2)./porosity(1:nv)-RXN(:,4);
% RXN_RATES(:,gas)=RXN(:,1)./porosity(1:nv)+RXN(:,4);
%RXN(:,gases)=RXN(:,gases)./repmat(porosity(1:nv),[1 ngc]);

drhodt=zeros(size(rho));
drhodt(1:nv,:)=3*(rhour2_r(1:nv,:)-rhour2_r(2:nv+1,:) ...
    +RXN_RATES)./repmat(r_ev(2:nv+1).^3-r_ev(1:nv).^3,[1 nc]);


%Evaluation of the change of mass at the frontier

%drhodt(end,solids)=fun(a_Rrxn,intervals,r_ev(end),solids);
drhodt(end,bio)=polyval(a_Rrxn(end,:,bio),r_ev(end));
drhodt(end,cha)=polyval(a_Rrxn(end,:,cha),r_ev(end));

%drhodt(end,gou)=polyval(a_rxn(end,:,2)/porosity(end)-a_rxn(end,:,4),r_ev(end))...
%drhodt(end,gou)=polyval(a_Rrxn(end,:,gou),r_ev(end))...
%drhodt(end,gou)=...
drhodt(end,gou)=polyval(a_Rrxn(end,:,gou),r_ev(end))...
        +(2*K_r(end)*dPdr_r(end)*rho(end,gou)/r_ev(end) ...
        +dKdr_r(end)*dPdr_r(end)*rho(end,gou)...
        +K_r(end)*d2Pdr2_r(end)*rho(end,gou)...
        +K_r(end)*dPdr_r(end)*drhodr_r(end,gou))/param.mu;

%drhodt(end,gas)=polyval(a_rxn(end,:,1)/porosity(end)+a_rxn(end,:,4),r_ev(end))...
%drhodt(end,gas)=polyval(a_Rrxn(end,:,gas),r_ev(end))...
%drhodt(end,gas)=...
drhodt(end,gas)=polyval(a_Rrxn(end,:,gas),r_ev(end))...
        +(2*K_r(end)*dPdr_r(end)*rho(end,gas)/r_ev(end) ...
        +dKdr_r(end)*dPdr_r(end)*rho(end,gas)...
        +K_r(end)*d2Pdr2_r(end)*rho(end,gas)...
        +K_r(end)*dPdr_r(end)*drhodr_r(end,gas))/param.mu;
    
dydt=(drhodt(end,gas)*(sum(rho(end,gases))) ...
     -sum(drhodt(end,gases))*rho(end,gas))/sum(rho(end,gases))^2;
dydt=(drhodt(:,gas).*sum(rho(:,gases),2) ...
     -sum(drhodt(:,gases),2).*rho(:,gas))./sum(rho(:,gases),2).^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               ENERGY BALANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------
% % Convection term
% % ------------------------------------------------------------
%thermal part

a_dPdr_dTdr=multipleConvolution(a_dPdr ,a_dTdr);
a_rhoCvgases_dPdr_dTdr=multipleConvolution(a_dPdr_dTdr, ...
    sum(a_rhoCv(:,:,gases),3) );
a_rhoCvgases_dPdr_dTdr=multipleConvolution(a_K, a_rhoCvgases_dPdr_dTdr)...
                                /param.mu;
RHOCVdPdrdTdr=volumeIntegral(a_rhoCvgases_dPdr_dTdr,intervals,r);
RHOCVdPdrdTdr=[RHOCVdPdrdTdr(1:nv-1);sum(RHOCVdPdrdTdr(nv:nv+1))];

%Pressure part
a_dPdr2=multipleConvolution(a_dPdr,a_dPdr);
a_KdPdr2=multipleConvolution(a_K,a_dPdr2);
KDPDR2=volumeIntegral(a_KdPdr2,intervals,r);
KDPDR2=[KDPDR2(1:nv-1);sum(KDPDR2,1)];
P_DIVP=(r_ev(2:nv+1).^2.*K_r(2:nv+1).*P_r(2:nv+1).*dPdr_r(2:nv+1) ...
       -r_ev(1:nv).^2.*K_r(1:nv).*P_r(1:nv).*dPdr_r(1:nv) ...
    -KDPDR2)/param.mu;

% % ------------------------------------------------------------
% % Conduction term
% % ------------------------------------------------------------

DivLambdadTdr=r_ev.^2.*lambda_r.*dTdr_r;
% 
% % ------------------------------------------------------------
% % Energy of reaction term
% % ------------------------------------------------------------
a_Urxn_rxnrate=zeros([np (size(a_rxn,2)-1)+(size(a_Urxn,2)-1)+1 nrxn]);
for i=1:nrxn
    a_Urxn_rxnrate(:,:,i)=multipleConvolution(a_Urxn(:,:,i),a_rxn(:,:,i));
end
SUM_URXN_RXNRATE=volumeIntegral(sum(a_Urxn_rxnrate,3),intervals,r);
SUM_URXN_RXNRATE=[SUM_URXN_RXNRATE(1:nv-1); sum(SUM_URXN_RXNRATE(nv:nv+1))];

% % ------------------------------------------------------------
% % Heat capacity of the volume
% % ------------------------------------------------------------


%SUM_RHOCVaux= volumeIntegral(sum(a_rhoCv,3),intervals,r);
%SUM_RHOCV=[SUM_RHOCVaux(1:nv-1); sum(SUM_RHOCVaux(nv:nv+1))];

SUM_RHOCVsolids= volumeIntegral(sum(a_rhoCv(:,:,solids),3),intervals,r);
SUM_RHOCVgases= volumeIntegral(sum(a_rhoCv(:,:,gases),3),intervals,r);
SUM_RHOCV=[SUM_RHOCVsolids(1:nv-1)+...
    SUM_RHOCVgases(1:nv-1);...
    sum(SUM_RHOCVsolids(nv:nv+1))+...
    sum(SUM_RHOCVgases(nv:nv+1)) ];

dTdt=zeros(size(T));
dTdt(1:nv)=3*(RHOCVdPdrdTdr+P_DIVP+DivLambdadTdr(2:nv+1)...
    -DivLambdadTdr(1:nv)-SUM_URXN_RXNRATE)./SUM_RHOCV;
    
dPdt=(drhodt(:,gases)*(1./W(gases)).*T...
    +rho(:,gases)*(1./W(gases)).*dTdt)*Rgas./porosity...
    +rho(:,gases)*(1./W(gases))*Rgas.*T*(1-param.porosity0)...
    /param.rhob0./porosity.^2.*sum(drhodt(:,solids),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Displacement of the center of mass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
%dCMdt=( int_rw^re[drhodt*r^3]dr * int_rw^re[rho*r^2]dr - int_rw^re[drhodt*r^2]dr * int_rw^re[rho*r^3]dr   ) 
%       ---------------------------------------------------------------------------------------------------   
%                                       int_rw^re[rho*r^2]dr 

% int_rw^re[drhodt*r^3]dr=[r^3 K dPdr rho(gases)|_rw^re]/mu - int_rw^re[r^2rho(gases) K dPdr/mu] dr
a_KdPdr=multipleConvolution(a_K,a_dPdr);
a_KdPdrrho=multipleConvolution(a_KdPdr,sum(a_rho(:,:,gases),3));
KdPdrrho_r=evaluatePolynomials(a_KdPdrrho,intervals,r_ev);
r3KdPdrrho_r=r_ev.^3.*KdPdrrho_r;
R2RHOKdPdr=volumeIntegral(a_KdPdrrho,intervals,r);
R2RHOKdPdr=[R2RHOKdPdr(1:nv-1);sum(R2RHOKdPdr(nv+1))];
dRHOdtR3=(r3KdPdrrho_r(2:nv+1)-r3KdPdrrho_r(1:nv) -R2RHOKdPdr)/param.mu;

% int_rw^re[drhodt*r^2]dr =[r^2 K dPdr rho(gases)|_rw^re]/mu
r2KdPdrrho_r=r_ev.^2.*KdPdrrho_r;
dRHOdtR2=(r2KdPdrrho_r(2:nv+1)-r2KdPdrrho_r(1:nv))/param.mu;

%int_rw^re[rho*r^2]dr
RHOR2=volumeIntegral(sum(a_rho,3),intervals,r);
RHOR2=[RHOR2(1:nv-1);sum(RHOR2(nv:nv+1))];
    
%int_rw^re[rho*r^2]dr
a_rhor=[sum(a_rho,3) zeros(np,1)];
RHOR3=volumeIntegral(a_rhor,intervals,r);
RHOR3=[RHOR3(1:nv-1);sum(RHOR3(nv:nv+1))];  

dCMdt= (dRHOdtR3.*RHOR2 - dRHOdtR2.*RHOR3)./(RHOR2.^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Displacement of the center of thermal mass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% int_rw^re[drhoCvdt*r^3]dr=[r^3 K dPdr Cvrho(gases)|_rw^re]/mu - int_rw^re[r^2rho(gases) K dPdr/mu] dr
a_KdPdrrhoCv=multipleConvolution(a_KdPdr,sum(a_rhoCv(:,:,gases),3));
KdPdrrhoCv_r=evaluatePolynomials(a_KdPdrrhoCv,intervals,r_ev);
r3KdPdrrhoCv_r=r_ev.^3.*KdPdrrhoCv_r;

R2RHOCVKdPdr=volumeIntegral(a_KdPdrrhoCv,intervals,r);
R2RHOCVKdPdr=[R2RHOCVKdPdr(1:nv-1);sum(R2RHOCVKdPdr(nv+1))];
a_CvRrxn=0;
for i=1:nc
    a_CvRrxn= a_CvRrxn + multipleConvolution(a_Cv(:,:,i),a_Rrxn(:,:,i));
end

R3CVRRXN=volumeIntegral([a_CvRrxn zeros(np,1)],intervals,r);
R3CVRRXN=[R3CVRRXN(1:nv-1);sum(R3CVRRXN(nv:nv+1))];
dRHOdtCVR3=(R3CVRRXN+r3KdPdrrhoCv_r(2:nv+1)-r3KdPdrrhoCv_r(1:nv) -R2RHOKdPdr)/param.mu;

% int_rw^re[Cv drhodt*r^2]dr =[r^2 K dPdr rho(gases) Cv |_rw^re]/mu
r2KdPdrrhoCv_r=r_ev.^2.*KdPdrrhoCv_r;
R2CVRRXN=volumeIntegral(a_CvRrxn ,intervals,r);
R2CVRRXN=[R2CVRRXN(1:nv-1);sum(R2CVRRXN(nv:nv+1))];
dRHOdtCVR2=(r2KdPdrrho_r(2:nv+1)-r2KdPdrrho_r(1:nv))/param.mu+ R2CVRRXN;

%int_rw^re[rho*r^2]dr
RHOCVR2=volumeIntegral(sum(a_rhoCv,3),intervals,r);
RHOCVR2=[RHOCVR2(1:nv-1);sum(RHOCVR2(nv:nv+1))];
    
%int_rw^re[rho*r^2]dr
a_rhocvr=[sum(a_rhoCv,3) zeros(np,1)];
RHOCVR3=volumeIntegral(a_rhocvr,intervals,r);
RHOCVR3=[RHOCVR3(1:nv-1);sum(RHOCVR3(nv:nv+1))];  

dCMTdt= (dRHOdtCVR3.*RHOCVR2 - dRHOdtCVR2.*RHOCVR3)./(RHOCVR2.^2);


%--------------------------------------------------
% Return balances
%--------------------------------------------------

dXdt=zeros(size(X));
dXdt(solidnodes)=drhodt(:,solids);
%dXdt(gounodes)=drhodt(1:nv,gou);
%dXdt(gasnodes)=drhodt(1:nv,gas);
dXdt(ygasnode)=dydt;
dXdt(Pnodes)=dPdt(1:nv);
dXdt(Tnodes)=dTdt(1:nv);
%dXdt(CMnodes)=dCMdt;



end


