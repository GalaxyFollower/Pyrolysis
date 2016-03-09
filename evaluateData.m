function [yield]=evaluateData(t,r,X,param, graph)

n=length(r);
nv=param.nv;%number of volumes
nt=length(t);

%nodes indices
%nodes indices
bionodes=1:nv+1;
charnodes=max(bionodes)+1:max(bionodes)+nv+1;
solidnodes=[bionodes charnodes];
%gounodes=max(charnodes)+1:max(charnodes)+nv+1;
%gasnodes=max(gounodes)+1:max(gounodes)+nv+1;
ygasnode=max(charnodes)+1:max(charnodes)+nv+1;
%gasesnodes=[gounodes gasnodes ];
Pnodes=max(ygasnode)+1:max(ygasnode)+nv;
Tnodes=max(Pnodes)+1:max(Pnodes)+nv;
CMnodes=max(Tnodes)+1:max(Tnodes)+nv;

%Component indices
bio=1; %biomass
cha=2;
gou=3;%tar=goudron (tar is a matlab function)
gas=4;
gases=param.nsc+1:param.nc;
solids=1:param.nsc;

rho1=X(:,bionodes);
rho2=X(:,charnodes);
%rho3=X(:,gounodes);
ygas=X(:,ygasnode);
%rho4=X(:,gasnodes);
P=X(:,Pnodes);
T=[X(:,Tnodes) zeros([size(X,1) 1])];
Rgas=param.Rgas;
W=param.W;
porosity=1-(1-param.porosity0)*(rho1(:,:)+rho2(:,:))/param.rhob0;
P(:,nv+1)=param.Pinf;
rhotot=porosity.*1./(ygas/W(gas) +(1-ygas)/W(gou)).*P/Rgas./T;
rho3=(1-ygas).*rhotot;
rho4=ygas.*rhotot;
CM=X(:,CMnodes);

for i=1:size(X,1)
    clc 
    fprintf('data for T= %3.0f, %2.2f%% completed\n',param.Tinf,i*100/size(X,1))
    rho=[rho1(i,:)' rho2(i,:)' [rho3(i,:)'] [rho4(i,:)']];
    r_ev=[r(1:nv);CM(i,end);r(nv+1)];
    [T(i,end) rhogasesend]=limitTemperature(CM(i,:)',T(i,1:nv)',rho, ...
        ygas(i),r_ev,param);
    rho3(i,end)=rhogasesend(1);
    rho4(i,end)=rhogasesend(2);
end     






dPdr=zeros(nt,1);
for i=1:nt
    [a_P, intervals] =getPolynomials(CM(i,:)',P(i,:)',r,param.unit,param.order,param.deriv2);
    
    a_dPdr=differentiatePolynomials(a_P);
    dPdr(i)=evaluatePolynomials(a_dPdr,intervals,r(end));
end
K=(param.K(bio)-param.K(cha))*rho1(:,end)/param.rhob0+param.K(cha);
u=-K/param.mu.*dPdr;
porosity=1-(1-param.porosity0)*(rho1(:,:)+rho2(:,:))/param.rhob0;
M3out=porosity(:,end).*rho3(:,end).*u*4*pi*param.r^2;
M4out=porosity(:,end).*rho4(:,end).*u*4*pi*param.r^2;

M40=param.porosity0*rho4(1,1)*(4*pi*param.r^3/3);

[a_rho1, intervals] =getPolynomials(CM(i,:)',rho1(end,:)',r,param.unit,...
    param.order,param.deriv2);
M1tot=4*pi*sum(volumeIntegral(a_rho1,intervals,intervals));

[a_rho2, intervals] =getPolynomials(CM(i,:)',rho2(end,:)',r,param.unit,...
    param.order,param.deriv2);
M2tot=4*pi*sum(volumeIntegral(a_rho2,intervals,intervals));

[a_rho3, intervals] =getPolynomials(CM(i,:)',(porosity(end,:).*rho3(end,:))',r,param.unit,...
    param.order,param.deriv2);
M3tot=4*pi*sum(volumeIntegral(a_rho3,intervals,intervals));

[a_rho4, intervals] =getPolynomials(CM(i,:)',(porosity(end,:).*rho4(end,:))',r,param.unit,...
    param.order,param.deriv2);
M4tot=4*pi*sum(volumeIntegral(a_rho4,intervals,intervals));




M3exits=trapz(t,M3out);
M4exits=trapz(t,M4out);
M0=param.rhob0*(4*pi*param.r^3/3);
yield=[M1tot; M2tot;M3exits;M4exits]/M0;



%CM=X(:,5*(nv+1)+1:6*(nv));
if(graph)
figure(1)
hdl=subplot(2,2,1);
gimmiedagraph(hdl,t,rho1,{'\rho_B  (kg\cdotm^{-3})','\rho_B  (kg\cdotm^{-3})'},false)
hdl=subplot(2,2,2);
gimmiedagraph(hdl,t,rho2,{'Time (s)','\rho_C  (kg\cdotm^{-3})'},false)
hdl=subplot(2,2,3);
gimmiedagraph(hdl,t,rho3,{'Time (s)','\rho_T  (kg\cdotm^{-3})'},true)
hdl=subplot(2,2,4);
gimmiedagraph(hdl,t,rho4,{'Time (s)','\rho_G  (kg\cdotm^{-3})'},true)
gases=param.nsc+1:param.nc;

figure(2)
hdl=subplot(2,2,1);
gimmiedagraph(hdl,t,P-param.Pinf,{'Time (s)','Pressure (Pa)'},true)
%plot(t,P)
hdl=subplot(2,2,2);
gimmiedagraph(hdl,t,rho1./(rho1+rho2),{'Time (s)',{'Fraction of biomass';'in solid phase ( - )'}},false)
hdl=subplot(2,2,3);
gimmiedagraph(hdl,t,T,{'Time (s)','Temperature (K)'},false)

figure(5)
plot(t,[M3out M4out])


end
end

function gimmiedagraph(hdl,X,Y,axisnames,semi_log)
    if(semi_log)
        semilogx(hdl,X,Y)
    else
        plot(hdl,X,Y)
    end
    xlabel(axisnames{1})
    ylabel(axisnames{2})
       
end