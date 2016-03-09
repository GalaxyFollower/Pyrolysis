%% Function that estimates the polynomial constants that describe the rates
%of reaction.
%   a_Rrxn: represent the constants for the rate of reaction of each
%       component
%   a_rxn: The speed of each reaction.
function [ a_rxn, a_Rrxn] = rxnRate(a_rho,T,r,CM,param)
%Component indices
bio=1; %biomass
cha=2;
gou=3;
gas=4;

k=exp(-(1./(param.Rgas*T))*param.kinetics.E).*repmat(param.kinetics.A,[param.nv+1 1]);
%CMx=CM(1:param.nv-1);
%rx=[r(1:param.nv-1);r(end)];

a_k=zeros([param.np param.order+1 4]);
a_rxn=zeros([param.np 2*param.order+1 4]);
x=0:param.r/(1000):param.r;


for i=1:4
    
    [a_k(:,:,i) intervals]=getPolynomials(CM,k(:,i),r,param.unit,param.order);
    if(i<=3)
        a_rxn(:,:,i)=multipleConvolution(a_k(:,:,i),a_rho(:,:,bio));    
    else
        a_rxn(:,:,i)=multipleConvolution(a_k(:,:,i),a_rho(:,:,gou));    
    end
%     
%     if(i<=3)
%         [a_rxn(:,:,i),intervals]=getPolynomials(CM,k(:,i).*rho(:,bio),r,param.unit,param.order);    
%     else
%         [a_rxn(:,:,i),intervals]=getPolynomials(CM,k(:,i).*rho(:,gou),r,param.unit,param.order);    
%     end
%     
  
end

a_Rrxn(:,:,bio)=-sum(a_rxn(:,:,1:3),3);
a_Rrxn(:,:,cha)=a_rxn(:,:,3);
a_Rrxn(:,:,gou)=a_rxn(:,:,2)-a_rxn(:,:,4);
a_Rrxn(:,:,gas)=a_rxn(:,:,1)+a_rxn(:,:,4);
    
end