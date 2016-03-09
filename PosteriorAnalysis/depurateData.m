%Takes a saved workspace from a simulationsaved under a 'name(.mat)' file
%, erases the useless variables and depurates the results into usable 
%variables before saving the results in a the 'saveFile(.mat)' file.
function succes=depurateData(name,j)
load(name)
%Clenning up

clear X deriv2 hdl i n odeoptions order P0 plotem t T t0 T0 tf Tinf
clear vtermult X0
ri=ri;
nt=length(time);
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
%gounodes=max(charnodes)+1:max(charnodes)+nv+1;
%gasnodes=max(gounodes)+1:max(gounodes)+nv+1;
ygasnode=max(charnodes)+1:max(charnodes)+nv+1;
%gasesnodes=[gounodes gasnodes ];
Pnodes=max(ygasnode)+1:max(ygasnode)+nv;
Tnodes=max(Pnodes)+1:max(Pnodes)+nv;
CMnodes=max(Tnodes)+1:max(Tnodes)+nv;

rho(:,:,bio)=XX(:,bionodes);
rho(:,:,cha)=XX(:,charnodes);
porosity=1-(1-param.porosity0)*sum(rho(:,:,solids),3)/param.rhob0;
ygas=XX(:,ygasnode);
P=[XX(:,Pnodes) zeros(nt,1)];
T=[XX(:,Tnodes) zeros(nt,1)];
rhotot=(P+param.Pinf).*porosity/param.Rgas./T./(ygas/param.W(gas)+(1-ygas)/param.W(gou));

rho(:,:,gou)=(1-ygas).*rhotot;
rho(:,:,gas)=(ygas).*rhotot;

CM=CM0;
%CM=XX(:,CMnodes);


for i=1:nt
	clc 
    fprintf('files analysed %2.2f%%\n',j)
    fprintf('data for %s, %2.2f%% completed\n',name,i*100/nt)
    rhoi=reshape(rho(i,:,:),[param.nv+1 param.nc]);
    [T(i,end), rhoend]=limitTemperature(CM,T(i,:)',rhoi, ygas(i,end),ri,param);
    rho(i,end,gases)=reshape(rhoend,[1 1 param.ngc ]);
end
clear XX ygas ygasnode

if(succes)
    save(['Simulations\Successful\' strrep(name,'Simulations\','')]);
else
    save(['Simulations\Failed\' strrep(name,'Simulations\','')]);
end


end




