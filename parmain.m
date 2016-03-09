
try
   matlabpool open local
catch err
    
end
    
%%
addpath('polynomialOperations')
addpath('AuxiliarScripts')
addpath('StarchStructure')
addpath('PosteriorAnalysis')
%[1]:Biomass->gas
%[2]:Biomass->tar
%[3]:Biomass->char
%[4]:tar->gas

kinetics(1).name=['Chan'];
kinetics(1).A=[1.3e8 2e8 1.08e7 0];
kinetics(1).E=[140e3 133e3 121e3 0];

kinetics(2).name=['Wangenaar'];
kinetics(2).A=[1.11E11 9.28E9 3.05E7 0];
kinetics(2).E=[177E3 149e3 125e3 0];

kinetics(3).name=['DiBlasi'];
kinetics(3).A=[4.4E9 1.1E10 3.3E6 0];
kinetics(3).E=[153E3 148e3 112e3 0];

kinetics(4).name=['Font'];
kinetics(4).A=[1.52E7 5.85E6 2.98E3 0];
kinetics(4).E=[139E3 119e3 73e3 0];

kinetics=kinetics(4);


%Enthalpies of reaction [J/kg]
for i=1:length(kinetics)
    kinetics(i).Hrxn=[538e3 0 -2E6 538e3];%[Milosavljevic; 1996]
end

r=[ 2.5E-6
    5E-6
    10E-6
    15E-6
    25E-6
    35E-6
    50E-6]; 

Tinf=[(750:50:1000)';1100;1200];
Tinf=900;
nvert=[0.001 0.25 1 5 10];
%nvert=1;
i=fullfact([length(kinetics) length(Tinf) length(r)  length(nvert)]);
i=i(randperm(size(i,1)),:);


yield=zeros(3,length(Tinf));
time=zeros(length(Tinf),1);
succes=zeros(length(Tinf),1);

parfor j=1:size(i,1)
	name=['Simulations\' kinetics(i(j,1)).name '_Tinf=' num2str(Tinf(i(j,2)))...
    '_r=' num2str(r(i(j,3))) '_nvter=' num2str(nvert(i(j,4))) '.mat'];
    if ~(exist(name)==2)

        changing(Tinf(i(j,2)),kinetics(i(j,1)),r(i(j,3)),nvert(i(j,4)));

            
    elseif deleteFile(name)
    %         changing(Tinf(i(j,2)),kinetics(i(j,1)),r(i(j,3)),nvert(i(j,4)));        
    end

end
try    
    matlabpool close
catch err
   close parpool
end

