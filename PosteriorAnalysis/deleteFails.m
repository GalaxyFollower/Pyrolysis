
kinetics(1).name=['Chan et al'];
kinetics(1).A=[1.3e8 2e8 1.08e7 0];
kinetics(1).E=[140e3 133e3 121e3 0];


%Enthalpies of reaction [J/kg]
kinetics(:).Hrxn=[538e3 0 -2E6 0 ] ;%[Milosavljevic; 1996]

r=[ 1; 5; 10; 15; 25; 35; 50]*1E-6; 
%Tinf=(750:50:1200)';
nvert=[1;5; 10];

i=fullfact([length(kinetics) length(Tinf) length(r)  length(nvert)]);
for j=1:size(i,1)
	name=[kinetics(i(j,1)).name '_Tinf=' num2str(Tinf(i(j,2)))...
    '_r=' num2str(r(i(j,3))) '_nvter=' num2str(nvert(i(j,4))) '.mat'];
    if (exist(name)==2)
        deleteFile(name)
    end            
end