clc 
clear
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

for i=1:length(kinetics)
    kinetics(i).Hrxn=[538e3 0 -2E6 0 538e3] ;%[Milosavljevic; 1996]
end

%kinetics=kinetics(4);
r=[2.5E-6
   5E-6 
   10E-6 
   15E-6 
   25E-6 
   35E-6 
   50E-6];
Tinf=(750:50:1200)';
nvert=[0.001;0.25;1;5;10];
%nvert=1;
i=fullfact([length(kinetics) length(Tinf) length(r)  length(nvert)]);



fi=numel(dir('*.mat'));
list=boolean(zeros(length(i),4));

for j=1:length(i)
 	name{j}=[kinetics(i(j,1)).name '_Tinf=' num2str(Tinf(i(j,2)))...
    '_r=' num2str(r(i(j,3))) '_nvter=' num2str(nvert(i(j,4))) '.mat'];
end



parfor j=1:size(i,1)
    succes=[];
	if (exist(['Simulations\' name{j} ])==2 || ...
        exist(['Simulations\Successful\' name{j}])==2 ||...
        exist(['Simulations\Failed\' name{j}])==2 )
        if(exist(['Simulations\Successful\' name{j}])==2)
           succes=true;
        elseif(exist(['Simulations\Failed\' name{j}])==2) 
           succes=false;
         else
            try
                succes=depurateData(['Simulations\' name{j}],j);   
                
            catch err                
                
            end
         end
         
        if succes 
            list(j,:)=[1 0 0 0];
        elseif ~succes
            list(j,:)=[0 1 0 0];
        else
            list(j,:)=[0 0 1 0];
        end

    else
              list(j,:)=[0 0 0 1];
    end
    
    
end
ngood=0;
nbad=0;
nleft=0;
nFckdUp=0;
for j=1:length(i)
    if(list(j,1))
        ngood=ngood+1;
        listOfGoodOnes{ngood,1}=kinetics(i(j,1)).name;
        listOfGoodOnes{ngood,2}=Tinf(i(j,2));
        listOfGoodOnes{ngood,3}=r(i(j,3));
        listOfGoodOnes{ngood,4}=nvert(i(j,4));
        listOfGoodOnes{ngood,5}=['Simulations\Successful\' name{j}];
    end
    if(list(j,2))
        nbad=nbad+1;            
        listOfBadOnes{nbad,1}=kinetics(i(j,1)).name;
        listOfBadOnes{nbad,2}=Tinf(i(j,2));
        listOfBadOnes{nbad,3}=r(i(j,3));
        listOfBadOnes{nbad,4}=nvert(i(j,4));
        listOfBadOnes{nbad,5}=['Simulations\Failed\' name{j}];
    end
    if(list(j,3))
        nFckdUp=nFckdUp+1;
        listOfFckdUpOnes{nFckdUp,1}=kinetics(i(j,1)).name;
        listOfFckdUpOnes{nFckdUp,2}=Tinf(i(j,2));
        listOfFckdUpOnes{nFckdUp,3}=r(i(j,3));
        listOfFckdUpOnes{nFckdUp,4}=nvert(i(j,4));
        listOfFckdUpOnes{nFckdUp,5}=name{j};
    end
    if(list(j,4))
        nleft=nleft+1;            
        listOfLeftOnes{nleft,1}=kinetics(i(j,1)).name;
        listOfLeftOnes{nleft,2}=Tinf(i(j,2));
        listOfLeftOnes{nleft,3}=r(i(j,3));
        listOfLeftOnes{nleft,4}=nvert(i(j,4));
        listOfLeftOnes{nleft,5}=name{j};
    end
end    
    


saved=false;
while ~saved
    try
        if(exist('listOfGoodOnes')||exist('listOfBadOnes'))
            xlswrite('Results.xls',[listOfGoodOnes;listOfBadOnes],'all');
        end
        if(exist('listOfGoodOnes'))
            xlswrite('Results.xls',listOfGoodOnes,'Good');
        end
        if(exist('listOfBadOnes'))
            xlswrite('Results.xls',listOfBadOnes,'Bad');
        end
        if(exist('listOfLeftOnes'))
            xlswrite('Results.xls',listOfLeftOnes,'left');
        end
        saved=true;
    catch err
        err
        fprintf('Close Results.xls you ass, and then come back and press a fckin key if U dare \n')
        pause
    end
end


    
    

