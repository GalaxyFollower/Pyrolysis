function keepemGoing
    kinetics(1).name=['Chan et al'];
    kinetics(1).A=[1.3e8 2e8 1.08e7 0];
    kinetics(1).E=[140e3 133e3 121e3 0];


    %Enthalpies of reaction [J/kg]
    kinetics(:).Hrxn=[538e3 0 -2E6 0 ] ;%[Milosavljevic; 1996]

    r=[ 1; 5; 10; 15; 25; 35; 50]*1E-6; 
    Tinf=(750:50:1200)';
    nvert=[1;5; 10];
    
    ngoodones=0;
    i=fullfact([length(kinetics) length(Tinf) length(r)  length(nvert)]);
    for j=1:size(i,1)
        name=[kinetics(i(j,1)).name '_Tinf=' num2str(Tinf(i(j,2)))...
        '_r=' num2str(r(i(j,3))) '_nvter=' num2str(nvert(i(j,4))) '.mat'];
        if (exist(['Failed/' name])==2)
            fprintf(['\n\n' name '\n'])
            plotParticulesData(['Failed/' name] )
            
            goon=[];
            while isempty(goon)
                answer=input('Should Go on?\n>','s');
                if (strcmp(answer,'n') || strcmp(answer,'y') )
                    goon=strcmp(answer,'y');
                end                
            end
            haveitgoing(['Failed/' name],goon)
            if(goon)
                ngoodones=ngoodones+1;
                listOfGoodOnes{ngoodones}={name};
            else
                eraseit=[];
                while isempty(eraseit)
                	answer=input('Erase it?\n>','s');
                    if (strcmp(answer,'n') || strcmp(answer,'y') )
                       eraseit=strcmp(answer,'y');
                    end                
                end
                if(eraseit)
                    delete(['Failed/' name] );
                    delete(name)                    
                end
            end
        end


    end
    try
        save('ListOfGoodOnes.mat','listOfGoodOnes')
    catch err
        err
    end
    
    
end

function haveitgoing(name,goon)
    load(name)
    shouldgoon=goon;
    save(name)
end