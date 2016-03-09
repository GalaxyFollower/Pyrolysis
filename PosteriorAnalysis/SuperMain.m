clc
clear
SortFiles
Completelist=[listOfGoodOnes;listOfBadOnes];
%organise the list
[~,i]=sort(cell2mat(Completelist(:,4)));
Completelist=Completelist(i,:);
[~,i]=sort(cell2mat(Completelist(:,2)));
Completelist=Completelist(i,:);
[~,i]=sort(cell2mat(Completelist(:,3)));
Completelist=Completelist(i,:);
kineticsindex=[];
for i=1:length(kinetics)
    kineticsindex= [kineticsindex;find(cellfun(@(x) strcmp(x,kinetics(i).name)...
        , Completelist(:,1)))];
end
Completelist=Completelist(kineticsindex,:);

[time, yield]=compareExploitedData(Completelist);
name=['AllExploiteddata.mat'];
save(name);
%%
for i=1:length(kinetics)
    FinalGraphs(name,kinetics(i).name,true);
end
