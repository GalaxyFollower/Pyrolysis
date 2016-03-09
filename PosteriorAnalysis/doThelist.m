function doThelist(listOfGoodOnes)

try
    poolobj = parpool;
catch
    matlabpool open
end

gottaproblem=zeros(length(listOfGoodOnes),1);
parfor i=1:length(listOfGoodOnes)
    try
        goFarther(char(listOfGoodOnes{i}))
    catch err        
        gottaproblem(i)=true;
    end
end

try
    delete(poolobj);
catch err
    matlabpool close
end

try
    fprintf('list of fckups\n')
    listOfGoodOnes{boolean(gottaproblem)}
catch err
    err
end

