function deleted=deleteFile(name)
load(name)
if ~succes
    delete(name)
end 
deleted=~succes;