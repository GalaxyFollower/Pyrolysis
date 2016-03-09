function [hdl ] = getFigureHdl(name)
%returns the handle for a figure with the name given, if the figure doen't
%exist, it creates one figure with that name,
hdl=findobj('name',name);

if(isempty(hdl))
    hdl=figure('name',name);
else
    figure(hdl);
end


