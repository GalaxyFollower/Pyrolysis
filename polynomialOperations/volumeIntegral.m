function [Y, aint]=volumeIntegral(a,I,x)
[np q]=size(a);
aint=[a./repmat((q+2:-1:3),[np 1]) zeros([np 3])];

n=length(x);
center=mean([x(2:n) x(1:n-1)],2);

Y=zeros([n-1 1]);
for i=1:np
    j=find(center>I(i) & center<=I(i+1));
    if(isempty(j)), 
    error(['The integration limits are not in accordance with the' ...
        'validity of the interval validity of each polynomial']);
    end
        
    Y(j)=polyval(aint(i,:),x(j+1))-polyval(aint(i,:),x(j));
end

