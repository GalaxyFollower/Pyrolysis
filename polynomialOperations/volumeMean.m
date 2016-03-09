function [volmean]=volumeMean(a,X,CM,I)

[n m]=size(a);

% intY= sum(a.*(xi(2:n).^3-xi(1:n-1).^3)/3+b.*(xi(2:n).^2-xi(1:n-1).^2)/2+c.*(xi(2:n)-xi(1:n-1)));
%Calculation of the integral over the volume of the sphere
q=m-1;

A=[a.*repmat(1./(q+3:-1:3),[n 1]) zeros(n,3)];

n=size(CM,1);

volmean=zeros([n 1]);

for i=1:length(I)-1
	j=find(CM>=I(i) & CM<=I(i+1));
    volmean(j)=3*(polyval(A(i,:),X(j+1))-polyval(A(i,:),X(j)))./(X(j+1).^3-X(j).^3);    
end


end

