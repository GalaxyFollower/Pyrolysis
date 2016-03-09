%Evaluate the values of the polynomiae with characteristic constants a
% between the intervals X at their center of mass cm 
function [y]=evaluatethem(a,X,upwind)

%selecting constants according to the direcction of the speed

[n m]=size(a);
repmatupwind=repmat(upwind(2:n),[1 m]);

al=[a(1,:);a(1:n-1,:).*repmatupwind + a(2:n,:).*~repmatupwind;a(end,:)] ;

XX=ones(size(X));
for i=2:1:m
    XX=[ X.^(i-1) XX];
end
y=sum(al.*XX,2);
end