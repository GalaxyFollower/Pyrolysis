%Evaluate the values of the polynomiae with characteristic constants a
% at the intervals X
function [y]=evaluatethemsimply(a,X)

[n m]=size(a);
al=[a(1,:);a;a(end,:)] ;
n=length(X);

XX=ones(size(X));
for i=2:1:m
    XX=[ X.^(i-1) XX];
end
y=sum(al.*XX,2);