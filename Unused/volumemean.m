function [volmean]=volumemean(a,X,upwind)


n=length(X);
xi=X;

[n m]=size(a);
a=[a zeros([n 2])];
m=m+2;
repmatupwind=repmat(upwind(2:n),[1 m]);
a=[a(1,:);a(1:n-1,:).*repmatupwind + a(2:n,:).*~repmatupwind;a(end,:)] ;

% intY= sum(a.*(xi(2:n).^3-xi(1:n-1).^3)/3+b.*(xi(2:n).^2-xi(1:n-1).^2)/2+c.*(xi(2:n)-xi(1:n-1)));
%Calculation of the integral over the volume of the sphere
[m n]=size(a);
XX=[];
for i=1:n
    XX=[ xi(2:m+1).^i/i-xi(1:m).^i/i XX];
end
volmean=3*sum(a.*XX,2)./( xi(2:m+1).^3-xi(1:m).^3);

