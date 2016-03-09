%Integrates a series of polynomia with characteristic constants, al over
%over the range of X used to calcuate "a"
function intY=integratethem(a,X,upwind)
%calculation of the integral
n=length(X);
xi=X;%[X(1);(X(1:n-1)+X(2:n))/2;X(end)];
[n m]=size(a);
repmatupwind=repmat(upwind(2:n),[1 m]);
a=[a(1,:);a(1:n-1,:).*repmatupwind + a(2:n,:).*~repmatupwind;a(end,:)] ;
[m n]=size(a);
m=length(xi);
XX=[];
for i=1:n
    XX=[ xi(2:m).^i/i-xi(1:m-1).^i/i XX];
end
intY=sum(a.*XX,2);
end