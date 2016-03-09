clear 
clc

order=3;
X=(0:1E-7:1E-6)';
n=length(X);
Y=300*(ones(size(X)));
Y(end)=450;
W=[0; mean([X(2:n) X(1:n-1)],2);X(end)+range(X)/2];
[a1, intervals1]=getPolynomials(X,Y, W,1, order,true);
[a2, intervals2]=getPolynomials(X*1E6,Y, W*1E6,1, order,true);
[a3, intervals3]=getPolynomials(X,Y, W,1E6, order,true);

[m n]=size(a1);
XX=ones(size(X));
WW=ones(size(W));
unit=1E6;
Units=ones([m 1]);
for i=1:order
    XX=[X.^i XX];
    WW=[W.^i WW];
    Units=[unit^i*ones([m 1]) Units];
end 