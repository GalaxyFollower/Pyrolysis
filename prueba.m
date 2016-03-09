clear
clc
X=(-3:0.1:3)';
Y=normpdf(X);
u=ones(size(X));
figure(1)
plot(X,Y,'*')
hold on

n=length(Y);
upwind=(u(1:n-1)>0 ) & ( u(2:n)>0 | u(1:n-1)>abs(u(2:n)));

[a]=polynomiame(X,Y);
[cm al]=centerOfMass(a,X,upwind);
[y]=evaluatethem(a,X,cm,upwind);
plot(cm,y,'og');

[y]=evaluatethem(a,X,cm,~upwind);
plot(cm,y,'.r');

intY=sum(integratethem(a,X,~upwind))
[volint]=volumeintegral(X,Y,upwind)

plot(cm,volint,'sk')
hold off
