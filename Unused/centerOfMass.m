function [cm al]=centerOfMass(a,X,upwind)

%selecting constants according to the direcction of the speed

[n m]=size(a);
repmatupwind=repmat(upwind(2:n),[1 m]);

al=[a(1,:);a(1:n-1,:).*repmatupwind + a(2:n,:).*~repmatupwind;a(end,:)] ;
n=length(X);

intDown=[];
intUp=[];
for i=m:-1:1
    intUp=[intUp (X(2:n).^(i+3)-X(1:n-1).^(i+3))/(i+3) ];
    intDown=[intDown (X(2:n).^(i+2)-X(1:n-1).^(i+2))/(i+2) ];
end
cm=sum(al.*intUp,2)./sum(al.*intDown,2);