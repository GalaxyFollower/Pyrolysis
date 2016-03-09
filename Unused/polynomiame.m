function [a]=polynomiame(X,Y, order)
%vertical vectors
X=X(:);
Y=Y(:);

    n=length(Y);
if(order==2)
    %dividing vectors

    Ye=Y(1:n-2);
    Yp=Y(2:n-1);
    Yw=Y(3:n);

    Xe=X(1:n-2);
    Xp=X(2:n-1);
    Xw=X(3:n);

    %calculation of polinomial constants

    a(:,1)=Ye./((Xe-Xp).*(Xe-Xw))+Yp./((Xp-Xe).*(Xp-Xw))+Yw./((Xw-Xe).*(Xw-Xp));
    a(:,2)=-(Ye.*(Xp+Xw)./((Xe-Xp).*(Xe-Xw))+Yp.*(Xe+Xw)./((Xp-Xe).*(Xp-Xw))+Yw.*(Xe+Xp)./((Xw-Xe).*(Xw-Xp)));
    a(:,3)=Ye.*Xp.*Xw./((Xe-Xp).*(Xe-Xw))+Yp.*Xe.*Xw./((Xp-Xe).*(Xp-Xw))+Yw.*Xe.*Xp./((Xw-Xe).*(Xw-Xp));
elseif(order==1)
    Ye=Y(1:n-1);
    Yp=Y(2:n);
    Xe=X(1:n-1);
    Xp=X(2:n);
    
    a(:,1)=(Yp-Ye)./(Xp-Xe);
    a(:,2)=(Ye.*Xp-Yp.*Xe)./(Xp-Xe);
    
end
end
