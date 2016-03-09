function Y=fixLast(X,Y, Yend, W )
n=length(Y);


[a, I]=getPolynomials(X(1:n-3),Y(1:n-2), [W(1:n-3); X(n-2)],1, 3,true);
I(end)=W(n-1);
dadx=differentiatePolynomials(a);
d2adx2=differentiatePolynomials(dadx);
b=[ Yend
    evaluatePolynomials(a,I,W(n-1))
    evaluatePolynomials(dadx,I,W(n-1))
    evaluatePolynomials(d2adx2,I,W(n-1))
    ];

A=[W(n)^3 W(n)^2 W(n) 1
   W(n-1)^3 W(n-1)^2 W(n-1) 1
 3*W(n-1)^2 2*W(n-1) 1 0
 6*W(n-1) 2 0 0 ];
c=(A\b)';
Y(end)=Yend;
Y(n-1)= polyval(c,X(n-1));
end