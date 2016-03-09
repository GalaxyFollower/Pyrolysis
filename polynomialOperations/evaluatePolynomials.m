function f=evaluatePolynomials(A,I,x)


f=zeros(size(x));
f(:)=nan;
for i=1:length(I)-1
	j=x>=I(i) & x<=I(i+1);
    f(j)=polyval(A(i,:),x(j));
end

