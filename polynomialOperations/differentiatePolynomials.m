%returns the coefficients of the polynomial resulting from the derivation
%of the polynpmial withcoefficiets a
function dadx=differentiatePolynomials(a)
[m n]=size(a);
dadx=a(:,1:n-1).*repmat([n-1:-1:1],[m 1]);