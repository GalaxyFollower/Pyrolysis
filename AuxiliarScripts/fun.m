% function for ealuating various sets of polynomials for various variables.
function Y=fun(a,intervals,x,i)
    k=1;
    Y=zeros([size(x,1) length(i)]);
    for j=i
        Y(:,k)=evaluatePolynomials(a(:,:,j),intervals,x);
        k=k+1;
    end
end