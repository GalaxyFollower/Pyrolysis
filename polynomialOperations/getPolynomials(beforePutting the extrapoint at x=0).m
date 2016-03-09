function [a, intervals]=getPolynomials(X,Y, W,unit, order,second)

if(second)
    [a, intervals]=getPoly2_2ndOrder(X,Y, W,unit, order);
else
    [a,interval]=getPoly2_1stOrder(X,Y, W, order);
end

end
   
function [a, intervals]=getPoly2_2ndOrder(X,Y, W, unit, order)
%putting one point before w(0) that mirror X(1)

%     X=[W(1)-abs(W(1)-X(1));X];
%     Y=[Y(1);Y];

n=length(X);
m=length(W);
X=X*unit;
W=W*unit;

%verify that the number of points is in accord with the order of the
%polynomials imposed
if(mod(n-1,order-2)~=0)
    error(['If ''n'' is the number of volumes (' num2str(n) ') and ''q'' the order of the polynomials (' num2str(order) '), then ''n-1'' must be a multiple of ''q-2'' ']);
end 
np=(n-1)/(order-2)+1;%number of polynomials

%Initialize the matrices
A=zeros(np*(order+1));
b=zeros([size(A,1) 1]);
intervals=zeros([np+1 1]);

XX=ones(size(X));
WW=ones(size(W));
Units=ones([np 1]);
for i=1:order
    XX=[X.^i XX];
    WW=[W.^i WW];
    Units=[unit^i*ones([np 1]) Units];
end 

DW=[repmat((order:-1:1),[m 1]).*WW(:,2:end) zeros(size(W))];
D2W=[repmat((order:-1:2).*(order-1:-1:1),[m 1]).*WW(:,3:end) zeros([length(W) 2])];


line=1; %count the lines writen in the matrix
j=1;%count used points
for i=1:np
    if(i==1)
        %dP1dx(W0)=0
        A(line,1:order+1)=DW(1,:);        
        intervals(1)=W(1);
        line=line+1;
        %P1(x1)=y1 extrapoint at the begining
        A(line,1:(order+1))=XX(1,:); 
        b(line)=Y(1);
        line=line+1;
        j=j+1;
        %Pi(xj)=yj
        A(line:line+order-3,(i-1)*(order+1)+1:(i)*(order+1))=XX(j:j+order-3,:); 
        b(line:line+order-3)=Y(j:j+order-3);
        line=line+order-2;
        j=j+order-2;
        %Pi(Wi)=Pi+1(Wi)

        k=W>X(j-1) & W<X(j);
        A(line,(i-1)*(order+1)+1:i*(order+1))=WW(k,:);
        A(line,(i)*(order+1)+1:(i+1)*(order+1))=-WW(k,:);
        intervals(i+1)=W(k);
        line=line+1;

        %dP1dx(Wi)=dP2dx(Wi)
        A(line,(i-1)*(order+1)+1:i*(order+1))=DW(k,:);
        A(line,(i)*(order+1)+1:(i+1)*(order+1))=-DW(k,:);
        line=line+1;

        %d2P1dx2(Wi)=d2P2dx2(Wi)
        A(line,(i-1)*(order+1)+1:i*(order+1))=D2W(k,:);
        A(line,(i)*(order+1)+1:(i+1)*(order+1))=-D2W(k,:);
        line=line+1;
    elseif(i<np-1)
        %Pi(xj)=yj
        A(line:line+order-3,(i-1)*(order+1)+1:(i)*(order+1))=XX(j:j+order-3,:); 
        b(line:line+order-3)=Y(j:j+order-3);
        line=line+order-2;
        j=j+order-2;
        
        %Pi(Wi)=Pi+1(Wi)

        k=W>X(j-1) & W<X(j);

        A(line,(i-1)*(order+1)+1:i*(order+1))=WW(k,:);
        A(line,(i)*(order+1)+1:(i+1)*(order+1))=-WW(k,:);
        intervals(i+1)=W(k);
        line=line+1;

        %dP1dx(Wi)=dP2dx(Wi)
        A(line,(i-1)*(order+1)+1:i*(order+1))=DW(k,:);
        A(line,(i)*(order+1)+1:(i+1)*(order+1))=-DW(k,:);
        line=line+1;

        %d2P1dx2(Wi)=d2P2dx2(Wi)
        A(line,(i-1)*(order+1)+1:i*(order+1))=D2W(k,:);
        A(line,(i)*(order+1)+1:(i+1)*(order+1))=-D2W(k,:);
        line=line+1;

    elseif(i==np)
%         %finish with 1st order polynomial
%         A(line:line+order-2,(i-1)*(order+1)+1:(i)*(order+1)-2)=eye(order-1);
%         A(end,i*(order+1)-1:i*(order+1))=WW(end,order:order+1);
%         b(end)=Y(end);
%         intervals(end)=W(end);  
        
        %finish with 2nd order polynomial
        A(line:line+order-3,(i-1)*(order+1)+1:(i)*(order+1)-3)=eye(order-2);
        A(end,i*(order+1)-2:i*(order+1))=WW(end,order-1:order+1);
        b(end)=Y(end);
        intervals(end)=W(end);  
    else
%         %Extra point at the end
%         %P1(x1)=y1
%         A(line,(i-1)*(order+1)+1:(i)*(order+1))=XX(j,:); 
%         b(line)=Y(j);
%         line=line+1;
%         j=j+1;
        
        %Pi(xj)=yj
        A(line:line+order-3,(i-1)*(order+1)+1:(i)*(order+1))=XX(j:j+order-3,:); 
        b(line:line+order-3)=Y(j:j+order-3);
        line=line+order-2;
        j=j+order-2;
        
%         %finish with 1st order polynomial
%         %Pi(Wi)=Pi+1(Wi)
%         %k=W>X(j-1) & W<X(j);
%         A(line,(i-1)*(order+1)+1:i*(order+1))=WW(np+1,:);
%         A(line,(i)*(order+1)+order:(i+1)*(order+1))=-WW(np+1,order:order+1);
%         intervals(i+1)=W(np+1);
%         line=line+1;
% 
%         %dP1dx(Wi)=dP2dx(Wi)
%         A(line,(i-1)*(order+1)+1:i*(order+1))=DW(np+1,:);
%         A(line,(i)*(order+1)+order:(i+1)*(order+1))=-DW(np+1,order:order+1);
%         line=line+1;
        
        %finish with 2nd order polynomial
        %Pi(Wi)=Pi+1(Wi)
        %k=W>X(j-1) & W<X(j);
        A(line,(i-1)*(order+1)+1:i*(order+1))=WW(np+1,:);
        A(line,(i)*(order+1)+order-1:(i+1)*(order+1))=-WW(np+1,order-1:order+1);
        intervals(i+1)=W(np+1);
        line=line+1;      

        %dP1dx(Wi)=dP2dx(Wi) 
        A(line,(i-1)*(order+1)+1:i*(order+1))=DW(np+1,:);
        A(line,(i)*(order+1)+order-1:(i+1)*(order+1))=-DW(np+1,order-1:order+1);
        line=line+1;

        %d2P1dx2(Wi)=d2P2dx2(Wi)
        A(line,(i-1)*(order+1)+1:i*(order+1))=D2W(np+1,:);
        A(line,(i)*(order+1)+order-1:(i+1)*(order+1))=-D2W(np+1,order-1:order+1);
        line=line+1;
         
    end   
end
acnst=zeros(size(b));
% if (range(Y)==0)
%    acnst(order+1:order+1:end)=Y(1);
% else
%     acnst=A\b;
% end
acnst=A\b;
%acnst=findActualValues(A,b,order);
%acnst=A\b;
ncnst=length(acnst);
%a=[reshape(acnst(1:ncnst-3),order+1,np-1)'; [ zeros(1,2) acnst(ncnst-2:end)'] ];
a=reshape(acnst,[order+1 np])';
a=a.*Units;
intervals=intervals/unit;
end

%function to invert matrix and to correct equation systemes for which the
%function is rather a constant more than a polynomial
function x=findActualValues(A,b,order)
    
    
    x=A\b;
    
%     %are the last constants similar to each other?
%     n=length(x);
%     a_nplus1_similar=abs(range(x(order+1:order+1:end))/mean(x(order+1:order+1:end)))<1e-4 ...
%         | isnan(range(x(order+1:order+1:end))/mean(x(order+1:order+1:end)));
%     %are the rest of the constants practically zero?
%     j=setdiff(1:n-3,order+1:order+1:n-3);
% 
%     a_i_zero= abs(mean(x(j)))==0 | (range(x(j)) & abs(mean(x(j)))<1e-5);
% 
%         
%     if(a_nplus1_similar && a_i_zero)
%         y=x;
%         y([order+1:order+1:n-3 n])=1;
%         y(setdiff(1:n,[order+1:order+1:n-3 n]))=0;
%         fun1=@(val,A,y,b)(1e10*sum(abs(A*(y*val)-b)));
%         fun_pre=10e20*sum(abs(A*x-b));
%         val=fminsearch(@(val)fun1(val,A,y,b),round(1E5*mean(x(order+1:order+1:end)))/1E5);
%         %val=fminsearch(@(val)fun1(val,A,y,b),mean(x(order+1:order+1:end)));
%         fun_post=fun1(val,A,y,b);
%         if (fun_post<fun_pre),x=val*y;end 
%     end
    
    
end

function a=getPoly2_1stOrder(X,Y, W, order)
n=length(X);


%verify that the number of points is in accord with the order of the
%polynomials imposed
if(mod(n-1,order-1)~=0)
    error(['If ''n'' is the number of points provided (' num2str(n) ') and ''q'' the order of the polynomials (' num2str(order) '), then ''n-1'' must be a multiple of ''q-1'' ']);
end

%Initialize the matrices
A=zeros((order+1)*(n-1)/(order-1));
b=zeros([size(A,1) 1]);


XX=ones(size(X));
II=ones(size(W));
DX=[ ones(size(W)) zeros(size(W)) ];
for i=1:order
    XX=[X.^i XX ];
    II=[W.^i II ];
    
    if(i>1)
        DX=[i*W.^(i-1) DX] ;
    end
end 

p=1;%polynomial counter
line=1; %count the lines writen in the matrix
j=1;%count intervals
for i=1:n
       
    if(i==1)
        %dP1_dx(xi_)=0\n'
        if( W(1)>=X(1))
            error('the first interval must be less than the first value of X (usually 0)')
        end
        %fprintf('dP1_dx=0\n')
        A(line,1:order+1)=DX(j,:);
        j=j+1;
        line=line+1;
    end
    %Pp(xi)=yi
    %fprintf(['P' num2str(p) '(x' num2str(i) ')=y'  num2str(i) '\n' ])
    A(line,(p-1)*(order +1)+1:(p)*(order +1))=XX(i,:);
    b(line)=Y(i);
    line=line+1;
    if mod(i,order-1)==0 && i<n-1
        %Pp(inter_j)=Pp+1(inter_j)
        %fprintf(['P' num2str(p) '(x)=P' num2str(p+1) '(x)\n'])
        A(line,(p-1)*(order +1)+1:p*(order +1))=II(j,:);
        A(line,(p)*(order +1)+1:(p+1)*(order +1))=-II(j,:);
        line=line +1;
       %dPp_dx(inter_j)=dPp+1_dx(inter_j)
       %fprintf(['dP' num2str(p) '_dx(x)=dP' num2str(p+1) '_dx(x)\n'])
        A(line,(p-1)*(order +1)+1:p*(order +1))=DX(j,:);
        A(line,(p)*(order +1)+1:(p+1)*(order +1))=-DX(j,:);
        line=line+1;
        
        j=j+1;
        
        p=p+1;
    end
    
end

a=reshape(inv(A)*b,order+1,p)';
end
