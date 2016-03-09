function [y x]=rungekutta(fun,s,y0,h)
a=s(1);
b=s(2);

x(1)=a;
y(1,:)=y0';
i=1;
while(x<b)
    k1=fun(x(i),y(i,:)');
    k2=fun(x(i)+0.5*h,y(i,:)'+0.5*k1*h);
    k3=fun(x(i)+0.5*h,y(i,:)'+0.5*k2*h);
    k4=fun(x(i)+h,y(i,:)'+k3*h);
    y(i+1,:)=y(i,:)+1/6*(k1+2*k2+2*k3+k4)'*h;
    x(i+1)=x(i)+h;
    i=i+1;
end
%[min_difference, position] =min(abs(b-x)');
%y=y(position);
