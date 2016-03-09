figure
rhomean=200:10:1400;
T=1000;
n=[0.001 0.25 1 5 10];
param.vter_mult=n;
param.r=2E-6;

h_conv=zeros(length(rhomean),length(n));
for i=1:length(rhomean)
    
    param.rho_mean=rhomean(i);
    [h_conv(i,:) Re]=hconv(T,param);
    
end

plot(rhomean,h_conv)