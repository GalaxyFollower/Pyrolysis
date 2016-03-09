clc 
clear
load param

kinetics(1).name=['Chan'];
kinetics(1).A=[1.3e8 2e8 1.08e7 0];
kinetics(1).E=[140e3 133e3 121e3 0];

kinetics(2).name=['Wangenaar'];
kinetics(2).A=[1.11E11 9.28E9 3.05E7 0];
kinetics(2).E=[177E3 149e3 125e3 0];

kinetics(3).name=['DiBlasi'];
kinetics(3).A=[4.4E9 1.1E10 3.3E6 0];
kinetics(3).E=[153E3 148e3 112e3 0];

kinetics(4).name=['Font'];
kinetics(4).A=[1.52E7 5.85E6 2.98E3 0];
kinetics(4).E=[139E3 119e3 73e3 0];


dp=unique([1E-6;exp(log(1E-6):0.01:log(500E-6))';500E-6]);
T=(300:20:10000)';
Rgas=param.Rgas;
figure(1)
clf

[TT, DP]=meshgrid(T,dp);
[m,n]=size(TT);
Cp=polyval(param.Cp(1,:),TT);
t_int=param.rhob0*DP.^2.*Cp/(6*param.lambda(1));
t_ext= real(param.rhob0*Cp.*(TT-param.T0).*DP./(6*(hconv(TT,param,DP).*(TT-param.T0)+param.stefBoltz*param.emis*(TT.^4-param.T0^4))));
for i=1:4
    subplot(2, 2,i)
    kinetics(i).t_pyro=1./sum(repmat(kinetics(i).A,size(T)).*exp(-repmat(kinetics(i).E,size(T) )./(Rgas*repmat(T,[1 4]))),2);
    kinetics(i).t_tot=t_int+t_ext+repmat(kinetics(i).t_pyro',[m 1]);
    [C,h]=contourf(TT,DP,kinetics(i).t_tot*1000,'showtext','on','fill','off','levellist',sort([10.^(-2:7) 3*10.^(-2:7)])');
    ax=get(h,'parent');
    xlim([700 10000])
    set(ax,'xscale','log','yscale','log','ytick',[(1:9)*1E-6 (1:10)*1E-5],'yticklabel',[(1:9) (10:10:100)])
    
%     semilogy(T,kinetics(i).tpyro)
%     hold on
end

% ax=get(1,'children')
% set(ax,'xgrid','on','ygrid','on', 'ytick',10.^(-6:0))
% ylim([1E-6 1]
% xlim([750 5000])set(ax,'xscale','log','yscale','log','yticklabel',[(1:9)*1E-6 (1:10)*1E-5],'yticklabel',[(1:9) (10:100)])
% hold off
% 
% 
% 
