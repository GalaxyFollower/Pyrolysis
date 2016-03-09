function plotParticulesData(fileName)

load(fileName)


subindex='BCTG';

s=min(5,nv);%number of lines to draw
cord=[CM;r(end)];
plotem=param.r-exp(linspace(min(log(cord)),max(log(cord)),s));
j=zeros(5,1);

for i=1:s
   
   [~, k]= sort(abs(plotem(i)-CM));
   OK=false;
   l=1;
   while(~OK)
    OK=isempty(find(j==k(l), 1));
 
    j(s-i+1)=k(l);
    l=l+1;
   end
    
   
end


name=cellfun(@(x)strrep(x,' ',''),mat2cell([repmat('r=',[s 1]) ...
    num2str(round(cord(j)*2E7)/20) repmat('µm',[s 1])],repmat([1],[s 1])),...
    'UniformOutput',false);


   %------------------gases------------------------- 
   
    k=find(rho(:,1,gou)>1E-10,1);
    minrhog=10^ceil(log10( rho(k,1,gou)));
    minrhos=10^ceil(log10( rho(k,1,cha)));
    minj=find(rho(:,1,gou)>0,1);
    mint=time(k);
    maxt=10^(ceil(log10(max(time))*10)/10);
    orders=ceil(log10(mint)):ceil(log10(maxt)) ;

    maxrhog= 10^(log10(ceil(max(max(rho(:,:,gou)))) ));
    maxrhos= 10^(log10(ceil(max(max(1.5*rho(:,:,bio))))));
    rhogorders=ceil(log10(minrhog)):ceil(log10(maxrhog)) ;
    rhosorders=ceil(log10(minrhos)):ceil(log10(maxrhos)) ;
    
    xlab=cellfun(@(x)strrep(x,' ',''),mat2cell([repmat('$10^{',[length(orders) 1]) ...
    num2str(orders') repmat('}$',[length(orders) 1])],ones([length(orders) 1])),...
    'UniformOutput',false);
    ylab=cellfun(@(x)strrep(x,' ',''),mat2cell([repmat('$10^{',[length(rhogorders) 1]) ...
    num2str(rhogorders') repmat('}$',[length(rhogorders) 1])],ones([length(rhogorders) 1])),...
    'UniformOutput',false);

    fig=getFigureHdl(['gases vs time' ] );
    figure(fig)
    clf
    set(fig,'units','normalized') 
    set(fig,'position',[0.38    0.08    0.6000    0.82]);
    
    ax(1)=subplot(4,2,1);    
    ax(2)=subplot(4,2,2);
    ax(3)=subplot(4,2,3);    
    ax(4)=subplot(4,2,4);
    
    %gases
    set(ax(1),'LineStyleOrder','-|--|:|-.|-','ColorOrder',[repmat([0 0 0],...
        [s 1]); repmat([0.4 0.4 0.4],[s 1])] ,'fontsize',20, 'box','on',...
        'xgrid','on','ygrid','on');
    hold(ax(1),'on')
    plot(ax(1),time,[rho(:,j,gou)  rho(:,j,gas) ],'linewidth',2)    
    hold(ax(1),'off')
    
    set(ax(2),'ColorOrder',[0 0 0], 'LineStyleOrder','-|--|:|-.|-',...
        'fontsize',20, 'box','on','xtick',10.^(orders),'xticklabel',xlab,...
        'ytick',10.^(rhogorders),'yticklabel',ylab,'xgrid','on','ygrid',...
        'on','yscale','log','xscale','log','minorgridlinestyle','-',...
        'minorgridalpha',0.05);
    
    hold(ax(2),'on')
    loglog(ax(2),time,rho(:,j,gou),'linewidth',2) 
    set(ax(2),'ColorOrder',[0.4 0.4 0.4], 'LineStyleOrder','-|--|:|-.|-')
    loglog(ax(2),time,rho(:,j,gas),'linewidth',2) 
    hold(ax(2),'off')
    %solids
    set(ax(3),'LineStyleOrder','-|--|:|-.','ColorOrder',[repmat([0 0 0],...
        [s 1]); repmat([0.4 0.4 0.4],[s 1])] ,'fontsize',20, 'box','on',...
        'xgrid','on','ygrid','on');
    hold(ax(3),'on')
    plot(ax(3),time,[rho(:,j,bio)  rho(:,j,cha) ],'linewidth',2)    
    hold(ax(3),'off')
    xlab=get(ax(3),'xticklabel');
    xlab=cellfun(@(x)strrep(x,' ',''),mat2cell([repmat(' ',...
        [length(xlab) 1])],ones([length(xlab) 1])),...
    'UniformOutput',false);
    set(ax(3),'xticklabel',xlab);
    
    xlab=cellfun(@(x)strrep(x,' ',''),mat2cell([repmat(' ',...
        [length(orders) 1])],ones([length(orders) 1])),...
    'UniformOutput',false);
    ylab=cellfun(@(x)strrep(x,' ',''),mat2cell([repmat('$10^{',[length(rhosorders) 1]) ...
    num2str(rhosorders') repmat('}$',[length(rhosorders) 1])],ones([length(rhosorders) 1])),...
    'UniformOutput',false);

    set(ax(4),'ColorOrder',[0 0 0], 'LineStyleOrder','-|--|:|-.|-',...
        'fontsize',20, 'box','on','xtick',10.^(orders),'xticklabel',xlab,...
        'ytick',10.^(rhosorders),'yticklabel',ylab,'xgrid','on','ygrid',...
        'on','yscale','log','xscale','log','minorgridlinestyle','-',...
        'minorgridalpha',0.05);
    
    hold(ax(4),'on')
    loglog(ax(4),time,rho(:,j,bio),'linewidth',2) 
    set(ax(4),'ColorOrder',[0.4 0.4 0.4], 'LineStyleOrder','-|--|:|-.|-')
    loglog(ax(4),time,rho(:,j,cha),'linewidth',2) 
    hold(ax(4),'off')
    
    
    ylabel(ax(1),'$\mathsf{\rho_{_{Tar}} ,\rho_{_{Gas}}~{\scriptstyle\left[~\frac{kg}{m^3}~\right]}}$',...
        'fontsize',30,'interpreter','latex')
    ylabel(ax(3),'$\mathsf{\rho_{_{Biomass}} ,\rho_{_{char}}~{\scriptstyle\left[~\frac{kg}{m^3}~\right]}}$',...
        'fontsize',30,'interpreter','latex')
    
    ypos= get(get(ax(1),'ylabel'),'position');
    ypos([1 3])=[-0.1 -0.5];
    
    xlabel(ax(2),['$\mathsf{time~{\scriptstyle\left[~s~\right]}}$'],'fontsize',30,...
        'interpreter','latex')
    xlabel(ax(1),['$\mathsf{time~{\scriptstyle\left[~s~\right]}}$'],'fontsize',30,...
        'interpreter','latex')
    
    set(ax(1),'position',[0.1 0.09 0.41 0.43])
    set(ax(2),'position',[0.57 0.09 0.41 0.43])
    set(ax(3),'position',[0.1 0.55 0.41 0.43])
    set(ax(4),'position',[0.57 0.55 0.41 0.43])
    
    
    set(get(ax(1),'xlabel'),'units','normalized','position',[ 0.5 -0.07 1])  
    set(get(ax(2),'xlabel'),'units','normalized','position',[ 0.5 -0.07 1])  
    set(get(ax(1),'ylabel'),'units','normalized','position',[ -0.11 0.5 1])
    set(get(ax(3),'ylabel'),'units','normalized','position',[ -0.11 0.5 1])
    
    xlim(ax(2),[mint maxt])
    ylim(ax(2),[minrhog maxrhog])
    ylim(ax(4),[minrhos maxrhos])
    ylim(ax(3),[0 1500])
    xlim(ax(4),[mint maxt])
    %ylim(ax(2),[minrho maxrho])
    %latexfigure(fig,['rho' subindex(i) 't_FontTinf1000r50micronvter1'],'pdf');


%---------------yield -----------------------

    fig=getFigureHdl('Yield profile');
    
    clf
    set(fig,'units','normalized')
    set(fig,'position',[0.0448    0.5074    0.3000    0.4139]);
    ax=axes;
    set(ax,'LineStyleOrder','-|--|:|-.','ColorOrder',[0 0 0],...
      'fontsize',20, 'box','on','xgrid','on','ygrid','on');
  
    hold on
    plot(time, yield,'linewidth',2);
    hold off
    
    ylabel(['$\mathsf{\rho_{_' num2str(i) '}~(\si{\kilogram\per\cubic\metre})}$'],'fontsize',30,'interpreter','latex')
    ypos= get(get(ax,'ylabel'),'position');
    ypos([1 3])=[-0.0827 -1];
    xlabel(['$\mathsf{time~(\si{\second})}$'],'fontsize',30,'interpreter','latex')
    set(ax,'position',[.18 .18 .81 .75])
    set(get(ax,'xlabel'),'units','normalized','position',[ 0.5 -0.08]);  
    set(get(ax,'ylabel'),'units','normalized','position',[ -0.1 0.5 0])
        

%---------plot temperature-------------------------
[minT minj]=min(T);
minj=max(minj);
mint=10^floor(log10(time(minj))-1);
maxt=10^(ceil(log10(max(time))*10)/10);
orders=ceil(log10(mint)):ceil(log10(maxt)) ;


fig=getFigureHdl('State Variables');
figure(fig)
clf
  set(fig,'units','normalized')
set(fig,'position',[0.0448   0 0.5   1]);
ax1=subplot(2,1,1);
ylab=300:100:round(max(max(T))/100+1)*100;
xlab=cellfun(@(x)strrep(x,' ',''),mat2cell([repmat(' ',[s 1])],repmat([1],[s 1])),...
    'UniformOutput',false);
set(ax1,'LineStyleOrder','-|--|:|-.','ColorOrder',[0 0 0],...
      'xscale','log','fontsize',20,  'ytick',ylab,'xtick',10.^(orders),'xticklabel',xlab,...
      'box','on','xgrid','on','ygrid','on','minorgridlinestyle','-',...
        'minorgridalpha',0.05);
  hold on
semilogx(time,T(:,[j(1:s-1); nv+1] ),'linewidth',1.5)

  hold off
ylabel(['$\mathsf{Temperature~{\scriptstyle(K)}}$'],'fontsize',30,'units','normalized','position',[-0.055 0.5 0])
legend([name(1:s-1);['r=' num2str(round(r(nv+2)*1E7)/10) ' µm']],'Fontsize',16,'Fontsize',16,'location','best')
xlim([mint maxt])
ylim([300 round(max(max(T))/50+1)*50])

%---------plot Pressure-------------------------

minP=min(P(minj:end,max(j(1:s-1))));
maxP =max(max(P));
minP=floor(minP*10^(-floor(log10(minP))+1))/10^(-floor(log10(minP))+1);
maxP=ceil(maxP*10^(-floor(log10(maxP))))/10^(-floor(log10(maxP)));
porders=(floor(log10(minP)):ceil(log10(maxP)));

xlab=cellfun(@(x)strrep(x,' ',''),mat2cell([repmat('$10^{',[length(orders) 1]) ...
    num2str(orders') repmat('}$',[length(orders) 1])],ones([length(orders) 1])),...
    'UniformOutput',false);
ylab=cellfun(@(x)strrep(x,' ',''),mat2cell([repmat('$10^{',[length(porders) 1]) ...
    num2str(porders') repmat('}$',[length(porders) 1])],ones([length(porders) 1])),...
    'UniformOutput',false);

ax2=subplot(2,1,2);
set(ax2,'LineStyleOrder','-|--|:|-.','ColorOrder',[0 0 0],...
      'xscale','log','yscale','log','fontsize',20, 'xtick',10.^(orders),...
      'xticklabel',xlab,'ytick', 10.^porders,'yticklabels', ylab,...
      'box','on','xgrid','on','ygrid','on','minorgridlinestyle','-',...
      'minorgridalpha',0.05);
  hold on
loglog(time,P(:,j),'linewidth',1.5)
hold off

ylabel(['$\mathsf{Overpresure~{\scriptstyle(Pa)}}$'],'fontsize',30,'units','normalized','position',[-0.055 0.5 0])
xlabel(['$\mathsf{time~(s)}$'],'fontsize',30,'units','normalized','position',[0.5 -0.06])

legend(name,'Fontsize',16,'location','best')
xlim([mint maxt])

set(ax1,'position',[0.11 0.54 0.87 0.44])
set(ax2,'position',[0.11 0.08 0.87 0.44])
ylim(ax2,[minP maxP])



end

