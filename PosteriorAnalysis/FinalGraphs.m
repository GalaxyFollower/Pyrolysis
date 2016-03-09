function FinalGraphs(name,kinName,savefig)
global ylabels FS LEL yst
load(name)


%Completelist=listOfGoodOnes;
Ti=cell2mat(Completelist(:,2));
ri=cell2mat(Completelist(:,3));

nvteri=cell2mat(Completelist(:,4));


Tinf=unique(Ti);
r=unique(ri);
nvter=unique(nvteri);
LEL=0.0783; % vol/vol
yst=0.232;

linestyle={'-';':';'--';'-.'};
colors=get(gca,'colororder');
colors=[0 0 0];
linestyle={'-'};
Pinf=1E5;%Pa
Rgas=8.314;%J/(mol*K)
Wgas=44;%g/mol

%ylabels{1}=['\textsf{$\mathsf{\frac{C_{Dust}}{y_{gas}}}$ \mbox{ } \mbox{ } \mbox{ }}'
%            '\textsf{$\left[\mathsf{\frac{g}{m^3}}\right]$\mbox{ }\mbox{ }}         '];
ylabels{1}='\begin{tabular}{@{}c@{}}$\mathsf{~\frac{C_{dust}}{y_{gas}}~}$~\\${}_{_{\mathsf{\left[\frac{g}{m^3}\right]}}}$~\end{tabular}';
ylabels{1}='';
ylabels{2}=['\textsf{$\mathsf{C_{dust}^{LEL}}$}'];
ylabels{3}=['\textsf{$\mathsf{C_{dust}^{y_{st}}}$}'];







if savefig
    FS.legend=14;
    FS.y=36;
    FS.x=24;
    FS.ax=20;
    FS.aux=20;
    FS.lw=3;
    FS.lbl=16;
else
    FS.legend=9;
    FS.y=16;
    FS.x=11;
    FS.ax=10;
    FS.aux=11;
    FS.lw=1;
    FS.lbl=8;
end
if savefig
    position=[0.1,0.01,0.56,0.88];
else
    position=[0,0,1,1];
end
%%
plotname=['Concentration for LEL vs time constant r'];
fignum=1;
fig=getFigureHdl([plotname ', part ' num2str(fignum)]);
figure(fig)
clf;

for i=1:length(r)

    if savefig
        clf
        hdl=axes;
    else
        try
            hdl=subplot(2,3,i-(fignum-1)*6);
        catch err
           set(fig,'name' ,[plotname ', part ' num2str(fignum)]);
           fignum=fignum+1;
           fig=getFigureHdl([plotname ', part ' num2str(fignum)]);
           hdl=subplot(2,3,i-(fignum-1)*6);
        end
    end
    set(fig, 'Units', 'normalized', 'Position', position);
        
    j=find(cell2mat(Completelist(:,3))==r(i) & ...
           cell2mat(Completelist(:,4))==1 &...
           cellfun(@(x)strcmp(x,kinName),Completelist(:,1)) &...
           cell2mat(Completelist(:,2))~=1050 & ...
           cell2mat(Completelist(:,2))~=1150  ...
       );
    mint=1;
    maxt=0;
    series=[];
    p=0;
    for k=1:length(j)
        x=yield{j(k)};
        Cdust=Pinf*Wgas./x(:,4)/Rgas/Ti(j(k));
        t=time{j(k)};
        nseries=find(Tinf==Ti(j(k)));
        mint=min([mint t(find(Cdust<1E9,1))]);
        maxt=max([maxt max(time{j(k)})]);
        
        %definition of series labeling

        ideltat=[find(Cdust<1E9,1) length(Cdust)];
        deltalabel=0.04*range(log(Cdust(ideltat)));
        try
            [~,m]=min(abs(log(t)-(log(t(ideltat(2))^0.6*t(ideltat(1))^0.4))));            
            [~,minf]=min(abs(log(Cdust)-log(Cdust(m))-...
            deltalabel)); 
            [~,msup]=min(abs(log(Cdust)-log(Cdust(m))+...
            deltalabel));
            series{p+1}.y=(Cdust(m-1:m+1));
            series{p+1}.x=(t(m-1:m+1));
            series{p+1}.name=[ '$' num2str(Ti(j(k))) '~\mathrm{K}$' ];
        catch err
            
            
        end

       % Cdust(minf:msup)=nan;
        
        if max(t)>mint && min(Cdust)<1E9
            if p==0

                [ax,hlines] = plotyyy(t,Cdust,[0],[0],[0],[0],ylabels);
                delete(hlines(2:3));
                set(hlines(1),'color', colors(mod(nseries,size(colors,1))+1,:),...
                    'LineStyle',linestyle{mod(nseries,length(linestyle))+1},...
                    'Linewidth',FS.lw);

            else
             loglog(t,Cdust,...
                 'color', colors(mod(nseries,size(colors,1))+1,:),...
                 'LineStyle',linestyle{mod(nseries,length(linestyle))+1},...
                 'Linewidth',FS.lw);
            end
            p=p+1;

            
        end
        hold on
      
    end
    if savefig
        graphtitle=[kinName 'Dp' num2str(2*r(i)*1E6) 'microm'];
    else
        graphtitle=['D_p=' num2str(2*r(i)*1E6) ' ' char(181) 'm'];
    end
    try
        makeItPretty(ax,LEL,graphtitle,series,[mint maxt],savefig,false)
    catch err
       
    end
end


%%
plotname=['Concentration for LEL vs time constant T'];
fignum=1;
fig=getFigureHdl([plotname ', part ' num2str(fignum)]);
figure(fig)
clf;

for i=1:length(Tinf)

    if savefig
        clf
        hdl=axes;
    else
        try
            hdl=subplot(2,3,i-(fignum-1)*6);
        catch err
           set(fig,'name' ,[plotname ', part ' num2str(fignum)]);
           fignum=fignum+1;
           fig=getFigureHdl([plotname ', part ' num2str(fignum)]);
           hdl=subplot(2,3,i-(fignum-1)*6);
        end
    end
    set(fig, 'Units', 'normalized', 'Position', position);
        
    j=find(cell2mat(Completelist(:,2))==Tinf(i) & ...
           cell2mat(Completelist(:,4))==1 &...
           cellfun(@(x)strcmp(x,kinName),Completelist(:,1)));
    mint=1;
    maxt=0;
    series=[];
    p=0;
    for k=1:length(j)
        x=yield{j(k)};
        Cdust=Pinf*Wgas./x(:,4)/Rgas/Ti(j(k));
        t=time{j(k)};
        nseries=find(r==ri(j(k)));
        mint=min([mint t(find(Cdust<1E9,1))]);
        maxt=max([maxt max(time{j(k)})]);
        
        %definition of series labeling

        ideltat=[find(Cdust<1E9,1) length(Cdust)];
        deltalabel=0.04*range(log(Cdust(ideltat)));
        try
            [~,m]=min(abs(log(t)-(log(t(ideltat(2))^0.2*t(ideltat(1))^0.8))));            
            [~,minf]=min(abs(log(Cdust)-log(Cdust(m))-...
            deltalabel)); 
            [~,msup]=min(abs(log(Cdust)-log(Cdust(m))+...
            deltalabel));
            series{p+1}.y=(Cdust(m-1:m+1));
            series{p+1}.x=(t(m-1:m+1));
            if savefig
                series{p+1}.name=[ '$' num2str(ri(j(k))*2E6) '~\mathrm{\upmu m}$' ];
            else
                series{p+1}.name=[ '$' num2str(ri(j(k))*2E6) '~\mathrm{\mu m}$' ];
            end
        catch err
            
            
        end

       % Cdust(minf:msup)=nan;
        
        if max(t)>mint && min(Cdust)<1E9
            if p==0

                [ax,hlines] = plotyyy(t,Cdust,[0],[0],[0],[0],ylabels);
                delete(hlines(2:3));
                set(hlines(1),'color', colors(mod(nseries,size(colors,1))+1,:),...
                    'LineStyle',linestyle{mod(nseries,length(linestyle))+1},...
                    'Linewidth',FS.lw);

            else
             loglog(t,Cdust,...
                 'color', colors(mod(nseries,size(colors,1))+1,:),...
                 'LineStyle',linestyle{mod(nseries,length(linestyle))+1},...
                 'Linewidth',FS.lw);
            end
            p=p+1;

            
        end
        hold on
      
    end
    if savefig
        graphtitle=[kinName 'Tinf' num2str(Tinf(i)) 'K'];
    else
        graphtitle=['Tinf=' num2str(Tinf(i)) ' K'];
    end
    try
        makeItPretty(ax,LEL,graphtitle,series,[mint maxt],savefig,false)
    catch err
       
    end
end



    
%%
plotname=['Concentration for LEL vs time variation of nvert'...
    '; D_p= ' num2str(30) char(181) 'm' ];
fignum=1;
fig=getFigureHdl([plotname ', part ' num2str(fignum)]);
figure(fig)
clf;
linestyle={'-.';'-';':';'--'};
colors=[0 0 0;0 0 0;0 0 0; 0 0 0; 0.4 0.4 0.4; 0.4 0.4 0.4; 0.4 0.4 0.4];
for i=1:length(Tinf)-4
    if savefig
        clf
        hdl=axes;
    else
        try
            hdl=subplot(2,3,i-(fignum-1)*6);
        catch err
           set(fig,'name' ,[plotname ', part ' num2str(fignum)]);
           fignum=fignum+1;
           fig=getFigureHdl([plotname ', part ' num2str(fignum)]);
           figure(fig)
           hdl=subplot(2,3,i-(fignum-1)*6);
        end
    end
    set (fig, 'Units', 'normalized', 'Position', position);

    j=find(cell2mat(Completelist(:,3))>=1.499e-5 & ...
        cell2mat(Completelist(:,3))<=1.501e-5  & ...
        cell2mat(Completelist(:,2))==Tinf(i) & ...
        cellfun(@(x)strcmp(x,kinName),Completelist(:,1)));

    maxt=0;
 
    series=[];
    mint=1;
    p=0;

    for k=1:length(j)
        x=yield{j(k)};
        Cdust=Pinf*Wgas./x(:,4)/Rgas/Tinf(i);
        t=time{j(k)};
        mint=min([mint t(find(Cdust<1E9,1))]);
        maxt=max([maxt max(time{j(k)})]);
        nseries=find(nvter==nvteri(j(k)));
        %definition of series labeling
        if(nvteri(j(k))==1)
           series{p+1}.name='$\mathrm{Re}(\sigma_\mathrm{p})$';
           seriesnames{p+1}=series{p+1}.name;
        else
            series{p+1}.name=['$\mathrm{Re}(' num2str(nvteri(j(k))) '\sigma_\mathrm{p})$'];
            seriesnames{p+1}=series{p+1}.name;
        end
        ideltat=[find(Cdust<1E9,1) length(Cdust)];
        deltalabel=0.03*range(log(Cdust(ideltat)));
        
         [~,m]=min(abs(log(t)-(log(t(ideltat(2))^0.6*t(ideltat(1))^0.4))));        
        [~,minf]=min(abs(log(Cdust)-log(Cdust(m))-...
            deltalabel)); 
        [~,msup]=min(abs(log(Cdust)-log(Cdust(m))+...
            deltalabel));
        series{p+1}.y=(Cdust(m-1:m+1));
        series{p+1}.x=(t(m-1:m+1));
       % Cdust(minf:msup)=nan;
        
        
        if max(t)>mint && min(Cdust)<1E9
            if p==0

                [ax,hlines] = plotyyy(t,Cdust,[0],[0],[0],[0],ylabels);
                delete(hlines(2:3));
                set(hlines(1),'color', colors(mod(nseries,size(colors,1))+1,:),...
                    'LineStyle',linestyle{mod(nseries,length(linestyle))+1},...
                    'Linewidth',FS.lw);

            else
             loglog(t,Cdust,...
                 'color', colors(mod(nseries,size(colors,1))+1,:),...
                 'LineStyle',linestyle{mod(nseries,length(linestyle))+1},...
                 'Linewidth',FS.lw);
            end
            p=p+1;

            
        end
        hold on             
    end
    if savefig
        graphtitle=[kinName 'Tinf' num2str(Tinf(i)) 'KnvterDp100microm' ];
    else
        graphtitle=['T_{\infty}=' num2str(Tinf(i)) 'K' ];
    end
    try
        makeItPretty(ax,LEL,graphtitle,seriesnames,[mint maxt],savefig,true)
    catch
    end
end

%%
plotname=['Concentration for LEL vs time variation of nvert'...
    '; Tinf= ' num2str(900) 'K' ];
fignum=1;
fig=getFigureHdl([plotname ', part ' num2str(fignum)]);
figure(fig)
set (fig, 'Units', 'normalized', 'Position', [0,0,0.5,1]);
clf;

for i=1:length(r)
    if savefig
        clf
        hdl=axes;
    else
        try
            hdl=subplot(2,3,i-(fignum-1)*6);
        catch err
           set(fig,'name' ,[plotname ', part ' num2str(fignum)]);
           fignum=fignum+1;
           fig=getFigureHdl([plotname ', part ' num2str(fignum)]);
           figure(fig)
           hdl=subplot(2,3,i-(fignum-1)*6);
        end
    end
    set (fig, 'Units', 'normalized', 'Position', position);

    j=find(cell2mat(Completelist(:,3))==r(i) & ...
           cell2mat(Completelist(:,2))==950 &...
           cellfun(@(x)strcmp(x,kinName),Completelist(:,1)));
    mint=1;
 
    maxt=0;
    series=[];
    seriesnames=[];
    p=0;
         
    for k=1:length(j)
        x=yield{j(k)};
        Cdust=Pinf*Wgas./x(:,4)/Rgas/Ti(j(k));
        t=time{j(k)};
        nseries=find(nvter==nvteri(j(k)));
        mint=min([mint t(find(Cdust<1E9,1))]);
        maxt=max([maxt max(time{j(k)})]);
        %definition of series labeling
        if(nvteri(j(k))==1)
            series{p+1}.name='$\mathrm{Re}(\sigma_\mathrm{p})$';
            seriesnames{p+1}=series{p+1}.name
        else
            series{p+1}.name=['$\mathrm{Re}(' num2str(nvteri(j(k))) '\sigma_\mathrm{p})$'];
            seriesnames{p+1}=series{p+1}.name
        end
        ideltat=[find(Cdust<1E9,1) length(Cdust)];
        deltalabel=0.03*range(log(Cdust(ideltat)));
        
        [~,m]=min(abs(log(t)-(log(t(ideltat(2))^0.6*t(ideltat(1))^0.4))));          
        [~,minf]=min(abs(log(Cdust)-log(Cdust(m))-...
            deltalabel)); 
        [~,msup]=min(abs(log(Cdust)-log(Cdust(m))+...
            deltalabel));
        series{p+1}.y=(Cdust(m-1:m+1));
        series{p+1}.x=(t(m-1:m+1));
        %Cdust(minf:msup)=nan;
        if max(t)>mint && min(Cdust)<1E9
            if p==0

                [ax,hlines] = plotyyy(t,Cdust,[0],[0],[0],[0],ylabels);
                delete(hlines(2:3));
                set(hlines(1),'color', colors(mod(nseries,size(colors,1))+1,:),...
                    'LineStyle',linestyle{mod(nseries,length(linestyle))+1},...
                    'Linewidth',FS.lw);

            else
             loglog(t,Cdust,...
                 'color', colors(mod(nseries,size(colors,1))+1,:),...
                 'LineStyle',linestyle{mod(nseries,length(linestyle))+1},...
                 'Linewidth',FS.lw);
            end
            p=p+1;
            
        end

        hold on
        
        
        
    end
    

    
    if savefig
        graphtitle=[kinName 'Dp' num2str(2*r(i)*1E6) 'micromNvterTinf900K'];
    else
        graphtitle=['D_p=' num2str(2*r(i)*1E6) ' ' char(181) 'm'];
    end
    try
    
        makeItPretty(ax,LEL,graphtitle,seriesnames,[mint maxt],savefig,true)

    catch
    end
    
end
end

function makeItPretty(ax,LEL,graphtitle,series,limt,savefig,makelegend)
    
    global ylabels FS yst 
    warning off
    fig=get(ax(1),'parent');
    if savefig
        position=[0.07 0.12 0.78 0.8
                  0.07 0.12 0.78 0.8
                  0.07 0.12 0.86 0.8];
    else
        position(1,:)=get(ax(1),'position');
        position(2,:)=get(ax(2),'position');
        position(3,:)=get(ax(3),'position');
    end
    xlim(ax(1),limt);        
    ylim(ax(1),[1E3 1e9]);

        
    ylim(ax(2),[1E3 1e9]*LEL);
    ylim(ax(3),[1E3 1e9]*yst);
    

    xorder=floor(log10(limt(1))):ceil(log10(limt(2)));
    yorder=2:8;

    ylabel(ax(1), ylabels{1},'interpreter','latex',...
        'HorizontalAlignment','center', 'VerticalAlignment','middle',...
        'units','normalized','position',[-.13 0.5 0],'rotation',0, ...
        'fontsize',FS.y,'fontweight','bold')
    
    
    xlabel('$\mathsf{Time~[~s~]}$','fontsize', FS.x,'interpreter',...
        'latex','units','normalized','position',[0.5 -0.08 0] ) 
    set(ax,'fontsize',FS.ax)
    grid on
    set(ax(1), 'GridLineStyle','-','MinorGridLineStyle', '-',...
    'YMinorTick','on',...
    'YScale','log',...
    'YTick',10.^yorder,...
    'XMinorTick','on',...
    'XScale','log',...+
    'XTick',10.^(xorder),...
    'TickLabelInterpreter','none',...    'yticklabel',yticks,...'xticklabel',xticks,...
    'linewidth',1.5,...
    'MinorGridAlpha',0.05,...
    'GridAlpha',0.25,...
    'MinorGridColor',[0 0 0],...
    'GridColor',[0 0 0],...
    'box','on',...
    'position',position(1,:));
    
    set(ax(2:3),'YMinorTick','on',...
    'YScale','log','linewidth',1.5,...
    'YTick',10.^(2:10),'Ticklabelinterpreter','none');
    %set(ax(2),'Ycolor',[0 0.2 0.6])
    %set(ax(3),'Ycolor',[0.6 0 0])
    ylab=get(ax,'ylabel');
    set(ax(2),'position',position(2,:))
    set(ax(3),'position',position(3,:))
    for a=2:3
        yl=get(ax(a),'ylim');
        xl=get(ax(a),'xlim');
        set(ylab{a},'position',[xl(2) 3.5*yl(2) 0]);
    end

   
    
    hold off
    
	yauxlabels=get(ax(2:3),'ylabel');
    for a=1:2
        set(yauxlabels{a},'fontsize',FS.aux);
    end
    set(ax,'box','off','fontsize',FS.ax)
    
    
    if(savefig)
        set(fig,'children',ax(3:-1:1))
        if(makelegend)
            legend(ax(1),series,'fontsize',FS.legend,'interpreter','latex',...
            'location','best')
        else
            labelthem(series,ax(1),FS)   
        end
        done=false;
        while ~done
            try
                figname=['Figures\' strrep(strrep(graphtitle,'Wangenaar','Wagenaar'),' ','')];
%                 latexfigure(fig,figname,'png');
%                 FileWait( [figname,'.png'] );
%                 latexfigure(fig,figname,'eps');
%                 FileWait( [figname,'.eps'] );
 

                 latexfigure(fig,figname,'pdf');
                 FileWait( [figname,'.pdf'] );
                done=true;
            catch err
                fprintf('\n\nI could not save %s files, I will try again\n\n',...
                    figname);
            end
            
        end
    else
        title(graphtitle); 
         if(makelegend)
            legend(ax(1),series,'fontsize',FS.legend,'interpreter','latex',...
            'location','best')
        else
            labelthem(series,ax(1),FS)   
        end
    end
    
    warning off;
    
    
end
    
function labelthem(series,hdl,FS)
    K=1;
               
    for i=1:length(series)
        orientation=180/pi*atan(K*(log(series{i}.y(3)/series{i}.y(1)))/...
            (log(series{i}.x(3)/series{i}.x(1))));
        text(series{i}.x(2),series{i}.y(2),[series{i}.name],'Rotation',...
            orientation,'horizontalalignment','center','fontsize',FS.lbl,...
            'verticalalignment','top','interpreter','latex','parent',hdl)
    end
end

% Busy waits for a file to finish being created. This is necessary because
% on some platforms the file isn't available immediately after performing a
% print.
  function FileWait(filename)
    counter = 0;
    WAIT=true;
    while WAIT
        if ~exist(filename,'file')
            pause(0.05);
            counter = counter + 1;
        else 
            WAIT=false;
        end
        if counter>1200
            error('The file can not be saved\n');
        end
    end
  end
