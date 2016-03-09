function FinalGraphsYield(name,kinName,savefig)
global ylabels FS LEL yst
load(name)


%Completelist=listOfGoodOnes;
Ti=cell2mat(Completelist(:,2));
ri=cell2mat(Completelist(:,3));

nvteri=cell2mat(Completelist(:,4));


Tinf=unique(Ti);
r=unique(ri);
nvter=unique(nvteri);
LEL=0.09; % vol/vol
yst=0.332;

linestyle={'-';':';'--';'-.'};
colors=get(gca,'colororder');
colors=[0 0 0];
linestyle={'-'};
Pinf=1E5;%Pa
Rgas=8.314;%J/(mol*K)
Wgas=44;%g/mol

%ylabels{1}=['\textsf{$\mathsf{\frac{C_{Dust}}{y_{gas}}}$ \mbox{ } \mbox{ } \mbox{ }}'
%            '\textsf{$\left[\mathsf{\frac{g}{m^3}}\right]$\mbox{ }\mbox{ }}         '];





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
    position=[1.1,0.25,0.5,0.55];
else
    position=[0,0,1,1];
end
%%
plotname=['Bio yield vs time constant r'];
fignum=1;
fig=getFigureHdl([plotname ', part ' num2str(fignum)]);
figure(fig)
clf;
 ylabels{1}='$\mathsf{Biomass~consumption~[~-~]}$';
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
           cell2mat(Completelist(:,2))~=1150  );
    mint=100;
    maxt=0;
    maxy=0;
    series=[];
    p=0;
    for k=1:length(j)
        x=yield{j(k)};        
        t=time{j(k)};
        nseries=find(Tinf==Ti(j(k)));
        

       l=find(x(:,1)<max(x(:,1))/2,1);
       if(~isempty(l))
            mint=min(mint,t(find(x(:,1)<0.99,1)));
            maxt=max(maxt,max(t));
       
            maxy=max(maxy,max(x(:,1)));
           series{p+1}.x(1)=t(l-10);
           series{p+1}.x(2)=t(l);
           series{p+1}.x(3)=t(l+10);
           
           series{p+1}.y(1)=x(l-10,1);
           series{p+1}.y(2)=x(l,1);
           series{p+1}.y(3)=x(l+10,1);
           
           series{p+1}.name=[ '$' num2str(Ti(j(k))) '~\mathrm{K}$' ];
           semilogx(t,1-x(:,1),...
             'color', colors(mod(nseries,size(colors,1))+1,:),...
             'LineStyle',linestyle{mod(nseries,length(linestyle))+1},...
             'Linewidth',FS.lw);
             p=p+1;
             hold on
       end



       
    end
    
    if savefig
        graphtitle=['yield\Bioyield' kinName 'Dp' num2str(2*r(i)*1E6) 'microm'];
    else
        graphtitle=['D_p=' num2str(2*r(i)*1E6) ' ' char(181) 'm'];
    end
    try
        makeItPretty(hdl,graphtitle,series,[mint maxt],ceil(maxy*10)/10 ,savefig,false)
    catch err
        err
    end
end



%%
plotname=['gas yield vs time constant r'];
fignum=1;
fig=getFigureHdl([plotname ', part ' num2str(fignum)]);
figure(fig)
clf;
ylabels{1}='$\mathsf{Gas~yield~[~-~]}$';
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
           cell2mat(Completelist(:,2))~=1150  );
    mint=100;
    maxt=0;
    maxy=0;
    series=[];
    p=0;
    for k=1:length(j)
        x=yield{j(k)};        
        t=time{j(k)};
        nseries=find(Tinf==Ti(j(k)));
        l=find(x(:,4)>max(x(:,4))/2,1);
       if(~isempty(l))
        mint=min(mint,t(find(x(:,1)<0.99,1)));
        maxt=max(maxt,max(t));
       
       maxy=max(maxy,max(x(:,4)));
       
           series{p+1}.y(1)=x(l-10,4);
           series{p+1}.y(2)=x(l,4);
           series{p+1}.y(3)=x(l+10,4);

           series{p+1}.x(1)=t(l-10);
           series{p+1}.x(2)=t(l);
           series{p+1}.x(3)=t(l+10);
           series{p+1}.name=[ '$' num2str(Ti(j(k))) '~\mathrm{K}$' ];
            semilogx(t,x(:,4),...
             'color', colors(mod(nseries,size(colors,1))+1,:),...
             'LineStyle',linestyle{mod(nseries,length(linestyle))+1},...
             'Linewidth',FS.lw);
        end
        hold on
        p=p+1;
    end
    
    if savefig
        graphtitle=['yield\Gasyield' kinName 'Dp' num2str(2*r(i)*1E6) 'microm'];
    else
        graphtitle=['D_p=' num2str(2*r(i)*1E6) ' ' char(181) 'm'];
    end
    try
        makeItPretty(hdl,graphtitle,series,[mint maxt],ceil(maxy*10)/10 ,savefig,false)
    catch err
    end
end


%%
plotname=['tar yield vs time constant r'];
fignum=1;
fig=getFigureHdl([plotname ', part ' num2str(fignum)]);
figure(fig)
clf;
ylabels{1}='$\mathsf{Tar~yield~[~-~]}$';
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
           cell2mat(Completelist(:,2))~=1150  );
    mint=100;
    maxt=0;
    maxy=0;
    series=[];
    p=0;
    for k=1:length(j)
        x=yield{j(k)};        
        t=time{j(k)};
        nseries=find(Tinf==Ti(j(k)));
    
        l=find(x(:,3)>max(x(:,3))/2,1);
       if(~isempty(l)) 
        mint=min(mint,t(find(x(:,1)<0.99,1)));
        maxt=max(maxt,max(t));
       
       maxy=max(maxy,max(x(:,3)));
       

       series{p+1}.y(1)=x(l-10,3);
       series{p+1}.y(2)=x(l,3);
       series{p+1}.y(3)=x(l+10,3);
       
       series{p+1}.x(1)=t(l-10);
       series{p+1}.x(2)=t(l);
       series{p+1}.x(3)=t(l+10);
       series{p+1}.name=[ '$' num2str(Ti(j(k))) '~\mathrm{K}$' ];
        semilogx(t,x(:,3),...
         'color', colors(mod(nseries,size(colors,1))+1,:),...
         'LineStyle',linestyle{mod(nseries,length(linestyle))+1},...
         'Linewidth',FS.lw);
       end
        hold on
        p=p+1;
    end
    
    if savefig
        graphtitle=['yield\Taryield' kinName 'Dp' num2str(2*r(i)*1E6) 'microm'];
    else
        graphtitle=['D_p=' num2str(2*r(i)*1E6) ' ' char(181) 'm'];
    end
    try
        makeItPretty(hdl,graphtitle,series,[mint maxt],ceil(maxy*10)/10 ,savefig,false)
    catch err
    end
end

%%
plotname=['char yield vs time constant r'];
fignum=1;
fig=getFigureHdl([plotname ', part ' num2str(fignum)]);
figure(fig)
clf;
ylabels{1}='$\mathsf{Char~yield~[~-~]}$';
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
        
    j=find(cell2mat(Completelist(:,3))>=0.95*r(i) & ...
           cell2mat(Completelist(:,3))<=1.05*r(i) & ...
           cell2mat(Completelist(:,4))==1 &...
           cellfun(@(x)strcmp(x,kinName),Completelist(:,1)) &...
           cell2mat(Completelist(:,2))~=1050 & ...
           cell2mat(Completelist(:,2))~=1150  );
    mint=100;
    maxt=0;
    maxy=0;
    series=[];
    p=0;
    for k=1:length(j)
        x=yield{j(k)};        
        t=time{j(k)};
        nseries=find(Tinf==Ti(j(k)));
        
        l=find(x(:,3)>max(x(:,3))/2,1);
         
       if(~isempty(l)) 
            mint=min(mint,t(find(x(:,1)<0.99,1)));
       
        
        maxt=max(maxt,max(t));
       
       maxy=max(maxy,max(x(:,2)));
       l=find(x(:,2)>max(x(:,2))/2,1);
       series{p+1}.y(1)=x(l-10,2);
       series{p+1}.y(2)=x(l,2);
       series{p+1}.y(3)=x(l+10,2);
       
       series{p+1}.x(1)=t(l-10);
       series{p+1}.x(2)=t(l);
       series{p+1}.x(3)=t(l+10);
       series{p+1}.name=[ '$' num2str(Ti(j(k))) '~\mathrm{K}$' ];
        semilogx(t,x(:,2),...
         'color', colors(mod(nseries,size(colors,1))+1,:),...
         'LineStyle',linestyle{mod(nseries,length(linestyle))+1},...
         'Linewidth',FS.lw);
       end
        hold on
        p=p+1;
    end
    
    if savefig
        graphtitle=['yield\Charyield' kinName 'Dp' num2str(2*r(i)*1E6) 'microm'];
    else
        graphtitle=['D_p=' num2str(2*r(i)*1E6) ' ' char(181) 'm'];
    end
    try
        makeItPretty(hdl,graphtitle,series,[mint maxt],ceil(maxy*10)/10 ,savefig,false)
    catch err
    end
end

%%
plotname=['Porosity vs time constant r'];
fignum=1;
fig=getFigureHdl([plotname ', part ' num2str(fignum)]);
figure(fig)
clf;
 ylabels{1}='$\mathsf{Average~porosity~[~-~]}$';
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
           cell2mat(Completelist(:,2))~=1150  );
    mint=100;
    maxt=0;
    maxy=0;
    series=[];
    p=0;
    for k=1:length(j)
        x=yield{j(k)};        
        t=time{j(k)};
        nseries=find(Tinf==Ti(j(k)));
        load(Completelist{j(k),5},'param');
       rhog0=param.porosity0*param.P0*param.W(4)./(param.Rgas*param.T0);
       porosity=1-(1-param.porosity0)/param.rhob0*(param.rhob0+rhog0)*(sum(x(:,1:2),2));
        l=find(porosity(:,1)>max(porosity(:,1))/2,1);             
       if(~isempty(l) && length(x)>10 ) 
        mint=min(mint,t(find(x(:,1)<0.99,1)));
        maxt=max(maxt,max(t));
       
       maxy=max(maxy,max(x(:,1)));
       
       series{p+1}.y(1)=porosity(l-10,1);
       series{p+1}.y(2)=porosity(l,1);
       series{p+1}.y(3)=porosity(l+10,1);
       
       series{p+1}.x(1)=t(l-10);
       series{p+1}.x(2)=t(l);
       series{p+1}.x(3)=t(l+10);
       series{p+1}.name=[ '$' num2str(Ti(j(k))) '~\mathrm{K}$' ];
        semilogx(t,porosity(:,1),...
         'color', colors(mod(nseries,size(colors,1))+1,:),...
         'LineStyle',linestyle{mod(nseries,length(linestyle))+1},...
         'Linewidth',FS.lw);
       end
        hold on
        p=p+1;
    end
    
    if savefig
        graphtitle=['yield\Porosity' kinName 'Dp' num2str(2*r(i)*1E6) 'microm'];
    else
        graphtitle=['D_p=' num2str(2*r(i)*1E6) ' ' char(181) 'm'];
    end
    try
        makeItPretty(hdl,graphtitle,series,[mint maxt],ceil(maxy*10)/10 ,savefig,false)
    catch err
    end
end

end

function makeItPretty(ax,graphtitle,series,limt,maxy,savefig,makelegend)
    
    global ylabels FS LEL yst 
    warning off
    fig=get(ax(1),'parent');
    if savefig
        position=[0.11 0.14 0.87 0.82
                  0.14 0.14 0.737 0.785
                  0.16 0.14 0.785 0.785];
    else
        position(1,:)=get(ax(1),'position');
    end
    xlim(ax(1),limt);        
    ylim(ax(1),[0 maxy]);

    

    xorder=floor(log10(limt(1))):ceil(log10(limt(2)));
    yorder=0:1;

    ylabel(ax(1), ylabels{1},'interpreter','latex',...
        'HorizontalAlignment','center', 'VerticalAlignment','middle',...
        'rotation',90, 'fontsize',FS.y,'fontweight','bold')
    
    
    xlabel('$\mathsf{Time~[~s~]}$','fontsize', FS.x,'interpreter','latex') 
    set(ax,'fontsize',FS.ax)
    grid on
    set(ax(1), 'GridLineStyle','-','MinorGridLineStyle', '-',...
    'YMinorTick','on',...    %'YScale','log',...    %'YTick',10.^yorder,...
    'XMinorTick','on',...
    'XScale','log',...
    'XTick',10.^(xorder),...
    'TickLabelInterpreter','none',...    'yticklabel',yticks,...'xticklabel',xticks,...
    'linewidth',1.5,...
    'MinorGridAlpha',0.05,...
    'GridAlpha',0.25,...
    'MinorGridColor',[0 0 0],...
    'GridColor',[0 0 0],...
    'box','on',...
    'position',position(1,:),...
    'fontsize',FS.ax);
    

    %set(ax(2),'Ycolor',[0 0.2 0.6])
    %set(ax(3),'Ycolor',[0.6 0 0])
    ylab=get(ax,'ylabel');
    set(ylab,'units','normalized' )
    set(ylab ,'position',[-0.095 0.5  0 ])
    xlab=get(ax,'xlabel');
    set(xlab,'units','normalized' )
    set(xlab ,'position',[0.5 -0.08  0 ])
    
    hold off


    
    
    if(savefig)

        if(makelegend)
            legend(ax(1),series,'fontsize',FS.legend,'interpreter','latex',...
            'location','best')
        else
            labelthem(series,ax(1),FS)   
        end
        done=false;
        while ~done
            try
                graphtitle=strrep(graphtitle,'Wangenaar','Wagenaar');
                figname=['Figures\' strrep(graphtitle,' ','')];
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
    K=25;
               
    for i=1:length(series)
        orientation=abs(180/pi*atan(K*((series{i}.y(3)-series{i}.y(1)))/...
            (log10(series{i}.x(3)/series{i}.x(1)))));
        text(series{i}.x(2),series{i}.y(2),[series{i}.name],'Rotation',...
            orientation,'horizontalalignment','center','fontsize',FS.lbl,...
            'verticalalignment','bottom','interpreter','latex','parent',hdl)
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
