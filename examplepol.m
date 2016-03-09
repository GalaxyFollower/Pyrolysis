clear 
clc
addpath('polynomialOperations')
addpath('AuxiliarScripts')
n=10;
order=3;
X=(0:n)';

W=(0.75*(X(2:n+1).^4-X(1:n).^4)./( X(2:n+1).^3-X(1:n).^3));
cord=[W;X(end)];
%X=[mean([W(2:n+1) W(1:n)],2)];

%Y=[X;W(n)];


%Y=300*(ones(length(W),1 ));



%Y=1000*exp([X;n]/5);
%Y=normpdf(cord,0,1)
%Y=1000*exp(-cord/5);
Y=8*exp(-cord/10).*cos(cord*5);
Y=8*exp(-cord/5).*cos(cord*1.5);
%Y(end)=0;
%Y(n+1)=2;
[a intervals]=getPolynomials(W,Y,X,1,order, true);




%% Make a graph

x=(X(1):0.001:X(end))';

xx=ones(size(x));
for i=1:order
    xx=[x.^i xx];
end
y=zeros(size(x,1),size(a,1));
for i=1:size(a,1)
    y(:,i)=xx*a(i,:)';
end
f=evaluatePolynomials(a,intervals,x);

fig=figure(11);
clf
%subplot(3,1,1)
minf=min(f);
maxf=max(f);

minf=minf-0.1*(maxf-minf);
maxf=maxf+0.1*(maxf-minf);




for i=1:n+1
    plot([intervals(i) intervals(i)],[2*minf 2*maxf],':','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
end

ax=get(fig,'children');
xlbl=cellfun(@(x)strrep(x,' ',''),mat2cell([repmat('$r_{',n+1,1) num2str((0:n)') repmat('}$',n+1,1)],ones(n+1,1)),'uniformoutput',false);
xlbl{n}='$r_n$';
xlbl{n-1}='$r_{n-1}$';


set(ax,'ytick',[],'xtick',(0:n+1)', 'xticklabel',xlbl,'ticklabelinterpreter','tex','fontsize',16)
ylabel('$y(r)$','interpreter','latex','rotation',0,'fontsize',24);
xlabel('$r$','interpreter','latex','fontsize',30);
cmlbl=cellfun(@(x)strrep(x,' ',''),mat2cell([repmat('$( \overline{r}_{',n,1) num2str((1:n)') repmat('}, \overline{y}_{',n,1) num2str((1:n)') repmat('})$',n,1)],ones(n,1)),'uniformoutput',false);
cmlbl{n-1}='$(\overline{r}_{n-1},\overline{y}_{n-1})$'; 
cmlbl{n}='$(\overline{r}_{n},\overline{y}_{n})$'; 
cmlbl{n+1}='$(r_n,y_n)$'; 
    ylim([minf maxf])

for i=[1:4 n-1:n+1]
    text(cord(i),Y(i),cmlbl{i},'interpreter','latex','fontsize',16);
end
for i=[1:4 n-2:n+1]
    y=polyval(a(i,:),x);
    plot(x,y,'--k','linewidth',1.5)
end

 
plot(x,f,'k','LineWidth',3)
plot(cord,Y,'ok','markersize',10,'linewidth',3,'markerfacecolor',[1 1 1])


if(minf~=maxf)
    ylim([minf maxf])
end
breakxaxis(ax,[3.9 8.1])
    
hold off

