clear 
clc
addpath('polynomialOperations')
addpath('AuxiliarScripts')
n=3;
order=3;
X=(0:3/(n+1):3)';


W=(0.75*(X(2:n+1).^4-X(1:n).^4)./( X(2:n+1).^3-X(1:n).^3));
cord=[W;X(end)];

Y=8*exp(-cord/8).*cos(cord*1.8);
Y=normpdf(cord,0,1)
%Y(end)=0;
%Y(n+1)=2;
[a intervals]=getPolynomials(W,Y,X,1,order);

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



subplot(3,1,1)
minf=min(f);
maxf=max(f);

minf=minf-0.1*(maxf-minf);
maxf=maxf+0.1*(maxf-minf);




for i=1:length(intervals)
    plot([intervals(i) intervals(i)],[minf maxf],':','color',[0.5 0.5 0.5],'linewidth',0.5);
    hold on
end
% xlbl=cellfun(@(x)strrep(x,' ',''),mat2cell([repmat('r_{',n+1,1) num2str((0:n)') repmat('}',n+1,1)],ones(n+1,1)),'uniformoutput',false);
% set(ax,'xtick',(0:n+1)', 'xticklabel',xlbl)

plot(x,f,'k','LineWidth',6)
%plot(x,y,'--','linewidth',2)
plot(cord,Y,'or','markersize',10,'linewidth',3)
if(minf~=maxf)
    ylim([minf maxf])
end
    
hold off

%% Differantiate and plot

dadx=differentiatePolynomials(a);

%plot

xx=ones(size(x));
for i=1:order-1
    xx=[x.^i xx];
end
y=zeros(size(x,1),size(dadx,1));
for i=1:size(dadx,1)
    y(:,i)=xx*dadx(i,:)';
end
f=evaluatePolynomials(dadx,intervals,x);

figure(1)
subplot(3,1,2)
minf=min(f);
maxf=max(f);
minf=minf-0.1*(maxf-minf);
maxf=maxf+0.1*(maxf-minf);


plot(x,f,'k','LineWidth',6)
hold on
%plot(x,y,'--','linewidth',2)

for i=1:length(intervals)
    plot([intervals(i) intervals(i)],[minf maxf],':','color',[0.5 0.5 0.5],'linewidth',0.5)
end
if(minf~=maxf)
    ylim([minf maxf])
end
hold off

%% Second Derivative and plot

d2adx2=differentiatePolynomials(dadx);

%plot

xx=ones(size(x));
for i=1:order-2
    xx=[x.^i xx];
end
y=zeros(size(x,1),size(d2adx2,1));
for i=1:size(d2adx2,1)
    y(:,i)=xx*d2adx2(i,:)';
end
f=evaluatePolynomials(d2adx2,intervals,x);

figure(1)
subplot(3,1,3)
minf=min(f);
maxf=max(f);
minf=minf-0.1*(maxf-minf);
maxf=maxf+0.1*(maxf-minf);


plot(x,f,'k','LineWidth',6)
hold on
%plot(x,y,'--','linewidth',2)

for i=1:length(intervals)
    plot([intervals(i) intervals(i)],[minf maxf],':','color',[0.5 0.5 0.5],'linewidth',0.5)
end
if(minf~=maxf)
    ylim([minf maxf])
end
hold off