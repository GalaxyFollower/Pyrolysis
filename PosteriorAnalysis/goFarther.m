function goFarther(name)
load(name)
delete(['Failed/' name]);

succes=true;
t0=param.tf;
tf=max(param.tf*1.5,0.5);
X0=[rho(end,:,bio)'
    rho(end,:,cha)'
    rho(end,1:nv,gou)'
    rho(end,1:nv,gas)'
    rho(end,end,gas)/sum(rho(end,end,gases))
    P(end,1:nv)'
    T(end,1:nv)'
    CM0
    ];
figure(3)
for i=1:4
    hdl(i)=subplot(2,2,i);
end
figure(4)
for i=1:4
    hdl(i+4)=subplot(2,2,i);
end
figure(5)
for i=1:4
    hdl(i+8)=subplot(2,2,i);
end

odeoptions=odeset('Stats','on','OutputFCn', ...
    @(t,X,init)plotDataAtStep(t,ri,X,param,init,hdl,false),'NonNegative',1);
tic
while mean(X0(1:nv+1))>0.01*param.rhob0
    try
        [t X]=ode15s(@(t,X)massAndEnergyBalances(t,X,ri,param),...
            [ t0 tf],X0,odeoptions);
        
        time=[time;t(2:end)];
        XX=[XX;X(2:end,:)];
        if(t(end)<tf)
            succes=false;
        end
        t0=tf;
        tf=1.5*tf;
        X0=X(end,:)';
        
    catch err
        succes=false;
    end
    if ~succes
        X0(:)=0;
    end
    
end

i=find(mean(XX(:,1:nv+1),2)<.01*param.rhob0);
if isempty(i)
    i=length(time);
end

time=time(1:i(1));
XX=XX(1:i(1),:);

param.tf=time(end);
param.duration=toc;
save(['BrokenOnes\' name]);
