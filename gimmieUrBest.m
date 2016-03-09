%script to extend the integration to a further time

tf=input('How far do you want to go? gimmie the time in seconds pal\n>');

param.tf=tf;

t0=t;
X0=X;

tic
extratime=toc;
param.duration=param.duration+extratime;
[t X]=ode15s(@(t,X)massAndEnergyBalances(t,X,ri*1e-6,param),[ t0(end) tf],X0(end,:)',odeoptions);
fprintf('Elapsed time is %f seconds.\n%f total\n\n',extratime, param.duration)
X=[X0;X];
t=[t0;t];
plotData(t,ri*1E-6,X,param)
