

Re=exp(0:0.01:20)';
Pr=(50:50:400);

RE=repmat(Re,[1 length(Pr)]);
PR=repmat(Pr,[length(Re) 1]);

figure(10)

re=1<=Re & Re<=7E4;
pr=0.6<=Pr & Pr<=400;
Nu1=2+0.6*PR(re,pr).^(1/3).*RE(re,pr).^(1/2);
plot(Re(re),Nu1,'-')
hold on


re=10<=Re & Re<=1E6;
pr=0.6<=Pr & Pr<=380;
Nu2=2+1.3*PR(re,pr).^(0.15)+0.66*PR(re,pr).^0.31.*RE(re,pr).^(1/2);
plot(Re(re),Nu2,':')


Nu3=2+0.03*PR.^0.33.*RE.^0.54+0.35*PR.^0.36.*RE.^0.58;
plot(Re,Nu3,'--')


figure(11)
kT_e=(0.3:0.05:90)';
A=1.16145;
B=0.14874;
C=0.52487;
D=0.7732;
E=2.16178;
F=2.43787;
Omega=A*kT_e.^-B+C*exp(-D*kT_e)+E*exp(-F*kT_e);






hold off