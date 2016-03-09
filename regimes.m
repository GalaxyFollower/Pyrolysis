%%%
%finding the regimes as functions of the diameter and temperature ...
%according to Di benedeto, modelling the effects...
%%%
clc
clear
load param
%Da=1
figure(1)
clf

kinetics(1).name=['Wangenaar'];
kinetics(1).A=[1.11E11 9.28E9 3.05E7 0];
kinetics(1).E=[177E3 149e3 125e3 0];

kinetics(2).name=['DiBlasi'];
kinetics(2).A=[4.4E9 1.1E10 3.3E6 0];
kinetics(2).E=[153E3 148e3 112e3 0];

kinetics(3).name=['Chan'];
kinetics(3).A=[1.3e8 2e8 1.08e7 0];
kinetics(3).E=[140e3 133e3 121e3 0];

kinetics(4).name=['Font'];
kinetics(4).A=[1.52E7 5.85E6 2.98E3 0];
kinetics(4).E=[139E3 119e3 73e3 0];

c='bgrk';

for kin=1:4
rp=@(T)sum(repmat(kinetics(kin).A,[length(T) 1]).*exp(-repmat(kinetics(kin).E,[length(T) 1])./(param.Rgas*repmat(T,[ 1 4]))),2)*param.rhob0;
Cp=@(T)(polyval(param.Cp(1,:),T)+polyval(param.Cp(1,:),param.T0))/2;
Cp=@(T)(polyval(param.Cp(1,:),param.T0));
dp=unique([1E-6;exp(log(1E-6):0.01:log(500E-6))';500E-6]);
T=(300:20:10000)';
Rgas=param.Rgas;


Da=@(T,d)rp(T)*abs(T-param.T0)*Cp(T)*d/( 6*(hconv(T,param,d/2)*abs(T-param.T0)+ param.emis*param.stefBoltz*abs(T.^4 -param.T0.^4)));
Bi=@(T,d)d.*( hconv(T,param,d/2).*abs(T-param.T0)+param.emis*param.stefBoltz.*abs(T.^4 -param.T0.^4))./(2*param.lambda(1)*abs(T -param.T0));
Th=@(T,d) rp(T).*Cp(T).*d.^2/(12*param.lambda(1));

T=[400:20:580  600:839 840:20:1680 1700:1763 1764:0.1:1765 1766:1770  1775:5:2295    2300:100:2900 3000:1000:5000]';
d_Da=zeros(size(T));
d_Da10=zeros(size(T));
d_Da01=zeros(size(T));
d_Bi=zeros(size(T));
d_Bi10=zeros(size(T));
d_Bi01=zeros(size(T));

d_Th=sqrt(12*param.lambda(1)./rp(T)./Cp(T));
d_Th10=sqrt(120*param.lambda(1)./rp(T)./Cp(T));
d_Th01=sqrt(1.2*param.lambda(1)./rp(T)./Cp(T));


dDa=@(T,d)(6*(hconv(T,param,d/2)*abs(T-param.T0)+param.emis*param.stefBoltz.*abs(T.^4 -param.T0.^4)) )/(rp(T)*abs(T-param.T0)*Cp(T));
d_Da(1)= dDa(T(1),d_Th(1));
d_Da10(1)=10*dDa(T(1),d_Th10(1));
d_Da01(1)=0.1*dDa(T(1),d_Th01(1));
TOL=1E-4;
while abs(dDa(T(1),d_Da(1)) - d_Da(1) )/d_Da(1)>TOL || ...
        abs(10*dDa(T(1),d_Da10(1)) - d_Da10(1) )/d_Da10(1)>TOL || ...
        abs(0.1*dDa(T(1),d_Da01(1)) - d_Da01(1) )/d_Da01(1)>TOL
    d_Da(1)=dDa(T(1),d_Da(1));    
    d_Da10(1)=10*dDa(T(1),d_Da10(1));
    d_Da01(1)=0.1*dDa(T(1),d_Da01(1));
end

for i=2:length(T)
	  
    d_Da(i)=dDa(T(i),d_Da(i-1));
    d_Da10(i)=10*dDa(T(i),d_Da10(i-1));
    d_Da01(i)=0.1*dDa(T(i),d_Da01(i-1));
    
    it=0;
    while ( abs(dDa(T(i),d_Da(i)) - d_Da(i) )/d_Da(i)>TOL || ...
          abs(10*dDa(T(i),d_Da10(i)) - d_Da10(i) )/d_Da10(i)>TOL || ...
          abs(0.1*dDa(T(i),d_Da01(i)) - d_Da01(i) )/d_Da01(i)>TOL ) && it<100
   
        d_Da(i)=dDa(T(i),d_Da(i));  
        d_Da10(i)=10*dDa(T(i),d_Da10(i));
        d_Da01(i)=0.1*dDa(T(i),d_Da01(i));
        it=it+1;

    end
   
end



%%
dBi=@(T,d)2*param.lambda(1)*abs(T-param.T0)./(hconv(T,param,d/2)*abs(T-param.T0)+param.emis*param.stefBoltz.*abs(T.^4 -param.T0.^4) );
d_Bi(1)=dBi(T(1),0.0001);
d_Bi10(1)=dBi(T(1),0.0001);
%d_Bi01(1)=dBi(T(1),0.0001);


while abs(dBi(T(1),d_Bi(1)) - d_Bi(1) )/d_Bi(1)>TOL || ...
        abs(10*dBi(T(1),d_Bi10(1)) - d_Bi10(1) )/d_Bi10(1)>TOL %|| ...        abs(0.1*dBi(T(1),d_Bi01(1)) - d_Bi01(1) )/d_Bi01(1)>TOL
    d_Bi(1)=dBi(T(1),d_Bi(1));    
    d_Bi10(1)=10*dBi(T(1),d_Bi10(1));
    %d_Bi01(1)=0.1*dBi(T(1),d_Bi01(1));
end

for i=2:length(T)
	    
    d_Bi(i)=dBi(T(i),d_Bi(i-1));
    d_Bi10(i)=10*dBi(T(i),d_Bi10(i-1));
   % d_Bi05(i)=0.5*dBi(T(i),d_Bi05(i-1));
    
    it=0;
    while ( abs(dBi(T(i),d_Bi(i)) - d_Bi(i) )/d_Bi(i)>TOL || ...
          abs(10*dBi(T(i),d_Bi10(i)) - d_Bi10(i) )/d_Bi10(i)>TOL  ...||          abs(0.5*dBi(T(i),d_Bi05(i)) - d_Bi05(i) )/d_Bi05(i)>TOL
          ) && it<100
   
        d_Bi(i)=dBi(T(i),d_Bi(i));  
        d_Bi10(i)=10*dBi(T(i),d_Bi10(i));
        %d_Bi01(i)=0.5*dBi(T(i),d_Bi05(i));
        it=it+1;

    end
   
end


iDa1=Bi(T,d_Da)>1;
d_Da1=d_Da;
d_Da1(iDa1)=nan;
hold on
semilogy(T,d_Da1,['-' c(kin)] ,'linewidth',2)
iDa=Bi(T,d_Da10)>1;
d_Da10(iDa)=nan;
semilogy(T,d_Da10,[':' c(kin)] ,'linewidth',2)
iDa=Bi(T,d_Da01)>1;
d_Da01(iDa)=nan;
semilogy(T,d_Da01,[':' c(kin)] ,'linewidth',2)




iBi=Th(T,d_Bi)<=1;
semilogy(T,d_Bi,[':' c(kin)],'linewidth',2)
d_Bi(iBi)=nan;

semilogy(T,d_Bi10,[':' c(kin)],'linewidth',2)
%semilogy(T,d_Bi05,[':' c(kin)],'linewidth',2)

        d_Bi1=d_Bi;
        iBi1=Th(T,d_Bi1)>1;
        d_Bi1(iBi1)=nan;
        iBi=find(~(Th(T,d_Bi)<1));
        d_Bi1(max(iBi1))=d_Da(max(find(~iDa1)));
        semilogy(T,d_Bi1,['-' c(kin)],'linewidth',2)
        semilogy(T,d_Bi10,[':' c(kin)],'linewidth',2)

        d_Bi1=d_Bi;
        iBi=find(Th(T,d_Bi)<1);
        d_Bi1(iBi)=nan;
        iBi=find(~(Th(T,d_Bi)<1));
        d_Bi1(max(find(iBi1)))=d_Da(max(find(~iDa1)));
        semilogy(T,d_Bi1,['-' c(kin)],'linewidth',2)



iTh=Bi(T,d_Th)<1;
d_Th(iTh)=nan;
semilogy(T,d_Th,['-' c(kin)],'linewidth',2)
iTh=Bi(T,d_Th10)<=1;
d_Th10(iTh)=nan;
semilogy(T,d_Th10,[':' c(kin)],'linewidth',2)

iTh=Bi(T,d_Th01)<=1;
d_Th01(iTh)=nan;
semilogy(T,d_Th01,[':' c(kin)],'linewidth',2)




ylim([1e-7 1]);
xlim([300 2000]);
ax=get(1,'children');
ylabel('Particle diametre, d_p( m )')
xlabel('External Temperature, T_\infty( K )')
set(ax,'fontsize',16,'ytick',10.^(-7:0),'xtick',400:200:2000,'xgrid','on','ygrid','on','linewidth',0.5,'minorgridlinestyle','-');
annotation('textbox',[0.2 0.3 0.2 0], 'String','Process limited by pyrolysis reaction','fontsize',16,'linestyle','none','horizontalalignment','center')
annotation('textbox',[0.45 0.55 0.2 0], 'String','Process limited by external heat transfer','fontsize',16,'linestyle','none','horizontalalignment','center')
annotation('textbox',[0.65 0.75 0.2 0], 'String','Process limited by internal heat transfer','fontsize',16,'linestyle','none','horizontalalignment','center')
%hold off
end
%%
lines=get(get(1,'children'),'children');
Font.solid=[];
Font.dashed=[];

Chan.solid=[];
Chan.dashed=[];

DiBlasi.solid=[];
DiBlasi.dashed=[];

Wagenaar.solid=[];
Wagenaar.dashed=[];

Solid=[];
Dashed=[];

for i=1:length(lines)
    switch( true)
        case isequal(get(lines(i),'color'),[0 0 0])
            if strcmp(get(lines(i),'linestyle'),'-')
                Font.solid=[Font.solid ;lines(i)];
                Solid=[Solid;lines(i)];
            else
                Font.dashed=[Font.dashed ;lines(i)];
                Dashed=[Dashed;lines(i)];
            end
        case isequal(get(lines(i),'color'),[1 0 0])
            if strcmp(get(lines(i),'linestyle'),'-')
                Chan.solid=[Chan.solid ;lines(i)];
                Solid=[Solid;lines(i)];
            else
                Chan.dashed=[Chan.dashed ;lines(i)];
                Dashed=[Dashed;lines(i)];
            end
            
        case isequal(get(lines(i),'color'),[0 0 1])
            
            if strcmp(get(lines(i),'linestyle'),'-')
                Wagenaar.solid=[Wagenaar.solid ;lines(i)];
                Solid=[Solid;lines(i)];
            else
                Wagenaar.dashed=[Wagenaar.dashed ;lines(i)];
                Dashed=[Dashed;lines(i)];
            end
        otherwise
            if strcmp(get(lines(i),'linestyle'),'-')
                DiBlasi.solid=[DiBlasi.solid ;lines(i)];
                Solid=[Solid;lines(i)];
            else
                DiBlasi.dashed=[DiBlasi.dashed ;lines(i)];
                Dashed=[Dashed;lines(i)];
            end
    
   end
end
get(DiBlasi.solid,'color')
set(DiBlasi.solid,'color',[0 0.5 0])
set(DiBlasi.dashed,'color',[0 0.5 0])
