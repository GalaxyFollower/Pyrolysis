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

kinetics(1).name=['Wagenaar'];
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

c='kkkk';
dp=unique([1E-7;exp(log(1E-7):0.01:log(.1))';0.1]);

T=[400:20:580  600:839 840:20:1740 1760:1763 1764:0.1:1765 1766:1770  1775:5:2295    2300:100:2900 3000:200:5000]';

Cp=@(T)(polyval(param.Cp(1,:),param.T0));
Rgas=param.Rgas;
[TT, DP]=meshgrid(T,dp);
[m,n]=size(TT);
Cp=@(T)(polyval(param.Cp(1,:),T)+polyval(param.Cp(1,:),param.T0))/2;
t_int=param.rhob0*DP.^2.*Cp(TT)/(6*param.lambda(1));
t_ext= real(param.rhob0*Cp(TT).*(TT-param.T0).*DP./(12*(hconv(TT,param,DP).*(TT-param.T0)+param.stefBoltz*param.emis*(TT.^4-param.T0^4))));
    
for kin=1:4
    rp=@(T)sum(repmat(kinetics(kin).A,[length(T) 1]).*exp(-repmat(kinetics(kin).E,[length(T) 1])./(param.Rgas*repmat(T,[ 1 4]))),2)*param.rhob0;
    Da=@(T,d)rp(T)*abs(T-param.T0)*Cp(T)*d/(6*( hconv(T,param,d/2)*abs(T-param.T0)+ param.emis*param.stefBoltz*abs(T.^4 -param.T0.^4)));
    Bi=@(T,d)d.*( hconv(T,param,d/2).*abs(T-param.T0)+param.emis*param.stefBoltz.*abs(T.^4 -param.T0.^4))./(2*param.lambda(1)*abs(T -param.T0));
    Th=@(T,d) rp(T).*Cp(T).*d.^2/(12*param.lambda(1));
    
    %----------------------------------------------------------------------
    %                   TOTAL TIME CONTOURS
    %----------------------------------------------------------------------

        kinetics(kin).t_pyro=1./sum(repmat(kinetics(kin).A,size(T)).*exp(-repmat(kinetics(kin).E,size(T) )./(Rgas*repmat(T,[1 4]))),2);
        kinetics(kin).t_tot=t_int+t_ext+repmat(kinetics(kin).t_pyro',[m 1]);
        Contours=10.^(round(-9:7));

        %labels on the right
            [C,h]=contourf(TT(:,T>2500 & T<3000),DP(:,T>2500 & T<3000),...
                kinetics(kin).t_tot(:,T>2500 & T<3000),'showtext','off',...
                'fill','off','levellist',Contours,'labelspacing',100);
            
            set(get(h,'parent'),'gridlinestyle','-');
            hold on
            tl=clabel(C,Contours,'fontsize',10);
            arrayfun(@(i)set(tl(i),'string',['10^{ '...
                num2str(round(log10(str2double(tl(i).String)))) '} s'],...
             'horizontalalignment','left','verticalalignment',...
                'bottom', 'position',...
                [ tl(i).Position(1:2) 1]),2:2:length(tl))
                delete(tl(1:2:length(tl)))
        %labels on the bottom left
            [C,h]=contourf(TT(dp>1E-7 & dp<7E-6,T>410 & T<1600),...
                DP(dp>1E-7 & dp<7E-6,T>410 &  T<1600),...
                kinetics(kin).t_tot(dp>1E-7 & dp<7E-6,T>410 & T<1600),...
                'showtext','off','fill','off','levellist',Contours,...
                'labelspacing',100);
            tl=clabel(C,Contours,'fontsize',12,'rotation',90);
            arrayfun(@(i)set(tl(i),'string',['10^{ ' ...
                num2str(round(log10(str2double(tl(i).String)))) '} s'],...
             'horizontalalignment','left','verticalalignment',...
                'bottom', 'position',...
                [ tl(i).Position(1:2) 1]),2:2:length(tl))
                delete(tl(2))
                set(tl(6),'horizontalalignment','right')
                delete(tl(1:2:length(tl)))
        contourf(TT,DP,kinetics(kin).t_tot,'showtext','off','fill','off','levellist',Contours,'labelspacing',100);
        ax=get(h,'parent');
        set(ax,'xscale','log','yscale','log')
        
%   ---        
        rp=@(T)sum(repmat(kinetics(kin).A,[length(T) 1]).*exp(-repmat(kinetics(kin).E,[length(T) 1])./(param.Rgas*repmat(T,[ 1 4]))),2)*param.rhob0;
        Cp=@(T)(polyval(param.Cp(1,:),T)+polyval(param.Cp(1,:),param.T0))/2;
    
    %----------------------------------------------------------------------
    %                   THIELE NUMBER
    %----------------------------------------------------------------------
        d_Th=sqrt(12*param.lambda(1)./rp(T)./Cp(T));
        d_Th10=sqrt(120*param.lambda(1)./rp(T)./Cp(T));
        d_Th01=sqrt(1.2*param.lambda(1)./rp(T)./Cp(T));

        iTh=Bi(T,d_Th)<1;
        d_Th(iTh)=nan;
        semilogy(T,d_Th,['-' c(kin)],'linewidth',3)
        
        iTh=Bi(T,d_Th10)<=1;
        d_Th10(iTh)=nan;
        semilogy(T,d_Th10,'-','color',[0.5 0.5 0.5],'linewidth',2)

        iTh=Bi(T,d_Th01)<=1;
        d_Th01(iTh)=nan;
        semilogy(T,d_Th01,'-','color',[0.5 0.5 0.5],'linewidth',2)


    
    %----------------------------------------------------------------------
    %                   DAMKOLER NUMBER
    %----------------------------------------------------------------------
        Tda=T(T>600 & T<1800);    
        d_Da=zeros(size(Tda));
        d_Da10=zeros(size(Tda));
        d_Da01=zeros(size(Tda));

        dDa=@(T,d)(6*hconv(T,param,d/2)*abs(T-param.T0)+param.emis*param.stefBoltz.*abs(T.^4 -param.T0.^4) )/(rp(T)*abs(T-param.T0)*Cp(T));
        d_Da(1)= dDa(Tda(1),d_Th(1));
        d_Da10(1)=10*dDa(Tda(1),d_Th10(1));
        d_Da01(1)=0.1*dDa(Tda(1),d_Th01(1));

        TOL=1E-4;

        for i=2:length(Tda)

            d_Da(i)=dDa(Tda(i),d_Da(i-1));
            d_Da10(i)=10*dDa(Tda(i),d_Da10(i-1));
            d_Da01(i)=0.1*dDa(Tda(i),d_Da01(i-1));

            it=0;
            while ( abs(dDa(Tda(i),d_Da(i)) - d_Da(i) )/d_Da(i)>TOL || ...
                  abs(10*dDa(Tda(i),d_Da10(i)) - d_Da10(i) )/d_Da10(i)>TOL || ...
                  abs(0.1*dDa(Tda(i),d_Da01(i)) - d_Da01(i) )/d_Da01(i)>TOL ) && it<100

                d_Da(i)=dDa(Tda(i),d_Da(i));  
                d_Da10(i)=10*dDa(Tda(i),d_Da10(i));
                d_Da01(i)=0.1*dDa(Tda(i),d_Da01(i));
                it=it+1;
            end

        end

        iDa=Bi(Tda,d_Da)>1;
        d_Da(iDa)=nan;
        hold on
        semilogy(Tda,d_Da,['-' c(kin)] ,'linewidth',3)
        iDa10=Bi(Tda,d_Da10)>1;
        d_Da10(iDa10)=nan;
        semilogy(Tda,d_Da10,'-','color',[0.5 0.5 0.5],'linewidth',2)
        iDa01=Bi(Tda,d_Da01)>1;
        d_Da01(iDa01)=nan;
        semilogy(Tda,d_Da01,'-','color',[0.5 0.5 0.5],'linewidth',2)


    %----------------------------------------------------------------------
    %                   BIOT NUMBER
    %----------------------------------------------------------------------

        Tbi=T(T<2500);
        d_Bi=zeros(size(Tbi));
        d_Bi10=zeros(size(Tbi));
        
        
        dBi=@(T,d)2*param.lambda(1)*abs(T-param.T0)./(hconv(T,param,d/2)*abs(T-param.T0)+param.emis*param.stefBoltz.*abs(T.^4 -param.T0.^4) );
        d_Bi(1)=dBi(Tbi(1),0.0001);
        d_Bi10(1)=dBi(Tbi(1),0.0001);

        while abs(dBi(Tbi(1),d_Bi(1)) - d_Bi(1) )/d_Bi(1)>TOL || ...
                abs(10*dBi(Tbi(1),d_Bi10(1)) - d_Bi10(1) )/d_Bi10(1)>TOL %|| ...        abs(0.1*dBi(T(1),d_Bi01(1)) - d_Bi01(1) )/d_Bi01(1)>TOL
            d_Bi(1)=dBi(Tbi(1),d_Bi(1));    
            d_Bi10(1)=10*dBi(Tbi(1),d_Bi10(1));

        end

        for i=2:length(Tbi)

            d_Bi(i)=dBi(Tbi(i),d_Bi(i-1));
            d_Bi10(i)=10*dBi(Tbi(i),d_Bi10(i-1));


            it=0;
            while ( abs(dBi(Tbi(i),d_Bi(i)) - d_Bi(i) )/d_Bi(i)>TOL || ...
                  abs(10*dBi(Tbi(i),d_Bi10(i)) - d_Bi10(i) )/d_Bi10(i)>TOL  ...||          abs(0.5*dBi(T(i),d_Bi05(i)) - d_Bi05(i) )/d_Bi05(i)>TOL
                  ) && it<100

                d_Bi(i)=dBi(Tbi(i),d_Bi(i));  
                d_Bi10(i)=10*dBi(Tbi(i),d_Bi10(i));

                it=it+1;

            end

        end
        

        
        d_Bi1=d_Bi;
        iBi1=Th(Tbi,d_Bi1)>1;
        d_Bi1(iBi1)=nan;
        iBi=find(~(Th(Tbi,d_Bi)<1));
        d_Bi1(max(iBi))=d_Da(max(find(~iDa)));
        semilogy(Tbi,d_Bi1,'-','color',[0.5 0.5 0.5],'linewidth',2)
        semilogy(Tbi,d_Bi10,'-','color',[0.5 0.5 0.5],'linewidth',2)

        d_Bi1=d_Bi;
        iBi=find(Th(Tbi,d_Bi)<1);
        d_Bi1(iBi)=nan;
        iBi=find(~(Th(Tbi,d_Bi)<1));
        d_Bi1(max(iBi))=d_Da(max(find(~iDa)));
        semilogy(Tbi,d_Bi1,'-k','linewidth',3)
                
        
    %----------------------------------------------------------------------
    %                   FIGURE DETAILS
    %----------------------------------------------------------------------

        ylim([1e-7 .1]);
        xlim([400 5000]);
        ax=get(1,'children');
        ylabel('Particle diametre, d_p( m )')
        xlabel('External Temperature, T_\infty( K )')
        scl=[400:100:1000 1100:100:10000];
        lbl=[400:200:1000 1400 2000:1000:10000];
        j=ismember(scl,lbl);        
        ylbl=arrayfun(@(x)num2str(x), scl,'uniformoutput',false);
        for i=find(~j), ylbl{i}=' '; end
        set(ax,'fontsize',16,'ytick',10.^(-7:0),'xtick',scl,'xticklabel',ylbl ,'xgrid','on','ygrid','on','linewidth',0.5,'gridlinestyle','none');
        
%         annotation('textbox',[0.2 0.3 0.2 0], 'String','Process limited by pyrolysis reaction','fontsize',16,'linestyle','none','horizontalalignment','center')
%         annotation('textbox',[0.45 0.55 0.2 0], 'String','Process limited by external heat transfer','fontsize',16,'linestyle','none','horizontalalignment','center')
%         annotation('textbox',[0.65 0.75 0.2 0], 'String','Process limited by internal heat transfer','fontsize',16,'linestyle','none','horizontalalignment','center')
        hold off
        set(1,'position',[2251 34 762 650]);
        set(get(1,'children'),'position',[0.13 0.1169 0.8327 0.8548])

       export_fig(1,['regimetimes' kinetics(kin).name],'-pdf','-transparent')
end
