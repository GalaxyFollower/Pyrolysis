function giveup= plotDataAtStep(t,r_ev,X,param,init,hdl,plothem)
labels='BCTG';

giveup=false;

if (isempty(init)|| strcmp(init,'init'))
    reInterpretInputs
    giveup=max(rho(solids))>1.1*param.rhob0 ||... 
    max(max(rho(:,bio)))<0.01*param.rhob0 || range(T)>param.Tinf || toc>7200;
    msg=['[' param.kinetics.name ']' ...       
         ' Tinf:' num2str(param.Tinf) ...
         ' r:' num2str(param.r) ...
         ' nvert:' num2str(param.vter_mult) ...
         '\n\t\t Elapsed time: ' num2str(toc) 's. simtime: ' num2str((t)) ...
         ' rhob: ' num2str(mean(X(bionodes))) ...
         ' T='  num2str(mean(X(Tnodes))) '\n\n'];
    fprintf(msg);
end


if ((isempty(init)|| strcmp(init,'init')) && plothem)
% ------------------------------------------------------
% Reinterpretation of inputs
%------------------------------------------------------------
X=X(:,end);
t=t(end);
    reInterpretInputs
    
    dXdt=massAndEnergyBalances(t,X,r_ev,param);
    drhodt=zeros(nv+1,param.nc);
    drhodt(:,solids)=reshape(dXdt([solidnodes]),[nv+1 2]);

    
%------------------------------------------------------------
% Principal Polynomials
%------------------------------------------------------------
    calculateAllPolynomials
    for i=1:nc
        RXNaux(:,:,i)=volumeIntegral(a_Rrxn(:,:,i),intervals,r);    
    end


    RXN_RATES=[RXNaux(1:nv-1,:);sum(RXNaux(nv:nv+1,:),1)]; 
    
    %pressure=rho(:,gases)*(1./W(gases))*Rgas.*T./porosity-param.Pinf;
    if(giveup)
       giveup 
    end
    
    cla(hdl)
      if plothem  
            x=(0:r(end)/101:r(end))';
            for i=1:param.nc

               f=evaluatePolynomials(a_rho(:,:,i),intervals,x);
               gimmiedagraph(hdl(i),cord,rho(:,i),'*',{'radius (m)', ...
                   ['\rho_' labels(i)  '(kg\cdotm^{-3})']},'tex')
               hold(hdl(i),'on')
               gimmiedagraph(hdl(i),x,f,'-',{'radius (m)', ...
                   ['\rho_' labels(i)  '(kg\cdotm^{-3})']},'tex')
               hold(hdl(i),'off')
            end

            for i=1:param.nc

              % a_drhodt=getPolynomials(CM,drhodt(:,i),r,unit,order,param.deriv2);
               f=evaluatePolynomials(a_drhodr(:,:,i),intervals,x);


               gimmiedagraph(hdl(i+4),r_ev,drhodr_r(:,i),'*',{'radius (m)', ...
                   ['$\frac{\partial \rho_{_' labels(i)  '}}{\partial r}$']},'latex')
               hold(hdl(i+4),'on')
               gimmiedagraph(hdl(i+4),x,f,'-',{'radius (m)', ...
                   ['$\frac{\partial \rho_{_' labels(i)  '}}{\partial r}$']},'latex')
               hold(hdl(i+4),'off')
            end        



            f=evaluatePolynomials(a_T,intervals,x);
            gimmiedagraph(hdl(9),cord,T,'*',{'radius (m)','T (K)'},'tex')
    %         f=evaluatePolynomials(a_K,intervals,x);
    %         gimmiedagraph(hdl(9),r_ev,K_r,'*',{'radius (m)','K'},'tex')
            hold(hdl(9),'on')
            %set(hdl(9),'ylim',[1E-8 1E-4])
            %set(hdl(9),'YScale','log')
            gimmiedagraph(hdl(9),x,f,'-',{'radius (m)','T (K)'},'tex')
            hold(hdl(9),'off')


            f=evaluatePolynomials(a_P,intervals,x);
            gimmiedagraph(hdl(10),cord,P,'*',{'radius (m)','P-Patm (Pa)'},'tex')          
            hold(hdl(10),'on')
            gimmiedagraph(hdl(10),x,f,'-',{'radius (m)','P (Pa)'},'tex')
            hold(hdl(10),'off')

            f=evaluatePolynomials(a_dTdr,intervals,x);
            gimmiedagraph(hdl(11),r_ev,dTdr_r,'*',{'radius (m)','$\frac{\partial T}{\partial r}$'},'latex')
    %         f=evaluatePolynomials(a_dKdr,intervals,x);
    %         gimmiedagraph(hdl(11),r_ev,dKdr_r,'*',{'radius (m)','$\frac{\partial K}{\partial r}$'},'latex')
            hold(hdl(11),'on')
            gimmiedagraph(hdl(11),x,f,'-',{'radius (m)','$\frac{\partial T}{\partial r}$'},'latex')
            hold(hdl(11),'off')
            
            
            f=evaluatePolynomials(a_dPdr,intervals,x);
            gimmiedagraph(hdl(12),r_ev,dPdr_r,'*',{'radius (m)','$\frac{\partial P}{\partial r}$'},'latex')
            hold(hdl(12),'on')
            gimmiedagraph(hdl(12),x,f,'-',{'radius (m)','$\frac{\partial P}{\partial r}$'},'latex')
            hold(hdl(12),'off')
            pause(0.0000001)
            
            CM0=(0.75*(r_ev(2:nv+1).^4-r_ev(1:nv).^4)./( r_ev(2:nv+1).^3-r_ev(1:nv).^3));
            plot(hdl(13),[r_ev r_ev ]',[-ones(nv+1,1)  ones(nv+1,1)]','-k');
            hold(hdl(13),'on')
            plot(hdl(13),[ CM0 CM0]',[-ones(nv,1)  ones(nv,1)]',':k');
            plot(hdl(13),[ CM CM]',[-ones(nv,1)  ones(nv,1)]',':b');
            plot(hdl(13),[ CMT CMT]',[-ones(nv,1)  ones(nv,1)]',':b');
            plot(hdl(13),CM0,zeros(nv,1),'ok');
            plot(hdl(13),CM0,zeros(nv,1),'.k');
            
            plot(hdl(13),CM,0.5*ones(nv,1),'ob');
            plot(hdl(13),CM,0.5*ones(nv,1),'.b');
            
            plot(hdl(13),CMT,-0.5*ones(nv,1),'or');
            plot(hdl(13),CMT,-0.5*ones(nv,1),'.r');
            hold(hdl(13),'off')

      end
  
        
        
        
end



end

function gimmiedagraph(hdl,X,Y,specs,axisnames,interpreter)

    plot(hdl,X,Y,specs) 
    xlabel(hdl,axisnames{1})
    ylabel(hdl,axisnames{2},'interpreter',interpreter)
    
    
end