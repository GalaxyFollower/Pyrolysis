function giveup= plotPreDataAtStep(t,r_ev,X,param,init,hdl,plothem)
labels='BCTG';

giveup=false;

if (isempty(init)|| strcmp(init,'init'))
    reInterpretInputs
    giveup=max(rho(solids))>1.1*param.rhob0 ||... 
    max(max(rho(:,solids)))<0.01*param.rhob0 || range(T)>param.Tinf ||...
    T(end)>500||toc>7200 ;
    msg=['[' param.kinetics.name ']' ...       
         ' Tinf:' num2str(param.Tinf) ...
         ' r:' num2str(param.r) ...
         ' nvert:' num2str(param.vter_mult) ...
         '\n\t\t Elapsed time: ' num2str(toc) 's. simtime: ' num2str(t) ...
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
    
    dXdt=massAndEnergyBalances(t,X,r,param);
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
               f=evaluatePolynomials(a_Rrxn,intervals,x);


               gimmiedagraph(hdl(i+4),CM,RXN_RATES(:,i),'*',{'radius (m)', ...
                   ['$\dot R_{_' labels(i)  '}$']},'latex')
               hold(hdl(i+4),'on')
%                gimmiedagraph(hdl(i+4),x,f,'-',{'radius (m)', ...
%                    ['$\frac{\dot R_{_' labels(i)  '}}{\partial t}$']},'latex')
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
            gimmiedagraph(hdl(10),cord,P-param.Pinf,'*',{'radius (m)','P-Patm (Pa)'},'tex')          
            hold(hdl(10),'on')
            gimmiedagraph(hdl(10),x,f-param.Pinf,'-',{'radius (m)','P (Pa)'},'tex')
            hold(hdl(10),'off')

            f=evaluatePolynomials(a_lambda,intervals,x);
            gimmiedagraph(hdl(11),r_ev,lambda_r,'*',{'radius (m)','$\lambda_\mathrm{eff}$'},'latex')
    %         f=evaluatePolynomials(a_dKdr,intervals,x);
    %         gimmiedagraph(hdl(11),r_ev,dKdr_r,'*',{'radius (m)','$\frac{\partial K}{\partial r}$'},'latex')
            hold(hdl(11),'on')
            gimmiedagraph(hdl(11),x,f,'-',{'radius (m)','$\lambda_\mathrm{eff}$'},'latex')
            hold(hdl(11),'off')
            
            [a_pre, intervals] =getPolynomials(CM,pressure-P,r_ev,unit,order,param.deriv2);
            f=evaluatePolynomials(a_pre,intervals,x);
            gimmiedagraph(hdl(12),r_ev,pressure-P,'*',{'radius (m)','$\frac{\partial P}{\partial r}$'},'latex')
            hold(hdl(12),'on')
            gimmiedagraph(hdl(12),x,f,'-',{'radius (m)','$\frac{\partial P}{\partial r}$'},'latex')
            hold(hdl(12),'off')
            pause(0.0000001)

      end
  
        
        
        
end



end

function gimmiedagraph(hdl,X,Y,specs,axisnames,interpreter)

    plot(hdl,X,Y,specs) 
    xlabel(hdl,axisnames{1})
    ylabel(hdl,axisnames{2},'interpreter',interpreter)
    
    
end