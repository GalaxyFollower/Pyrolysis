
%Internal energy of reaction
function deltaUrxn=Urxn(T,param)
    bio=1; %biomass
    cha=2;
    tar=3;
    gas=4;
    
    Cv=param.Cv;
    a_U=integratePolynomials(Cv);
    Hrxn=param.kinetics.Hrxn;
    
    deltaUrxn(:,1)=polyval(a_U(gas,:)-a_U(bio,:),T) ...
              -polyval(a_U(gas,:)-a_U(bio,:),ones(size(T))*param.Tref)...
              +Hrxn(1);

    deltaUrxn(:,2)=polyval(a_U(tar,:)-a_U(bio,:),T) ...
              -polyval(a_U(tar,:)-a_U(bio,:),ones(size(T))*param.Tref)...
              +Hrxn(2);
           
    deltaUrxn(:,3)=polyval(a_U(cha,:)-a_U(bio,:),T) ...
              -polyval(a_U(cha,:)-a_U(bio,:),ones(size(T))*param.Tref)...
              +Hrxn(3);
          
    deltaUrxn(:,4)=polyval(a_U(gas,:)-a_U(tar,:),T) ...
              -polyval(a_U(gas,:)-a_U(tar,:),ones(size(T))*param.Tref)...
              +Hrxn(4);
end

