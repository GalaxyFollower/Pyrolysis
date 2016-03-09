function minsqdiff=objectiveFunction(X,param)
e_k=X(1);
sigma=X(2);
Tcond=[100	200	300	400	500	600	700	800	1500];
lambda_exp=[0.00934	0.0181	0.0263	0.0338	0.0407	0.0469	0.0524	0.0573	0.1];

Tvisc=[100	200	300	400	500	600	700	800	900 1000 1100 1200 1300];
mu_exp=[0.00000711	0.00001325	0.00001846	0.00002301	0.00002701	0.00003058	0.00003388	0.00003698	0.00003981	0.00004244	0.0000449	0.0000473	0.0000496];

T=union(Tcond,Tvisc);
[mu_the, lambda_the]=hconv(T,e_k,sigma, param);

[C,ia,ib] = intersect(T,Tvisc);

mu_the=mu_the(ia);
 
[C,ia,ib] = intersect(T,Tcond) ;
lambda_the=lambda_the(ia);
 
minsqdiff=(sum(((mu_exp-mu_the)./mu_exp).^2)...
    +sum(((lambda_exp-lambda_the)./lambda_exp).^2))*1000;
    


