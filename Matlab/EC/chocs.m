function [ M2,sigma,tau_p,tau_Rho,tau_T,tau_Pi] = chocs(M1,theta)
% Calcul d'une onde de choc droite ou oblique
%   Airiau, avril 2012
% 
if  M1 < 1 
    disp('M amont < 1 : faux')
    return
end 

coef=pi/180;
if theta == 0
   disp('ONDE DE CHOC DROITE')
   sigma = 90;
else 
    disp('ONDE DE CHOC OBLIQUE')
    % calcul de sigma
	fprintf('theta            = %10f ° \n',theta);
    sigma=valeur_sigma(M1,theta);
end
fprintf('sigma            = %10f ° \n',sigma);
Mn1=M1*sin(sigma*coef);
fprintf('Mach amont       = %10f, Mach normal amont = %10f \n',M1,Mn1)

% calcul du mach normal aval
Mn2=Mach_aval(Mn1);
M2=Mn2/sin((sigma-theta)*coef);
fprintf('Mach aval        = %10f, Mach normal aval  = %10f \n',M2,Mn2)

tau_p=P2surP1(Mn1);
tau_Rho=Rho2surRho1(Mn1);  
tau_T=tau_p/tau_Rho;
fprintf('PAval/PAmont     = %10f\n',tau_p);
fprintf('RhoAval/RhoAmont = %10f \n',tau_Rho);
fprintf('TAval/TAmont     = %10f \n',tau_T);
tau_Pi=Pi2surPi1(tau_p,tau_Rho);
fprintf('PiAval/PiAmont   = %10f  \n',tau_Pi);

end

