function [ M2,sigma,tau_p,tau_Rho,tau_T,tau_Pi] = solve_shock(M1,theta)
% to calcul normal and oblique shock waves
%   Airiau, avril 2012
% 
if  M1 < 1 
    disp('M upstream < 1 : error')
    return
end 

coef=pi/180;
if theta == 0
   disp('NORMAL shock wave')
   sigma = 90;
else 
    disp('OBLIQUE shock wave')
    % calculus of the shock angle sigma
	fprintf('theta            = %10f ° \n',theta);
    sigma=shock_angle(M1,theta);
end
fprintf('sigma            = %10f ° \n',sigma);
Mn1=M1*sin(sigma*coef);
fprintf('Upstream Mach    = %10f, Normal upstream Mach    = %10f \n',M1,Mn1)

% downstream normal Mach
Mn2=downstream_normal_Mach(Mn1);
M2=Mn2/sin((sigma-theta)*coef);
fprintf('Downstream Mach  = %10f, Normal downstream Mach  = %10f \n',M2,Mn2)

tau_p=P2overP1(Mn1);
tau_Rho=Rho2overRho1(Mn1);  
tau_T=tau_p/tau_Rho;
fprintf('Pressure jump    = %10f\n',tau_p);
fprintf('density jump     = %10f \n',tau_Rho);
fprintf('Temperature jump = %10f \n',tau_T);
tau_Pi=Pi2overPi1(tau_p,tau_Rho);
fprintf('Pi jump          = %10f  \n',tau_Pi);

end

