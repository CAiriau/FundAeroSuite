function [tau_Pi] = Pi2surPi1(tau_p,tau_Rho)
% Rapport des pressions isentropiques à travers un choc
% Airiau, Avril 2012
global gam
x=-1/(gam-1);
y=-gam*x;      
tau_Pi= tau_p^x * tau_Rho^y;
end

 