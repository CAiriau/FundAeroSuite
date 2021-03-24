function [tau_Pi] = Pi2overPi1(tau_p,tau_Rho)
% Total pressure jump across a shock wave 
% Airiau, Avril 2012
global gam
x=-1/(gam-1);
y=-gam*x;      
tau_Pi= tau_p^x * tau_Rho^y;
end