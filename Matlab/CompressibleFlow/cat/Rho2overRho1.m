function [ tau_Rho] = Rho2overRho1(Mn)
% density ratio across a shock wave
% C. Airiau, avril 2012

global gam
tau_Rho= 1./( 2./((gam+1).*Mn.^2)+ (gam-1)/(gam+1) );
end