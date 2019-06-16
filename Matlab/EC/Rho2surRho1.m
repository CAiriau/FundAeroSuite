function [ tau_Rho] = Rho2surRho1(Mn)
% Saut de masse volumique à travers un choc
% C. Airiau, avril 2012

global gam
tau_Rho= 1./( 2./((gam+1).*Mn.^2)+ (gam-1)/(gam+1) );
end

