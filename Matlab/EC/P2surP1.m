function [tau_p] =P2surP1(Mn)
% Saut de pression à travers un choc
%  Airiau, avril 2012  
global gam
tau_p= 2*gam/(gam+1)* Mn^2- (gam-1)/(gam+1);
end  

