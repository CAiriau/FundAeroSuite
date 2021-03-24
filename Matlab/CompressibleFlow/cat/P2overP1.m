function [tau_p] =P2overP1(Mn)
% Pressure jump across a normal shock
%  Airiau, avril 2012  
global gam
tau_p= 2*gam/(gam+1)* Mn^2- (gam-1)/(gam+1);
end  

