function [Mn2]=Mach_aval(Mn)
% Mach aval pour un choc
% C. Airiau avril 2012
% Mn: mach normal amont
% M_aval : mach normal aval
global gam
Mn2= sqrt((1+ 0.5*(gam-1)* Mn.^2)./(gam*Mn.^2-0.5*(gam-1)));
end 

