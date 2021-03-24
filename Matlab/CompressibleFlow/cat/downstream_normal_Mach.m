function [Mn2]=downstream_normal_Mach(Mn)
% Normal downstream Mach number in  shock wave
% C. Airiau avril 2012
% Mn: upstream normal Mach number
global gam
Mn2= sqrt((1+ 0.5*(gam-1)* Mn.^2)./(gam*Mn.^2-0.5*(gam-1)));
end 

