function [F]=Rayleigh(Mach,opt)
% Rayleigh function
% opt= 1 : the function, opt= 2 : its derivative

global gam

switch opt
case 1
    F=2.0*(gam+1.0)*Mach^2/(1.0+gam*Mach^2)^2*(1.0+(gam-1.0)/2.0*Mach^2);
otherwise
    F=-4.0*Mach*(gam+1.0)*(Mach^2-1.0)/(1.0+gam*Mach^2)^3;
end 

end 
