function [F]=Fanno(Mach,opt)
% fonction de Fanno et sa dérivée
% opt= 1 : la fonction, opt= 2 : sa dérivée

global gam

omega=1.+(gam-1.)/2.*Mach^2;
switch opt
case 1
    F= (1.-Mach^2)/(gam*Mach^2)+(gam+1)/(2.*gam)*log(0.5*(gam+1.)*Mach^2/omega);
otherwise
    F=-2.*(1.-Mach^2)/(gam*Mach^3*omega);
end 

end 
