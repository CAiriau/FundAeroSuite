function [omega] = Prandtl_Meyer(Mach)
% Prandtl-Meyer function (angle in degrees) as a function of the Mach number
% C.Airiau, avril 2012
global gam
x=Mach.^2-1;
coef=180/pi;
omega=coef*(sqrt((gam+1)/(gam-1))*atan(sqrt((gam-1)/(gam+1).*x))-atan(sqrt(x)));
end

