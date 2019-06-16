function [omega] = valeur_omega(Mach)
% calcul de la fonction omega (angle en degrés) en fonction du Mach
% C.Airiau, avril 2012
global gam
x=Mach.^2-1;
coef=180/pi;
omega=coef*(sqrt((gam+1)/(gam-1))*atan(sqrt((gam-1)/(gam+1).*x))-atan(sqrt(x)));
end

