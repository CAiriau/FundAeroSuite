function [theta,Sigma] = unitary_downstream_Mach_curve(Mach)
% curve of the unitary downstream Mach number in oblique shock
% C Airiau, avril 2012

global gam
coef=180/pi;

c=2*(gam+1);
a=-gam*c;
b=Mach.^2*(gam+1)^2+gam^2-2*gam-3;
Delt=b.*b-4.*a.*c;
x=-(b+sqrt(Delt))./(2.*a);
Sigma=asin(sqrt(x)./Mach);
inv_tan_theta=((gam+1)/2.*Mach.^2./ (Mach.^2.*sin(Sigma).^2-1) -1 ) .*tan(Sigma);
theta=atan(1./inv_tan_theta)*coef;
Sigma=Sigma*coef;
return
end