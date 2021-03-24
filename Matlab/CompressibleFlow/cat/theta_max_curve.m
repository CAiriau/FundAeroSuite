function [theta,Sigma] = theta_max_curve(Mach)
%  curve of the maximal wall deviation angle for oblique shock 
% C Airiau, avril 2012

global gam
coef=180/pi;

a=-4*gam.*Mach.^4;
b=8*Mach.^2+(6*gam-2).*Mach.^4;
c=4+(2-2*gam).*Mach.^4+(2*gam-6).*Mach.^2;
Delt=b.*b-4.*a.*c;
x=(-b+sqrt(Delt))./(2.*a);
Sigma=acos(sqrt(x));
inv_tan_theta=((gam+1)/2.*Mach.^2./ (Mach.^2.*sin(Sigma).^2-1) -1 ) .*tan(Sigma);
theta=atan(1./inv_tan_theta)*coef;

if length(Mach) == 1 
    Mn=Mach.*sin(Sigma);
    f=1./Rho2overRho1(Mn);
    tan_theta=f.*tan(Sigma);
    theta2=(Sigma-atan(tan_theta))*coef;
    %fprintf('theta = %f, theta2 =  %f \n',theta,theta2);
end
Sigma=Sigma*coef;
return
end
