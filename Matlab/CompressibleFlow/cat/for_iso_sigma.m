function [theta,Sigma]=for_iso_sigma(Mach)
% fast calculus of theta w.r.t. sigma
% C Airiau, Avril 2012
global gam
coef=180/pi;
n=500;
Sigma_min=asin(1/Mach); Sigma_max=pi/2;
dSigma=(Sigma_max-Sigma_min)/(n-1);
Sigma=Sigma_min:dSigma:Sigma_max;
inv_tan_theta=((gam+1)/2*Mach^2./ (Mach^2*sin(Sigma).^2-1) -1 ) .*tan(Sigma);
theta=atan(1./inv_tan_theta)*coef;
Sigma=Sigma*coef;
end