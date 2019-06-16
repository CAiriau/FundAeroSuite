function [sigma] = valeur_sigma(M1, theta_d)

% ------------------------------------------
% Usage : sigma = valeur_sigma(M1, theta)
%
% avec M1 = nombre de Mach incident
% et theta = angle de deviation (en degres)
% (sigma est calcule en degres)
% ------------------------------------------

if nargin ~= 2
  help valeur_sigma
  return
end

global gam;


theta = theta_d/180*pi;

sigma = fzero(@(sigma)f(sigma, M1, theta), pi/4, optimset('TolX', 1e-8, 'TolFun', 1e-8));

sigma = sigma/pi*180;
%fprintf('sigma = %5f degrés \n',sigma)
 
function y = f(sigma, M1, theta)

global gam;

y = tan(sigma-theta)./tan(sigma) - 2/(gam+1)/(M1*sin(sigma))^2 - (gam-1)/(gam+1);

 
