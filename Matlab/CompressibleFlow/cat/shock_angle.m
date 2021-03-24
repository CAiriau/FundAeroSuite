function [sigma] = shock_angle(M1, theta_d)
%  get the shock angle in degrees for oblique shock  wave
% Usage : sigma = shock_angle(M1, theta)
% 
% Args:
%     M1 : upstream Mach number
%     theta_d : wall deviation angle in degrees
% Returns:
%   sigma,  shock angle  :math:`\sigma` in degree

if nargin ~= 2
  help shock_angle
  return
end

global gam;


theta = theta_d/180*pi;
sigma = fzero(@(sigma)f(sigma, M1, theta), pi/4, optimset('TolX', 1e-8, 'TolFun', 1e-8));
sigma = sigma/pi*180;
%fprintf('sigma = %5f degrés \n',sigma)

function y = f(sigma, M1, theta)
% function to solve to get oblique shock angle
global gam;
y = tan(sigma-theta)./tan(sigma) - 2/(gam+1)/(M1*sin(sigma))^2 - (gam-1)/(gam+1);

