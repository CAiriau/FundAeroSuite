function [a,b]=rotation(i,j,ul,wl,th)
% rotation of the base

%%*******************************************************
% Original coding obtained from J. Katz, in fortran 77
% and found in :
% Low-Speed Aerodynamics,Joseph Katz et Allen Plotkin
% Second Edition,Cambridge Aerospace Series, 2001
% 
% translation in Matlab by Christophe Airiau, 2015
% free use of the matlab code, under J. Katz agreement.
%%*******************************************************
%   return velocity to global ref. frame
u=ul*cos(th(j))-wl*sin(th(j));
w=ul*sin(th(j))+wl*cos(th(j));
%    a(i,j) is the component of velocity induced in the
%    direction normal to panel i by panel j at the ith 
%    colocation point

%   a(i,j) is the influence coeff. defined by the 
%   tangency condition.  b(i,j) is the induced local 
%   tangential velocity to be used in cp calculation.
a=-u*sin(th(i))+w*cos(th(i));
b=u*cos(th(i))+w*sin(th(i));
% We found a= ul, b=wl ... bizarre ?
end
