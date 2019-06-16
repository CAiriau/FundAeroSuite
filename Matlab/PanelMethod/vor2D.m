function [u,w]=vor2D(X,Z,Gamma)
% calculates influence of vortex
%

%%*******************************************************
% Original coding obtained from J. Katz, in fortran 77
% and found in :
% Low-Speed Aerodynamics,Joseph Katz et Allen Plotkin
% Second Edition,Cambridge Aerospace Series, 2001
% 
% translation in Matlab by Christophe Airiau, 2015
% free use of the matlab code, under J. Katz agreement.
%%*******************************************************
R=sqrt(X^2+Z^2);
if R >= 0.001
    v=Gamma/(2*pi*R);
    u= v*Z/R;
    w=-v*X/R;
else
    u=0;w=0;
end
end
