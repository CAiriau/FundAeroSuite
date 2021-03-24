function [V]=wake(Pt,itloc,qw,gam)
%  calculates semi wake induced velocity at point (x,y,z) 
%  at t=itloc*dt, in the inertial frame of reference
%
% Args:
%   Pt  : point coordinates
%   qw  : panel points coordinates
%   gam : panel circulation
%   itloc : index of wortex location 
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
global j_panel
V=zeros(1,3);
i1=itloc-1;
for i=1:i1
    for j=1:j_panel
    [V1]=vortex(Pt,qw(i  ,j,  1:3),qw(i+1,j,  1:3),gam(i,j));
    [V2]=vortex(Pt,qw(i+1,j,  1:3),qw(i+1,j+1,1:3),gam(i,j));
    [V3]=vortex(Pt,qw(i+1,j+1,1:3),qw(i  ,j+1,1:3),gam(i,j));
    [V4]=vortex(Pt,qw(i  ,j+1,1:3),qw(i  ,j,  1:3),gam(i,j));
    V=V+V1+V2+V3+V4;
    end
end
