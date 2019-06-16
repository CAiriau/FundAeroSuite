function [S,c,b,AR,qf,ds,qc]=grid(alpha,b,ch,dxw,x,i_panel,j_panel)
%  x(1) - is root l.e., x(2) tip l.e., x(3) tip t.e., and x(4) is root t.e.
%  i_panel: no. of chordwise boxes,  j_panel: no. of spanwise boxes 
%  alpha: angle of attack in degrees 
%  b    : semi-span
%  ch   : ground clearance
%  dxw  : position of the wake panel
%
%  qf   : Grid points (panel)
%  qc   : Collocation points

%%*******************************************************
% Original coding obtained from J. Katz, in fortran 77
% and found in :
% Low-Speed Aerodynamics,Joseph Katz et Allen Plotkin
% Second Edition,Cambridge Aerospace Series, 2001
% 
% translation in Matlab by Christophe Airiau, 2015
% free use of the matlab code, under J. Katz agreement.
%%*******************************************************
ib1=i_panel+1;
ib2=i_panel+2;
jb1=j_panel+1;

Cos=cos(alpha/180*pi);
Sin=sin(alpha/180*pi);
%        wing fixed vortices location   ( qf(i,j,(x,y,z))...) 
dy=b/j_panel;
for j=1:jb1
    yle=dy*(j-1);
    xle=x(1)+(x(2)-x(1))*yle/b;
    xte=x(4)+(x(3)-x(4))*yle/b;
    % xle and xte are l.e. and t.e. x-coordinates
    dx=(xte-xle)/i_panel;
    for  i=1:ib1
        qf(i,j,1)=(xle+dx*(i-0.75))*Cos;
        qf(i,j,2)=yle;
        qf(i,j,3)=-qf(i,j,1)*Sin+ch;
    end
    %  wake far field points
    qf(ib2,j,1)=xte+dxw;
    qf(ib2,j,2)=qf(ib1,j,2);
    qf(ib2,j,3)=qf(ib1,j,3);
end
%
%        wing collocation points
%
for j=1:j_panel
    for i=1:i_panel
        qc(i,j,1:3)=(qf(i,j,1:3)+qf(i,j+1,1:3)+qf(i+1,j+1,1:3)+qf(i+1,j,1:3))/4;
        %
        %        computation of normal vectors
        %
        [ds(i,j,1:3),ds(i,j,4)]= panel(qf(i,j,1:3),qf(i+1,j,1:3),qf(i,j+1,1:3),qf(i+1,j+1,1:3));
    end
end
% ds(i,j,1:3) : panel normal direction
% ds(i,j,4)   : panel surface

%    b -is semi span, c -av. chord, s - area 
S=0.5*(x(3)-x(2)+x(4)-x(1))*b;  % wing  area
c=S/b;                          % wing mean chord
AR=2*b^2/S;                     % wing Aspect Ratio
end
