function [V]=vortex(Pt,Pt1_in,Pt2_in,Gamma)
%        it  calculates the induced velocity (u,v,w) at a point
%        (x,y,z) due to a vortex element vith strength gamma per unit length
%        pointing to the direction (x2,y2,z2)-(x1,y1,z1).

%%*******************************************************
% Original coding obtained from J. Katz, in fortran 77
% and found in :
% Low-Speed Aerodynamics,Joseph Katz et Allen Plotkin
% Second Edition,Cambridge Aerospace Series, 2001
% 
% translation in Matlab by Christophe Airiau, 2015
% free use of the matlab code, under J. Katz agreement.
%%*******************************************************
rcut=1.0e-10;
V=zeros(1,3);
%fprintf('vortex: size Pt = %d %d  \n',size(Pt));
%        calculation of r1 x r2
Pt1=reshape(Pt1_in,1,3);
Pt2=reshape(Pt2_in,1,3);
P1=Pt(:)-Pt1(:);P2=Pt(:)-Pt2(:);
P1xP2=cross(P1,P2);
%        calculation of (r1 x r2 )**2
square=dot(P1xP2,P1xP2);
%        calculation of r0(r1/r(r1)-r2/r(r2))
r1=sqrt(dot(P1,P1));
r2=sqrt(dot(P2,P2));
if (r1 < rcut) || (r2 < rcut) || (square < rcut)
%        when point (x,y,z) lyes on vortex element; its induced velocity is
    V=[0 0 0];
else
    P12=Pt2-Pt1;
    coef=Gamma/(4.0*pi*square)*(dot(P12,P1)/r1-dot(P12,P2)/r2);
    V=coef*P1xP2';
end
%fprintf('vortex: size V = %d %d  \n',size(V));
end
