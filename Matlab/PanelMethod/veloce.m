function [V]=veloce(Pt_in,itloc,ch,Sign,alpha,sx,sz,qf,qw,sno,cso,gam)
%     subroutine veloce calculates induced velocities due to the wing
%     and its wakes in a point (x,y,z) given in the inertial frame of
%     reference.

%%*******************************************************
% Original coding obtained from J. Katz, in fortran 77
% and found in :
% Low-Speed Aerodynamics,Joseph Katz et Allen Plotkin
% Second Edition,Cambridge Aerospace Series, 2001
% 
% translation in Matlab by Christophe Airiau, 2015
% free use of the matlab code, under J. Katz agreement.
%%*******************************************************
%
global j_panel
Pt=reshape(Pt_in,1,3);
Cos=cos(alpha);Sin=sin(alpha);
x1=(Pt(1)-sx)*Cos+(Pt(3)-sz)*Sin;
y1=Pt(2);
z1=-(Pt(1)-sx)*Sin+(Pt(3)-sz)*Cos;
[V1]=wake(Pt,itloc,qw,gam);
[V2]=wake([Pt(1) -Pt(2) Pt(3)],itloc,qw,gam);
[V3]=wing0([x1 y1 z1]  ,gam,qf,sno,cso,Sign);
[V4]=wing0([x1 -y1 z1] ,gam,qf,sno,cso,Sign);
u33=Cos*(V3(1)+V4(1))-Sin*(V3(3)+V4(3));
w33=Sin*(V3(1)+V4(1))+Cos*(V3(3)+V4(3));
%     influence of mirror image
if ch >100.0
    V5=zeros(1,3); V6=zeros(1,3); V7=zeros(1,3); V8=zeros(1,3); u77=0;w77=0;
else
    x2=(Pt(1)-sx)*Cos+(-Pt(3)-sz)*Sin;
    z2=-(Pt(1)-sx)*Sin+(-Pt(3)-sz)*Cos;
    [V5]=wake([Pt(1) Pt(2) -Pt(3)],itloc,qw,gam);
    [V6]=wake([Pt(1) -Pt(2) -Pt(3)],itloc,qw,gam);
    [V7]=wing0([x2 y1 z2] ,gam,qf,sno,cso,Sign);
    [V8]=wing0([x2 -y1 z2],gam,qf,sno,cso,Sign);
    u77=Cos*(V7(1)+V8(1))-Sin*(V7(3)+V8(3));
    w77=Sin*(V7(1)+V8(1))+Cos*(V7(3)+V8(3));
end
%     velocities measured in inertial frame                             
V(1)=V1(1)+V2(1)+u33  +V5(1)   +V6(1)   +u77;
V(2)=V1(2)-V2(2)+V3(2)+V5(2)-V4(2)-V6(2)+V7(2)-V8(2);
V(3)=V1(3)+V2(3)+w33  -V5(3)   -V6(3)   -w77;
end
