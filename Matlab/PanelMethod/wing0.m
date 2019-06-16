function [a1,V] = wing0(Pt,Gamma,qf,sno,cso,Sign)
%
% calculates semi wing induced velocity at a point (x,y,z) due to wing
% vorticity distribution gama(i,j)  in a wing fixed coordinate system

%%*******************************************************
% Original coding obtained from J. Katz, in fortran 77
% and found in :
% Low-Speed Aerodynamics,Joseph Katz et Allen Plotkin
% Second Edition,Cambridge Aerospace Series, 2001
% 
% translation in Matlab by Christophe Airiau, 2015
% free use of the matlab code, under J. Katz agreement.
%%*******************************************************
global i_panel j_panel
V=zeros(1,3);V0=zeros(1,3);
N=size(qf);
%   fprintf('wing0: i_panel = %d, j_panel=%d \n',i_panel,j_panel);
%   fprintf('wing0: size Gamma = %d %d \n',size(Gamma));
%  fprintf('wing0: size qf = %d %d %d \n',N);

for i=1:i_panel
    for j=1:j_panel
        [V1]=vortex(Pt,qf(i,j,:)      ,qf(i,j+1,:)    ,Gamma(i,j));
        [V2]=vortex(Pt,qf(i,j+1,1:3)  ,qf(i+1,j+1,1:3),Gamma(i,j));
        [V3]=vortex(Pt,qf(i+1,j+1,1:3),qf(i+1,j,1:3)  ,Gamma(i,j));
        [V4]=vortex(Pt,qf(i+1,j,1:3)  ,qf(i,j,1:3)    ,Gamma(i,j));
        V0=V2+V4+V1+V3;
        if Sign >= 1 
            a1(i,j)=V0(1)*sno(i)-V0(3)*cso(i);
        else
            a1(i,j)=V0(1)*sno(i)+V0(3)*cso(i);
        end
        V=V+V0;
    end
end
