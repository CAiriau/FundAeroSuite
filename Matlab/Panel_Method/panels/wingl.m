function [V] = wingl(Pt,Gamma,qf)
%     calculates induced velocity at a point (x,y,z) due to longitudinal
%     vorticity distribution gamax(i,j) only(semi-span), in a wing fixed
%     coordinate system + (t.e. unsteady vortex).
%     ** serves for induced drag calculation only **


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
V=zeros(1,3);
for i=1:i_panel
    for j=1:j_panel
        [V2]=vortex(Pt,qf(i,j+1,1:3),qf(i+1,j+1,1:3),Gamma(i,j));
        [V4]=vortex(Pt,qf(i+1,j,1:3),qf(i,j,1:3),Gamma(i,j));
        V=V2+V4;
    end
end

%     add influence of latest unsteady wake element:
i=i_panel;
for j=1:j_panel
    [V3]=vortex(Pt,qf(i+1,j+1,1:3),qf(i+1,j,1:3),Gamma(i,j));
    V=V+V3;
end

end
