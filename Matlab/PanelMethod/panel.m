function [Normal,A]=panel(Pt1,Pt2,Pt3,Pt4)
%  calculation of panel area  = A and normal vector =Normal.


%%*******************************************************
% Original coding obtained from J. Katz, in fortran 77
% and found in :
% Low-Speed Aerodynamics,Joseph Katz et Allen Plotkin
% Second Edition,Cambridge Aerospace Series, 2001
% 
% translation in Matlab by Christophe Airiau, 2015
% free use of the matlab code, under J. Katz agreement.
%%*******************************************************

rot=cross(Pt2-Pt3,Pt4-Pt1);    % vec M_3M_2 x vec M_1M_4
norm=sqrt(dot(rot,rot));
Normal=rot/norm;        % normal unit vector

% Panel area :

rot1=cross(Pt2-Pt1,Pt4-Pt1);    % vec M_1M_2 x vec M_1M_4
rot2=cross(Pt4-Pt1,Pt3-Pt1);    % vec M_1M_4 x vec M_1M_3
A=0.5*(sqrt(dot(rot1,rot1))+sqrt(dot(rot2,rot2)));
end
