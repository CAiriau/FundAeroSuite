function [u,w]=downwash(x,z,i1,i2,vortic)
%   calculates downwash induced by it-1 wake vortices
%
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
u=0;w=0;
for i=i1:i2
    [u1,w1]=vor2D(x-vortic(i,1),z-vortic(i,2),vortic(i,3));
    u=u+u1;w=w+w1;
end

end
