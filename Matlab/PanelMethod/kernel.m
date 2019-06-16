function [x,z,x2,z2,r1,r2,th1,th2]=kernel(i,j,co,pt1,pt2,th)
% kernel to calculate influence coefficients
% convert colocation point to local panel coords.

%%*******************************************************
% Original coding obtained from J. Katz, in fortran 77
% and found in :
% Low-Speed Aerodynamics,Joseph Katz et Allen Plotkin
% Second Edition,Cambridge Aerospace Series, 2001
% 
% translation in Matlab by Christophe Airiau, 2015
% free use of the matlab code, under J. Katz agreement.
%%*******************************************************
xt=co(i,1)-pt1(j,1);
zt=co(i,2)-pt1(j,2);
x2t=pt2(j,1)-pt1(j,1);
z2t=pt2(j,2)-pt1(j,2);
Cos=cos(th(j));Sin=sin(th(j));
x=xt*Cos+zt*Sin;
z=-xt*Sin+zt*Cos;
x2=x2t*Cos+z2t*Sin;
%z2=-x2t*Sin+z2t*Cos;  verified : it is = 0
%fprintf('z2= %f\n',z2)
z2=0;
%   find r1, r2, th1, th2
r1=sqrt(x^2+z^2);
r2=sqrt((x-x2)^2+z^2);
th1=atan2(z,x);
th2=atan2(z,x-x2);
end

