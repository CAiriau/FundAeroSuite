function [n,m,pt1,pt2,co,th,al]=initialization(x,y,alpha)
% initialization for all functions
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
m=length(x)-1;                  %  number of panels
n=m+1; 
ept(:,1)=x; ept(:,2)=y;         % profile coordinates
if nargin == 3
    % alpha = angle of attack in degrees
    al=alpha/180.*pi;
    fprintf('angle of attack = %4.2f °\n',alpha);
else
    disp('no angle of attack required');
    al=0.;
end

% convert paneling to clockwise
for i=1:m+1
     ep(i,:)=ept(n-i+1,:);
end
% coordinates of panel end points
for i=1:m
     pt1(i,:)=ep(i,:);
     pt2(i,:)=ep(i+1,:);
end 
% find panel angles th(j)
for i=1:m
     th(i)=atan2(pt2(i,2)-pt1(i,2),pt2(i,1)-pt1(i,1));
end
% colocation points
for i=1:m
     co(i,:)=(pt2(i,:)-pt1(i,:))/2+pt1(i,:);
end

fid = fopen('coordinates.dat','w');
fprintf(fid,'# Pt1, Pt2, Co, i \n');
for i = 1:m
    fprintf(fid,'%10.3f %10.3f \t | \t %10.3f %10.3f \t | \t %10.3f %10.3f  \t|\t %10.3f  %d \n',pt1(i,1:2),pt2(i,1:2),co(i,1:2),th(i)/pi*180,i);
    fprintf('%10.4f %10.4f \t | \t %10.4f %10.4f \t | \t %10.4f %10.4f  \t|\t %10.4f  %d \n',pt1(i,1:2),pt2(i,1:2),co(i,1:2),th(i)/pi*180,i);
end
fclose(fid);

disp('End initialization')
end

