function [sc,Cp,Cl]= Constant_Strength_Doublet(alpha,x,y)
%    program no. 3: constant strength doublet
%       ----------------------------------------
%    this program finds the pressure distribution on an arbitrary airfoil 
%       by representing the surface as a finite number of doublet panels with 
%       const. strength (neumann b.c., program by steven yon, 1989).
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



[n,m,pt1,pt2,co,th,al]=initialization(x,y,alpha);
a=zeros(n,n+1);


% influence coefficients
for i=1:m
     for j=1:m
         [x,z,x2,z2,r1,r2,th1,th2]=kernel(i,j,co,pt1,pt2,th);
         %    compute the velocity induced at the ith 
         %    colocation point by the jth panel
          if i==j 
            ul=0; wl=-1/(pi*x);
          else
            ul=1/(2*pi)*(z/(r1^2)-z/(r2^2));
            wl=-1/(2*pi)*(x/(r1^2)-(x-x2)/(r2^2));
          end
         [a(i,j),b(i,j)]=rotation(i,j,ul,wl,th);
    end 
    %    include the influence of the wake panel
   r=sqrt((co(i,1)-pt2(m,1))^2+(co(i,2)-pt2(m,2))^2);
   u=1/(2*pi)*(co(i,2)/(r^2));
   w=-1/(2*pi)*(co(i,1)-pt2(m,1))/(r^2);

    
    a(i,n)=-u*sin(th(i))+w*cos(th(i));
    b(i,n)=u*cos(th(i))+w*sin(th(i));
    a(i,n+1)=cos(al)*sin(th(i))-sin(al)*cos(th(i));
end 

%   prepare the matrix for solution by providing 
%   a kutta condition
a(n,1:n+1)=0;
a(n,1)=-1;
a(n,m)=1;
a(n,n)=-1;
%   solve for the solution vector of doublet strengths
RHS=a(:,n+1);
gam=a(:,1:n)\RHS;

%    convert doublet strengths into tangential 
%    velocities along the airfoil surface and cp's 
%    on each of the panels
temp=b*gam;
for i=1:m
    if (i ~= 1) && (i ~= m) 
        r=sqrt((co(i+1,1)-co(i-1,1))^2+(co(i+1,2)-co(i-1,2))^2);
        vloc(i)=(gam(i+1)-gam(i-1))/r;
    elseif i == 1 
        r=sqrt((co(2,1)-co(1,1))^2 +(co(2,2)-co(1,2))^2);
        vloc(i)=(gam(2)-gam(1))/r;
    elseif i == m
        r=sqrt((co(m,1)-co(m-1,1))^2+(co(m,2)-co(m-1,2))^2);
        vloc(i)=(gam(m)-gam(m-1))/r;
    end
end 
vel=cos(al)*cos(th)+sin(al)*sin(th)+temp'+vloc/2;
Cp=1-vel.^2;
sc=co(:,1);             % x coordinates
plot_Cp(sc,Cp);
Cl= gam(m+1);
fprintf('lift coefficient  Cl= %6.4f\n',Cl);
end


