function [sc,Cp,Cl]=Constant_Strength_Source(x,y)
%    program no. 2: constant strength source
%     -------------------------------------
%    
%    this program finds the pressure distribution on an arbitrary airfoil 
%    by representing the surface as a finite number of source panels with 
%    const. strength (alpha=0, neumann b.c., program by steven yon, 1989).

%%*******************************************************
% Original coding obtained from J. Katz, in fortran 77
% and found in :
% Low-Speed Aerodynamics,Joseph Katz et Allen Plotkin
% Second Edition,Cambridge Aerospace Series, 2001
% 
% translation in Matlab by Christophe Airiau, 2015
% free use of the matlab code, under J. Katz agreement.
%%*******************************************************

[n,m,pt1,pt2,co,th]=initialization(x,y);

% influence coefficients
for i=1:m
    for j=1:m
        [x,z,x2,z2,r1,r2,th1,th2]=kernel(i,j,co,pt1,pt2,th);
        % compute velocity in local referece frame
        if i==j 
          ul=0; wl=0.5;
        else
          ul=1/(2*pi)*log(r1/r2);
          wl=1/(2*pi)*(th2-th1);
        end 
        [a(i,j),b(i,j)]=rotation(i,j,ul,wl,th);
    end 
    a(i,n)=sin(th(i));
end
%   solve for the solution vector of source strengths
RHS=a(:,n);
gam=a(:,1:n-1)\RHS;

%   convert source strengths into tangential
%   velocities along the airfoil surface and cp's
%   on each of the panels
  
vel=b*gam;
Cp=1-(vel+cos(th')).^2; % pressure coefficient
sc=co(:,1);             % x coordinates
Cl=0;                   % lift coefficient
plot_Cp(sc,Cp);         % Cp symmetrical ...
disp('lift coefficient=0');
save_data('Constant_S_Source',sc,Cp)
end 


