function [sc,Cp,Cl]=Constant_Strength_Doublet_Potential(alpha,x,y)
%     program no. 7: constant strength doublet potential
%      --------------------------------------------------
%
%     this program finds the pressure distribution on an arbitrary airfoil 
%      by representing the surface as a finite number of doublet panels with 
%      constant strength (dirichlet b.c., program by steven yon, 1989).
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

%     establish influence coefficients
for i=1:m
    for j=1:m
         [x,z,x2,z2,r1,r2,th1,th2]=kernel(i,j,co,pt1,pt2,th);

         % save panel lengths for lift coeff. calc.
        if i==1
             dl(j)=x2;
        end 

        %     compute influence coefs. a(i,j)
        if i==j
            a(i,j)=0.5;
        else
            a(i,j)=-1/(2*pi)*(th2-th1);
        end
    end
    %   add wake influence
    xw=co(i,1)-pt2(m,1);
    zw=co(i,2)-pt2(m,2);
    dthw=-atan(zw/xw);
    a(i,n)=-1/(2*pi)*dthw;
    a(i,n+1)=co(i,1)*cos(al)+co(i,2)*sin(al);
end
%     add an explicit kutta condition
a(n,1)=-1;
a(n,m)=1;
a(n,n)=-1;

%     solve for the solution vector of doublet strengths

RHS=a(:,n+1);
gam=a(:,1:n)\RHS;

%     convert doublet strengths into tangential 
%     velocities along the airfoil surface and 
%     cp's on each of the panels.

for  i=1:m-1
  r=(dl(i)+dl(i+1))/2;
  vel(i)=(gam(i+1)-gam(i))/r;
  Cp(i)=1-vel(i)^2;
  sc(i)=pt2(i,1);
end
Cl= gam(m+1);
plot_Cp(sc,Cp);
fprintf('lift coefficient  Cl= %6.4f\n',Cl);
end


