function [sc,Cp,Cl]=Quadratic_Strength_Doublet_Potential(alpha,x,y,x_internal)
%     program no. 10: quadratic strength doublet potential
%      ----------------------------------------------------
%
%     this program finds the pressure distribution on an arbitrary airfoil 
%      by representing the surface as a finite number of doublet panels with 
%      quadratic strength (dirichlet b.c., program by steven yon, 1989).


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

%     establish location of additional colocation point
co(m+1,1)=x_internal;
co(m+1,2)=0.0;
 
%     establish influence coefficients
for i=1:m+1
    for j=1:m
         [x,z,x2,z2,r1,r2,th1,th2]=kernel(i,j,co,pt1,pt2,th);
         % save panel lengths for lift coeff. calc.
         if i==1
             dl(j)=x2;
         end 

        %     compute the influence coefficients in 
        %     the 'b' matrix (unreduced).
    if i==j
        b(i,j,1)=0.5;
        b(i,j,2)=0.5*x;
        b(i,j,3)=0.5*x^2;
    else
        b(i,j,1)=-1/(2*pi)*(th2-th1);
        b(i,j,2)=-1/(2*pi)*(x*(th2-th1)+z*log(r2/r1));
        b(i,j,3)=1/(2*pi)*((x^2-z^2)*(th1-th2)-2*x*z*log(r2/r1)-z*x2);
    end
    end
end

%     add doublet gradient condition
b(m+2,1,1)=0;
b(m+2,1,2)=1;
b(m+2,1,3)=0;
b(m+2,m,1)=0;
b(m+2,m,2)=1;
b(m+2,m,3)=2*dl(m);

%     add kutta condition
b(m+3,1,1)=-1;
b(m+3,1,2)=0;
b(m+3,1,3)=0;
b(m+3,m,1)=1;
b(m+3,m,2)=dl(m);
b(m+3,m,3)=dl(m)^2;

%     back substitute the 'b' matrix with the 
%     regression formula to get the components
%     of the 'a' matrix.
for i=1:m+3
    for j=m-1:-1:1
        b(i,j,1)=b(i,j,1)+b(i,j+1,1);
        b(i,j,2)=b(i,j,2)+b(i,j+1,1)*dl(j)+b(i,j+1,2);
        b(i,j,3)=b(i,j,3)+b(i,j+1,1)*dl(j)^2+2*b(i,j+1,2)*dl(j);
    end
    a(i,1)=b(i,1,1);
    a(i,2)=b(i,1,2);
    for  j=1:m
        a(i,j+2)=b(i,j,3);
    end
end
%
%    add influence of wake as a(i,m+3)
%    and rhs as a(i,m+4)
for i=1:m+1
    xw=co(i,1)-pt2(m,1);
    zw=co(i,2)-pt2(m,2);
    dthw=-atan(zw/xw);

    a(i,m+3)=-1/(2*pi)*dthw;
    a(i,m+4)=(co(i,1)*cos(al)+co(i,2)*sin(al));

end

%     complete kutta cond. by adding wake coeff and rhs
%     to rows m+2 and m+3
a(m+2,m+3)=0;
a(m+2,m+4)=0;
a(m+3,m+3)=-1;
a(m+3,m+4)=0;
%     solve for the solution vector of doublet strengths
RHS=a(:,m+4);
gam=a(:,1:m+3)\RHS;

%     convert doublet strengths, linear correction
%     and quadratic correction into tangential 
%     velocities along the airfoil surface and cp's 
%     on each of the panels.

%
%     forward substitute using the solution vector to 
%     get the doublet strength parameters 
%     for each panel.
u1(1)=gam(1);
a1(1)=gam(2);
for i=3:m+2
    b1(i-2)=gam(i);
end

for i=1:m-1
    u1(i+1)=u1(i)+a1(i)*dl(i)+b1(i)*(dl(i))^2;
    a1(i+1)=a1(i)+2*b1(i)*dl(i);
end

%     the derivative of the doublet strength is the 
%     surface speed along each panel.

vel(1:m)=a1(1:m)+b1(1:m).*dl(1:m);
Cp=1-vel(1:m).^2;
sc=co(1:m,1);

Cl= gam(m+3);
plot_Cp(sc,Cp);
fprintf('lift coefficient  Cl= %6.4f\n',Cl);
end

