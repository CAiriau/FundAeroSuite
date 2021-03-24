function [sc,Cp,Cl]=Linear_Strength_Doublet_Potential(alpha,x,y)
%    program no. 9: linear strength doublet potential
%       ---------------------------------------------
%
%    this program finds the pressure distribution on an arbitrary airfoil 
%    by representing the surface as a finite number of doublet panels with 
%    linear strength (dirichlet b.c., program by steven yon, 1989).

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

        %    compute the potential components as   functions of r, th
        if i==j
            ph1=-0.5*(x/x2-1);
            ph2=0.5*x/x2;
        else
            ph1=1/(2*pi)*(x/x2*(th2-th1)+z/x2*log(r2/r1)-(th2-th1));
            ph2=-1/(2*pi)*(x/x2*(th2-th1)+z/x2*log(r2/r1));
        end

        %    compute the coefficients in the influence matrix
        if j==1 
            a(i,1)=ph1;
            holda=ph2;
        elseif j==m
            a(i,m)=holda+ph1;
            a(i,m+1)=ph2;
        else
            a(i,j)=holda+ph1;
            holda=ph2;
        end
    end
    %    add influence of wake as a(m+2)
    xw=co(i,1)-pt2(m,1);
    zw=co(i,2)-pt2(m,2);
    dthw=-atan(zw/xw);
    a(i,m+2)=-1/(2*pi)*dthw;
    a(i,m+3)=co(i,1)*cos(al)+co(i,2)*sin(al);

end
%    add the doublet gradient condition 
a(m+1,1)=-1;
a(m+1,2)=1;
a(m+1,m)=1;
a(m+1,m+1)=-1;

%    add the kutta condition
a(m+2,1)=-1;
a(m+2,m+1)=1;
a(m+2,m+2)=-1;
%     solve for the solution vector of doublet strengths

RHS=a(:,m+3);
gam=a(:,1:m+2)\RHS;

%    convert doublet strengths into tangential 
%    velocities along the airfoil surface and 
%    cp's on each of the panels.

for i=1:m-1
    r=(dl(i)+dl(i+1))/2;
    t1=(gam(i)+gam(i+1))/2;
    t2=(gam(i+1)+gam(i+2))/2;
    vel(i)=(t2-t1)/r;
    Cp(i)=1-vel(i)^2;
    sc(i)=pt2(i,1);
end

Cl= gam(m+2);
plot_Cp(sc,Cp);
fprintf('lift coefficient  Cl= %6.4f\n',Cl);
end


