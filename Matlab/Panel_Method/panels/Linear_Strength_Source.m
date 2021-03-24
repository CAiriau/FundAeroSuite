function [sc,Cp,Cl]= Linear_Strength_Source(alpha,x,y,x_wake)
%     program no. 5: linear strength source
%       -------------------------------------
%
%     this program finds the pressure distribution on an arbitrary airfoil 
%     by representing the surface as a finite number of source panels with 
%     linear strength (alpha=0, neumann b.c., program by steven yon, 1989).

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

th(m+1)=0;
co(m+1,1)=x_wake;
co(m+1,2)=0;
%     establish influence coefficients
for i=1:m+1
    for j=1:m
        [x,z,x2,z2,r1,r2,th1,th2]=kernel(i,j,co,pt1,pt2,th);
        %     compute velocity components as functions of 
        %     sigma1 and sigma2. these velocities are in 
        %     the jth reference frame.
        if i==j 
            u1l=1/(2*pi);
            u2l=-1/(2*pi);
            w1l=-0.5*(x-x2)/x2;
            w2l= 0.5*x/x2;
        else
            w1l=-(z*log(r2/r1)+x*(th2-th1)-x2*(th2-th1))/(2*pi*x2);
            w2l=(z*log(r2/r1)+x*(th2-th1))/(2*pi*x2);
            u1l=((x2-z*(th2-th1))-x*log(r1/r2)+x2*log(r1/r2))/(2*pi*x2);
            u2l=-((x2-z*(th2-th1))-x*log(r1/r2))/(2*pi*x2);
        end
        %     transform the local velocities into the global 
        %     reference frame.

        u1=u1l*cos(th(j))-w1l*sin(th(j));
        u2=u2l*cos(th(j))-w2l*sin(th(j));
        w1=u1l*sin(th(j))+w1l*cos(th(j));
        w2=u2l*sin(th(j))+w2l*cos(th(j));
        %     compute the coefficients of sigma in the 
        %     influence matrix
        if j==1
            a(i,1)=-u1*sin(th(i))+w1*cos(th(i));
            holda=-u2*sin(th(i))+w2*cos(th(i));
            b(i,1)=u1*cos(th(i))+w1*sin(th(i));
            holdb=u2*cos(th(i))+w2*sin(th(i));
        elseif  j==m
            a(i,m)=-u1*sin(th(i))+w1*cos(th(i))+holda;
            a(i,n)=-u2*sin(th(i))+w2*cos(th(i));
            b(i,m)=u1*cos(th(i))+w1*sin(th(i))+holdb;
            b(i,n)=u2*cos(th(i))+w2*sin(th(i));
        else
            a(i,j)=-u1*sin(th(i))+w1*cos(th(i))+holda;
            holda=-u2*sin(th(i))+w2*cos(th(i));
            b(i,j)=u1*cos(th(i))+w1*sin(th(i))+holdb;
            holdb=u2*cos(th(i))+w2*sin(th(i));
        end
    end
    a(i,n+1)=sin(th(i));
end
%     solve for the solution vector of source strengths
RHS=a(:,n+1);
gam=a(:,1:n)\RHS;
%     convert source strengths into tangential 
%     velocities along the airfoil surface and cp's 
%     on each of the panels.
vel=b*gam;
v=vel'+cos(al)*cos(th)+sin(al)*sin(th);

ans=input('smooth the velocity distribution yes=1 \n ?');

if ans == 1 
    for i=2:m
        Cp(i)=1-((v(i)+v(i-1))/2)^2;
        sc(i)=pt1(i,1);
    end
    sc(1)=1;Cp(1)=1;
else
    Cp(1:m)=1-v(1:m).^2;
    Cp=Cp';
    sc=co(1:m,1);             % x coordinates
end 
Cl=0;
plot_Cp(sc,Cp);
fprintf('lift coefficient  Cl= %6.4f\n',Cl);
end 


