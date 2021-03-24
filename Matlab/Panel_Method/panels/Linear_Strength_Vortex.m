function [sc,Cp,Cl]=Linear_Strength_Vortex(alpha,x,y)
%     program no. 6: linear strength vortex
%         -------------------------------------
%
%     this program finds the pressure distribution on an arbitrary airfoil 
%     by representing the surface as a finite number of vortex panels with 
%     linear strength (neumann b.c., program by steven yon, 1989).

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

          %     compute velocity componants as functions of 
          %     gamma1 and gamma2. these velocities are in 
          %     the jth reference frame.
          if i==j 
               u1l=-0.5*(x-x2)/x2;
               u2l=0.5*x/x2;
               w1l=-1/(2*pi);
               w2l=1/(2*pi);
          else
               u1l=-(z*log(r2/r1)+x*(th2-th1)-x2*(th2-th1))/(2*pi*x2);
               u2l=(z*log(r2/r1)+x*(th2-th1))/(2*pi*x2);
               w1l=-((x2-z*(th2-th1))-x*log(r1/r2)+x2*log(r1/r2))/(2*pi*x2);
               w2l=((x2-z*(th2-th1))-x*log(r1/r2))/(2*pi*x2);
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
     a(i,n+1)=cos(al)*sin(th(i))-sin(al)*cos(th(i));
end
%     add the kutta condition
a(n,1)=1;
a(n,n)=1;

%     solve for the solution vector of source strengths
RHS=a(:,n+1);
gam=a(:,1:n)\RHS;

%     convert vortex strengths into tangential 
%     velocities along the airfoil surface and cp's 
%     on each of the panels.

Cl=0;
vel=b*gam;
v=vel'+cos(al)*cos(th)+sin(al)*sin(th);
Cl=sum(v.*dl);
Cp=1-v.^2;
sc=co(:,1);             % x coordinates
plot_Cp(sc,Cp);
fprintf('lift coefficient  Cl= %6.4f\n',Cl);

end



