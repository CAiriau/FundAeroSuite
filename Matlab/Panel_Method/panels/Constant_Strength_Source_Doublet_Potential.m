function [sc,Cp,Cl]=Constant_Strength_Source_Doublet_Potential(alpha,x,y)
%     program no. 8: constant strength source/doublet potential
%       ---------------------------------------------------------
%
%     this program finds the pressure distribution on an arbitrary airfoil 
%     by representing the surface as a finite number of source/doublet panels
%     with constant strength (dirichlet b.c., program by steven yon, 1989).

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
fprintf('i=1, x,y= %f,%f i=N, x,y=  %f,%f\n',x(2),y(2),x(n-1),y(n-1));
%     establish source strengths (sigma=v dot n)
for i=1:m
    sigma(i)=cos(al)*sin(th(i))-sin(al)*cos(th(i));
end 
fprintf('sum sigma_j = %f \n',sum(sigma))
%     establish influence coefficients
for i=1:m
    temp=0;
    for j=1:m
        [x,z,x2,z2,r1,r2,th1,th2]=kernel(i,j,co,pt1,pt2,th);
        % save panel lengths for lift coeff. calc.
        if i==1
            dl(j)=x2;
        end 

        %     compute the doublet influence coefficients
        if i==j
            a(i,j)=0.5;
        else
            a(i,j)=-1/(2*pi)*(th2-th1);
        end

        %     compute the source influence coeff's and add 
        %     them up to give the rhs
        if i==j
            temp=temp+sigma(j)/pi*x*log(r1);
            utmp=1/pi*x*log(r1);
            fprintf('xk = %f \t rk = %f  \t u=%f \t i= %d\n',x,r1,utmp,i);
        else
            temp=temp+sigma(j)/(2*pi)*(x*log(r1) -(x-x2)*log(r2)+z*(th2-th1));
        end
    end
%   add wake influence
    xw=co(i,1)-pt2(m,1);
    zw=co(i,2)-pt2(m,2);
    dthw=-atan(zw/xw);
    a(i,n)=-1/(2*pi)*dthw;
    a(i,n+1)=temp;
end

fprintf('Pt wake 2= %f %f\n',pt2(m,:))
fprintf('Pt wake 1= %f %f\n',pt1(m,:))
%     add an explicit kutta condition
a(n,1:n+1)=0;
a(n,1)=-1;
a(n,m)=1;
a(n,n)=-1;
%     solve for the solution vector of doublet strengths
disp('RHS')
RHS=a(:,n+1)
disp('gamma')
gam=a(:,1:n)\RHS
fprintf('size of mu = %d \n',length(gam))
size(co)

%     convert doublet strengths into tangential 
%     velocities along the airfoil surface and cp's 
%     on each panel.
phi=co(:,1)*cos(al)+co(:,2)*sin(al)+gam(1:m);
for  i=1:m-1
    r=(dl(i+1)+dl(i))/2;
    vel(i)=(phi(i)-phi(i+1))/r;
    Cp(i)=1-vel(i)^2;
    sc(i)=pt2(i,1);
end

Cl= gam(m+1);
plot_Cp(sc,Cp);
fprintf('lift coefficient  Cl= %6.4f\n',Cl);
save_data('Constant_Source_Doublet_Potential',sc,Cp);
plot_Source(co(:,1),sigma);

h=figure();
plot(co(:,1),gam(1:m),'r-')
title('doublet distribution');
xlabel('x/c');
ylabel('mu');
size(co)
size(gam)
save_data('Doublet_distribution',co(:,1),gam(1:m));
end


