function [sc,Cp,Cl]= Constant_Strength_Vortex(alpha,x,y)
%    program no. 4: constant strength vortex
%     ---------------------------------------
%
%    this program finds the pressure distribution on an arbitrary airfoil 
%    by representing the surface as a finite number of vortex panels with 
%    const. strength (neumann b.c., program by steven yon, 1989).

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
        %    save panel lengths for later use
        if i==1
            dl(j)=x2;
        end 
        if i==j 
          ul=0.5; wl=0;
        else
          ul=1/(2*pi)*(th2-th1);
          wl=1/(2*pi)*log(r2/r1);
        end
        [a(i,j),b(i,j)]=rotation(i,j,ul,wl,th);
    end 
    a(i,n)=cos(al)*sin(th(i))-sin(al)*cos(th(i));
end 

%    replace equation m/4 with a kutta condition

a(m/4,:)=0;
a(m/4,1)=1;
a(m/4,m)=1;
%    solve for the solution vector of vortex strengths

%   solve for the solution vector of doublet strengths
RHS=a(:,n);
gam=a(:,1:n-1)\RHS;
%    convert source strengths into tangential 
%    velocities along the airfoil surface and cp's 
%    on each of the panels

Cl=0;
temp=b*gam;
vel=temp'+cos(al)*cos(th)+sin(al)*sin(th);
Cl=sum(vel.*dl);

ans=input('smooth the velocity distribution yes=1 \n ?');

if ans == 1 
    for i=2:m
        Cp(i)=1-((vel(i)+vel(i-1))/2)^2;
        sc(i)=pt2(i-1,1);
    end
    sc(1)=1;Cp(1)=1;
else
    Cp(1:m)=1-vel(1:m).^2;
    sc=co(:,1);             % x coordinates
end 
plot_Cp(sc,Cp);
fprintf('lift coefficient  Cl= %6.4f\n',Cl);

end


