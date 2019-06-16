function [out]=Discrete_Vortex_Method(n)
%	program no. 15: discrete vortex method (thin wing, elliptic camber)
%       -------------------------------------------------------------------
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
global name

% ATTENTION : A LA FIN ce n'est pas le Cp mais 2x Cp  = dL/dx/q

%     input data
%

%
% ==========
% input data
% ==========
%
line='====================================';
%n=10;                   % number of panels
c=1;                    % root chord
epsilon=0.1*c;
alpha_1=0;             % Angle of attack
rho=1;                  % density
Vinf=1;                 % velocity m/s
alpha=pi/180*alpha_1;
uinf=cos(alpha)*Vinf;
winf=sin(alpha)*Vinf;
q=0.5*rho*Vinf^2;



% ===========================
% grid generation (n panels)
% ===========================

dx=c/n;
for i=1:n
% collocation point
    xc(i) = c/n*(i-0.25);
    zc(i) = 4.*epsilon*xc(i)/c*(1.-xc(i)/c);
% vortex point
    x(i)  = c/n*(i-0.75);
    z(i)  = 4.*epsilon*x(i)/c*(1.-x(i)/c);
% normal at collocation point; n=(enx,enz)
    detadx=4.*epsilon/c*(1.-2.*xc(i)/c);
    sq=sqrt(1+detadx^2);
    enx(i)= -detadx/sq;
    enz(i)=  1./sq;
end

% plot the profile
xp=[0:0.01:1];
zp=4.*epsilon.*xp.*(1-xp);
figure()
plot(xp,zp)
axis equal;
axis equal;
title('Profile');
 
% ======================
% influence coefficients
% ======================
for i=1:n
    for j=1:n
      [u,w]=vor2D(xc(i)-x(j),zc(i)-z(j),1.0);
      a(i,j)=u*enx(i)+w*enz(i);
    end
%   the rhs vector 
    rhs(i)=-uinf*enx(i)-winf*enz(i);
end


% ===============================================                                                                    
% solution of the problem: rhs(i)=a(i,j)*gamma(i)
% =============================================== 
%fprintf('size %s : %d %d \n','a',size(a));
%fprintf('size %s : %d %d \n','rhs',size(rhs));
gam=a\rhs';
%
figure()
plot(x/c,gam);
title('Circulation versus x/c');
% ==================
% aerodynamic loads
% ==================

bl=0.0;
dL=rho*Vinf*gam;
Cp=dL./(dx*q);     % lower surface Cp^-
temp=32*epsilon/c.*sqrt(x/c.*(1-x/c));
% Cp_ref is the analytic solution
Cp_ref=4.*sqrt((c-x)./x)*alpha+temp;
Lift=sum(dL);
Cl=Lift/(q*c);
Cl_ref=2.*pi*(alpha+2*epsilon/c);
%     Cl_ref, Cp_ref - are the exact solutions

fprintf('Lift                       = %12.8f \n',Lift);
fprintf('Cl                         = %12.8f \n',Cl);
fprintf('Cl   (exact)               = %12.8f \n',Cl_ref);

figure()
plot(x,Cp_ref,'k-',x,Cp,'b--');
title([ 'Cp : ' name]);
legend('Exact','Panel','Location','northeast');

end
