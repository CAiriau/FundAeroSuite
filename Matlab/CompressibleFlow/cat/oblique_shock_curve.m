function oblique_shock_curve(Mach)
% Shock curves for a given Mach number 
% C Airiau, avril 2012

global gam
close all;
disp('There are many ways to plot the shock curves');
coef=180/pi;
n=500;
Sigma_min=asin(1/Mach); Sigma_max=pi/2;
dSigma=(Sigma_max-Sigma_min)/(n-1);
Sigma=Sigma_min:dSigma:Sigma_max;

Mn=Mach*sin(Sigma);
f=1./Rho2overRho1(Mn);
u=cos(Sigma).^2+f.*sin(Sigma).^2;
v=sin(Sigma).*cos(Sigma).*(1-f);
theta=atan(v./u);
fprintf('maximal deviation angle = %10f °\n',max(theta)*coef);
figure();
h1=plot(theta*coef,Sigma*coef);
title('Shock curve (method 1)');
xlabel('Wall angle deviation \theta (°)');
ylabel('Shock angle (°)');


% Second possibility
inv_tan_theta=((gam+1)/2*Mach^2./ (Mach^2*sin(Sigma).^2-1) -1 ) .*tan(Sigma);
theta1=atan(1./inv_tan_theta);
figure();
h1=plot(theta1*coef,Sigma*coef,'ro',theta*coef,Sigma*coef,'k-');
title('Shock curve (method 2)');
xlabel('Wall angle deviation \theta (°)');
ylabel('Shock angle (°)');
legend('met. 2','met. 1');

% Third  possibility
tan_theta=f.*tan(Sigma);
theta2=Sigma-atan(tan_theta);
figure();
h1=plot(theta2*coef,Sigma*coef,'ro',theta*coef,Sigma*coef,'k-');
title('Shock curve (method 3)');
xlabel('Wall angle deviation \theta (°)');
ylabel('Shock angle (°)');
legend('met. 3','met. 1')


figure();
Mn2=downstream_normal_Mach(Mn);
M2=Mn2./sin(Sigma-theta);
h2=plot(u,v);
title(['hodograph plane : BUSEMANN polar, M_{upstream}=' num2str(Mach)]);
xlabel('u');ylabel('v');
for i=1:50:500
    text(u(i),v(i),num2str(theta(i)*coef));
end
text(0.45,0.1,'Deviation angle  \theta in °');

figure();
plot(theta*coef,u,theta*coef,v)
title('u and v function of \theta, behind the shock');
xlabel('\theta (°)');
ylabel('u,v');
legend('u','v')

figure()
plot(theta*coef,M2);
title('Downstream Mach as a function of \theta');
xlabel('\theta (°)');
ylabel('Downstream Mach Number');
[val, i]=min(abs(M2-1));
theta_1=theta(i)*coef;
fprintf('Deviation angle for a downstream Mach number equal to 1 = %10f ° \n',theta_1);

end

