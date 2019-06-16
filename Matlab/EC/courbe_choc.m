function courbe_choc(Mach)
% courbe de choc pour un Mach donné
% C Airiau, avril 2012

global gam
close all;
disp('Il y a plusieurs façon de tracer les courbes de choc');
coef=180/pi;
n=500;
Sigma_min=asin(1/Mach); Sigma_max=pi/2;
dSigma=(Sigma_max-Sigma_min)/(n-1);
Sigma=Sigma_min:dSigma:Sigma_max;

Mn=Mach*sin(Sigma);
f=1./Rho2surRho1(Mn);
u=cos(Sigma).^2+f.*sin(Sigma).^2;
v=sin(Sigma).*cos(Sigma).*(1-f);
theta=atan(v./u);
fprintf('Angle de deviation maximale = %10f °\n',max(theta)*coef);
figure();
h1=plot(theta*coef,Sigma*coef);
title('Courbe de choc (méthode 1)');
xlabel('Deviation de la paroi \theta (°)');
ylabel('Angle du choc (°)');


% Seconde possibilite
inv_tan_theta=((gam+1)/2*Mach^2./ (Mach^2*sin(Sigma).^2-1) -1 ) .*tan(Sigma);
theta1=atan(1./inv_tan_theta);
figure();
h1=plot(theta1*coef,Sigma*coef,'ro',theta*coef,Sigma*coef,'k-');
title('Courbe de choc (méthode 2)');
xlabel('Deviation de la paroi \theta (°)');
ylabel('Angle du choc (°)');
legend('mét. 2','mét. 1');

% Troisième  possibilite
tan_theta=f.*tan(Sigma);
theta2=Sigma-atan(tan_theta);
figure();
h1=plot(theta2*coef,Sigma*coef,'ro',theta*coef,Sigma*coef,'k-');
title('Courbe de choc (méthode 3)');
xlabel('Deviation de la paroi \theta (°)');
ylabel('Angle du choc (°)');
legend('mét. 3','mét. 1');


figure();
Mn2=Mach_aval(Mn);
M2=Mn2./sin(Sigma-theta);
h2=plot(u,v);
title(['Plan de l''hodographe: POLAIRE DE BUSEMANN, M_{amont}=' num2str(Mach)]);
xlabel('u');ylabel('v');
for i=1:50:500
    text(u(i),v(i),num2str(theta(i)*coef));
end
text(0.45,0.1,'Angle de déviation \theta indiqué en °');

figure();
plot(theta*coef,u,theta*coef,v)
title('u et v fonction de \theta, derrière le choc');
xlabel('\theta (°)');
ylabel('u,v');
legend('u','v')

figure()
plot(theta*coef,M2);
title('Mach aval en fonction de \theta');
xlabel('\theta (°)');
ylabel('Mach aval');
[val, i]=min(abs(M2-1));
theta_1=theta(i)*coef;
fprintf('Deviation pour un Mach aval de 1 = %10f ° \n',theta_1);


end

