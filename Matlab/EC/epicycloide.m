function epicycloide
% tracé de l'épicycloide
% C. Airiau, avril 2012
global gam

Mstar=sqrt((gam+1)/(gam-1));
fprintf('M*               = %10f\n',Mstar);
Mach=1:0.1:50;
[omega]=omega_star(Mach,Mstar);
Ms=sqrt((gam+1).*Mach.^2./(2+(gam-1).*Mach.^2));
us=Ms.*cos(omega); vs=Ms.*sin(omega);
ur=Mstar*cos(omega); vr=Mstar*sin(omega);
figure()
h=plot(us,vs,'b-',ur,vr,'k-','LineWidth',2);
xlabel('u/a');ylabel('v/a');title('Epicycloide');
title('Plan de l''hodographe, epicycloide, Ecoulement isentropique')
hleg=legend('(u_a,v_a) = f(Mach)','limite');
set(hleg,'Location','SouthWest');
for i=1:5:50
  text(us(i),vs(i),['M= ' num2str(Mach(i))]);
end


end

function [omega]=omega_star(Mach,Mstar)
x=sqrt(Mach.^2-1);
omega=Mstar*atan(x/Mstar)-atan(x);
end
