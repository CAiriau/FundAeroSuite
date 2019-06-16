function y = f1(Mach,omega)

global gam;
x=Mach^2-1;
coef=180/pi;
y=coef*(sqrt((gam+1)/(gam-1))*atan(sqrt((gam-1)/(gam+1)*x))-atan(sqrt(x)))-omega;
