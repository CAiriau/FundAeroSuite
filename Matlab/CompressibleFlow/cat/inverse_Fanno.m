function Mach=inverse_Fanno(x_init,RHS)
% solve the  Mach number for a given Fanno function value
% with a Newton method
% if x_init < 1 : subsonic case
% if x_init > 1 : supersonic case

tol_max=1.e-12;
iter_max=25;
x0=x_init;
tol=1.0;
iter=0;
Mach=0.;

% Analytical derivative
y0=Fanno(x0,1)-RHS;

disp('Interted Fanno function')
fprintf('iter, \t tol,\t Mach \t  , err Mach, \t err F \n');

while  (iter <= iter_max) && (tol >= tol_max)
    iter=iter+1;
    dx2=-y0/Fanno(x0,2);
    x2=x0+dx2;
    y2=Fanno(x2,1)-RHS;
    errx=abs(dx2);
    erry=abs(y2);
    x0=x2; y0=y2;
    tol=erry;
    fprintf( '%3d  \t %7.4f \t %7.4f \t %7.4f \t %7.4f   \n',iter,tol,x0,errx,erry)
end

if iter > iter_max
    disp('not convergence : stop')
    return
end

Mach=x0;
if Mach < 0.d0
    disp('Convergence toward a bad value')
    disp('Modify the initial condition of the Mach from the estimated value ')
    return
end
disp('')