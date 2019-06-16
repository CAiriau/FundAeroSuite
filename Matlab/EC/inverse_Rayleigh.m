function Mach=inverse_Rayleigh(x_init,RHS)
% calcul du Mach connaissant sa fonction de Rayleigh
%
% if x_init < 1 : subsonic case
% if x_init > 1 : supersonic case

tol_max=1.e-12;
iter_max=25;
x0=x_init;
tol=1.0;
iter=0;
Mach=0.;

%calcul de la dérivée analytique
y0=Rayleigh(x0,1)-RHS;


fprintf('iter, \t tol,\t Mach \t  , err Mach, \t err F \n');

while  (iter <= iter_max) && (tol >= tol_max)
    iter=iter+1;
    dx2=-y0/Rayleigh(x0,2);
    x2=x0+dx2;
    y2=Rayleigh(x2,1)-RHS;
    errx=abs(dx2);
    erry=abs(y2);
    x0=x2; y0=y2;
    tol=erry;
    fprintf( '%3d  \t %7.4f \t %7.4f \t %7.4f \t %7.4f \t %7.4f  \n',iter,tol,x0,errx,erry)
end
 
if iter > iter_max
    disp('pas de convergence : stop')
    return
end
 

Mach=x0;
if Mach < 0.d0
    disp('Convergence vers une mauvaise valeur')
    disp('Il faut rapprocher la condition initiale (le Mach) de la solution')
    return
end

  