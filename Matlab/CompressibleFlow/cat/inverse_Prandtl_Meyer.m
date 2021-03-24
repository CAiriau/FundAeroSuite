function [Mach] = inverse_Prandtl_Meyer(omega,M_init,method)
% For a given value of the Prandtl-Meyer function, find the Mach number with 
% a Newton method  
% C. Airiau, avril 2012
% only the firt argument is mandatory
switch nargin

case 1
    M_init=2
    method=1
case 2
    fprintf('change of the initial value of the  Mach number for omega \n')
    method=1
case 3
    fprintf('Newton method\n')
end


if method == 1
    % to display details uncomment the line below and  comment the uncommented line
    %Mach = fzero(@(Mach)f1(Mach,omega), M_init, optimset('TolX', 1e-4, 'TolFun', 1e-4,'Display','Iter'));
    %Mach = fzero(@(Mach)f1(Mach,omega), M_init, optimset('TolX', 1e-4, 'TolFun', 1e-4));
    Mach = fsolve(@(Mach)f1(Mach,omega), M_init, optimset('TolX', 1e-4, 'TolFun', 1e-4));
else
    % In case of problem of the previous algorithm choose the Newton algorithm below
    iterMax=50;erreur_Max=1e-6;
    iter=0;
    erreur=1;
    Mach=M_init;
    
    while erreur > erreur_Max && iter < iterMax
        dMach=-f1(Mach,omega)/df1(Mach);
        erreur=abs(dMach);
        Mach=Mach+dMach;
        %fprintf('i= %d \t e= %f \t Mach= %f \n',iter,erreur,Mach)
        iter=iter+1;
    end 
    fprintf('i= %d \t e= %f \t Mach= %f \n',iter-1,erreur,Mach)

end
fprintf('Mach = %2.9e \n',Mach);
end

function y = f1(Mach,omega)
% Prandtl-Meyer function in degrees
global gam;
x=Mach^2-1;
coef=180/pi;
y=coef*(sqrt((gam+1)/(gam-1))*atan(sqrt((gam-1)/(gam+1)*x))-atan(sqrt(x)))-omega;
end

function y=df1(Mach)
% d omega/ d Mach en degrés
global gam;
coef=180/pi;
y=2*sqrt(Mach^2-1)/(2+(gam-1)*Mach^2)*coef;
end


