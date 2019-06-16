function [Mach] = valeur_Mach(omega,M_init,method)
% connaissant la valeur de omega, trouver le nombre de Mach correspondant
% C. Airiau, avril 2012
% Seul le premier argument est obligatoire.
switch nargin

case 1
    M_init=2
    method=1
case 2
    fprintf('modification de la valeur initiale du Mach pour Omega \n')
    method=1
case 3
    fprintf('Méthode de Newton\n')
end


if method == 1
    % Pour afficher les détails décommenter la ligne ci dessous et commenter la ligne suivante
    %Mach = fzero(@(Mach)f1(Mach,omega), M_init, optimset('TolX', 1e-4, 'TolFun', 1e-4,'Display','Iter'));
    %Mach = fzero(@(Mach)f1(Mach,omega), M_init, optimset('TolX', 1e-4, 'TolFun', 1e-4));
    Mach = fsolve(@(Mach)f1(Mach,omega), M_init, optimset('TolX', 1e-4, 'TolFun', 1e-4));
else
    % En cas de difficulté de l'algorithme précédent, prendre l'algorithme de Newton...
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
% fonction omega en degrés
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


