function [M2,omega2] = evolution_isentropique(M1,Theta2,option_char)
% Méthode des caractéristiques
% C Airiau, avril 2012

Theta1=0;
omega1=valeur_omega(M1);
fprintf('| Déviation |    = %10f °  \n',abs(Theta2));
if option_char == 0
   omega2 = omega1 + Theta2-Theta1;
else
   omega2 = omega1 + Theta1 - Theta2 ;
end
M2=valeur_Mach(omega2);
[r_T1,r_P1,r_Rho1]=rapports_isentropiques(M1);
[r_T2,r_P2,r_Rho2]=rapports_isentropiques(M2);
tau_P=r_P2/r_P1;tau_T=r_T2/r_T1;tau_Rho=r_Rho2/r_Rho1;

fprintf('omega_Amont      = %10f °  \n',omega1);
fprintf('omega_Aval       = %10f °  \n',omega2);
fprintf('Mach aval M2     = %10f \n',M2);
fprintf('PAval/PAmont     = %10f\n',tau_P);
fprintf('RhoAval/RhoAmont = %10f \n',tau_Rho);
fprintf('TAval/TAmont     = %10f \n',tau_T);

end

