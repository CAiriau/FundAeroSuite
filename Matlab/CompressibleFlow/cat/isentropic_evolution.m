function [M2,omega2] = isentropic_evolution(M1,Theta2,option_char)
% Characteristic method
% C Airiau, avril 2012

Theta1=0;
omega1=Prandtl_Meyer(M1);
fprintf('| Deviation |    = %10f °  \n',abs(Theta2));
if option_char == 0
   omega2 = omega1 + Theta2-Theta1;
else
   omega2 = omega1 + Theta1 - Theta2 ;
end
M2=inverse_Prandtl_Meyer(omega2);
[r_T1,r_P1,r_Rho1]=isentropic_ratios(M1);
[r_T2,r_P2,r_Rho2]=isentropic_ratios(M2);
tau_P=r_P2/r_P1;tau_T=r_T2/r_T1;tau_Rho=r_Rho2/r_Rho1;

fprintf('Upstream omega    = %10f °  \n',omega1);
fprintf('Downstream omega  = %10f °  \n',omega2);
fprintf('Downstream Mach   = %10f \n',M2);
fprintf('P_down / P_up     = %10f\n',tau_P);
fprintf('Rho_down / Rho_up = %10f \n',tau_Rho);
fprintf('T_down / T_up     = %10f \n',tau_T);

end

