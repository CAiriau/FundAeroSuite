function [q]=Rayleigh_flow(Ma)
% valeurs pour résoudre l'écoulement de Rayleigh
% Rapports quantité/quantité critique
global gam
q.Mach=Ma;

alpha=(gam+1.d0)/2.d0;
beta=gam/(gam-1.d0);
omega=1.0+(gam-1.0)/2.0*Ma^2;
q.P=(gam+1.0)/(1.0+gam*Ma^2);
q.T=Ma^2*q.P^2;
q.rho=1.0/(Ma^2*q.P);
q.u=1.0/q.rho;
q.Ti=q.T*omega/alpha;
q.Pi=q.P*(omega/alpha)^beta;
q.s=beta*log(q.T)-log(q.P);
fprintf('T/T*                               :  %15.13f\n', q.T)
fprintf('P/P*                               :  %15.13f\n', q.P)
fprintf('rho/rho*                           :  %15.13f\n', q.rho)
fprintf('Pi/Pi*                             :  %15.13f\n', q.Pi)
fprintf('Ti/Ti*                             :  %15.13f\n', q.Ti)
fprintf('(s-s*)/r                           :  %15.13f\n', q.s)
fprintf('(s-s*)/Cp                          :  %15.13f\n', q.s/beta)
end  