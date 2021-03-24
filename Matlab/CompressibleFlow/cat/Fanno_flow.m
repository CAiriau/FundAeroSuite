function [q]=Fanno_flow(Ma)
% value to solve Fanno problem
% Ratios w.r.t. critical quantities 

global gam

q.Mach=Ma;

alpha=(gam+1.d0)/2.d0;
beta=gam/(gam-1.d0);

omega=1.0+(gam-1.0)/2.0*Ma^2;
tmp=0.5*(gam+1.)/omega;
q.T=tmp;
q.P=1./Ma*sqrt(tmp);
alpha=-(gam+1.)/(2.*(gam-1.));
q.rho=1./(Ma*sqrt(tmp));
q.Pi=1.0/Ma*tmp^alpha;

fprintf('T/T*                               :  %15.13f\n', q.T)
fprintf('P/P*                               :  %15.13f\n', q.P)
fprintf('rho/rho*                           :  %15.13f\n', q.rho)
fprintf('Pi/Pi*                             :  %15.13f\n', q.Pi)
end  