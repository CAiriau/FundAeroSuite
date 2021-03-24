function [TsurTi,PsurPi,RhosurRhoi]=isentropic_ratios(Mach)
% isentropic ratio T/T_i, P/P_i, rho/rho_i

global gam

TisurT=1+(gam-1)/2*Mach^2;
TsurTi=1/TisurT;
x=-gam/(gam-1);
y=-1/(gam-1);
PsurPi=(TisurT)^x;
RhosurRhoi=(TisurT)^y;
end