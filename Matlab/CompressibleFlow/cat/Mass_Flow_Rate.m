function [qm]=Mass_Flow_Rate(S,Pi,Ti,Mach)
% Mass flow rate in Fanno's problem
global gam
global r
qm=S*sqrt(gam/r)*Pi/sqrt(Ti)*qm_Mach(Mach);
end