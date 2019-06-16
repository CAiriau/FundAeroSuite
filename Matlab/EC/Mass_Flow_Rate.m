function [qm]=Mass_Flow_Rate(S,Pi,Ti,Mach)
% débit massique dans le problème de Fanno
global gam
global r
qm=S*sqrt(gam/r)*Pi/sqrt(Ti)*qm_Mach(Mach);
end


 