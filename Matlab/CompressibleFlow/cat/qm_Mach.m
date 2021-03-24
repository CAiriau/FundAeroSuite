function [res]=qm_Mach(Mach) 
% part of the Mach dependence in the mass flow rate formula
global gam
res=Mach*(1.0+(gam-1.0)/2.0*Mach^2)^(-(gam+1.)/2./(gam-1.));
end 
