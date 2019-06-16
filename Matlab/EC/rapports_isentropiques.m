function [TsurTi,PsurPi,RhosurRhoi]=rapports_isentropiques(Mach);
% rapport des quantités par rapport aux valeurs isentropiques
%  C. Airiau
global gam

TisurT=1+(gam-1)/2*Mach^2;
TsurTi=1/TisurT;
x=-gam/(gam-1);
y=-1/(gam-1);
PsurPi=(TisurT)^x;
RhosurRhoi=(TisurT)^y;
end