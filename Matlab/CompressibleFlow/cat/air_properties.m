function  air_properties()
% C. Airiau, avril 2012
%  Air properties

global gam
global r
global Cp
global Cv

M = 28.965338e-3;
R = 8.3144621;
r = R/M;
Cp = gam*r/(gam-1);
Cv = r/(gam-1);
end

