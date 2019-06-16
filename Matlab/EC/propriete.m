function  propriete
%function [r,Cp,Cv] = propriete(gam)
% C. Airiau, avril 2012
%  propriétés de l'air
 
global gam
global r
global Cp
global Cv

M=28.965338e-3;
R = 8.3144621;
r=R/M;
Cp=gam*r/(gam-1);
Cv=r/(gam-1);
end

