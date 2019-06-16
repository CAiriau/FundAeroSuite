function [s,Cp]=panel_method(option,alpha)
% Main program, choose the option
% alpha is optionnaly the angle of attack in degrees.
%
%%*******************************************************
% Original coding obtained from J. Katz, in fortran 77
% and found in :
% Low-Speed Aerodynamics,Joseph Katz et Allen Plotkin
% Second Edition,Cambridge Aerospace Series, 2001
% 
% translation in Matlab by Christophe Airiau, 2015
% free use of the matlab code, under J. Katz agreement.
%%*******************************************************
% Use : 3 cases
%  main()       
%  main(n)              % n is the option (< 16)
%  main(n,alpha)        % alpha is the AoA in degrees

close all
global name

line='===========================================================';
fprintf('%s \n PANEL METHOD for airfoil  and wing: various distribution,  \n %s\n',line,line);
fprintf('Adapted in matlab by Christophe Airiau\n');
fprintf('From J. Katz routines in f77 \n %s\n',line);

switch nargin
case 0 
    option=menu();
    alpha=0.;  % angle of Attack in °
case 1
    alpha=0.;  % angle of Attack in °
case 2
    sprintf(' option : %i \n',option);
end

sprintf(' angle of attack : %5.2f °\n ',alpha);

switch option
case 1
    [s,Y,Cp_ref]=grid_generator(true);
case 2  
    [s,Y,Cp_ref]=grid_generator(true);
    [sc,Cp,Cl]=Constant_Strength_Source(s,Y);
    name='Constant Strength Source';
case 3  
    [s,Y,Cp_ref]=grid_generator(true);
    [sc,Cp,Cl]=Constant_Strength_Doublet(alpha,s,Y);
    name='Constant Strength Doublet';
case 4  
    [s,Y,Cp_ref]=grid_generator(true);
    [sc,Cp,Cl]=Constant_Strength_Vortex(alpha,s,Y);
    name='Constant Strength Vortex';
case 5  
    alpha=0; x_wake=1.2;
    [s,Y,Cp_ref]=grid_generator(true);
    [sc,Cp,Cl]=Linear_Strength_Source(alpha,s,Y,x_wake);
    name='Linear Strength Source';
case 6  
    [s,Y,Cp_ref]=grid_generator(true);
    [sc,Cp,Cl]=Linear_Strength_Vortex(alpha,s,Y);
    name='Linear Strength Vortex';
case 7  
    [s,Y,Cp_ref]=grid_generator(true);
    [sc,Cp,Cl]=Constant_Strength_Doublet_Potential(alpha,s,Y);
    name='Constant Strength Doublet Potential';
case 8  
    [s,Y,Cp_ref]=grid_generator(true);
    [sc,Cp,Cl]=Constant_Strength_Source_Doublet_Potential(alpha,s,Y);
    name='Constant Strength Source Doublet Potential';
case 9  
    [s,Y,Cp_ref]=grid_generator(true);
    [sc,Cp,Cl]=Linear_Strength_Doublet_Potential(alpha,s,Y);
    name='Linear Strength Doublet Potential';
case 10  
    x_internal=0.98758;
    [s,Y,Cp_ref]=grid_generator(true);
    [sc,Cp,Cl]=Quadratic_Strength_Doublet_Potential(alpha,s,Y,x_internal);
    name='Quadratic Strength Doublet Potential';
case 11  
    Influence_Coefficients();
    name='Influence Coefficient';
case 12  
    Rectangular_Lifting_Surface();
    name='Rectangular Lifting Surface';
case 13  
    [time,Cl,Cd,Clt,gam1]=lumped_vortex();
    name='Lumped Vortex';
case 14
    Unsteady_Rectangular_Lifting_Surface();
    name='Unsteady Rectangular Lifting Surface';
case 15
    Discrete_Vortex_Method(30);
    name='Discrete Vortex Method';
case 16
    [s,Y,Cp_ref]=grid_generator([0.1,0,10,181,1,2]);
    [sc,Cp,Cl]=Constant_Strength_Vortex(alpha,s,Y);
    figure();
    plot(s,Cp_ref,'k-',sc,Cp,'b--');
    title([ 'Cp : ' name]);
    xlim([0 1]);
    ylim([-3  1]);
    legend('Ref','Current','Location','southeast');

end
if option > 1 && option < 11
    figure();
    plot(s,Cp_ref,'k-',sc,Cp,'b--');
    title([ 'Cp : ' name]);
    legend('Ref','Current','Location','southeast');
end
disp('Normal end of execution');
end
