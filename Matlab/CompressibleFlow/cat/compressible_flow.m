function compressible_flow(choice)
% Main function to use Compressible Aerodynamic Tool
% C. Airiau, 2018 - 2020 
% 
% global variables
global gam
global r
global Cp
global Cv

close all;
gam = 1.4 ;  % gamma 

format long;
% air properties :
air_properties

option_number=21
% Menu 
display_menu;

fprintf('number of arguments : %i \n',nargin);
if nargin == 1 && choice != 0
	option = choice;
else
	option=input('enter your choice ( integer ) : ? ');
end

fprintf('Your choice is  %i  \n',option)
hline();
if option > option_number
	fprintf('bad choice, the option must be  < %i \n',option_number+1)
	return
end

switch option
    case 1
        Mach=input('Enter the upstream Mach number : ');
        solve_shock(Mach,0);
        
    case 2
        Mach=input('Enter the upstream Mach number : ');
        Theta=input('Enter the wall deviation angle in °: ( theta > 0 ) : ');
        solve_shock(Mach,Theta);
        
    case 3
        Mach=input('Enter the upstream Mach number : ');
        Theta=input('Enter the deviation angle in ° ( take care of the sign ): ');
        option_char=input('Characteristic type : [ C+ ] : 0  or [ C- ] : 1, your choice : ');
        isentropic_evolution(Mach,Theta,option_char);
        
    case 4
        Mach=input('Enter the Mach number : ');
        omega=Prandtl_Meyer(Mach);
        fprintf('omega            = %10f °\n',omega);
    case 5
        omega=input('Enter omega in ° ( omega > 0 ) : ');
        Mach=inverse_Prandtl_Meyer(omega);
        fprintf('Mach             = %10f \n',Mach)
    case 6
        Mach=input('Enter the Mach number : ');
        [TsurTi,PsurPi,RhosurRhoi]=isentropic_ratios(Mach);
        fprintf('T/Ti             = %10f\n',TsurTi);
        fprintf('P/Pi             = %10f\n',PsurPi);
        fprintf('Rho/Rhoi         = %10f\n',RhosurRhoi);
    case 7
        Mach=1:0.05:5;
        [omega]=Prandtl_Meyer(Mach);h=plot(Mach,omega); 
        xlabel('Mach');ylabel('omega in °');title('Omega function of Mach');
    case 8
        table_omega;
    case 9
        %Mach=input('Enter the Mach number : ');
        epicycloide;
    case 10
        Mach=input('Enter the Mach number : ');
        oblique_shock_curve(Mach);
    case 11
        disp('Option 1 : 2D - iso ');
        disp('       2 : 3D - iso ');
        disp('       3 : unitary M_2 ');
        option=input('Enter the option choice :'); 

        if option==3
            iso_mach(2,1);
        end
        iso_mach(2,option);     
    case 12
        Mach=input('Enter the Mach number : ');
        disp ('Point of the shock curve with  M_2 = 1')
        [theta,Sigma] = unitary_downstream_Mach_curve(Mach);
        fprintf('theta            = %10f °\n',theta);
        fprintf('sigma            = %10f °\n',Sigma);
        disp('Point of the shock curve with theta max')
        [theta,Sigma] = theta_max_curve(Mach);
        fprintf('theta            = %10f °\n',theta);
        fprintf('sigma            = %10f °\n',Sigma);
    case 13
		beta1=20; beta2=10;
		Theta1=beta1;M0=2.2;Theta2=0;Theta3=-beta2;P0=10000;T0=250;
		e_sur_l=0.003/2;
        fprintf('e / l = %8.6f \n ',e_sur_l);
		V0=M0*sqrt(gam*r*T0);
        
		fprintf('V_0 = %8.4f m/s \n',V0);
		[r_T0,r_P0,r_Rho0]=isentropic_ratios(M0);
		Pi0=P0/r_P0;Ti0=T0/r_T0;
		fprintf('P_i_0 = %8.4f Pa, \t T_i_0 = %8.4f K \n',Pi0,Ti0);

		fprintf('\n first oblique shock  \n \n')

        [M1,sigma1,tau_p1,tau_Rho,tau_T,tau_Pi] =solve_shock(M0,Theta1);
		[r_T1,r_P1,r_Rho1]=isentropic_ratios(M1);
		fprintf('M1    = %8.4f \n',M1);
		fprintf('P1/P0 = %8.4f \n',tau_p1);

		fprintf('\n first expansion \n \n')

        [M2,omega2]=isentropic_evolution(M1,Theta2-Theta1,1);
		[r_T2,r_P2,r_Rho2]=isentropic_ratios(M2);
		tau_p2=r_P2/r_P1;
		tau_p20=tau_p2*tau_p1;
		fprintf('P2/P1 = %8.4f \n',tau_p2);
		fprintf('P2/P0 = %8.4f \n',tau_p20);

		fprintf('\n second expansion \n \n')

        [M3,omega3]=isentropic_evolution(M2,Theta3-Theta2,1);
		[r_T3,r_P3,r_Rho3]=isentropic_ratios(M3);
		tau_p3=r_P3/r_P2;
		tau_p30=tau_p3*tau_p20;
		fprintf('P3/P2 = %8.4f \n',tau_p3);
		fprintf('P3/P0 = %8.4f \n',tau_p30);

		fprintf('\n second oblique shock \n \n')

        [M5,sigma5,tau_p5,tau_Rho,tau_T,tau_Pi] =solve_shock(M3,beta2);
		tau_p50=tau_p5*tau_p30;
		fprintf('P5/P3 = %8.4f \n',tau_p5);
		fprintf('P5/P0 = %8.4f \n',tau_p50);
		
		Cd=2/(gam*M0^2)*e_sur_l*(tau_p1-tau_p30);
		
        beta1=beta1*pi/180; beta2=beta2*pi/180;  
        tmp=1-tau_p1*e_sur_l*cot(beta1)-tau_p30*e_sur_l*cot(beta2);
        tmp=tmp-tau_p20*(1-e_sur_l*(cot(beta1)+cot(beta2)));
		Cl=2/(gam*M0^2)*tmp;
        f=Cl/Cd;
        fprintf('Cl = %12.10f, \t C_d = %12.10f, \t f= %7.4f \n',Cl,Cd,Cl/Cd);

    case 14
        Mach=input('Enter the Mach number : ');
        F=Rayleigh(Mach,1);
        fprintf('F_Rayleigh = %15.13f \n',F)
    
    case 15
        F=input('Enter the Rayleigh function value  : ');
        Mach=input('Enter the estimated Mach number (  subsonic M <1,  supersonic M >1) :')
        Mach=inverse_Rayleigh(Mach,F)
        fprintf('Mach  = %15.13f \n',Mach)
    case 16
        Mach=input('Enter the Mach number : ');
        Rayleigh_flow(Mach);
    case 17
        plot_Rayleigh_entropy();
    case 18
        Mach=input('Enter the Mach number : ');
        F=Fanno(Mach,1);
        fprintf('F_Fanno = %15.13f \n',F)
    case 19
        F=input('Enter the Fanno function value : ');
        Mach=input('Enter the estimated Mach number (  subsonic M <1,  supersonic M >1) :')
        Mach=inverse_Fanno(Mach,F)
        fprintf('Mach  = %15.13f \n',Mach)
    case 20
        Mach=input('Enter the Mach number : ');
        Fanno_flow(Mach);
    case 21
        Mach=input('Enter the Mach number      : ');
        Pi=input('Entrer the total pressure    Pi : ');
        Ti=input('Entrer the total temperature Ti : ');
        S=input('Enter the section in m^2 S       : ');
        fprintf('qm =  %15.13f \n',Mass_Flow_Rate(S,Pi,Ti,Mach))

    otherwise
        disp('bad choice')
end        

end
