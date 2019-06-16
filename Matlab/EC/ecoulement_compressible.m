function ecoulement_compressible(choix)
% Aérodynamique compressible
% C. Airiau, 2018 
% 
% Variables globales :
global gam
global r
global Cp
global Cv

close all;
gam = 1.4 ;  % gamma 
 
format long;
% proriétés de l'air :
propriete

nombre_option=21
% Menu
disp(' *********************************************')
disp('          ECOULEMENTS COMPRESSIBLES')
disp('           (C. Airiau, 2012 - 2015)')
disp(' *********************************************')
disp(' ')
disp('***** Menu du programme *****')
disp(' ')
disp(' 1 - Ondes de choc droites')
disp(' 2 - Ondes de choc obliques')
disp(' 3 - Compression ou détente isentropique')
disp(' 4 - Omega fonction du Mach (une valeur)')
disp(' 5 - Mach fonction de Omega (une valeur)')
disp(' 6 - Rapport des grandeurs isentropiques pour une valeur de Mach')    
disp(' 7 - Figure Omega fonction du Mach')    
disp(' 8 - Table Omega fonction du Mach')
disp(' 9 - Epicycloïde')    
disp('10 - Courbes de choc et Polaire de choc de Busemann')
disp('11 - Courbe de choc multiples en 2D ou 3D, courbe M2=1 et theta max ')
disp('12 - Choc : Mach aval unitaire et theta max pour une valeur du Mach amont')
disp('13 - Exercice: choc, 2 détentes + choc')
disp('14 - Fonction de Rayleigh')
disp('15 - Fonction inverse de Rayleigh')
disp('16 - Problème de Rayleigh')
disp('17 - Courbes de Rayleigh')
disp('18 - Fonction de Fanno')
disp('19 - Fonction inverse de Fanno')
disp('20 - Problème Fanno')
disp('21 - débit Fanno')





fprintf('number of arguments : %i \n',nargin);
if nargin == 1
	option = choix;
else
	option=input('entrer un entier en fonction de votre choix : ? ');
end

fprintf('votre choix est %i  \n',option)
hline();
if option > nombre_option
	disp('mauvais choix, il faut un nombre < 14')
	return
end

switch option
    case 1
        Mach=input('Entrer le nombre de Mach amont : ');
        chocs(Mach,0);
        
    case 2
        Mach=input('Entrer le nombre de Mach amont : ');
        Theta=input('Entrer l''angle de la paroi en °: (theta >0) : ');
        chocs(Mach,Theta);
         
    case 3
        Mach=input('Entrer le nombre de Mach amont : ');
        Theta=input('Entrer l''angle de deviation en degrés (attention au signe): ');
        option_char=input('type de caractéristiques : [ C+ ] : 0  ou [ C- ] : 1, choix : ');
        evolution_isentropique(Mach,Theta,option_char);
        
	case 4
         Mach=input('Entrer le nombre de Mach : ');
         omega=valeur_omega(Mach);
         fprintf('Omega            = %10f °\n',omega);
    case 5
        omega=input('Entrer omega en degrés (omega > 0) : ');
        Mach=valeur_Mach(omega);
        fprintf('Mach             = %10f \n',Mach)
    case 6
        Mach=input('Entrer le nombre de Mach : ');
        [TsurTi,PsurPi,RhosurRhoi]=rapports_isentropiques(Mach);
        fprintf('T/Ti             = %10f\n',TsurTi);
        fprintf('P/Pi             = %10f\n',PsurPi);
        fprintf('Rho/Rhoi         = %10f\n',RhosurRhoi);
    case 7
        Mach=1:0.05:5;
        [omega]=valeur_omega(Mach);h=plot(Mach,omega); 
        xlabel('Mach');ylabel('omega en degrés');title('Omega fonction du Mach');
    case 8
        table_omega;
    case 9
        %Mach=input('Entrer le nombre de Mach : ');
        epicycloide;
    case 10
        Mach=input('Entrer le nombre de Mach : ');
        courbe_choc(Mach);
    case 11
        disp('Option 1 : iso en 2D');
        disp('       2 : iso en 3D');
        disp('       3 : M2 unitaire ');
        option=input('Entrer l''option '); 

        if option==3
          iso_mach(2,1);
        end
        iso_mach(2,option);     
    case 12
        Mach=input('Entrer le nombre de Mach : ');
         disp ('Point de la courbe de choc avec M2 = 1')
        [theta,Sigma] = courbe_Mach_aval_unitaire(Mach);
        fprintf('Theta            = %10f °\n',theta);
        fprintf('Sigma            = %10f °\n',Sigma);
        disp('Point de la courbe de choc avec theta max')
        [theta,Sigma] = courbe_theta_max(Mach);
        fprintf('Theta            = %10f °\n',theta);
        fprintf('Sigma            = %10f °\n',Sigma);
    case 13
		beta1=20; beta2=10;
		Theta1=beta1;M0=2.2;Theta2=0;Theta3=-beta2;P0=10000;T0=250;
		e_sur_l=0.003/2;
        fprintf('e/l = %8.6f \n ',e_sur_l);
		V0=M0*sqrt(gam*r*T0);
        
		fprintf('V_0 = %8.4f m/s \n',V0);
		[r_T0,r_P0,r_Rho0]=rapports_isentropiques(M0);
		Pi0=P0/r_P0;Ti0=T0/r_T0;
		fprintf('P_i_0 = %8.4f Pa, \t T_i_0 = %8.4f K \n',Pi0,Ti0);

		fprintf('\n premier choc  \n \n')

        [M1,sigma1,tau_p1,tau_Rho,tau_T,tau_Pi] =chocs(M0,Theta1);
		[r_T1,r_P1,r_Rho1]=rapports_isentropiques(M1);
		fprintf('M1    = %8.4f \n',M1);
		fprintf('P1/P0 = %8.4f \n',tau_p1);

		fprintf('\n première détente \n \n')

        [M2,omega2]=evolution_isentropique(M1,Theta2-Theta1,1);
		[r_T2,r_P2,r_Rho2]=rapports_isentropiques(M2);
		tau_p2=r_P2/r_P1;
		tau_p20=tau_p2*tau_p1;
		fprintf('P2/P1 = %8.4f \n',tau_p2);
		fprintf('P2/P0 = %8.4f \n',tau_p20);

		fprintf('\n seconde détente \n \n')

        [M3,omega3]=evolution_isentropique(M2,Theta3-Theta2,1);
		[r_T3,r_P3,r_Rho3]=rapports_isentropiques(M3);
		tau_p3=r_P3/r_P2;
		tau_p30=tau_p3*tau_p20;
		fprintf('P3/P2 = %8.4f \n',tau_p3);
		fprintf('P3/P0 = %8.4f \n',tau_p30);

		fprintf('\n second choc \n \n')

        [M5,sigma5,tau_p5,tau_Rho,tau_T,tau_Pi] =chocs(M3,beta2);
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
        Mach=input('Entrer le nombre de Mach : ');
        F=Rayleigh(Mach,1);
        fprintf('F_Rayleigh = %15.13f \n',F)
    
    case 15
        F=input('Entrer la valeur de la fonction de Rayleigh : ');
        Mach=input('Entrer le Mach estimé (  subsonique M <1,  supersonique M >1) :')
        Mach=inverse_Rayleigh(Mach,F)
        fprintf('Mach  = %15.13f \n',Mach)
    case 16
        Mach=input('Entrer le nombre de Mach : ');
        Rayleigh_flow(Mach);
    case 17
        plot_Rayleigh_entropie();
    case 18
        Mach=input('Entrer le nombre de Mach : ');
        F=Fanno(Mach,1);
        fprintf('F_Fanno = %15.13f \n',F)
    case 19
        F=input('Entrer la valeur de la fonction de Fanno : ');
        Mach=input('Entrer le Mach estimé (  subsonique M <1,  supersonique M >1) :')
        Mach=inverse_Fanno(Mach,F)
        fprintf('Mach  = %15.13f \n',Mach)
    case 20
        Mach=input('Entrer le nombre de Mach : ');
        Fanno_flow(Mach);
    case 21
        Mach=input('Entrer le nombre de Mach      : ');
        Pi=input('Entrer la pression totale Pi    : ');
        Ti=input('Entrer la température  arrêt Ti : ');
        S=input('Entrer la section en m^2 S       : ');
        fprintf('qm =  %15.13f \n',Mass_Flow_Rate(S,Pi,Ti,Mach))

    otherwise
        disp('mauvais choix')
end        

end
