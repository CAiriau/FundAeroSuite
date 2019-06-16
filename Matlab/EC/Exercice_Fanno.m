function Exercice_Fanno()
% exercice fait dans le M2 MSME.
% C. AIRIAU, Octobre 2017

% Variables globales :
global gam
global r
global Cp
global Cv

close all;
gam = 1.4 ;  % gamma 
 
format long;
% proriétés de l'air :
propriete;

%**********************
% CAS SUBSONIQUE
%**********************
% Question 1 
disp('Question 1, subsonique')
D=0.1;L=20;Cf=0.005;
M1=0.24;
F1=Fanno(M1,1);
L_critique=F1*D/(4*Cf);
printf('Longueur critique       : %f\n',L_critique)
F2=F1-4/D*Cf*L;
M2=inverse_Fanno(0.4,F2);
disp('')
printf('M2                      : %f\n',M2)
q1=Fanno_flow(M1);
q2=Fanno_flow(M2);
disp('')
printf('P2/P1                    : %f \n',q2.P/q1.P)
printf('Pi2/Pi1-1 (en pourcent)  : %f\n',(q2.Pi/q1.Pi-1)*100)

disp('')
disp('Question 2, subsonique')
% conduite rectangulaire
h=0.02;w=0.01;S=h*w;c=2*(h+w);D=4*S/c;
printf('Dh                      : %f\n',D)
M1=0.5;
F1=Fanno(M1,1);
L_critique=F1*D/(4*Cf);
printf('Longueur critique       : %f\n',L_critique)
L=2.16*L_critique;
F1=L*4*Cf/D;
M1p=inverse_Fanno(0.4,F1);
printf('M1 nouveau               : %f\n',M1p)
q1=Fanno_flow(M1p);
disp('')
printf('P2/P1                    : %f \n',1/q1.P)
printf('variation de d�bit relative: %f en pourcent\n',100*(qm_Mach(M1p)/qm_Mach(M1)-1))


%**********************
% CAS SUPERSONIQUE
%**********************
% Question 1 
disp('Question 1, supersonique')
L=40;D=1;Cf=0.005;M1=2.6;
F1=Fanno(M1,1);
L_critique=F1*D/(4*Cf);
printf('Longueur critique       : %f\n',L_critique)
disp('Question 2, supersonique')

[test]=calcul_longueur(M1,1.5,L,D,Cf);
[test]=calcul_longueur(M1,2.5,L,D,Cf);
[test]=calcul_longueur(M1,2.128,L,D,Cf);
 
end


function [test]=calcul_longueur(M1,M1p,L,D,Cf)
  F1=Fanno(M1,1);
  L1=D/(4*Cf)*(F1-Fanno(M1p,1));
  printf('Longueur supersonique   : %f\n',L1)
  M2=Mach_aval(M1p);
  printf('Mach aval au choc       : %f\n',M2)
  L2=Fanno(M2,1)*D/(4*Cf);
  printf('Longueur subsonique     : %f\n',L2)
  printf('Longueur totale         : %f\n',L1+L2)
  test=false
  if abs(L1+L2-L) <= 1e-2
    disp('Vous avez trouve la solution')
    test=true;
  elseif L1+L2 < L 
    disp('La longueur est inferieure a celle visee')
    disp('Augmenter M1p')
  else  L1+L2 > L 
    disp('La longueur est superieure a celle visee')
    disp('Diminuer M1p')
  end
  
end