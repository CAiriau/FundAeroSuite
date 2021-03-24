function Exercise_Fanno()
% exercise done in M2 MSME.
% C. AIRIAU, Octobre 2017

% global variable:
global gam
global r
global Cp
global Cv

close all;
gam = 1.4 ;  % gamma 

format long;

air_properties;

%**********************
% SUBSONIC CASE
%**********************
% Question 1 
disp('Question 1, subsonic')
D=0.1;L=20;Cf=0.005;
M1=0.24;
F1=Fanno(M1,1);
L_critique=F1*D/(4*Cf);
printf('critical length       : %f\n',L_critique)
F2=F1-4/D*Cf*L;
M2=inverse_Fanno(0.4,F2);
disp('')
printf('M2                      : %f\n',M2)
q1=Fanno_flow(M1);
q2=Fanno_flow(M2);
disp('')
printf('P2/P1                    : %f \n',q2.P/q1.P)
printf('Pi2/Pi1-1 (in percent)  : %f\n',(q2.Pi/q1.Pi-1)*100)

disp('')
disp('Question 2, subsonic')
% conduite rectangulaire
h=0.02;w=0.01;S=h*w;c=2*(h+w);D=4*S/c;
printf('Dh                      : %f\n',D)
M1=0.5;
F1=Fanno(M1,1);
L_critique=F1*D/(4*Cf);
printf('critical length       : %f\n',L_critique)
L=2.16*L_critique;
F1=L*4*Cf/D;
M1p=inverse_Fanno(0.4,F1);
printf('M1 new               : %f\n',M1p)
q1=Fanno_flow(M1p);
disp('')
printf('P2/P1                    : %f \n',1/q1.P)
printf('relative variation of the mass flow rate : %f in percent\n',100*(qm_Mach(M1p)/qm_Mach(M1)-1))


%**********************
% SUPERSONIC CASE
%**********************
% Question 1 
disp('Question 1, supersonic')
L=40;D=1;Cf=0.005;M1=2.6;
F1=Fanno(M1,1);
L_critique=F1*D/(4*Cf);
printf('critical length     : %f\n',L_critique)
disp('Question 2, supersonic')

[test]=calculate_length(M1,1.5,L,D,Cf);
[test]=calculate_length(M1,2.5,L,D,Cf);
[test]=calculate_length(M1,2.128,L,D,Cf);

end

function [test]=calculate_length(M1,M1p,L,D,Cf)
  F1=Fanno(M1,1);
  L1=D/(4*Cf)*(F1-Fanno(M1p,1));
  printf('supersonic length   : %f\n',L1)
  M2=downstream_normal_Mach(M1p);
  printf('Downstream Mach of the normal shock       : %f\n',M2)
  L2=Fanno(M2,1)*D/(4*Cf);
  printf('subsonic length     : %f\n',L2)
  printf('total length        : %f\n',L1+L2)
  test=false
  if abs(L1+L2-L) <= 1e-2
    disp('You have the solution')
    test=true;
  elseif L1+L2 < L 
    disp('The lenght is lower than the target')
    disp('Increase M1p')
  else  L1+L2 > L 
    disp('The lenght is greater than the target')
    disp('Decrease M1p')
  end
  
end