!    This source code is a part of the FundAeroSuite
!    Copyright (C) <2018>  < Christophe Airiau >
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Affero General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!    @author : Christophe.airiau@imft.fr
!

 
!> MODULE :   Correction des exercices du chapitre 10
!! 
!! @author C. Airiau
!!
module exercices_chap10

use mod_EC
use Constante

contains

!******************************************************************
subroutine Exercice_10_1
!******************************************************************
!>@details Calcul d'un choc droit
implicit none
real(kind=8)        :: M0,rho0,P0,Ti0,Pi0,rhoi0
real(kind=8)        :: rT0,rP0,rr0,U0,a0,T0
real(kind=8)        :: M1,rho1,P1,Pi1,U1,a1,T1
real(kind=8)        :: rT1,rP1,rr1,rPi1,delta_s

if (fichier) open(6,form='formatted',file='Solution_10-1.out')
write(6,110) 'Exercice 10.1 : Choc droit'
M0=3.d0;rho0=1.25;P0=1.0d5

rT0=T_Ti(M0);
T0=P0/(r*rho0);Ti0=T0/rT0;
rP0=P_Pi(M0);Pi0=P0/rP0;
rr0=rho_rhoi(M0);rhoi0=rho0/rr0;
a0=sqrt(gam*r*T0); U0=M0*a0

write(6,*) 'T0/Ti0                  = ',rT0
write(6,*) 'P0/Pi0                  = ',rP0
write(6,*) 'T0                      = ',T0
write(6,*) 'rho0                    = ',rho0
write(6,*) 'Ti0                     = ',Ti0
write(6,*) 'Pi0                     = ',Pi0
write(6,*) 'rhoi0                   = ',rhoi0
write(6,*) 'a0                      = ',a0
write(6,*) 'U0                      = ',U0

write(6,100) 'Grandeur en aval du choc droit'

M1=Mach_Aval(M0);
rP1=P2_P1(M0);P1=P0*rP1;
rr1=rho2_rho1(M0);rho1=rho0*rr1;
rPi1=pi2_Pi1(M0);Pi1=Pi0*rPi1;
write(6,*) 'Pi1/Pi0                 = ',rPi1
write(6,*) 'P1/P0                   = ',rP1
write(6,*) 'rho1/rho0               = ',rr1
rT1=T_Ti(M1);T1=rT1*Ti0;
a1=sqrt(gam*r*T1); U1=M1*a1
delta_s=Saut_Entropie(Pi0,Pi1)
write(6,*) 'M1                      = ',M1
write(6,*) 'T1                      = ',T1
write(6,*) 'P1                      = ',P1
write(6,*) 'rho1                    = ',rho1
write(6,*) 'Pi1                     = ',Pi1
write(6,*) 'a1                      = ',a1
write(6,*) 'U1                      = ',U1
write(6,*) 's1-s0                   = ',delta_s


100  format(/,50('#'),/,3x,a,/,50('#'),/)
110  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_10_1

!******************************************************************
subroutine Exercice_10_2
!******************************************************************
!>@details Calcul pratique d'un choc droit

implicit none
real(kind=8)        :: M0,P0,Ti0,Pi0
real(kind=8)        :: rT0,rP0,U0,a0,T0
real(kind=8)        :: M1,P1,Pi1,U1,a1,T1
real(kind=8)        :: rT1,rP1,rPi1,delta_s,rPi0
real(kind=8)        :: Cp

if (fichier) open(6,form='formatted',file='Solution_10-2.out')
write(6,110) 'Exercice 10.2 : Calcul pratique d un choc droit'
Cp=gam*r/(gam-1.d0)
U0=800.d0;P0=1.0d5;T0=300.d0;

a0=sqrt(gam*r*T0); M0=U0/a0

rT0=T_Ti(M0);Ti0=T0/rT0;
rP0=P_Pi(M0);Pi0=P0/rP0;

write(6,*) 'a0                      = ',a0
write(6,*) 'M0                      = ',M0
write(6,*) 'Ti0                     = ',Ti0
write(6,*) 'Pi0                     = ',Pi0
write(6,*) 'T0/Ti0                  = ',rT0
write(6,*) 'Pi0/Po                  = ',1.d0/rP0


write(6,100) 'Grandeur en aval du choc droit'

M1=Mach_Aval(M0);
rP1=P2_P1(M0);P1=P0*rP1;
rPi1=pi2_Pi1(M0);Pi1=Pi0*rPi1;
rT1=T_Ti(M1);T1=rT1*Ti0;
a1=sqrt(gam*r*T1); U1=M1*a1
delta_s=Saut_Entropie(Pi0,Pi1)

write(6,*) 'Pi1/Pi0                 = ',rPi1
write(6,*) 'T1/Ti1                  = ',rT1
write(6,*) 'P1/P0                   = ',rP1
write(6,*) 'M1                      = ',M1
write(6,*) 'T1                      = ',T1
write(6,*) 'P1                      = ',P1
write(6,*) 'Pi1                     = ',Pi1
write(6,*) 'a1                      = ',a1
write(6,*) 'U1                      = ',U1
write(6,*) 's1-s0                   = ',delta_s


write(6,100) 'Ecoulement isentropique de U0 à U1'
! Cp T1 + U1**2/2 = Cp T0+U0**2/2 = a**2/(gam-1)+U**2/2
! calcul de a1 
! meme Ti, meme vitesse donc meme a1 et meme Mach et meme température
! seule la pression change mais Pi est conservée
a1=sqrt(a0**2+(U0**2-U1**2)*(gam-1)/2)
M1=U1/a1; T1=a1**2/(gam*r);
rPi1=P_Pi(M1);
write(6,*) 'P1/Pi0                  = ',rPi1
rPi0=P_Pi(M0);
write(6,*) 'P0/Pi0                  = ',rPi0
P1=P0*rPi1/rPi0;

write(6,*) 'P1/P0                   = ',P1/P0
write(6,*) 'a1                      = ',a1
write(6,*) 'M1                      = ',M1
write(6,*) 'T1                      = ',T1
write(6,*) 'P1                      = ',P1




100  format(/,50('#'),/,3x,a,/,50('#'),/)
110  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_10_2




!******************************************************************
subroutine Exo_Laval
!******************************************************************
!>@details  un exercice sur la tuyère de Laval
implicit none

real(kind=8)        :: Rinlet=0.4d0,Rexit=0.12d0,Rcol=0.09d0
real(kind=8)        :: Ti0=350.d0,Pi0,Pa=101325.d0
real(kind=8)        :: Me,omega,rapport_section,Rcrit,Mcol,Scol,Pexit
real(kind=8)        :: Minlet,Tinlet,ainlet,Uinlet,rhoinlet,Pinlet
real(kind=8)        :: Qm,Pi2
real(kind=8)        :: Rchoc,Mchoc,Maval

if (fichier) open(6,form='formatted',file='Solution_Laval.out')
write(6,110) 'Exercice Laval'

write(6,100) 'Question 2-1'
Pi0=103000d0
write(6,*) 'Me                      = ',inverse_P_Pi(Pa/Pi0)


write(6,100) 'Question 2-2'
Me=0.3d0
omega=P_Pi(Me); 
write(6,*) 'P/Pi                    = ',omega
Pi0=Pa/omega;   
write(6,*) 'Pi_0                    = ',Pi0
rapport_section= A_Acritic(Me)
Rcrit= Rexit/sqrt(rapport_section)
write(6,*) 'A/Ac                    = ',rapport_section
write(6,*) 'R_crit                  = ',Rcrit

write(6,*) 'Col :'
rapport_section=(Rcol/Rcrit)**2
Mcol=inverse_A_Acritic(0.5d0,rapport_section,1,1,45,1.d-12) 
write(6,*) 'Acol/Ac                 = ',rapport_section
write(6,*) ' Mcol                   = ',Mcol


write(6,*) 'Entrée :'

rapport_section=(Rinlet/Rcrit)**2
Minlet=inverse_A_Acritic(0.5d0,rapport_section,1,1,45,1.d-12) 
write(6,*) 'Ainlet/Ac               = ',rapport_section
write(6,*) 'Minlet                  = ',Minlet

omega=T_Ti(Minlet)
write(6,*) 'T/Ti                    = ',omega
Tinlet=Ti0*omega; 
write(6,*) 'Tinlet                  = ',Tinlet
ainlet=sqrt(gam*r*Tinlet);
write(6,*) 'a_inlet                 = ',ainlet
Uinlet=Minlet*ainlet;
write(6,*) 'Uinlet                  = ',Uinlet

omega=P_Pi(Minlet)
write(6,*) 'P/Pi                    = ',omega
Pinlet=Pi0*omega; 
write(6,*) 'Pinlet                  = ',Pinlet
rhoinlet=Pinlet/(r*Tinlet);
write(6,*) 'rho_inlet               = ',rhoinlet
Qm=rhoinlet*Uinlet*val_pi*Rinlet**2
write(6,*) 'Débit massique Qm       = ',Qm

write(6,100) 'Question 2-3'
Rcrit=Rcol
rapport_section=(Rexit/Rcrit)**2
Me=inverse_A_Acritic(0.5d0,rapport_section,1,1,45,1.d-12) 
omega=P_Pi(Me); 
write(6,*) 'P/Pi                    =',omega
Pi0=Pa/omega;   
write(6,*) 'Pi_O                    =',Pi0

write(6,100) 'Question 3-1'

Me=2.d0
omega=P_Pi(Me); 
write(6,*) 'P/Pi                    =',omega
Pi0=Pa/omega;   
write(6,*) 'Pi_O                    =',Pi0

write(6,100) 'Question 3-2'
Scol=val_pi*Rcol**2
write(6,*) 'Qm                      = ',Sonic_Mass_flow (Pi0,Ti0,Scol) 

write(6,100) 'Question 4-1'
write(6,*) ' Si Pi > Pi0 cas précédent et/ou'
write(6,*) 'Pexit > Pa'

write(6,100) 'Question 4-2'

Me=2.d0  ! pour avoir P2/P1 = 4.5
!vérification
write(6,*) 'Me                      = ',Me
write(6,*) 'P2/P1 à travers le choc = ',P2_P1(Me)
Pexit=Pa/P2_P1(Me)
write(6,*) 'P exit                  = ',Pexit
omega=P_Pi(Me); 
write(6,*) 'P/Pi                    =',omega
Pi0=Pexit/omega;   
write(6,*) 'Pi_O                    =',Pi0

write(6,100) 'Question 4-3'
write(6,*) 'Mach exit avant le choc = ', Me
write(6,*) 'Mach exit après le choc = ', Mach_aval(Me)


write(6,100) 'Question 4-4'
write(6,*) ' On diminue de 20% Pi0 '
Pi0=0.80d0*Pi0
Rchoc=Rcol+ 0.90d0*(Rexit-Rcol)
write(6,*)  'supposons que le choc se trouve dans la section de rayon =  ',Rchoc
rapport_section=(Rchoc/Rcol)**2
Mchoc=inverse_A_Acritic(1.5d0,rapport_section,1,1,45,1.d-12) 
write(6,*) 'Mach amont              = ', Mchoc
Maval=Mach_aval(Mchoc)
write(6,*) 'Mach aval               = ', Maval
Rcrit=Rchoc/A_Acritic(Maval)
write(6,*) 'Nouvelle section critique= ', Rcrit
rapport_section=(Rexit/Rcrit)**2
Me=inverse_A_Acritic(0.5d0,rapport_section,1,1,45,1.d-12) 

rapport_section=(Rexit/Rcrit)**2
Me=inverse_A_Acritic(0.5d0,rapport_section,1,1,45,1.d-12) 

Pi2=Pi0*Pi2_Pi1(Mchoc);
write(6,*) 'Pi derrière le choc    = ',Pi2
Pexit=Pi2*P_Pi(Me)
write(6,*) 'Pexit                  = ', Pexit
write(6,*) 'Pexit /Pa              = ', Pexit/Pa

if (Pexit/Pa .ge.1.d0) then 
        write(*,*) 'il faut augmenter Rchoc'
        write(*,*) 'le choc se trouve plus proche de la sortie'
else 
        write(*,*) 'il faut diminuer Rchoc'
        write(*,*) 'le choc se trouve plus proche du col'
end if




100  format(/,50('#'),/,3x,a,/,50('#'),/)
110  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine exo_Laval


!******************************************************************
subroutine Exercice_10_3
!******************************************************************
!>@details Ecoulement de Fanno subsonique
! Livre Aérodynamique Fondamentale, exercice 10.3
implicit none
real(kind=8)                ::  cf=0.008d0
real(kind=8)                ::  Dh, M1,M2,RHS,L,Lcrit,F1,F2
real(kind=8)                ::  r_T1,r_P1,r_rho1,r_Pi1
real(kind=8)                ::  r_T2,r_P2,r_rho2,r_Pi2
real(kind=8)                ::  haut,larg,S,c
real(kind=8)                ::  L1crit,fM1,fM1_new,M1_new



if (fichier) open(6,form='formatted',file='Solution_10-3.out')
write(6,110) 'Exercice Fanno 10.3'

!************************************************
write(6,100) 'Ecoulement subsonique à section circulaire'
!************************************************

M1=0.30d0; L=15.d0;Dh=0.1d0
! calcul de L
write(6,*) 'M1                      = ',M1
write(6,*) 'L (m)                   = ',L
write(6,*) 'Dh (m)                  = ',Dh
F1=Fanno(M1)
RHS=4.d0*Cf*L/Dh;
Lcrit=Sonic_Length(M1)*Dh/(4*Cf)
write(6,*) 'F(M1)                   = ',F1
write(6,*) 'RHS                     = ',RHS
write(6,*) 'Lcrit                   = ',Lcrit

call sonic_ratio_Fanno(M1,r_T1,r_P1,r_rho1,r_Pi1)

F2=F1-RHS
if (Lcrit.gt. L) then
    write(6,*) ' Ecoulement totalement subsonique'
    M2=inverse_Fanno(0.1d0,F2)
    write(6,*) 'F(M2)                   = ',F2
    write(6,*) 'M2                      = ',M2
else
    write(6,*) 'Lcrit < L :  ==> nouveau calcul à faire'
end if

call sonic_ratio_Fanno(M2,r_T2,r_P2,r_rho2,r_Pi2)

write(6,*) 'P2/P1                   = ',r_P2/r_P1
write(6,*) '(P_i2-P_i1)/P_i1        = ',r_Pi2/r_Pi1-1.d0

!************************************************
write(6,100) 'Ecoulement subsonique à section rectangulaire'
!************************************************
M1=0.5d0
! si section rectangulaire
larg=0.02d0;haut=0.01d0;
S=larg*haut;c=2.d0*(larg+haut)
Dh=4.d0*S/c
 ! si section cylindrique
Dh=0.2d0
write(6,*) 'M1                      = ',M1
write(6,*) 'L (m)                   = ',L
write(6,*) 'Dh (m)                  = ',Dh

F1=Fanno(M1)
write(6,*) 'F(M1)                   = ',F1
write(6,*) 'F(M1)                   = ',Sonic_Length(M1)
L1crit=Sonic_Length(M1)*Dh/(4.d0*cf)
write(6,*) 'L1 critique             = ',L1crit
!L=2.16d0*L1crit
L=20.d0
write(6,*) 'L2                      = ',L
F2=L/Dh*4.d0*cf
write(6,*) 'F2                      = ',F2
M1_new=inverse_Fanno(0.1d0,F2)
write(6,*) 'Nouvel M1               = ',M1_new
call sonic_ratio_Fanno(M1_new,r_T1,r_P1,r_rho1,r_Pi1)
! le point 2 est le point critique
write(6,*) 'P2/P1                   = ',1.d0/r_P1
! Ti et Pi sont inchangés
fM1=qm_Mach(M1)
write(6,*) 'Qm(M1) propto           = ',fM1
fM1_new=qm_Mach(M1_new)
write(6,*) 'Qm(M1_new) propto       = ',fM1_new

write(6,*) 'Variation relative qm   = ',fM1_new/fM1-1.d0

100  format(/,50('#'),/,3x,a,/,50('#'),/)
110  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine exercice_10_3


!******************************************************************
subroutine Exercice_10_4
!>@details Ecoulement de Fanno supersonique
!******************************************************************
! Livre Aérodynamique Fondamentale, exercice 10.4
implicit none
real(kind=8)                ::  cf=0.008d0
real(kind=8)                ::  Dh, M1,M2,L,F1
real(kind=8)                ::  L1crit
real(kind=8)                ::  Mx,L2,L2crit,L1



if (fichier) open(6,form='formatted',file='Solution_10-4.out')
write(6,110) 'Exercice Fanno 10.4'

!************************************************
write(6,110) 'Ecoulement supersonique choké, exercice 10.4'
!************************************************
Dh=1.d0; L=30.d0*Dh; 
M1=3.0d0
write(6,*) 'M1                      = ',M1
write(6,*) 'L (m)                   = ',L
write(6,*) 'Dh (m)                  = ',Dh

F1=Fanno(M1)
write(6,*) 'F(M1)                   = ',F1
L1crit=Sonic_Length(M1)*Dh/(4.d0*cf)
write(6,*) 'L1 critique             = ',L1crit
if (L1crit.gt. L) then
    write(6,*) ' Ecoulement totalement supersonique, sans choc'
    stop
else
    write(6,*) 'Lcrit < L :  ==> Ecoulement supersonique avec choc'
end if

Mx=1.5d0
write(6,*) 'Mx                      = ',Mx
L1=Dh/(4.d0*cf)*(Fanno(M1)-Fanno(Mx))
write(6,*) 'L1                      = ',L1
L2=L-L1;
M2=Mach_aval(Mx);
write(6,*) 'Mach aval               = ',M2
L2crit=Sonic_Length(M2)*Dh/(4.d0*cf)
write(6,*) 'L2 critique             = ',L2crit
write(6,*) 'L1+L2crit               = ',L2crit+L1
write(6,*) 
Mx=2.50d0
write(6,*) 'Mx                      = ',Mx
L1=Dh/(4.d0*cf)*(Fanno(M1)-Fanno(Mx))
write(6,*) 'L1                      = ',L1
L2=L-L1;
M2=Mach_aval(Mx);
write(6,*) 'Mach aval               = ',M2
L2crit=Sonic_Length(M2)*Dh/(4.d0*cf)
write(6,*) 'L2 critique             = ',L2crit
write(6,*) 'L1+L2crit               = ',L2crit+L1
write(6,*) 
Mx=2.306d0
!write(*,*) 'entrer Mx = '
!read(*,*) Mx
write(6,*) 'Mx                      = ',Mx
L1=Dh/(4.d0*cf)*(Fanno(M1)-Fanno(Mx))
write(6,*) 'L1                      = ',L1
L2=L-L1;
M2=Mach_aval(Mx);
write(6,*) 'Mach aval               = ',M2
L2crit=Sonic_Length(M2)*Dh/(4.d0*cf)
write(6,*) 'L2 critique             = ',L2crit
write(6,*) 'L1+L2crit               = ',L2crit+L1

110  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine exercice_10_4

!******************************************************************
subroutine Exercice_10_5
!******************************************************************
!>@details Ecoulement de Rayleigh en sub et supersonique
 
real(kind=8)        :: Ti0,T0,P0,Ra0,hi0,M0
real(kind=8)        :: Tic,Tc,Pc,qmax
real(kind=8)        :: rTic,rTi0,rT,rP,Cp
 
 
if (fichier) open(6,form='formatted',file='Solution_10-5.out')
!************************************************
write(6,100) 'Ecoulement de Rayleigh, exercice 10.5'
!************************************************


P0=1.d5         ! Pression en Pa en entrée
T0=300.d0       ! température statique en entrée
M0=2.d0         ! Mach en entrée
Cp=gam*r/(gam-1.d0)

call rayleigh_flow_case(2.d0)
call rayleigh_flow_case(0.2d0)


100  format(/,50('*'),/,3x,a,/,50('*'),/)

if (fichier) close(6)
contains

!================================
subroutine rayleigh_flow_case(M0)
!================================
!>@details écoulement de Rayleigh : résolution pour un Mach donné
implicit none
real(kind=8),intent(in) :: M0
real(kind=8)            :: M_init,M0prime,rp0,P0prime,T0prime,rT0


write(6,110) 'Q1: quantité de chaleur  qmax'

rTi0=T_Ti(M0)
Ti0=T0/rTi0
hi0=Cp*Ti0
Ra0=Rayleigh(M0,1)
Tic=Ti0/Ra0       
qmax=(1.d0/Ra0-1.d0)*hi0

write(6,*) 'M amont                  = ',M0
write(6,*) 'P amont                  = ',P0
write(6,*) 'T amont                  = ',T0
write(6,*) 'T0/Ti0                   = ',rTi0
write(6,*) 'Ti0                      = ',Ti0
write(6,*) 'Ra0                      = ',Ra0
write(6,*) 'Tic/Ti0 =1/Ra0           = ',1.d0/Ra0
write(6,*) 'Tic                      = ',Tic
write(6,*) 'qmax  (kJ)               = ',qmax/1000.d0

write(6,110) 'Q2: Pc et Tc'


rP=(1.d0+gam*M0**2)/(gam+1.d0)  ! Pc/P0  
rT=(rP/M0)**2                   ! Tc/T0
Tc=rT*T0
Pc=rP*P0 

write(6,*) 'Pc/P0                    = ',rP
write(6,*) 'Tc/T0                    = ',rT
write(6,*) 'Tc                       = ',Tc
write(6,*) 'Pc                       = ',Pc

if (M0.gt.1.d0) then
    write(6,110) 'Q3: qmax=2 qmax, supersonique'
    M_init=2d0
else
    write(6,110) 'Q3: qmax=2 qmax, subsonique'
    M_init=0.5d0
end if
qmax=2.d0*qmax
rTic=1.d0+qmax/hi0
Ra0=1.d0/rTic
M0prime=inverse_Rayleigh(M_init,Ra0)

write(6,*) 'qmax (kJ)                = ',qmax/1000.d0
write(6,*) 'Tic/Ti0                  = ',rTic
write(6,*) 'Ra0                      = ',Ra0
write(6,*) 'M amont (nouveau)        = ',M0prime

rp0=P_Pi(M0prime)/P_Pi(M0)
P0prime=rp0*P0
rT0=T_Ti(M0prime)/T_Ti(M0)
T0prime=T0*rT0
rP=(1.d0+gam*M0prime**2)/(gam+1.d0)  ! Pc/P0  
rT=(rP/M0prime)**2                   ! Tc/T0
Tc=rT*T0prime
Pc=rP*P0prime 
write(6,*) 'T/Ti0                    = ',T_Ti(M0prime)
write(6,*) 'Tic                      = ',Ti0/Ra0       
write(6,*) 'M0,M0prime               = ',M0,M0prime
write(6,*) 'P0prime/P0               = ',rP0
write(6,*) 'T0prime/T0               = ',rT0
write(6,*) 'Pc/P0                    = ',rP
write(6,*) 'Tc/T0                    = ',rT
write(6,*) 'Tc                       = ',Tc
write(6,*) 'Pc                       = ',Pc
write(6,*) 'P0 prime                 = ',P0prime
write(6,*) 'T0 prime                 = ',T0prime



110  format(/,50('#'),/,3x,a,/,50('#'),/)
end subroutine rayleigh_flow_case


end subroutine Exercice_10_5



!******************************************************************
subroutine Exercice_10_6
!******************************************************************
!>@details Ecoulement de Rayleigh + combustion

 
real(kind=8)        :: Ti0,T0,P0,Ra0,hi0,M0,a0,u0,rho0,S0
real(kind=8)        :: Ti1,Ra1,qm,qmax,q,rho1,Pi0,Pi1
real(kind=8)        :: Cp,rT,rTi,rP,T1,P1,u1,rPi,ru 
 
 
if (fichier) open(6,form='formatted',file='Solution_10-6.out')
!************************************************
write(6,100) 'Ecoulement de Rayleigh + Combustion, exercice 10.6'
!************************************************

u0=80.d0         ! vitesse à l'entrée
P0=1.4d4         ! Pression en Pa en entrée
T0=290.d0        ! température statique en entrée
q=800.d3         ! Dégagement de chaleur q
Cp=gam*r/(gam-1.d0)
S0=1.d0

write(6,*) 'u0                       = ',u0
write(6,*) 'P0                       = ',P0
write(6,*) 'T0                       = ',T0
write(6,*) 'q                        = ',q


write(6,110) 'Q1: Ti0 et Ti1'

a0=sqrt(gam*r*T0)
M0=u0/a0
rT=T_Ti(M0)
Ti0=T0/rT
Ra0=Rayleigh(M0,1)
hi0=Cp*Ti0
rTi=1.d0+q/hi0
Ti1=rTi*Ti0
! ou bien  Ti1=Ti0+q/Cp

write(6,*) 'a0                       = ',a0  
write(6,*) 'M0                       = ',M0
write(6,*) 'T0/Ti0                   = ',rT
write(6,*) 'Ti0                      = ',Ti0
write(6,*) 'Ra0                      = ',Ra0
write(6,*) 'Ti1/Ti0                  = ',rTi
write(6,*) 'Ti1                      = ',Ti1,Ti0+q/Cp
  
write(6,110) 'Q2: M1 et P1'

Ra1=Ra0*rTi
M1=inverse_Rayleigh(0.5d0,Ra1)
rP=(1.d0+gam*M0**2)/(1.d0+gam*M1**2)
rT=(M1/M0 *rP)**2
P1=rP*P0
T1=rT*T0

write(6,*) 'Ra1                      = ',Ra1
write(6,*) 'M1                       = ',M1
write(6,*) 'P1/P0                    = ',rP
write(6,*) 'T1/T0                    = ',rT
write(6,*) 'P1                       = ',P1
write(6,*) 'T1                       = ',T1



write(6,110) 'Q3: q_m'

rho0=P0/(r*T0)
qm=rho0*u0*S0
write(6,*) 'rho_0                    = ',rho0
write(6,*) 'qm                       = ',qm

!********************************************* 
write(6,110) 'Q4: q_max'
!********************************************* 

qmax=(1.d0/Ra0-1.d0)*hi0
write(6,*) 'qmax                     = ',qmax
Ti1=Ti0/Ra0
write(6,*) 'Ti1                      = ',Ti1

!********************************************* 
write(6,110) 'Q5: q = 1600'
!********************************************* 

M1=1.d0
qmax=1600.d3
rTi=1.d0+qmax/hi0
Ra0=1.d0/rTi
M0=inverse_Rayleigh(0.5d0,Ra0)
Ti1=rTi*Ti0

write(6,*) 'Ti0                      = ',Ti0
write(6,*) 'Ti1/Ti0                  = ',rTi
write(6,*) 'Ra0                      = ',Ra0
write(6,*) 'M0                       = ',M0
write(6,*) 'Ti1                      = ',Ti1

!********************************************* 
write(6,110) 'Q6: u0,T0,rho0'
!********************************************* 
 
T1=Ti1*T_Ti(M1)                   ! M1=1
u1=M1*sqrt(gam*r*T1)

write(6,*) 'T1                       = ',T1
write(6,*) 'u1                       = ',u1
write(6,*) 

T0=Ti0*T_Ti(M0)                   
u0=M0*sqrt(gam*r*T0)
rho0=P0/(r*T0)

write(6,*) 'Ti0                      = ',Ti0
write(6,*) 'T0                       = ',T0
write(6,*) 'P0                       = ',P0
write(6,*) 'rho0                     = ',rho0
write(6,*) 'u0                       = ',u0
write(6,*) 'débit                    = ',rho0*u0




write(6,*)  

rP=(1.d0+gam*M0**2)/(gam+1.d0)  ! P1/P0 (critique) 
ru=rP*(M1/M0)**2
rT=(M1/M0 *rP)**2
P1=rP*P0
rho1=P1/(r*T1) 

write(6,*) 'P1/P0                    = ',rP
write(6,*) 'T1/T0                    = ',rT,T1/T0
write(6,*) 'P1                       = ',P1
write(6,*) 'rho1                     = ',rho1
write(6,*) 'débit                    = ',rho1*u1
write(6,*) 

rPi=rP*P_Pi(M0)/P_Pi(1.d0)
Pi0=P0/P_Pi(M0)
Pi1=rPi*Pi0


write(6,*) 'Pi0                      = ',Pi0
write(6,*) 'Pi1/Pi0                  = ',rPi 
write(6,*) 'Pi1                      = ',Pi1,P1/P_Pi(M1)
write(6,*) 'vérification des vitesses = ',u1/u0,ru,rho0/rho1



100  format(/,50('*'),/,3x,a,/,50('*'),/)
110  format(/,50('#'),/,3x,a,/,50('#'),/)
if (fichier) close(6)
end subroutine Exercice_10_6
!******************************************************************
subroutine Exercice_10_7
!******************************************************************
!>@details Pitot subsonique
real(kind=8)        :: M0,rho0,P0,Ti0,Pi0,rP0,alpha,a0,u0,T0

if (fichier) open(6,form='formatted',file='Solution_10-7.out')
!************************************************
write(6,100) 'Tube de Pitot en subsonique, exercice 10.7'
!************************************************
P0=26400.d0; Pi0=41000.d0;Ti0=253.d0;
rP0=P0/Pi0
alpha=(1-gam)/gam;
M0=sqrt(2.d0/(gam-1.d0)*(rP0**alpha -1.d0))
write(6,*) 'Mach M0                 = ',M0
T0=Ti0*T_Ti(M0)
write(6,*) 'T0                      = ',T0
a0=sqrt(gam*r*T0);
u0=M0*a0;
write(6,*) 'u0                      = ',u0
rho0=P0/(r*T0)
write(6,*) 'rho0                    = ',rho0


100  format(/,50('*'),/,3x,a,/,50('*'),/)

if (fichier) close(6)
end subroutine Exercice_10_7


!******************************************************************
subroutine Exercice_10_8
!******************************************************************
!>@details Pitot supersonique
real(kind=8)        :: M0,rho0,P0,Ti0,Pi0,rP0,alpha,a0,u0,T0,rPi0
real(kind=8)        :: M2,Pi2,P2,Ti2,rp2,rT0

if (fichier) open(6,form='formatted',file='Solution_10-8.out')
!************************************************
write(6,100) 'Tube de Pitot en supersonique, exercice 10.8'
!************************************************
P2=3.6d4; Pi2=4.51d4;Ti2=390.d0;
rP2=P2/Pi2
alpha=(1-gam)/gam;
!M2=sqrt(2.d0/(gam-1.d0)*(rP2**alpha -1.d0))
M2=Inverse_P_Pi(rP2)
write(6,*) 'Mach M2                 = ',M2
M0=Mach_Aval(M2) ! la relation inverse est égale à elle meme !
write(6,*) 'Mach M0                 = ',M0
Ti0=Ti2
rT0=T_Ti(M0)
write(6,*) 'T0/Ti0                  = ',rT0
T0=Ti0*rT0
write(6,*) 'T0                      = ',T0
a0=sqrt(gam*r*T0);
write(6,*) 'a0                      = ',a0
u0=M0*a0;
write(6,*) 'u0                      = ',u0

rPi0=Pi2_Pi1(M0)
write(6,*) 'Pi1/Pi0                 = ',rPi0
Pi0=Pi2/rPi0
write(6,*) 'Pi0                     = ',Pi0
rP0=P_Pi(M0)
write(6,*) 'P0/Pi0                  = ',rP0
P0=P_Pi(M0)*Pi0
write(6,*) 'P0                      = ',P0
rho0=P0/(r*T0)
write(6,*) 'rho0                    = ',rho0

100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_10_8

!******************************************************************
subroutine Exercice_10_9
!******************************************************************
!>@details Régimes dans une tuyere

real(kind=8)        :: Ti0,Pi0,rPs
real(kind=8)        :: M2,rP2,rPi2,P2,Pi2,Sc2
real(kind=8)        :: Sc,Ss,Ms,Ps,rapport_section,qm
real(kind=8)        :: S1,rP1,M1,Pa,sigma2,Mns,tmp,theta
real(kind=8)        :: Om2,Oms,deltaTheta

if (fichier) open(6,form='formatted',file='Solution_10-9.out')
!************************************************
write(6,100) 'Régime dans une tuyère, exercice 10.9'
!************************************************

Pi0=1.d6      ! Pression en bar
Ti0=600.d0
Sc=5.d-4        ! section critique en cm^2
Ss=20.d-4       ! section de sortie en cm^2


rapport_section=Ss/Sc
write(6,110) 'Q1 : Régime C'
Ms= inverse_A_Acritic(2.d0,rapport_section,1,1,45,1.d-12) 
write(6,*) 'Ms                      = ',Ms
rPs=P_Pi(Ms);
write(6,*) 'Ps/Pis                  = ',rPs
Ps=Pi0*rPs
write(6,*) 'Ps                      = ',Ps

write(6,110) 'Q2: débit massique qm'
qm=Sonic_Mass_flow (Pi0,Ti0,Sc)
write(6,*) 'qm                      = ',qm


write(6,110) 'Q3: Régime B subsonique'
Ms= inverse_A_Acritic(0.5d0,rapport_section,1,1,45,1.d-12)
write(6,*) 'Ms                      = ',Ms
rPs=P_Pi(Ms);
write(6,*) 'Ps/Pis                  = ',rPs
Ps=Pi0*rPs
write(6,*) 'Ps                      = ',Ps

write(6,110) 'Q4: Régime D supersonique choc en sortie'
Ms= inverse_A_Acritic(2.d0,rapport_section,1,1,45,1.d-12) 
M2=Mach_Aval(Ms)
rP2=P_Pi(M2);
rPi2=Pi2_Pi1(Ms);
Pi2=Pi0*rPi2;
P2=Pi2*rP2
 
write(6,*) 'Mach amont               = ',Ms
write(6,*) 'Pression amont           = ',P_Pi(Ms)*Pi0
write(6,*) 'Mach aval                = ',M2
write(6,*) 'P/Pi aval                = ',rP2
write(6,*) 'P aval/P amont           = ',P2_P1(Ms)
write(6,*) 'Pi aval/Pi amont         = ',rPi2
write(6,*) 'P sortie                 = ',P2
write(6,*) 'Pi sortie                = ',Pi2
 

write(6,110) 'Q5: Régime D prime supersonique dans la tuyere '
write(6,*) ' entrer S1 en cm2'
!read*,S1
S1=S1*1.d-4
S1=12.2d-4
M1= inverse_A_Acritic(1.5d0,S1/Sc,1,1,45,1.d-12)
M2= Mach_Aval(M1)
write(6,*) 'Section amont            = ',S1
write(6,*) 'Mach amont               = ',M1
write(6,*) 'Mach aval                = ',M2
rP1=P2_P1(M1)
write(6,*) 'P aval/P amont           = ',rP1
rapport_section=A_Acritic(M2)
Sc2=S1/rapport_section
write(6,*) 'Nouvelle Section Critique= ',Sc2
M2= inverse_A_Acritic(0.5d0,Ss/Sc2,1,1,45,1.d-12)
write(6,*) 'Mach sortie              = ',M2
rPi2=Pi2_Pi1(M1);
write(6,*) 'Pi aval/Pi amont         = ',rPi2

Pa=5.d5
Ps=P_Pi(M2)*rPi2*Pi0/Pa
write(6,*) ' = 1 ?                   = ',Ps

write(6,110) 'Q5: Pa= 100 000, choc oblique, tuyère sous détendue'
P2=1.d5
! on a besoin des résultats de la question 1 :
write(6,*) 'rappel des résultats de la question 1 :'
rapport_section=Ss/Sc
Ms= inverse_A_Acritic(2.d0,rapport_section,1,1,45,1.d-12) 
write(6,*) 'Ms                      = ',Ms
Ps=Pi0*P_Pi(Ms)
write(6,*) 'Ps                      = ',Ps
rP2=P2/Ps
write(6,*) 'P aval/P amont           = ',rP2
Mns=Inverse_P2_P1(rP2)
write(6,*) 'Mn sortie amont          = ',Mns
sigma2=asin(Mns/Ms)
write(6,*) 'Angle de choc            = ',sigma2/val_pi*180.d0

tmp=(2.d0+(gam-1.d0)*Mns**2)/((gam+1.d0)*Mns**2)*tan(sigma2)
theta=(sigma2-atan(tmp))/val_pi*180.d0
write(6,*) 'Angle de deviation        = ',theta,atan(tmp)/val_pi*180.d0

write(6,110) 'Q6: 20 000 Pa, tuyère sur détendue '
P2=2.d4;
rP2=P2/Ps
M2=Inverse_P_Pi(P2/Pi0)
write(6,*) 'Ms                       = ',Ms
write(6,*) 'P aval/P amont           = ',rP2
write(6,*) 'M aval détente           = ',M2
Om2=omega(M2);Oms=omega(Ms)
deltaTheta=Om2-Oms
write(6,*) 'omega sortie  amont      = ',Oms/val_pi*180.d0
write(6,*) 'omega sortie  aval       = ',Om2/val_pi*180.d0
write(6,*) 'delta theta              = ',deltaTheta/val_pi*180.d0


100  format(/,50('*'),/,3x,a,/,50('*'),/)
110  format(/,50('#'),/,3x,a,/,50('#'),/)
if (fichier) close(6)
end subroutine Exercice_10_9

end module exercices_chap10


