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


!******************************************************************
subroutine TD1_N7
!******************************************************************
implicit none

real(kind=8)        :: Rinlet=0.4d0,Rexit=0.1d0,Rcol=0.08d0
real(kind=8)        :: Ti0=350.d0,Pi0,Pa=101325.d0
real(kind=8)        :: Me,omega,rapport_section,Rcrit,Mcol,Scol,Pexit,Sexit
real(kind=8)        :: Minlet,Tinlet,ainlet,Uinlet,rhoinlet,Pinlet
real(kind=8)        :: Qm,Pi2
real(kind=8)        :: Rchoc,Mchoc,Maval

!**********
!EXERCICE 2
!**********

write(6,100)'Question 2-1'
Pi0=103000d0
write(6,*)'Me                      = ',inverse_P_Pi(Pa/Pi0)


write(6,100)'Question 2-2'
Me=0.3d0
omega=P_Pi(Me); 
write(6,*)'P/Pi                    = ',omega
Pi0=Pa/omega;   
write(6,*)'Pi_0                    = ',Pi0
rapport_section= A_Acritic(Me)
Rcrit= Rexit/sqrt(rapport_section)
write(6,*)'A/Ac                    = ',rapport_section
write(6,*)'R_crit                  = ',Rcrit

write(6,*)'Col :'
rapport_section=(Rcol/Rcrit)**2
Mcol=inverse_A_Acritic(0.5d0,rapport_section,1,1,45,1.d-12) 
write(6,*)'Acol/Ac                 = ',rapport_section
write(6,*)' Mcol                   = ',Mcol


write(6,*)'Entrée :'

rapport_section=(Rinlet/Rcrit)**2
Minlet=inverse_A_Acritic(0.5d0,rapport_section,1,1,45,1.d-12) 
write(6,*)'Ainlet/Ac               = ',rapport_section
write(6,*)'Minlet                  = ',Minlet

omega=T_Ti(Minlet)
write(6,*)'T/Ti                    = ',omega
Tinlet=Ti0*omega; 
write(6,*)'Tinlet                  = ',Tinlet
ainlet=sqrt(gam*r*Tinlet);
write(6,*)'a_inlet                 = ',ainlet
Uinlet=Minlet*ainlet;
write(6,*)'Uinlet                  = ',Uinlet

omega=P_Pi(Minlet)
write(6,*)'P/Pi                    = ',omega
Pinlet=Pi0*omega; 
write(6,*)'Pinlet                  = ',Pinlet
rhoinlet=Pinlet/(r*Tinlet);
write(6,*)'rho_inlet               = ',rhoinlet
Qm=rhoinlet*Uinlet*val_pi*Rinlet**2
write(6,*)'Débit massique Qm       = ',Qm

write(6,100)'Question 2-3'
Rcrit=Rcol
rapport_section=(Rexit/Rcrit)**2
Me=inverse_A_Acritic(0.5d0,rapport_section,1,1,45,1.d-12) 
omega=P_Pi(Me); 
write(6,*)'P/Pi                    =',omega
Pi0=Pa/omega;   
write(6,*)'Pi_O                    =',Pi0

!**********
!EXERCICE 3
!**********


write(6,100)'Question 3-1'

Me=2.d0
omega=P_Pi(Me); 
write(6,*)'P/Pi                    =',omega
Pi0=Pa/omega;   
write(6,*)'Pi_O                    =',Pi0
rapport_section= A_Acritic(Me)
write(6,*)'Se/Scrit                =',rapport_section
Scol=val_pi*Rcol**2
Sexit=Scol*rapport_section
Rexit=Rcol*sqrt(rapport_section)
write(6,*)'Se                      =',Sexit
write(6,*)'Re                      =',Rexit


write(6,100)'Question 3-2'
write(6,*)'Qm                      = ',Sonic_Mass_flow (Pi0,Ti0,Scol) 

!**********
!EXERCICE 4
!**********


write(6,100)'Question 4-1'
write(6,*)' Si Pi > Pi0 cas précédent et/ou'
write(6,*)'Pexit > Pa'

write(6,100)'Question 4-2'

Me=2.d0  ! pour avoir P2/P1 = 4.5
!vérification
write(6,*)'Me                      = ',Me
write(6,*)'P2/P1 à travers le choc = ',P2_P1(Me)
Pexit=Pa/P2_P1(Me)
write(6,*)'P exit                  = ',Pexit
omega=P_Pi(Me); 
write(6,*)'P/Pi                    =',omega
Pi0=Pexit/omega;   
write(6,*)'Pi_O                    =',Pi0

write(6,100)'Question 4-3'
write(6,*)'Mach exit avant le choc = ', Me
write(6,*)'Mach exit après le choc = ', Mach_aval(Me)


write(6,100)'Question 4-4'
write(6,*)' On diminue de 20% Pi0 '
Pi0=0.80d0*Pi0
Rchoc=Rcol+ 0.90d0*(Rexit-Rcol)
write(6,*) 'supposons que le choc se trouve dans la section de rayon =  ',Rchoc
rapport_section=(Rchoc/Rcol)**2
Mchoc=inverse_A_Acritic(1.5d0,rapport_section,1,1,45,1.d-12) 
write(6,*)'Mach amont              = ', Mchoc
Maval=Mach_aval(Mchoc)
write(6,*)'Mach aval               = ', Maval
Rcrit=Rchoc/A_Acritic(Maval)
write(6,*)'Nouvelle section critique= ', Rcrit
rapport_section=(Rexit/Rcrit)**2
Me=inverse_A_Acritic(0.5d0,rapport_section,1,1,45,1.d-12) 

rapport_section=(Rexit/Rcrit)**2
Me=inverse_A_Acritic(0.5d0,rapport_section,1,1,45,1.d-12) 

Pi2=Pi0*Pi2_Pi1(Mchoc);
write(6,*)'Pi derrière le choc    = ',Pi2
Pexit=Pi2*P_Pi(Me)
write(6,*)'Pexit                  = ', Pexit
write(6,*)'Pexit /Pa              = ', Pexit/Pa

if (Pexit/Pa .ge.1.d0) then 
        write(6,*) 'il faut augmenter Rchoc'
        write(6,*) 'le choc se trouve plus proche de la sortie'
else 
        write(6,*) 'il faut diminuer Rchoc'
        write(6,*) 'le choc se trouve plus proche du col'
end if




100  format(/,50('#'),/,3x,a,/,50('#'),/)
end subroutine TD1_N7


!******************************************************************
subroutine TD2_N7
!******************************************************************
implicit none
real(kind=8)                :: cf=0.005d0
real(kind=8)                ::  Dh, M1,M2,RHS,L,Lcrit,F1,F2
real(kind=8)                ::  r_T1,r_P1,r_rho1,r_Pi1
real(kind=8)                ::  r_T2,r_P2,r_rho2,r_Pi2
real(kind=8)                ::  haut,larg,S,c
real(kind=8)                ::  L1crit,fM1,fM1_new,M1_new
real(kind=8)                ::  Mx,L2,L2crit,L1



!************************************************
write(6,100)'Question 3-1'
!************************************************

M1=0.24d0; L=20.d0;Dh=0.1d0
! calcul de L
F1=Fanno(M1)
RHS=4.d0*Cf*L/Dh;
Lcrit=Sonic_Length(M1)*Dh/(4*Cf)
write(6,*)'F(M1)                   = ',F1
write(6,*)'RHS                     = ',RHS
write(6,*)'Lcrit                   = ',Lcrit

call sonic_ratio_Fanno(M1,r_T1,r_P1,r_rho1,r_Pi1)

F2=F1-RHS
if (Lcrit.gt. L) then
    write(6,*)' Ecoulement totalement subsonique'
    M2=inverse_Fanno(0.1d0,F2)
    write(6,*)'F(M2)                   = ',F2
    write(6,*)'M2                      = ',M2
else
    write(6,*)'Lcrit < L :  ==> nouveau calcul à faire'
end if

call sonic_ratio_Fanno(M2,r_T2,r_P2,r_rho2,r_Pi2)

write(6,*)'P2/P1                   = ',r_P2/r_P1
write(6,*)'(P_i2-P_i1)/P_i1        = ',r_Pi2/r_Pi1-1.d0

!************************************************
write(6,100)'Question 3-2'
!************************************************
M1=0.5d0
larg=0.02d0;haut=0.01d0;
S=larg*haut;c=2.d0*(larg+haut)
Dh=4.d0*S/c
write(6,*)'Dh                      = ',Dh
L1crit=Sonic_Length(M1)*Dh/(4.d0*cf)
write(6,*)'L1 critique             = ',L1crit
L=2.16d0*L1crit
write(6,*)'L2                      = ',L
F2=L/Dh*4.d0*cf
write(6,*)'F2                      = ',F2
M1_new=inverse_Fanno(0.1d0,F2)
write(6,*)'Nouvel M1               = ',M1_new
call sonic_ratio_Fanno(M1_new,r_T1,r_P1,r_rho1,r_Pi1)
! le point 2 est le point critique
write(6,*)'P2/P1                   = ',1.d0/r_P1
! Ti et Pi sont inchangés
fM1=qm_Mach(M1)
write(6,*)'Qm(M1) propto           = ',fM1
fM1_new=qm_Mach(M1_new)
write(6,*)'Qm(M1_new) propto       = ',fM1_new

write(6,*)'Variation relative qm   = ',fM1_new/fM1-1.d0

!************************************************
write(6,100)'Question 4-1'
!************************************************
Dh=1.d0; L=40.d0*Dh; 
M1=2.6d0
L1crit=Sonic_Length(M1)*Dh/(4.d0*cf)
write(6,*)'L1 critique             = ',L1crit
if (L1crit.gt. L) then
    write(6,*)' Ecoulement totalement supersonique, sans choc'
    stop
else
    write(6,*)'Lcrit < L :  ==> Ecoulement supersonique avec choc'
end if

Mx=1.5d0
write(6,*)'Mx                      = ',Mx
L1=Dh/(4.d0*cf)*(Fanno(M1)-Fanno(Mx))
write(6,*)'L1                      = ',L1
L2=L-L1;
M2=Mach_aval(Mx);
write(6,*)'Mach aval               = ',M2
L2crit=Sonic_Length(M2)*Dh/(4.d0*cf)
write(6,*)'L2 critique             = ',L2crit
write(6,*)'L1+L2crit               = ',L2crit+L1
write(6,*)
Mx=2.50d0
write(6,*)'Mx                      = ',Mx
L1=Dh/(4.d0*cf)*(Fanno(M1)-Fanno(Mx))
write(6,*)'L1                      = ',L1
L2=L-L1;
M2=Mach_aval(Mx);
write(6,*)'Mach aval               = ',M2
L2crit=Sonic_Length(M2)*Dh/(4.d0*cf)
write(6,*)'L2 critique             = ',L2crit
write(6,*)'L1+L2crit               = ',L2crit+L1
write(6,*)
Mx=2.128d0
write(6,*)'Mx                      = ',Mx
L1=Dh/(4.d0*cf)*(Fanno(M1)-Fanno(Mx))
write(6,*)'L1                      = ',L1
L2=L-L1;
M2=Mach_aval(Mx);
write(6,*)'Mach aval               = ',M2
L2crit=Sonic_Length(M2)*Dh/(4.d0*cf)
write(6,*)'L2 critique             = ',L2crit
write(6,*)'L1+L2crit               = ',L2crit+L1

100  format(/,50('#'),/,3x,a,/,50('#'),/)

end subroutine TD2_N7

!**************************
subroutine examen_Fev2017
!**************************
use Constante
implicit none
real(kind=8)                        :: P0,M0,tetap,Ti0
integer                             :: k

if (fichier) open(6,form='formatted',file='Solution_N7_Fev_2017.out')
write(6,100)'Examen de rattrapage, Fev 2017'

nz=4

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),zone_amont(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(type_caracteristique(0:nz))
allocate(longueur(0:nz),profil(0:nz),CL_zone(0:nz),CD_zone(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;zone_amont=0;Mn=0.d0;
com='';type_caracteristique='';

M0=3.5d0;p0=5000.d0; Ti0=300.d0;tetap=15.d0
print*,'data ok'
! définitions des zones
! zone 0 : amont
k=0
type_zone(k)    = 'Uniforme'
M(k)            = M0
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
zone_amont(k)   = 0             ! si elle même, pas de correction d'incidence
longueur(k)     = 1.d0
profil(k)       = 1
! zone 1 :  choc oblique
k=1
type_zone(k)    = 'Choc'
theta(k)        =  tetap*deg2rad
zone_amont(k)   = 0
longueur(k)     = 1.d0
profil(k)       = 1
! zone 2 : détente
k=2
type_zone(k)    = 'Detente'
theta(k)        = 0.d0
type_caracteristique(k)='-'
zone_amont(k)   = 1
longueur(k)     = 1.d0
profil(k)       = 1
! zone 3 : détente
k=3
type_zone(k)    = 'Detente'
theta(k)        =  -tetap*deg2rad
type_caracteristique(k)='-'
zone_amont(k)   = 2
longueur(k)     = 1.d0
profil(k)       = 1
! zone 4 :  choc oblique
k=4
type_zone(k)    = 'Choc'
theta(k)        =  0.d0
zone_amont(k)   = 3
longueur(k)     = 1.d0
profil(k)       = 1
call solve_case
call sorties

call coefficients_aerodynamiques  ! calcul des coefficients aérodynamiques


deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,zone_amont)
deallocate(longueur,profil,CL_zone,CD_zone)

100  format(/,50('*'),/,3x,a,/,50('*'),/)
end subroutine examen_Fev2017


