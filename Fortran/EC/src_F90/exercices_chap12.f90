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

 
!> MODULE :   Correction des exercices du chapitre 12
!! 
!!
!!
module exercices_chap12

use mod_EC
use Constante

contains

!******************************************************************
subroutine Exercice_12_3
!******************************************************************
!>@details application du tube à choc
implicit none
real(kind=8)        :: p1,p4,T1,T4,P41,P21,P52,P51
real(kind=8)        :: a1,u2,a2,w2,a4,u3,pente3
real(kind=8)        :: Mc,Wc,alpha,tmp,MR,WR,M2
real(kind=8)        :: L1,Lc,L2,temps2,temps3,tempsA
real(kind=8)        :: ech=1.d0

if (fichier) open(6,form='formatted',file='Solution_12-3.out')
write(6,100)'Exercice 12.3 : Tube à choc'

T1=290.d0;T4=290.d0;P1=1.d0;P4=100.d0;
L1=10.d0;Lc=8.d0;L2=5.d0
P41=P4/P1
P21=3.d0

!********************************************* 
write(6,110)'Q2  : P2/P1 ? '
!********************************************* 

call recherche_p2_p1(P41,P21)

!===============================
write(6,110)'Choc 1 vers 2 '
!===============================

print*,'rapport P4/P1            = ',P41
print*,'rapport P2/P1 =P3/P1     = ',P21
print*,'rapport P3/P4            = ',P21/P41
! Mach amont au choc droit
Mc=Inverse_P2_P1(P21)
print*,'M1=Mc (Mach amont)       = ',Mc
a1=sqrt(gam*r*T1)
print*,'a1                       = ',a1
Wc=a1*Mc;
print*,'Wc                       = ',Wc
print*
w2=vitesse_aval(a1,Mc)
print*,'w2=w3 (repere absolu)    = ',w2
print*
tmp=1.d0+2.d0*(gam-1.d0)/(gam+1)**2*(Mc**2-1.d0)*(gam+1.d0/Mc**2)
print*,'T2/T1                    = ',tmp
print*,'a2/a1                    = ',sqrt(tmp)
a2=a1*sqrt(tmp)
print*,'a2                       = ',a2
print*,'M2=w2/a2                 = ',w2/a2
u2=Wc-w2
print*,'(Wc-w2)/a2 = u2/a2       = ',u2/a2
print*,'u2=Wc-w2 (repere mobile) = ',u2


!======================================
write(6,110)'Choc reflechi   2 vers 5 '
!======================================

alpha=Mc/(Mc**2-1.d0)*sqrt(tmp)
alpha=1.d0/alpha
MR=(alpha+sqrt(alpha**2+4.d0))/2.d0
print*,'MR                       = ',MR
print*
M2=2.d0/(gam+1.d0)*(MR-1/MR)
print*,'M2  (repere absolu)      = ',M2
WR=a2*MR-w2
print*,'WR = a2 MR - w2          = ',WR
print*,'WR papier                = ',a1*(1+2*(gam-1)/(gam+1)*(Mc**2-1))/Mc 
tmp=(1+(gam-1)/(gam+1)*(Mc**2-1))/(1+2*gam/(gam+1)*(Mc**2-1))
print*,'M2(papier, repere mobile,choc amont)   = ',sqrt(tmp)
print*
P52=p2_p1(MR)
print*,'rapport P5/P2            = ',P52
! article
alpha=(gam+1.d0)/(gam-1.d0)
P52=((alpha+2.d0)*P21-1)/(P21+alpha)
print*,'rapport P5/P2 article    = ',P52
P52=((3*gam-1)*P21+1-gam)/((gam-1)*P21+1+gam)
print*,'rapport P5/P2 article    = ',P52

P51=P52*P21
print*,'rapport P5/P1            = ',P51

P51=(2*gam*MC**2+1-gam)/(gam+1)*(-2*(gam-1)+MC**2*(3*gam-1))/&
     (2+MC**2*(gam-1))
print*,'rapport P5/P1 article    = ',P51

write(6,110)'Temps sur le capteur '

tempsA=L1/Wc
temps2=Lc/Wc
temps3=tempsA+(L1-Lc)/WR
print*,'axe en t en ms, echelle   = ',ech
print*,'Temps A   (ms)           = ',TempsA*ech
print*,'Temps 2   (ms)           = ',Temps2*ech
print*,'Temps 3   (ms)           = ',Temps3*ech
print*,'Temps 2 / Temps A        = ',Temps2/TempsA
print*,'Temps 3 / Temps A        = ',Temps3/TempsA

100  format(/,50('*'),/,3x,a,/,50('*'),/)
110  format(/,50('#'),/,3x,a,/,50('#'),/)

!======================================
write(6,110)'détentes zone 4 et 3 '
!======================================
u3=w2;a4=sqrt(r*gam*T4);

print*,'a4,a4^-1(°)               = ',a4,atan(ech/a4)*rad2deg
pente3=a4-(gam-1.d0)/2.d0*u3
print*,'a4- (gam-1)/2 u3          = ',pente3
print*,'(a4- (gam-1)/2 u3)^-1 (°) = ',atan(ech/pente3)*rad2deg
print*,'Ligne de glissement'
print*,'u3, u3^-1 (°)             = ',u3,atan(ech/u3)*rad2deg
print*,'1/Wc (°)                  = ',atan(ech/Wc)*rad2deg
print*,'1/Wr (°)                  = ',atan(ech/WR)*rad2deg
!
! Pour tracer le graphique x-t à l'échelle
open(1,form='formatted',file='points.dat')
print*,'Point pour le trace'
write(1,*)0,0
write(1,*)L1,tempsA
write(1,*)Lc,temps3  
write(1,*)2*L1/3,tempsA+L1/3/WR
write(1,*)
write(1,*)0,0
write(1,*)L1/2,L1/u3
write(1,*)
write(1,*)0,0
write(1,*)-L2,L2/a4
write(1,*)
write(1,*)0,0
write(1,*)-L2,L2/pente3
close(1)


contains
    function vitesse_aval(a_amont,Mc) result(res)
    !>@details
    ! A travers un choc droit
    real(kind=8),intent(in) ::a_amont,Mc
    real(kind=8)            :: res
    res=a_amont*2.d0/(gam+1)*(Mc-1d0/Mc)
    end function vitesse_aval

    !*********************************
    function p4_p1(x) result(rapport)
    !*********************************
    !>@details
    real(kind=8),intent(in) :: x   ! p2/p1
    real(kind=8)            :: rapport,Num,DeNom,rT
    rT=T4/T1
    Num=(gam-1.d0)*rT*(x-1.d0)
    Denom=2.d0*gam*((gam+1.d0)*x+gam-1.d0)
    rapport=x*(1.d0-Num/sqrt(Denom))**(2.d0*gam/(1.d0-gam))
    end function p4_p1

    !********************************************************
    subroutine recherche_p2_p1(cible,x)
    !********************************************************
    !>@details recherche d'un rapport de pression par une méthode de Newton
    implicit none
    real(kind=8),intent(inout)  :: x
    real(kind=8),intent(in)     :: cible
    real(kind=8),parameter      :: dx=0.01d0
    real(kind=8),parameter      :: erreur_max=1.d-4
    integer,parameter           :: iterMax=20
    integer                     :: iter
    real(kind=8)                :: gradient,variation,valeur

    write(6,100) "Méthode de Newton"
    valeur=p4_p1(x)
    iter=0;
    variation=1.d0
    show=.false.                ! pour ne pas afficher le détail des calculs 
    do while ((abs(variation).gt.erreur_max).and.(iter.lt.iterMax))
        gradient=(p4_p1(x+dx)-valeur)/dx
        variation=-(valeur-cible)/gradient
        x=x+variation
        valeur= p4_p1(x)
        iter=iter+1
        write(*,120)'ITER = ', iter,variation,x
    end do

    if (iter.ne.iterMax) then
        print*,'convergence'
    else
        stop 'pas de convergence'
    end if
    !valeur= p4_p1(x)
    
100  format(/,50('.'),/,3x,a,/,50('.'),/)
120  format(3x,a6,2x,i2,3x,2(f13.6,3x))
    end subroutine recherche_p2_p1

end subroutine Exercice_12_3

end module exercices_chap12
