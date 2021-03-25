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


!> MODULE :    Solution of chapter 11 exercises
!! 
!!
!!
module exercices_chap12

use mod_cat
use Constante

contains

!******************************************************************
subroutine Exercice_12_3
!******************************************************************
!>@details Shock tube application  
implicit none
real(kind=8)        :: p1,p4,T1,T4,P41,P21,P52,P51
real(kind=8)        :: a1,u2,a2,w2,a4,u3,pente3
real(kind=8)        :: Mc,Wc,alpha,tmp,MR,WR,M2
real(kind=8)        :: L1,Lc,L2,time2,time3,timeA
real(kind=8)        :: ech=1.d0

if (fichier) open(6,form='formatted',file='Solution_12-3.out')
write(6,100)'Exercice 12.3 : Shock tube application'

T1=290.d0;T4=290.d0;P1=1.d0;P4=100.d0;
L1=10.d0;Lc=8.d0;L2=5.d0
P41=P4/P1
P21=3.d0

!********************************************* 
write(6,110)'Q2  : P2/P1 ? '
!********************************************* 

call solve_p2_p1(P41,P21)

!===============================
write(6,110)'Shock 1 to 2 '
!===============================

print*,'ratio P4/P1              = ',P41
print*,'ratio P2/P1 =P3/P1       = ',P21
print*,'ratio P3/P4              = ',P21/P41
! upstream Mach of the normal shock
Mc=Inverse_P2_P1(P21)
print*,'M1=Mc (Mach upstream)     = ',Mc
a1=sqrt(gam*r*T1)
print*,'a1                       = ',a1
Wc=a1*Mc;
print*,'Wc                       = ',Wc
print*
w2=velocity_downstr(a1,Mc)
print*,'w2=w3 (absolute base)    = ',w2
print*
tmp=1.d0+2.d0*(gam-1.d0)/(gam+1)**2*(Mc**2-1.d0)*(gam+1.d0/Mc**2)
print*,'T2/T1                    = ',tmp
print*,'a2/a1                    = ',sqrt(tmp)
a2=a1*sqrt(tmp)
print*,'a2                       = ',a2
print*,'M2=w2/a2                 = ',w2/a2
u2=Wc-w2
print*,'(Wc-w2)/a2 = u2/a2       = ',u2/a2
print*,'u2=Wc-w2 (mobile base)   = ',u2

!======================================
write(6,110)'Reflected shock  2 to 5 '
!======================================

alpha=Mc/(Mc**2-1.d0)*sqrt(tmp)
alpha=1.d0/alpha
MR=(alpha+sqrt(alpha**2+4.d0))/2.d0
print*,'MR                       = ',MR
print*
M2=2.d0/(gam+1.d0)*(MR-1/MR)
print*,'M2  (absolute base)      = ',M2
WR=a2*MR-w2
print*,'WR = a2 MR - w2          = ',WR
print*,'WR papier                = ',a1*(1+2*(gam-1)/(gam+1)*(Mc**2-1))/Mc 
tmp=(1+(gam-1)/(gam+1)*(Mc**2-1))/(1+2*gam/(gam+1)*(Mc**2-1))
print*,'M2(paper, mobile base, upstream shock)   = ',sqrt(tmp)
print*
P52=p2_p1(MR)
print*,'ratio P5/P2            = ',P52
! paper
alpha=(gam+1.d0)/(gam-1.d0)
P52=((alpha+2.d0)*P21-1)/(P21+alpha)
print*,'ratio P5/P2 paper      = ',P52
P52=((3*gam-1)*P21+1-gam)/((gam-1)*P21+1+gam)
print*,'ratio P5/P2 paper      = ',P52

P51=P52*P21
print*,'ratio P5/P1            = ',P51

P51=(2*gam*MC**2+1-gam)/(gam+1)*(-2*(gam-1)+MC**2*(3*gam-1))/&
     (2+MC**2*(gam-1))
print*,'ratio P5/P1 paper      = ',P51

write(6,110)'Time from the sensor '

timeA=L1/Wc
time2=Lc/Wc
time3=timeA+(L1-Lc)/WR
print*,'time t is ms, scale      = ',ech
print*,'Time  A   (ms)           = ',timeA*ech
print*,'Time  2   (ms)           = ',time2*ech
print*,'Time  3   (ms)           = ',time3*ech
print*,'Time  2 / Time A         = ',time2/timeA
print*,'Time  3 / Time A         = ',time3/timeA

100  format(/,50('*'),/,3x,a,/,50('*'),/)
110  format(/,50('#'),/,3x,a,/,50('#'),/)

!======================================
write(6,110)'expansions zone 4 and 3 '
!======================================
u3=w2;a4=sqrt(r*gam*T4);

print*,'a4,a4^-1(°)               = ',a4,atan(ech/a4)*rad2deg
pente3=a4-(gam-1.d0)/2.d0*u3
print*,'a4- (gam-1)/2 u3          = ',pente3
print*,'(a4- (gam-1)/2 u3)^-1 (°) = ',atan(ech/pente3)*rad2deg
print*,'Isobar line'
print*,'u3, u3^-1 (°)             = ',u3,atan(ech/u3)*rad2deg
print*,'1/Wc (°)                  = ',atan(ech/Wc)*rad2deg
print*,'1/Wr (°)                  = ',atan(ech/WR)*rad2deg
!
! to plot in x-t with a right scaling
open(1,form='formatted',file='points.dat')
print*,'Point for plotting'
write(1,*)0,0
write(1,*)L1,timeA
write(1,*)Lc,time3  
write(1,*)2*L1/3,timeA+L1/3/WR
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
    function velocity_downstr(a_amont,Mc) result(res)
    !>@details
    ! Across the normal shock
    real(kind=8),intent(in) ::a_amont,Mc
    real(kind=8)            :: res
    res=a_amont*2.d0/(gam+1)*(Mc-1d0/Mc)
    end function velocity_downstr

    !*********************************
    function p4_p1(x) result(ratio)
    !*********************************
    !>@details pressure ration p_4 / p_1
    real(kind=8),intent(in) :: x   ! p_2/p_1
    real(kind=8)            :: ratio,Num,DeNom,rT
    rT=T4/T1
    Num=(gam-1.d0)*rT*(x-1.d0)
    Denom=2.d0*gam*((gam+1.d0)*x+gam-1.d0)
    ratio=x*(1.d0-Num/sqrt(Denom))**(2.d0*gam/(1.d0-gam))
    end function p4_p1

    !********************************************************
    subroutine solve_p2_p1(cible,x)
    !********************************************************
    !>@details search for a given pressure ratio with a  Newton method
    implicit none
    real(kind=8),intent(inout)  :: x
    real(kind=8),intent(in)     :: cible
    real(kind=8),parameter      :: dx=0.01d0
    real(kind=8),parameter      :: erreur_max=1.d-4
    integer,parameter           :: iterMax=20
    integer                     :: iter
    real(kind=8)                :: gradient,variation,valeur

    write(6,100) "Newton method"
    valeur=p4_p1(x)
    iter=0;
    variation=1.d0
    show=.false.
    do while ((abs(variation).gt.erreur_max).and.(iter.lt.iterMax))
        gradient=(p4_p1(x+dx)-valeur)/dx
        variation=-(valeur-cible)/gradient
        x=x+variation
        valeur= p4_p1(x)
        iter=iter+1
        write(*,120)'ITER = ', iter,variation,x
    end do

    if (iter.ne.iterMax) then
        print*,'convergence reached'
    else
        stop 'no convergence'
    end if
    !valeur= p4_p1(x)
    
100  format(/,50('.'),/,3x,a,/,50('.'),/)
120  format(3x,a6,2x,i2,3x,2(f13.6,3x))
    end subroutine solve_p2_p1

end subroutine Exercice_12_3

end module exercices_chap12
