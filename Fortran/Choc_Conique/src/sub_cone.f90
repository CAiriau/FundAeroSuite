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

!*********************
function correction(k)
!*********************
real(kind=8)    :: correction
integer         :: k
integer,parameter   :: k_lim=10
real(kind=8),parameter ::coef_init=0.03d0

if (k.gt.k_lim) then
    correction=1.d0
else
    correction=coef_init+(1.d0-coef_init)/(real(k_lim-1,kind=8))**2*(real(k-1,kind=8))**2
end if

end function correction

!*******************
subroutine read_data
!*******************
open(1,form='formatted',file='cone.in')
read(1,*); 
read(1,*) option_Mach
read(1,*)n_Mach_tmp
read(1,*)delta_Mach
read(1,*)Mach_initial
read(1,*)Mach_final
read(1,*)n_points
read(1,*)dteta
read(1,*)opt_Newton
read(1,*)nc
read(1,*)error_max
read(1,*)opt_velocity
dteta=dteta/coef
close(1)
end subroutine read_data

!***************************
subroutine open_output_files
!***************************

open(80,form='formatted',file='deviation.out')
write(80,110)
110 format('#',2x,'Theta paroi en °',2x,'Angle Choc en °', &
            5x,'Mach amont')

open(81,form='formatted',file='theta_cone_max.out')
write(81,100)
100 format('#',5x,'Mach amont',5x,'Angle Choc en °',4x,'ThetaMax en °', &
         2x,'Inverse Mach cone',8x,'Kp')

open(82,form='formatted',file='Mach_aval_unitaire.out')
write(82,120)
120 format('#',2x,'Mach amont',6x,'Angle Choc en °',2x,&
        'Theta Mach = 1 en °',2x,'Inverse Mach',9x,'Kp')

open(83,form='formatted',file='Mach_cone_unitaire.out')
write(83,150)
150 format('#',2x,'Mach amont',6x,'Angle Choc en °',2x,&
        'Theta Mc = 1 en °',2x,'Mach aval',9x,'Kp')

! OK 
open(90,form='formatted',file='Kp.out')
write(90,130)
130 format('#',2x,'Theta paroi en °',9x,'Kp',11x,'Mach amont')

! OK
open(85,form='formatted',file='Inverse_Mach_Aval.out')
write(85,140)
140 format('#',2x,'Theta paroi en °',2x,&
                  'Inverse Mach Paroi',2x,'Mach amont')

end subroutine open_output_files 

!************************
subroutine set_Mach_table
!************************
select case (option_Mach)
case (1)
    n_Mach_Max=25
    allocate(table_Mach(N_Mach_Max))
    table_Mach=(/ 1.01d0, 1.05d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, &
                  1.7d0,  1.8d0, 1.9d0, 2.d0,  2.2d0, 2.4d0, 2.6d0, 2.8d0, &
                  3.0d0,  3.2d0, 3.4d0, 3.6d0, 3.8d0, 4.0d0, 4.5d0, 5.0d0, &
                  10.d0,  50.d0 /)
case DEFAULT
    if (Mach_initial.eq.0.d0) Mach_initial=2.05d0
    if (Mach_final.eq.0.d0)   Mach_final=3.5d0
    if (n_Mach_tmp.ne.0) then
        n_Mach_Max=n_Mach_tmp
        dMach=(Mach_final-Mach_initial)/dfloat(n_Mach_Max-1)
    else
        dMach=delta_Mach
        write(*,*)((Mach_final-Mach_initial)/dMach)
        n_Mach_Max=int((Mach_final-Mach_initial)/dMach)+2
    end if
    allocate(table_Mach(N_Mach_Max))
    table_Mach(1)=Mach_initial
    do i=2,n_Mach_Max
        table_Mach(i)=table_Mach(i-1)+dMach
    end do
end select
write(*,*) 'nombre de Mach : ', n_Mach_Max
write(*,'(f10.6)') table_Mach
write(*,'(i5)') int(1000*table_Mach)
end subroutine set_Mach_table

!****************************
FUNCTION carac(N_in,i_enter)
!****************************
! 
!  conversion d'un nombre entier en une chaine de caractère
!  exemple : 123 --> '123'
!  N_in : nombre de chiffres
!  i_enter : nombre entier à convertir

IMPLICIT NONE
INTEGER, INTENT(IN) :: i_enter,N_in
CHARACTER(LEN=N_in) :: carac
INTEGER             :: i,k,n

IF (N_in.le.1) THEN
        WRITE(*,*) 'Problème avec la fonction carac'
        STOP
ENDIF

n=N_in-1
i=i_enter
DO k=n,0,-1
        carac(n-k+1:n-k+1)=CHAR(48+INT(i/10**k))
        i=MOD(i,10**k)
END DO
END FUNCTION carac

!********************
function Rpchoc(mach)
!********************
!     saut de pression du choc droit
real (kind=8)::  Rpchoc,mach
Rpchoc= 2*gam/(gam+1)* mach**2- (gam-1)/(gam+1)
end function Rpchoc

!********************
function mach1(mach)
!********************
!     Mach normal aval d'un choc droit plan
real (kind=8):: mach,mach1
mach1= sqrt((1+ 0.5*(gam-1)* mach**2)/(gam*mach**2-0.5*(gam-1)))
end function mach1


!**********************************************************
subroutine integration_O1(teta_init,Angle_Paroi,V1,VP1,VS1)
!*********************************************************
implicit none
real(kind=8)        :: teta_init,Angle_Paroi
real(kind=8)        :: V1,VP1,VS1
real(kind=8)        :: teta,dteta_tmp,teta_old,VP1_old,VS1_old,V1_old,alpha

V1_old=V1; VP1_old=VP1; VS1_old=VS1
teta=teta_init
teta_old=teta
write(3,110)k, Angle_choc*coef
dteta_tmp=dteta
do  while (teta.ge.0.d0)
    if (VP1.gt.-0.05d0) dteta_tmp=dteta*0.1d0
    VP1=VP1-VS1*dteta_tmp
    V1=V1-VP1*dteta_tmp
    VS1=calcul_VS1(teta,V1,VP1)
    write(3,100) teta*coef,V1,VP1,VS1,Angle_Choc*coef,Mach0
    if (VP1.gt.0.d0) then                           ! V_theta est <0 et = 0 à la paroi
         alpha=-VP1_old/(VP1-VP1_old)
         Angle_Paroi=alpha*teta+(1.d0-alpha)*teta_old
         V1=alpha*V1+(1.d0-alpha)*V1_old
         VS1=alpha*VS1+(1.d0-alpha)*VS1_old
         if (opt_velocity.eq.1) write(3,100) Angle_paroi*coef,V1,0.d0,VS1,Angle_Choc*coef,Mach0
         write(48,*) (Angle_Paroi-teta)/teta
         !Angle_Paroi=teta                           ! le cone est atteint
         exit
    endif
    teta_old=teta
    V1_old=V1; VP1_old=VP1; VS1_old=VS1
    teta=teta-dteta_tmp
end do
100 format(6(e20.13,3x))
110 format('# theta',i4,/,'#',7x,f12.4,/,'#')
end subroutine integration_O1

!*******************************
function calcul_VS1(teta,V1,VP1)
!*******************************
implicit none
real(kind=8), intent(in)        :: teta,V1,VP1
real(kind=8)                    :: calcul_VS1
real(kind=8)                    :: num,den
num=(gam-1.d0)/2.d0*(1.d0-V1**2-VP1**2)*(2.d0*V1+VP1/(tan(teta))) -V1*VP1**2
den=(gam-1.d0)/2.d0*(1.d0-V1**2-VP1**2)-VP1**2
calcul_VS1=-num/den            ! dérivée seconde / theta de Vr
end function calcul_VS1

!***********************************
function deviation(Mach0,Angle_choc)
!***********************************
! deviation derriere un choc oblique
implicit none
real(kind=8),intent(in)         :: Mach0,Angle_choc
real(kind=8)                    :: deviation
real(kind=8)                    :: num,den
num=2.d0/tan(Angle_choc)*((Mach0*sin(Angle_choc))**2-1.d0)
den=Mach0**2.*(gam+cos(2.*Angle_choc))+2.d0
deviation=atan(num/den)
end function deviation

!**************************************************
function calcul_Mach_Normal_Aval(Mach_Normal_Amont)
!**************************************************
! Mach normal aval derrière un choc oblique
implicit none
real(kind=8),intent(in)     :: Mach_Normal_Amont
real(kind=8)                :: calcul_Mach_Normal_Aval
real(kind=8)                :: num,den
num=Mach_Normal_Amont**2.+2.d0/(gam-1.d0)
den=2.d0*gam/(gam-1.d0)*Mach_Normal_Amont**2.-1.d0
calcul_Mach_Normal_Aval=sqrt(num/den)
end function calcul_Mach_Normal_Aval

!****************************************************************
function coefficient_pression(Mach0,Mach_cone,Mach_normal_Amont)
!****************************************************************
implicit none
real(kind=8),intent(in) :: Mach0,Mach_cone,Mach_normal_Amont
real(kind=8)            :: coefficient_pression
real(kind=8)            :: a,b,c,d,rpi

! calcul du coefficient de pression Kp
a=1.+(gam-1.d0)/2.d0*Mach_Cone**2                                       ! Ti/T  pour l'aval
b=1.+(gam-1.d0)/2.d0*Mach0**2                                           ! Ti/T  pour l'amont
c=2.d0*gam/(gam+1.d0)*Mach_Normal_Amont**2-(gam-1.d0)/(gam+1.d0)        ! P_aval/P_amont
d=((gam-1)*Mach_Normal_Amont**2+2.0)/((gam+1)*Mach_Normal_Amont**2)     ! rho_amont/rho_aval
rpi=(a/b*d)**(-gam/(gam-1.d0))*c**(-1.d0/(gam-1.d0)) 
!     rpi :  P aval/ P amont= P_aval/Pi_aval* Pi aval/Pi_amont*Pi_amont/P_amont
coefficient_pression=2./gam/Mach0/Mach0*(rpi-1.d0)                        ! coefficient de pression fonction de P_aval/P_amont
end function coefficient_pression

