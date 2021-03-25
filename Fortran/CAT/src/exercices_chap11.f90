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


!> MODULE :   Solution of chapter 11 exercises
!! 
!! @author C. Airiau
!!

module exercices_chap11

use mod_cat
use Constante

contains


!******************************************************************
subroutine Exercice_11_1
!******************************************************************
!>@details flat plate in incidence

implicit none
integer             :: k

if (fichier) open(6,form='formatted',file='Solution_11-1.out')
write(6,100)'Exercice 11.1 : flat plate in incidence'

!*********************
! parameters
!*********************
nz=4                        ! 5 zones with the upstream zone
alpha=0.d0                  ! incidence of the plate (the slope are changed below)
! table allocation for state in the different zones:
allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),upstr_zone(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(longueur(0:nz),profil(0:nz),CL_zone(0:nz),CD_zone(0:nz))

! zone 0 : upstream
k=0
type_zone(k)    = 'Uniform'
M(k)            = 2.0d0
Ti(k)           = 300.d0
P(k)            = 1.d0
theta(k)        = 0.d0
upstr_zone(k)   = 0     ! no need to correct the incidence
! zone 1 : upper_wall
k=1
type_zone(k)    = 'Expansion'
theta(k)        = -alpha
upstr_zone(k)   = 0
longueur(k)     = 1.d0
profil(k)       = 1
! zone 2 : lower_wall
k=2
type_zone(k)    = 'Shock'
theta(k)        = -alpha
upstr_zone(k)   = 0
longueur(k)     = 1.d0
profil(k)       = -1
! zone 3 : upper_wall, wake
k=3
type_zone(k)    = 'Shock'
theta(k)        = 0.d0
upstr_zone(k)   = 1
! zone 4 : lower_wall,wake
k=4
type_zone(k)    = 'Expansion'
theta(k)        = 0.d0
upstr_zone(k)   = 2
theta=theta*deg2rad

!********************************************* 
write(6,110)'Q2: solutions for each zone'
!********************************************* 

call solve_case              
call display_outputs

!****************************************************
write(6,110)'Q3: lift and pressure drag coefficients '
!****************************************************

call aerodynamic_coefficients  

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,upstr_zone)
deallocate(longueur,profil,CL_zone,CD_zone)

110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_11_1

!
!******************************************************************
subroutine Exercice_11_2
!******************************************************************
!>@details diamond airfoil  
implicit none

if (fichier) open(6,form='formatted',file='Solution_11-2.out')
write(6,100)'Exercice 11.2 :  diamond airfoil '

alpha=0.d0 
call diamond_airfoil
alpha=8.d0
call diamond_airfoil
alpha=10.d0
call diamond_airfoil
alpha=12.d0
call diamond_airfoil

100  format(/,50('*'),/,3x,a,/,50('*'),/)

if (fichier) close(6)
end subroutine Exercice_11_2

!*************************
subroutine diamond_airfoil
!*************************
!>@details diamond airfoil definition

integer             :: k
real(kind=8)        :: phi
!*********************
! inputs
!*********************

nz=6                        ! 7 zones with the upstream zone
phi=10.d0
write(6,120)phi,alpha
120 format(/,'# DIAMOND AIRFOIL : 1/2 summit angle  : ',f8.2,' alpha = ',f8.2,/)
allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),upstr_zone(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(longueur(0:nz),profil(0:nz),CL_zone(0:nz),CD_zone(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;upstr_zone=0;Mn=0.d0;longueur=0.d0;profil=0;CL_zone=0.d0;CD_zone=0.d0
! zone 0 : upstream
k=0
type_zone(k)    = 'Uniform'
M(k)            = 2.0d0
Ti(k)           = 300.d0
P(k)            = 1.d0
theta(k)        = 0.d0
upstr_zone(k)   = 0  

! zone 1 : upper_wall
k=1
type_zone(k)    = 'Shock'
theta(k)        = phi
upstr_zone(k)   = 0
longueur(k)     = 0.5d0
profil(k)       = 1
! zone 2 : upper_wall
k=2
type_zone(k)    = 'Expansion'
theta(k)        = -phi
upstr_zone(k)   = 1
longueur(k)     = 0.5d0
profil(k)       = 1
! zone 3 : lower_wall
k=3
type_zone(k)    = 'Shock'
theta(k)        = -phi
upstr_zone(k)   = 0
longueur(k)     = 0.5d0
profil(k)       = -1
! zone 4 : lower_wall
k=4
type_zone(k)    = 'Expansion'
theta(k)        = phi
upstr_zone(k)   = 3
longueur(k)     = 0.5d0
profil(k)       = -1
! zone 5 : upper_wall wake
k=5
type_zone(k)    = 'Shock'
theta(k)        = 0.d0
upstr_zone(k)   = 2
longueur(k)     = 0.d0
profil(k)       = 0
! zone 6 : lower_wall wake
k=6
type_zone(k)    = 'Shock'
theta(k)        = 0.d0
upstr_zone(k)   = 4
longueur(k)     = 0.d0
profil(k)       = 0

if (alpha.lt.phi) then
    type_zone(1:6)    = (/'Shock       ', 'Expansion   ','Shock       ','Expansion   ','Shock       ','Shock       '/)
else if (alpha.gt.phi) then
    type_zone(1:6)    = (/'Expansion   ','Expansion   ','Shock       ','Expansion   ','Shock       ','Expansion   '/)
else   ! alpha=phi
    type_zone(1:6)    = (/'Nothing     ','Expansion   ','Shock       ','Expansion   ','Shock       ','Nothing     '/)
end if

do k=1,nz
    if ((upstr_zone(k).ne.k).and.(profil(k).ne.0)) theta(k)=theta(k)-alpha
    write(6,*)' k = ', k, ' theta = ', theta(k)
end do
theta=theta*deg2rad

!********************************************* 
write(6,110)'Q1: solution for each zone'
!********************************************* 

call solve_case
call display_outputs

!****************************************************
write(6,110)'Q2: lift and pressure drag coefficients'
!****************************************************

call aerodynamic_coefficients

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,upstr_zone)
deallocate(longueur,profil,CL_zone,CD_zone)

110  format(/,50('#'),/,3x,a,/,50('#'),/)

end subroutine diamond_airfoil



!
!******************************************************************
subroutine Exercice_11_3
!******************************************************************
!>@details shock-expansion method
implicit none

if (fichier) open(6,form='formatted',file='Solution_11-3.out')
write(6,100)'Exercice 11.3 : shock-expansion method on a lenticular airfoil'

alpha=0.d0                  
call shock_expansion_method('ShockExpansion_0.dat')

alpha=4.d0                  ! airfoil incidence
call shock_expansion_method('ShockExpansion_4.dat')
100  format(/,50('*'),/,3x,a,/,50('*'),/)

if (fichier) close(6)
end subroutine Exercice_11_3


!*******************************
subroutine shock_expansion_method(filename)
!*******************************
!>@details shock-expansion method for a given incidence
implicit none
integer,parameter                   :: n_panel=10
integer                             :: k
real(kind=8),dimension(0:n_panel)   :: x,ye,yi
real(kind=8)                        :: dx,theta0,beta,xloc,Kpe,Kpi
character(len=*),intent(in)         :: filename
integer,parameter                   :: x_option=1   ! 1 is the best choice

write(6,*) "option x_option                    : ",x_option
write(6,*) "number of panels on the upper_wall : ",n_panel

if (x_option.eq.0) then
    write(6,*) "the panel slop is the local slope of the "
    write(6,*) "first point of the  upstream panel, given with the derivative of the aifoil equation"
else
    write(6,*)"The panel slope is given with the straight line slope"
end if
!
! I use the existing subroutines,
write(6,100) alpha
chord=1.d0; thickness=0.12d0
nz=2*n_panel        ! the upstream zone is the first zone, 
                    ! then there are n_panel+1 zones on the airfoil

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),upstr_zone(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(longueur(0:nz),profil(0:nz),CL_zone(0:nz),CD_zone(0:nz))

! point distribution along the chord                         
x(0)=0.;dx=chord/real(n_panel,kind=8)
ye(0)=lenticular_airfoil(x(0),0,1)
yi(0)=ye(0)
do k=1,n_panel
    x(k)=x(k-1)+dx
    ye(k)=lenticular_airfoil(x(k),0,1)
    yi(k)=-ye(k)            !symmetrical airfoil
end do

open(1,form='formatted',file='lenticular_airfoil.dat')
write(1,*)'# lenticular airfoil, x, y_upper_wall, y_lower_wall'
write(1,'(3(f12.7,3x))')(x(k),ye(k),yi(k),k=0,n_panel)
close(1)

theta0=atan(lenticular_airfoil(0.d0,1,1))
write(6,*)'Theta in x=0             = ',theta0*rad2deg

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;upstr_zone=0;Mn=0.d0;longueur=0.d0;profil=0;CL_zone=0.d0;CD_zone=0.d0
! default initialisation:
type_zone='Expansion'
!  upper_wall / lower_wall
profil(1:n_panel)=1;
profil(n_panel+1:2*n_panel)=-1
longueur=1.d0/real(n_panel,kind=8)

! zone 0 : upstream
k=0
type_zone(k)    = 'Uniform'
M(k)            = 2.0d0
Ti(k)           = 300.d0
P(k)            = 1.d0
theta(k)        = 0.d0
upstr_zone(k)   = 0              

! zone 1 : upper_wall
k=1
type_zone(k)    = 'Shock'
theta(k)        = theta0
upstr_zone(k)   = 0
!
write(6,110)'calcul du choc oblique de tête avec le vrai angle'
call oblique_shock(upstr_zone(k),k)

!  the  slope are defined from the straight panel lines (segments)
do k=1,n_panel
    if (x_option.eq.1) then
        theta(k)=atan((ye(k)-ye(k-1))/dx)      !  segments
    else
        theta(k)=atan(lenticular_airfoil(x(k-1),1,1))  ! angle at the segment start
    end if
    upstr_zone(k)=k-1
end do

write(6,110)'Head oblique shock wave with the slope theta0'
theta(n_panel+1:2*n_panel)=-theta(1:n_panel)
theta(1:nz)=theta(1:nz)-alpha*deg2rad
! particular case, the first point is from the lower_wall
upstr_zone(n_panel+1)=0
type_zone(n_panel+1)='Shock'
do k=n_panel+2,2*n_panel
    upstr_zone(k)=k-1
end do
open(1,form='formatted',file='data_ShockExpansionMethod.dat')
do k=0,nz
    write(1,*) 'index         = ',k
    write(1,*) 'zone type     = ',type_zone(k)
    write(1,*) 'theta (°)     = ',theta(k)*rad2deg
    write(1,*) 'upstream zone = ',upstr_zone(k)
    write(1,*) 'length        = ',longueur(k)
    write(1,*) 'airfoil       = ',profil(k)
    write(1,*)
end do
close(1)

call solve_case  
call display_outputs
call aerodynamic_coefficients
beta=sqrt(M(0)**2-1.d0)
open(1,form='formatted',file=filename)
write(1,210) alpha
! for the linearized theory, the slope at calculated on the median points
do k=1,n_panel
    if (x_option.eq.1) then
        xloc=0.5d0*(x(k)+x(k-1))
    else
        xloc=x(k-1)
    end if
    Kpe=2.d0*(atan(lenticular_airfoil(xloc,1,1))-alpha*deg2rad)/beta
    Kpi=-2.d0*(atan(lenticular_airfoil(xloc,1,-1))-alpha*deg2rad)/beta
    write(1,200)xloc ,CL_zone(k),CL_zone(k+n_panel),&
                theta(k)*rad2deg,theta(k+n_panel)*rad2deg, Kpe,Kpi
end do
close(1)

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,upstr_zone)
deallocate(longueur,profil,CL_zone,CD_zone)

100 format(/,'#',50('-'),/,3x,'Incidence = ',f10.2,'°',/,'#',50('-'),/)
110 format(/,50('#'),/,3x,a,/,50('#'),/)
210 format('#  Incidence = ',f8.2,' °',/, &
    '#',5x, 'x', 13x,'Kp^+',10x,'Kp^-',7x,'theta^+ (°)',3x,'theta^- (°)', &
        5x,'Kp^+ Ack',6x,'Kp^- Ack')
200 format(7(f12.6,2x))
!*******************************
end subroutine shock_expansion_method
!*******************************


!************************************************
function lenticular_airfoil(x,opt,surface)  result(y)
!************************************************
!>@details  lenticular airfoil definition
!> or airfoil slope  pente w.r.t. x
!
! surface = 1 : upper_wall, -1 : lower_wall
! opt = 0 : airfoil camber line  y(x),   1 : derivative dy/dx

implicit none
real(kind=8)                :: y        !< y coordinate
real(kind=8),intent(in)     :: x        !< x coordinate
integer,intent(in)          :: opt      !< 0: y(x), 1:  dy(x)/dx
integer,intent(in)          :: surface  !< 1: upper_wall, -1: lower_wall
real(kind=8)                :: eta      
integer,parameter           :: airfoil_type=1 !< 1 : parabola, 2 : diamond
real(kind=8)                :: angle

eta=x/chord

select case (airfoil_type)
case (1)
    if (opt.eq.1) then
        y=2.d0*thickness/chord*(1.d0-2.d0*eta)  ! derivative
    else
        y=2.d0*thickness*eta*(1.d0-eta)         ! airfoil equation
    end if

case(2)
    angle= thickness/chord
    if (opt.eq.1) then
        if (eta.le.0.5) then
            y=angle   ! derivative
        else
            y=-angle
        end if
    else
        if (eta.le.0.5) then
            y=angle*eta         ! airfoil equation
        else
            y=angle * (1-eta)
        end if
    end if
end select

if (surface.eq.-1) y=-y ! if lower_wall

end function lenticular_airfoil
!
!******************************************************************
subroutine Exercice_11_4
!******************************************************************
!>@details Oblique shock in a channel, reflection

implicit none
real(kind=8)                        :: Pi0,P0,HsurL,M0,Ti0
integer                             :: k

if (fichier) open(6,form='formatted',file='Solution_11-4.out')
write(6,100)'Exercice 11.4 : Shock reflection in a channel '
nz=2
alpha=5.d0   ! wall slope
M0=1.5d0
Ti0=500.d0
Pi0=1.0d6
P0=Pi0*P_Pi(M0);

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),upstr_zone(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))

!********************************************* 
write(6,110)'Q1,  Q2 et Q3 : theta= 5°'
!********************************************* 
! zone 0 : upstream
k=0
type_zone(k)    = 'Uniform'
M(k)            = M0
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
upstr_zone(k)   = 0             
! zone 1 : fisrt shock
k=1
type_zone(k)    = 'Shock'
theta(k)        = alpha*deg2rad
upstr_zone(k)   = 0
! zone 2 : second shock
k=2
type_zone(k)    = 'Shock'
theta(k)        = 0.d0
upstr_zone(k)   = 1

call solve_case
call display_outputs

!********************************************* 
write(6,110)'Q4  : h/L'
!********************************************* 
alpha=alpha*deg2rad
HsurL=tan(sigma(1))*(tan(alpha)+tan(sigma(2)-alpha) )/ &
        (tan(sigma(1))+tan(sigma(2)-alpha)) 

write(6,*)'sigma_1 and sigma_2 : ', sigma(1)*rad2deg,sigma(2)*rad2deg
write(6,*)'H/L and L/H         : ',HsurL,1.d0/HsurL

!********************************************* 
write(6,110)'Q5 : theta= 10°'
!********************************************* 
alpha=10.d0  
theta(1)        = alpha*deg2rad

call solve_case

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,upstr_zone)

110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_11_4

!******************************************************************
subroutine Exercice_11_5
!******************************************************************
!>@details Interaction of 2 oblique shocks
implicit none
real(kind=8)                        :: Pi0,P0,M0,Ti0,theta1,theta2,theta_g
integer                             :: k

if (fichier) open(6,form='formatted',file='Solution_11-5.out')
write(6,100)'Exercice 11.5 :  Interaction of 2 oblique shocks'
nz=4
M0=2.d0
Ti0=500.d0
Pi0=1.0d6
P0=Pi0*P_Pi(M0);
theta1=5.d0;
theta2=15.d0;
theta_g=theta2
allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),upstr_zone(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))

! zone 0 : upstream
k=0
type_zone(k)    = 'Uniform'
M(k)            = M0
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
upstr_zone(k)   = 0              
! zone 1 : first shock
k=1
type_zone(k)    = 'Shock'
theta(k)        = theta1*deg2rad
upstr_zone(k)   = 0
! zone 2 : second shock
k=2
type_zone(k)    = 'Shock'
theta(k)        =  theta2*deg2rad
upstr_zone(k)   = 1
! zone 3 : shock-shock reflection
k=3
type_zone(k)    = 'Nothing'
theta(k)        =  theta_g*deg2rad
upstr_zone(k)   = 2
! zone 4 : 
k=4
type_zone(k)    = 'Shock'
theta(k)        = theta_g*deg2rad
upstr_zone(k)   = 0

!********************************************* 
write(6,110)'Q1, Q2 et Q3'
!********************************************* 
call solve_case
call display_outputs

!********************************************* 
write(6,110)'Q1, Q2 et Q3 '
!********************************************* 
theta_g=theta2
theta(3)        =  (theta_g-0.1d0)*deg2rad
theta(4)        =  theta(3)
type_zone(3)    = 'Shock'
call solve_case
call display_outputs

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,upstr_zone)

110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_11_5

!******************************************************************
subroutine Exercice_11_6
!******************************************************************
!>@details backward facing step (simple program)

implicit none

if (fichier) open(6,form='formatted',file='Solution_11-6.out')
write(6,100)'Exercice 11.6 : backward facing ste'

nz=2
allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))

Ti(0)=400; M(0)=2.d0; theta(0)=0.d0*deg2rad;P(0)=10000.d0;

!********************************************* 
write(6,110)'Q2 à  Q4 '
!********************************************* 
tau(1)=0.5  ! pression ratio from 0 to 1
write(*,200) 0
call solve_uniform_zone(0)
write(*,200) 1
p(1)=tau(1)*p(0);pi(1)=pi(0);Ti(1)=Ti(0)
call solve_isobarline(0,1)

!********************************************* 
write(6,110)'Q5  : oblique shock'
!********************************************* 
write(*,200) 2
theta(2)=0.d0; 
call oblique_shock(1,2)
deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn)

110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
200 format(50('*'),/,'ZONE : ', i1,/,50('*'))
if (fichier) close(6)
end subroutine Exercice_11_6


!******************************************************************
subroutine Exercice_11_61
!******************************************************************
!>@details backward facing step (2nd program)
implicit none
real(kind=8)                        :: P0,M0,Ti0,taug
integer                             :: k

if (fichier) open(6,form='formatted',file='Solution_11-61.out')
write(6,100)'Exercice 11.6 : backward facing step (2nd program)'

nz=2

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),upstr_zone(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(type_characteristic(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;upstr_zone=0;Mn=0.d0;
com='';type_characteristic='';

Ti0=400.d0; M0=2.d0; P0=10000.d0;taug=0.5d0

! zone 0 : upstream
k=0
type_zone(k)    = 'Uniform'
M(k)            = M0
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
upstr_zone(k)   = 0       
! zone 1 :  isobar line
k=1
type_zone(k)    = 'Isobarline'
P(k)            = taug*P0
type_characteristic(k)='-'
upstr_zone(k)   = 0
! zone 2 : oblique shock
k=2
type_zone(k)    = 'Shock'
theta(k)        =  0.d0
upstr_zone(k)   = 1

!********************************************* 
write(6,110)'Q2 à Q5  : '
!********************************************* 
call solve_case
call display_outputs

write(6,*)'V1/V0                          = ',V(1)/V(0)
write(6,*)'V2/V0                          = ',V(2)/V(0)
write(6,*)'recirculation length           = ',1.d0/abs(tan(theta(1)))

deallocate(M,P,Ti,theta,V,com,type_zone,upstr_zone,a,tau,T,Pi,sigma,Om,Mn)
deallocate(type_characteristic)

110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_11_61


!******************************************************************
subroutine Exercice_11_7
!******************************************************************
!>@details interaction shock - isobar line with fluid at rest
implicit none
real(kind=8)                        :: P0,M0,Ti0,tetap
integer                             :: k

if (fichier) open(6,form='formatted',file='Solution_11-7.out')
write(6,100)'Exercice 11.7 : interaction shock - isobar line with fluid at rest'

nz=2

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),upstr_zone(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(type_characteristic(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;upstr_zone=0;Mn=0.d0;
com='';type_characteristic='';

Ti0=400.d0; M0=2.d0; P0=100000.d0;tetap=10.d0

! zone 0 : upstream
k=0
type_zone(k)    = 'Uniform'
M(k)            = M0
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
upstr_zone(k)   = 0           
! zone 1 :  oblique shock
k=1
type_zone(k)    = 'Shock'
theta(k)        =  tetap*deg2rad
upstr_zone(k)   = 0
! zone 2 :  isobar line
k=2
type_zone(k)    = 'Isobarline'
P(k)            = P0
type_characteristic(k)='+'
upstr_zone(k)   = 1

!********************************************* 
write(6,110)'Q1 à Q4  : '
!********************************************* 
call solve_case
call display_outputs

write(6,*)'V1/V0                          = ',V(1)/V(0)
write(6,*)'V2/V0                          = ',V(2)/V(0)

!********************************************* 
write(6,110)'Q5  : '
!********************************************* 
M(0)            = 3.d0
theta(1)        = 20.d0*deg2rad 
call solve_case
call display_outputs

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,upstr_zone)
deallocate(type_characteristic)
110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_11_7


!******************************************************************
subroutine Exercice_11_8
!******************************************************************
!>@details interaction expansion - isobar line with fluid at rest 

implicit none
real(kind=8)                        :: P0,M0,tetap,T0
integer                             :: k

if (fichier) open(6,form='formatted',file='Solution_11-8.out')
write(6,100)'Exercice 11.8 : interaction expansion - isobar line with fluid at rest'

nz=2

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),upstr_zone(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(type_characteristic(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;upstr_zone=0;Mn=0.d0;
com='';type_characteristic='';
M0=2.5d0; T0=600.d0;  P0=100000.d0;tetap=-12.d0

! zone 0 : upstream
k=0
type_zone(k)    = 'Uniform'
M(k)            = M0
T(k)            = T0
P(k)            = P0
theta(k)        = 0.d0
upstr_zone(k)   = 0 
! zone 1 :  expansion C-
k=1
type_zone(k)    = 'Expansion'
theta(k)        =  tetap*deg2rad
type_characteristic(k)='-'
upstr_zone(k)   = 0
! zone 2 : isobar line
k=2
type_zone(k)    = 'Isobarline'
P(k)            = P0
type_characteristic(k)='+'
upstr_zone(k)   = 1

!********************************************* 
write(6,110)'Q1 à Q5  : '
!********************************************* 
call solve_case
call display_outputs

write(6,*)'V1/V0                          = ',V(1)/V(0)
write(6,*)'V2/V0                          = ',V(2)/V(0)

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,upstr_zone)
deallocate(type_characteristic)
110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_11_8

!******************************************************************
subroutine Exercice_11_9
!******************************************************************
!>@details interaction expansion - wall
implicit none
real(kind=8)                        :: P0,M0,tetap,T0
integer                             :: k

if (fichier) open(6,form='formatted',file='Solution_11-9.out')
write(6,100)'Exercice 11.9 : interaction expansion - wall'

nz=2

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),upstr_zone(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(type_characteristic(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;upstr_zone=0;Mn=0.d0;
com='';type_characteristic='';
M0=3.d0; T0=500.d0;  P0=200000.d0;tetap=-10.d0

! zone 0 : upstream
k=0
type_zone(k)    = 'Uniform'
M(k)            = M0
T(k)            = T0
P(k)            = P0
theta(k)        = 0.d0
upstr_zone(k)   = 0 
! zone 1 :  expansion C-
k=1
type_zone(k)    = 'Expansion'
theta(k)        =  tetap*deg2rad
type_characteristic(k)='-'
upstr_zone(k)   = 0
! zone 2 : expansion
k=2
type_zone(k)    = 'Expansion'
theta(k)        =  0.d0
type_characteristic(k)='+'
upstr_zone(k)   = 1

!********************************************* 
write(6,110)'Q1 à Q5  : '
!********************************************* 
call solve_case
call display_outputs

write(6,*)'V1/V0                          = ',V(1)/V(0)
write(6,*)'V2/V0                          = ',V(2)/V(0)

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,upstr_zone)
deallocate(type_characteristic)
110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_11_9

!******************************************************************
subroutine Exercice_11_10
!******************************************************************
!>@details Two upstream flow  separated by a isobar line  - oblique shock interaction 
implicit none
real(kind=8)                        :: P0,M0,tetap,Ti0,M1,theta3
integer                             :: k

if (fichier) open(6,form='formatted',file='Solution_11-10.out')
write(6,100)'Exercice 11.10 : Two upstream flow  separated by a isobar line  - oblique shock interaction'

nz=4

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),upstr_zone(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(type_characteristic(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;upstr_zone=0;Mn=0.d0;
com='';type_characteristic='';
M0=4.d0; Ti0=500.d0;  P0=10000.d0;tetap=10.d0; M1=2.d0
theta3=tetap+1.d0

! zone 0 : upstream
k=0
type_zone(k)    = 'Uniform'
M(k)            = M0
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
upstr_zone(k)   = 0  
! zone 1 : upstream
k=1
type_zone(k)    = 'Uniform'
M(k)            = M1
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
upstr_zone(k)   = 1    
! zone 2 :  oblique shock
k=2
type_zone(k)    = 'Shock'
theta(k)        =  tetap*deg2rad
upstr_zone(k)   = 0
! zone 3 : expansion
k=3
type_zone(k)    = 'Expansion'
theta(k)        =  theta3*deg2rad
type_characteristic(k)='+'
upstr_zone(k)   = 2
! zone 4 :  oblique shock
k=4
type_zone(k)    = 'Shock'
theta(k)        =  theta3*deg2rad
upstr_zone(k)   = 1


!********************************************* 
write(6,110)'Q1  : cas a '
!********************************************* 
call solve_case
call display_outputs

write(6,*) 'p4            = ',P(4)
write(6,*) 'p3            = ',P(3)
print*,'p4/p3  =',p(4)/p(3)

!********************************************* 
write(6,110)'angle of the isobar line'
!********************************************* 
call solve_isobarline_angle(theta3,p(4)/p(3),'A')

call solve_case
call display_outputs

!********************************************* 
write(6,110)'Q2  : cas b '
!********************************************* 

theta3=tetap-1.d0
! zone 3 : oblique shock
k=3
type_zone(k)    = 'Shock'
theta(k)        =  theta3*deg2rad
upstr_zone(k)   = 2
theta(4)=theta(3)
!
M(0)=M1;M(1)=M0
call solve_case
call display_outputs

write(6,*) 'p4            = ',P(4)
write(6,*) 'p3            = ',P(3)
print*,'p4/p3  =',p(4)/p(3)

!********************************************* 
write(6,110)'angle of the isobar line '
!********************************************* 
call solve_isobarline_angle(theta3,p(4)/p(3),'B')

call solve_case
call display_outputs

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,upstr_zone)
deallocate(type_characteristic)
110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)

contains

    !********************************************************
    subroutine solve_isobarline_angle(theta_in,ratio_in,cas)
    !********************************************************
!>@details  deviation angle of a isobar line with a Newton method 
    implicit none
    character(len=1),intent(in) :: cas
    real(kind=8),intent(inout)  :: theta_in
    real(kind=8),intent(in)     :: ratio_in
    real(kind=8),parameter      :: dtheta=0.1d0
    real(kind=8),parameter      :: erreur_max=1.d-2
    integer,parameter           :: iterMax=20
    integer                     :: iter
    real(kind=8)                :: gradient,variation,ratio
    ratio=ratio_in
    iter=0;
    variation=1.d0
    show=.false.                ! .true. to display details
    do while ((abs(variation).gt.erreur_max).and.(iter.lt.iterMax))
        gradient=(isobar_line_equilibrium(theta_in+dtheta,cas)-ratio)/dtheta
        variation=-(ratio-1.d0)/gradient
        theta_in=theta_in+variation
        ratio= isobar_line_equilibrium(theta_in,cas)
        iter=iter+1
        write(*,*)'ITER = ', iter,variation,theta_in
    end do

    if (iter.ne.iterMax) then
        print*,'convergence reached'
    else
        stop 'no convergence'
    end if
    show=.true. 
    end subroutine solve_isobarline_angle

    !*****************************************************
    function  isobar_line_equilibrium(thetae,cas) result(rapport)
    !*****************************************************
!>@details solve the deviation angle of the isobar line to get pressure equilibrium
!! with a Newton method
    implicit none
    real(kind=8), intent(in)    :: thetae
    character(len=1),intent(in) :: cas
    real(kind=8)                :: rapport

    theta(3)        =  thetae*deg2rad
    theta(4)        =  theta(3)

    ! zone 3
    if (cas.eq.'A') then
        Om(3)=Om(2)+theta(3)-theta(2); M(3)=invomega(Om(3))
        tau(3)=P_Pi(M(3))/P_Pi(M(2))*tau(2)
    else
        call Newton_shock_angle(sigma(3),abs(theta(3)-theta(2)),M(2))
        Mn(2)=M(2)*sin(sigma(3))
        Mn(3)=mach_downstr(Mn(2))
        M(3) = Mn(3)/sin(sigma(3)-abs(theta(3)-theta(2)))
        tau(3)=P2_P1(Mn(2))*tau(2)
    end if
    ! zone 4
    call Newton_shock_angle(sigma(4),abs(theta(4)),M(1))
    Mn(1)=M(1)*sin(sigma(4))
    Mn(4)=mach_downstr(Mn(1))
    M(4) = Mn(4)/sin(sigma(4)-abs(theta(4)))
    tau(4)=P2_P1(Mn(1))

    write(6,100),tau(3),tau(4),tau(4)/tau(3)
    100 format('#', 10x,'p3/p0 = ',f8.4,3x,'p4/p0 = ',f8.4,3x,'p4/p3 = ',f8.4)
    rapport=tau(4)/tau(3)
    end function isobar_line_equilibrium

end subroutine Exercice_11_10


!******************************************************************
subroutine Exercice_11_11
!******************************************************************
!>@detailsInteraction of shock waves in a nozzle outlet

 
real(kind=8)        :: Ti0,M0,P0,Pa
real(kind=8)        :: tetap,val
integer             :: k

if (fichier) open(6,form='formatted',file='Solution_11-11.out')
!****************************************************************************
write(6,100) "Interaction of shock waves in a nozzle outlet, exercice 11.11"
!****************************************************************************
nz=3

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),upstr_zone(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(type_characteristic(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;upstr_zone=0;Mn=0.d0;
com='';type_characteristic='';

M0=2.d0; Ti0=500.d0;  P0=10000.d0;tetap=12.5d0; 
Pa=0.2d0/P_Pi(M0)*P0;

! zone 0 : upstream
k=0
type_zone(k)    = 'Uniform'
M(k)            = M0
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
upstr_zone(k)   = 0   
! zone 1 :  oblique shock
k=1
type_zone(k)    = 'Shock'
upstr_zone(k)   = 0
! zone 2 :  oblique shock
k=2
type_zone(k)    = 'Shock'
theta(k)        = 0.d0
upstr_zone(k)   = 1

!********************************************* 
write(6,110)'Q1  : cas a '
!********************************************* 
! zone 3 : isobar line
k=3
type_zone(k)    = 'Isobarline'
P(k)            = Pa
type_characteristic(k)='+'
upstr_zone(k)   = 2

theta(1)=-search_deviation_shock_isobarline(M0,Pa/P0)

call solve_case
call display_outputs

write(6,*) 'Verification Pa/Pi0  =',Pa/Pi(0)

!********************************************* 
write(6,110)'Q2  : cas b '
!********************************************* 
! delete some data
theta(3)=0.d0
Om(1:3)=0.d0
type_characteristic=''

! zone 3 :  normal shock
k=3
type_zone(k)    = 'Normal shock'
upstr_zone(k)   = 1
! Adding to zone 2
sigma(2)        = 88.d0*deg2rad     ! initialisation to search the second solution
! Adding to zone 1
theta(1)        =  tetap*deg2rad    ! approximation of the solution
!Pi0=1.d6      ! Presure in Pa

call solve_uniform_zone(0)
call normal_shock(0,3)
val=solve_normal_shock(tetap,1)
call solve_deviation_angle(tetap,val)

!call solve_case
!call display_outputs

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,upstr_zone)
deallocate(type_characteristic)
110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
contains

    !********************************************************
    subroutine solve_deviation_angle(theta_in,ratio_in)
    !********************************************************
!>@details  deviation angle for shock interaction with a Newton method
!! 
    implicit none
    real(kind=8),intent(inout)  :: theta_in
    real(kind=8),intent(in)     :: ratio_in
    real(kind=8),parameter      :: dtheta=0.01d0
    real(kind=8),parameter      :: erreur_max=1.d-2
    integer,parameter           :: iterMax=20
    integer                     :: iter
    real(kind=8)                :: gradient,variation,ratio

    write(6,100) "Newton method"
    ratio=ratio_in
    iter=0;
    variation=1.d0
    show=.false.               
    do while ((abs(variation).gt.erreur_max).and.(iter.lt.iterMax))
        gradient=(solve_normal_shock(theta_in+dtheta,0)-ratio)/dtheta
        variation=-(ratio-1.d0)/gradient
        theta_in=theta_in+variation
        ratio= solve_normal_shock(theta_in,0)
        iter=iter+1
        write(*,120)'ITER = ', iter,variation,theta_in
    end do

    if (iter.ne.iterMax) then
        print*,'convergence reached'
    else
        stop 'no convergence'
    end if
    show=.true.
    ratio= solve_normal_shock(theta_in,1)
100  format(/,50('.'),/,3x,a,/,50('.'),/)
120  format(/,80('='),/,3x,a6,2x,i2,3x,2(f13.6,3x),/,80('='),/)
    end subroutine solve_deviation_angle


    !************************************************************
    function search_deviation_shock_isobarline(M_upstream,ratio_in) result(theta_out)
    !************************************************************
!>@details function used in the Newton method
    implicit none
    real(kind=8),intent(in)     :: ratio_in,M_upstream
    real(kind=8)                :: theta_out
    real(kind=8)                :: Mn_upstream,sigma_downstr,tmp

    write(6,100) "Deviation angle of an isobar line after a oblique shock"
    Mn_upstream=Inverse_P2_P1(ratio_in)
    write(6,*) 'Mn upstream                 = ',Mn_upstream
    sigma_downstr=asin(Mn_upstream/M_upstream)
    write(6,*) 'Shock angle            = ',sigma_downstr*rad2deg
    tmp=(2.d0+(gam-1.d0)*Mn_upstream**2)/((gam+1.d0)*Mn_upstream**2)*tan(sigma_downstr)
    theta_out=(sigma_downstr-atan(tmp))
    write(6,*) 'Deviation angle       = ',theta_out*rad2deg

100  format(/,50('.'),/,3x,a,/,50('.'),/)
    end function search_deviation_shock_isobarline

    !************************************************************
    function solve_normal_shock(deviation,opt) result(valeur)
    !************************************************************

    implicit none
    real(kind=8), intent(in)    :: deviation
    integer, intent(in)         :: opt
    real(kind=8)                :: valeur

    theta(1)=deviation*deg2rad
    call oblique_shock(0,1)
    call oblique_shock(1,2)
    valeur=tau(2)/tau(3)
    if (opt.eq.1) call display_outputs

    end function solve_normal_shock

end subroutine Exercice_11_11


!******************************************************************
subroutine Exercice_11_12
!******************************************************************
!>@details tuyère : Interaction of a centered expansion in a nozzle outle
 
real(kind=8)        :: Ti0,M0,P0,Pa
integer             :: k

if (fichier) open(6,form='formatted',file='Solution_11-12.out')
!****************************************************************************
write(6,100) "Interaction of a centered expansion in a nozzle outle, exercice 11.12"
!****************************************************************************
nz=3

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),upstr_zone(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(type_characteristic(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;upstr_zone=0;Mn=0.d0;
com='';type_characteristic='';
M0=2.d0; Ti0=500.d0;  P0=10000.d0
Pa=0.1d0/P_Pi(M0)*P0;

! zone 0 : upstream
k=0
type_zone(k)    = 'Uniform'
M(k)            = M0
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
upstr_zone(k)   = 0 
! zone 1 : upstream
k=1
type_zone(k)    = 'Isobarline'
P(k)            = Pa
type_characteristic(k)='+'
upstr_zone(k)   = 0
! zone 2 : expansion
k=2
type_zone(k)    = 'Expansion'
theta(k)        =  0.d0
type_characteristic(k)='-'
upstr_zone(k)   = 1
! zone 3 : upstream
k=3
type_zone(k)    = 'Isobarline'
P(k)            = Pa
type_characteristic(k)='+'
upstr_zone(k)   = 2

call solve_case
call display_outputs

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,upstr_zone)
deallocate(type_characteristic)
!110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
end subroutine Exercice_11_12

end module exercices_chap11
