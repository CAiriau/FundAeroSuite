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

 
!> MODULE :   Correction des exercices du chapitre 11
!! 
!! @author C. Airiau
!!

module exercices_chap11

use mod_EC
use Constante

contains


!******************************************************************
subroutine Exercice_11_1
!******************************************************************
!>@details plaque en incidence

implicit none
integer             :: k

if (fichier) open(6,form='formatted',file='Solution_11-1.out')
write(6,100)'Exercice 11.1 : Plaque en incidence'

!*********************
! Entrée des données 
!*********************
nz=4                        ! 5 zones avec l'amont
alpha=0.d0                  ! incidence de la plaque (ici on change les pentes directements)
! allocation des tableaux donnant l'état dans les différentes zones :
allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),zone_amont(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(longueur(0:nz),profil(0:nz),CL_zone(0:nz),CD_zone(0:nz))

! zone 0 : amont
k=0
type_zone(k)    = 'Uniforme'
M(k)            = 2.0d0
Ti(k)           = 300.d0
P(k)            = 1.d0
theta(k)        = 0.d0
zone_amont(k)   = 0             ! si elle même, pas de correction d'incidence
! zone 1 : extrados
k=1
type_zone(k)    = 'Detente'
theta(k)        = -alpha
zone_amont(k)   = 0
longueur(k)     = 1.d0
profil(k)       = 1
! zone 2 : intrados
k=2
type_zone(k)    = 'Choc'
theta(k)        = -alpha
zone_amont(k)   = 0
longueur(k)     = 1.d0
profil(k)       = -1
! zone 3 : extrados, sillage
k=3
type_zone(k)    = 'Choc'
theta(k)        = 0.d0
zone_amont(k)   = 1
! zone 4 : intrados,sillage
k=4
type_zone(k)    = 'Detente'
theta(k)        = 0.d0
zone_amont(k)   = 2
theta=theta*deg2rad

!********************************************* 
write(6,110)'Q2: solutions par zones'
!********************************************* 

call solve_case             ! calcul de toutes les zones 
call sorties                ! formattage des sorties


!****************************************************
write(6,110)'Q3: Coefficients de portance et de trainée'
!****************************************************

call coefficients_aerodynamiques  ! calcul des coefficients aérodynamiques

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,zone_amont)
deallocate(longueur,profil,CL_zone,CD_zone)



110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_11_1

!
!******************************************************************
subroutine Exercice_11_2
!******************************************************************
!>@details profil losangique  
implicit none

if (fichier) open(6,form='formatted',file='Solution_11-2.out')
write(6,100)'Exercice 11.2 : profil losangique'

alpha=0.d0                  ! incidence de la plaque (ici on change les pentes directements)
call profil_losange
alpha=8.d0
call profil_losange
alpha=10.d0
call profil_losange
alpha=12.d0
call profil_losange

100  format(/,50('*'),/,3x,a,/,50('*'),/)

if (fichier) close(6)
end subroutine Exercice_11_2

!*************************
subroutine profil_losange
!*************************
!>@details définition du profil en losange

integer             :: k
real(kind=8)        :: phi
!*********************
! Entrée des données 
!*********************

nz=6                        ! 7 zones avec l'amont
phi=10.d0
write(6,120)phi,alpha
120 format(/,'# PROFIL LOSANGIQUE : 1/2 angle au sommet : ',f8.2,' alpha = ',f8.2,/)
! allocation des tableaux donnant l'état dans les différentes zones :
allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),zone_amont(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(longueur(0:nz),profil(0:nz),CL_zone(0:nz),CD_zone(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;zone_amont=0;Mn=0.d0;longueur=0.d0;profil=0;CL_zone=0.d0;CD_zone=0.d0
! zone 0 : amont
k=0
type_zone(k)    = 'Uniforme'
M(k)            = 2.0d0
Ti(k)           = 300.d0
P(k)            = 1.d0
theta(k)        = 0.d0
zone_amont(k)   = 0             ! si elle même, pas de correction d'incidence

! zone 1 : extrados
k=1
type_zone(k)    = 'Choc'
theta(k)        = phi
zone_amont(k)   = 0
longueur(k)     = 0.5d0
profil(k)       = 1
! zone 2 : extrados
k=2
type_zone(k)    = 'Detente'
theta(k)        = -phi
zone_amont(k)   = 1
longueur(k)     = 0.5d0
profil(k)       = 1
! zone 3 : intrados
k=3
type_zone(k)    = 'Choc'
theta(k)        = -phi
zone_amont(k)   = 0
longueur(k)     = 0.5d0
profil(k)       = -1
! zone 4 : intrados
k=4
type_zone(k)    = 'Detente'
theta(k)        = phi
zone_amont(k)   = 3
longueur(k)     = 0.5d0
profil(k)       = -1
! zone 5 : extrados sillage
k=5
type_zone(k)    = 'Choc'
theta(k)        = 0.d0
zone_amont(k)   = 2
longueur(k)     = 0.d0
profil(k)       = 0
! zone 6 : intrados sillage
k=6
type_zone(k)    = 'Choc'
theta(k)        = 0.d0
zone_amont(k)   = 4
longueur(k)     = 0.d0
profil(k)       = 0

if (alpha.lt.phi) then
    type_zone(1:6)    = (/'Choc      ','Detente   ','Choc      ','Detente   ','Choc      ','Choc      '/)
else if (alpha.gt.phi) then
    type_zone(1:6)    = (/'Detente   ','Detente   ','Choc      ','Detente   ','Choc      ','Detente   '/)
else   ! alpha=phi
    type_zone(1:6)    = (/'Rien      ','Detente   ','Choc      ','Detente   ','Choc      ','Rien      '/)
end if

do k=1,nz
    if ((zone_amont(k).ne.k).and.(profil(k).ne.0)) theta(k)=theta(k)-alpha
   write(6,*)' k = ', k, ' theta = ', theta(k)
end do
theta=theta*deg2rad

!********************************************* 
write(6,110)'Q1: solutions par zones'
!********************************************* 

call solve_case             ! calcul de toutes les zones 
call sorties                ! formattage des sorties


!****************************************************
write(6,110)'Q2: Coefficients de portance et de trainée'
!****************************************************

call coefficients_aerodynamiques  ! calcul des coefficients aérodynamiques

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,zone_amont)
deallocate(longueur,profil,CL_zone,CD_zone)



110  format(/,50('#'),/,3x,a,/,50('#'),/)

end subroutine profil_losange



!
!******************************************************************
subroutine Exercice_11_3
!******************************************************************
!>@details methode choc detente
implicit none

if (fichier) open(6,form='formatted',file='Solution_11-3.out')
write(6,100)'Exercice 11.3 : Méthode choc-détente sur un profil lenticulaire'

alpha=0.d0                  ! incidence du profil
call methode_choc_detente('ChocDetente_0.dat')

alpha=4.d0                  ! incidence du profil
call methode_choc_detente('ChocDetente_4.dat')
100  format(/,50('*'),/,3x,a,/,50('*'),/)

if (fichier) close(6)
end subroutine Exercice_11_3


!*******************************
subroutine methode_choc_detente(filename)
!*******************************
!>@details methode pour une incidence donnée
implicit none
integer,parameter                   :: n_panel=10
integer                             :: k
real(kind=8),dimension(0:n_panel)   :: x,ye,yi
real(kind=8)                        :: dx,theta0,beta,xloc,Kpe,Kpi
character(len=*),intent(in)         :: filename
integer,parameter                   :: x_option=1   ! 1 est le meilleur choix

write(6,*) "option x_option                   : ",x_option
write(6,*) "nombre de panneaux sur l'extrados : ",n_panel

if (x_option.eq.0) then
    write(6,*) "la pente du panneau est calculée à partir de la pente locale"
    write(6,*) "du premier point amont du panneau avec la dérivée de l'équation du profil"
else
    write(6,*)"la pente du panneau est calculée à partir du segment "
end if
!

! j'utilise les subroutines déjà programmées, cela peut paraitre plus 
! complexe, mais c'est plus clair d'un point de vue méthodologie


write(6,100) alpha

chord=1.d0; thickness=0.12d0
nz=2*n_panel          ! la première zone est la zone amont, 
                      ! puis il y a n_panel+1 zones sur le profil

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),zone_amont(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(longueur(0:nz),profil(0:nz),CL_zone(0:nz),CD_zone(0:nz))

! répartition des points le long de la corde                        
x(0)=0.;dx=chord/real(n_panel,kind=8)
ye(0)=profil_lenticulaire(x(0),0,1)
yi(0)=ye(0)
do k=1,n_panel
    x(k)=x(k-1)+dx
    ye(k)=profil_lenticulaire(x(k),0,1)
    yi(k)=-ye(k)            ! car le profil est symétrique
end do

open(1,form='formatted',file='profil_lenticulaire.dat')
write(1,*)'# profil lenticulaire, x, y_extrados, y_intrados'
write(1,'(3(f12.7,3x))')(x(k),ye(k),yi(k),k=0,n_panel)
close(1)

theta0=atan(profil_lenticulaire(0.d0,1,1))
write(6,*)'Theta en x=0             = ',theta0*rad2deg


sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;zone_amont=0;Mn=0.d0;longueur=0.d0;profil=0;CL_zone=0.d0;CD_zone=0.d0
! initialisation par défaut:
type_zone='Detente'
! délimitation extrados/intrados
profil(1:n_panel)=1;
profil(n_panel+1:2*n_panel)=-1
longueur=1.d0/real(n_panel,kind=8)


! zone 0 : amont
k=0
type_zone(k)    = 'Uniforme'
M(k)            = 2.0d0
Ti(k)           = 300.d0
P(k)            = 1.d0
theta(k)        = 0.d0
zone_amont(k)   = 0             ! si elle même, pas de correction d'incidence


! zone 1 : extrados
k=1
type_zone(k)    = 'Choc'
theta(k)        = theta0
zone_amont(k)   = 0
!
write(6,110)'calcul du choc oblique de tête avec le vrai angle'
call choc_oblique(zone_amont(k),k)

! maintenant je fais tous les calculs à partie des angles définis par les segments.
do k=1,n_panel
    if (x_option.eq.1) then
        theta(k)=atan((ye(k)-ye(k-1))/dx)      ! calculé sur les segments
    else
        theta(k)=atan(profil_lenticulaire(x(k-1),1,1))  ! calculé sur l'angle de pente du début du segment .
    end if
    zone_amont(k)=k-1
end do

write(6,110)'calcul du choc oblique de tête pour une pente de theta0 recalculé éventuellement'
theta(n_panel+1:2*n_panel)=-theta(1:n_panel)
theta(1:nz)=theta(1:nz)-alpha*deg2rad
! cas particulier du premier point sur l'intrados
zone_amont(n_panel+1)=0
type_zone(n_panel+1)='Choc'
do k=n_panel+2,2*n_panel
    zone_amont(k)=k-1
end do
open(1,form='formatted',file='data_MethodeChocDetente.dat')
do k=0,nz
    write(1,*) 'indice     = ',k
    write(1,*) 'type       = ',type_zone(k)
    write(1,*) 'theta (°)  = ',theta(k)*rad2deg
    write(1,*) 'zone amont = ',zone_amont(k)
    write(1,*) 'longueur   = ',longueur(k)
    write(1,*) 'profil     = ',profil(k)
    write(1,*)
end do
close(1)

! les données sont entrées, on calcule les zones
call solve_case             ! calcul de toutes les zones 
call sorties                ! formattage des sorties

call coefficients_aerodynamiques  ! calcul des coefficients aérodynamiques
beta=sqrt(M(0)**2-1.d0)
open(1,form='formatted',file=filename)
write(1,210) alpha
! pour la théorie linéarisée les pentes aux points médiants doivent être recalculés.
do k=1,n_panel
    if (x_option.eq.1) then
        xloc=0.5d0*(x(k)+x(k-1))
    else
        xloc=x(k-1)
    end if
    Kpe=2.d0*(atan(profil_lenticulaire(xloc,1,1))-alpha*deg2rad)/beta
    Kpi=-2.d0*(atan(profil_lenticulaire(xloc,1,-1))-alpha*deg2rad)/beta
    write(1,200)xloc ,CL_zone(k),CL_zone(k+n_panel),&
               theta(k)*rad2deg,theta(k+n_panel)*rad2deg, Kpe,Kpi
end do
close(1)


deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,zone_amont)
deallocate(longueur,profil,CL_zone,CD_zone)


100 format(/,'#',50('-'),/,3x,'Incidence = ',f10.2,'°',/,'#',50('-'),/)
110 format(/,50('#'),/,3x,a,/,50('#'),/)
210 format('#  incidence = ',f8.2,' °',/, &
    '#',5x, 'x', 13x,'Kp_e',10x,'Kp_i',7x,'theta_e (°)',3x,'theta_i (°)', &
        5x,'Kp_e Ack',6x,'Kp_i Ack')
200 format(7(f12.6,2x))
!*******************************
end subroutine methode_choc_detente
!*******************************



!************************************************
function profil_lenticulaire(x,opt,surface)  result(y)
!************************************************
!>@details  définition du profil lenticulaire
! calcul de la pente du profil en fonction de x
!
! surface = 1 : extrados, -1 : intrados
! opt = 0 : profil y(x),   1 : dérivée profil dy/dx

implicit none

real(kind=8)                :: y
real(kind=8),intent(in)     :: x
integer,intent(in)          :: opt      ! 0: y(x), 1:  dy(x)/dx
integer,intent(in)          :: surface  ! 1: extrados, -1: intrados
real(kind=8)                :: eta
integer,parameter           :: type_profil=1 ! 1 : parabolique, 2 : losangique
real(kind=8)                :: angle


eta=x/chord

select case (type_profil)
case (1)
    if (opt.eq.1) then
         y=2.d0*thickness/chord*(1.d0-2.d0*eta)  ! dérivée
    else
         y=2.d0*thickness*eta*(1.d0-eta)         ! equation du profil
    end if

case(2)
    angle= thickness/chord
    if (opt.eq.1) then
         if (eta.le.0.5) then
            y=angle  ! dérivée
         else
            y=-angle
         end if
    else
         if (eta.le.0.5) then
            y=angle*eta         ! equation du profil
         else
            y=angle * (1-eta)
         end if
    end if
end select


if (surface.eq.-1) y=-y                     ! si intrados

end function profil_lenticulaire
!
!******************************************************************
subroutine Exercice_11_4
!******************************************************************
!>@details Choc dans un canal

implicit none
real(kind=8)                        :: Pi0,P0,HsurL,M0,Ti0
integer                             :: k

if (fichier) open(6,form='formatted',file='Solution_11-4.out')
write(6,100)'Exercice 11.4 : Réflexion de choc dans un canal'
nz=2
alpha=5.d0                  ! pente de la paroi
M0=1.5d0
Ti0=500.d0
Pi0=1.0d6
P0=Pi0*P_Pi(M0);
! allocation des tableaux donnant l'état dans les différentes zones :
allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),zone_amont(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))



!********************************************* 
write(6,110)'Q1,  Q2 et Q3 : theta= 5°'
!********************************************* 
! zone 0 : amont
k=0
type_zone(k)    = 'Uniforme'
M(k)            = M0
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
zone_amont(k)   = 0             ! si elle même, pas de correction d'incidence
! zone 1 : premier choc
k=1
type_zone(k)    = 'Choc'
theta(k)        = alpha*deg2rad
zone_amont(k)   = 0
! zone 2 : second choc
k=2
type_zone(k)    = 'Choc'
theta(k)        = 0.d0
zone_amont(k)   = 1

call solve_case
call sorties

!********************************************* 
write(6,110)'Q4  : h/L'
!********************************************* 
alpha=alpha*deg2rad
HsurL=tan(sigma(1))*(tan(alpha)+tan(sigma(2)-alpha) )/ &
      (tan(sigma(1))+tan(sigma(2)-alpha)) 

write(6,*)'sigma_1 et sigma_2', sigma(1)*rad2deg,sigma(2)*rad2deg

write(6,*)'H/L et L/H              = ',HsurL,1.d0/HsurL

!********************************************* 
write(6,110)'Q5 : theta= 10°'
!********************************************* 

alpha=10.d0                  ! pente de la paroi
theta(1)        = alpha*deg2rad

call solve_case


deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,zone_amont)

110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_11_4

!******************************************************************
subroutine Exercice_11_5
!******************************************************************
!>@details Interaction de 2 chocs obliques
implicit none
real(kind=8)                        :: Pi0,P0,M0,Ti0,theta1,theta2,theta_g
integer                             :: k

if (fichier) open(6,form='formatted',file='Solution_11-5.out')
write(6,100)'Exercice 11.5 : Interaction de 2 chocs obliques'
nz=4
M0=2.d0
Ti0=500.d0
Pi0=1.0d6
P0=Pi0*P_Pi(M0);
theta1=5.d0;
theta2=15.d0;
theta_g=theta2
! allocation des tableaux donnant l'état dans les différentes zones :
allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),zone_amont(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))

! définitions des zones
! zone 0 : amont
k=0
type_zone(k)    = 'Uniforme'
M(k)            = M0
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
zone_amont(k)   = 0             ! si elle même, pas de correction d'incidence
! zone 1 : premier choc
k=1
type_zone(k)    = 'Choc'
theta(k)        = theta1*deg2rad
zone_amont(k)   = 0
! zone 2 : second choc
k=2
type_zone(k)    = 'Choc'
theta(k)        =  theta2*deg2rad
zone_amont(k)   = 1
! zone 3 : réflexion du choc choc
k=3
type_zone(k)    = 'Rien'
theta(k)        =  theta_g*deg2rad
zone_amont(k)   = 2
! zone 4 : 
k=4
type_zone(k)    = 'Choc'
theta(k)        = theta_g*deg2rad
zone_amont(k)   = 0


!********************************************* 
write(6,110)'Q1, Q2 et Q3 simultanément '
!********************************************* 
call solve_case
call sorties

!********************************************* 
write(6,110)'Q1, Q2 et Q3 simultanément '
!********************************************* 
theta_g=theta2
theta(3)        =  (theta_g-0.1d0)*deg2rad
theta(4)        =  theta(3)
type_zone(3)    = 'Choc'
call solve_case
call sorties


deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,zone_amont)

110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_11_5

!******************************************************************
subroutine Exercice_11_6
!******************************************************************
!>@details Marche descendante (programmation basique)

implicit none

if (fichier) open(6,form='formatted',file='Solution_11-6.out')
write(6,100)'Exercice 11.6 : Marche descendante'

nz=2
allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))

Ti(0)=400; M(0)=2.d0; theta(0)=0.d0*deg2rad;P(0)=10000.d0;

!********************************************* 
write(6,110)'Q2 à  Q4 '
!********************************************* 
tau(1)=0.5  ! rapport de pression entre 0 et 1
write(*,200) 0
call calcul_zone_uniforme(0)
write(*,200) 1
p(1)=tau(1)*p(0);pi(1)=pi(0);Ti(1)=Ti(0)
call calcul_ligne_glissement(0,1)

!********************************************* 
write(6,110)'Q5  : choc oblique'
!********************************************* 
write(*,200) 2
theta(2)=0.d0; 
call choc_oblique(1,2)
deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn)

110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
200 format(50('*'),/,'ZONE : ', i1,/,50('*'))
if (fichier) close(6)
end subroutine Exercice_11_6


!******************************************************************
subroutine Exercice_11_6bis
!******************************************************************
!>@details Marche descendante (2nd programmation)
implicit none
real(kind=8)                        :: P0,M0,Ti0,taug
integer                             :: k

if (fichier) open(6,form='formatted',file='Solution_11-6bis.out')
write(6,100)'Exercice 11.6 : Marche descendante, seconde programmation'

nz=2

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),zone_amont(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(type_caracteristique(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;zone_amont=0;Mn=0.d0;
com='';type_caracteristique='';

Ti0=400.d0; M0=2.d0; P0=10000.d0;taug=0.5d0

! définitions des zones
! zone 0 : amont
k=0
type_zone(k)    = 'Uniforme'
M(k)            = M0
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
zone_amont(k)   = 0             ! si elle même, pas de correction d'incidence
! zone 1 :  ligne de glissement
k=1
type_zone(k)    = 'Glissement'
P(k)            = taug*P0
type_caracteristique(k)='-'
zone_amont(k)   = 0
! définitions des zones
! zone 2 : choc oblique
k=2
type_zone(k)    = 'Choc'
theta(k)        =  0.d0
zone_amont(k)   = 1

!********************************************* 
write(6,110)'Q2 à Q5  : '
!********************************************* 
call solve_case
call sorties

write(6,*)'V1/V0                          = ',V(1)/V(0)
write(6,*)'V2/V0                          = ',V(2)/V(0)
write(6,*)'Longueur de la recirculation   = ',1.d0/abs(tan(theta(1)))


deallocate(M,P,Ti,theta,V,com,type_zone,zone_amont,a,tau,T,Pi,sigma,Om,Mn)
deallocate(type_caracteristique)

110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_11_6bis


!******************************************************************
subroutine Exercice_11_7
!******************************************************************
!>@details interaction choc - ligne isobare avec fluide au repos
implicit none
real(kind=8)                        :: P0,M0,Ti0,tetap
integer                             :: k

if (fichier) open(6,form='formatted',file='Solution_11-7.out')
write(6,100)'Exercice 11.7 : interaction choc oblique / ligne isobare avec fluide au repos'

nz=2

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),zone_amont(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(type_caracteristique(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;zone_amont=0;Mn=0.d0;
com='';type_caracteristique='';

Ti0=400.d0; M0=2.d0; P0=100000.d0;tetap=10.d0

! définitions des zones
! zone 0 : amont
k=0
type_zone(k)    = 'Uniforme'
M(k)            = M0
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
zone_amont(k)   = 0             ! si elle même, pas de correction d'incidence
! zone 1 :  choc oblique
k=1
type_zone(k)    = 'Choc'
theta(k)        =  tetap*deg2rad
zone_amont(k)   = 0
! zone 2 : choc oblique
k=2
type_zone(k)    = 'Glissement'
P(k)            = P0
type_caracteristique(k)='+'
zone_amont(k)   = 1

!********************************************* 
write(6,110)'Q1 à Q4  : '
!********************************************* 
call solve_case
call sorties

write(6,*)'V1/V0                          = ',V(1)/V(0)
write(6,*)'V2/V0                          = ',V(2)/V(0)

!********************************************* 
write(6,110)'Q5  : '
!********************************************* 
M(0)            = 3.d0
theta(1)        = 20.d0*deg2rad 
call solve_case
call sorties


deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,zone_amont)
deallocate(type_caracteristique)
110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_11_7


!******************************************************************
subroutine Exercice_11_8
!******************************************************************
!>@details interaction détente - ligne isobare avec fluide au repos

implicit none
real(kind=8)                        :: P0,M0,tetap,T0
integer                             :: k

if (fichier) open(6,form='formatted',file='Solution_11-8.out')
write(6,100)'Exercice 11.8 : interaction détente / ligne isobare avec fluide au repos'

nz=2

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),zone_amont(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(type_caracteristique(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;zone_amont=0;Mn=0.d0;
com='';type_caracteristique='';
M0=2.5d0; T0=600.d0;  P0=100000.d0;tetap=-12.d0

! définitions des zones
! zone 0 : amont
k=0
type_zone(k)    = 'Uniforme'
M(k)            = M0
T(k)            = T0
P(k)            = P0
theta(k)        = 0.d0
zone_amont(k)   = 0             ! si elle même, pas de correction d'incidence
! zone 1 :  detente C-
k=1
type_zone(k)    = 'Detente'
theta(k)        =  tetap*deg2rad
type_caracteristique(k)='-'
zone_amont(k)   = 0
! zone 2 : ligne isobare
k=2
type_zone(k)    = 'Glissement'
P(k)            = P0
type_caracteristique(k)='+'
zone_amont(k)   = 1

!********************************************* 
write(6,110)'Q1 à Q5  : '
!********************************************* 
call solve_case
call sorties

write(6,*)'V1/V0                          = ',V(1)/V(0)
write(6,*)'V2/V0                          = ',V(2)/V(0)

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,zone_amont)
deallocate(type_caracteristique)
110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_11_8

!******************************************************************
subroutine Exercice_11_9
!******************************************************************
!>@details interaction détente - paroi
implicit none
real(kind=8)                        :: P0,M0,tetap,T0
integer                             :: k

if (fichier) open(6,form='formatted',file='Solution_11-9.out')
write(6,100)'Exercice 11.9 : interaction détente / paroi'

nz=2

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),zone_amont(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(type_caracteristique(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;zone_amont=0;Mn=0.d0;
com='';type_caracteristique='';
M0=3.d0; T0=500.d0;  P0=200000.d0;tetap=-10.d0

! définitions des zones
! zone 0 : amont
k=0
type_zone(k)    = 'Uniforme'
M(k)            = M0
T(k)            = T0
P(k)            = P0
theta(k)        = 0.d0
zone_amont(k)   = 0             ! si elle même, pas de correction d'incidence
! zone 1 :  detente C-
k=1
type_zone(k)    = 'Detente'
theta(k)        =  tetap*deg2rad
type_caracteristique(k)='-'
zone_amont(k)   = 0
! zone 2 : ligne isobare
k=2
type_zone(k)    = 'Detente'
theta(k)        =  0.d0
type_caracteristique(k)='+'
zone_amont(k)   = 1

!********************************************* 
write(6,110)'Q1 à Q5  : '
!********************************************* 
call solve_case
call sorties

write(6,*)'V1/V0                          = ',V(1)/V(0)
write(6,*)'V2/V0                          = ',V(2)/V(0)

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,zone_amont)
deallocate(type_caracteristique)
110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)
end subroutine Exercice_11_9

!******************************************************************
subroutine Exercice_11_10
!******************************************************************
!>@details deux écoulements amont + Ligne de glissement + choc oblique
implicit none
real(kind=8)                        :: P0,M0,tetap,Ti0,M1,theta3
integer                             :: k

if (fichier) open(6,form='formatted',file='Solution_11-10.out')
write(6,100)'Exercice 11.10 : ligne de glissment amont + choc oblique'

nz=4

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),zone_amont(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(type_caracteristique(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;zone_amont=0;Mn=0.d0;
com='';type_caracteristique='';

M0=4.d0; Ti0=500.d0;  P0=10000.d0;tetap=10.d0; M1=2.d0

theta3=tetap+1.d0


! définitions des zones
! zone 0 : amont
k=0
type_zone(k)    = 'Uniforme'
M(k)            = M0
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
zone_amont(k)   = 0             ! si elle même, pas de correction d'incidence
! zone 1 : amont
k=1
type_zone(k)    = 'Uniforme'
M(k)            = M1
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
zone_amont(k)   = 1             ! si elle même, pas de correction d'incidence
! zone 2 :  choc oblique
k=2
type_zone(k)    = 'Choc'
theta(k)        =  tetap*deg2rad
zone_amont(k)   = 0
! zone 3 : détente
k=3
type_zone(k)    = 'Detente'
theta(k)        =  theta3*deg2rad
type_caracteristique(k)='+'
zone_amont(k)   = 2
! zone 4 :  choc oblique
k=4
type_zone(k)    = 'Choc'
theta(k)        =  theta3*deg2rad
zone_amont(k)   = 1


!********************************************* 
write(6,110)'Q1  : cas a '
!********************************************* 
call solve_case
call sorties

write(6,*) 'p4            = ',P(4)
write(6,*) 'p3            = ',P(3)
print*,'p4/p3  =',p(4)/p(3)

!********************************************* 
write(6,110)'Recherche angle ligne de glissement '
!********************************************* 
call recherche_angle_glissement(theta3,p(4)/p(3),'A')

call solve_case
call sorties

!********************************************* 
write(6,110)'Q2  : cas b '
!********************************************* 

theta3=tetap-1.d0
! zone 3 : détente
k=3
type_zone(k)    = 'Choc'
theta(k)        =  theta3*deg2rad
zone_amont(k)   = 2
theta(4)=theta(3)
!
M(0)=M1;M(1)=M0
call solve_case
call sorties

write(6,*) 'p4            = ',P(4)
write(6,*) 'p3            = ',P(3)
print*,'p4/p3  =',p(4)/p(3)

!********************************************* 
write(6,110)'Recherche angle ligne de glissement '
!********************************************* 
call recherche_angle_glissement(theta3,p(4)/p(3),'B')

call solve_case
call sorties




deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,zone_amont)
deallocate(type_caracteristique)
110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
if (fichier) close(6)

contains

    !********************************************************
    subroutine recherche_angle_glissement(theta_in,ratio_in,cas)
    !********************************************************
!>@details calcul de l'angle de déviation d'une ligne de glissement par 
!!  une méthode de Newton
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
    show=.false.                ! pour ne pas afficher le détail des calculs 
    do while ((abs(variation).gt.erreur_max).and.(iter.lt.iterMax))
        gradient=(equilibre_glissement(theta_in+dtheta,cas)-ratio)/dtheta
        variation=-(ratio-1.d0)/gradient
        theta_in=theta_in+variation
        ratio= equilibre_glissement(theta_in,cas)
        iter=iter+1
        write(*,*)'ITER = ', iter,variation,theta_in
    end do

    if (iter.ne.iterMax) then
        print*,'convergence'
    else
        stop 'pas de convergence'
    end if
    show=.true.             ! pour revenir à l'option d'avant la subroutine
    end subroutine recherche_angle_glissement

    !*****************************************************
    function  equilibre_glissement(thetae,cas) result(rapport)
    !*****************************************************
!>@details calcul de la déviation pour obtenir l'équilibre d'une ligne de glissement
!! fonction utilisée dans la méthode de Newton
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
        call Newton_Angle_Choc(sigma(3),abs(theta(3)-theta(2)),M(2))
        Mn(2)=M(2)*sin(sigma(3))
        Mn(3)=Mach_Aval(Mn(2))
        M(3) = Mn(3)/sin(sigma(3)-abs(theta(3)-theta(2)))
        tau(3)=P2_P1(Mn(2))*tau(2)
    end if
    ! zone 4
    call Newton_Angle_Choc(sigma(4),abs(theta(4)),M(1))
    Mn(1)=M(1)*sin(sigma(4))
    Mn(4)=Mach_Aval(Mn(1))
    M(4) = Mn(4)/sin(sigma(4)-abs(theta(4)))
    tau(4)=P2_P1(Mn(1))

    write(6,100),tau(3),tau(4),tau(4)/tau(3)
    100 format('#', 10x,'p3/p0 = ',f8.4,3x,'p4/p0 = ',f8.4,3x,'p4/p3 = ',f8.4)
    rapport=tau(4)/tau(3)
    end function equilibre_glissement

end subroutine Exercice_11_10


!******************************************************************
subroutine Exercice_11_11
!******************************************************************
!>@details tuyère : interaction d'ondes de chocs en sortie

 
real(kind=8)        :: Ti0,M0,P0,Pa
real(kind=8)        :: tetap,val
integer             :: k

if (fichier) open(6,form='formatted',file='Solution_11-11.out')
!****************************************************************************
write(6,100) "Interactions de chocs à la sortie d'une tuyère, exercice 11.11"
!****************************************************************************
nz=3

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),zone_amont(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(type_caracteristique(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;zone_amont=0;Mn=0.d0;
com='';type_caracteristique='';

M0=2.d0; Ti0=500.d0;  P0=10000.d0;tetap=12.5d0; 
Pa=0.2d0/P_Pi(M0)*P0;

! définitions des zones
! zone 0 : amont
k=0
type_zone(k)    = 'Uniforme'
M(k)            = M0
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
zone_amont(k)   = 0             ! si elle même, pas de correction d'incidence
! zone 1 :  choc oblique
k=1
type_zone(k)    = 'Choc'
zone_amont(k)   = 0
! zone 2 :  choc oblique
k=2
type_zone(k)    = 'Choc'
theta(k)        = 0.d0
zone_amont(k)   = 1


!********************************************* 
write(6,110)'Q1  : cas a '
!********************************************* 
! zone 3 : ligne isobare
k=3
type_zone(k)    = 'Glissement'
P(k)            = Pa
type_caracteristique(k)='+'
zone_amont(k)   = 2

theta(1)=-cherche_deviation_choc_glissement(M0,Pa/P0)

call solve_case
call sorties

write(6,*) 'Verification Pa/Pi0  =',Pa/Pi(0)

!********************************************* 
write(6,110)'Q2  : cas b '
!********************************************* 
! effacement de quelques données
theta(3)=0.d0
Om(1:3)=0.d0


type_caracteristique=''

! zone 3 :  choc droit
k=3
type_zone(k)    = 'Choc droit'
zone_amont(k)   = 1
! Additif zone 2
sigma(2)        = 88.d0*deg2rad     ! initialisation pour chercher la seconde solution 

! Additif zone 1
theta(1)        =  tetap*deg2rad    ! approximation de la solution
!Pi0=1.d6      ! Pression en bar

call calcul_zone_uniforme(0)
call choc_droit(0,3)
!  do while (tetap.gt.0)
!      val=solution_choc_droit(tetap,1)
!      print*, 'valeur finale = ',val
!      print*,'entrer theta en degrés (négatif pour sortir)'
!      read*,tetap
!  end do
!           
!      print*,'entrer theta en degrés pour le test final'
!      read*,tetap

val=solution_choc_droit(tetap,1)
call recherche_angle_deviation(tetap,val)

!call solve_case
!call sorties

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,zone_amont)
deallocate(type_caracteristique)
110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
contains

    !********************************************************
    subroutine recherche_angle_deviation(theta_in,ratio_in)
    !********************************************************
!>@details  calcul de l'angle de déviation pour l'interaction de chocs
!! par une méthode de Newton
    implicit none
    real(kind=8),intent(inout)  :: theta_in
    real(kind=8),intent(in)     :: ratio_in
    real(kind=8),parameter      :: dtheta=0.01d0
    real(kind=8),parameter      :: erreur_max=1.d-2
    integer,parameter           :: iterMax=20
    integer                     :: iter
    real(kind=8)                :: gradient,variation,ratio

    write(6,100) "Méthode de Newton"
    ratio=ratio_in
    iter=0;
    variation=1.d0
    show=.false.                ! pour ne pas afficher le détail des calculs 
    do while ((abs(variation).gt.erreur_max).and.(iter.lt.iterMax))
        gradient=(solution_choc_droit(theta_in+dtheta,0)-ratio)/dtheta
        variation=-(ratio-1.d0)/gradient
        theta_in=theta_in+variation
        ratio= solution_choc_droit(theta_in,0)
        iter=iter+1
        write(*,120)'ITER = ', iter,variation,theta_in
    end do

    if (iter.ne.iterMax) then
        print*,'convergence'
    else
        stop 'pas de convergence'
    end if
    show=.true.             ! pour revenir à l'option d'avant la subroutine
    ratio= solution_choc_droit(theta_in,1)
100  format(/,50('.'),/,3x,a,/,50('.'),/)
120  format(/,80('='),/,3x,a6,2x,i2,3x,2(f13.6,3x),/,80('='),/)
    end subroutine recherche_angle_deviation


    !************************************************************
    function cherche_deviation_choc_glissement(M_amont,ratio_in) result(theta_out)
    !************************************************************
!>@details fonction utilisée par la méthode de Newton
    implicit none
    real(kind=8),intent(in)     :: ratio_in,M_amont
    real(kind=8)                :: theta_out
    real(kind=8)                :: Mn_amont,sigma_aval,tmp

    write(6,100) "Angle de déviation d'une ligne de glissement derrière un choc oblique"
    Mn_amont=Inverse_P2_P1(ratio_in)
    write(6,*) 'Mn amont                 = ',Mn_amont
    sigma_aval=asin(Mn_amont/M_amont)
    write(6,*) 'Angle de choc            = ',sigma_aval*rad2deg
    tmp=(2.d0+(gam-1.d0)*Mn_amont**2)/((gam+1.d0)*Mn_amont**2)*tan(sigma_aval)
    theta_out=(sigma_aval-atan(tmp))
    write(6,*) 'Angle de deviation       = ',theta_out*rad2deg

100  format(/,50('.'),/,3x,a,/,50('.'),/)
    end function cherche_deviation_choc_glissement

    !************************************************************
    function solution_choc_droit(deviation,opt) result(valeur)
    !************************************************************

    implicit none
    real(kind=8), intent(in)    :: deviation
    integer, intent(in)         :: opt
    real(kind=8)                :: valeur

    theta(1)=deviation*deg2rad
    call choc_oblique(0,1)
    call choc_oblique(1,2)
    valeur=tau(2)/tau(3)
    if (opt.eq.1) call sorties

    end function solution_choc_droit

end subroutine Exercice_11_11


!******************************************************************
subroutine Exercice_11_12
!******************************************************************
!>@details tuyère : interaction d'un faisceau de détente en sortie
 
real(kind=8)        :: Ti0,M0,P0,Pa
integer             :: k

if (fichier) open(6,form='formatted',file='Solution_11-12.out')
!****************************************************************************
write(6,100) "Interactions d'un faisceau de détente à la sortie d'une tuyère, exercice 11.12"
!****************************************************************************
nz=3

allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),zone_amont(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
allocate(type_caracteristique(0:nz))

sigma=0.d0; theta=0.d0;Om=0.d0;M=0.d0;P=0.d0;Ti=0.d0;V=0.d0
a=0.d0;tau=0.d0;T=0.d0;Pi=0.d0;zone_amont=0;Mn=0.d0;
com='';type_caracteristique='';

M0=2.d0; Ti0=500.d0;  P0=10000.d0
Pa=0.1d0/P_Pi(M0)*P0;

! définitions des zones
! zone 0 : amont
k=0
type_zone(k)    = 'Uniforme'
M(k)            = M0
Ti(k)           = Ti0
P(k)            = P0
theta(k)        = 0.d0
zone_amont(k)   = 0             ! si elle même, pas de correction d'incidence
! zone 1 : amont
k=1
type_zone(k)    = 'Glissement'
P(k)            = Pa
type_caracteristique(k)='+'
zone_amont(k)   = 0
! zone 2 : détente
k=2
type_zone(k)    = 'Detente'
theta(k)        =  0.d0
type_caracteristique(k)='-'
zone_amont(k)   = 1
! zone 3 : amont
k=3
type_zone(k)    = 'Glissement'
P(k)            = Pa
type_caracteristique(k)='+'
zone_amont(k)   = 2


call solve_case
call sorties

deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,zone_amont)
deallocate(type_caracteristique)
!110  format(/,50('#'),/,3x,a,/,50('#'),/)
100  format(/,50('*'),/,3x,a,/,50('*'),/)
end subroutine Exercice_11_12

end module exercices_chap11
