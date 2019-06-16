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

 program busemann
! Choc oblique
! Christophe Airiau, mars 2011
! trace la polaire de busemann ainsi
!  que la courbe de choc sigma fonction de theta
!  aide : 2 parametres à entrer : M1 et i1
!         - si i1 = 0, calcul de la polaire pour M1 = x1
!         - si i1 > 0 : calcul de la polaire pour les valeurs du tableau entre
!              int(M1) et int(M1)+i1-1  soit i1 valeurs.
!  par exemple: 1.01 , 30 va donner les courbes pour tous les nombres de Mach du tableau.
!  
!   en sortie, 3 types de fichier  :
!           - courbes_Mxxx.dat  : contient des sorties pour chaque Mach
!           - mach_unitaire.dat : la ligne où M_2 = 1 (ligne sonique)
!           - theta_max.dat     : contient les nombres de Mach où on a la déviation maximale
!
!! UTILISATION :
!!  ./chocs_obliques
!!  2 paramètres à entrer :
!!      - si M1,0      : calcul du choc pour la valeur de M1
!!      _ si i_debut,n : calcul des chocs pour les valeurs de M1 du tableau M1(i_debut,i_debut+n-1)
!!
implicit none

real (kind=8)                       :: M1,sigma,sigma_min,sigma_max,dsigma,value_min
real (kind=8)                       :: us,vs,Mn1,Mn2,M2,theta,M2min,value,theta2,sigma2
real (kind=8)                       :: theta_max
integer, parameter                  :: n=1000,npt_max=40
real (kind=8), dimension(npt_max)   :: Mach1
real (kind=8),parameter             :: pi = 4.d0*atan(1.d0), coef = 180.d0/pi

integer i,npt,k,ind1,option,nmin
character(len=50)                  :: filename

 
Mach1(1:15)=(/1.01,1.03,1.1, 1.15, 1.20, 1.25,  1.3,  1.35, 1.4,  1.45, 1.5, 1.6,  1.7, 1.8, 2. /)
Mach1(16:30)=(/2.2,2.4,2.6,2.8,3.,3.2,3.4,3.6,3.8,4.,4.5,5.,6.,10.,100./)
npt=30

write(*,180) Mach1(1:npt)
write(*,*) 'Nombre de courbes de chocs possibles : ',npt

  180 format('Mach numbers M1 : ',/,3(10(F8.2,2x),/))

call menus
!  open(3,form='formatted',file='theta_max.dat',status='old', access='append')
!  open(2,form='formatted',file='m2.dat',status='old', access='append')
open(4,form='formatted',file='mach_unitaire.dat' )
write(4,170)
  170  format ('# theta (deg)',3x,'sigma (deg)',8x,'M1',9x,'M2 unitaire')
open(3,form='formatted',file='theta_max.dat' )
write(3,160)
  160  format ('# theta (deg)',3x,'sigma (deg)',8x,'M1',12x ,'M2')
open(2,form='formatted',file='M2.dat' )

do k=ind1,npt
      M1=Mach1(k)
      sigma_min = asin (1./M1)
      sigma_max = 0.5*pi

      dsigma=(sigma_max-sigma_min)/float(n-1)

      filename='courbes_M'//carac(3,int(100*M1))//'.dat'
      open(1,form='formatted',file=trim(filename))
      write(1,100)M1
      100 format(/,'# Polaire de Busemann pour M = ',f11.4,/,'#',4x,'theta',10x,'sigma', 10x,'u2/V1', &
            10x,'v2/V1',10x,'Mn2', 10x, 'M2', 10X, 'Mn1', /,'#', /)
      sigma=sigma_min
      value_min=10.
      theta_max= 0.
      do i=1,n
        call calcul(sigma,M1,us,vs,Mn1,Mn2,theta,M2)
        write(1,110) coef*theta,sigma*coef, us,vs,Mn2,M2,Mn1
        value=abs(M2-1.)

        ! recherche de la valeur de M2 le plus proche de 1.
        if (value.lt.value_min) then
          value_min=value; theta2=theta*coef; sigma2=sigma*coef
          M2min=M2; nmin=i
        end if

        !  angle de deviation maximale
        if (theta.gt.theta_max) then
          theta_max=theta; sigma_max=sigma
        end if

       ! valeur à la nouvelle itération
        sigma=sigma+dsigma
      end do

      write(*,*) 'sigma  =', sigma2,' theta = ',theta2, 'M2 = ',M2min, &
          'valeur ', value_min, '  nmin = ',nmin
      write(2,110) theta2,sigma2,M1

      write(*,*) 'theta_max = ', theta_max*coef, ' sigma_max = ', sigma_max*coef
      write(3,110) theta_max*coef, sigma_max*coef,M1,M2
      write(4,110) theta2, sigma2,M1,M2min
end do

 110 format(8(e12.5,3x))
 close(1); close(2); close(3); close(4)


 contains

!***************
subroutine menus
!***************
implicit none

 write(*,90)
 90 format(50('-'),/,10x, 'Polaire de Busemann et de choc',/,& 
           20x, 'C. Airiau, 2018',/,50('-'))
 write(*,60)
 60 format('Aide : on peut tracer les polaires pour des valeurs dans le tableau ci-dessus uniquement',/, &
          4x,'Il faut entrer 2 paramètres sous la forme Mach, 0    ou    i_debut,n',/,&
          8x,'par exemple  3 , 0 => un seul Mach est calculé pour M1=3',/, &
          8x,'              3 , 5 => les Mach  du tableau  Mach1(3:3+5-1) sont calculés.'/, &
          4x , 'si le second paramètre est 0, le calcul se fait pour toute valeur du Mach donnée',/,&
          4x , 'sinon n est le nombre de Mach calculé')
          
 write(*,*) 'Entrer les deux paramètres M1, n   ou   i_debut,i_fin  :'
 read(*,*) M1,option
 if (option.eq.0) then
   ind1=1
   Mach1(1)=abs(M1);npt=1
 else
   ind1=int(M1); npt=ind1+option-1
   write(*,180) Mach1(ind1:npt)
 end if
 
  180 format('Mach numbers M1 : ',/,3(10(F8.2,2x),/))
end subroutine menus

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


 end program busemann

!*****************************************************
 subroutine  calcul(sigma,M1,us,vs,Mn1,Mn2,theta,M2)
!*****************************************************
 implicit none
 real (kind=8)::sigma, M1,M2,Mn1,Mn2,theta, f
  real (kind=8)::coss,sins,us,vs
 real (kind=8), parameter ::gamma=1.4
 coss=cos(sigma); sins=sin(sigma)
 Mn1= M1*sins
 ! temp= rho1/rho2
 f =  (gamma-1)/(gamma+1) +2./((gamma+1)*Mn1**2)
 us=coss**2+sins**2*f
 vs=sins*coss*(1-f)
 theta=atan(vs/us)
! Mn2 fonction de Mn1
 Mn2=sqrt((1+(gamma-1)/2*Mn1**2)/(gamma*Mn1**2-(gamma-1)/2.))
 M2=Mn2/ sin(sigma-theta)
 end subroutine calcul


