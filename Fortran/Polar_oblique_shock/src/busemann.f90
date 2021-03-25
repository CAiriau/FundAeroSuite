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
!------------------------------------------------------------------------------
! TITLE         : BUSEMANN POLAR
! PROJECT       : FundAeroSuite
! MODULE        : busemann
! URL           : 
! AFFILIATION   : Paul Sabatier University, Toulouse
! DATE          : 2011 - 2020
! REVISION      : 2020 V2
!> @author
!> Christophe Airiau
!
! DESCRIPTION:
!>  generate Busemann polar for oblique shock wave to create a database for plotting
!------------------------------------------------------------------------------

!> @details:
!! generate data to plot:  
!! * Busememann polar
!! * shock wave versus deviation angle theta
!! 
!! Parameters:
!!       Mach number (to choose in a given list) 
!!   
!! Ouputs :
!!     Three file types:
!!      * data_Machxyz.dat : all data for a given Mach (xyz/100 = Mach)   
!!      * unitary_downstream_Mach.dat  : sonic line (downstream Mach number = 1)
!!      * theta_max.dat     : Mach number at maximal deviation angle
!!    
!! Note:
!!  1 : refers to upstream
!!  2 : refers to downstream
program busemann

implicit none

real (kind=8)                       :: M1,sigma,sigma_min,sigma_max,dsigma,value_min
real (kind=8)                       :: us,vs,Mn1,Mn2,M2,theta,M2min,value,theta2,sigma2
real (kind=8)                       :: theta_max
integer, parameter                  :: n=1000    ! must be large for accuracy (related to shock angle step)
integer, parameter                  :: npt_max=40
real (kind=8), dimension(npt_max)   :: Mach1
real (kind=8),parameter             :: pi = 4.d0*atan(1.d0), coef = 180.d0/pi

integer i,npt,k,ind1,option,nmin
character(len=50)                  :: filename

Mach1(1:15)=(/1.01, 1.03, 1.1, 1.15, 1.20, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7, 1.8, 2. /)
Mach1(16:30)=(/2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4., 4.5, 5., 6., 10., 100./)
npt=30

write(*,180) Mach1(1:npt)
write(*,*) 'Number of Mach number available : ',npt

180 format('upstream Mach numbers M1 : ',/,3(10(F8.2,2x),/))

call menus

! open files

! file status can be changed in some uses.
! open(3,form='formatted',file='theta_max.dat',status='old', access='append')
! open(2,form='formatted',file='m2.dat',status='old', access='append')
open(4,form='formatted',file='unitary_downstream_Mach.dat' )
write(4,170)
170  format ('# theta (deg)',3x,'sigma (deg)',8x,'M1',9x,'unitary M2')
open(3,form='formatted',file='theta_max.dat' )
write(3,160)
160  format ('# theta (deg)',3x,'sigma (deg)',8x,'M1',12x ,'M2')

write(*,130) 
130 format(5x,'sigma',10x,'theta',13x,'M2',11x, 'M2 - 1',3x, 'n_min', 3x,'theta_max',6x, 'sigma_max')

do k=ind1,npt
    M1=Mach1(k)
    sigma_min = asin (1./M1)
    sigma_max = 0.5*pi
    dsigma=(sigma_max-sigma_min)/float(n-1)

    ! one generic file per Mach number 
    filename='data_M'//carac(int(log10(100*M1)+1),int(100*M1+0.001))//'.dat'
    open(1,form='formatted',file=trim(filename))
    write(1,100)M1
    100 format(/,'# Busemann polar for M = ',f11.4,/,'#',4x,'theta',10x,'sigma', 10x,'u2/V1', &
          10x,'v2/V1',10x,'Mn2', 10x, 'M2', 10X, 'Mn1', /,'#', /)
    sigma=sigma_min
    value_min=10.
    theta_max= 0.
    do i=1,n     ! loop on the shock angle
        call calcul(sigma,M1,us,vs,Mn1,Mn2,theta,M2)
        write(1,110) coef*theta,sigma*coef, us,vs,Mn2,M2,Mn1
        value=abs(M2-1.)

        ! look for the angle value wich provides the nearest M2 number from 1
        if (value.lt.value_min) then
          value_min=value; theta2=theta*coef; sigma2=sigma*coef
          M2min=M2; nmin=i
        end if

        !  maximal deviation angle
        if (theta.gt.theta_max) then
          theta_max=theta; sigma_max=sigma
        end if

        ! shock angle at the next iteration
        sigma=sigma+dsigma
    end do

    write(*,120) sigma2,theta2,M2min,value_min,nmin, theta_max*coef, sigma_max*coef
    write(3,110) theta_max*coef, sigma_max*coef,M1,M2
    write(4,110) theta2, sigma2,M1, M2min
end do

110 format(8(e12.5,3x))
120 format(4(f12.5,3x),i4,2(f12.5,3x))
close(1); close(2); close(3); close(4)

contains

!> @details
!! menu and calculate the Mach numbers to use
!***************
subroutine menus
!***************


  implicit none
  real(kind=8) :: error
  real(kind=8), parameter :: error_min = 0.001
  integer      :: i

  write(*,90)
  90 format(50('-'),/,10x, 'Busemann polar and shock curve',/,& 
            20x, 'C. Airiau, 2011-2020',/,50('-'))
  write(*,60)
  60 format('Help : only the Mach given in the table below are available to get data for the polar',/, &
            4x,'2 parameters are asked, and must written as  < Mach, 0 >  or more generally < i_debut, n >',/,&
            8x,' for instance  3 , 1  =>   single Mach = 3',/, &
            8x,'               3 , 4 =>  first Mach = 3 and the next 3 Mach are the next one in the table'/, &
            8x,' for instance  3 , 0  =>  Mach = 3 to all the next values in the Mach table',/, &
            4x , 'n = 0  : all the next value',/,&
            4x , 'n <> 0 : number of Mach calculated')
            
  write(*,*) 'Enter the two parameters  M1, n  (example 2.2, 2 ):'
  read(*,*) M1, option

  error=1
  i=1
  error = abs(M1-Mach1(i))
  do while ((error.gt.error_min).and.(i.le.npt-1)) 
    i = i +1
    error = abs(M1-Mach1(i))
    !write(6,*)'error = ', error, 'for M =', Mach1(i)
  end do 

  if  (error.gt.error_min) then
    stop 'Mach number not found in the Mach table'
  else 
    ind1 = i
    write(6,*) 'Mach found in the table : ', Mach1(i)
  end if

  if ((option.ne.0).and.(ind1+option-1.le.npt)) then 
    npt=ind1+option-1
  end if
  write(*,180) Mach1(ind1:npt)
  180 format('Mach numbers M1 : ',/,3(10(F8.2,2x),/))
end subroutine menus

!> @details
!! convert a integer to a string chain
!!
!!  exemple : 123 --> '123'
!****************************
function carac(n_in,i_enter)
!****************************

    implicit none
    integer, intent(in) :: i_enter   !> integer number to convert
    integer, intent(in) :: n_in      !> number of digits
    character(len=n_in) :: carac
    integer             :: i, k, n
  
    if (n_in .le. 1) then
      stop  'problem in  carac function'
    endif
  
    n = n_in-1
    i = i_enter
    do k = n, 0, -1
      carac(n-k+1:n-k+1) = char(48+int(i/10**k))
      i = mod(i,10**k)
    end do
  end function carac


end program busemann

!< details
!! oblique shock state
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


