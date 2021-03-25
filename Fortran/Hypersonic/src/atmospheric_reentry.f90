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

!>@details solve parameters of an atmospheric reentry
!
!> @author C. Airiau, august 2015; d'après A. Giovannini
program atmospheric_reentry

implicit none
integer             :: i
integer,parameter   :: n=100
!calcul du mach
real(kind=8),dimension(n)   :: z,a,visc, temp,press,dens
real(kind=8),dimension(n)   :: Velocity,Mach
real(kind=8),dimension(n)   :: flux                     ! heat flux
real(kind=8),dimension(n)   :: rey                      ! Reynolds number
real(kind=8),dimension(n)   :: q                        ! dynamical pressure

real(kind=8)                :: beta=0.14d0              ! beta parameter
real(kind=8)                :: vinf=7.8d0               ! reentry velocity in  km/s
real(kind=8)                :: rho=1.225d0              ! air density (without importance)
real(kind=8)                :: mcx=200.d0               ! parameter  m x Cx
real(kind=8),parameter      :: gamma=0.5d0              ! gamma for the gas
real(kind=8)                :: Lref=1.d0                ! reference length
real(kind=8)                :: dum1,dum2,dum3,dum4
! Need atmospheric data
open(1,form='formatted',file='Atmosphere.dat')
read(1,*);read(1,*)
do i=1,44
    read(1,*)z(i),dum1,dum2,dum3,temp(i),press(i),dens(i), a(i),dum4,visc(i)
end do
close(1)

! velocity 

do i=1,44
    Velocity(i)=vinf*exp(-rho/2./0.14/mCx/gamma*1000*exp(-beta*z(i)))
    Mach(i)=Velocity(i)/a(i)*1000
    rey(i)=Velocity(i)*Lref/visc(i)
    flux(i)=0.094*sqrt(dens(i)/rho)*Velocity(i)**3.
    ! velocity in  km/s
    q(i)=0.5d0*dens(i)*(1000*Velocity(i))**2
end do

dum1=maxval(flux); flux=flux/dum1
dum1=maxval(q); q=q/dum1
open(1,form='formatted',file='reentry.dat')
do i=1,44
    write(1,300) z(i),Velocity(i),Mach(i),rey(i), flux(i), q(i)
end do
300    format(6(1x,e12.5))
close(1)
end program atmospheric_reentry
