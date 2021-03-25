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
!> get the pressure coefficient in hypersonic flow
!! see the chapter 14 of the book "Aerodynamique Fondamentale"
!! @author C. Airiau, aout 2015
program  hyperso

implicit none
real(kind=8),parameter      :: gam=1.4
real(kind=8),parameter      :: gam4=(gam+1)/4.d0
integer,parameter           :: n=201
real(kind=8),dimension(3,2) :: K0
integer                     :: i,j
character(len=8)            :: filename
real(kind=8)                :: dk,k,Kp

K0(1,1)=0.01d0;   K0(1,2)=10.d0
K0(2,1)=1.0d0;    K0(2,2)=10.d0
K0(3,1)=0.01d0;   K0(3,2)=1.d0

do i=1,3
    select case (i)
    case (1)
        filename='eq27.dat'
    case (2)
        filename='eq29.dat'
    case (3)
        filename='eq35.dat'
    end select

    dk=(log10(K0(i,2))-log10(K0(i,1)))/real(n-1,kind=8)

    open(1,form='formatted',file=trim(filename))
    do j=1,n
        k=10**(dk*real(j-1,kind=8)+log10(K0(i,1)))
        select case (i)
        case (1)
            Kp=2.d0*(gam4 + sqrt(gam4**2+1.d0/k**2))
        case (2)
            Kp=gam+1.d0
        case (3)
            Kp=2.d0/k
        end select
        write(1,'(2(e15.8,3x))')k,Kp
    end do
    close(1)
end do

end program hyperso


