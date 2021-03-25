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

!  piston method : integral calculus by solving an ode
!  C. Airiau, Janvier 2016, d'après A. Giovannini
program piston_method

implicit none
integer                                 :: n
real(kind=8),allocatable,dimension(:)   :: f,g,h,eta,s
real(kind=8),parameter                  :: gam=1.4d0
real(kind=8),parameter                  :: eps=1.d-4    ! convergence error
real(kind=8)                            :: eta_step,som,c0,df,dg,dh,som1
real(kind=8)                            :: m,alpha
integer                                 :: i, i_max, k
!  flow parameter m  < 1, 
!  k = 0 : 2D, k=1, conical flow
! initialisation
eta_step=0.000001d0
n=int(1.d0/eta_step)+1
print*,' maximal number of point  = ',n
allocate(f(n),g(n),h(n),eta(n),s(n))
eta(1)=1.d0
f(1)=2.d0/(gam+1.d0)
g(1)=f(1)
h(1)=(gam+1.d0)/(gam-1.d0)

print *,'enter the value m with  0 < m < 1'
read(*,*) m
alpha=(m-1.d0)/m
print *,'enter the value k (0 or 1)'
read(*,*)k
select case (k)
case(0)
    print*, '2D flow'
case(1) 
    print*, 'conical flow'
case default
    stop 'bad value for k'
end select
! order 1 integration (Euler scheme)

print*,'first method in  piston.dat'
i=1
s(i)= eta(i)**k*eta_step*(h(i)*f(i)*f(i)/2.d0+1.d0/(gam-1.d0)*g(i))
do  i=2,n
    eta(i)=eta(i-1)-eta_step
    c0=eta(i-1)-f(i-1)
    df=-(alpha*f(i-1)*h(i-1)*h(i-1)*c0/gam/g(i-1) &
        + float(k)*f(i-1)*h(i-1)/eta(i)+2.d0*alpha*h(i-1)/gam)/&
        (h(i-1)-c0*c0*h(i-1)*h(i-1)/g(i-1)/gam)
    f(i)=f(i-1)-df*eta_step
    c0=eta(i)-f(i)
    dg=h(i-1)*(c0*df-alpha*f(i))
    g(i)=g(i-1)-dg*eta_step
    dh=h(i-1)/gam*(dg/g(i)-2.d0*alpha/c0)
    h(i)=h(i-1)-dh*eta_step

    if(c0.lt.eps)  then
        exit
    else if (eta(i).lt.0) then
        print*,'i = ',i,' eta(i) = ',eta(i)
        stop 'problem'
    else if (i.ge.n)  then
        print*,'i = ',i,' eta(i) = ',eta(i)
        stop 'i = n, increase n or the step in eta'
    end if
    s(i)= eta(i)**k*eta_step*(h(i)*f(i)*f(i)/2.d0+1.d0/(gam-1.d0)*g(i))
end do
i_max=i
som=0.5d0*(s(1)+s(i_max))+sum(s(2:i_max-1))

call outputs(som,'piston.dat')


print*,'second method in piston1.dat'
i_max=1
do  i=2,n
    eta(i)=eta(i-1)-eta_step
    call edo(k,gam,alpha,eta(i-1),f(i-1),g(i-1),h(i-1),df,dg,dh)
    f(i)=f(i-1)-df*eta_step
    g(i)=g(i-1)-dg*eta_step
    h(i)=h(i-1)-dh*eta_step
    c0=eta(i)-f(i)
    if(c0.lt.eps)  then
        exit
    else if (eta(i).lt.0) then
        print*,'i = ',i,' eta(i) = ',eta(i)
        stop 'problem'
    else if (i.ge.n)  then
        print*,'i = ',i,' eta(i) = ',eta(i)
        stop 'i = n, increase n or the step in eta'
    end if
    s(i)= eta(i)**k*eta_step*(h(i)*f(i)*f(i)/2.d0+1.d0/(gam-1.d0)*g(i))
end do
i_max=i
som1=0.5d0*(s(1)+s(i_max))+sum(s(2:i_max-1))
call outputs(som1,'piston.dat')
print *,'error in the integral  = ',abs(1.d0-som/som1)

deallocate(f,g,h,eta,s)

contains

!******************************
subroutine outputs(som,filename)
!******************************
!> @brief ouputs on screen and in a file

implicit none
real(kind=8),intent(in)     :: som
character(len=*),intent(in) ::filename
real(kind=8)                :: coor

print *,'stop wall'
print *,'c0 conver. = ', c0
print *,'eta(im)    = ',eta(i_max)
print *,'f(im)      = ',f(i_max)
print *,'g(im)      = ',g(i_max)
print *,'h(im)      = ',h(i_max)
print *,'il         = ',i_max
print *,'integral   = ',som
open(1,form='formatted',file=trim(filename))
write(1,200)
do  i=1,i_max
    coor=(eta(i)-eta(i_max))/(1.-eta(i_max))
    write(1,100) coor,f(i)/f(1),g(i)/g(1),h(i)/h(1)
end do
close(1)

100 format(4(e15.8,2x))
200 format('#',50('*'),/,'#',3x,'eta*',10x,'f/f1',10x,'g/g1',10x,'h/h1')
end subroutine outputs

end program piston_method

!> ODE to solve the heat flux q
!*********************************************
subroutine edo(k,gam,alpha,eta,f,g,h,df,dg,dh)
!*********************************************
!
implicit none
integer, intent(in)             :: k
real(kind=8), intent(in)        :: gam,alpha
real(kind=8), intent(in)        :: eta,f,g,h
real(kind=8), intent(out)       :: df,dg,dh
real(kind=8)                    :: den,c0,a1,a2

c0=eta-f
a1=real(k,kind=8)*f*h/eta
a2=-alpha*f
den=c0**2 * h - g * gam
df= (-a2 * c0 * h** 2 + a1 * g * gam + 2.d0 * alpha * g * h)/(den*h)
dg=g*(a1 * c0 * gam - a2 * gam * h + 2.d0 * alpha * c0 * h)  /den
dh=h * (a1 * c0** 2 - a2 * c0 * h + 2.d0 * alpha * g) / (den*c0)

end subroutine edo
