!------------------------------------------------------------------------------
! TITLE         : NACA AIRFOIL GEOMETRY
! PROJECT       : FundAeroSuite
! MODULE        : naca_airfoil_design  
! URL           : 
! AFFILIATION   : Paul Sabatier University, Toulouse
! DATE          : 2016 - 2020
! REVISION      : 2020 V2
!> @author
!> Christophe Airiau
!
! DESCRIPTION:
!>  get geometry of NACA ' and 5 digits  
!!
!! use of a module to define the airfoil type
!!
!! x distribution : regular grid or  regular angle theta  
!------------------------------------------------------------------------------

module naca45_airfoil

type airfoil
integer             :: num
character(len=15)   :: name
character(len=5)    :: digits
real(kind=8)        :: tmax                          ! maximun of thickness
real(kind=8)        :: chord=1.d0                    ! chord
integer             :: n                             ! number of points per side
integer             :: grid=1                        ! 1 : dx = cste, 0 : d\theta=cste.
logical             :: sym=.false.                   ! for symmetrical airfoil
logical             :: TE=.true.
real(kind=8),dimension(:),allocatable   :: x,xu,xl   ! u : upper side (extrados), l: lower size (intrados)
real(kind=8),dimension(:),allocatable   :: yt,yu,yl,yc,dy,theta
end type airfoil

type(airfoil)               :: q
real(kind=8)                :: pi=4.d0*datan(1.d0)
integer                     :: n
integer,allocatable,dimension(:) :: ntmp

contains

!> @brief
!! main subroutine to test different airfil
!**********************
subroutine run_airfoil
!**********************
call title
! symmetrical 4 digits 
call NACA45(0012,101,0)
call SaveAirfoil; call deallocate_airfoil
! non symmetrical 4 digits 
call NACA45(4412,101,0)
call SaveAirfoil; call deallocate_airfoil
! 5 digits, single camber line
call NACA45(24012,101,0)
call SaveAirfoil; call deallocate_airfoil
! 5 digits, double camber line
call NACA45(25112,101,0)
call SaveAirfoil; call deallocate_airfoil

print*,'end of run_airfoil'

end subroutine run_airfoil

!> @brief to deallocate table from airfoil type
!******************************
subroutine deallocate_airfoil
!******************************
!print*,allocated(q%xu),allocated(q%xl),allocated(q%yu)
!print*,allocated(q%yl),allocated(q%yc),allocated(q%yt)
deallocate(q%xu,q%xl,q%yu,q%yl,q%yc,q%yt)
if (allocated(q%x))     deallocate(q%x)
if (allocated(q%dy))    deallocate(q%dy)
if (allocated(q%theta)) deallocate(q%theta)
if (allocated(ntmp))    deallocate(ntmp)
end subroutine deallocate_airfoil

!>@brief define the parameters for the NACA5 airfoil
!************************************************************
subroutine select_parameter_NACA5(ntmp,camber_line,P,M,K1,K21)
!************************************************************
implicit none
integer, intent(in)             :: camber_line
integer,dimension(5),intent(in) :: ntmp
real(kind=8),intent(out)        :: M,P,K1,K21
real(kind=8),dimension(20),parameter :: tab1=(/  &
                    210d0,0.05d0,0.0580d0 ,361.40d0, &
                    220d0,0.10d0,0.1260d0 ,51.640d0, &
                    230d0,0.15d0,0.2025d0 ,15.957d0, &
                    240d0,0.20d0,0.2900d0 , 6.643d0, &
                    250d0,0.25d0,0.3910d0 , 3.230d0 /)
real(kind=8),dimension(4,5),parameter :: Coef1a=reshape(tab1,(/4,5/))
real(kind=8),dimension(5,4),parameter :: Coef1=transpose(Coef1a)
real(kind=8),dimension(20),parameter  :: tab2=(/  & 
221d0,0.10d0    ,0.130d0   ,51.990d0  ,0.000764d0 ,&
231d0,0.15d0    ,0.217d0   ,15.793d0  ,0.00677d0  ,&
241d0,0.20d0    ,0.318d0   , 6.520d0  ,0.0303d0   ,&
251d0,0.25d0    ,0.441d0   , 3.191d0  ,0.1355d0   /)
real(kind=8),dimension(5,4),parameter  :: Coef2a=reshape(tab2,(/5,4/))
real(kind=8),dimension(4,5),parameter :: Coef2=transpose(Coef2a)
integer     :: i


if (ntmp(3).eq.0) then
    K21=0.d0;
    print 100,(Coef1(i,:),i=1,5)
    i=1;
    do while ((camber_line.ne.Coef1(i,1)).and.(i.le.5))
        i=i+1
    end do
    if (i.le.5) then
        print*,'parameters matched for a single camber line'
        P=Coef1(i,2); M=Coef1(i,3); K1=Coef1(i,4)
    else
        print*,' profile not in the database for a single camber line'; stop
    end if

else
    print 110,(Coef2(i,:),i=1,4)
    i=1;
    do while ((camber_line.ne.Coef2(i,1)).and.(i.le.4))
        i=i+1
    end do
    if (i.le.4) then
        print*,'parameters matched for a double camber line'
        P=Coef2(i,2); M=Coef2(i,3); K1=Coef2(i,4); K21=Coef2(i,5)
    else
        print*,' profile not in the database for a double camber line'; stop
    end if
end if
print 120,'P = ', P,'M = ',M,'K1 = ',K1,'K2/K1 = ',K21

100 format(5(4(f12.6,2x),/))
110 format(4(5(f12.6,2x),/))
120 format(4(a6,2x,f12.6,2x))


end subroutine select_parameter_NACA5


!> @brief calculate NACA 4 and 5  geometry
!******************************************
subroutine NACA45(num,nsize,x_distribution)
!*****************************************
implicit none
! the airfoil data are in q
integer,intent(in)                  :: num,nsize,x_distribution
integer                             :: i
real(kind=8)                        :: M,P,T,K1,K21
real(kind=8)                        :: eta,c
integer                             :: n_digits  ! number of digits
integer                             :: camber_line ! NACA 5 parameter
real(kind=8),dimension(0:4)         :: a=(/ 0.2969d0,-0.1260d0,-0.3516d0,0.2843d0,-0.1015d0 /)
integer,parameter,dimension(4)      :: k=(/1,2,3,4/)
q%grid=x_distribution
q%num=num;q%n=nsize;
c=q%chord                               ! chord


n_digits=airfoil_type()
call allocate_vectors(nsize,n_digits)

ntmp=extract_digits(q%num,n_digits) ! digit extraction
select case (n_digits)
case(4)
    M = float(ntmp(4)) / 100.d0                !
    P = float(ntmp(3)) / 10.d0                 !
case(5)
    camber_line=q%num/100
    print*,' camber line number : ',camber_line
    call select_parameter_NACA5(ntmp,camber_line,P,M,K1,K21)
case default
    print*,'error, not a NACA 4 or 5 digits'; stop
end select


T = FLOAT(ntmp(2)*10 + ntmp(1)) / 100.d0   ! thickness
write(*,100),M,P,T
100 format ('# M = ', f4.2, 3x,'P= ',f4.2,3x,'T= ',f4.2)
n=q%n

if (.not.q%TE) a(4)=-0.1036d0

call grid
! thickness equation
do i=1,n
    !q%yt(i)=T/20.d0*(a(0)*sqrt(q%x(i))+a(1)*q%x(i)+a(2)*q%x(i)**2&
    !        +a(3)*q%x(i)**3+a(4)*q%x(i)**4
    eta=q%x(i)/c
    q%yt(i)=T/0.2d0*c*(a(0)*sqrt(eta)+sum(a(k)*eta**k))
end do

q%yc=0.d0

select case (n_digits)

case(4)         ! NACA 4 digits
    if (.not.q%sym) then
        ! cambered line equation
        do i=1,n
            eta=q%x(i)/c
            if (eta.le.P) then
                q%yc(i)=c*M/P**2*eta*(2.d0*P-eta)
                q%dy(i)=2.d0*M/P**2*(P-eta)
                q%theta(i)=atan(q%dy(i))
            else
                q%yc(i)=c*M*(1-eta)/(1.d0-P)**2*(1.d0-2.d0*P+eta)
                q%dy(i)=2.d0*M/(1.d0-P)**2*(P-eta)
                q%theta(i)=atan(q%dy(i))
            end if
        end do
        call transformation
    else
        ! symmetrical airfoil
            q%xu=q%x; q%xl=q%x; q%yu=q%yt; q%yl=-q%yt
    end if
case(5)         ! NACA 5 digits
    if (ntmp(3).eq.0) then
        do i=1,n
            eta=q%x(i)/c
            if (eta.le.M) then
                q%yc(i)=c*K1/6.d0*(eta**3 -3.d0*M*eta**2+M**2*(3.d0-M) * eta)
                q%dy(i)=K1/6.d0*(3.d0*eta**2-6.d0*M*eta +M**2 *(3.d0-M))
            else
                q%yc(i)=c*K1/6.d0*M**3*(1-eta)
                q%dy(i)=-K1/6.d0*M**3
            end if
            q%theta(i)=atan(q%dy(i))
        end do
    else
        do i=1,n
            eta=q%x(i)/c
            if (eta.le.M) then
                q%yc(i)=c*K1/6.d0* ( (eta-M)**3 -K21*(1.d0-M)**3*eta-M**3 * eta+M**3)
                q%dy(i)=K1/6.d0*(3.d0*(eta-M)**2-K21*(1.d0-M)**3-M**3)
            else
                q%yc(i)=c*K1/6.d0*(K21*(eta-M)**3-K21*(1.d0-M)**3*eta-M**3 * eta+M**3)
                q%dy(i+1:n)=K1/6.d0*(3.d0*K21*(eta-M)**2-K21*(1.d0-M)**3-M**3 )
            end if
            q%theta(i)=atan(q%dy(i))
        end do

    end if
    call transformation
end select
end subroutine NACA45

!> @brief transformation to take into account of airfoil thickness
!*************************
subroutine transformation
!*************************
implicit none
q%xu=q%x-q%yt*sin(q%theta)
q%xl=q%x+q%yt*sin(q%theta)
q%yu=q%yc+q%yt*cos(q%theta)
q%yl=q%yc-q%yt*cos(q%theta)
end subroutine transformation


!> @brief define if the airfoil is NACA 4 or NACA5
!********************************
function airfoil_type()  result(nd)
!********************************
implicit none
integer         :: nd
if (q%num.ge.10000) then
    print*,' NACA 5 digits';nd=5;
else
    print*,' NACA 4 digits'; nd=4;
end if

q%digits=carac(nd,q%num)      ! numbers to string 
q%name='NACA'//q%digits             ! airfoil name
write(*,*) 'Airfoil name :',q%name
if (q%num.lt.100) then 
    q%sym=.true.; print*,q%name,' is symmetrical '
else
    q%sym=.false.;print*,q%name,' is non symmetrical '
end if
end function airfoil_type

!> @brief allocation of table in airfoil type
!**************************************
subroutine allocate_vectors(n,n_digits)
!***************************************
implicit none
integer, intent(in)     :: n,n_digits
!print*,allocated(q%x),allocated(q%yc),allocated(q%yt)
allocate(ntmp(n_digits))
allocate(q%x(n),q%yt(n),q%yc(n))
allocate(q%xu(n),q%xl(n),q%yu(n),q%yl(n))
allocate(q%dy(n),q%theta(n))
end subroutine allocate_vectors

!>@brief define the x distribution points
!****************
subroutine grid
!****************
q%x(1)=0.d0;q%x(n)=q%chord;

if (q%grid.eq.1) then
    dx=(q%x(n)-q%x(1))/real(n,kind=8)
    do i=2,n-1; q%x(i)=q%x(i-1)+dx; end do
else
    dtheta=pi/real(n,kind=8); q%theta(1)=0.d0
    do i=2,n-1; q%theta(i)=q%theta(i-1)+dtheta; end do
    q%x=0.5d0*(1.d0-cos(q%theta))
end if

end subroutine grid

!> @brief get digit from a integer
!******************************
function extract_digits(m,n) result(d)
!******************************
integer,intent(in)     :: m
integer,intent(in)     :: n
integer,dimension(n)   :: d
integer                :: i,itmp,expo
itmp=m
do i=n,1,-1
    expo=10**(i-1)
    d(i)=int(itmp/expo)
    itmp=itmp-d(i)*expo
end do

end function extract_digits

!> @brief set the title
!******************************
subroutine title
!******************************
write(*,100)
100 format (50('*'),/,'#',10x,' Airfoil Design ',/,50('*'))
end subroutine title

!> @brief save the Airfoil into a file with its name
!******************************
subroutine SaveAirfoil
!******************************
    implicit none
    integer         :: i
    
    open(1,form='formatted',file=trim(q%name)//'.dat')
    write(1,100) q%name
    write(1,110)
    do i=1,q%n-1; write(1,200) q%xu(i),q%yu(i); end do
    write(1,*)
    do i=1,q%n-1; write(1,200) q%xl(i),q%yl(i); end do
    close(1)
    open(1,form='formatted',file=trim(q%name)//'_carac.dat')
    write(1,100) q%name
    write(1,120)
    do i=1,q%n-1
        write(1,200) q%x(i),q%yc(i),q%yt(i),q%dy(i)
    end do
    close(1)
    
    print*,'Airfoil saved'
    100 format('#',50('*'),/,'#',4x,a10,/,'#',50('*'))
    110 format('#',4x, 'x/c',10x,'y/c')
    120 format('#',4x, 'x/c',7x,'yc/c', 4x,'yt/c')
    200 format(5(e12.5,2x))
    
end subroutine SaveAirfoil
    

!> @details
!!  conversion integer to a string chain  
!!  exemple : 123 --> '123'
!!  N_in : number of digits
!!  i_enter : integer to convert
!****************************
function carac(n_in,i_enter)
!****************************
    implicit none
    integer, intent(in) :: i_enter,n_in
    character(len=n_in) :: carac
    integer             :: i,k,n
    if (n_in.le.1) then
        stop 'problem in carac function'
    endif
    n=n_in-1
    i=i_enter
    do k=n,0,-1
        carac(n-k+1:n-k+1)=char(48+int(i/10**k))
        i=mod(i,10**k)
    end do
end function carac
end module naca45_airfoil

!> @brief main program to test NACA airfoil, an example of use
!=======================
PROGRAM Airfoil_Design
!=======================

use naca45_airfoil
implicit none

call run_airfoil
print*,' normal end of execution of Airfoil_Design'

END PROGRAM Airfoil_Design
