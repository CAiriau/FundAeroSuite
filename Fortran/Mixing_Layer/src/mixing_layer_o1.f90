!...............................................................
!>@details
!! ........ PROGRAMM TO SOLVE AUTOSIMILAR MIXING LAYER ......... 
!!
!! **LAMINAR and TURBULENT flow**
!! 
!! Governing equation :
!!              \f$ \displaystyle \frac{1}{2}~f~f''+ f'''~=~0 \f$ 
!!   
!!
!! **Heat transfer** problem is also solved, but not validated for turbulent flow
!!
!! 
!!
!!   \f$ \begin{array}{lll}
!!   {\rm variable-function | regime}  \qquad& {\rm   laminar} \qquad   & {\rm   turbulent} \vspace{0.3cm}  \\ \hline \hline
!!   \eta          &   \eta_{\ell}      &   \eta_t \vspace{0.3cm}\\
!!   f(\eta)       &  f_{\ell}          &   f_t    \vspace{0.3cm}\\
!!   U_{ref}       &   U_1              & \displaystyle \frac{U_1+U_2}{2}
!!  \end{array} \f$
!!
!! We get the same results by a change of variable and functions:
!!
!!  \f$\eta_{\ell}~=~\alpha~\eta_t\f$ and  \f$ f_{\ell} = f_t~/~\alpha\f$ with
!!   \f$\displaystyle \alpha^2~=~\frac{2}{1+\lambda}\f$ and 
!!  \f$\displaystyle \lambda~=~\frac{U_2}{U_1}\f$
!!
!! Methodology:         
!!  * shooting method
!!  * first order Euler scheme for the integration
!!  * for boundary conditions, use \f$ \lambda \f$ parameters      
!! .............................................................

!> details
!! general module to solve the autosimilar profile of a mixing layer
module mod_mixing_layer
!      eta  : normal direction coordinate 
!      deta : grid step in this direction 
!      option_f0  = 0 : calculus with the value of f(0)
!      option_f0  = 1 : calculus with iteration on  f(0)

integer,parameter       :: il=120000       ! maximal size of the grid (depend on eta_max and deta)
real(kind=8),parameter  :: const=0.5d0     ! constante in the main ODE.
real(kind=8),parameter  :: prt=0.5         ! Prandtl number
integer                 :: option_f0 = 1   ! Newton algorithm to get f0, otherwise f0=0
real(kind=8), parameter :: eta_max = 20    ! height of the grid in + or - eta
real(kind=8)            :: deta = 0.002    ! step in eta
real(kind=8)            :: lambda          ! parameter related to infinite upstream velocity
real(kind=8)            :: coef=1.         ! 2/(1+lambda) for turbulent flow, 1 for laminar flow
integer                 :: turbulence=1    ! by default = 1 - True for turbulent flow
real(kind=8)            :: delta_u1=1      ! velocity at + infinity, default value = 1
real(kind=8)            :: delta_u2=0      ! velocity at - infinity, default value = 0
real(kind=8)            :: alpha, alpha2   ! scale factor and its square
real(kind=8)            :: u_ref           ! reference velocity
real(kind=8)            :: reduced_uinf_negative, reduced_uinf_positive
logical                 :: heat_transfer_problem=.false. 
integer                 :: imax            ! maximun size of eta lower or eta upper 


contains

!> @details
!! main subroutine to read the parameter and solve the problem
subroutine solve_mixing_layer
    
    implicit none
    real(kind=8), dimension(0:2*il+2) :: f,u,g,eta, g1, g2
    real(kind=8)            :: s
    real(kind=8)            :: u0,f0,testf0
    real(kind=8)            :: Vinf_positive,Vinf_negative,ured
    real(kind=8)            :: ratio_velocity, scale_factor
    integer                 :: i,k ,k1,k2
    real(kind=8)            :: x_zero,u_zero,g_zero, scale
    integer,parameter       :: dk=100
    character(len=11), parameter :: output_file = 'outputs.out' 

    write(6,105)                ! title
    
    call read_parameters        ! inputs from a file

    lambda=delta_u2/delta_u1
    imax=int(eta_max/abs(deta)) 
    reduced_uinf_positive=2.d0/(lambda+1.d0)
    reduced_uinf_negative=2.d0*lambda/(lambda+1.d0)
    alpha2=2.d0/(lambda+1.d0)    
    alpha=sqrt(alpha2)

    if (turbulence.eq.1) then
        coef= alpha2
        u_ref = (delta_u2+delta_u1)/2
        print*,'regime          : turbulent'
        scale_factor=1.
        scale=1. 
    else
        coef=1.
        u_ref = delta_u1
        scale_factor = sqrt(2.0)
        scale        = 2.0
        print*,'regime          : laminar'
    endif

    print*,'main parameters '
    print*,'delta_u1        : ', delta_u1 
    print*,'delta_u2        : ', delta_u2
    print*,'lambda = u2/u1  : ', lambda
    print*,'deta            : ', deta
    print*,'imax            : ', imax 
    print*, 'U (inf)  / Um  : ', reduced_uinf_positive
    print*, 'U (-inf) / Um  : ', reduced_uinf_negative
    print*, 'c (coef)       : ', coef
    print*, 'U ref          : ', u_ref
    print*, 'alpha          : ', alpha
    print*, 'alpha^2        : ', alpha2
    print*, 'scale          : ', scale

    if (imax.gt.il) then
        write(*,*) 'il < imax, imax = ',imax
        stop "modify deta or il"
    end if

    if (option_f0.eq.0) then
        f0 = 0
    else
        f0 = 0.1
    end if

    open(25,form='formatted',file='convergence.out')
    write(25,200)

    !
    ! ......... integration of the velocity profile .............
    !
    call  shooting_f0(k,k1,k2,u0,s,testf0,f0,eta,f,u,g)
    !
    !.....................  OUTPUTS .............................
    !
    102 format ('#',4x,'eta',5x,'f',7x,' df/deta',7x,'d2f/deta^2', &
                9x,'g1',9x,'g2',7x,'reduced U')
    100 format(i6,3x,9(e15.8,1x))
    105 format(50('*'),/,'*',4x,'Autosimilar shear layer',/, &
        ' LAMINAR / TURBULENT FLOW ',/,50('*'),/)
    200 format('#  n',7x,'error s',7x,'error u0',10x,'s',13x,'u0',11x,&
        'eta1',12x,'eta2',10x,'error F0')

    call matlab(k,s,u0,eta,f,u,g)      

    ! ........ integration of the temperature profile ..................
    !
    !         g'''/Prt + f g" = 0,  g'(+infty ) = 1, g'(-infty = 0)
    !        change of function :  g1 = g', g2 = g",  and  g3 = g'''
    ! ..................................................................
    !
    if (heat_transfer_problem) then
        call shooting__heat_transfer_problem(k,k1,k2,eta,f,g1,g2)
    else
        g1=0.d0;g2=0.d0
    end if
    !
    ! search intersection with y=0 axis or the line of  v = 0
    !
    call zeros(k,eta,f,u,g,x_zero,u_zero,g_zero)

    Vinf_negative = (eta(0)*u(0)-f(0))/scale
    Vinf_positive = (eta(k)*u(k)-f(k))/scale
    ratio_velocity=Vinf_negative/Vinf_positive

    open(unit=20,file=output_file)
    write(20,102)

    ! k is a cut-off on the grid

    do i=0,k,dk
        ured=(u(i)- reduced_uinf_negative)/(reduced_uinf_positive- reduced_uinf_negative)    ! reduced nondimensional velocity
        write(20,100)i,eta(i),f(i),u(i),g(i),g1(i),g2(i), &
                    (eta(i)*u(i)-f(i))/scale, ured
    end do

    write(*,*) 'eta = ',eta(0),' ==> Vinf . C(x) = ',Vinf_negative
    write(*,*) 'eta = ',eta(k),' ==> Vinf . C(x) = ',Vinf_positive

    write(*,*) 'Vinf ratio and inverse : ', ratio_velocity, 1.d0/ratio_velocity
    write(*,*) 'eta() ', eta(k1),'==> U() =',u(k1),'==> f() =',f(k1),'k1 = ',k1
    write(*,*) 'eta() ', eta(k2),'==> U() =',u(k2),'==> f() =',f(k2),'k2 = ',k2
    write(*,*) 'eta() ', eta(k),'==> U() =',u(k),'==> f() =',f(k),'k2 = ',k
    close(20)
    write(*,*)"Book table :"
    write(*,300) string("f(0)"),f(k2)
    write(*,300) string("f'(0)"),u(k2)
    write(*,300) string("eta_i"),x_zero
    write(*,300) string("f'(eta_i)"),u_zero
    write(*,300) string("f^(2)(eta_i)"),g_zero
    write(*,300) string("delta_omega"),coef*(1.d0-lambda)/g_zero
    write(*,300) string("v ( infinity)"),Vinf_positive 
    write(*,300) string("v (-infinity)"),Vinf_negative 


    open(1,form='formatted',file='positive.dat')
    write(1,110)
    do i=k2,k,dk
        write(1,100)i,eta(i),f(i),u(i),g(i),g1(i),g2(i),(eta(i)*u(i)-f(i))/scale
    end do
    close(1)

    open(1,form='formatted',file='negative.dat')
    write(1,110)
    do i=0,k2,dk
        write(1,100)i,eta(i),f(i),u(i),g(i),g1(i),g2(i),(eta(i)*u(i)-f(i))/scale
    end do
    close(1)
    close(25)
    write(6,*) "normal end of execution"

    300 format('#',5x,a30,2x,':',2x,f12.4)
    110 format ('#',3x,'i',3x,'eta',3x,'f',3x,'u',5x,'g',3x,'g1',3x,'g2',3x,'v/vref')
end subroutine solve_mixing_layer

function string(in)
    character(len=*)    :: in
    character(len=30)   :: string
    integer             :: n
    n=len(in)
    string='..............................'   
    string(1:n) = in
    string(n+1:n+2) = '  '
     
end function string

!> @details
!! read the parameters from the input file **mixing_layer.in**
subroutine read_parameters
    implicit none
    logical               :: file_exist
    inquire ( file = 'mixing_layer.in', exist = file_exist )
    if (file_exist) then
        write(*,*) 'the input file has been found'
        open(1,form='formatted',file='mixing_layer.in')
        read(1,*)
        read(1,*) turbulence
        read(1,*) delta_u1
        read(1,*) delta_u2
        read(1,*) deta
        
        close(1)
    else
        write(*,*) 'input file not found'
        write(*,*) 'use of default parameters'
    end if

end subroutine read_parameters
! .............................................................................
!> @details
!! the boundary condition s is solved with a Newton algorihtm (Shooting method)
!! for the fluid equation
!! s = f''(O)
!! two target conditions  at + and - infinity
!!   f  :  autosimilar function
!!
!!   u  : df/d eta
!!
!!   g  : d2 f / d eta 2
subroutine shooting_(k,k1,k2,u0,s,testf0,f0,eta,f,u,g)
! to get a very goo accuracy take imax large and epseta small
    implicit none
    integer, intent(out) :: k   !> total number of points used to define the mixing layer
    real(kind=8), dimension(0:2*il+2) :: f,u,g,eta
    real(kind=8), dimension (0:il) :: f1,u1,g1
    real(kind=8), dimension (0:il) :: f2,u2,g2
    integer      :: i,k2,k1,iter_max,n,nstep
    real(kind=8) :: tests,testu,epss,epsu
    real(kind=8) :: gu1,gu2,gs1,gs2,uold,s0
    real(kind=8) :: du1, du2,detu,dets,du,ds,det
    real(kind=8) :: s,u0,f0,deta1,deta2
    real(kind=8) :: testf0,eta1,eta2,vplus,vmoins

    ! initialisation of Newton method
    tests=1.
    testu=1.
    iter_max=100
    epss=1e-7
    epsu=1e-7
    testf0=1.d0
    !
    !- Loop on  n: calculul of  s  
    !
    !   u value at infinity      
    deta1 = deta
    deta2 = -deta    
    u0=(delta_u1+delta_u2)/2.
    s=0.1 
    n=1      
    
    do while (((tests.gt.epss).or.(testu.gt.epsu)).and.(n.lt.iter_max))

        call integration(deta1,u0,f0,f1,u1,g1,k1,s,gs1,gu1)
        print*,'result y > 0 ', k1,s,u1(k1)        
        call integration(deta2,u0,f0,f2,u2,g2,k2,s,gs2,gu2)
        print*,'result y < 0 ', k2,s,u2(k2)
        ! @ infinity   u1 = delta_u1
        ! @ -infinity  u2 = delta_u2
        ! du1=delta_u1-u1(k1)
        ! du2=delta_u2-u2(k2)
        du1=coef-u1(k1)
        du2=lambda*coef-u2(k2)
        det=gu1*gs2-gu2*gs1
        detu=du1*gs2-gs1*du2
        dets=gu1*du2-du1*gu2
        du=detu/det
        ds=dets/det
        s0=s
        uold=u0
        s=s+ds
        u0=u0+du
        tests=abs(ds)
        testu=abs(du)
        eta1 = float(k1)*deta1
        eta2 = float(k2)*deta2
        ! vmoins = delta_u2/delta_u1*(eta2*u2(k2)-f2(k2))  
        ! vplus =  eta1*u1(k1)-f1(k1)
        vmoins = lambda*coef*eta2-f2(k2)  
        vplus =  coef*eta1-f1(k1)
        testf0=(lambda*vmoins+vplus)**2
        write(25,110) n,tests,testu,s,u0,eta1,eta2,testf0
    end do
    if (n.eq.iter_max) then
        stop 'no convergence in shooting_'
    end if
    110   format(i5,2x,10(e12.5,3x))
    ! end n loop for Newton 

    nstep=1
    write(6, *) 's = ', s, 'test V_infinity = ', testf0
    k=-1
    ! lower part solution
    do i=k2,1,-nstep
        k=k+1
        eta(k)=-float(i)*deta
        f(k)=f2(i)
        u(k)=u2(i)
        g(k)=g2(i)
        !write(22,100)eta(k),f(k),u(k),g(k)
    end do
    ! upper part solution
    do i=0,k1,nstep
        k=k+1
        eta(k)=float(i)*deta
        f(k)=f1(i)
        u(k)=u1(i)
        g(k)=g1(i)
        if (i.eq.0) then
            print*,'@ eta = 0, f, u, du/deta : ',eta(k), f(k),u(k),g(k)
        end if 
        !write(22,100)eta(k),f(k),u(k),g(k)  
    end do
    
    print*, 'number of grid points : ', k
    print*, 'k2 (eta > 0) = ', k2, ' k1 (eta <0) =  ',k1,' k1+k2 = ',k2+k1
    !100    format(5(e12.5,2x))
end  subroutine shooting_


! ......................................................................
!> @details 
!! integration of the ODE with an Euler order 1 Finite Difference Scheme
!
subroutine integration(deta_layer,u0,f0,f,u,g,k,s,gu,fu)
!   f  : auto-similar function
!   u  : df/d eta
!   g  : d2 f / d eta 2
!   t  : d3 f / d eta 3
    implicit none
    real(kind=8), dimension (0:il), intent(out) :: f,u,g
    integer      :: i,k
    real(kind=8) :: testeta 
    real(kind=8) :: epseta = 0.5d-5 
    real(kind=8) :: t,deta_layer,f0
    real(kind=8) :: gf,gg,gu,gt,s,u0,ff,fu,fg,ft

    testeta=1.    
    ! Initial conditions
    ! system (2)
    f(0)=f0
    u(0)=u0
    t=0.
    g(0)=s
        
    ! gradient of the primal system  with respect to s=f''(0)          
    ! system (3)
    gf=0.
    gu=0.
    gt=0.          
    gg=1.
        
    ! gradient of the primal system with respect to u0=f'(0)          
    ! system (4)
    ff=0.
    fu=1
    ft=0.          
    fg=0
    !
    !- Loop on i: system (2), (3) and (4) are solved simultaneously
    !
    k=0
    do i=1,imax
        if (testeta .lt. epseta) exit
        ! system (2)
        f(i)=f(i-1)+deta_layer*u(i-1)
        u(i)=u(i-1)+deta_layer*g(i-1)
        g(i)=g(i-1)+deta_layer*t
        t= -const*f(i)*g(i)   
        !
        ! system (3)
        !   
        gf=gf+deta_layer*gu
        gu=gu+deta_layer*gg
        gg=gg+deta_layer*gt
        gt= -const*(f(i)*gg+g(i)*gf)  
        !
        ! system (4)   
        ff=ff+deta_layer*fu
        fu=fu+deta_layer*fg
        fg=fg+deta_layer*ft
        ft= -const*(f(i)*fg+g(i)*ff)  
        ! test to exit of the loop at the mixing layer thickness
        testeta=abs(g(i))
        !  write(*,*) i,f(i),u(i),g(i),gf,gu,gg,s 
        k=k+1
    end do
end subroutine integration


! ....................................................
!> save the profile in a matlab file,
!! 1 point every 500 points is kept
subroutine matlab(k,s,u0,eta,f,u,g) 

    implicit none  
    real(kind=8), dimension(0:2*il+2) ::  f,u,g,eta
    integer         :: k,nstep,n,i
    real(kind=8)    :: s,u0

    open(1,form='formatted',file='profile.m')
    write(1,100)
    100    format('function [eta,f,u,g]=profile')
    write(1,*)'% Generated in fortran from turbulent_mixing.f90 '
    write(1,*)'% Christophe Airiau, mars 2006 - mai 2020'
    write(1,*)'% data of the solution'
    write(1,200) delta_u1,delta_u2
    write(1,201) u0,s

    200 format('% U1 = ',2x,f7.5 ,4x,'U2=',2x, f7.5)
    201 format('% U(0) = ',2x,f7.5,4x,'dU/d eta(0) = ',2x,f7.5)
    nstep = int(k/500)
    if (nstep.eq.0) then
        nstep = 1
    end if 
    write(1,101)
    write(1,*)'% eta grid '
    n=0
    do i=0,k,nstep
        n=n+1
        write(1,102)n,eta(i)
        end do
    write(1,*)'% function f(eta)'
    n=0
    do i=0,k,nstep
        n=n+1
        write(1,103)n,f(i)
        end do
    write(1,*)'% velocity u(eta)'
    n=0
    do i=0,k,nstep
        n=n+1
        write(1,104)n,u(i)
        end do
    write(1,*)'% g(eta)'
    n=0
    do i=0,k,nstep
        n=n+1
        write(1,105)n,g(i)
    end do
    close(1)
    101 format ('%',4x,'eta',5x,'f',7x,' df/deta',7x,'d2f/deta^2')
    102 format('eta(',i5,') = ',e12.5,' ;')  
    103 format('f(',i5,') = ',e12.5,' ;') 
    104 format('u(',i5,') = ',e12.5,' ;') 
    105 format('g(',i5,') = ',e12.5,' ;')   
end subroutine matlab


!.....................................................
!> @details
!! first part of the shooting method :
!! solve f0.
subroutine shooting_f0(k,k1,k2,u0,s,testf0,f0, eta,f,u,g)

    implicit none
    real(kind=8), intent(inout)  :: f0
    integer, intent(out)      :: k,k1,k2
    
    real(kind=8), intent(out) :: u0
    real(kind=8), dimension(0:2*il+2) ::  f,u,g,eta
    
    real(kind=8) ::  s
    real(kind=8) :: testf0,testf0g
    real(kind=8) :: df0,f0g,dtest,f0corr,test,testrep
    integer      :: i,imax

    open(24,form='formatted',file='F0.out')
    write(24,200)

    if (option_f0.eq.0) then 
        print*,"direct => shooting_"
        call shooting_(k,k1,k2,u0,s,testf0,f0,eta,f,u,g) 
        return
    else
        test = 1.
        i = 0
        imax = 50
        testf0 =1.
        df0 = 1.d-6
        testrep = 1.d-12
        do while ((i.le.imax).and.(abs(testf0).gt.testrep))  
            print*,' iteration sur f(0) '
            call shooting_(k,k1,k2,u0,s, testf0,f0,eta,f,u,g) 
            f0g = f0+df0
            call shooting_(k,k1,k2,u0,s, testf0g,f0g,eta,f,u,g)
            dtest = (testf0g-testf0)/df0
            f0corr = -testf0/dtest
            test = abs(f0corr/f0)
            write(24,100)i,f0,test,testf0
            i = i+1
            f0 = f0+f0corr
        end do
    end if
    100 format(i4,3(e12.5,3x))
    200 format('#   it',5x,"F0",10x,"err F0",10x,"Test V infinity")
    close(24)
end subroutine shooting_f0


! ...................................................
!> @details
!! integration of the ODE for the heat transfer problem

!! **methodology, notations:**
!!
!! g1n = g1(eta_N), eta_N = eta(k1+k2),
!! g10 = g1(eta_0), eta_0 = eta(0)
!!
!! from eta = 0 to eta(k2) => lower layer from  1 to k2
!!
!! from k2 to k=k1+k2 => upper layder but eta = 0, eta(k2+1) 
!!
!! a = g1(0) and b = g2(0)
!!
!! search a and b such that  g1n = 1. and g10 = 0.
!! Need to calculate the ODE gradients with respect to a and b 
!! refer as dag1,dag2 et dbg1 dbg2

subroutine int_heat_transfer_problem(k,k1,k2,a,b,eta,f,g1,g2, g1n,g10,dag1n,dbg1n,dag10,dbg10)

    implicit none

    real(kind=8), dimension(0:2*il+2) ::  f,g1,g2,eta
    real(kind=8) ::  deta_tmp,dag1,dag2,dag3,dbg1,dbg2,dbg3,g1n,g10
    real(kind=8) ::   dag1n,dbg1n,dag10,dbg10
    real(kind=8) ::  a,b,g3
    integer      :: k,k1,k2,i
    ! a = 0.1
    ! b = 0.1
    print*,'K1,k2, int th',k1,k2,k
    print*,eta(0),eta(k2),eta(k)
    ! condition @ eta = 0 
    ! primal systeme 
    g1(k2) = a
    g2(k2) = b
    g3 = 0.
    ! integration of the upper layer
    ! gradient wrt / a
    dag1 = 1
    dag2 = 0
    dag3=0.
    ! gradient wrt / b
    dbg1 = 0
    dbg2 = 1
    dbg3=0.
    
    do i = k2+1,k
        deta_tmp = eta(i)-eta(i-1)
        ! primal system      
        g1(i)=g1(i-1)+deta_tmp*g2(i-1)
        g2(i)=g2(i-1)+deta_tmp*g3
        g3=- Prt*f(i)*g2(i)
        ! gradient wrt / a     
        dag1=dag1+deta_tmp*dag2
        dag2=dag2+deta_tmp*dag3
        dag3=- Prt*f(i)*dag2
        ! gradient wrt / b       
        dbg1=dbg1+deta_tmp*dbg2
        dbg2=dbg2+deta_tmp*dbg3
        dbg3=- Prt*f(i)*dbg2
    
    end do
    g1n = g1(k)
    dag1n = dag1
    dbg1n =dbg1     
    print*,'upper layer @ ', k2+1,i 
    ! integration of the upper layer
    ! gradient wrt / a
    dag1 = 1
    dag2 = 0
    dag3=0.
    ! gradient wrt / b
    dbg1 = 0
    dbg2 = 1
    dbg3=0. 
    do i = k2-1,0,-1
        deta_tmp = eta(i)-eta(i+1)
        
        g1(i)=g1(i+1)+deta_tmp*g2(i+1)
        g2(i)=g2(i+1)+deta_tmp*g3
        g3=- Prt*f(i)*g2(i)   
        ! gradient wrt / a
        dag1=dag1+deta_tmp*dag2
        dag2=dag2+deta_tmp*dag3
        dag3=- Prt*f(i)*dag2  
        ! gradient wrt / b 
        dbg1=dbg1+deta_tmp*dbg2
        dbg2=dbg2+deta_tmp*dbg3
        dbg3=- Prt*f(i)*dbg2              
    end do       
    print*,'lower layer @ ', k2,i +1
    g10=g1(0)
    dag10 = dag1
    dbg10 =dbg1             
end subroutine int_heat_transfer_problem

!* ........................................
!> @brief
!! main subroutine to solve heat transfer
!! not validated for turbulent flow
subroutine shooting__heat_transfer_problem(k,k1,k2,eta,f,g1,g2)
!
! g1n = g1(eta_N),   eta_N = eta(k1+k2)
! g10 = g1(eta_0),   eta_0 = eta(0)
! eta = 0  @ eta(k2), 
! lower layer : from 1 to k2
! upper layer from k2 to k=k1+k2
!      
! @ eta = 0, eta(k2+1) 
!       a = g1(0) et b = g2(0)
! search a and b such that  g1n = 1. and g10 = 0.
! Need to calculate the ODE gradients with respect to a and b 
! refer as dag1,dag2 et dbg1 dbg2

    implicit none
    real(kind=8), dimension(0:2*il+2) ::  f,g1,g2,eta
    real(kind=8) ::   g1n,g10
    real(kind=8) ::   dag1n,dbg1n,dag10,dbg10
    real(kind=8) ::  a,b,testa,testb
    real(kind=8) ::  test,det_a,det_b,da,db,det,eps
    integer k,k1,k2,i,itmax

    !  initial values for a and b (g1 and g2 @ eta = 0.)
    a=0.1
    b=0.1
    print*, 'k1,k2 , th', k1,k2
    eps= 0.001
    itmax = 4
    i=0
    test = 1.

    do while  ( (i.lt.itmax).and.(test.gt.eps))
        i=i+1
        call int_heat_transfer_problem(k,k1,k2,a,b,eta,f,g1,g2, g1n,g10,dag1n,dbg1n,dag10,dbg10)   
        det =  dag1n*dbg10-dag10*dbg1n
        det_a =   (1.-g1n)*dbg10+ g10*dbg1n
        det_b = -dag1n*g10-dag10*(1-g1n)
        da = det_a/det
        db = det_b/det
        testa = abs(da/a)
        testb = abs(db/b)
        write(24,100)i,a,b,testa,testb,g1n,g10
        100   format (i3,2x,' a = ',f12.5,'  b= ',f12.5, &
            ' testa ',e12.5,' testb ',e12.5, ' G1n ',f12.5,' G10 ',f12.5)   
        if (testa.lt.testb) then
            test = testb
        else
            test  = testa
        end if
        a = a + da
        b = b + db                
    end do    

end subroutine shooting__heat_transfer_problem


! ..............................................................
!> @details
!! solve f=0, and give the values of x,u and g at this location
subroutine zeros(n,x,f,u,g,x_zero,u_zero,g_zero)

    implicit none
    integer,intent(in)                     :: n
    real(kind=8),dimension(n),intent(in)   :: x,f,u,g
    real(kind=8),dimension(n)              :: grad
    integer                                :: nbre_zeros,i
    real(kind=8)                           :: alpha
    real(kind=8),intent(out)               :: x_zero,u_zero,g_zero
    open(13,form='formatted',file='essai.out')
    x_zero=0.d0;u_zero=0.d0;g_zero=0.d0
    nbre_zeros=0
    print*,'seach the median line'
    grad(1:n-1)=f(1:n-1)*f(2:n)
    do i=1,n-1
        if (grad(i).lt.0.d0) then
            nbre_zeros=nbre_zeros+1
            alpha=-f(i)/(f(i+1)-f(i))
            x_zero=(1-alpha)*x(i)+alpha*x(i+1)
            u_zero=(1-alpha)*u(i)+alpha*u(i+1)
            g_zero=(1-alpha)*g(i)+alpha*g(i+1)
            print*, 'a zero has been found',x(i),' and ', x(i+1), ' @ ',x_zero
            write(*,100) i,i+1
            write(*,200) x(i),x(i+1)
            write(*,200) f(i),f(i+1)
            write(*,200) u(i),u(i+1)
            write(*,200) g(i),g(i+1)
            exit
        end if
    end do

    if (i.ge.n-1) then
        write(*,*) 'no zero found'
    end if
    write(*,110)'x zero',x_zero
    write(*,110)'u zero',u_zero
    write(*,110)'g zero',g_zero
    close(13)
    110 format('#', a10,4x,e12.5)
    100 format ('#', 5x,'|',6x,'i',6x,'|',6x,'i+1',/,10x,i8,7x,i8)
    200 format(6x,e12.5,3x,e12.5)
end subroutine zeros

end module mod_mixing_layer


!> @details
!! program to test the module
program main_mixing_layer

use mod_mixing_layer
call solve_mixing_layer

end program main_mixing_layer
