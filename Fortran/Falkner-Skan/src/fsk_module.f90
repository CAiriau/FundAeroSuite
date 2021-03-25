!------------------------------------------------------------------------------
! TITLE         : FALKNER6SKAN FLOW + CONTINUATION APPROACH
! PROJECT       : FundAeroSuite
! MODULE        : 
! AFFILIATION   : Paul Sabatier University, Toulouse
! DATE          : 2010 - 2020
! REVISION      : 2020 V2
!> @author
!> Christophe Airiau
!
! DESCRIPTION:
!>   Solve the Falkner-Skan equation with very high accuracy
!!
!!   * Get the boundary profile
!!
!!   * Get the boundary layer characteristic
!! 
!!   * Able to  treat multiple value of the m parameter by continuation  
!!
!!   * main equation : 
!!    
!!      - \f$\displaystyle f'''+ f~f''+ \beta(1-f'^2) = 0  \f$
!!      - the reference thickness (lenght) is \f$\displaystyle  \delta = \sqrt{\frac{\nu x}{Ue} }  \sqrt{\frac{2}{m+1}} \f$
!!      - details are in "Couche limite laminaire", Cousteix, Cepadues edition, 1987
!!      - in our book "Aérodynamique fondamentale", Giovannini, Airiau, 2016 is solved:
!!        \f$\displaystyle f'''+\frac{m+1}{2} f f''+ m (1-f'^2)=0  \f$ with \f$\displaystyle \delta = \sqrt{\frac{\nu x}{Ue} } \f$
!!      - There is a scale factor between the solution of the both ODE
!!      - in the program, the variable "xk" is this scale factor, and it is written in "save. characteristics"         
!!      - The boundary layer profiles are saved with the Cousteix formulation.  
!!
!!   * Methodogies 
!!      - Runge Kutta fourth order integration
!!      - double integrated Newton method with derivative of the systems (2 parameters)               
!!      - Continuation method to follow beta or m in the intervall of resolution called 'Arclength Continuation'
!!  * Run 
!!      - parameters can be entered in the input file **fsk_continuation.in**
!!      - or in the terminal by answering to the question
!!  * Ouputs
!!      - BL_characteristic.dat
!!      - profile00000000.out
!!      - Convergence.dat
!!      - output.out
!!  * Notations
!!      -  f'= u and  f''= g
!!
!******************************************************************
MODULE  fsk_module
!> @brief main program with parameters
IMPLICIT NONE

REAL(KIND=8),PARAMETER  :: tolerance_s=1.d-9     ! tolerance for s = f(1)
REAL(KIND=8),PARAMETER  :: tolerance_g=1.d-5     ! tolerance for g = g(infinity)
REAL(KIND=8),PARAMETER  :: CutOff_U=1.d-4        ! tolerance for g = g(infinity) ou t(infinity)
REAL(KIND=8),PARAMETER  :: CutOff_end_BL=1.d-2   ! to define  delta_0.99
INTEGER,PARAMETER       :: m= 3                  ! system order  +1 (why + 1 ?)
REAL(KIND=8)            :: deta                  ! step in eta
REAL(KIND=8)            :: eta_max = 10.d0       ! eta_max, default value, change later
INTEGER                 :: imax
INTEGER                 :: size_max              ! size of the vectors f, u, g
REAL(KIND=8)            :: dL                    ! length of the arc
INTEGER                 :: ind_end_BL
INTEGER                 :: mmax=10,lmax=30       ! maximal number of the loops, default value

REAL(KIND=8)            :: s_init, s0, s1, ds 
REAL(KIND=8)            :: beta0, dbeta,beta1


INTEGER                 :: cas,flag
INTEGER                 :: flag_continuation=1  ! default value to not run continuation if not asked

INTEGER                 :: type_schema=1         ! 0 : order 1, 1 : RK 4

REAL(KIND=8),DIMENSION(:,:), ALLOCATABLE ::  y
INTEGER,PARAMETER         :: nr=8
REAL(KIND=8),DIMENSION(nr):: results

CONTAINS

SUBROUTINE Main_Continuation(read_option)
integer, intent(in) :: read_option !< = 1 : test to read parameters in file
REAL(KIND=8)            :: L1 
LOGICAL                 :: divergence

write(*,110)
110 format(50('*'),/,'* FSK flow, Continuation arclength method', /,50('*'),/)

!  output files
OPEN(14,FORM='FORMATTED',FILE='Convergence.dat')
OPEN(13,FORM='FORMATTED',FILE='BL_characteristics.dat')
WRITE(13,200)
OPEN(20,form='formatted', file = 'output.out' )
write(20,'(a1,6x,a4,10x,a1,10x,a)')'#','beta','s','m'

IF (read_option.eq.1) THEN
    CALL read_data(flag)     ! read in the input file, flag is an output.
ELSE 
    flag = 0
END IF
IF (flag.eq.0) THEN
    ! if the input file is not found, the program will run with default parameters.
    beta0 = 0.d0
    dbeta =-0.005d0
    deta  = 2.0d-4      ! step in eta
    dL    = 0.005d0     ! length of the arc
    s_init = 0.5d0
    print*, 'Use default parameters'
END IF
imax=int (eta_max / deta)-1
size_max=imax+1 
write(*,'(a,5x,i7)') 'Size of the vectors f, u, g', size_max

ALLOCATE (y(0:size_max,m))
IF (type_schema.eq.1) THEN
    PRINT*,'Integration scheme:  RK 4'
ELSE
    PRINT*,'Integration scheme : Euler order 1'
END IF

! calculus of the first convergence point named M0, of the function F(s,beta)=0
! .
! Calculus of  M0(beta=0 ; s= s_Blasius  or  another point for another beta)
!-----------------------------------------
WRITE ( * , '(a,/,a)' ) 'calcul de M0', '************'
cas=0
write(*,100) beta0,m_factor(beta0)
CALL Solve_FSK(beta0, s_init, s0) 

WRITE(20,103) beta0,s0,m_factor(beta0)

if (flag_continuation.eq.1) then
    CLOSE(13);CLOSE(14); CLOSE(20)
    DEALLOCATE (y)
    stop 'normal end of execution after a single B.L. profile'
end if
!
! Calculus of  M1(beta= beta_init+dbeta, s1)
!-------------------------------------------

cas=1
WRITE ( * , '(a,/,a)' ) 'calcul de M1', '************'
beta1=beta0+dbeta
write(*,100) beta1,m_factor(beta1)
CALL Solve_FSK (beta1, s0, s1) 
WRITE(20,103) beta1,s1,m_factor(beta1)

! data : ds, L1 
!------------------

ds = s1 - s0
L1 = sqrt (dbeta * dbeta + ds * ds) 
WRITE ( *, '(a,e12.6)' ) 'arc length L1 = ', L1
WRITE ( *, * )

! stop 'beta1 : ok'
!
! continuation method
!------------------------
CALL Continuation (s1, beta1, s0, beta0, divergence) 
CLOSE(20);CLOSE(13);CLOSE(14)
DEALLOCATE (y)

WRITE ( * , * ) 'normal end of execution after the continuation process' 

100 FORMAT('# beta =',2x,f10.5,4x,'m =',2x,f10.5)
103 FORMAT (5(e15.8,1x)) 
200 FORMAT('#     beta',10x,'m',13x,'H',12x,'s',10x,'d1/d',9x,'d2/d',&
    9x,'Cf',10x,'xk',9x,'Eta Max',7x,'Lambda',6x,'End_BL')

! ************************************************
! end of main program
! ************************************************

END SUBROUTINE Main_Continuation


!> @brief beta to m factor 
!**********************
FUNCTION m_factor(beta)
!**********************
    REAL(KIND=8) :: m_factor, beta
    m_factor=beta/(2.d0-beta)
END FUNCTION m_factor

!> @brief m to beta factor
!**********************
FUNCTION beta_factor(mloc)
!**********************
    REAL(KIND=8) :: beta_factor, mloc
    beta_factor=2.d0*mloc/(mloc+1.d0)
END FUNCTION beta_factor

!> @details
!! if flag = 1 : parameters are read in the input file
!!
!*************************
SUBROUTINE read_data(flag)
!*************************
IMPLICIT NONE
LOGICAL               :: file_exist
INTEGER, INTENT(OUT)  :: flag

INQUIRE ( file = 'fsk_continuation.in', exist = file_exist )
IF (file_exist) then
    WRITE(*,*) 'The input file has been found'
    flag=1
    OPEN(1,FORM='formatted',FILE='fsk_continuation.in')
    READ(1,*)
    READ(1,*) flag
    READ(1,*) flag_continuation
    READ(1,*) beta0
    READ(1,*) dbeta
    READ(1,*) deta
    READ(1,*) dL  
    READ(1,*) s_init
    READ(1,*) Lmax
    READ(1,*) eta_max
    READ(1,*) type_schema
    CLOSE(1)
    select case (flag)
    case(0,1)
        print*,'flag = ',flag
    case default
        flag=0
    end select
ELSE
    WRITE(*,*) 'Input file not found'
    WRITE(*,*) 'Use of default parameters'
    flag=0
END IF
END SUBROUTINE read_data

!> @brief
!! use of Runge-Kutta algorithm to solve primal system (main equation)
!**************************
SUBROUTINE rk4(beta,dt,tin,y0,dy0,yout,dyout)
!**************************
IMPLICIT NONE
REAL(KIND=8),INTENT(IN)                     :: dt,beta,tin
REAL(KIND=8),DIMENSION(m),INTENT(OUT)       :: yout,dyout
REAL(KIND=8),DIMENSION(m),INTENT(IN)        :: y0,dy0

REAL(KIND=8),DIMENSION(m)                   :: dydt,ddydt
REAL(KIND=8),DIMENSION(m)                   :: k1,k2,k3,k4,dk1,dk2,dk3,dk4
REAL(KIND=8),DIMENSION(m)                   :: y_tmp1,y_tmp2,y_tmp3,dy_tmp1,dy_tmp2,dy_tmp3

! t : theta, y
! First sub-step
CALL ode(beta,y0,dy0,tin,dydt,ddydt)
k1= dt*dydt; y_tmp1=y0+0.5d0*k1
dk1= dt*ddydt; dy_tmp1=dy0+0.5d0*dk1
! Second sub-step
CALL ode(beta,y_tmp1,dy_tmp1,tin+0.5d0*dt,dydt,ddydt)
k2= dt*dydt; y_tmp2=y0+0.5d0*k2
dk2= dt*ddydt; dy_tmp2=dy0+0.5d0*dk2
! Third sub-step
CALL ode(beta,y_tmp2,dy_tmp2,tin+0.5d0*dt,dydt,ddydt)
k3= dt*dydt; y_tmp3=y0+k3
dk3= dt*ddydt; dy_tmp3=dy0+dk3
! Fourth sub-step
CALL ode(beta,y_tmp3,dy_tmp3,tin+dt,dydt,ddydt)
k4= dt*dydt; yout=y0+(k1+2.0d0*(k2+k3)+k4)/6.0d0
dk4= dt*ddydt; dyout=dy0+(dk1+2.0d0*(dk2+dk3)+dk4)/6.0d0
END SUBROUTINE rk4

!>  @details
!! use of Runge-Kutta algorithm to solve the system for continuation 
!! (derivative of primal system w.r.t. beta)
!*************************************************************
SUBROUTINE rk4_cont(beta,dt,tin,y0,dy0,dy0b,yout,dyout,dyoutb)
!*************************************************************
IMPLICIT NONE
REAL(KIND=8),INTENT(IN)                     :: dt,beta,tin
REAL(KIND=8),DIMENSION(m),INTENT(OUT)       :: yout,dyout,dyoutb
REAL(KIND=8),DIMENSION(m),INTENT(IN)        :: y0,dy0,dy0b

REAL(KIND=8),DIMENSION(m)                   :: dydt,ddydt,dydtb
REAL(KIND=8),DIMENSION(m)                   :: k1,k2,k3,k4,dk1,dk2,dk3,dk4,dk1b,dk2b,dk3b,dk4b
REAL(KIND=8),DIMENSION(m)                   :: y_tmp1,y_tmp2,y_tmp3,dy_tmp1,dy_tmp2,dy_tmp3
REAL(KIND=8),DIMENSION(m)                   :: dy_tmp1b,dy_tmp2b,dy_tmp3b

! t : theta, y
! First sub-step
CALL ode_cont(beta,y0,dy0,dy0b,tin,dydt,ddydt,dydtb)
k1= dt*dydt; y_tmp1=y0+0.5d0*k1;
dk1= dt*ddydt; dy_tmp1=dy0+0.5d0*dk1
dk1b= dt*dydtb; dy_tmp1b=dy0b+0.5d0*dk1b
! Second sub-step
CALL ode_cont(beta,y_tmp1,dy_tmp1,dy_tmp1b,tin+0.5d0*dt,dydt,ddydt,dydtb)
k2= dt*dydt; y_tmp2=y0+0.5d0*k2
dk2= dt*ddydt; dy_tmp2=dy0+0.5d0*dk2
dk2b= dt*dydtb; dy_tmp2b=dy0b+0.5d0*dk2b
! Third sub-step
CALL ode_cont(beta,y_tmp2,dy_tmp2,dy_tmp2b,tin+0.5d0*dt,dydt,ddydt,dydtb)
k3= dt*dydt; y_tmp3=y0+k3
dk3= dt*ddydt; dy_tmp3=dy0+dk3
dk3b= dt*dydtb; dy_tmp3b=dy0b+dk3b
! Fourth sub-step
CALL ode_cont(beta,y_tmp3,dy_tmp3,dy_tmp3b,tin+dt,dydt,ddydt,dydtb)
k4= dt*dydt; yout=y0+(k1+2.0d0*(k2+k3)+k4)/6.0d0
dk4= dt*ddydt; dyout=dy0+(dk1+2.0d0*(dk2+dk3)+dk4)/6.0d0
dk4b= dt*dydtb; dyoutb=dy0b+(dk1b+2.0d0*(dk2b+dk3b)+dk4b)/6.0d0
END SUBROUTINE rk4_cont

!> @details
!! ODE of the primal equation and it derivatives w.r.t. s parameter
!************************************************************
SUBROUTINE ode(beta,x,x_d,tloc,dxdt,dxdt_d) 
!************************************************************
IMPLICIT NONE
REAL(KIND=8),INTENT(IN)                  :: tloc,beta
REAL(KIND=8),DIMENSION(m),INTENT(IN)     :: x,x_d
REAL(KIND=8),DIMENSION(m),INTENT(OUT)    :: dxdt,dxdt_d
REAL(KIND=8)                             :: tmp
tmp=tloc   ! unuseful, to keep tloc in the function argument only
dxdt=0.d0;dxdt_d=0.d0
dxdt(1) = x(2)                      !f' = u
dxdt(2) = x(3)                      !u' = g
dxdt(3) = - x(1)*x(3)-beta *(1.d0-x(2)*x(2))    ! g' = - f f'' -beta (1-f'^2) 

dxdt_d(1) = x_d(2)                      !df' = du
dxdt_d(2) = x_d(3)                      !du' = dg
dxdt_d(3) = -x_d(1)*x(3)-x(1)*x_d(3)+ 2.d0 *beta *x_d(2)* x(2)   ! t' = - f df'' - df f'' -2 beta f' df' 
END SUBROUTINE  ode

!> @details
!! ODE of the primal equation and it derivatives w.r.t. s and beta
!! for the continuation approach
!************************************************************
SUBROUTINE ode_cont(beta,x,x_d,x_b,tloc,dxdt,dxdt_d,dxdt_b) 
!************************************************************
IMPLICIT NONE
REAL(KIND=8),INTENT(IN)                  :: tloc,beta
REAL(KIND=8),DIMENSION(m),INTENT(IN)     :: x,x_d,x_b
REAL(KIND=8),DIMENSION(m),INTENT(OUT)    :: dxdt,dxdt_d,dxdt_b
REAL(KIND=8)                             :: tmp
tmp=tloc   ! unuseful, to keep tloc in the function argument only
dxdt=0.d0;dxdt_d=0.d0
dxdt(1) = x(2)                      !f' = u
dxdt(2) = x(3)                      !u' = g
dxdt(3) = - x(1)*x(3)-beta *(1.d0-x(2)*x(2))    ! g' = - f f'' -beta (1-f'^2) 

dxdt_d(1) = x_d(2)                      !df' = du
dxdt_d(2) = x_d(3)                      !du' = dg
dxdt_d(3) = -x_d(1)*x(3)-x(1)*x_d(3)+ 2.d0 *beta *x_d(2)* x(2)   ! t' = - f df'' - df f'' -2 beta f' df' 
dxdt_b(1) = x_b(2)                      !df' = du
dxdt_b(2) = x_b(3)                      !du' = dg
dxdt_b(3) = -x_b(1)*x(3)-x(1)*x_b(3)+ 2.d0 *beta *x_b(2)* x(2)+x(2)**2-1.d0
END SUBROUTINE  ode_cont

!> @details
!! main subroutine to solve the Falkner-Skan equation
!! with a Newton algorithm to determine s : boundary condition for eta \to \infty
!**********************************
SUBROUTINE Solve_FSK(beta, s_in, s_out) 
!**********************************
IMPLICIT NONE
INTEGER      :: i, k, n, nmax
REAL(KIND=8) :: test_eta, test_s, beta, s_in, s_out, s 
REAL(KIND=8) :: s_old,eta_loc
REAL(KIND=8),DIMENSION(0:size_max,m) :: dy
REAL(KIND=8),DIMENSION(m) :: grady,gradys
LOGICAL      :: test
!
!      y(:,1)  : auto similar function
!      y(:,2)  : df/d eta
!      y(:,3)  : d2 f / d eta 2
!
!-INITIALISATION 
!---------------
!    To get a very good accuracy, we have to set a very large value to
!      imax (between 30 000 to 100 000) and 
!      tolerance_g, the tolerance very small
test=.true.
WRITE(14,*)'# beta = ',beta
test_eta = 1.d0; test_s = 1.d0; 
nmax = 100
s = s_in

!- START of the n LOOP : solve s
!--------------------------------
DO n = 1, nmax 

    IF (test_s.lt.tolerance_s)  exit
    y=0.d0;dy=0.d0
    y(0,1:m) =(/0.d0,0.d0,s/)
    dy(0,1:m)=(/0.d0,0.d0,1.d0/)
    test_eta = 1.d0
!
!- START of the i LOOP: solve simultaneouslty the ODE and its gradient / s
!------------------------------------------------------------------------
    k = 0
    DO i = 1,imax 
        eta_loc=dfloat(i)*deta
        IF (eta_loc.gt.eta_max+deta) exit
        if (type_schema.eq.1) then
            CALL rk4(beta,deta,eta_loc,y(i-1,1:m),dy(i-1,1:m),y(i,1:m),dy(i,1:m))
        else
            CALL ode(beta,y(i-1,1:m),dy(i-1,1:m),eta_loc,grady,gradys)
            y(i,:)=y(i-1,:)+deta*grady
            dy(i,:)=dy(i-1,:)+deta*gradys
        end if
        test_eta = abs(y(i,m)) 
        if (test.and.(abs(test_eta).le.tolerance_g)) then
            write(*,*) 'tolerance on u" reached for  i = ',i
            write(*,*) 'eta                            = ',eta_loc 
            test=.false.
        end if
        k = k + 1 
    END DO 

!- END of i LOOP
!----------------

!           ...... Calculus of s(n+1) .......................

    s_old = s
    print*,' s= , k ',s,k,i
    s = s - (y(k,2) - 1.d0) / dy(k,2)

!           ...... CONVERGENCE TEST for s ...................

    test_s = abs (s - s_old)
    WRITE(14,200) test_s, s_old
    if (isnan(s)) then
        write(*,*) 's is not a number, s_old=',s_old
        stop 'divergence'
    end if
    if (abs(y(i,2)).gt.10.d0) then
        write(*,*) 'explosion for i= ',i,' beta = ', beta
        stop 'not a number'
    end if

ENDDO

if (n.eq.nmax) stop 'Divergence'

!-END of n LOOP
! --------------

WRITE (* ,100 )s,test_s,deta*float(i)
s_out = s
print*," convergence : s = ",s
CALL save_profile(k,cas,beta)
CALL save_characteristics(k,beta)

100 format('Convergence, s = ',f11.8,3x,'test = ',e12.5,' etamax ',e12.5,/)
200 FORMAT(2(e15.8,3x))

END SUBROUTINE Solve_FSK

!> @details
!! main subroutine to perform the continuation algorithm
!********************************************************************
SUBROUTINE Continuation (s1, beta1, s_init, beta_init,divergence)
!********************************************************************

IMPLICIT NONE 
REAL(KIND=8),PARAMETER      :: tolerance_f=1.0d-5    ! tolerance on f

INTEGER  i, k, k_beta, j, p
LOGICAL test_L 

REAL(KIND=8),INTENT(IN) :: s_init,beta_init,s1,beta1
LOGICAL, INTENT(INOUT)  :: divergence

REAL(KIND=8) :: delta_L
REAL(KIND=8) :: beta, beta_arc 
REAL(KIND=8) :: s, s_arc
REAL(KIND=8) :: s_new, s_old, beta_new, beta_old 

REAL(KIND=8) :: test_eta, test_f1, test_f2 

REAL(KIND=8) :: delta_s, delta_beta 
REAL(KIND=8) :: pu,  qu, eta_loc
REAL(KIND=8) :: p1, q1, f1, f2, det,test_s 

REAL(KIND=8),DIMENSION(0:size_max,m) :: dy   ! dy/ds
REAL(KIND=8),DIMENSION(0:size_max,m) :: dyb  ! dy/dbeta
REAL(KIND=8),DIMENSION(m) :: grady,gradys,gradyb


!      f  : auto similar function
!      u  : df/d eta
!      g  : d^2 f / d eta^2
!      t  : d^3 f / d eta^3

!      pf : p=df/ds
!      pu : dp/d eta
!      pg : d^2 p/d eta^2
!      pt : d^3 p/d eta^3

!      qf : q=df/dbeta
!      qu : dq/d eta
!      qg : d^2 q/d eta^2
!      qt : d^3 q/d eta^3

!
!-INITIALISATION
!---------------
divergence=.false.
delta_L  = 0.d0
s_arc    = 0.d0
beta_arc = 0.d0
beta_old = beta_init
s_old    = s_init
beta_new = beta1
s_new    = s1
test_L   = .false.

!- START of L LOOP
!==================

k_beta = 0
WRITE(*,120)
DO j = 1, lmax

    IF (test_L) exit
    ! initialisation
    ! --------------

    delta_s    = s_new - s_old
    delta_beta = beta_new - beta_old
    delta_L    = dsqrt ( delta_s**2 + delta_beta**2)
    s_arc      = s_new + (dL / delta_L) * delta_s
    beta_arc   = beta_new + (dL / delta_L) * delta_beta

    !- CALCULUS of  q1= df2/dbeta and  p1=df2/ds
    !-------------------------------------------

    q1 = delta_beta / (delta_L * dL)
    p1 = delta_s / (delta_L * dL)

    !- CALCULUS of df1/dbeta=qu  and df1/ds=pu
    !------------------------------------------

    test_f1 = 1.d0
    test_f2 = 1.d0
    s       = s_arc
    beta    = beta_arc
    if (isnan(beta)) then
            WRITE(*,*) 'Not a number for beta =',beta
            return
    end if

    WRITE(14,*)'# beta = ',beta
    DO p = 1, mmax

        IF ( (test_f1.lt.tolerance_f) .and. (test_f2.lt.tolerance_f) ) exit

        y(0,1:m) =(/0.d0,0.d0,s/)
        dy(0,1:m)=(/0.d0,0.d0,1.d0/)
        dyb(0,1:m)=(/0.d0,0.d0,0.d0/)

        !- START of i LOOP :solve simultaneouslty the system  (f), (p) and (q)  
        !---------------------------------------------------------------------   

        test_eta = 1.d0
        k = 0
        DO i = 1, imax
            eta_loc=float(i)*deta
            IF (eta_loc.gt.eta_max) exit
            if (type_schema.eq.1) then
                CALL rk4_cont(beta,deta,eta_loc,y(i-1,1:m),dy(i-1,1:m),dyb(i-1,1:m),y(i,1:m),dy(i,1:m),dyb(i,1:m))
            else
                CALL ode_cont(beta,y(i-1,1:m),dy(i-1,1:m),dyb(i-1,1:m),eta_loc,grady,gradys,gradyb)
                y(i,:)=y(i-1,:)+deta*grady
                dy(i,:)=dy(i-1,:)+deta*gradys
                dyb(i,:)=dyb(i-1,:)+deta*gradyb
            end if
        !           ......TEST OF CONVERGENCE OF THE BOUNDARY LAYER FRONTIER ....   
            test_eta = abs(y(i,3)) 
            if (abs(y(i,2)).gt.10.d0) then
                write(*,*) 'explosion for i = ',i,' beta = ', beta
                stop 'not a number'
            end if
            if (isnan(test_eta)) then
                write(*,*) 'NaN i= ',i,y(i,2),y(i-1,2),y(i-2,2),dfloat(i)*deta
                stop 'not a number'
            end if
            k = k + 1 
        enddo 
        !if ((k.lt.100).or.(k.ge.imax)) then
        if (k.lt.100) then
            divergence=.true.
            write(*,*) 'divergence in continuation'
            write(*,*) 'k = ',k
            return
        end if

        !- END OF i LOOP
        !----------------
        pu=dy(k,2)
        qu=dyb(k,2)

        f2 = delta_s / delta_L * (s - s_new) / dL + delta_beta / delta_L * (beta - beta_new) / dL - 1.d0
        f1 = y(k,2) - 1.d0
        test_f1 = abs (f1)
        test_f2 = abs (f2)

        det = (q1 * pu) - (p1 * qu)
        if (abs(det).lt.1e-7) stop 'NaN'
        test_s=- ( (q1 * f1) - (qu * f2) ) / det
        beta = beta - ( (pu * f2) - (p1 * f1) ) / det
        WRITE(14,200) test_s, s,det
        if (isnan(s)) then
            WRITE(*,*) 'Not a number ',det,y(k,2),k
            return
        end if
        s = s +test_s

    enddo 

    if (mod(k_beta,5).eq.0)  CALL save_profile(k,k_beta,beta)
    CALL save_characteristics(k,beta)

    !- NEW COEFFICIENTs s and beta
    !------------------------------

    WRITE( * ,  110 ) k_beta, beta, s, m_factor(beta),float(k)*deta,float(i)*deta,test_eta,k,i
    WRITE(20, 103) beta, s, m_factor(beta)

103 FORMAT (5(e15.8,1x)) 
110 FORMAT(3x,i4,6(3x,e12.5),2(3x,i8))
120 FORMAT('#    k',7x,'beta',13x,'s',12x,'m')
200 FORMAT(3(e15.8,3x))

    k_beta = k_beta + 1 
    s_old = s_new 
    beta_old = beta_new 
    s_new = s 
    beta_new = beta 
enddo 

!- END OF L LOOP
!=======================

END SUBROUTINE Continuation

!> @details
!! conversion integer to string chain,
!!
!! example : 123 --> '123'
!****************************
function charac(N_in,i_enter)
!****************************
! implemented by C Airiau, march 2010
! from Burkardt subroutines
!  
! written on Nin digits

integer, intent(in):: i_enter,N_in
character(len=N_in):: charac
integer            :: i,k,n
! n=int(log10(float(i)))
if (N_in.le.1) then 
    stop 'problem ni charact function'
endif
n=N_in-1
i=i_enter
do k=n,0,-1
    charac(n-k+1:n-k+1)=char(48+int(i/10**k))
    i=mod(i,10**k)
end do
return 
end function charac

!> @details
!! the B.L. profile are save in a output file
!*********************************
SUBROUTINE SAVE_PROFILE(k,ind,beta)
!*********************************

IMPLICIT NONE
INTEGER, INTENT(IN) :: k,ind
INTEGER             :: i
REAL(KIND=8)        :: beta
LOGICAL             :: test
test=.false.
write(*,*) 'save profile case ',ind, 'beta= ',beta
OPEN(1,FORM='formatted',FILE='profile'//charac(8,ind)//'.out')
write(1,100) beta
do i=0,k
    if (mod(i,100).eq.0) then
        write(1,'(5(e12.5,2x))') deta*real(i,kind=8),y(i,:)
    end if
    if (abs(y(i,2)-1.d0).le.CutOff_U) exit
    if ((.not.test).and.(abs(y(i,2)-1.d0).le.CutOff_end_BL)) then
        ind_end_BL=i; test=.true.
        write(*,*) 'end of Boundaray layer at eta= ',float(ind_end_BL)*deta, &
            'index = ',ind_end_BL
    end if
end do

CLOSE(1)

print*,'B.L. profile has been saved successfully'
100 format('# beta = ',f10.5,/,'# ',3x,'eta',12x,'f',12x,'u',12x,'g')
END SUBROUTINE SAVE_PROFILE

!> @details
!! calculate the usual boundary layer characteristics
SUBROUTINE save_characteristics(k,beta)

IMPLICIT NONE
INTEGER, INTENT(IN)     :: k
REAL(KIND=8),INTENT(IN) :: beta
INTEGER                 :: i
REAL(KIND=8)            :: delta1,delta2,H,xm,MaxEta,xk,s,Lambda

print*,'Boundary layer data are saved in a output file'
s=y(0,3)
xm=m_factor(beta)
MaxEta=ind_end_BL*deta
xk=sqrt(2.d0/(xm+1.d0))
delta1=0.d0
delta2=0.d0
do i=0,k-1
    delta1 =delta1+(1.-y(i,2)+1-y(i+1,2))*deta*0.5
    delta2 =delta2+(y(i,2)*(1.-y(i,2))+y(i+1,2)*(1-y(i+1,2)))*deta*0.5
end do
H=delta1/delta2
!         delta = sqrt(nu x/ue)
!         delta1 : displacement thickness
!         delta2 : momentum thickness	 
Lambda=xm*(delta2*xk)**2
WRITE(13,200)beta,xm,H,s,delta1*xk,delta2*xk,s/xk,xk,MaxEta*xk,Lambda,ind_end_BL

results(1:nr)=(/ s, h, delta1*xk, delta2*xk, xk, MaxEta, MaxEta/xk, Lambda /)

200 FORMAT (10(f11.5,2x),3x,i6)
END SUBROUTINE save_characteristics

END MODULE fsk_module

!< @details
!! main test program for Falkner-Skan boundary layer

PROGRAM test_fsk
! read_option = 1 : parameters read in input file
!             = 0 : use default parameters defined in the module
use fsk_module
integer read_option

read_option = 1
call Main_Continuation(read_option)

END PROGRAM test_fsk