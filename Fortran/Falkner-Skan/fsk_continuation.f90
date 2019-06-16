PROGRAM Main_Continuation
!*********** PROGRAMME DE RESOLUTION DE L'EQUATION DE FALKNER-SKAN****
!
!       f'''+ ff''+beta(1-f'^2)=0
!       Attention l'épaisseur de référence est delta = sqrt(nu x / Ue ) * sqrt(2/(m+1))
!       les calculs dans dans le "Cousteix" : "Couche limite laminaire"
!
!       dans notre livre 'Aérodynamique fondamentale" on résout
!                  f'''+(m+1)/2 f f''+m (1-f'^2)=0 avec delta = sqrt(nu x /Ue)
!                  il y a un facteur d'échelle entre les deux. 
!                  Dans le programme ci-dessous  il s'appelle "xk" dans save_characteristics
!       les profils de vitesse sont sortis dans la formulation initiale, celle de Cousteix
!
!                  
! déc 2017 : ajout du RK4 dans tout le programme.
!................................................................
!..     Ce programme intègre l'equation de Falkner-Skan (FSK)
!..     par la methode de continuation.
!..     Entrees: beta parametre de l'eq de FSK
!..              "nom" du fichier d'ecriture de eta, f, f', f''
!                f'=u;  f''=g
!..     Sorties: fichier "sortie3"
!******************************************************************
!
! une méthode de continuation permet de tracer la solution pour
! les valeurs de beta variables
!  Arclength Continuation
!........................... PROGRAMME PRINCIPAL ...................
IMPLICIT NONE

REAL(KIND=8),PARAMETER  :: tolerance_s=1.d-9     ! tolerance for s = f(1)
REAL(KIND=8),PARAMETER  :: tolerance_g=1.d-5     ! tolerance for g = g(infinity)
REAL(KIND=8),PARAMETER  :: CutOff_U=1.d-4        ! tolerance for g = g(infinity) ou t(infinity)
REAL(KIND=8),PARAMETER  :: CutOff_end_BL=1.d-2   ! pour delta_0.99
INTEGER,PARAMETER       :: m= 3                  ! ordre du systeme +1
REAL(KIND=8)            :: deta                  ! step in eta
REAL(KIND=8)            :: eta_max = 10.d0       ! eta_max, par défaut 
INTEGER                 :: imax
INTEGER                 :: size_max              ! size of the vectors
REAL(KIND=8)            :: dL                    ! length of the arc
INTEGER                 :: ind_end_BL
INTEGER                 :: mmax=10,lmax=30       ! maximun size of the loops, valeurs par défaut

REAL(KIND=8)            :: s_init, s0, s1, ds 
REAL(KIND=8)            :: beta0, dbeta,beta1
REAL(KIND=8)            :: L1 
INTEGER                 :: i, k 
INTEGER                 :: cas,flag
INTEGER                 :: flag_continuation=1
LOGICAL                 :: divergence
INTEGER                 :: type_schema=1         ! 0 : ordre 1, 1 : RK 4

 
REAL(KIND=8),DIMENSION(:,:), ALLOCATABLE ::  y
INTEGER,PARAMETER         :: nr=8
REAL(KIND=8),DIMENSION(nr):: results

write(*,110)
110 format(50('*'),/,'* FSK flow, Continuation arclength method', &
          /,50('*'),/)

!  fichiers de sortie généraux
OPEN(14,FORM='FORMATTED',FILE='Convergence.dat')
OPEN(13,FORM='FORMATTED',FILE='BL_characteristics.dat')
WRITE(13,200)
OPEN(20,form='formatted', file = 'output.out' )
write(20,'(a1,6x,a4,10x,a1,10x,a)')'#','beta','s','m'

! valeurs par défaut
CALL read_data(flag)
IF (flag.eq.0) THEN
    ! le fichier n'existant pas on lit les données par défaut
    beta0 = 0.d0
    dbeta =-0.005d0
    deta  = 2.0d-4      ! step in eta
    dL    = 0.005d0     ! length of the arc
    s_init = 0.5d0
END IF
imax=int (eta_max / deta)-1; size_max=imax+1 
write(*,'(a,5x,i7)') 'Size of the vectors',size_max
 
ALLOCATE (y(0:size_max,m))
IF (type_schema.eq.1) THEN
    PRINT*,'Intégration RK 4'
ELSE
    PRINT*,'Intégration : schéma O1 Euler'
END IF

 
! calcul du premier point de convergence de la fonction F(s,beta)=0
!.
! Calcul de M0(beta=0 ; s= s blasius  ou un autre point)
!------------------------------------
WRITE ( * , '(a,/,a)' ) 'calcul de M0', '************'
cas=0
write(*,100) beta0,m_factor(beta0)
CALL Solve_FSK(beta0, s_init, s0) 

WRITE(20,103) beta0,s0,m_factor(beta0)

if (flag_continuation.eq.1) then
    CLOSE(13);CLOSE(14); CLOSE(20)
    DEALLOCATE (y)
    stop 'normal end after a single profile'
end if
! Calcul de M1(beta= beta_init+dbeta, s1)
!----------------------------------------

cas=1
WRITE ( * , '(a,/,a)' ) 'calcul de M1', '************'
beta1=beta0+dbeta
write(*,100) beta1,m_factor(beta1)
CALL Solve_FSK (beta1, s0, s1) 
WRITE(20,103) beta1,s1,m_factor(beta1)

! donnees calculees
!------------------

ds = s1 - s0
L1 = sqrt (dbeta * dbeta + ds * ds) 
WRITE ( *, '(a,e12.6)' ) 'L1 = ', L1
WRITE ( *, * )

!stop 'beta1 : ok'
! Methode de continuation
!------------------------
CALL Continuation (L1, s1, beta1, s0, beta0, divergence) 
CLOSE(20);CLOSE(13);CLOSE(14)
DEALLOCATE (y)
 
WRITE ( * , * ) 'normal end of execution' 

100 FORMAT('# beta =',2x,f10.5,4x,'m =',2x,f10.5)
103 FORMAT (5(e15.8,1x)) 
200 FORMAT('#     beta',10x,'m',13x,'H',12x,'s',10x,'d1/d',9x,'d2/d',&
    9x,'Cf',10x,'xk',9x,'Eta Max',7x,'Lambda',6x,'End_BL')

! ************************************************
! end of main program
! ************************************************

CONTAINS

!**********************
FUNCTION m_factor(beta)
!**********************
    REAL(KIND=8) :: m_factor, beta
    m_factor=beta/(2.d0-beta)
END FUNCTION m_factor

!**********************
FUNCTION beta_factor(mloc)
!**********************
    REAL(KIND=8) :: beta_factor, mloc
    beta_factor=2.d0*mloc/(mloc+1.d0)
END FUNCTION beta_factor

!*************************
SUBROUTINE read_data(flag)
!*************************
IMPLICIT NONE
LOGICAL               :: file_exist
INTEGER, INTENT(OUT)  :: flag

INQUIRE ( file = 'fsk_continuation.in', exist = file_exist )
IF (file_exist) then
    WRITE(*,*) 'Le fichier d''entrée a été trouvé'
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
    WRITE(*,*) 'Le fichier d''entrée n''a pas été trouvé'
    WRITE(*,*) 'Lecture des paramètres par défaut'
    flag=0
END IF

END SUBROUTINE read_data


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

!**************************
SUBROUTINE rk4_cont(beta,dt,tin,y0,dy0,dy0b,yout,dyout,dyoutb)
!**************************
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


 
!************************************************************
SUBROUTINE ode(beta,x,x_d,tloc,dxdt,dxdt_d) 
!************************************************************
! ODE du problème
IMPLICIT NONE
REAL(KIND=8),INTENT(IN)                  :: tloc,beta
REAL(KIND=8),DIMENSION(m),INTENT(IN)     :: x,x_d
REAL(KIND=8),DIMENSION(m),INTENT(OUT)    :: dxdt,dxdt_d
dxdt=0.d0;dxdt_d=0.d0
dxdt(1) = x(2)                      !f' = u
dxdt(2) = x(3)                      !u' = g
dxdt(3) = - x(1)*x(3)-beta *(1.d0-x(2)*x(2))    ! g' = - f f'' -beta (1-f'^2) 
  
dxdt_d(1) = x_d(2)                      !df' = du
dxdt_d(2) = x_d(3)                      !du' = dg
dxdt_d(3) = -x_d(1)*x(3)-x(1)*x_d(3)+ 2.d0 *beta *x_d(2)* x(2)   ! t' = - f df'' - df f'' -2 beta f' df' 
END SUBROUTINE  ode

!************************************************************
SUBROUTINE ode_cont(beta,x,x_d,x_b,tloc,dxdt,dxdt_d,dxdt_b) 
!************************************************************
! ODE du problème
IMPLICIT NONE
REAL(KIND=8),INTENT(IN)                  :: tloc,beta
REAL(KIND=8),DIMENSION(m),INTENT(IN)     :: x,x_d,x_b
REAL(KIND=8),DIMENSION(m),INTENT(OUT)    :: dxdt,dxdt_d,dxdt_b
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


!**********************************
SUBROUTINE Solve_FSK(beta, s_in, s_out) 
!**********************************

IMPLICIT NONE
!     - la dicretisation est realisee par un schema de Runge-Kutta d'ordre 4
!     - la condition a la limite infinie s est determinee par
!       la methode de Newton

INTEGER      :: i, k, n, nmax
REAL(KIND=8) :: test_eta, test_s, beta, s_in, s_out, s 
REAL(KIND=8) :: s_old,eta_loc
REAL(KIND=8),DIMENSION(0:size_max,m) :: dy
REAL(KIND=8),DIMENSION(m) :: grady,gradys
LOGICAL      :: test
!
!      y(:,1)  :  fonction de similutude
!      y(:,2)  : df/d eta
!      y(:,3)  : d2 f / d eta 2
!
!-INITIALISATION DU CALCUL
!-------------------------
!      pour avoir une tres bonne precision, il faut prendre
!      imax tres grand et tolerance_g tres petit

test=.true.
WRITE(14,*)'# beta = ',beta
test_eta = 1.d0; test_s = 1.d0; 
nmax = 100
s = s_in

!- DEBUT BOUCLE SUR n: calcul de s
!--------------------------------
DO n = 1, nmax 
   
    IF (test_s.lt.tolerance_s)  exit
    y=0.d0;dy=0.d0
    y(0,1:m) =(/0.d0,0.d0,s/)
    dy(0,1:m)=(/0.d0,0.d0,1.d0/)
    test_eta = 1.d0
!
!- DEBUT BOUCLE SUR i: resolution simultannee de l'ODE et de son gradient
!--------------------------------------------------------------------
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
              write(*,*) 'tolerance sur u" atteinte pour i = ',i
              write(*,*) 'eta                              = ',eta_loc 
              test=.false.
          end if
          k = k + 1 
      END DO 
  

!-FIN DE LA BOUCLE SUR i
!-----------------------

!           ......CALCUL DE s(n+1)...................................

       s_old = s
       print*,' s= , k ',s,k,i
       s = s - (y(k,2) - 1.d0) / dy(k,2)

!           ......VARIABLE DE CONVERGENCE SUR s.....................

       test_s = abs (s - s_old)
       WRITE(14,200) test_s, s_old
       if (isnan(s)) then
           write(*,*) 's is not a number, s_old=',s_old
           stop 'divergence'
       end if
       if (abs(y(i,2)).gt.10.d0) then
            write(*,*) 'explosion pour i= ',i,'beta = ', beta
            stop 'not a number'
      end if

ENDDO

if (n.eq.nmax) stop 'Divergence'

!-FIN DE LA BOUCLE SUR n
! ----------------------

WRITE (* ,100 )s,test_s,deta*float(i)
s_out = s
CALL save_profile(k,cas,beta)
CALL save_characteristics(k,beta)

100 format('Convergence, s = ',f11.8,3x,'test = ',e12.5,' etamax ',e12.5,/)
200 FORMAT(2(e15.8,3x))

END SUBROUTINE Solve_FSK


!********************************************************************
SUBROUTINE Continuation (L1, s1, beta1, s_init, beta_init,divergence)
!********************************************************************

IMPLICIT NONE
REAL(KIND=8),PARAMETER :: tolerance_f=1.0d-5    ! tolerance on f

INTEGER  i, k, k_beta, j, n, nmax, p
LOGICAL test_L 

REAL(KIND=8),INTENT(IN) :: s_init,beta_init,s1,L1,beta1
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

CHARACTER(15) fichier 

!      f  :  fonction de similutude
!      u  : df/d eta
!      g  : d2 f / d eta 2
!      t  : d3 f / d eta 3

!      pf : p=df/ds
!      pu : dp/d eta
!      pg : d2 p/d eta 2
!      pt : d3 p/d eta 3

!      qf : q=df/dbeta
!      qu : dq/d eta
!      qg : d2 q/d eta 2
!      qt : d3 q/d eta 3


!
!-INITIALISATION DU CALCUL
!-------------------------
divergence=.false.
delta_L  = 0.d0
s_arc    = 0.d0
beta_arc = 0.d0
beta_old = beta_init
s_old    = s_init
beta_new = beta1
s_new    = s1
test_L   = .false.

!- DEBUT BOUCLE SUR L
!====================

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


    !- CALCUL DE q1= df2/dbeta et de p1=df2/ds
    !_______________________________________


    q1 = delta_beta / (delta_L * dL)
    p1 = delta_s / (delta_L * dL)


    !- CALCUL DE df1/dbeta=qu  et de df1/ds=pu
    !----------------------------------------

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

        !- DEBUT BOUCLE SUR i:resolution simultannee des systemes (f),(p)et(q)  
        !--------------------------------------------------------------------   

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
        !           ......VARIABLE DE SORTIE DE LA COUCHE LIMITE.............   
            test_eta = abs(y(i,3)) 
            if (abs(y(i,2)).gt.10.d0) then
                write(*,*) 'explosion pour i= ',i,'beta = ', beta
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

        !-FIN DE LA BOUCLE SUR i
        !-----------------------
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

    !- NOUVEAUX COEFF S ET BETA
    !---------------------------

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

!-FIN DE LA BOUCLE SUR n
! ----------------------

enddo 

!-FIN DE LA BOUCLE SUR L
!=======================

END SUBROUTINE Continuation

!****************************
function charac(N_in,i_enter)
!****************************
!  implemented by C Airiau, march 2010
!  from Burkardt subroutines
!
! conversion of an integer to a character
!  example : 123 --> '123'
! written on Nin digits

integer, intent(in):: i_enter,N_in

character(len=N_in):: charac
integer            :: i,k,n

! n=int(log10(float(i)))
if (N_in.le.1) then
        write(*,*) 'problem ni charact function'
        stop
endif

n=N_in-1
i=i_enter
do k=n,0,-1
        charac(n-k+1:n-k+1)=char(48+int(i/10**k))
        i=mod(i,10**k)
        !     write(*,*) 'k', k,'<<',carac
end do
return 
end function charac

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
    write(1,'(5(e12.5,2x))') deta*real(i,kind=8),y(i,:)
    !if (abs(y(i,2)-1.d0).le.CutOff_U) exit
    if ((.not.test).and.(abs(y(i,2)-1.d0).le.CutOff_end_BL)) then
        ind_end_BL=i; test=.true.
        write(*,*) 'end of Boundaray layer at eta= ',float(ind_end_BL)*deta, &
            'index = ',ind_end_BL
    end if
end do
CLOSE(1)

print*,'fichier des profils écrit'
100 format('# beta = ',f10.5,/,'# ',3x,'eta',12x,'f',12x,'u',12x,'g', &
       12x,'t')
END SUBROUTINE SAVE_PROFILE

SUBROUTINE save_characteristics(k,beta)
IMPLICIT NONE
INTEGER, INTENT(IN)     :: k
REAL(KIND=8),INTENT(IN) :: beta
INTEGER                 :: i
REAL(KIND=8)            :: delta1,delta2,H,xm,MaxEta,xk,s,Lambda

print*,'sauvegarde des données couche limite'
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
!         delta1 : epaisseur de deplacement
!         delta2 : epaisseur de quantite de mouvement	 
Lambda=xm*(delta2*xk)**2
WRITE(13,200)beta,xm,H,s,delta1*xk,delta2*xk,s/xk,xk,MaxEta*xk,Lambda,ind_end_BL

results(1:nr)=(/ s, h, delta1*xk, delta2*xk, xk, MaxEta, MaxEta/xk, Lambda /)

200 FORMAT (10(f11.5,2x),3x,i6)
END SUBROUTINE save_characteristics

END PROGRAM Main_Continuation


