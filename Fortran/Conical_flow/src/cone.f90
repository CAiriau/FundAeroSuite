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
!
!------------------------------------------------------------------------------
! TITLE         : CONICAL FLOW SOLVER
! PROJECT       : FundAeroSuite
! MODULE        : cone
! URL           : 
! AFFILIATION   : Paul Sabatier University, Toulouse
! DATE          : 2016 - 2020
! REVISION      : 2020 V2
!> @author
!> Christophe Airiau
!
! DESCRIPTION:
!>  Solve the shock over a cone, create a database
!------------------------------------------------------------------------------
!!
! wall_angle    : cone angle
!
! shock_angle   : shock angle
! theta         : velocity deviation angle
! M0            : upstream Mach at infinity
!
program  main_conical_flow
!> @details     
!! supersonic conical flow : study of the shock wave
!!
!!    approach inspired by LASSALINE's paper 
!!   
implicit none
real(kind=8),parameter  :: gam=1.4d0, pi = 4.d0*atan(1.d0)
real(kind=8),parameter  :: coef=180.d0/pi, epsil=1.0d-7
integer                 :: option_Mach
integer                 :: opt_Newton
integer                 :: opt_velocity
integer                 :: n_Mach_Max, n_Mach_tmp
real(kind=8),dimension(:),allocatable:: table_Mach
integer                 :: n_Mach, nc
real(kind=8)            :: Mach0, Mach_initial, Mach_final, dMach

real(kind=8)            :: error_max
real(kind=8)            :: delta_Mach
real(kind=8)            :: shock_angle_init
real(kind=8)            :: wall_angle                  
real(kind=8)            :: shock_angle
real(kind=8)            :: Inv_Mach_Cone, Kp, Mach_Cone
real(kind=8)            :: upstream_normal_mach, downstream_normal_mach, downstream_mach
real(kind=8)            :: V0, V1, VP1, VS1
real(kind=8)            :: tetamax
real(kind=8)            :: downstream_mach_1_old
! when the Mach number on the cone is unity Mc = 1
real(kind=8)            :: shock_angle_MC1, wall_angle_MC1, downstream_mach_Mc1
! when the downstream Mach number is  unity
real(kind=8)            :: shock_angle_downstream_M1, wall_angle_downstream_M1, Mach_cone_downstream_M1
real(kind=8)            :: Kp_downstream_M1, Inv_Mach_cone_downstream_M1, Kp_Mc_unitaire
real(kind=8)            :: alpha
real(kind=8)            :: dteta = 1.d-5          ! parameter of accuracy for the wall calculation
real(kind=8)            :: shock_angle_step, Angle_deviation
real(kind=8)            :: shock_angle_old, Inv_Mach_cone_old, downstream_mach_old
real(kind=8)            :: wall_angle_old, Kp_old, Mach_cone_old, Angle_cone_old
integer                 :: flag_M1_unity, flag_Mc_unity, i, k, flag_theta_max
integer                 :: n_points = 201           ! number of point on the shock angle curve, wall angle theta
character(len=20)       :: filename, filename1

call read_data
call set_Mach_table
call open_output_files

n_Mach = n_Mach_Max

! initialisation 
shock_angle_old = 0.d0
Inv_Mach_cone_old = 10.d0
Kp_old = 0.d0
Mach_cone_old = 10.d0
wall_angle_old = shock_angle_old
downstream_mach_old = 0.d0

do i = 1, n_Mach
    Mach0 = table_Mach(i)
    filename = 'Mach_'//carac(5,int(1000.0*(Mach0+epsil)))//'.dat'
    filename1 = 'velocity_'//carac(5,int(1000.0*(Mach0+epsil)))//'.dat'
    open(3, form='formatted', file=trim(filename1))
    open(2, form='formatted', file=trim(filename)); write(2,120); print *,'Mach0 = ', Mach0
    flag_M1_unity = 1; flag_Mc_unity = 1;flag_theta_max = 1
    shock_angle_init = asin(1.d0/Mach0)
    write(2,100)0.d0, shock_angle_init*coef, 0.d0, 1.d0-1.d0/Mach0, Mach0   ! wrong value, to change in the output file
    write(*,'(a40,f15.8)')  'angle minimal shock angle   ', shock_angle_init*coef
    wall_angle = 0.d0                     ! initial cone angle 
    write(80,100)0.0, shock_angle_init*coef, Mach0
    shock_angle = shock_angle_init                              ! initial shock angle 
    shock_angle_step = (pi/2-shock_angle)/dfloat(n_points-1)    ! n_points for the shock curve
    ! loop on the shock angle of   arcsing( 1/Mach_upstream)  to  pi/2 
    downstream_mach_1_old = 10.d0                               ! initialisation at each loop
    Angle_cone_old = 0.20d0*shock_angle
    wall_angle_old = 0.0d0
    shock_angle_old = shock_angle
    
    do k = 1, n_points-2
        shock_angle = shock_angle + shock_angle_step*correction(k)  ! next shock angle
        Angle_deviation = deviation(Mach0,shock_angle)
        upstream_normal_mach = Mach0*sin(shock_angle)
        downstream_normal_mach = eval_downstream_normal_mach(upstream_normal_mach)
        downstream_mach = downstream_normal_mach/sin(shock_angle-Angle_deviation)

        !  V components in spherical coordinates 
        ! V1 = V_r, VP1 = V'_r=V_theta, VS1 = V''_r
        V0 = 1.d0/(sqrt(1.d0+2.d0/((gam-1.d0)*downstream_mach**2)))   ! V1/VL 
                                                                ! VL : reduced velocity
        V1 = V0*cos(shock_angle-Angle_deviation)                ! reduced radial velocity Vr
        VP1 = -V0*sin(shock_angle-Angle_deviation)              ! derivative  of reduced radial velocity V'r = V_theta
        VS1 = calcul_VS1(shock_angle,V1,VP1)                    ! initial value before ODE integration   of 
        !                                                        from theta= Sigma to theta = wall_angle where V_theta=0
        ! we have the initial value of the shock angle 
        ! loop where angle decrease from the shock to the cone wall 
        ! radial velocity integration to solve the problem 
        ! order 4 integration

        select case (opt_Newton)
        case(0)
            call integration_O1(shock_angle, wall_angle, V1, VP1, VS1)
        case(1)
            call solve_angle_cone2(shock_angle, wall_angle, V1, VP1, VS1)
            !write(51,*)wall_angle*coef,shock_angle*coef
        case(2)
            wall_angle = Angle_cone_old
            print*,'Newton : order 4'
            call solve_angle_cone(shock_angle, wall_angle, V1, VP1, VS1)
            !write(51,*)wall_angle*coef,shock_angle*coef
            print*, wall_angle*coef, shock_angle*coef
        end select
        
        ! calculus of   1-1/Mc
        Mach_Cone = sqrt(2.d0/(gam-1.d0)*V1*V1/(1.d0-V1*V1))
        Inv_Mach_Cone = 1.-1./Mach_Cone                          ! cone surface parameter
        Kp = pressure_coefficient(Mach0, Mach_cone, upstream_normal_mach)

        write(85,100) wall_angle*coef, Inv_Mach_Cone, Mach0
        write(80,100) wall_angle*coef, shock_angle*coef, Mach0
        write(90,100) wall_angle*coef, Kp, Mach0

        ! maximal deviation and loop break 	
        if ((wall_angle-wall_angle_old .lt. 0.d0) .and. (flag_theta_max .eq. 1)) then
            flag_theta_max = 0; tetamax = wall_angle
            write(81,100) Mach0, shock_angle*coef, tetamax*coef, Inv_Mach_Cone, Kp
        endif
    
        ! locus of  unitary downstream Mach
        !===========================
        ! unitary  downstream  Mach 
        !===========================

        if ((downstream_mach .le. 1.d0) .and. (flag_M1_unity .eq. 1)) then
            alpha = -(downstream_mach_old-1.d0)/(downstream_mach-downstream_mach_old)
            flag_M1_unity = 0
            Mach_cone_downstream_M1 = alpha*Mach_cone + (1.d0-alpha)*Mach_cone_old
            wall_angle_downstream_M1 = alpha*wall_angle + (1.d0-alpha)*wall_angle_old
            shock_angle_downstream_M1 = alpha*shock_angle + (1.d0-alpha)*shock_angle_old
            Kp_downstream_M1 = alpha*Kp + (1.d0-alpha)*Kp_old
            Inv_Mach_Cone_downstream_M1 = 1.d0-1.d0/Mach_cone_downstream_M1
            write(*,'(a,2(e12.5,2x))')'error on downstream Mach 1 ',&
                        abs(1.d0-wall_angle_downstream_M1/wall_angle), & 
                        abs(1.d0-Mach_cone_downstream_M1/Mach_cone) 
            write(82,100) Mach0, shock_angle_downstream_M1*coef, wall_angle_downstream_M1*coef,&
                            Inv_Mach_Cone_downstream_M1, Kp_downstream_M1
        endif

        !====================
        ! Mach cone unitary
        !====================

        if ((Mach_Cone .le. 1.d0) .and. (flag_Mc_unity .eq. 1)) then
            alpha = -inv_Mach_cone_old / (inv_Mach_cone-inv_Mach_cone_old)
            shock_angle_Mc1 = alpha*shock_angle + (1.d0-alpha)*shock_angle_old
            wall_angle_Mc1 = alpha*wall_angle + (1.d0-alpha)*wall_angle_old
            downstream_mach_Mc1 = alpha*downstream_mach + (1.d0-alpha)*downstream_mach_old
            Kp_Mc_unitaire = alpha*Kp + (1.d0-alpha)*Kp_old
            print*, 'error on  Mc1 ', abs(1.d0-wall_angle_Mc1/wall_angle)
            write(83,100) table_Mach(i), shock_angle_Mc1*coef, wall_angle_Mc1*coef, downstream_mach_Mc1, Kp_Mc_unitaire
            flag_Mc_unity = 0
        end if

        if (flag_theta_max .ne. 0) then
            write(2,100) wall_angle*coef, shock_angle*coef, Kp, Inv_Mach_Cone, Mach0, downstream_mach
        end if
        Inv_Mach_Cone_old = Inv_Mach_Cone
        shock_angle_old = shock_angle
        Kp_old = Kp
        wall_angle_old = wall_angle
        Mach_cone_old = Mach_Cone
        downstream_mach_old = downstream_mach

    end do

    write(80,100)0.0d0, 90.d0, Mach0
    write(80,*)'#'; write(85,*)'#'
    ! Kp=-2./gam/Mach0/Mach0             ! pressure  fonction de P_downstream/P_upstream
    ! non physical value for the last point:
    ! write(2,100)0.d0,90.d0,Kp, 1.d0-1.d0/mach1(Mach0),Mach0   ! wrong value ?
    close(2); close(3)
end do
close(90); close(80); close(81); close(82)
print*,'normal end of execution'
deallocate(table_Mach)
100 format(10(f15.8,3x))
120 format('#',2x,'Theta wall in deg.',2x,'shock angle in deg. ',2x,'Kp',11x,'Inverse Mach',7x,'upstream Mach',&
       3x,'downstream Mach')

contains


!==================
! RK 4 Integration
!==================


!> @details
!!
!!  ODE system to solve with Runge-Kutta scheme
!!
!! t  : theta,  y =(U_r, U_theta=dU_r/dtheta), 
!! y(1)= U_r,  y(2)= d U_r / dtheta = U_theta
!!
!! system to solve
!!
!! d y1/dt = y2
!!
!! d y2/dt = (y2**2*y1-(gam-1)/2(1-y1**2-y2**2)*(2*y1+y2/tant(theta))/ 
!!          ( (gam-1)/2*(1-y1**2-y2**2)-y2**2
!!
!
!******************************
subroutine rhs(theta, y, dydt)
!******************************
    implicit none
    real(kind=8),intent(in)                 :: theta  !< deviation angle
    real(kind=8),dimension(2), intent(in)  :: y       !< state
    real(kind=8),dimension(2), intent(out) :: dydt    !< state derivative 
    real(kind=8)                           :: num,den,y1,y2
    y1 = y(1)
    y2 = y(2)
    num = -(gam-1.d0)/2.d0*(1.d0-y1**2-y2**2)*(2.d0*y1+y2/(tan(theta)))+y1*y2**2
    den = (gam-1.d0)/2.d0*(1.d0-y1**2-y2**2)-y2**2
    dydt(1) = y2
    dydt(2) = num/den
end subroutine rhs

!> @details 
!! order 4 Runge Kutta algorithm 
!! 
!! t : theta, y(1)=V_r, y(2)=U_theta=d V_r/d theta
!! 
!******************************
subroutine rk4(dt, t, y, dydt)
!******************************
    implicit none
    real(kind=8),intent(in)        :: dt !< time step
    real(kind=8),intent(in)        :: t  !< time
    real(kind=8),dimension(2), intent(inout) :: y       !< state
    real(kind=8),dimension(2), intent(out) :: dydt      !< state derivative 
    real(kind=8),dimension(2)      :: k1,k2,k3,k4,y_tmp

    ! First sub-step
    call rhs(t,y,dydt)
    k1= dt*dydt; y_tmp=y+0.5d0*k1
    ! Second sub-step
    call rhs(t,y_tmp,dydt)
    k2= dt*dydt; y_tmp=y+0.5d0*k2
    ! Third sub-step
    call rhs(t,y_tmp,dydt)
    k3= dt*dydt; y_tmp=y+k3
    ! Fourth sub-step
    call rhs(t,y_tmp,dydt)
    k4= dt*dydt; y=y+(k1+2.0d0*(k2+k3)+k4)/6.0d0
    !call rhs(t,y,dydt)
end subroutine rk4

!> @details
!! solve tangential velocity
!************************************************************
subroutine U_theta(sauver, theta_init, thetac, y0, y, dy_dt)
!***********************************************************
    implicit none
    logical,intent(in)                     :: sauver   !< option to save in a file
    real(kind=8),dimension(2),intent(in)   :: y0       !< initial condition
    real(kind=8),intent(in)                :: theta_init !< initial condition for the deviation angle theta 
    real(kind=8),intent(in)                :: thetac     !< cone angle
    real(kind=8),dimension(2),intent(inout):: y          !< state y(:,1)=U_r, y(:,2)=U_theta
    real(kind=8),dimension(2),intent(inout):: dy_dt      !< state derivative

    integer,parameter                      :: nc_max=1000
    real(kind=8),dimension(nc_max)         :: t
    integer,parameter                      :: iter_max=30
    real(kind=8)                           :: dt
    integer                                :: i

    if (sauver) write(3,110)k, shock_angle*coef
    dy_dt = 0.d0
    ! on can modify the integration step
    ! nc=int((abs(thetac-theta_init)/dteta))+1
    ! if (nc.gt.nc_max) nc=nc_max
    dt = (thetac-theta_init)/real(nc-1,kind=8)
    print*,'U_theta :     dt =', dt*coef, '  nc = ',nc
    do i =1, nc-1
        t(i) = theta_init + dt*real(i-1,kind=8)
    end do
    t(nc) = thetac
    y = y0
    do i = 1, nc
        call rk4(dt, t(i), y, dy_dt)
        if (sauver .and. (opt_velocity .eq. 1)) then
            write(3,100) t(i)*coef, y(1), y(2), dy_dt(2), shock_angle*coef,Mach0
        end if
    end do
    100 format(6(e20.13,3x))
    110 format('# theta',i4,/,'#',7x,f12.4,/,'#')
end subroutine U_theta

!> @details 
!! to calculate the cone angle, last point in the loop process 
!!   
!! Newton method    + RK4
!****************************************************************************
subroutine solve_angle_cone(theta_init, theta_cone, Ur, dUr_dtheta, dU_theta)
!****************************************************************************

    implicit none
    real(kind=8),intent(inout)  :: theta_init                
    real(kind=8),intent(inout)  :: Ur,dUr_dtheta,dU_theta    
    real(kind=8),intent(inout)  :: theta_cone                
    integer                     :: iter
    real(kind=8)                :: error,delta_theta
    integer,parameter           :: iter_max=20
    real(kind=8),dimension(2)   :: y0,y,dy_dt
    print*,'solve_angle_cone, error_max =', error_max
    y0(1) = Ur 
    y0(2) = dUr_dtheta
    error = 1.d0
    iter = 0
    !theta_cone=0.25d0*theta_init
    if (theta_cone.le.0.d0) then 
        print*, 'theta_cone ', theta_cone*coef, 'theta_init = ', theta_init*coef
        stop 'error theta_cone'
    end if
    do while ((error .gt. error_max) .and. (iter.le. iter_max)) 
        iter=iter+1
        print*,'Newton : iteration = ', iter, 'theta c = ', theta_cone*coef
        call U_theta(.false., theta_init, theta_cone, y0, y, dy_dt)
        error = abs(y(2))
        delta_theta = -y(2)/dy_dt(2)
        write(50,*) iter, theta_cone*coef, error, delta_theta*coef
        theta_cone = theta_cone+delta_theta
        if (theta_cone .le. 0.d0) then 
            print*, 'theta_cone ', theta_cone*coef, 'theta_init = ', theta_init*coef
            stop 'errorr theta_cone loop'
        end if
    end do

    if (iter .eq. iter_max) then
        stop 'no convergence in Newton_iteration procedure'
    else
        print*,'convergence in Newton'
    end if
    call U_theta(.true., theta_init, theta_cone, y0, y, dy_dt)
    Ur = y(1)
    dUr_dtheta = y(2)
    dU_theta = dy_dt(2)
    write(50,*) 'Cone angle in deg ', theta_cone*coef, ' with ', iter, ' iterations'

end subroutine solve_angle_cone

!> details
!! to calculate the cone angle, last point in the loop process 
!!
!! old version
!! Newton method
!************************************************************************
subroutine solve_angle_cone2(theta_init,theta_cone,Ur,dUr_dtheta,dU_theta)
!************************************************************************

    implicit none
    real(kind=8),intent(inout)  :: theta_init
    real(kind=8),intent(inout)  :: Ur, dUr_dtheta, dU_theta
    real(kind=8),intent(inout)  :: theta_cone                 !< initial and final angle
    integer                     :: iter
    real(kind=8)                :: error,teta,teta_old,alpha
    integer,parameter           :: iter_max=20
    real(kind=8),dimension(2)   :: y0,y,dy_dt,y_old,dy_dt_old

    y0(1) = Ur
    y0(2) = dUr_dtheta
    error = 1.d0
    iter = 0
    write(3,110)k, shock_angle*coef
    teta = theta_init  
    teta_old = teta
    y_old = y0
    dy_dt_old = 0.d0
    y = y0
    do while (teta .ge. 0.d0)
        call rk4(-dteta, teta, y, dy_dt)
        if (opt_velocity .eq. 1)  write(3,100) teta*coef, y(1), y(2), dy_dt(2), shock_angle*coef,Mach0
        if (y(2) .gt. 0.d0) then                           ! V_theta is <0 and = 0 @ wall
            alpha = -y_old(2)/(y(2)-y_old(2))
            y = alpha*y+(1.d0-alpha)*y_old
            dy_dt = alpha*dy_dt+(1.d0-alpha)*dy_dt_old
            theta_cone = alpha*teta+(1.d0-alpha)*teta_old
            write(3,100) theta_cone*coef, y(1), y(2), dy_dt(2), shock_angle*coef,Mach0
            !wall_angle=teta                           ! the cone is reached
            exit
        endif
        y_old = y
        dy_dt_old = dy_dt
        teta_old = teta
        teta = teta-dteta
    end do
    Ur = y(1)
    dUr_dtheta = y(2)
    dU_theta = dy_dt(2)

    100 format(6(e20.13,3x))
    110 format('# theta',i4,/,'#',7x,f12.4,/,'#')
end subroutine solve_angle_cone2

include 'sub_cone.f90'

end program main_conical_flow
