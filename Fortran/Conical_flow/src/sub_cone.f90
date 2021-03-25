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

!*********************
function correction(k)
!*********************
!> brief use to calculate a new shock angle

    real(kind=8)    :: correction
    integer         :: k
    integer,parameter   :: k_lim=10
    real(kind=8),parameter :: coef_init=0.03d0

    if (k.gt.k_lim) then
        correction = 1.d0
    else
        correction = coef_init+(1.d0-coef_init)/(real(k_lim-1,kind=8))**2*(real(k-1,kind=8))**2
    end if

end function correction

!> brief
!! to read the input parameters in the input file *cone.in*
!*******************
subroutine read_data
!*******************
    open(1, form='formatted', file='cone.in')
    read(1, *)
    read(1, *) option_Mach
    read(1, *) n_Mach_tmp
    read(1, *) delta_Mach
    read(1, *) Mach_initial
    read(1, *) Mach_final
    read(1, *) n_points
    read(1, *) dteta
    read(1, *) opt_Newton
    read(1, *) nc
    read(1, *) error_max
    read(1, *) opt_velocity
    dteta = dteta/coef
    close(1)
end subroutine read_data

!***************************
subroutine open_output_files
!***************************
!> @brief open all the output file used to save data

    open(80, form='formatted', file='deviation.out')
    write(80, 110)
    110 format('#', 2x, 'theta wall in  °', 2x, 'shock angle in °', &
                5x, 'upstream Mach')

    open(81, form='formatted', file='theta_cone_max.out')
    write(81, 100)
    100 format('#', 5x, 'upstream Mach', 2x,'shock angle in °', 4x, 'thetaMax in °', &
                2x, 'Inverse cone Mach', 8x, 'Kp')

    open(82, form='formatted', file='downstream_unitary_mach.out')
    write(82,120)
    120 format('#', 2x, 'upstream Mach ', 3x, 'shock angle in °',2x,&
            'Theta Mach = 1 in °',2x,'Inverse Mach',9x,'Kp')

    open(83, form='formatted', file='unitary_cone_Mach.out')
    write(83, 150)
    150 format('#', 2x, 'upstream Mach', 6x, 'shock angle in °', 2x,&
            'Theta Mc = 1 in °', 2x, 'downstream Mach', 9x, 'Kp')

    open(90, form='formatted', file='Kp.out')
    write(90, 130)
    130 format('#', 2x,'Theta wall in °', 9x, 'Kp', 11x, 'upstream Mach')

    open(85, form='formatted', file='Inverse_downstream_mach.out')
    write(85, 140)
    140 format('#', 2x, 'Theta wall in °',2x,&
                        'Inverse Mach wall',2x,'upstream Mach')

end subroutine open_output_files 

!************************
subroutine set_Mach_table
!************************
!> @brief to fix the Mach range for the simulations

    select case (option_Mach)
    case (1)
        n_Mach_Max = 25
        allocate(table_Mach(N_Mach_Max))
        table_Mach = (/ 1.01d0, 1.05d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, &
                        1.7d0,  1.8d0, 1.9d0, 2.d0,  2.2d0, 2.4d0, 2.6d0, 2.8d0, &
                        3.0d0,  3.2d0, 3.4d0, 3.6d0, 3.8d0, 4.0d0, 4.5d0, 5.0d0, &
                        10.d0,  50.d0 /)
    case DEFAULT
        if (Mach_initial .eq. 0.d0) Mach_initial = 2.05d0
        if (Mach_final .eq. 0.d0)   Mach_final = 3.5d0
        if (n_Mach_tmp .ne. 0) then
            n_Mach_Max = n_Mach_tmp
            dMach = (Mach_final-Mach_initial)/dfloat(n_Mach_Max-1)
        else
            dMach = delta_Mach
            write(*, *)((Mach_final-Mach_initial)/dMach)
            n_Mach_Max = int((Mach_final-Mach_initial)/dMach)+2
        end if
        allocate(table_Mach(N_Mach_Max))
        table_Mach(1) = Mach_initial
        do i = 2, n_Mach_Max
            table_Mach(i) = table_Mach(i-1)+dMach
        end do
    end select
    write(*,*) 'number of Mach : ', n_Mach_Max
    write(*, '(f10.6)') table_Mach
    write(*, '(i5)') int(1000*table_Mach)
end subroutine set_Mach_table

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
            write(*,*) 'problem in  carac function'
            stop
    endif

    n = n_in-1
    i = i_enter
    do k = n, 0, -1
            carac(n-k+1:n-k+1) = char(48+int(i/10**k))
            i = mod(i,10**k)
    end do
end function carac

!********************
function Rpchoc(mach)
!********************
!> @brief 
!! pressure jump across a normal plane shock    
    real (kind=8)::  Rpchoc,mach
    Rpchoc= 2*gam/(gam+1)* mach**2- (gam-1)/(gam+1)
end function Rpchoc

!********************
function mach1(mach)
!********************
!> brief
!!    Mach normal downstream a normal plane shock
    real (kind=8):: mach,mach1
    mach1= sqrt((1+ 0.5*(gam-1)* mach**2)/(gam*mach**2-0.5*(gam-1)))
end function mach1

!*************************************************************
subroutine integration_O1(teta_init, wall_angle, V1, VP1, VS1)
!*************************************************************
!> @brief order 1 Euler scheme to integrate the ODE 
    implicit none
    real(kind=8)        :: teta_init, wall_angle
    real(kind=8)        :: V1, VP1, VS1
    real(kind=8)        :: teta, dteta_tmp, teta_old, VP1_old, VS1_old, V1_old, alpha

    V1_old = V1; VP1_old = VP1; VS1_old = VS1
    teta = teta_init
    teta_old = teta
    write(3,110) k, shock_angle*coef
    dteta_tmp = dteta
    do  while (teta .ge. 0.d0)
        if (VP1 .gt. -0.05d0) dteta_tmp = dteta*0.1d0
        VP1 = VP1-VS1*dteta_tmp
        V1 = V1-VP1*dteta_tmp
        VS1 = calcul_VS1(teta,V1,VP1)
        write(3,100) teta*coef, V1, VP1, VS1, shock_angle*coef, Mach0
        if (VP1 .gt. 0.d0) then                           ! V_theta  <0 and = 0 @ wall
            alpha = -VP1_old/(VP1-VP1_old)
            wall_angle = alpha*teta+(1.d0-alpha)*teta_old
            V1 = alpha*V1+(1.d0-alpha)*V1_old
            VS1 = alpha*VS1+(1.d0-alpha)*VS1_old
            if (opt_velocity .eq. 1) write(3,100) wall_angle*coef, V1, 0.d0, VS1, shock_angle*coef, Mach0
            write(48,*) (wall_angle-teta)/teta
            !wall_angle=teta                           ! the cone is reached
            exit
        endif
        teta_old = teta
        V1_old = V1; VP1_old = VP1; VS1_old = VS1
        teta = teta-dteta_tmp
    end do
    100 format(6(e20.13, 3x))
    110 format('# theta', i4, /, '#', 7x, f12.4, /, '#')
end subroutine integration_O1

!*******************************
function calcul_VS1(teta,V1,VP1)
!*******************************
!> @brief a part of the unknown in the ODE system, seconde derivative of Vr
    implicit none
    real(kind=8), intent(in)        :: teta, V1, VP1
    real(kind=8)                    :: calcul_VS1
    real(kind=8)                    :: num, den
    num = (gam-1.d0)/2.d0*(1.d0-V1**2-VP1**2)*(2.d0*V1+VP1/(tan(teta))) -V1*VP1**2
    den = (gam-1.d0)/2.d0*(1.d0-V1**2-VP1**2)-VP1**2
    calcul_VS1 = -num/den            ! second derivative w.r.t. theta of  Vr
end function calcul_VS1

!***********************************
function deviation(Mach0,shock_angle)
!***********************************
!> @brief deviation angle behind  oblique shock
    implicit none
    real(kind=8),intent(in)         :: Mach0,shock_angle
    real(kind=8)                    :: deviation
    real(kind=8)                    :: num, den
    num = 2.d0/tan(shock_angle)*((Mach0*sin(shock_angle))**2-1.d0)
    den = Mach0**2.*(gam+cos(2.*shock_angle))+2.d0
    deviation = atan(num/den)
end function deviation

!**************************************************
function eval_downstream_normal_mach(upstream_normal_mach)
!**************************************************
!> @brief downstream normal Mach behind an oblique shock
    implicit none
    real(kind=8),intent(in)     :: upstream_normal_mach
    real(kind=8)                :: eval_downstream_normal_mach
    real(kind=8)                :: num, den
    num = upstream_normal_mach**2.+2.d0/(gam-1.d0)
    den = 2.d0*gam/(gam-1.d0)*upstream_normal_mach**2.-1.d0
    eval_downstream_normal_mach = sqrt(num/den)
end function eval_downstream_normal_mach

!****************************************************************
function pressure_coefficient(Mach0, Mach_cone, upstream_normal_mach)
!****************************************************************
!> @brief wall pressure coefficient on the cone
    implicit none
    real(kind=8),intent(in) :: Mach0, Mach_cone, upstream_normal_mach
    real(kind=8)            :: pressure_coefficient
    real(kind=8)            :: a, b, c, d, rpi

a = 1.+(gam-1.d0)/2.d0*Mach_Cone**2                                         ! Ti/T  downstream
b = 1.+(gam-1.d0)/2.d0*Mach0**2                                             ! Ti/T  upstream
c = 2.d0*gam/(gam+1.d0)*upstream_normal_mach**2-(gam-1.d0)/(gam+1.d0)       ! P_downstream/P_upstream
d = ((gam-1)*upstream_normal_mach**2+2.0)/((gam+1)*upstream_normal_mach**2) ! rho_upstream/rho_downstream
rpi = (a/b*d)**(-gam/(gam-1.d0))*c**(-1.d0/(gam-1.d0)) 
!     rpi :  P downstream/ P upstream= P_downstream/Pi_downstream* Pi downstream/Pi_upstream*Pi_upstream/P_upstream
pressure_coefficient = 2./gam/Mach0/Mach0*(rpi-1.d0)                        ! Kp as a function of  P_downstream/P_upstream
end function pressure_coefficient

