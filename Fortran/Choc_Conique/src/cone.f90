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

program  ecoulement_conique
!    ecoulement  supersonique conique, onde de chocs
!    ceci est inspire du papier de LASSALINE 
!    sur l angle de choc  le Mach sur le cone et la pression
!    sont constantes

! Angle_paroi  : angle du cone
! Angle_choc  : angle de choc
! theta  : déviation de la vitesse
! M0     : Mach amont
!
implicit none
real(kind=8),parameter  :: gam=1.4d0,pi=4.d0*atan(1.d0)
real(kind=8),parameter  :: coef=180.d0/pi,epsil=1.0d-7
integer                 :: option_Mach
integer                 :: opt_Newton
integer                 :: opt_velocity
integer                 :: n_Mach_Max,n_Mach_tmp
real(kind=8),dimension(:),allocatable::table_Mach
integer                 :: n_Mach,nc
real(kind=8)            :: Mach0,Mach_initial,Mach_final,dMach

real(kind=8)            :: error_max
real(kind=8)            :: delta_Mach
real(kind=8)            :: Angle_choc_init
real(kind=8)            :: Angle_Paroi                  ! angle de la paroi (du cone)
real(kind=8)            :: Angle_choc                   ! angle du choc
real(kind=8)            :: Inv_Mach_Cone,Kp,Mach_Cone
real(kind=8)            :: Mach_Normal_Amont,Mach_Normal_Aval,Mach_Aval
real(kind=8)            :: V0,V1,VP1,VS1
real(kind=8)            :: tetamax
real(kind=8)            :: Mach_aval_1_old
! pour le Mach cone unitaire
real(kind=8)            :: Angle_choc_MC1,Angle_paroi_MC1,Mach_aval_Mc1
! pour le Mach aval unitaire
real(kind=8)            :: Angle_choc_aval_M1,Angle_paroi_aval_M1,Mach_cone_aval_M1
real(kind=8)            :: Kp_aval_M1,Inv_Mach_cone_Aval_M1,Kp_Mc_unitaire
real(kind=8)            :: alpha
real(kind=8)            :: dteta=1.d-5          ! paramètre de précision sur theta paroi
real(kind=8)            :: Angle_choc_step,Angle_deviation
real(kind=8)            :: Angle_choc_old,Inv_Mach_cone_old,Mach_aval_old
real(kind=8)            :: Angle_paroi_old,Kp_old,Mach_cone_old,Angle_cone_old
integer                 :: flag_M1_unite, flag_Mc_unite, i,k,flag_theta_max
integer                 :: n_points=201           ! nombre de points sur la courbe Angle Choc, theta paroi
character(len=20)       :: filename,filename1

call read_data
call set_Mach_table
call open_output_files

n_Mach=n_Mach_Max

! initialisation pas trop stupide
Angle_choc_old=0.d0
Inv_Mach_cone_old=10.d0
Kp_old=0.d0
Mach_cone_old=10.d0
Angle_paroi_old=Angle_choc_old
Mach_aval_old=0.d0

do i=1,n_Mach
    Mach0=table_Mach(i)
    filename='Mach_'//carac(5,int(1000.0*(Mach0+epsil)))//'.dat'
    filename1='velocity_'//carac(5,int(1000.0*(Mach0+epsil)))//'.dat'
    open(3,form='formatted',file=trim(filename1))
    open(2,form='formatted',file=trim(filename)); write(2,120); print *,'Mach0=',Mach0
    flag_M1_unite=1; flag_Mc_unite=1;flag_theta_max=1
    Angle_choc_init=asin(1.d0/Mach0)
    write(2,100)0.d0,Angle_choc_init*coef,0.d0,1.d0-1.d0/Mach0,Mach0   ! pas les bonnes valeurs
    write(*,'(a40,f15.8)')  'angle minimal de choc ',Angle_choc_init*coef
    Angle_Paroi=0.d0                     ! angle du cone initial
    write(80,100)0.0,Angle_choc_init*coef,Mach0
    Angle_choc=Angle_choc_init                              ! angle de choc initial
    Angle_choc_step=(pi/2-Angle_choc)/dfloat(n_points-1)    ! n_points sur la courbe de choc
    ! boucle sur l'angle de choc  de arcsing( 1/Mach_amont) à pi/2 
    Mach_aval_1_old=10.d0                                  ! initialisation à chaque boucle
    Angle_cone_old=0.20d0*Angle_choc
    Angle_paroi_old=0.0d0
    Angle_choc_old=Angle_choc
    
    do k=1,n_points-2
        Angle_choc=Angle_choc + Angle_choc_step*correction(k)  ! prochain angle de choc
        Angle_deviation=deviation(Mach0,Angle_choc)
        Mach_Normal_Amont=Mach0*sin(Angle_choc)
        Mach_Normal_Aval=calcul_Mach_Normal_Aval(Mach_Normal_Amont)
        Mach_Aval=Mach_Normal_Aval/sin(Angle_choc-Angle_deviation)

        ! calcul des composantes de V en coordonnees spheriques
        ! V1 = V_r, VP1 = V'_r=V_theta, VS1 = V''_r
        V0=1.d0/(sqrt(1.d0+2.d0/((gam-1.d0)*Mach_Aval**2)))   ! V1/VL 
                                                              ! VL : vitesse réduite
        V1=V0*cos(Angle_choc-Angle_deviation)                 ! vitesse réduite radiale Vr
        VP1=-V0*sin(Angle_choc-Angle_deviation)               ! dérivée de la vitesse  réduite radiale V'r = V_theta

        VS1=calcul_VS1(Angle_choc,V1,VP1)                     ! valeur initiale de l'ODE à intégrer
        !                                                     de theta= Sigma à theta = Angle_paroi où V_theta=0
        ! on a les valeurs initiales au choc 
        ! boucle du choc vers le cone
        ! intégration de la vitesse radiale pour résoudre le problème
        ! schéma à l'ordre 1

        select case (opt_Newton)
        case(0)
            call integration_O1(Angle_choc,Angle_Paroi,V1,VP1,VS1)
        case(1)
            call solve_angle_cone2(Angle_choc,Angle_Paroi,V1,VP1,VS1)
            !write(51,*)Angle_Paroi*coef,Angle_choc*coef
        case(2)
            Angle_paroi=Angle_cone_old
            print*,'Newton : ordre 4'
            call solve_angle_cone(Angle_choc,Angle_Paroi,V1,VP1,VS1)
            !write(51,*)Angle_Paroi*coef,Angle_choc*coef
            print*,Angle_Paroi*coef,Angle_choc*coef
        end select
        
        ! calcul  1-1/Mc
        Mach_Cone=sqrt(2.d0/(gam-1.d0)*V1*V1/(1.d0-V1*V1))
        Inv_Mach_Cone=1.-1./Mach_Cone                                           ! parametre de surface du cone
        Kp=coefficient_pression(Mach0,Mach_cone,Mach_normal_Amont)

        write(85,100) Angle_Paroi*coef,Inv_Mach_Cone,Mach0
        write(80,100) Angle_Paroi*coef,Angle_choc*coef,Mach0
        write(90,100) Angle_Paroi*coef,Kp,Mach0

        ! déviation maximale et arret des calculs 	
        if ((Angle_Paroi-Angle_Paroi_old.lt.0.d0).and.(flag_theta_max.eq.1)) then
            flag_theta_max=0; tetamax=Angle_Paroi
            write(81,100) Mach0,Angle_choc*coef,tetamax*coef,Inv_Mach_Cone,Kp
        endif
       
        ! lieu du Mach aval unitaire
        !====================
        ! Mach aval unitaire
        !====================
        if ((Mach_aval.le.1.d0).and.(flag_M1_unite.eq.1)) then
            alpha=-(Mach_aval_old-1.d0)/(Mach_aval-Mach_aval_old)
            flag_M1_unite=0
            Mach_cone_aval_M1   = alpha*Mach_cone  +(1.d0-alpha)*Mach_cone_old
            Angle_paroi_aval_M1 = alpha*Angle_paroi+(1.d0-alpha)*Angle_paroi_old
            Angle_choc_aval_M1  = alpha*Angle_choc +(1.d0-alpha)*Angle_choc_old
            Kp_aval_M1          = alpha*Kp         +(1.d0-alpha)*Kp_old
            Inv_Mach_Cone_aval_M1 = 1.d0-1.d0/Mach_cone_aval_M1
            write(*,'(a,2(e12.5,2x))')'erreur sur Mach aval 1 ',&
                     abs(1.d0-Angle_paroi_aval_M1/Angle_paroi), & 
                     abs(1.d0-Mach_cone_aval_M1/Mach_cone) 
            write(82,100) Mach0,Angle_choc_aval_M1*coef,Angle_paroi_aval_M1*coef,&
                          Inv_Mach_Cone_aval_M1,Kp_aval_M1
        endif

        !====================
        ! Mach cone unitaire
        !====================

        if ((Mach_Cone.le.1.d0).and.(flag_Mc_unite.eq.1)) then
            ! il faudrait faire une interpolation linéaire
           alpha=-inv_Mach_cone_old/(inv_Mach_cone-inv_Mach_cone_old)
           Angle_choc_Mc1 =alpha*Angle_choc +(1.d0-alpha)*Angle_choc_old
           Angle_paroi_Mc1=alpha*Angle_paroi+(1.d0-alpha)*Angle_paroi_old
           Mach_aval_Mc1  =alpha*Mach_aval  +(1.d0-alpha)*Mach_aval_old
           Kp_Mc_unitaire = alpha*Kp        +(1.d0-alpha)*Kp_old
           print*,'erreur sur Mc1 ',abs(1.d0-Angle_paroi_Mc1/Angle_paroi)
           write(83,100) table_Mach(i),Angle_choc_Mc1*coef,Angle_paroi_Mc1*coef, Mach_Aval_Mc1,Kp_Mc_unitaire
            flag_Mc_unite=0
        end if

        if (flag_theta_max.ne.0) then
            write(2,100)Angle_Paroi*coef,Angle_choc*coef,Kp,Inv_Mach_Cone,Mach0,Mach_aval
        end if
        Inv_Mach_Cone_old=Inv_Mach_Cone
        Angle_choc_old=Angle_Choc
        Kp_old=Kp
        Angle_paroi_old=Angle_paroi
        Mach_cone_old=Mach_Cone
        Mach_aval_old=Mach_aval

    end do
    write(80,100)0.0d0,90.d0,Mach0
    write(80,*)'#'; write(85,*)'#'
    ! Kp=-2./gam/Mach0/Mach0                          ! coefficient de pression fonction de P_aval/P_amont
    ! je ne trace pas le dernier point en haut, la solution n'est pas physique:
    !write(2,100)0.d0,90.d0,Kp, 1.d0-1.d0/mach1(Mach0),Mach0   ! pas les bonnes valeurs. ?
    close(2); close(3)
end do
close(90);close(80);close(81);close(82)
print*,'fin normale du programme '
deallocate(table_Mach)
100 format(10(f15.8,3x))
120 format('#',2x,'Theta paroi en °',2x,'Angle Choc en °',8x,'Kp',11x,'Inverse Mach',7x,'Mach amont',&
       7x,'Mach aval')

contains

include 'sub_cone.f90'

!==================
! RK 4 Intégration
!==================

!****************************
subroutine rhs(theta,y,dydt)
!****************************
implicit none
! t  : theta,  y =(U_r, U_theta=dU_r/dtheta)
! y(1)= U_r,  y(2)= d U_r / dtheta = U_theta
! system to solve
! d y1/dt = y2
! d y2/dt = (y2**2*y1-(gam-1)/2(1-y1**2-y2**2)*(2*y1+y2/tant(theta))/ &
!           ( (gam-1)/2*(1-y1**2-y2**2)-y2**2
real(kind=8),intent(in)        :: theta
real(kind=8),dimension(2)      :: dydt,y
real(kind=8)                   :: num,den,y1,y2
y1=y(1);y2=y(2)
num=-(gam-1.d0)/2.d0*(1.d0-y1**2-y2**2)*(2.d0*y1+y2/(tan(theta)))+y1*y2**2
den=(gam-1.d0)/2.d0*(1.d0-y1**2-y2**2)-y2**2
dydt(1)=y2
dydt(2)=num/den
end subroutine rhs

!**************************
subroutine rk4(dt,t,y,dydt)
!**************************
implicit none
real(kind=8),dimension(2)      :: y,dydt
real(kind=8),intent(in)        :: dt,t
real(kind=8),dimension(2)      :: k1,k2,k3,k4,y_tmp
! t : theta, y(1)=V_r, y(2)=U_theta=d V_r/d theta
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

!***********************************************
subroutine U_theta(sauver,theta_init,thetac,y0,y,dy_dt)
!***********************************************
implicit none
logical,intent(in)                     :: sauver
real(kind=8),dimension(2),intent(in)   :: y0
real(kind=8),intent(in)                :: theta_init,thetac
real(kind=8),dimension(2),intent(inout):: y,dy_dt  !y(:,1)=U_r, y(:,2)=U_theta
integer,parameter                      :: nc_max=1000
real(kind=8),dimension(nc_max)         :: t
integer,parameter                      :: iter_max=30
real(kind=8)                           :: dt
integer                                :: i
if (sauver) write(3,110)k, Angle_choc*coef
dy_dt=0.d0
! on peut modifier le pas d'intégration
!nc=int((abs(thetac-theta_init)/dteta))+1
!if (nc.gt.nc_max) nc=nc_max
dt=(thetac-theta_init)/real(nc-1,kind=8)
print*,'U_theta :     dt =',dt*coef,'  nc = ',nc
do i=1,nc-1
    t(i)=theta_init+dt*real(i-1,kind=8)
end do
t(nc)=thetac
y=y0
do i=1,nc
    call rk4(dt,t(i),y,dy_dt)
    if (sauver.and.(opt_velocity.eq.1)) write(3,100) t(i)*coef,y(1),y(2),dy_dt(2),Angle_Choc*coef,Mach0
end do
100 format(6(e20.13,3x))
110 format('# theta',i4,/,'#',7x,f12.4,/,'#')
end subroutine U_theta

!************************************************************************
subroutine solve_angle_cone(theta_init,theta_cone,Ur,dUr_dtheta,dU_theta)
!************************************************************************

implicit none
real(kind=8),intent(inout)  :: theta_init
real(kind=8),intent(inout)  :: Ur,dUr_dtheta,dU_theta
real(kind=8),intent(inout)  :: theta_cone       !angle initial et final
integer                     :: iter
real(kind=8)                :: error,delta_theta
integer,parameter           :: iter_max=20
real(kind=8),dimension(2)   :: y0,y,dy_dt
print*,'solve_angle_cone, error_max =',error_max
y0(1)=Ur; y0(2)=dUr_dtheta
error=1.d0
iter=0
!theta_cone=0.25d0*theta_init
if (theta_cone.le.0.d0) then 
    print*, 'theta_cone ', theta_cone*coef, 'theta_init = ', theta_init*coef
    stop 'erreur theta_cone'
end if
do while ((error.gt.error_max).and.(iter.le.iter_max)) 
    iter=iter+1; print*,'Newton : iteration = ',iter,'theta c = ',theta_cone*coef
    call U_theta(.false.,theta_init,theta_cone,y0,y,dy_dt)
    error=abs(y(2))
    delta_theta=-y(2)/dy_dt(2)
    write(50,*) iter,theta_cone*coef,error, delta_theta*coef
    theta_cone=theta_cone+delta_theta
    if (theta_cone.le.0.d0) then 
        print*, 'theta_cone ', theta_cone*coef, 'theta_init = ', theta_init*coef
        stop 'erreur theta_cone boucle'
    end if

end do

if (iter.eq.iter_max) then
    stop 'no convergence in Newton_iteration procedure'
else
    print*,'convergence dans NR'
end if
call U_theta(.true.,theta_init,theta_cone,y0,y,dy_dt)
Ur=y(1);dUr_dtheta=y(2);dU_theta=dy_dt(2)
write(50,*) 'solution : en degrés ', theta_cone*coef, ' avec ', iter, ' iterations'

end subroutine solve_angle_cone

!************************************************************************
subroutine solve_angle_cone2(theta_init,theta_cone,Ur,dUr_dtheta,dU_theta)
!************************************************************************

implicit none
real(kind=8),intent(inout)  :: theta_init
real(kind=8),intent(inout)  :: Ur,dUr_dtheta,dU_theta
real(kind=8),intent(inout)  :: theta_cone       !angle initial et final
integer                     :: iter
real(kind=8)                :: error,teta,teta_old,alpha
integer,parameter           :: iter_max=20
real(kind=8),dimension(2)   :: y0,y,dy_dt,y_old,dy_dt_old

y0(1)=Ur; y0(2)=dUr_dtheta
error=1.d0
iter=0
write(3,110)k, Angle_choc*coef
teta=theta_init;  
teta_old=teta; y_old=y0; dy_dt_old=0.d0
y=y0
do while (teta.ge.0.d0)
    call rk4(-dteta,teta,y,dy_dt)
    if (opt_velocity.eq.1)  write(3,100) teta*coef,y(1),y(2),dy_dt(2),Angle_Choc*coef,Mach0
    if (y(2).gt.0.d0) then                           ! V_theta est <0 et = 0 à la paroi
         alpha=-y_old(2)/(y(2)-y_old(2))
         y=alpha*y+(1.d0-alpha)*y_old
         dy_dt=alpha*dy_dt+(1.d0-alpha)*dy_dt_old
         theta_cone=alpha*teta+(1.d0-alpha)*teta_old
         write(3,100) theta_cone*coef,y(1),y(2),dy_dt(2),Angle_Choc*coef,Mach0
         !Angle_Paroi=teta                           ! le cone est atteint
         exit
    endif
    y_old=y; dy_dt_old=dy_dt
    teta_old=teta
    teta=teta-dteta
end do
Ur=y(1);dUr_dtheta=y(2);dU_theta=dy_dt(2)

100 format(6(e20.13,3x))
110 format('# theta',i4,/,'#',7x,f12.4,/,'#')
end subroutine solve_angle_cone2




end program ecoulement_conique
