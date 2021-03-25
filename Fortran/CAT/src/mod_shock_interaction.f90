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


!>@details
!! MODULE : to treat  interaction between to shock waves
!! 
!! @author C. Airiau
!!

module mod_shock_interaction


contains

!********************************************************
subroutine main_shock_interaction(cas)
!********************************************************
!>@details
!! C. Airiau, Oct 2015
!! in (theta, p) plane : shock curve iso Mach,
!!                           and expansion curve at iso Mach
!!
!! cas 1         :: shock -  Mach (choice)
!!
!! cas 2         :: choc - given Mach 
!!
!! cas 3         :: isentropic expansion
!!
!! cas 4         :: 2 oblique shocks + isobar line
!!
!! cas 5         :: intersection of 2 curves
!!
!! cas 6         :: intersection Shock / isobar line
!!
!! cas 7         :: characteristic lines in the (u,v) plane
!********************************************************
implicit none
integer,parameter               :: n = 10001
real(kind=8),parameter          :: gam=1.4d0,val_pi=4.d0*atan(1.d0)
real(kind=8),parameter          :: rad2deg=180.d0/val_pi
real(kind=8),parameter          :: deg2rad=1.d0/rad2deg
real(kind=8),dimension(n)       :: p_out,tetap_out
real(kind=8),dimension(:),allocatable  :: p2,t2,mach2
real(kind=8)                    :: p0,M0,teta_ref,P_downstr,M_downstr,p_val,t_val
real(kind=8)                    :: M1,M2,teta2,P_tmp,M_tmp
integer                         :: cas,cas_out
integer,parameter               :: casMax=11
character(len=1)                :: car
character(len=40),dimension(casMax) :: liste_menu

liste_menu(1)="Curves of simple interaction"
liste_menu(2)="Three successive shock"
liste_menu(3)="Expansion"
liste_menu(4)="Two shocks + isobar line"
liste_menu(5)="Intersection of 2 curves"
liste_menu(6)="Interaction Shock/Isobar line"
liste_menu(7)="Charactacteristic in (u,v) map"
liste_menu(8)="isobar line + shock - expansion"
liste_menu(9)="isobar line + shock - shock     "
liste_menu(10)="nozzle + recompression : case A "
liste_menu(11)="nozzle + recompression : case B "
!
select case (cas)
case (0)
    call menu(cas_out)
case (1,2,3,4,5)
     cas_out=cas
case default
    stop 'bad argument in main_shock_interaction'
end select

write(6,200)liste_menu(cas_out)

select case (cas_out)
case(1)
   write(6,*) 'enter upstream pressure'
    read(*,*) p0
   write(6,*) 'enter upstream pressure'
    read(5,*) M0
    write(6,*)'angle in degrees ? =>'
    read*,teta_ref
    call shock_interaction(M0,p0,teta_ref,P_downstr,M_downstr,p_out,tetap_out)
case(2)
    p0=1.d0;M0=4.d0;teta_ref=10.d0
    call shock_interaction(M0,p0,teta_ref,P_downstr,M_downstr,p_out,tetap_out)
    M0=M_downstr;p0=P_downstr
    call shock_interaction(M0,p0,teta_ref,P_downstr,M_downstr,p_out,tetap_out)
    M0=M_downstr;p0=P_downstr
    call shock_interaction(M0,p0,teta_ref,P_downstr,M_downstr,p_out,tetap_out)
case(3)
    write(6,*)'Expansion'
    write(6,*) 'enter upstream pressure'
    read(*,*) p0
    write(6,*) 'enter upstream pressure'
    read(5,*) M0
    call expansion(M0,P0,p_out,tetap_out)
case(4)
    !p0=1.d0;M0=1.4d0;teta_ref=4.d0
    ! zone 1
    p0=1.d0;M0=2.d0;teta_ref=5.d0
    call shock_interaction(M0,p0,teta_ref,P_downstr,M_downstr,p_out,tetap_out)
    ! zone 2 et 3
    M0=M_downstr;p0=P_downstr
    call shock_interaction(M0,p0,2.d0*teta_ref,P_downstr,M_downstr,p_out,tetap_out)
    ! zone 4
    M0=2.d0;p0=1.d0;
    call shock_interaction(M0,p0,3.d0*teta_ref,P_downstr,M_downstr,p_out,tetap_out)

case(5)
    allocate(t2(n),p2(n))
    M1=2.d0;p0=1;teta_ref=10.d0
    call shock_interaction(M1,p0,teta_ref,P_downstr,M_downstr,p_out,tetap_out)
    M2=4.d0;p0=1;teta_ref=10.d0
    call shock_interaction(M2,p0,teta_ref,P_downstr,M_downstr,p2,t2)
    call curve_interaction(M1,M2,tetap_out,p_out,t2,p2,p_val,t_val)
    deallocate(t2,p2)

case(6)
    allocate(t2(n),p2(n),mach2(n))
    p0=1.d0;M0=2.d0;teta_ref=10.d0
    call shock_interaction(M0,p0,teta_ref,P_downstr,M_downstr,p_out,tetap_out)
    call characteristic(M_downstr,P_downstr,teta_ref,'+',M0,0.25d0,mach2,p2,t2)

    p0=1.d0;M0=3.d0;teta_ref=20.d0
    call shock_interaction(M0,p0,teta_ref,P_downstr,M_downstr,p_out,tetap_out)
    write(6,*) 'upstream Mach     expansion ',M_downstr
    write(6,*) 'upstream pressure expansion ',P_downstr
    write(6,*) 'upstream teta     expansion ',teta_ref
    call characteristic(M_downstr,P_downstr,teta_ref,'+',M0,0.25d0,mach2,p2,t2)
    deallocate(t2,p2,mach2)
case(7)
    allocate(t2(n),p2(n),mach2(n))
    p0=1.d0;M0=2.5d0;teta_ref=0.d0
    write(6,*) 'enter M0, P0, teta0, + ou - => ' 
    read(*,*) M0,p0,teta_ref,car
    call characteristic(M0,P0,teta_ref,car,M0,0.5d0,mach2,p2,t2)
    deallocate(t2,p2,mach2)
case(8)
    allocate(t2(n),p2(n),mach2(n))
    p0=1.d0;M0=4.d0;teta_ref=10.d0;M1=2.d0
    call shock_interaction(M0,p0,teta_ref,P_downstr,M_downstr,p_out,tetap_out)
    call characteristic(M_downstr,P_downstr,teta_ref,'+',M_downstr+0.5,0.25d0,mach2,p2,t2)
    teta_ref=12.810385
    call shock_interaction(M1,p0,teta_ref,P_downstr,M_downstr,p_out,tetap_out)
    deallocate(t2,p2,mach2)
case(9)
    allocate(t2(n),p2(n),mach2(n))
    p0=1.d0;M0=2.d0;teta_ref=10.d0;M1=4.d0
    call shock_interaction(M0,p0,teta_ref,P_downstr,M_downstr,p_out,tetap_out)
    teta2=7.1594113
    M_tmp=M_downstr;P_tmp=P_downstr
    call shock_interaction(M_tmp,P_tmp,abs(teta_ref-teta2),P_downstr,M_downstr,p_out,tetap_out)
    call shock_interaction(M1,p0,teta2,P_downstr,M_downstr,p_out,tetap_out)
    deallocate(t2,p2,mach2)
case(10)
    allocate(t2(n),p2(n),mach2(n))
    p0=1.d0;M0=2.d0;teta_ref=8.3101883465242228d0;
    call shock_interaction(M0,p0,teta_ref,P_downstr,M_downstr,p_out,tetap_out)
    teta2=0.d0
    M_tmp=M_downstr;P_tmp=P_downstr
    call shock_interaction(M_tmp,P_tmp,abs(teta_ref-teta2),P_downstr,M_downstr,p_out,tetap_out)
    call characteristic(M_downstr,P_downstr,teta2,'+',M0,0.25d0,mach2,p2,t2)
    deallocate(t2,p2,mach2)
case(11)
    allocate(t2(n),p2(n),mach2(n))
    p0=1.d0;M0=2.d0;teta_ref= 12.6272;
    print*,"enter the deciation angle in degrees (12.6272)  :"
    read*,teta_ref
    call shock_interaction(M0,p0,teta_ref,P_downstr,M_downstr,p_out,tetap_out)
    teta2=0.d0
    M_tmp=M_downstr;P_tmp=P_downstr
    call shock_interaction(M_tmp,P_tmp,abs(teta_ref-teta2),P_downstr,M_downstr,p_out,tetap_out)
    deallocate(t2,p2,mach2)

end select

200 format('#',50('*'),/,'#',10x,a30,/,'#',50('*'),/)

contains

    !************************
    subroutine menu(cas)
    !************************
    !>@details  display the menu
    implicit none
    integer,intent(out)     :: cas
    integer                 :: k
    write(6,200)
    do while ((cas.lt.1).or.(cas.gt.casMax))
        do k=1,casMax
            write(6,110) liste_menu(k),k
        end do
        write(6,100)casMax
        read(*,*) cas
    end do
    200 format(50('*'),/,5x,"Diagramm of  Pressure/deviation",/,&
            8x,"for the study of shock interactions",/,10x,'C. Airiau, Oct. 2016',/,&
            50('*'),/)
    110 format(a40,' : [',i2,'] ')
    100 format('Entrer an integer between  1 and ',i2,'   : ')
    end subroutine menu

    !*********************************************
    subroutine expansion(M0,P0,M_out,P_out)
    !*********************************************
    !>@details isentropic expansion
    ! C Airiau, oct. 2015, from Giovannini, 2015
    implicit none
    real(kind=8),intent(in)         :: P0,M0
    real(kind=8),dimension(n),intent(out) :: M_out,P_out
    real(kind=8)                    :: dMach,M_init,M_end,delta=0.25d0
    real(kind=8)                    :: PM0,mu0,mu,Mach,d,d0,P
    real(kind=8)                    :: teta0,teta
    !character(len=30)               :: set_filename
    integer                         :: i

    open(1,form='formatted',file=trim(set_filename('mach_teta',M0,3)) )
    write(1,200)
    teta0=0.d0
    M_init=M0-delta*M0; M_end=M0+delta*M0
    dMach=(M_end-M_init)/real(n-1,kind=8)
    mu0=asin(1.d0/M0)*rad2deg
    PM0=omega(M0)
    d0=1.d0+(gam-1.d0)/2.d0*M0**2
    do i=1,n
        Mach=M_init+dMach*real(i-1,kind=8)
        mu=asin(1.d0/Mach)*rad2deg
        teta=teta0+PM0-omega(Mach)
        d=1.d0+(gam-1.d0)/2.d0*Mach**2
        P=P0*(d/d0)**(-gam/(gam-1.d0))
        write(1,100) Mach, teta, P, mu
        M_out(i)=Mach;P_out(i)=P
    end do
    close(1)
    write(6,'(a,5x,a)')'See the file :',trim(set_filename('mach_teta',M0,3))
    100 format(4(e12.5,2x))
    200 format('#',3x,'Mach',10x,'theta',11x,'P',12x,'mu')
    end subroutine expansion

    !******************************************************************************
    subroutine characteristic(M0,P0,teta0,type_onde,M_final,delta,M_out,P_out,theta_out)
    !*****************************************************************************
    ! C Airiau, oct. 2015
    !>@details   to plot characteristic lines
    implicit none
    real(kind=8),intent(in)         :: P0,M0,teta0,M_final,delta
    character(len=1),intent(in)     :: type_onde
    real(kind=8),dimension(n),intent(out) :: M_out,P_out,theta_out
    real(kind=8)                    :: dMach,M_init,M_end
    real(kind=8)                    :: PM0,mu0,mu,Mach,d,d0,P
    real(kind=8)                    :: teta
    real(kind=8)                    :: epsil,Vtilde
    !character(len=30)               :: set_filename
    integer                         :: i

    if (type_onde=='+') then
        epsil=-1.d0
    else
        epsil=1.d0
    end if

    open(1,form='formatted',file=trim(set_filename('mach_teta_expansion_',M0,3) ))
    write(1,200)

    M_init=M0*(1.d0-delta);
    M_end=M_final*(1.d0+delta)
    dMach=(M_end-M_init)/real(n-1,kind=8)
    mu0=asin(1.d0/M0)*rad2deg
    PM0=omega(M0)
    d0=1.d0+(gam-1.d0)/2.d0*M0**2
    do i=1,n
        Mach=M_init+dMach*real(i-1,kind=8)
        mu=asin(1.d0/Mach)*rad2deg
        teta=teta0+epsil*(PM0-omega(Mach))
        d=1.d0+(gam-1.d0)/2.d0*Mach**2
        P=P0*(d/d0)**(-gam/(gam-1.d0))
        Vtilde=Mach*sqrt((gam+1.d0)/2.d0/d)
        theta_out(i)=teta
        write(1,100) Mach, teta, P, Vtilde,Vtilde*cos(teta*deg2rad),Vtilde*sin(teta*deg2rad),mu
        M_out(i)=Mach;P_out(i)=P
    end do
    close(1)
    write(6,'(a,5x,a)')'See the file :',trim(set_filename('mach_teta',M0,3))
    write(6,*) 'interval of  theta = ',theta_out(1),theta_out(n)
    100 format(7(e12.5,2x))
    200 format('#',3x,'Mach',10x,'theta',11x,'P',12x,'V/Vc',12x,'u/Vc',12x,'v/Vc',12x,'mu')
    end subroutine characteristic

    !*****************************
    function omega(M) result(res)
    !*****************************
    !>@details Prandtl-Meyer function in degrees
    implicit none
    real(kind=8)                    :: M,tmp,beta,res
    real(kind=8),parameter          :: c1=(gam+1.d0)/(gam-1.d0)

    if (M.lt.1.d0) then
        write(6,*)'error in omega for  Mach < 1   => M = ',M; stop
    end if
    beta=sqrt(M**2-1.d0)
    tmp=sqrt(c1)*atan(beta/sqrt(c1))-atan(beta)
    res=tmp*rad2deg
    end function omega

    !************************************************    
    function P_Pi(mach)
    !************************************************    
    !>@details    statique pressure ratio / isentropic pressure 
    implicit none
    real (kind=8)::  P_Pi,mach,y
    y=-gam/(gam-1.d0)
    P_Pi= (1.d0+ 0.5d0*(gam-1.d0)* mach**2)**y
    end function P_Pi


    !********************************************************
    subroutine shock_interaction(M0,p0,teta_ref,P_downstr,M_downstr,p,tetap)
    !********************************************************
    ! C Airiau, oct. 2015, from Giovannini, 2015
    !
    !>@details Pitot Tube relationship
    implicit none

    real(kind=8),intent(in)         :: p0,M0,teta_ref
    real(kind=8),intent(out)        :: P_downstr,M_downstr
    real(kind=8),dimension(n),intent(out)   :: p,tetap

    real(kind=8),parameter          :: eps=1.d-9
    real(kind=8),dimension(n)       :: M1,sigma
    real(kind=8)                    :: Mx1,Mx0,dMx,Mx_star,Mn0,Mn1,M0_star,den,anum
    real(kind=8)                    :: teta,My_star,teta_max,coef,tmp
    integer                         :: i
    logical,parameter               :: display=.true.
    !character(len=30)               :: set_filename

    write(6,300)M0,P0,teta_ref

    M0_star=M0*sqrt((gam+1.d0)/(2.0+(gam-1.d0)*M0**2))
    write(6,'(2(a15,f7.4,3x))')'M amont  =',M0, 'Mach M0^* = ',M0_star
    ! incident shock
    Mx0=1/M0_star+eps; Mx1=M0_star-eps; dMx=(Mx1-Mx0)/real(n-1,kind=8)

    do i=1,n
        Mx_star=Mx0+real(i-1,kind=8)*dMx
        anum=(M0_star-Mx_star)**2*(Mx_star*M0_star-1.0d0)
        den= 2.d0*M0_star**2/(gam+1.d0)-Mx_star*M0_star+1.0d0
        My_star=sqrt(anum/den)
        teta=atan(My_star/Mx_star)
        sigma(i)=atan((M0_star-Mx_star)/My_star)
        Mn0=M0*sin(sigma(i))
        p(i)=p0*(2.d0*gam*Mn0**2+1.d0-gam)/(gam+1.d0)
        tetap(i)=teta*rad2deg
        Mn1=sqrt((2.d0+(gam-1.d0)*Mn0**2)/(2.d0*gam*Mn0**2+1.d0-gam))
        M1(i)=Mn1/sin(sigma(i)-abs(teta))
    end do

    !
    !  display_outputs
    !
    open(2,form='formatted',file=trim(set_filename('theta_pression',M0,3)) )
    open(1,form='formatted',file=trim(set_filename('mach_pression',M0,3)) )
    write(1,*)'#  M0, M1, P, theta, sigma, P/Pi0'
    write(2,*)'#  theta+, theta-, P, P/Pi0'
    do i=n,1,-1
        tmp=p(i)*P_Pi(M0)
        write(1,*)M0,M1(i),p(i),tetap(i),sigma(i)*rad2deg,tmp
        write(2,*)tetap(i),-tetap(i),p(i),tmp
    end do
    close(1);close(2)
    teta_max=maxval(tetap)
    write(6,'(a40,f7.4)') 'Maximal deviation angle = ',teta_max

    if (teta_ref.ne.0.d0) then
        write(6,100)'searched theta   = ',teta_ref
        do i=n,1,-1
            if (tetap(i)-teta_ref.gt.0.d0) exit
        end do
        if (i.le.1) then
            write(6,*) 'bad value of teta '
        else
            write(6,'(a40,i5)' ) 'i = ',i
            coef=(teta_ref-tetap(i))/(tetap(i+1)-tetap(i))
            P_downstr=coef*p(i+1)+(1.d0-coef)*p(i)
            M_downstr=coef*M1(i+1)+(1.d0-coef)*M1(i)
            if (display) then
                write(6,100)'tetap(i)        = ',tetap(i)
                write(6,100)'tetap(i+1)      = ',tetap(i+1)
                write(6,100)'coef =          = ',coef
                write(6,100)'p(i)            = ',p(i)
                write(6,100)'p(i+1)          = ',p(i+1)
                write(6,100)'M1(i)           = ',M1(i)
                write(6,100)'M1(i+1)         = ',M1(i+1)
            end if
            write(6,100)'p downstream        = ',P_downstr
            write(6,100)'Mach downstream     = ',M_downstr
        endif
    else
        p_downstr=p0
        M_downstr=M0
    end if

    100 format(a40,f10.4)
    300 format('Upstream Mach = ',f7.4,3x,'Upstream pressure  = ',f7.4,3x,&
            'Deviation angle in degrees = ',f7.4)
    end subroutine shock_interaction

    !***********************************************************
    subroutine curve_interaction(M1,M2,t1,p1,t2,p2,p_out,t_out)
    !************************************************************
    !>@details  to verify curve_intersection  
    implicit none
    real(kind=8),dimension(n),intent(in):: t1,t2,p1,p2 ! pressure et theta angle
    real(kind=8),intent(in)             :: M1,M2
    real(kind=8),intent(out)            :: p_out,t_out
    logical                             :: error,trouve
    integer                             :: i,j,indMax,m,counter
    real(kind=8)                        :: thetaMax
    real(kind=8),dimension(n)           :: x1,x2,y1,y2 
    real(kind=8)                        :: x,y
    write(6,'(a,f10.4,a,f10.4)') 'Curve intersections for  Mach numbers of  ',M1,' and ', M2

    call Maximun(n,t1,thetaMax,indMax)
    write(6,'(a,f10.4,3x,a,f12.5,a,i5)') 'For M = ',M1,'theta Max = ',thetaMax,' at the index i= ',indMax
    call Maximun(n,t2,thetaMax,indMax)
    write(6,'(a,f10.4,3x,a,f12.5,a,i5)') 'For M = ',M2,'theta Max = ',thetaMax,' at the index i= ',indMax
    ! look for intersection
    error=.false.;  trouve=.false.
    p_out=0.d0;t_out=0.d0
    counter=0
    call flip(n,t1,y1);call flip(n,t2,y2);call flip(n,p1,x1);call flip(n,p2,x2)

    ! write(6,*)x1(1),x1(n),x2(1),x2(n)
    if ( (x1(n)-x2(1))*(x1(n)-x2(n)).lt.0.d0) then
        write(6,*)'A solution exists'
    else
        write(6,*)"No intersection exists"
        error=.true.; return
    end if

    m=1
    ! I am looking for the maximal  of the x index (x2(m)>x1(n)) , to avoid unuseful loop
    do while (x2(m).lt.x1(n))
        m=m+1
    end do
    !write(6,*)'last interesting x2 index : m = ',m

    do j=1,m-1
        do i=1,n-1
            call line_intersection(x1(i:i+1),y1(i:i+1),x2(j:j+1),y2(j:j+1),x,y,trouve)
            if (trouve) then
                write(6,*)"a solution if found"
                write(6,210)M1,M2, i,i+1,j,j+1
                write(6,100) 'theta: ',y1(i:i+1),y2(j:j+1),'pressure: ',x2(i:i+1),x2(j:j+1)
                write(6,200) y,x
                100 format(a10,2f12.5,2x,'|',2f12.5,/,a10,2f12.5,2x,'|',2f12.5)
                200 format('theta approx/ = ',f12.5,2x,'P approx. = ',f12.5)
                210 format('Table of values : ',/,10x,'M1 = ',f10.3,11x,'| M2 = ', f10.3,&
                    /,'Indices : ',3x,i5,7x,i5,6x,'|',5x,i5,7x,i5)
                counter=counter+1
                trouve=.false.
            end if
        end do
    end do

    if (counter.eq.0) then
        write(6,*)"no intersection"
    end if

    end subroutine curve_interaction


end subroutine main_shock_interaction

!*******************
subroutine flip(n,u,v)
!*******************
!>@details inversionof the vector components
implicit none
integer,intent(in)                     :: n
real(kind=8),dimension(n),intent(in)   :: u
real(kind=8),dimension(n),intent(out)  :: v
integer                                :: i
do i=n,1,-1
    v(i)=u(n-i+1)
end do

end subroutine flip


!********************************
subroutine Maximun(n,u,Umax,indMax)
!********************************
!>@details search maximal value of the vector and its index
implicit none
integer,intent(in)                  :: n
real(kind=8),dimension(n),intent(in):: u
real(kind=8),intent(out)            :: Umax
integer,intent(out)                 :: indMax
integer,dimension(1)                :: ind_tmp

Umax=maxval(u); ind_tmp=maxloc(u); indMax=ind_tmp(1)
end subroutine Maximun


!*****************************************************
subroutine line_intersection(x1,y1,x2,y2,x,y,test)
!*****************************************************
!>@details intersection of two straitgh lines
implicit none
real(kind=8),dimension(2),intent(in)    :: x1,y1,x2,y2
real(kind=8),intent(inout)              :: x,y
logical,intent(inout)                   :: test
real(kind=8)                            :: a1,b1,a2,b2,tmpx,tmpy
x=0.d0;y=0.d0
! equation of the first line
a1=(y1(2)-y1(1))/(x1(2)-x1(1))
b1=y1(1)-a1*x1(1)

! equation of the second line
a2=(y2(2)-y2(1))/(x2(2)-x2(1))
b2=y2(1)-a2*x2(1)

! intersection point
x=-(b2-b1)/(a2-a1)
y=a1*x+b1

test=.false.
tmpx=(x-x1(2))*(x-x1(1))
tmpy=(x-x2(2))*(x-x2(1))
!tmpy=(y-y1(2))*(y-y1(1))
if ((tmpx.le.0.d0).and.(tmpy.le.0.d0)) test=.true.
! if the point is located in the intervall of the x1 or y1.
end subroutine line_intersection



!****************************************************
function set_filename(root,value,digit) result(filename)
!****************************************************
!>@details
!! set a filename including the value of a parameter
!! root      : root of the file name
!! value     : value of the parameter as integer to avoid comma or dot
!! digit     : number of digit,    

implicit none

character(len=*),intent(in) :: root
integer,intent(in)          :: digit
real(kind=8)                :: value
real(kind=8),parameter      :: epsil=1.0d-7
character(len=30)            :: filename
real(kind=8)                :: m
integer,dimension(8) :: values


call date_and_time(VALUES=values)
! character(len=digit)        :: carac
m=10.d0**(digit-1)
filename=trim(root)//trim(carac(digit,int(m*(value+epsil))))//'_'//carac(4,values(8))//'.dat'
!write(6,*)' new filename ',filename
end function set_filename

!****************************
FUNCTION carac(N_in,i_enter)
!****************************
!>@details
!
!!  to convert an integer on a string
!  exemple : 123 --> '123'
!  N_in : number of digits
!  i_enter : integer

IMPLICIT NONE
INTEGER, INTENT(IN) :: i_enter,N_in
CHARACTER(LEN=N_in) :: carac
INTEGER             :: i,k,n

IF (N_in.le.1) THEN
    write(6,*) 'Problem in the charac function'
    STOP
ENDIF
n=N_in-1
i=i_enter
DO k=n,0,-1
    carac(n-k+1:n-k+1)=CHAR(48+INT(i/10**k))
    i=MOD(i,10**k)
END DO
END FUNCTION carac


end module mod_shock_interaction


! PROGRAM test
!     use mod_shock_interaction

!     call main_shock_interaction(0)
! end PROGRAM test
