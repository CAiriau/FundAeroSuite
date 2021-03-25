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

 
!> MODULE Constante : contient les constantes des problèmes
!======================
module Constante
!======================
! définition des constantes du programme
    real(kind=8),parameter        :: gam=1.4d0    !<   \f$ \displaystyle \gamma =   \frac{C_p}{C_v}  \f$ 
    real(kind=8),parameter        :: val_pi=4.d0*atan(1.d0)
    real(kind=8),parameter        :: rad2deg=180.d0/val_pi, deg2rad=val_pi/180.d0
    real(kind=8),parameter        :: r=287.058d0  !<   \f$ \displaystyle r =   C_p-C_v  \f$ 
    ! grandeur de référence pour le calcul de l'atmosphère standard 1976
    real(kind=8),parameter        :: P_atm=101325.d0    !< pression de reference
    real(kind=8),parameter        :: T_atm=288.15d0     !< temperature de reference
    real(kind=8),parameter        :: rho_atm=P_atm/(T_atm*r)  !< masse volumique de refence
    logical,parameter             :: fichier=.true. !< mettre true pour une sortie dans un fichier, false pour une sortie à l ecran.
end Module Constante


 
!> MODULE type_def : definitions of fortran types
!======================
Module type_def
!======================
!> type to define the physical aerodynamic quantities (few uses)
type state
real(kind=8)        :: M,T,P,rho,u,Ti,Pi,s
end type state

!> type to define input data, for the menu only
type input
integer             :: n
character(len=100)  :: title=''
end type

end Module type_def

!> MODULE mod_CAT  : it contains all the subroutine used to solve compressible aerodynamic problems
!!
!!         
!!
!! HISTORT :
!!
!!      - May 2001  : first subroutine for oblique shock waves and P.M. expansion   
!!      - March 2010 : translation from f77 to f90
!!      - May 2014 to August 2015 : improvements, configuration read, automatisation of the run
!!      - December 2016 to july 2017 : 
!!                  - include exercises of the book "Exercices et Problèmes d'aérodynamique fondamentale"
!!                  - add of shock interaction p
!!      - April 2018  : 
!!                      - the exercise to solve are now entered into an input file 
!!                      - add some modules , 
!!                      - comment implementation for doxygen
!!      - April 2020 : 
!!                      - translation into english, improvement of comments for doxygen
!!                      - EC has been renamed CAT for Compressible-Aero-Tool
!!
!! NOTATIONS
!!      - for shock waves :  upstream 1, downstream 2, index n : normal
!!      - isentropic or stagnation values : index i
!!
!! USE of CAT with command line into a shell terminal, several ways :
!!       1. CAT
!!       2. CAT opt= n                              # n : option
!!       3. CAT opt= n M= 3 teta= 10 A/Ac= 2.5      # n : option
!!
!>               depending on the choice (opt) M: Mach, teta : angle, ans : A/Acritique
!!
!>               M corresponds as wall to   'cas' for option 10 :
!>                   CAT opt= 10 M= 1  (for case 1)
!!
!>               space after = is very important
!!       4. use a script shell ./x1.run_CAT.sh 
!!         
!!           Caution:  the option is defined in the file  CAT.in
!!       5. A new mode has been implemented where the configuration is written in the file cas.in
!!          
!!
!! @author  Christophe Airiau
!! \note
!!  - Le program is run from the file  main_CAT.f90
!!  - The modules can be implemented in any  code following main_CAT.f90 example. 
!!
!
!==================================================================================
module mod_CAT
!==================================================================================
!
use Constante
use type_def
use mod_read
use atmos

implicit none

integer, parameter                :: ndim=600                   ! number of value for omega(M) table (P.M. function)
integer, parameter                :: menu_size=40               ! menu size
real(kind=8)                      :: Mn1,M1,M2,Mn2              ! Mach numbers
real(kind=8)                      :: om2,om1                    ! omega values (P.M. function)
real (kind=8), dimension(ndim)    :: xm,ym
real(kind=8)                      :: sigmac,teta                ! shock and deviation angle
real(kind=8)                      :: alpha=0.d0                 ! incidence angle
real(kind=8)                      :: ans                        ! answer value
integer                           :: i,j,ichoix
logical                           :: ask=.true.

integer                                     :: nz                     ! number maximal of zone to solve in a compressible flow
real(kind=8),dimension(:),allocatable       :: M,Ti,P,V               ! Mach, stagnation temperature, pressure, velocity
real(kind=8),dimension(:),allocatable       :: a,tau,T,Pi             ! sound velocity, pressure ratio, temperature, isentropic pressure
real(kind=8),dimension(:),allocatable       :: sigma, theta,Om,Mn     ! shock angle, flow deviation, Omega, normal Mach number
character(len=36),dimension(:),allocatable  :: com                    ! comments
character(len=12),dimension(:),allocatable  :: type_zone              ! type of zone
character(len=1),dimension(:),allocatable   :: type_characteristic    ! type of  characteristic line
integer,dimension(:),allocatable            :: upstr_zone             ! index of the upstream zone
integer,dimension(:),allocatable            :: profil                 ! -1 : upper_wall, 1: lower_wall, 0 : no airfoil
real(kind=8),dimension(:),allocatable       :: longueur               ! zone length in  % of chord
real(kind=8),dimension(:),allocatable       :: CL_zone,CD_zone        ! CL and CD in a given zone
real(kind=8)                                :: Cl,Cd                  ! lift and pressure drag coefficient
real(kind=8)                                :: tmp1,tmp2,tmp3,tmp4,F
character(len=80),dimension(0:menu_size)    :: menu_title=''
type(state)                                 :: q
type(input),dimension(0:menu_size)          :: menu 
real(kind=8)                                :: chord                  ! airfoil chord
real(kind=8)                                :: thickness              ! airfoil thickness parameter
real(kind=8)                                :: Kp                     ! pressure coefficient
integer                                     :: N_exo                  ! number of exercise to solve
character(len=4),allocatable,dimension(:)   :: liste_exo              ! list of exercises
logical                                     :: show=.true.

contains

!********************
subroutine solve_case
!********************
!> @details
!! SUBROUTINE to calculate the physical state as a function of the zone type
!! 
! 
implicit none
integer         :: k
write(6,*)'solve case'
do k=0,nz
    write(6,100) k
    select case (trim(type_zone(k)))
    case ('Uniform')
        call solve_uniform_zone(k)
    case ('Shock')
        call oblique_shock(upstr_zone(k),k)
    case ('Normal shock')
        call normal_shock(upstr_zone(k),k)
    case ('Expansion')
        call isentropic_expansion(upstr_zone(k),k)
    case ('Isobarline')
        call solve_isobarline(upstr_zone(k),k)
    case ('Nothing')
        call nothing(upstr_zone(k),k)
    case default
        write(6,*)'the zone is not found ...',k,' : ',trim(type_zone(k))
        stop
    end select
end do

100 format(50('*'),/,'ZONE : ', i3,/,50('*'))

end subroutine solve_case

!*********************
subroutine read_case
!*********************
!> @details
!! SUBROUTINE to read the flow configuration in the input file **cas.in**
!! 
! 
implicit none
integer                 :: i,k

open(1,form='formatted',file='cas.in')
read(1,*);read(1,*);read(1,*)
read(1,*) nz
read(1,*) alpha
allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),upstr_zone(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
do i=0,nz
    read(1,*)
    read(1,*) k
    read(1,*) type_zone(k)
    write(6,*) "zone n° ",k,trim(type_zone(k))
    select case (trim(type_zone(k)))
    case ('Uniform')
        read(1,*) M(k)
        read(1,*) Ti(k)
        read(1,*) P(k)
        read(1,*) theta(k)
        read(1,*) upstr_zone(k)
        if (upstr_zone(k).ne.k) theta(k)=theta(k)-alpha
    case ('Shock','Expansion','Nothing')
        read(1,*) theta(k)
        read(1,*) upstr_zone(k)
        if (upstr_zone(k).ne.k) theta(k)=theta(k)-alpha
    case ('Normal shock')
        read(1,*) upstr_zone(k)
    case ('Isobarline')
        read(1,*) theta(k)
        read(1,*) P(k)
        read(1,*) upstr_zone(k)
    end select
end do
theta=theta*deg2rad
close(1)

end subroutine read_case

!************************
subroutine set_menu
!************************
!> @details
!! SUBROUTINE display the menu options and the index of the option
!! 
! 
implicit none
integer     :: k
write(6,300)
300  format(60('*'),/,10x,'Paul Sabatier University, Toulouse, France',/,   &
    10X,'Compressible-Aero-Tool (CAT)',/,/, &
    20X,' the program is a part of  FundAeroSuite ',/,/, &
    2X, ' codes for education purpose only '  &
    ,/, 15X,'author : Christophe Airiau, 2015 - 2020',/,60('*'),/)

! k defines the order in the menu display
!
k=0 ; menu(k)%title='NORMAL SHOCK WAVE                    ';menu(k)%n=00
k=1 ; menu(k)%title='OBLIQUE SHOCK WAVE                   ';menu(k)%n=01
k=2 ; menu(k)%title='ISENTROPIC EXPANSION AND COMPRESSION ';menu(k)%n=02
k=3 ; menu(k)%title='Prandtl-Meyer omega(Mach) function   ';menu(k)%n=03
k=4 ; menu(k)%title='MACH NUMBER from a given omega angle ';menu(k)%n=04
k=5 ; menu(k)%title='ISENTROPIC EVOLUTION                 ';menu(k)%n=05
k=6 ; menu(k)%title='Shock interactions, (p,theta) plane  ';menu(k)%n=10
k=13; menu(k)%title='A / A critical for a given Mach      ';menu(k)%n=13
k=14; menu(k)%title='Mach for a given A / A critical      ';menu(k)%n=14 
k=16; menu(k)%title='Fanno problem  F=f(Mach)             ';menu(k)%n=16
k=17; menu(k)%title='Inverted Fanno  problem M=Fanno(F)   ';menu(k)%n=17
k=19; menu(k)%title='Rayleigh problem                     ';menu(k)%n=19
!
k=20; menu(k)%title='PLOT of Rayleigh function            ';menu(k)%n=20
k=21; menu(k)%title='PLOT of Fanno function               ';menu(k)%n=18
!
k=34; menu(k)%title='Rayleigh function  f(Mach)           ';menu(k)%n=21
k=35; menu(k)%title='inverted of Rayleigh function        ';menu(k)%n=22
!
k=22; menu(k)%title='EPICYCLOIDE PLOT                     ';menu(k)%n=08
k=23; menu(k)%title='Omega function of Mach PLOT          ';menu(k)%n=06
k=24; menu(k)%title='A / A critical = f(Mach) PLOT        ';menu(k)%n=15
k=25; menu(k)%title=' '                                    ;menu(k)%n=25   ! not used
k=26; menu(k)%title=' '                                    ;menu(k)%n=26   ! not used
k=27; menu(k)%title='Livre Aerodynamique Fondamentale     ';menu(k)%n=27
k=28; menu(k)%title='ENSEEIHT TD 1'                        ;menu(k)%n=28    
k=29; menu(k)%title='ENSEEIHT TD 2'                        ;menu(k)%n=29   
k=30; menu(k)%title='TABLES of the book                   ';menu(k)%n=12
k=31; menu(k)%title='Shock TABLES                         ';menu(k)%n=09
k=32; menu(k)%title='Prandtl-Meyer omega function TABLE   ';menu(k)%n=07
k=33; menu(k)%title='Personal TABLE                       ';menu(k)%n=23
k=36; menu(k)%title='Standard ATMOSPHERE                  ';menu(k)%n=24
k=39; menu(k)%title='May  2014 exam   (example)           ';menu(k)%n=90
k=38; menu(k)%title='June 2014 exam (example)             ';menu(k)%n=91
k=menu_size; menu(k)%title='configuration in input file'   ;menu(k)%n=99

end subroutine set_menu


!************************
subroutine display_menu
!************************
!> @details
!! SUBROUTINE to display the menu in the standard output (terminal or file)
! 
implicit none
integer     :: i
do i=0,menu_size
    if (menu(i)%title.ne.'') write(6,302) menu(i)%title,menu(i)%n
    if (mod(i+1,10).eq.0) write(6,100)
end do
302 format(3x,'=>',4x,a40,2x,'[ ',i2,' ]')
100 format(50('.'))
end subroutine display_menu

!***************
subroutine help
!***************
!> @details
!! SUBROUTINE to display an help on the standard output 
! 

write(6,100)

100 format('===>  Aide : ',/, 'Name of the program : CAT',/,     &
    'Make a choice, in command line, for instance:',/, &
    'CAT',/,                                                  &
    'CAT opt= n                    with n the option index ',/,        &
    'CAT opt= 1 M= 2.0 angle= 10.0 angle in deg',/,        &
    'CAT opt= 14 M= 0.5 A/Ac= 3.5  (for subsonic solution)',/,    &
    'CAT opt= 14 M= 1.5 A/Ac= 3.5  (for supersonic solution)',/,  &
    'CAT opt= 17 M= 0.5 F= 0.4     (for subsonic solution)',/,  &
    'CAT opt= 17 M= 1.5 F= 0.4     (for supersonic solution)',/,  &
    'CAT opt= 10 M= 1              (shock interaction, M is the case index )',/,  &
    'CAT h                         display help',/,           &
    /,' Do not forget the space after = sign in the command line !')
write(6,*)
write(6,'(a,/)') 'task choice  : '
call display_menu
end subroutine help

!*************************************
subroutine solve_uniform_zone(upstr)
!*************************************
    !> @details
    !! SUBROUTINE to calculate the state of an uniform zone
    !! 
    !! Inputs : Mach, P, Ti
    !!
    !! display_outputs : Pi, T, a, V, tau=P_downstream/P_upstream
    implicit none
    integer         :: upstr
    write(6,*) 'indice upstr         = ',upstr
    write(6,*) 'uniform zone : Mach  = ',M(upstr)
    com(upstr)='uniform upstream'
    if (Ti(upstr).ne.0.d0) then
        print*,'get temperature T'
        T(upstr)=T_Ti(M(upstr))*Ti(upstr)
    elseif (T(upstr).ne.0.d0) then
        print*,'get temperature Ti'
        Ti(upstr)=T(upstr)/T_Ti(M(upstr))
    else
        stop 'Problem with the uniform zone and T'
    end if
    write(6,100) 'T/Ti = ',T_Ti(M(upstr)),' T = ',T(upstr)
    write(6,100) 'Ti = ',Ti(upstr)
    a(upstr)=sqrt(gam*r*T(upstr));V(upstr)=M(upstr)*a(upstr);
    write(6,100) 'a = ',a(upstr),'V= ',V(upstr)
    if (P(upstr).ne.0.d0) then
        Pi(upstr)=P(upstr)/P_Pi(M(upstr))
    elseif (Pi(upstr).ne.0.d0) then
        P(upstr)=Pi(upstr)*P_Pi(M(upstr))
    else
        stop 'Problem with the uniform zone and P'
    end if

    write(6,100) 'P/Pi = ',P_Pi(M(upstr)),' Pi = ',Pi(upstr)
    write(6,100) 'P = ',P(upstr)
    write(6,100) 'ac = ',critical_velocity(Ti(upstr))
    write(6,100) 'V/ac = ',V(upstr)/critical_velocity(Ti(upstr))
    tau(upstr)=1.d0
    100 format(2(a20,3x,f15.4,5x))
end subroutine solve_uniform_zone

!**************************
subroutine nothing(upstr,downstr)
!**************************
    !> @details
    !! SUBROUTINE to copy a upstream uniform zone in another uniform zone
    implicit none
    integer         :: upstr,downstr
    write(6,*) 'uniform zone : Mach  = ',M(upstr) ;
    write(6,*) 'nothing  happens'
    write(6,*) 'state identical to zone ',upstr
    com(downstr)='state equal to the previous zone'
    M(downstr)=M(upstr)
    P(downstr)=P(upstr)
    tau(downstr)=tau(upstr)
    Pi(downstr)=Pi(upstr)
    theta(downstr)=theta(upstr)
    T(downstr)=T(upstr)
    Ti(downstr)=Ti(upstr)
    a(downstr)=a(upstr);V(downstr)=V(upstr);
    write(6,110)'M',downstr,M(downstr)
    write(6,110)'Ti',downstr,Ti(downstr)
    write(6,110) 'T/Ti',downstr,T_Ti(M(downstr))
    write(6,110) 'P/Pi',downstr,P_Pi(M(downstr))
    write(6,110) 'T',downstr,T(downstr)
    write(6,110) 'a',downstr,a(downstr)
    write(6,110) 'V',downstr,V(downstr)
    write(6,110) 'ac = ',critical_velocity(Ti(downstr))
    write(6,110) 'V/ac = ',V(downstr)/critical_velocity(Ti(downstr))

    110 format(a50,1x,'(',i3,')',':',3x,f12.4)
end subroutine nothing

!*****************************************
subroutine isentropic_expansion(upstr,downstr)
!*****************************************
!> @details
!! SUBROUTINE to solve a zone   
!! with an isentropic expansion or a compression
!! with the characteristic method
    implicit none
    integer         :: upstr,downstr
    real(kind=8)    :: deviation
    write(6,*) 'centered isentropic expansion : Mach = ',M(upstr)
    deviation=theta(downstr)-theta(upstr)
    write(6,100) 'deviaton angle in deg. ',deviation*rad2deg
    com(downstr)='isentropic expansion'
    Om(upstr)=omega(M(upstr))
    write(6,110)'Omega upstream in deg.',upstr,om(upstr)*rad2deg
    if (allocated(type_characteristic)) then
        if (type_characteristic(downstr)=='+') then
            write(6,*) 'outgoing wave C+'
            Om(downstr)=Om(upstr)+theta(downstr)-theta(upstr)
        else  
            write(6,*) 'incoming wave C-'
            Om(downstr)=Om(upstr)+theta(upstr)-theta(downstr)
        end if
    else
        write(6,*) 'the wave type is not given in type_characteristic'
        write(6,*) 'downstream Mach number is set from an expansion only'
        Om(downstr)=Om(upstr)+abs(deviation)
    end if

    M(downstr)=invomega(Om(downstr))
    write(6,110)'Omega upstream in deg.',upstr,om(upstr)*rad2deg
    write(6,110)'Omega downstream in deg.',downstr,om(downstr)*rad2deg
    write(6,110)'downstream Mach number',downstr,M(downstr)
    tau(downstr)=P_Pi(M(downstr))/P_Pi(M(upstr))
    write(6,110)'P/Pi',downstr,P_Pi(M(downstr))
    write(6,110) 'tau (downstream/upstream)',downstr,tau(downstr)
    P(downstr)=P(upstr)*tau(downstr)
    Pi(downstr)=Pi(upstr)
    Ti(downstr)=Ti(upstr)
    T(downstr)=T_Ti(M(downstr))*Ti(downstr)
    a(downstr)=sqrt(gam*r*T(downstr));V(downstr)=M(downstr)*a(downstr);

    write(6,110) 'T/Ti',downstr,T_Ti(M(downstr))
    write(6,110) 'T',downstr,T(downstr)
    write(6,110) 'a',downstr,a(downstr)
    write(6,110) 'V',downstr,V(downstr)
    write(6,110) 'Pressure P ',downstr,P(downstr)
    write(6,110) 'tau = P(k)/P(0)',downstr,P(downstr)/P(0)
    tau(downstr)=P(downstr)/P(0)

    write(6,120)'downstream Mach ',M(downstr),'downstream theta in deg.',theta(downstr)*rad2deg
    write(6,100) 'ac = ',critical_velocity(Ti(downstr))
    write(6,100) 'V/ac = ',V(downstr)/critical_velocity(Ti(downstr))

    100 format(a50,1x,':',3x,f12.4)
    110 format(a50,1x,'(',i3,')',':',3x,f12.4)
    120 format(2(a20,3x,f15.4,5x))

end subroutine isentropic_expansion
!***********************************
subroutine oblique_shock(upstr,downstr)
!***********************************
!> @details
!! SUBROUTINE to solve an oblique shock wave
    implicit none
    integer         :: upstr,downstr
    real(kind=8)    :: deviation

    write(6,*) 'oblique shock wave between zone ',upstr,' and zone ', downstr
    deviation=abs(theta(downstr)-theta(upstr))
    com(downstr)= 'oblique shock'
    write(6,100) 'upstream Mach',M(upstr)
    write(6,100) 'deviaton  angle in deg.',deviation*rad2deg
    call Newton_shock_angle(sigma(downstr),deviation,M(upstr))
    if (sigma(downstr).gt.-1.0d0) then
        write(6,100)'SHOCK ANGLE in deg. ',sigma(downstr)*rad2deg
        Mn(upstr)=M(upstr)*sin(sigma(downstr))
        write(6,110)'upstream NORMAL MACH Mn',upstr,Mn(upstr)
        Mn(downstr)=mach_downstr(Mn(upstr))
        M(downstr) = Mn(downstr)/sin(sigma(downstr)-deviation)
        write(6,110)'downstream NORMAL MACH Mn',downstr,Mn(downstr)
        write(6,110)'downstream MACH  M',downstr,M(downstr)
        write(6,110)'temperature ratio T/Ti'       ,downstr, T_Ti(M(downstr))
        write(6,120)'density ratio rho ',downstr,upstr,rho2_rho1(Mn(upstr))
        write(6,120)'pressure ratio P',downstr,upstr,P2_P1(Mn(upstr))
        write(6,120)'total pressure ratio Pi',downstr,upstr,Pi2_Pi1(Mn(upstr))

        tau(downstr)=P2_P1(Mn(upstr))*P(upstr)/P(0)
        Pi(downstr)=Pi(upstr)*Pi2_Pi1(Mn(upstr))
        P(downstr)=P(upstr)*P2_P1(Mn(upstr))
        write(6,110)'P/Pi',downstr,P_Pi(M(downstr))
        write(6,110)'tau= P(k)/P(0) =  ',downstr,tau(downstr)
        write(6,110)'Pi',downstr,Pi(downstr)
        write(6,110)'P',downstr,P(downstr)

        Ti(downstr)=Ti(upstr)
        write(6,110)'Ti',downstr,Ti(downstr)
        T(downstr)=T_Ti(M(downstr))*Ti(downstr)
        a(downstr)=sqrt(gam*r*T(downstr));V(downstr)=M(downstr)*a(downstr);
        write(6,110) 'T/Ti',downstr,T_Ti(M(downstr))
        write(6,110) 'T',downstr,T(downstr)
        write(6,110) 'a',downstr,a(downstr)
        write(6,110) 'V',downstr,V(downstr)
        write(6,100) 'ac = ',critical_velocity(Ti(downstr))
        write(6,100) 'V/ac = ',V(downstr)/critical_velocity(Ti(downstr))
    else
        write(6,*)' the shock is detached ! '
    end if

    100 format(a57,1x,':',3x,f12.4)
    110 format(a50,1x,'(',i3,')',':',3x,f12.4)
    120 format(a50,1x,'(',i3,' / ',i3,')',':',3x,f12.4)

end subroutine oblique_shock


!***********************************
subroutine normal_shock(upstr,downstr)
!***********************************
!> @details
!! SUBROUTINE to solve a normal shock wave
    implicit none
    integer         :: upstr,downstr
    sigma(downstr)=val_pi/2.d0
    write(6,*) 'normal shock bewteen zone ',upstr,' and zone ', downstr
    com(downstr)= 'normal shock'
    write(6,100) 'upstream Mach ',M(upstr)
    M(downstr)=mach_downstr(M(upstr))
    write(6,110)'downstream MACH M',downstr,M(downstr)
    write(6,110)'temperature ratio T/Ti',downstr, T_Ti(M(downstr))
    write(6,120)'density ratio rho ',downstr,upstr,rho2_rho1(M(upstr))
    write(6,120)'pressure ratio P',downstr,upstr,P2_P1(M(upstr))
    write(6,120)'total pressure ratio Pi',downstr,upstr,Pi2_Pi1(M(upstr))

    tau(downstr)=P2_P1(M(upstr))*P(upstr)/P(0)
    Pi(downstr)=Pi(upstr)*Pi2_Pi1(M(upstr))
    P(downstr)=P(upstr)*P2_P1(M(upstr))
    write(6,110)'P/Pi',downstr,P_Pi(M(downstr))
    write(6,110)'tau= P(k)/P(0) =  ',downstr,tau(downstr)
    write(6,110)'Pi',downstr,Pi(downstr)
    write(6,110)'P',downstr,P(downstr)

    Ti(downstr)=Ti(upstr)
    write(6,110)'Ti',downstr,Ti(downstr)
    T(downstr)=T_Ti(M(downstr))*Ti(downstr)
    a(downstr)=sqrt(gam*r*T(downstr));V(downstr)=M(downstr)*a(downstr);
    write(6,110) 'T/Ti',downstr,T_Ti(M(downstr))
    write(6,110) 'T',downstr,T(downstr)
    write(6,110) 'a',downstr,a(downstr)
    write(6,110) 'V',downstr,V(downstr)
    write(6,100) 'ac = ',critical_velocity(Ti(downstr))
    write(6,100) 'V/ac = ',V(downstr)/critical_velocity(Ti(downstr))
    
    100 format(a57,1x,':',3x,f12.4)
    110 format(a50,1x,'(',i3,')',':',3x,f12.4)
    120 format(a50,1x,'(',i3,' / ',i3,')',':',3x,f12.4)

end subroutine normal_shock


!*********************************************
subroutine solve_isobarline(upstr,downstr)
!*********************************************
!> @details
!! SUBROUTINE to solve a zone with a isobar line as a frontier
    implicit none
    integer         :: upstr,downstr
    real(kind=8)    :: theta_montante,theta_descendante,deviation
    write(6,*) 'isobar line solution'
    com(downstr)='isobar line'
    tau(downstr)=P(downstr)/P(upstr)
    M(downstr)=dsqrt(2.d0/(gam-1.d0)*( 1.d0/T_Ti(M(upstr))*tau(downstr)**(1.d0/gam-1.d0)-1.d0))
    om(upstr)=omega(M(upstr));om(downstr)=omega(M(downstr))
    deviation=abs(om(downstr)-om(upstr))
    Ti(downstr)=Ti(upstr)
    Pi(downstr)=P(downstr)/P_Pi(M(downstr))
    T(downstr)=T_Ti(M(downstr))*Ti(downstr)
    a(downstr)=sqrt(gam*r*T(downstr));V(downstr)=M(downstr)*a(downstr);

    write(6,110)'upstream MACH',upstr,M(upstr)
    write(6,110)'upstream omega in deg.',upstr,om(upstr)*rad2deg
    write(6,110)'downstream omega in deg. ',downstr,om(downstr)*rad2deg
    write(6,110)'downstream MACH',downstr,M(downstr)
    write(6,100)'deviaton angle in deg.',deviation*rad2deg
    theta_montante=theta(upstr)+om(downstr)-om(upstr)
    write(6,100)'theta_g, outgoing wave C+',theta_montante*rad2deg
    theta_descendante=theta(upstr)-om(downstr)+om(upstr)
    write(6,100)'theta_g, incoming wave C-',theta_descendante*rad2deg 
    write(6,110)'P(k)/P(upstream) =  ',downstr,tau(downstr)
    write(6,110)'Pi',downstr,Pi(downstr)
    write(6,110)'P',downstr,P(downstr)
    write(6,110)'P/Pi',downstr,P_Pi(M(downstr))
    write(6,110) 'T/Ti',downstr,T_Ti(M(downstr))
    write(6,110)'Ti',downstr,Ti(downstr)
    write(6,110) 'T',downstr,T(downstr)
    write(6,110) 'a',downstr,a(downstr)
    write(6,110) 'V',downstr,V(downstr)
    write(6,110)'tau= P(k)/P(0) =  ',downstr,tau(downstr)
    write(6,100) 'ac = ',critical_velocity(Ti(downstr))
    write(6,100) 'V/ac = ',V(downstr)/critical_velocity(Ti(downstr))

    if (allocated(type_characteristic)) then
        if (type_characteristic(downstr)=='+') then
            write(6,*) 'outgoing wave C+'
            theta(downstr)=theta_montante
        else  
            write(6,*) 'incoming wave C-'
            theta(downstr)=theta_descendante
        end if
    else
        write(6,*) 'the wave type is not given in type_characteristic'
        write(6,*) 'downstream theta is set arbitrarily to an ongoing wave C+'
        theta(downstr)=theta_montante
    end if
    write(6,120)'downstream Mach ',M(downstr),'downstream theta in deg',theta(downstr)*rad2deg
    tau(downstr)=P(downstr)/P(0)
    120 format(2(a20,3x,f15.4,5x))

    100 format(a57,1x,':',3x,f12.4)
    110 format(a50,1x,'(',i3,')',':',3x,f12.4)
end subroutine solve_isobarline


!***********************
subroutine calcul_Cl_Cd
!***********************
!> @details
!! SUBROUTINE to evaluate the Cl and Cd  coefficient as a function of the pressure ratio
    implicit none
    real(kind=8)    :: const
    write(6,*) 'Get aerodynamic coefficients for isocele triangle airfoil'
    const=-2.d0/(gam*M(0)**2)
    ! isocele triangle airfoil 
    Cl=const*(0.5d0/cos(theta(0))*(tau(1)+cos(theta(2))*tau(2))-cos(theta(3))*tau(3))
    Cd=const*(0.5d0/cos(theta(0))*sin(theta(2))*tau(2)-sin(theta(3))*tau(3))
    write(6,100) 'Cl ',Cl
    write(6,100) 'Cd ',Cd

    100 format(a50,1x,':',3x,f12.4)

end subroutine calcul_Cl_Cd

!*************************************
subroutine aerodynamic_coefficients
!*************************************
!> @details
!! SUBROUTINE to calculate aerodynamic coefficients
    implicit none
    real(kind=8)    :: coeff,tmp,beta
    integer         :: k
    coeff=2.d0/(gam*M(0)**2)
    write(6,*) ' aerodynamic coefficients'
    CL_zone=0.d0;CD_zone=0.d0
    write(6,100)
    do k=0,nz
        select case (profil(k))
        case(1)
            !write(6,*)'zone upper_wall :',k,theta(k)*rad2deg
            tmp=coeff*(tau(k)-1.d0)
            CL_zone(k)=tmp
            CD_zone(k)=tmp*tan(theta(k))
        case(-1)
            !write(6,*)'zone lower_wall :',k,theta(k)*rad2deg
            tmp=coeff*(tau(k)-1.d0)
            CL_zone(k)=tmp
            CD_zone(k)=tmp*tan(theta(k))
        case default
            longueur(k)=0.d0
        end select
        if (profil(k).ne.0.d0) then
            write(6,110)k,CL_zone(k),theta(k)*rad2deg,tau(k),CD_zone(k),CL_zone(k)*theta(k),&
                        longueur(k),profil(k)
            !write(6,*)'zone                                   = ', k
            !write(6,*)'Kp, theta en degrés, tau               = ',CL_zone(k),theta(k)*rad2deg,tau(k)
            !write(6,*)'Kp tan(theta) (nonlinear, linear)      = ',CD_zone(k),CL_zone(k)*theta(k)
        end if
    end do
    ! exact formula by the pressure force integration on the both directions
    write(6,120)'Total CL by Kp projection on panels           = ',-sum(CL_zone*longueur*profil)
    write(6,120)'Total CD by Kp projection on panels           = ', sum(CD_zone*longueur*profil)
    beta=sqrt(M(0)**2-1)
    ! Ackeret linear theory 
    ! only effect of incidence : 
    write(6,120)'CL linear theory : 4 alpha/beta               = ',4.d0*alpha*deg2rad/beta
    write(6,120)'CD linear theory : 4 alpha^2/beta             = ',4.d0*(alpha*deg2rad)**2/beta, ' for the flat plate ( CD0 = 0)' 
    ! thickness and incidence effect (incidence is included into theta)
    write(6,120)'CD theorie linéaire : int theta^2 dx (pente)  = ',2.d0/beta*sum(longueur*tan(theta)**2), 'for any body'
    write(6,120)'CD theorie linéaire : int theta^2 dx (angle)  = ',2.d0/beta*sum(longueur*theta**2), 'for any body'
    !            le tan theta in the formula above, and to respect linear theory assumption is given by the slope value, the theta angle itself
    !            it is different from the real tangente value!

    100 format('#', 1x,'zone', 4x,'Kp',9x,'theta (°)',8x,'P/P0',5x,'Kp tan(theta)',5x,'Kp theta' &
            3x,'length',3x,'upper_wall/lower_wall')
    110 format(i3,6(f12.5,2x),i2)
    120 format('# ',a50,2x,f12.5,2x,a)
end subroutine aerodynamic_coefficients


!********************
subroutine display_outputs
!********************
!> @details
!! SUBROUTINE to display outputs in the standard output - terminal or file
    implicit none
    write(6,90)
    do i=0,nz
        if (sigma(i).gt.-1.d0) then
            write(6,100) i, M(i),tau(i),P(i),Pi(i),theta(i)*rad2deg,sigma(i)*rad2deg,Om(i)*rad2deg,com(i),T(i),V(i),a(i)
        else
            write(6,110) i
        end if
    end do

    90 format(/,'#',80('*'),/,"# BILAN",/,'#',80('*'),/,'# i ',7x,'Mach',11x,'P/P0',10x,'P',13x,'Pi',14x,'theta',10x,&
        'sigma',10x,'Omega',14x,'Comments',25x,'T',13x,'V',13x,'a')
    100 format(i4,2x,7(f12.4,3x),a40,3(f12.4,3x))
    110 format(i3, 15x, ' ==> detached shock, no result')

end subroutine display_outputs


!
! I - Problem with shocks
!
!************************************************
subroutine Newton_shock_angle(sigma,teta,M0)
!************************************************
!> @details
!! Newton method to get the oblique shock angle 
    use Constante
    implicit none
    real(kind=8)              :: sigma,teta,M0
    real (kind=8)             :: er0,ds,f0,f1,dsigma,erreur
    integer i
    if (sigma.eq.0.d0) then
        write(6,*) 'classical initialization'
        sigma=asin(1.d0/M0)
    end if
    er0=1d-6
    ds=1d-2*deg2rad
    erreur=1.d0
    i=1
    if (show) write(6,101)
    101 format(50('-'),/, 'iter',4x,'Sigma',8x,'dsigma',10x,'error')
    100 format(i3,f12.5,2x,2(f14.9,2x))             
    do while ( (abs(erreur).gt.er0).and.(i.le.20))
        f0=shock_angle(sigma,teta,M0)
        f1=shock_angle(sigma+ds,teta,M0)
        dsigma=- ds * f0/(f1-f0)
        !      write(6,*)*,f0,f1,i
        erreur=dsigma/sigma
        if (show) write(6,100)i,sigma*rad2deg,dsigma*rad2deg,erreur
        sigma=sigma+dsigma
        i=i+1
    end do

    if (i.gt.20) then
        write(6,*)'no convergence'
        write(6,*)'  deviation angle is greater than'
        write(6,*)'  the angle to get an attached shock => detached shock ' 
        sigma=-1.d0
    else if (sigma.lt.0) then
        write(6,*) 'No convergence, angle < 0 :'
        write(6,*) 'Increase the upstream Mach '
        write(6,*) "or decrease the deviation angle"
        stop "error to calculate the oblique shock angle"
        sigma=-1.d0
    else if (sigma*rad2deg.gt.90.d0) then
        write(6,*) 'No convergence, angle > 90 :'
        write(6,*) 'Increase the upstream Mach'
        write(6,*) "or decrease the deviation angle"
        sigma=-1.d0
    else
        if (show) write(6,*) 'convergence on sigma'
    endif             
    if (sigma < 0) stop 'no  convergence on sigma'
    if (show) write(6,102)
    102 format(50('-'))  
end subroutine Newton_shock_angle


!************************************************       
function shock_angle(sigma,teta,M0) result(res)
!************************************************       
!> @details
!!     function to calcule oblique shock angle (density ratio) : f
    use Constante 
    implicit none 
    real (kind=8):: res, sigma,teta,M0
    res= tan (sigma-teta)/tan(sigma) -2.d0/(gam+1)/M0**2/(sin(sigma))**2 - (gam-1.d0)/(gam+1.d0)  
end function shock_angle


    !************************************************       
    function dfsigma(sigma,teta,M0)
!************************************************       
!> @details
!! derivative of the function shock_angle / sigma
    use Constante
    implicit none
    real (kind=8)::dfsigma, sigma,teta,M0,f1,f2
    f1= tan(sigma)/(cos(sigma-teta))**2 - tan (sigma-teta)/(cos(sigma))**2
    f2=f1/(tan(sigma))**2  +4.d0/(gam+1.d0)/M0**2    * cos(sigma)/(sin(sigma))**3
    dfsigma=f1+f2
end function dfsigma


!************************************************       
function Pi2_Pi1(mach)
!************************************************       
!> @details
!!     stagnation pressure jump across a shock
    use Constante
    implicit none
    real(kind=8)::  Pi2_Pi1,mach,x,y
    x=-1.d0/(gam-1.d0)
    y=-gam*x
    Pi2_Pi1= P2_P1(mach)**x * rho2_rho1(mach)**y
end function Pi2_Pi1

!************************************************       
function P2_P1(mach)
!************************************************       
!> @details
!!      static pressure jump across a shock     
    use Constante
    implicit none
    real (kind=8)::  P2_P1,mach
    P2_P1= 2.d0*gam/(gam+1.d0)* mach**2- (gam-1.d0)/(gam+1.d0)
end function P2_P1

!************************************************       
function Inverse_P2_P1(rapport) result(Mach)
!************************************************       
!> @details
!!     inverted function, as input the pressure jump across a shock,  
!!      as output : the normal upstream Mach number    
    use Constante
    implicit none
    real (kind=8)::  rapport,Mach
    Mach= sqrt(1.d0/(2.d0*gam)* ((gam+1.d0)* rapport+ (gam-1.d0)))
end function Inverse_P2_P1      


!************************************************       
function rho2_rho1(mach)
!************************************************       
!> @details
!!     density jump across a shock    
    use Constante
    implicit none
    real (kind=8):: rho2_rho1,mach
    rho2_rho1= 1.d0/( 2.d0/((gam+1.d0)* mach**2)+ (gam-1.d0)/(gam+1.d0) )
end function rho2_rho1


!************************************************       
function mach_downstr(mach)
!************************************************       
!> @details
!!     downstream normal Mach across a shock     
    use Constante
    implicit none
    real (kind=8):: mach,mach_downstr
    mach_downstr= sqrt((1.d0+ 0.5d0*(gam-1.d0)* mach**2)/(gam*mach**2-0.5d0*(gam-1.d0)))
end function mach_downstr


!************************************************       
function entropy_variation(Pi_upstr,Pi_Aval)
!************************************************
!> @details
!! entropy jump across a shock      
    use Constante
    implicit none
    real (kind=8):: Pi_upstr,Pi_Aval,entropy_variation
    entropy_variation=-r*log(Pi_Aval/Pi_upstr)
end function entropy_variation

    !
! II- Isentropic problems 
!

!************************************************    
function inverse_P_Pi(ratio) result(Mach)
!************************************************    
!> @details
!!     Mach as a function of the ratio static pressure / isentropic pressure
    use Constante
    implicit none
    real (kind=8):: Mach, ratio, y
    y=-(gam-1.d0)/gam
    Mach=sqrt((ratio**y-1.d0)*2.d0/(gam-1.d0))
end function inverse_P_Pi


!************************************************    
function P_Pi(mach)
!************************************************    
!> @details
!!     ratio static pressure / isentropic pressure 
    use Constante
    implicit none
    real (kind=8)::  P_Pi,mach,y
    y=-gam/(gam-1.d0)
    P_Pi= (1.d0+ 0.5d0*(gam-1.d0)* mach**2)**y
end function P_Pi


!************************************************       
function rho_rhoi(mach)
!************************************************       
!> @details
!!     ratio density / isentropic density
    use Constante
    implicit none
    real (kind=8)::  rho_rhoi,mach,y
    y=-1.d0/(gam-1.d0)
    rho_rhoi= (1.d0+ 0.5d0*(gam-1.d0)* mach**2)**y
    end function rho_rhoi


!************************************************       
function T_Ti(mach)
!************************************************       
!> @details
!!     ratio temperature / isentropic temperature 
    use Constante
    implicit none
    real(kind=8)::  T_Ti,mach
    T_Ti= 1.d0/ (1.d0+ 0.5d0*(gam-1.d0)* mach**2)
end function T_Ti


!************************************************       
function critical_velocity(Ti)
!************************************************       
!> @details
!!     critical velocity as a function of the isentropic temperature
    use Constante
    implicit none
    real(kind=8)::  critical_velocity,Ti
    critical_velocity=sqrt( 2*gam*r/(gam+1)*Ti)
end function critical_velocity


!************************************************
function omega(M)
!************************************************
!> @details
!!    Prandtl-Meyer (omega) function w.r.t. Mach number, in radians 
!!    two formulas are possible (cf comments )
    use Constante
    implicit none
    real (kind=8):: M,omega,c,beta
!   real (kind=8):: demi_pi,mu,omega0
    if (M.lt.1.d0) then
        write(6,*)'error in omega function :  Mach < 1   => M = ',M
        stop
    end if             
    ! demi_pi=val_pi/2.d0
    c=sqrt((gam+1.d0)/(gam-1.d0))
    beta=sqrt(M**2-1.d0)
    ! mu=asin(1.d0/M)
    ! omega0=mu+c*atan(1.d0/c/tan(mu))-demi_pi
    omega=c*atan(beta/c)-atan(beta)
    !      write(6,*)*,'OMEGA',omega0,omega,omega0/omega
    end function omega


!*************************************************  
function invomega(angle)
!*************************************************  
!> @details
!!    inverted omega function : Mach number is solved for a omega given 
!!  with a Newton method
    use Constante
    implicit none
    real (kind=8)::invomega,angle
    real (kind=8):: omega0,omega1,er0, erreur,mach0,dm,dmach
    integer i
    er0=1d-6
    if (show) write(6,*)'value of omega in deg.                     :  ', angle*rad2deg
    
    dmach=1d-3
    mach0=2.d0
    i=1
    erreur=1.d0

    do while ((erreur.gt.er0).and.(i.le.20))
        omega0=omega(mach0)-angle
        omega1=omega(mach0+dmach)-angle
        dm=- dmach * omega0/(omega1-omega0)
        erreur=abs(dm/mach0)
        if (show) write(6,*)'i = ',i,'error = ',erreur,'dm =',dm,'dM= ',mach0
        mach0=mach0+dm
        i=i+1
    end do

    if (i.gt.20) then
        write(6,*)'no convergence'
        stop
    else
        if (show) write(6,*)'convergence reached for Mach, iter              :  ', i
        invomega=mach0
    endif             
end function invomega


!**********************************************
subroutine omega_plot
!**********************************************
!> @details
!! to plot  omega (Mach)
    use Constante
    implicit none
    real(kind=8):: M1,M2,DM
    integer, parameter ::npt=1000
    integer i
    open(1,form='formatted',file='Prandtl_Meyer_function.dat')
    M1=1.d0;  M2=15.d0
    DM=(M2-M1)/float(npt-1)
    write(1,100) 
    100 format('#',10x,'function omega (Mach), in deg.')
    do i=1,npt
        write(1,'(2(e12.5,3X))')M1,omega(M1)*rad2deg
        M1=M1+DM
    end do   
    close(1)
end subroutine omega_plot     


!*******************
function S_sur_Sc(M)
!*******************
!> @details
!! ratio section / critical section  w.r.t. Mach number
    use Constante
    implicit none
    real(kind=8)            :: S_sur_Sc,M
    real(kind=8),parameter  :: c1=(gam-1.d0)/(gam+1.d0)
    real(kind=8),parameter  :: c2=0.5d0/c1
    S_sur_Sc=1.d0/M*(2.d0/(gam+1.d0)+c1*M**2)**c2
end function


!************************************************        
subroutine  epicycloide
!************************************************        
!> @details
!
!! epicycloide  plot in the isentropic expansion problem
!
    use Constante
    implicit none

    ! epicylcoide in critical coordinates
    integer, parameter          :: npt=10000
    real(kind=8), parameter     :: M_inf=1000.d0
    real(kind=8)                :: mach,dmach,mach_star,theta,theta_max,M_max
    integer                     :: i

    M_max=sqrt((gam+1.d0)/(gam-1.d0))
    open(1,form='formatted',file='epicycloide.dat')
    write(1,90)
    mach=1.d0
    dmach=(M_inf-mach)/float(npt-1)
    do i=1,npt
        call calcul_epi(mach,theta,mach_star)
        write(1,100) mach, theta*rad2deg,mach_star,mach_star*cos(theta),mach_star*sin(theta),&
        M_max*cos(theta), M_max*sin(theta),cos(theta), sin(theta)
        mach =mach+dmach
    end do
    100 format(9(e12.5,3x))
    90  format('#',4x,'Mach',9x,'theta', 9x,'Mach*',12x, 'u/a*',11x,'v/a*',&
    9x,'Mmax cos', 7x, 'Mmax sin', 9x, 'cos',11x, 'sin')
    close(1)

    theta_max= (M_max-1) *90.d0
    write(6,*) 'M* max                                    :  ', M_max
    write(6,*) 'theta_max                                 :  ', theta_max
end subroutine epicycloide


!**********************************
subroutine calcul_epi(x,omega,mach)
!**********************************
!> @details
!! get  the  Mach M* as a function of  M (to verify or M w.r.t. M* ?)
    use Constante
    implicit none
    real(kind=8):: x,omega, mach
    real(kind=8):: c1,c2
    c2=(gam-1.d0)/(gam+1.d0)
    c1=sqrt(1.d0/c2)
    omega = c1*atan(sqrt(c2*(x*x-1.d0)))-atan(sqrt(x*x-1.d0))
    mach=sqrt((gam+1.d0)*x*x/(2.d0+(gam-1.d0)*x*x))
end subroutine calcul_epi


!
! III-  tables
!

!**********************************
subroutine shock_table
!**********************************
!> @details
!! shock table printed in a output file, latex output
!
implicit none
integer                     ::n,i,j,iloc,step
real(kind=8),allocatable,dimension(:) ::m1,m2,rp,rpi,rti
real(kind=8)                ::dm,m1i,m1f
integer                     ::k=3

    m1i=1.00d0          !  initial Mach
    m1f=2.00d0          ! final Mach
    dm=0.010d0
    n=int((m1f-m1i)/dm)+2
    n=k*(int(n/k)+1)
    allocate(m1(n),m2(n),rp(n),rpi(n),rti(n))

    write(6,100) m1i,m1f,dm,n
    100 format('Shock table',/,'initial Mach  : ',f5.2,3x,'final Mach: ',f5.2,3x, &
        'Mach step : ',f5.3,3x,'Number of  points :',i5)
    m1(1)=m1i
    do i=1,n-1
        m2(i)=mach_downstr(m1(i)); rp(i)=P2_P1(m1(i));rpi(i)=P_Pi(m1(i))
        rti(i)=T_Ti(m1(i))
        m1(i+1)=m1(i)+dm
    end do

    open(1,form='formatted',file='shock_table.out')
    open(2,form='formatted',file='shock_table.tex')
    open(3,form='formatted',file='shock_table_vert.tex')
    write(1,110)
    write(2,111);write(3,111)
    111 format('\tiny{\begin{tabular}{',3(5('|c'),'|')'}',&
            '\hline')

    write(3,150)'&';write(3,150)'&';write(3,150)'\\ \hline'
    write(2,150)'&';write(2,150)'&';write(2,150)'\\ \hline'

    110 format(3(3x,'M1',4x,'M2',4x,'P2/P1 '))
    150 format('$Mn_1$&','$Mn_2$&','$P_2/P_1$&$P_1/P_{i1}$&$T_1/T_{i1}$',a)
    do i=1,n-k-1,k
        write(1,120)(m1(i+j),m2(i+j), rp(i+j),rpi(i+j),rti(i+j),j=0,k-1)
        write(2,122)(m1(i+j),m2(i+j), rp(i+j),rpi(i+j),rti(i+j),j=0,k-1)
        120 format(3('|',f5.3,'|',f6.4,'|',f7.4,'|',f6.4,'|',f6.4))
        122 format(2(f5.3,'&',f6.4,'&',f7.4,'&',f6.4,'&',f6.4,'&'),  & 
                    f5.3,'&',f6.4,'&',f7.4,'&',f6.4,'&',f6.4,'\\')
    end do
    step=int(n/k)-1
    write(6,*) 'step                                      :  ',step
    do i=1,step
        iloc=i
        write(3,122)(m1(iloc+j*step),m2(iloc+j*step), rp(iloc+j*step),rpi(iloc+j*step),rti(iloc+j*step),j=0,k-1)
    end do
    write(2,132);write(3,132)
    132 format('\hline \end{tabular}}')
    close(1);close(2);close(3)
end subroutine shock_table


!**********************
subroutine tables_livre
!**********************
!> @details
!! Table as there are displayed in the book
    use Constante
    implicit none
    real(kind=8)            :: dM,M_init,M_final,Mach
    integer                 :: nptM     ! number of points
    integer,dimension(4)    :: npt       
    integer                 :: nc       ! number of columns
    integer                 :: nl       ! number of rows
    integer                 :: reste,i,nl1,k,j,i1,i2,i_tmp
    real(kind=8),allocatable,dimension(:,:) :: table(:,:)
    real(kind=8),dimension(4)   :: dM_tmp,M_init_tmp,M_final_tmp
    logical                 :: test
    !=====================================
    write(6,*) 'table for isentropic subsonic flow'
    !=====================================

    nc = 2; k=6
    write(6,*) 'number of columns                         :  ', nc
    M_init=0.0d0; M_final=1.0d0; dM=0.01d0
    nptM=int((M_final-M_init)/dM)+1
    write(6,*) 'number of points                          : ', nptM
    allocate(table(nptM,k))
    do i=1,nptM
        Mach=M_init+dM*real(i-1,kind=8)
        table(i,1)=Mach
        table(i,2)=T_Ti(Mach);table(i,3)=P_Pi(Mach)
        table(i,4)=rho_rhoi(Mach);table(i,5)=S_sur_Sc(Mach)
        table(i,6)=Fanno(Mach);
    end do
    nl=nptM/nc
    reste=mod(nptM,nc)
    if (reste.ne.0) then
        write(6,*) 'one column will not be filled'
        write(6,*) 'It remains ',reste,' values'
        nl1=nl+1
    else
        nl1=nl
    end if

    !
    !   2 column : 1 file
    !
    write(6,*)' number of lines : ',nl,nptM,nptM/nc
    open(10,form='formatted',file='subsonic.tex')
    write(10,110)
    ! uniquement pour deux colonnes
    write(10,130)'&';write(10,130)'\\ \hline \hline'
    write(6,*) 'nl1                                       : ', nl1
    do i=1,nl1-1
        write(10,100) table(i,:),table(i+nl1,:)
        if (mod(i,5).eq.0) write(10,*) '\hline'
    end do
    write(6,*) 'mod                                       : ',mod(reste,2)
    if (mod(reste,2).eq.1) then
        write(6,*)' last line'
        write(10,101) table(nl1,:)
    end if
    write(10,120)
    close(10)

    !
    ! 1 column, 2 filess
    !
    open(10,form='formatted',file='subsonic1.tex')
    write(10,111); write(10,130)'\\ \hline \hline'
    do i=1,nptM/2
        write(10,102) table(i,:)
        if (mod(i,5).eq.0) write(10,*) '\hline'
    end do
    write(10,120)
    close(10)

    open(10,form='formatted',file='subsonic2.tex')
    write(10,111); write(10,130)'\\ \hline \hline'
    do i=nptM/2+1,nptM
        write(10,102) table(i,:)
        if (mod(i,5).eq.0) write(10,*) '\hline'
    end do
    write(10,120)
    close(10)

    deallocate(table)


!===================================================
write(6,*) 'table for isentropic supersonic flows'
!===================================================
    k=7
    M_init=1.00d0; M_final=5.0d0; dM=0.05d0
    nptM=int((M_final-M_init)/dM)+1
    write(6,*) 'number of points                          : ', nptM
    allocate(table(nptM,k))
    do i=1,nptM
        Mach=M_init+dM*real(i-1,kind=8)
        table(i,1)=Mach
        table(i,2)=T_Ti(Mach);table(i,3)=10.d0*P_Pi(Mach)
        table(i,4)=10.d0*rho_rhoi(Mach);table(i,5)=S_sur_Sc(Mach)
        table(i,6)=asin(1.d0/Mach)*rad2deg;table(i,7)=omega(Mach)*rad2deg
    end do
    nl=nptM/nc
    reste=mod(nptM,nc)
    if (reste.ne.0) then
        write(6,*) 'one column will not be filled'
        write(6,*) 'It remains ',reste,' values'
        nl1=nl+1
    else
        nl1=nl
    end if
    write(6,*)' number of lines                          : ',nl,nptM,nptM/nc
    !
    !   2 colonnes : 1 fichier
    !
    open(10,form='formatted',file='supersonic.tex')
    write(10,210)
    ! uniquement pour deux colonnes
    write(10,230)'&';write(10,230)'\\ \hline \hline'
    write(6,*) 'nl1                                       : ', nl1
    do i=1,nl1-1
        write(10,200) table(i,:),table(i+nl1,:)
        if (mod(i,5).eq.0) write(10,*) '\hline'
    end do
    write(6,*) 'mod                                       : ',mod(reste,2)
    if (mod(reste,2).eq.1) then
        write(6,*)' last line'
        write(10,201) table(nl1,:)
    end if
    write(10,220)
    close(10)
    deallocate(table)

    !
    ! 1 column, several file
    !

    dM_tmp=(/0.02d0, 0.05d0, 1.0d0, 0.0d0/)
    M_init_tmp =(/1.d0, 3.d0,  5.d0, 0.d0/)
    M_final_tmp=(/3.d0, 5.d0, 30.d0, 0.d0/)

    npt(4)=0
    do j=1,3
        write(6,*)'j, ',j,' step ', ((M_final_tmp(j)-M_init_tmp(j))/dM_tmp(j)), M_final_tmp(j),M_init_tmp(j),dM_tmp(j)
        npt(j)=int((M_final_tmp(j)-M_init_tmp(j))/dM_tmp(j))
    end do

    nptM=sum(npt)+1; 

    k=6
    allocate(table(nptM,k))
    write(6,*) 'npt ',npt(:), ' total = ', sum(npt)
    i1=1
    table(i1,1)=M_init_tmp(1)
    write(6,*)'initial Mach                               : ',table(1,1)
    do j=1,3
        i2=i1+npt(j)
        write(30,*)' Mach init and final and step',M_init_tmp(j), M_final_tmp(j),dM_tmp(j)
        do i=i1+1,i2
            table(i,1)=table(i-1,1)+dM_tmp(j)
            write(30,'(i4,3x,f12.3,5x,i3)')i,table(i,1),i-i1
        end do
        i1=i2
    end do

    do i=1,nptM
        Mach=table(i,1)
        table(i,2)=T_Ti(Mach);    table(i,3)=P_Pi(Mach)
        table(i,4)=rho_rhoi(Mach);  table(i,5)=S_sur_Sc(Mach)
        table(i,6)=Fanno(Mach); 
    end do

    do j=1,4
        open(10,form='formatted',file='supersonic'//charac(2,j)//'.tex')
        write(10,111); write(10,130)'\\ \hline \hline'
        do i=1,40
            i_tmp=i+(j-1)*40
            test=(table(i_tmp,5).lt.100.d0)
            if (minval(table(i_tmp,2:4)).gt.0.1d0) then
                if (test) then
                    write(10,103) table(i+(j-1)*40,1:k)
                else
                    write(10,104) table(i+(j-1)*40,1:k)
                end if
            else
                if (test) then
                    write(10,105) table(i+(j-1)*40,1:k)
                else
                    write(10,106) table(i+(j-1)*40,1:k)
                end if
            end if
            if (mod(i,5).eq.0) write(10,*) '\hline'
        end do
        write(10,120)
        close(10)
    end do
    deallocate(table)

    !==================================================
    write(6,*) 'table for normal and oblique shock '
    !==================================================
    k=6



    M_init=1.00d0; M_final=5.0d0; dM=0.05d0
    nptM=int((M_final-M_init)/dM)+1
    write(6,*) 'number of points : ', nptM
    allocate(table(nptM,k))
    do i=1,nptM
        Mach=M_init+dM*real(i-1,kind=8)
        table(i,1)=Mach                     ! upstream normal Mach
        table(i,2)=mach_downstr(Mach)       ! downstream normal Mach
        table(i,3)=P2_P1(Mach)              ! static pressure ratio
        table(i,4)=T_Ti(table(i,2))/T_Ti(Mach)  ! temperature ratio
        table(i,5)=rho2_rho1(Mach)           ! density ratio
        table(i,6)=Pi2_Pi1(Mach)            ! isentropic pressure ratio
    end do
    nl=nptM/nc
    reste=mod(nptM,nc)
    if (reste.ne.0) then
        write(6,*) 'one column will not be filled'
        write(6,*) 'It remains ',reste,' values'
        nl1=nl+1
    else
        nl1=nl
    end if
    write(6,*)' number of lines : ',nl,nptM,nptM/nc
    open(10,form='formatted',file='normal_shock.tex')
    write(10,310)
    ! for two columns
    write(10,330)'&';write(10,330)'\\ \hline \hline'
    write(6,*) 'nl1                                       : ',nl1
    do i=1,nl1-1
        write(10,300) table(i,:),table(i+nl1,:)
        if (mod(i,5).eq.0) write(10,*) '\hline'
    end do
    write(6,*) 'mod                                       : ',mod(reste,2)
    if (mod(reste,2).eq.1) then
        write(6,*)' last line'
        write(10,301) table(nl1,:)
    end if
    write(10,320)
    close(10)
    deallocate(table)

    !
    !  1 column, several files
    !
    k=6
    nptM=sum(npt)+1; 
    allocate(table(nptM,k))
    write(6,*) 'SHOCKS : npt ',npt(:), ' total = ', sum(npt)
    i1=1
    table(i1,1)=M_init_tmp(1)
    write(6,*)'initial mach  : ',table(1,1)
    do j=1,3
        i2=i1+npt(j)
        do i=i1+1,i2; table(i,1)=table(i-1,1)+dM_tmp(j); end do
        i1=i2
    end do

    do i=1,nptM
        Mach=table(i,1)
        table(i,2)=mach_downstr(Mach)               
        table(i,3)=P2_P1(Mach)              
        table(i,4)=T_Ti(table(i,2))/T_Ti(Mach)   
        table(i,5)=rho2_rho1(Mach)            
        table(i,6)=Pi2_Pi1(Mach)
    end do
    do j=1,4
        open(10,form='formatted',file='shocks'//charac(2,j)//'.tex')
        write(10,311); write(10,330)'\\ \hline \hline'

        do i=1,40
            i_tmp=i+(j-1)*40
            test=(table(i_tmp,3).lt.100.d0)
            if (minval(table(i_tmp,2:k)).gt.0.1d0) then
                if (test) then
                    write(10,203) table(i+(j-1)*40,1:k)
                else
                    write(10,204) table(i+(j-1)*40,1:k)
                end if
            else
                if (test) then
                    write(10,205) table(i+(j-1)*40,1:k)
                else
                    write(10,206) table(i+(j-1)*40,1:k)
                end if
            end if
            if (mod(i,5).eq.0) write(10,*) '\hline'
        end do
        write(10,120)
        close(10)
    end do

    deallocate(table)

    !
    ! Table of Prandtl Meyer function
    !
    write(6,*)' PRANDTL-MEYER TABLE'
    dM_tmp=(/0.01d0, 0.02d0, 0.5d0, 0.0d0/)
    M_init_tmp =(/1.d0, 3.5d0,  5.d0, 0.d0/)
    M_final_tmp=(/3.5d0, 5.d0, 30.d0, 0.d0/)
    npt(4)=0
    do j=1,3
        write(6,*)'j, ',j,' step ', ((M_final_tmp(j)-M_init_tmp(j))/dM_tmp(j)), M_final_tmp(j),M_init_tmp(j),dM_tmp(j)
        npt(j)=int((M_final_tmp(j)-M_init_tmp(j))/dM_tmp(j))
    end do

    nptM=sum(npt)+1; 
    write(6,*) 'number of points                          : ', nptM

    k=3
    allocate(table(nptM,k))
    write(6,*) 'npt ',npt(:), ' total = ', sum(npt)
    i1=1
    table(i1,1)=M_init_tmp(1)
    write(6,*) 'inital mach                              : ',table(1,1)
    do j=1,3
        i2=i1+npt(j)
        write(31,*)' initial,  final and step Mach',M_init_tmp(j), M_final_tmp(j),dM_tmp(j)
        do i=i1+1,i2
            table(i,1)=table(i-1,1)+dM_tmp(j)
            write(31,'(i4,3x,f12.3,5x,i3)')i,table(i,1),i-i1
        end do
        i1=i2
    end do

    do i=1,nptM
        Mach=table(i,1)
        table(i,2)=omega(Mach)*rad2deg              
        table(i,3)=asin(1.d0/Mach)*rad2deg
    end do

    do j=1,3      ! number of pages, 3 columns
        open(10,form='formatted',file='omega'//charac(2,j)//'.tex')
        write(10,411); 
        write(10,430) ' & '; write(10,430) ' & '; write(10,430)'\\ \hline \hline'

        do i=1,40           ! for every rows
            i_tmp=i+(j-1)*120
            write(10,403) table(i_tmp,1:k),table(i_tmp+40,1:k),table(i_tmp+80,1:k)
            if (mod(i,5).eq.0) write(10,*) '\hline'
        end do
        write(10,120)
        close(10)
    end do

    deallocate(table)


    !*****************************************************************************
    !               FORMATS
    !*****************************************************************************
    100 format((f5.2,' & ',f6.4,' & ',f7.4,' & ',f7.4,' & ',f7.4,' & ', f7.4), &
                f5.3,' & ',f6.4,' & ',f7.4,' & ',f7.4,' & ',f7.4,' & ', f7.4,'\\')
    101 format((f5.3,' & ',f6.4,' & ',f7.4,' & ',f7.4,' & ',f7.4,' & '), &
                5x,  ' & ',6x ,' & ',7x, ' & ',7x, ' & \\')
    102 format(f5.2,' & ',f7.5,' & ',f7.5,' & ',f7.5,' & ',f9.4,' & ',f9.4, '\\')
    103 format(f5.2,' & ',f7.5,' & ',f7.5,' & ',f7.5,' & ',f9.4,' & ',f9.4, '\\')
    104 format(f5.2,' & ',f7.5,' & ',f7.5,' & ',f7.5,' & ',e12.5,' & ',e12.5, '\\')
    105 format(f5.2,' & ',es12.4,' & ',es12.4,' & ',es12.4,' & ',f9.4,' & ',f9.4, '\\')
    106 format(f5.2,' & ',es12.4,' & ',es12.4,' & ',es12.4,' & ',es12.4,' & ',es12.4, '\\')
    110 format('{\small \begin{tabular}{',2(6('|c'),'|')'}', '\hline')
    111 format('{\small \begin{tabular}{',(6('|c'),'|')'}', '\hline')
    120 format('\hline \end{tabular}}')
    130 format('$M$ &  $T/T_i$  &  $P/P_i$ & $\rho / \rho_i$ & $S / S_c$ & $Fa$',a)
    200 format((f5.2,' & ',f7.4,' & ',f7.4,' & ',f7.4,' & ',3(f7.4,' & ')), &
                f5.3,' & ',f7.4,' & ',f7.4,' & ',f7.4,' & ',f7.4,' & ', &
                f9.2,'&',f9.2,'\\')
    201 format((f5.3,' & ',f6.4,' & ',f7.5,' & ',f7.5,' & ',3(f7.5,' & ')), &
                5x,  ' & ',6x ,' & ',7x, ' & ',7x, ' & ',7x,' & ', 7x, ' &   \\')
    203 format(f5.2,' & ',f7.5,3(' & ',f9.4),' & ',f7.5, '\\')
    204 format(f5.2,' & ',f7.5,3(' & ',e12.5),' & ',f7.5, '\\')
    205 format(f5.2,' & ',es12.4,3(' & ',f9.4),' & ',es12.4, '\\')
    206 format(f5.2,' & ',es12.4,3(' & ',es12.4),' & ',es12.4, '\\')
    210 format('{\small \begin{tabular}{',2(7('|c'),'|')'}', '\hline')
    220 format('\hline \end{tabular}}')
    230 format('$M$ &  $T/T_i$  &  $10  P/P_i$ & $10 \rho / \rho_i$ & $S / S_c$',&
                '&  $\mu (^\circ) $ & $\omega (^\circ)$',a)

    ! normal and oblique shocks
    300 format(f5.2,' & ',f6.4,' & ',f7.4,' & ',f7.4,' & ',2(f7.4,' & '), &
                f5.3,' & ',f6.4,' & ',f7.4,' & ',f7.4,' & ',f7.4,' & ', &
                f7.4,'\\')
    301 format((f5.3,' & ',f6.4,' & ',f7.4,' & ',f7.4,' & ',2(f7.4,' & ')), &
                5x,  ' & ',6x ,' & ',7x, ' & ',7x, ' & ',7x,' & \\')
    310 format('{\small \begin{tabular}{',2(6('|c'),'|')'}', '\hline')
    311 format('{\small \begin{tabular}{',6('| c'),'| } \hline')
    320 format('\hline \end{tabular}} \hline')
    330 format('$Mn_0$ &  $Mn_1$  &  $P_1/P_0$ & $ T_1 / T_0$ & $\rho_1/ \rho_0$',&
                '&  $pi_1/pi_0 $',a)
    411 format('{\small \begin{tabular}{',3('| c | c | c | '),' } \hline')
    430 format('$M_0$ &  $\omega \ (^\circ)$  &  $\mu \ (^\circ)$ ',a)
    403 format(2(f5.2,' & ',f10.2,' & ',f10.2, ' & '), &
            f5.2,' & ',f10.2,' & ',f10.2,2x,' \\')

end subroutine tables_livre

!*****************************
function charac(N_in,i_enter)
!*****************************
!> @details
!!  implemented by C Airiau, march 2010
!!  from Burkardt subroutines
!
!! conversion of an integer to a character
!  example : 123 --> '123'
! written on Nin digits
    implicit none
    integer, intent(in)::i_enter,N_in

    character(len=N_in)::charac
    integer i,k,n
    ! n=int(log10(float(i)))
    if (N_in.le.1) then
            write(6,*) 'problem in charact function'
            stop
    endif

    n=N_in-1
    i=i_enter
    do k=n,0,-1
            charac(n-k+1:n-k+1)=char(48+int(i/10**k))
            i=mod(i,10**k)
            !     write(6,*) 'k', k,'<<',carac
    end do
end function charac


!
! IV- critical section , Laval nozzle
!
! 
!*********************************
function Sonic_Mass_flow (pi_tmp,Ti_tmp,Sc) result(res)
!*********************************
!> @details
!! Sonic Mass flow
    use Constante
    implicit none
    real(kind=8)    :: pi_tmp,Ti_tmp,Sc,res
    real(kind=8)    :: alpha,Cste
    alpha=(gam+1.d0)/(2.d0*(gam-1.d0))
    Cste=sqrt(gam/r)*(2.d0/(gam+1.d0))**alpha
    res=Cste*Sc*pi_tmp/sqrt(Ti_tmp)
end function Sonic_Mass_flow


!*********************************
function A_Acritic(Ma) result(res)
!*********************************
!> @details
!! Isentropic relation for area ratio   A/A*
    use Constante
    implicit none
    real(kind=8)    :: Ma,res
    real(kind=8)    :: alpha,omega,Cste
    omega=1.d0+(gam-1.d0)/2.d0*Ma**2
    alpha=(gam+1.d0)/(2.d0*(gam-1.d0))
    Cste=(2.d0/(gam+1.d0))**alpha
    res=Cste/Ma*omega**alpha;
    !A_Acritic=(1.d0/Ma)*((1.d0+(gam-1.d0)/2.d0*Ma**2)/((gam+1.d0)/2.d0))**((gam+1.d0)/(2.d0*(gam-1.d0)));
end function A_Acritic

!************************
function der_A_Acritic(Ma) result(res)
!************************
    !> @details
    !! Isentropic relation for area ratio   A/A*
    !! first derivative w.r.t. Mach
    use Constante
    implicit none
    real(kind=8)    :: Ma,res
    real(kind=8)    :: alpha,omega,Cste
    omega=1.d0+(gam-1.d0)/2.d0*Ma**2
    alpha=(gam+1.d0)/(2.d0*(gam-1.d0))
    Cste=(2.d0/(gam+1.d0))**alpha
    res=Cste*omega**(alpha-1.d0)* (1.d0-1.d0/Ma**2);
end function der_A_Acritic

!**************************************
function der2_A_Acritic(Ma) result(res)
!**************************************
    !> @details
    !! Isentropic relation for area ratio   A/A*
    !! second derivative wrt Mach
    use Constante
    implicit none
    real(kind=8)    :: Ma,res
    real(kind=8)    :: alpha,omega,Cste
    omega=1.d0+(gam-1.d0)/2.d0*Ma**2
    alpha=(gam+1.d0)/(2.d0*(gam-1.d0))
    Cste=(2.d0/(gam+1.d0))**alpha
    res=- Cste*omega**(alpha-2.d0)* (Ma**4*(gam-3)+Ma**2*(5.d0-3.d0*gam-4.d0))/(2.d0*Ma**3);
end function der2_A_Acritic

!*********************************************************************************************
function inverse_A_Acritic(Ma_ref,C_ref,opt_newton,opt_function,iter_max,tol_max) result(res)
!*********************************************************************************************
!> @details
!! inverted function of  A/A_crit=F(Mach) to get a Mach number for A/A_crit given
    implicit none

    real(kind=8),intent(in)     :: Ma_ref,C_ref,tol_max
    integer,intent(in)          :: opt_newton,opt_function,iter_max
    character(len=3)            :: opt
    real(kind=8)                :: res,M_init,M_end,const
    const=C_ref
    if (const.lt.1.d0) then
        stop 'A /A critical < 1 , no solution'
    else if (C_ref-1.d0.lt.1e-8) then
        res=1.d0; return
    end if

    if (Ma_ref.le.1.d0) then
        opt='sub'; M_init=0.01d0; M_end=1.d0
    else
        opt='sup';M_init=2.d0;M_end=6.d0
    end if

    select case (opt_newton)
        case(0)
            res=newton(opt_function,const,M_init,iter_max,tol_max)
        case default
            res=newton_new(opt_function,const,M_init,iter_max,tol_max)
    end select

end function inverse_A_Acritic

!****************************************************
FUNCTION newton(opt,const,x_init,iter_max,tol_max) result(res)
!****************************************************
!> @details
!! Newton method with numerical derivative for  Mach = f(A/A_crit)
    IMPLICIT NONE
    INTEGER,INTENT(IN)      :: opt
    REAL(KIND=8),INTENT(IN) :: const     
    REAL(KIND=8),INTENT(IN) :: x_init    
    REAL(KIND=8),INTENT(IN) :: tol_max   
    INTEGER,INTENT(IN)      :: iter_max  
    REAL(KIND=8)            :: res

    REAL(KIND=8)        :: y0,dx2,x0,dx
    REAL(KIND=8)        :: x1,x2,y1,y2,errx,erry
    INTEGER             :: iter
    REAL(KIND=8),PARAMETER :: epsil=1d-9   
    REAL(KIND=8)        :: tol
    x0=x_init
    tol=1.d0
    iter=0

    ! calculus of the numerical derivative
    !y0=(A_Acritic(x0)-const)**2
    y0=f_analytique(opt,0,x0,const)
    !write(6,*)'x0=',x0,' y0=',y0,'cons= ',const,'A/A* =',A_Acritic(x0)
    if (x0.ne.0.d0) then 
        dx=x0*epsil
    else
        dx=epsil
    end if

    ! PRINT'(a,4x,10(a,12x))','iter', 'tol',' x ', 'errx','erry'
    DO WHILE ((iter.le.iter_max).and.(tol.ge.tol_max))
        iter=iter+1
        x1=x0+dx
        ! y1=(A_Acritic(x1)-const)**2
        y1=f_analytique(opt,0,x1,const)
        dx2=-y0*dx/(y1-y0)
        x2=x0+dx2
        !y2=(A_Acritic(x2)-const)**2
        y2=f_analytique(opt,0,x2,const)
        errx=abs(dx2)
        erry=abs(y2)
        x0=x2; y0=y2
        tol=erry
        ! PRINT'(i4,5(e23.15,2x))',iter,tol,x0,errx,erry
    END DO

    !write(6,*)'Newton , iter :',iter
    IF (iter.gt.iter_max) THEN
        write(6,*) 'no convergence : stop'
        STOP
    END IF

    !write(6,*) 'Solution : x=  ',x0, 'f(x) = ', y0
    res=x0
END FUNCTION newton


!****************************************************
FUNCTION newton_new(opt,const,x_init,iter_max,tol_max) result(res)
!****************************************************
!> @details
!! Newton method with analytical derivative for  Mach = f(A/A_crit)
!!  
    IMPLICIT NONE
    INTEGER,INTENT(IN)      :: opt
    REAL(KIND=8),INTENT(IN) :: const   
    REAL(KIND=8),INTENT(IN) :: x_init    
    REAL(KIND=8),INTENT(IN) :: tol_max  
    INTEGER,INTENT(IN)      :: iter_max  
    REAL(KIND=8)            :: res

    REAL(KIND=8)        :: y0,dx2,x0
    REAL(KIND=8)        :: x2,y2,errx,erry
    INTEGER             :: iter
    REAL(KIND=8)        :: tol
    x0=x_init
    tol=1.d0
    iter=0

    ! analytical derivative
    y0=f_analytique(opt,0,x0,const)
    DO WHILE ((iter.le.iter_max).and.(tol.ge.tol_max))
        iter=iter+1
        dx2=-y0/f_analytique(opt,1,x0,const)
        x2=x0+dx2
        y2=f_analytique(opt,0,x2,const)
        errx=abs(dx2)
        erry=abs(y2)
        x0=x2; y0=y2
        tol=erry
        ! PRINT'(i4,5(e23.15,2x))',iter,tol,x0,errx,erry
    END DO
    !write(6,*)'Newton new , iter :',iter

    IF (iter.gt.iter_max) THEN
        write(6,*) 'no convergence : stop'
        STOP
    END IF

    !write(6,*) 'Solution : x=  ',x0, 'f(x) = ', y0
    res=x0
END FUNCTION newton_new


!***********************************
function f_analytique(opt,f_type,x,const)  result(f)
!***********************************
!> @details  function A/A_crit=F(Mach), f_type = 0 : function, f_type = 1 : derivative
!!

    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN)    :: x,const
    INTEGER,INTENT(IN)          :: opt,f_type
    REAL(KIND=8)                :: f,tmp

    SELECT CASE(opt)
    CASE(0)
        if (f_type.eq.0) then
            f=(A_Acritic(x)-const)**2
        else
            f=2.d0*der_A_Acritic(x)*(A_Acritic(x)-const)
        end if
    CASE(1)
        if (f_type.eq.0) then
            f=A_Acritic(x)-const
        else
            f=der_A_Acritic(x)
        end if
    CASE(2)
        if (f_type.eq.0) then
            f=(1.d0/const-1.d0/A_Acritic(x))**2
        else
            tmp=A_Acritic(x)
            f=2.d0*der_A_Acritic(x)/tmp**2*(1.d0/const-1.d0/tmp)
        end if
    CASE DEFAULT
        f=0.d0
    END SELECT
end function f_analytique


!************************
subroutine plot_A_Acritic
!************************
!> @details  save function A/A_critical in a file
!!
    implicit none
    integer,parameter           :: n=201
    real(kind=8)                :: Ma,dMa,tmp
    integer                     :: i
    open(1,form='formatted',file='AoverAc.dat')
    write(1,*)'# M,    F=A/A*,     A*/A,  dF/dM,    d^2F/dM^2'
    dMa=5.d0/real(n-1,kind=8)
    do i=1,n
        Ma=dMa*real(i,kind=8)
        tmp=A_Acritic(Ma)
        write(1,'(5(e15.8,3x))')Ma,tmp,1.d0/tmp,der_A_Acritic(Ma),der2_A_Acritic(Ma)
    end do
    close(1)
end subroutine plot_A_Acritic


!*********************************
function der_Fanno(Ma) result(res)
!*********************************
!> @details
!! derivative of the Fanno function w.r.t.  Mach number
    use Constante
    implicit none
    real(kind=8)    :: Ma,res
    real(kind=8)    :: omega
    omega=1.d0+(gam-1.d0)/2.d0*Ma**2
    res=-2.d0*(1.d0-Ma**2)/(gam*Ma**3*omega)
end function der_Fanno


!*********************************
function Fanno(Ma) result(res)
!*********************************
!> @details
!! Fanno function w.r.t. Mach number 
    use Constante
    implicit none
    real(kind=8)    :: Ma,res
    real(kind=8)    :: omega
    omega=1.d0+(gam-1.d0)/2.d0*Ma**2
    !res=-(gam+1.d0)/(2.d0*gam)*log(Ma**2/omega)-1.d0/(gam*Ma**2)
    res=(1.d0-Ma**2)/(gam*Ma**2)+(gam+1)/(2.d0*gam)*log(0.5d0*(gam+1.d0)*Ma**2/omega)
end function Fanno


!*********************************
function entropy_Fanno(Ma) result(res)
!*********************************
!> @details
!! s=f(Mach) for the Fanno problem
    use Constante
    implicit none
    real(kind=8)    :: Ma,res
    real(kind=8)    :: omega,tmp
    omega=1.d0+(gam-1.d0)/2.d0*Ma**2
    tmp=r*(gam+1)/(2.d0*(gam-1))*log((gam+1)/2.d0)
    res=r*(log(Ma)-(gam+1)*log(omega)/(2*(gam-1)))+tmp
end function entropy_Fanno


!*********************************
function der_entropy_Fanno(Ma) result(res)
!*********************************
!> @details
!! ds/dMach = d f(Mach) / d Mach for the Fanno problem
    use Constante
    implicit none
    real(kind=8)    :: Ma,res
    real(kind=8)    :: omega
    omega=1.d0+(gam-1.d0)/2.d0*Ma**2
    res=r*(1.d0-Ma**2)/(Ma*omega)
end function der_entropy_Fanno


!*********************************
function Sonic_Length(Ma) result(res)
!*********************************
!> @details
!! Fanno : get the sonic length
    use Constante
    implicit none
    real(kind=8)    :: Ma,res
    real(kind=8)    :: omega
    omega=1.d0+(gam-1.d0)/2.d0*Ma**2
    res=(1.d0-Ma**2)/(gam*Ma**2)+(gam+1)/(2.d0*gam)*log(0.5d0*(gam+1.d0)*Ma**2/omega)
end function Sonic_Length


!************************
subroutine plot_Fanno
!************************
!> @details
!! save Fanno function in a file 
    implicit none
    integer,parameter           :: n=201
    real(kind=8)                :: Ma,dMa,tmp,Cp,Ma_init=0.1,Ma_end=5
    integer                     :: i
    open(1,form='formatted',file='Fanno.dat')
    write(1,*)'# M,    F,      dF/dM,    s/Cp    ds/dM/Cp          4 fmoy L*/D'
    dMa=(Ma_end-Ma_init)/real(n-1,kind=8)
    Cp=gam*r/(gam-1.d0)
    do i=1,n
        Ma=dMa*real(i,kind=8)+Ma_init
        tmp=Fanno(Ma)
        write(1,'(6(e15.8,3x))')Ma,Fanno(Ma),der_Fanno(Ma),entropy_Fanno(Ma)/Cp, &
                der_entropy_Fanno(Ma)/Cp,Sonic_Length(Ma)
    end do
    close(1)
end subroutine plot_Fanno

!**************************************************
subroutine sonic_ratio_Fanno(Ma,r_T,r_P,r_rho,r_Pi)
!**************************************************
!> @details
!! display ratio of  T, P, rho and Pi with respect to critical values
    implicit none
    real(kind=8),intent(in)             :: Ma
    real(kind=8),intent(out)            :: r_T,r_P,r_rho,r_Pi
    real(kind=8)                        :: omega,tmp,alpha
    omega=1.d0+(gam-1.d0)/2.d0*Ma**2
    tmp=0.5d0*(gam+1.d0)/omega
    r_T=tmp
    r_P=1.d0/Ma*sqrt(tmp)
    alpha=-(gam+1.d0)/(2.d0*(gam-1.d0))
    r_rho=1.d0/(Ma*sqrt(tmp))
    r_Pi=1.d0/Ma*tmp**alpha
        write(6,*)'T/T*                                       :  ', r_T
        write(6,*)'P/P*                                       :  ', r_P
        write(6,*)'rho/rho*                                   :  ', r_rho
        write(6,*)'Pi/Pi*                                     :  ', r_Pi
end subroutine sonic_ratio_Fanno


!**************************************************
function inverse_Fanno(x_init,RHS)  result(Mach)
!**************************************************
!> @details 
!! inverted Fanno function : Mach = f(Fanno)
! if x_init < 1 : subsonic case
! if x_init > 1 : supersonic case
!
! use analytical derivative
    implicit none
    real(kind=8),intent(in)  :: RHS     !< Fanno function value
    real(kind=8),intent(in)  :: x_init  !< first guess of Mach number
    real(kind=8)             :: Mach
    real(kind=8),parameter   :: tol_max=1.d-12
    integer,parameter        :: iter_max=25

    REAL(KIND=8)            :: y0,dx2,x0
    REAL(KIND=8)            :: x2,y2,errx,erry
    INTEGER                 :: iter
    REAL(KIND=8)            :: tol
    x0=x_init
    tol=1.d0
    iter=0

    y0=Fanno(x0)-RHS
    !write(6,*)'x0=',x0,' y0=',y0,'cons= ',const,'A/A* =',A_Acritic(x0)
    write(6,'(a,9x,10(a,20x))') ' iter', 'tol',' M  ', 'err Mach','err F'
    DO WHILE ((iter.le.iter_max).and.(tol.ge.tol_max))
        iter=iter+1
        dx2=-y0/der_Fanno(x0)
        x2=x0+dx2
        y2=Fanno(x2)-RHS
        errx=abs(dx2)
        erry=abs(y2)
        x0=x2; y0=y2
        tol=erry
        write(6,'(i4,5(e23.15,2x))'),iter,tol,x0,errx,erry
    END DO
    !write(6,*)'Newton new , iter :',iter

    IF (iter.gt.iter_max) THEN
        write(6,*) 'no convergence : stop'
        STOP
    END IF

    !write(6,*) 'Solution : x=  ',x0, 'f(x) = ', y0
    Mach=x0
    if (Mach.lt.0.d0) then
        write(6,*)'convergence to a bad value'
        write(6,*)'Change the initial guessed Mach number'
        stop
    end if

end function inverse_Fanno


!*********************************
function qm_Mach(Mach) result(res)
!*********************************
!> @details
!! part of the Mach dependence in the mass flow rate formula
    implicit none
    real(kind=8)    :: Mach,res
    res=Mach*(1.d0+(gam-1.d0)/2.d0*Mach**2)**(-(gam+1.d0)/2.d0/(gam-1.d0))
end function qm_Mach

!************************************
function Mass_Flow_Rate(S,Pi,Ti,Mach) result(res)
!************************************
!> @details
!! mass flow rate in a 1D channel
    implicit none
    real(kind=8),intent(in) :: S, Pi,Ti,Mach
    real(kind=8)            :: res
    res=S*sqrt(gam/r)*Pi/sqrt(Ti)*qm_Mach(Mach)
end function Mass_Flow_Rate


!**********************************************************
function Rayleigh(M,opt) result(F)
!**********************************************************
!> @details
!! fonction de Rayleigh et sa dérivée
!> opt = 1: Rayleigh function
!! opt = 2: Derivative of Rayleigh function
    implicit none
    real(kind=8),intent(in)         :: M    !< Mach number
    integer,intent(in)              :: opt  !< = 1 : function, = 2 : derivative
    real(kind=8)                    :: F

    select case (opt)
    case(1)
        F=2.d0*(gam+1.d0)*M**2/(1.d0+gam*M**2)**2*(1.d0+(gam-1.d0)/2.d0*M**2)
    case default
        F=-4.d0*M*(gam+1.d0)*(M**2-1.d0)/(1.d0+gam*M**2)**3
    end select
end function Rayleigh


!**************************************************
function inverse_Rayleigh(x_init,RHS)  result(Mach)
!**************************************************
!> @details
!! inverted  Rayleigh function  Mach = f(Rayleigh function)
!! if x_init < 1 : subsonic case; 
!! if x_init > 1 : supersonic case
    implicit none
    real(kind=8),intent(in)  :: RHS     !< Rayleigh function value
    real(kind=8),intent(in)  :: x_init  !< first guess of Mach number
    real(kind=8)             :: Mach
    real(kind=8),parameter   :: tol_max=1.d-12
    integer,parameter        :: iter_max=25

    REAL(KIND=8)            :: y0,dx2,x0
    REAL(KIND=8)            :: x2,y2,errx,erry
    INTEGER                 :: iter
    REAL(KIND=8)            :: tol
    x0=x_init
    tol=1.d0
    iter=0

    !use of analytical derivative in Newton method
    y0=Rayleigh(x0,1)-RHS
    !write(6,*)'x0=',x0,' y0=',y0,'cons= ',const,'A/A* =',A_Acritic(x0)
    write(6,'(a,9x,10(a,20x))'),' iter', 'tol',' M  ', 'err Mach','err F'
    DO WHILE ((iter.le.iter_max).and.(tol.ge.tol_max))
        iter=iter+1
        dx2=-y0/Rayleigh(x0,2)
        x2=x0+dx2
        y2=Rayleigh(x2,1)-RHS
        errx=abs(dx2)
        erry=abs(y2)
        x0=x2; y0=y2
        tol=erry
        write(6,'(i4,5(e23.15,2x))'),iter,tol,x0,errx,erry
    END DO
    !write(6,*)'Newton new , iter :',iter

    IF (iter.gt.iter_max) THEN
        write(6,*) 'no convergence : stop'
        STOP
    END IF

    !write(6,*) 'Solution : x=  ',x0, 'f(x) = ', y0
    Mach=x0
    if (Mach.lt.0.d0) then
        write(6,*)'Convergence to a bad value'
        write(6,*)'change the initial Mach number '
        stop
    end if
end function inverse_Rayleigh



!**********************************************************
subroutine Rayleigh_flow(q)
!**********************************************************
!> @details
!! display the ratios to critical state for the Rayleigh flow
    implicit none
    type(state),intent(inout)           :: q
    real(kind=8),parameter              :: alpha=(gam+1.d0)/2.d0
    real(kind=8),parameter              :: beta=gam/(gam-1.d0)
    real(kind=8)                        :: omega,Ma
    Ma=q%M
    omega=1.d0+(gam-1.d0)/2.d0*Ma**2
    q%P=(gam+1.d0)/(1.d0+gam*Ma**2)
    q%T=Ma**2*q%P**2
    q%rho=1.d0/(Ma**2*q%P)
    q%u=1.d0/q%rho
    q%Ti=q%T*omega/alpha
    q%Pi=q%P*(omega/alpha)**beta
    q%s=beta*log(q%T)-log(q%P)
    write(6,*)'T/T*                                       :  ', q%T
    write(6,*)'P/P*                                       :  ', q%P
    write(6,*)'rho/rho*                                   :  ', q%rho
    write(6,*)'Pi/Pi*                                     :  ', q%Pi
    write(6,*)'Ti/Ti*                                     :  ', q%Ti
    write(6,*)'(s-s*)/r                                   :  ', q%s
    write(6,*)'(s-s*)/Cp                                  :  ', q%s/beta
end subroutine Rayleigh_flow


!**************************************************
subroutine plot_Rayleigh_entropy
!**************************************************
!> @details
!! (s-s_c)/Cp function of Ma
    implicit none
    integer,parameter           :: n=981
    real(kind=8)                :: Ma,dMa,Ma_init=0.1,Ma_end=5
    real(kind=8)                :: r_T,r_P,r_u,r_rho,r_Ti,r_Pi
    real(kind=8)                :: tmp,omega
    real(kind=8),parameter      :: alpha=(gam+1.d0)/2.d0
    real(kind=8),parameter      :: beta=gam/(gam-1.d0)
    integer                     :: i
    open(1,form='formatted',file='Rayleigh.dat')
    write(1,100)
    dMa=(Ma_end-Ma_init)/real(n-1,kind=8)


    do i=1,n
        Ma=dMa*real(i-1,kind=8)+Ma_init
        omega=1.d0+(gam-1.d0)/2.d0*Ma**2
        r_P=(gam+1.d0)/(1.d0+gam*Ma**2)
        r_T=Ma**2*r_P**2
        r_rho=1.d0/(Ma**2*r_P)
        r_u=1.d0/r_rho
        r_Ti=r_T*omega/alpha
        r_Pi=r_P*(omega/alpha)**beta
        tmp=2.d0*log(Ma)+alpha*log(r_P)
        write(1,'(8(e15.8,3x))')Ma,tmp,r_P,r_T,r_u,r_Ti,r_Pi,r_rho
    end do
    close(1)
    100 format('#',5x, 'M',14x,'(s-sc)/Cp',12x,'P/Pc',13x,'T/TC',14x,'u/uc',14x,'Ti/Tic',12x,&
            'Pi/Pic',12x,'rho/rhoc')

end subroutine plot_Rayleigh_entropy


!********************************************************************
include 'tables_perso.f90'
!********************************************************************

end module mod_CAT


