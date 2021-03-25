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
!! MODULE :  base to select the exercises or problem to solve
!! 
!! @author C. Airiau
!!
module mod_etude
use mod_cat
use mod_shock_interaction
use exercices_chap10
use exercices_chap11
use exercices_chap12
contains


!***************************
subroutine liste_exercices()
!***************************
!>@details 
!! list of exercises with solution
!!
!! the list must be included in the file  "exercise_choice.in"
!!
!! the solution is then printed into an output file whose name refers to the exercise reference

implicit none
integer                              :: chap,exo
logical                              :: file_exists
character(len=18),parameter          :: fichier='exercise_choice.in'
integer                              :: nc,i

print*,'Choice of the running exercise'
INQUIRE(FILE=fichier,exist=file_exists)
if (file_exists) then
    open(20, form='formatted',file='comment_exos.out')
    open(1,form='formatted',file=fichier)
    read(1,*);read(1,*);read(1,*)
    read(1,*) nc
    print*,"number of exercises which will be solved : ",nc
    do i=1,nc
        read(1,*) chap,exo
        write(20,*)'Chapter :', chap,' exercise: ',exo
        if (exo.lt. 10) then
            write(20,100)chap,exo
        else
            write(20,101)chap,exo
        endif
        100 format('output in the file : ','Solution_',i2,'-',i1,'.out')
        101 format('output in the file : ','Solution_',i2,'-',i2,'.out')
        
        select case(chap)
        case(10)
            select case(exo)
                case(1)
                    call Exercice_10_1    ! Normal shock wave / Calcul d'un choc droit
                case(2)
                    call Exercice_10_2    ! Normal shock wave calculus / Calcul pratique d'un choc droit
                case(10)
                    call Exercice_10_10   ! Laval Nozzle  / Tuyère de Laval
                case(3)
                    call Exercice_10_3    ! Subsonic Fanno's flow / Ecoulement de Fanno subsonique
                case(4)
                    call Exercice_10_4    ! Supersonic Fanno's flow / Ecoulement de Fanno supersonique
                case(5)
                    call Exercice_10_5    ! Rayleigh's flow / Ecoulement de Rayleigh
                case(6)
                    call Exercice_10_6    ! Rayleigh's flow with combustion /  Ecoulement de Rayleigh + combustion
                case(7)
                    call Exercice_10_7    ! Pitot tube in subsonic regime / Pitot subsonique
                case(8)
                    call Exercice_10_8    !  Pitot tube in supersonic regime / Pitot supersonique
                case(9)
                    call Exercice_10_9    ! Regimes in nozzles  / Régimes dans une tuyere
                case default
                    write(20,*) "bad reference to exercise"
            end select

        case(11)
            select case(exo)
                case(1)
                    call Exercice_11_1    ! Flat plate in incidence / plaque en incidence
                case(2)
                    call Exercice_11_2    ! Diamond airfoil /profil losangique   
                case(3)
                    call Exercice_11_3    ! methode choc expansion
                case(4)
                    call Exercice_11_4    ! Shocks in a channel / Shock dans un canal
                case(5)
                    call Exercice_11_5    ! Two obliques shock wave interactions /Interaction de 2 chocs obliques
                case(6)
                    call Exercice_11_6    ! Backward facing step  (simple program) / Marche descendante (programmation basique)
                case(61)
                    call Exercice_11_61   ! Backward facing step (2nd program) / Marrche descendante (2nd programmation)
                case(7)
                    call Exercice_11_7    ! Shock - isobar line interactions of a fluid at rest  / interaction choc - ligne isobare avec fluide au repos
                case(8)
                    call Exercice_11_8    ! Expansion - isobar line interactions of a fluid at rest / interaction détente - ligne isobare avec fluide au repos
                case(9)
                    call Exercice_11_9    ! Expansion - wall interaction /interaction détente - paroi
                case(10)
                    call Exercice_11_10   ! 2 upstream flow  separated by a isobar line  - oblique shock interaction /deux écoulements amont + Ligne de glissement + choc oblique
                case(11)
                    call Exercice_11_11   ! Interaction of shock waves in a nozzle outlet / tuyère : interaction d'ondes de chocs en sortie
                case(12)
                    call Exercice_11_12   ! Interaction of a centered expansion in a nozzle outlet / tuyère : interaction d'un faisceau de détente en sortie
                case default
                    write(20,*) "bad reference to exercise"
            end select
        case(12)
            select case(exo)
                case(3)
                    call Exercice_12_3     ! Shock tube application / application du tube à choc.
                case default
                    write(20,*)"bad reference to exercise"
            end select
        end select
    end do

else
    print*,'file', fichier," not found"
end if
close(20)
call system("more comment_exos.out")

end subroutine liste_exercices

!*****************
subroutine run_CAT
!*****************
!> @details
!! SUBROUTINE to run the study with the choice of options
!! They can be enter on a terminal, by command line or in a input file
! 
implicit none
real(kind=8) :: sigma_atm,delta_atm,theta_atm,alt
character(len=50) :: text

call driver_aerodynamics(ask,ichoix,M1,teta,ans,F)
call set_menu
if (ichoix== -2) then
    call help
    stop 'end of execution'
else if (ichoix == -1) then 
    call display_menu
    write(6,*) 'enter your choice:'; read*,ichoix 
else 
    write(6,*) 'menu choice : ' , ichoix
    write(6,*) menu_title(ichoix)
end if

select case (ichoix)
case(0)

    write(6,*)
    if (ask) then 
        call question('Enter upstream Mach number',M1)
    else
        call display('upstream Mach number',M1)
    end if
    M2=mach_downstr(M1)
    call display('downstream Mach M2',M2)
    call display('Ratio T1 /  Ti', T_Ti(M1))
    call display('Ratio T2 /  Ti', T_Ti(M2))
    call display('Density ratio rho2 / rho1',rho2_rho1(M1))
    call display('Pressure ratio p2 / p1',P2_P1(M1))
    call display('Total pressure ratio pi2 / pi1',Pi2_Pi1(M1))

case(1)

    write(6,*)
    if (ask) then
        call question('Enter upstream Mach number',M1)
        call question('Enter the deviation angle in deg.',teta)
    else
        call display('upstream Mach number',M1)
        call display('deviation angle in deg.',teta)
    end if
    teta=teta*deg2rad
    call Newton_shock_angle(sigmac,teta,M1)
    call display('Shock angle in deg.',sigmac*rad2deg)
    Mn1=M1*sin(sigmac)
    call display('Upstream normal Mach   Mn1',Mn1)
    Mn2=mach_downstr(Mn1)
    M2 = Mn2/sin(sigmac-teta)
    call display('Downstream normal Mach Mn2', Mn2)
    call display('downstream Mach M2', M2)
    call display('temperature ratio T1 / Ti', T_Ti(M1))
    call display('temperature ratio T2 /  Ti', T_Ti(M2))
    call display('Density ratio rho2 / rho1', rho2_rho1(Mn1))
    call display('Pressure ratio p2 / p1', P2_P1(Mn1))
    call display('Total pressure ratio pi2 / pi1', Pi2_Pi1(Mn1))
    Kp=2.d0/(gam*M1**2)*(P2_P1(Mn1)-1.d0)
    call display('Kp', Kp)

case(2)
    if (ask) then
        call question('Enter upstream Mach number',M1)
        text='teta > 0 => expansion, teta < 0 => compression '
        call question_help('Enter the deviation angle in deg. ',text, teta)
    else
        call display('upstream Mach number',M1)
        call display(' deviation angle in deg.',teta)
    end if
    teta=teta*deg2rad
    om1=omega(M1)
    om2=om1+teta
    M2=invomega(om2)
    call display('Upstream P.M. function omega in deg.', om1*rad2deg)
    call display('Downstream P.M. function omega in deg.', om2*rad2deg)
    call display('Downstream Mach', M2)
    call display('Pressure ratio', P_Pi(M2)/P_Pi(M1))

case(3)
    if (ask) then
        call question('Enter Mach number',M1)
    else
        call display('upstream Mach number',M1)
    end if
    om1=omega(M1)
    call display('P.M. function omega in deg.', om1*rad2deg)
    
case(4)
    if (ask) then
        call question('Enter P.M. function omega in deg.',teta)
        call display('omega angle',teta)
    end if
    om1=teta
    om1=om1*deg2rad
    M1=invomega(om1)
    call display('Mach number',M1)

case(5)
    if (ask) then
        call question('Enter Mach number : ',M1)
    else
        call display('upstream Mach number',M1)
    end if
    call display('Temperature ratio T / Ti', T_Ti(M1))
    call display('Pressure ratio P / Pi ', P_Pi(M1))
    call display('Density ratio rho / rhoi ', rho_rhoi(M1))

case(6)
    call omega_plot 

case(7)
    !     output of omega(Mach)
    open(1,form='formatted',file='omega.dat')
    xm(1)=1.
    do i=1,ndim
        ym(i)=omega(xm(i))*rad2deg
        xm(i+1)=xm(i)+0.005d0
    end do

    do i=1,ndim,5
    write(1,200)(xm(i-1+j),ym(i-1+j),j=1,5)
    if (abs(xm(i)-int(xm(i))-0.975d0).lt.1d-6)  write(1,*)
    end do
    200 format(4(f6.3,' & ', f8.3,' & '),f6.3, ' & ', f8.3,'\\')
    close(1)

case(8)
    call epicycloide

case(9)
    call shock_table

case(10)
    if (ask) then
        M1=0
    else
        call display('choosen case',M1)
    end if
    ! Here M1 is not the Mach number but the index of the case!
    ! If we don't not set to 0.
    call main_shock_interaction(int(M1))

case(90)
        call examen_mai_2014
case(91)
        call examen_juin_2014
case(12)
        call tables_livre
case(13)
    write(6,*)
    if (ask) then
        call question('Enter the Mach number',M1)
    else
        call display('Mach number ',M1)
    end if
    call display('A / A critical', A_Acritic(M1))
case(14)
    write(6,*)
    if (ask) then
        call question('Enter A / A critical',ans)
        call question('enter initial Mach guess (supersonic > 1, subsonic < 1)',M1)
    else
        call display('A / A_critical',ans)
        if (M1.le.1) then
            write(6,*)'subsonic solution'
        else
            write(6,*)'supersonic solution'
        end if
    end if
    M2=inverse_A_Acritic(M1,ans,1,1,45,1.d-12) 
    call display('corresponding Mach number', M2)
case(15)
    write(6,*) ' The A/Acritical function is in the file AoverAc.dat'
    call plot_A_Acritic
case(16)
    write(6,*) menu_title(ichoix)
    if (ask) then
        call question('Enter the Mach number',M1)
    else
        call display('Mach',M1)
    end if
    call display('Fanno function', Fanno(M1))
    call display('entropy S(M) - S_crit', entropy_Fanno(M1))
    call display(' 4 f_moy L* /D ', Sonic_Length(M1))
    call sonic_ratio_Fanno(M1,tmp1,tmp2,tmp3,tmp4)

case(17)
    write(6,*) menu_title(ichoix)
    if (ask) then
        call question('Enter the FANNO function value',F)
        call question('Enter the estimated Mach  (  subsonic M <1,  supersonic M >1)',M1)
        if (M1.eq.1.d0) stop 'bad value of estimated Mach'
    else
        call display('Fanno function',F)
        if (M1.le.1.d0) then
            write(6,*)'subsonic solution'
        else
            write(6,*)'supersonic solution'
        end if

    end if

    M2=inverse_Fanno(M1,F) 
    call display('corresponding Mach', M2)
case(18)
    write(6,*) menu_title(ichoix)
    call plot_Fanno
case(19)
    write(6,*) menu_title(ichoix)
    if (ask) then
        call question('Enter the Mach number',M1)
    else
        call display('Mach number',M1)
    end if
    q%M=M1;
    call Rayleigh_flow(q)
case(20)
    write(6,*) menu_title(ichoix)
    call plot_Rayleigh_entropy
case(21)
    write(6,*) menu_title(ichoix)
    if (ask) then
        call question('Enter the Mach',M1)
    else
        call display('Mach number',M1)
    end if
    call display('Rayleigh function', Rayleigh(M1,1))
case(22)
    write(6,*) menu_title(ichoix)
    if (ask) then
        call question('Enter the Rayleigh function : ',F)
        call question('Enter the estimated Mach  (  subsonic M <1,  supersonic M >1) :',M1)
        if (M1.eq.1.d0) stop 'bad value of estimated Mach'
    else
        call display('Rayleigh function',F)
        if (M1.le.1.d0) then
            write(6,*)'subsonic solution'
        else
            write(6,*)'supersonic solution'
        end if
    end if
    if ((M1.gt.1.d0).and.(F.lt.1.d0-1.d0/gam**2)) stop 'Mauvaise valeur de F'
    M2=inverse_Rayleigh(M1,F) 
    call display('corresponding Mach number', M2)

case(23)
    write(6,*) menu_title(ichoix)
    call tables_perso

case(24)
    write(6,*) menu_title(ichoix)
    call question('Enter  the altitude in km',alt)
    call atmosphere(alt,sigma_atm,delta_atm,theta_atm)
    call display('density rho', rho_atm*sigma_atm)
    call display('pressure p ', P_atm*delta_atm)
    call display('temperature T', T_atm*theta_atm)

case(27)
    write(6,*)'Exercise of the book AERODYNAMIQUE FONDAMENTALE Chap. 10, 11 et 12'
    !include 'liste_exercices.f90'
    call liste_exercices()

case(28)
    write(6,*) 'ENSEEIHT TD 1'
    call TD1_N7

case(29)
    write(6,*) 'ENSEEIHT TD 2'
    call TD2_N7

case(30)
    write(6,*) 'correction N7, fevrier 2017'
    call examen_Fev2017()

case(99) 
        call read_case
        call solve_case
        call display_outputs
        deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,upstr_zone)
end select

write(*,*) 'Normal end of the CAT program'

end subroutine run_CAT
!  end of the main program


!******************************
subroutine question(text, valeur)
!******************************
!> @details
!! SUBROUTINE formatted question with an expected value
    character(len=*), intent(in) :: text !< content of the question
    real (kind=8), intent(out)   :: valeur    !< user value, real number
    write(6,100) trim(text)
    read*,valeur
    100 format(a,' ? ')
end subroutine question

!******************************
subroutine question_help(text, help, valeur)
    !******************************
    !> @details
    !! SUBROUTINE formatted question with an expected value, and an help
        character(len=*), intent(in) :: text    !< content of the question
        character(len=*), intent(in) :: help    !< help of the question
        real (kind=8), intent(out)   :: valeur  !< user value, real number
        write(6,100) trim(text)
        write(6,101) trim(help)
        read*,valeur
        100 format(a,' ? ')
        101 format("help : ",a)
    end subroutine question_help

!**********************************
subroutine display(comment, valeur)
!*********************************
    !> @details
    !! SUBROUTINE formatted display of a number with its meaning 
    character(len=*), intent(in) :: comment !< content of the comment
    real (kind=8), intent(in)    :: valeur    !< user value to display, real number
    character(len=45),parameter :: line=' ..........................................'
    n=len(comment)
    write(6,100)adjustl(comment//line(1:45-n)),valeur
    100 format(a45,' : ',f15.6)
end subroutine display



include 'exercices_N7.f90'
include 'examens.f90'

end module mod_etude
