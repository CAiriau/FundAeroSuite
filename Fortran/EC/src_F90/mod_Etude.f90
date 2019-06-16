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
!! MODULE :  résolution des exercices
!! 
!! @author C. Airiau
!!
module mod_etude
use mod_EC
use mod_interaction_choc
use exercices_chap10
use exercices_chap11
use exercices_chap12
contains




!***************************
subroutine liste_exercices()
!***************************
!>@details 
!! liste des exercices corrigés.
!!
!! il faut mettre la liste des exercices choisis dans le fichie "choix_exercices.in"
!!
!! la sortie d'effectue dans un fichier qui porte le nom de l'exercice
 
implicit none
integer                              :: chap,exo
logical                              :: file_exists
character(len=18),parameter          :: fichier='choix_exercices.in'
integer                              :: nc,i

print*,'Choix des exercices en cours'
INQUIRE(FILE=fichier,exist=file_exists)
if (file_exists) then
    open(20, form='formatted',file='commentaires_exos.out')
    open(1,form='formatted',file=fichier)
    read(1,*);read(1,*);read(1,*)
    read(1,*) nc
    print*,"nombre d'exercices résolus : ",nc
    do i=1,nc
        read(1,*) chap,exo
        write(20,*)'Chapitre :', chap,' exercice: ',exo
        if (exo.lt. 10) then
            write(20,100)chap,exo
        else
            write(20,101)chap,exo
        endif
        100 format('sortie dans le fichier : ','Solution_',i2,'-',i1,'.out')
        101 format('sortie dans le fichier : ','Solution_',i2,'-',i2,'.out')
        
        select case(chap)
        case(10)
            select case(exo)
                case(1)
                    call Exercice_10_1    ! Calcul d'un choc droit
                case(2)
                    call Exercice_10_2    ! Calcul pratique d'un choc droit
                case(30)
                    call exo_Laval
                case(3)
                    call Exercice_10_3    ! Ecoulement de Fanno subsonique
                case(4)
                    call Exercice_10_4    ! Ecoulement de Fanno supersonique
                case(5)
                    call Exercice_10_5    ! Ecoulement de Rayleigh
                case(6)
                    call Exercice_10_6    ! Ecoulement de Rayleigh + combustion
                case(7)
                    call Exercice_10_7    ! Pitot subsonique
                case(8)
                    call Exercice_10_8    ! Pitot supersonique
                case(9)
                    call Exercice_10_9    ! Régimes dans une tuyere
                case default
                    write(20,*) "cet exercice  n'existe pas"
            end select

        case(11)
            select case(exo)
                case(1)
                    call Exercice_11_1    ! plaque en incidence
                case(2)
                    call Exercice_11_2    ! profil losangique   
                case(3)
                    call Exercice_11_3    ! methode choc detente
                case(4)
                    call Exercice_11_4    ! Choc dans un canal
                case(5)
                    call Exercice_11_5    ! Interaction de 2 chocs obliques
                case(6)
                    call Exercice_11_6    ! Marche descendante (programmation basique)
                case(61)
                    call Exercice_11_6bis ! Marche descendante (2nd programmation)
                case(7)
                    call Exercice_11_7    ! interaction choc - ligne isobare avec fluide au repos
                case(8)
                    call Exercice_11_8    ! interaction détente - ligne isobare avec fluide au repos
                case(9)
                    call Exercice_11_9    ! interaction détente - paroi
                case(10)
                    call Exercice_11_10   ! deux écoulements amont + Ligne de glissement + choc oblique
                case(11)
                    call Exercice_11_11   ! tuyère : interaction d'ondes de chocs en sortie
                case(12)
                    call Exercice_11_12   ! tuyère : interaction d'un faisceau de détente en sortie
                case default
                    write(20,*) "cet exercice  n'existe pas"
            end select
        case(12)
            select case(exo)
                case(3)
                    call Exercice_12_3     ! application du tube à choc.
                case default
                    write(20,*) "cet exercice  n'existe pas"
            end select
        end select
    end do

else
    print*,'le fichier ', fichier," n'a pas été trouvé"
end if
close(20)
call system("more commentaires_exos.out")

end subroutine liste_exercices

!*****************
subroutine run_EC
!*****************
!> @details
!! SUBROUTINE permettant de lancer l'étude, avec le choix des options
!! soit à entrer en lignes, soit au clavier, soit dans un fichier
! 
implicit none
real(kind=8)                                :: sigma_atm,delta_atm,theta_atm,alt
call driver_aerodynamique(ask,ichoix,M1,teta,ans,F)
call set_menu
if (ichoix== -2) then
    call help
    stop 'fin du programme'
else if (ichoix == -1) then 
    call display_menu
    write(6,*) 'entrer votre choix:'; read*,ichoix 
else 
    write(6,*) 'choix du menu : ' , ichoix
    write(6,*) menu_title(ichoix)
end if

select case (ichoix)
case(0)

     write(6,*)
     if (ask) then 
        write(6,*)'Entrer le MACH amont :'; read*,M1
     else
        write(6,*)'MACH amont                                 :  ',M1
     end if
     M2=Mach_Aval(M1)
     write(6,*)'MACH AVAL  M2                              :  ', M2
     write(6,*)'RAPPORT T1/Ti                              :  ', T_Ti(M1)
     write(6,*)'RAPPORT T2/Ti                              :  ', T_Ti(M2)
     write(6,*)'RAPPORT DES MASSES VOLUMIQUES rho2/rho1    :  ', rho2_rho1(M1)
     write(6,*)'RAPPORT DES PRESSIONS p2/p1                :  ', P2_P1(M1)
     write(6,*)'RAPPORT DES PRESSIONS TOTALES pi2/pi1      :  ', Pi2_Pi1(M1)

case(1)

     write(6,*)
     if (ask) then
        write(6,*)'Entrer le MACH amont : '; read*,M1
        write(6,*)'Entrer l''angle deviation de l''ecoulement en degres :'
        read*,teta
     else
        write(6,*)'MACH amont                                 :  ',M1
        write(6,*)'angle deviation de l''ecoulement en degres :  ',teta
     end if
     teta=teta*deg2rad
     call Newton_Angle_Choc(sigmac,teta,M1)
     write(6,*)'ANGLE DU CHOC en degres                    :  ',sigmac*rad2deg
     Mn1=M1*sin(sigmac)
     write(6,*)'MACH NORMAL AMONT Mn1                      :  ',Mn1
     Mn2=Mach_Aval(Mn1)
     M2 = Mn2/sin(sigmac-teta)
     write(6,*)'MACH NORMAL AVAL  Mn2                      :  ', Mn2
     write(6,*)'MACH AVAL  M2                              :  ', M2
     write(6,*)'RAPPORT T1/Ti                              :  ', T_Ti(M1)
     write(6,*)'RAPPORT T2/Ti                              :  ', T_Ti(M2)
     write(6,*)'RAPPORT DES MASSES VOLUMIQUES rho2/rho1    :  ', rho2_rho1(Mn1)
     write(6,*)'RAPPORT DES PRESSIONS p2/p1                :  ', P2_P1(Mn1)
     write(6,*)'RAPPORT DES PRESSIONS TOTALES pi2/pi1      :  ', Pi2_Pi1(Mn1)
     Kp=2.d0/(gam*M1**2)*(P2_P1(Mn1)-1.d0)
     write(6,*)'Kp                                         :  ', Kp
     
 case(2)
     if (ask) then
        write(6,*)'Entrer le MACH amont :'; read*,M1
        write(6,*)'Entrer l''angle deviation de l''ecoulement en degres : '
        write(6,*)'avec teta > 0 : detente, teta < compression '
        read*,teta
     else
        write(6,*)'MACH amont                                 :  ',M1
        write(6,*)'angle deviation de l''ecoulement en degres :  ',teta
     end if
     teta=teta*deg2rad
     om1=omega(M1)
     om2=om1+teta
     M2=invomega(om2)
     write(6,*)'OMEGA amont en degres                      :  ', om1*rad2deg
     write(6,*)'OMEGA aval  en degres                      :  ', om2*rad2deg
     write(6,*)'MACH aval                                  :  ', M2
     write(6,*)'Rapport des pressions                      :  ', P_Pi(M2)/P_Pi(M1)

 case(3)
     if (ask) then
        write(6,*)'Entrer le MACH :'; read*,M1
     else
        write(6,*)'MACH amont                                 :  ',M1
     end if
     om1=omega(M1)
     write(6,*) 'OMEGA en degres                           ;  ', om1*rad2deg
     
 case(4)
     if (ask) then
        write(6,*)'Entrer OMEGA en degres : '; read*,teta
        write(6,*)'angle omega                                :  ',teta
     end if
     om1=teta
     om1=om1*deg2rad
     M1=invomega(om1)
     write(6,*)'MACH                                       :  ', M1
     
   
 case(5)
    if (ask) then
        write(6,*)'Entrer le MACH : '; read*,M1
    else
        write(6,*)'MACH amont                                 :  ',M1
    end if
    write(6,*)'RAPPORT T/Ti                                :  ', T_Ti(M1)
    write(6,*)'RAPPORT P/Pi                                :  ', P_Pi(M1)
    write(6,*)'RAPPORT rho/rhoi                            :  ', rho_rhoi(M1)

case(6)
   call courbe_omega 

    
case(7)
     !     sortie de la fonction omega(Mach)
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
     call table_choc
case(10)

    if (ask) then
        M1=0
    else
        write(6,*)'cas choisi                                 :  ',M1
    end if
    ! attention M1 n'est pas le nombre de Mach mais le cas considéré!
    ! si on ne sait pas on met 0.
    call main_interaction_choc(int(M1))

case(90)
        call examen_mai_2014
case(91)
        call examen_juin_2014
case(12)
        call tables_livre
case(13)
     write(6,*)
     if (ask) then
        write(6,*)'Entrer le MACH  :'; read*,M1
     else
        write(6,*)'MACH                                       :  ',M1
     end if
     write(6,*)'A / A critique                             :  ', A_Acritic(M1)
case(14)
     write(6,*)
     if (ask) then
        write(6,*)'Entrer A / A critique  :'; read*,ans
        write(6,*)'type de solution :  subsonique (oui : 1, non : autre)'
        read*,M1
     else
        write(6,*)'A / A_critique                             :  ',ans
        if (M1.le.1) then
           write(6,*)'Solution subsonique'
        else
           write(6,*)'Solution supersonique'
        end if
     end if
     M2=inverse_A_Acritic(M1,ans,1,1,45,1.d-12) 
     write(6,*)'Mach correspondant                         :  ', M2
case(15)
    write(6,*) 'la courbe A/Acritique est dans AoverAc.dat'
    call plot_A_Acritic
case(16)
    write(6,*) menu_title(ichoix)
    if (ask) then
        write(6,*)'Entrer le MACH  :'; read*,M1
    else
        write(6,*)'MACH                                       :  ',M1
    end if
    write(6,*)'Fonction de Fanno                          :  ', Fanno(M1)
    write(6,*)'Entropie S(M) - S_crit                     :  ', Entropie_Fanno(M1)
    write(6,*)' 4 f_moy L* /D                             :  ', Sonic_Length(M1)
    call sonic_ratio_Fanno(M1,tmp1,tmp2,tmp3,tmp4)

case(17)
    write(6,*) menu_title(ichoix)
    if (ask) then
        write(6,*)'Entrer la fonction de FANNO  : '; read*,F
        write(6,*)'Entrer le Mach estimé (  subsonique M <1,  supersonique M >1) :'
        read*,M1
        if (M1.eq.1.d0) stop 'Mauvaise valeur du Mach estime'
    else
        write(6,*)'Fonction de Fanno                          :  ',F
        if (M1.le.1.d0) then
           write(6,*)'Solution subsonique'
        else
           write(6,*)'Solution supersonique'
        end if

    end if

     M2=inverse_Fanno(M1,F) 
     write(6,*)'Mach correspondant                         :  ', M2
case(18)
    write(6,*) menu_title(ichoix)
    call plot_Fanno
case(19)
    write(6,*) menu_title(ichoix)
    if (ask) then
        write(6,*)'Entrer le MACH  :'; read*,M1
    else
        write(6,*)'MACH                                       :  ',M1
    end if
    q%M=M1;
    call Rayleigh_flow(q)
case(20)
    write(6,*) menu_title(ichoix)
    call plot_Rayleigh_entropie
case(21)
    write(6,*) menu_title(ichoix)
    if (ask) then
        write(6,*)'Entrer le MACH  :'; read*,M1
    else
        write(6,*)'MACH                                       :  ',M1
    end if
    write(6,*)'Fonction de Rayleigh                       :  ', Rayleigh(M1,1)
case(22)
     write(6,*) menu_title(ichoix)
    if (ask) then
        write(6,*)'Entrer la fonction de RAYLEIGH  : '; read*,F
        write(6,*)'Entrer le Mach estimé (  subsonique M <1,  supersonique M >1) :'
        read*,M1
        if (M1.eq.1.d0) stop 'Mauvaise valeur du Mach estime'
    else
        write(6,*)'Fonction de Rayleigh                     :  ',F
        if (M1.le.1.d0) then
           write(6,*)'Solution subsonique'
        else
           write(6,*)'Solution supersonique'
           
        end if
    end if
    if ((M1.gt.1.d0).and.(F.lt.1.d0-1.d0/gam**2)) stop 'Mauvaise valeur de F'
    M2=inverse_Rayleigh(M1,F) 
    write(6,*)'Mach correspondant                         :  ', M2

case(23)
     write(6,*) menu_title(ichoix)
     call tables_perso
case(24)

     write(6,*) menu_title(ichoix)
     write(6,*)'Entrer  altitude en km: '; read*,alt
     call atmosphere(alt,sigma_atm,delta_atm,theta_atm)
    write(6,*)' rho                                       :  ', rho_atm*sigma_atm
    write(6,*)' P                                         :  ', P_atm*delta_atm
    write(6,*)' T                                         :  ', T_atm*theta_atm
case(27)
    write(6,*)'Exercice du livre AERODYNAMIQUE FONDAMENTALE Chap. 10, 11 et 12'
    !include 'liste_exercices.f90'
    call liste_exercices()
   case(28)
    write(6,*) 'TD 2'
    call TD2_N7
case(29)
    write(6,*) 'TD 1'
    call TD1_N7
case(30)
    write(6,*) 'correction N7, fevrier 2017'
    call examen_Fev2017
case(99) 
        call read_case
        call solve_case
        call sorties
        deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn,zone_amont)
end select

write(*,*) 'FIN DU PROGRAMME'

end subroutine run_EC
!  FIN DU PROGRAMME PRINCIPAL


include 'exercices_N7.f90'
include 'examens.f90'
end module mod_etude
