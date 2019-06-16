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
    real(kind=8),parameter        :: P_atm=101325.d0    ! pression de reference
    real(kind=8),parameter        :: T_atm=288.15d0     ! temperature de reference
    real(kind=8),parameter        :: rho_atm=P_atm/(T_atm*r)  ! masse volumique de refence
    logical,parameter             :: fichier=.true. ! mettre true pour une sortie dans un fichier, false pour une sortie à l ecran.
end Module Constante


 
!> MODULE type_def : contient des définitions 
!======================
Module type_def
!======================
!> type pour introduire un état aérodynamique (peu utilisé finalement).
type state
real(kind=8)        :: M,T,P,rho,u,Ti,Pi,s
end type state

!> type pour définir une option d entree.
type input
integer             :: n
character(len=100)  :: title=''
end type

end Module type_def

!> MODULE mod_EC  : contient l'ensemble des procédures utilisées pour l'aérodynamique 
!!                  compressible.
!!
!!                 calcul des ondes de chocs ou des detentes de Prandtl-Meyer
!!
!! HISTORIQUE :
!!
!!      - Mai 2001  : premier programme sur les chocs et les détentes  
!!      - Mars 2010 : passage du f77 au f90
!!      - Mai 2014 à Aout 2015 : améliorations, lecture de configuration, automatisation
!!      - Décembre 2016 à juillet 2017 : exercices du livre "Exercices et Problèmes d'aérodynamique fondamentale"
!!                                problème de l'interaction de chocs ajouté
!!      - Avril 2018  : 
!!                      - lecture des exercices à résoudre dans un fichier,
!!                      - ajouts de module, 
!!                      - prise en charge des commentaires avec doxygen
!!
!! NOTATIONS
!!      - pour les chocs :  amont 1, aval 2, indice n : normal
!!      - grandeurs isentropiques ou d'arrêt: incide i
!!
!! UTILISATION DU LOGICIEL dans un terminal : plusieurs possibililités
!!       1. EC
!!       2. EC opt= n                              # n : choix de l'option
!!       3. EC opt= n M= 3 teta= 10 A/Ac= 2.5     # n : choix de l'option
!!
!>               en fonction du choix (opt) M: Mach, teta : angle, ans : A/Acritique
!!
!>               M correspond aussi à 'cas' pour l'option 10 :
!>                   EC opt= 10 M= 1  (pour le cas 1)
!!
!>               Il faut absolument respecter l'espace après le signe =
!!       4. Lancement automatique par un script shell ./x1.run_EC.sh 
!!         
!!           Attention l'option est définie directement dans EC.in
!!       5. Il existe un mode où le programme lit la configuration dans le fichier cas.in
!!          et fait tous les calculs.
!!
!! @author  Christophe Airiau
!! \note
!!  - Le programme est lancé à partir du fichier main_EC.f90
!!  - On peut écrire un fichier fortran contenant juste l'instruction run_EC() avec l'utilisation
!!          du module mod_etude
!!
 
!
!==================================================================================
module mod_EC
!==================================================================================
!
use Constante
use type_def
use mod_read
use atmos

implicit none



integer, parameter                :: ndim=600                   ! nombre de valeurs pour omega(M)
integer, parameter                :: menu_size=40               ! taille du menu
real(kind=8)                      :: Mn1,M1,M2,Mn2              ! nombre de Mach
real(kind=8)                      :: om2,om1                    ! fonction omega
real (kind=8), dimension(ndim)    :: xm,ym
real(kind=8)                      :: sigmac,teta                ! angle de choc et déviation
real(kind=8)                      :: alpha=0.d0                 ! angle d'incidence
! real(kind=8), external            :: Angle_Choc,Mach_Aval       ! f : angle de choc, Mach aval
!real(kind=8), external            :: rho2_rho1,P2_P1, Pi2_Pi1   ! rapports à travers le choc
!real(kind=8), external            :: P_Pi,T_Ti,rho_rhoi         ! rapports isentropiques
! real(kind=8), external            :: omega,invomega             ! fonction omega
! real(kind=8), external            :: A_Acritic,inverse_A_Acritic! A/Acritique et sa réciproque
real(kind=8)                      :: ans                        ! answer value
integer                           :: i,j,ichoix
logical                           :: ask=.true.

integer                                     :: nz                     ! nombre de zone maximal
real(kind=8),dimension(:),allocatable       :: M,Ti,P,V               ! Mach, température d'arrêt, pression, vitesse
real(kind=8),dimension(:),allocatable       :: a,tau,T,Pi             ! Célérité du son, rapport de pression, température, pression d'arrêt
!real(kind=8),dimension(:),allocatable       :: tau_Ref                ! P_k/P_ref
real(kind=8),dimension(:),allocatable       :: sigma, theta,Om,Mn     ! angle de choc, déviation, Omega, Mach normal
character(len=30),dimension(:),allocatable  :: com                    ! commentaires
character(len=10),dimension(:),allocatable  :: type_zone              ! type de zone
character(len=1),dimension(:),allocatable   :: type_caracteristique   ! type de caractéristique
integer,dimension(:),allocatable            :: zone_amont             ! index de la zone amont
integer,dimension(:),allocatable            :: profil                 ! -1 : extrados, 1: intrados, 0 : pas profil
real(kind=8),dimension(:),allocatable       :: longueur               ! longueur de la zone en % de corde
real(kind=8),dimension(:),allocatable       :: CL_zone,CD_zone        !  CL et CD de la zone considérée
real(kind=8)                                :: Cl,Cd                  ! coefficient de portance et de traînée
real(kind=8)                                :: tmp1,tmp2,tmp3,tmp4,F
character(len=80),dimension(0:menu_size)    :: menu_title=''
type(state)                                 :: q
type(input),dimension(0:menu_size)          :: menu 
real(kind=8)                                :: chord                  ! corde d'un profil
real(kind=8)                                :: thickness              ! paramètre d'épaisseur
real(kind=8)                                :: Kp                     ! coefficient de pression
integer                                     :: N_exo                  ! nombre d'exercices à résoudre
character(len=4),allocatable,dimension(:)   :: liste_exo              ! liste des exercices
logical                                     :: show=.true.

contains





!********************
subroutine solve_case
!********************
!> @details
!! SUBROUTINE permettant de calculer l'état aérodynamique
!! en fonction du type de zone
! 
implicit none
integer         :: k
write(6,*)'solve case '
do k=0,nz
    write(6,100) k
    select case (trim(type_zone(k)))
    case ('Uniforme')
        call calcul_zone_uniforme(k)
    case ('Choc')
        call choc_oblique(zone_amont(k),k)
    case ('Choc droit')
        call choc_droit(zone_amont(k),k)
    case ('Detente')
        call detente_isentropique(zone_amont(k),k)
    case ('Glissement')
        call calcul_ligne_glissement(zone_amont(k),k)
    case ('Rien')
        call rien(zone_amont(k),k)
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
!! SUBROUTINE permettant de lire la configuration de l'écoulement
!! dans le fichier cas.in
! 
implicit none
integer                 :: i,k

open(1,form='formatted',file='cas.in')
read(1,*);read(1,*);read(1,*)
read(1,*) nz
read(1,*) alpha
allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz),zone_amont(nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
do i=0,nz
    read(1,*)
    read(1,*) k
    read(1,*) type_zone(k)
    write(6,*) "zone lue ",k,trim(type_zone(k))
    select case (trim(type_zone(k)))
    case ('Uniforme')
        read(1,*) M(k)
        read(1,*) Ti(k)
        read(1,*) P(k)
        read(1,*) theta(k)
        read(1,*) zone_amont(k)
        if (zone_amont(k).ne.k) theta(k)=theta(k)-alpha
    case ('Choc','Detente','Rien')
        read(1,*) theta(k)
        read(1,*) zone_amont(k)
        if (zone_amont(k).ne.k) theta(k)=theta(k)-alpha
    case ('Choc droit')
        read(1,*) zone_amont(k)
    case ('Glissement')
        read(1,*) theta(k)
        read(1,*) P(k)
        read(1,*) zone_amont(k)
    end select
end do
theta=theta*deg2rad
close(1)

end subroutine read_case

!************************
subroutine set_menu
!************************
!> @details
!! SUBROUTINE décrivant le menu du programme
!! et le numéro des options
! 
implicit none
integer     :: k
write(6,300)
 300  format(60('*'),/,10x,'UNIVERSITE PAUL SABATIER, TOULOUSE',/,   &
       10X,'MASTER  DE MECANIQUE ET ENERGETIQUE, AERODYNAMIQUE',/,/, &
       10X,'ENSEEIHT,  FORMATION APPRENTISSAGE, Eclt. Compr.',/,/, &
      2X, 'CALCUL DES ECOULEMENTS COMPRESSIBLES, CHOCS ET DETENTES, ...'  &
       ,/, 15X,'auteur : Christophe Airiau, Aout  2015',/,60('*'),/)
     
k=0 ; menu(k)%title='ONDES DE CHOC DROITES                ';menu(k)%n=00
k=1 ; menu(k)%title='ONDES DE CHOC OBLIQUES               ';menu(k)%n=01
k=2 ; menu(k)%title='DETENTE ou COMPRESSION ISENTROPIQUE  ';menu(k)%n=02
k=3 ; menu(k)%title='Valeur de OMEGA pour un MACH donne   ';menu(k)%n=03
k=4 ; menu(k)%title='Valeur du MACH pour un OMEGA donne   ';menu(k)%n=04
k=5 ; menu(k)%title='EVOLUTION ISENTROPIQUE               ';menu(k)%n=05
k=6 ; menu(k)%title='Interactions de chocs, plan (p,theta)';menu(k)%n=10
k=13; menu(k)%title='A / A critique pour un Mach          ';menu(k)%n=13
k=14; menu(k)%title='Mach pour A / A critique             ';menu(k)%n=14 
k=16; menu(k)%title='Probleme de Fanno F=f(Mach)          ';menu(k)%n=16
k=17; menu(k)%title='Probleme de Fanno inverse M=Fanno(F) ';menu(k)%n=17
k=19; menu(k)%title='Probleme de Rayleigh                 ';menu(k)%n=19
!
k=20; menu(k)%title='Plot Courbes de Rayleigh             ';menu(k)%n=20
k=21; menu(k)%title='Plot courbe de Fanno                 ';menu(k)%n=18

k=34; menu(k)%title='Fonction de Rayleigh  f(Mach)        ';menu(k)%n=21
k=35; menu(k)%title='Fonction de Rayleigh inverse         ';menu(k)%n=22

k=22; menu(k)%title='Plot EPICYCLOIDE                     ';menu(k)%n=08
k=23; menu(k)%title='Plot OMEGA fonction du MACH          ';menu(k)%n=06
k=24; menu(k)%title='Plot A / A critique = f(Mach)        ';menu(k)%n=15
k=25; menu(k)%title=' '                                    ;menu(k)%n=25
k=26; menu(k)%title=' '                                    ;menu(k)%n=26
k=27; menu(k)%title='Livre Aerodynamique Fondamentale     ';menu(k)%n=27
k=28; menu(k)%title=' '                                    ;menu(k)%n=28    ! TD 2 N7 
k=29; menu(k)%title=' '                                    ;menu(k)%n=29    ! TD 1 N7
k=30; menu(k)%title='Tables pour le livre                 ';menu(k)%n=12
k=31; menu(k)%title='TABLE DE CHOC                        ';menu(k)%n=09
k=32; menu(k)%title='TABLE DE LA FONCTION OMEGA           ';menu(k)%n=07
k=33; menu(k)%title='Tables personnalisees                ';menu(k)%n=23
k=36; menu(k)%title='Atmosphere standard                  ';menu(k)%n=24
k=39; menu(k)%title='EXAMEN Mai  2014                     ';menu(k)%n=90
k=38; menu(k)%title='EXAMEN Juin 2014                     ';menu(k)%n=91
k=menu_size; menu(k)%title='Cas dans un fichier'           ;menu(k)%n=99

end subroutine set_menu


!************************
subroutine display_menu
!************************
!> @details
!! SUBROUTINE permettant d'afficher le menu
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
!! SUBROUTINE permettant d'afficher une aide pour l'utilisation du logiciel
! 

write(6,100)

100 format('===>  Aide : ',/, 'Nom du programme : EC',/,     &
    'Faire au choix, en ligne de commande  et par exemple:',/, &
    'EC',/,                                                  &
    'EC opt= n                   avec n l''option',/,        &
    'EC opt= 1 M= 2.0 angle= 10.0  angle en degrés',/,        &
    'EC opt= 14 M= 0.5 A/Ac= 3.5 (pour une solution subsonique)',/,    &
    'EC opt= 14 M= 1.5 A/Ac= 3.5 (pour une solution supersonique)',/,  &
    'EC opt= 17 M= 0.5 F= 0.4  (pour une solution subsonique)',/,  &
    'EC opt= 17 M= 1.5 F= 0.4  (pour une solution supersonique)',/,  &
    'EC opt= 10 M= 1          (interaction choc, M est le cas )',/,  &
    'EC h                         pour l''aide',/,           &
    /,' Il faut respecter l''espace entre le nom et la valeur des options')
write(6,*)
write(6,'(a,/)') 'choix des tâches : '
call display_menu
end subroutine help

!*************************************
subroutine calcul_zone_uniforme(amont)
!*************************************
    !> @details
    !! SUBROUTINE permettant de calculer l'état aérodynamique d'une zone uniforme
    !! 
    !! entrées : Mach, P, Ti
    !!
    !! sorties : Pi, T, a, V, tau=P_aval/P_amont
    implicit none
    integer         :: amont
    write(6,*) 'indice amont         = ',amont
    write(6,*) 'Zone uniforme : Mach = ',M(amont)
    com(amont)='uniforme amont'
    if (Ti(amont).ne.0.d0) then
        print*,'calcul de T'
        T(amont)=T_Ti(M(amont))*Ti(amont)
    elseif (T(amont).ne.0.d0) then
        print*,'calcul de Ti'
        Ti(amont)=T(amont)/T_Ti(M(amont))
    else
        stop 'Probleme avec la zone uniforme et T'
    end if
    write(6,100) 'T/Ti = ',T_Ti(M(amont)),' T = ',T(amont)
    write(6,100) 'Ti = ',Ti(amont)
    a(amont)=sqrt(gam*r*T(amont));V(amont)=M(amont)*a(amont);
    write(6,100) 'a = ',a(amont),'V= ',V(amont)
    if (P(amont).ne.0.d0) then
        Pi(amont)=P(amont)/P_Pi(M(amont))
    elseif (Pi(amont).ne.0.d0) then
        P(amont)=Pi(amont)*P_Pi(M(amont))
    else
        stop 'Probleme avec la zone uniforme et P'
    end if

    write(6,100) 'P/Pi = ',P_Pi(M(amont)),' Pi = ',Pi(amont)
    write(6,100) 'P = ',P(amont)
    write(6,100) 'ac = ',critical_velocity(Ti(amont))
    write(6,100) 'V/ac = ',V(amont)/critical_velocity(Ti(amont))
    tau(amont)=1.d0
    100 format(2(a20,3x,f15.4,5x))
end subroutine calcul_zone_uniforme
 
!**************************
subroutine rien(amont,aval)
!**************************
    !> @details
    !! SUBROUTINE permettant de recopier une zone où il ne se passe rien
    implicit none
    integer         :: amont,aval
    write(6,*) 'Zone uniforme : Mach = ',M(amont) ;
    write(6,*) 'il ne se passe rien'
    write(6,*) ' comme dans la zone ',amont
    com(aval)='Comme la zone précédente'
    M(aval)=M(amont)
    P(aval)=P(amont)
    tau(aval)=tau(amont)
    Pi(aval)=Pi(amont)
    theta(aval)=theta(amont)
    T(aval)=T(amont)
    Ti(aval)=Ti(amont)
    a(aval)=a(amont);V(aval)=V(amont);
    write(6,110)'M',aval,M(aval)
    write(6,110)'Ti',aval,Ti(aval)
    write(6,110) 'T/Ti',aval,T_Ti(M(aval))
    write(6,110) 'P/Pi',aval,P_Pi(M(aval))
    write(6,110) 'T',aval,T(aval)
    write(6,110) 'a',aval,a(aval)
    write(6,110) 'V',aval,V(aval)
    write(6,110) 'ac = ',critical_velocity(Ti(aval))
    write(6,110) 'V/ac = ',V(aval)/critical_velocity(Ti(aval))

    
    110 format(a50,1x,'(',i3,')',':',3x,f12.4)
end subroutine rien

!*****************************************
subroutine detente_isentropique(amont,aval)
!*****************************************
!> @details
!! SUBROUTINE permettant de calculer une zone 
!! avec une compression/détente  isentropique
!!  en fonction du type de caractéristique
    implicit none
    integer         :: amont,aval
    real(kind=8)    :: deviation
    write(6,*) 'détente centrée isentropique : Mach = ',M(amont)
    deviation=theta(aval)-theta(amont)
    write(6,100) 'déviaton en degrés',deviation*rad2deg
    com(aval)='détente isentropique'
    Om(amont)=omega(M(amont))
    write(6,110)'OMEGA amont en degres',amont,om(amont)*rad2deg
    if (allocated(type_caracteristique)) then
        if (type_caracteristique(aval)=='+') then
            write(6,*) 'Onde montante'
            Om(aval)=Om(amont)+theta(aval)-theta(amont)
        else  
            write(6,*) 'Onde descendante'
            Om(aval)=Om(amont)+theta(amont)-theta(aval)
        end if
    else
        write(6,*) 'le type d onde n est pas precise avec type_caracteristique'
        write(6,*) 'on fixe Mach(aval) avec nécessairement une détente '
        Om(aval)=Om(amont)+abs(deviation)
    end if

    M(aval)=invomega(Om(aval))
    write(6,110)'OMEGA amont en degres',amont,om(amont)*rad2deg
    write(6,110)'OMEGA amont en degres',aval,om(aval)*rad2deg
    write(6,110)'MACH aval',aval,M(aval)
    tau(aval)=P_Pi(M(aval))/P_Pi(M(amont))
    write(6,110)'P/Pi',aval,P_Pi(M(aval))
    write(6,110) 'tau (aval/amont)',aval,tau(aval)
    P(aval)=P(amont)*tau(aval)
    Pi(aval)=Pi(amont)
    Ti(aval)=Ti(amont)
    T(aval)=T_Ti(M(aval))*Ti(aval)
    a(aval)=sqrt(gam*r*T(aval));V(aval)=M(aval)*a(aval);

    write(6,110) 'T/Ti',aval,T_Ti(M(aval))
    write(6,110) 'T',aval,T(aval)
    write(6,110) 'a',aval,a(aval)
    write(6,110) 'V',aval,V(aval)
    write(6,110) 'Pression P ',aval,P(aval)
    write(6,110) 'tau = P(k)/P(0)',aval,P(aval)/P(0)
    tau(aval)=P(aval)/P(0)

    write(6,120)'Mach aval ',M(aval),'theta(aval) en degres',theta(aval)*rad2deg
    write(6,100) 'ac = ',critical_velocity(Ti(aval))
    write(6,100) 'V/ac = ',V(aval)/critical_velocity(Ti(aval))

    100 format(a50,1x,':',3x,f12.4)
    110 format(a50,1x,'(',i3,')',':',3x,f12.4)
    120 format(2(a20,3x,f15.4,5x))

end subroutine detente_isentropique
!***********************************
subroutine choc_oblique(amont,aval)
!***********************************
!> @details
!! SUBROUTINE permettant de calculer un choc oblique
    implicit none
    integer         :: amont,aval
    real(kind=8)    :: deviation

    write(6,*) 'choc oblique entre la zone ',amont,' et la zone ', aval
    deviation=abs(theta(aval)-theta(amont))
    com(aval)= 'choc oblique'
    write(6,100) 'Mach amont',M(amont)
    write(6,100) 'déviaton en degrés',deviation*rad2deg
    call Newton_Angle_Choc(sigma(aval),deviation,M(amont))
    if (sigma(aval).gt.-1.0d0) then
        write(6,100)'ANGLE DU CHOC en degres',sigma(aval)*rad2deg
        Mn(amont)=M(amont)*sin(sigma(aval))
        write(6,110)'MACH NORMAL AMONT Mn',amont,Mn(amont)
        Mn(aval)=Mach_Aval(Mn(amont))
        M(aval) = Mn(aval)/sin(sigma(aval)-deviation)
        write(6,110)'MACH NORMAL AVAL  Mn',aval,Mn(aval)
        write(6,110)'MACH AVAL  M'        ,aval,M(aval)
        write(6,110)'RAPPORT T/Ti'       ,aval, T_Ti(M(aval))
        write(6,120)'RAPPORT DES MASSES VOLUMIQUES rho ',aval,amont,rho2_rho1(Mn(amont))
        write(6,120)'RAPPORT DES PRESSIONS P',aval,amont,P2_P1(Mn(amont))
        write(6,120)'RAPPORT DES PRESSIONS TOTALES Pi',aval,amont,Pi2_Pi1(Mn(amont))

        tau(aval)=P2_P1(Mn(amont))*P(amont)/P(0)
        Pi(aval)=Pi(amont)*Pi2_Pi1(Mn(amont))
        P(aval)=P(amont)*P2_P1(Mn(amont))
        write(6,110)'P/Pi',aval,P_Pi(M(aval))
        write(6,110)'tau= P(k)/P(0) =  ',aval,tau(aval)
        write(6,110)'Pi',aval,Pi(aval)
        write(6,110)'P',aval,P(aval)

        Ti(aval)=Ti(amont)
        write(6,110)'Ti',aval,Ti(aval)
        T(aval)=T_Ti(M(aval))*Ti(aval)
        a(aval)=sqrt(gam*r*T(aval));V(aval)=M(aval)*a(aval);
        write(6,110) 'T/Ti',aval,T_Ti(M(aval))
        write(6,110) 'T',aval,T(aval)
        write(6,110) 'a',aval,a(aval)
        write(6,110) 'V',aval,V(aval)
        write(6,100) 'ac = ',critical_velocity(Ti(aval))
        write(6,100) 'V/ac = ',V(aval)/critical_velocity(Ti(aval))
    else
        write(6,*)' le choc est détaché '
    end if

    100 format(a57,1x,':',3x,f12.4)
    110 format(a50,1x,'(',i3,')',':',3x,f12.4)
    120 format(a50,1x,'(',i3,' / ',i3,')',':',3x,f12.4)

end subroutine choc_oblique


!***********************************
subroutine choc_droit(amont,aval)
!***********************************
!> @details
!! SUBROUTINE permettant de calculer un choc droit
    implicit none
    integer         :: amont,aval
    sigma(aval)=val_pi/2.d0
    write(6,*) 'choc droit entre la zone ',amont,' et la zone ', aval
    com(aval)= 'choc droit'
    write(6,100) 'Mach amont',M(amont)
    M(aval)=Mach_Aval(M(amont))
    write(6,110)'MACH AVAL  M'        ,aval,M(aval)
    write(6,110)'RAPPORT T/Ti'        ,aval, T_Ti(M(aval))
    write(6,120)'RAPPORT DES MASSES VOLUMIQUES rho ',aval,amont,rho2_rho1(M(amont))
    write(6,120)'RAPPORT DES PRESSIONS P',aval,amont,P2_P1(M(amont))
    write(6,120)'RAPPORT DES PRESSIONS TOTALES Pi',aval,amont,Pi2_Pi1(M(amont))

    tau(aval)=P2_P1(M(amont))*P(amont)/P(0)
    Pi(aval)=Pi(amont)*Pi2_Pi1(M(amont))
    P(aval)=P(amont)*P2_P1(M(amont))
    write(6,110)'P/Pi',aval,P_Pi(M(aval))
    write(6,110)'tau= P(k)/P(0) =  ',aval,tau(aval)
    write(6,110)'Pi',aval,Pi(aval)
    write(6,110)'P',aval,P(aval)

    Ti(aval)=Ti(amont)
    write(6,110)'Ti',aval,Ti(aval)
    T(aval)=T_Ti(M(aval))*Ti(aval)
    a(aval)=sqrt(gam*r*T(aval));V(aval)=M(aval)*a(aval);
    write(6,110) 'T/Ti',aval,T_Ti(M(aval))
    write(6,110) 'T',aval,T(aval)
    write(6,110) 'a',aval,a(aval)
    write(6,110) 'V',aval,V(aval)
    write(6,100) 'ac = ',critical_velocity(Ti(aval))
    write(6,100) 'V/ac = ',V(aval)/critical_velocity(Ti(aval))
    
    100 format(a57,1x,':',3x,f12.4)
    110 format(a50,1x,'(',i3,')',':',3x,f12.4)
    120 format(a50,1x,'(',i3,' / ',i3,')',':',3x,f12.4)

end subroutine choc_droit


!*********************************************
subroutine calcul_ligne_glissement(amont,aval)
!*********************************************
!> @details
!! SUBROUTINE permettant de  calculer une zone délimitée par une ligne de glissement
    implicit none
    integer         :: amont,aval
    real(kind=8)    :: theta_montante,theta_descendante,deviation
    write(6,*) 'calcul d une ligne isobare (ligne de glissement)'
    com(aval)='Ligne de glissement isobare'
    tau(aval)=P(aval)/P(amont)
    M(aval)=dsqrt(2.d0/(gam-1.d0)*( 1.d0/T_Ti(M(amont))*tau(aval)**(1.d0/gam-1.d0)-1.d0))
    om(amont)=omega(M(amont));om(aval)=omega(M(aval))
    deviation=abs(om(aval)-om(amont))
    Ti(aval)=Ti(amont)
    Pi(aval)=P(aval)/P_Pi(M(aval))
    T(aval)=T_Ti(M(aval))*Ti(aval)
    a(aval)=sqrt(gam*r*T(aval));V(aval)=M(aval)*a(aval);

    write(6,110)'MACH amont',amont,M(amont)
    write(6,110)'OMEGA amont en degres',amont,om(amont)*rad2deg
    write(6,110)'OMEGA amont en degres',aval,om(aval)*rad2deg
    write(6,110)'MACH aval',aval,M(aval)
    write(6,100)'deviaton en degres',deviation*rad2deg
    theta_montante=theta(amont)+om(aval)-om(amont)
    write(6,100)'theta_g, onde montante',theta_montante*rad2deg
    theta_descendante=theta(amont)-om(aval)+om(amont)
    write(6,100)'theta_g, onde descendante',theta_descendante*rad2deg 
    write(6,110)'P(k)/P(amont) =  ',aval,tau(aval)
    write(6,110)'Pi',aval,Pi(aval)
    write(6,110)'P',aval,P(aval)
    write(6,110)'P/Pi',aval,P_Pi(M(aval))
    write(6,110) 'T/Ti',aval,T_Ti(M(aval))
    write(6,110)'Ti',aval,Ti(aval)
    write(6,110) 'T',aval,T(aval)
    write(6,110) 'a',aval,a(aval)
    write(6,110) 'V',aval,V(aval)
    write(6,110)'tau= P(k)/P(0) =  ',aval,tau(aval)
    write(6,100) 'ac = ',critical_velocity(Ti(aval))
    write(6,100) 'V/ac = ',V(aval)/critical_velocity(Ti(aval))

    if (allocated(type_caracteristique)) then
        if (type_caracteristique(aval)=='+') then
            write(6,*) 'Onde montante'
            theta(aval)=theta_montante
        else  
            write(6,*) 'Onde descendante'
            theta(aval)=theta_descendante
        end if
    else
        write(6,*) 'le type d onde n est pas precise avec type_caracteristique'
        write(6,*) 'on fixe theta(aval) avec l onde montante arbitrairement '
        theta(aval)=theta_montante
    end if
    write(6,120)'Mach aval ',M(aval),'theta(aval) en degres',theta(aval)*rad2deg
    tau(aval)=P(aval)/P(0)
    120 format(2(a20,3x,f15.4,5x))

    100 format(a57,1x,':',3x,f12.4)
    110 format(a50,1x,'(',i3,')',':',3x,f12.4)
end subroutine calcul_ligne_glissement


!***********************
subroutine calcul_Cl_Cd
!***********************
!> @details
!! SUBROUTINE permettant de  calculer le Cl et le Cd en fonction des rapports de pression
    implicit none
    real(kind=8)    :: const
    write(6,*) 'Calcul des coefficients aérodynamiques'
    const=-2.d0/(gam*M(0)**2)
    ! profil triangle isocèle en incidence
    Cl=const*(0.5d0/cos(theta(0))*(tau(1)+cos(theta(2))*tau(2))-cos(theta(3))*tau(3))
    Cd=const*(0.5d0/cos(theta(0))*sin(theta(2))*tau(2)-sin(theta(3))*tau(3))
    write(6,100) 'Cl ',Cl
    write(6,100) 'Cd ',Cd

    100 format(a50,1x,':',3x,f12.4)

end subroutine calcul_Cl_Cd

!*************************************
subroutine coefficients_aerodynamiques
!*************************************
!> @details
!! SUBROUTINE permettant de  calculer les coefficients aérodynamiques
   implicit none
   real(kind=8)    :: coeff,tmp,beta
   integer         :: k
   coeff=2.d0/(gam*M(0)**2)
   write(6,*) 'coefficients aérodynamiques'
   CL_zone=0.d0;CD_zone=0.d0
   write(6,100)
   do k=0,nz
     select case (profil(k))
    case(1)
          !write(6,*)'zone extrados :',k,theta(k)*rad2deg
          tmp=coeff*(tau(k)-1.d0)
          CL_zone(k)=tmp
          CD_zone(k)=tmp*tan(theta(k))
     case(-1)
          !write(6,*)'zone intrados :',k,theta(k)*rad2deg
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
         !write(6,*)'Kp, theta en degrés, tau              = ',CL_zone(k),theta(k)*rad2deg,tau(k)
         !write(6,*)'Kp tan(theta) (nonlineaire, linéaire)  = ',CD_zone(k),CL_zone(k)*theta(k)
      end if
   end do
  ! formule exacte pour l'intégration des forces de pression sur les deux directions
  write(6,120)'CL Total par projection du Kp des panneaux    = ',-sum(CL_zone*longueur*profil)
  write(6,120)'CD Total par projection du Kp des panneaux    = ', sum(CD_zone*longueur*profil)
  beta=sqrt(M(0)**2-1)
  ! théorie linéaire d'Ackeret 
  ! effet d'incidence seulement : 
  write(6,120)'CL theorie linéaire : 4 alpha/beta            = ',4.d0*alpha*deg2rad/beta
  write(6,120)'CD theorie linéaire : 4 alpha^2/beta          = ',4.d0*(alpha*deg2rad)**2/beta, ' pour la plaque plane (pas de CD0)' 
  ! effet d'épaisseur et d'incidence (dans theta on a pris en compte l'incidence
  write(6,120)'CD theorie linéaire : int theta^2 dx (pente)  = ',2.d0/beta*sum(longueur*tan(theta)**2), 'pour un corps quelconque'
  write(6,120)'CD theorie linéaire : int theta^2 dx (angle)  = ',2.d0/beta*sum(longueur*theta**2), 'pour un corps quelconque'
  !            le tan theta dans la formule au dessus est pour respecté l'esprit de la théorie linéarisée où
  ! on assimile l'angle à la pente ... Sinon y'a des petites différences
   
   100 format('#', 1x,'zone', 4x,'Kp',9x,'theta (°)',8x,'P/P0',5x,'Kp tan(theta)',5x,'Kp theta' &
            3x,'longueur',3x,'extrados/intrados')
   110 format(i3,6(f12.5,2x),i2)
   120 format('# ',a50,2x,f12.5,2x,a)
end subroutine coefficients_aerodynamiques


!********************
subroutine sorties
!********************
!> @details
!! SUBROUTINE permettant d'afficher les sorties à l'écran
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
        'sigma',10x,'Omega',14x,'Commentaire',25x,'T',13x,'V',13x,'a')
 100 format(i4,2x,7(f12.4,3x),a40,3(f12.4,3x))
 110 format(i3, 15x, ' ==> Choc détaché, pas de résultats')

end subroutine sorties


!
! I - Problème avec chocs
!
!************************************************
 subroutine Newton_Angle_Choc(sigma,teta,M0)
!************************************************
!> @details
!! méthode de Newton pour déterminer l'angle de choc oblique
      use Constante
      implicit none
      real(kind=8)              :: sigma,teta,M0
      real (kind=8)             :: er0,ds,f0,f1,dsigma,erreur
!      real (kind=8),external    :: Angle_Choc
      integer i
      if (sigma.eq.0.d0) then
          write(6,*) 'initialisation classique'
          sigma=asin(1.d0/M0)
      end if
!     if (fichier) then
!         init=0.d0
!     else
!         print*,'entrer initialisation : 0 ou 1'
!         read*,init
!     end if
      er0=1d-6
      ds=1d-2*deg2rad
      erreur=1.d0
      i=1
      if (show) write(6,101)
  101 format(50('-'),/, 'iter',4x,'Sigma',8x,'dsigma',10x,'erreur')
  100 format(i3,f12.5,2x,2(f14.9,2x))             
      do while ( (abs(erreur).gt.er0).and.(i.le.20))
          f0=Angle_Choc(sigma,teta,M0)
          f1=Angle_Choc(sigma+ds,teta,M0)
          dsigma=- ds * f0/(f1-f0)
    !      write(6,*)*,f0,f1,i
          erreur=dsigma/sigma
          if (show) write(6,100)i,sigma*rad2deg,dsigma*rad2deg,erreur
          sigma=sigma+dsigma
          i=i+1
      end do

      if (i.gt.20) then
        write(6,*)'pas de convergence'
        write(6,*)'  angle de déviation  est supérieur'
        write(6,*)'   à celui détachant le choc' 
        sigma=-1.d0
      else if (sigma.lt.0) then
        write(6,*) 'Non convergence, angle négatif :'
        write(6,*) 'Il faut augmenter le Mach amont'
        write(6,*) "ou diminuer angle de déviation"
        stop "erreur sur angle de choc"
        sigma=-1.d0
      else if (sigma*rad2deg.gt.90.d0) then
        write(6,*) 'Non convergence, angle > 90 :'
        write(6,*) 'Il faut augmenter le Mach amont'
        write(6,*) "ou diminuer angle de déviation"
        sigma=-1.d0
      else
          if (show) write(6,*) 'convergence sur sigma'
      endif             
       if (sigma < 0) stop 'pas de convergence sur sigma'
      if (show) write(6,102)
  102 format(50('-'))  
      return
      end subroutine Newton_Angle_Choc
      
 
!************************************************       
      function Angle_Choc(sigma,teta,M0) result(res)
!************************************************       
!> @details
!!     fonction donnant l'angle du choc : f
      use Constante 
      implicit none 
      real (kind=8):: res, sigma,teta,M0
      res= tan (sigma-teta)/tan(sigma) -2.d0/(gam+1)/M0**2/(sin(sigma))**2 - (gam-1.d0)/(gam+1.d0)  
      end function Angle_Choc
!************************************************       
      function dfsigma(sigma,teta,M0)
!************************************************       
!> @details
!! dérivée de la fonction f / sigma
      use Constante
      implicit none
!      gradient de la function pr�c�dente      
      real (kind=8)::dfsigma, sigma,teta,M0,f1,f2

      f1= tan(sigma)/(cos(sigma-teta))**2 - tan (sigma-teta)/(cos(sigma))**2
      f2=f1/(tan(sigma))**2  +4.d0/(gam+1.d0)/M0**2    * cos(sigma)/(sin(sigma))**3
      dfsigma=f1+f2
      end function dfsigma
!************************************************       
      function Pi2_Pi1(mach)
!************************************************       
!> @details
!!     saut de pression d'arret à travers un choc
      use Constante
      implicit none
      real(kind=8)::  Pi2_Pi1,mach,x,y
      !real(kind=8)::  rho2_rho1,P2_P1
      x=-1.d0/(gam-1.d0)
      y=-gam*x
      Pi2_Pi1= P2_P1(mach)**x * rho2_rho1(mach)**y
      end function Pi2_Pi1
          
!************************************************       
      function P2_P1(mach)
!************************************************       
!> @details
!!     saut de pression à travers un choc      
      use Constante
      implicit none
      real (kind=8)::  P2_P1,mach
      P2_P1= 2.d0*gam/(gam+1.d0)* mach**2- (gam-1.d0)/(gam+1.d0)
      end function P2_P1
!************************************************       
      function Inverse_P2_P1(rapport) result(Mach)
!************************************************       
!> @details
!!     Fonction inverse, connaissant le saut de pression à travers un choc
!!       on retrouve le Mach normal amont      
      use Constante
      implicit none
      real (kind=8)::  rapport,Mach
      Mach= sqrt(1.d0/(2.d0*gam)* ((gam+1.d0)* rapport+ (gam-1.d0)))
      end function Inverse_P2_P1      
!************************************************       
      function rho2_rho1(mach)
!************************************************       
!> @details
!!     saut de masse volumique à travers un choc     
      use Constante
      implicit none
      real (kind=8):: rho2_rho1,mach
      rho2_rho1= 1.d0/( 2.d0/((gam+1.d0)* mach**2)+ (gam-1.d0)/(gam+1.d0) )
      end function rho2_rho1
      
!************************************************       
      function Mach_Aval(mach)
!************************************************       
!> @details
!!     Mach normal aval à travers un choc      
      use Constante
      implicit none
      real (kind=8):: mach,Mach_Aval
      Mach_Aval= sqrt((1.d0+ 0.5d0*(gam-1.d0)* mach**2)/(gam*mach**2-0.5d0*(gam-1.d0)))
      end function Mach_Aval

!************************************************       
      function Saut_Entropie(Pi_Amont,Pi_Aval)
!************************************************
!> @details
!! calcul du saut d'entropie à travers un choc       
      use Constante
      implicit none
      real (kind=8):: Pi_Amont,Pi_Aval,Saut_Entropie
      Saut_Entropie=-r*log(Pi_Aval/Pi_Amont)

      end function Saut_Entropie
!
! II- Problème isentropique
!
!************************************************    
      function inverse_P_Pi(ratio) result(Mach)
!************************************************    
!> @details
!!     rapport de pression statique sur la pression d'arrêt
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
!!     rapport de pression statique sur la pression d'arrêt
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
!!     rapport de masse volumique sur la masse volumique à l'arrêt
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
!!     rapport de température sur la température d'arrêt
      use Constante
      implicit none
      real(kind=8)::  T_Ti,mach
      T_Ti= 1.d0/ (1.d0+ 0.5d0*(gam-1.d0)* mach**2)
      end function T_Ti
!************************************************       
      function critical_velocity(Ti)
!************************************************       
!> @details
!!     rapport de température sur la température d'arrêt
      use Constante
      implicit none
      real(kind=8)::  critical_velocity,Ti
      critical_velocity=sqrt( 2*gam*r/(gam+1)*Ti)
      end function critical_velocity
!************************************************
      function omega(M)
!************************************************
!> @details
!!    fonction omega en fonction du Mach, en radians
!!   il y a deux formules pour le calcul de omega (cf commentaires )
      use Constante
      implicit none
      real (kind=8):: M,omega,c,beta
!      real (kind=8):: demi_pi,mu,omega0
      if (M.lt.1.d0) then
        write(6,*)'erreur dans omega Mach < 1   => M = ',M
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
!!    inverse de la fonction omega : M en fonction de omega
      use Constante
      implicit none
      !real (kind=8):: omega
      real (kind=8)::invomega,angle
      real (kind=8):: omega0,omega1,er0, erreur,mach0,dm,dmach
       
      integer i
      er0=1d-6
      if (show) write(6,*)'valeur cherchee en deg                     :  ', angle*rad2deg
      
      dmach=1d-3
      mach0=2.d0
      i=1
      erreur=1.d0

      do while ((erreur.gt.er0).and.(i.le.20))
        omega0=omega(mach0)-angle
        omega1=omega(mach0+dmach)-angle
        dm=- dmach * omega0/(omega1-omega0)
  !      write(6,*)*,omega0,omega1,i
        erreur=abs(dm/mach0)
        if (show) write(6,*)'i = ',i,'erreur = ',erreur,'dm =',dm,'dM= ',mach0
        mach0=mach0+dm
        i=i+1
      end do

      if (i.gt.20) then
        write(6,*)'pas de convergence'
        stop
      else
          if (show) write(6,*)'convergence sur le Mach, iter              :  ', i
        invomega=mach0
      endif             
      end function invomega

 !**********************************************
 subroutine courbe_omega
 !**********************************************
 !> @details
 !! fonction omega (Mach)
      use Constante
      implicit none
      real(kind=8):: M1,M2,DM
!     real(kind=8),external:: omega
      integer, parameter ::npt=1000
      integer i
      open(1,form='formatted',file='fonction_omega.dat')
      M1=1.d0;  M2=15.d0
      DM=(M2-M1)/float(npt-1)
      write(1,100) 
      100 format('#',10x,'fonction omega (Mach), en degres')
      do i=1,npt
        write(1,'(2(e12.5,3X))')M1,omega(M1)*rad2deg
        M1=M1+DM
      end do   
      close(1)
 end subroutine courbe_omega     

!*******************
function S_sur_Sc(M)
!*******************
!> @details
!! rapport de la section à la section critique en fonction du Mach
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
!! trace le l'epicycloide pour le probleme de détente
!
!
use Constante
implicit none

! epicylcoide en coordonnees critiques

integer, parameter          :: npt=10000
real(kind=8), parameter     :: M_inf=1000.d0
real(kind=8)                :: mach,dmach,mach_star,theta,theta_max,M_max
integer                     :: i

M_max=sqrt((gam+1.d0)/(gam-1.d0))

 open(1,form='formatted',file='epi.dat')
 write(1,90)
 mach=1.d0
 dmach=(M_inf-mach)/float(npt-1)
 do i=1,npt
   call calcul(mach,theta,mach_star)
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
subroutine calcul(x,omega,mach)
!**********************************
!> @details
!! calcul du Mach M* en fonction de M (ou inversement à vérifier)
use Constante
implicit none
real(kind=8):: x,omega, mach
real(kind=8):: c1,c2
 
 c2=(gam-1.d0)/(gam+1.d0)
 c1=sqrt(1.d0/c2)
 omega = c1*atan(sqrt(c2*(x*x-1.d0)))-atan(sqrt(x*x-1.d0))
! write(6,*) x,atan(sqrt(c2*(x*x-1))), atan(sqrt(x*x-1))
 mach=sqrt((gam+1.d0)*x*x/(2.d0+(gam-1.d0)*x*x))
end subroutine calcul


!
! III-  les tables
!

!**********************************
subroutine table_choc
!**********************************
!> @details
!! table de chocs dans un fichier de sortie
!
implicit none
integer                     ::n,i,j,iloc,step
real(kind=8),allocatable,dimension(:) ::m1,m2,rp,rpi,rti
real(kind=8)                ::dm,m1i,m1f
! real (kind=8), external     ::Mach_Aval,P2_P1,P_Pi,T_Ti
integer                     ::k=3

m1i=1.00d0          ! mach initial
m1f=2.00d0          ! mach final
dm=0.010d0
n=int((m1f-m1i)/dm)+2
n=k*(int(n/k)+1)
allocate(m1(n),m2(n),rp(n),rpi(n),rti(n))

write(6,100) m1i,m1f,dm,n
100 format('Table de choc',/,'Mach initial : ',f5.2,3x,'Mach final : ',f5.2,3x, &
    'pas en Mach : ',f5.3,3x,'Nombre de points :',i5)
m1(1)=m1i
do i=1,n-1
    m2(i)=Mach_Aval(m1(i)); rp(i)=P2_P1(m1(i));rpi(i)=P_Pi(m1(i))
    rti(i)=T_Ti(m1(i))
    m1(i+1)=m1(i)+dm
end do

open(1,form='formatted',file='table_choc.out')
open(2,form='formatted',file='table_choc.tex')
open(3,form='formatted',file='table_choc_vert.tex')
write(1,110)
write(2,111);write(3,111)
111 format('\tiny{\begin{tabular}{',3(5('|c'),'|')'}',&
        '\hline')
!write(2,112);
!write(3,112)

write(3,150)'&';write(3,150)'&';write(3,150)'\\ \hline'
write(2,150)'&';write(2,150)'&';write(2,150)'\\ \hline'

110 format(3(3x,'M1',4x,'M2',4x,'P2/P1 '))
!112 format(3('$Mn_1$&','$Mn_2$&','$P_2/P_1$&$P_1/P_{i1}$&$T_1/T_{i1}$&'),'\\ \hline')
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
    !iloc=(i-1)*k+1;write(6,*)*,'iloc = ',iloc
    iloc=i
    write(3,122)(m1(iloc+j*step),m2(iloc+j*step), rp(iloc+j*step),rpi(iloc+j*step),rti(iloc+j*step),j=0,k-1)
end do
write(2,132);write(3,132)
132 format('\hline \end{tabular}}')
close(1);close(2);close(3)
end subroutine table_choc

!**********************
subroutine tables_livre
!**********************
!> @details
!! Table telles que sorties dans le livre
use Constante
implicit none
real(kind=8)            :: dM,M_init,M_final,Mach
integer                 :: nptM     ! nombre de points
integer,dimension(4)    :: npt      ! nombre de points
integer                 :: nc       ! nombre de colonnes
integer                 :: nl       ! nombre de lignes
integer                 :: reste,i,nl1,k,j,i1,i2,i_tmp
real(kind=8),allocatable,dimension(:,:) :: table(:,:)
!real(kind=8), external            :: T_Ti,P_Pi,rho_rhoi,S_sur_Sc
!real(kind=8),external            :: P2_P1,rho2_rho1,Pi2_Pi1,Mach_Aval,omega
real(kind=8),dimension(4)   :: dM_tmp,M_init_tmp,M_final_tmp
!character(len=2)        :: charac
logical                 :: test
!=====================================
write(6,*) 'table pour le subsonique'
!=====================================

nc = 2; k=6
write(6,*) 'nombre de colonnes choisi                 :  ', nc
M_init=0.0d0; M_final=1.0d0; dM=0.01d0
nptM=int((M_final-M_init)/dM)+1
write(6,*) 'nombre de points                          : ', nptM
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
    write(6,*) 'une des colonnes ne sera pas pleine'
    write(6,*) 'il y aura ',reste,' valeurs uniquement'
    nl1=nl+1
else
    nl1=nl
end if

!
!   2 colonnes : 1 fichier
!
write(6,*)' nombre de lignes : ',nl,nptM,nptM/nc
open(10,form='formatted',file='subsonique.tex')
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
    write(6,*)' dernière ligne'
    write(10,101) table(nl1,:)
end if
write(10,120)
close(10)

!
! 1 colonne, 2 fichiers
!
open(10,form='formatted',file='subsonique1.tex')
write(10,111); write(10,130)'\\ \hline \hline'
do i=1,nptM/2
    write(10,102) table(i,:)
    if (mod(i,5).eq.0) write(10,*) '\hline'
end do
write(10,120)
close(10)

open(10,form='formatted',file='subsonique2.tex')
write(10,111); write(10,130)'\\ \hline \hline'
do i=nptM/2+1,nptM
    write(10,102) table(i,:)
    if (mod(i,5).eq.0) write(10,*) '\hline'
end do
write(10,120)
close(10)



deallocate(table)


!===================================================
write(6,*) 'table pour le supersonique isentropique'
!===================================================
k=7
M_init=1.00d0; M_final=5.0d0; dM=0.05d0
nptM=int((M_final-M_init)/dM)+1
write(6,*) 'nombre de points                          : ', nptM
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
    write(6,*) 'une des colonnes ne sera pas pleine'
    write(6,*) 'il y aura ',reste,' valeurs uniquement'
    nl1=nl+1
else
    nl1=nl
end if
write(6,*)' nombre de lignes                          : ',nl,nptM,nptM/nc
!
!   2 colonnes : 1 fichier
!
open(10,form='formatted',file='supersonique.tex')
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
    write(6,*)' dernière ligne'
    write(10,201) table(nl1,:)
end if
write(10,220)
close(10)
deallocate(table)


!
! 1 colonne, plusieurs fichiers
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
write(6,*)'mach initial                               : ',table(1,1)
do j=1,3
    i2=i1+npt(j)
    write(30,*)' Mach init et final et pas',M_init_tmp(j), M_final_tmp(j),dM_tmp(j)
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
    open(10,form='formatted',file='supersonique'//charac(2,j)//'.tex')
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
write(6,*) 'table pour les chocs droits ou obliques '
!==================================================
k=6



M_init=1.00d0; M_final=5.0d0; dM=0.05d0
nptM=int((M_final-M_init)/dM)+1
write(6,*) 'nombre de points : ', nptM
allocate(table(nptM,k))
do i=1,nptM
    Mach=M_init+dM*real(i-1,kind=8)
    table(i,1)=Mach                     ! mach normal amont
    table(i,2)=Mach_Aval(Mach)              ! mach normal aval
    table(i,3)=P2_P1(Mach)             ! rapport des pressions statiques
    table(i,4)=T_Ti(table(i,2))/T_Ti(Mach)  ! rapport des températures
    table(i,5)=rho2_rho1(Mach)           ! rapport des masses volumiques
    table(i,6)=Pi2_Pi1(Mach)            ! rapport des pressions isentropiques
end do
nl=nptM/nc
reste=mod(nptM,nc)
if (reste.ne.0) then
    write(6,*) 'une des colonnes ne sera pas pleine'
    write(6,*) 'il y aura ',reste,' valeurs uniquement'
    nl1=nl+1
else
    nl1=nl
end if
write(6,*)' nombre de lignes : ',nl,nptM,nptM/nc
open(10,form='formatted',file='chocs_droits.tex')
write(10,310)
! uniquement pour deux colonnes
write(10,330)'&';write(10,330)'\\ \hline \hline'
write(6,*) 'nl1                                       : ',nl1
do i=1,nl1-1
    write(10,300) table(i,:),table(i+nl1,:)
    if (mod(i,5).eq.0) write(10,*) '\hline'
end do
write(6,*) 'mod                                       : ',mod(reste,2)
if (mod(reste,2).eq.1) then
    write(6,*)' dernière ligne'
    write(10,301) table(nl1,:)
end if
write(10,320)
close(10)
deallocate(table)

!
!  1 colonnes, plusieurs fichiers
!
k=6
nptM=sum(npt)+1; 
allocate(table(nptM,k))
write(6,*) 'CHOCS: npt ',npt(:), ' total = ', sum(npt)
i1=1
table(i1,1)=M_init_tmp(1)
write(6,*)'mach initial : ',table(1,1)
do j=1,3
    i2=i1+npt(j)
    do i=i1+1,i2; table(i,1)=table(i-1,1)+dM_tmp(j); end do
    i1=i2
end do

do i=1,nptM
    Mach=table(i,1)
    table(i,2)=Mach_Aval(Mach)              ! mach normal aval
    table(i,3)=P2_P1(Mach)             ! rapport des pressions statiques
    table(i,4)=T_Ti(table(i,2))/T_Ti(Mach)  ! rapport des températures
    table(i,5)=rho2_rho1(Mach)           ! rapport des masses volumiques
    table(i,6)=Pi2_Pi1(Mach)            ! rapport des pressions isentropiques
end do
do j=1,4
    open(10,form='formatted',file='chocs'//charac(2,j)//'.tex')
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
! Table pour Prandtl Meyer
!
write(6,*)' TABLE DE PRANDTL-MEYER'
dM_tmp=(/0.01d0, 0.02d0, 0.5d0, 0.0d0/)
M_init_tmp =(/1.d0, 3.5d0,  5.d0, 0.d0/)
M_final_tmp=(/3.5d0, 5.d0, 30.d0, 0.d0/)
npt(4)=0
do j=1,3
    write(6,*)'j, ',j,' step ', ((M_final_tmp(j)-M_init_tmp(j))/dM_tmp(j)), M_final_tmp(j),M_init_tmp(j),dM_tmp(j)
    npt(j)=int((M_final_tmp(j)-M_init_tmp(j))/dM_tmp(j))
end do

nptM=sum(npt)+1; 
write(6,*) 'nombre de points                          : ', nptM


k=3
allocate(table(nptM,k))
write(6,*) 'npt ',npt(:), ' total = ', sum(npt)
i1=1
table(i1,1)=M_init_tmp(1)
write(6,*) 'mach initial                              : ',table(1,1)
do j=1,3
    i2=i1+npt(j)
    write(31,*)' Mach init et final et pas',M_init_tmp(j), M_final_tmp(j),dM_tmp(j)
    do i=i1+1,i2
        table(i,1)=table(i-1,1)+dM_tmp(j)
        write(31,'(i4,3x,f12.3,5x,i3)')i,table(i,1),i-i1
    end do
    i1=i2
end do

do i=1,nptM
    Mach=table(i,1)
    table(i,2)=omega(Mach)*rad2deg             ! mach normal aval
    table(i,3)=asin(1.d0/Mach)*rad2deg
end do

do j=1,3      ! nombres de pages, 3 colonnes
    open(10,form='formatted',file='omega'//charac(2,j)//'.tex')
    write(10,411); 
    write(10,430) ' & '; write(10,430) ' & '; write(10,430)'\\ \hline \hline'

    do i=1,40           ! chaque ligne du tableau
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

! chocs droits ou obliques
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

! table de mach, omega, mu

! supersonique

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
! IV- Section critique, Tuyere de Laval
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
    !! first derivative wrt Mach
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
    !! fonction réciproque de A/A_crit=F(Mach)
    implicit none

    real(kind=8),intent(in)     :: Ma_ref,C_ref,tol_max
    integer,intent(in)          :: opt_newton,opt_function,iter_max
    character(len=3)            :: opt
    real(kind=8)                :: res,M_init,M_end,const
!     real(kind=8),external       ::newton,newton_new 
    const=C_ref
    if (const.lt.1.d0) then
         stop 'A /A critique < 1 , pas de solution'
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
!! Méthode de Newton avec la dérivée numérique 
!! pour la fonction réciproque de A/A_crit=F(Mach)
IMPLICIT NONE
INTEGER,INTENT(IN)      :: opt
REAL(KIND=8),INTENT(IN) :: const  ! valeur recherchée
REAL(KIND=8),INTENT(IN) :: x_init   ! proche du zero cherché
REAL(KIND=8),INTENT(IN) :: tol_max  ! tolerance maximale
INTEGER,INTENT(IN)      :: iter_max ! nombre d'itérations maximales
REAL(KIND=8)            :: res

REAL(KIND=8)        :: y0,dx2,x0,dx
REAL(KIND=8)        :: x1,x2,y1,y2,errx,erry
INTEGER             :: iter
REAL(KIND=8),PARAMETER :: epsil=1d-9  !pourcentage du pas
REAL(KIND=8)        :: tol
! REAL(KIND=8),EXTERNAL:: f_analytique
x0=x_init
tol=1.d0
iter=0

     
    ! calcul de la dérivée numérique

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
          write(6,*) 'pas de convergence : stop'
          STOP
END IF

!write(6,*) 'Solution : x=  ',x0, 'f(x) = ', y0
res=x0
END FUNCTION newton

!****************************************************
FUNCTION newton_new(opt,const,x_init,iter_max,tol_max) result(res)
!****************************************************
!> @details
!! Méthode de Newton avec la dérivée numérique 
!! pour la fonction réciproque de A/A_crit=F(Mach)
IMPLICIT NONE
INTEGER,INTENT(IN)      :: opt
REAL(KIND=8),INTENT(IN) :: const  ! valeur recherchée
REAL(KIND=8),INTENT(IN) :: x_init   ! proche du zero cherché
REAL(KIND=8),INTENT(IN) :: tol_max  ! tolerance maximale
INTEGER,INTENT(IN)      :: iter_max ! nombre d'itérations maximales
REAL(KIND=8)            :: res
! REAL(KIND=8),EXTERNAL:: f_analytique

REAL(KIND=8)        :: y0,dx2,x0
REAL(KIND=8)        :: x2,y2,errx,erry
INTEGER             :: iter
REAL(KIND=8)        :: tol
x0=x_init
tol=1.d0
iter=0

! calcul de la dérivée analytique

!y0=(A_Acritic(x0)-const)**2
y0=f_analytique(opt,0,x0,const)
!write(6,*)'x0=',x0,' y0=',y0,'cons= ',const,'A/A* =',A_Acritic(x0)
! PRINT'(a,4x,10(a,12x))','iter', 'tol',' x ', 'errx','erry'
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
          write(6,*) 'pas de convergence : stop'
          STOP
END IF

!write(6,*) 'Solution : x=  ',x0, 'f(x) = ', y0
res=x0
END FUNCTION newton_new
!***********************************
function f_analytique(opt,f_type,x,const)  result(f)
!***********************************
!> @details
!!f_type = 0 : function
!!f_type = 1 : derivative
!! A/A_crit=F(Mach)

IMPLICIT NONE
REAL(KIND=8), INTENT(IN)    :: x,const
INTEGER,INTENT(IN)          :: opt,f_type
REAL(KIND=8)                :: f,tmp
! REAL(KIND=8),EXTERNAL       ::der_A_Acritic,A_Acritic 
 
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
implicit none
integer,parameter           :: n=201
real(kind=8)                :: Ma,dMa,tmp
! real(kind=8),external       :: der_A_Acritic,der2_A_Acritic,A_Acritic
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
!! dérivée de la fonction de Fanno / Mach
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
!! fonction de Fanno 
use Constante
implicit none
real(kind=8)    :: Ma,res
real(kind=8)    :: omega
omega=1.d0+(gam-1.d0)/2.d0*Ma**2
!res=-(gam+1.d0)/(2.d0*gam)*log(Ma**2/omega)-1.d0/(gam*Ma**2)
res=(1.d0-Ma**2)/(gam*Ma**2)+(gam+1)/(2.d0*gam)*log(0.5d0*(gam+1.d0)*Ma**2/omega)
end function Fanno


!*********************************
function Entropie_Fanno(Ma) result(res)
!*********************************
!> @details
!! s=f(Mach) pour le problème de Fanno
use Constante
implicit none
real(kind=8)    :: Ma,res
real(kind=8)    :: omega,tmp
omega=1.d0+(gam-1.d0)/2.d0*Ma**2
tmp=r*(gam+1)/(2.d0*(gam-1))*log((gam+1)/2.d0)
res=r*(log(Ma)-(gam+1)*log(omega)/(2*(gam-1)))+tmp
end function Entropie_Fanno

!*********************************
function der_Entropie_Fanno(Ma) result(res)
!*********************************
!> @details
!! ds/dMach=df(Mach)/dMach pour le problème de Fanno
use Constante
implicit none
real(kind=8)    :: Ma,res
real(kind=8)    :: omega
omega=1.d0+(gam-1.d0)/2.d0*Ma**2
res=r*(1.d0-Ma**2)/(Ma*omega)
end function der_Entropie_Fanno

!*********************************
function Sonic_Length(Ma) result(res)
!*********************************
!> @details
!! Fanno : Calcul de la longueur sonique
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
    write(1,'(6(e15.8,3x))')Ma,Fanno(Ma),der_Fanno(Ma),Entropie_Fanno(Ma)/Cp, &
         der_Entropie_Fanno(Ma)/Cp,Sonic_Length(Ma)
end do
close(1)
end subroutine plot_Fanno

!**************************************************
subroutine sonic_ratio_Fanno(Ma,r_T,r_P,r_rho,r_Pi)
!**************************************************
!> @details
!! rapport de T, P,rho et  pi par rapport aux grandeurs critiques
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
!! fonction réciproque de la fonction de Fanno
! if x_init < 1 : subsonic case
! if x_init > 1 : supersonic case
implicit none
real(kind=8),intent(in)  :: RHS,x_init
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

! calcul de la dérivée analytique

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
          write(6,*) 'pas de convergence : stop'
          STOP
END IF

!write(6,*) 'Solution : x=  ',x0, 'f(x) = ', y0
Mach=x0
if (Mach.lt.0.d0) then
     write(6,*)'Convergence vers une mauvaise valeur'
     write(6,*)'Il faut rapprocher la condition initiale (le Mach) de la solution'
     stop
end if
 
end function inverse_Fanno


!*********************************
function qm_Mach(Mach) result(res)
!*********************************
!> @details
!! part of the Mach dependance in the mass flow rate formula
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
implicit none
real(kind=8),intent(in)         :: M
integer,intent(in)              :: opt
real(kind=8)                    :: F
!> opt = 1: Rayleigh function
!! opt = 2: Derivative of Rayleigh function

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
!! réciproque de la fonction de Rayleigh
!! if x_init < 1 : subsonic case
!! if x_init > 1 : supersonic case
implicit none
real(kind=8),intent(in)  :: RHS,x_init
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

! calcul de la dérivée analytique

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
          write(6,*) 'pas de convergence : stop'
          STOP
END IF

!write(6,*) 'Solution : x=  ',x0, 'f(x) = ', y0
Mach=x0
if (Mach.lt.0.d0) then
     write(6,*)'Convergence vers une mauvaise valeur'
     write(6,*)'Il faut rapprocher la condition initiale (le Mach) de la solution'
     stop
end if
 
end function inverse_Rayleigh



!**********************************************************
subroutine Rayleigh_flow(q)
!**********************************************************
!> @details
!! ratios to critical state for the Rayleigh flow
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
subroutine plot_Rayleigh_entropie
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

end subroutine plot_Rayleigh_entropie


!********************************************************************
include 'tables_perso.f90'
!********************************************************************

end module mod_EC


 



