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

!**********************
subroutine tables_perso
!**********************
use Constante
implicit none
real(kind=8)            :: dM,M_init,M_final,Mach
integer                 :: nptM     ! nombre de points
integer,dimension(4)    :: npt      ! nombre de points
integer                 :: nc       ! nombre de colonnes
integer                 :: nl       ! nombre de lignes
integer                 :: nl_super ! nombre de lignes par pages pour le supersonique
integer                 :: reste,i,nl1,k,j,i1,i2,i_tmp
integer                 :: nt       ! nombre de cas de tables pour le supersonique
integer                 :: nf       ! nombre de fichiers pour les tables en supersonique
real(kind=8),allocatable,dimension(:,:) :: table(:,:)
!real(kind=8), external            :: T_Ti,P_Pi,rho_rhoi,S_sur_Sc
!real(kind=8),external            :: P2_P1,rho2_rho1,Pi2_Pi1,mach_downstr,omega
real(kind=8),dimension(4)   :: dM_tmp,M_init_tmp,M_final_tmp
!character(len=2)        :: charac
logical                 :: test
logical,parameter       :: opt_super_table=.false.

!  opt_super_table  = .true.    # on fixe le pas et l'intervalle de Mach pour chaque table
!  opt_super_table  = .false.   # on fixe la dimension du tableau, le pas et le Mach initial

write(6,*) 'tables personnelles'

!=====================================
write(6,*) 'table pour le subsonique'
!=====================================

nc = 2;                     ! nombre de colonnes dans le fichier
k  = 6;                     ! nombre de paramètres de display_outputs

write(6,*) 'nombre de colonnes choisi                 :  ', nc
M_init=0.0d0; M_final=1.0d0; dM=0.01d0          ! paramètre de la table pour le subsonique
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

! ************************************************
!   SUBSONIQUE : 2 colonnes : 1 fichier
! ************************************************
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

! ************************************************
!   SUBSONIQUE : 1 colonne : 2 fichiers
! ************************************************
!
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
k=7             ! nombre de display_outputs
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

! ******************************************************
!   SUPERSONIQUE ISENTROPIQUE  : 2 colonnes : 1 fichier
! ******************************************************
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


! **************************************************************
!   SUPERSONIQUE ISENTROPIQUE  :  1 colonne : plusieurs fichiers
! **************************************************************
!

if (opt_super_table) then 
! valeurs originales
    nt=3
    nf=4
    dM_tmp=(/0.02d0, 0.05d0, 1.0d0, 0.0d0/)
    M_init_tmp =(/1.d0, 3.d0,  5.d0, 0.d0/)
    M_final_tmp=(/3.d0, 5.d0, 30.d0, 0.d0/)
    nl_super=40
    npt(4)=0
    do j=1,nt
        write(6,*)'j, ',j,' step ', ((M_final_tmp(j)-M_init_tmp(j))/dM_tmp(j)), M_final_tmp(j),M_init_tmp(j),dM_tmp(j)
        npt(j)=int((M_final_tmp(j)-M_init_tmp(j))/dM_tmp(j))
    end do
else
    nt=4
    nf=4
    nl_super=50
    dM_tmp=(/0.005d0, 0.005d0, 0.005d0, 0.01d0/)
    M_init_tmp =(/1.d0, 1.25d0,  1.5d0, 1.75d0/)
    M_final_tmp=(/0d0, 0.d0, 0.d0, 0.d0/)
    do j=1,nt
        M_final_tmp(j)=M_init_tmp(j)+float(nl_super-1)*dM_tmp(j)
        npt(j)=nl_super
        write(6,*)'j, ',j,' step ', ((M_final_tmp(j)-M_init_tmp(j))/dM_tmp(j)), M_final_tmp(j),M_init_tmp(j),dM_tmp(j)
    end do

end if

nptM=sum(npt)+1; 

k=6         !  nombres de sortie
allocate(table(nptM,k))
write(6,*) 'npt ',npt(:), ' total = ', sum(npt)


if (opt_super_table) then 
    i1=1
    table(i1,1)=M_init_tmp(1)
    write(6,*)'mach initial                               : ',table(1,1)
    do j=1,nt
        i2=i1+npt(j)
        write(30,*)' Mach init et final et pas',M_init_tmp(j), M_final_tmp(j),dM_tmp(j)
        do i=i1+1,i2
            table(i,1)=table(i-1,1)+dM_tmp(j)
            write(30,'(i4,3x,f12.3,5x,i3)')i,table(i,1),i-i1
        end do
        i1=i2
    end do
else
    i1=1
    do j=1,nt
        write(6,*) 'i1= ',i1
        do i=1,nl_super
            table(i1+i-1,1)=M_init_tmp(j)+float(i-1)*dM_tmp(j)
        end do
        i1=i1+nl_super
    end do
end if


do i=1,nptM
    Mach=table(i,1)
    table(i,2)=T_Ti(Mach);    table(i,3)=P_Pi(Mach)
    table(i,4)=rho_rhoi(Mach);  table(i,5)=S_sur_Sc(Mach)
    table(i,6)=Fanno(Mach); 
end do

do j=1,nf
    open(10,form='formatted',file='supersonique'//charac(2,j)//'.tex')
    write(10,111); write(10,130)'\\ \hline \hline'
    do i=1,nl_super
        i_tmp=i+(j-1)*nl_super
        test=(table(i_tmp,5).lt.100.d0)
        if (minval(table(i_tmp,2:4)).gt.0.1d0) then
            if (test) then
                write(10,103) table(i+(j-1)*nl_super,1:k)
            else
                write(10,104) table(i+(j-1)*nl_super,1:k)
            end if
        else
            if (test) then
                write(10,105) table(i+(j-1)*nl_super,1:k)
            else
                write(10,106) table(i+(j-1)*nl_super,1:k)
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
    table(i,2)=mach_downstr(Mach)              ! mach normal downstr
    table(i,3)=P2_P1(Mach)             ! rapport des pressions statiques
    table(i,4)=T_Ti(table(i,2))/T_Ti(Mach)  ! rapport des températures
    table(i,5)=rho2_rho1(Mach)           ! rapport des masses volumiques
    table(i,6)=1.d0-Pi2_Pi1(Mach)            ! rapport des pressions isentropiques
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

if (opt_super_table) then 
    i1=1
    table(i1,1)=M_init_tmp(1)
    write(6,*)'mach initial                               : ',table(1,1)
    do j=1,nt
        i2=i1+npt(j)
        do i=i1+1,i2
            table(i,1)=table(i-1,1)+dM_tmp(j)
        end do
        i1=i2
    end do
else
    i1=1
    do j=1,nt
        write(6,*) 'i1= ',i1
        do i=1,nl_super
            table(i1+i-1,1)=M_init_tmp(j)+float(i-1)*dM_tmp(j)
        end do
        i1=i1+nl_super
    end do
end if

do i=1,nptM
    Mach=table(i,1)
    table(i,2)=mach_downstr(Mach)              ! mach normal downstr
    table(i,3)=P2_P1(Mach)             ! rapport des pressions statiques
    table(i,4)=T_Ti(table(i,2))/T_Ti(Mach)  ! rapport des températures
    table(i,5)=rho2_rho1(Mach)           ! rapport des masses volumiques
    table(i,6)=1.d0-Pi2_Pi1(Mach)            ! 1-rapport des pressions isentropiques
end do

do j=1,nf
    open(10,form='formatted',file='chocs'//charac(2,j)//'.tex')
    write(10,311); write(10,330)'\\ \hline \hline'

    do i=1,nl_super
        i_tmp=i+(j-1)*nl_super
        test=(table(i_tmp,3).lt.100.d0)
        if (minval(table(i_tmp,2:k)).gt.0.1d0) then
            if (test) then
                write(10,203) table(i+(j-1)*nl_super,1:k)
            else
                write(10,204) table(i+(j-1)*nl_super,1:k)
            end if
        else
            if (test) then
                write(10,205) table(i+(j-1)*nl_super,1:k)
            else
                write(10,206) table(i+(j-1)*nl_super,1:k)
            end if
        end if
        if (mod(i,5).eq.0) write(10,*) '\hline'
    end do
    write(10,120)
    close(10)
end do

deallocate(table)

!*****************************************
! Table pour Prandtl Meyer
!*****************************************

write(6,*)' TABLE DE PRANDTL-MEYER'

if (opt_super_table) then 

    write(6,*) 'cas de référence'
    dM_tmp=(/0.01d0, 0.02d0, 0.5d0, 0.0d0/)
    M_init_tmp =(/1.d0, 3.5d0,  5.d0, 0.d0/)
    M_final_tmp=(/3.5d0, 5.d0, 30.d0, 0.d0/)
    nt=3
    nl_super=40
    nl=nl_super
    npt(4)=0
    do j=1,nt
        npt(j)=int((M_final_tmp(j)-M_init_tmp(j))/dM_tmp(j))
        write(6,*)'j, ',j,' step ', ((M_final_tmp(j)-M_init_tmp(j))/dM_tmp(j)), M_init_tmp(j), &
                                    M_final_tmp(j),dM_tmp(j)
    end do

else

    dM_tmp     =(/0.005d0, 0.005d0, 0.01d0, 0.2d0/)
    M_init_tmp =(/1.d0   , 1.75d0 , 2.5d0  , 4.0d0/)
    M_final_tmp=(/0d0    , 0.d0   , 0.d0   , 0.d0/)
    nt=4
    nl_super=3*50
    nl=50
    do j=1,nt
        M_final_tmp(j)=M_init_tmp(j)+float(nl_super-1)*dM_tmp(j)
        npt(j)=nl_super
        write(6,*)'j, ',j,' step ', ((M_final_tmp(j)-M_init_tmp(j))/dM_tmp(j)), M_final_tmp(j), &
                                M_init_tmp(j),dM_tmp(j)
    end do

end if

nptM=sum(npt)+1; 
write(6,*) 'nombre de points                          : ', nptM
k=3
allocate(table(nptM,k))
write(6,*) 'npt ',npt(:), ' total = ', sum(npt)

if (opt_super_table) then 
    i1=1
    table(i1,1)=M_init_tmp(1)
    write(6,*)'mach initial                               : ',table(1,1)
    do j=1,nt
        i2=i1+npt(j)
        write(31,*)' Mach init et final et pas',M_init_tmp(j), M_final_tmp(j),dM_tmp(j)
        do i=i1+1,i2
            table(i,1)=table(i-1,1)+dM_tmp(j)
            write(31,'(i4,3x,f12.3,5x,i3)')i,table(i,1),i-i1
        end do
        i1=i2
    end do
else
    i1=1
    do j=1,nt
        write(6,*) 'i1= ',i1
        do i=1,nl_super
            table(i1+i-1,1)=M_init_tmp(j)+float(i-1)*dM_tmp(j)
        end do
        i1=i1+nl_super
    end do
end if

do i=1,nptM
    Mach=table(i,1)
    table(i,2)=omega(Mach)*rad2deg             ! mach normal downstr
    table(i,3)=asin(1.d0/Mach)*rad2deg
end do

write(6,*) ' nt = ',nt,'  nl_super =', nl_super

do j=1,nt      ! nombres de pages, 3 colonnes
    open(10,form='formatted',file='omega'//charac(2,j)//'.tex')
    write(10,411); 
    write(10,430) ' & '; write(10,430) ' & '; write(10,430)'\\ \hline \hline'

    do i=1,nl           ! chaque ligne du tableau
        i_tmp=i+(j-1)*3*nl
        write(10,403) table(i_tmp,1:k),table(i_tmp+nl,1:k),table(i_tmp+2*nl,1:k)
        if (mod(i,5).eq.0) write(10,*) '\hline'
    end do
    write(10,120)
    close(10)
end do
deallocate(table)

!*****************************************************************************
!               FORMATS
!*****************************************************************************
100 format((f6.3,' & ',f6.4,' & ',f7.4,' & ',f7.4,' & ',f7.4,' & ', f7.4), &
            f6.3,' & ',f6.4,' & ',f7.4,' & ',f7.4,' & ',f7.4,' & ', f7.4,'\\')
101 format((f6.3,' & ',f6.4,' & ',f7.4,' & ',f7.4,' & ',f7.4,' & '), &
                5x,  ' & ',6x ,' & ',7x, ' & ',7x, ' & \\')
102 format(f6.3,' & ',f7.5,' & ',f7.5,' & ',f7.5,' & ',f9.4,' & ',f9.4, '\\')
103 format(f6.3,' & ',f7.5,' & ',f7.5,' & ',f7.5,' & ',f9.4,' & ',f9.4, '\\')
104 format(f6.3,' & ',f7.5,' & ',f7.5,' & ',f7.5,' & ',e12.5,' & ',e12.5, '\\')
105 format(f6.3,' & ',es12.4,' & ',es12.4,' & ',es12.4,' & ',f9.4,' & ',f9.4, '\\')
106 format(f6.3,' & ',es12.4,' & ',es12.4,' & ',es12.4,' & ',es12.4,' & ',es12.4, '\\')
110 format('{\small \begin{tabular}{',2(6('|c'),'|')'}', '\hline')
111 format('{\small \begin{tabular}{',(6('|c'),'|')'}', '\hline')
120 format('\hline \end{tabular}}')
130 format('$M$ &  $T/T_i$  &  $P/P_i$ & $\rho / \rho_i$ & $S / S_c$ & $Fa$',a)
200 format((f5.2,' & ',f7.4,' & ',f7.4,' & ',f7.4,' & ',3(f7.4,' & ')), &
            f5.3,' & ',f7.4,' & ',f7.4,' & ',f7.4,' & ',f7.4,' & ', &
            f9.2,'&',f9.2,'\\')
201 format((f5.3,' & ',f6.4,' & ',f7.5,' & ',f7.5,' & ',3(f7.5,' & ')), &
            5x,  ' & ',6x ,' & ',7x, ' & ',7x, ' & ',7x,' & ', 7x, ' &   \\')
203 format(f6.3,' & ',f7.5,3(' & ',f9.4),' & ',f7.5, '\\')
204 format(f6.3,' & ',f7.5,3(' & ',e12.5),' & ',f7.5, '\\')
205 format(f6.3,' & ',es12.4,3(' & ',f9.4),' & ',es12.4, '\\')
206 format(f6.3,' & ',es12.4,3(' & ',es12.4),' & ',es12.4, '\\')
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
330 format('$Mn_1$ &  $Mn_2$  &  $P_2/P_1$ & $ T_2 / T_1$ & $\rho_2/ \rho_1$',&
            '&  $1-P_{i2}/P_{i1} $',a)
411 format('{\small \begin{tabular}{',3('| c | c | c | '),' } \hline')
430 format('$M_0$ &  $\omega \ (^\circ)$  &  $\mu \ (^\circ)$ ',a)
403 format(2(f6.3,' & ',f10.2,' & ',f10.2, ' & '), &
            f6.3,' & ',f10.2,' & ',f10.2,2x,' \\')
end subroutine tables_perso


