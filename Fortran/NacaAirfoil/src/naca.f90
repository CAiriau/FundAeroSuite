!------------------------------------------------------------------------------
! TITLE         : NACA AIRFOIL GEOMETRY
! PROJECT       : FundAeroSuite
! MODULE        : naca_airfoil_design (version with Drela's subroutines)
! URL           : 
! AFFILIATION   : Paul Sabatier University, Toulouse
! DATE          : 2016 - 2020
! REVISION      : 2020 V2
!> @author
!> Christophe Airiau
!
! DESCRIPTION:
!>  get geometry of NACA ' and 5 digits (main subroutine from Drela f77)
!!
!! subroutine translated from f77 to f90 in 2016
!------------------------------------------------------------------------------
program  test_profil_naca
implicit none
integer,parameter               :: n = 201     ! nombre de points sur l'extrados et sur l'intrados
real(kind=8),dimension(n)       :: x,yt,yc
real(kind=8),dimension(2*n)     :: xb,yb
integer                         :: nb,i
character(len=20)               :: NomProfil='NACA-'

! test avec 4 digits 

! profil non symétrique
call NACA4(4412,x,yt,yc,n,xb,yb,nb,NomProfil)
print*,'nb                  = ', nb
print*,'Nom du profil       = ',trim(NomProfil)
call save_data

! profil symétrique
call NACA4(0012,x,yt,yc,n,xb,yb,nb,NomProfil)
print*,'nb                  = ', nb
print*,'Nom du profil       = ',trim(NomProfil)
call save_data

! profil  5 digits
call NACA5(23021,x,yt,yc,n,xb,yb,nb,NomProfil)
print*,'nb                  = ', nb
print*,'Nom du profil       = ',trim(NomProfil)
call save_data


contains

!*******************
subroutine save_data
!*******************
open(1,form='formatted',file=trim(NomProfil)//'.dat')
write(1,100) NomProfil
write(1,110)
do i=1,nb
    write(1,200) xb(i),yb(i)
end do
close(1)
open(1,form='formatted',file=trim(NomProfil)//'_carac.dat')
write(1,100) NomProfil
write(1,120)
do i=1,n
    write(1,200) x(i),yc(i),yt(i)
end do
close(1)


100 format('#',50('*'),/,'#',4x,a10,/,'#',50('*'))
110 format('#',4x, 'x/c',10x,'y/c')
120 format('#',4x, 'x/c',7x,'cambrure/c', 4x,'epaisseur/c')
200 format(3(e12.5,2x))

end subroutine save_data

end program test_profil_naca


!***********************************************************************
!    Module:  naca.f
!    traduit en f90, C. Airiau, Janv. 2016
! 
!    Copyright (C) 2000 Mark Drela 
! 
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!***********************************************************************

!> @brief
!!  NACA four digits calculator
!!  @author Mark Drela, f77 -> f90 : C. Airiau
!**********************************************************
SUBROUTINE NACA4(IDES,XX,YT,YC,NSIDE,XB,YB,NB,Profile_Name)
!**********************************************************
IMPLICIT NONE
INTEGER,INTENT(IN)              :: NSIDE,IDES
REAL(kind=8),dimension(NSIDE)   :: XX, YT, YC
REAL(kind=8),dimension(2*NSIDE) :: XB,YB
REAL(kind=8)                    :: M,P,T
CHARACTER(len=*)                :: Profile_Name
CHARACTER(len=10),parameter     :: DIGITS='0123456789'
REAL(kind=8)                    :: ANP
REAL(kind=8),parameter          :: AN=1.5d0 !---- TE point bunching parameter
INTEGER                         :: N1,N2,N3,N4, I,IB,NB
REAL(kind=8)                    :: FRAC

!
N4 =  IDES                             / 1000
N3 = (IDES - N4*1000                 ) / 100
N2 = (IDES - N4*1000 - N3*100        ) / 10
N1 = (IDES - N4*1000 - N3*100 - N2*10)
!
M = FLOAT(N4) / 100.d0
P = FLOAT(N3) / 10.d0
T = FLOAT(N2*10 + N1) / 100.d0
!
ANP = AN + 1.d0
DO  I=1, NSIDE
    FRAC = FLOAT(I-1)/FLOAT(NSIDE-1)
    IF(I.EQ.NSIDE) THEN
        XX(I) = 1.0d0
    ELSE
        XX(I) = 1.0d0 - ANP*FRAC*(1.0d0-FRAC)**AN - (1.0d0-FRAC)**ANP
    ENDIF
    YT(I) = ( 0.29690d0*SQRT(XX(I)) -0.12600d0*XX(I) -0.35160d0*XX(I)**2&
          + 0.28430d0*XX(I)**3 - 0.10150d0*XX(I)**4) * T / 0.20d0
    IF(XX(I).LT.P) THEN
        YC(I) = M/P**2 * (2.0d0*P*XX(I) - XX(I)**2)
    ELSE
        YC(I) = M/(1.0d0-P)**2 * ((1.0d0-2.0d0*P) + 2.0d0*P*XX(I)-XX(I)**2)
    ENDIF
END DO
!
IB = 0
DO I=NSIDE, 1, -1
    IB = IB + 1
    XB(IB) = XX(I)
    YB(IB) = YC(I) + YT(I)
END DO
DO  I=2, NSIDE
    IB = IB + 1
    XB(IB) = XX(I)
    YB(IB) = YC(I) - YT(I)
END DO
NB = IB
!
Profile_Name = 'NACA-'
Profile_Name(6:9) = DIGITS(N4+1:N4+1)// DIGITS(N3+1:N3+1)// DIGITS(N2+1:N2+1)// DIGITS(N1+1:N1+1)
END SUBROUTINE NACA4

!> @brief
!!  NACA five digits calculator
!!  @author Mark Drela, f77 -> f90 : C. Airiau
!**********************************************************
SUBROUTINE NACA5(IDES,XX,YT,YC,NSIDE,XB,YB,NB,Profile_Name)
!**********************************************************

IMPLICIT NONE
INTEGER,INTENT(IN)              :: IDES
INTEGER,INTENT(IN)              :: NSIDE
REAL(kind=8),dimension(NSIDE)   :: XX, YT, YC
REAL(kind=8),dimension(2*NSIDE) :: XB,YB
REAL(kind=8)                    :: M,P,T,C      ! P unused, just for airfoil characteristic
CHARACTER(len=*)                :: Profile_Name
CHARACTER(len=10),parameter     :: DIGITS='0123456789'
REAL(kind=8)                    :: ANP
REAL(kind=8),parameter          :: AN=1.5d0 !---- TE point bunching parameter
INTEGER                         :: N1,N2,N3,N4,N5,N543, I,IB,NB
REAL(kind=8)                    :: FRAC
!
N5 =  IDES                                        / 10000
N4 = (IDES - N5*10000                           ) / 1000
N3 = (IDES - N5*10000 - N4*1000                 ) / 100
N2 = (IDES - N5*10000 - N4*1000 - N3*100        ) / 10
N1 = (IDES - N5*10000 - N4*1000 - N3*100 - N2*10)
!
N543 = 100*N5 + 10*N4 + N3
!
IF (N543 .EQ. 210) THEN
    !     P = 0.05
    M = 0.0580d0
    C = 361.4d0
ELSE IF (N543 .EQ. 220) THEN
    !     P = 0.10
    M = 0.1260d0
    C = 51.64d0
ELSE IF (N543 .EQ. 230) THEN
    !     P = 0.15
    M = 0.2025d0
    C = 15.957d0
ELSE IF (N543 .EQ. 240) THEN
    !     P = 0.20
    M = 0.2900d0
    C = 6.643d0
ELSE IF (N543 .EQ. 250) THEN
    !     P = 0.25
    M = 0.3910d0
    C = 3.230d0
ELSE
    WRITE(*,*) 'Illegal 5-digit designation'
    WRITE(*,*) 'First three digits must be 210, 220, ... 250'
    RETURN
ENDIF
!
T = FLOAT(N2*10 + N1) / 100.0
!
ANP = AN + 1.0d0
DO I=1, NSIDE
    FRAC = FLOAT(I-1)/FLOAT(NSIDE-1)
    IF(I.EQ.NSIDE) THEN
        XX(I) = 1.0d0
    ELSE
        XX(I) = 1.0d0 - ANP*FRAC*(1.0d0-FRAC)**AN - (1.0d0-FRAC)**ANP
    ENDIF
    !
    YT(I) = ( 0.29690d0*SQRT(XX(I))- 0.12600d0*XX(I)-0.35160d0*XX(I)**2&
          + 0.28430d0*XX(I)**3 - 0.10150d0*XX(I)**4) * T / 0.20d0
    IF(XX(I).LT.M) THEN
        YC(I) = (C/6.0d0) * (XX(I)**3 - 3.0d0*M*XX(I)**2 + M*M*(3.0d0-M)*XX(I))
    ELSE
        YC(I) = (C/6.0d0) * M**3 * (1.0d0 - XX(I))
    ENDIF
END DO
!
IB = 0
DO  I=NSIDE, 1, -1
    IB = IB + 1
    XB(IB) = XX(I)
    YB(IB) = YC(I) + YT(I)
END DO
DO  I=2, NSIDE
    IB = IB + 1
    XB(IB) = XX(I)
    YB(IB) = YC(I) - YT(I)
END DO
NB = IB
!
Profile_Name = 'NACA-'
Profile_Name(6:10) = DIGITS(N5+1:N5+1)//DIGITS(N4+1:N4+1)//DIGITS(N3+1:N3+1) &
                   // DIGITS(N2+1:N2+1)// DIGITS(N1+1:N1+1)
!
END SUBROUTINE NACA5
