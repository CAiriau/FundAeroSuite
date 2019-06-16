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


!**************************
subroutine examen_mai_2014
!**************************
! triangle isocèle en incidence
nz=5
allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
! Mai 2014 
Ti(0)=280;P(0)=10000.d0; M(0)=3.d0; theta(0)=10.d0*deg2rad
! question f i)
write(6,100) 0
call calcul_zone_uniforme(0)
! question f ii)
! zone 1 : il ne se passe rien P1=P0, etc ...
write(6,100) 1
call rien (0,1)
theta(2)=2.d0*theta(0)
! question f iii) zone 2 : détente isentropique de 2 theta
write(6,100) 2
call detente_isentropique(1,2)
!
! question f iv) choc oblique en zone 3
write(6,100) 3
theta(3)=theta(0)
call choc_oblique(0,3)
! question f v) choc oblique en zone 4
write(6,100) 4
theta(4)=theta(2); 
call choc_oblique(2,4)
! question f V) zone 5 : détente isentropique de  theta
write(6,100) 5
theta(5)=theta(3)
call detente_isentropique(3,5)
!
! question f Vi) Cl Cd
call calcul_Cl_Cd

call sorties
deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn)


100 format(50('*'),/,'ZONE : ', i1,/,50('*'))
end subroutine examen_mai_2014



!**************************
subroutine examen_juin_2014
!**************************
! marche descendante
nz=2
allocate(M(0:nz),P(0:nz),Ti(0:nz),theta(0:nz),V(0:nz))
allocate(type_zone(0:nz),com(0:nz))
allocate(a(0:nz),tau(0:nz),T(0:nz),Pi(0:nz))
allocate(sigma(0:nz),Om(0:nz),Mn(0:nz))
! Mai 2014, marche descendante
Ti(0)=400; M(0)=2.d0; theta(0)=0.d0*deg2rad;P(0)=10000.d0;
write(6,*) 'Question 2.1'
tau(1)=0.5  ! rapport de pression entre 0 et 1
write(6,100) 0
call calcul_zone_uniforme(0)
write(6,100) 1
p(1)=tau(1)*p(0);pi(1)=pi(0);Ti(1)=Ti(0)
call calcul_ligne_glissement(0,1)
call calcul_zone_uniforme(1)
write(6,*) 'Question 2.2'
write(6,100) 2
theta(2)=0.d0; 
call choc_oblique(1,2)
deallocate(M,P,Ti,theta,V,com,type_zone,a,tau,T,Pi,sigma,Om,Mn)

100 format(50('*'),/,'ZONE : ', i1,/,50('*'))
end subroutine examen_juin_2014

