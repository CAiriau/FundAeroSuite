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
!! MODULE : concerne les entrées en lignes de commandes
!! @author C. Airiau
!!
!! DATE :  May 2014
!!
!! 1) read the inline argument of the command line to modify program parameters values
!!
!! 2) affect value to parameters depending on the parameter type.

module mod_read
    ! 
    ! see the driver at the end of the file

    !

implicit none
!  private
!  save
!  public:: save_par  
!  public:: nparameters,Narg
!  public:: par

interface save_par
        module procedure save_par_real,save_par_integer,save_par_character
end interface save_par

interface s_to_val
        module procedure s_to_integer, s_to_real
end interface s_to_val

!>@details
!! paramètres sur une ligne de commande
type param
    character(len=1)    :: t='-'            ! type : i,r,c
    character(len=15)   :: nom='empty'      ! name
    integer             :: ent=0            ! integer value
    real(kind=8)        :: ree=0.d0         ! real value
    character(len=15)   :: car='-'          ! character value
end type param

integer,parameter               :: n_par_max=20
type(param),dimension(n_par_max):: par
integer                         :: nparameters
integer                         :: Narg

contains


!**************************************
subroutine save_par_real(i,typ,nom,val)
!**************************************
    integer,          intent(in) :: i
    character(len=1), intent(in) :: typ
    character(len=*),intent(in) :: nom
    real(kind=8),     intent(in) :: val
    par(i)%t=typ ; par(i)%nom=nom ; par(i)%ree=val
end subroutine save_par_real

!*****************************************
subroutine save_par_integer(i,typ,nom,val)
!*****************************************
    integer,          intent(in) :: i
    character(len=1), intent(in) :: typ
    character(len=*),intent(in) :: nom
    integer,          intent(in) :: val
    par(i)%t=typ ; par(i)%nom=nom ; par(i)%ent=val
end subroutine save_par_integer

!*******************************************
subroutine save_par_character(i,typ,nom,val)
!*******************************************
    integer,          intent(in) :: i
    character(len=1), intent(in) :: typ
    character(len=*),intent(in) :: nom,val
    par(i)%t=typ ; par(i)%nom=nom ; par(i)%car=val
end subroutine save_par_character



!================================================
subroutine s_to_integer ( s, ival, ierror, last )
!================================================
!>@details
!*****************************************************************************80
!! S_TO_I4 reads an I4 from a string.
!    John Burkardt
!  Parameters:
!    Input, character ( len = * ) S, a string to be examined.
!    Output, integer IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!    Output, integer LAST, the last character of S used 
!    to make IVAL.
!
implicit none

character           :: c
integer             :: i, ierror,isgn,istate,ival,last
character(len=*):: s

ierror = 0; istate = 0; isgn = 1; ival = 0

do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1; isgn = -1
      else if ( c == '+' ) then
        istate = 1; isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2; ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1; return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2; ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1; return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival; last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival; last = len_trim ( s )
  else
    ierror = 1; last = 0
  end if

end subroutine s_to_integer

!============================================
subroutine s_to_real ( s, r, ierror, lchar )
!============================================
!>@details
!*****************************************************************************80
!! S_TO_R8 reads an R8 from a string.
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Examples:
!
!    S                 R
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2D-1'           0.2
!    '23.45'           23.45
!    '-4.2D+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal number.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) R, the value that was read from the string.
!
!    Output, integer IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
implicit none

character         :: c
!logical           :: ch_eqi
integer           :: ierror,ihave,isgn,iterm,jbot,jsgn,jtop,lchar,nchar,ndig
real(kind=8)      :: r,rbot,rexp,rtop
real(kind=8),parameter::ten=10.D0
character ( len = * ) s
character, parameter :: TAB = char ( 9 )

nchar = len_trim ( s )
ierror = 0
r = 0.0D+00
lchar = - 1; isgn = 1; rtop = 0.0D+00; rbot = 1.0D+00
jsgn = 1; jtop = 0; jbot = 1; ihave = 1; iterm = 0

do

    lchar = lchar + 1; c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        lchar = lchar + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = - 1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = - 1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit2 ( c, ndig )

      if ( ihave == 3 ) then
        rtop = ten * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = ten * rtop + real ( ndig, kind = 8 )
        rbot = ten * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 .or. nchar <= lchar + 1 ) then
      exit
    end if

end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
if ( iterm /= 1 .and. lchar + 1 == nchar ) then
    lchar = nchar
end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave
    return

end if
!
!  Number seems OK.  Form it.
!
if ( jtop == 0 ) then
    rexp = 1.0D+00
else

    if ( jbot == 1 ) then
      rexp = ten**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = ten**rexp
    end if

end if

r = isgn * rexp * rtop / rbot

end subroutine s_to_real 

!============================================
subroutine ch_to_digit2 ( c, digit )
!============================================

!*****************************************************************************80
!>@details
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the corresponding integer value.  
!    If C was 'illegal', then DIGIT is -1.
!
implicit none

character   :: c
integer     :: digit

if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
    digit = ichar ( c ) - 48
else if ( c == ' ' ) then
    digit = 0
else
    digit = -1
end if

end subroutine ch_to_digit2 

!************************
subroutine ch_cap ( c )
!************************
!>@details
!
!! CH_CAP capitalizes a single character.
!
!  Modified: C. Airiau, august 2015
!  Author:   John Burkardt
!  Parameters:   Input/output, character C, the character to capitalize.
!
implicit none

character (len=1)   :: c
integer             :: itemp
itemp = ichar ( c )
if ( 97 <= itemp .and. itemp <= 122 )  c = char ( itemp - 32 )
end subroutine ch_cap

!**************************
function ch_eqi ( c1, c2 )
!**************************
!>@details
!*****************************************************************************80
!> CH_EQI is a case insensitive comparison of two characters for equality.
!  Examples:    CH_EQI ( 'A', 'a' ) is .TRUE.
!  Modified:   august 2015 by C. Airiau
!  Author:   John Burkardt
!  Parameters:
!    Input, character C1, C2, the characters to compare.
!    Output, logical CH_EQI, the result of the comparison.
!
implicit none

character(len=1),intent(in)   :: c1,c2
character(len=1)              :: c1_cap, c2_cap
logical ch_eqi

c1_cap = c1; c2_cap = c2

call ch_cap ( c1_cap ); call ch_cap ( c2_cap )

if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
else
    ch_eqi = .false.
end if

end function ch_eqi


!**********************
subroutine get_argument
!**********************
!>@details : pour obtenir un argument du clavier
character(len=15),allocatable,dimension(:)  :: arg
integer             :: ierror,info,i,k
logical             :: found

Narg = command_argument_count()
if (Narg.le.1) then
   write(6,*) 'no argument'; return
endif
allocate(arg(Narg))
!write(6,*) 'number of arguments : ',Narg
i=0
do  while (i.lt.Narg)
    found=.false.
    i=i+1
    call get_command_argument(i,arg(i))
    do k=1,nparameters
        if (arg(i).eq.trim(par(k)%nom)) then
            !write(6,*) 'argument in the list : N° ',i,'  arg = ',arg(i)
            i=i+1; found=.true.
            call get_command_argument(i,arg(i))
            select case (par(k)%t)
                case ('i')
                    call s_to_val(arg(i),par(k)%ent,ierror,info)
                case ('r')
                    call s_to_val(arg(i),par(k)%ree,ierror,info)
                case ('c')
                    par(k)%car=arg(i)
            end select
            cycle
        end if
    end do

    if (.not.found) then
           write(6,*) 'argument not in the list : N° ',i,'  arg = ',arg(i)
    end if
end do

deallocate(arg)
end subroutine get_argument

!**********************
subroutine get_argument2
!**********************
!>@details 
!! pour obtenir un argument du clavier
implicit none
character(len=15),allocatable,dimension(:)  :: arg
integer             :: ierror,info,i

ierror=0
Narg = command_argument_count()
if (Narg==0) then
   write(6,*) 'no argument'; return
endif
allocate(arg(Narg))
!write(6,*) 'number of arguments : ',Narg

i=0;
do  while ((i.lt.Narg).and.(ierror.eq.0))
    i=i+1
    call get_command_argument(i,arg(i))
    call s_to_val(arg(i),par(i)%ent,ierror,info)
end do
deallocate(arg)

if (ierror.eq.1)  stop'bad argument in commande line'

end subroutine get_argument2

 
!*******************************
subroutine driver_test_argument
!*******************************
!>@details
!! driver pour tester l'entrer des arguments en ligne de commande
implicit none
integer i

write(6,*) 'Procedure to get argument and affect the values'

nparameters=4       ! number of interesting parameters
! name and initial values
par(1)=param('c','cvar',0,0.d0, 'contenu')
par(2)=param('r','xvar',0,  3.141592, '')
par(3)=param('i','ivar',-1, 0, '')

write(6,*) 'before test'
do i=1,nparameters;write(6,*) par(i); end do

call get_argument

write(6,*) 'after test'
do i=1,nparameters;write(6,*) par(i); end do

end subroutine driver_test_argument


!*************************************************************
subroutine driver_aerodynamique(ask,opt_read,Mach,theta,ans,F)
!************************************************************
!>@details
!! driver spécifique pour le code EC 
implicit none
integer,optional        :: opt_read
! integer                 :: i
logical,optional        :: ask
real(kind=8),optional   :: Mach,theta,ans,F

ask=.true.
!write(6,*) 'Procedure to get argument and affect the values'

nparameters=6       ! number of interesting parameters
! name and initial values
par(1)=param('i','opt=',0,0.d0, '')
par(2)=param('r','M=',0,  2.0, '')
par(3)=param('r','angle=',0,  10.0, '')
par(4)=param('r','A/Ac=',0,  0.0, '')
par(5)=param('r','F=',0,  0.0, '')
par(nparameters)=param('c','h',0,  0.0, 'Aide')

!write(6,*) 'options par défaut'
!do i=1,nparameters;write(6,*) par(i); end do

call get_argument

!  write(6,*) 'Après commande :'
!   do i=1,nparameters;write(6,100) par(i); end do
!   100 format(a1,3x,a,3x,i3,3x,e12.5,2x,a)
!write(6,*) 'nombre argument: ', Narg
select case (Narg)
case (0)
    opt_read=-1
case(1)
    opt_read=-2
case(2)
    if (par(1)%ent.ge.0) then
        opt_read=par(1)%ent
    else
        opt_read=-1
    end if
case default
    opt_read=par(1)%ent
    Mach=par(2)%ree
    theta=par(3)%ree
    ans=par(4)%ree
    F=par(5)%ree
    ask=.false.
end select

!write(6,*) 'opt_read= ',opt_read
!write(6,*) par(1)


end subroutine driver_aerodynamique

end module mod_read
