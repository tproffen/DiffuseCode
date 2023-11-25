module atom_line_mod
!-
!   Contains the information on an atom line in the input structure file
!+
!
use precision_mod
!
implicit none
!
private
!
public  at_number
public  at_style
public  at_vals
public  AT_COMMA
public  AT_XYZ
public  AT_XYZB
public  AT_XYZBP
public  AT_XYZBPMMOS
public  AT_XYZBPMMOSU
public  atom_line_inter
public  atom_get_size
public  atom_line_get_style
public  read_atom_line
public  atom_alloc
public  atom_dealloc
public   get_atom_werte
public  atom_verify_param
public  atom_get_type
!
integer, parameter :: AT_MAXP = 22
!
character(len=8), dimension(22), PARAMETER :: at_cnam = (/ &
               'X       ', 'Y       ', 'Z       ', 'BISO    ', 'PROPERTY',  &
               'MOLENO  ', 'MOLEAT  ', 'OCC     ', 'ST      ', 'SH      ',  &
               'SK      ', 'SL      ',                                      &
               'U11     ', 'U22     ', 'U33     ', 'U23     ', 'U13     ',  &
               'U12     ',                                                  &
               'MM      ', 'MV      ', 'MU      ', 'MW      '               &
             /)
character(len=8), dimension(22)      :: at_param
integer         , dimension(AT_MAXP) :: at_look  ! at_param(i) uses at_cnam(at_look(i))
integer         , dimension(AT_MAXP) :: at_kool  ! at_cnam(j) is used by at_param(at_kool(j))
logical         , dimension(0:AT_MAXP) :: at_user  ! User provided parameters on atom line
logical                              :: at_number ! Atom lines contain numbers only
integer                              :: at_style ! Aktual style number
integer                              :: at_vals  ! Aktual values present 
integer                              :: at_ianz
integer                              :: at_natoms ! Current number of atoms
integer, parameter :: AT_COMMA     =  1          ! Comma delimited style, full flexibility
integer, parameter :: AT_XYZ       =  2          ! Just coordinates as pure numbers, no comma
integer, parameter :: AT_XYZB      =  3          ! Just coordinates, Biso   numbers, no comma
integer, parameter :: AT_XYZBP     =  4          ! Just coordinates, Biso, Property, no comma
integer, parameter :: AT_XYZBPMMOS =  5          ! All the way to Surface in DISCUS sequence
integer, parameter :: AT_XYZBPMMOSU=  6          ! All the way to Uij   e in DISCUS sequence
!
integer, parameter :: COL_X        =  1
integer, parameter :: COL_Y        =  2
integer, parameter :: COL_Z        =  3
integer, parameter :: COL_BISO     =  4
integer, parameter :: COL_PROPERTY =  5
integer, parameter :: COL_MOLENO   =  6
integer, parameter :: COL_MOLEAT   =  7
integer, parameter :: COL_OCC      =  8
integer, parameter :: COL_SURFT    =  9
integer, parameter :: COL_SURFH    = 10
integer, parameter :: COL_SURFK    = 11
integer, parameter :: COL_SURFL    = 12
integer, parameter :: COL_U11      = 13
integer, parameter :: COL_U22      = 14
integer, parameter :: COL_U33      = 15
integer, parameter :: COL_U23      = 16
integer, parameter :: COL_U13      = 17
integer, parameter :: COL_U12      = 18
integer, parameter :: COL_MM       = 19
integer, parameter :: COL_MU       = 20
integer, parameter :: COL_MV       = 21
integer, parameter :: COL_MW       = 22
!
character(len=4)  , dimension(:  ), allocatable :: at_names   ! All atom names
real(kind=PREC_DP), dimension(:,:), allocatable :: at_values  ! All atom coordinates etc
!
contains
!
!*******************************************************************************
!
!-
!  Subroutines to handle the 'atom' line, and to interpret the actual atom lines
!+
!
!*******************************************************************************
!
subroutine atom_line_inter(lline, llength)
!-
!  Determine the style of the 'atom' line, determine location of parameters etc.
!+
!
use get_params_mod
use precision_mod
use string_convert_mod
!
implicit none
!
character(len=*), intent(inout) :: lline      ! The atom line
integer         , intent(inout) :: llength    ! Its length
!
character(len=len(lline)) :: line
integer                   :: length
character(len=PREC_STRING), dimension(AT_MAXP) :: cpara
integer                   , dimension(AT_MAXP) :: lpara
integer :: ianz
integer :: i, j
integer :: inew
!
line = lline
length = llength
call get_params (line, ianz, cpara, lpara, AT_MAXP, length)
at_param = ' '              ! Ensure pristine state
at_ianz  = 0
at_look  = 0
at_user = .false.
!
if(ianz==0) then ! Pre 5.17.2 style, no params
   at_style = -1  ! Undetermined old style
   at_ianz = 4    ! At least x,y,z,Biso
   at_param(1) = 'X'
   at_param(2) = 'Y'
   at_param(3) = 'Z'
   at_param(4) = 'BISO'
   at_param(5:) = ' '
   do i=1, AT_MAXP
      at_look(i) = i
      at_kool(i) = i
   enddo
   at_user(1:4) = .TRUE.
   at_user(5: ) = .FALSE.
else
   at_style = AT_COMMA
   do i=1, ianz
      call do_cap(cpara(i))
      at_param(i) = cpara(i)(1:MIN(LEN(at_param),lpara(i)))
   enddo
   at_ianz = ianz
   inew    = 0
   loop_used: do i=1, at_ianz
      loop_names: do j=1,AT_MAXP
         if(at_param(i)==at_cnam(j)) then
            at_look(i) = j
            at_kool(j) = i
            at_user(j) = .TRUE.   ! User provide atom parameter J
            cycle loop_used
         endif
      enddo loop_names
   enddo loop_used
!         inew = inew + 1
!         at_look(inew)  = j
!         at_kool(j)     = inew
!         at_param(inew) = at_cnam(j)
!         at_user(j) = .FALSE.   ! User did not provide atom parameter J
  if(    line=='     x,              y,              z,             Biso,    Property,  MoleNo,  MoleAt,   Occ,     St,  Sh,  Sk,  Sl,  U11,       U22,       U33,       U23,       U13,       U12') then
    at_vals = AT_XYZBPMMOSU
  elseif(line=='     x,              y,              z,             Biso,    Property,  MoleNo,  MoleAt,   Occ,     St,  Sh,  Sk,  Sl') then
    at_vals = AT_XYZBPMMOS
  elseif(line=='     x,              y,              z,             Biso,    Property') then
    at_vals = AT_XYZBP
  elseif(line=='     x,              y,              z,             Biso'             ) then
    at_vals = AT_XYZB
  elseif(line=='     x,              y,              z'                               ) then
    at_vals = AT_XYZ
  endif
endif
at_natoms = 0
!
end subroutine atom_line_inter
!
!*******************************************************************************
!
subroutine atom_get_size(infile, nlines)
!-
! quickly determine number of lines in input file
!+
!
use envir_mod
use errlist_mod
use precision_mod
!
implicit none
!
character(len=*), intent(in)  :: infile
integer         , intent(out) :: nlines
!
integer, parameter :: IRD = 63
character(len=PREC_STRING) :: line
!character(len=PREC_STRING) :: tfile
!integer                    :: tfile_l
integer :: ios       ! I/O status flag
!
open(unit=ird, file=infile, status='old', iostat=ios)
if(ios/=0) then
   ier_num = -2
   ier_typ =ER_IO
   ier_msg(1) = 'Cannot open structure file '
   ier_msg(2) = infile(1:min(40,len_trim(infile)))
   return
endif
ios = 0
nlines = 0
loop_main: do
  read(ird,'(a)', iostat=ios) line
  if(is_iostat_end(ios)) exit loop_main
  nlines = nlines + 1
enddo loop_main
close(unit=ird)
return
!
!tfile = tmp_dir(1:tmp_dir_l) // '/discus_atom.lines'
!tfile_l = tmp_dir_l + 17
!!
!line = 'cat ' // infile(1:len_trim(infile)) // ' | wc -l > ' // tfile(1:tfile_l)
!call execute_command_line(line)
!!
!open(IRD, file=tfile(1:tfile_l), status='old')
!read(IRD, *, iostat=ios) nlines
!close(IRD)
!!
!line = ' rm -f ' // tfile(1:tfile_l)
!call execute_command_line(line)
!!
!at_natoms = 0
!
end subroutine atom_get_size
!
!*******************************************************************************
!
subroutine atom_line_get_style(line, ibl, length, MAXW, werte)
!-
!  Determine the style of the first line in the input file, only called if the
!  'atom' line is empty
!+
!
use errlist_mod
use precision_mod
!
implicit none
!
character(len=*)                     , intent(inout) :: line      ! The atom line
integer                              , intent(in)    :: ibl       ! Position of blank space
integer                              , intent(in)    :: length    ! Its length
integer                              , INTENT(in)    :: MAXW
real(kind=PREC_DP), dimension(1:MAXW), INTENT(out)   :: werte
!
integer :: j      ! Dummy index
integer :: ios       ! I/O status flag
!
at_user = .false.
read(line(ibl:length), *, iostat=ios) (werte(j), j = 1, 5)   !Try to read 5 params
if(.not.is_iostat_end(ios)) then        ! five paramters
   at_vals = AT_XYZBP
   at_user(1:5) = .true.
   if(index(line(ibl:length),',')>0) at_style=AT_COMMA
else
   read(line(ibl:length), *, iostat=ios) (werte(j), j = 1, 4)   !Try to read 4 params
   if(.not.is_iostat_end(ios)) then        ! four paramters
      at_vals = AT_XYZB
      at_user(1:4) = .true.
      if(index(line(ibl:length),',')>0) at_style=AT_COMMA
   else
      read(line(ibl:length), *, iostat=ios) (werte(j), j = 1, 3)   !Try to read 3 params
      if(.not.is_iostat_end(ios)) then        ! four paramters
         at_vals = AT_XYZ
         at_user(1:3) = .true.
         if(index(line(ibl:length),',')>0) at_style=AT_COMMA
      else
         ier_num = -49
         ier_typ = ER_APPL
         ier_msg(1) = 'Error reading first atom line '
      endif
   endif
endif
!
end subroutine atom_line_get_style
!
!*******************************************************************************
!
subroutine atom_alloc(nlines)
!-
! Allocate array atom_values
!+
!
implicit none
!
integer, intent(in) :: nlines
!
if(allocated(at_values)) deallocate(at_values)
if(allocated(at_names )) deallocate(at_names )
allocate(at_values(AT_MAXP, nlines))
allocate(at_names (         nlines))
!
end subroutine atom_alloc
!
!*******************************************************************************
!
subroutine atom_dealloc
!-
! Deallocate array atom_values
!+
!
implicit none
!
if(allocated(at_values)) deallocate(at_values)
if(allocated(at_names )) deallocate(at_names )
!
end subroutine atom_dealloc
!
!*******************************************************************************
!
subroutine read_atom_line(lline, ibl,llength, cr_natoms, MAXW, werte)
!-
!  reads a line from the cell file/structure file                    
!+
!
use errlist_mod
use ber_params_mod
use get_params_mod
use precision_mod
!
implicit none
!
character(len=*)                     , intent(inout) ::lline      ! The atom line
character(len=PREC_STRING)                           :: line      ! The atom line
integer                              , intent(in)    :: ibl       ! Its length
integer                              , intent(inout) ::llength    ! Its length
integer                                              :: length    ! Its length
integer                              , intent(in)    :: cr_natoms
integer                              , INTENT(in)    :: MAXW
real(kind=PREC_DP), dimension(1:MAXW), INTENT(out)   :: werte
!
character(len=12) :: cform
character(len=1 ) :: csurf
character(len=max(PREC_STRING,len(line))), dimension(MAXW) :: cpara
integer                                  , dimension(MAXW) :: lpara
character(len=max(PREC_STRING,len(line))), dimension(1   ) :: ccpara
integer                                  , dimension(1   ) :: llpara
integer                                  , dimension(6   ) :: iwerte
integer :: laenge ! length of significant string
integer :: ianz   ! Number of comma delimited values in the line
integer :: iianz   ! Number of comma delimited values in the line
integer :: ios    ! I/O status flag
integer :: i      ! Dummy index
integer :: j      ! Dummy index
!real(kind=PREC_DP), dimension(1:MAXW) :: wwerte
real(kind=PREC_DP), dimension(1:1   ) :: wwerte
!
line = lline
length = llength
werte    = 0.0D0          ! Initialize values
werte(5) = 1.0D0
werte(8) = 1.0D0
!
if_number: if(at_number) then                                ! Numbers only in DISCUS format
   if(at_vals==AT_XYZBPMMOSU) then                           ! Line with X,y,Z,B,Prop, Mole, mole, Occ, Surface
      read(line(5:length), 44445, &
      iostat=ios) werte(1:4), iwerte(1:3), werte(8), csurf, iwerte(4:6), werte(13:18)
44445 FORMAT(3(1x,f14.6,1x ),4x,f10.6,1x ,i8, 1x , I8, 1x , I8,2x  , F10.6,2x  ,A1,3(2x  ,I3),6(1x,f10.6))
      if(ios/=0) then
         ier_num = -49 ! 
         ier_typ = ER_APPL
         ier_msg(1) = 'File format specified as numbers only'
         ier_msg(2) = 'Check atom line starting with '
         ier_msg(3)(1:40) = line(1:40)
         return
      endif
      werte( 5: 7) = real(iwerte(1:3),kind=PREC_DP)
      select case(csurf)
         case('_')
            werte(9) = 0.0D0
         case('P')
            werte(9) = 1.0D0
         case('S')
            werte(9) = 2.0D0
         case('Y')
            werte(9) = 3.0D0
         case('E')
            werte(9) = 4.0D0
         case('C')
            werte(9) = 5.0D0
         case('L')
            werte(9) = 6.0D0
         case('T')
            werte(9) = 7.0D0
         case default
      end select
      werte(10:12) = real(iwerte(4:6),kind=PREC_DP)
   elseif(at_vals==AT_XYZBPMMOS) then                           ! Line with X,y,Z,B,Prop, Mole, mole, Occ, Surface
      read(line(5:length), 44444, &
      iostat=ios) werte(1:4), iwerte(1:3), werte(8), csurf, iwerte(4:6)
44444 FORMAT(3(1x,f14.6,1x ),4x,f10.6,1x ,i8, 1x , I8, 1x , I8,2x  , F10.6,2x  ,A1,3(2x  ,I3))
      if(ios/=0) then
         ier_num = -49 ! 
         ier_typ = ER_APPL
         ier_msg(1) = 'File format specified as numbers only'
         ier_msg(2) = 'Check atom line starting with '
         ier_msg(3)(1:40) = line(1:40)
         return
      endif
      werte( 5: 7) = real(iwerte(1:3),kind=PREC_DP)
      select case(csurf)
         case('_')
            werte(9) = 0.0D0
         case('P')
            werte(9) = 1.0D0
         case('S')
            werte(9) = 2.0D0
         case('Y')
            werte(9) = 3.0D0
         case('E')
            werte(9) = 4.0D0
         case('C')
            werte(9) = 5.0D0
         case('L')
            werte(9) = 6.0D0
         case('T')
            werte(9) = 7.0D0
         case default
      end select
      werte(10:12) = real(iwerte(4:6),kind=PREC_DP)
   endif
else if_number                                               ! Might contain variables
!
if_style: if(at_style==AT_COMMA) then             ! x,y,z,B, ...
   cpara((COL_PROPERTY)) = '1.0D0'
   cpara((COL_OCC ))     = '1.0D0'
!
   laenge = length - ibl + 1
   call get_params(line(ibl:length), ianz, cpara, lpara, MAXW, laenge)
!
!  The line has comma separated parameters, compare to expectation from 'atom' line
   got_params: if(ier_num == 0) then
      if(at_user((COL_SURFT))) then  ! User specfied surface type
         select case(cpara((COL_SURFT)))
         case('_')
            cpara((COL_SURFT)) = '0.0D0'
         case('P')
            cpara((COL_SURFT)) = '1.0D0'
         case('S')
            cpara((COL_SURFT)) = '2.0D0'
         case('Y')
            cpara((COL_SURFT)) = '3.0D0'
         case('E')
            cpara((COL_SURFT)) = '4.0D0'
         case('C')
            cpara((COL_SURFT)) = '5.0D0'
         case('L')
            cpara((COL_SURFT)) = '6.0D0'
         case('T')
            cpara((COL_SURFT)) = '7.0D0'
         case default
            cpara((COL_SURFT)) = '0.0D0'
         end select
      endif
      iianz     = 1
      do i=1, AT_MAXP
         if(at_user(at_look(i))) then      ! User specified this on 'atom x,y,z, ...' line
            if(index(cpara(i)(1:lpara(i)),'.')==0) then
               cform = '(*)'
            else
               write(cform,'(a2,i4.4,a1,i4.4,a1)') '(f',lpara(i), '.',lpara(i)-index(cpara(i)(1:lpara(i)),'.'),')'
            endif
            read(cpara((i)),cform, iostat=ios) werte(at_look(i))
            if(ios/=0) then
               ccpara(1) = cpara((i))
               llpara(1) = lpara((i))
               call ber_params(iianz, ccpara, llpara, wwerte, iianz)
               if(ier_num/=0) exit if_style
               werte(at_look(i)) = wwerte(1)
            endif
         endif
      enddo
   else
      exit if_style
   endif got_params
else                             if_style
    if(at_vals==AT_XYZB) then              ! x y z B no comma
   read(line(ibl:length), *, iostat=ios) (werte(j), j = 1, 4)
   elseif(at_vals==AT_XYZ) then!   if_style         ! x y z   no comma
   read(line(ibl:length), *, iostat=ios) (werte(j), j = 1, 3)
   elseif(at_vals==AT_XYZBP) then! if_style         ! x y z B P   no comma
   read(line(ibl:length), *, iostat=ios) (werte(j), j = 1, 5)
   endif
!else                             if_style
endif   if_style
endif if_number
at_natoms = at_natoms + 1
at_names(at_natoms) = line(1:ibl-1)
at_values(:,at_natoms) = werte
!write(*,*) ' AT_LINE ', line(1:len_trim(line))
!write(*,'(a,19f8.3)') ' Werte ', werte
!
if(ier_num/=0) then
  ier_msg (1) = 'Error reading parameters for'
  ier_msg (2) = 'coordinates for atom '//line (1:ibl)
  write(ier_msg(3), '(''Atom Nr. '',i4)') cr_natoms + 1
endif
!
end subroutine read_atom_line
!
!*******************************************************************************
!
subroutine get_atom_werte(natom, MAXP, werte)
!
use precision_mod
use errlist_mod
!
implicit none
!
integer, intent(in) :: natom
integer, intent(in) :: MAXP
real(kind=PREC_DP), dimension(MAXP), intent(out) :: werte
!
werte(1:AT_MAXP) = at_values(1:MAXP,natom)
ier_num = 0
ier_typ = 0
!
end subroutine get_atom_werte
!
!*******************************************************************************
!
subroutine atom_verify_param
!-
! Check that all user parameters on the "atoms" line are valid names
!+
!
use errlist_mod
!
implicit none
!
integer :: i   ! Loop index
integer :: j   ! Loop index
!
do_param: do i= 1, at_ianz
   do_allowed: do j=1, AT_MAXP
      if(at_param(i)==at_cnam(j)) cycle do_param
   enddo do_allowed
   ier_num = -189
   ier_typ = ER_APPL
   ier_msg(3) = 'Illegal : ' // at_param(i)(1:len_trim(at_param(i)))
   return
enddo do_param
end subroutine atom_verify_param
!
!*******************************************************************************
!
function atom_get_type(MAXSCAT, nscat_start, nscat, MAXMASK,      &
                       at_lis, at_dw, at_occ,        &
                       nw_name, nw_dw, nw_occ, mask) &
         result(is_type)
!-
!  Determine atom type from name, charge, Bval, (occ) 
!+
!
use precision_mod
use string_convert_mod
!
implicit none
!
integer, intent(in) :: MAXSCAT
integer, intent(in) :: nscat_start
integer, intent(in) :: nscat
integer, intent(in) :: MAXMASK
character(len=4)  , dimension(nscat_start:MAXSCAT), intent(in) :: at_lis
real(kind=PREC_DP), dimension(nscat_start:MAXSCAT), intent(in) :: at_dw
real(kind=PREC_DP), dimension(nscat_start:MAXSCAT), intent(in) :: at_occ
character(len=*)  ,                       intent(inout) :: nw_name
real(kind=PREC_DP),                       intent(in) :: nw_dw
real(kind=PREC_DP),                       intent(in) :: nw_occ
logical           , dimension(0:MAXMASK), intent(in) :: mask
!
real(kind=PREC_DP), parameter :: TOL = 1.0D-5
integer :: is_type
!
character(len=4) :: local_name
character(len=4) :: local_full
integer :: i   ! Loop index
integer :: j   ! Loop index
!
is_type = -1
!
if(.not.any(mask(1:))) return       ! No masking at all; all are new atom types
!
j = len_trim(nw_name)
if(nw_name(j:j)=='+' .or. nw_name(j:j)=='-')  j = j - 2 ! Shorten name
local_name = nw_name(1:j)
local_full = nw_name
call do_cap(local_name)       ! Capitalized, pure atom name
call do_cap(local_full)       ! Capitalized, full atom name
!
do_search: do i=nscat_start, nscat
   if(mask(4)) then           ! name, charge, Bval and occ are required
      if(at_lis(i)==local_full .and. abs(at_dw(i)-nw_dw)<TOL .and. abs(at_occ(i)-nw_occ)<TOL) then
         is_type = i
         exit do_search
      else
         cycle do_search
      endif
   elseif(mask(3)) then       ! name, charge, Bval are required
      if(at_lis(i)==local_full .and. abs(at_dw(i)-nw_dw)<TOL) then
         is_type = i
         exit do_search
      else
         cycle do_search
      endif
   elseif(mask(2)) then       ! name, charge
      if(at_lis(i)==local_full) then
         is_type = i
         exit do_search
      else
         cycle do_search
      endif
   elseif(mask(1)) then       ! name
      if(at_lis(i)(1:j)==local_name(1:j)) then
         is_type = i
         nw_name = local_name
         exit do_search
      else
         cycle do_search
      endif
   else
      is_type = i
      exit do_search
   endif
enddo do_search
if(is_type==-1 .and. .not.mask(2)) then   ! New atom type and ignore charge
   nw_name = local_name
endif
!
end function atom_get_type
!
!*******************************************************************************
!
end module atom_line_mod
