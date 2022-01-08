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
public  at_style
public  atom_line_inter
public  atom_get_size
public  atom_line_get_style
public  read_atom_line
public  atom_alloc
public  atom_dealloc
public   get_atom_werte
!
integer, parameter :: AT_MAXP = 16
!
character(len=8), dimension(16), PARAMETER :: at_cnam = (/ &
               'X       ', 'Y       ', 'Z       ', 'BISO    ', 'PROPERTY',  &
               'MOLENO  ', 'MOLEAT  ', 'OCC     ', 'ST      ', 'SH      ',  &
               'SK      ', 'SL      ', 'MM      ', 'MV      ', 'MU      ',  &
               'MW      '                                                   &
             /)
character(len=8), dimension(16)      :: at_param
integer         , dimension(AT_MAXP) :: at_look  ! at_param(i) uses at_cnam(at_look(i))
integer         , dimension(AT_MAXP) :: at_kool  ! at_cnam(j) is used by at_param(at_kool(j))
logical         , dimension(0:AT_MAXP) :: at_user  ! User provided parameters on atom line
integer                              :: at_style ! Aktual style number
integer                              :: at_ianz
integer                              :: at_natoms ! Current number of atoms
integer, parameter :: AT_COMMA     =  1          ! Comma delimited style, full flexibility
integer, parameter :: AT_XYZ       =  2          ! Just coordinates as pure numbers, no comma
integer, parameter :: AT_XYZB      =  3          ! Just coordinates, Biso   numbers, no comma
integer, parameter :: AT_XYZBP     =  4          ! Just coordinates, Biso, Property, no comma
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
integer, parameter :: COL_MM       = 13
integer, parameter :: COL_MU       = 14
integer, parameter :: COL_MV       = 15
integer, parameter :: COL_MW       = 16
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
   at_user(4: ) = .FALSE.
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
!write(*,*) ' New at ', inew
!         at_look(inew)  = j
!         at_kool(j)     = inew
!         at_param(inew) = at_cnam(j)
!         at_user(j) = .FALSE.   ! User did not provide atom parameter J
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
use precision_mod
!
implicit none
!
character(len=*), intent(in)  :: infile
integer         , intent(out) :: nlines
!
integer, parameter :: IRD = 63
character(len=PREC_STRING) :: line
character(len=PREC_STRING) :: tfile
integer                    :: tfile_l
integer :: ios       ! I/O status flag
!
tfile = tmp_dir(1:tmp_dir_l) // '/discus_atom.lines'
tfile_l = tmp_dir_l + 17
!
line = 'cat ' // infile(1:len_trim(infile)) // ' | wc -l > ' // tfile(1:tfile_l)
call execute_command_line(line)
!
open(IRD, file=tfile(1:tfile_l), status='old')
read(IRD, *, iostat=ios) nlines
close(IRD)
!
line = ' rm -f ' // tfile(1:tfile_l)
call execute_command_line(line)
!
at_natoms = 0
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
read(line(ibl:length), *, iostat=ios) (werte(j), j = 1, 5)   !Try to read 5 params
if(.not.is_iostat_end(ios)) then        ! five paramters
   at_style = AT_XYZBP
else
   read(line(ibl:length), *, iostat=ios) (werte(j), j = 1, 4)   !Try to read 4 params
   if(.not.is_iostat_end(ios)) then        ! four paramters
      at_style = AT_XYZB
   else
      read(line(ibl:length), *, iostat=ios) (werte(j), j = 1, 3)   !Try to read 3 params
      if(.not.is_iostat_end(ios)) then        ! four paramters
         at_style = AT_XYZ
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
character(len=max(PREC_STRING,len(line))), dimension(MAXW) :: cpara
integer                                  , dimension(MAXW) :: lpara
character(len=max(PREC_STRING,len(line))), dimension(1   ) :: ccpara
integer                                  , dimension(1   ) :: llpara
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
            read(cpara((i)),*, iostat=ios) werte(at_look(i))
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
elseif(at_style==AT_XYZB) then              ! x y z B no comma
   read(line(ibl:length), *, iostat=ios) (werte(j), j = 1, 4)
elseif(at_style==AT_XYZ) then    if_style         ! x y z   no comma
   read(line(ibl:length), *, iostat=ios) (werte(j), j = 1, 3)
elseif(at_style==AT_XYZBP) then  if_style         ! x y z B P   no comma
   read(line(ibl:length), *, iostat=ios) (werte(j), j = 1, 5)
else                             if_style
endif   if_style
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
end module atom_line_mod
