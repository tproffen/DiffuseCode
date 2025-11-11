module perioditize_mod
!-
!  Attempt to rearrange the structure as periodic structure in the standard
!  DISCUS sequence
!
use precision_mod
!
private
public perioditize_menu
!public estimate_ncells
!public estimate_ncatom
!
integer :: pdt_ncells     ! Number of cells in periodic crystal volume
integer :: pdt_nsite = 0  ! Number of sites in an averaged unit cell
integer :: pdt_usite = 0  ! Number of sites in an averaged unit cell, defined by user with 'set site'
integer :: pdt_asite = 0  ! Achieved sites either from set site or find_aver
!integer           , dimension(:,:,:), allocatable :: pdt_cell     ! Atoms per cell
logical           , dimension(:    ), allocatable :: pdt_lsite    ! This site has been defined
integer           , dimension(:,  :), allocatable :: pdt_itype    ! Atom types at each site
integer           , dimension(3)                  :: pdt_ilow = 0    ! Unit cell dimensions in periodic
integer           , dimension(3)                  :: pdt_ihig = 0    ! low and high inidces
logical                                           :: pdt_usr_ncell = .FALSE.  ! True if user fixed number of cells
logical                                           :: pdt_usr_nsite = .FALSE.  ! True if user fixed number of sites
logical                                           :: pdt_usr_atom  = .FALSE.  ! True if user fixed atoms at sites
real(kind=PREC_DP), dimension(3, 2)               :: pdt_dims ! Crystal dimensions
real(kind=PREC_DP), dimension(:,:)  , allocatable :: pdt_pos  ! Atom positions in average cell
!
contains
!
!*******************************************************************************
!
subroutine perioditize_menu
!-
!  Main menu routine for perioditizing
!+
!
use discus_kdo_common_mod
!
use calc_expr_mod
use chem_mod
use class_macro_internal
use doact_mod
use errlist_mod
use lib_length
use lib_errlist_func
use lib_macro_func
use macro_mod
!use macro_internal_mod
use precision_mod
use prompt_mod
use str_comp_mod
use sup_mod
!
implicit none
!
character(len=PREC_STRING) :: zeile        ! Input line
character(len=PREC_STRING) :: line         ! Input line
character(len=6          ) :: befehl       ! The actual command
character(len=len(prompt)) :: orig_prompt
!
integer :: indxg
integer :: length            ! Input line length
integer :: lbef              ! Command length
integer :: lp                ! Line length
logical :: lend              ! if true terminate the menu
logical :: success
!
call no_error
if(.not.allocated(pdt_pos))   then
   allocate(pdt_pos  (1:3, 1:1))
   pdt_pos = 0.0
endif
if(.not.allocated(pdt_itype)) then
   allocate(pdt_itype( 0:1, 1:1))
   pdt_itype(0,1) =  1
   pdt_itype(1,1) = -1
endif
!
orig_prompt = prompt
prompt = prompt (1:len_str (prompt) ) //'/period'
lend   = .FALSE.
!
loop_menu: do
   if(ier_num == 0) then
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt)  ! get new command
   endif
   if(lend) exit loop_menu
!
   IF(ier_num /= 0) THEN   ! Error in previous command or in get_cmd
      CALL errlist
      IF (ier_sta.ne.ER_S_LIVE) THEN
         IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
            IF(sprompt /= prompt ) THEN
               ier_num = -10
               ier_typ = ER_COMM
               ier_msg(1) = ' Error occured in perioditize menu'
               prompt_status = PROMPT_ON 
               prompt = orig_prompt
               RETURN
            ELSE
               IF(lmacro_close) THEN
                  CALL macro_close(-1)
                  prompt_status = PROMPT_ON 
               ENDIF 
            ENDIF 
         ENDIF 
         IF (lblock) THEN 
            ier_num = - 11 
            ier_typ = ER_COMM 
            prompt_status = PROMPT_ON 
            prompt = orig_prompt
            RETURN 
         ENDIF 
         CALL no_error 
         lmakro_error = .FALSE.
         sprompt = ' '
      ENDIF 
      cycle loop_menu
   ENDIF 
!
!  Normal menu commands
!
   IF (line == ' '      .or. line(1:1) == '#' .or. &
       line == char(13) .or. line(1:1) == '!'        ) cycle loop_menu
!                                                                       
!     ----search for "="                                                
!                                                                       
   indxg = index (line, '=') 
   IF (indxg.ne.0.AND..NOT. (str_comp (befehl, 'echo',   2, lbef, 4) ) &
                 .AND..NOT. (str_comp (befehl, 'system', 2, lbef, 6) )    &
                 .AND..NOT. (str_comp (befehl, 'help',   2, lbef, 4) .OR. &
                             str_comp (befehl, '?   ',   2, lbef, 4) )    &
                 .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
!     ---evaluate an expression and assign the value to a variabble   
!                                                                       
      CALL do_math (line, indxg, length) 
      cycle loop_menu
   endif
!
! - Test common menu entries
!
   CALL discus_kdo_common(befehl, lbef, line, length, zeile, lp, 'peri' , &
                          lend, success)
   if(lend) exit loop_menu
   if(success) cycle loop_menu
!                                                                       
!     ----run symmetry 'rese'                                            
!                                                                       
   if_rese: IF(str_comp (befehl, 'reset', 2, lbef, 5) ) THEN 
      call perioditize_reset
      cycle loop_menu
   endif if_rese
!-------------------------------------------------------------------------------
!-------- Perioditize specific commands if no common command was found
!-------------------------------------------------------------------------------
!
!  IF( cr_nscat > SYM_MAXSCAT .or. mole_num_type > SYM_MAXSCAT) THEN
!     nscat = max ( cr_nscat, mole_num_type)
!     nsite = max ( nsite, cr_ncatoms, SYM_MAXSITE)
!     CALL alloc_symmetry ( nscat, nsite )
!     IF ( ier_num < 0 ) THEN
!        RETURN
!     ENDIF
!  ENDIF
!                                                                       
!     ----exit menu        'exit'                                     
!                                                                       
   if(str_comp (befehl, 'run', 2, lbef, 3) ) then
      call perioditize_run
      call perioditize_reset
   elseif(str_comp (befehl, 'show', 2, lbef, 4) ) then
      call perioditize_show
   elseif(str_comp (befehl, 'set', 2, lbef, 3) ) then
      call perioditize_set(zeile, lp)
   else
      ier_num = -8
      ier_typ = ER_COMM
   endif
   if(lend) exit loop_menu
enddo loop_menu
!
prompt = orig_prompt
!
end subroutine perioditize_menu
!
!*******************************************************************************
!
subroutine perioditize_set(zeile, lp)
!-
!  Set user defined parameters
!+
use crystal_mod
use errlist_mod
use get_params_mod
use precision_mod
use str_comp_mod
!
implicit none
!
character(len=*), intent(inout) :: zeile      ! input line
integer         , intent(inout) :: lp         ! input length
!
integer            :: maxw                        ! Actual parameter number
integer, parameter :: MIN_PARA =  4               ! Minimum parameter numbers
character(len=PREC_STRING), dimension(MAX(MIN_PARA,MAXSCAT+1)) :: cpara
integer                   , dimension(MAX(MIN_PARA,MAXSCAT+1)) :: lpara
real(kind=PREC_DP)        , dimension(MAX(MIN_PARA,MAXSCAT+1)) :: werte
integer :: ianz          ! Number of parameters
!
!
!
MAXW = MAX(MIN_PARA,MAXSCAT+1)
!
call get_params (zeile, ianz, cpara, lpara, maxw, lp)
if(ier_num/=0) return
!
!call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  &
!                  ncalc, oname, loname, opara, lopara, lpresent, owerte)
!
if(str_comp (cpara(1), 'ncells', 2, lpara(1), 6) ) then
   call perioditize_set_ncells(ianz, cpara, lpara, werte, maxw)
elseif(str_comp (cpara(1), 'natoms', 2, lpara(1), 6) ) then
   call perioditize_set_natoms(ianz, cpara, lpara, werte, maxw)
elseif(str_comp (cpara(1), 'site', 2, lpara(2), 4) ) then
   call perioditize_set_site  (ianz, cpara, lpara, maxw, .false.)
elseif(str_comp (cpara(1), 'wyckoff', 2, lpara(2), 7) ) then
   call perioditize_set_site  (ianz, cpara, lpara, maxw, .true.)
endif
!
end subroutine perioditize_set
!
!*******************************************************************************
!
subroutine perioditize_set_ncells(ianz, cpara, lpara, werte, maxw)
!-
!  Set the number of atoms per unit cell
!+
!
use ber_params_mod
use errlist_mod
use precision_mod
!
implicit none
!
integer :: MAXW          ! max number of parameters
integer :: ianz          ! Number of parameters
integer :: i             ! Dummy loop parameter
character(len=PREC_STRING), dimension(MAXW), intent(inout)     :: cpara
integer                   , dimension(MAXW), intent(inout)     :: lpara
real(kind=PREC_DP)        , dimension(MAXW), intent(inout)     :: werte
!
cpara(1) = '0'
lpara(1) = 1
!
if(cpara(2) == 'auto') then
   pdt_ilow = 1
   pdt_ihig = 1
   pdt_usr_ncell = .FALSE.                 ! Determine sites automatically
else
   cpara(1) = '0'
   lpara(1) = 1
   call ber_params(ianz, cpara, lpara, werte, maxw)
   if(ier_num/=0) return
   pdt_ncells = 1
   do i=1, 3
      if(mod(nint(werte(i+1)),2)==0) then    ! Even number of unit cells
         pdt_ilow(i) = -nint(werte(i+1)/2.)
         pdt_ihig(i) = pdt_ilow(i) + nint(werte(i+1)) - 1
      else
         pdt_ilow(i) = -nint((werte(i+1)-1.)/2.)
         pdt_ihig(i) = pdt_ilow(i) + nint(werte(i+1)) - 1
      endif
      pdt_ncells = pdt_ncells * nint(werte(i+1))
   enddo
   pdt_usr_ncell = .TRUE.                  ! User fixed number of sites
endif
!write(*,'(a)') ' User defined unit cell numbers'
!do i=1, 3
!  write(*,'(a,2i5)') ' Unit cells ', pdt_ilow(i), pdt_ihig(i)
!enddo
!
end subroutine perioditize_set_ncells
!
!*******************************************************************************
!
subroutine perioditize_set_natoms(ianz, cpara, lpara, werte, maxw)
!-
!  Set the number of atoms per unit cell
!+
!
use allocate_generic
use ber_params_mod
use errlist_mod
use precision_mod
!
implicit none
!
integer :: MAXW          ! max number of parameters
integer :: ianz          ! Number of parameters
character(len=PREC_STRING), dimension(MAXW), intent(inout)     :: cpara
integer                   , dimension(MAXW), intent(inout)     :: lpara
real(kind=PREC_DP)        , dimension(MAXW), intent(inout)     :: werte
integer :: osite
integer :: nsite
integer :: all_status
!
cpara(1) = '0'
lpara(1) = 1
!
if(cpara(2) == 'auto') then
   pdt_nsite = 0
   pdt_usr_nsite = .FALSE.                 ! Determine sites automatically
else
   cpara(1) = '0'
   lpara(1) = 1
   call ber_params(ianz, cpara, lpara, werte, maxw)
   if(ier_num/=0) return
   pdt_nsite = nint(werte(2))
   pdt_usr_nsite = .TRUE.                  ! User fixed number of sites
endif
!write(*,'(a, i6, l2)') ' User set number sites ', pdt_nsite, pdt_usr_nsite
if(pdt_nsite>=ubound(pdt_pos,2)) then 
   osite = ubound(pdt_pos, 2)
   nsite = max(pdt_nsite, pdt_usite, pdt_asite) + 5
!write(*,*) 'PRE REALLOC ', ubound(pdt_pos, 2), osite, nsite, pdt_usite
   call alloc_arr(pdt_pos, 1, 3, 1, nsite, all_status, 0.0D0  )
   call alloc_arr(pdt_lsite, 1, nsite, all_status, .false.)
endif
end subroutine perioditize_set_natoms
!
!*******************************************************************************
!
subroutine perioditize_set_site(ianz, cpara, lpara, maxw, lwyckoff)
!-
!  Set explicit positions expected in the average unit cell
!+
!
use crystal_mod
use spcgr_apply, only:symmetry_make_wyckoff
!
use allocate_generic
use ber_params_mod
use berechne_mod
use errlist_mod
use get_iscat_mod
use get_params_mod
use precision_mod
use take_param_mod
!
implicit none
!
integer :: MAXW          ! max number of parameters
integer :: ianz          ! Number of parameters
character(len=PREC_STRING), dimension(MAXW), intent(inout)     :: cpara
integer                   , dimension(MAXW), intent(inout)     :: lpara
!real(kind=PREC_DP)        , dimension(MAXW), intent(inout)     :: werte
logical                                    , intent(in)        :: lwyckoff
!
!
integer, parameter :: MAXU= 3                   ! Max params for pos:[}
integer, parameter :: NNMAX= 192                ! Max Wyckoff positions
character(len=PREC_STRING) :: line
character(len=PREC_STRING), dimension(MAXSCAT) :: ccpara
integer                   , dimension(MAXSCAT) :: llpara
real(kind=PREC_DP)        , dimension(MAXSCAT) :: wwerte
real(kind=PREC_DP)        , dimension(MAXU)    :: uwerte
real(kind=PREC_DP)        , dimension(4   )    :: v_pos    ! Symmetry resulting position
real(kind=PREC_DP)        , dimension(3,192 )  :: wyck     ! Symmetry resulting position
integer :: iianz          ! Number of parameters
integer :: jjanz          ! Number of parameters
integer :: maxww          ! Number of parameters
integer :: next           ! next site number
integer :: length
integer :: osite
integer :: nsite
integer :: ntype 
integer :: all_status
integer :: istart
integer :: ifinish
logical :: success
integer :: i, j           ! Dummy loop variable
!
!
integer, parameter :: NOPTIONAL = 3

integer, parameter :: O_NUMBER  = 1
integer, parameter :: O_ATOM    = 2
integer, parameter :: O_POS     = 3
character(LEN=   6), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 0 ! Number of values to calculate 
!
data oname  / 'number', 'atoms ', 'pos'   /
data loname /  6      ,  5      ,  3      /
opara  =  (/ 'next   ', 'all    ', '[0,0,0]' /)   ! Always provide fresh default values
lopara =  (/  4       ,  3       , 7         /)
owerte =  (/  -1.0    ,  0.0     , 0.00      /)
!
MAXWW = MAXSCAT
!
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num/=0) return
cpara(1) = '0'
lpara(1) = 1
next = -1
!
if(cpara(2) == 'auto') then
   pdt_nsite     = 0                       ! No sites found automatically yet
   pdt_asite     = 0                       ! No sites found automatically yet
   pdt_usr_atom  = .FALSE.                 ! Determine sites automatically
   pdt_usite     = 0                       ! No sites fixed by user
   return
else
   if(lpresent(O_ATOM)) then               ! atom:[] is present
      if(opara(O_ATOM)=='all' ) then       ! Allow all atom ptypes at this site
         pdt_usr_atom = .TRUE.             ! User did flag atom positions
      else
         if(opara(O_ATOM)(1:1)=='[' .and. opara(O_ATOM)(lopara(O_ATOM):lopara(O_ATOM))==']') then
            line = opara(O_ATOM)(2:lopara(O_ATOM)-1)
            length = lopara(O_ATOM)-2
         else
            line = opara(O_ATOM)
            length = lopara(O_ATOM)
         endif
         call get_params(line, iianz, ccpara, llpara, maxwW, length)
         if(ier_num/=0) return
         call get_iscat (iianz, ccpara, llpara, wwerte, MAXWW, .FALSE. )   ! No new atom types
         if(ier_num/=0) return
      endif
   else
      ier_num = -6
      ier_typ = ER_COMM
      ier_msg(1) = '''atoms:'' parameter is missing'
      return
   endif
   if(lpresent(O_POS )) then               ! atom:[] is present
      line = opara(O_POS)
      length = lopara(O_POS)
      jjanz = 3
      call get_optional_multi(MAXU, line, length, uwerte, jjanz)
      if(ier_num/=0) return
      if(uwerte(1)<0.0 .or. uwerte(1)>1.0 .or.                                  &
         uwerte(2)<0.0 .or. uwerte(2)>1.0 .or.                                  &
         uwerte(3)<0.0 .or. uwerte(3)>1.0       ) then
         ier_num = -50
         ier_typ = ER_FORT
         ier_msg(1) = 'At least on of the coordinates outside [0,1]'
         return
      endif
   else
      ier_num = -6
      ier_typ = ER_COMM
      ier_msg(1) = '''pos:'' parameter is missing'
      return
   endif
   if(lpresent(O_NUMBER)) then               ! number: is present
      if(opara(O_NUMBER)=='next') then
         next = -1
      else
         line   = '(' // opara(O_NUMBER)(1:lopara(O_NUMBER)) // ')'
         length = lopara(O_NUMBER) + 2
         next = nint(berechne(line, length))
         if(ier_num/=0) return
         if(next<=0) then
            ier_num = -50
            ier_typ = ER_FORT
            ier_msg(1) = 'number:<value> must be > 1'
            return
         endif
      endif
   else
      next = -1
   endif
endif
!
pdt_usr_atom = .TRUE.                ! User did flag atom positions
!
if(lwyckoff) then                    ! Generate Wyckoff orbit
  wyck(:,1) = uwerte(1:3)
  call symmetry_make_wyckoff(NNMAX, wyck, ifinish)
  istart = 1
  next = -1
else                                 ! Just this site
  istart = 1
  ifinish = 1
  wyck(:,1) = uwerte(1:3)
endif
!
do j=istart, ifinish
!
   v_pos(1:3) = wyck(1:3,j)
!  if(j>1) then                      ! Need to apply symmetry operation
!     
!     v_pos = matmul(spc_mat(:,:,j),u_pos)
!     v_pos(1) = v_pos(1) - REAL(INT(v_pos(1)), PREC_DP)
!     if(v_pos(1)< 0.0) v_pos(1) = v_pos(1) + 1.0D0
!     v_pos(2) = v_pos(2) - REAL(INT(v_pos(2)), PREC_DP)
!     if(v_pos(2)< 0.0) v_pos(2) = v_pos(2) + 1.0D0
!     v_pos(3) = v_pos(3) - REAL(INT(v_pos(3)), PREC_DP)
!     if(v_pos(3)< 0.0) v_pos(3) = v_pos(3) + 1.0D0
!  endif
!
   
if(next == -1) then
   success = .false.
   search: do i=1, pdt_usite         ! Search for a free site
      if(.not.pdt_lsite(i)) then     ! this site has not been defined
         pdt_usite = i
         pdt_lsite(pdt_usite) = .true.
         success = .true.
         exit search
      endif
   enddo search
   if(.not.success) then             ! No free site found, make next
      pdt_usite = pdt_usite + 1
   endif
else
   pdt_usite = iabs(next)            ! User defined explicit site number
endif
pdt_asite = max(pdt_asite, pdt_usite)
if(pdt_usite>=ubound(pdt_pos,2)) then 
   osite = ubound(pdt_pos, 2)
   nsite = max(pdt_nsite, pdt_usite, pdt_asite) + 5
!write(*,*) 'PRE REALLOC ', ubound(pdt_pos, 2), osite, nsite, pdt_usite
   call alloc_arr(pdt_pos, 1, 3, 1, nsite, all_status, 0.0D0  )
   call alloc_arr(pdt_lsite, 1, nsite, all_status, .false.)
endif
if(pdt_usite>=ubound(pdt_itype,2) .or. iianz>ubound(pdt_itype,1)) then
   nsite = ubound(pdt_pos, 2)
   ntype = max(iianz, ubound(pdt_itype,1))
   call alloc_arr(pdt_itype, 0, ntype, 1, nsite, all_status, -2)
endif
!
pdt_pos(:,pdt_usite) = v_pos(1:3)
pdt_itype(1:iianz, pdt_usite) = nint(wwerte(1:iianz))
pdt_itype(0      , pdt_usite) = iianz
pdt_lsite(pdt_usite) = .true.
enddo
!
!write(*,*) 'SITE: ', pdt_usite, pdt_itype(0:iianz,pdt_usite), pdt_pos(1:jjanz,pdt_usite)
pdt_nsite = max(pdt_nsite, pdt_usite, pdt_asite)
!
end subroutine perioditize_set_site
!
!*******************************************************************************
!
subroutine perioditize_run
!-
! Perform the actual perioditization
!+
use crystal_mod
use dis_estimate_mod
!
use errlist_mod
use precision_mod
!
implicit none
!
integer, dimension(3) :: ncell_dummy
real(kind=PREC_DP) :: aver
real(kind=PREC_DP) :: sigma
integer           , dimension(:)    , allocatable :: at_site  ! Atom is at this site
!
!integer, dimension(3)              :: ncell_test
!real(kind=PREC_DP), dimension(3, 2):: test_dims
!integer           , dimension(3)   :: test_ilow        ! Unit cell dimensions in periodic
!integer           , dimension(3)   :: test_ihig        ! low and high inidce
!integer                            :: test_ncells     ! Number of cells in periodic crystal volume
!
aver = 0.0
sigma = 0.0
if(.not. (pdt_usr_ncell) ) then            ! User did not define number of atoms
   call estimate_ncells(ncell_dummy, pdt_dims, pdt_ilow, pdt_ihig, pdt_ncells)
   pdt_usr_ncell = .false.
!else                                       ! User did define number of unit cells, ch
!   call estimate_ncells(ncell_test, test_dims, test_ilow, test_ihig, test_ncells)
endif
!
if(.not. (pdt_usr_nsite) ) then            ! User did not define number of sites in unit cell
   call estimate_ncatom(aver, sigma, pdt_ilow, pdt_ihig, pdt_ncells)
endif
!
if(.not. (pdt_usr_atom) ) then             ! User did not specify unit cell content
   if(pdt_usr_nsite) aver = pdt_nsite      ! User specified number of atoms per cell
   allocate(at_site(cr_natoms))
   call find_average(aver, pdt_dims, pdt_ncells, pdt_nsite, pdt_asite, &
           pdt_itype, pdt_pos, cr_natoms, at_site)
   if(allocated(at_site)) deallocate(at_site)
   if(ier_num /= 0) return
else                                       ! User specified (partial) unit cell content
   if(pdt_usr_nsite) aver = pdt_nsite      ! User specified number of atoms per cell
   if(pdt_asite*pdt_ncells <cr_natoms) then ! Number of sites does not match
      call find_average(aver, pdt_dims, pdt_ncells, pdt_nsite, pdt_asite, &
           pdt_itype, pdt_pos, cr_natoms, at_site)
      if(ier_num /= 0) return
   elseif(pdt_nsite*pdt_ncells >2.*cr_natoms) then ! Number of sites does not match
      ier_num = -179
      ier_typ = ER_APPL
      ier_msg(1) = 'Expecting too many atoms '
      write(ier_msg(2),'(a,i9)') 'Expecting : ', pdt_nsite*pdt_ncells
      write(ier_msg(3),'(a,i9)') 'In crystal: ', cr_natoms
      return
   endif
endif
!
call perioditize_show
!
call map_to_aver
!
if(allocated(pdt_pos)) deallocate(pdt_pos)
!
cr_flags(:,1) = .true.       ! Crystal is supercell
cr_flags(1,2) = .false.      ! Not an asymmetric unit
cr_flags(2,2) = .false.      ! Not an asymmetric unit and Certain
cr_flags(:,3:5) = .true.     ! Periodic boundaries may be applied
cr_flags(:,6)   = .true.     ! Crystal is homogeneous and Certain about this
!
end subroutine perioditize_run
!
!*******************************************************************************
!
subroutine map_to_aver
!-
!  Map the current structure to an averaged and periodic crystal
!+
!
use crystal_mod
use celltoindex_mod
use chem_mod
use discus_allocate_appl_mod , only:alloc_crystal_scat, alloc_crystal_nmax
use metric_mod
use molecule_mod
!
use errlist_mod
!
implicit none
!
logical, parameter :: lspace =.TRUE.
real(kind=PREC_DP) :: EPS = 2.5
!logical :: lperiod               ! Test with periodic boundary conditions
integer :: i,j,k,l,m,n, nn       ! dummy counters
integer :: nprior                ! Number of atoms in crystal
integer :: natoms                ! Number of atoms in new crystal
integer :: ifail                 ! Atom type in case of failure
integer :: n_max                 ! MAX atom numbers for allocate
real(kind=PREC_DP) :: dmin
real(kind=PREC_DP) :: dist
real(kind=PREC_DP), dimension(3) :: u             ! fractional coordinates of atom
real(kind=PREC_DP), dimension(3) :: v             ! fractional coordinates of atom
integer, dimension(3) :: shift                    ! Integer shift vector
integer, dimension(3) :: icell                    ! Old atom is in this new cell
!
integer, dimension(  :), allocatable ::  tmp_iscat  ! (  1:NMAX)  !Atom type 0 to cr_nscat
integer, dimension(  :), allocatable ::  tmp_prop   ! (  1:NMAX)  !Property flag
integer, dimension(  :), allocatable ::  tmp_mole   ! (  1:NMAX)  !Atom is in this molecule
integer, dimension(:,:), allocatable ::  tmp_surf   ! (  1:NMAX)  !Atom is on this surface 
real(kind=PREC_DP)   , dimension(:,:), allocatable ::  tmp_magn   ! (  1:NMAX)  !Magnetic moment 
real(kind=PREC_DP)   , dimension(  :), allocatable ::  tmp_valu   ! (  1:NMAX)  !Magnetic moment 
real(kind=PREC_DP)   , dimension(:,:), allocatable ::  tmp_pos    ! (3,1:NMAX)  !Atom coordinates
!           
eps = 2.5
nprior = cr_natoms
if(pdt_usr_nsite) then
!write(*,*) ' pdt_usr_nsite is ', pdt_usr_nsite, pdt_nsite, pdt_ncells
   natoms = pdt_nsite * pdt_ncells
else
   natoms = pdt_nsite*pdt_ncells
!write(*,*) ' pdt_usr_nsite is ', pdt_usr_nsite, pdt_nsite, pdt_ncells
endif
!

do i=1, 3
   if(mod(pdt_ihig(i)-pdt_ilow(i) + 1,2)==0) then    ! Even number of unit cells
      pdt_ilow(i) = -(pdt_ihig(i)-pdt_ilow(i) + 1)/2
      pdt_ihig(i) =  (pdt_ihig(i)-pdt_ilow(i) + 1)/2 - 1
!        pdt_ihig(i) =  (pdt_ilow(i) + nint(werte(i+1)) - 1
   else
      pdt_ilow(i) = -(pdt_ihig(i)-pdt_ilow(i)    )/2
      pdt_ihig(i) =  (pdt_ihig(i)-pdt_ilow(i)    )/2
!        pdt_ilow(i) = -nint((werte(i+1)-1.)/2.)
!        pdt_ihig(i) = pdt_ilow(i) + nint(werte(i+1)) - 1
   endif
!      pdt_ncells = pdt_ncells * nint(werte(i+1))
enddo
!write(*,*) ' NATOMS ', natoms
!do l=1, pdt_nsite
!write(*,*) ' SITE   ', pdt_pos(:,l)
!enddo
!write(*,*) ' low  ', pdt_ilow
!write(*,*) ' high ', pdt_ihig
!
allocate(tmp_iscat(     1:natoms))
allocate(tmp_prop (     1:natoms))
allocate(tmp_mole (     1:natoms))
allocate(tmp_surf (0:3, 1:natoms))
allocate(tmp_magn (0:3, 1:natoms))
allocate(tmp_valu (     1:natoms))
allocate(tmp_pos  (3  , 1:natoms))
!
tmp_iscat  = -1                        ! Initialize with an error flag
!
!  Build the periodic crystal
!
m = 0
do k=pdt_ilow(3), pdt_ihig(3)
   do j=pdt_ilow(2), pdt_ihig(2)
      do i=pdt_ilow(1), pdt_ihig(1)
         do l=1, pdt_nsite
            m = m + 1
            tmp_pos(1,m) = real(i, kind=PREC_DP) + pdt_pos(1,l)
            tmp_pos(2,m) = real(j, kind=PREC_DP) + pdt_pos(2,l)
            tmp_pos(3,m) = real(k, kind=PREC_DP) + pdt_pos(3,l)
         enddo
      enddo
   enddo
enddo
!write(*,*) 'MIN/max x ', minval(tmp_pos(1,:)), maxval(tmp_pos(1,:))
!write(*,*) 'MIN/max y ', minval(tmp_pos(2,:)), maxval(tmp_pos(2,:))
!write(*,*) 'MIN/max z ', minval(tmp_pos(3,:)), maxval(tmp_pos(3,:))
!do m=1,pdt_nsite
!write(*,*) ' LL ', tmp_pos(:,m)
!enddo
!do m=natoms+1-pdt_nsite, natoms
!write(*,*) ' TR ', tmp_pos(:,m)
!enddo
!write(*,*) ' CL ', cr_dim(:,1)
!write(*,*) ' CT ', cr_dim(:,2), cr_dim(:,2)-cr_dim(:,1)
!write(*,*) 'SHFT', (nint(pdt_ilow(l)-cr_dim(l,1)), l=1,3)
!
! if necessary shift the old crystal
!
shift = nint(pdt_ilow(:)-cr_dim(:,1))
!write(*,*) ' SHIFT IN PERIOD ', shift
!write(*,*) ' CDIM  IN PERIOD ', cr_dim(:,1), cr_dim(:,2)
!write(*,*) ' PDT_ILOW PERIOD ', pdt_ilow(:), pdt_ihig(:)
if(any(shift/=0)) then
  v = real(shift,kind=PREC_DP)
  do i=1, cr_natoms
    cr_pos(:,i) = cr_pos(:,i) + v(:)
  enddo
endif
!
! Map the old crystal onto the new one
!
cr_icc = pdt_ihig - pdt_ilow + 1
!write(*,*) ' cr_icc ', cr_icc
!do m=1,5
!   icell = floor(cr_pos(:, m)+0.15) - pdt_ilow + 1
!   if(icell(1)<1          ) icell(1) = icell(1) + cr_icc(1)
!   if(icell(2)<1          ) icell(2) = icell(2) + cr_icc(2)
!   if(icell(3)<1          ) icell(3) = icell(3) + cr_icc(3)
!   call celltoindex(icell, 1, j)
!write(*,'(a, i9, 3f12.6, 3i4, i9)') ' ATOM ', m, cr_pos(:,m), icell, j
!enddo
!read(*,*) n
n = 0
nn = 0
do m=1, nprior
   icell = floor(cr_pos(:, m) + 0.15) - pdt_ilow + 1
   u = cr_pos(:, m)
   if(icell(1)<0          ) then
      icell(1) = icell(1) + cr_icc(1)
      u(1)     = u(1)     + cr_icc(1)
   elseif(icell(1)==0         ) then
      icell(1) = 1
   elseif(icell(1)>cr_icc(1)  ) then
      icell(1) = icell(1) - cr_icc(1)
      u(1)     = u(1)     - cr_icc(1)
   endif
   if(icell(2)<0          ) then
      icell(2) = icell(2) + cr_icc(2)
      u(2)     = u(2)     + cr_icc(2)
   elseif(icell(2)==0         ) then
      icell(2) = 1
   elseif(icell(2)>cr_icc(2)  ) then
      icell(2) = icell(2) - cr_icc(2)
      u(2)     = u(2)     - cr_icc(2)
   endif
   if(icell(3)<0          ) then
      icell(3) = icell(3) + cr_icc(3)
      u(3)     = u(3)     + cr_icc(3)
   elseif(icell(3)==0         ) then
      icell(3) = 1
   elseif(icell(3)>cr_icc(3)  ) then
      icell(3) = icell(3) - cr_icc(3)
      u(3)     = u(3)     - cr_icc(3)
   endif
   call celltoindex(icell, 1, j)
   dmin = 1.0E12
   do i=j, j+pdt_nsite-1    ! Loop over sites in this cell
      v = tmp_pos(:,i)
      dist = do_blen(lspace, u           , v)                ! Calculate distance
      if(dist < dmin) then                 ! Found new shortest distance
         dmin = dist
         k    = i                          ! Archive site
      endif
   enddo
      if(dmin < EPS) then                     ! Found close atom
      if(tmp_iscat(   k) == -1) then
         tmp_iscat(   k) = cr_iscat(1, m)
         tmp_prop (   k) = cr_prop (   m)
         tmp_mole (   k) = cr_mole (   m)
         tmp_surf (:, k) = cr_surf (:, m)
         tmp_magn (:, k) = cr_magn (:, m)
         tmp_valu (   k) = cr_valu (   m)
         tmp_pos  (:, k) = cr_pos  (:, m)
         cr_iscat(1,m) = -1                     ! Flag as done
           n = n + 1
         do j=mole_off(cr_mole(m)), mole_off(cr_mole(m))+mole_len(cr_mole(m))
            if(mole_cont(j) == m) then
               mole_cont(j) = k           ! Change atom number in Molecule
            endif
         enddo
      else
         nn = nn + 1
      endif
      else
!write(*,'(a, 2i9, 7f12.6, 3i3)') ' Failed ', m, cr_iscat(1, m), cr_pos(:,m), u, dmin, icell
!read(*,*) k
      endif
enddo
!write(*,'(a, i9, a, 2i9)') ' Mapped ', n, ' out of old/new ', nprior, natoms
!j=0
!do i=1, nprior
!  if(cr_iscat(1,i)>=0) j=j+1
!enddo
!write(*,'(a, i9)') ' mismatched old ', j
!j=0
!do i=1, natoms
!  if(tmp_iscat(i)< 0) j=j+1
!enddo
!write(*,'(a, i9, i9)') ' mismatched new ', j, nn
!read(*,*) m
! use : k index of closest atom
!Qn = 0
!Qsoften: do l=1,2
!Q   loop_new: do i=1, natoms
!Q      if(tmp_iscat(i)>=0) cycle loop_new          ! Already occupied
!Q      k = 0
!Q      dmin = 1.0E12
!Q      dist = 1.0E12
!Q      v = tmp_pos(:,i)
!Q      loop_old:do m=1, nprior
!Q         if( cr_iscat(1,m) <0) cycle loop_old          ! Already transfered
!Q         u = cr_pos(:, m)
!Q         dist = do_blen(lspace, u, v)                ! Calculate distance
!Q!if(            m==1713) then !abs(cr_pos(1,m)+9.2)<0.1 .and. abs(cr_pos(2,m)+10.0)<0.1) then
!Q!if(m<6) then
!Q!  write(*,*) ' ATOM ', cr_pos(:,m), cr_iscat(1,m), dist, dmin
!Q!  write(*,*) ' ATOM ', cr_pos(:,m)+20, cr_iscat(1,m), dist, dmin
!Q!endif
!Q         if(dist < dmin) then                 ! Found new shortest distance
!Q            dmin = dist
!Q            k    = m                          ! Archive site
!Q         endif
!Q         lperiod = .false.                    ! Test if periodic boundary conditions may be needed
!Q         do nn = 1, 3
!Q            if(u(nn)<pdt_ilow(nn)) then
!Q               u(nn) = u(nn) + pdt_ihig(nn) - pdt_ilow(nn) + 1
!Q               lperiod = .true.
!Q            endif
!Q         enddo 
!Q         if(lperiod) then    ! Test with periodic boundary conditions
!Q            dist = do_blen(lspace, u, v)                ! Calculate distance
!Q            if(dist < dmin) then                 ! Found new shortest distance
!Q               dmin = dist
!Q               k    = m                          ! Archive site
!Q            endif
!Q         endif
!Q      enddo loop_old
!Q      if(dmin < EPS) then                     ! Found close atom
!Q         tmp_iscat(   i) = cr_iscat(1, k)
!Q         tmp_prop (   i) = cr_prop (   k)
!Q         tmp_mole (   i) = cr_mole (   k)
!Q         tmp_surf (:, i) = cr_surf (:, k)
!Q         tmp_magn (:, i) = cr_magn (:, k)
!Q         tmp_pos  (:, i) = cr_pos  (:, k)
!Q         cr_iscat(1,k) = -1                     ! Flag as done
!Q           n = n + 1
!Q         do j=mole_off(cr_mole(k)), mole_off(cr_mole(k))+mole_len(cr_mole(k))
!Q            if(mole_cont(j) == k) then
!Q               mole_cont(j) = i           ! Change atom number in Molecule
!Q            endif
!Q         enddo
!Q      else
!Q      endif
!Q   enddo loop_new
!Q   if(n==nprior) exit soften
!Q   eps = eps * 2.0
!Qenddo soften
!
cr_iscat = 0
cr_prop  = 0
cr_mole  = 0
cr_surf  = 0
cr_magn  = 0.0
cr_valu  = 0.0
cr_pos   = 0.0
cr_iscat(1,1:natoms) = tmp_iscat(  1:natoms)
cr_prop (  1:natoms) = tmp_prop (  1:natoms)
cr_mole (  1:natoms) = tmp_mole (  1:natoms)
cr_surf (:,1:natoms) = tmp_surf (:,1:natoms)
cr_magn (:,1:natoms) = tmp_magn (:,1:natoms)
cr_valu (  1:natoms) = tmp_valu (  1:natoms)
cr_pos  (:,1:natoms) = tmp_pos  (:,1:natoms)
cr_natoms = natoms
!
cr_icc (1) = pdt_ihig(1) - pdt_ilow(1) + 1
cr_icc (2) = pdt_ihig(2) - pdt_ilow(2) + 1
cr_icc (3) = pdt_ihig(3) - pdt_ilow(3) + 1
cr_ncatoms = pdt_nsite
chem_purge = .FALSE.     ! Crystal dimension should allow periodic boundary
if(cr_icc(1)>1) chem_period(1) = .TRUE.
if(cr_icc(2)>1) chem_period(2) = .TRUE.
if(cr_icc(3)>1) chem_period(3) = .TRUE.
chem_quick = .TRUE.
!
if(n /= nprior) then
   ier_num =  -178
   ier_typ = ER_APPL
   write(ier_msg(1),'(a,i8,a,i8)') 'Map ', n, ' out of ', nprior
   write(ier_msg(2),'(a        )') 'Non matched flagged as ''FAIL'' '
   ifail = cr_nscat + 1
   if(ifail> MAXSCAT) then
      n_max = max(cr_natoms, NMAX)
      call alloc_crystal_scat ( ifail )
      call alloc_crystal_nmax ( n_max )
   endif
! 
   do i=1, nprior
      if(cr_iscat(1,i)<0) then
         cr_iscat(1,i) = ifail
      endif
   enddo
   cr_at_lis(ifail) = 'FAIL'
endif
!
if(natoms>nprior) then
   do i=1, natoms
      if(cr_iscat(1,i)<0)  then
         cr_iscat(1,i)  = 0
         cr_prop(i)   = 0
         cr_mole(i)   = 0
         cr_surf(:,i) = 0
         cr_magn(:,i) = 0.0
      endif
   enddo
endif
!
deallocate(tmp_iscat)
deallocate(tmp_prop )
deallocate(tmp_mole )
deallocate(tmp_surf )
deallocate(tmp_magn )
deallocate(tmp_valu )
deallocate(tmp_pos  )
!
end subroutine map_to_aver
!
!*******************************************************************************
!
subroutine perioditize_show
!-
!  Show current setting for the Perioditize menu
!+
use prompt_mod
!
implicit none
!
integer :: i   ! Dummy loop index
!
write(output_io,*)
write(output_io,'(a)') ' Perioditize menu'
write(output_io,*)
if(pdt_usr_nsite) then
   write(output_io,'(a, i5)')                                                   &
      ' User defined sites per unit cell        :', pdt_nsite
else
   write(output_io,'(a    )')                                                   &
      ' Number of sites determined automatically'
   write(output_io,'(a, i5)')                                                   &
      ' Number of sites estimated               ', pdt_nsite
endif
if(pdt_usr_atom) then
   do i=1, pdt_usite
      write(output_io, '(a,3(f12.3,2x))')                                       &
      ' User defined site position              :', pdt_pos  (1:3,i)
      write(output_io, '(a,4x,20i4)')                                           &
      ' User defined site content               :', pdt_itype(1:pdt_itype(0,i),i)
   enddo
else
   write(output_io, '(a)')                                                      &
      ' Sites determined automatically'
   if(pdt_nsite>0) then
   do i=1, pdt_nsite
      write(output_io, '(a,3(f12.3,2x))')                                       &
      ' estimated    site position              :', pdt_pos  (1:3,i)
      write(output_io, '(a,4x,20i4)')                                           &
      ' Estimated    site content               :', pdt_itype(1:pdt_itype(0,i),i)
   enddo
   endif
endif
if(pdt_usr_ncell) then
   write(output_io,'(a,3i4)')                                                   &
      ' User defined number of unit cells       :', pdt_ihig-pdt_ilow + 1
else
   write(output_io,'(a, 3i4)')                                                  &
      ' Unit cells determined automatically     :', pdt_ihig-pdt_ilow + 1 
endif
write(output_io,*)
!
end subroutine perioditize_show
!
!*******************************************************************************
!
subroutine perioditize_reset
!
!
implicit none
!
pdt_ncells = 0
pdt_asite  = 0
pdt_nsite  = 0
pdt_usite  = 0
pdt_ilow   = 0
pdt_ihig   = 0
pdt_dims   = 0.0D0
pdt_usr_nsite = .FALSE.
pdt_usr_ncell = .FALSE.
pdt_usr_atom  = .FALSE.
!
if(allocated(pdt_pos)) deallocate(pdt_pos)
if(allocated(pdt_itype)) deallocate(pdt_itype)
!
allocate(pdt_pos  (1:3, 1:1))
allocate(pdt_itype(0:1, 1:1))
pdt_pos   = 0.0
pdt_itype(0,1) =  1
pdt_itype(1,1) = -1
!
end subroutine perioditize_reset
!
!*******************************************************************************
!
end module perioditize_mod
