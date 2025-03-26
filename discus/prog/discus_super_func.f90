module superspace_func_mod
!-
!  Menus/ Functions for superspace
!+
!
contains
!
!*******************************************************************************
!
subroutine super_menu
!-
!  Main superspace menu
!+
!
use discus_kdo_common_mod
!
use calc_expr_mod
use class_macro_internal
use doact_mod
use errlist_mod
use kdo_all_mod
use lib_errlist_func
use lib_length
use lib_macro_func
use prompt_mod
use str_comp_mod
use sup_mod
!
implicit none
!
CHARACTER(len=11)          :: befehl
CHARACTER(LEN=LEN(prompt)) :: orig_prompt
CHARACTER(LEN=PREC_STRING) :: line, zeile
INTEGER :: lp, length, lbef
INTEGER :: indxg
logical :: success
logical :: lend   ! User typed exit command
!
call no_error
!
orig_prompt = prompt
prompt = prompt (1:len_str(prompt) ) //'/super'
lend = .false.
!
loop_main: do while (.not.lend)                 ! Main super loop
!
   if_error: IF (ier_num.ne.0) THEN
      CALL errlist
      IF (ier_sta.ne.ER_S_LIVE) THEN
         IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
            IF(sprompt /= prompt ) THEN
               ier_num = -10
               ier_typ = ER_COMM
               ier_msg(1) = ' Error occured in exp2pdf menu'
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
   ENDIF if_error
!
   CALL get_cmd(line, length, befehl, lbef, zeile, lp, prompt)
   IF (ier_num /= 0) cycle loop_main

   IF(line == ' '      .or. line(1:1) == '#' .or. &
      line == char(13) .or. line(1:1) == '!'        ) cycle loop_main
!                                                                       
!     ----search for "="                                                
!                                                                       
   indxg = index (line, '=')
   IF (indxg.ne.0                                              &
        .AND..NOT. (str_comp (befehl, 'echo',   2, lbef, 4) )    &
        .AND..NOT. (str_comp (befehl, 'system', 2, lbef, 6) )    &
        .AND..NOT. (str_comp (befehl, 'help',   2, lbef, 4) .OR. &
                    str_comp (befehl, '?   ',   2, lbef, 4) )    &
        .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
!     ------evaluate an expression and assign the value to a variable   
!                                                                       
            CALL do_math(line, indxg, length)
      cycle loop_main
   endif
!
!
! - Test common menu entries
!
   CALL discus_kdo_common(befehl, lbef, line, length, zeile, lp, 'super' , &
                          lend, success)
   if(success) cycle loop_main

!--SUPER   specific commands
!
!  --- exit
!
   if(str_comp(befehl, 'exit', 2, lbef, 4)) then
      lend = .TRUE.
!
!  Assign properties to dummy atoms
!
   elseif(str_comp(befehl, 'assign', 2, lbef, 6)) then
      call super_set(zeile, lp)
!
!  Reset 
!
   elseif(str_comp(befehl, 'reset', 2, lbef, 5)) then
      call super_reset
!
!  Perform actual task
!
   elseif(str_comp(befehl, 'run', 2, lbef, 3)) then
      call super_run
!
!  Show current settings
!
   elseif(str_comp(befehl, 'show', 2, lbef, 4)) then
      call super_show
!
!  Perform Set parameters
!
   elseif(str_comp(befehl, 'set', 2, lbef, 3)) then
      call super_set(zeile, lp)
!                                                                       
!     ----Unknown exp2pdf command, try general commands
!                                                                       
   else
      call kdo_all(befehl, lbef, zeile, lp)
   endif
!
enddo loop_main
!
prompt = orig_prompt
!
end subroutine super_menu
!
!*******************************************************************************
!
subroutine super_set (zeile, lp)
!-
!  Set parameters
!+
!
use crystal_mod
use discus_allocate_appl_mod
use get_iscat_mod
use superspace_mod
!use modify_mod
!
use ber_params_mod
use errlist_mod
use get_params_mod
use lib_metric_mod
use precision_mod
use string_convert_mod
use str_comp_mod
use take_param_mod
!
implicit none
!
character(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp
!
integer, parameter :: MAXW  = 6
integer, parameter :: MAXWW = 6
!logical, parameter :: LOLD = .TRUE.
!
character(len=4) :: at_new
character(len=PREC_STRING) :: string
character(LEN=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
real(kind=PREC_DP)        , dimension(MAXW) :: werte   ! Calculated values
!
character(LEN=PREC_STRING), dimension(10        ) :: ccpara
integer                   , dimension(10        ) :: llpara
real(kind=PREC_DP)        , dimension(10        ) :: wwerte   ! Calculated values
real(kind=PREC_DP)        , dimension(3         ) :: v        ! A vector
real(kind=PREC_DP) :: vv ! Vector length
real(kind=PREC_DP) :: up_up ! 'UP' length crenel
integer :: ianz, iianz
integer :: i1, i2, j     ! Dummy indices
!integer :: isite         ! Current site
logical, dimension(:), allocatable :: l_isite
logical :: lselect
!
integer, parameter :: NOPTIONAL = 18
integer, parameter :: O_CURRENT = 1
integer, parameter :: O_GROUP   = 2
integer, parameter :: O_UP      = 3
integer, parameter :: O_PHASE   = 4
integer, parameter :: O_SITE    = 5
integer, parameter :: O_PSEUDO  = 6
integer, parameter :: O_CHAR    = 7
integer, parameter :: O_FUNC    = 8
integer, parameter :: O_DISP    = 9
integer, parameter :: O_REPLACE = 10
integer, parameter :: O_FILE    = 11
integer, parameter :: O_QVEC    = 12
integer, parameter :: O_AMP     = 13
integer, parameter :: O_VEC     = 14
integer, parameter :: O_NUMBER  = 15
integer, parameter :: O_PROB    = 16
integer, parameter :: O_OLD     = 17
integer, parameter :: O_NEW     = 18
character(LEN=  12), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 2 ! Number of values to calculate 
!
data oname  / 'current', 'group', 'up'   , 'phase' , 'site' , 'pseudo', 'character', 'function', 'displacement',    &
              'replace', 'file' , 'qvec'  , 'amp'      , 'vec'     , 'number'      ,    &
              'prob'   , 'old'  , 'new'                                             /
data loname /  7    ,  5         ,  2,           5,       4,       6,      9,           8,          12,               &
               7    ,  4         ,  4       ,  3   ,  3   , 6   ,  4, 3, 3/
opara  =  (/ '1.0000        ', '1.0000        ', '0.5000        ', '0.0000        ', '0.0000        ', &
             'VOID          ', 'displa        ', &
             'sine          ', '[0,0,0,0]     ', '[VOID,0.0,1.0]', 'internal      ', &
             '[1,0,0]       ', '0.0000        ', '[1,0,0]       ', '1             ', &
             '[0.0,1.0]     ', 'VOID          ', 'VOID          '                   /)   ! Always provide fresh default values
lopara =  (/  6         ,  6         ,  6         ,  6         ,  6,           4,           6, &
              6         ,  9         ,    &
              14        ,  9         ,  7,           6,           7         ,  1         ,    &
              9         ,  4,           4                                           /)
owerte =  (/  1.0, 1.0, 0.50, 1.0, 1.0, 0.0,      0.0,      0.0    ,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0  /)
!
call get_params(zeile, ianz, cpara, lpara, maxw, lp)
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
!
if(lpresent(O_FILE)) then   ! Set file:
   sup_file = opara(O_FILE)
   return
endif
!
if(lpresent(O_CURRENT)) then   ! Set current wave:
   sup_current = nint(owerte(O_CURRENT))
   if(sup_current>sup_nwaves) then
      call alloc_super(cr_ncatoms, cr_nscat, sup_current, sup_ngroups)
      sup_nwaves = sup_current
   endif
!
   if(lpresent(O_GROUP)) then   ! Set group
      j = nint(owerte(O_GROUP))
      if(j>sup_ngroups) then
         i1 = max(1, j, sup_ngroups)
         call alloc_super(cr_ncatoms, cr_nscat, sup_current, i1)
      endif
      sup_group(sup_current) = j
   else
      sup_group(sup_current) = 1
   endif
!
   return
endif
!
!
if(lpresent(O_QVEC)) then   ! Set qvec:
   string=opara(O_QVEC)
   j = lopara(O_QVEC)
   call get_optional_multi(MAXWW, string, j, werte, iianz)
   if(ier_num/=0) return
   sup_qvec(:, sup_current) = werte(1:3)
   return
endif
!
cond_site: if(lpresent(O_SITE)) then
   if(allocated(l_isite)) deallocate(l_isite)
   allocate(l_isite(1:cr_ncatoms))
!
   if(cr_ncatoms>ubound(sup_char,1 ) .or.  &
      cr_nscat  >ubound(sup_old,1 )      ) then
      call alloc_super(cr_ncatoms, cr_nscat, sup_nwaves, sup_ngroups)
   endif
   j = lopara(O_SITE)
   if(opara(O_SITE)(1:1)=='[') then    ! Multiple sites
      i1 = 2
      i2 = lopara(O_SITE) - 1
   else
      i1 = 1
      i2 = lopara(O_SITE)
   endif
   string = opara(O_SITE)(i1:i2)
   j = i2-i1+1
   call get_params(string, iianz, ccpara, llpara, MAXWW, j)
   if(ier_num/=0) return
   l_isite = .false.
   if(str_comp(string, 'all', 3, j, 3)) then
      l_isite = .true.
   else
      call ber_params(iianz, ccpara, llpara, wwerte, maxww)
      if(ier_num/=0) return
      do j=1, iianz
         l_isite(nint(wwerte(j))) = .true.
      enddo
   endif
!  isite = nint(owerte(O_SITE))
!
   if(lpresent(O_PSEUDO)) then   ! assign type:   must always be present
      string = opara(O_PSEUDO)(2:lopara(O_PSEUDO)-1)
      j = lopara(O_PSEUDO) - 2
      call get_params(string, iianz, ccpara, llpara, MAXWW, j)
      do j=1, cr_ncatoms
         if(l_isite(j)) then
            do i1=1, iianz
               sup_atom(i1, j, sup_current, sup_group(sup_current)) = ccpara(i1)(1:llpara(i1))
               call do_cap(sup_atom(i1, j, sup_current, sup_group(sup_current)))
            enddo
         endif
      enddo
   endif
!
!  CHARACTER
!
   if(lpresent(O_CHAR)) then     ! "character:"
      if(str_comp(opara(O_CHAR), 'displacement', 3, lopara(O_CHAR), 12)) then
         i1 = SUP_DISP
      elseif(str_comp(opara(O_CHAR), 'density', 3, lopara(O_CHAR), 7)) then
         i1 = SUP_DENS
      endif
!     sup_char = 0
      do j=1, cr_ncatoms
         if(l_isite(j)) then
            sup_char(j, sup_current) = i1
         endif
      enddo
   endif
!
!  FUNCTION
!
   if(lpresent(O_FUNC)) then
      up_up = 0.50_PREC_DP
      if(str_comp(opara(O_FUNC), 'sine', 3, lopara(O_FUNC), 4)) then
         i1 = SUP_SINE
      elseif(str_comp(opara(O_FUNC), 'crenel', 3, lopara(O_FUNC), 7)) then
         i1 = SUP_CREN 
         if(lpresent(O_UP)) then
            up_up = owerte(O_UP)
         else
            up_up = 0.50_PREC_DP
         endif
      endif
!     sup_func = 0
      do j=1, cr_ncatoms
         if(l_isite(j)) then
            sup_func(j, sup_current) = i1
            sup_func_p(j, sup_current) = up_up
         endif
      enddo
   endif
!
!  DISPLACEMENT VECTOR
!  displacement:[amp:0.20, vec:[1.0,1.0,0.0], phase:0.0, number:1]
!
   if(lpresent(O_DISP)) then
      string = opara(O_DISP)(2:lopara(O_DISP)-1)
      i1 = len_trim(string)
      call get_params(string, ianz, cpara, lpara, maxw, i1)
      call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                        oname, loname, opara, lopara, lpresent, owerte)
!     Get amplitude part
      if(lpresent(O_AMP)) then
         ccpara(1) = opara(O_AMP)
         llpara(1) = lopara(O_AMP)
      else
         ccpara(1) = cpara(1)
         llpara(1) = lpara(1)
      endif
!     Get phase part
      if(lpresent(O_PHASE)) then
         ccpara(5) = opara(O_PHASE)
         llpara(5) = lopara(O_PHASE)
      else
         if(ianz>=5) then
            ccpara(5) = cpara(5)
            llpara(5) = lpara(5)
         else
            ccpara(5) = '1.0'
            llpara(5) =  3
         endif
      endif
!     Get number part
      if(lpresent(O_NUMBER)) then
         ccpara(6) = opara(O_NUMBER)
         llpara(6) = lopara(O_NUMBER)
      else
         if(ianz==6) then
            ccpara(6) = cpara(6)
            llpara(6) = lpara(6)
         else
            ccpara(6) = '1.0'
            llpara(6) =  3
         endif
      endif
!     get displacment vector
      if(lpresent(O_VEC)) then
         string = opara(O_VEC)(2:lopara(O_VEC)-1)
         i1 = len_trim(string)
         call get_params(string, ianz, cpara, lpara, maxw, i1)
         ccpara(2:4) = cpara(1:3)
         llpara(2:4) = lpara(1:3)
      else
         ccpara(2:4) = cpara(2:4)
         llpara(2:4) = lpara(2:4)
      endif
!     Calculate corresponding values
      iianz = 6
      call ber_params(iianz, ccpara, llpara, wwerte, maxww)
      if(ier_num/=0) return
      i1 = nint(wwerte(6))
      do j=1,cr_ncatoms
         if(l_isite(j)) then
            sup_phase(i1,j, sup_current) = wwerte(5)
            sup_ampl(1:4, i1, j, sup_current) = wwerte(1:4)
            v = wwerte(2:4)
            vv = lib_blen(cr_gten, v)
            if(vv>0.0) then
               sup_ampl(2:4, i1, j, sup_current) = sup_ampl(2:4, i1, j, sup_current)*sup_ampl(1, i1, j, sup_current)/vv
            endif
         end if
      enddo
   endif
!
!  REPLACEMENT 
!  replace:[old:[C], new:[Si], prob:[0.05, 0.999], phase:0.00]
!  
   if(lpresent(O_REPLACE)) then
      string = opara(O_REPLACE)(2:lopara(O_REPLACE)-1)
      j = lopara(O_REPLACE) - 2
      call get_params(string, ianz, cpara, lpara, maxw, j)
      call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                        oname, loname, opara, lopara, lpresent, owerte)
      at_new = cpara(1)(1:lpara(1))
      i1     = len_trim(at_new)
      ccpara(1) = '0.0'
      llpara(1) = 3
      if(lpresent(O_OLD)) then
         string = opara(O_OLD)(2:lopara(O_OLD)-1)
         i1 = len_trim(string)
         call get_params(string, ianz, cpara, lpara, maxw, i1)
         lselect = .FALSE.
         call get_iscat(ianz, cpara, lpara, werte, maxw, lselect)
         if(nint(werte(1))==-1) then
            do j=1,cr_ncatoms
               if(l_isite(j)) then
                  sup_old(:,j, sup_current) = .TRUE.
               endif
            enddo
         else
            do j=1,cr_ncatoms
               if(l_isite(j)) then
                  sup_old(:,j, sup_current) = .false.
                  do i1=1, ianz
                     sup_old(nint(werte(i1)),j, sup_current) = .TRUE.
                  enddo
               endif
            enddo
         endif
      else
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'Mandatory ''old:[]'' is missing'
         return
      endif
      if(lpresent(O_NEW)) then
         string = ' '
         string(1:lopara(O_NEW)-2) = opara(O_NEW)(2:lopara(O_NEW)-1)
         i1 = len_trim(string)
         call get_params(string, ianz, cpara, lpara, maxw, i1)
         lselect = .TRUE.
         call get_iscat(ianz, cpara, lpara, werte, maxw, lselect)
         if(cr_ncatoms>ubound(sup_char,1) .or.  &
            cr_nscat  >ubound(sup_old,1)      ) then
            call alloc_super(cr_ncatoms, cr_nscat, sup_nwaves, sup_ngroups)
         endif
         if(nint(werte(1))==-1) then
            do j=1,cr_ncatoms
               if(l_isite(j)) then
                  do i1=0, cr_nscat
                     sup_new(i1,j, sup_current) = i1
                  enddo
               endif
            enddo
         else
            do j=1,cr_ncatoms
               if(l_isite(j)) then
                  sup_new(:,j, sup_current) = -1
                  sup_irepl(j, sup_current) = ianz
                  do i1=1, ianz
                     sup_new(ianz-1,j, sup_current) = nint(werte(i1))
                  enddo
               endif
            enddo
         endif
      else
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'Mandatory ''new:[]'' is missing'
         return
      endif
      if(lpresent(O_PROB)) then
         string = opara(O_PROB)(2:lopara(O_PROB)-1)
         i1 = len_trim(string)
         call get_params(string, ianz, cpara, lpara, maxw, i1)
         ccpara(2:3) = cpara(1:2)
         llpara(2:3) = lpara(1:2)
      else
         ccpara(2:3) = cpara(2:3)
         llpara(2:3) = lpara(2:3)
      endif
#
      if(lpresent(O_PHASE)) then
         ccpara(4) = opara(O_PHASE)
         llpara(4) = lopara(O_PHASE)
      else
         if(ianz>=4) then
            ccpara(4) = cpara(5)
            llpara(4) = lpara(5)
         else
            ccpara(4) = '0.0'
            llpara(4) =  3
         endif
      endif
      iianz = 4
      call ber_params(iianz, ccpara, llpara, wwerte, maxww)
      if(ier_num/=0) return
      loop_atom: do j=1, cr_ncatoms
         if(l_isite(j)) then
!           sup_irepl(j) = nint(werte(2))
!           sup_repl (j) = at_new
            sup_prob (1,j, sup_current) = wwerte(2)
            sup_prob (2,j, sup_current) = wwerte(3)
            sup_phase(1,j, sup_current) = wwerte(4)
         endif
      enddo loop_atom
   endif
!else
!   ier_num = -6
!   ier_typ = ER_COMM
!   return
!endif
else cond_site
   ier_num = -6
   ier_typ = ER_COMM
   return
endif cond_site
!
if(allocated(l_isite)) deallocate(l_isite)
!
end subroutine super_set
!
!*******************************************************************************
!
subroutine super_run 
!-
!  Run actual modulation
!+
!
use celltoindex_mod
use crystal_mod
use read_internal_mod
use superspace_mod
use super_waves_mod
!
use lib_functions_mod
!
use precision_mod
use wink_mod
!
implicit none
!
interface
  real(kind=PREC_DP    ) function sup_fun_gen (amp, amp0, arg1, arg2)
  USE precision_mod
  real(kind=PREC_DP), intent(in) :: amp
  real(kind=PREC_DP), intent(in) :: amp0
  real(kind=PREC_DP), intent(in) :: arg1
  real(kind=PREC_DP), intent(in) :: arg2
  end function sup_fun_gen 
end interface
!
!interface
!  real(kind=kind(1.0_D0)) function sup_fun_sine(amp, amp0, arg1)
!  USE precision_mod
!  real(kind=PREC_DP), intent(in) :: amp
!  real(kind=PREC_DP), intent(in) :: amp0
!  real(kind=PREC_DP), intent(in) :: arg1
!  real(kind=PREC_DP), intent(in) :: arg2
!  end function sup_fun_sine
!end interface
!
!interface
!  real(kind=PREC_DP    ) function sup_fun_cren(amp, amp0, arg1, arg2)
!  USE precision_mod
!  real(kind=PREC_DP), intent(in) :: amp
!  real(kind=PREC_DP), intent(in) :: amp0
!  real(kind=PREC_DP), intent(in) :: arg1
!  real(kind=PREC_DP), intent(in) :: arg2
!  end function sup_fun_cren
!end interface
!
PROCEDURE(sup_fun_gen   )   , POINTER :: p_sub_function   => NULL()
!
real(kind=PREC_DP), parameter :: TOL = 0.00001_PREC_DP
!
integer  :: dim_natoms     ! Total number of atoms 
integer  :: dim_ncatoms    ! Number of atoms per uni cell
integer  :: dim_nscat      ! Number of chemical species
integer  :: dim_nanis      ! Number of ADPs
integer  :: dim_n_mole     ! Number of molecules
integer  :: dim_n_type     ! Number of molecule types
integer  :: dim_n_atom     ! Total number of atoms in molecules
integer  :: natoms     ! Total number of atoms 
integer  :: ncatoms    ! Number of atoms per unit cell
integer  :: nscat      ! Number of chemical species
integer  :: nanis      ! Number of ADPs
integer  :: n_mole     ! Number of molecules
integer  :: n_type     ! Number of molecule types
integer  :: n_atom     ! Total number of atoms in molecules
!
integer  :: current  ! Current wave
integer, dimension(sup_ngroups) :: group_iscat
integer  :: iatom  ! Dummy atom number from "super file"
integer  :: katom  ! Dummy atom number from "super file"
integer  :: jatom  ! Dummy atom number from  crystal structure
real(kind=PREC_DP), dimension(3) :: arg
real(kind=PREC_DP), dimension(3) :: vector
!
integer               :: jj   ! Counter atoms per cell
integer, dimension(3) :: ic   ! unit cell of crystal atom
integer               :: is   ! Site of crystal atom
!integer               :: ip   ! Pseudo atom counter
integer               :: ii   ! Diplacement number
integer               :: ig    ! group counter
integer               :: istep ! Number of crystal atoms per pseudo atom
real(kind=PREC_DP)    :: xstep ! Number of crystal atoms per pseudo atom
!
character(len=4), dimension(:), allocatable :: rd_at_lis
REAL(kind=PREC_DP)  , DIMENSION(3)  :: rd_cr_pos
INTEGER             , dimension(3)  :: rd_cr_iscat
INTEGER                             :: rd_cr_prop
INTEGER             , DIMENSION(0:3):: rd_cr_surf
REAL(kind=PREC_DP)  , DIMENSION(0:3):: rd_cr_magn
INTEGER                             :: rd_cr_mole
INTEGER                             :: rd_cr_moleatom
!
real(kind=PREC_DP), dimension(:,:,:), allocatable :: rep_ampl
real(kind=PREC_DP)                              :: r1
real(kind=PREC_DP)                              :: phase
real(kind=PREC_DP)                              :: prob
real(kind=PREC_DP), dimension(3) :: v
!
p_sub_function => sup_fun_sine          ! Default to sine function
!
! For replacement waves, calculate amplitude and average probability
!
allocate(rep_ampl(2, cr_ncatoms, sup_nwaves))
current = 1
do iatom=1, cr_ncatoms
  do current=1, sup_nwaves
      rep_ampl(1,iatom, current) = (sup_prob(2,iatom, current) - sup_prob(1,iatom, current))*0.5_PREC_DP ! Amplitude
      rep_ampl(2,iatom, current) = (sup_prob(2,iatom, current) + sup_prob(1,iatom, current))*0.5_PREC_DP ! Average
   enddo
enddo
!
!  Read list of pseudo atoms
!
call readstru_size_int(sup_file, dim_natoms, dim_ncatoms,      &
     dim_nscat, dim_nanis, dim_n_mole, dim_n_type, dim_n_atom, &
     natoms, ncatoms, nscat, nanis, n_mole, n_type, n_atom)
!
! Get pseudo atom names 
allocate(rd_at_lis(0:dim_nscat))
call stru_get_atlis(sup_file, dim_nscat, rd_at_lis)
!
! Determine steps:  Crystal atoms / Pseudo atoms / N_groups
istep = 1
xstep = real(cr_natoms)/real(natoms)
if(abs(xstep-real(nint(xstep)))<TOL) then
   istep = cr_natoms/natoms!*sup_ngroups
endif
!write(*,*) ' STEPS ', istep, cr_natoms, natoms, dim_natoms, sup_nwaves, sup_ngroups
!write(*,*) ' STEPS ', istep, cr_ncatoms, dim_ncatoms, dim_natoms, sup_nwaves, sup_ngroups
!write(*,*) ' G1    ', sup_group(1)
!write(*,*) ' G1    ', sup_group(2)
!write(*,*) ' GROUP 1:2, 1, 1, : ', sup_atom(1:2, 1, 1, sup_group(1))
!write(*,*) ' GROUP 1:2, 1, 2, : ', sup_atom(1:2, 1, 2, sup_group(2))
!write(*,*) ' GROUP 1:2, 2, 1, : ', sup_atom(1:2, 2, 1, sup_group(1))
!write(*,*) ' GROUP 1:2, 2, 2, : ', sup_atom(1:2, 2, 2, sup_group(2))
!
!  Main loop
!
loop_main: do iatom=1, dim_natoms, sup_ngroups
!
!  Read "pseudo atoms" from "super file"
!
   do ig=1, sup_ngroups
      katom = iatom + ig - 1
      call struc_read_one_atom_internal(sup_file, katom,  &
                 rd_cr_pos, rd_cr_iscat, rd_cr_prop, rd_cr_surf, &
                 rd_cr_magn, rd_cr_mole, rd_cr_moleatom )
      group_iscat(ig) = rd_cr_iscat(1)
   enddo
   loop_cell: do jj=1, istep
      jatom = (iatom-1)*istep + jj
      call indextocell(jatom, ic, is)
      if(all(sup_atom(1,is,:,:)==' ')) cycle loop_cell   ! nothing defined for this site
      v = cr_pos(:,jatom) !- cr_dim0(:,1)
!
      vector = 0.0_PREC_DP
      prob   = 1.0_PREC_DP
!
      loop_current: do current=1, sup_nwaves   ! Accumulate vectors / probabilities
         if(sup_atom(1,is,current,sup_group(current))==' ') cycle loop_current   ! nothing defined for this current wave
            phase = 0.0_PREC_DP
!GROUP   if(sup_atom(1,is, current, sup_group(sup_current))==rd_at_lis(rd_cr_iscat(1))) then
         if(sup_atom(1,is, current, sup_group(current))==rd_at_lis(group_iscat(sup_group(current)))) then
            phase = 0.0_PREC_DP
!GROUP   elseif(sup_atom(2,is, current, sup_group(sup_current))==rd_at_lis(rd_cr_iscat(1))) then
         elseif(sup_atom(2,is, current, sup_group(current))==rd_at_lis(group_iscat(sup_group(current)))) then
            phase = 0.5_PREC_DP
         endif
         do ii=1,3
         arg(ii) = (v(1)*sup_qvec(1, current) + v(2)*sup_qvec(2, current) + v(3)*sup_qvec(3, current) + &
                  sup_phase(ii,is, current)+ phase                             )!*zpi
         enddo
!if(jatom<12) then
! write(*,'(3i5,12f8.3)') iatom, jatom, is, arg(1), sup_phase(1,is,current), sin(arg(1)), sup_ampl(2:3, 1,is,current), &
! arg(2), sup_phase(2,is,current), sin(arg(2)), sup_ampl(2:3, 2,is,current), v(1:2)
!endif
!if(iatom<=4) then
!if(iatom<=200 .and. ic(2)==1 .and. ic(3)==1) then
!write(*,'(a, 3i4)') 'LOOP CELL', jj, current, sup_group(current)
!write(*,'(a,i7,3f8.2, i3,1x, 2a5, i3, 2a6, 2(2x,3f8.3))') ' PSEUDO ',iatom, rd_cr_pos, &
!  group_iscat(current), rd_at_lis(group_iscat(sup_group(current))),&
!  ' HOST ', jatom, sup_atom(:,is,current, sup_group(current)), v, sup_qvec(:, current)
!write(*,'(a,i7, 2(2a4,l3),3f10.4)')      ' sup_at ', is, &
!  sup_atom(1,is,current, (sup_group(current))), rd_at_lis(group_iscat(sup_group(current))), &
!  sup_atom(1,is,current, (sup_group(current)))==rd_at_lis(group_iscat(sup_group(current))), &
!  sup_atom(2,is,current, (sup_group(current))), rd_at_lis(group_iscat(sup_group(current))), &
!  sup_atom(2,is,current, (sup_group(current)))==rd_at_lis(group_iscat(sup_group(current))), &
!  phase, sup_phase(1,is,current), arg(1)
!endif
         select case(sup_func(is, current))
            case(SUP_SINE)
               p_sub_function => sup_fun_sine
            case(SUP_CREN)
               p_sub_function => sup_fun_cren
            case default
               p_sub_function => sup_fun_sine
         end select
         if(sup_char(is, current)==SUP_DISP) then     ! Displacement wave
            do ii=1, 3
               vector = vector + sup_ampl(2:4,ii,is, current) * & !sin(arg(ii))
                                 p_sub_function(1.0_PREC_DP, 0.0_PREC_DP, arg(1), sup_func_p(is,current))
            enddo
         elseif(sup_char(is, current)==SUP_DENS) then
            prob = prob * p_sub_function(rep_ampl(1,is, current), rep_ampl(2,is,current), arg(1), sup_func_p(is,current))
!if(iatom<= 12.or. iatom > dim_natoms-12) then
!write(*,'(i5,i9,i5, 7f8.3, 2f20.16, l3, f8.3)') iatom, jatom, is, arg(1), cr_pos(:,jatom) , sup_qvec(:, current), &
!frac(arg(1)), frac( (frac(arg(1))+1.0D0)), frac( (frac(arg(1))+1.0D0))<sup_func_p(is,current), sup_func_p(is,current)
!write(*,'(i5,i9,i5, 7f8.3)') iatom, jatom, is, arg(1), v
!write(*, '(a,4f8.4, 2i8, l2)')  ' PROB ', prob, rep_ampl(1,is, current) , rep_ampl(2,is, current), &
!p_sub_function (rep_ampl(1,is, current), rep_ampl(2,is,current), arg(1), sup_func_p(is,current)), jatom, is, &
!sup_old(cr_iscat(1,jatom),is, current)
!endif
         endif
      enddo loop_current
!
      loop_current2: do current=1, sup_nwaves   ! Performe shift/ replacement
         if(sup_char(is, current)==SUP_DISP) then     ! Displacement wave
            cr_pos(:,jatom) = cr_pos(:,jatom) + vector
         elseif(sup_char(is, current)==SUP_DENS) then
            cond_old: if(sup_old(cr_iscat(1,jatom),is, current)) then   ! Proper old atom type
               call random_number(r1)
               if(r1<prob) then
                  call random_number(r1)
                  cr_iscat(1,jatom) = sup_new(int(r1)*sup_irepl(is, current),is, current)
               endif
            endif cond_old
         endif
      enddo loop_current2
   enddo loop_cell
enddo loop_main
!
deallocate(rep_ampl)
!
end subroutine super_run 
!
!*******************************************************************************
!
subroutine super_show
!-
! Show super space setting
!
use crystal_mod
use metric_mod
use superspace_mod
!
use param_mod
use prompt_mod
!
implicit none
!
character(len=PREC_STRING) :: string

character(len=12), dimension(1:2) :: cchar
character(len=12), dimension(1:4) :: cfunc
integer :: current  ! Current wave
integer :: is   ! Site number
integer :: iv   ! Displacement vector number
data cchar /'Displacement', &
            'Substitution'/
data cfunc /'sine        ', 'crenel      ', 'triangle    ', 'User        '/
!
write(output_io,*)
write(output_io, '(a )') ' Super space setting'
write(output_io, '(2a)') ' Super space structure  : ', sup_file(1:len_trim(sup_file))
loop_current: do current=1, sup_nwaves
   cond_show:if(any(sup_func(:, current)>0)) then
      write(string,'(3(f12.6:,'',''))') sup_qvec(:,current)
      is = len_trim(string)
      call d2r(string, is, .false.)
   write(output_io,*)
   write(output_io, '(a,i4, i4)')  ' Current modulation     : ', current, sup_group(current)
   write(output_io, '(a,2(3f10.5:, ''; ''))')  ' Satellite vector; Wave : ', sup_qvec(:,current), res_para(4:6)
do is = 1, cr_ncatoms
   cond_any:if(sup_func(is, current)>0) then
   if(sup_func(is, current)==SUP_SINE) then   ! Sine function
      write(output_io,'(a, i4,('', '',a12),'', '',a4)')  ' Site                   : ', &
            is, cchar(sup_char(is, current)), cfunc(sup_func(is, current))
   elseif(sup_func(is, current)==SUP_CREN) then   ! Crenel function
      write(output_io,'(a, i4,('', '',a12),'', '',a6,a6, 2f10.5)')  ' Site                   : ', &
            is, cchar(sup_char(is, current)), cfunc(sup_func(is, current)), ' up +-', &
            sup_func_p(is, current)!, 1.0_PREC_DP-sup_func_p(is, current)
      write(output_io,'(a, 2(3f10.5:, ''; ''))')    ' Up / down vectors      : ', res_para(4:6)*sup_func_p(is, current), &
            res_para(4:6)*(1.0_PREC_DP-sup_func_p(is, current))
   endif
   write(output_io,'(a, 4a4)') ' Pseudo atoms           : ', sup_atom(1:2, is, current, sup_group(current))
   if(sup_char(is, current)==0) then   ! Displacement wave
      do iv=1, 3
        if(sup_ampl(1,iv,is, current)>0.0) then
            write(output_io,'(a,i4, a, f10.5,a,3f10.5,a,f10.5)') &
            ' Ampl, Vector, Phase    : ', iv, ' ; ', sup_ampl(1,iv,is, current), ' ; ', &
            sup_ampl(2:4,iv,is, current), ' ; ', sup_phase (iv,is, current)
        endif
      enddo
   else                       ! Replacement = Density wave
      string = ' '
      do iv = 0, cr_nscat
         if(sup_old(iv,is, current)) then
             string = string(1:len_trim(string)) // cr_at_lis(iv)(1:4) //', '
         endif
      enddo
      iv = len_trim(string)
      if(string(iv:iv)==',') string(iv:iv)=' '
      write(output_io, '(a, a)') &
         ' Old atom types         : ', string(1:len_trim(string))
      string = ' '
      do iv = 0, sup_irepl(is, current)-1
         if(sup_new(iv,is, current)/=-1) then
         string = string(1:len_trim(string)) // cr_at_lis(sup_new(iv,is, current))(1:4) //', '
         endif
      enddo
      iv = len_trim(string)
      if(string(iv:iv)==',') string(iv:iv)=' '
      write(output_io, '(a, a)') &
         ' New atom types         : ', string(1:len_trim(string))
      write(output_io,'(a, 3f10.5)') &
         ' Low, high, Phase       : ', sup_prob(1,is, current), sup_prob(2,is, current), sup_phase(1,is, current)
   endif
   endif cond_any
enddo
   endif cond_show
enddo loop_current
!
end subroutine super_show
!
!*******************************************************************************
!
subroutine super_reset
!-
!  Reset all parameters
!+
use superspace_mod
use discus_allocate_appl_mod
!
integer :: is, ic, iw, ig
is = 1
ic = 1
iw = 1
ig = 1
!
call alloc_super( is, ic, iw, ig)
!
sup_nwaves = 1

end subroutine super_reset
!
!*******************************************************************************
!
end module superspace_func_mod
