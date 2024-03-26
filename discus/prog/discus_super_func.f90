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
                  CALL macro_close
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
!  Perform actual takso dummy atoms
!
   elseif(str_comp(befehl, 'run', 2, lbef, 3)) then
      call super_run
!
!  Perform actual takso dummy atoms
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
use str_comp_mod
use take_param_mod
!
implicit none
!
character(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp
!
integer, parameter :: MAXW  = 4
integer, parameter :: MAXWW = 3
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
real(kind=PREC_DP) :: vv ! Vectro length
integer :: ianz, iianz
integer :: i1, i2, j     ! Dummy indices
!integer :: isite         ! Current site
logical, dimension(1:10        ) :: l_isite
logical :: lselect
!
integer, parameter :: NOPTIONAL = 9
integer, parameter :: O_PHASE   = 1
integer, parameter :: O_SITE    = 2
integer, parameter :: O_TYPE    = 3
integer, parameter :: O_CHAR    = 4
integer, parameter :: O_FUNC    = 5
integer, parameter :: O_DISP    = 6
integer, parameter :: O_PROB    = 7
integer, parameter :: O_FILE    = 8
integer, parameter :: O_QVEC    = 9
character(LEN=  12), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 0 ! Number of values to calculate 
!
data oname  / 'phase'  , 'site' , 'type', 'character', 'function', 'displacement',    &
              'replace', 'file'     , 'qvec'                       /
data loname /  5,       4,       4,      9,           8,          12,               &
               7    ,  4         ,  4                           /
opara  =  (/ '0.0000        ', '0.0000        ', '0.0000        ', 'displa        ', &
             'sine          ', '[0,0,0,0]     ',                                     &
             '[VOID,0.0,1.0]', 'internal      ', '[1,0,0]       '          /)   ! Always provide fresh default values
lopara =  (/  6         ,  6,           6,           6,           6         ,  9         ,    &
              14        ,  8         ,  7                   /)
owerte =  (/  1.0, 1.0, 0.0,      0.0,      0.0    ,  0.0, 0.0, 0.0, 0.0      /)
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
if(lpresent(O_QVEC)) then   ! Set qvec:
   string=opara(O_QVEC)
   j = lopara(O_QVEC)
   call get_optional_multi(MAXWW, string, j, werte, iianz)
   if(ier_num/=0) return
   sup_qvec = werte(1:3)
   return
endif
!
cond_site: if(lpresent(O_SITE)) then
   if(cr_ncatoms>ubound(sup_char,1)) then
      call alloc_super(cr_ncatoms)
   endif
   if(opara(O_SITE)(1:1)=='[') then 
      i1 = 2
   else
      i1 = 1
   endif
   j = lopara(O_SITE)
   if(opara(O_SITE)(j:j)==']') then 
      i2 = j-1
   else
      i2 = j
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
!if(lpresent(O_TYPE)) then   ! assign type:   must always be present
!   where(l_isite)
!      sup_atom = opara(O_TYPE)(1:lopara(O_TYPE))
!   else where
!      sup_atom = ' '
!   end where
!
!  CHARACTER
!
   if(lpresent(O_CHAR)) then     ! "character:"
      if(str_comp(opara(O_CHAR), 'displacement', 3, lopara(O_CHAR), 12)) then
         i1 = SUP_DISP
      elseif(str_comp(opara(O_CHAR), 'density', 3, lopara(O_CHAR), 7)) then
         i1 = SUP_DENS
      endif
      sup_char = 0
      do j=1, cr_ncatoms
         if(l_isite(j)) then
            sup_char(j) = i1
         endif
      enddo
   endif
!
!  FUNCTION
!
   if(lpresent(O_FUNC)) then
      if(str_comp(opara(O_FUNC), 'sine', 3, lopara(O_FUNC), 4)) then
         i1 = SUP_SINE
      endif
      sup_func = 0
      do j=1, cr_ncatoms
         if(l_isite(j)) then
            sup_func(j) = i1
         endif
      enddo
   endif
!
!  DISPLACEMENT VECTOR
!
   if(lpresent(O_DISP)) then
      string = opara(O_DISP)
      j = lopara(O_DISP)
      call get_optional_multi(MAXW, string, j, werte, iianz)
      if(ier_num/=0) return
      do j=1,cr_ncatoms
         if(l_isite(j)) then
            sup_ampl(1:4,j) = werte(1:4)
            v = werte(2:4)
            vv = lib_blen(cr_gten, v)
            if(vv>0.0) then
               sup_ampl(2:4, j) = sup_ampl(2:4, j)*sup_ampl(1,j)/vv
            endif
         else
            sup_ampl(:, j) = 0.0_PREC_DP
         end if
      enddo
   endif
!
!  PHASE
!
   if(lpresent(O_PHASE)) then
      sup_phase = 0
      do j=1, cr_ncatoms
         if(l_isite(j)) then
            sup_phase(j) = owerte(O_PHASE)
         endif
      enddo
   endif
!
!  REPLACEMENT 
!  
   if(lpresent(O_PROB)) then
      string = opara(O_PROB)(2:lopara(O_PROB)-1)
      j = lopara(O_PROB) - 2
      call get_params(string, iianz, ccpara, llpara, MAXWW, j)
      at_new = ccpara(1)(1:llpara(1))
      i1     = len_trim(at_new)
      ccpara(1) = '0.0'
      llpara(1) = 3
      call ber_params(iianz, ccpara, llpara, wwerte, maxww)
      if(ier_num/=0) return
      lselect = .true.
      ccpara = ' '
      llpara = 0
      ccpara(1) = cr_at_lis(1)
      llpara(1) = len_trim(ccpara(1))
      ccpara(2) = at_new
      llpara(2) = len_trim(ccpara(2))
      iianz = 2
      call get_iscat(iianz, ccpara, llpara, werte, maxw, lselect)
!      CALL atom_select(at_new, i1, 0, cr_nscat, &
!                           wv_latom,  &
!                           wv_lsite, 1, cr_ncatoms,                     &
!                           lselect    ,        .true., .true.)
      loop_atom: do j=1, cr_ncatoms
         if(l_isite(j)) then
            sup_irepl(j) = nint(werte(2))
            sup_repl (j) = at_new
            sup_prob (1,j) = wwerte(2)
            sup_prob (2,j) = wwerte(3)
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
!
use precision_mod
use wink_mod
!
implicit none
!
!
integer  :: dim_natoms     ! Total number of atoms 
integer  :: dim_ncatoms    ! Number of atoms per uni cell
integer  :: dim_nscat      ! Number of chemical species
integer  :: dim_nanis      ! Number of ADPs
integer  :: dim_n_mole     ! Number of molecules
integer  :: dim_n_type     ! Number of molecule types
integer  :: dim_n_atom     ! Total number of atoms in molecules
integer  :: natoms     ! Total number of atoms 
integer  :: ncatoms    ! Number of atoms per uni cell
integer  :: nscat      ! Number of chemical species
integer  :: nanis      ! Number of ADPs
integer  :: n_mole     ! Number of molecules
integer  :: n_type     ! Number of molecule types
integer  :: n_atom     ! Total number of atoms in molecules
!
integer  :: iatom  ! Dummy atom number from "super file"
real(kind=PREC_DP) :: arg
real(kind=PREC_DP), dimension(3) :: vector
!
integer, dimension(3) :: ic
integer               :: is

!
REAL(kind=PREC_DP)  , DIMENSION(3)  :: rd_cr_pos
INTEGER             , dimension(3)  :: rd_cr_iscat
INTEGER                             :: rd_cr_prop
INTEGER             , DIMENSION(0:3):: rd_cr_surf
REAL(kind=PREC_DP)  , DIMENSION(0:3):: rd_cr_magn
INTEGER                             :: rd_cr_mole
INTEGER                             :: rd_cr_moleatom
!
real(kind=PREC_DP), dimension(:,:), allocatable :: rep_ampl
real(kind=PREC_DP)                              :: r1
real(kind=PREC_DP)                              :: prob
!
allocate(rep_ampl(2, cr_ncatoms))
do iatom=1, cr_ncatoms
   rep_ampl(1,iatom) = (sup_prob(2,iatom) + sup_prob(1,iatom))*0.5_PREC_DP
   rep_ampl(2,iatom) = (sup_prob(2,iatom) - sup_prob(1,iatom))*0.5_PREC_DP
enddo
!
call readstru_size_int(sup_file, dim_natoms, dim_ncatoms,      &
     dim_nscat, dim_nanis, dim_n_mole, dim_n_type, dim_n_atom, &
     natoms, ncatoms, nscat, nanis, n_mole, n_type, n_atom)
!
loop_main: do iatom=1, dim_natoms
   call indextocell(iatom, ic, is)
!
!  Read "pseudo atoms" from "super file"
!
   call struc_read_one_atom_internal(sup_file, iatom,  &
              rd_cr_pos, rd_cr_iscat, rd_cr_prop, rd_cr_surf, &
              rd_cr_magn, rd_cr_mole, rd_cr_moleatom )
   arg = (cr_pos(1,iatom)*sup_qvec(1) + cr_pos(2,iatom)*sup_qvec(2) +   &
          cr_pos(3,iatom)*sup_qvec(3) +                                 &
          sup_phase(is)+ real(0.5*mod(rd_cr_iscat(1)-1,2)) )*zpi
   if(sup_char(is)==SUP_DISP) then     ! Displacement wave
      vector = sup_ampl(2:4,is) * sin(arg)
      cr_pos(:,iatom) = cr_pos(:,iatom) + vector
   elseif(sup_char(is)==SUP_DENS) then
      prob = rep_ampl(1,is) + rep_ampl(2,is)*sin(arg)
      call random_number(r1)
      if(r1>prob) then
         cr_iscat(1,iatom) = sup_irepl(is)
      endif
   endif
enddo loop_main
!
deallocate(rep_ampl)
!
end subroutine super_run 
!
!*******************************************************************************
!
end module superspace_func_mod
