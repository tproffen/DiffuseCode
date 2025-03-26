MODULE discus_export
!
CONTAINS
!
SUBROUTINE do_export(line, lp)
!
! Exports a structure in different formats
!
use crystal_mod
use update_cr_dim_mod
!
use errlist_mod
use get_params_mod
use precision_mod
use take_param_mod
use str_comp_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT)  :: line
INTEGER         , INTENT(INOUT)  :: lp
!
INTEGER, PARAMETER :: MAXW = 5
CHARACTER(LEN=    PREC_STRING           ), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
INTEGER  :: ianz
integer  :: ncycle
logical  :: lforce
INTEGER, PARAMETER :: NOPTIONAL = 5
INTEGER, PARAMETER :: O_RMCVS   = 1
INTEGER, PARAMETER :: O_CYCLE   = 2
INTEGER, PARAMETER :: O_SPCGR   = 3
INTEGER, PARAMETER :: O_SITE    = 4
INTEGER, PARAMETER :: O_FORCE   = 5
CHARACTER(LEN=   7), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=    PREC_STRING           ), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 2 ! Number of values to calculate 
!
INTEGER  :: rmcversion
integer  :: scatty_site
!
!
DATA oname  / 'version', 'cycle  ', 'spcgr  ', 'site   ', 'force  ' /
DATA loname /  7       ,  5       ,  5       ,  4       ,  5        /
!
opara  =  (/ '6.00000', '20     ','P1     ', 'average', 'no     ' /)   ! Always provide fresh default values
lopara =  (/  7       ,  2       ,  7      ,  7       ,  2        /)
owerte =  (/  6.0D0   ,  20.0D0  ,  0.0D0  ,  0.0D0   , 0.0D0     /)
!
CALL get_params (line, ianz, cpara, lpara, MAXW, lp)
IF (ier_num.ne.0) THEN
   RETURN
ENDIF
if(cr_natoms<1) then
   ier_num = -198     ! Empty structure
   ier_typ = ER_APPL
   ier_msg(1) = 'No export performed as structure is emty'
   return
endif
!
CALL update_cr_dim
!                                                                       
IF (ianz.ge.1) THEN
   IF (str_comp (cpara (1) , 'cif', 2, lpara (1) , 3) ) THEN
      CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                        oname, loname, opara, lopara, lpresent, owerte)
      IF(ier_num/=0) RETURN
!
      IF (ianz >= 2) THEN
         CALL del_params (1, ianz, cpara, lpara, maxw)
         IF (ier_num.ne.0) RETURN
         ianz = ianz + 1
         cpara(ianz) = opara(O_SPCGR)
         lpara(ianz) = lopara(O_SPCGR)
         CALL discus2cif (ianz, cpara, lpara, MAXW)
      ELSE
         ier_num = - 6
         ier_typ = ER_COMM
      ENDIF
   ELSEIF (str_comp (cpara (1) , 'shelx', 2, lpara (1) , 5) ) THEN
      CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                        oname, loname, opara, lopara, lpresent, owerte)
      IF(ier_num/=0) RETURN
!
      IF (ianz >= 2) THEN
         CALL del_params (1, ianz, cpara, lpara, maxw)
         IF (ier_num.ne.0) RETURN
         ncycle = NINT(owerte(O_CYCLE))
         lforce = opara(O_FORCE)=='yes'
         CALL discus2ins (ianz, cpara, lpara, MAXW, ncycle, lforce)
      ELSE
         ier_num = - 6
         ier_typ = ER_COMM
      ENDIF
   ELSEIF (str_comp (cpara (1) , 'rmcprofile', 2, lpara (1) , 10) ) THEN
      CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                        oname, loname, opara, lopara, lpresent, owerte)
      IF(ier_num/=0) RETURN
!
      IF (ianz >= 2) THEN
         CALL del_params (1, ianz, cpara, lpara, maxw)
         IF (ier_num.ne.0) RETURN
         rmcversion = NINT(owerte(O_RMCVS))
         CALL discus2rmc6f (ianz, cpara, lpara, MAXW, rmcversion)
      ELSE
         ier_num = - 6
         ier_typ = ER_COMM
      ENDIF
   ELSEIF (str_comp (cpara (1) , 'scatty', 2, lpara (1) , 6) ) THEN
      CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                        oname, loname, opara, lopara, lpresent, owerte)
      IF(ier_num/=0) RETURN
!
      IF (ianz >= 2) THEN
         CALL del_params (1, ianz, cpara, lpara, maxw)
         IF (ier_num.ne.0) RETURN
         if(opara(O_SITE)(1:lopara(O_SITE))=='average') then
            scatty_site = 0
         elseif(opara(O_SITE)(1:lopara(O_SITE))=='indiv') then
            scatty_site = 1
         else
            ier_num = - 6
            ier_typ = ER_COMM
            ier_msg(1) = 'site parameter must be ''average'' or ''indiv'' '
            return
         endif
         CALL discus2scatty (ianz, cpara, lpara, MAXW, scatty_site)
      ELSE
         ier_num = - 6
         ier_typ = ER_COMM
      ENDIF
   ELSEIF (str_comp (cpara (1) , 'vasp', 2, lpara (1) , 4)  .OR.  &
           str_comp (cpara (1) , 'poscar', 2, lpara (1) , 6) ) THEN
      IF (ianz >= 2) THEN
         CALL del_params (1, ianz, cpara, lpara, maxw)
         IF (ier_num.ne.0) RETURN
         CALL discus2poscar (ianz, cpara, lpara, MAXW)
      ELSEIF(ianz==1) THEN
         ianz = 1
         cpara(1) = 'POSCAR'
         lpara(1) = 6
         CALL discus2poscar (ianz, cpara, lpara, MAXW)
      ELSE
         ier_num = - 6
         ier_typ = ER_COMM
      ENDIF
   ELSE
      ier_num = - 86
      ier_typ = ER_APPL
   ENDIF
ELSE
   ier_num = - 6
   ier_typ = ER_COMM
ENDIF
!
END SUBROUTINE do_export
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus2cif (ianz, cpara, lpara, MAXW)
!
USE crystal_mod
USE discus_allocate_appl_mod
USE discus_plot_mod
USE discus_plot_export_mod
USE modify_mod
USE molecule_mod
!
USE build_name_mod
USE errlist_mod
USE precision_mod
USE take_param_mod
USE support_mod
!
IMPLICIT NONE
!
INTEGER            ,                  INTENT(IN)    :: MAXW
INTEGER            ,                  INTENT(INOUT) :: ianz
CHARACTER(LEN=*   ), DIMENSION(MAXW), INTENT(INOUT) :: cpara
INTEGER            , DIMENSION(MAXW), INTENT(INOUT) :: lpara
!
INTEGER, PARAMETER :: IWR = 35
LOGICAL, PARAMETER :: lold = .FALSE.
!
CHARACTER(LEN=    PREC_STRING            )      :: ofile = ' '
CHARACTER(LEN=    PREC_STRING            )      :: zeile = ' '
CHARACTER(LEN=16)        :: do_spcgr = 'P1'
INTEGER                  :: lp
INTEGER                  :: nscat = 1
INTEGER                  :: nsite = 1
INTEGER                  :: nline = 1
REAL(KIND=PREC_DP)   , DIMENSION(MAXW) :: werte
!
do_spcgr = cpara(ianz)(1:LEN(do_spcgr))
ianz = ianz -1
!
CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1)
IF (ier_num.ne.0) THEN
   RETURN
ENDIF
ofile = cpara (1)
IF(ofile(lpara(1)-3:lpara(1)) /= '.cif') ofile = cpara (1) (1:lpara (1) ) //'.cif'
CALL oeffne (IWR, ofile, 'unknown')
IF (ier_num.ne.0) THEN
   CLOSE(IWR)
   RETURN
ENDIF
!
!
!     Allocate the necessary arrays
!
IF ( cr_nscat > PL_MAXSCAT .or. mole_num_type > PL_MAXSCAT .OR. & 
     MAXSCAT > PL_MAXSCAT ) THEN
   nscat = max ( cr_nscat, mole_num_type, MAXSCAT)
   nsite = max ( cr_ncatoms, PL_MAXSITE, MAXSCAT)
   nline =       PL_MAXLINE
   CALL alloc_plot ( nscat, nsite, nline )
   IF(ier_num < 0) THEN
     RETURN
   ENDIF
ENDIF
zeile = 'all'
lp    = 3
CALL atom_select (zeile, lp, 0, PL_MAXSCAT, pl_latom, &
                  pl_lsite, 0, PL_MAXSITE,            &
                  pl_sel_atom, lold, .TRUE.)
!
CALL plot_cif(IWR, .FALSE., do_spcgr)
!
CLOSE(IWR)
!
END SUBROUTINE discus2cif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus2ins (ianz, cpara, lpara, MAXW, ncycle, lforce)
!
USE chem_aver_mod
USE class_internal
USE crystal_mod
USE generate_mod
USE diffuse_mod
USE discus_save_mod
USE molecule_mod
!USE modify_mod
use prep_anis_mod
USE prop_para_func
USE prop_para_mod
use spcgr_apply,  only:setup_lattice, get_symmetry_matrices
USE structur
USE save_menu !, ONLY: save_internal
USE wyckoff_mod
!
USE build_name_mod
use charact_mod
USE errlist_mod
use lib_errlist_func
use element_data_mod
USE param_mod
USE precision_mod
USE support_mod
!
IMPLICIT NONE
!
INTEGER            ,                  INTENT(IN)    :: MAXW
INTEGER            ,                  INTENT(INOUT) :: ianz
CHARACTER(LEN=*   ), DIMENSION(MAXW), INTENT(INOUT) :: cpara
INTEGER            , DIMENSION(MAXW), INTENT(INOUT) :: lpara
integer            ,                  intent(in)    :: ncycle
logical            ,                  intent(in)    :: lforce
!
INTEGER, PARAMETER :: IWR = 35
integer, parameter :: MAXMASK = 4
real(kind=PREC_DP), parameter :: TOL = 0.001
!
CHARACTER(LEN=    PREC_STRING            )      :: ofile = ' '
CHARACTER(LEN=    PREC_STRING            )      :: line = ' '
CHARACTER(LEN=    PREC_STRING            )      :: names = ' '
CHARACTER(LEN=    PREC_STRING            )      :: units = ' '
CHARACTER(LEN=19)        :: origfile = ' '
CHARACTER(LEN= 5)        :: befehl   = 'lcell'
integer :: is_cond   ! return condition from element guess work
INTEGER                  :: i, j, k, l,i1, i2
INTEGER                  :: length
INTEGER                  :: lbef
INTEGER                  :: lattice = 1
INTEGER                  :: unique_n
INTEGER                  :: shelx_n
real(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: true_occ
real(kind=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: n_atoms 
character(len=4)                             :: atom_name
CHARACTER (LEN=2), DIMENSION(:), ALLOCATABLE :: unique_names
CHARACTER (LEN=4), DIMENSION(:), ALLOCATABLE :: shelx_names
integer           ,DIMENSION(:), ALLOCATABLE :: shelx_types 
real(kind=PREC_DP),DIMENSION(:), ALLOCATABLE :: unique_n_atoms 
logical                  :: lout = .FALSE.
logical                  :: l_not_full = .TRUE.
LOGICAL                  :: orig_OK =.FALSE.
logical, dimension(0:MAXMASK) :: uni_mask
REAL(KIND=PREC_DP)       :: z_unit
REAL(KIND=PREC_DP)   , DIMENSION(MAXW) :: werte
REAL(KIND=PREC_DP)   , DIMENSION(3), PARAMETER :: NULL = (/0.00, 0.00, 0.00/)
uni_mask(0)   = .true.
uni_mask(1:3) = .true.
uni_mask(4)   = .false.
!
if(.not.lforce .and. maxval(cr_icc)>1) then
  ier_num = -185
  ier_typ = ER_APPL
  ier_msg(1) = 'To save as SHELX instruction file the '
  ier_msg(2) = 'must consist of a single unit cell only!'
  return
endif
CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1)
IF (ier_num.ne.0) THEN
   RETURN
ENDIF
ofile = cpara (1)
IF(ofile(lpara(1)-3:lpara(1)) /= '.ins') ofile = cpara (1) (1:lpara (1) ) //'.ins'
CALL oeffne (IWR, ofile, 'unknown')
IF (ier_num.ne.0) THEN
   CLOSE(IWR)
   RETURN
ENDIF
!
! Test if centrosymmetric and origin at 0,0,0. 
! Generator no. 6 must be present and be the last operator except for 
! centering operators
!
IF(.NOT.cr_acentric) THEN                      ! Is centrosymmetric
   orig_OK = .FALSE.
   is_null: DO j=1,generspcgr(0,cr_spcgrno)
      IF(generspcgr(j,cr_spcgrno) == 6 ) THEN
         orig_OK = .TRUE.                      ! Found generator No. 6
         EXIT is_null
      ENDIF
   ENDDO is_null
   IF(.NOT. orig_OK) THEN
      ier_num = -148
      ier_typ = ER_APPL
      ier_msg(1) = 'Spacegroup is centrosymmetric but the center'
      ier_msg(2) = 'of symmetry is not at the origin. '
      ier_msg(3) = 'Transform to setting with 1bar at origin'
      CLOSE(IWR)
      RETURN
   ENDIF
ELSE
   orig_OK = .TRUE.   ! Acentric space group origin is always OK
ENDIF
IF(cr_spcgr_set(1:1)=='P') lattice = 1    ! Requires the alternative setting 
IF(cr_spcgr_set(1:1)=='I') lattice = 2    ! to get the correct centering in the
                                          ! orthorhombic space groups
IF(cr_spcgr_set(1:1)=='R') then
   if(abs(cr_win(3)-120.0_PREC_DP)<TOL) then
      lattice = 3    ! rhombohedral space group, hexagonal setting
   else
      lattice = 1    ! rhombohedral space group, rhombohedral setting
   endif
endif
IF(cr_spcgr_set(1:1)=='F') lattice = 4
IF(cr_spcgr_set(1:1)=='A') lattice = 5
IF(cr_spcgr_set(1:1)=='B') lattice = 6
IF(cr_spcgr_set(1:1)=='C') lattice = 7
!
! Temporarily save the original structure
!
call save_store_setting
sav_w_gene = .true.
sav_w_ncell= .true.
sav_w_symm = .true.
sav_w_scat = .true.
sav_w_adp  = .true.
sav_w_occ  = .true.
sav_w_surf = .true.
sav_w_magn = .FALSE.   ! MAGNETIC_WORK
sav_w_mole = .true.
sav_w_obje = .true.
sav_w_doma = .true.
sav_w_prop = .true.
origfile  = 'internal.shelx_orig'
line       = 'ignore, all'          ! Ignore all properties
length     = 11
CALL property_select(line, length, sav_sel_prop)
CALL save_internal(origfile)
ALLOCATE(n_atoms(2,1:cr_nscat))
allocate(true_occ(1:ubound(cr_occ,1)))
true_occ(1:cr_nscat) = cr_occ(1:cr_nscat)
CALL chem_elem(.FALSE.)
DO i=1, cr_nscat
   n_atoms(1,i) = (res_para(i+1)*cr_natoms)*cr_occ(i)
ENDDO
!
! Read structure as unit cell file and expand to full crystal
!
befehl = 'lcell'
lbef   = 5
cpara(1) = origfile
lpara(1) = 19
CALL do_readcell(befehl, lbef, ianz, MAXW, cpara, lpara, .TRUE., 1.0D-5, 0, .TRUE., .FALSE., MAXMASK, uni_mask)
CALL chem_elem(.FALSE.)
z_unit = 192.0
DO i=1, cr_nscat
   n_atoms(2,i) = (res_para(i+1)*cr_natoms*true_occ(i))
   IF(n_atoms(1,i)>0 .AND. n_atoms(2,i)>0) THEN
      z_unit = MIN(z_unit, REAL(n_atoms(2,i)/n_atoms(1,i), kind=PREC_DP))
   ENDIF
ENDDO
deallocate(true_occ)
!
! Restore original structure
!
CALL do_readstru(MAXMASK, origfile, .FALSE., uni_mask, l_not_full)
lout = .false.
CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps, &
     cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat,       &
     cr_fmat, cr_cartesian,                                      &
     cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
CALL get_symmetry_matrices 
l_not_full = .false.
call prep_anis(cr_natoms, l_not_full)
!
call save_restore_setting
CALL store_remove_single(origfile, ier_num)
IF(ier_num/=0) THEN
   call save_restore_setting
   ier_typ = ER_APPL
   ier_msg(1) = 'Could not remove temporary internal storage'
   ier_msg(2) = 'in shelx export'
   ier_msg(3) = 'Please document and report'
   RETURN
ENDIF
!
! Reduce names to unique chemical element names
!
ALLOCATE(unique_n_atoms(1:cr_nscat))
ALLOCATE(unique_names(1:cr_nscat))
unique_n = 1
!atom_name= cr_at_lis(1)(1:2)
!i1 = iachar(  atom_name(2:2))
!if(.not. ( (a<=i1 .and. i1<=z) .or. (aa<=i1 .and. i1<=zz)))   atom_name(2:2) = ' '
call guess_element(atom_name, is_cond, cr_at_lis(1), scat_equ=cr_scat_equ(1),   &
                                                     scat_equ_name=cr_at_equ(1))
unique_names(1)   = atom_name(1:2)       ! Chemical element name only
unique_n_atoms(1) = n_atoms(2,1)
loop_types: DO i=2, cr_nscat
   DO j=1,unique_n
!     atom_name= cr_at_lis(i)(1:2)
!     i1 = iachar(  atom_name(2:2))
!     if(.not. ( (a<=i1 .and. i1<=z) .or. (aa<=i1 .and. i1<=zz)))   atom_name(2:2) = ' '
      call guess_element(atom_name, is_cond, cr_at_lis(i), scat_equ=cr_scat_equ(i),   &
                                                           scat_equ_name=cr_at_equ(i))
      if(unique_names(j) == atom_name(1:2)) then  ! Previous name
         unique_n_atoms(j) = unique_n_atoms(j) + n_atoms(2,i)
         CYCLE loop_types
      ENDIF
   ENDDO
!  atom_name= cr_at_lis(i)(1:2)
!  i1 = iachar(  atom_name(2:2))
!  if(.not. ( (a<=i1 .and. i1<=z) .or. (aa<=i1 .and. i1<=zz)))   atom_name(2:2) = ' '
   call guess_element(atom_name, is_cond, cr_at_lis(i), scat_equ=cr_scat_equ(i),   &
                                                        scat_equ_name=cr_at_equ(i))
   unique_n = unique_n + 1
   unique_names(unique_n) = atom_name(1:2)
   unique_n_atoms(unique_n) = n_atoms(2,i)
ENDDO loop_types
!
names = ' '
units = ' '
DO i=1,unique_n                                ! Prepare SFAC / UNIT statement
   i1 = (i-1) * 3 + 1
   i2 = (i  ) * 3
   names(i1:i2) = unique_names(i)(1:2)
   i1 = (i-1) * 9 + 1
   i2 = (i  ) * 9
   WRITE(units(i1:i2),'(1x,f8.2)') unique_n_atoms(i)
ENDDO
!
! Create unique Shelx names
!
ALLOCATE(shelx_names(1:cr_natoms))
allocate(shelx_types(1:cr_natoms))
shelx_names(:) = ' '
shelx_types(:) = 0
shelx_n = 0
!!loop_shelx_o: DO i=1,cr_natoms
!!   l = shelx_n
!!   DO j=l,1,-1
!!      atom_name= cr_at_lis(cr_iscat(1,i))(1:2)
!!      i1 = iachar(  atom_name(2:2))
!!      if(.not. ( (a<=i1 .and. i1<=z) .or. (aa<=i1 .and. i1<=zz)))   atom_name(2:2) = ' '
!!!     IF(cr_at_lis(cr_iscat(1,i))(1:2) == shelx_names(j)(1:2)) THEN
!!      IF(atom_name               (1:2) == shelx_names(j)(1:2)) THEN
!!         shelx_n = shelx_n + 1
!!!        shelx_names(i)(1:2) = cr_at_lis(cr_iscat(1,i))(1:2)
!!         shelx_names(i)(1:2) = atom_name(1:2)
!!         shelx_names(i)(3:3) = shelx_names(j)(3:3)
 !!        i1 = iachar(shelx_names(i)(2:2))
!!         if(.not. ( (a<=i1 .and. i1<=z) .or. (aa<=i1 .and. i1<=zz)))  shelx_names(i)(2:2) = '0'
!!         shelx_types(shelx_n)      = cr_iscat(1,i)
!!         read(shelx_names(j)(3:4),'(i2)') k
!!         write(shelx_names(i)(3:4),'(i2.2)') k + 1
!!         CYCLE loop_shelx_o
!!      ENDIF
!!   ENDDO
!!   shelx_n = shelx_n + 1
!!   shelx_names(shelx_n)(1:2) = cr_at_lis(cr_iscat(1,i))(1:2)
!!   shelx_names(shelx_n)(3:4) = '01'
!!   i1 = iachar(shelx_names(shelx_n)(2:2))
!!   if(.not. ( (a<=i1 .and. i1<=z) .or. (aa<=i1 .and. i1<=zz)))  shelx_names(shelx_n)(2:2) = '0'
!!   shelx_types(shelx_n)      = cr_iscat(1,i)
!!ENDDO loop_shelx_o
shelx_names(:) = ' '
shelx_types(:) = 0
shelx_n = 0
loop_shelx: do i=1,cr_natoms
   l = shelx_n
   do j=l,1,-1                              ! Count previous names backwards
      atom_name = cr_at_lis(cr_iscat(1,i))
      if(atom_name == shelx_names(j)) then  ! Found previous name, augment a new name
         if(len_trim(atom_name)< 4) then    ! Previous atom name is short enough
            shelx_n = shelx_n + 1
            shelx_names(shelx_n) = atom_name   ! New atom name is identical but:
            shelx_types(shelx_n) = cr_iscat(1,i)
            read(shelx_names(j)(4:4),'(i1)') k
            write(shelx_names(shelx_n)(4:4),'(i1)') k+1 ! We increment a number at digit 4
            cycle loop_shelx
         else                               ! Error previous atom name is 4 characters as well!
            ier_num = 10
            ier_typ = ER_APPL
            ier_msg(1) = 'Atom ' //atom_name//' occurs multiple times in unit cell'
            ier_msg(2) = 'Edit SHELXL instruction file to avoid error'
            call errlist
         endif
      endif
   enddo
   shelx_n = shelx_n + 1
   shelx_names(shelx_n) = cr_at_lis(cr_iscat(1,i))
   shelx_types(shelx_n)      = cr_iscat(1,i)
enddo loop_shelx
!
DO i=1,shelx_n
   DO j=1,4
      IF(shelx_names(i)(j:j) == ' ' .and. len_trim(shelx_names(i))>j) shelx_names(i)(j:j) = '0'
   ENDDO
ENDDO
!
!  Start writing the actual file 
!
IF(cr_syst==4 .AND. cr_iset/=1) THEN
   WRITE(IWR, 1010) cr_name(1:LEN_TRIM(cr_name)), cr_spcgr(1:LEN_TRIM(cr_spcgr)), cr_set, cr_spcgr_set
ELSE
   WRITE(IWR, 1000) cr_name(1:LEN_TRIM(cr_name)), cr_spcgr(1:LEN_TRIM(cr_spcgr))
ENDIF
WRITE(IWR, 1100) rlambda, cr_a0(:), cr_win(:)
WRITE(IWR, 1200) z_unit, NULL, NULL
IF(cr_acentric) THEN
   WRITE(IWR, 1310) -1*IABS(lattice)
ELSE
   WRITE(IWR, 1300) IABS(lattice)
ENDIF
!
! Write Symmetry matrices
j = spc_n
!
IF(cr_spcgr(1:1)=='P') THEN
   j = spc_n                              ! Need all symmetry operations
ELSEIF(cr_spcgr(1:1)=='A' .OR. cr_spcgr(1:1)=='B' .OR.         &
       cr_spcgr(1:1)=='C' .OR. cr_spcgr(1:1)=='I'      ) THEN
   j =spc_n / 2                           ! Only need the first half symmetry operations
ELSEIF(cr_spcgr(1:1)=='F' ) THEN
   j =spc_n / 4                           ! Only need the first quarter symmetry operations
ELSEIF(cr_spcgr(1:1)=='R' ) THEN
   if(abs(cr_win(3)-120.0_PREC_DP)<TOL) then
   j =spc_n / 3                           ! Only need the first third symmetry operations
   else
      j = spc_n
   endif
ENDIF
IF(.NOT.cr_acentric) j = j / 2            ! Centrosymmetric need the first half only
DO i=2, j                                 ! Omit identity x,y,z
   WRITE(IWR, 1400) spc_xyz(i)(1:LEN_TRIM(spc_xyz(i)))
ENDDO
!
WRITE(IWR, 1500) names(1:LEN_TRIM(names)) ! SFAC
do i=1, unique_n
   write(IWR,'(a6,a2,1x,3f9.5)') 'DISP $',unique_names(i)(1:2), &
   cr_delfr(shelx_types(i)), cr_delfi(shelx_types(i)), 0.0_PREC_DP
enddo
WRITE(IWR, 1600) units(1:LEN_TRIM(units)) ! UNIT
WRITE(IWR, 1700) ncycle      ! L.S.
WRITE(IWR, 1800)             ! BOND
WRITE(IWR, 1900)             ! BOND $H
WRITE(IWR, 2000)             ! FMAP 2
WRITE(IWR, 2100)             ! PLAN 20
WRITE(IWR, 2200)             ! ACTA
WRITE(IWR, '(a)') 'LIST 6'   ! LIST 6
WRITE(IWR, 2300)             ! WGHT
if(diff_exti > 0.0_PREC_DP) write(IWR, '(a4,1x,f9.5)') 'EXTI',diff_exti ! EXTI
WRITE(IWR, 2400)             ! FVAR
!
!  First write atoms that are not part of a molecule
!
DO i=1, cr_natoms
   IF(cr_mole(i)==0) THEN
      CALL shelx_write_atom(IWR,i, unique_n, unique_names, cr_natoms, shelx_names)
   ENDIF
ENDDO
!
!  Now write any molecule
!
DO i=1,mole_num_mole
   WRITE(IWR,2450) i
   DO j = 1, mole_len (i)
      k = mole_cont (mole_off (i) + j)
      CALL shelx_write_atom(IWR,k, unique_n, unique_names, cr_natoms, shelx_names)
   ENDDO
ENDDO
WRITE(IWR, 2600)             ! HKLF 4
WRITE(IWR, 2700)             ! END 
!
DEALLOCATE(n_atoms)
DEALLOCATE(unique_n_atoms)
DEALLOCATE(unique_names)
DEALLOCATE(shelx_names)
DEALLOCATE(shelx_types)
CLOSE(IWR)
!
!
1010 FORMAT('TITL ',a,' in ',a,' setting: ',a3,' == ',a16)
1000 FORMAT('TITL ',a,' in ',a)
1100 FORMAT('CELL ',f7.5,3f11.6,3f9.4)
1200 FORMAT('ZERR ',f7.2,3f11.6,3f9.4)
1300 FORMAT('LATT ', i1)
1310 FORMAT('LATT ', i2)
1400 FORMAT('SYMM ', a )
1500 FORMAT('SFAC ', a )
1600 FORMAT('UNIT ', a )
1700 FORMAT('L.S. ', i6)
1800 FORMAT('BOND '    )
1900 FORMAT('BOND $H ' )
2000 FORMAT('FMAP 2'   )
2100 FORMAT('PLAN 20'  )
2200 FORMAT('ACTA '    )
2300 FORMAT('WGHT    0.000000')
2400 FORMAT('FVAR 1.00000')
2450 FORMAT('MOLE ',i3 )
2600 FORMAT('HKLF 4'   )
2700 FORMAT('END   '   )
!
!
END SUBROUTINE discus2ins
!
!*******************************************************************************
!
SUBROUTINE shelx_write_atom(IWR, i, unique_n, unique_names, natoms, shelx_names)
!
USE crystal_mod
USE param_mod
USE spcgr_apply, ONLY: get_wyckoff
!
use charact_mod
USE wink_mod
!
IMPLICIT NONE
!
INTEGER                              , INTENT(IN) :: IWR
INTEGER                              , INTENT(IN) :: i
INTEGER                              , INTENT(IN) :: unique_n
CHARACTER(LEN=*), DIMENSION(unique_n), INTENT(IN) :: unique_names
INTEGER                              , INTENT(IN) :: natoms
CHARACTER(LEN=*), DIMENSION(natoms  ), INTENT(IN) :: shelx_names
!
character(len=4) :: atom_name
real(kind=PREC_DP), parameter :: TOL = 2.0D-5
INTEGER :: stype
INTEGER :: j, i1
REAL(KIND=PREC_DP) :: occup,biso
REAL(KIND=PREC_DP), DIMENSION(3) :: vec
!
   stype = 0
   loop_stype: DO j=1,unique_n
      if(cr_scat_equ(cr_iscat(1,i))) then
         if(cr_at_equ(cr_iscat(1,i))(1:2) == unique_names(j)(1:2)) then
         stype = j
         exit loop_stype
         endif
      else
      atom_name = cr_at_lis(cr_iscat(1,i))(1:2)
      i1 = iachar(  atom_name(2:2))
      if(.not. ( (a<=i1 .and. i1<=z) .or. (aa<=i1 .and. i1<=zz)))   atom_name(2:2) = ' '
      IF(atom_name               (1:2) == unique_names(j)(1:2))  THEN
!     IF(cr_at_lis(cr_iscat(1,i))(1:2) == unique_names(j)(1:2))  THEN
         stype = j
         EXIT loop_stype
      ENDIF
      endif
   ENDDO loop_stype
   vec(:) = cr_pos(:,i)
   CALL get_wyckoff(vec,.FALSE.,1)
   occup = 10.000 + REAL(res_para(1)/res_para(3))*cr_occ(cr_iscat(1,i))
   biso = cr_dw(cr_iscat(1,i))/8./REAL(pi**2)
if(abs(cr_prin(4,1,cr_iscat(3,i))-cr_prin(4,2,cr_iscat(3,i)))>TOL  .or.  &
   abs(cr_prin(4,1,cr_iscat(3,i))-cr_prin(4,3,cr_iscat(3,i)))>TOL      ) then
   write(IWR, 2510) shelx_names(i)(1:4), stype, cr_pos(:,i), occup, cr_anis_full(1:2,cr_iscat(3,i))
   write(IWR, 2511) cr_anis_full(3:6,cr_iscat(3,i))
else
   WRITE(IWR,2500) shelx_names(i), stype, cr_pos(:,i), occup,biso
endif
!
2500 FORMAT(a4,1x,i2,3(f12.6),f12.5,f11.5)
2510 FORMAT(a4,1x,i2,3(f12.6),f12.5,2f11.5,' =')
2511 format(5x,4f11.5)
!
END SUBROUTINE shelx_write_atom
!
!*******************************************************************************
!
SUBROUTINE discus2rmc6f (ianz, cpara, lpara, MAXW, rmcversion)
!
!  Exports the current structure in RMC6F format
!
USE crystal_mod
use chem_mod
USE celltoindex_mod
USE modify_func_mod
!
USE charact_mod
USE build_name_mod
USE errlist_mod
USE param_mod
USE precision_mod
USE prompt_mod
USE times_mod
USE trig_degree_mod
USE support_mod
!
IMPLICIT NONE
!
INTEGER            ,                  INTENT(IN)    :: MAXW
INTEGER            ,                  INTENT(INOUT) :: ianz
CHARACTER(LEN=*   ), DIMENSION(MAXW), INTENT(INOUT) :: cpara
INTEGER            , DIMENSION(MAXW), INTENT(INOUT) :: lpara
INTEGER                             , INTENT(IN)    :: rmcversion     ! Version is 6 or 7
!
INTEGER, PARAMETER :: IWR = 35
!
CHARACTER(LEN=    PREC_STRING            )      :: ofile  = ' '
CHARACTER(LEN=    PREC_STRING            )      :: line   = ' '
CHARACTER(LEN=    PREC_STRING            )      :: string = ' '
CHARACTER(LEN=9), DIMENSION(:), ALLOCATABLE :: atom_names
INTEGER         , DIMENSION(:), ALLOCATABLE :: atom_number
INTEGER                  :: iatom, iscat, is, i, j, ii, ll
INTEGER                  :: ntypes      ! Actual atom types      to be written to file
INTEGER                  :: natoms      ! Actual number of atoms to be written to file
INTEGER                  :: nrealatoms      ! Actual number of real atoms to be written to file
INTEGER, DIMENSION(3)    :: icell 
integer, dimension(:,:,:), allocatable :: isites
integer              , DIMENSION(3)    :: scalef
REAL(KIND=PREC_DP)   , DIMENSION(3)    :: shift
REAL(KIND=PREC_DP)   , DIMENSION(3)    :: posit    ! Position to write, restricted to [0,1[
REAL(KIND=PREC_DP)   , DIMENSION(MAXW) :: werte
!
! Build the output file name
!
CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1)
IF (ier_num.ne.0) THEN
   RETURN
ENDIF
ofile = cpara (1)
if(rmcversion==7) then
   IF(.NOT.(ofile(lpara(1)-4:lpara(1)) == '.rmc7')) THEN
      ofile = cpara (1) (1:lpara (1) ) //'.rmc7'
   ENDIF
elseif(rmcversion==6) then
   IF(.NOT.(ofile(lpara(1)-5:lpara(1)) == '.rmc6f')) THEN
      ofile = cpara (1) (1:lpara (1) ) //'.rmc6f'
   ENDIF
endif
CALL oeffne (IWR, ofile, 'unknown')
IF (ier_num.ne.0) THEN
   CLOSE(IWR)
   RETURN
ENDIF
!
! Make RMCprofile atom names
!
ALLOCATE(atom_names(0:MAXSCAT))
ALLOCATE(atom_number(0:MAXSCAT))
atom_names(:)  = ' '
atom_number(:) = 0
IF(rmcversion==6) THEN
   DO i=0,MAXSCAT
      atom_names(i)(1:2) = cr_at_lis(i)(1:2)
      IF(IACHAR(atom_names(i)(2:2))>=aa .AND. IACHAR(atom_names(i)(2:2))<=zz) THEN
         atom_names(i)(2:2) = ACHAR(IACHAR(atom_names(i)(2:2)) - aa + a)
      ELSEIF(.NOT. (IACHAR(atom_names(i)(2:2))>=a .AND. IACHAR(atom_names(i)(2:2))<=z)) THEN
         atom_names(i)(2:2) = ' '
      ENDIF
      IF(cr_at_lis(i)(1:2)=='H.') atom_names(i)='H   '
      IF(cr_at_lis(i)(1:2)=='D ') atom_names(i)='H   '
      IF(cr_at_lis(i)(1:4)=='VOID') atom_names(i)(1:4)='Va  '
   ENDDO
ELSE
   DO i=0,MAXSCAT
      atom_names(i)(1:4) = cr_at_lis(i)(1:4)
      IF(IACHAR(atom_names(i)(2:2))>=aa .AND. IACHAR(atom_names(i)(2:2))<=zz) THEN
         atom_names(i)(2:2) = ACHAR(IACHAR(atom_names(i)(2:2)) - aa + a)
      ENDIF
      IF(cr_at_lis(i)(1:2)=='H.' ) atom_names(i)(1:4)='H   '
      IF(cr_at_lis(i)(1:2)=='D ' ) atom_names(i)(1:4)='2H  '
      IF(cr_at_lis(i)(1:3)=='D1-') atom_names(i)(1:4)='2H  '
      IF(cr_at_lis(i)(1:4)=='VOID') atom_names(i)(1:4)='Va  '
      ii = 1
      DO j=0,i-1
         IF(atom_names(i)(1:4) == atom_names(j)(1:4)) THEN
            ii = ii + 1
         ENDIF
      ENDDO
      IF(ii<10) THEN
         WRITE(atom_names(i)(5:9),'(A2,I1,A2)') ' [',ii,'] '
      ELSEIF(ii<100) THEN
         WRITE(atom_names(i)(5:9),'(A2,I2,A1)') ' [',ii,']'
      ENDIF
   ENDDO
ENDIF
!
! Count atoms that will be written to file
!
natoms = 0
nrealatoms = 0
ntypes = 0
count: DO iatom=1,cr_natoms
   IF (check_select_status (iatom, .true., cr_prop (iatom),  cr_sel_prop) ) THEN
      natoms = natoms + 1
      IF(cr_at_lis(cr_iscat(1,iatom))/='VOID') THEN
            nrealatoms = nrealatoms + 1
      endif
      find_entry: DO j=0,MAXSCAT
         IF(atom_names(cr_iscat(1,iatom))(1:4)==atom_names(j)(1:4)) THEN
            atom_number(j) = atom_number(j) + 1
            EXIT find_entry
         ENDIF
      ENDDO find_entry
   ENDIF
ENDDO count
DO i=0,MAXSCAT
   IF(atom_number(i) > 0 ) ntypes = ntypes + 1
ENDDO
!write(*,'(a,3f8.4)') ' DIMENSIONS x ', cr_dim(1,1), cr_dim(1,2), cr_dim(1,2)-cr_dim(1,1)
!write(*,'(a,3f8.4)') ' DIMENSIONS y ', cr_dim(2,1), cr_dim(2,2), cr_dim(2,2)-cr_dim(2,1)
!write(*,'(a,3f8.4)') ' DIMENSIONS z ', cr_dim(3,1), cr_dim(3,2), cr_dim(3,2)-cr_dim(3,1)
!write(*,'(a,3i8  )') ' CR_ICC       ', cr_icc
DO i=1,3
   IF(NINT(cr_dim(i,2)-cr_dim(i,1))-(cr_dim(i,2)-cr_dim(i,1))== 0.000D0) THEN
      scalef(i) = MAX(1,NINT((cr_dim(i,2)-cr_dim(i,1)))+1)
   ELSE
      if(cr_dim(i,2)-cr_dim(i,1) - int(cr_dim(i,2)-cr_dim(i,1))>0.50D0) then
         scalef(i) =  INT((cr_dim(i,2)-cr_dim(i,1))) + 1
      else
         scalef(i) =  INT((cr_dim(i,2)-cr_dim(i,1)))
      endif
   ENDIF
ENDDO
!
! Write the actual file
!
CALL datum
IF(rmcversion==6) THEN
   WRITE(IWR, '(a)')   '(Version 6f format configuration file)'
else
   WRITE(IWR, '(a)')   '(Version 7 format configuration file)'
endif
WRITE(IWR, '(a,a,a)') '(Generated by DISCUS Version ',version,')'
WRITE(IWR, '(a)')   'Metadata owner:     DISCUS'
WRITE(IWR, 1100 )    int_date(3), int_date(2), int_date(1)
WRITE(IWR, '(a,a)') 'Metadata comment:   ',cr_name(1:LEN_TRIM(cr_name))
WRITE(IWR, 1110 )   ntypes
line   = 'Atom types present:         '
string = 'Number of each atom type:   '
ii     = 28
ll     = 28
DO i=0, cr_nscat
   IF(atom_number(i) > 0) THEN
      line = line(1:ii) // ' ' // atom_names(i)(1:4)
      ii = LEN_TRIM(line)
      WRITE(string(ll+1:ll+13),'(I12,1x)') atom_number(i)
      ll = ll + 13
   ENDIF
ENDDO
WRITE(IWR, '(a)')   line(1:LEN_TRIM(line))
WRITE(IWR, '(a)')   string(1:LEN_TRIM(string))
WRITE(IWR, '(a)')   'Number of moves generated:                  0'
WRITE(IWR, '(a)')   'Number of moves tried:                      0'
WRITE(IWR, '(a)')   'Number of moves accepted:                   0'
WRITE(IWR, '(a)')   'Number of prior configuration saves: 0'
WRITE(IWR, 1200 )   natoms
!WRITE(IWR, 1300 )   REAL(nrealatoms, kind=PREC_DP)/(REAL(cr_icc(1)*cr_icc(2)*cr_icc(3), kind=PREC_DP)*cr_v)
WRITE(IWR, 1300 )   REAL(nrealatoms, kind=PREC_DP)/(REAL(scalef(1)*scalef(2)*scalef(3), kind=PREC_DP)*cr_v)
WRITE(IWR, 1400 )   scalef(:)
!WRITE(IWR, 1400 )   cr_icc(:)
!WRITE(IWR, 1500 )   (cr_a0 (i)*REAL(cr_icc(i)),i=1,3), cr_win(:)
WRITE(IWR, 1500 )   (cr_a0 (i)*REAL(scalef(i), kind=PREC_DP),i=1,3), cr_win(:)
WRITE(IWR, '(a)')   'Lattice vectors (Ang):'
IF(cr_syst==cr_ortho .OR. cr_syst==cr_tetragonal .OR. cr_syst==cr_cubic) THEN
   WRITE(IWR, 1600 )   cr_a0(1)*REAL(scalef(1),kind=PREC_DP), 0.000000D0,                            0.000000D0
   WRITE(IWR, 1600 )   0.000000D0,                            cr_a0(2)*REAL(scalef(2),kind=PREC_DP), 0.000000D0
   WRITE(IWR, 1600 )   0.000000D0,                            0.000000D0,                            cr_a0(3)*REAL(scalef(3),kind=PREC_DP)
ELSE
!  WRITE(IWR, 1600 )   cr_a0(1)*REAL(cr_icc(1)), 0.000000, 0.00000
!  WRITE(IWR, 1600 )   cr_a0(2)*REAL(cr_icc(2))*cosd(cr_win(3))                     , &
!                      cr_a0(2)*REAL(cr_icc(2))*sind(cr_win(3)), 0.000000
!  WRITE(IWR, 1600 )   cr_a0(3)*REAL(cr_icc(3))*cosd(cr_win(2))                     , &
!                      cr_a0(3)*REAL(cr_icc(3))*sind(cr_win(2))*cosd(cr_win(1))     , &
!                      cr_v    *REAL(cr_icc(3))/(cr_a0(1)*cr_a0(2)*sind(cr_win(3)))
   WRITE(IWR, 1600 )   cr_a0(1)*REAL(scalef(1)), 0.000000D0, 0.00000D0
   WRITE(IWR, 1600 )   cr_a0(2)*REAL(scalef(2), KIND=PREC_DP)*cosd(cr_win(3))                     , &
                       cr_a0(2)*REAL(scalef(2), KIND=PREC_DP)*sind(cr_win(3)), 0.000000D0
   WRITE(IWR, 1600 )   cr_a0(3)*REAL(scalef(3), KIND=PREC_DP)*cosd(cr_win(2))                     , &
                       cr_a0(3)*REAL(scalef(3), KIND=PREC_DP)*sind(cr_win(2))*cosd(cr_win(1))     , &
                       cr_v    *REAL(scalef(3), KIND=PREC_DP)/(cr_a0(1)*cr_a0(2)*sind(cr_win(3)))
ENDIF
WRITE(IWR, '(a)')   'Atoms:'
!
! Finished header, loop over all atom types and all atoms
! Property induced selection rules hold, VOIDS are never written
!
if(.not.chem_quick) then                                 ! Nonperiodic strucure 
   allocate(isites(scalef(1),scalef(2),scalef(3)))       ! we need to set the sites
   isites = 0
endif
shift(:) = REAL(NINT(cr_dim(:,1)))
!shift = cr_dim(:,1)
ii       = 0
DO iscat=0,MAXSCAT
   DO iatom=1,cr_natoms
      IF(cr_iscat(1,iatom)==iscat) THEN
!        IF(cr_at_lis(cr_iscat(1,iatom))/='VOID') THEN
            IF (check_select_status (iatom, .true., cr_prop (iatom),  cr_sel_prop) ) THEN
               ii = ii + 1                               ! increment atom number for RMCprofile sequence
               if(chem_quick) then                       ! Structure is periodic
                  CALL indextocell (iatom, icell, is)
               else                                      ! Nonperiodic cell, use coordinates
                  icell(:) = int(cr_pos(:, iatom)-cr_dim(:,1)) + 1 
                  if(icell(1)>scalef(1)) icell(1) = icell(1)-scalef(1)
                  if(icell(2)>scalef(2)) icell(2) = icell(2)-scalef(2)
                  if(icell(3)>scalef(3)) icell(3) = icell(3)-scalef(3)
                  isites(icell(1),icell(2),icell(3)) = isites(icell(1),icell(2),icell(3)) + 1
                  is = isites(icell(1),icell(2),icell(3))
               endif
               posit = (cr_pos(:,iatom)-shift(:))/REAL(scalef(:), kind=PREC_DP)
               if(posit(1)>=1.0D0) posit(1) = posit(1)-1.0D0
               if(posit(2)>=1.0D0) posit(2) = posit(2)-1.0D0
               if(posit(3)>=1.0D0) posit(3) = posit(3)-1.0D0
               if(posit(1)< 0.0D0) posit(1) = posit(1)+1.0D0
               if(posit(2)< 0.0D0) posit(2) = posit(2)+1.0D0
               if(posit(3)< 0.0D0) posit(3) = posit(3)+1.0D0
               IF(rmcversion==6) THEN
                  WRITE(IWR, 1700) ii, atom_names(cr_iscat(1,iatom))(1:2),   &
                     posit, &
                     is, icell(1)-1, icell(2)-1, icell(3)-1
!                    ((cr_pos(i,iatom)-shift(i))/REAL(scalef(i)),i=1,3),  &
               ELSE
                  WRITE(IWR, 1800) ii, atom_names(cr_iscat(1,iatom))(1:9),   &
                     posit, &
                     is, icell(1)-1, icell(2)-1, icell(3)-1
!                    ((cr_pos(i,iatom)-shift(i))/REAL(scalef(i)),i=1,3),  &
               ENDIF
            ENDIF
!        ENDIF
      ENDIF
   ENDDO
ENDDO
!
CLOSE(IWR)
DEALLOCATE(atom_names)
DEALLOCATE(atom_number)
if(allocated(isites)) deallocate(isites)
!
1100 FORMAT('Metadate date:      ',I2.2,'-',I2.2,'-',I4)
1110 FORMAT('Number of types of atoms: ', I27)
1200 FORMAT('Number of atoms: ', I27)
1300 FORMAT('Number density (Ang^-3): ', F24.15)
1400 FORMAT('Supercell dimensions:             ',3I4)
!1501 FORMAT('Cell (Ang/deg):',6(1x,F20.15))
1500 FORMAT('Cell (Ang/deg):',6(1x,F10.5,'0000000000'))
1600 FORMAT(3(1x,F10.5,'0000000000'))
1700 FORMAT(i6,3x,a2,3(1x,f9.6,'000000000'),i6,3i4)
1800 FORMAT(i6,3x,a9,3(1x,f9.6,'000000000'),i6,3i4)
!
END SUBROUTINE discus2rmc6f
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus2poscar (ianz, cpara, lpara, MAXW)
!
!  Exports the current structure in VASP POSCAR format
!
USE crystal_mod
USE celltoindex_mod
USE modify_func_mod
!
USE charact_mod
USE build_name_mod
USE errlist_mod
USE param_mod
USE precision_mod
USE prompt_mod
USE times_mod
USE trig_degree_mod
USE support_mod
!
IMPLICIT NONE
!
INTEGER            ,                  INTENT(IN)    :: MAXW
INTEGER            ,                  INTENT(INOUT) :: ianz
CHARACTER(LEN=*   ), DIMENSION(MAXW), INTENT(INOUT) :: cpara
INTEGER            , DIMENSION(MAXW), INTENT(INOUT) :: lpara
!
INTEGER, PARAMETER :: IWR = 35
!
CHARACTER(LEN=    PREC_STRING            )      :: ofile  = ' '
CHARACTER(LEN=    PREC_STRING            )      :: string = ' '
CHARACTER(LEN=9), DIMENSION(:), ALLOCATABLE :: atom_names
INTEGER         , DIMENSION(:), ALLOCATABLE :: atom_number
INTEGER                  :: iatom, iscat, i, j, ii, ll
INTEGER                  :: ntypes      ! Actual atom types      to be written to file
INTEGER                  :: natoms      ! Actual number of atoms to be written to file
REAL(KIND=PREC_DP)   , DIMENSION(3)    :: shift
REAL(KIND=PREC_DP)   , DIMENSION(MAXW) :: werte
!
! Build the output file name
!
CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1)
IF (ier_num.ne.0) THEN
   RETURN
ENDIF
ofile = cpara (1)
CALL oeffne (IWR, ofile, 'unknown')
IF (ier_num.ne.0) THEN
   CLOSE(IWR)
   RETURN
ENDIF
!
! Make VASP atom names
!
ALLOCATE(atom_names(0:MAXSCAT))
ALLOCATE(atom_number(0:MAXSCAT))
atom_names(:)  = ' '
atom_number(:) = 0
   DO i=0,MAXSCAT
      atom_names(i)(1:4) = cr_at_lis(i)(1:4)
      IF(IACHAR(atom_names(i)(2:2))>=aa .AND. IACHAR(atom_names(i)(2:2))<=zz) THEN
         atom_names(i)(2:2) = ACHAR(IACHAR(atom_names(i)(2:2)) - aa + a)
      ENDIF
      ii = 1
      DO j=0,i-1
         IF(atom_names(i)(1:4) == atom_names(j)(1:4)) THEN
            ii = ii + 1
         ENDIF
      ENDDO
      IF(ii<10) THEN
         WRITE(atom_names(i)(5:9),'(A2,I1,A2)') ' [',ii,'] '
      ELSEIF(ii<100) THEN
         WRITE(atom_names(i)(5:9),'(A2,I2,A1)') ' [',ii,']'
      ENDIF
   ENDDO
!
! Count atoms that will be written to file
!
natoms = 0
ntypes = 0
count: DO iatom=1,cr_natoms
   IF(cr_at_lis(cr_iscat(1,iatom))/='VOID') THEN
      IF (check_select_status (iatom, .true., cr_prop (iatom),  cr_sel_prop) ) THEN
         natoms = natoms + 1
         find_entry: DO j=1,MAXSCAT
            IF(atom_names(cr_iscat(1,iatom))(1:4)==atom_names(j)(1:4)) THEN
               atom_number(j) = atom_number(j) + 1
               EXIT find_entry
            ENDIF
         ENDDO find_entry
      ENDIF
   ENDIF
ENDDO count
DO i=1,MAXSCAT
   IF(atom_number(i) > 0 ) ntypes = ntypes + 1
ENDDO
!
! Write the actual file
!
CALL datum
WRITE(IWR, '(a,a,a,a)') cr_name(1:LEN_TRIM(cr_name)), ' (Generated by DISCUS Version ',version,')'
WRITE(IWR, '(a)')   '  1.0000'
IF(cr_syst==cr_ortho .OR. cr_syst==cr_tetragonal .OR. cr_syst==cr_cubic) THEN
   WRITE(IWR, 1600 )   cr_a0(1)*REAL(cr_icc(1)), 0.000000,                  0.000000
   WRITE(IWR, 1600 )   0.000000,                  cr_a0(2)*REAL(cr_icc(2)), 0.000000
   WRITE(IWR, 1600 )   0.000000,                  0.000000,                  cr_a0(3)*REAL(cr_icc(3))
ELSE
   WRITE(IWR, 1600 )   cr_a0(1)*REAL(cr_icc(1)), 0.000000, 0.00000
   WRITE(IWR, 1600 )   cr_a0(2)*REAL(cr_icc(2))*cosd(cr_win(3))                     , &
                       cr_a0(2)*REAL(cr_icc(2))*sind(cr_win(3)), 0.000000
   WRITE(IWR, 1600 )   cr_a0(3)*REAL(cr_icc(3))*cosd(cr_win(2))                     , &
                       cr_a0(3)*REAL(cr_icc(3))*sind(cr_win(2))*cosd(cr_win(1))     , &
                       cr_v    *REAL(cr_icc(3))/(cr_a0(1)*cr_a0(2)*sind(cr_win(3)))
ENDIF
!
string = '   '
ll     =  3
DO i=1, MAXSCAT
   IF(atom_number(i) > 0) THEN
      WRITE(string(ll+1:ll+13),'(I12,1x)') atom_number(i)
      ll = ll + 13
   ENDIF
ENDDO
WRITE(IWR, '(a)')   string(1:LEN_TRIM(string))
WRITE(IWR, '(a)')   'Direct'
!
!
! Finished header, loop over all atom types and all atoms
! Property induced selection rules hold, VOIDS are never written
!
shift(:) = REAL(NINT(cr_dim(:,1)))
DO iscat=1,MAXSCAT
   DO iatom=1,cr_natoms
      IF(cr_iscat(1,iatom)==iscat) THEN
         IF(cr_at_lis(cr_iscat(1,iatom))/='VOID') THEN
            IF (check_select_status (iatom, .true., cr_prop (iatom),  cr_sel_prop) ) THEN
                  WRITE(IWR, 1700) ((cr_pos(i,iatom)-shift(i))/REAL(cr_icc(i)),i=1,3)
            ENDIF
         ENDIF
      ENDIF
   ENDDO
ENDDO
!
CLOSE(IWR)
DEALLOCATE(atom_names)
!
1600 FORMAT(3(1x,F13.7))
1700 FORMAT(3(1x,f10.7))
!
END SUBROUTINE discus2poscar
!
SUBROUTINE discus2scatty(ianz, cpara, lpara, MAXW, scatty_site)
!
! Write the current structure in format for SCATTY (Joe Paddison)
! 
USE crystal_mod
USE atom_name
USE celltoindex_mod
USE chem_mod
USE chem_aver_mod
!
USE build_name_mod
USE errlist_mod
USE string_convert_mod
USE precision_mod
USE support_mod
!use take_param_mod
!
IMPLICIT NONE
!
INTEGER            ,                  INTENT(IN)    :: MAXW
INTEGER            ,                  INTENT(INOUT) :: ianz
CHARACTER(LEN=*   ), DIMENSION(MAXW), INTENT(INOUT) :: cpara
INTEGER            , DIMENSION(MAXW), INTENT(INOUT) :: lpara
integer                             , intent(in)    :: scatty_site
!
INTEGER, PARAMETER :: IWR = 35
!LOGICAL, PARAMETER :: lold = .FALSE.
LOGICAL, PARAMETER :: lout = .FALSE.
LOGICAL, PARAMETER :: lsite= .TRUE.
!
CHARACTER(LEN=    PREC_STRING            )      :: ofile = ' '
CHARACTER(LEN=    PREC_STRING            )      :: line  = ' '
CHARACTER(LEN=    PREC_STRING            )      :: string= ' '
CHARACTER(LEN=   4)      :: at_name_i = ' '
!INTEGER                  :: lp
INTEGER                  :: i, k, ia
INTEGER                  :: isite
INTEGER, DIMENSION(3)    :: icell
!INTEGER                  :: nscat
REAL(KIND=PREC_DP)   , DIMENSION(MAXW) :: werte
REAL(KIND=PREC_DP)   , DIMENSION(3   ) :: vec
!
CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1)
IF (ier_num.ne.0) THEN
   RETURN
ENDIF
ofile = cpara (1)
IF(lpara(1)>13) THEN
   IF(ofile(lpara(1)-12:lpara(1)) /= '_atoms_01.txt') ofile = cpara (1) (1:lpara (1) ) //'_atoms_01.txt'
ELSE
   ofile = cpara (1) (1:lpara (1) ) //'_atoms_01.txt'
ENDIF
CALL oeffne (IWR, ofile, 'unknown')
IF (ier_num.ne.0) THEN
   CLOSE(IWR)
   RETURN
ENDIF
!
! Determine average structure
!
CALL chem_aver(lout, lsite)
!
! Ready to write the structure
!
WRITE(IWR, '(a,a)') 'TITLE ',cr_name
!
ia = cr_icc(1)*cr_icc(2)*cr_icc(3)
if(scatty_site==0) then                   !use average sites 
   WRITE(IWR, '(a, 6(2x,f12.5))') 'CELL', cr_a0(:), cr_win(:)
   WRITE(IWR, '(a, 3(3x,i6))')    'BOX ', cr_icc(:)
   DO i=1, cr_ncatoms
      WRITE(IWR, '(a,3(1x,f15.12))') 'SITE ',chem_ave_pos(:,i)
   ENDDO
   DO i=1, cr_ncatoms
      line = 'OCC '
      types: DO k=1, chem_ave_n(i)
         string = ' '
         at_name_i = cr_at_lis(chem_ave_iscat(k,i))
         IF(at_name_i=='VOID') CYCLE types
         CALL do_cap(at_name_i(1:1))
         CALL do_low(at_name_i(2:4))
         WRITE(string,'(a4,1x,f8.4)') at_name_i, chem_ave_bese(i,k)/REAL(ia)
         line = line(1:LEN_TRIM(line))// ' ' // string(1:11)
      ENDDO types
      WRITE(IWR, '(a)') line(1:LEN_TRIM(line))
   ENDDO
!
! Write atom list
!
   atoms: DO i=1, cr_natoms
      at_name_i = cr_at_lis(cr_iscat(1,i))
      IF(at_name_i=='VOID') CYCLE atoms
      CALL do_cap(at_name_i(1:1))
      CALL do_low(at_name_i(2:4))
      icell = 0
      isite = 0
      CALL indextocell(i, icell, isite)
      WRITE(IWR, '(a,4i12,3(1x,g26.12e3),1x,a4)') 'ATOM ', &
           isite, icell(1)-1, icell(2)-1, icell(3)-1, &
           (cr_pos(1,i) - chem_ave_pos(1,isite))  - &
      NINT((cr_pos(1,i) - chem_ave_pos(1,isite)))  ,&
           (cr_pos(2,i) - chem_ave_pos(2,isite))  - &
      NINT((cr_pos(2,i) - chem_ave_pos(2,isite)))  ,&
           (cr_pos(3,i) - chem_ave_pos(3,isite))  - &
      NINT((cr_pos(3,i) - chem_ave_pos(3,isite))),  &
      at_name_i(1:4)
   ENDDO atoms
else                       ! Use all atom positions
   WRITE(IWR, '(a, 6(2x,f12.5))') 'CELL', cr_icc*cr_a0, cr_win
   WRITE(IWR, '(a, 3(3x,i6))')    'BOX ', 1,1,1
   loop_atoms3:DO i=1, cr_natoms
      at_name_i = cr_at_lis(cr_iscat(1,i))
      IF(at_name_i=='VOID') CYCLE loop_atoms3
      vec = cr_pos(:,i) - real(int(cr_pos(:,i)),kind=PREC_DP)
      if(vec(1)<0.0D0) vec(1:) = vec(1) + 1.0D0
      if(vec(2)<0.0D0) vec(2:) = vec(2) + 1.0D0
      if(vec(3)<0.0D0) vec(3:) = vec(3) + 1.0D0
      WRITE(IWR, '(a,3(1x,f15.12))') 'SITE ',vec
   ENDDO loop_atoms3
   loop_atoms4: DO i=1, cr_natoms
      line = 'OCC '
      string = ' '
      at_name_i = cr_at_lis(cr_iscat(1,i))
      IF(at_name_i=='VOID') CYCLE loop_atoms4
      CALL do_cap(at_name_i(1:1))
      CALL do_low(at_name_i(2:4))
      WRITE(string,'(a4,1x,f8.4)') at_name_i, 1.0D0
      line = line(1:LEN_TRIM(line))// ' ' // string(1:11)
      WRITE(IWR, '(a)') line(1:LEN_TRIM(line))
   ENDDO loop_atoms4
!
! Write atom list
!
   k = 0
   atoms2: DO i=1, cr_natoms
      at_name_i = cr_at_lis(cr_iscat(1,i))
      IF(at_name_i=='VOID') CYCLE atoms2
      CALL do_cap(at_name_i(1:1))
      CALL do_low(at_name_i(2:4))
      k = k + 1
      icell = 0
      isite = 0
      CALL indextocell(i, icell, isite)
      WRITE(IWR, '(a,4i12,3(1x,g26.12e3),1x,a4)') 'ATOM ', &
           k, 0, 0, 0, &
           0.0D0, 0.0D0, 0.0D0, at_name_i(1:4)
!          (cr_pos(1,i) - chem_ave_pos(1,isite))  - &
!     NINT((cr_pos(1,i) - chem_ave_pos(1,isite)))  ,&
!          (cr_pos(2,i) - chem_ave_pos(2,isite))  - &
!     NINT((cr_pos(2,i) - chem_ave_pos(2,isite)))  ,&
!          (cr_pos(3,i) - chem_ave_pos(3,isite))  - &
!     NINT((cr_pos(3,i) - chem_ave_pos(3,isite))),  &
!     at_name_i(1:4)
   ENDDO atoms2
endif
!
CLOSE(IWR)
!
END SUBROUTINE discus2scatty
!
END MODULE discus_export
