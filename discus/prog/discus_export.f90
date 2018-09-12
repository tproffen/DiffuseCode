MODULE discus_export
!
CONTAINS
!
SUBROUTINE do_export(line, lp)
!
! Exports a structure in different formats
!
USE errlist_mod
USE get_params_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT)  :: line
INTEGER         , INTENT(INOUT)  :: lp
!
INTEGER, PARAMETER :: MAXW = 5
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
INTEGER  :: ianz
INTEGER, PARAMETER :: NOPTIONAL = 1
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 1 ! Number of values to calculate 
!
INTEGER  :: rmcversion
!
LOGICAL str_comp
!
DATA oname  / 'version'                    /
DATA loname /  7                           /
!
opara  =  (/ '6.0000' /)   ! Always provide fresh default values
lopara =  (/  6       /)
owerte =  (/  6.0     /)
!
CALL get_params (line, ianz, cpara, lpara, MAXW, lp)
IF (ier_num.ne.0) then
   RETURN
ENDIF
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, owerte)
!
!                                                                       
IF (ianz.ge.1) then
   IF (str_comp (cpara (1) , 'cif', 2, lpara (1) , 3) ) then
      IF (ianz.eq.2) then
         CALL del_params (1, ianz, cpara, lpara, maxw)
         IF (ier_num.ne.0) return
         CALL discus2cif (ianz, cpara, lpara, MAXW)
      ELSE
         ier_num = - 6
         ier_typ = ER_COMM
      ENDIF
   ELSEIF (str_comp (cpara (1) , 'shelx', 2, lpara (1) , 5) ) then
      IF (ianz.eq.2) then
         CALL del_params (1, ianz, cpara, lpara, maxw)
         IF (ier_num.ne.0) return
         CALL discus2ins (ianz, cpara, lpara, MAXW)
      ELSE
         ier_num = - 6
         ier_typ = ER_COMM
      ENDIF
   ELSEIF (str_comp (cpara (1) , 'rmcprofile', 2, lpara (1) , 10) ) then
      IF (ianz.eq.2) then
         CALL del_params (1, ianz, cpara, lpara, maxw)
         IF (ier_num.ne.0) RETURN
         rmcversion = NINT(owerte(1))
         CALL discus2rmc6f (ianz, cpara, lpara, MAXW, rmcversion)
      ELSE
         ier_num = - 6
         ier_typ = ER_COMM
      ENDIF
   ELSEIF (str_comp (cpara (1) , 'vasp', 2, lpara (1) , 4)  .OR.  &
           str_comp (cpara (1) , 'poscar', 2, lpara (1) , 6) ) then
      IF (ianz.eq.2) then
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
!
IMPLICIT NONE
!
INTEGER            ,                  INTENT(IN)    :: MAXW
INTEGER            ,                  INTENT(INOUT) :: ianz
CHARACTER(LEN=1024), DIMENSION(MAXW), INTENT(INOUT) :: cpara
INTEGER            , DIMENSION(MAXW), INTENT(INOUT) :: lpara
!
INTEGER, PARAMETER :: IWR = 35
LOGICAL, PARAMETER :: lold = .FALSE.
!
CHARACTER(LEN=1024)      :: ofile = ' '
CHARACTER(LEN=1024)      :: zeile = ' '
INTEGER                  :: lp
INTEGER                  :: nscat
REAL   , DIMENSION(MAXW) :: werte
!
CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1)
IF (ier_num.ne.0) then
   RETURN
ENDIF
ofile = cpara (1)
IF(ofile(lpara(1)-3:lpara(1)) /= '.cif') ofile = cpara (1) (1:lpara (1) ) //'.cif'
CALL oeffne (IWR, ofile, 'unknown')
IF (ier_num.ne.0) then
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
   CALL alloc_plot ( nscat )
   IF(ier_num < 0) THEN
     RETURN
   ENDIF
ENDIF
zeile = 'all'
lp    = 3
CALL atom_select (zeile, lp, 0, PL_MAXSCAT, pl_latom, &
                  pl_sel_atom, lold, .TRUE.)
!
CALL plot_cif(IWR, .FALSE.)
!
CLOSE(IWR)
!
END SUBROUTINE discus2cif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus2ins (ianz, cpara, lpara, MAXW)
!
USE chem_aver_mod
USE class_internal
USE crystal_mod
USE generate_mod
USE diffuse_mod
USE discus_save_mod
USE molecule_mod
!USE modify_mod
USE prop_para_func
USE prop_para_mod
USE structur
USE save_menu, ONLY: save_internal
USE wyckoff_mod
!
USE build_name_mod
USE errlist_mod
USE param_mod
!
IMPLICIT NONE
!
INTEGER            ,                  INTENT(IN)    :: MAXW
INTEGER            ,                  INTENT(INOUT) :: ianz
CHARACTER(LEN=1024), DIMENSION(MAXW), INTENT(INOUT) :: cpara
INTEGER            , DIMENSION(MAXW), INTENT(INOUT) :: lpara
!
INTEGER, PARAMETER :: IWR = 35
!
CHARACTER(LEN=1024)      :: ofile = ' '
CHARACTER(LEN=1024)      :: line = ' '
CHARACTER(LEN=1024)      :: names = ' '
CHARACTER(LEN=1024)      :: units = ' '
CHARACTER(LEN=19)        :: origfile = ' '
CHARACTER(LEN= 5)        :: befehl   = 'lcell'
INTEGER                  :: i, j, k, l,i1, i2
INTEGER                  :: length
INTEGER                  :: lbef
INTEGER                  :: lattice = 1
INTEGER                  :: unique_n
INTEGER                  :: shelx_n
INTEGER, DIMENSION(:,:), ALLOCATABLE :: n_atoms 
CHARACTER (LEN=2), DIMENSION(:), ALLOCATABLE :: unique_names
CHARACTER (LEN=4), DIMENSION(:), ALLOCATABLE :: shelx_names
INTEGER          , DIMENSION(:), ALLOCATABLE :: unique_n_atoms 
LOGICAL                  :: orig_OK =.FALSE.
REAL                     :: z_unit
REAL   , DIMENSION(MAXW) :: werte
REAL   , DIMENSION(3), PARAMETER :: NULL = (/0.00, 0.00, 0.00/)
!
CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1)
IF (ier_num.ne.0) then
   RETURN
ENDIF
ofile = cpara (1)
IF(ofile(lpara(1)-3:lpara(1)) /= '.ins') ofile = cpara (1) (1:lpara (1) ) //'.ins'
CALL oeffne (IWR, ofile, 'unknown')
IF (ier_num.ne.0) then
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
IF(cr_spcgr(1:1)=='P') lattice = 1
IF(cr_spcgr(1:1)=='I') lattice = 2
IF(cr_spcgr(1:1)=='R') lattice = 3
IF(cr_spcgr(1:1)=='F') lattice = 4
IF(cr_spcgr(1:1)=='A') lattice = 5
IF(cr_spcgr(1:1)=='B') lattice = 6
IF(cr_spcgr(1:1)=='C') lattice = 7
!
! Temporarily save the original structure
!
origfile  = 'internal.shelx_orig'
line       = 'ignore, all'          ! Ignore all properties
length     = 11
CALL property_select(line, length, sav_sel_prop)
CALL save_internal(origfile)
ALLOCATE(n_atoms(2,1:cr_nscat))
CALL chem_elem(.FALSE.)
DO i=1, cr_nscat
   n_atoms(1,i) = NINT(res_para(i+1)*cr_natoms)
ENDDO
!
! Read structure as unit cell file and expand to full crystal
!
befehl = 'lcell'
lbef   = 5
cpara(1) = origfile
lpara(1) = 19
CALL do_readcell(befehl, lbef, ianz, MAXW, cpara, lpara, .TRUE., 1.0E-5, 0)
CALL chem_elem(.FALSE.)
z_unit = 192.0
DO i=1, cr_nscat
   n_atoms(2,i) = NINT(res_para(i+1)*cr_natoms)
   IF(n_atoms(1,i)>0 .AND. n_atoms(2,i)>0) THEN
      z_unit = MIN(z_unit, REAL(n_atoms(2,i)/n_atoms(1,i)))
   ENDIF
ENDDO
!
! Restore original structure
!
CALL do_readstru(origfile)
CALL store_remove_single(origfile, ier_num)
IF(ier_num/=0) THEN
   ier_typ = ER_APPL
   ier_msg(1) = 'Could not remove temporary internal storage'
   ier_msg(2) = 'in shelx export'
   ier_msg(3) = 'Please document and report'
   RETURN
ENDIF
!
! Reduce names to unique atom names
!
ALLOCATE(unique_n_atoms(1:cr_nscat))
ALLOCATE(unique_names(1:cr_nscat))
unique_n = 1
unique_names(1)   = cr_at_lis(1)(1:2)
unique_n_atoms(1) = n_atoms(2,1)
loop_types: DO i=2, cr_nscat
   DO j=1,unique_n
      IF(unique_names(j) == cr_at_lis(i)(1:2)) THEN  ! Previous name
         unique_n_atoms(j) = unique_n_atoms(j) + n_atoms(2,i)
         CYCLE loop_types
      ENDIF
   ENDDO
   unique_n = unique_n + 1
   unique_names(unique_n) = cr_at_lis(i)(1:2)
   unique_n_atoms(unique_n) = n_atoms(2,i)
ENDDO loop_types
!
names = ' '
units = ' '
DO i=1,unique_n                                ! Prepare SFAC / UNIT statement
   i1 = (i-1) * 3 + 1
   i2 = (i  ) * 3
   names(i1:i2) = unique_names(i)(1:2)
   i1 = (i-1) * 6 + 1
   i2 = (i  ) * 6
   WRITE(units(i1:i2),'(1x,i5)') unique_n_atoms(i)
ENDDO
!
! Create unique Shelx names
!
ALLOCATE(shelx_names(1:cr_natoms))
shelx_names(:) = ' '
shelx_n = 0
loop_shelx: DO i=1,cr_natoms
   l = shelx_n
   DO j=l,1,-1
      IF(cr_at_lis(cr_iscat(i))(1:2) == shelx_names(j)(1:2)) THEN
         shelx_n = shelx_n + 1
         shelx_names(i)(1:2) = cr_at_lis(cr_iscat(i))(1:2)
         shelx_names(i)(3:3) = shelx_names(j)(3:3)
         IF(shelx_names(j)(4:4) == 'Z') THEN
            IF(shelx_names(j)(3:3) == 'x') THEN
               shelx_names(i)(3:4) = 'AA'
            ELSE
               shelx_names(i)(3:3) = ACHAR(IACHAR(shelx_names(j)(3:3))+1)
               shelx_names(i)(4:4) = 'A'
            ENDIF
         ELSE
            shelx_names(i)(4:4) = ACHAR(IACHAR(shelx_names(j)(4:4))+1)
         ENDIF
         CYCLE loop_shelx
      ENDIF
   ENDDO
   shelx_n = shelx_n + 1
   shelx_names(shelx_n)(1:2) = cr_at_lis(cr_iscat(i))(1:2)
   shelx_names(shelx_n)(3:4) = 'xA'
ENDDO loop_shelx
!
DO i=1,shelx_n
   IF(shelx_names(i)(3:4) == 'xA') THEN
      shelx_names(i)(3:4) = '  '
      IF(shelx_names(i)(2:2) == 'x') shelx_names(i)(2:2) = ' '
   ELSE
      DO j=1,4
         IF(shelx_names(i)(j:j) == ' ') shelx_names(i)(j:j) = 'x'
      ENDDO
   ENDIF
ENDDO
!
!  Start writing the actual file 
!
WRITE(IWR, 1000) cr_name(1:LEN_TRIM(cr_name)), cr_spcgr(1:LEN_TRIM(cr_spcgr))
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
   j =spc_n / 3                           ! Only need the first third symmetry operations
ENDIF
IF(.NOT.cr_acentric) j = j / 2            ! Centrosymmetric need the first half only
DO i=2, j                                 ! Omit identity x,y,z
   WRITE(IWR, 1400) spc_xyz(i)(1:LEN_TRIM(spc_xyz(i)))
ENDDO
!
WRITE(IWR, 1500) names(1:LEN_TRIM(names)) ! SFAC
WRITE(IWR, 1600) units(1:LEN_TRIM(units)) ! UNIT
WRITE(IWR, 1700)             ! L.S.
WRITE(IWR, 1800)             ! BOND
WRITE(IWR, 1900)             ! BOND $H
WRITE(IWR, 2000)             ! FMAP 2
WRITE(IWR, 2100)             ! PLAN 20
WRITE(IWR, 2200)             ! ACTA
WRITE(IWR, 2300)             ! WGHT
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
CLOSE(IWR)
!
1000 FORMAT('TITL ',a,' in ',a)
1100 FORMAT('CELL ',f7.5,3f11.6,3f9.4)
1200 FORMAT('ZERR ',f7.2,3f11.6,3f9.4)
1300 FORMAT('LATT ', i1)
1310 FORMAT('LATT ', i2)
1400 FORMAT('SYMM ', a )
1500 FORMAT('SFAC ', a )
1600 FORMAT('UNIT ', a )
1700 FORMAT('L.S. 100' )
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
END SUBROUTINE discus2ins
!
SUBROUTINE shelx_write_atom(IWR, i, unique_n, unique_names, natoms, shelx_names)
!
USE crystal_mod
USE param_mod
USE spcgr_apply, ONLY: get_wyckoff
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
INTEGER :: stype
INTEGER :: j
REAL               :: occup,biso
REAL, DIMENSION(3) :: vec
!
   stype = 0
   loop_stype: DO j=1,unique_n
      IF(cr_at_lis(cr_iscat(i))(1:2) == unique_names(j)(1:2))  THEN
         stype = j
         EXIT loop_stype
      ENDIF
   ENDDO loop_stype
   vec(:) = cr_pos(:,i)
   CALL get_wyckoff(vec,.FALSE.,1)
   occup = 10.000 + REAL(res_para(1)/res_para(3))*cr_occ(cr_iscat(i))
   biso = cr_dw(cr_iscat(i))/8./REAL(pi**2)
   WRITE(IWR,2500) shelx_names(i), stype, cr_pos(:,i), occup,biso
!
2500 FORMAT(a4,1x,i2,3(f12.6),f12.5,f11.5)
!
END SUBROUTINE shelx_write_atom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus2rmc6f (ianz, cpara, lpara, MAXW, rmcversion)
!
!  Exports the current structure in RMC6F format
!
USE crystal_mod
USE celltoindex_mod
USE modify_func_mod
!
USE charact_mod
USE build_name_mod
USE errlist_mod
USE param_mod
USE prompt_mod
USE times_mod
USE trig_degree_mod
!
IMPLICIT NONE
!
INTEGER            ,                  INTENT(IN)    :: MAXW
INTEGER            ,                  INTENT(INOUT) :: ianz
CHARACTER(LEN=1024), DIMENSION(MAXW), INTENT(INOUT) :: cpara
INTEGER            , DIMENSION(MAXW), INTENT(INOUT) :: lpara
INTEGER                             , INTENT(IN)    :: rmcversion     ! Version is 6 or 7
!
INTEGER, PARAMETER :: IWR = 35
!
CHARACTER(LEN=1024)      :: ofile  = ' '
CHARACTER(LEN=1024)      :: line   = ' '
CHARACTER(LEN=1024)      :: string = ' '
CHARACTER(LEN=9), DIMENSION(:), ALLOCATABLE :: atom_names
INTEGER         , DIMENSION(:), ALLOCATABLE :: atom_number
INTEGER                  :: iatom, iscat, is, i, j, ii, ll
INTEGER                  :: ntypes      ! Actual atom types      to be written to file
INTEGER                  :: natoms      ! Actual number of atoms to be written to file
INTEGER, DIMENSION(3)    :: icell 
REAL   , DIMENSION(3)    :: shift
REAL   , DIMENSION(MAXW) :: werte
!
! Build the output file name
!
CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1)
IF (ier_num.ne.0) then
   RETURN
ENDIF
ofile = cpara (1)
IF(.NOT.(ofile(lpara(1)-5:lpara(1)) == '.rmc6f')) THEN
  ofile = cpara (1) (1:lpara (1) ) //'.rmc6f'
ENDIF
CALL oeffne (IWR, ofile, 'unknown')
IF (ier_num.ne.0) then
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
ntypes = 0
count: DO iatom=1,cr_natoms
   IF(cr_at_lis(cr_iscat(iatom))/='VOID') THEN
      IF (check_select_status (iatom, .true., cr_prop (iatom),  cr_sel_prop) ) THEN
         natoms = natoms + 1
         find_entry: DO j=1,MAXSCAT
            IF(atom_names(cr_iscat(iatom))(1:4)==atom_names(j)(1:4)) THEN
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
WRITE(IWR, '(a)')   '(Version 6f format configuration file)'
WRITE(IWR, '(a,a,a)') '(Generated by DISCUS Version ',version,')'
WRITE(IWR, '(a)')   'Metadata owner:     DISCUS'
WRITE(IWR, 1100 )    int_date(3), int_date(2), int_date(1)
WRITE(IWR, '(a,a)') 'Metadata comment:   ',cr_name(1:LEN_TRIM(cr_name))
WRITE(IWR, 1110 )   ntypes
line   = 'Atom types present:         '
string = 'Number of each atom type:   '
ii     = 28
ll     = 28
DO i=1, MAXSCAT
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
WRITE(IWR, 1300 )   FLOAT(natoms)/(FLOAT(cr_icc(1)*cr_icc(2)*cr_icc(3))*cr_v)
WRITE(IWR, 1400 )   cr_icc(:)
WRITE(IWR, 1500 )   (cr_a0 (i)*FLOAT(cr_icc(i)),i=1,3), cr_win(:)
WRITE(IWR, '(a)')   'Lattice vectors (Ang):'
IF(cr_syst==cr_ortho .OR. cr_syst==cr_tetragonal .OR. cr_syst==cr_cubic) THEN
   WRITE(IWR, 1600 )   cr_a0(1)*FLOAT(cr_icc(1)), 0.000000,                  0.000000
   WRITE(IWR, 1600 )   0.000000,                  cr_a0(2)*FLOAT(cr_icc(2)), 0.000000
   WRITE(IWR, 1600 )   0.000000,                  0.000000,                  cr_a0(3)*FLOAT(cr_icc(3))
ELSE
   WRITE(IWR, 1600 )   cr_a0(1)*FLOAT(cr_icc(1)), 0.000000, 0.00000
   WRITE(IWR, 1600 )   cr_a0(2)*FLOAT(cr_icc(2))*cosd(cr_win(3))                     , &
                       cr_a0(2)*FLOAT(cr_icc(2))*sind(cr_win(3)), 0.000000
   WRITE(IWR, 1600 )   cr_a0(3)*FLOAT(cr_icc(3))*cosd(cr_win(2))                     , &
                       cr_a0(3)*FLOAT(cr_icc(3))*sind(cr_win(2))*cosd(cr_win(1))     , &
                       cr_v    *FLOAT(cr_icc(3))/(cr_a0(1)*cr_a0(2)*sind(cr_win(3)))
ENDIF
WRITE(IWR, '(a)')   'Atoms:'
!
! Finished header, loop over all atom types and all atoms
! Property induced selection rules hold, VOIDS are never written
!
shift(:) = FLOAT(NINT(cr_dim(:,1)))
ii       = 0
DO iscat=1,MAXSCAT
   DO iatom=1,cr_natoms
      IF(cr_iscat(iatom)==iscat) THEN
         IF(cr_at_lis(cr_iscat(iatom))/='VOID') THEN
            IF (check_select_status (iatom, .true., cr_prop (iatom),  cr_sel_prop) ) THEN
               ii = ii + 1                               ! increment atom number for RMCprofile sequence
               CALL indextocell (iatom, icell, is)
               IF(rmcversion==6) THEN
                  WRITE(IWR, 1700) ii, atom_names(cr_iscat(iatom))(1:2),   &
                     ((cr_pos(i,iatom)-shift(i))/FLOAT(cr_icc(i)),i=1,3),  &
                     is, icell(1)-1, icell(2)-1, icell(3)-1
               ELSE
                  WRITE(IWR, 1800) ii, atom_names(cr_iscat(iatom))(1:9),   &
                     ((cr_pos(i,iatom)-shift(i))/FLOAT(cr_icc(i)),i=1,3),  &
                     is, icell(1)-1, icell(2)-1, icell(3)-1
               ENDIF
            ENDIF
         ENDIF
      ENDIF
   ENDDO
ENDDO
!
CLOSE(IWR)
DEALLOCATE(atom_names)
!
1100 FORMAT('Metadate date:      ',I2.2,'-',I2.2,'-',I4)
1110 FORMAT('Number of types of atoms: ', I27)
1200 FORMAT('Number of atoms: ', I27)
1300 FORMAT('Number density (Ang^-3): ', F24.15)
1400 FORMAT('Supercell dimensions:             ',3I4)
1500 FORMAT('Cell (Ang/deg):',6(1x,F20.15))
1600 FORMAT(3(1x,F20.15))
1700 FORMAT(i6,3x,a2,3(1x,f18.15),i6,3i4)
1800 FORMAT(i6,3x,a9,3(1x,f18.15),i6,3i4)
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
USE prompt_mod
USE times_mod
USE trig_degree_mod
!
IMPLICIT NONE
!
INTEGER            ,                  INTENT(IN)    :: MAXW
INTEGER            ,                  INTENT(INOUT) :: ianz
CHARACTER(LEN=1024), DIMENSION(MAXW), INTENT(INOUT) :: cpara
INTEGER            , DIMENSION(MAXW), INTENT(INOUT) :: lpara
!
INTEGER, PARAMETER :: IWR = 35
!
CHARACTER(LEN=1024)      :: ofile  = ' '
CHARACTER(LEN=1024)      :: line   = ' '
CHARACTER(LEN=1024)      :: string = ' '
CHARACTER(LEN=9), DIMENSION(:), ALLOCATABLE :: atom_names
INTEGER         , DIMENSION(:), ALLOCATABLE :: atom_number
INTEGER                  :: iatom, iscat, is, i, j, ii, ll
INTEGER                  :: ntypes      ! Actual atom types      to be written to file
INTEGER                  :: natoms      ! Actual number of atoms to be written to file
INTEGER, DIMENSION(3)    :: icell 
REAL   , DIMENSION(3)    :: shift
REAL   , DIMENSION(MAXW) :: werte
!
! Build the output file name
!
CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1)
IF (ier_num.ne.0) then
   RETURN
ENDIF
ofile = cpara (1)
CALL oeffne (IWR, ofile, 'unknown')
IF (ier_num.ne.0) then
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
   IF(cr_at_lis(cr_iscat(iatom))/='VOID') THEN
      IF (check_select_status (iatom, .true., cr_prop (iatom),  cr_sel_prop) ) THEN
         natoms = natoms + 1
         find_entry: DO j=1,MAXSCAT
            IF(atom_names(cr_iscat(iatom))(1:4)==atom_names(j)(1:4)) THEN
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
   WRITE(IWR, 1600 )   cr_a0(1)*FLOAT(cr_icc(1)), 0.000000,                  0.000000
   WRITE(IWR, 1600 )   0.000000,                  cr_a0(2)*FLOAT(cr_icc(2)), 0.000000
   WRITE(IWR, 1600 )   0.000000,                  0.000000,                  cr_a0(3)*FLOAT(cr_icc(3))
ELSE
   WRITE(IWR, 1600 )   cr_a0(1)*FLOAT(cr_icc(1)), 0.000000, 0.00000
   WRITE(IWR, 1600 )   cr_a0(2)*FLOAT(cr_icc(2))*cosd(cr_win(3))                     , &
                       cr_a0(2)*FLOAT(cr_icc(2))*sind(cr_win(3)), 0.000000
   WRITE(IWR, 1600 )   cr_a0(3)*FLOAT(cr_icc(3))*cosd(cr_win(2))                     , &
                       cr_a0(3)*FLOAT(cr_icc(3))*sind(cr_win(2))*cosd(cr_win(1))     , &
                       cr_v    *FLOAT(cr_icc(3))/(cr_a0(1)*cr_a0(2)*sind(cr_win(3)))
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
shift(:) = FLOAT(NINT(cr_dim(:,1)))
DO iscat=1,MAXSCAT
   DO iatom=1,cr_natoms
      IF(cr_iscat(iatom)==iscat) THEN
         IF(cr_at_lis(cr_iscat(iatom))/='VOID') THEN
            IF (check_select_status (iatom, .true., cr_prop (iatom),  cr_sel_prop) ) THEN
                  WRITE(IWR, 1700) ((cr_pos(i,iatom)-shift(i))/FLOAT(cr_icc(i)),i=1,3)
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
!
END MODULE discus_export
