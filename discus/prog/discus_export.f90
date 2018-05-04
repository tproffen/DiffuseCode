MODULE discus_export
!
CONTAINS
!
SUBROUTINE do_export(line, lp)
!
! Exports a structure in different formats
!
USE errlist_mod
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
!
LOGICAL str_comp
!
CALL get_params (line, ianz, cpara, lpara, MAXW, lp)
IF (ier_num.ne.0) then
   RETURN
ENDIF
!                                                                       
IF (ianz.ge.1) then
   IF (str_comp (cpara (1) , 'shelx', 2, lpara (1) , 5) ) then
      IF (ianz.eq.2) then
         CALL del_params (1, ianz, cpara, lpara, maxw)
         IF (ier_num.ne.0) return
         CALL discus2ins (ianz, cpara, lpara, MAXW)
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
SUBROUTINE discus2ins (ianz, cpara, lpara, MAXW)
!
USE chem_aver_mod
USE crystal_mod
USE generate_mod
USE diffuse_mod
USE discus_save_mod
USE molecule_mod
USE modify_mod
USE prop_para_mod
USE structur
USE save_menu, ONLY: save_internal
USE wyckoff_mod
!
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
         IF(shelx_names(j)(4:4) == 'Z') THEN
            shelx_names(3:4) = 'AA'
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
   IF(shelx_names(i)(4:4) == 'A') THEN
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
END MODULE discus_export
