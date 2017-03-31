MODULE molecule_func_mod
!
USE conn_mod
USE molecule_mod
!
IMPLICIT NONE
PRIVATE
PUBLIC  :: do_molecularize
PUBLIC  :: molecularize_sub
PUBLIC  :: molecularize_numbers
!
!
LOGICAL, DIMENSION(:  ), ALLOCATABLE :: t_list    ! Temporary list of all atoms in the molecule
!
CONTAINS
!
SUBROUTINE do_molecularize (line, length)
!
USE crystal_mod
!
USE errlist_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: line
INTEGER         , INTENT(IN) :: length
!
INTEGER            , PARAMETER       :: MAXW= 20
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara ! 
INTEGER            , DIMENSION(MAXW) :: lpara ! 
REAL               , DIMENSION(MAXW) :: werte ! 
INTEGER                              :: ianz
INTEGER                              :: istart, ifinish
INTEGER                              :: natom
INTEGER, DIMENSION(:), ALLOCATABLE   :: iatoms
INTEGER                              :: itype
INTEGER                              :: i
REAL                                 :: biso
!
LOGICAL :: str_comp
!
CALL get_params (line, ianz, cpara, lpara, MAXW, length)
!      CALL ber_params (ianz, cpara, lpara, werte, maxw)
!call molecularize_numbers(nint(werte(1)),nint(werte(2)), 1, 0.0)
IF(ier_num == 0) THEN
   IF(IANZ > 0 ) THEN
      IF(str_comp(cpara (1) , 'conn', 1, lpara (1) , 4)) THEN
         CALL del_params (1, ianz, cpara, lpara, maxw)
         CALL ber_params (ianz, cpara, lpara, werte, maxw)
         IF(ier_num == 0) THEN
            IF(NINT(werte(1)) == -1 .AND. ianz==1) THEN
               istart  = 1
               ifinish = cr_natoms
            ELSEIF(NINT(werte(1)) > 0) THEN
               istart  = NINT(werte(1))
               ifinish = istart
            ELSE
               ier_num = -6
               ier_typ = ER_FORT
            ENDIF
            IF(ier_num == 0) THEN
               natom = ianz
               ALLOCATE(iatoms(1:natom))
               iatoms(1:natom) = NINT(werte(1:ianz))
               DO i=istart,ifinish
                  iatoms(1) = i
                  CALL molecularize(natom,iatoms)
               ENDDO
               DEALLOCATE(iatoms)
            ENDIF
         ENDIF
      ELSEIF(str_comp(cpara (1) , 'range', 1, lpara (1) , 5)) THEN
         CALL del_params (1, ianz, cpara, lpara, maxw)
         CALL ber_params (ianz, cpara, lpara, werte, maxw)
         IF(ier_num == 0) THEN
            IF(ianz == 4) THEN
            istart  = NINT(werte(1))
            ifinish = NINT(werte(2))
            itype   = NINT(werte(3))
            biso    =      werte(4)
            CALL molecularize_numbers(istart, ifinish, itype, biso)
            ELSE
               ier_num = -6
               ier_typ = ER_FORT
               ier_msg(1) = 'Molecularize command needs parameters:'
               ier_msg(2) = ' range, <from> , <to>, <type>, <biso>'
            ENDIF
         ENDIF
      ENDIF
   ELSE
      ier_num = -6
      ier_typ = ER_FORT
      ier_msg(1) = 'Molecularize command needs parameters:'
      ier_msg(2) = 'conn, <central> [, <exclude>...]'
      ier_msg(3) = 'range, <from> , <to>, <type>, <biso>'
   ENDIF
ENDIF
!
END SUBROUTINE do_molecularize
!
SUBROUTINE molecularize(natom,iatoms)
!
! Groups atoms starting with the first into a new molecule with new molecule type
! Atoms no. 2 etc and their connectivity are excluded.
!
USE discus_allocate_appl_mod
USE crystal_mod
USE prop_para_mod
USE errlist_mod
IMPLICIT NONE
!
INTEGER,                                  INTENT(IN)  :: natom
INTEGER, DIMENSION(1:natom),              INTENT(IN)  :: iatoms
!
INTEGER   :: jatom
!
INTEGER   :: katom
INTEGER   :: n_new,n_atom, n_type, n_mole
INTEGER   :: i, j
INTEGER   :: imole
!
jatom   = iatoms(1)                        ! Set central atom number
IF(.NOT. btest(cr_prop(jatom),1)) THEN      ! Atom is not yet inside a molecule
!
   ALLOCATE(t_list(1:cr_natoms))              ! Make temporary list of atoms already in the molecule
   t_list(:)     = .FALSE.
   t_list(jatom) = .TRUE.
!
!  Accumulate all neighbors of the current atom into temp_list
!
   IF(natom==1) THEN
      CALL mol_add_conn(jatom)
   ELSE
      CALL mol_add_excl(jatom, natom, iatoms)
   ENDIF
   n_new = 0
   DO i=1, UBOUND(t_list,1)
      IF(t_list(i)) THEN
         n_new = n_new + 1
      ENDIF
   ENDDO
!
   n_mole = mole_num_mole + 1   ! Reserve new molecule
   n_type = mole_num_type + 1   ! Reserve new molecule type
   n_atom = mole_off(mole_num_mole) + mole_len(mole_num_mole) + n_new

   IF ( n_mole > MOLE_MAX_MOLE  .or.  &
        n_type > MOLE_MAX_TYPE  .or.  &
        n_atom > MOLE_MAX_ATOM      ) THEN
      n_mole = MAX(n_mole,MOLE_MAX_MOLE)
      n_type = MAX(n_type,MOLE_MAX_TYPE)
      n_atom = MAX(n_atom,MOLE_MAX_ATOM)
      CALL alloc_molecule(1, 1, n_mole, n_type, n_atom)
      IF ( ier_num /= 0) THEN
         RETURN
      ENDIF
   ENDIF
   mole_num_mole = mole_num_mole + 1
   mole_num_type = mole_num_type + 1
   mole_len(mole_num_mole)  = n_new
   mole_off(mole_num_mole)  = mole_off(mole_num_mole-1)+mole_len(mole_num_mole-1)
   mole_type(mole_num_mole) = n_type
   mole_char(mole_num_mole) = MOLE_ATOM
   mole_file(mole_num_mole) = ' '
   mole_biso(mole_num_mole) = 0.0
   mole_num_atom = mole_off(mole_num_mole) + mole_len(mole_num_mole)
   j = 0
   DO i=1,UBOUND(t_list,1)
      IF(t_list(i)) THEN
         j = j + 1
         mole_cont(mole_off(mole_num_mole)+j) = i
         cr_prop(i) = ibset(cr_prop(i),PROP_MOLECULE)
      ENDIF
   ENDDO
!
   DEALLOCATE(t_list)
!
ELSE                                          ! Atom is already inside a molecule
!
   ALLOCATE(t_list(1:cr_natoms))              ! Make temporary list of atoms already in the molecule
   t_list(:)     = .FALSE.
   t_list(jatom) = .TRUE.
!
!  Accumulate all neighbors of the current atom into temp_list
!
   IF(natom==1) THEN
      CALL mol_add_conn(jatom)
   ELSE
      CALL mol_add_excl(jatom, natom, iatoms)
   ENDIF
   imole = cr_mole(jatom)                    ! Atom is in this molecule
   n_new = 0
   DO i=1, UBOUND(t_list,1)
      IF(t_list(i) .AND. cr_mole(i)/=imole) THEN
         n_new = n_new + 1
      ENDIF
   ENDDO
!
   CALL molecule_shift(n_new, imole+1)
   j = mole_len(imole)
   DO i=1,UBOUND(t_list,1)
      IF(t_list(i).AND.cr_mole(i)/=imole) THEN
         j = j + 1
         mole_cont(mole_off(imole)+j) = i
         cr_prop(i) = ibset(cr_prop(i),PROP_MOLECULE)
         cr_mole(i) = imole
      ENDIF
   ENDDO
   mole_len(imole) = mole_len(imole) + n_new
!
   DEALLOCATE(t_list)
!
ENDIF
!
END SUBROUTINE molecularize
!
!
SUBROUTINE molecularize_sub(jatom,natom,iatoms,max_sub,n_sub, sub_list)
!
! Groups atoms starting with the first into a new molecule with new molecule type
! Atoms no. 2 etc and their connectivity are excluded.
!
USE discus_allocate_appl_mod
USE crystal_mod
USE errlist_mod
IMPLICIT NONE
!
INTEGER,                                  INTENT(IN)  :: jatom
INTEGER,                                  INTENT(IN)  :: natom
INTEGER, DIMENSION(1:natom),              INTENT(IN)  :: iatoms
INTEGER,                                  INTENT(IN)  :: max_sub
INTEGER,                                  INTENT(OUT) :: n_sub
INTEGER, DIMENSION(1:max_sub),            INTENT(OUT) :: sub_list
!
INTEGER   :: i
!
!ALLOCATE(t_list(1:max_sub))              ! Make temporary list of atoms already in the molecule
ALLOCATE(t_list(1:cr_natoms))              ! Make temporary list of atoms already in the molecule
t_list(:)     = .FALSE.
t_list(jatom) = .TRUE.
!
!  Accumulate all neighbors of the current atom into temp_list
!
IF(natom==0) THEN
   CALL mol_add_conn(jatom)
ELSE
   CALL mol_add_excl(jatom, natom, iatoms)
ENDIF
n_sub = 0
DO i=1, UBOUND(t_list,1)
   IF(t_list(i)) THEN
      n_sub = n_sub + 1
      sub_list(n_sub) = i
   ENDIF
ENDDO
!write(*,*) ' Centr ', jatom
!write(*,*) ' Excl  ', natom,' : ', iatoms(:)
!write(*,*) ' GROUP ', n_sub
!write(*,*) ' GROUP ', sub_list(:)
!
!write(*,*) ' STATUS ', t_list
DEALLOCATE(t_list)
!
END SUBROUTINE molecularize_sub
!
!
RECURSIVE SUBROUTINE mol_add_conn(katom)
!
USE crystal_mod
IMPLICIT NONE
!
INTEGER , INTENT(IN) :: katom
!
INTEGER   :: ktype
INTEGER   :: ino_max  ! Number of connectivities around central atom
INTEGER   :: ino
INTEGER                              :: c_natoms  ! Number of atoms connected
CHARACTER(LEN=256)                   :: c_name    ! Connectivity name
INTEGER, DIMENSION(:  ), ALLOCATABLE :: c_list    ! List of atoms connected to current
INTEGER, DIMENSION(:,:), ALLOCATABLE :: c_offs    ! Offsets for atoms connected to current
!
INTEGER :: latom
INTEGER :: i
!
!  Accumulate all neighbors of the current atom into temp_neig
!
ktype = cr_iscat(katom)
ino_max = get_connectivity_numbers(ktype)  ! Get number of connectivity definitions
DO ino = 1, ino_max                        ! Loop over all connectivity types
   CALL get_connectivity_list(katom, ktype, ino, c_list, c_offs, c_natoms)
   IF(c_natoms>0) THEN
   DO i=1,c_natoms                         ! Loop over all neighbors
      IF(.NOT.t_list(c_list(i))) THEN      ! Not yet in the molecule
         t_list(c_list(i)) = .TRUE.        ! Add to temporary list
         latom = c_list(i)                 ! Now search current neighbor for further ones
!        Apply periodic boundary shifts to keep molecule contiguous
         IF(c_offs(1,i)/=0 .OR. c_offs(2,i)/=0 .OR. c_offs(3,i)/=0 ) THEN
            cr_pos(:,latom) = cr_pos(:,latom) + c_offs(:,i)
            CALL conn_update(latom, FLOAT(c_offs(:,i)))
         ENDIF
         CALL mol_add_conn(latom)
      ENDIF
   ENDDO
   DEALLOCATE(c_list)                      ! Remove temporary list of atoms
   DEALLOCATE(c_offs)                      ! Remove temporary offsets
   ENDIF
ENDDO
END SUBROUTINE mol_add_conn
!
RECURSIVE SUBROUTINE mol_add_excl(katom, n_excl, excl)
!
USE crystal_mod
IMPLICIT NONE
!
INTEGER , INTENT(IN) :: katom
INTEGER , INTENT(IN) :: n_excl
INTEGER , DIMENSION(1:n_excl), INTENT(IN) :: excl
!
INTEGER   :: ktype
INTEGER   :: ino_max  ! Number of connectivities around central atom
INTEGER   :: ino
INTEGER                              :: c_natoms  ! Number of atoms connected
CHARACTER(LEN=256)                   :: c_name    ! Connectivity name
INTEGER, DIMENSION(:  ), ALLOCATABLE :: c_list    ! List of atoms connected to current
INTEGER, DIMENSION(:,:), ALLOCATABLE :: c_offs    ! Offsets for atoms connected to current
!
INTEGER :: latom
INTEGER :: i, j
!
!  Accumulate all neighbors of the current atom into temp_neig
!
ktype = cr_iscat(katom)
ino_max = get_connectivity_numbers(ktype)  ! Get number of connectivity definitions
!write(*,*) ' EXCL START  ', katom, ktype, ino_max
DO ino = 1, ino_max                        ! Loop over all connectivity types
   CALL get_connectivity_list(katom, ktype, ino, c_list, c_offs, c_natoms)
!write(*,*) ' EXCL CENTRAL', katom, c_list(1:c_natoms)
!write(*,*) ' EXCL T_LIST ', t_list
!write(*,*) ' EXCL excl   ', n_excl,' : ', excl(:)
   search: DO i=1,c_natoms                         ! Loop over all neighbors
      DO j=1,n_excl
         IF(c_list(i)==excl(j)) CYCLE search
      ENDDO
!write(*,*) ' EXCL NOT EX ', i, c_list(i), t_list(c_list(i))
      IF(.NOT.t_list(c_list(i))) THEN      ! Not yet in the molecule
         t_list(c_list(i)) = .TRUE.        ! Add to temporary list
         latom = c_list(i)                 ! Now search current neighbor for further ones
!        Apply periodic boundary shifts to keep molecule contiguous
         IF(c_offs(1,i)/=0 .OR. c_offs(2,i)/=0 .OR. c_offs(3,i)/=0 ) THEN
            cr_pos(:,latom) = cr_pos(:,latom) + c_offs(:,i)
            CALL conn_update(latom, FLOAT(c_offs(:,i)))
         ENDIF
!write(*,*) '      search ', latom
!write(*,*) ' INTE T_LIST ', t_list
         CALL mol_add_excl(latom, n_excl, excl)
      ENDIF
   ENDDO search
!write(*,*) ' ENDE T_LIST ', t_list
   DEALLOCATE(c_list)                      ! Remove temporary list of atoms
   DEALLOCATE(c_offs)                      ! Remove temporary offsets
ENDDO
END SUBROUTINE mol_add_excl
!
SUBROUTINE molecularize_numbers(istart,ifinish, new_type, biso)
!
! Groups atoms istart to ifinish into a molecule
! with molecule type new_type. 
!
USE discus_allocate_appl_mod
USE crystal_mod
USE prop_para_mod
USE errlist_mod
use discus_show_menu
IMPLICIT NONE
!
INTEGER, INTENT(IN)  :: istart
INTEGER, INTENT(IN)  :: ifinish
INTEGER, INTENT(IN)  :: new_type
REAL   , INTENT(IN)  :: biso
!
INTEGER   :: n_new,n_atom, n_type, n_mole
INTEGER   :: i, j
character(len=1):: dummy
!
IF(istart  < 1 .OR. istart  > cr_natoms .OR.  &
   ifinish < 1 .OR. ifinish > cr_natoms     ) THEN
   ier_num = -137
   ier_typ = ER_APPL
   WRITE(ier_msg(1),'(a4,i8,a13,2i8)') 'n[1]',cr_natoms,' lower/upper ',istart,ifinish
   RETURN
ENDIF
DO i=istart,ifinish
   IF(btest(cr_prop(i),1)) THEN      ! Atom is already inside a molecule
      ier_num = -138
      ier_typ = ER_APPL
      WRITE(ier_msg(1),'(a9,i8)') 'Atom no. ',i
      RETURN
   ENDIF
ENDDO
!
n_new = ifinish - istart + 1            ! no of atoms in new molecule
!
   n_mole = mole_num_mole + 1           ! Reserve new molecule
   n_type = MAX(mole_num_type,new_type) ! If necessary new molecule type number
                                        ! (offset+length) of previous + current length
   n_atom = mole_off(mole_num_mole) + mole_len(mole_num_mole) + n_new

   IF ( n_mole > MOLE_MAX_MOLE  .or.  &
        n_type > MOLE_MAX_TYPE  .or.  &
        n_atom > MOLE_MAX_ATOM      ) THEN
      n_mole = INT(MAX(n_mole,MOLE_MAX_MOLE)*1.25)        ! Add 25%
      n_type = MAX(n_type,MOLE_MAX_TYPE)
      n_atom = INT(MAX(n_atom+n_new,MOLE_MAX_ATOM)*1.25)  ! reserve space for two more molecules+25%
      CALL alloc_molecule(1, 1, n_mole, n_type, n_atom)
      IF ( ier_num /= 0) THEN
         RETURN
      ENDIF
   ENDIF
   mole_num_mole = mole_num_mole + 1
   mole_num_type = MAX(mole_num_type, new_type)
   mole_len(mole_num_mole)  = n_new
   mole_off(mole_num_mole)  = mole_off(mole_num_mole-1)+mole_len(mole_num_mole-1)
   mole_type(mole_num_mole) = n_type
   mole_char(mole_num_mole) = MOLE_ATOM
   mole_file(mole_num_mole) = ' '
   mole_biso(mole_num_mole) = biso
   mole_num_atom = mole_off(mole_num_mole) + n_new
   j = 0
   DO i=istart,ifinish
         j = j + 1
         mole_cont(mole_off(mole_num_mole)+j) = i
         cr_prop(i) = ibset(cr_prop(i),PROP_MOLECULE)
         cr_mole(i) = mole_num_mole
   ENDDO
!
END SUBROUTINE molecularize_numbers
!
SUBROUTINE molecule_shift(nshift, istart)
!
!   Shifts the atoms in molecule istart and up by nshift atoms further up
!
USE discus_allocate_appl_mod
USE molecule_mod
USE errlist_mod
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: nshift
INTEGER, INTENT(IN) :: istart
!
INTEGER :: i,j,k
INTEGER :: n_mole
INTEGER :: n_type
INTEGER :: n_atom
!
n_atom = MAX(mole_off(mole_num_mole)+mole_len(mole_num_mole)+ nshift, MOLE_MAX_ATOM)
IF(n_atom > MOLE_MAX_ATOM) THEN
   n_mole = MOLE_MAX_MOLE
   n_type = MOLE_MAX_TYPE
   n_atom = INT(MAX(n_atom, MOLE_MAX_ATOM)*1.25)  ! reserve space for two more molecules+25%
   CALL alloc_molecule(1, 1, n_mole, n_type, n_atom)
   IF ( ier_num /= 0) THEN
      RETURN
   ENDIF
ENDIF
!
DO i=mole_num_mole, istart, -1
   DO j=mole_len(i), 1, -1
      k = mole_off(i) + j
      mole_cont(k+nshift) = mole_cont(k)
   ENDDO
   mole_off(i) = mole_off(i) + nshift
ENDDO
IF(istart <= mole_num_mole) THEN
   DO j=1,mole_len(istart)
      mole_cont(mole_off(istart-1)+mole_len(istart-1) + j) = 0
   ENDDO
ENDIF
mole_num_atom = mole_num_atom + nshift
!
END SUBROUTINE molecule_shift
!
!
END MODULE molecule_func_mod
