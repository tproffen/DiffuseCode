MODULE do_molecule_alloc
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE molecule_set_array_size(n_gene, n_symm, n_mole, n_type, n_atom)
!
USE discus_allocate_appl_mod
USE molecule_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: n_gene
INTEGER, INTENT(IN) :: n_symm
INTEGER, INTENT(IN) :: n_mole
INTEGER, INTENT(IN) :: n_type
INTEGER, INTENT(IN) :: n_atom
!
INTEGER :: i_gene
INTEGER :: i_symm
INTEGER :: i_mole
INTEGER :: i_type
INTEGER :: i_atom
!
LOGICAL :: lalloc
!
i_gene = MAX(1,n_gene)
i_symm = MAX(1,n_symm)
i_mole = MAX(1,n_mole)
i_type = MAX(1,n_type)
i_atom = MAX(1,n_atom)
!
IF(n_gene > MOLE_MAX_GENE) THEN
   i_gene = MOLE_MAX_GENE +  5
   lalloc = .true.
ENDIF
!
IF(n_symm > MOLE_MAX_SYMM) THEN
   i_symm = MOLE_MAX_SYMM +  5
   lalloc = .true.
ENDIF
!
IF(n_mole > MOLE_MAX_MOLE) THEN
   i_mole = MOLE_MAX_MOLE + 20
   lalloc = .true.
ENDIF
!
IF(n_type > MOLE_MAX_TYPE) THEN
   i_type = MOLE_MAX_TYPE +  5
   lalloc = .true.
ENDIF
!
IF(n_atom > MOLE_MAX_ATOM) THEN
   i_atom = MOLE_MAX_ATOM + 200
   lalloc = .true.
ENDIF
!
IF(lalloc) THEN

   CALL alloc_molecule(i_gene, i_symm, i_mole, i_type, i_atom)
ENDIF
!
mole_gene_n   = MAX(1,mole_gene_n, n_gene)
mole_symm_n   = MAX(1,mole_symm_n, n_symm)
mole_num_mole = MAX(1,mole_num_mole, n_mole)
mole_num_type = MAX(1,mole_num_type, n_type)
mole_num_atom = MAX(1,mole_num_atom, n_atom)
!
END SUBROUTINE molecule_set_array_size
!
!*******************************************************************************
!
SUBROUTINE molecule_copy_prop(from,to)
!
USE molecule_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: from
INTEGER, INTENT(IN) :: to
!
mole_char(to)  = mole_char(from)
mole_dens(to)  = mole_dens(from)
mole_file(to)  = mole_file(from)
mole_fuzzy(to) = mole_fuzzy(from)
!
END SUBROUTINE molecule_copy_prop
!
END MODULE do_molecule_alloc
