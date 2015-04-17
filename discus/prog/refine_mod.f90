MODULE refine_mod
!
USE crystal_mod
USE rmc_mod
!
PRIVATE
PUBLIC refine_alloc
PUBLIC refine_dealloc
PUBLIC refine_adapt_move
PUBLIC refine_adapt_lattice
PUBLIC refine_adapt_pdf_sd
!
PUBLIC ref_maxlatt
PUBLIC ref_maxpdfsd
!
INTEGER, PARAMETER                   :: REF_MV_INT=100
INTEGER, DIMENSION(:,:), ALLOCATABLE :: ref_move_success
INTEGER, DIMENSION(  :), ALLOCATABLE :: ref_move_success_n
INTEGER, PARAMETER                   :: REF_LAT_INT=30
INTEGER, DIMENSION(:  ), ALLOCATABLE :: ref_latt_success
INTEGER                              :: ref_latt_success_n
INTEGER, PARAMETER                   :: REF_PDF_SD =30     !PDF Scale Density
INTEGER, DIMENSION(:  ), ALLOCATABLE :: ref_pdfs_success
INTEGER                              :: ref_pdfs_success_n
INTEGER                              :: isuccess
REAL   , DIMENSION(2)                :: ref_maxlatt
REAL   , DIMENSION(2)                :: ref_maxpdfsd
!
CONTAINS
!
SUBROUTINE refine_alloc
!
!  allocates the refinement arrays
!
IF(.not.ALLOCATED(ref_move_success)) THEN
   ALLOCATE(ref_move_success( 0:REF_MV_INT-1,0:MAXSCAT))
   ALLOCATE(ref_move_success_n (             0:MAXSCAT))
   ref_move_success   = 0
   ref_move_success_n = 0
ENDIF
!
IF(.not.ALLOCATED(ref_latt_success)) THEN
   ALLOCATE(ref_latt_success( 0:REF_LAT_INT-1))
   ref_latt_success   = 0
   ref_latt_success_n = 0
ENDIF
ref_maxlatt(1) = 0.003
ref_maxlatt(2) = 0.030
!
IF(.not.ALLOCATED(ref_pdfs_success)) THEN
   ALLOCATE(ref_pdfs_success( 0:REF_LAT_INT-1))
   ref_pdfs_success   = 0
   ref_pdfs_success_n = 0
ENDIF
ref_maxpdfsd(1) = 0.002
ref_maxpdfsd(2) = 0.0004
!
END SUBROUTINE refine_alloc
!
SUBROUTINE refine_dealloc
!
!  deallocates the refinement arrays
!
IF(ALLOCATED(ref_move_success)) THEN
   DEALLOCATE(ref_move_success  )
   DEALLOCATE(ref_move_success_n)
ENDIF
!
IF(ALLOCATED(ref_latt_success)) THEN
   DEALLOCATE(ref_latt_success)
ENDIF
!
IF(ALLOCATED(ref_pdfs_success)) THEN
   DEALLOCATE(ref_pdfs_success)
ENDIF
!
END SUBROUTINE refine_dealloc
!
SUBROUTINE refine_adapt_move(natoms, isel, spin )
!
!  Adapts the sigma for position moves
!
IMPLICIT NONE
!
INTEGER,                           INTENT(IN) :: natoms
INTEGER, DIMENSION(1:rmc_max_atom),INTENT(IN) :: isel
INTEGER,                           INTENT(IN) :: spin
!

INTEGER  :: i,j,k    ! dummy loop variables
!
DO i = 1, natoms 
   ref_move_success_n (cr_iscat(isel(i))) =   &
   MOD(ref_move_success_n (cr_iscat(isel(i))) + 1, REF_MV_INT)
   j = ref_move_success_n (cr_iscat(isel(i)))
   ref_move_success(j,cr_iscat(isel(i))) = spin
!write(*,*) ' ADAP ', j, ref_move_success(j,cr_iscat(isel(i))), rmc_maxmove(1,j)
   IF(MOD(j, REF_MV_INT/10)==REF_MV_INT/10-1) THEN
      isuccess = 0
      DO k=0,REF_MV_INT-1
         isuccess = isuccess + ref_move_success(k,cr_iscat(isel(i)))
      ENDDO
!write(*,*) ' IN ADAPTATION +', isuccess, REF_MV_INT/5
      IF(isuccess < REF_MV_INT/5) THEN
         rmc_maxmove(:,cr_iscat(isel(i))) = rmc_maxmove(:,cr_iscat(isel(i))) * 0.85
      ELSEIF(isuccess > REF_MV_INT/5) THEN
         rmc_maxmove(:,cr_iscat(isel(i))) = rmc_maxmove(:,cr_iscat(isel(i))) / 0.85
      ENDIF
   ENDIF
ENDDO
!
END SUBROUTINE refine_adapt_move
!
SUBROUTINE refine_adapt_lattice(spin )
!
!  Adapts the sigma for position moves
!
IMPLICIT NONE
!
INTEGER,                           INTENT(IN) :: spin
!

INTEGER  :: j,k    ! dummy loop variables
!
   ref_latt_success_n =   &
   MOD(ref_latt_success_n + 1, REF_LAT_INT)
   j = ref_latt_success_n 
   ref_latt_success(j) = spin
!write(*,*) ' ADAP ', j, ref_latt_success(j), ref_maxlatt(1)
   IF(MOD(j, REF_LAT_INT/10)==REF_LAT_INT/10-1) THEN
      isuccess = 0
      DO k=0,REF_LAT_INT-1
         isuccess = isuccess + ref_latt_success(k)
      ENDDO
!write(*,*) ' IN ADAPTATION +', isuccess, REF_LAT_INT/5, (isuccess < REF_LAT_INT/5)
      IF(isuccess < REF_LAT_INT/5) THEN
         ref_maxlatt(:) = ref_maxlatt(:) * 0.85
      ELSEIF(isuccess > REF_LAT_INT/5) THEN
         ref_maxlatt(:) = ref_maxlatt(:) / 0.85
      ENDIF
   ENDIF
!
END SUBROUTINE refine_adapt_lattice
!
SUBROUTINE refine_adapt_pdf_sd(spin )
!
!  Adapts the sigma for PDF scale and density
!
IMPLICIT NONE
!
INTEGER,                           INTENT(IN) :: spin
!

INTEGER  :: j,k    ! dummy loop variables
!
   ref_pdfs_success_n =   &
   MOD(ref_pdfs_success_n + 1, REF_PDF_SD)
   j = ref_pdfs_success_n 
   ref_pdfs_success(j) = spin
!write(*,*) ' ADAP ', j, ref_pdfs_success(j), ref_maxpdfsd
   IF(MOD(j, REF_PDF_SD/10)==REF_PDF_SD/10-1) THEN
      isuccess = 0
      DO k=0,REF_PDF_SD-1
         isuccess = isuccess + ref_pdfs_success(k)
      ENDDO
!write(*,*) ' IN ADAPTATION +', isuccess, REF_PDF_SD/5, (isuccess < REF_PDF_SD/5)
      IF(isuccess < REF_PDF_SD/5) THEN
         ref_maxpdfsd(:) = ref_maxpdfsd(:) * 0.85
      ELSEIF(isuccess > REF_PDF_SD/5) THEN
         ref_maxpdfsd(:) = ref_maxpdfsd(:) / 0.85
      ENDIF
   ENDIF
!
END SUBROUTINE refine_adapt_pdf_sd
!
END MODULE refine_mod
