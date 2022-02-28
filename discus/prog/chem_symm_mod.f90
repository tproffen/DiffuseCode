MODULE chem_symm_mod
!
!  Generates all vectors that are symmetrically equivalent 
!  to those that had been provided by the user via a
!  set vec, <no>, isite, jsite, dx, dy, dz
!
USE errlist_mod
!
CONTAINS
!
SUBROUTINE chem_symm(zeile, lp)
!
USE discus_allocate_appl_mod 
USE chem_aver_mod
USE crystal_mod
USE chem_mod
USE wyckoff_mod
!
USE blanks_mod
USE build_name_mod
USE errlist_mod
USE get_params_mod
USE precision_mod
USE support_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: zeile
INTEGER         , INTENT(INOUT) :: lp
!
INTEGER, PARAMETER    :: MAXW = 20
INTEGER, PARAMETER    :: IWR = 91
REAL(kind=PREC_DP)   , PARAMETER    :: EPS =  0.001
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(MAXW) :: cpara
INTEGER,             DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP),  DIMENSION(MAXW) :: werte
!
!CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile)))   :: outfile
INTEGER               :: i, j
INTEGER               :: ianz
INTEGER               :: isym
INTEGER               :: iv, ivv
INTEGER               :: ic
INTEGER               :: iv_max
INTEGER               :: isite, jsite
INTEGER               :: n_vec, n_cor, n_old
INTEGER, DIMENSION(3) :: vec
LOGICAL, DIMENSION(:),ALLOCATABLE :: old_vector
LOGICAL               :: lisite, ljsite
LOGICAL               :: lisnew
LOGICAL               :: lout
REAL(kind=PREC_DP)   , DIMENSION(4) :: ri
REAL(kind=PREC_DP)   , DIMENSION(4) :: rj
REAL(kind=PREC_DP)   , DIMENSION(4) :: ri_n
REAL(kind=PREC_DP)   , DIMENSION(4) :: rj_n
REAL(kind=PREC_DP)   , DIMENSION(4) :: ri_f
REAL(kind=PREC_DP)   , DIMENSION(4) :: rj_f
!
ri(:) = 0.0
rj(:) = 0.0
ri(4) = 1.0
rj(4) = 1.0
lout = .false.
iv_max = 0
!
CALL chem_aver(.FALSE., .TRUE.)        ! Ensure we have the average structure
!
IF(lp>0) THEN                           ! output file parameter given
   CALL get_params(zeile, ianz, cpara, lpara, MAXW, lp)
   IF(ier_num==0) THEN
      CALL do_build_name(ianz, cpara, lpara, werte, MAXW, 1)
!   outfile = zeile(1:lp)
      CALL oeffne(IWR, cpara(1), 'unknown')
      IF(ier_num/=0) RETURN
      lout = .TRUE.
   ELSE
      RETURN
   ENDIF
ENDIF
!
DO iv = 1, CHEM_MAX_VEC    ! Loop over all defined vectors
   IF(chem_cvec(1,iv) /= -9999 ) THEN   ! Vector has been defined
      iv_max = iv                       ! Log highest vector number
!     If file was defined, write output file
      IF(lout) WRITE(IWR,1000) iv, chem_cvec(:,iv), iv
   ENDIF
ENDDO
!
!  Create a list of all previously used correlation vectors, 
!  only these will be copied via symmetry expansion
!
ALLOCATE(old_vector(1:CHEM_MAX_VEC))   ! space for old vectors
old_vector(:) = .false.
DO iv = 1, CHEM_MAX_VEC                ! Loop over all defined vectors
   IF(chem_cvec(1,iv) /= -9999 ) THEN
      old_vector(iv) = .true.          ! 
   ENDIF
ENDDO
n_old = CHEM_MAX_VEC
!
ivv = iv_max
!
symmetry: DO isym = 2, spc_n   !Loop over all space group symmetry operations
!   IF(lout) WRITE(IWR,1100)
!
   old_vectors: DO iv = 1, n_old           ! Loop over all defined vectors
      IF(old_vector(iv)           ) THEN   ! Vector has been defined
         DO i=1,3                          ! copy sites from average unit cell
            ri(i) = chem_ave_pos(i,chem_cvec(1,iv))
            rj(i) = chem_ave_pos(i,chem_cvec(2,iv)) + REAL(chem_cvec(2+i,iv))
         ENDDO
         ri_n(:) = 0.0                     ! initialize new vectors
         rj_n(:) = 0.0
         DO i=1,4                          ! apply space group symmetry
            DO j=1,4
              ri_n(i) = ri_n(i) + spc_mat(i,j,isym)*ri(j)
              rj_n(i) = rj_n(i) + spc_mat(i,j,isym)*rj(j)
            ENDDO
         ENDDO
!write(*,*) ' old sites ', ri(1:3), rj(1:3)
         DO i=1,3                          ! convert sites into first unit cell [0:1]
            ri_f(i) = ri_n(i) - REAL(INT(ri_n(i))) + 1.0 ! create a copy for site i
            ri_f(i) = ri_f(i) - REAL(INT(ri_f(i))) 
            rj_f(i) = rj_n(i) - REAL(INT(rj_n(i))) + 1.0 ! create a copy for site j
            rj_f(i) = rj_f(i) - REAL(INT(rj_f(i))) 
            vec(i)  = nint(rj_n(i)-rj_f(i)) - nint(ri_n(i)-ri_f(i))
         ENDDO
!write(*,*) ' new sites ', ri_f(1:3), rj_f(1:3)
!read(*,*) i
!
!        To identify the new sites within the unit cell, we have to compare the new
!        position to all atoms in the average unit cell. Multiple hits flag an error
!
         isite = 0
         jsite = 0
         lisite = .false.
         ljsite = .false.
         DO i=1, cr_ncatoms                ! Identify new site in unit cell
            IF(ABS(ri_f(1)-chem_ave_pos(1,i))<EPS .AND. &
               ABS(ri_f(2)-chem_ave_pos(2,i))<EPS .AND. &
               ABS(ri_f(3)-chem_ave_pos(3,i))<EPS      ) THEN
               isite = i
               IF(lisite) THEN
                  ier_num = -30
                  ier_typ = ER_CHEM
                  ier_msg(1) = 'There are several sites in the unit cell with'
                  ier_msg(2) = '(almost) identical coordinates. Symmetry ex-  '
                  ier_msg(3) = 'pansion restricted to file, flagged as -site'
                  isite = -isite
               ENDIF
               lisite = .true.
            ENDIF
            IF(ABS(rj_f(1)-chem_ave_pos(1,i))<EPS .AND. &
               ABS(rj_f(2)-chem_ave_pos(2,i))<EPS .AND. &
               ABS(rj_f(3)-chem_ave_pos(3,i))<EPS      ) THEN
               jsite = i
               IF(ljsite) THEN
                  ier_num = -30
                  ier_typ = ER_CHEM
                  ier_msg(1) = 'There are several sites in the unit cell with'
                  ier_msg(2) = '(almost) identical coordinates. Symmetry ex-'
                  ier_msg(3) = 'pansion restricted to file, flagged as -site'
                  jsite = -jsite
               ENDIF
               ljsite = .true.
            ENDIF
         ENDDO
!
         lisnew = .TRUE.
         search: DO i=1, ivv
            IF(chem_cvec(1,i)==isite .AND. chem_cvec(2,i)==jsite .AND.  &
               chem_cvec(3,i)==vec(1).AND. chem_cvec(4,i)==vec(2).AND.  &
               chem_cvec(5,i)==vec(3) ) THEN
               lisnew = .FALSE.
               EXIT search
           ENDIF
         ENDDO search
!write(*,*) ' NEW VECTOR ', isite, jsite, vec, iv, lisnew
         IF(lisnew) THEN
         ivv = ivv + 1
         IF(lout) WRITE(IWR,1000) ivv, isite, jsite, vec, iv
!
! update list of defined correlation vectors
! and enter these new vectors into neighborhood definitions
!   Loops over all correlations and checks if it uses vectors
!
         IF(isite>0 .AND. jsite>0) THEN    !no error
            IF (ivv > CHEM_MAX_VEC) THEN
               n_vec = CHEM_MAX_VEC + 10   ! Increment size by 10
               n_cor = CHEM_MAX_COR
               CALL alloc_chem_vec ( n_vec , n_cor )
               IF (ier_num /= 0) RETURN
            ENDIF
            chem_cvec(1,ivv) = isite
            chem_cvec(2,ivv) = jsite
            chem_cvec(3:5,ivv) = vec(1:3)
            DO ic=1,chem_ncor
               IF(chem_ctyp (ic) .eq.CHEM_VEC) THEN
is_used:          DO i=1, chem_nvec (ic)
                     IF(chem_use_vec(i,ic)==iv) THEN
                        chem_use_vec(chem_nvec(ic)+1,ic) = ivv
                        chem_nvec(ic) = chem_nvec(ic) + 1
                        EXIT is_used
                     ENDIF
                  ENDDO is_used
               ENDIF
            ENDDO
         ENDIF
         ENDIF
!
      ENDIF ! old_vector(iv) == .true.
   ENDDO old_vectors
ENDDO symmetry
!
WRITE(IWR,'(a)') 'set neig, rese' 
DO ic=1,chem_ncor
   IF(chem_ctyp(ic) == CHEM_VEC) THEN
      WRITE(IWR,'(a,50(a2,i3))') 'set neig, vec', (', ',chem_use_vec(i,ic),i=1, chem_nvec(ic))
   ELSEIF(chem_ctyp(ic) == CHEM_CON) THEN
      WRITE(IWR,'(a,50(a2,i3))') 'set neig, con', (', ',chem_use_con(i,ic),i=1, chem_ncon(ic))
   ELSEIF(chem_ctyp(ic) == CHEM_ANG) THEN
      WRITE(IWR,'(a,50(a2,i3))') 'set neig, ang', (', ',chem_use_win(i,ic),i=1, chem_nwin(ic))
   ELSEIF(chem_ctyp(ic) == CHEM_RANGE) THEN
      WRITE(IWR,'(a,50(a2,i3))') 'set neig, ra ', (', ',chem_use_ran(i,ic),i=1, chem_nran(ic))
   ELSEIF(chem_ctyp(ic) == CHEM_ENVIR) THEN
      WRITE(IWR,'(a,50(a2,i3))') 'set neig, env', (', ',chem_use_env(i,ic),i=1, chem_nenv(ic))
!  ELSEIF(chem_ctyp(ic) == CHEM_DIST ) THEN
!     WRITE(IWR,'(a,50(a2,i3))') 'set neig, dis', (', ',chem_use_env(i,ic),i=1, chem_nenv(ic))
!  ELSEIF(chem_ctyp(ic) == CHEM_DIR  ) THEN
!     WRITE(IWR,'(a, 6(a2,f9.4))') 'set neig, dir', (', ',chem_dir(1:3,1,ic), chem_dir(1:3,2,ic)
   ENDIF
   IF(ic<chem_ncor) WRITE(IWR,'(a)') 'set neig, add'
ENDDO
!
IF(lout) CLOSE (IWR)
DEALLOCATE(old_vector)
!
1000 FORMAT('set vec, ', 5(i3,', '),i3,'  ! symm. equiv. to vector ', i3) 
1100 FORMAT('!',/,'!')
!
END SUBROUTINE chem_symm
!
END MODULE chem_symm_mod
