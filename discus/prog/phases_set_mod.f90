MODULE phases_set_mod
!
!  Set values for the phases
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE phases_set(zeile,  lcomm)
!
USE crystal_mod
USE discus_allocate_appl_mod
USE phases_mod
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE precision_mod
USE take_param_mod
!
IMPLICIT NONE
!                                                                       
CHARACTER(LEN=*), INTENT(INOUT) :: zeile
INTEGER         , INTENT(INOUT) :: lcomm
!
INTEGER, PARAMETER :: MAXW = 7
!                                                                       
INTEGER, PARAMETER :: NOPTIONAL = 3
INTEGER, PARAMETER :: O_CURRENT = 1                  ! Current phase number
INTEGER, PARAMETER :: O_WGHT    = 2                  ! Weight fraction
INTEGER, PARAMETER :: O_MODE    = 3                  ! Single / multiple
CHARACTER(LEN=   8), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 2 ! Number of values to calculate 
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
INTEGER             :: ianz
!
DATA oname  / 'current ', 'fraction', 'mode    ' /
DATA loname /  8        ,  8        ,  8         /
!
INTEGER :: n_pha
INTEGER :: n_pts
INTEGER :: n_scat
!
!                                                                       
opara  =  (/ '1.0000', '1.0000', 'single' /)   ! Always provide fresh default values
lopara =  (/  6      ,  6      ,  6       /)
owerte =  (/  1.0    ,  1.0    ,  0.0     /)
!
CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm)
err_para: IF (ier_num.eq.0) THEN
   CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                     oname, loname, opara, lopara, lpresent, owerte)
   err_opti: IF (ier_num.eq.0) THEN
      IF(lpresent(O_MODE)) THEN                ! mode: is present
         IF(opara(O_MODE)=='single') THEN
            pha_multi = .FALSE.
            pha_n     = 1
            pha_curr  = 1
            n_pha  = 1
            n_pts  = 1
            n_scat = 1
            CALL alloc_phases(n_pha, n_pts, n_scat)
            pha_frac(:)    = 0.0
            pha_frac(1)    = 0.0
            pha_nscat(1)   = 0.0
         ELSEIF(opara(O_MODE)=='multiple') THEN
            pha_multi = .TRUE.
            IF(lpresent(O_CURRENT) .AND. lpresent(O_WGHT)) THEN
               pha_curr = NINT(owerte(O_CURRENT))
               IF(pha_curr>PHA_MAXPHA) THEN
                  n_pha  = pha_curr
                  n_pts  = PHA_MAXPTS
                  n_scat = PHA_MAXSCAT
                  CALL alloc_phases(n_pha, n_pts, n_scat)
               ENDIF
               pha_frac(pha_curr) = owerte(O_WGHT)
               pha_nscat(pha_curr) = 0.0
               pha_n = MAX(pha_n, pha_curr)
            ELSEIF((     lpresent(O_CURRENT) .AND. .NOT.lpresent(O_WGHT)) .OR. &
                   (.NOT.lpresent(O_CURRENT) .AND.      lpresent(O_WGHT))      &
                  ) THEN
               ier_num = -6
               ier_typ = ER_APPL
               ier_msg(1) = '''current'' and ''fraction'' must both be present'
               ier_msg(2) = 'or both be absent'
            ENDIF
         ELSE
            ier_num = -6
            ier_typ = ER_COMM
            ier_msg(1) = 'Mode requires ''single'' or ''multiple'''
         ENDIF
      ELSE                                     ! mode: is absent
         IF(pha_multi) THEN                    ! We are in multiple mode
            IF(lpresent(O_CURRENT) .AND. lpresent(O_WGHT)) THEN
               pha_curr = NINT(owerte(O_CURRENT))
               IF(pha_curr>PHA_MAXPHA) THEN
                  n_pha  = pha_curr
                  n_pts  = PHA_MAXPTS
                  n_scat = PHA_MAXSCAT
                  CALL alloc_phases(n_pha, n_pts, n_scat)
               ENDIF
               pha_frac(pha_curr) = owerte(O_WGHT)
               pha_nscat(pha_curr) = 0.0
               pha_n = MAX(pha_n, pha_curr)
            ELSEIF((     lpresent(O_CURRENT) .AND. .NOT.lpresent(O_WGHT)) .OR. &
                   (.NOT.lpresent(O_CURRENT) .AND.      lpresent(O_WGHT))      &
                  ) THEN
               ier_num = -6
               ier_typ = ER_APPL
               ier_msg(1) = '''current'' and ''fraction'' must both be present'
               ier_msg(2) = 'or both be absent'
            ENDIF
         ELSE                                  ! We are in single mode
            IF(lpresent(O_CURRENT) .OR. lpresent(O_WGHT)) THEN
               ier_num = -6
               ier_typ = ER_APPL
               ier_msg(1) = 'Phases are in single mode, activate with'
               ier_msg(2) = 'mode:multiple                           '
            ENDIF
         ENDIF
      ENDIF                                    ! mode: is present
   ENDIF err_opti
ENDIF err_para
!
END SUBROUTINE phases_set
!
!*******************************************************************************
!
SUBROUTINE phases_place
!-
! Place the current powder pattern into the phases entry
!+
!
USE crystal_mod
USE crystal_task_mod
USE chem_mod
USE debye_mod
USE diffuse_mod
USE pdf_mod
USE phases_mod
USE powder_mod
USE powder_tables_mod
!
USE precision_mod
!
IMPLICIT NONE
!
integer :: i,j
INTEGER :: k,iscat
INTEGER :: npkt
REAL( KIND(0.0D0))             :: signum
!
npkt = NINT((pow_qmax-pow_qmin)/pow_deltaq) + 1
stack: IF(pow_four_mode==POW_FOURIER) THEN     ! Standard Fourier, not Stacking fault mode
!write(*,*) ' PHASES ', pha_multi, pha_n, pha_curr
!write(*,*) ' entry  ', lbound(pha_entry), ubound(pha_entry)
!write(*,*) ' frac   ', lbound(pha_frac ), ubound(pha_frac )
!write(*,*) ' powder ', lbound(pha_powder), ubound(pha_powder), num(1)*num(2)
!write(*,*) ' adp    ', lbound(pha_adp   ), ubound(pha_adp   ), cr_nscat
   CALL crystal_calc_mass
!write(*,*) ' pha_frac   ', lbound(pha_frac), ubound(pha_frac), pha_frac
!write(*,*) ' pha_weight ', lbound(pha_weight), ubound(pha_weight), pha_weight
   pha_weight(pha_curr) = cr_mass            ! Mass in multiples of u for this phase
   if(pow_l_partial) then
      loop_partial: do i=1,cr_nscat
         if(pow_do_partial(i,i)) then
            pha_nscat(pha_curr) = 1
            exit loop_partial
         else
            do j=i,cr_nscat
              if(pow_do_partial(i,j)) then
                 pha_nscat(pha_curr) = 2
                 exit loop_partial
              endif
            enddo
         endif
      enddo loop_partial
   else
      pha_nscat(pha_curr)  = cr_nscat        ! Number of atom types for this phase
   end if
   pha_calc (pha_curr)  = pow_four_type      ! Fourier type Complete / Debye
!write(*,*) ' Crystal mass ', cr_mass
!write(*,*) ' SET FORM ', num(1), num(2), num(1)*num(2), npkt
   DO iscat = 1, cr_nscat
      DO k=0, npkt
         signum = 1.0D0
         IF(REAL(cfact_pure(1, iscat))< 0.0D0) signum = -1.0D0
         pha_form(k, iscat, pha_curr) = SQRT( DBLE (       cfact_pure (powder_istl (k), iscat)   *   &
                                                    CONJG (cfact_pure (powder_istl (k), iscat) )  )  &
                                            )                                                        &
                                        * signum
      ENDDO
      pha_adp  (iscat, pha_curr) = cr_dw(iscat)
      pha_occ  (iscat, pha_curr) = cr_occ(iscat)
      pha_niscat(iscat, pha_curr) = cr_niscat(iscat)    ! cr_niscat is the number of atoms of type iscat
!write(*,*) pha_adp  (iscat, pha_curr), pha_occ  (iscat, pha_curr), &
!                                      pha_niscat(iscat, pha_curr)
   ENDDO
   pha_nreal(pha_curr)  = SUM(pha_niscat(1:pha_nscat(pha_curr),pha_curr)* &
                              pha_occ(1:pha_nscat(pha_curr),pha_curr)     &
                             )  ! number of real atoms in phase 
   IF(chem_quick) THEN
      pha_ncreal(pha_curr) = pha_nreal(pha_curr)/(cr_icc(1)*cr_icc(2)*cr_icc(3))
   ELSE
      pha_ncreal(pha_curr) = cr_ncatoms
   ENDIF
ENDIF stack
!write(*,*) ' IN PHA_set    nreal, ', pha_nreal(pha_curr), pha_ncreal(pha_curr), cr_nscat
!write(*,*) ' PHASES_SET ', pow_conv(0), pow_conv(1), pow_conv(npkt)
!write(*,*) ' PHASES_SET ', minval(pow_conv), maxval(pow_conv)
!
! Place powder pattern into appropriate phase entry
!
IF (pow_four_type.eq.POW_COMPL) THEN                 ! Complete powder patterm, normalizer is 1
!write(*,*) ' COPY ', npkt
!open(66, file='phases.place', status='unknown')
   DO k=0, npkt
      pha_powder(k,pha_curr) = pow_conv(k)                          &
                    /(REAL(pha_nreal(pha_curr), KIND=PREC_DP))**2 * &
                      REAL(pha_ncreal(pha_curr) , KIND=PREC_DP)   / &
                      REAL(cr_v, KIND=PREC_DP)                    * &
             20.00D0  ! Needs to be verified where this 20 comes from
!  pow_tmp = pow_tmp/(REAL(pow_nreal, KIND=PREC_DP))**2 * &
!                     REAL(pow_ncreal , KIND=PREC_DP)/REAL(cr_v, KIND=PREC_DP) * &
!            20.00D0  ! Needs to be verified where this 20 comes from
!q = pow_qmin + (k-1)*pow_deltaq
!  write(66, '(3F20.6)') q, pow_conv(k), pha_powder(k,pha_curr)
   ENDDO
!close(66)
ELSE                                                 ! Complete powder patterm, normalizer is nreal
!write(*,*) ' COPY ', npkt, num(1)*num(2), num(1:2)
!open(66, file='phases.place', status='unknown')
   DO k=0, npkt
      pha_powder(k,pha_curr) = pow_conv(k) / REAL(pha_nreal(pha_curr))
!write(66, '(3F20.6)') pow_qmin + (k-1)*pow_deltaq, pow_conv(k), pha_powder(k,pha_curr)
   ENDDO
!close(66)
ENDIF
!write(ofile,'(a,i1.1)') 'phases.place.',pha_curr
!open(66, file=ofile, status='unknown')
!   DO k=0, npkt
!write(66, '(3F20.6)') pow_qmin + (k-1)*pow_deltaq, pow_conv(k), pha_powder(k,pha_curr)
!  ENDDO
!close(66)
!
! Adjust powder pattern if clin or cquad are non-zero
!
! Currently the sine fourier is not reversible need to find the issue at hand
!IF(pdf_clin_a/=0.0 .OR. pdf_cquad_a/=0.0) THEN
!  CALL phases_corr(npkt)
!ENDIF
!
END SUBROUTINE phases_place
!
!*******************************************************************************
!
SUBROUTINE phases_average(xmin, xdel, npkt)
!-
!  Average the powder intensities and form factors for all phases
!
!  The required scale is the (intended weight fraction) / (current weight fraction)
!  current weight fraction = (pha_weight) / (grand weight)
!+
!
USE crystal_mod
USE debye_mod
USE diffuse_mod
USE pdf_mod
USE phases_mod
USE powder_mod
!
USE errlist_mod
USE precision_mod
USE trig_degree_mod
USE wink_mod
!
IMPLICIT NONE
!
REAL(KIND=PREC_SP), INTENT(IN) :: xmin      ! Smallest Q-value
REAL(KIND=PREC_SP), INTENT(IN) :: xdel      ! Q-step
INTEGER           , INTENT(IN) :: npkt      ! Number of points in powder pattern
!
REAL(KIND=PREC_SP), PARAMETER :: EPS = 0.001  ! A small fraction
INTEGER            :: i,j, k
INTEGER            :: empty                 ! number of empty phases with no atoms
REAL(KIND=PREC_SP) :: weight                ! Current grand weight
REAL(KIND=PREC_SP) :: fractions             ! Current grand fractions
REAL(KIND=PREC_SP) :: q                     ! Temporary number of real atoms per phase
REAL(KIND=PREC_SP) :: ttheta                ! 2Theta
!

pow_f2aver(:)    = 0.0D0
pow_faver2(:)    = 0.0D0
pow_fu    (:)    = 0.0D0
pow_f2    (:,:)  = 0.0D0
pow_u2aver       = 0.0
!
weight = SUM(pha_weight(1:pha_n))           ! Calculate grand weight
fractions = SUM(pha_frac(1:pha_n))          ! Sum up all phase fractions
empty = 0
DO i=1,pha_n
   IF(pha_nreal(i)==0.0) THEN
      empty = empty+1
   ENDIF
ENDDO
!write(*,*) ' pha_frac   ', lbound(pha_frac), ubound(pha_frac), pha_frac
!write(*,*) ' pha_weight ', lbound(pha_weight), ubound(pha_weight), pha_weight
!write(*,*) ' pha_n      ', pha_n, empty, pha_nreal(1:pha_n)
IF(ABS(1.-fractions) > EPS) THEN
   ier_num = -169
   ier_typ = ER_APPL
   WRITE(ier_msg(1),'(a,f6.3)') 'Fractions add up to ',fractions
   RETURN
ENDIF
!
! Scale = (intended fraction)/(current weight fraction) / number of phases
pha_scale = pha_frac * weight / pha_weight / (pha_n-empty)  
!
! Add all form factors into pow_faver2
!
!write(*,*) ' PHA_AVERAGE ', pha_n, pha_nscat(1), pha_niscat(:,1), pha_nreal(1)
!write(*,*) ' xmin, xdel, npkt', xmin, xdel, npkt
!
DO i=1,pha_n                                ! Sum over all phases
   IF(pha_nreal(i)>0.0) THEN
   DO j=1,pha_nscat(i)                      ! Sum over the entries for each phase
!write(*,*) ' ', pha_form(1,j,i), pha_niscat(j,i),pha_occ(j,i), pha_nreal(i), pha_scale(i), cr_dw(j)
      DO k = 0, npkt
         q = (k*xdel + xmin)
         pow_faver2(k) = pow_faver2(k) + &
            pha_form(k,j,i)   *          &
            REAL(pha_niscat(j,i)*pha_occ(j,i)/pha_nreal(i)*pha_scale(i),KIND=PREC_DP)
         pow_f2aver(k) = pow_f2aver(k) + &
            pha_form(k,j,i)**2*          &
            REAL(pha_niscat(j,i)*pha_occ(j,i)/pha_nreal(i)*pha_scale(i),KIND=PREC_DP)
         pow_fu(k) = pow_fu(k) +                                               &
            pha_form(k,j,i)**2*                                                &
            REAL(pha_occ(j,i)*pha_scale(i)*EXP(-q**2/8./PI**2*pha_adp(j,i)),KIND=PREC_DP)  &
           *REAL(pha_niscat(j,i)/pha_nreal(i))                                
      ENDDO
!write(*,*) 'u2aver ', i,j, pha_adp(j,i), pha_niscat(j,i), pha_nreal(i), pha_occ(j,i), pha_frac(i),pha_scale(i)
      pow_u2aver = pow_u2aver + (pha_adp(j,i)) * pha_niscat(j,i)/pha_nreal(i)*pha_occ(j,i)*pha_scale(i)
   ENDDO
   ENDIF
ENDDO
pow_faver2(:) = pow_faver2(:)**2
pow_u2aver    = pow_u2aver /8./REAL(PI)**2
!write(*,*) ' U2aver ', pow_u2aver, pow_u2aver*8*PI**2
!write(*,*) ' f2aver ', pow_f2aver(0), pow_f2aver(1)
!write(*,*) ' faver2 ', pow_faver2(0), pow_faver2(1)
!write(*,*) ' Weight ', pha_weight(1:pha_n), weight
!write(*,*) ' Fract  ', pha_frac(1:pha_n), fractions
!write(*,*) ' Scale  ', pha_scale(1:pha_n)
!write(*,*) ' Wcor   ', (pha_weight(i)*pha_scale(i),i=1,pha_n)
!open(77,file='POWDER/phases.faver', status='unknown')
!do k=0, npkt
!  q = ((k)*xdel + xmin)
!  write(77, '(8F12.6)') q, pow_fu(k), pha_form(k,1:pha_nscat(1),1)**2, pow_faver2(k) 
!enddo
!close(77)
!
! Add all powder pattern, each is multiplied with pha_scale
! Distinguish if calculation mode was COMPLETE or DEBYE
! Make separate results for intensities and S(Q) [F(Q) and PDF will work with S(Q)]
!
pow_conv = 0.0
pow_sq   = 0.0
DO i=1,pha_n
   IF(pha_nreal(i)>0.0) THEN
   IF(pha_calc(i) == POW_COMPL) THEN           ! Complete calculation mode
!
!     Prepare intensity output
!
      DO k=0, npkt
         q = (k*xdel + xmin)
         ttheta = 2.*asind ( REAL(q / 2. /REAL(zpi) *rlambda ))
         pow_conv(k) = pow_conv(k) + pha_scale(i)*pha_powder(k,i) / q**2 &
                                    * polarisation(ttheta)
      ENDDO
!
!     Prepare S(Q)      output
!
      DO k=0, npkt
         q = (k*xdel + xmin)
         pow_sq(k) = pow_sq(k) + pha_scale(i)* pha_powder(k,i) / q**2
      ENDDO
!
   ELSEIF(pha_calc(i) == POW_DEBYE) THEN           ! Complete calculation mode
!
!     Prepare intensity output
!
      DO k=0, npkt
         q = (k*xdel + xmin)
         ttheta = 2.*asind ( REAL(q / 2. /REAL(zpi) *rlambda ))
         pow_conv(k) = pow_conv(k) + pha_scale(i)*pha_powder(k,i) &
                                    * polarisation(ttheta)
      ENDDO
!
!     Prepare S(Q)      output
!
      DO k=0, npkt
         pow_sq(k) = pow_sq(k) + pha_scale(i)* pha_powder(k,i)
      ENDDO
!
   ENDIF
   ENDIF
ENDDO
!open(87,file='POWDER/multi_average.conv1', status='unknown')
!do k= 0, npkt
!q = ((k)*xdel + xmin)
!write(87,'(2(f16.6,2x))') q,      pow_conv(k), pow_sq(k)
!enddo
!close(87)
!
!write(*,*) ' IN PHASES ', deb_conv .OR. .NOT.ldbw, deb_conv , .NOT.ldbw
IF(deb_conv .OR. .NOT.ldbw) THEN              ! DEBYE was done with convolution of ADP
!write(*,*) ' DEBYE was     done  with conv '
   DO k=0, npkt
      pow_sq(k) = pow_sq(k) / REAL(pow_faver2(k))                               &
                  + 1.0 - pow_f2aver(k)/pow_faver2(k)
!        ypl(j) =  (ypl(j)/REAL(pow_faver2(j))/normalizer    &
!                  + 1.0 - & !exp(-q**2*pow_u2aver*0.25)      * &
!                   pow_f2aver(j)/pow_faver2(j))
   ENDDO
ELSE
!write(*,*) ' DEBYE was not done  with conv '
! .NOT. (pdf_clin_a/=0.0 .OR. pdf_cquad_a/=0.0), &
! pdf_clin_a/=0.0, pdf_cquad_a/=0.0
   IF(.NOT. (pdf_clin_a/=0.0 .OR. pdf_cquad_a/=0.0)) THEN
   DO k=0, npkt
      q = ((k)*xdel + xmin)
      pow_conv(k) = pow_conv(k) +                                               &
                    (+ 1.0 - exp(-q**2*pow_u2aver))*pow_faver2(k)
   ENDDO
   ENDIF
   DO k=0, npkt
      pow_sq(k) = 1.0 + (pow_sq(k) - pow_fu(k))/pow_faver2(k) 
!        ypl(j) = 1.0 + (ypl(j) -pow_fu(j)                     )/pow_faver2(j)
   ENDDO
ENDIF
!
!open(87,file='POWDER/multi_average.conv', status='unknown')
!do k= 0, npkt
!q = ((k)*xdel + xmin)
!write(87,'(2(f16.6,2x))') q,      pow_conv(k), pow_sq(k)
!enddo
!close(87)
!write(*,*) ' PHASES_FIN ', pha_powder(0, 1), pha_powder(1, 1), pha_powder(4300, 1), pha_powder(npkt, 1)
!write(*,*) ' PHASES_FIN ', pow_conv(0), pow_conv(1), pow_conv(4300), pow_conv(npkt)
!write(*,*) ' PHASES_FIN ', pow_sq(0), pow_sq(1), pow_sq(4300), pow_sq(npkt)
!
END SUBROUTINE phases_average
!
!*******************************************************************************
!
SUBROUTINE phases_corr(npkt)
!-
! correct the powder pattern in case of correlated motion pdf_clin / pdf_cquad
!
! Preparatory steps:
!   Needs faver2, f2aver, u2aver  for this phase
!
! Necessary steps are First forward, then backward:
!
!   Divide by an artificial Debye-Waller term
!   Calculate S(Q)/F(Q)/PDF
!   Convolute PDF with distacne dependent Gaussian
!   Transform back to F(Q) ==> S(Q) == Inte
!+
USE crystal_mod
USE debye_mod
USE diffuse_mod
USE output_mod
USE pdf_mod
USE phases_mod
USE powder_mod
USE powder_fft_mod
!
USE precision_mod
USE wink_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: npkt
!
INTEGER            :: npkt_fft    ! number of points in powder pattern for Fast Fourier
INTEGER            :: j,k, ii
INTEGER            :: npkt_pdf
INTEGER            :: all_status
REAL(KIND=PREC_DP) :: q
REAL(KIND=PREC_DP) :: u2aver                           ! Average squared displacement
REAL(KIND=PREC_DP) :: u2aver_scale = 2.00D0   ! Scale to multiply <u^2> if conversion
!                     ! with corrlin_corrquad is needed. This increases the calculated
!                     ! intensity actually by more than the damping by <u^2>, in order
!                     ! to sharpen the distances sufficiently for later broadening
REAL(KIND=PREC_DP) :: rmin, rmax, rstep
REAL(KIND=PREC_SP) :: sigma   ! sigma for PDF Corrlin correction

REAL(KIND=PREC_SP), DIMENSION(:), ALLOCATABLE :: fq    ! Temporary array, might do inplace replacement later
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: f2aver! Temporary array, might do inplace replacement later
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: faver2! Temporary array, might do inplace replacement later
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: fu    ! Temporary array, might do inplace replacement later
REAL(KIND=PREC_SP), DIMENSION(:), ALLOCATABLE :: xq    ! Temporary array, might do inplace replacement later
REAL(KIND=PREC_SP), DIMENSION(:), ALLOCATABLE :: xfour ! Temporary array, might do inplace replacement later
REAL(KIND=PREC_SP), DIMENSION(:), ALLOCATABLE :: yfour ! Temporary array, might do inplace replacement later
!
IF(deb_conv .OR. .NOT.ldbw) THEN              ! DEBYE was done with convolution of ADP
   RETURN
ENDIF
!
npkt_fft = 2**18
!
ALLOCATE(fq    (0:POW_MAXPKT),stat = all_status)
ALLOCATE(faver2(0:POW_MAXPKT),stat = all_status)
ALLOCATE(f2aver(0:POW_MAXPKT),stat = all_status)
ALLOCATE(fu    (0:POW_MAXPKT),stat = all_status)
ALLOCATE(xq    (0:POW_MAXPKT),stat = all_status)
!
DO k=0, npkt
   q = k*pow_deltaq + pow_qmin
   xq(k) = q
ENDDO
!
! Prepare average displacement factor
!
!write(*,*) ' NPKT ', npkt, pha_form(1,1,1), pha_form(npkt,1,1)
u2aver = 0
DO j=1,pha_nscat(pha_curr)                      ! Sum over the entries for current phase
   DO k=0, npkt
      q = k*pow_deltaq + pow_qmin
      faver2(k) = faver2(k) + &
         pha_form(k,j,pha_curr)   *          &
         REAL(pha_niscat(j,pha_curr)*pha_occ(j,pha_curr)/pha_nreal(pha_curr),KIND=PREC_DP)
      f2aver(k) = f2aver(k) + &
         pha_form(k,j,pha_curr)**2*          &
         REAL(pha_niscat(j,pha_curr)*pha_occ(j,pha_curr)/pha_nreal(pha_curr),KIND=PREC_DP)
      fu(k) = fu(k) +                                               &
         pha_form(k,j,pha_curr)**2*                                                &
         REAL(pha_niscat(j,pha_curr)/pha_nreal(pha_curr))*                                &
         REAL(pha_occ(j,pha_curr)*EXP(-q**2/8./PI**2*cr_dw(j)),KIND=PREC_DP)
   ENDDO
   u2aver = u2aver + (cr_dw(j)) * pha_niscat(j,pha_curr)/pha_nreal(pha_curr)*pha_occ(j,pha_curr)
ENDDO
faver2 = faver2**2
u2aver = u2aver /8./REAL(PI)**2
!
! Prepare S(Q), initially copy powder pattern and divide by Debye-Waller term
!
!     k = npkt
!     q = (k-1)*pow_deltaq + pow_qmin
!
!open(77,file='POWDER/prae_corrlin.INTE',status='unknown')
!DO k =1,npkt
!      q = (k-1)*pow_deltaq + pow_qmin
!write(77,'(2(2x,G17.7E3))') q, pha_powder(k,pha_curr)
!enddo
!close(77)
!
!write(*,*) ' EXPONENT PH ? ', u2aver*u2aver_scale, 0.5*(u2aver*u2aver_scale)*q**2, EXP(-0.5*(u2aver*u2aver_scale)*q**2)
!read(*,*) k
IF(pha_calc(pha_curr) == POW_COMPL) THEN           ! Complete calculation mode
   DO k=0, npkt
      q = k*pow_deltaq + pow_qmin
      fq(k) = pha_powder(k, pha_curr)/q**2 /EXP(-0.5*(u2aver*u2aver_scale)*q**2)
   ENDDO
ELSEIF(pha_calc(pha_curr) == POW_DEBYE) THEN       ! DEBYE calculation mode
   DO k=0, npkt
      q = k*pow_deltaq + pow_qmin
      fq(k) = pha_powder(k, pha_curr)      /EXP(-0.5*(u2aver*u2aver_scale)*q**2)
   ENDDO
ENDIF
!
!open(87,file='POWDER/prae_corrlin.INTERMEDIATE', status='unknown')
!DO k=1, npkt
!      q = (k-1)*pow_deltaq + pow_qmin
!write(87,'(5F16.6)') q, fq(k)
!enddo
!close(87)
IF(deb_conv .OR. .NOT.ldbw) THEN              ! DEBYE was done with convolution of ADP
   DO k=0, npkt
      q = k*pow_deltaq + pow_qmin
      fq(k) = (fq(k) / REAL(faver2(k))                              &
                  + 0.0 - f2aver(k)/faver2(k) ) * q
   ENDDO
ELSE
!write(*,*) ' DO CORRECTION FOR CLIN '
   DO k=0, npkt
      q = k*pow_deltaq + pow_qmin
      fq(k) = ((fq(k) - f2aver(k))/faver2(k) ) * q
   ENDDO
ENDIF
!open(87,file='phases.clin_fq1', status='unknown')
!DO k=1, npkt
!      q = (k-1)*pow_deltaq + pow_qmin
!write(87,'(5F16.6)') q, faver2(k), f2aver(k), fu(k) , fq(k)
!enddo
!close(87)
!
! Set limits for PDF calculation
!
rmin = 0.01D0
rmax = 400.0D0
rstep = 0.01D0
npkt_pdf = INT((rmax-rmin)/rstep) + 1
!IF(out_user_limits) THEN
!   rmin     = 0.50 ! out_user_values(1)
!   rmax     = out_user_values(2)*1.25
!   rstep    = out_user_values(3)
!   npkt_pdf = NINT((out_user_values(2)*1.25-0.5000000000000000)/out_user_values(3)) + 1
!!        npkt_pdf = NINT((out_user_values(2)*1.25-out_user_values(1))/out_user_values(3)) + 1
!ELSE
!   rmin     = pdf_rminu
!   rmax     = pdf_rmaxu*1.25
!   rstep    = pdf_deltaru
!   npkt_pdf = NINT((rmax-rmin)/pdf_deltaru) + 1
!ENDIF
rstep = REAL((rmax-rmin)/(npkt_pdf-1), KIND=PREC_DP)
!
ALLOCATE(xfour(1:npkt_pdf))
ALLOCATE(yfour(1:npkt_pdf))
!
zero_last1: DO ii = npkt,2,-1
   IF(fq(ii)*fq(ii-1)<0) THEN   ! Different signs , found zero point
      fq(ii) = 0.0
      EXIT zero_last1
   ELSEIF(fq(ii)==0.0) THEN       ! Found exact zero point
      EXIT zero_last1
   ELSE
      fq(ii) = 0.0
   ENDIF
ENDDO zero_last1
!
!open(87,file='phases.clin_fq', status='unknown')
!DO k=1, npkt !_pdf
!write(87,'(5F16.6)') xq(k), fq(k) 
!enddo
!close(87)
!write(*,*) ' START FFT'
CALL fft_fq(npkt, xq, fq, pow_qmin, pow_qmax, pow_deltaq, rmin, rmax, rstep, &
            npkt_fft, npkt_pdf, xfour, yfour)
!open(87,file='phases.clin_grcalc', status='unknown')
!DO k=1, npkt_pdf
!write(87,'(5F16.6)') xfour(k), yfour(k) 
!enddo
!close(87)
!write(*,*) ' FINISHED FFT'
sigma = 2.0*(u2aver*u2aver_scale)              ! TO BE REPLACED BY ATOMIC B VALUE
!CALL powder_conv_corrlin(yfour, REAL(rmin),REAL(rmax), REAL(rstep),   &
!                         sigma, pdf_clin_a, pdf_cquad_a, pdf_rcut, pow_width,   &
!                         npkt_pdf)
!open(77,file='POWDER/post_corrlin.PDF',status='unknown')
!DO ii=1,npkt_pdf
!write(77,'(2(2x,G17.7E3))') xfour(ii), yfour(ii)
!enddo
!close(77)
!
! NOw start the process back to the original powder intensity
!
!write(*,*) ' GOING BACKWARDS '
CALL fft_fq(npkt_pdf, xfour, yfour, rmin, rmax, rstep, &
            pow_qmin, pow_qmax, pow_deltaq, npkt_fft, npkt, xq, fq)
!open(77,file='POWDER/post_corrlin.FQ',status='unknown')
!DO ii=1,npkt
!write(77,'(2(2x,G17.7E3))') xq(ii), fq(ii)
!enddo
!write(*,*) ' FINISHED FFT'
!close(77)
DO k=1, npkt
   fq(k) = fq(k)/xq(k)
ENDDO
If(xq(0)> 0.0) THEN
  fq(0) = fq(0)/xq(0)
ELSE
  fq(0) = fq(1)
ENDIF
fq = fq + 1.0
!
!open(77,file='POWDER/post_corrlin.SQ',status='unknown')
!DO ii=1,npkt
!write(77,'(2(2x,G17.7E3))') xq(ii), fq(ii)
!enddo
!close(77)
!write(*,*) ' FINISHED POST_SQ'
DO k=0, npkt
   q = k*pow_deltaq + pow_qmin
!  fq(k) = ((fq(k) - f2aver(k))/faver2(k) ) * q
   fq(k) = (fq(k)  - &
                    (+ 1.0 - exp(-q**2*u2aver)))*faver2(k)
ENDDO
!open(87,file='POWDER/post_corrlin.INTERMEDIATE', status='unknown')
!DO k=1, npkt
!      q = (k-1)*pow_deltaq + pow_qmin
!write(87,'(5F16.6)') q, fq(k)
!enddo
!close(87)
IF(pha_calc(pha_curr) == POW_COMPL) THEN           ! Complete calculation mode
   DO k=0, npkt
      q = k*pow_deltaq + pow_qmin
!     fq(k) = pha_powder(k, pha_curr)/q**2 /EXP(-0.5*(u2aver*u2aver_scale)*q**2)
      pha_powder(k, pha_curr) = fq(k) * q**2 * EXP(-0.5*(u2aver*u2aver_scale)*q**2)
   ENDDO
ELSEIF(pha_calc(pha_curr) == POW_DEBYE) THEN       ! DEBYE calculation mode
   DO k=0, npkt
      q = k*pow_deltaq + pow_qmin
!     fq(k) = pha_powder(k, pha_curr)      /EXP(-0.5*(u2aver*u2aver_scale)*q**2)
      pha_powder(k, pha_curr) = fq(k)        * EXP(-0.5*(u2aver*u2aver_scale)*q**2)
   ENDDO
ENDIF
!open(77,file='POWDER/post_corrlin.INTE',status='unknown')
!DO ii=1,npkt
!write(77,'(2(2x,G17.7E3))') xq(ii), pha_powder(ii,pha_curr)
!enddo
!close(77)
!
DEALLOCATE(fq)
DEALLOCATE(faver2)
DEALLOCATE(f2aver)
DEALLOCATE(fu)
DEALLOCATE(xq)
DEALLOCATE(xfour)
DEALLOCATE(yfour)
!
END SUBROUTINE phases_corr
!
!*******************************************************************************
!
REAL function lorentz (ttheta, flag_fq) 
!+                                                                      
!-                                                                      
USE discus_config_mod 
USE powder_mod 
USE trig_degree_mod
!                                                                       
IMPLICIT none 
!                                                                       
REAL   , INTENT(IN) :: ttheta 
INTEGER, INTENT(IN) :: flag_fq
!                                                                       
!                                                                       
lorentz = 1.0
      
IF (pow_four_type.eq.POW_DEBYE) THEN 
   lorentz = 1.0 
ELSE 
   IF(flag_fq==0) THEN
      IF (pow_lp.eq.POW_LP_BRAGG) THEN 
         lorentz = 0.5 / sind (0.5 * ttheta) / sind (ttheta) 
      ELSEIF (pow_lp.eq.POW_LP_NEUT) THEN 
         lorentz = 0.5 / sind (0.5 * ttheta) / sind (ttheta) 
      ELSEIF (pow_lp.eq.POW_LP_NONE) THEN 
         lorentz = 1.0 
      ELSEIF (pow_lp.eq.POW_LP_SYNC) THEN 
         lorentz = 0.5 / sind (0.5 * ttheta) / sind (ttheta) 
      ENDIF 
   ELSEIF(flag_fq==1) THEN 
      lorentz = 0.5 / sind (0.5 * ttheta) / sind (ttheta) 
   ENDIF 
ENDIF 
!                                                                       
END FUNCTION lorentz                          
!
!*****7*****************************************************************
!
REAL FUNCTION polarisation (ttheta) 
!+                                                                      
!-                                                                      
USE discus_config_mod 
USE powder_mod 
USE trig_degree_mod
!                                                                       
IMPLICIT none 
!                                                                       
REAL ttheta 
!
!
!                                                                       
polarisation = 1.0
      
IF (pow_lp.eq.POW_LP_BRAGG) THEN 
   polarisation = (1. + (cosd(ttheta))**2 * pow_lp_fac)/(1. + pow_lp_fac)
ELSEIF (pow_lp.eq.POW_LP_NEUT) THEN 
   polarisation = 1.0 
ELSEIF (pow_lp.eq.POW_LP_NONE) THEN 
   polarisation = 1.0 
ELSEIF (pow_lp.eq.POW_LP_SYNC) THEN 
   polarisation = pow_lp_fac + (1. - pow_lp_fac) * (cosd(ttheta))**2 * pow_lp_cos
ENDIF 
!                                                                       
END FUNCTION polarisation                     
!
!*****7*****************************************************************
!
REAL FUNCTION lorentz_pol (ttheta) 
!+                                                                      
!-                                                                      
USE discus_config_mod 
USE powder_mod 
USE trig_degree_mod
!                                                                       
IMPLICIT none 
!                                                                       
REAL ttheta 
!                                                                       
!     REAL sind, cosd 
!                                                                       
IF (pow_four_type.eq.POW_DEBYE) THEN 
   lorentz_pol = 1.0 
ELSE
   lorentz_pol = (1-pow_lp_fac+pow_lp_fac*(cosd(pow_lp_ang))**2*(cosd(ttheta))**2)/ &
                 (2.*(sind(0.5*ttheta))**2*cosd(0.5*ttheta))
ENDIF
!
END FUNCTION lorentz_pol
!
!*****7*****************************************************************
!
SUBROUTINE phases_reset
!
USE discus_allocate_appl_mod
USE phases_mod
!
IMPLICIT NONE
!
CALL alloc_phases(1, 1, 1)
!
pha_multi  = .FALSE.
pha_n      = 1
pha_curr   = 1
pha_nscat  = 1    ! No of atom types for each phase (      i)
pha_calc   = 0    ! Calc mode Comp/Debye each phase (      i)
pha_frac   = 0.0  ! weight fraction User intent     (      i)
pha_weight = 0.0  ! Weight          for each phase   (      i)
pha_scale  = 0.0  ! Scale temp>frac for each phase   (      i)
pha_nreal  = 0.0  ! Number real atoms at     phase   (      i)
pha_ncreal = 0.0  ! Number real atoms /unit call at phase   (      i)
pha_powder = 0.0  ! Powder pattern for each phase   (q,    i)
pha_form   = 0.0D0  ! Form factors   for each phase   (q,is,i)
pha_adp    = 0.0  ! B-values ADP   for each phase   (is,   i)
pha_occ    = 0.0  ! Occupancies    for each phase   (is ,  i)
pha_niscat = 0
!
END SUBROUTINE phases_reset
!
!*****7*****************************************************************
!
END MODULE phases_set_mod
