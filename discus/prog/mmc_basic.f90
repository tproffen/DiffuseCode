module mmc_basic_mod
!
!  Basic procedures for mmc, that are common to atom and molecule
!
!  mmc_test_multi                      ! Test if energy has been lowered
!  mmc_correlations                    ! Accumulate all correlations
!  mmc_correlations_occ                ! Occupational correlation
!  mmc_correlations_uni                ! Unidirectional occupational corrrelation
!  mmc_correlations_group              ! Group-wsie correlations
!  index2angles                        ! helper routine to convert indices fro angular correlations
!
contains
!
!*****7*****************************************************************
!
subroutine mmc_initial(old_chem_period, itry, igen, iacc_good, iacc_neut, iacc_bad, &
           done, loop)
!-
!  Perform the initial mmc run setup
!+
!
use crystal_mod
use chem_mod
use discus_allocate_appl_mod
use mmc_mod
use molecule_mod
!
implicit none
!
LOGICAL, DIMENSION(3), INTENT(OUT)  :: old_chem_period    ! Store current periodic boundary settings
INTEGER(KIND=PREC_INT_LARGE), INTENT(OUT) :: itry         ! Number tired moves
INTEGER(KIND=PREC_INT_LARGE), INTENT(OUT) :: igen         ! number generated moves
INTEGER                     , INTENT(out) :: iacc_good    ! Number accepted good
INTEGER                     , INTENT(out) :: iacc_neut    ! Number accepted neutral
INTEGER                     , INTENT(out) :: iacc_bad     ! Number accepted bad
!NTEGER       ,dimension(3) , intent(out) :: lbeg
logical                     , intent(out) :: done
logical                     , intent(out) :: loop
!
integer :: i                 ! Dummy index
INTEGER :: n_corr     ! Dummy for angle allocation
INTEGER :: n_scat     ! Dummy for angle allocation
INTEGER :: n_site     ! Dummy for angle allocation
INTEGER :: n_mole     ! Dummy for angle allocation
INTEGER :: n_angles   ! Dummy for angle allocation
!
old_chem_period = chem_period
DO i=1,3
  IF(cr_icc(i)==1) chem_period(i) = .FALSE.
ENDDO
!
! Basic allocation
!
n_corr = MAX(CHEM_MAX_COR,MMC_MAX_CORR)
n_scat = MAX(MAXSCAT, MMC_MAX_SCAT, 3)    ! 3 is needed for 'group'
n_site = MAX(MAXSCAT, MMC_MAX_SITE)
n_mole = MOLE_MAX_TYPE
call alloc_mmc      (CHEM_MAX_COR, MC_N_ENERGY, MAXSCAT, n_site)
call alloc_mmc_pid  (CHEM_MAX_COR, MC_N_ENERGY, MAXSCAT, n_site)
IF (mmc_n_angles == 0          .OR.     &
    mmc_n_angles >= MMC_MAX_ANGLES) THEN 
   n_angles = max(mmc_n_angles+20, int(MMC_MAX_ANGLES*1.025))
   call alloc_mmc_angle (CHEM_MAX_COR,n_angles)
ENDIF
IF(CHEM_MAX_COR>UBOUND(mmc_buck_a,1) .OR. MAXSCAT>UBOUND(mmc_buck_a,2)) THEN
   call alloc_mmc_buck( CHEM_MAX_COR, MAXSCAT )
ENDIF
IF(CHEM_MAX_COR>UBOUND(mmc_len_a ,1) .OR. MAXSCAT>UBOUND(mmc_len_a ,2)) THEN
   call alloc_mmc_lenn( CHEM_MAX_COR, MAXSCAT )
ENDIF
IF(CHEM_MAX_COR>UBOUND(mmc_rep_a ,1) .OR. MAXSCAT>UBOUND(mmc_rep_a ,2)) THEN
   CALL alloc_mmc_rep ( CHEM_MAX_COR, MAXSCAT )
ENDIF
!
mmc_h_nfeed = 0                                   ! Clear overall feedback run counter
IF(ALLOCATED(mmc_h_diff  )) DEALLOCATE(mmc_h_diff)
IF(ALLOCATED(mmc_h_targ  )) DEALLOCATE(mmc_h_targ)
IF(ALLOCATED(mmc_h_aver  )) DEALLOCATE(mmc_h_aver)
IF(ALLOCATED(mmc_h_aver_r)) DEALLOCATE(mmc_h_aver_r)
IF(ALLOCATED(mmc_h_maxd  )) DEALLOCATE(mmc_h_maxd)
ALLOCATE(mmc_h_diff  (mmc_h_number+10,0:MMC_H_NNNN-1))
ALLOCATE(mmc_h_targ  (mmc_h_number+10))
ALLOCATE(mmc_h_aver  (mmc_h_number+10))
ALLOCATE(mmc_h_aver_r(mmc_h_number+10))
ALLOCATE(mmc_h_maxd  (mmc_h_number+10, 0:MMC_H_NNNN-1))
mmc_h_diff   = 0.0
mmc_h_targ   = 0.0
mmc_h_aver   = 0.0
mmc_h_aver_r = 0.0
mmc_h_maxd   = 0.0
!
!
n_e_av_p = 0                                 ! Initialize counters for energy changes
n_e_av_m = 0    !
n_e_av_z = 0    !
e_aver_p = 0.0  !
e_aver_m = 0.0  !
!                                                                       
!------ reset some counters                                             
!                                                                       
igen      = 0 
itry      = 0 
iacc_good = 0 
iacc_neut = 0 
iacc_bad = 0 
loop = .TRUE. 
done = .TRUE. 
!lbeg(1) = -1 
!lbeg(2) =  0 
!lbeg(3) =  0 
mmc_pid_diff = 0.0
mmc_pid_inte = 0.0
mmc_pid_deri = 0.0
mmc_pre_corr = 0.0
mmc_pid_pid  = 0.0
mmc_pid_pid(1,1) = 0.30
mmc_pid_pid(2,1) = 0.05   !0.06
mmc_pid_pid(3,1) =-0.10
mmc_pid_pid_n = 0
mmc_pid_change = 0.0
!
end subroutine mmc_initial
!
!*****7*****************************************************************
!
SUBROUTINE mmc_test_multi (iacc_good, iacc_neut, iacc_bad, &
           e_new, e_old, laccept)                                                          
!+                                                                      
!     Tests performed MC move                                           
!-                                                                      
USE mc_mod 
USE mmc_mod 
USE lib_random_func
USE random_mod
!                                                                       
IMPLICIT none 
!SAVE
!                                                                       
INTEGER, INTENT(INOUT) :: iacc_good, iacc_neut, iacc_bad 
REAL(kind=PREC_DP),    INTENT(IN)    :: e_old (0:MC_N_ENERGY) 
REAL(kind=PREC_DP),    INTENT(IN)    :: e_new (0:MC_N_ENERGY) 
LOGICAL, INTENT(OUT)   :: laccept 
!                                                                       
INTEGER :: i 
REAL(kind=PREC_DP)    :: e_del 
real(kind=PREC_DP)    :: e_neu
REAL(kind=PREC_DP)    :: e_ran 
REAL(kind=PREC_DP)    :: e_delta 
REAL(kind=PREC_DP)    :: r1
!                                                                       
!                                                                       
e_del = 0.0
DO i = 1, MC_N_ENERGY 
   IF (mmc_cor_energy(0, i) ) THEN 
      e_delta = NINT((e_new (i) - e_old (i) ) *1.0E4)*1.0E-4
      e_del = e_del + e_delta 
   ENDIF 
ENDDO 
!write(*,*) 'OLD ', -(sqrt(abs(e_old(7))))
!write(*,*) 'new ', -(sqrt(abs(e_new(7))))
!write(*,*) 'del ',(e_new(i)-e_old(i),i=7,7) , mo_kt !1,MC_N_ENERGY)
if(mo_kt>0) then
IF (e_del <  0) THEN 
   laccept = .TRUE. 
ELSE 
   IF (mo_kt <  1.0e-10) THEN 
      laccept = .FALSE. 
   ELSE 
      e_ran = exp ( - e_del / mo_kt) 
      e_ran = e_ran / (1 + e_ran) 
      CALL RANDOM_NUMBER(r1)
      laccept = (e_ran > r1          ) 
   ENDIF 
ENDIF
else
   laccept = .FALSE.
   if(e_new(7)>0.0) then                    ! Positive energy, 
      e_neu = sqrt(e_new(7))
   else
      e_neu = - sqrt(abs(e_new(7)))
   endif
!     if(e_del < 0) then
!        laccept = .TRUE.
!     endif
!  else
!     e_ran = exp ( - e_del / mo_kt) 
      e_ran = exp(-(e_neu + mmc_depth(1,7,1,1)             )/abs(mo_kt))
!      e_ran = exp(-(e_neu + mmc_depth(1,7,1,1) - abs(mo_kt))/abs(mo_kt))
!!!      e_ran =      exp ( - (mmc_depth(1,7,1,1)-sqrt(abs(e_new(7))) -abs(mo_kt)) / abs(mo_kt)) 
!write(*,*) e_ran , exp ( - (mmc_depth(1,7,1,1)-sqrt(abs(e_new(7))) -abs(mo_kt)) / abs(mo_kt)) , &
!(mmc_depth(1,7,1,1)-sqrt(abs(e_new(7))) - abs(mo_kt) ), &
!(mmc_depth(1,7,1,1)-sqrt(abs(e_new(7))) - abs(mo_kt) )/ abs(mo_kt)
      e_ran = e_ran / (1 + e_ran) 
      CALL RANDOM_NUMBER(r1)
      laccept = (e_ran > r1          ) 
!write(*,*) e_ran, r1, laccept
!  endif
endif
!                                                                       
IF (laccept) THEN 
   IF (e_del <  0.0) THEN 
      iacc_good = iacc_good+1 
!write(*,*) ' ACCEPTED ?', laccept, ' Good'
   ELSEIF(e_del==0) THEN
      iacc_neut = iacc_neut + 1
!write(*,*) ' ACCEPTED ?', laccept, ' Neutral'
   ELSE 
      iacc_bad = iacc_bad+1 
!write(*,*) ' ACCEPTED ?', laccept, ' Bad '
   ENDIF 
   DO i = 1, MC_N_ENERGY 
      IF (mmc_cor_energy(0, i) ) THEN 
         e_delta = NINT((e_new (i) - e_old (i) ) *1.0E4)*1.0E-4
         IF (e_delta > 0.0) THEN 
            e_aver_p (i) = e_aver_p (i) + e_delta 
            n_e_av_p (i) = n_e_av_p (i) + 1 
         ELSEIF (e_delta < 0.0) THEN 
            e_aver_m (i) = e_aver_m (i) + e_delta 
            n_e_av_m (i) = n_e_av_m (i) + 1 
         ELSE 
            n_e_av_z (i) = n_e_av_z (i) + 1 
         ENDIF 
      ENDIF 
   ENDDO
ENDIF 
!                                                                       
END SUBROUTINE mmc_test_multi                 
!
!*****7*************************************************************************
!
SUBROUTINE mmc_correlations (lout, rel_cycl, done, lfinished, lfeed, maxdev) 
!-                                                                      
!     Determines the achieved correlations                              
!                                                                       
!+                                                                      
USE crystal_mod 
USE chem_mod 
!USE chem_menu
USE chem_aver_mod
USE chem_neig_multi_mod
USE atom_env_mod
USE celltoindex_mod
USE metric_mod
USE mc_mod 
USE mmc_mod 
!
USE debug_mod 
USE errlist_mod 
use precision_mod
use param_mod
USE prompt_mod 
!
IMPLICIT none 
!                                                                       
LOGICAL , INTENT(IN) :: lout        ! Flag for output yes/no
REAL(kind=PREC_DP)    , INTENT(IN) :: rel_cycl ! Relative progress along cycles
LOGICAL , INTENT(INOUT) :: done     ! MMC is converged/ stagnates
LOGICAL , INTENT(IN)    :: lfinished ! MMC is finished
LOGICAL , INTENT(IN)    :: lfeed     ! Perform feedback algorithm
real(kind=PREC_DP), dimension(2), intent(inout) :: maxdev
! 
!                                                                       
CHARACTER(LEN=30) :: energy_name (0:MC_N_ENERGY) 
!                                                                       
INTEGER :: ic, je, ic_a 
INTEGER :: is, js, ls 
INTEGER :: iis, jjs, lls, iic, kk 
INTEGER :: i, j, k, l 
INTEGER :: icent 
!
LOGICAL :: searching 
LOGICAL   :: lfirst = .TRUE.  ! Flag to write output only at first instance
!                                                                       
INTEGER :: ncent 
REAL(KIND=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE :: patom ! (3, 0:MAX_ATOM_ENV, MMC_MAX_CENT) 
INTEGER           , DIMENSION(  :,:), ALLOCATABLE :: iatom ! (0:MAX_ATOM_ENV, MMC_MAX_CENT) 
LOGICAL           , DIMENSION(  :,:), ALLOCATABLE :: tatom ! (0:MAX_ATOM_ENV, MMC_MAX_CENT) 
INTEGER           , DIMENSION(    :), ALLOCATABLE :: natom ! ( MMC_MAX_CENT) 
!                                                                       
INTEGER           , DIMENSION(:,:), ALLOCATABLE :: bl_anz ! (0:DEF_maxscat, 0:DEF_maxscat) 
REAL(KIND=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: bl_sum ! (0:DEF_maxscat, 0:DEF_maxscat) 
REAL(KIND=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: bl_s2  ! (0:DEF_maxscat, 0:DEF_maxscat) 
REAL(kind=PREC_DP) :: u (3), v (3), d (3) 
REAL(kind=PREC_DP) :: dist 
real(kind=PREC_DP) :: change    ! Changes in (target-achieved) from feedback to feedback
real(kind=PREC_DP), dimension(4), save :: conv_val ! Convergence values
!
INTEGER :: n_cn
INTEGER, DIMENSION(:), ALLOCATABLE :: ncentral
INTEGER, DIMENSION(:,:), ALLOCATABLE :: p_cn
!                                                                       
INTEGER nneigh 
REAL(kind=PREC_DP) :: damp = 1.0
!                                                                       
REAL(kind=PREC_DP) :: wi, wis 
REAL(kind=PREC_DP) :: divider 
!                                                                       
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: ba_sum ! (CHEM_MAX_COR * MMC_MAX_ANGLES) 
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: ba_s2  ! (CHEM_MAX_COR * MMC_MAX_ANGLES) 
INTEGER           , DIMENSION(:), ALLOCATABLE :: ba_anz ! (CHEM_MAX_COR * MMC_MAX_ANGLES) 
!                                                                       
INTEGER :: icc (3), jcc (3) 
REAL(kind=PREC_DP) :: idir (3), jdir (3), disi (3), disj (3) 
REAL(kind=PREC_DP) :: rdi=1.0, rdj=1.0, dpi=1.0, dpj
INTEGER           , DIMENSION(:,:), ALLOCATABLE :: xnn !(0:maxscat, 0:maxscat) 
REAL(KIND=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: xij !(0:maxscat, 0:maxscat) 
REAL(KIND=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: xi2 !(0:maxscat, 0:maxscat) 
REAL(KIND=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: xj2 !(0:maxscat, 0:maxscat) 
!
!                                                                       
DATA energy_name / 'none', 'Chemical correlation    ', 'Displacement correlation', &
     'Distance correlation (Hooke)', 'Angular      correlation', &
     'Vector       correlation',     'Distance     correlation', &
     'Lennard Jones potential ',     'Buckingham    potential ', & 
     'Repulsive     potential ',     'Coordination number     ', &
     'Unidirectional corr     ',     'Groupwise    correlation' /         
!
!lout = lout_in .and. mmc_out_feed    ! Combine system and user settings
!
mmc_pid_pid(:,2) = 0.0    ! Reset all average (diff, inte, deriv) terms
mmc_pid_pid_n    = 0      ! Reset counter of contributing correlations
damp = 0.01 + 0.99*exp(-4.0*rel_cycl)
damp = 0.01 + 0.99*exp(-4.0*rel_cycl**2)
!damp = 1.0
!
ALLOCATE(patom(3, 0:MAX_ATOM_ENV, MMC_MAX_CENT))
ALLOCATE(iatom(0:MAX_ATOM_ENV, MMC_MAX_CENT))
ALLOCATE(tatom(0:MAX_ATOM_ENV, MMC_MAX_CENT))
ALLOCATE(natom(MMC_MAX_CENT))
ALLOCATE( xnn (0:MAXSCAT, 0:MAXSCAT) )
ALLOCATE( xij (0:MAXSCAT, 0:MAXSCAT) )
ALLOCATE( xi2 (0:MAXSCAT, 0:MAXSCAT) )
ALLOCATE( xj2 (0:MAXSCAT, 0:MAXSCAT) )
ALLOCATE( ba_sum (CHEM_MAX_COR * MMC_MAX_ANGLES) )
ALLOCATE( ba_s2  (CHEM_MAX_COR * MMC_MAX_ANGLES) )
ALLOCATE( ba_anz (CHEM_MAX_COR * MMC_MAX_ANGLES) )
ALLOCATE( bl_anz (0:MAXSCAT, 0:MAXSCAT) )
ALLOCATE( bl_sum (0:MAXSCAT, 0:MAXSCAT) )
ALLOCATE( bl_s2 (0:MAXSCAT, 0:MAXSCAT) )
ALLOCATE( p_cn(0:MAXSCAT, 0:MAXSCAT))
!ALLOCATE( pneig(0:MAXSCAT, 0:MAXSCAT, 1:CHEM_MAX_COR) )
!
IF(.NOT.ALLOCATED(mmc_h_diff) ) THEN
   ALLOCATE(mmc_h_diff  (1:mmc_h_number + 10,0:MMC_H_NNNN-1))
   ALLOCATE(mmc_h_targ  (mmc_h_number+10))
   ALLOCATE(mmc_h_aver  (mmc_h_number+10))
   ALLOCATE(mmc_h_aver_r(mmc_h_number+10))
   ALLOCATE(mmc_h_maxd  (mmc_h_number+10,0:MMC_H_NNNN-1))
   mmc_h_diff  = 0.0
   mmc_h_targ  = 0.0
   mmc_h_aver  = 0.0
   mmc_h_aver_r= 0.0
   mmc_h_maxd  = 0.0
   mmc_h_index =  -1
   mmc_h_ncycl =  0
ENDIF
mmc_h_ctarg = 0          ! Start with no targets in history
mmc_h_index = MOD(mmc_h_index+1,MMC_H_NNNN)   ! increment current index
mmc_m_index = max(mmc_m_index, mmc_h_index)   ! Record maximum number of feedback cycles so far
mmc_h_nfeed = mmc_h_nfeed + 1                 ! Increment overall feedback counter
!                                                                       
!------ Write title line                                                
!                                                                       
IF (lout) THEN 
   WRITE (output_io, 410) 
ENDIF 
!                                                                       
!     Get the average structure for the distance energies               
!                                                                       
!     IF (mmc_cor_energy (0, MC_DISP)    .OR.mmc_cor_energy (0, MC_SPRING) &
!     .OR.mmc_cor_energy (0, MC_LENNARD) .OR.mmc_cor_energy (0, MC_BUCKING)&
!     .OR.mmc_cor_energy (0,MC_REPULSIVE) ) THEN                                                
IF(     mmc_cor_energy (0, MC_DISP)                                  &
   .OR. mmc_move_prob(MC_MOVE_SWDISP) > 0                            &
   .OR. mmc_move_prob(MC_MOVE_INVDISP) > 0                           &
   ) THEN
   CALL chem_aver (.FALSE., .TRUE.) 
ENDIF 
!                                                                       
!     Reset all achieved correlations                                   
!                                                                       
!
      DO ic = 1, CHEM_MAX_COR 
         DO je = 1, MC_N_ENERGY 
            DO is = - 1, MAXSCAT 
               DO js = - 1, MAXSCAT 
                  mmc_ach_corr (ic, je, is, js) = 0.0 
                  mmc_ach_sigm (ic, je, is, js) = 0.0 
               ENDDO 
            ENDDO 
         ENDDO 
      ENDDO 
      DO i = 1, CHEM_MAX_COR * MMC_MAX_ANGLES 
         mmc_ach_angl (i) = 0.0 
         mmc_ang_sigm (i) = 0.0 
      ENDDO 
ALLOCATE(ncentral(0:MAXSCAT))
!                                                                       
!     Loop over all correlations                                        
!                                                                       
main_corr: DO ic = 1, chem_ncor 
   ncentral(:) = 0
!  IF (lout) THEN 
!     WRITE (output_io, * ) 
!  ENDIF 
   DO is = 0, MAXSCAT 
      DO js = 0, MAXSCAT 
         mmc_pneig (is, js,ic) = 0 
         p_cn  (is, js) = 0 
         bl_sum (is, js) = 0.0 
         bl_s2 (is, js) = 0.0 
         bl_anz (is, js) = 0 
         xnn (is, js) = 0 
         xij (is, js) = 0 
         xi2 (is, js) = 0 
         xj2 (is, js) = 0 
      ENDDO 
   ENDDO 
   DO i = 1, CHEM_MAX_COR * MMC_MAX_ANGLES 
      ba_sum (i) = 0.0 
      ba_s2 (i) = 0.0 
      ba_anz (i) = 0 
   ENDDO 
   IF (mmc_cor_energy (0, MC_DISP) ) THEN 
      DO i = 1, 3 
         idir (i) = chem_dir (i, 1, ic) 
         jdir (i) = chem_dir (i, 2, ic) 
      ENDDO 
!                                                                       
!------ calculate correlations                                          
!                                                                       
      rdi = skalpro (idir, idir, cr_gten) 
      rdj = skalpro (jdir, jdir, cr_gten) 
      IF (rdi > 0.0) rdi = sqrt (rdi) 
      IF (rdj > 0.0) rdj = sqrt (rdj) 
   ENDIF 
!                                                                       
!     -- Loop over all atoms                                            
!                                                                       
main_atoms: DO i = 1, cr_natoms 
      is = cr_iscat (1,i) 
      ncentral(is) = ncentral(is) + 1
      CALL chem_neighbour_multi (i, ic, iatom, patom, tatom, natom, ncent, MAX_ATOM_ENV, MMC_MAX_CENT)
      IF (ier_num /= 0) THEN
         RETURN 
      ENDIF
!                                                                       
!------ ---- In case of Displacement correlation, calculate             
!            displacement of central atom                                                 
!                                                                       
is_mc_disp: IF (mmc_cor_energy (ic, MC_DISP) ) THEN 
         CALL indextocell (i, icc, is) 
         DO j = 1, 3 
            disi (j) = cr_pos (j, i) - chem_ave_pos (j, is) - &
                       REAL(icc (j) - 1) - cr_dim0 (j, 1)                                          
         ENDDO 
!                                                                       
      IF (chem_ldall (ic) ) THEN 
         DO j = 1, 3 
            jdir (j) = disi (j) 
         ENDDO 
         rdj = skalpro (jdir, jdir, cr_gten) 
         IF (rdj > 0.0) THEN 
            rdj = sqrt (rdj) 
         ELSE 
            rdj = 1.0 
         ENDIF 
         dpi = 1.0 
      ELSE 
         dpi = skalpro (disi, idir, cr_gten) / rdi 
      ENDIF 
   ENDIF is_mc_disp
!                                                                       
!     ---- Loop over all centers                                        
!                                                                       
main_cent: DO icent = 1, ncent 
!write(*,'(3(2i5,a)20i5)') i, is, ' : ', icent, ncent, ' : ', natom(icent), iatom(0,icent), ' : ', iatom(1:natom(icent), icent)
!                                                                       
!     ------- Since the loop is over all atoms, do central atoms only   
!                                                                       
is_cent:  IF (i == iatom (0, icent) ) THEN 
is_energy:  IF (mmc_cor_energy (ic, MC_OCC)        .OR. &
                mmc_cor_energy (ic, MC_UNI)        .OR. &
                mmc_cor_energy (ic, MC_GROUP)      .OR. &
                mmc_cor_energy (ic, MC_DISP)       .OR. &
                mmc_cor_energy (ic, MC_SPRING)     .OR. &
                mmc_cor_energy (ic, MC_LENNARD)    .OR. &
                mmc_cor_energy (ic, MC_REPULSIVE)  .OR. &
                mmc_cor_energy (ic, MC_COORDNUM )  .OR. &
                mmc_cor_energy (ic, MC_BUCKING)         ) THEN    
!                                                                       
!     ---------- Loop over all neighbours                               
!                                                                       
loop_neig:  DO j = 1, natom (icent) 
               IF(.NOT.tatom(j,icent)) CYCLE loop_neig        ! Current atom is neighbour to this 'neighbor'
               js = cr_iscat (1,iatom (j, icent) ) 
               DO k = 1, 3 
                  u (k) = patom (k, 0, icent) 
               ENDDO 
!                                                                       
!     --------- Accumulate values for all Energies                      
!                                                                       
               IF (mmc_cor_energy (ic, MC_OCC)  ) THEN
!                                                                       
!     ----------- Chemical correlation, add number of atom pairs        
!                                                                       
                  IF(mmc_pair(ic, MC_OCC,is,js) /=0) THEN
!write(*,'(5i5)') i, is, j, iatom(j, icent), js
                  mmc_pneig (is, js, ic) = mmc_pneig (is, js, ic) + 1 
                  ENDIF
               ENDIF 
               IF (mmc_cor_energy (ic, MC_UNI) ) THEN 
!                                                                       
!     ----------- Chemical correlation, add number of atom pairs        
!                                                                       
                  IF(mmc_pair(ic, MC_UNI,is,js) /=0) THEN
                  mmc_pneig (is, js, ic) = mmc_pneig (is, js, ic) + 1 
                  ENDIF
               ENDIF 
               IF(mmc_cor_energy(ic, MC_GROUP) ) THEN 
!                                                                       
!     ----------- Chemical correlation, add number of atom pairs        
!                                                                       
!write(*,*) is, js
!write(*,*) mmc_left (ic, MC_GROUP,:), ' | ', &
!mmc_left (ic, MC_GROUP,is),  &
!mmc_left (ic, MC_GROUP,js)

!write(*,*) mmc_right(ic, MC_GROUP,:), ' | ',  &
!mmc_right(ic, MC_GROUP,is), &
!mmc_right(ic, MC_GROUP,js)

                  IF(mmc_left (ic, MC_GROUP,is)/=0) THEN                    ! (A) => 
                     IF(    mmc_right(ic, MC_GROUP,js)/=0) THEN                ! (A) => (B) 
                        mmc_pneig (1 , 2, ic ) = mmc_pneig (1 , 2, ic ) + 1 
!write(*,*) ' (A) => (B)'
                     ELSEIF(mmc_left (ic, MC_GROUP,js)/=0) THEN                ! (A) => (A)
                        mmc_pneig (1 , 1, ic ) = mmc_pneig (1 , 1, ic ) + 1 
!write(*,*) ' (A) => (A)'
                     ELSE                                                      ! (A) => (other)
                        mmc_pneig (1 , 3, ic ) = mmc_pneig (1 , 3, ic ) + 1 
!write(*,*) ' (A) => (o)'
                     ENDIF
                  ELSEIF(mmc_right(ic, MC_GROUP,is)/=0) THEN                ! (B) =>
                     IF(    mmc_left (ic, MC_GROUP,js)/=0) THEN                ! (B) => (A)
                        mmc_pneig (2 , 1, ic ) = mmc_pneig (2 , 1, ic ) + 1 
!write(*,*) ' (B) => (A)'
                     ELSEIF(mmc_right(ic, MC_GROUP,js)/=0) THEN                ! (B) => (B) 
                        mmc_pneig (2 , 2, ic ) = mmc_pneig (2 , 2, ic ) + 1 
!write(*,*) ' (B) => (B)'
                     ELSE                                                      ! (B) => (other)
                        mmc_pneig (2 , 3, ic ) = mmc_pneig (2 , 3, ic ) + 1 
!write(*,*) ' (B) => (o)'
                     ENDIF
                  ELSE
                     mmc_pneig (3 , 3, ic ) = mmc_pneig (3 , 3, ic ) + 1       ! (other) => (other)
!write(*,*) ' (o) => (o)'
                  ENDIF
               ENDIF 
!read(*,*) k
!
!
!--- Coordination number, ic
!
               IF(mmc_cor_energy(ic, MC_COORDNUM)) THEN
                  IF(mmc_pair(ic,MC_COORDNUM, is,js)==-1) THEN
                     p_cn (is, js) = p_cn (is, js) + 1 
                  ENDIF
               ENDIF 
               IF (mmc_cor_energy (ic, MC_SPRING)    .OR.     &
                   mmc_cor_energy (ic, MC_LENNARD)   .OR.     &
                   mmc_cor_energy (ic, MC_REPULSIVE) .OR.     &
                   mmc_cor_energy (ic, MC_BUCKING)     ) THEN      
                  DO k = 1, 3 
                     v (k) = patom (k, j, icent) 
                     d (k) = v (k) - u (k) 
                  ENDDO 
!if(ic==3) THEN
!write(*,*) ' should add for atom pair is,js', i,is,js
!endif
                  dist = do_blen (.TRUE., u, v) 
                  js   = cr_iscat (1,iatom (j, icent) ) 
                  bl_sum (is, js) = bl_sum (is, js) + dist 
                  bl_s2  (is, js) = bl_s2 (is, js) + dist**2 
                  bl_anz (is, js) = bl_anz (is, js) + 1 
                  mmc_pneig (is, js, ic) = mmc_pneig (is, js, ic) + 1 
               ENDIF 
               IF (mmc_cor_energy (ic, MC_DISP) ) THEN 
                  CALL indextocell (iatom (j, icent), jcc, js) 
                  DO k = 1, 3 
                     disj (k) = cr_pos (k, iatom (j, icent) ) - chem_ave_pos (&
                     k, js) - REAL(jcc (k) - 1) - cr_dim0 (k, 1)            
                  ENDDO 
                  dpj = skalpro (disj, jdir, cr_gten) / rdj 
                  xij (is, js) = xij (is, js) + dpi * dpj 
                  xi2 (is, js) = xi2 (is, js) + dpi**2 
                  xj2 (is, js) = xj2 (is, js) + dpj**2 
                  xnn (is, js) = xnn (is, js) + 1 
               ENDIF 
            ENDDO loop_neig
!                     j ! Loop over all neighbours                      
         ENDIF is_energy
         IF (mmc_cor_energy (ic, MC_ANGLE) ) THEN 
!                                                                       
!     ---------- Angular Correlations                                   
!                                                                       
            is = cr_iscat (1,i) 
            DO k = 1, 3 
               u (k) = patom (k, 0, icent) 
            ENDDO 
!                                                                       
!     ---------- Double loop over all neighbours                        
!                                                                       
            DO j = 1, natom (icent) - 1 
               js = cr_iscat (1,iatom (j, icent) ) 
               DO k = 1, 3 
                  v (k) = patom (k, j, icent) 
               ENDDO 
               DO l = j + 1, natom (icent) 
                  ls = cr_iscat (1,iatom (l, icent) ) 
!                                                                       
!     -------------- Find proper entry in correlation table             
!                                                                       
                  k = 0 
                  ic_a = 0 
                  searching = .TRUE. 
                  DO WHILE (searching.AND.k <= mmc_n_angles) 
                     k = k + 1 
                     CALL index2angles (mmc_angles (k), iic, kk, iis, jjs, lls,  &
                                          MAXSCAT)
                     IF (iic == ic) THEN 
                        IF (iis == is.OR.iis ==  - 1) THEN 
                           IF (jjs ==  - 1.OR.jjs == min (js, ls) ) THEN 
                              IF (lls ==  - 1.OR.lls == max (js, ls) ) THEN 
                                 searching = .FALSE. 
                                 ic_a = k 
                              ENDIF 
                           ENDIF 
                        ENDIF 
                     ENDIF 
                  ENDDO 
!                         ! Find proper entry in correlation table      
                  IF (.NOT.searching) THEN 
                     DO k = 1, 3 
                        d (k) = patom (k, l, icent) 
                     ENDDO 
                     wi = do_bang (.TRUE., v, u, d) 
                     wis = mmc_target_angl (ic_a) 
                     IF (wis <= 90.) THEN 
                        IF (wi > 1.5 * wis) THEN 
                           wi = mod (wi + wis / 2., wis) + wis / 2. 
                        ENDIF 
                  ENDIF 
                  ba_sum (ic_a) = ba_sum (ic_a) + wi 
                  ba_s2 (ic_a) = ba_s2 (ic_a) + wi**2 
                  ba_anz (ic_a) = ba_anz (ic_a) + 1 
               ENDIF 
            ENDDO 
         ENDDO 
!                      ! j     ! Double loop over neighbours                   
      ENDIF 
      ENDIF is_cent    !       ! center atoms only                                   
   ENDDO main_cent     ! icent ! Loop over centers                               
ENDDO main_atoms       ! i     ! Loop over all atoms                                   
!                                                                       
!------ -- Summ up all energies, write output                           
!                                                                       
!     ----- Chemical correlation                                        
!                                                                       
IF (mmc_cor_energy (0, MC_OCC)) THEN
   CALL mmc_correlations_occ(ic, mmc_pneig, rel_cycl, damp, lout, lfeed,        &
        max(MAXSCAT,3), CHEM_MAX_COR, maxdev)
ENDIF
!                                                                       
!     ----- Unidirectional Chemical correlation                                        
!                                                                       
IF (mmc_cor_energy (0, MC_UNI)) THEN
   CALL mmc_correlations_uni(ic, mmc_pneig, rel_cycl, damp, lout, lfeed,        &
        max(MAXSCAT,3), CHEM_MAX_COR, maxdev)
ENDIF
!                                                                       
!     ----- Group wise correlations
!                                                                       
IF (mmc_cor_energy (0, MC_GROUP)) THEN
CALL mmc_correlations_group(ic, mmc_pneig, rel_cycl, damp, lout, lfeed,         &
     max(MAXSCAT,3), CHEM_MAX_COR, maxdev)
ENDIF
!
!  Coordination number
!
   je = MC_COORDNUM
   n_cn = 0
   cn_pair: DO is = 0, cr_nscat 
      DO js = 0 , cr_nscat 
         IF(mmc_pair(ic, MC_COORDNUM, is, js) /=  0 ) THEN 
            n_cn = n_cn + p_cn(is, js)
         ENDIF
      ENDDO
   ENDDO cn_pair
!write(*,*) ' p_cn 1 ', p_cn(1,1:cr_nscat), n_cn, ncentral(1)
!write(*,*) ' p_cn 2 ', p_cn(2,1:cr_nscat), n_cn, ncentral(2)
!write(*,*) ' p_cn 3 ', p_cn(3,1:cr_nscat), n_cn, ncentral(3)
   cn_out: DO is = 0, cr_nscat 
      DO js = 0 , cr_nscat 
         IF(mmc_pair(ic, MC_COORDNUM, is, js) /=  0 ) THEN 
            mmc_ach_corr(ic, je, is, js) = REAL(n_cn        )/REAL(ncentral(is))
!           Feedback mechanism                                      
!           mmc_depth (ic, MC_COORDNUM, 0, 0) = mmc_depth (ic, MC_COORDNUM, 0, 0) - &
!           mmc_cfac (ic, MC_COORDNUM) * (mmc_target_corr (ic, MC_COORDNUM, is,js)- &
!                                            mmc_ach_corr (ic, MC_COORDNUM, is, js) ) / 2. &
!           *ABS(mmc_target_corr (ic, MC_OCC, is, js)) &
!           * damp
            mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
            mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_corr(ic, je, is, js) - mmc_ach_corr(ic, je, is, js)
            mmc_h_targ(mmc_h_ctarg)              = mmc_target_corr(ic, je, is, js)
            IF(lout) THEN
               IF(mmc_target_corr(ic, je, is, js)/=0.0) THEN
                  divider = mmc_target_corr(ic, je, is, js)
               ELSE
                  divider = 1.0
               ENDIF
               WRITE(output_io, 3150) ic, cr_at_lis(is), cr_at_lis(js),         &
                  mmc_target_corr(ic, je, is, js),                              &
                  mmc_ach_corr(ic, je, is, js),                                 &
                  mmc_target_corr(ic, je, is, js) - mmc_ach_corr(ic,je, is, js),&
                 (mmc_target_corr(ic, je, is, js) - mmc_ach_corr(ic,je, is, js))/divider,&
                  ncentral(is)
            ENDIF
         ENDIF
      ENDDO
   ENDDO cn_out
!                                                                       
!     ----- correlation of displacements                                
!                                                                       
disp_pair: DO is = 0, cr_nscat 
         DO js = is, cr_nscat 
         IF (mmc_pair (ic, MC_DISP, is, js) == -1 ) THEN 
         je = MC_DISP 
         xnn (is, js) = xnn (is, js) + xnn (js, is) 
         xij (is, js) = xij (is, js) + xij (js, is) 
         xi2 (is, js) = xi2 (is, js) + xi2 (js, is) 
         xj2 (is, js) = xj2 (is, js) + xj2 (js, is) 
         IF (xnn (is, js)  /= 0) THEN 
            xij (is, js) = xij (is, js) / REAL(xnn (is, js) ) 
            xi2 (is, js) = xi2 (is, js) / REAL(xnn (is, js) ) 
            xj2 (is, js) = xj2 (is, js) / REAL(xnn (is, js) ) 
!                                                                       
            IF (xi2 (is, js)  /= 0.AND.xj2 (is, js)  /= 0.0) THEN 
               mmc_ach_corr (ic, je, is, js) = xij (is, js) / sqrt (xi2 &
               (is, js) * xj2 (is, js) )                                
               mmc_ach_corr (ic, je, js, is) = mmc_ach_corr (ic, je, is,&
               js)                                                      
!               Feedback mechanism                                      
               mmc_depth (ic, MC_DISP, 0, 0) = mmc_depth (ic, MC_DISP, 0,0) - &
               mmc_cfac  (ic, MC_DISP) * (mmc_target_corr (ic, MC_DISP, is, js) - &
               mmc_ach_corr (ic, MC_DISP, is, js) ) / 2. &
                    *ABS(mmc_target_corr (ic, MC_DISP, is, js)) &
                    * damp
            ELSE 
               mmc_ach_corr (ic, je, is, js) = 0.0 
               mmc_ach_corr (ic, je, js, is) = 0.0 
            ENDIF 
         ELSE 
            mmc_ach_corr (ic, je, is, js) = 0.0 
            mmc_ach_corr (ic, je, js, is) = 0.0 
         ENDIF 
         mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
         mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_corr(ic, je, is, js) - mmc_ach_corr(ic, je, is, js)
            mmc_h_targ(mmc_h_ctarg)              = mmc_target_corr(ic, je, is, js)
!                                                                       
         IF (lout) THEN 
               IF(mmc_target_corr(ic, je, is, js)/=0.0) THEN
                  divider = mmc_target_corr(ic, je, is, js)
               ELSE
                  divider = 1.0
               ENDIF
            WRITE (output_io, 3200) ic, cr_at_lis (is), cr_at_lis (js), &
               mmc_target_corr (ic, je, is, js), mmc_ach_corr(ic, je, is,js), &
               mmc_target_corr (ic, je, is, js) - mmc_ach_corr (ic, je, is, js), &
               (mmc_target_corr(ic, je, is, js) - mmc_ach_corr(ic,je, is, js))/divider,&
               nneigh                                         
         ENDIF 
         ENDIF 
         ENDDO 
         ENDDO disp_pair
!                                                                       
!                                                                       
!     -- Loop over all atom pairs to do Hooke potential                 
!                                                                       
spri_pair: DO is = 0, cr_nscat 
         DO js = is, cr_nscat 
         IF (mmc_pair (ic, MC_SPRING, is, js) == -1 ) THEN 
         je = MC_SPRING 
!                                                                       
!     ----- Spring                                                      
!                                                                       
         IF (bl_anz (is, js)  /= 0.OR.bl_anz (js, is)  /= 0) THEN 
            mmc_ach_corr (ic, je, is, js) = (bl_sum (is, js) + bl_sum ( &
            js, is) ) / (bl_anz (is, js) + bl_anz (js, is) )            
            mmc_ach_sigm (ic, je, is, js) = (bl_s2 (is, js) + bl_s2 (js,&
            is) ) / (bl_anz (is, js) + bl_anz (js, is) )                
            mmc_ach_sigm (ic, je, is, js) = (mmc_ach_sigm (ic, je, is,  &
            js) - (mmc_ach_corr (ic, je, is, js) **2) )                 
            IF (mmc_ach_sigm (ic, je, is, js)  > 0) THEN 
               mmc_ach_sigm (ic, je, is, js) = sqrt (mmc_ach_sigm (ic,  &
               je, is, js) )                                            
            ELSE 
               mmc_ach_sigm (ic, je, is, js) = 0.0 
            ENDIF 
         mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
         mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_corr (ic, je, is, js) - mmc_ach_corr(ic, je, is, js)
            mmc_h_targ(mmc_h_ctarg)              = mmc_target_corr(ic, je, is, js)
!                                                                       
            IF (lout) THEN 
               IF(mmc_target_corr(ic, je, is, js)/=0.0) THEN
                  divider = mmc_target_corr(ic, je, is, js)
               ELSE
                  divider = 1.0
               ENDIF
               WRITE (output_io, 3300) ic, cr_at_lis (is),  cr_at_lis (js),   &
                  mmc_target_corr (ic, je, is, js),  mmc_ach_corr (ic, je, is, js),  &
                  mmc_ach_sigm (ic, je, is, js),                                     &
                  mmc_target_corr (ic, je, is, js)  - mmc_ach_corr (ic, je, is, js), &
                 (mmc_target_corr(ic, je, is, js) - mmc_ach_corr(ic,je, is, js))/divider,&
                  bl_anz (is, js)  + bl_anz (js, is) 
            ENDIF 
         ENDIF 
         ENDIF 
         ENDDO 
         ENDDO  spri_pair
!                                                                       
!     -- Loop over all defined angle correlations                       
!                                                                       
angl_pair: IF (mmc_cor_energy (ic, MC_ANGLE) ) THEN 
         je = MC_ANGLE 
         DO k = 1, mmc_n_angles 
         CALL index2angles (mmc_angles (k), iic, kk, iis, jjs, lls,     &
         MAXSCAT)
         IF (ba_anz (k)  /= 0) THEN 
            mmc_ach_angl (k) = (ba_sum (k) ) / (ba_anz (k) ) 
            mmc_ang_sigm (k) = (ba_s2 (k) ) / (ba_anz (k) ) 
            mmc_ang_sigm (k) = (mmc_ang_sigm (k) - (mmc_ach_angl (k) ** &
            2) )                                                        
            IF (mmc_ang_sigm (k)  > 0) THEN 
               mmc_ang_sigm (k) = sqrt (mmc_ang_sigm (k) ) 
            ELSE 
               mmc_ang_sigm (k) = 0.0 
            ENDIF 
            mmc_ach_corr (ic, je, jjs, lls) = mmc_ach_angl (k) 
            mmc_ach_sigm (ic, je, jjs, lls) = mmc_ang_sigm (k) 
         mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
         mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_angl (k               ) - mmc_ach_angl(k               )
            mmc_h_targ(mmc_h_ctarg)              = mmc_target_angl(k               )
!                                                                       
            IF (lout) THEN 
               IF (iic == ic) THEN 
               IF(mmc_target_angl(k               )/=0.0) THEN
                  divider = mmc_target_angl(k)
               ELSE
                  divider = 1.0
               ENDIF
               WRITE (output_io, 3400) ic, cr_at_lis (iis),  cr_at_lis (jjs), &
                  cr_at_lis (lls),  mmc_target_angl (k),  mmc_ach_angl (k),   &
                  mmc_ang_sigm (k),                                           &
                  mmc_target_angl (k)  - mmc_ach_angl (k),                    &
                 (mmc_target_angl (k)  - mmc_ach_angl (k))/divider,           &
                  ba_anz (k)  + ba_anz (k)                                                        
                  IF (iis ==  - 1) THEN 
                     IF (jjs ==  - 1) THEN 
                        IF (lls ==  - 1) THEN 
                           searching = .FALSE. 
                           ic_a = k 
                        ELSE 
                        ENDIF 
                     ENDIF 
                  ENDIF 
               ENDIF 
            ENDIF 
         ELSE 
            IF (lout) THEN 
               IF (iic == ic) THEN 
      WRITE (output_io, 3410) ic, cr_at_lis (iis),  cr_at_lis (jjs),  cr&
     &_at_lis (lls),  mmc_target_angl (k),  0.0, 0                           
                  ENDIF 
               ENDIF 
            ENDIF 
            ENDDO 
         ENDIF angl_pair
!                                                                       
!     -- Loop over all atom pairs to do Lennard Jones potential         
!                                                                       
lenn_pair: DO is = 0, cr_nscat 
      DO js = is, cr_nscat 
      IF (mmc_pair (ic, MC_LENNARD, is, js) == -1 ) THEN 
         je = MC_LENNARD 
!                                                                       
!     ----- Lennard                                                     
!                                                                       
         IF (bl_anz (is, js)  /= 0.OR.bl_anz (js, is)  /= 0) THEN 
            mmc_ach_corr (ic, je, is, js) = (bl_sum (is, js) + bl_sum ( &
            js, is) ) / (bl_anz (is, js) + bl_anz (js, is) )            
            mmc_ach_sigm (ic, je, is, js) = (bl_s2 (is, js) + bl_s2 (js,&
            is) ) / (bl_anz (is, js) + bl_anz (js, is) )                
            mmc_ach_sigm (ic, je, is, js) = (mmc_ach_sigm (ic, je, is,  &
            js) - (mmc_ach_corr (ic, je, is, js) **2) )                 
            IF (mmc_ach_sigm (ic, je, is, js)  > 0) THEN 
               mmc_ach_sigm (ic, je, is, js) = sqrt (mmc_ach_sigm (ic,  &
               je, is, js) )                                            
            ELSE 
               mmc_ach_sigm (ic, je, is, js) = 0.0 
            ENDIF 
         mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
         mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_corr (ic, je, is, js) - mmc_ach_corr(ic, je, is, js)
            mmc_h_targ(mmc_h_ctarg)              = mmc_target_corr(ic, je, is, js)
            IF (lout) THEN 
               IF(mmc_target_corr(ic, je, is, js)/=0.0) THEN
                  divider = mmc_target_corr(ic, je, is, js)
               ELSE
                  divider = 1.0
               ENDIF
               WRITE (output_io, 3700) ic, cr_at_lis (is),  cr_at_lis (js), &
                   mmc_target_corr(ic, je, is, js),  mmc_ach_corr(ic, je, is, js), &
                   mmc_ach_sigm(ic, je, is, js),                                   &
                   mmc_target_corr(ic, je, is, js) - mmc_ach_corr (ic, je, is, js),&
                  (mmc_target_corr(ic, je, is, js) - mmc_ach_corr (ic, je, is, js))/divider,&
                   bl_anz (is, js)  + bl_anz (js, is) 
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
      ENDDO lenn_pair
!                                                                       
!     -- Loop over all atom pairs to do Repulsive     potential         
!                                                                       
repu_pair: DO is = 0, cr_nscat 
      DO js = is, cr_nscat 
      IF (mmc_pair (ic, MC_REPULSIVE, is, js) == -1 ) THEN 
!write(*,*) is,js,ic, mmc_pair (ic, MC_REPULSIVE, is, js), bl_anz (is, js), bl_anz (js, is)
         je = MC_REPULSIVE 
!                                                                       
!     ----- REPULSIVE                                                     
!                                                                       
         IF (bl_anz (is, js)  /= 0.OR.bl_anz (js, is)  /= 0) THEN 
            mmc_ach_corr (ic, je, is, js) = (bl_sum (is, js) + bl_sum (js, is) ) &
                                          / (bl_anz (is, js) + bl_anz (js, is) )            
            mmc_ach_sigm (ic, je, is, js) = (bl_s2  (is, js) + bl_s2  (js, is) ) &
                                          / (bl_anz (is, js) + bl_anz (js, is) )                
            mmc_ach_sigm (ic, je, is, js) = (mmc_ach_sigm (ic, je, is, js)       &
                                          - (mmc_ach_corr (ic, je, is, js) **2) )                 
            IF (mmc_ach_sigm (ic, je, is, js)  > 0) THEN 
               mmc_ach_sigm (ic, je, is, js) = sqrt(mmc_ach_sigm(ic, je, is, js) )                                            
            ELSE 
               mmc_ach_sigm (ic, je, is, js) = 0.0 
            ENDIF 
         mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
         mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_corr (ic, je, is, js)  - mmc_ach_corr(ic, je, is, js)
            mmc_h_targ(mmc_h_ctarg)              = mmc_target_corr(ic, je, is, js)
            IF (lout) THEN 
               IF(mmc_target_corr(ic, je, is, js)/=0.0) THEN
                  divider = mmc_target_corr(ic, je, is, js)
               ELSE
                  divider = 1.0
               ENDIF
               WRITE (output_io, 3900) ic, cr_at_lis (is),  cr_at_lis (js),       &
               mmc_target_corr (ic, je, is, js),  mmc_ach_corr (ic, je, is, js),  &
               mmc_ach_sigm (ic, je, is, js),                                     &
               mmc_target_corr (ic, je, is, js)  - mmc_ach_corr (ic, je, is, js), &
              (mmc_target_corr (ic, je, is, js)  - mmc_ach_corr (ic, je, is, js))/divider, &
               bl_anz (is, js)  + bl_anz (js, is) 
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
         ENDDO repu_pair      
!                                                                       
!     -- Loop over all atom pairs to do Buckingham potential            
!                                                                       
   buck_pair: DO is = 0, cr_nscat 
      DO js = is, cr_nscat 
         IF (mmc_pair (ic, MC_BUCKING, is, js) == -1 ) THEN 
            je = MC_BUCKING 
!                                                                       
!     ----- Buckingham                                                  
!                                                                       
            IF (bl_anz (is, js)  /= 0.OR.bl_anz (js, is)  /= 0) THEN 
               mmc_ach_corr (ic, je, is, js) =              &
                     (bl_sum (is, js) + bl_sum (js, is) ) / &
                     (bl_anz (is, js) + bl_anz (js, is) )            
               mmc_ach_sigm (ic, je, is, js) =              &
                     (bl_s2  (is, js) + bl_s2  (js, is) ) / &
                     (bl_anz (is, js) + bl_anz (js, is) )                
               mmc_ach_sigm (ic, je, is, js) =              &
                     (mmc_ach_sigm (ic, je, is, js) -       &
                     (mmc_ach_corr (ic, je, is, js) **2) )                 
               IF (mmc_ach_sigm (ic, je, is, js)  > 0) THEN 
                  mmc_ach_sigm (ic, je, is, js) =           &
                      sqrt (mmc_ach_sigm (ic, je, is, js) )                                            
               ELSE 
                  mmc_ach_sigm (ic, je, is, js) = 0.0 
               ENDIF 
         mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
         mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_corr (ic, je, is, js) - mmc_ach_corr(ic, je, is, js)
            mmc_h_targ(mmc_h_ctarg)              = mmc_target_corr(ic, je, is, js)
               IF (lout) THEN 
               IF(mmc_target_corr(ic, je, is, js)/=0.0) THEN
                  divider = mmc_target_corr(ic, je, is, js)
               ELSE
                  divider = 1.0
               ENDIF
                  WRITE (output_io, 2100) ic, cr_at_lis(is), cr_at_lis(js),  &
                     mmc_ach_corr(ic, je, is, js), mmc_ach_sigm(ic, je,is, js), &
                     bl_anz (is, js) + bl_anz (js, is)               
               ENDIF 
            ENDIF 
         ENDIF 
      ENDDO 
   ENDDO buck_pair
ENDDO main_corr
!
if(lfeed) then
!          mmc_pid_pid(1,2) = mmc_pid_pid(1,2) + sum(mmc_pid_deri(ic, MC_OCC, is, js, :)) ! Derivat. PID ==> P
!          mmc_pid_pid(2,2) = mmc_pid_pid(2,2) + mmc_pid_diff(ic, MC_OCC, is, js)         ! Integral PID ==> I
!          mmc_pid_pid(3,2) = mmc_pid_pid(3,2) + change                                   ! Change       ==> D
!          mmc_pid_pid_n    = mmc_pid_pid_n + 1
   if(mmc_pid_pid_n > 0 ) then
      mmc_pid_pid(1,2) = mmc_pid_pid(1,2) / mmc_pid_pid_n
      mmc_pid_pid(2,2) = mmc_pid_pid(2,2) / mmc_pid_pid_n
      mmc_pid_pid(3,2) = mmc_pid_pid(3,2) / mmc_pid_pid_n
   endif
!          mmc_pid_pid_n    = mmc_pid_pid_n + 1
   if(mmc_pid_pid(1,2)>0.0) then
      mmc_pid_pid(1,1) = min(mmc_pid_pid(1,1) * 1.05D0, 1.000D0)
   else
      mmc_pid_pid(1,1) = max(mmc_pid_pid(1,1) * 0.90D0, 0.001D0)
   endif
   if(mmc_pid_pid(2,2) > 0) then
      mmc_pid_pid(2,1) = min(mmc_pid_pid(2,1) * 1.025D0, 1.000D0)
   else
      mmc_pid_pid(2,1) = max(mmc_pid_pid(2,1) * 0.950D0, 0.001D0)
   endif
   if(abs(mmc_pid_pid(3,2)) > mmc_pid_change) then
      mmc_pid_pid(3,1) = min(mmc_pid_pid(3,1) * 0.950D0, 1.000D0)
   else
      if(mmc_pid_pid(3,1)>0.0) then
         mmc_pid_pid(3,1) = max(mmc_pid_pid(3,1) * 1.025D0, 0.001D0)
      elseif(mmc_pid_pid(3,1)<0.0) then
         mmc_pid_pid(3,1) = min(mmc_pid_pid(3,1) * 1.025D0,-0.001D0)
      endif
   endif
   mmc_pid_change = abs(mmc_pid_pid(3,2))
endif
!
!
DEALLOCATE(ncentral)
DEALLOCATE(patom)
DEALLOCATE(iatom)
DEALLOCATE(tatom)
DEALLOCATE(natom)
DEALLOCATE(xnn)
DEALLOCATE(xij)
DEALLOCATE(xi2)
DEALLOCATE(xj2)
DEALLOCATE(ba_sum)
DEALLOCATE(ba_s2)
DEALLOCATE(ba_anz)
DEALLOCATE( bl_anz)
DEALLOCATE( bl_sum)
DEALLOCATE( bl_s2)
DEALLOCATE( p_cn)
!DEALLOCATE( pneig)
!
!  Check for convergence
!
done = .FALSE.
IF(mmc_h_stop) THEN     ! Apply convergence criteria
   IF(.NOT.lfinished .AND.mmc_h_ctarg>0) THEN
!
      mmc_h_aver  = 0.0                     ! Initialize average changes in difference (target -achieved)
      mmc_h_aver_r= 0.0                     ! Initialize maximum changes in difference (target -achieved)
      if(mmc_m_index>0) then                ! We need at least feedback cycles 0 and 1
         conv_val = 0.0
         do i = 1, mmc_h_number             ! Loop over all targets
            do is = 0,mmc_m_index -1        ! Loop over all feedback cycles but one
               j = mod(is+1, MMC_H_NNNN)   ! Next cycle after is
               change = (mmc_h_diff(i,j) - mmc_h_diff(i,is))
               if(mmc_h_targ(i) /= 0.0)  change = change / mmc_h_targ(i)   ! Relative change
               mmc_h_aver(i) = mmc_h_aver(i) + change                      ! Accumulate all changes feedback to feedback
               mmc_h_aver_r(i) = max(mmc_h_aver_r(i), abs(change))             ! Determine maximum absolute change
            enddo
            mmc_h_aver = mmc_h_aver / mmc_m_index
            if(mmc_h_targ(i)/=0.0) then
               conv_val(2) = max(conv_val(2), abs(mmc_h_diff(i,mmc_h_index))/abs(mmc_h_targ(i)))
            endif
         enddo
         conv_val(1) = maxval(abs(mmc_h_diff(1:mmc_h_number,mmc_h_index)))   ! Maximum difference (target-achieved)
!        conv_val(2) = maxval(abs(mmc_h_diff(1:mmc_h_number,mmc_h_index))/abs(mmc_h_targ(1:mmc_h_number)), &
!                             mmc_h_targ(1:mmc_h_number) > 0.0)              ! Maximum difference (target-achieved)/target
         conv_val(3) = maxval(mmc_h_aver(1:mmc_h_number))                    ! Average change between feedbacks
         conv_val(4) = maxval(mmc_h_aver_r(1:mmc_h_number))                    ! Maximum change between feedbacks
      else
         conv_val = 1.0e20                  ! No convergence yet
      endif
!write(output_io, '(a,5i5)    ') ' INDEX               ', mmc_h_index, mmc_m_index, mmc_h_number, MMC_H_NNNN, mmc_h_nfeed
!do is = 0, mmc_m_index                                     ! Only for actual feedback cycles
!  WRITE(output_io,'(a,i3,10F8.4)') 'Target-Achieved      ', is, abs(mmc_h_diff(1:mmc_h_number,is)) ,  &
!                                                            maxval(abs(mmc_h_diff(1:mmc_h_number,is))),  &
!                                                            conv_val(1), mmc_h_conv_m
!enddo
!do is = 0, mmc_m_index                                     ! Only for actual feedback cycles
!  change = 0
!  do i = 1, mmc_h_number
!    change = max(change, abs(mmc_h_diff(i,is)/mmc_h_targ(i)))
!  enddo
!  WRITE(output_io,'(a,i3,10F8.4)') 'Target-Achieved rel  ', is, (abs(mmc_h_diff(i,is)/mmc_h_targ(i)),i=1,mmc_h_number), &
!                                                            change, conv_val(2), mmc_h_conv_r
!enddo
!write(output_io,'(a,1x,10F8.4)') 'Average changes        ', mmc_h_aver(1:mmc_h_number),  conv_val(3), mmc_h_conv_a
!write(output_io,'(a,1x,10F8.4)') 'Maximum changes        ', mmc_h_aver_r(1:mmc_h_number),conv_val(4), mmc_h_conv_c
!write(output_io,'(a,4l3)')       'Convergence             ', &
!          conv_val(1) < mmc_h_conv_m ,     conv_val(2) < mmc_h_conv_r ,         &
!          conv_val(3) < mmc_h_conv_a ,     conv_val(4) < mmc_h_conv_c 
maxdev(1) = conv_val(1)
maxdev(2) = conv_val(2)
!write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
!  DO i=1, mmc_h_number                     ! Loop over all targets
!     j = MOD(mmc_h_index+2,MMC_H_NNNN)     ! Increment Feedback number
!     mmc_h_aver(i) = mmc_h_diff(i,0)       ! initialize average deviations over MMC_H_NNNN feedbacks
!     DO is = 1, MMC_H_NNNN-1
!        mmc_h_aver(i) = mmc_h_aver(i) + mmc_h_diff(i,is)
!     ENDDO
!     mmc_h_aver  (i) = mmc_h_aver(i)/MMC_H_NNNN
!     mmc_h_aver_r(i) = mmc_h_aver(i)/MMC_H_NNNN
!     mmc_h_maxd(i,mmc_h_index) = (mmc_h_diff(i,mmc_h_index) - mmc_h_diff(i,j))
!     IF(mmc_h_targ(i) /= 0.0) THEN         ! Normalize relative to target value
!        mmc_h_aver_r(i) = ABS(mmc_h_aver_r(i)/mmc_h_targ(i))
!        mmc_h_maxd(i,mmc_h_index) = (mmc_h_maxd(i,mmc_h_index)/mmc_h_targ(i))
!     ENDIF
!  ENDDO
!mmc_h_ctarg
   ENDIF
!write(*,*) mmc_h_nfeed , mmc_h_nfeed > 2
   done = mmc_h_nfeed >= MMC_H_NNNN                                   .and.     &
          conv_val(1) < mmc_h_conv_m .and. conv_val(2) < mmc_h_conv_r .and.     &
          conv_val(3) < mmc_h_conv_a .and. conv_val(4) < mmc_h_conv_c 
!         MAXVAL(    mmc_h_aver  (1:mmc_h_ctarg))    < mmc_h_conv_m .AND.     &
!         MAXVAL(    mmc_h_aver_r(1:mmc_h_ctarg))    < mmc_h_conv_r .AND.     &
!         MAXVAL(ABS(mmc_h_maxd  (1:mmc_h_ctarg,:))) < mmc_h_conv_c .AND.     &
!         ABS(SUM(mmc_h_maxd(1:mmc_h_ctarg,0:MMC_H_NNNN-1)))/               &
!            (REAL(mmc_h_ctarg*MMC_H_NNNN))   < mmc_h_conv_a
!if(mmc_h_ctarg>0) then
!WRITE(output_io,'(a,10F8.4)') ' Average deviations    ', mmc_h_aver  (1:mmc_h_ctarg), mmc_h_conv_m
!WRITE(output_io,'(a,10F8.4)') ' Average rel. deviation', mmc_h_aver_r(1:mmc_h_ctarg), mmc_h_conv_r
!DO is = 0, mmc_m_index                                     ! Only for actual feedback cycles
!  WRITE(output_io,'(a,10F8.4)') ' Differences           ', (mmc_h_maxd(1:mmc_h_ctarg,is))
!enddo
!WRITE(output_io,'(a,10F8.4)') ' Stagnant    largest   ', MAXVAL(ABS(mmc_h_maxd(1:mmc_h_ctarg,:))), &
! mmc_h_conv_c
!WRITE(output_io,'(a,10F8.4)') ' Stagnant    average   ', &
!SUM(mmc_h_maxd(1:mmc_h_ctarg,0:MMC_H_NNNN-1))/(REAL(mmc_h_ctarg*MMC_H_NNNN)), mmc_h_conv_a
!write(output_io,* ) ' DONE, finished ', done,lfinished, rel_cycl, mmc_h_ctarg
!endif
   IF(lout) then
   if(lfinished) THEN
      WRITE(output_io,*)
      IF(done) THEN
         WRITE(output_io,'(a,f6.0,a,i3,a)') ' Convergence reached after ',      &
            rel_cycl*100., ' % of user cycles',                                 &
            MMC_H_NNNN, ' Feedbacks averaged'
         WRITE(output_io,'(a,f12.4,a,f12.4)') ' Largest          difference ',   &
            conv_val(1), ' < ', mmc_h_conv_m
         WRITE(output_io,'(a,f12.4,a,f12.4)') ' Largest relative difference ',   &
            conv_val(2), ' < ', mmc_h_conv_r
         WRITE(output_io,'(a,f12.4,a,f12.4)') ' Stagnant: average changes   ',   &
            conv_val(3), ' < ', mmc_h_conv_a
         WRITE(output_io,'(a,f12.4,a,f12.4)') ' Stagnant: largest changes   ',   &
            conv_val(4), ' < ', mmc_h_conv_c
      ELSEIF(rel_cycl==1.0) THEN
         WRITE(output_io,*) 'Maximum Cycles reached '
         WRITE(output_io,'(a,f12.4,a,f12.4)') ' Largest          difference ',   &
            conv_val(1), ' ? ', mmc_h_conv_m
         WRITE(output_io,'(a,f12.4,a,f12.4)') ' Largest relative difference ',   &
            conv_val(2), ' ? ', mmc_h_conv_r
         WRITE(output_io,'(a,f12.4,a,f12.4)') ' Stagnant: average changes   ',   &
            conv_val(3), ' ? ', mmc_h_conv_a
         WRITE(output_io,'(a,f12.4,a,f12.4)') ' Stagnant: largest changes   ',   &
            conv_val(4), ' ? ', mmc_h_conv_c
      ENDIF
      else
         if(mmc_h_log) then
         WRITE(output_io,'(a,f6.0,a,i3,a)') ' Convergence status  after ',      &
            rel_cycl*100., ' % of user cycles',                                 &
            MMC_H_NNNN, ' Feedbacks averaged'
         WRITE(output_io,'(a,f12.4,a,f12.4)') ' Largest          difference ',   &
            conv_val(1), ' < ', mmc_h_conv_m
         WRITE(output_io,'(a,f12.4,a,f12.4)') ' Largest relative difference ',   &
            conv_val(2), ' < ', mmc_h_conv_r
         WRITE(output_io,'(a,f12.4,a,f12.4)') ' Stagnant: average changes   ',   &
            conv_val(3), ' < ', mmc_h_conv_a
         WRITE(output_io,'(a,f12.4,a,f12.4)') ' Stagnant: largest changes   ',   &
            conv_val(4), ' < ', mmc_h_conv_c
         endif
      ENDIF
   ENDIF
ENDIF
!! PRINT FOR res_para
if(lfinished) THEN
res_para(0) = 0
loop_cor: DO ic = 1, CHEM_MAX_COR 
   loop_ener: DO je = 1, MC_N_ENERGY 
      cond_ener: if(mmc_cor_energy (ic, je)) then   ! This energy is set
         lfirst = .true.
         loop_is: DO is =  -1, MAXSCAT 
            loop_js: DO js  =  -1, MAXSCAT 
               cond_pair: if(mmc_pair(ic, je, is, js)>-3 .and. mmc_pair(ic, je, is, js)<0 ) then
                  if((je==MC_OCC      .and. js> is  .and. is>=0  .and. lfirst   ) .or. & ! 1 Chemical correlation
                     (je==MC_DISP     .and. js>=is  .and. is>=0                 ) .or. & ! 2 Displacement correlation
                     (je==MC_SPRING   .and. mmc_ach_corr(ic, je, is, js)/=0.0D0 ) .or. & ! 3 Hooke
                     (je==MC_LENNARD  .and. mmc_ach_corr(ic, je, is, js)/=0.0D0 ) .or. & ! 7 Lennard Jones
                     (je==MC_BUCKING  .and. mmc_ach_corr(ic, je, is, js)/=0.0D0 ) .or. & ! 8 Buckingham
                     (je==MC_REPULSIVE.and. mmc_ach_corr(ic, je, is, js)/=0.0D0 ) .or. & ! 9 Repulsive
                     (je==MC_COORDNUM .and. js> is  .and. is>=0  .and. lfirst   ) .or. & !10 Coordination number
                     (je==MC_UNI      .and. js> is  .and. is>=0  .and. lfirst   ) .or. & !11 Unidirectional 
                     (je==MC_GROUP    .and. is/=js .and. is>0 .and. js>0 .and.         & !12 Group correlations
                                            mmc_cor_energy(0, MC_GROUP)   .and.        &
                                            mmc_ach_corr(ic, je, is, js)/=0.0D0 .and.  & 
                                            mmc_target_corr(ic, je, is, js)/=0.0D0)    & !12 Group Correlations
                    )  then 
!
! write(*,'(4i5, f10.4, i5)') ic, je,    (is),    (js), mmc_ach_corr (ic, je, is, js), &
!mmc_pair(ic, je, is, js)
                    res_para(0) = res_para(0) + 1
                    res_para(nint(res_para(0))) = mmc_ach_corr(ic, je, is, js)
                    lfirst = .false.
!
                  endif
               endif cond_pair
            ENDDO loop_js
         ENDDO loop_is
         cond_angl: if(je==MC_ANGLE) then                                               ! 4 Angular correlation
                     do k = 1, mmc_n_angles
!write(*,'(4i5, f10.4, i9)') ic, je, 0,0, mmc_ach_angl(k), mmc_pair(ic, je, -1, -1)
               res_para(0) = res_para(0) + 1
               res_para(nint(res_para(0))) = mmc_ach_angl(k)
            enddo
         endif cond_angl
      endif cond_ener
   ENDDO loop_ener
ENDDO loop_cor
endif
!                                                                       
  410 FORMAT ( 45x,'Correlations/',/                                    &
     &   ' Neig.- Energy-',7x,'Atoms',11x,'Target',2x,'Distance/',4x,   &
     &    'Sigma',5x,'Diff',6x,'Diff/',5x,'Number',/                               &
     &    ' Def.   Type    central  Neighbors',13x,'Angle'              &
     &   ,27x,'Target  of pairs')                                               
 2100 FORMAT (1x,i3,3x,a9,3x,a9,5x,f7.3,3x,f7.3,3x,i8) 
 3100 FORMAT (1x,i3,3x,'Occupancy',a5,3x,a5,      8x,2(f7.3,3x),        &
     &        10x,f7.3,3x,i8)
 3150 FORMAT (1x,i3,3x,'Coord.No.',a5,3x,a5,      8x,2(f7.3,3x),        &
     &        10x,f7.3,3x,f7.3,3x,i8)
 3200 FORMAT (1x,i3,3x,'Disp.Cor.',a5,3x,a5,      8x,2(f7.3,3x),        &
     &        10x,f7.3,3x,f7.3,3x,i8)                                           
 3300 FORMAT (1x,i3,3x,'Hooke    ',a5,3x,a5,      8x,5(f7.3,3x)         &
     &             ,i8)                                                 
 3400 FORMAT (1x,i3,3x,'Angle    ',a5,3x,a5,2x,a5,1x,5(f7.3,3x)         &
     &             ,i8)                                                 
 3410 FORMAT (1x,i3,3x,'Angle    ',a5,3x,a5,2x,a5,1x,  f7.3,23x         &
     &             ,f7.3,3x,i8)                                                 
 3700 FORMAT (1x,i3,3x,'Lennard  ',a5,3x,a5,      8x,5(f7.3,3x)         &
     &             ,i8)                                                 
 3900 FORMAT (1x,i3,3x,'Repuls.  ',a5,3x,a5,      8x,5(f7.3,3x)         &
     &             ,i8)                                                 
!                                                                       
END SUBROUTINE mmc_correlations               
!
!*******************************************************************************
!
SUBROUTINE mmc_correlations_occ(ic, pneig, rel_cycl, damp, lout, lfeed, MAXSCAT_L, &
                                MAX_COR, maxdev)
!
!     ----- Chemical correlation                                        
!
USE crystal_mod 
USE mc_mod 
USE mmc_mod 
!
USE prompt_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: ic
INTEGER, INTENT(IN) :: MAXSCAT_L
INTEGER, INTENT(IN) :: MAX_COR
INTEGER, DIMENSION(0:MAXSCAT_L, 0:MAXSCAT_L, 1:MAX_COR) , INTENT(IN) :: pneig
REAL(kind=PREC_DP)   , INTENT(IN) :: rel_cycl ! Relative progress along cycles
REAL(kind=PREC_DP)   , INTENT(IN) :: damp
LOGICAL, INTENT(IN) :: lout
LOGICAL, INTENT(IN) :: lfeed
real(kind=PREC_DP), dimension(2), intent(inout) :: maxdev
!
!
INTEGER :: pair11
INTEGER :: pair12
INTEGER :: pair21
INTEGER :: pair22
INTEGER :: is, js, je
INTEGER :: nneigh
LOGICAL :: lfirst
REAL(kind=PREC_DP) :: prob11, prob12, prob22 
REAL(PREC_DP) :: thet !, thet2
REAL(PREC_DP) :: divisor
REAL(PREC_DP) :: change
integer :: k
!
je = MC_OCC 
prob11 = 0.0
prob12 = 0.0
prob22 = 0.0
pair11 = 0
pair12 = 0
pair21 = 0
pair22 = 0
thet   = 0.0
k = 0
DO is = 0, cr_nscat 
!IF(MAxVAL(mmc_pair (ic, MC_OCC, is, 0:cr_nscat))>0 .OR.   &
!   minval(mmc_pair (ic, MC_OCC, is, 0:cr_nscat))<0      ) then
!write(*,*) 'PAIRS is= : ',is,' :: ',pneig(is, 0:cr_nscat, ic), ' ::', ic, ' :: ', mmc_pair (ic, MC_OCC, is, 0:cr_nscat)
!endif
   DO js =  0, cr_nscat 
      IF     (mmc_pair (ic, MC_OCC, is, js) == -1 ) THEN 
         pair12 = pair12 + pneig (is,js, ic)
      ELSEIF (mmc_pair (ic, MC_OCC, is, js) == -2 ) THEN 
         pair21 = pair21 + pneig (is,js, ic)
      ELSEIF (mmc_pair (ic, MC_OCC, is, js) == +1 ) THEN 
         pair11 = pair11 + pneig (is,js, ic)
      ELSEIF (mmc_pair (ic, MC_OCC, is, js) == +2 ) THEN 
         pair22 = pair22 + pneig (is,js, ic)
      ENDIF
   ENDDO
ENDDO
je = MC_OCC 
!write(*,'(a,4i8)') 'PAIRS ', pair11, pair22, pair12, pair21
!                                                                       
nneigh = pair11 + pair12 + pair21 + pair22 
IF (nneigh > 0.) THEN 
   prob11 =  pair11           / REAL(nneigh) 
   prob12 = (pair12 + pair21) / REAL(nneigh) 
   prob22 =  pair22           / REAL(nneigh) 
   thet = 0.5 * (2.0 * pair11 + pair12 + pair21) / REAL(nneigh)                                                     
!  thet2= 0.5 * (2.0*prob11 + prob12)
ENDIF 
!write(*,'(a,4i8,  i8)'      ) 'PAIRS ', pair11, pair12, pair21, pair22, nneigh
!write(*,'(a,3f7.3,3x,2f7.3)') 'PROBs ', prob11, prob22, prob12, thet, thet2 !0.5 * ((pair22 + pair11) + pair12 + pair21) / REAL(nneigh)
!write(*,'(a,2F7.3)'         ) 'corrs ', (prob11 - thet**2) / (thet * (1 - thet)) , &
!                                        (prob11 - thet2**2) / (thet2 * (1 - thet2))
!
lfirst = .TRUE.
corr_pair: DO is = 0, cr_nscat 
   DO js = is, cr_nscat 
      IF     (mmc_pair (ic, MC_OCC, is, js) /=  0 ) THEN 
        IF (thet /= 0.0.AND.thet /= 1.0) THEN 
           mmc_ach_corr (ic, je, is, js) = (prob11 - thet**2) /&
                                           (thet * (1 - thet) )
           mmc_ach_corr (ic, je, js, is) = (prob11 - thet**2) /&
                                           (thet * (1 - thet) )
        ELSE 
           IF((pair11>0 .OR. pair22>0) .AND. (pair12==0 .AND. pair21==0)) THEN
              mmc_ach_corr (ic, je, is, js) = 1.0 
              mmc_ach_corr (ic, je, js, is) = 1.0 
           ELSEIF((pair12>0 .OR. pair12>0) .AND. (pair11==0 .AND. pair22==0)) THEN
              mmc_ach_corr (ic, je, is, js) =-1.0 
              mmc_ach_corr (ic, je, js, is) =-1.0 
           ELSE 
              mmc_ach_corr (ic, je, is, js) = 0.0 
              mmc_ach_corr (ic, je, js, is) = 0.0 
           ENDIF 
        ENDIF 
!               Feedback mechanism                                      
        IF(mmc_target_corr (ic, MC_OCC, is, js) /= 0.0) THEN
           divisor = ABS(mmc_target_corr (ic, MC_OCC, is, js))
        ELSE
           divisor = 1.0
        ENDIF
        IF_FEED: IF(lfeed) THEN
           IF_lfeed: IF(mmc_lfeed(ic,MC_OCC)) THEN
           IF_RELC: IF(rel_cycl>0.0) THEN
              call calc_change_pid(ic, MC_OCC, is, js, damp, divisor, 1.0D0, maxdev, change)
!
!          change = mmc_cfac(ic, MC_OCC) * (mmc_target_corr(ic, MC_OCC, is, js) -      &
!                                           mmc_ach_corr   (ic, MC_OCC, is, js) ) / 0.1& !divisor &
!                  *ABS(mmc_target_corr(ic, MC_OCC, is, js)) * damp
!ENAGE
               if(abs(mmc_target_corr(ic, MC_OCC, is, js))<0.001D0) then
!write(*,*) ' Zerol ', change, change*0.0001, mmc_target_corr(ic, MC_OCC, is, js)
                  mmc_depth(ic, MC_OCC, is, js) = 0.0D0
                  mmc_depth(ic, MC_OCC, js, is) = mmc_depth (ic, MC_OCC, is, js)
                  change = 0.0D0
               elseif(abs(mmc_target_corr(ic, MC_OCC, is, js))<0.10) then
!write(*,*) ' Small ', change, change*0.005
                   change = change*0.005
               elseif(abs(mmc_target_corr(ic, MC_OCC, is, js))<0.2) then
!write(*,*) ' SMALL ', change, change*0.05
                   change = change*0.05
               elseIF(mmc_target_corr(ic, MC_OCC, is, js)*mmc_ach_corr(ic, MC_OCC, is, js)>=0.0  .AND. &
                  mmc_depth(ic, MC_OCC, is, js)*(mmc_depth(ic, MC_OCC, is, js)+change)<0.0         ) THEN
!write(*,*) ' CHANGE ', change, -mmc_depth(ic, MC_OCC, is, js)*0.100
!                  change = change*0.05
!                  change = -mmc_depth(ic, MC_OCC, is, js)*0.005
                   change = -mmc_depth(ic, MC_OCC, is, js)*0.100
               ENDIF
!
               mmc_depth(ic, MC_OCC, is, js) = mmc_depth (ic, MC_OCC, is, js) + change
               mmc_depth(ic, MC_OCC, js, is) = mmc_depth (ic, MC_OCC, is, js)
            ENDIF IF_RELC
            ENDIF IF_lfeed
         ENDIF IF_FEED
!                                                                       
!if(ic==1 .and. rel_cycl==0.0) then
!write(88,'(5g18.6e3)') rel_cycl, mmc_ach_corr(ic, je,1,2),mmc_depth(ic, MC_OCC, 1,2), change, damp
!endif
         IF (lout .AND. mmc_pair(ic,MC_OCC,is,js) < 0 .AND. lfirst) THEN
!if(ic==1 .and. lfeed) then
!write(88,'(5g18.6e3)') rel_cycl, mmc_ach_corr(ic, je,1,2),mmc_depth(ic, MC_OCC, 1,2), change, damp
!endif
            mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
            mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_corr(ic, je, is, js) - &
                                                   mmc_ach_corr   (ic, je, is, js)
            mmc_h_targ(mmc_h_ctarg)              = mmc_target_corr(ic, je, is, js)
            lfirst = .FALSE.
            WRITE (output_io, 3100) ic, cr_at_lis (is), cr_at_lis (js),         &
                mmc_target_corr(ic, je, is, js),                                &
                mmc_ach_corr   (ic, je, is, js),                                &
                mmc_target_corr(ic, je, is, js) - mmc_ach_corr (ic,je, is, js), &
               (mmc_target_corr(ic, je, is, js) - mmc_ach_corr (ic,je, is, js))/divisor, &
                nneigh! , &
!                mmc_depth(ic, MC_OCC, is, js), change
!write(*,'(a,f16.8)') ' depth ', mmc_depth(ic, MC_OCC, is,js)
!                                                                     
         ENDIF 
!                                                                       
      ENDIF 
   ENDDO 
ENDDO corr_pair
!
 3100 FORMAT (1x,i3,3x,'Occupancy',a5,3x,a5,      8x,2(f7.3,3x),        &
     &        10x,f7.3,3x,f7.3,3x,i8) !, 2f9.5)
!
END SUBROUTINE mmc_correlations_occ
!
!*****7*************************************************************************
!
SUBROUTINE mmc_correlation_write!(lcurrent)
!-                                                                      
!     Determines the achieved correlations                              
!                                                                       
!+                                                                      
USE crystal_mod 
USE chem_mod 
!USE chem_menu
!USE chem_aver_mod
!USE chem_neig_multi_mod
!USE atom_env_mod
!USE celltoindex_mod
!USE metric_mod
!USE mc_mod 
USE mmc_mod 
!
!USE debug_mod 
!USE errlist_mod 
use precision_mod
USE prompt_mod 
!
IMPLICIT none 
!
!logical, intent(in) :: lcurrent
!
integer :: ic
integer :: is, js, je
logical :: lfirst
real(kind=PREC_DP) :: divisor
!
write(output_io,'(a)') ' --- Initial multiple energy configuration ---'
write(output_io,*)
write(output_io,410)
!
main_corr: do ic = 1, chem_ncor 
!
!  OCCUPATION CORRELATION 
!
   je = MC_OCC 
   call mmc_correlation_write_corr(ic, je, 'Occupancy', .false.)
!
!  Unidirectional OCCUPATION CORRELATION 
!
   je = MC_UNI 
   call mmc_correlation_write_corr(ic, je, 'Unidirect', .false.)
!
!  Group OCCUPATION CORRELATION 
!
   je = MC_GROUP 
   lfirst = .TRUE.
   IF (mmc_cor_energy (0, MC_GROUP)) THEN
   is = MAXLOC(mmc_left( ic, MC_GROUP,:), 1) - 1
   js = MAXLOC(mmc_right(ic, MC_GROUP,:), 1) - 1
   IF(mmc_target_corr(ic, je, is, js) /= 0.0) then
      divisor = ABS(mmc_target_corr(ic, je, is, js))
   ELSE
      divisor = 1.0
   ENDIF
   write(output_io, 3100) ic, 'Group cor', cr_at_lis (is), cr_at_lis (js),         &
      mmc_target_corr(ic, je, is, js),                                &
      mmc_ini_corr   (ic, je, is, js),                                &
      mmc_target_corr(ic, je, is, js) - mmc_ini_corr (ic,je, is, js), &
      (mmc_target_corr(ic, je, is, js) - mmc_ini_corr (ic,je, is, js))/divisor
   endif
!
!  Unidirectional OCCUPATION CORRELATION 
!
   je = MC_COORDNUM
   call mmc_correlation_write_corr(ic, je, 'Coor num ', .false.)
!
!  Displacement correlations
!
   je = MC_DISP 
   call mmc_correlation_write_corr(ic, je, 'Disp.Cor.', .true.)
!
!  Angle potentials
!
   je = MC_ANGLE 
   call mmc_correlation_write_angl(ic, je, 'Angular  ', .true.)
!
!  Spring       correlations
!
   je = MC_SPRING 
   call mmc_correlation_write_pote(ic, je, 'Hooke    ', .true.)
!
!  Lennard Jones Potential s
!
   je = MC_LENNARD
   call mmc_correlation_write_pote(ic, je, 'Lennard  ', .true.)
!
!  Repulsive Potential
!
   je = MC_REPULSIVE
   call mmc_correlation_write_pote(ic, je, 'Repulsive', .true.)
!
enddo main_corr
!
write(output_io,*)
write(output_io,'(a)') ' ---------------------------------------------'
!
!                                                                       
  410 FORMAT ( 45x,'Correlations/',/                                    &
     &   ' Neig.- Energy-',7x,'Atoms',11x,'Target',2x,'Distance/',4x,   &
     &    'Sigma',5x,'Diff',6x,'Diff/',5x,'Number',/                               &
     &    ' Def.   Type    central  Neighbors',13x,'Angle'              &
     &   ,27x,'Target  of pairs')                                               
!
 3100 FORMAT (1x,i3,3x,a9,a5,3x,a5,      8x,2(f7.3,3x),        &
     &        10x,f7.3,3x,f7.3,3x,i8)
!
end subroutine mmc_correlation_write
!
!*******************************************************************************
!
subroutine mmc_correlation_write_corr(ic, je, title, luse_all)
!-
!  Write a specific correlation for chemical coordination
!+
!
use crystal_mod
use chem_mod
use mmc_mod
!
use precision_mod
use prompt_mod
!
implicit none
!
integer         , intent(in) :: ic         ! Correlations number
integer         , intent(in) :: je         ! Energy type
character(len=*), intent(in) :: title      ! Energy name
logical         , intent(in) :: luse_all ! Do not exist at first pair
!
integer :: is, js
logical :: lfirst
real(kind=PREC_DP) :: divisor
!
cond_ener: IF(mmc_cor_energy (ic, je) ) THEN 
   lfirst = .TRUE.
   loop_type: do is = 0, cr_nscat 
      do js = is, cr_nscat 
         if(mmc_pair(ic,je,is,js) < 0 .AND. lfirst) then
            if(mmc_target_corr (ic, je, is, js) /= 0.0) then
               divisor = ABS(mmc_target_corr (ic, je, is, js))
            else
               divisor = 1.0
            endif
            if     (mmc_pair (ic, je, is, js) /=  0 ) then 
            lfirst = luse_all
            write(output_io, 3100) ic, title      , cr_at_lis (is), cr_at_lis (js),         &
                mmc_target_corr(ic, je, is, js),                                &
                mmc_ini_corr   (ic, je, is, js),                                &
                mmc_target_corr(ic, je, is, js) - mmc_ini_corr (ic,je, is, js), &
               (mmc_target_corr(ic, je, is, js) - mmc_ini_corr (ic,je, is, js))/divisor
            endif
         endif
      enddo
   enddo loop_type
endif cond_ener
!
 3100 FORMAT (1x,i3,3x,a9,a5,3x,a5,      8x,2(f7.3,3x),        &
     &        10x,f7.3,3x,f7.3,3x,i8)
!
end subroutine mmc_correlation_write_corr
!
!*******************************************************************************
!
subroutine mmc_correlation_write_pote(ic, je, title, luse_all)
!-
!  Write a specific correlation for chemical coordination
!+
!
use crystal_mod
use chem_mod
use mmc_mod
!
use precision_mod
use prompt_mod
!
implicit none
!
integer         , intent(in) :: ic         ! Correlations number
integer         , intent(in) :: je         ! Energy type
character(len=*), intent(in) :: title      ! Energy name
logical         , intent(in) :: luse_all ! Do not exist at first pair
!
integer :: is, js
logical :: lfirst
real(kind=PREC_DP) :: divisor
!
cond_ener: IF(mmc_cor_energy (ic, je) ) THEN 
   lfirst = .TRUE.
   loop_type: do is = 0, cr_nscat 
      do js = is, cr_nscat 
         if(mmc_pair(ic,je,is,js) < 0 .AND. lfirst) then
            if(mmc_target_corr (ic, je, is, js) /= 0.0) then
               divisor = ABS(mmc_target_corr (ic, je, is, js))
            else
               divisor = 1.0
            endif
            if     (mmc_pair (ic, je, is, js) /=  0 ) then 
            lfirst = luse_all
            write(output_io, 3100) ic, title      , cr_at_lis (is), cr_at_lis (js),         &
                mmc_target_corr(ic, je, is, js),                                &
                mmc_ini_corr   (ic, je, is, js),                                &
                mmc_ini_sigm (ic, je, is, js),                                     &
                mmc_target_corr(ic, je, is, js) - mmc_ini_corr (ic,je, is, js), &
               (mmc_target_corr(ic, je, is, js) - mmc_ini_corr (ic,je, is, js))/divisor
            endif
         endif
      enddo
   enddo loop_type
endif cond_ener
!
 3100 FORMAT (1x,i3,3x,a9,a5,3x,a5,      8x,5(f7.3,3x)         &
     &             ,i8)                                                 
!
end subroutine mmc_correlation_write_pote
!
!*******************************************************************************
!
subroutine mmc_correlation_write_angl(ic, je, title, luse_all)
!-                                                                       
!     -- Loop over all defined angle correlations                       
!+                                                                       
!
use crystal_mod
use chem_mod
use mmc_mod
!
use precision_mod
use prompt_mod
!
implicit none
!
integer         , intent(in) :: ic         ! Correlations number
integer         , intent(in) :: je         ! Energy type
character(len=*), intent(in) :: title      ! Energy name
logical         , intent(in) :: luse_all ! Do not exist at first pair
!
integer :: k
integer :: iic, kk, iis, jjs, lls
real(kind=PREC_DP) :: divisor
!
cond_ener: IF(mmc_cor_energy (ic, je) ) THEN 
   DO k = 1, mmc_n_angles 
      CALL index2angles (mmc_angles (k), iic, kk, iis, jjs, lls, MAXSCAT)
      IF(mmc_target_angl(k               )/=0.0) THEN
         divisor = mmc_target_angl(k)
      ELSE
         divisor = 1.0
      ENDIF
      WRITE (output_io, 3400) ic, cr_at_lis (iis),  cr_at_lis (jjs), &
         cr_at_lis (lls),  mmc_target_angl (k),  mmc_ach_angl (k),   &
         mmc_ini_sang (k),                                           &
         mmc_target_angl (k)  - mmc_ach_angl (k),                    &
         (mmc_target_angl (k)  - mmc_ach_angl (k))/divisor
   enddo
endif cond_ener
!
 3400 FORMAT (1x,i3,3x,'Angle    ',a5,3x,a5,2x,a5,1x,5(f7.3,3x))
!
end subroutine mmc_correlation_write_angl
!
!*******************************************************************************
!
subroutine calc_change_pid(ic, ie, is, js, damp, divisor, fact, maxdev, change)
!-
!  Calculate depth change and new PID parameters
!+
use mmc_mod
!
use precision_mod
!
implicit none
!
integer,            intent(in)  :: ic
integer,            intent(in)  :: ie
integer,            intent(in)  :: is
integer,            intent(in)  :: js
real(kind=PREC_DP),            intent(in)  :: damp
real(kind=PREC_DP),               intent(in)    :: divisor
real(kind=PREC_DP),               intent(in)    :: fact
real(kind=PREC_DP), dimension(2), intent(inout) :: maxdev
real(kind=PREC_DP), intent(out) :: change
!
real(kind=PREC_DP) :: mmc_pid_summ
real(kind=PREC_DP) :: mmc_pid_dddd
!
mmc_pid_diff(ic, ie, is, js) = mmc_target_corr(ic, ie, is, js) -                &    ! Proportional
                                   mmc_ach_corr   (ic, ie, is, js) 
mmc_pid_inte(ic, ie, is, js, mmc_h_index) =                                     &    ! Integral
                                   mmc_pid_diff(ic, ie, is, js)
mmc_pid_summ = sum(mmc_pid_inte(ic, ie, is, js, :))
mmc_pid_deri(ic, ie, is, js, mmc_h_index) =                                     &    ! Derivative
                                   mmc_ach_corr   (ic, ie, is, js) -            &
                                   mmc_pre_corr   (ic, ie, is, js)
mmc_pid_dddd = sum(mmc_pid_deri(ic, ie, is, js, :))
mmc_pre_corr   (ic, ie, is, js) = mmc_ach_corr   (ic, ie, is, js)
change = -(mmc_pid_pid(1,1) * mmc_pid_diff(ic, ie, is, js) +                    &
           mmc_pid_pid(2,1) * mmc_pid_summ                 +                    &
           mmc_pid_pid(3,1) * mmc_pid_dddd                   ) * damp
!
!
mmc_pid_pid(1,2) = mmc_pid_pid(1,2) + sum(mmc_pid_deri(ic, ie, is, js, :))*fact ! Derivat. PID ==> P
mmc_pid_pid(2,2) = mmc_pid_pid(2,2) + mmc_pid_diff(ic, ie, is, js)        *fact ! Integral PID ==> I
mmc_pid_pid(3,2) = mmc_pid_pid(3,2) + change                              *fact ! Change       ==> D
mmc_pid_pid_n    = mmc_pid_pid_n + 1
maxdev(1) = max(maxdev(1), abs( mmc_target_corr(ic, ie, is, js) -               &
                                mmc_ach_corr   (ic, ie, is, js)))
maxdev(2) = max(maxdev(2), abs((mmc_target_corr(ic, ie, is, js) -               &
                                mmc_ach_corr   (ic, ie, is, js))/divisor))
!
end subroutine calc_change_pid
!
!*******************************************************************************
!
SUBROUTINE mmc_correlations_uni(ic, pneig, rel_cycl, damp, lout, lfeed, MAXSCAT_L, &
                                MAX_COR, maxdev)
!
!     ----- Chemical correlation                                        
!
USE crystal_mod 
USE mc_mod 
USE mmc_mod 
!
USE prompt_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: ic
INTEGER, INTENT(IN) :: MAXSCAT_L
INTEGER, INTENT(IN) :: MAX_COR
INTEGER, DIMENSION(0:MAXSCAT_L, 0:MAXSCAT_L, 1:MAX_COR) , INTENT(IN) :: pneig
REAL(kind=PREC_DP)   , INTENT(IN) :: rel_cycl ! Relative progress along cycles
REAL(kind=PREC_DP)   , INTENT(IN) :: damp
LOGICAL, INTENT(IN) :: lout
LOGICAL, INTENT(IN) :: lfeed
real(kind=PREC_DP), dimension(2), intent(inout) :: maxdev
!
INTEGER :: pair11
INTEGER :: pair12
INTEGER :: pair21
INTEGER :: pair22
INTEGER :: is, js, je
INTEGER :: nneigh
LOGICAL :: lfirst
      REAL(kind=PREC_DP) :: prob11=0.0, prob12, prob22 
REAL(PREC_DP) :: thet
REAL(PREC_DP) :: divisor
REAL(PREC_DP) :: change
REAL(PREC_DP) :: fact
!integer :: iwr, j=0
!
pair11 = 0
pair12 = 0
pair21 = 0
pair22 = 0
thet   = 0.0
DO is = 0, cr_nscat 
!IF(MAxVAL(mmc_pair (ic, MC_UNI, is, 0:cr_nscat))>0 .OR.   &
!   minval(mmc_pair (ic, MC_UNI, is, 0:cr_nscat))<0      ) then
!write(*,*) 'PAIRS is= : ',is,' :: ',pneig(is, 0:cr_nscat, ic), ' ::', ic, ' :: ', mmc_pair (ic, MC_UNI, is, 0:cr_nscat)
!endif
   DO js =  0, cr_nscat 
      IF     (mmc_pair (ic, MC_UNI, is, js) == -1 ) THEN 
         pair12 = pair12 + pneig (is,js, ic)
      ELSEIF (mmc_pair (ic, MC_UNI, is, js) == -2 ) THEN 
         pair21 = pair21 + pneig (is,js, ic)
      ELSEIF (mmc_pair (ic, MC_UNI, is, js) == +1 ) THEN 
         pair11 = pair11 + pneig (is,js, ic)
      ELSEIF (mmc_pair (ic, MC_UNI, is, js) == +2 ) THEN 
         pair22 = pair22 + pneig (is,js, ic)
      ENDIF
   ENDDO
ENDDO
je = MC_UNI 
!                                                                       
nneigh = pair11 + pair12 + pair21 + pair22 
IF (nneigh > 0.) THEN 
   prob11 =  pair11           / REAL(nneigh) 
   prob12 = (pair12 + pair21) / REAL(nneigh) 
   prob22 =  pair22           / REAL(nneigh) 
   thet = 0.5 * (2.0 * pair11 + pair12 + pair21) / REAL(nneigh)                                                     
!           thet = 0.5 * ((pair22 + pair11) + pair12 + pair21) / REAL(nneigh)                                                     
ENDIF 
lfirst = .TRUE.
corr_pair: DO is = 0, cr_nscat 
   DO js = 0 , cr_nscat 
      IF     (mmc_pair (ic, MC_UNI, is, js) /=  0 ) THEN 
        IF (thet /= 0.0.AND.thet /= 1.0) THEN 
           mmc_ach_corr (ic, je, is, js) = (prob11 - thet**2) /&
                                           (thet * (1 - thet) )
           mmc_ach_corr (ic, je, js, is) = (prob11 - thet**2) /&
                                           (thet * (1 - thet) )
        ELSE 
           IF((pair11>0 .OR. pair22>0) .AND. (pair12==0 .AND. pair21==0)) THEN
              mmc_ach_corr (ic, je, is, js) = 1.0 
              mmc_ach_corr (ic, je, js, is) = 1.0 
           ELSEIF((pair12>0 .OR. pair12>0) .AND. (pair11==0 .AND. pair22==0)) THEN
              mmc_ach_corr (ic, je, is, js) =-1.0 
              mmc_ach_corr (ic, je, js, is) =-1.0 
           ELSE 
              mmc_ach_corr (ic, je, is, js) = 0.0 
              mmc_ach_corr (ic, je, js, is) = 0.0 
           ENDIF 
        ENDIF 
!               Feedback mechanism                                      
        IF(mmc_target_corr (ic, MC_UNI, is, js) /= 0.0) THEN
           divisor = ABS(mmc_target_corr (ic, MC_UNI, is, js))
        ELSE
           divisor = 1.0
        ENDIF
        if(mmc_target_corr(ic, MC_UNI, is, js)<0.0) then
          fact = -1.0          ! Adjust shift for PID parameters
        else
          fact =  1.0
        endif
        IF_FEED: IF(lfeed) THEN
           IF_lfeed: IF(mmc_lfeed(ic,MC_UNI)) THEN
        IF_RELC: IF(rel_cycl>0.0) THEN
              call calc_change_pid(ic, MC_UNI, is, js, damp, divisor, fact, maxdev,  change)
           change = -change
!CHANGE    change= mmc_cfac (ic, MC_UNI) * (mmc_target_corr(ic, MC_UNI, is, js)-     &
!CHANGE                                     mmc_ach_corr   (ic, MC_UNI, is, js) )/2. &
!CHANGE           *ABS(mmc_target_corr (ic, MC_UNI, is, js)) * damp
!CHANGE    IF(mmc_target_corr(ic, MC_UNI, is, js)*mmc_ach_corr(ic, MC_UNI, is, js)>=0.0  .AND. &
!CHANGE       mmc_depth(ic, MC_UNI, is, js)*(mmc_depth(ic, MC_UNI, is, js)-change)<0.0         ) THEN
!CHANGE        change = -mmc_depth(ic, MC_UNI, is, js)*0.005
!CHANGE    ENDIF
!          mmc_depth (ic, MC_UNI, is, js) = mmc_depth      (ic, MC_UNI, is, js)    - &
!                  mmc_cfac (ic, MC_UNI) * (mmc_target_corr(ic, MC_UNI, is, js)-     &
!                                           mmc_ach_corr   (ic, MC_UNI, is, js) )/2. &
!                 *ABS(mmc_target_corr (ic, MC_UNI, is, js)) * damp
           mmc_depth(ic, MC_UNI, is, js) = mmc_depth (ic, MC_UNI, is, js) - change
           mmc_depth(ic, MC_UNI, js, is) = mmc_depth (ic, MC_UNI, is, js)
!          mmc_depth(ic, MC_UNI,  0,  0) = mmc_depth (ic, MC_UNI, is, js)
!if(is==1 .and. js==1) then
!if(ic==1) j=j+1
!IF(is==js.AND.is>0) then
!iwr = 90+ic
!write(iwr,'(1(i6), 1x, 7(f12.6,1x),2l2)') j, mmc_depth (ic, MC_UNI, is, js), & 
!mmc_target_corr(ic, MC_UNI, is, js), mmc_ach_corr   (ic, MC_UNI, is, js), &
!mmc_target_corr(ic, MC_UNI, is, js)- mmc_ach_corr   (ic, MC_UNI, is, js), &
! (mmc_target_corr(ic, MC_UNI, is, js)-                               &
!  mmc_ach_corr   (ic, MC_UNI, is, js) ) / 2., damp, -change !          , &
!endif
        ENDIF IF_RELC
        ENDIF IF_lFEED
        ENDIF IF_FEED
!                                                                       
        IF (lout .AND. mmc_pair(ic,MC_UNI,is,js) < 0 .AND. lfirst) THEN
         mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
         mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_corr(ic, je, is, js) - &
                                                mmc_ach_corr   (ic, je, is, js)
            mmc_h_targ(mmc_h_ctarg)           = mmc_target_corr(ic, je, is, js)
           lfirst = .FALSE.
           WRITE (output_io, 3100) ic, cr_at_lis (is), cr_at_lis (js),         &
               mmc_target_corr (ic, je, is, js),                               &
               mmc_ach_corr (ic, je, is, js),                                  &
               mmc_target_corr (ic, je, is, js) - mmc_ach_corr (ic,je, is, js),&
              (mmc_target_corr (ic, je, is, js) - mmc_ach_corr (ic,je, is, js))/divisor,&
               nneigh
!write(*,*) ' DEPTH ', mmc_depth(ic, MC_UNI, 1:2, 1:2), change
        ENDIF 
!                                                                       
      ENDIF 
!
      mmc_ach_corr (ic, je,  0,  0) = mmc_ach_corr (ic, je, is, js)
   ENDDO 
ENDDO corr_pair
!
 3100 FORMAT (1x,i3,3x,'Unidirect',a5,3x,a5,      8x,2(f7.3,3x),        &
     &        10x,f7.3,3x,f7.3,3x,i8)
!
END SUBROUTINE mmc_correlations_uni
!
!*******************************************************************************
!
SUBROUTINE mmc_correlations_group(ic, pneig, rel_cycl, damp, lout, lfeed, MAXSCAT_L, &
                                MAX_COR, maxdev)
!
!     ----- Chemical correlation                                        
!
USE crystal_mod 
USE mc_mod 
USE mmc_mod 
!
USE prompt_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: ic
INTEGER, INTENT(IN) :: MAXSCAT_L
INTEGER, INTENT(IN) :: MAX_COR
INTEGER, DIMENSION(0:MAXSCAT_L, 0:MAXSCAT_L, 1:MAX_COR) , INTENT(IN) :: pneig
REAL(kind=PREC_DP)   , INTENT(IN) :: rel_cycl ! Relative progress along cycles
REAL(kind=PREC_DP)   , INTENT(IN) :: damp
LOGICAL, INTENT(IN) :: lout
LOGICAL, INTENT(IN) :: lfeed
real(kind=PREC_DP), dimension(2), intent(inout) :: maxdev
!
INTEGER :: pair11
INTEGER :: pair12
INTEGER :: pair21
INTEGER :: pair22
INTEGER :: is, js, je
INTEGER :: nneigh
LOGICAL :: lfirst
      REAL(kind=PREC_DP) :: prob11=0.0, prob12, prob22 
REAL(PREC_DP) :: achieved
REAL(PREC_DP) :: target_corr
REAL(PREC_DP) :: depth    
REAL(PREC_DP) :: thet
REAL(PREC_DP) :: divisor
REAL(PREC_DP) :: change
real(kind=PREC_DP) :: fact
!integer :: iwr, j=0
!write(*,*) ' IN GROUP CORRELATIONS ', ic, lout, lfeed
!
pair11 = 0
pair12 = 0
pair21 = 0
pair22 = 0
thet   = 0.0
pair11 = pneig(1, 1, ic)
pair12 = pneig(1, 2, ic)
pair21 = pneig(2, 1, ic)
pair22 = pneig(2, 2, ic)
!
je = MC_GROUP 
!                                                                       
nneigh = pair11 + pair12 + pair21 + pair22 + pneig(1,3,ic) + pneig(2,3,ic) + pneig(3,3,ic) 
!write(*,*) ' LEFT     ', mmc_left( ic, MC_GROUP,:)
!write(*,*) ' RIGHT    ', mmc_right(ic, MC_GROUP,:)
!write(*,*) ' PAIRS 1,:', pneig(1,1:3, ic) 
!write(*,*) ' PAIRS 2,:', pneig(2,1:3, ic) 
!write(*,*) ' PAIRS 3,:', pneig(3,1:3, ic) , sum(pneig)
IF (nneigh > 0.) THEN 
   prob11 =  pair11           / REAL(nneigh) 
   prob12 = (pair12 + pair21) / REAL(nneigh) 
   prob22 =  pair22           / REAL(nneigh) 
   thet = 0.5 * (2.0 * pair11 + pair12 + pair21) / REAL(nneigh)                                                     
!           thet = 0.5 * ((pair22 + pair11) + pair12 + pair21) / REAL(nneigh)                                                     
ENDIF 
IF (thet /= 0.0.AND.thet /= 1.0) THEN 
   achieved = (prob11 - thet**2) /&
              (thet * (1 - thet) )
ELSE 
   IF((pair11>0 .OR. pair22>0) .AND. (pair12==0 .AND. pair21==0)) THEN
      achieved = 1.0 
   ELSEIF((pair12>0 .OR. pair12>0) .AND. (pair11==0 .AND. pair22==0)) THEN
      achieved =-1.0 
   ELSE 
      achieved = 0.0 
   ENDIF 
ENDIF 
lfirst = .TRUE.
is = MAXLOC(mmc_left( ic, MC_GROUP,:), 1) - 1
js = MAXLOC(mmc_right(ic, MC_GROUP,:), 1) - 1
target_corr = mmc_target_corr (ic, MC_GROUP, is, js)
IF(target_corr /= 0.0) then
   divisor = ABS(target_corr)
ELSE
   divisor = 1.0
ENDIF
if(target_corr<0.0) then
  fact = -1.0          ! Adjust shift for PID parameters
else
  fact = +1.0
endif
!write(*,*) ' Achieved ', achieved, prob11, prob12, prob22, thet, target_corr
mmc_ach_corr(ic, MC_GROUP, is, js) = achieved
!write(*,*) ' IS, JS ', is, js, divisor
IF_FEED: IF(lfeed) THEN
   IF_lfeed: IF(mmc_lfeed(ic, MC_GROUP)) THEN
   IF_RELC: IF(rel_cycl>0.0) THEN
      CFAC: IF(mmc_cfac(ic, MC_GROUP)>0.0) THEN
!write(*,*) ' GROUP ', is, js
              call calc_change_pid(ic, MC_GROUP, is, js,  damp, divisor, fact, maxdev,  change)
!     change = mmc_cfac (ic, MC_GROUP) * (target_corr- achieved)/2.             &
!                                      * ABS(target_corr) * damp
!write(*,*) '         CHANGE ', mmc_cfac (ic, MC_GROUP) , (target_corr- achieved)/2.,            &
!                                        ABS(target_corr) , damp, change
       change = -change
!SIG  IF(target_corr*achieved>=0.0  .AND. &
!SIG          mmc_depth(ic, MC_GROUP, is, js)*(mmc_depth(ic, MC_GROUP, is, js)-change)<0.0         ) THEN
!SIG     change = -mmc_depth(ic, MC_GROUP, is, js)*0.005
!SIG  ENDIF
      mmc_depth(ic, MC_GROUP, is, js) = mmc_depth (ic, MC_GROUP, is, js) - change
      mmc_depth(ic, MC_GROUP, js, is) = mmc_depth (ic, MC_GROUP, is, js)
      depth = mmc_depth(ic, MC_GROUP, is, js)
!     mmc_depth_def = ABS(depth)

      DO is = 0, cr_nscat 
         DO js = 0 , cr_nscat 
            IF(mmc_left(ic, MC_GROUP, is) /= 0  .AND. mmc_left(ic, MC_GROUP, js) /= 0 ) THEN
               mmc_depth   (ic, MC_GROUP, is, js) = depth
               mmc_ach_corr(ic, MC_GROUP, is, js) = achieved
            ENDIF
            IF(mmc_right(ic, MC_GROUP, is) /= 0 .AND. mmc_right(ic, MC_GROUP, js) /= 0 ) THEN
               mmc_depth   (ic, MC_GROUP, is, js) = depth
               mmc_ach_corr(ic, MC_GROUP, is, js) = achieved
            ENDIF
            IF(mmc_left(ic, MC_GROUP, is) /= 0  .AND. mmc_right(ic, MC_GROUP, js) /= 0 ) THEN
               mmc_depth   (ic, MC_GROUP, is, js) = depth
               mmc_ach_corr(ic, MC_GROUP, is, js) = achieved
            ENDIF
            IF(mmc_right(ic, MC_GROUP, is) == 0  .AND. mmc_left(ic, MC_GROUP, js) /= 0 ) THEN
               mmc_depth   (ic, MC_GROUP, is, js) = depth
               mmc_ach_corr(ic, MC_GROUP, is, js) = achieved
            ENDIF
         ENDDO
      ENDDO
      ENDIF CFAC
   ENDIF IF_RELC
   ENDIF IF_lfeed
ENDIF IF_FEED
!
IF(lout .AND. lfirst) THEN
   mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
   mmc_h_diff(mmc_h_ctarg, mmc_h_index) = target_corr - achieved
   mmc_h_targ(mmc_h_ctarg)              = target_corr
   lfirst = .FALSE.
   WRITE(output_io, 3100) ic, 'GR1 ', 'GR2 ', target_corr, achieved ,          &
         target_corr - achieved, (target_corr - achieved)/divisor, pneig(1,2,ic) !nneigh
ENDIF
!write(*,*) ' DEPTH , change ', mmc_depth(ic, MC_GROUP, 1:2, 1:2), -change
!     ENDIF 
!do is=0,cr_nscat
!write(*,'(a,12f7.2)') ' Depth  ',mmc_depth(ic, MC_GROUP, is, :)
!enddo
!write(*,*) ' DEPTH ', mmc_depth_def(ic)
!read(*,*) is
!
 3100 FORMAT (1x,i3,3x,'Group cor',a5,3x,a5,      8x,2(f7.3,3x),        &
     &        10x,f7.3,3x,f7.3,3x,i8)
!
END SUBROUTINE mmc_correlations_group
!
!*****7*****************************************************************
!
SUBROUTINE index2angles (ind, ic, nr, is, js, ls, MAXSCAT)
!-                                                                      
!     Calculates the correlation number, and the angle triplet from     
!     a unique number that was determined by angles2index               
!+                                                                      
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: X_ANGLES =  200  ! This needs work, may not be unique!
INTEGER, PARAMETER :: X_SCAT   =   50  ! This needs work, may not be unique!
!
INTEGER, INTENT(IN)  :: ind 
INTEGER, INTENT(OUT) :: ic 
INTEGER, INTENT(OUT) :: nr 
INTEGER, INTENT(OUT) :: is 
INTEGER, INTENT(OUT) :: js 
INTEGER, INTENT(OUT) :: ls 
INTEGER, INTENT(IN)  :: MAXSCAT 
!                                                                       
INTEGER i 
INTEGER na 
!                                                                       
i = ind 
na = MAXSCAT + 2 
na = X_SCAT
!                                                                       
ic = i / (      X_ANGLES * na * na * na) + 1 
i = i - (ic - 1) *       X_ANGLES * na * na * na 
nr = i / (na * na * na) + 1 
i = i - (nr - 1) * na * na * na 
is = i / (na * na) - 1 
i = i - (is + 1) * na * na 
js = i / (na) - 1 
i = i - (js + 1) * na 
ls = i - 2 
!                                                                       
END SUBROUTINE index2angles                   
!
!*****7*****************************************************************
!
end module mmc_basic_mod
