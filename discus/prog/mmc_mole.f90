module mmc_mole
!-
!  Multiple monte Carlo procedues for molecules
!+
implicit none
!
!*******************************************************************************
!
contains
!
subroutine mmc_run_multi_mole(lout_feed)
!-
!  Main MMC routine for molecules
!+
!
use chem_mod
use mc_mod
use mmc_mod
use mmc_basic_mod
use symm_sup_mod
USE atom_env_mod
!
!$ use omp_lib
use parallel_mod
!
use precision_mod
use prompt_mod
!
implicit none
!
logical, intent(in) :: lout_feed   ! Print feedback output to screen true/false
LOGICAL             :: lfeed       ! Use Feedback algorithm
!
integer :: iacc_good
integer :: iacc_neut
integer :: iacc_bad
integer(KIND=PREC_INT_LARGE) :: itry
integer(KIND=PREC_INT_LARGE) :: igen
integer :: tid
integer :: nthreads
integer(KIND=PREC_INT_LARGE) :: imodulus
logical                      :: lserial    ! serial calculation if TRUE
LOGICAL                      :: done
logical                      :: lfinished
logical                      :: lout
logical                      :: loop
LOGICAL, DIMENSION(3)        :: old_chem_period     ! Original periodic boundary conditions
REAL(kind=PREC_DP)                         :: rel_cycl
real(kind=PREC_DP), dimension(2) :: maxdev = (/ 0.0, 0.0/)
!
! Privare variables for parallel computing
!
integer           , dimension(2)             :: isel         ! selected molecule numbers
integer                                      :: nmoles       ! number of selected molecules 1/2
integer                                      :: is_move      ! number of selected molecules 1/2
INTEGER           , DIMENSION(:,:)    , allocatable :: iatom
REAL(kind=PREC_DP), DIMENSION(:, :, :), allocatable :: patom
INTEGER           , DIMENSION(  :)    , allocatable :: natom
LOGICAL           , DIMENSION(:,:)    , allocatable :: tatom
real(kind=PREC_DP), dimension(:,:,:), allocatable   :: orig_pos  ! original atom positions
real(kind=PREC_DP), dimension(3)                    :: disp  ! Displacement vector
real(kind=PREC_DP), dimension(3, 3)                 :: matr  ! Rotation matrix
INTEGER                                      :: ncent
logical                                      :: laccept      ! Accept move T/F
real(kind=PREC_DP), dimension(0:MC_N_ENERGY) :: e_old        ! old energies
real(kind=PREC_DP), dimension(0:MC_N_ENERGY) :: e_new        ! new energies
!
ALLOCATE(patom(3, 0:MAX_ATOM_ENV, MMC_MAX_CENT))
ALLOCATE(iatom(   0:MAX_ATOM_ENV, MMC_MAX_CENT))
ALLOCATE(tatom(   0:MAX_ATOM_ENV, MMC_MAX_CENT))
ALLOCATE(natom(                   MMC_MAX_CENT))
!
call symm_setup      ! Do once with default, to avoid allocation that may happen in parallel
!
call mmc_initial(old_chem_period, itry, igen, iacc_good, iacc_neut, iacc_bad, &
           done, loop)
lserial = .FALSE.
!
!iacc_good = 0
!iacc_neut = 0
!iacc_bad  = 0
!itry = 0
nthreads  = 1
imodulus  = max(1_PREC_INT_LARGE, mo_feed/nthreads)
!done = .false.
rel_cycl = 0.0
lfeed = .false.
!                                                                       
!     Initialize the different energies                                 
!                                                                       
lout      = .FALSE.
lfinished = .FALSE.
lfeed     = .FALSE.
CALL mmc_correlations (lout, 0.0D0, done, lfinished, lfeed, maxdev)
!
tid = 0
nthreads = 1
IF(.NOT.lserial .AND. par_omp_use) THEN
!$OMP PARALLEL PRIVATE(tid)
!$   tid = OMP_GET_THREAD_NUM()
!$   IF (tid == 0) THEN
!$      IF(par_omp_maxthreads == -1) THEN
!$         nthreads = OMP_GET_NUM_THREADS()
!$      ELSE
!$         nthreads = MAX(1,MIN(par_omp_maxthreads, OMP_GET_NUM_THREADS()))
!$      ENDIF
!$   END IF
!$OMP END PARALLEL
ENDIF
!
call symm_store                                    ! save symmetry settings
!
cond_parallel: if(nthreads > 1) then
   !$OMP PARALLEL PRIVATE(tid, isel, nmoles, is_move, &
   !$OMP          iatom, patom, tatom, natom,ncent, &
   !$OMP          orig_pos, disp, matr, &
   !$OMP          laccept, e_old, e_new)   &
   !$OMP          SHARED(done)
   !$   tid = OMP_GET_THREAD_NUM()
   !$OMP DO SCHEDULE(DYNAMIC, mo_cyc/nthreads/32)
!
   parallel_loop: DO itry=1, mo_cyc     ! Do mmc in parallel
      IF(done) CYCLE parallel_loop    ! Quickly cycle to end if an error occured
      call mmc_run_loop_mole(tid, nthreads, &
           igen, itry, &
           isel, nmoles, is_move, &
           iatom, patom, natom, tatom, ncent, &
           orig_pos, disp, matr, &
           laccept, e_old, e_new, &
           iacc_good, iacc_neut, iacc_bad, &
           rel_cycl, lout_feed, lfeed, imodulus, done, &
           MAX_ATOM_ENV, MMC_MAX_CENT, MMC_MAX_ATOM)
   enddo parallel_loop
   !$OMP END DO NOWAIT
   !$ IF(tid==0) igen = igen*nthreads
   !$OMP END PARALLEL
   itry = igen
else cond_parallel
   tid = 0
   serial_loop: do itry=1, mo_cyc
      call mmc_run_loop_mole(tid, nthreads, &
           igen, itry, &
           isel, nmoles, is_move, &
           iatom, patom, natom, tatom, ncent, &
           orig_pos, disp, matr, &
           laccept, e_old, e_new, &
           iacc_good, iacc_neut, iacc_bad, &
           rel_cycl, lout_feed, lfeed, imodulus, done, &
           MAX_ATOM_ENV, MMC_MAX_CENT, MMC_MAX_ATOM)
      if(done) exit serial_loop
   end do serial_loop
endif cond_parallel
!
!                                                                       
!------ Loop finished                                                   
!                                                                       
IF(lout_feed) THEN
   WRITE (output_io, 3000)
   WRITE (output_io, 2000) itry, iacc_good, iacc_neut, iacc_bad
ENDIF
lfinished = .TRUE.
lfeed     = .FALSE.   ! no feedback algorithm
CALL mmc_correlations (lout_feed, rel_cycl, done, lfinished, lfeed, maxdev)
!
call symm_restore                                  ! restore symmetry settings
chem_period = old_chem_period                      ! Restore boundary conditions
deALLOCATE(patom)
deALLOCATE(iatom)
deALLOCATE(tatom)
deALLOCATE(natom)
!
 2000 FORMAT (/,' try: ',I10,' acc: (g/n/b): ',I8,        &
     &          ' / ',I8,' / ',I8,'  MC moves ')
 3000 FORMAT (/,' --- Final multiple energy configuration ---')
!
end subroutine mmc_run_multi_mole
!
!*******************************************************************************
!
subroutine mmc_run_loop_mole(tid, nthreads, &
           igen, itry, &
           isel, nmoles, is_move, &
           iatom, patom, natom, tatom, ncent, &
           orig_pos, disp, matr, &
           laccept, e_old, e_new, &
           iacc_good, iacc_neut, iacc_bad, &
           rel_cycl, lout_feed, lfeed, imodulus, done, &
           MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
!-
!  The main running loop for molecule MMC
!+
use chem_mod
use mc_mod
use mmc_mod
use mmc_basic_mod
!
use errlist_mod
use precision_mod
use prompt_mod
!
implicit none
!
INTEGER                           , INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER                           , INTENT(IN) :: MMC_MAX_CENT_L
INTEGER                           , INTENT(IN) :: MMC_MAX_ATOM_L
integer                           , intent(in)    :: tid      ! Thread identifier
integer                           , intent(in)    :: nthreads ! number of threads
integer(KIND=PREC_INT_LARGE)      , intent(inout) :: igen
integer(KIND=PREC_INT_LARGE)      , intent(IN   ) :: itry
integer           , dimension(2)            , intent(inout) :: isel         ! selected molecule numbers
integer                                     , intent(inout) :: nmoles       ! number of selected molecules 1/2
INTEGER, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: iatom
REAL(kind=PREC_DP)   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: patom
INTEGER, DIMENSION(     MMC_MAX_CENT_L)                 , INTENT(INOUT) :: natom
LOGICAL, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: tatom
INTEGER                                                 , INTENT(INOUT) :: ncent
real(kind=PREC_DP), dimension(:,:,:), allocatable , intent(inout) :: orig_pos  ! original atom positions
real(kind=PREC_DP), dimension(3)   , intent(inout) :: disp  ! Displacement vector
real(kind=PREC_DP), dimension(3, 3), intent(inout) :: matr  ! Rotation matrix
logical                                     , intent(inout) :: laccept      ! Accept move T/F
real(kind=PREC_DP), dimension(0:MC_N_ENERGY), intent(inout) :: e_old        ! old energies
real(kind=PREC_DP), dimension(0:MC_N_ENERGY), intent(inout) :: e_new        ! new energies
integer                           , intent(INOUT) :: iacc_good
integer                           , intent(INOUT) :: iacc_neut
integer                           , intent(INOUT) :: iacc_bad
REAL(kind=PREC_DP)                              , INTENT(INOUT) :: rel_cycl
LOGICAL                           , INTENT(IN )   :: lout_feed
LOGICAL                           , INTENT(IN )   :: lfeed       ! Use feedback algorithm  T/F
INTEGER(KIND=PREC_INT_LARGE)      , INTENT(IN)    :: imodulus
LOGICAL                           , INTENT(INOUT) :: done
!
integer                                      :: is_move      ! Type of move displacement or switch
!
real(kind=PREC_DP), dimension(2) :: maxdev = (/ 0.0, 0.0/)
!
if(done) return          ! Quickly cycle to end if an error occured
if(tid==0) then
   igen = igen + 1
endif
call mmc_select_mole(tid, nthreads, isel, nmoles, is_move)                                   ! Select move shift/rotate
!write(*,*) ' SELECTED MOLE ', isel, 'NM ',nmoles, 'MOVE ',is_move
!write(*,*) ' ENERGIES ', mmc_cor_energy ( :, MC_LENNARD)
!write(*,*) ' mmc_cor_energy ', lbound(mmc_cor_energy), ubound(mmc_cor_energy), &
!' | ', mmc_cor_energy(1, MC_LENNARD)
call mmc_energies_mole(isel, nmoles, iatom, patom, natom, tatom, ncent, &
     chem_ncor, CHEM_MAX_COR, MC_N_ENERGY, mmc_cor_energy, is_move, e_old, &
         MAX_ATOM_ENV_L, MMC_MAX_CENT_L)  ! calculate old energies
!
!   !$OMP CRITICAL
call mmc_modify_mole(tid, nthreads, isel, nmoles, is_move, orig_pos, disp, matr)             ! Modify molecule(s)
!   !$OMP end critical
call mmc_energies_mole(isel, nmoles, iatom, patom, natom, tatom, ncent, &
     chem_ncor, CHEM_MAX_COR,  MC_N_ENERGY, mmc_cor_energy, is_move, e_new, &
         MAX_ATOM_ENV_L, MMC_MAX_CENT_L)  ! calculate new energies
!
call mmc_test_multi (iacc_good, iacc_neut, iacc_bad, e_new, e_old, laccept)                   ! Do we have a good move ?
if(.not. laccept) then
   call mmc_unmodify_mole(tid, isel, nmoles, orig_pos)
endif
!
!$OMP CRITICAL
IF(tid==0) THEN
   if(mod(igen, imodulus)==0) then
      IF(lout_feed) WRITE(output_io, 2000) itry*nthreads, &
                iacc_good, iacc_neut, iacc_bad
!
      rel_cycl = REAL(igen)/REAL(mo_cyc)*REAL(NTHREADS)
      CALL mmc_correlations(lout_feed, rel_cycl, done, .FALSE., lfeed, maxdev)

   endif
!
   IF (igen > mo_cyc/nthreads) THEN
      done = .TRUE.
   ENDIF
endif
!$OMP end critical
!
2000 FORMAT (/,' Gen: ',I10,' acc: (g/n/b): ',I8,        &
     &          ' / ',I8,' / ',I8,'  MC moves ')
!
end subroutine mmc_run_loop_mole
!
!*******************************************************************************
!
subroutine mmc_select_mole(tid, nthreads, isel, nmoles, is_move)
!-
!  Select a molecule / a pair of molecules
!+
use mc_mod
use mmc_mod
use molecule_mod
!
use precision_mod
use random_mod
!
implicit none
!
integer              , intent(in ) :: tid           ! Thread ID
integer              , intent(in ) :: nthreads      ! number threads
integer, dimension(2), intent(out) :: isel          ! selected molecule numbers
integer              , intent(out) :: nmoles        ! number of selected molecules 1/2
integer              , intent(out) :: is_move       ! Type of move displacement or switch
!
integer, parameter :: MAXTRY = 10000                ! Maximum number of searches for a molecule
!
integer :: i                                        ! dummy loop index
integer :: ntry                                     ! Number of trials so far
real(kind=PREC_DP) :: r1                            ! a random number
!
ntry = 0
loop_search: do                                          ! Search until a valid molecule is found
   ntry = ntry + 1
   if(ntry > MAXTRY) exit loop_search
!
   call random_number(r1)
   i = 1
   do while (r1 >= mmc_move_cprob(i) .and. i <  MC_N_MOVE)
      i = i + 1                                         ! Increment until move type is found
   end do
   is_move = i
   isel    = 0                                          ! No molecules selected
!
   cond_moves: if (is_move == MC_MOVE_DISP .or.                                  &
                   is_move == MC_MOVE_ROTATE    ) then  ! Choose proper moves type
!
      call random_number(r1)
      isel(1) = int(r1*mole_num_mole/nthreads) + tid*mole_num_mole/nthreads + 1               ! Choose random molecule
      nmoles  = 1
   end if cond_moves                                    ! end Choose proper moves type
!
end do loop_search                                      ! end Search until a valid molecule is found
!
end subroutine mmc_select_mole
!
!*******************************************************************************
!
subroutine mmc_modify_mole(tid, nthreads, isel, nmoles, is_move, orig_pos,      &
           disp, matr)
!-
!  Modifys the molecule(s) by a move or rotation
!+
use crystal_mod
use mmc_mod
use molecule_mod
!use symm_mod
!use symm_sup_mod
!
use metric_mod
!
!use symm_menu
!
use lib_random_func
!use matrix_mod
!
implicit none
!
integer              , intent(in) :: tid             ! thread ID
integer              , intent(in) :: nthreads        ! number of threads
integer, dimension(2), intent(in) :: isel            ! selected molecule numbers
integer              , intent(in) :: nmoles          ! number of selected molecules 1/2
integer              , intent(in) :: is_move         ! Type of move displacement or switch
real(kind=PREC_DP), dimension(:,:,:), allocatable, intent(OUT) :: orig_pos  ! original atom positions
real(kind=PREC_DP), dimension(3)   , intent(inout) :: disp  ! Displacement vector
real(kind=PREC_DP), dimension(3, 3), intent(inout) :: matr  ! Rotation matrix
!
integer :: i
!real(kind=PREC_DP), dimension(3, 3) :: invm  ! Rotation matrix
real(kind=PREC_DP), dimension(3)    :: u                     ! Dummy vector 
real(kind=PREC_DP)                  :: sigma                 ! Dummy vector 
!
allocate(orig_pos(3, mole_len(isel(1)), 0:nthreads-1))
!
do i=1,mole_len(isel(1))
  orig_pos(:,i,tid) = cr_pos(:,mole_cont(mole_off(isel(1))+i))
! write(*,*) ' BLEN PRIOR ', do_blen(.TRUE., cr_pos(:,mole_cont(mole_off(isel(1))+i)) , &
!                                            cr_pos(:,mole_cont(mole_off(isel(1))+1)))
enddo
!
cond_moves: if(is_move == MC_MOVE_DISP) then
   do i=1,3
     disp(i) = gasdev(DBLE(mo_maxmove_mole(i, mole_type(isel(1)) )))
   enddo
   do i=1, mole_len(isel(1))                            ! Loop over all atoms in molecule
      cr_pos(:,mole_cont(mole_off(isel(1))+i)) = &
      cr_pos(:,mole_cont(mole_off(isel(1))+i)) + disp
! write(*,*) ' BLEN SHIFT ', do_blen(.TRUE., cr_pos(:,mole_cont(mole_off(isel(1))+i)) , &
!                                            cr_pos(:,mole_cont(mole_off(isel(1))+1)))
   enddo
!
elseif (is_move == MC_MOVE_ROTATE) then cond_moves      ! Rotate molecule
!
   sigma = DBLE(mo_maxrota_mole(1, mole_type(isel(1)) ))
!write(*,*) ' ISEL  ', isel(1), mole_type(isel(1))
!write(*,*) ' SIGMA ', sigma, DBLE(mo_maxrota_mole(1, mole_type(isel(1)) ))
   call symm_setup_mmc(matr, sigma)
!  call matinv3(matr, invm)
!  write(*,'(a,5f17.10)') 'MATR ',matr(1,:), sqrt(matr(1,1)**2+matr(1,2)**2+matr(1,3)**2), sym_angle
!  write(*,'(a,4f17.10)') 'MATR ',matr(2,:), sqrt(matr(2,1)**2+matr(2,2)**2+matr(2,3)**2)
!  write(*,'(a,4f17.10)') 'MATR ',matr(3,:), sqrt(matr(3,1)**2+matr(3,2)**2+matr(3,3)**2)
!  write(*,'(a,3f17.10)') 'RTEN ',cr_rten(1,:)
!  write(*,'(a,3f17.10)') 'RTEN ',cr_rten(2,:)
!  write(*,'(a,3f17.10)') 'RTEN ',cr_rten(3,:)
!  do i=1,3
!  write(*,'(a,3f17.10)') 'EPS  ',cr_eps (1,:,i)
!  write(*,'(a,3f17.10)') 'EPS  ',cr_eps (2,:,i)
!  write(*,'(a,3f17.10)') 'EPS  ',cr_eps (3,:,i)
!  enddo
!  do i=1,3
!  write(*,'(a,3f17.10)') 'REPS ',cr_reps(1,:,i)
!  write(*,'(a,3f17.10)') 'REPS ',cr_reps(2,:,i)
!  write(*,'(a,3f17.10)') 'REPS ',cr_reps(3,:,i)
!  enddo
!  write(*,*) ' Determinant rot ', det3(matr), det3(invm)
   do i=1, mole_len(isel(1))                            ! Loop over all atoms in molecule
      u = cr_pos(:,mole_cont(mole_off(isel(1))+i)) - cr_pos(:,mole_cont(mole_off(isel(1))+1)) ! subtract origin
      cr_pos(:,mole_cont(mole_off(isel(1))+i)) = matmul(matr, u) + cr_pos(:,mole_cont(mole_off(isel(1))+1))
! write(*,*) ' BLEN ROTO  ', do_blen(.TRUE., cr_pos(:,mole_cont(mole_off(isel(1))+i)) , &
!                                            cr_pos(:,mole_cont(mole_off(isel(1))+1)))
   enddo
end if cond_moves
!
end subroutine mmc_modify_mole
!
!*******************************************************************************
!
subroutine mmc_unmodify_mole(tid, isel, nmoles, orig_pos)
!-
!  Modifys the molecule(s) by a move or rotation
!+
use crystal_mod
use mmc_mod
use molecule_mod
!
use metric_mod
!
use lib_random_func
!
implicit none
!
integer              , intent(in) :: tid             ! thread ID
integer, dimension(2), intent(in) :: isel            ! selected molecule numbers
integer              , intent(in) :: nmoles          ! number of selected molecules 1/2
real(kind=PREC_DP), dimension(:,:,:), allocatable, intent(INOUT) :: orig_pos  ! original atom positions
!
integer :: i
!
do i=1,mole_len(isel(1))
   cr_pos(:,mole_cont(mole_off(isel(1))+i)) = orig_pos(:,i, tid)   ! Restore original position
! write(*,*) ' BLEN RESTO ', do_blen(.TRUE., cr_pos(:,mole_cont(mole_off(isel(1))+i)) , &
!                                            cr_pos(:,mole_cont(mole_off(isel(1))+1)))
enddo
!read(*,*) i
!
deallocate(orig_pos)
!
end subroutine mmc_unmodify_mole
!
!*******************************************************************************
!
subroutine mmc_energies_mole(isel, nmoles, iatom, patom, natom, tatom, ncent, &
     n_corr, MAXCOR, N_ENERGY, mmc_cor_energy, is_move, e_cur, &
         MAX_ATOM_ENV_L, MMC_MAX_CENT_L)
!     call mmc_energies_mole(isel, nmoles, chem_ncor, MC_N_ENERGY, mmc_cor_energy, is_move, e_old)  ! calculate old energies
!-
!  Calculate energy state for current molecule(s)
!
use mc_mod
use precision_mod
!
implicit none
!
INTEGER                           , INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER                           , INTENT(IN) :: MMC_MAX_CENT_L
integer, dimension(2), intent(in ) :: isel            ! selected molecule numbers
integer              , intent(in ) :: nmoles          ! number of selected molecules 1/2
INTEGER           , DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: iatom
REAL(kind=PREC_DP), DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: patom
INTEGER           , DIMENSION(     MMC_MAX_CENT_L)                 , INTENT(INOUT) :: natom
LOGICAL           , DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: tatom
INTEGER                                                            , INTENT(INOUT) :: ncent
integer              , intent(in ) :: n_corr          ! Number of correlations defined by user
integer              , intent(in ) :: MAXCOR          ! MAXIMUM number correlations
integer              , intent(in ) :: N_ENERGY        ! Maximum number energies
logical, dimension(0:MAXCOR, 0:N_ENERGY), intent(in) :: mmc_cor_energy 
integer              , intent(in ) :: is_move         ! Type of move displacement or switch
real(kind=PREC_DP)   , dimension(0:N_ENERGY), intent(OUT) :: e_cur        ! current energies
!
integer :: im                                         ! Counter molecules
integer :: ic                                         ! Counter correlations/targets
!
e_cur = 0.0
!
loop_moles: do im = 1, nmoles                           ! Loop over all molecules
   loop_target: do ic = 1, n_corr                       ! Loop over all targets
!write(*,*) ' ENERGY ? ', im, isel(im), ic, mmc_cor_energy (ic, MC_LENNARD)
      if(mmc_cor_energy (ic, MC_LENNARD)) then          ! Target is Lennard-Jones
         e_cur(MC_LENNARD) = e_cur(MC_LENNARD) + mmc_energy_len_mole(isel(im), ic, &
                             iatom, patom, natom, tatom, ncent, &
         MAX_ATOM_ENV_L, MMC_MAX_CENT_L)
      endif
   enddo loop_target                                    ! END Loop over all targets
enddo loop_moles                                        ! END Loop over all molecules
!
end subroutine mmc_energies_mole
!
!*******************************************************************************
!
function mmc_energy_len_mole(im, ic, iatom, patom, natom, tatom, ncent, &
         MAX_ATOM_ENV_L, MMC_MAX_CENT_L) result(e_len)
!-
!  Calculate the lenard Jones energy at molecule im
!+
use crystal_mod
use chem_neig_multi_mod
use metric_mod
use mmc_mod
use molecule_mod
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP)  :: e_len
!
INTEGER                           , INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER                           , INTENT(IN) :: MMC_MAX_CENT_L
integer, intent(in) :: im                               ! Molecule number
integer, intent(in) :: ic                               ! Target   number
INTEGER           , DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: iatom
REAL(kind=PREC_DP), DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: patom
INTEGER           , DIMENSION(     MMC_MAX_CENT_L)                 , INTENT(INOUT) :: natom
LOGICAL           , DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: tatom
INTEGER                                                            , INTENT(INOUT) :: ncent
!
!integer, parameter ::  MAXW = 20    ! Maximum array size
!integer, parameter ::  MAX_CENT = 1    ! Maximum array size
!                                                                       
integer ::  jatom   ! Central atom no, get neighbours around jatom
!integer, dimension(   0:MAXW, MAX_CENT) :: iatom ! indices of neighbs
!real(kind=PREC_DP)   , dimension(3, 0:maxw, MAX_CENT) :: patom ! Coordinates
!logical, dimension(   0:maxw, MAX_CENT) :: tatom ! indices of neighbs
!integer, dimension(           MAX_CENT) :: natom ! no of neigh
!integer                                 :: ncent ! no of central atoms
integer :: ia                                    ! Loop index molecule content
integer :: j                                     ! Loop index neighbors
integer :: is, js                                ! Atom types 
!
real(kind=PREC_DP)               :: d
real(kind=PREC_DP), dimension(3) :: u
real(kind=PREC_DP), dimension(3) :: v
!              
!character(len=40) :: string
!do ia=1,2
!is = ia
!js = ia
!if(mmc_pair(ic, MC_LENNARD,is,js)/=0) then
!  write(string,'(a8,i2.2)') 'lennard.',ia
!open(8,file=string      , status='unknown')
!write(*,*) ' LENN ', ic, ia, mmc_len_a(ic, is, js), mmc_len_b(ic, is, js), &
!mmc_len_m(ic, is, js), mmc_len_n(ic, is, js)
!do j=1, 200
!  d=2.0 + j*0.05
!  write(8,'(f6.3,2x, g18.8e3)') d,  &
!  (mmc_len_a(ic, is, js)/d**mmc_len_m(ic, is, js) - mmc_len_b(ic, is, js)/d**mmc_len_n(ic, is, js))
!enddo
!if(mmc_pair(ic, MC_LENNARD,is,js)/=0) then
!close(8)
!endif
!enddo
!
e_len = 0.0
!
loop_atoms: do ia = 1, mole_len(im)                     ! Loop over all atoms in the molecule
   jatom = mole_cont(mole_off(im)+ia)                    ! Actual atom number i crystal
   is = cr_iscat(jatom)                                  ! Atom type at start atom
   call chem_neighbour_multi(jatom, ic, iatom, patom, tatom, natom, &
                                 ncent, MAX_ATOM_ENV_L, MMC_MAX_CENT_L)
   u = cr_pos(:, jatom)
   loop_neig: do j = 1, natom(1)
      js = cr_iscat(iatom(j,1))
      if(mmc_pair(ic, MC_LENNARD,is,js)/=0) then
         v = patom(:, j, 1)                                 ! patom contains offset
         d = do_blen (.TRUE., u, v)
         e_len = e_len + (mmc_len_a(ic, is, js)/d**mmc_len_m(ic, is, js) -  &
                          mmc_len_b(ic, is, js)/d**mmc_len_n(ic, is, js))
      endif
   end do loop_neig
end do loop_atoms
!write(*,*) ' ENERGY LEN ', e_len
!
end function mmc_energy_len_mole
!
!*******************************************************************************
!
subroutine symm_setup_mmc(matr, sigma)
!-
! Set up a symmetry matrix locally only
!+
use crystal_mod
!
use precision_mod
use lib_random_func
USE trig_degree_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3,3), intent(inout) :: matr
real(kind=PREC_DP)                , intent(in   ) :: sigma
!
integer :: i, j, k, l
integer                            :: sym_use 
integer                            :: sym_orig_type
integer                            :: sym_axis_type
integer                            :: sym_power
logical                            :: sym_orig_mol 
real(kind=PREC_DP)                 :: sym_angle
real(kind=PREC_DP), dimension(3)   :: sym_orig
real(kind=PREC_DP), dimension(3)   :: sym_uvw 
real(kind=PREC_DP), dimension(3)   :: sym_hkl 
real(kind=PREC_DP), dimension(4,4) :: sym_mat
!real(kind=PREC_DP), dimension(4,4) :: sym_rmat
!
REAL(kind=PREC_DP) :: uij
REAL(kind=PREC_DP) :: ctheta, stheta
REAL(kind=PREC_DP) :: length
REAL(kind=PREC_DP) :: sym_d (3), sym_r (3)
REAL(kind=PREC_DP) :: kron (3, 3)
DATA kron / 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0 /
!
sym_use = 0                                          ! Use explicit matrix
sym_orig      = 0.0D0                                ! Axis through 1st mole atom
sym_orig_type = 0                                    ! Axis through 1st mole atom
sym_orig_mol  = .true.
do i=1,3
   sym_uvw(i) = 2.0 * ran1(0) - 1.0                  ! Make random rotation axis
enddo
sym_hkl = matmul(real(cr_gten,KIND=PREC_DP), sym_uvw)
sym_angle = gasdev(sigma)                            ! random angle
sym_axis_type = 0                                    ! Symmetry axis defined by uvw
sym_power     = 1                                    ! Apply once
!
sym_mat  = 0.0d0
sym_mat (4, 4) = 1.0d0
ctheta = cosd (sym_angle)
stheta = sind (sym_angle)
!
length = sqrt(dot_product(sym_uvw, matmul(real(cr_gten, kind=PREC_DP), sym_uvw)))
sym_d = sym_uvw/length
length = sqrt(dot_product(sym_hkl, matmul(real(cr_rten, kind=PREC_DP), sym_hkl)))
sym_r = sym_hkl/length
!
DO i = 1, 3
   DO j = 1, 3
      sym_mat (i, j) = 0.0d0
      DO k = 1, 3
         DO l = 1, 3
            sym_mat(i, j) = sym_mat(i, j) + real(cr_rten(i, k), kind=PREC_DP) * &
                                            real(cr_eps(k, l,j),kind=PREC_DP) * sym_d(l)
         ENDDO
      ENDDO
      uij = sym_d (i) * sym_r (j)
      sym_mat(i, j) = sym_mat(i, j) * stheta + uij + (kron(i, j) - uij) * ctheta
   ENDDO
ENDDO
!
!write(*,*) ' sigma ', sigma
!write(*,*) ' angle ', sym_angle,  ctheta, stheta
!write(*,*) ' uvw   ', sym_uvw, sym_d
!write(*,*) ' hkl   ', sym_hkl, sym_r
!write(*,*) ' Matr  ', sym_mat(1,:)
!write(*,*) ' Matr  ', sym_mat(2,:)
!write(*,*) ' Matr  ', sym_mat(3,:)
matr = sym_mat(1:3, 1:3)
!
end subroutine symm_setup_mmc
!
!*******************************************************************************
!
end module mmc_mole
