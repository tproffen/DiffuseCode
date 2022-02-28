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
integer :: nthreads
integer(KIND=PREC_INT_LARGE) :: imodulus
LOGICAL                      :: done
logical                      :: lfinished
logical                      :: lout
logical                      :: loop
LOGICAL, DIMENSION(3)        :: old_chem_period     ! Original periodic boundary conditions
REAL(kind=PREC_DP)                         :: rel_cycl
real(kind=PREC_DP), dimension(2) :: maxdev = (/ 0.0, 0.0/)
!
call mmc_initial(old_chem_period, itry, igen, iacc_good, iacc_neut, iacc_bad, &
           done, loop)
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
!
call symm_store                                    ! save symmetry settings
serial_loop: do itry=1, mo_cyc
   call mmc_run_loop_mole(itry, nthreads, iacc_good, iacc_neut, iacc_bad, &
        rel_cycl, lout_feed, lfeed, imodulus, done)
   if(done) exit serial_loop
end do serial_loop
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
!
 2000 FORMAT (/,' try: ',I10,' acc: (g/n/b): ',I8,        &
     &          ' / ',I8,' / ',I8,'  MC moves ')
 3000 FORMAT (/,' --- Final multiple energy configuration ---')
!
end subroutine mmc_run_multi_mole
!
!*******************************************************************************
!
subroutine mmc_run_loop_mole(itry, nthreads, iacc_good, iacc_neut, iacc_bad, &
           rel_cycl, lout_feed, lfeed, imodulus, done)
!-
!  The main running loop for molecule MMC
!+
use chem_mod
use mc_mod
use mmc_mod
use mmc_basic_mod
!
use precision_mod
use prompt_mod
!
implicit none
!
integer(KIND=PREC_INT_LARGE)      , intent(IN   ) :: itry
integer                           , intent(INOUT) :: nthreads
integer                           , intent(INOUT) :: iacc_good
integer                           , intent(INOUT) :: iacc_neut
integer                           , intent(INOUT) :: iacc_bad
REAL(kind=PREC_DP)                              , INTENT(INOUT) :: rel_cycl
LOGICAL                           , INTENT(IN )   :: lout_feed
LOGICAL                           , INTENT(IN )   :: lfeed       ! Use feedback algorithm  T/F
INTEGER(KIND=PREC_INT_LARGE)      , INTENT(IN)    :: imodulus
LOGICAL                           , INTENT(INOUT) :: done
!
integer           , dimension(2)             :: isel         ! selected molecule numbers
integer                                      :: nmoles       ! number of selected molecules 1/2
integer                                      :: is_move      ! Type of move displacement or switch
logical                                      :: laccept      ! Accept move T/F
real(kind=PREC_DP), dimension(0:MC_N_ENERGY) :: e_old        ! old energies
real(kind=PREC_DP), dimension(0:MC_N_ENERGY) :: e_new        ! new energies
!
real(kind=PREC_DP), dimension(:,:), allocatable :: orig_pos  ! original atom positions
real(kind=PREC_DP), dimension(2) :: maxdev = (/ 0.0, 0.0/)
!
call mmc_select_mole(isel, nmoles, is_move)                                                   ! Select move shift/rotate
!write(*,*) ' SELECTED MOLE ', isel, 'NM ',nmoles, 'MOVE ',is_move
!write(*,*) ' ENERGIES ', mmc_cor_energy ( :, MC_LENNARD)
!write(*,*) ' mmc_cor_energy ', lbound(mmc_cor_energy), ubound(mmc_cor_energy), &
!' | ', mmc_cor_energy(1, MC_LENNARD)
call mmc_energies_mole(isel, nmoles, chem_ncor, CHEM_MAX_COR, MC_N_ENERGY, mmc_cor_energy, is_move, e_old)  ! calculate old energies
call mmc_modify_mole(isel, nmoles, is_move, orig_pos)                                         ! Modify molecule(s)
call mmc_energies_mole(isel, nmoles, chem_ncor, CHEM_MAX_COR,  MC_N_ENERGY, mmc_cor_energy, is_move, e_new)  ! calculate new energies
call mmc_test_multi (iacc_good, iacc_neut, iacc_bad, e_new, e_old, laccept)                   ! Do we have a good move ?
if(.not. laccept) then
   call mmc_unmodify_mole(isel, nmoles, orig_pos)
endif
!
if(mod(itry, imodulus)==0) then
   IF(lout_feed) WRITE(output_io, 2000) itry*nthreads, &
             iacc_good, iacc_neut, iacc_bad
!
   rel_cycl = REAL(itry)/REAL(mo_cyc)*REAL(NTHREADS)
   CALL mmc_correlations(lout_feed, rel_cycl, done, .FALSE., lfeed, maxdev)

endif
!
IF (itry > mo_cyc/nthreads) THEN
   done = .TRUE.
ENDIF
!
2000 FORMAT (/,' Gen: ',I10,' acc: (g/n/b): ',I8,        &
     &          ' / ',I8,' / ',I8,'  MC moves ')
!
end subroutine mmc_run_loop_mole
!
!*******************************************************************************
!
subroutine mmc_select_mole(isel, nmoles, is_move)
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
      isel(1) = int(r1*mole_num_mole) + 1               ! Choose random molecule
      nmoles  = 1
   end if cond_moves                                    ! end Choose proper moves type
!
end do loop_search                                      ! end Search until a valid molecule is found
!
end subroutine mmc_select_mole
!
!*******************************************************************************
!
subroutine mmc_modify_mole(isel, nmoles, is_move, orig_pos)
!-
!  Modifys the molecule(s) by a move or rotation
!+
use crystal_mod
use mmc_mod
use molecule_mod
use symm_mod
use symm_sup_mod
!
use metric_mod
!
use symm_menu
!
use lib_random_func
use matrix_mod
!
implicit none
!
integer, dimension(2), intent(in) :: isel            ! selected molecule numbers
integer              , intent(in) :: nmoles          ! number of selected molecules 1/2
integer              , intent(in) :: is_move         ! Type of move displacement or switch
real(kind=PREC_DP), dimension(:, :), allocatable, intent(OUT) :: orig_pos  ! original atom positions
!
integer :: i
real(kind=PREC_DP), dimension(3)    :: disp  ! Displacement vector
real(kind=PREC_DP), dimension(3, 3) :: matr  ! Rotation matrix
!real(kind=PREC_DP), dimension(3, 3) :: invm  ! Rotation matrix
real(kind=PREC_DP), dimension(3)    :: u                     ! Dummy vector 
!
allocate(orig_pos(3, mole_len(isel(1))))
!
do i=1,mole_len(isel(1))
  orig_pos(:,i) = cr_pos(:,mole_cont(mole_off(isel(1))+i))
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
   sym_use = 0                                          ! Use explicit matrix
   sym_orig      = 0.0D0                                ! Axis through 1st mole atom
   sym_orig_type = 0                                    ! Axis through 1st mole atom
   sym_orig_mol  = .true.
   do i=1,3
      sym_uvw(i) = 2.0 * ran1(0) - 1.0                  ! Make random rotation axis
   enddo
   sym_hkl = matmul(real(cr_gten,KIND=PREC_DP), sym_uvw)
   sym_angle = gasdev(DBLE(mo_maxrota_mole(1, mole_type(isel(1)) )))  ! random angle
   sym_axis_type = 0                                    ! Symmetry axis defined by uvw
   sym_power     = 1                                    ! Apply once
   call symm_setup
   matr = sym_mat(1:3, 1:3)
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
subroutine mmc_unmodify_mole(isel, nmoles, orig_pos)
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
integer, dimension(2), intent(in) :: isel            ! selected molecule numbers
integer              , intent(in) :: nmoles          ! number of selected molecules 1/2
real(kind=PREC_DP), dimension(:, :), allocatable, intent(INOUT) :: orig_pos  ! original atom positions
!
integer :: i
!
do i=1,mole_len(isel(1))
   cr_pos(:,mole_cont(mole_off(isel(1))+i)) = orig_pos(:,i)   ! Restore original position
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
subroutine mmc_energies_mole(isel, nmoles, n_corr, MAXCOR, N_ENERGY, mmc_cor_energy, is_move, e_cur)
!     call mmc_energies_mole(isel, nmoles, chem_ncor, MC_N_ENERGY, mmc_cor_energy, is_move, e_old)  ! calculate old energies
!-
!  Calculate energy state for current molecule(s)
!
use mc_mod
use precision_mod
!
implicit none
!
integer, dimension(2), intent(out) :: isel            ! selected molecule numbers
integer              , intent(out) :: nmoles          ! number of selected molecules 1/2
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
         e_cur(MC_LENNARD) = e_cur(MC_LENNARD) + mmc_energy_len_mole(isel(im), ic)
      endif
   enddo loop_target                                    ! END Loop over all targets
enddo loop_moles                                        ! END Loop over all molecules
!
end subroutine mmc_energies_mole
!
!*******************************************************************************
!
function mmc_energy_len_mole(im, ic) result(e_len)
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
integer, intent(in) :: im                               ! Molecule number
integer, intent(in) :: ic                               ! Target   number
!
integer, parameter ::  MAXW = 20    ! Maximum array size
integer, parameter ::  MAX_CENT = 1    ! Maximum array size
!                                                                       
integer ::  jatom   ! Central atom no, get neighbours around jatom
integer, dimension(   0:MAXW, MAX_CENT) :: iatom ! indices of neighbs
real(kind=PREC_DP)   , dimension(3, 0:maxw, MAX_CENT) :: patom ! Coordinates
logical, dimension(   0:maxw, MAX_CENT) :: tatom ! indices of neighbs
integer, dimension(           MAX_CENT) :: natom ! no of neigh
integer                                 :: ncent ! no of central atoms
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
                                 ncent, maxw, MAX_CENT)
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
end module mmc_mole
