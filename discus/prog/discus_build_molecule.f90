module build_molecule_mod
!-
!  Routines to build up molecules from symmetry information
!
contains
!
!*******************************************************************************
!
subroutine do_build_molecule(cr_natoms, cr_iscat, cr_pos, cr_mole, temp_nmax,   &
           temp_inmole, SPC_MAXX, spc_nn, spc_table_l)
!-
! build the molecule info into cr_mole
!
use molecule_mod
use spcgr_apply , only:get_wyckoff
use wyckoff_mod
!
use lib_functions_mod
use precision_mod
!
implicit none
!
integer                                    , intent(in)    :: cr_natoms
integer           , dimension(3, cr_natoms), intent(in)    :: cr_iscat
real(kind=PREC_DP), dimension(3, cr_natoms), intent(inout) :: cr_pos
integer           , dimension(   cr_natoms), intent(inout) :: cr_mole
integer                                    , intent(in)    :: temp_nmax
integer           , dimension(temp_nmax   ), intent(inout) :: temp_inmole
integer                                        , intent(in) :: SPC_MAXX
integer                                        , intent(in) :: spc_nn  
integer           , dimension(SPC_MAXX, SPC_MAXX), intent(in) :: spc_table_l
!
logical, parameter :: loutput = .false.                ! No output
integer, parameter :: mode = 0                         ! Wyckoff site determination mode
REAL(kind=PREC_DP), PARAMETER      :: EPS = 0.0001
!
integer                              :: i, j, k, l, m  ! Dummy loop indices
integer                              :: temp_n_mole    ! Temporary number of molecules
integer, dimension(:)  , allocatable :: temp_len       ! Temporay molecule lengths
integer, dimension(:)  , allocatable :: temp_type      ! Temporay molecule types
integer, dimension(:,:), allocatable :: temp_mole_cont ! Temporary molecule content
logical, dimension(:)  , allocatable :: temp_is_done   ! Atom is placed
logical               :: lident
!
real(kind=PREC_DP), dimension(3)     :: vec            ! Dummy vectors
real(kind=PREC_DP), dimension(3)     :: trans          ! Dummy vectors
real(kind=PREC_DP), dimension(4)     :: u, orig, copy  ! Dummy vectors
!
!
allocate(temp_len(1:cr_natoms))                        ! We have at most this many molecules
allocate(temp_type(1:cr_natoms))                       ! We have at most this many molecules
allocate(temp_is_done(1:cr_natoms))                    ! We have at most this many molecules
allocate(temp_mole_cont(1:cr_natoms,1:cr_natoms))      ! We have at most this many molecules with this many atoms
temp_len       = 0
temp_type      = 0
temp_mole_cont = 0
temp_is_done   = .false.
temp_n_mole    = 0
!write(*,'(a)') '  i    x      y      z   typ  s   A  M temp'
!do i=1, cr_natoms
!write(*,'(i3, 3f7.3, 3i4, 3i3)') i, cr_pos(:,i), cr_iscat(:,i), cr_mole(i), temp_inmole(i)
!enddo
!write(*,'(a,30i3)' ) ' temp_inmole ', temp_inmole(1:cr_natoms)
!write(*,'(a,30i3)' ) ' cr_mole     ', cr_mole(1:cr_natoms)
!
! 1st loop, collect atoms with positive entries in temp_inmole
!
loop_first: do i= 1, cr_natoms
   if(cr_mole(i)==0) then                              ! Atom does not belong to a molecule
      temp_is_done(i) = .true.
      cycle loop_first
   endif
   if(cr_mole(i)>0 .and. temp_inmole(i)>0) then        ! Atom belongs to a molecule
      if(temp_inmole(i)==1) then                       ! 1st atom, increase number of molecules
         temp_n_mole = temp_n_mole + 1                 ! increase number of molecules
         temp_type(temp_n_mole) = abs(cr_mole(i))
      endif
      temp_len(cr_mole(i)) = temp_len(cr_mole(i)) + 1       ! Incrementz molecule length
      temp_mole_cont(cr_mole(i), temp_len(cr_mole(i))) = i
      temp_is_done(i) = .true.
   endif
enddo loop_first
!write(*,*) 
!write(*,*) 'after loop ONE'
!do i=1, temp_n_mole
!write(*,'(a, i3, a, 30i3)') 'Molecule ', i, ' || ', temp_mole_cont(i,1:temp_len(i))
!enddo
!write(*,'(a,30l3)' ) ' temp_is_done', temp_is_done(1:cr_natoms)
!write(*,'(a,30i3)' ) ' temp_inmole ', temp_inmole(1:cr_natoms)
!write(*,'(a,30i3)' ) ' cr_mole     ', cr_mole(1:cr_natoms)
!
! Second loop, collect atoms in netative molecule numbers
!
loop_negative: do i = 1, cr_natoms                    ! Loop over all atoms
   if(temp_is_done(i)) cycle loop_negative            ! Atom does not belong to a molecule
   if(cr_mole(i)<0 .and. temp_inmole(i)==-1) then     ! First atom of new molecule
      temp_n_mole = temp_n_mole + 1                   ! increase number of molecules
      temp_type(temp_n_mole) = abs(cr_mole(i))        ! Set molecule type
      temp_len(temp_n_mole ) = temp_len(temp_n_mole) + 1       ! Increment molecule length
      temp_mole_cont(temp_n_mole, temp_len(temp_n_mole)) = i  ! Enter atom into molecule
      k = cr_mole(i)                                  ! Keep track of original molecule type
      cr_mole(i) = temp_n_mole
      temp_inmole(i) = 1
      temp_is_done(i) = .true.                              ! Atom  i is done
      loop_rem: do j = i+1, cr_natoms                           ! Loop over remaining atoms
         if(temp_is_done(j))  cycle loop_rem
         if(cr_mole(j)<0 .and. temp_inmole(j)<0) then
            if(k==cr_mole(j)) then
            if(cr_iscat(2,i) == cr_iscat(2,j)) then     ! Same symmetry as 1st atom in molecule
               temp_len(temp_n_mole) = temp_len(temp_n_mole) + 1       ! Increment molecule length
               temp_mole_cont(temp_n_mole, temp_len(temp_n_mole)) = j  ! Enter atom into molecule
               cr_mole(j) = temp_n_mole
               temp_inmole(j) = temp_len(temp_n_mole)
               temp_is_done(j) = .true.                              ! Atom  i is done
            endif
            endif
         endif
      enddo loop_rem
   endif
enddo loop_negative
!write(*,*) 
!write(*,*) 'after loop TWO'
!do i=1, temp_n_mole
!write(*,'(a, i3, a, 30i3)') 'Molecule ', i, ' || ', temp_mole_cont(i,1:temp_len(i))
!enddo
!write(*,'(a,30l3)' ) ' temp_is_done', temp_is_done(1:cr_natoms)
!write(*,'(a,30i3)' ) ' temp_inmole ', temp_inmole(1:cr_natoms)
!write(*,'(a,30i3)' ) ' cr_mole     ', cr_mole(1:cr_natoms)
!
! third loop, find all symmetrically equivalent atoms
!
u    = 0.0_PREC_DP
u(4) = 1.0_PREC_DP
orig(4) = 1.0_PREC_DP
!
loop_scnd: do j=1, temp_n_mole
!write(*,*)
!write(*,'(a,i3, 3f7.3)') ' At molecule ', j, cr_pos(:,temp_mole_cont(j,1))
   vec(1:3) = cr_pos(:,temp_mole_cont(j,1))
   call get_wyckoff(vec, loutput, mode, .false.)
!write(*,'(a, i3, a, 30i3)') 'Molecule Wyckoff ', wyc_n, ' || ', wyc_list(1:wyc_n)
!write(*,'(a, i3, a, 30l3)') 'Molecule Wyckoff ', wyc_n, ' || ', wyc_extra(1:wyc_n)
   loop_inmole: do i=1, ubound(temp_len, 1)
      if(i>temp_len(j)) exit loop_inmole
      u(1:3) = cr_pos(:,temp_mole_cont(j,i))
!write(*,'(a, 2i3, 3f7.3, i3)') ' At Atom ', i, temp_mole_cont(j,i), u(1:3), cr_iscat(1,temp_mole_cont(j,i))
      loop_sym: do k=2,wyc_n                   ! Loop over wyckoff symmetry
         copy   = matmul(wyc_mat(:,:,k), u)
         if(any(abs(u-copy)>EPS)) then         ! This atom is not copied onto itself
!write(*,'(a, 2i3, 3f7.3)') ' At copy ', k, temp_mole_cont(j,i), copy(1:3)
            loop_inner: do l=i+1, cr_natoms    ! Loop over all remaining atoms
!write(*,'(a,i3, 3f7.2, i3)') ' Testing Atom ', l, cr_pos(:,l), cr_iscat(1,l)
               if(cr_iscat(1,temp_mole_cont(j,i))/=cr_iscat(1,l)) cycle loop_inner   ! Different atom type 
               if(temp_is_done(l)) cycle loop_inner     ! Atom is done alrady
!write(*,'(a,i3, 3f7.2, i3)') ' Testing Atom ', l, cr_pos(:,l), cr_iscat(1,l)
               orig(1:3) = cr_pos(:,l)
!
               lident = .true.
               trans = 0.0_PREC_DP
               do m = 1, 3
                  lident = lident.and.  &
                                 (    abs(orig(m) - copy(m) )      .lt.eps .OR.   &
                                 frac(abs(orig(m) - copy(m) )    ) .lt.eps)
                  trans(m) = orig(m) - copy(m)
               enddo
               if(lident) then           ! Atom l is a symmetrically equivalent atom to i
!write(*,'(a, i3, 3(a,3f7.3), l3, i3)') ' At imag ', l, ' || ', orig(1:3), ' || ', copy(1:3), ' || ', trans(1:3), lident, cr_iscat(1,l)
                  temp_len(j) = temp_len(j) + 1
                  temp_mole_cont(j, temp_len(j)) = l      ! Insert atom l into molecule
                  temp_inmole(l)  = temp_len(j)
                  temp_is_done(l) = .true.
                  cr_pos(:,l) = cr_pos(:,l) - trans(1:3)
                  cr_mole(l)  = j
                  temp_inmole(l) = temp_len(j)
               endif
            enddo loop_inner
         endif
      enddo loop_sym
   enddo loop_inmole
enddo loop_scnd
!write(*,*) 
!write(*,*) 'after loop THREE'
!do i=1, temp_n_mole
!write(*,'(a, i3, a, 30i3)') 'Molecule ', i, ' || ', temp_mole_cont(i,1:temp_len(i))
!enddo
!write(*,'(a,30l3)' ) ' temp_is_done', temp_is_done(1:cr_natoms)
!write(*,'(a,30i3)' ) ' temp_inmole ', temp_inmole(1:cr_natoms)
!write(*,'(a,30i3)' ) ' cr_mole     ', cr_mole(1:cr_natoms)
!
mole_num_mole = temp_n_mole
mole_off = 0
mole_cont(1:temp_len(1)) = temp_mole_cont(1,1:temp_len(1))
do i=1, mole_num_mole
   mole_type(i) = temp_type(i)
   mole_char(i) = mole_char(mole_type(i))
   mole_file(i) = mole_file(mole_type(i))
   mole_dens(i) = mole_dens(mole_type(i))
   mole_fuzzy(i) = mole_fuzzy(mole_type(i))

   mole_len(i) = temp_len(i)
enddo
do i=2,mole_num_mole
   mole_off(i) = mole_off(i-1) + mole_len(i-1)
   mole_cont(mole_off(i)+1:mole_off(i)+mole_len(i)) = temp_mole_cont(i,1:temp_len(i))
enddo
!
! Shift molecule with first atom in interval [0,1]
!
do i=1, mole_num_mole
   trans(1) = frac(cr_pos(1,mole_cont(mole_off(i)+1))+1.0D0) - cr_pos(1,mole_cont(mole_off(i)+1))
   trans(2) = frac(cr_pos(2,mole_cont(mole_off(i)+1))+1.0D0) - cr_pos(2,mole_cont(mole_off(i)+1))
   trans(3) = frac(cr_pos(3,mole_cont(mole_off(i)+1))+1.0D0) - cr_pos(3,mole_cont(mole_off(i)+1))
   trans = real(nint(trans), kind=PREC_DP)
   do j=1, mole_len(i)
      k = mole_cont(mole_off(i)+j)
      cr_pos(:,k) = cr_pos(:,k) + trans
   enddo
enddo
!
deallocate(temp_len)
deallocate(temp_type)
deallocate(temp_is_done)
deallocate(temp_mole_cont)
!
end subroutine do_build_molecule
!
!*******************************************************************************
!
end module build_molecule_mod
