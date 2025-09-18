module build_molecule_mod
!-
!  Routines to build up molecules from symmetry information
!
contains
!
!*******************************************************************************
!
subroutine do_build_molecule(cr_natoms, cr_iscat, cr_pos, cr_mole, temp_nmax,   &
           temp_inmole, SPC_MAX, spc_table)
!-
! build the molecule info into cr_mole
!+
!
use metric_mod
use molecule_mod
use spcgr_apply, only:get_wyckoff, first_mole
!
use errlist_mod
use param_mod
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
integer                                        , intent(in) :: SPC_MAX
integer           , dimension(SPC_MAX, SPC_MAX), intent(in) :: spc_table
!integer                                    , intent(inout)  :: mole_num_mole ! Number of molecules
!
logical, parameter :: lspace  = .true.
logical, parameter :: loutput = .false.
integer, parameter :: mode = 0
!
integer :: i, j, k, l   ! Dummy indices
!integer, dimension(1:cr_natoms) :: mole_len     ! Number of atoms in molecule(:)
logical, dimension(1:cr_natoms) :: tested       ! This atom needs to be tested )
integer, dimension(0:48) :: wyckoff_sym
integer, dimension(:), allocatable :: temp_len     ! Temporay molecule lengths
real(kind=PREC_DP), dimension(3) :: vec   !temporary vector
real(kind=PREC_DP), dimension(3) :: u     !temporary vector
real(kind=PREC_DP), dimension(3) :: v     !temporary vector
real(kind=PREC_DP)               :: dist  ! Distance to centrasl atom
!
tested = .true.
mole_len = 0
mole_num_mole = maxval(temp_inmole)       ! Determine current number of molecules
mole_num_mole = maxval(cr_mole(1:cr_natoms))       ! Determine current number of molecules
loop_init: do i=1, cr_natoms              ! Determine current molecule lengths
   if(cr_mole(i)>0 .and. temp_inmole(i)>0) then   ! Atom belongs to a molecule  
      mole_len(cr_mole(i)) = max(mole_len(cr_mole(i)), temp_inmole(i))
   endif
!write(*,'(i3, 3f6.3, 3i4, 3i3)') i, cr_pos(:,i), cr_iscat(:,i), cr_mole(i), temp_inmole(i)
enddo loop_init
!write(*,*) '#############################################################################'
!write(*,*) ' MOLE LEN ', mole_len(1:mole_num_mole)
!write(*,*) ' ATOMS    ', cr_natoms
!
loop_main: do i=1, cr_natoms 
!write(*,*) '??? atom ', i, cr_iscat(1,i), cr_mole(i), temp_inmole(i), cr_mole(i)>0, temp_inmole(i) == 1, mole_num_mole
   cond_mole: if(cr_mole(i)>0) then       ! Atom is in molecule
      if(temp_inmole(i) == 1) then        ! This is the first atom, determine wyckoff symmetry
         tested(i) = .false.
         vec=cr_pos(:,i)
         wyckoff_sym = 0
         call get_wyckoff(vec, loutput, mode)
         wyckoff_sym(0) = nint(res_para(2))
         wyckoff_sym(1:wyckoff_sym(0)) = nint(res_para(4:3+wyckoff_sym(0)))
!write(*,'(a,2i4, 16i4)') 'Atom ', i, cr_iscat(1,i), wyckoff_sym(1:wyckoff_sym(0))
         u = cr_pos(:,i)
!
         loop_remain: do j=i+1, cr_natoms              ! Search for remaining atoms in this molecule
            if(tested(j) .and. cr_mole(j)==cr_mole(i) .and. temp_inmole(j)>1) then
               v = cr_pos(:,j)
               dist = do_blen (lspace, u, v)
!write(*,*) ' FOUND FURTHER ATOM   ', j, cr_iscat(:,j), cr_mole(i), dist
               tested(j) = .false.        ! No need for further testing
               loop_inner: do k=j+1, cr_natoms
                  if(cr_iscat(2,k)==cr_iscat(2,i)) cycle loop_remain   ! This is a new Wyckoff group
                  if(tested(k) .and. cr_mole(k)==-cr_mole(i) .and. &
                     cr_iscat(1,j)==cr_iscat(1,k) .and. temp_inmole(k)<-1) then
                  do l=1,wyckoff_sym(0)
                     v = cr_pos(:,k)
                     if(cr_iscat(2,k)==wyckoff_sym(l)) then   ! Atoms belongs to Wyckoff group
!write(*,*) ' Testing  ', k, cr_iscat(1,k), cr_iscat(2,k), l, wyckoff_sym(l), mole_len(cr_mole(i)), &
!do_blen (lspace, u, v), abs(dist-do_blen (lspace, u, v))>0.01
                        mole_len(cr_mole(i)) = mole_len(cr_mole(i)) + 1
                        cr_mole(k) = cr_mole(i)
                        temp_inmole(k) = mole_len(cr_mole(i))
                        tested(k) = .false.
                        if(abs(dist-do_blen (lspace, u, v))>0.01) then    ! Atom needs a shift
                           call do_shift(u, v, dist, ier_num)
                           if(ier_num /= 0) then
                              ier_typ = ER_APPL
                              return
                           endif
                           cr_pos(:,k) = v
                        endif
                        cycle loop_inner
                     endif
                  enddo
                  endif
               enddo loop_inner
            endif
         enddo loop_remain
      endif
!  else cond_mole                         ! Atom is not inside a molecule
!     cycle loop_main
   endif cond_mole
enddo loop_main
!write(*,*) ' SECOND LOOP '
!
!  Second loop, search remainig atoms inside a molecule
!
loop_secnd: do i=1, cr_natoms
!write(*,*) '??? atom ', i, cr_iscat(1,i), cr_mole(i), temp_inmole(i), cr_mole(i)<0, temp_inmole(i) == 1, mole_num_mole
   cond_new: if(cr_mole(i)<0) then       ! Atom is in molecule, make new molecule
      mole_num_mole = mole_num_mole + 1
      mole_len(mole_num_mole) = 1
      j = -cr_mole(i)
      mole_type (mole_num_mole) = mole_type (j)
      mole_char (mole_num_mole) = mole_char (j)
      mole_file (mole_num_mole) = mole_file (j)
      mole_dens (mole_num_mole) = mole_dens (j)
      mole_fuzzy(mole_num_mole) = mole_fuzzy(j)
      if(temp_inmole(i) == -1) then      ! This is the first atom, determine wyckoff symmetry
         tested(i) = .false.
         vec=cr_pos(:,i)
         wyckoff_sym = 0
         call get_wyckoff(vec, loutput, mode)
         wyckoff_sym(0) = nint(res_para(2))
         wyckoff_sym(1:wyckoff_sym(0)) = nint(res_para(4:3+wyckoff_sym(0)))
         u = cr_pos(:,i)
!write(*,'(a, 20i4)') 'New atom ', i, cr_iscat(1,i), cr_mole(i), temp_inmole(i), wyckoff_sym(1:wyckoff_sym(0)), mole_num_mole
         loop_remain2: do j=i+1, cr_natoms              ! Search for remaining atoms in this molecule
            if(tested(j) .and. cr_mole(j)==cr_mole(i) .and. temp_inmole(j)<-1) then
               v = cr_pos(:,j)
               dist = do_blen (lspace, u, v)
!write(*,*) ' FOUND FURTHER ATOM   ', j, cr_iscat(1,i), cr_mole(i)
               mole_len(mole_num_mole) = mole_len(mole_num_mole) + 1
               tested(j) = .false.        ! No need for further testing
               cr_mole(j) = mole_num_mole
               temp_inmole(j) = -temp_inmole(j)
!write(*,*) ' found FURTHER ATOM   ', j, cr_iscat(1,i), cr_mole(i), temp_inmole(j), mole_num_mole, mole_len(mole_num_mole)
               loop_inner2: do k=j+1, cr_natoms
                  if(cr_iscat(2,k)==cr_iscat(2,i)) cycle loop_remain2   ! This is a new Wyckoff group
!                 if(tested(k) .and. cr_mole(k)== cr_mole(i) .and. temp_inmole(k)<-1) then
                  if(tested(k) .and. cr_mole(k)== cr_mole(i) .and. &
                     cr_iscat(1,j)==cr_iscat(1,k) .and. temp_inmole(k)<-1) then
                     v = cr_pos(:,k)
                  do l=1,wyckoff_sym(0)
                     if(cr_iscat(2,k)==spc_table(wyckoff_sym(l),cr_iscat(2,i))) then
!write(*,*) ' Testing  ', k, cr_iscat(2,k), l, wyckoff_sym(l), mole_len(abs(cr_mole(i))), &
!do_blen (lspace, u, v), abs(dist-do_blen (lspace, u, v))>0.01
                        mole_len(mole_num_mole) = mole_len(mole_num_mole) + 1
                        cr_mole(k) = mole_num_mole
                        temp_inmole(k) = mole_len(mole_num_mole)
                        tested(k) = .false.
                        if(abs(dist-do_blen (lspace, u, v))>0.01) then    ! Atom needs a shift
                           call do_shift(u, v, dist, ier_num)
                           if(ier_num /= 0) then
                              ier_typ = ER_APPL
                              return
                           endif
                           cr_pos(:,k) = v
                        endif
                        cycle loop_inner2
                     endif
                  enddo
                  endif
               enddo loop_inner2
            endif
         enddo loop_remain2
      endif
      cr_mole(i) = mole_num_mole
      temp_inmole(i) = -temp_inmole(i)
   endif cond_new
enddo loop_secnd
!write(*,*) 'FINISHED SECOND ', ier_num, ier_typ
!
allocate(temp_len(mole_num_mole))
temp_len = 0
!write(*,*) ' MOLE_LEN  ', mole_len(1:mole_num_mole)
mole_off = 0
do i=1, mole_num_mole
   mole_off(i) = mole_off(i-1) + mole_len(i-1)
enddo
!write(*,*) ' MOLE_OFF  ', mole_off(1:mole_num_mole)
do i=1, cr_natoms
   if(cr_mole(i)>0) then    ! This atom is in a molecule
      temp_len(cr_mole(i)) = temp_len(cr_mole(i)) + 1
      mole_cont(mole_off(cr_mole(i))+temp_len(cr_mole(i))) = i
   endif
enddo
!do l= mole_num_mole,1, -1
!j = l
!mole_num_curr = l
!write(*,*) ' AT MOLE ', l
!call first_mole(j)
!write(*,*) ' ERROR ?', ier_num, ier_typ
!enddo
!do i=1, mole_num_mole
!write(*,'(a,i4,a,20i4)') ' Mole Nr: ',i, ' : ', mole_cont(mole_off(i)+1:mole_off(i)+mole_len(i))
!enddo
!
end subroutine do_build_molecule
!
!*******************************************************************************
!
subroutine do_shift(central, neighbor, distance, ier_num)
!
!     Moves atoms by +- one unit cell to keep molecules concatenated    
!     The algorithm of this subroutine is based on the assumption,      
!     that the first atom is on the highest point of symmetry of the    
!     molecule, i.e. is not copied by any of the symmetry operations    
!     or generators that make up the molecule symmetry. This has the    
!     desired effect that the bond distance of all symmetry             
!     equivalent atoms to this first atom is constant and can be taken  
!     as a reference. If a symmetry operation of the space group        
!     moves an atom out of the current unit cell, one can move this     
!     atom back by integer unit cell vectors until the bond distance    
!     is correct again.                                                 
!     If the molecule does not include an atom on the highest point     
!     of the molecule symmetry, you must insert a "void" on this site.
!-
!
use metric_mod
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3), intent(in)    :: central   ! Central position
real(kind=PREC_DP), dimension(3), intent(inout) :: neighbor  ! Original neighbor position
real(kind=PREC_DP)              , intent(in)    :: distance  ! Intended distance
integer                         , intent(out)   :: ier_num
!
logical, parameter :: lspace = .true.
real(kind=PREC_DP), parameter :: TOL =  0.01 ! Tolerance
!
integer :: k1, k2, k3  ! Dummy indices
real(kind=PREC_DP), dimension(3) :: vec  ! temporary vector
real(kind=PREC_DP), dimension(3) :: v_min! Best solution
real(kind=PREC_DP)               :: d_min! temporary vector
real(kind=PREC_DP)               :: dd   ! temporary vector
!
ier_num = 0
d_min = 0.1
v_min = neighbor   ! Default to original position
loop_1: do k1 = 2, - 2, - 1
   vec(1) = neighbor(1) + real(k1, kind=PREC_DP)
   do k2 = 2, - 2, - 1
      vec(2) = neighbor(2) + real(k2, kind=PREC_DP)
      do k3 = 2, - 2, - 1
         vec(3) = neighbor(3) + real(k3, kind=PREC_DP)
         dd = do_blen(lspace, central, vec)
         if(abs(distance-dd)<d_min) then  ! Found a better solution
            d_min = abs(distance-dd)
            v_min = vec
            if(d_min<TOL) exit loop_1
         endif
      enddo
   enddo
enddo loop_1
!
if(d_min<TOL) then
   neighbor = v_min
else
!write(*,'(a, 3f6.3, 2x, 3f6.3, 2x, 2f6.3)') ' ERROR ', central , neighbor, distance, d_min
   ier_num = -83
endif
!
end subroutine do_shift
!
end module build_molecule_mod
