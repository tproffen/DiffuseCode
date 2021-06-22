MODULE do_find_mod
!
!*****7*****************************************************************
!
integer, DIMENSION(:)  , ALLOCATABLE :: grand     ! Initial long list of neighbors
integer, DIMENSION(:,:), ALLOCATABLE :: grand_o   ! Initial long list of neighbors, incl. offsets
integer                              :: grand_n   ! Number of neighbors found
!
private 
public  do_find_env
public  do_find_mol
public  do_find_nei_vect_many
public  do_find_neig_conn
public  do_find_neig_conn_all
!
contains
!
!*****7*****************************************************************
!
SUBROUTINE do_find_env (ianz, werte, maxw, x, rmin, rmax, fq, fp) 
!                                                                       
!     This routine finds all atoms around x with a minimal              
!     distance of rmin and a maximum distance of rmax.                  
!-                                                                      
USE discus_config_mod 
USE crystal_mod 
USE atom_env_mod 
USE check_bound_mod
USE check_blen_mod
USE celltoindex_mod
USE modify_func_mod
USE param_mod 
USE errlist_mod 
USE sorting_mod 
USE precision_mod
!                                                                       
IMPLICIT none 

!                                                                       
integer                            , intent(in) :: ianz
integer                            , intent(in) :: maxw 
REAL(KIND=PREC_DP), dimension(MAXW), intent(in) :: werte ! (maxw) 
REAL              , dimension(3)   , intent(in)  :: x     ! (3) 
REAL                               , intent(in)  :: rmin
REAL                               , intent(in)  :: rmax 
LOGICAL                            , intent(in)  :: fq
LOGICAL           , dimension(3)   , intent(in)  :: fp 
!                                                                       
integer :: i, j, k, ii 
integer :: ix, iy, iz
integer :: ix1, ix2, iy1,iy2, iz1,iz2
integer :: istart (3), iend (3), iii (3), cell (3), iatom 
REAL    :: offset (3), nooffset (3) 
LOGICAL :: ltype 
integer, DIMENSION(  :), ALLOCATABLE :: tmp_ind
integer, DIMENSION(  :), ALLOCATABLE :: tmp_env
REAL   , DIMENSION(:,:), ALLOCATABLE :: tmp_pos
REAL   , DIMENSION(  :), ALLOCATABLE :: tmp_dis
!                                                                       
DATA nooffset / 0.0, 0.0, 0.0 / 
!                                                                       
atom_env (0) = 0 
res_para (0) = 0 
!
IF(fq) then 
!                                                                       
!------ --quick version only looks at neighbouring unit cells           
!                                                                       
   DO i = 1, 3 
!     iii(i)    = int(x(i) - cr_dim0(i, 1) ) + 1 
      iii(i)    = int(x(i) - cr_dim(i, 1) ) + 1 
      istart(i) = iii(i) - 1 - int(rmax / cr_a0(i) )  -1
      iend(i)   = iii(i) + 1 + int(rmax / cr_a0(i) )  + 1
   ENDDO 
!                                                                       
   DO k = istart (3), iend (3) 
      DO j = istart (2), iend (2) 
         DO i = istart (1), iend (1) 
            cell (1) = i 
            cell (2) = j 
            cell (3) = k 
!                                                                       
            CALL check_bound (cell, offset, fp, ltype) 
            IF(ltype) then 
               DO ii = 1, cr_ncatoms 
                  CALL celltoindex (cell, ii, iatom) 
                    ltype = atom_allowed (iatom, werte, ianz, maxw) .and. &
                            check_select_status(iatom, .true., cr_prop (iatom),  cr_sel_prop)                    
                  IF (ltype) then 
                     CALL check_blen (x, iatom, rmin, rmax, offset) 
                     IF (ier_num.ne.0) return 
                  ENDIF 
               ENDDO 
            ENDIF 
!                                                                       
         ENDDO 
      ENDDO 
   ENDDO 
!                                                                       
!     --exact: loop over all atoms in crystal that have been added      
!                                                                       
   DO i = cr_ncatoms*cr_icc(1)*cr_icc(2)*cr_icc(3) + 1, cr_natoms
!                                                                       
      IF(fp(1) .or. fp(2) .or. fp(3)) then 
         ier_num = - 16 
         ier_typ = ER_CHEM 
         ier_msg(1) = 'Number of atoms in crystal is larger than'
         ier_msg(2) = 'atoms_per_unit_cell*number_of_unit_cells'
         ier_msg(3) = 'Use set crystal, 1, 1, 1,n[1] to correct'
         RETURN 
      ENDIF 
!                                                                       
      ltype = atom_allowed(i, werte, ianz, maxw)                &
      .and.check_select_status (i, .true., cr_prop (i), cr_sel_prop)    
      IF (ltype) then 
         CALL check_blen (x, i, rmin, rmax, nooffset) 
         IF (ier_num /= 0) return 
      ENDIF 
   ENDDO 
ELSE 
!                                                                       
!     --exact: loop over all atoms in crystal                           
!                                                                       
!!!      IF (fp (1) .or.fp (2) .or.fp (3) ) then 
!!!         ier_num = - 16 
!!!         ier_typ = ER_CHEM 
!!!         ier_msg(1) = 'Chem/exact mode is incompatible with '
!!!!        ier_msg(2) = 'periodic boundary conditions'
!!!         ier_msg(3) = ' '
!!!         RETURN 
!!!      ENDIF 
!                                                                       
!        IF(fp(1)) THEN
!           ix1 = -1
!           ix2 = 1
!        ELSE
!           ix1 = 0
!           ix2 = 0
!        ENDIF
!        IF(fp(2)) THEN
!           iy1 = -1
!           iy2 = 1
!        ELSE
!           iy1 = 0
!           iy2 = 0
!        ENDIF
!        IF(fp(3)) THEN
!           iz1 = -1
!           iz2 = 1
!        ELSE
!           iz1 = 0
!           iz2 = 0
!        ENDIF
   ix1 = 0
   ix2 = 0
   iy1 = 0
   iy2 = 0
   iz1 = 0
   iz2 = 0
   IF(fp(1)) THEN
      IF(x(1)       -cr_dim(1,1) < rmax*1.5) ix1 = -1
      IF(cr_dim(1,2)-x(1)        < rmax*1.5) ix2 =  1
   ENDIF
   IF(fp(2)) THEN
      IF(x(2)       -cr_dim(2,1) < rmax*1.5) iy1 = -1
      IF(cr_dim(2,2)-x(2)        < rmax*1.5) iy2 =  1
   ENDIF
   IF(fp(3)) THEN
      IF(x(3)       -cr_dim(3,1) < rmax*1.5) iz1 = -1
      IF(cr_dim(3,2)-x(3)        < rmax*1.5) iz2 =  1
   ENDIF
   DO i = 1, cr_natoms 
      ltype = atom_allowed (i, werte, ianz, maxw)                    &
             .and.check_select_status (i, .true., cr_prop (i), cr_sel_prop)    
      IF (ltype) then 
         DO ix = ix1, ix2, 1
            offset(1) = ix*cr_icc(1)
            DO iy = iy1, iy2, 1
               offset(2) = iy*cr_icc(2) 
               DO iz = iz1, iz2, 1
                     offset(3) = iz*cr_icc(3)
                  CALL check_blen (x, i, rmin, rmax, offset) 
               ENDDO
            ENDDO
         ENDDO
         IF(ier_num /= 0) RETURN
      ENDIF 
   ENDDO 
ENDIF 
!
!     Sort neighbors according to distance
!
IF(atom_env(0)>0) THEN
   ALLOCATE(tmp_ind(  1:atom_env(0)))
   ALLOCATE(tmp_env(  1:atom_env(0)))
   ALLOCATE(tmp_pos(3,0:atom_env(0)))
   ALLOCATE(tmp_dis(  1:atom_env(0)))
   tmp_env    = atom_env(1:atom_env(0))
   tmp_ind    = 0
   tmp_pos    = atom_pos
   tmp_dis    = atom_dis(1:atom_env(0))
   CALL indexx(atom_env(0),tmp_dis,tmp_ind)
   DO i=1,atom_env(0)
      atom_env(i)   = tmp_env(  tmp_ind(i))
      atom_pos(:,i) = tmp_pos(:,tmp_ind(i))
      atom_dis(i)   = tmp_dis(  tmp_ind(i))
   ENDDO
   DEALLOCATE(tmp_ind)
   DEALLOCATE(tmp_env)
   DEALLOCATE(tmp_pos)
   DEALLOCATE(tmp_dis)
ENDIF
!                                                                       
END SUBROUTINE do_find_env                    
!
!*****7*****************************************************************
!
      SUBROUTINE do_find_mol (ianz, werte, maxw, x, rmin, rmax) 
!                                                                       
!     This routine finds all molecules around x with a minimal          
!     distance of rmin and a maximum distance of rmax.                  
!-                                                                      
!                                                                       
      USE discus_config_mod 
      USE check_blen_mod
      USE crystal_mod 
      USE molecule_mod 
      USE mole_env_mod 
      USE errlist_mod 
      USE param_mod 
USE precision_mod
      IMPLICIT none 
       
!                                                                       
      integer ianz, maxw 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL x (3) 
      REAL rmin, rmax 
!                                                                       
      integer i, j 
!                                                                       
      REAL nooffset (3) 
      LOGICAL ltype 
!                                                                       
      DATA nooffset / 0.0, 0.0, 0.0 / 
!                                                                       
      mole_env (0) = 0 
      res_para (0) = 0 
!                                                                       
!     Exact loop over all molecules, no periodic boundary conditions    
!                                                                       
      DO i = 1, mole_num_mole 
      IF (ianz.eq. - 1) then 
         ltype = .true. 
      ELSE 
         ltype = .false. 
         DO j = 1, ianz 
         ltype = i.eq.nint (werte (j) ) 
         ENDDO 
      ENDIF 
      IF (ltype) then 
         CALL check_blen_mol (x, i, rmin, rmax, nooffset) 
      ENDIF 
      IF (ier_num.ne.0) return 
      ENDDO 
!DBG                                                                    
!DBG      do i=1,mole_env(0)                                            
!DBG        write (output_io,*) 'molecule ',mole_env(i)                 
!DBG      ENDDO                                                         
      END SUBROUTINE do_find_mol                    
!
!*****7*****************************************************************
!
subroutine do_find_nei_vect_many(ianz, werte, maxw, nvect, iatom, vect)
!-
!  Find the neighbor under the given vectors
!+
!
use crystal_mod
use atom_env_mod
use check_bound_mod
use chem_mod
use celltoindex_mod
USE modify_func_mod
!
use precision_mod
!
implicit none
!
integer                            , intent(in) :: ianz
integer                            , intent(in) :: MAXW 
REAL(KIND=PREC_DP), dimension(MAXW), intent(in) :: werte ! (maxw) 
integer                     , intent(in) :: nvect    ! Number of vectors
integer                     , intent(in) :: iatom    ! The central==start atom
integer, dimension(5, nvect), intent(in) :: vect     ! (isite, jsite, cx, cy, cz)
!
integer, dimension(3) :: icell                ! Unit cell number for start atom
integer, dimension(3) :: jcell                ! Unit cell number for end   atom
integer               :: isite                ! Site number for start atom
integer               :: jsite                ! Site number for end   atom
integer               :: jatom                ! Target atom number
integer               :: i                    ! Dummy index
logical               :: lok                  ! Flag if boundary is OK
logical               :: ltype                ! Flag atom is OK
real(kind=PREC_SP), dimension(3) :: offset    ! in case of periodic boundary conditions, the offset
!                      
atom_env(0) = 0
!ltype = atom_allowed (iatom, werte, ianz, maxw) .and. &
!        check_select_status(iatom, .true., cr_prop (iatom),  cr_sel_prop)                    
!if(.not. ltype) return                        ! Atom has incorrect type or properties
call indextocell(iatom, icell, isite)         ! Determine unit cell location
!
loop_vect: do i=1, nvect
   if(isite==vect(1,i)) then                       ! Start atom is at correct site
      jcell = icell + vect(3:5,i)                  ! Add unit cell shift to atom
      lok   = .false.
      call check_bound(jcell, offset, chem_period, lok)
      if(lok) then
         jsite = vect(2,i)                         ! Target site
         call celltoindex(jcell, jsite, jatom)     ! Determine unit cell location
         ltype = atom_allowed (jatom, werte, ianz, maxw) .and. &
                 check_select_status(jatom, .true., cr_prop (jatom),  cr_sel_prop)                    
         if(.not. ltype) cycle loop_vect
         atom_env(0)             = atom_env(0) + 1
         atom_env(atom_env(0))   = jatom
         atom_pos(:,atom_env(0)) = cr_pos(:,jatom) + offset
      endif
   endif
enddo loop_vect
!
end subroutine do_find_nei_vect_many
!
!*****7*****************************************************************
!
subroutine do_find_neig_conn(ianz, MAXW, werte) !iatom, iconn)
!-
!  Find the connectivity for atom iatom, use connectivitie stored in
!  Werte (2:), no search beyond these connectivities is performed
!+
!
use crystal_mod
use atom_env_mod
use conn_sup_mod
!
use errlist_mod
use precision_mod
!
implicit none
!
integer, intent(in) :: ianz                       ! Number of parameters
integer, intent(in) :: MAXW                       ! Array size
real(kind=PREC_DP), dimension(MAXW) :: werte      ! Values
integer :: iatom                                  ! Find connectivity for this atom
integer :: nconn                                  ! Central atom has this many connectivities
!
integer                              :: inoo      ! Find this connectivity number
integer                              :: c_natoms  ! Number of atoms connected
integer, DIMENSION(:  ), ALLOCATABLE :: c_list    ! List of atoms connected to current
integer, DIMENSION(:,:), ALLOCATABLE :: c_offs    ! Offsets for atoms connected to current
integer, DIMENSION(:)  , ALLOCATABLE :: list      ! List of connectivities to use
!
integer :: i, j, k ! Dummy index
integer :: itype   ! Atom type for iatom
!
atom_env(0) = 0
!
iatom = nint(werte(1))
itype = cr_iscat(iatom)
nconn = get_connectivity_numbers(itype)           ! Determine number of connectivities
if(nconn == 0) return                             ! No connectivities, nothing to do
allocate(list(1:nconn))                           ! List of connectivities to search
!
if(nint(werte(ianz))==-1) then                    ! Use all connectivities
   do i=1,nconn
      list(i) = i
   enddo
else                                              ! Use explicit list
   k = 0
   do i=2, ianz                                   ! Conn's start in werte(2)
      if(nint(werte(i))<=nconn .and. nint(werte(i))>0) then
         k = k + 1
         list(k) = nint(werte(i))
      else
         ier_num = -177
         ier_typ = ER_APPL
         write(ier_msg(1),'(a,i8)') 'Wrong number ', nint(werte(i))
         return
      endif
   enddo
   nconn = k                                      ! Use corrected length
endif
!
atom_env(0) = 0
!
loop_conns: do j=1, nconn                         ! Loop over all requested connectivities
   inoo  = list(j)
   call get_connectivity_list(iatom, itype, inoo, c_list, c_offs, c_natoms)
   if(ier_num /= 0) then
      atom_env(0) = 0
      if(allocated(list))   deallocate(list)
      if(allocated(c_list)) deallocate(c_list)
      if(allocated(c_offs)) deallocate(c_offs)
      return
   endif
!
!  Copy into atom_env
!
   do i=1, c_natoms
      atom_env(atom_env(0)+i)   = c_list(i)        ! Add atoms to end of list
      atom_pos(:,atom_env(0)+i) = cr_pos(:,c_list(i)) + c_offs(:,(i))
   enddo
   atom_env(0) = atom_env(0) + c_natoms            ! Update length
enddo loop_conns
!
if(allocated(list)) deallocate(list)
if(allocated(c_list)) deallocate(c_list)
if(allocated(c_offs)) deallocate(c_offs)
!
end subroutine do_find_neig_conn
!
!*****7*****************************************************************
!
subroutine do_find_neig_conn_all(ianz, MAXW, werte)
!-
!  Find the connectivity for atom iatom, use connectivities stored in
!  Werte (2:), do search beyond these connectivities as well
!+
!
use crystal_mod
use atom_env_mod
use conn_sup_mod
!
use errlist_mod
use precision_mod
!
implicit none
!
integer, intent(in) :: ianz                       ! Number of parameters
integer, intent(in) :: MAXW                       ! Array size
real(kind=PREC_DP), dimension(MAXW) :: werte      ! Values
integer :: iatom                                  ! Find connectivity for this atom
integer :: nconn                                  ! Central atom has this many connectivities
!
integer                              :: inoo      ! Find this connectivity number
!
integer :: i, j, k ! Dummy index
integer :: itype   ! Atom type for iatom
!
atom_env(0) = 0
!
iatom = nint(werte(1))
itype = cr_iscat(iatom)
nconn = get_connectivity_numbers(itype)           ! Determine number of connectivities
if(nconn == 0) return                             ! No connectivities, nothing to do
!
atom_env(0) = 0
allocate(grand(1:cr_natoms))                      ! Initial (very) long list of neighbors
allocate(grand_o(3, 1:cr_natoms))                 ! Initial (very) long list of neighbors, incl offset
grand   = 0
grand_o = 0
grand_n = 1                                       ! We start with one atom that has been found
grand(iatom) = 1                                  ! This atom is done, conn has been searched
!
!  Now find (recursive) neighbors
!
call do_find_neig_conn_rec(iatom)

j = 0
do i=1, cr_natoms
   if(grand(i)==1) then
      j = j + 1
      atom_env(j)     = i
      atom_pos(:,j) = cr_pos(:,i) + grand_o(:, i)
   endif
enddo
atom_env(0) = j
!
if(allocated(grand))    deallocate(grand)
if(allocated(grand_o))  deallocate(grand_o)
!
end subroutine do_find_neig_conn_all
!
!*****7*****************************************************************
!
recursive subroutine do_find_neig_conn_rec(iatom)
!-
!  Recursively find neighbors
!+
!
use crystal_mod
use conn_sup_mod
!
use errlist_mod
!
integer, intent(in) :: iatom
!
integer                              :: c_natoms  ! Number of atoms connected
integer, DIMENSION(:  ), ALLOCATABLE :: c_list    ! List of atoms connected to current
integer, DIMENSION(:,:), ALLOCATABLE :: c_offs    ! Offsets for atoms connected to current
integer :: itype              ! Central atom type
integer :: nconn              ! Number of connectivities this atom
integer :: ino                ! Current   connectivity   this atom
integer :: j, i, k
!
grand(iatom)= 1                   ! Mark this atom as done
itype = cr_iscat(iatom)
nconn = get_connectivity_numbers(itype)           ! Determine number of connectivities
if(nconn==0) return                               ! No connectivity, we are done
!
loop_conns: do j=1, nconn                         ! Loop over all requested connectivities
   ino = j
   call get_connectivity_list(iatom, itype, ino, c_list, c_offs, c_natoms)
   if(ier_num /= 0) then
      if(allocated(c_list)) deallocate(c_list)
      if(allocated(c_offs)) deallocate(c_offs)
      return
   endif
   do i=1, c_natoms
      k = c_list(i)
      if(k>0) then
         if(grand(k)==0) then                     ! Atom has not been searched
            grand(k)     = -1                     ! Mark as new neighbor
            grand_n      = grand_n + 1            ! New neighbor, increment total number
            grand_o(:,k) = grand_o(:,iatom) + c_offs(:,i)   ! Remember accumulated conn
            call do_find_neig_conn_rec(k)         ! Now search for neighbors of this neighbor
         endif
      endif
   enddo
enddo loop_conns
!
if(allocated(c_list)) deallocate(c_list)
if(allocated(c_offs)) deallocate(c_offs)
!
end subroutine do_find_neig_conn_rec
!
!*****7*****************************************************************
!
END MODULE do_find_mod
