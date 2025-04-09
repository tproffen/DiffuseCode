module dis_estimate_mod
!
! Procedures to estimate number of unit cells, average atoms per unit cell
!
private
!
public estimate_ncells
public estimate_ncatom
public find_average
!
contains
!
!*******************************************************************************
!
subroutine estimate_ncells(ncell_out, pdt_dims, pdt_ilow, pdt_ihig, pdt_ncells)
!-
!  Try to estimate the number of unit cells
!+
use crystal_mod
!
use precision_mod
!
implicit none
!
integer, dimension(3), intent(out) :: ncell_out
real(kind=PREC_DP), dimension(3, 2), intent(inout) :: pdt_dims
integer           , dimension(3)   , intent(out)   :: pdt_ilow        ! Unit cell dimensions in periodic
integer           , dimension(3)   , intent(out)   :: pdt_ihig        ! low and high inidce
integer                            , intent(out)   :: pdt_ncells     ! Number of cells in periodic crystal volume
!
integer :: i          ! Dummy counter
!
pdt_dims(:,1) = cr_pos(:,1)           ! Initialize the min/ max dimensions
pdt_dims(:,2) = cr_pos(:,1)
!
do i=1, cr_natoms
   pdt_dims(1,1) = min(pdt_dims(1,1), cr_pos(1,i))
   pdt_dims(1,2) = max(pdt_dims(1,2), cr_pos(1,i))
   pdt_dims(2,1) = min(pdt_dims(2,1), cr_pos(2,i))
   pdt_dims(2,2) = max(pdt_dims(2,2), cr_pos(2,i))
   pdt_dims(3,1) = min(pdt_dims(3,1), cr_pos(3,i))
   pdt_dims(3,2) = max(pdt_dims(3,2), cr_pos(3,i))
enddo
!
!write(*,'(a, f9.5, 2x, f9.5)') ' dims x ', pdt_dims(1,:)
!write(*,'(a, f9.5, 2x, f9.5)') ' dims y ', pdt_dims(2,:)
!write(*,'(a, f9.5, 2x, f9.5)') ' dims z ', pdt_dims(3,:)
!
pdt_ilow = 0
pdt_ihig = 0
pdt_ncells = 1
do i=1, 3
   if(nint(pdt_dims(i,2)-pdt_dims(i,1))>2) then
!  pdt_ilow(i) = nint(pdt_dims(i,1))! + 1
!  pdt_ihig(i) = nint(pdt_dims(i,2))! - 1
      pdt_ilow(i) = nint(pdt_dims(i,1))! + 1
      pdt_ihig(i) = pdt_ilow(i) + int(pdt_dims(i,2)-pdt_dims(i,1))! - 1
   else
      pdt_ilow(i) = nint(pdt_dims(i,1))
      pdt_ihig(i) = pdt_ilow(i) + int(pdt_dims(i,2)-pdt_dims(i,1))
   endif
   pdt_ncells = pdt_ncells * (pdt_ihig(i)-pdt_ilow(i)+1)
   ncell_out(i) = (pdt_ihig(i)-pdt_ilow(i)+1)
!write(*,'( i3,  i3, i3 )' ) i, pdt_ilow(i), pdt_ihig(i)
enddo
!!!pdt_usr_ncell = .false.              ! Unit cells were estimates automatically
!
end subroutine estimate_ncells
!
!*******************************************************************************
!
subroutine estimate_ncatom(aver, sigma, pdt_ilow, pdt_ihig, pdt_ncells)
!-
!  Try to estimate the number of atoms per unit cell
!+
use crystal_mod
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), intent(out) :: aver
real(kind=PREC_DP), intent(out) :: sigma
integer           , dimension(3)   , intent(out)   :: pdt_ilow        ! Unit cell dimensions in periodic
integer           , dimension(3)   , intent(out)   :: pdt_ihig        ! low and high inidce
integer                            , intent(out)   :: pdt_ncells     ! Number of cells in periodic crystal volume
!
integer :: i,j,k, l   ! Dummy counter
integer :: natoms     ! Total atoms in search window
integer :: maxn       ! Max number per unit cell
integer           , dimension(3)                  :: ixyz         ! Atom is in this unit cell
integer           , dimension(0:500)              :: counter
integer           , dimension(:,:,:), allocatable :: pdt_cell     ! Atoms per cell
!
allocate(pdt_cell(pdt_ilow(1):pdt_ihig(1), pdt_ilow(2):pdt_ihig(2),pdt_ilow(3):pdt_ihig(3)))
pdt_cell = 0
do i=1, cr_natoms
   ixyz(1) = int(cr_pos(1,i) + 0.05 - pdt_ilow(1)) + pdt_ilow(1)             ! add 0.5 to avoid loosing atoms at the low edge
   ixyz(2) = int(cr_pos(2,i) + 0.05 - pdt_ilow(2)) + pdt_ilow(2)
   ixyz(3) = int(cr_pos(3,i) + 0.05 - pdt_ilow(3)) + pdt_ilow(3)
   if(ixyz(1)>=pdt_ilow(1) .and. ixyz(1)<=pdt_ihig(1)   .and.          &
      ixyz(2)>=pdt_ilow(2) .and. ixyz(2)<=pdt_ihig(2)   .and.          &
      ixyz(3)>=pdt_ilow(3) .and. ixyz(3)<=pdt_ihig(3)         ) then
      pdt_cell(ixyz(1), ixyz(2), ixyz(3)) = pdt_cell(ixyz(1), ixyz(2), ixyz(3)) + 1
   endif
enddo
!
maxn = maxval(pdt_cell)
natoms = sum(pdt_cell)
pdt_ncells =  (pdt_ihig(1)-pdt_ilow(1)+1)*(pdt_ihig(2)-pdt_ilow(2)+1)*(pdt_ihig(3)-pdt_ilow(3)+1)
!write(*,*) ' Nunit ', minval(pdt_cell), maxval(pdt_cell), sum(pdt_cell),  &
!  (pdt_ihig(1)-pdt_ilow(1)+1)*(pdt_ihig(2)-pdt_ilow(2)+1)*(pdt_ihig(3)-pdt_ilow(3)+1)
counter = 0
do k=pdt_ilow(3), pdt_ihig(3)
   do j=pdt_ilow(2), pdt_ihig(2)
      do i=pdt_ilow(1), pdt_ihig(1)
         l = pdt_cell(i,j,k)
         counter(l) = counter(l) + 1
      enddo
   enddo
enddo
aver  = 0.0
sigma = 0.0
do i=0,maxn
  aver = aver + real(i*counter(i), kind=PREC_DP)
enddo
aver = aver/real(pdt_ncells, kind=PREC_DP)
do i=0,maxn
   sigma = sigma + counter(i)*(aver-i)**2
enddo
sigma = sqrt(sigma)/real(pdt_ncells*(pdt_ncells-1), kind=PREC_DP)
!
deallocate(pdt_cell)
!
end subroutine estimate_ncatom
!
!*******************************************************************************
!
subroutine find_average_(aver, pdt_dims, pdt_ncells, pdt_nsite, pdt_asite, &
           pdt_itype, pdt_pos, MAXATOM, at_site)
!-
!  Find the clusters that form the sites in the average unit cell
!+
!
use crystal_mod
use metric_mod
!
use errlist_mod
!
implicit none
!
real(kind=PREC_DP)                               , intent(in)  :: aver       ! Average number of atoms per unit cell
real(kind=PREC_DP), dimension(3, 2)              , intent(in)  :: pdt_dims   ! Crystal dimensions
integer                                          , intent(in)  :: pdt_ncells ! Number of cells in periodic crystal volume
integer                                          , intent(out) :: pdt_nsite  ! Number of sites in an averaged unit cell
integer                                          , intent(out) :: pdt_asite  ! Achieved sites either from set site or find_aver
integer           , dimension(:,  :), allocatable, intent(out) :: pdt_itype  ! Atom types at each site
real(kind=PREC_DP), dimension(:,:)  , allocatable, intent(out) :: pdt_pos    ! Atom positions in average cell
integer                                          , intent(in)  :: MAXATOM    ! Number of atoms in crystal
integer           , dimension(MAXATOM)           , intent(out) :: at_site    ! Atom is at this site
!
logical, parameter :: lspace =.TRUE.
integer :: i, j, ic, k                 ! Dummy counters
integer :: i1,i2,i3
integer, dimension(3) :: iic
real(kind=PREC_DP) :: eps              ! Distance to say we are in the same cluster
real(kind=PREC_DP), dimension(3,4*nint(aver)) :: pdt_temp      ! Positions of the average cell
real(kind=PREC_DP), dimension(3)              :: uvw           ! fractional coordinates of atom
real(kind=PREC_DP), dimension(3)              :: u,v           ! fractional coordinates of atom
real(kind=PREC_DP), dimension(3)              :: shift         ! Shift to make sure atom coordinates are poitive
real(kind=PREC_DP) :: dist
real(kind=PREC_DP) :: dmin
integer :: istart    ! Start index for atom loop
integer :: nred      ! Number of sites to remove
integer                                       :: temp_nsite    ! Initial number of sites
integer           , dimension(  4*nint(aver)) :: pdt_temp_n    ! Number of atoms per site
logical           , dimension(  4*nint(aver)) :: is_valid      ! This site is OK T/F
logical :: lsuccess ! Test variable
!
pdt_asite = 0
pdt_nsite = nint(aver)                 ! Initialize number of sites per unit cell
eps = 0.30                             ! Start with 0.5 Angstroem
pdt_temp = 0.0_PREC_DP
pdt_temp_n = 0
at_site    = 0                         ! Atoms not yet known
shift(:) = REAL(NINT(cr_dim(:,1))+1)
!
if(pdt_asite==0) then                  ! No sites yet
   pdt_asite = 1                              ! No sites yet
   pdt_temp(:, 1) =    (cr_pos(:,1) + 2*nint(abs(pdt_dims(:,1)))) - &
                    int(cr_pos(:,1) + 2*nint(abs(pdt_dims(:,1))))            ! Arbitrarily put atom one onto site 1
   pdt_temp_n(1)  = 1                                                        ! One atom at this site
   istart = 2
else
   pdt_temp(:,1:pdt_asite) = pdt_pos(:,1:pdt_asite)
   istart = 1
endif
!
loop_atoms: do i=istart, cr_natoms
   uvw(:) = (cr_pos(:,i)+shift(:)) - real(int(cr_pos(:,i)+shift(:)), PREC_DP)   ! Create fractional
   if(uvw(1) < 0.0E0) uvw(1) = uvw(1) + 1.0                 ! coordinates in 
   if(uvw(2) < 0.0E0) uvw(2) = uvw(2) + 1.0                 ! range [0,1]
   if(uvw(3) < 0.0E0) uvw(3) = uvw(3) + 1.0
   dmin = 1e12
   ic   = 0
   do j=1, pdt_asite                      ! Loop ovr all old sites
      u = pdt_temp(:,j)
      do i3=-1,1                           ! Do each axis at: -1, 0, +1 
      v(3) = uvw(3) + i3
      do i2=-1,1                           ! to avoid rounding errors at faces
      v(2) = uvw(2) + i2
      do i1=-1,1                           ! of the unit cell
      v(1) = uvw(1) + i1
      dist = do_blen(lspace, u, v  )
      if(dist < dmin) then                 ! Found new shortest distance
         dmin = dist
         ic   = j                          ! Archive site
         iic(1) = i1
         iic(2) = i2
         iic(3) = i3
      endif
      enddo
      enddo
      enddo
   enddo
      if(dmin < eps) then                  ! Atom belongs to cluster
         uvw = uvw + iic
         pdt_temp(:, ic) = (pdt_temp(:,ic)*pdt_temp_n(ic) + uvw) / (pdt_temp_n(ic) + 1)
         pdt_temp_n( ic) = pdt_temp_n( ic) + 1
         at_site(i)      = ic              ! Atom number i is on site ic
      else                                 ! To far , new cluster
         pdt_asite = pdt_asite + 1
         if(pdt_asite>ubound(pdt_temp, 2)) then
            ier_num = -184
            ier_typ = ER_APPL
            write(ier_msg(1),'(a,i6,a)') 'More than ',pdt_asite,' sites'
            ier_msg(2) = 'Structure appears too disordered to perioditize'
            return
         endif
         pdt_temp(:, pdt_asite) = uvw
         pdt_temp_n( pdt_asite) = 1
         at_site(i)      = pdt_asite       ! Atom number i is on site pdt_asite
      endif
enddo loop_atoms
!
temp_nsite = pdt_asite
do j=1, ubound(pdt_temp_n,1)
   if(pdt_temp_n(j)>0) then
      is_valid(j)  = .true.                             ! Initial estimate, all sites are OK
   else
      is_valid(j) = .false.
   endif
enddo
!
if(pdt_asite<pdt_nsite) then                   ! Too few sites, give up
   ier_num = -179
   ier_typ = ER_APPL
   ier_msg(1) = 'Found too few site '
   return
endif
!
nred = 0
reduce_equiv: do j=1, pdt_asite
   if(pdt_asite-nred == pdt_nsite) then        ! Correct site number
      if((pdt_asite-nred)*pdt_ncells==cr_natoms) then
         exit reduce_equiv                     ! Correct site number
      else
         ier_num = -179
         ier_typ = ER_APPL
         return
      endif
   endif
   if(.not. is_valid(j)) cycle reduce_equiv
!                                              ! Too many sites, delete least occupied
   uvw  = pdt_temp(:,j)                        ! Store the 'good' position
!
   do k = j+1, temp_nsite                        ! Loop over all original sites
   dmin = 1.e12
   ic   = 0
      if(is_valid(k)) then                     ! Use valid sites only
         u = pdt_temp(:,k)
      do i3=-1,1                           ! Do each axis at: -1, 0, +1 
      v(3) = uvw(3) + i3
      do i2=-1,1                           ! to avoid rounding errors at faces
      v(2) = uvw(2) + i2
      do i1=-1,1                           ! of the unit cell
      v(1) = uvw(1) + i1
         dist = do_blen(lspace, u, v)
         if(dist < dmin) then                 ! Found new shortest distance
            dmin = dist
            ic   = k                          ! Archive site
         endif
         enddo
         enddo
         enddo
      endif
   if(dmin<1.0_PREC_DP) then                  ! Found other site close by
      pdt_temp_n(j) = pdt_temp_n(j) + pdt_temp_n(ic)
      pdt_temp_n(ic) = 0
      is_valid(ic) = .false.
      nred = nred + 1
   endif
   enddo
enddo reduce_equiv
!
reduce: do
   if(pdt_asite-nred == pdt_nsite) then        ! Correct site number
      if((pdt_asite-nred)*pdt_ncells==cr_natoms) then
         exit reduce                           ! Correct site number
      else
         ier_num = -179
         ier_typ = ER_APPL
         return
      endif
   endif
!                                              ! Too many sites, delete least occupied
   j = minloc(pdt_temp_n, dim=1,mask=is_valid)            ! Find least occupied valid site
   is_valid(j) = .false.                       ! This site is no longer valid
   uvw  = pdt_temp(:,j)                        ! Store the 'bad' position
   dmin = 1.e12
   ic   = 0
   do k = 1, temp_nsite                        ! Loop over all original sites
      if(is_valid(k)) then                     ! Use valid sites only
         u = pdt_temp(:,k)
      do i3=-1,1                           ! Do each axis at: -1, 0, +1 
      v(3) = uvw(3) + i3
      do i2=-1,1                           ! to avoid rounding errors at faces
      v(2) = uvw(2) + i2
      do i1=-1,1                           ! of the unit cell
      v(1) = uvw(1) + i1
         dist = do_blen(lspace, u, uvw)
         if(dist < dmin) then                 ! Found new shortest distance
            dmin = dist
            ic   = k                          ! Archive site
         endif
         enddo
         enddo
         enddo
      endif
   enddo
   pdt_temp_n(ic) = pdt_temp_n(ic) + pdt_temp_n(j)
   pdt_temp_n(j) = 0
   nred = nred + 1
enddo reduce
!
!
if(allocated(pdt_pos)) deallocate(pdt_pos)
if(allocated(pdt_itype)) deallocate(pdt_itype)
pdt_asite = pdt_asite-nred                     ! Set final number of sites
allocate(pdt_pos(3,pdt_asite))                 ! Atom positions in the average unit cell
pdt_pos = 0.0_PREC_DP
allocate(pdt_itype(0:cr_nscat+1, 1:pdt_asite))
pdt_itype(0,:) =  0
pdt_itype(1,:) = -1
!
! Perform loop over all atoms to get more accurate estimate of positions
!
shift(:) = REAL(-NINT(cr_dim(:,1))+2)    ! Shift to make atom position positive
pdt_temp_n    = 0                              ! One atom at this site
at_site       = 0
istart = 1
eps = 1.50_PREC_DP   ! Start with 1.0 Angstroem
!write(*,*) ' SHIFT ', shift, cr_natoms
!
loop_fine: do i=istart, cr_natoms                    ! Loop over all atoms to fine tune the average positions
   uvw(:) = (cr_pos(:,i)+shift(:)) - real(int(cr_pos(:,i)+shift(:)), PREC_DP)   ! Create fractional
   if(uvw(1) < 0.0E0) uvw(1) = uvw(1) + 1.0                 ! coordinates in 
   if(uvw(2) < 0.0E0) uvw(2) = uvw(2) + 1.0                 ! range [0,1]
   if(uvw(3) < 0.0E0) uvw(3) = uvw(3) + 1.0
!
   dmin = 1.e12
   ic   = 0
   do j=1, pdt_asite
      u = pdt_temp(:,j)
      do i3=-1,1                           ! Do each axis at: -1, 0, +1 
         v(3) = uvw(3) + i3
         do i2=-1,1                           ! to avoid rounding errors at faces
         v(2) = uvw(2) + i2
         do i1=-1,1                           ! of the unit cell
         v(1) = uvw(1) + i1
         dist = do_blen(lspace, u, v  )
         if(dist < dmin) then                 ! Found new shortest distance
            dmin = dist
            ic   = j                          ! Archive site
            iic(1) = i1
            iic(2) = i2
            iic(3) = i3
         endif
         enddo
      enddo
      enddo
   enddo
   cond_found: if(dmin < eps) then                  ! Atom belongs to average site
      uvw = uvw + iic
      pdt_pos(:,ic) = pdt_pos(:,ic) + uvw
      pdt_temp_n(ic)     = pdt_temp_n(ic) + 1   ! Increment number atoms on this site
      at_site(i)         = ic
      lsuccess = .false.
      loop_itype: do j=1, pdt_itype(0,ic)
         if(pdt_itype(j,ic)==cr_iscat(1,i)) then
            lsuccess = .true.
            exit loop_itype
         endif
      enddo loop_itype
      if(.not. lsuccess) then
         pdt_itype(0,ic) = pdt_itype(0,ic) + 1
         pdt_itype(pdt_itype(0,ic),ic) = cr_iscat(1,i)
      endif
   endif cond_found
enddo loop_fine
!
do ic=1, pdt_asite
   if(pdt_temp_n(ic)>0) pdt_pos(:,ic) = pdt_pos(:,ic)/pdt_temp_n(ic)
enddo
!*******************************************************************************
!
k = 0
do j=1, pdt_asite
    k = k + 1
    do i=1,3
       if(pdt_pos(i,k)<-0.010_PREC_DP) pdt_pos(i,k) = pdt_pos(i,k) + 1.0_PREC_DP
       if(pdt_pos(i,k)>=0.990_PREC_DP) pdt_pos(i,k) = pdt_pos(i,k) - 1.0_PREC_DP
    enddo
enddo
end subroutine find_average
!
!*******************************************************************************
!
end module dis_estimate_mod
