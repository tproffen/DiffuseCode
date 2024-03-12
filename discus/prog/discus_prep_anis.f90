module prep_anis_mod
!-
! Prepare the anisotropic ADPs
!   Determine eigenvectors and Eigenvalues from U_ij
!   Apply symmetry to Eigenvectors for all sites
!+
private
public prep_anis
public anis2iso
public lookup_anis
public do_anis
public update_biso
!
!public calc_prin
!
interface lookup_anis
   module procedure lookup_anis_6,      & !  (U11, U22,U33, U23, U13, U12)
                    lookup_anis_3x3,    & !  (U_ij)
                    lookup_anis_6_prin, & !  !U11, U22,U33, U23, U13, U12), principal vectors
                    lookup_anis_1         !  Single Biso
end interface lookup_anis
!
interface calc_prin
!  module procedure calc_prin_6, calc_prin_3x3
   module procedure              calc_prin_3x3
end interface calc_prin
!
contains
!
!*******************************************************************************
!
subroutine prep_anis(natom, l_not_full)
!
use crystal_mod
use wyckoff_mod
!
use errlist_mod
use lib_metric_mod
use matrix_mod
use math_sup, only:eigen_value, test_eigen_value
use precision_mod
use wink_mod
use support_mod
!
implicit none
!
integer, intent(in) :: natom  ! Number of atoms to cycle 
logical, intent(in) :: l_not_full  ! Do not use old cr_anis_full (from internal)
!
real(kind=PREC_DP), parameter :: TOL = 1.0D-6
integer :: itype    ! Loop index atom types
integer :: iatom    ! Loop index atoms
integer :: iref     ! Reference atom, symmetrically equivalent atoms take this atom as model
integer :: j, l     ! Dummy loop indices
logical :: lold        ! RLookup previous ADPs
logical :: lanis       ! This atom has anisotropic ADP
logical :: lsuccess    ! Success in subroutines
logical :: is_invar    ! xx tensor is invariant under local Wyckoff point group
real(kind=PREC_DP), dimension(6)   :: cr_anis_use ! copy of cr_anis or cr_anis_full
real(kind=PREC_DP), dimension(3)   :: vec     ! Dummy vector
real(kind=PREC_DP), dimension(3)   :: ar_inv  ! ( 1/a*, 1/b*, 1/c*)
real(kind=PREC_DP), dimension(3,3) :: imat  ! Unit matrix
real(kind=PREC_DP), dimension(3,3) :: uij   ! U_ij as in SHELX, reference to ((a*).a, (b*).b, (c*))c
real(kind=PREC_DP), dimension(3,3) :: xx    ! Displacement tensor <DELTAx^i . DELTAx^j> refers to standard basis
real(kind=PREC_DP), dimension(3,3) :: ucij  ! Displacement tensor <DELTAx^i . DELTAx^j> at cartesian coordinates
real(kind=PREC_DP), dimension(3,3) :: symm_mat    ! Symmetry matrix
real(kind=PREC_DP), dimension(4,3) :: prin_new ! Symmetrical principal vectors
real(kind=PREC_DP) :: ss
real(kind=PREC_DP), dimension(0:5) :: zeit
!
data imat / 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0/
!
!write(*,*) ' PREP ', cr_iscat(1,1: natom)
ar_inv(1) = 1./cr_ar(1)
ar_inv(2) = 1./cr_ar(2)
ar_inv(3) = 1./cr_ar(3)
!
lold = .TRUE.
!
if(l_not_full) then         ! Use cr_anis (read from disk)
   if(allocated(cr_anis_full)) deallocate(cr_anis_full)
   if(allocated(cr_prin     )) deallocate(cr_prin     )
   allocate(cr_anis_full(6,    cr_ncatoms))
   allocate(cr_prin     (4, 3, cr_ncatoms))
   cr_anis_full = 0.0D0
   cr_prin      = 0.0D0
   cr_nanis = 0
   cr_is_anis = .false.
!write(*,*) ' USING CR_ANIS     ', ubound(cr_anis_full,2), natom
!else
!write(*,*) ' USING CR_ANIS_FULL', ubound(cr_anis_full,2), natom
!do j=1, ubound(cr_anis_full,2)
!   write(*,'(a,6f10.6)') 'OLD UIJ ', cr_anis_full(:,j)
!enddo
endif
uij          = 0.0D0
ucij         = 0.0D0
xx           = 0.0D0
!
! Step 0:  Determine number of different atoms for a given scattering type
!
j = 0
lsuccess= .false.
cr_ndiffer = 0
!cr_is_anis = .true.
iref = 1
!
do iatom=1,natom
   l = cr_iscat(1,iatom)        ! Current atom type
   cr_ndiffer(l) = cr_ndiffer(l) + 1
enddo
!write(*,*) ' PR 0'
!do j=1, cr_nanis
!write(*,'(a, 6f9.6)') 'PR 1 ', cr_anis_full(1:6,j)
!enddo
!
! Step 1:  Determine Eigenvalues and Eigenvectors 
ss = seknds (0.0)
zeit = 0.0D0
loop_atoms:do iatom=1,natom 
!ss = seknds (ss )
   j = 0
   lsuccess= .false.
   itype = cr_iscat(1,iatom)
   cr_ianis(iatom) = iatom
   vec = cr_pos(:, iatom)
   if(l_not_full) then        ! Use cr_anis
      cr_anis_use = cr_anis(:, itype)
   else                       ! Use (old) cr_anis_full
      itype = cr_iscat(3,iatom)
      cr_anis_use = cr_anis_full(:, itype)
   endif
!   write(*,*) ' NANIS ', cr_nanis, iatom
!   write(*,'(a,   5i3, 3f7.3, i4)') ' ATOM ', iatom, itype, cr_iscat(:,iatom), cr_pos(:,iatom), cr_iscat(2, iatom)
!   write(*,'(a,        6f20.16  )') ' UANI ', cr_anis_use !(:,itype)
!   write(*,*                      ) ' UANI ', cr_anis(:,itype)
!  if((l_not_full .and. maxval(abs(cr_anis(2:,itype)))> TOL) .or.  &            ! Elements 1 to 6 are given, Anisotropic atom type
!     (.not. l_not_full .and. maxval(abs(cr_anis_full(2:,itype)))> TOL)  ) then ! Elements 1 to 6 are given, Anisotropic atom type
   if(maxval(abs(cr_anis_use(2:)))>TOL) then   ! Elements 1 to 6 are given, Anisotropic atom type
!write(*,*) ' ANISOTROPIC VALUES '
      lanis =.true.
      uij(1,1) = cr_anis_use(1) !, itype)
      uij(2,2) = cr_anis_use(2) !, itype)
      uij(3,3) = cr_anis_use(3) !, itype)
      uij(2,3) = cr_anis_use(4) !, itype)
      uij(3,2) = cr_anis_use(4) !, itype)
      uij(1,3) = cr_anis_use(5) !, itype)
      uij(3,1) = cr_anis_use(5) !, itype)
      uij(1,2) = cr_anis_use(6) !, itype)
      uij(2,1) = cr_anis_use(6) !, itype)
!     cr_iscat(3,iatom) = cr_is_sym(iatom)
      cr_iscat(3,iatom) = cr_iscat(2,iatom)
      call uij_to_xx(uij, cr_ar, xx)      ! Transform Uij to XX tensor
      call xx_to_cart(xx, cr_eimat, ucij) ! Transform XX to cartesian basis
      cr_is_anis = .true.
   else                       ! Isotropic atom set cartesian UCij and transform back to UIJ
!write(*,*) '   ISOTROPIC VALUES ', cr_anis_use     , cr_dw(itype)
      lanis =.false.
      ucij      = 0.0D0
      if(cr_anis_use(1)> 0.0D0) then                ! Element U11 is the only one take as Uiso
         ucij(1,1) = cr_anis_use(1) !,itype)
         ucij(2,2) = cr_anis_use(1) !,itype)
         ucij(3,3) = cr_anis_use(1) !,itype)
      elseif(cr_dw(itype)>0.0_PREC_DP) then
         ucij(1,1) = cr_dw(itype)/8.0D0/pi**2    ! WORK !?!??!
         ucij(2,2) = cr_dw(itype)/8.0D0/pi**2
         ucij(3,3) = cr_dw(itype)/8.0D0/pi**2
      else                                       ! Default, make sure ADP parameters are not zero
         ucij(1,1) = 1.0/8.0D0/pi**2
         ucij(2,2) = 1.0/8.0D0/pi**2
         ucij(3,3) = 1.0/8.0D0/pi**2
      endif
      cr_iscat(3,iatom ) = 1
      call xx_to_cart(ucij, cr_emat, xx)   ! Transform XX from cartesian basis
      call uij_to_xx(xx, ar_inv, uij)      ! Transform Uij from XX tensor
   endif
!ss = seknds (ss )
!zeit(0) = zeit(0) + ss
!write(*,'(a, 3f10.6)') ' Uij 1: ', uij(1, 1:3)
!write(*,'(a, 3f10.6)') ' Uij 2: ', uij(2, 1:3)
!write(*,'(a, 3f10.6)') ' Uij 3: ', uij(3, 1:3)
!write(*,*) 
!write(*,'(a, 3f10.6)') ' XX   1: ', xx(1, 1:3)
!write(*,'(a, 3f10.6)') ' XX   2: ', xx(2, 1:3)
!write(*,'(a, 3f10.6)') ' XX   3: ', xx(3, 1:3)
!write(*,*) 
!write(*,'(a, 3f10.6)') ' UCij 1: ', ucij(1, 1:3)
!write(*,'(a, 3f10.6)') ' UCij 2: ', ucij(2, 1:3)
!write(*,'(a, 3f10.6)') ' UCij 3: ', ucij(3, 1:3)
!write(*,*)
!ss = seknds (0.0)
   if(cr_iscat(2,iatom)==1) then            ! Symmetrically first atom
!write(*,*) ' XX ', xx
      call test_symm(vec, xx, is_invar)            ! Test if xx is invariant under Wyckoff symmetry
!write(*,*) ' TESTE SYMM ', ier_num, ier_typ
      if(ier_num /= 0) then
         return
      endif
   endif
!ss = seknds (ss )
!zeit(1) = zeit(1) + ss
!
!  if(cr_is_sym(iatom)>1) then             ! Symmetrically equivalent atom
!if(cr_iscat(2,iatom)==1) then
!write(*,'(a,i3, a,6f10.5, 3i4, l2)') ' UIJ ', 1, ' | ', cr_anis(:,itype), cr_iscat(:, iatom), lanis
!endif
   if(cr_iscat(2,iatom)>1) then            ! Symmetrically equivalent atom
!write(*,'(a,i3, a,6f10.5, 3i4, l2)') ' UIJ ', 1, ' | ', cr_anis(:,itype), cr_iscat(:, iatom), lanis
      if(lanis) then                       ! Anisotropic ADP, apply symmetry to obtain new ADP, principal vectors
!        symm_mat = spc_mat(1:3,1:3, cr_is_sym(iatom))
         symm_mat = spc_mat(1:3,1:3, cr_iscat(2,iatom))
!        call anis_symm(cr_ncatoms, cr_nanis, iref , symm_mat, ubound(cr_anis_full,2), &
!             cr_anis_full, cr_prin, xx, uij, cr_emat, cr_eimat, ar_inv)
!ss = seknds (ss )
         call anis_symm(symm_mat, xx, cr_prin(:,:,iref), cr_emat, cr_eimat, ar_inv, uij, prin_new)
!ss = seknds (ss )
!zeit(2) = zeit(2) + ss
!             anis_full_new, prin_new, xx_new, uij_new)
!real(kind=PREC_DP), dimension(6) :: anis_full_new
         if(ier_num/=0) then
            return
         endif
!write(*,'(a, 3f5.1)') 'sym ', symm_mat(1,:)
!write(*,'(a, 3f5.1)') 'sym ', symm_mat(2,:)
!write(*,'(a, 3f5.1)') 'sym ', symm_mat(3,:)
!write(*,'(a,i3, a,6f10.5)') ' uij ', cr_iscat(2,iatom), ' | ', uij(1,1), uij(2,2), uij(3,3), uij(2,3), uij(1,3), uij(1,2)
!ss = seknds (ss )
         call lookup_anis(lold, cr_nanis, cr_anis_full, cr_prin, uij, ucij, .FALSE., j, lsuccess)   ! Check if old ADPs are reproduced
         if(.not.lsuccess) then
            cr_prin(:,:, j) = prin_new
         endif
!ss = seknds (ss )
!zeit(3) = zeit(3) + ss
!write(*,'(a,i3, a,4f10.5,3i6, l3)') ' prin', 1, ' | ', cr_prin(:,1,j), j, cr_nanis, ubound(cr_prin,3), lsuccess
!write(*,'(a,i3, a,6f10.5)') ' prin', 2, ' | ', cr_prin(:,2,j)
!write(*,'(a,i3, a,6f10.5)') ' prin', 3, ' | ', cr_prin(:,3,j)
!        call test_symm(vec, xx, is_invar)            ! Test if xx is invariant under Wyckoff symmetry
         if(ier_num /= 0) then
            return
         endif
         cr_iscat(3,iatom) = j            ! Assign ADP type
         cycle loop_atoms
      endif 
!  else
!     iref = iatom    ! Reference atom for ADP values 
   endif 
!
!ss = seknds (ss )
   call lookup_anis(lold, cr_nanis, cr_anis_full, cr_prin, uij, ucij, .FALSE., j, lsuccess) ! Check if previous ADPs exist with identical uij
!ss = seknds (ss )
!zeit(4) = zeit(4) + ss
   if(ier_num/=0) then
      return
   endif
   cr_iscat(3,iatom) = j            ! Assign ADP type
   iref = j           ! Reference atom for ADP values
   if(lsuccess) then
!     cr_dw(itype) = (cr_anis_full(1,j) + cr_anis_full(2,j) + cr_anis_full(3,j))/3.0_PREC_DP*8.0_PREC_DP*pi*pi
      cr_dw(cr_iscat(1,iatom)) = (cr_prin(4,   1,j) + cr_prin(4,   2,j) + cr_prin(4,   3,j))/3.0_PREC_DP*8.0_PREC_DP*pi*pi
      cycle loop_atoms
   endif
!  This atom has new UIJ, we need to calculate the principal vectors and Eigenvalues
!ss = seknds (ss )
   call calc_prin(cr_iscat(3,iatom), cr_nanis, ucij, ubound(cr_prin,3), cr_prin)
!ss = seknds (ss )
!zeit(4) = zeit(4) + ss
!  cr_dw(itype) = (cr_anis_full(1,j) + cr_anis_full(2,j) + cr_anis_full(3,j))/3.0_PREC_DP*8.0_PREC_DP*pi*pi
   cr_dw(cr_iscat(1,iatom)) = (cr_prin(4,   1,j) + cr_prin(4,   2,j) + cr_prin(4,   3,j))/3.0_PREC_DP*8.0_PREC_DP*pi*pi
   if(ier_num/=0) then
      return
   endif
!write(*,'(a,i3, a,4f10.5,i6)') ' PRIN', 1, ' | ', cr_prin(:,1,cr_nanis), iref
!write(*,'(a,i3, a,6f10.5)') ' PRIN', 2, ' | ', cr_prin(:,2,cr_nanis)
!write(*,'(a,i3, a,6f10.5)') ' PRIN', 3, ' | ', cr_prin(:,3,cr_nanis)
enddo loop_atoms
cr_is_anis = .FALSE.     ! Assume isotropic
do j=1, cr_nanis
   if((abs(cr_prin(4,1,j    )-cr_prin(4,2,j   ))>TOL  .or.  &
       abs(cr_prin(4,1,j    )-cr_prin(4,3,j   ))>TOL      )) then
      cr_is_anis =.TRUE.
   endif
!write(*,'(a, 6f20.16)') 'PR 2 ', cr_anis_full(1:6,j)
enddo
!
!cr_nanis = cr_ncatoms
!write(*,*) ' cr_anis_full ', ubound(cr_anis_full)
!write(*,*) ' cr_prin      ', ubound(cr_prin     )
!do j=1, cr_nanis
!write(*,'(i3, 6f10.6)') j, cr_anis_full(:,j)
!enddo
!write(*,*)
!do j=1, cr_nanis
!write(*,'(a,3f9.5,a,f9.5,a,f9.5)') 'E(1,0,0) => ', cr_prin(1:3,1,j),  ' VAL ', cr_prin(4,1,j)
!write(*,'(a,3f9.5,a,f9.5,a,f9.5)') 'E(0,1,0) => ', cr_prin(1:3,2,j),  ' VAL ', cr_prin(4,2,j)
!write(*,'(a,3f9.5,a,f9.5,a,f9.5)') 'E(0,0,1) => ', cr_prin(1:3,3,j),  ' VAL ', cr_prin(4,3,j)
!write(*,*)
!enddo
!write(*,*) ' DONE PREP_ANIS ', cr_is_anis
!write(*,*) ' ZEIT ', zeit
!read(*,*) j
!
!
end subroutine prep_anis
!
!*******************************************************************************
!
subroutine anis2iso(nanis, idim1 , anis_full, prin, ientry, emat, ar_inv)
!-
!  Replace anisotropic ADPs by equivalent isotropic values
!+
!
use errlist_mod
use precision_mod
!
implicit none
!
integer                                       , intent(in)    :: nanis  ! Number of ADPs
integer                                       , intent(in)    :: idim1    ! Dimension of anis_full, prin
real(kind=PREC_DP), dimension(1:6,1:idim1)    , intent(inout) :: anis_full  ! List of six UIJ values
real(kind=PREC_DP), dimension(1:4,1:3,1:idim1), intent(inout) :: prin    ! Principal vectors Eigenvalues
integer                                       , intent(in)    :: ientry  ! entry to be replaced
real(kind=PREC_DP), dimension(1:3,1:3)        , intent(in)    :: emat    ! Transformation of base to cartesian
real(kind=PREC_DP), dimension(1:3)            , intent(in)    :: ar_inv  ! 1/reciprocel lattice params
!
real(kind=PREC_DP), dimension(3,3) :: ucij  ! Uij cartesian
real(kind=PREC_DP), dimension(3,3) :: xx    ! Displacement tensor crystal space
real(kind=PREC_DP), dimension(3,3) :: uij   ! Uij crystal space
real(kind=PREC_DP)                 :: uiso  ! Eqivalent Uiso
!
if(ientry > nanis) then
   ier_num = -999
   ier_typ = ER_APPL
   write(ier_msg(1),'(a5,i3,a9,i3)') 'Entry', ientry, ' outside ',nanis
   return
endif
!
uiso = (prin(4,1,ientry) + prin(4,2,ientry) + prin(4,2,ientry))/3.0_PREC_DP
ucij      = 0.0_PREC_DP
ucij(1,1) = uiso
ucij(2,2) = uiso
ucij(3,3) = uiso

call xx_to_cart(ucij, emat, xx)   ! Transform XX from cartesian basis
call uij_to_xx(xx, ar_inv, uij)      ! Transform Uij from XX tensor
call calc_prin(ientry, nanis, uij, idim1, prin)
!
end subroutine anis2iso
!
!*******************************************************************************
!
subroutine update_biso(itype, biso)
!-
!  Update the ADP parameters for all atoms of type "itype" to use isotropic biso
!-
!
use crystal_mod
!
use precision_mod
!
implicit none
!
integer           , intent(in) :: itype
real(kind=PREC_DP), intent(in) :: biso
!
integer :: i
logical, dimension(:), allocatable :: is_uses_adp
logical, dimension(:), allocatable :: adp_used_is
integer                            :: new_adp
integer                            :: ientry
logical                            :: lsuccess
!
allocate(is_uses_adp(cr_nanis))     ! itype uses these ADPs
allocate(adp_used_is(0:cr_nscat))   ! This adp is used by these atom types
is_uses_adp = .FALSE.
adp_used_is = .FALSE.
loop_atoms1: do i=1, cr_natoms
   cond_iscat1: if(cr_iscat(1,i)==itype) then
      is_uses_adp(cr_iscat(3,i)) = .TRUE.
   endif cond_iscat1
enddo loop_atoms1
!
loop_atoms2: do i=1, cr_natoms
   if(cr_iscat(1,i)==itype) cycle loop_atoms2
   cond_iscat2: if(is_uses_adp(cr_iscat(3,i))) then
      adp_used_is(cr_iscat(1,i)) = .TRUE.
   endif cond_iscat2
enddo loop_atoms2
!
new_adp = cr_nanis + 1
cond_new: if(any(adp_used_is)) then
   new_adp = cr_nanis + 1
else cond_new
   do i=1, cr_nanis
      if(is_uses_adp(i)) then
         new_adp = i
         exit cond_new
      endif
   enddo
endif cond_new
call lookup_anis_1(.FALSE., cr_nanis, cr_anis_full, cr_prin, biso , cr_emat, cr_ar, &
                         ientry, lsuccess)
!
deallocate(is_uses_adp)
deallocate(adp_used_is)
!
end subroutine update_biso
!
!*******************************************************************************
!
subroutine lookup_anis_6_prin(lold, cr_nanis, cr_anis_full, cr_prin, uij_6, prin, ientry, lsuccess)
!-
!  Check if Uij are different from any cr_anis_full entry, if so set new entry
!  Version for Uij in SHELX style, Principal axes are provided
!  Copy the input principal matrix into cr_prin
!+
!
use allocate_generic
use precision_mod
!
implicit none
!
logical                                          , intent(in)    :: lold         ! Try to find old values
integer                                          , intent(inout) :: cr_nanis
real(kind=PREC_DP), dimension(:,:)  , allocatable, intent(inout) :: cr_anis_full
real(kind=PREC_DP), dimension(:,:,:), allocatable, intent(inout) :: cr_prin
real(kind=PREC_DP), dimension(6)                 , intent(in)    :: uij_6
real(kind=PREC_DP), dimension(4,3)               , intent(in)    :: prin
integer                                          , intent(out)   :: ientry
logical                                          , intent(out)   :: lsuccess
!
logical :: lcopy    ! Copy principal vectors
real(kind=PREC_DP), dimension(3,3) :: uij
real(kind=PREC_DP), dimension(3,3) :: ucij
!
uij(1,1) = uij_6(1)
uij(2,2) = uij_6(2)
uij(3,3) = uij_6(3)
uij(3,2) = uij_6(4)
uij(2,3) = uij_6(4)
uij(3,1) = uij_6(5)
uij(1,3) = uij_6(5)
uij(2,1) = uij_6(6)
uij(1,2) = uij_6(6)
!
ucij = 0.0_PREC_DP   ! Dummy cartesian matrix
!
lcopy = .false.
!
call lookup_anis_3x3(lold, cr_nanis, cr_anis_full, cr_prin, uij, ucij, lcopy, ientry, lsuccess)
cr_prin(:,:,ientry) = prin(:,:)
!
end subroutine lookup_anis_6_prin
!
!*******************************************************************************
!
subroutine lookup_anis_6(lold, cr_nanis, cr_anis_full, cr_prin, uij_6, cr_eimat, cr_ar, &
           ientry, lsuccess)
!-
!  Check if Uij are different from any cr_anis_full entry, if so set new entry
!  Version for Uij in SHELX style
!  The principal vectors will be calculated downstreams in lookup_anis_3x3
!+
!
use allocate_generic
use precision_mod
!
implicit none
!
logical                                          , intent(in)    :: lold         ! Try to find old values
integer                                          , intent(inout) :: cr_nanis
real(kind=PREC_DP), dimension(:,:)  , allocatable, intent(inout) :: cr_anis_full
real(kind=PREC_DP), dimension(:,:,:), allocatable, intent(inout) :: cr_prin
real(kind=PREC_DP), dimension(6)                 , intent(in)    :: uij_6
real(kind=PREC_DP), dimension(3,3)               , intent(in)    :: cr_eimat
real(kind=PREC_DP), dimension(3)                 , intent(in)    :: cr_ar
integer                                          , intent(out)   :: ientry
logical                                          , intent(out)   :: lsuccess
!
logical :: lcopy    ! Copy principal vectors
real(kind=PREC_DP), dimension(3,3) :: uij
real(kind=PREC_DP), dimension(3,3) :: ucij
real(kind=PREC_DP), dimension(3,3) :: xx
!
uij(1,1) = uij_6(1)
uij(2,2) = uij_6(2)
uij(3,3) = uij_6(3)
uij(3,2) = uij_6(4)
uij(2,3) = uij_6(4)
uij(3,1) = uij_6(5)
uij(1,3) = uij_6(5)
uij(2,1) = uij_6(6)
uij(1,2) = uij_6(6)
!
call uij_to_xx(uij, cr_ar, xx)      ! Transform Uij to XX tensor
call xx_to_cart(xx, cr_eimat, ucij) ! Transform XX to cartesian basis
!
lcopy = .true.
!
call lookup_anis_3x3(lold, cr_nanis, cr_anis_full, cr_prin, uij, ucij, lcopy, ientry, lsuccess)
!
end subroutine lookup_anis_6
!
!*******************************************************************************
!
subroutine lookup_anis_1(lold, cr_nanis, cr_anis_full, cr_prin, cr_dw, cr_emat, cr_ar, &
                         ientry, lsuccess)
!-
!  Check if Uij are different from any cr_anis_full entry, if so set new entry
!  Version for a single Biso value
!+
!
!use crystal_mod
!
use allocate_generic
use precision_mod
use wink_mod
!
implicit none
!
logical                                          , intent(in)    :: lold         ! Try to find old values
integer                                          , intent(inout) :: cr_nanis
real(kind=PREC_DP), dimension(:,:)  , allocatable, intent(inout) :: cr_anis_full
real(kind=PREC_DP), dimension(:,:,:), allocatable, intent(inout) :: cr_prin
real(kind=PREC_DP)                               , intent(in)    :: cr_dw
real(kind=PREC_DP), dimension(3,3)               , intent(in)    :: cr_emat  ! Transformation to cartesian
real(kind=PREC_DP), dimension(3)                 , intent(in)    :: cr_ar    ! Reciprocal lattice params
integer                                          , intent(out)   :: ientry
logical                                          , intent(out)   :: lsuccess
!
logical :: lcopy    ! Copy principal vectors
real(kind=PREC_DP), dimension(3)   :: ar_inv  ! 1./reciprocal lattice parameters
real(kind=PREC_DP), dimension(3,3) :: uij     ! Crystal space
real(kind=PREC_DP), dimension(3,3) :: ucij    ! Cartesian space
real(kind=PREC_DP), dimension(3,3) :: xx      ! Temporary tensor
! 
ar_inv(1) = 1./cr_ar(1)
ar_inv(2) = 1./cr_ar(2)
ar_inv(3) = 1./cr_ar(3)
!
ucij = 0.0_PREC_DP 
ucij(1,1) = cr_dw/8.0D0/pi**2    ! WORK !?!??!
ucij(2,2) = cr_dw/8.0D0/pi**2
ucij(3,3) = cr_dw/8.0D0/pi**2
!cr_iscat(3,iatom ) = 1
!
call xx_to_cart(ucij, cr_emat, xx)   ! Transform UCij to XX (cartesian => crystal)
call uij_to_xx(xx, ar_inv, uij)      ! Transform xx tensor into Uij
!
lcopy = .true.
!
call lookup_anis_3x3(lold, cr_nanis, cr_anis_full,  cr_prin, uij, ucij, lcopy, ientry, lsuccess)
!
end subroutine lookup_anis_1
!
!*******************************************************************************
!
subroutine lookup_anis_3x3(lold, cr_nanis, cr_anis_full,  cr_prin, uij, ucij, lcopy, ientry, lsuccess)
!-
!  Check if Uij are different from any cr_anis_full entry, if so set new entry
!  Version for Uij in as 3x3 matrix
!+
!
use allocate_generic
use errlist_mod
use precision_mod
!
implicit none
!
logical                                          , intent(in)    :: lold         ! Try to find old values
integer                                          , intent(inout) :: cr_nanis
real(kind=PREC_DP), dimension(:,:)  , allocatable, intent(inout) :: cr_anis_full
real(kind=PREC_DP), dimension(:,:,:), allocatable, intent(inout) :: cr_prin
real(kind=PREC_DP), dimension(3,3)               , intent(in)    :: uij       ! Crystal space
real(kind=PREC_DP), dimension(3,3)               , intent(in)    :: ucij      ! Cartesian space
logical                                          , intent(in)    :: lcopy    ! Copy principal vectors
integer                                          , intent(out)   :: ientry
logical                                          , intent(out)   :: lsuccess
!
real(kind=PREC_DP), parameter :: TOL = 1.0D-6
!
integer :: idim1
integer :: j 
integer :: all_status
!
ientry = -1
lsuccess = .FALSE.
!
if(lold) then
loop_nanis: do j=1, cr_nanis
   if(&
     (abs(cr_anis_full(1,j)-uij(1,1))<TOL .and.  &
      abs(cr_anis_full(2,j)-uij(2,2))<TOL .and.  &
      abs(cr_anis_full(3,j)-uij(3,3))<TOL .and.  &
      abs(cr_anis_full(4,j)-uij(3,2))<TOL .and.  &
      abs(cr_anis_full(5,j)-uij(3,1))<TOL .and.  &
      abs(cr_anis_full(6,j)-uij(2,1))<TOL      ) &
      .or.                                     &
     (abs(cr_anis_full(1,j)+uij(1,1))<TOL .and.  &
      abs(cr_anis_full(2,j)+uij(2,2))<TOL .and.  &
      abs(cr_anis_full(3,j)+uij(3,3))<TOL .and.  &
      abs(cr_anis_full(4,j)+uij(3,2))<TOL .and.  &
      abs(cr_anis_full(5,j)+uij(3,1))<TOL .and.  &
      abs(cr_anis_full(6,j)+uij(2,1))<TOL      ) &
      ) then  
      ientry = j
      lsuccess = .TRUE.
      exit loop_nanis
   endif
enddo loop_nanis
endif
!
!write(*,*) ' SUCCESS ', lsuccess, cr_nanis, lcopy
if(ientry==-1) then   ! new entry
   cr_nanis = cr_nanis + 1
   ientry   = cr_nanis
   if(cr_nanis>ubound(cr_anis_full,2)) then
      call alloc_arr(cr_anis_full, 1,6, 1, cr_nanis, all_status, 0.0D0      )
      call alloc_arr(cr_prin, 1, 4,1,3, 1, cr_nanis, all_status, 0.0D0      )
   endif
   cr_anis_full(1,ientry) = uij(1,1)
   cr_anis_full(2,ientry) = uij(2,2)
   cr_anis_full(3,ientry) = uij(3,3)
   cr_anis_full(4,ientry) = uij(3,2)
   cr_anis_full(5,ientry) = uij(3,1)
   cr_anis_full(6,ientry) = uij(2,1)
   do j=1, 6
!write(*,*) abs(cr_anis_full(j,ientry)), abs(cr_anis_full(j,ientry))<TOL
      if(abs(cr_anis_full(j,ientry))<TOL) then
         cr_anis_full(j,ientry) = 0.0D0
      endif
   enddo
   if(lcopy) then    ! Need to calculate new principal vectors and copy into global storage
      idim1 = cr_nanis
      call calc_prin_3x3(ientry, cr_nanis, ucij, idim1, cr_prin)
   endif
!else
!write(*,*) ' found old ', ientry
endif
!
end subroutine lookup_anis_3x3
!
!*******************************************************************************
!
!QQsubroutine calc_prin_6(ientry, nanis, uij_6, dmat, emat, gten, prin)
!QQ!-
!QQ! Calculate principal axis for Uij
!QQ!-
!QQ!
!QQuse errlist_mod
!QQuse lib_metric_mod
!QQuse math_sup, only:eigen_value, test_eigen_value
!QQuse precision_mod
!QQ!
!QQimplicit none
!QQ!
!QQinteger                               , intent(in)  :: ientry
!QQinteger                               , intent(in)  :: nanis
!QQreal(kind=PREC_DP), dimension(6)      , intent(in)  :: uij_6
!QQreal(kind=PREC_DP), dimension(3,3)    , intent(in)  :: dmat
!QQreal(kind=PREC_DP), dimension(3,3)    , intent(in)  :: emat
!QQreal(kind=PREC_DP), dimension(3,3)    , intent(in)  :: gten
!QQreal(kind=PREC_DP), dimension(4,3,nanis), intent(out) :: prin
!QQ!
!QQreal(kind=PREC_DP), dimension(3,3) :: uij
!QQ!
!QQuij(1,1) = uij_6(1)
!QQuij(2,2) = uij_6(2)
!QQuij(3,3) = uij_6(3)
!QQuij(3,2) = uij_6(4)
!QQuij(2,3) = uij_6(4)
!QQuij(3,1) = uij_6(5)
!QQuij(1,3) = uij_6(5)
!QQuij(2,1) = uij_6(6)
!QQuij(1,2) = uij_6(6)
!QQ!
!QQcall calc_prin_3x3(ientry, nanis, uij, prin)
!QQ!
!QQend subroutine calc_prin_6
!
!*******************************************************************************
!
subroutine calc_prin_3x3(ientry, nanis, ucij, idim1, prin)
!-
! Calculate principal axis for Uij
!-
!
!use crystal_mod
!
use errlist_mod
use lib_metric_mod
use math_sup, only:eigen_value, test_eigen_value
use precision_mod
!
implicit none
!
integer                               , intent(in)  :: ientry
integer                               , intent(in)  :: nanis
real(kind=PREC_DP), dimension(3,3)    , intent(in)  :: ucij    ! Displacement tensor at Cartesian
integer                               , intent(in)  :: idim1   ! Dimension prin
real(kind=PREC_DP), dimension(4,3,idim1), intent(out) :: prin
!
integer :: neigen
!integer :: i,j,k,l
real(kind=PREC_DP), dimension(3)   :: eigen_val   ! U_ij in cartesian coordinates
real(kind=PREC_DP), dimension(3,3) :: eigen_car   ! Eigenvectors in cartesian coordinates
real(kind=PREC_DP), DIMENSION(3,3)          ::      imat     = &
   reshape((/1.0D0,0.0D0,0.0D0, 0.0D0,1.0D0,0.0D0, 0.0D0,0.0D0,1.0D0/),SHAPE(imat))
!
call eigen_value(ucij, eigen_val, eigen_car, imat, neigen)
if(ier_num /=0 ) then
   ier_msg(1) = 'Cartesian UCij matrix has zero det.'
   ier_msg(2) = 'Check anisotropic ADP values '
   return
endif
prin  (1:3,1, ientry) = (eigen_car(:,1))
prin  (1:3,2, ientry) = (eigen_car(:,2))
prin  (1:3,3, ientry) = (eigen_car(:,3))
prin  (4,:  , ientry) = eigen_val
!write(*,*) ' DID CALC_PRIN_3x3 ', ier_num, ier_typ, ientry, nanis
!
end subroutine calc_prin_3x3
!
!*******************************************************************************
!
subroutine prep_anis_iscat
!-
!  Map the scattering curve etc onto the anis list
!+
use crystal_mod
!
use precision_mod
!
implicit none
!
integer i
!
if(cr_nanis==0) return
!
deallocate(cr_dw)
allocate(cr_dw(0:max(MAXSCAT,cr_nanis)))
cr_dw = 0.0D0
do i=1, cr_nanis
   cr_dw(i) = cr_ianis(i)
enddo
!
end subroutine prep_anis_iscat
!
!*******************************************************************************
!
subroutine uij_to_xx(uij, cr_ar, xx)   ! Transform Uij to XX tensor
!-
!  A simple routine that converts Uij to XX and vice versa
!  xx  = cr_ar * UIJ * cr_ar
!  UIJ = (1/cr_ar) * XX *  (1./cr_ar)
!+
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3,3), intent(in)  :: uij   ! Input matrix
real(kind=PREC_DP), dimension(3)  , intent(in)  :: cr_ar ! Elements are reciprocal lengths; or 1/reciprocal_lengths
real(kind=PREC_DP), dimension(3,3), intent(out) :: xx    ! Result matrix
!
integer :: j,k
!
xx = 0.0D0
do j=1, 3
  do k=1,3
    xx(j,k) = uij(j,k)*cr_ar(j)*cr_ar(k)
  enddo
enddo
!
end subroutine uij_to_xx
!
!*******************************************************************************
!
subroutine xx_to_cart(xx, eimat, ucij) ! Transform XX to cartesian basis
!-
! A simple 1 line routine that transforms XX to cartesian space or back
!  call xx_to_cart(xx, eimat, ucij)
!  ucij = matmul((eimat), matmul(xx,   transpose(eimat)))  ! With eimat into cartesian
!
!  call xx_to_cart(ucij, emat, xx)
!  xx   = matmul((emat) , matmul(ucij, transpose(emat)))   ! With emat into crystal space
!+
!
use precision_mod
!implicit none
!
real(kind=PREC_DP), dimension(3,3), intent(in)  :: xx    ! Input matrix
real(kind=PREC_DP), dimension(3,3), intent(in)  :: eimat ! transformation matrix
real(kind=PREC_DP), dimension(3,3), intent(out) :: ucij  ! Result matrix
!
ucij = matmul((eimat), matmul(xx, transpose(eimat)))
!
end subroutine xx_to_cart
!
!*******************************************************************************
!
!subroutine anis_symm(ncatoms, symm_mat, prin_in, xx_in, emat, eimat, ar_inv, &
!           anis_full_new, prin_new, xx_new, uij_new)
subroutine anis_symm(symm_mat, xx_in, prin_in, emat, eimat, ar_inv, uij_new, prin_new)
!-
!  Apply the symmetry matrix symm_mat to the current displacement tensor XX and the 
!  principal vectors
!+
!
use matrix_mod
use precision_mod
!
implicit none
!
!integer                           , intent(in)  :: ncatoms   ! Number of atoms in unit cell, 
real(kind=PREC_DP), dimension(3,3), intent(in)  :: symm_mat  ! Symmetry matrix
real(kind=PREC_DP), dimension(4,3)  , intent(in)  :: prin_in      ! Principal vectors cartesian
real(kind=PREC_DP), dimension(3,3), intent(in)  :: xx_in        ! Displacement tensor crystal space
real(kind=PREC_DP), dimension(3,3), intent(in)  :: emat      ! Transformation basis to Cartesian
real(kind=PREC_DP), dimension(3,3), intent(in)  :: eimat     ! Transformation basis to crystal
real(kind=PREC_DP), dimension(3)  , intent(in)  :: ar_inv    ! Reciprocal lattice parameters
!real(kind=PREC_DP), dimension(6)  , intent(out) :: anis_full_new ! List of UIJ
!real(kind=PREC_DP), dimension(4)  , intent(out) :: prin_new     ! Principal vectors cartesian
!real(kind=PREC_DP), dimension(3,3), intent(out) :: xx_new       ! Displacement tensor crystal space
real(kind=PREC_DP), dimension(3,3), intent(out) :: uij_new       ! Full UIJ matrix
real(kind=PREC_DP), dimension(4,3), intent(out) :: prin_new ! Symmetrical principal vectors
!
integer                            :: j
real(kind=PREC_DP), dimension(3)   :: v,w  ! Dummy vector
real(kind=PREC_DP), dimension(3,3) :: xx_new       ! Displacement tensor after symmetry operation
real(kind=PREC_DP), dimension(3,3) :: symm_imat    ! Inverse Symmetry matrix
!
symm_imat = 0.0_PREC_DP
call matinv(symm_mat, symm_imat)
if(ier_num/=0) then
   ier_msg(1) = 'Could not invert symmetry operation'
   return
endif
!
!write(*,*) ' SYMMETRY ', ncatoms, nanis, iref
!write(*,'(a,4f10.6)') ' OLD  1', prin(:,1,iref )
!write(*,'(a,4f10.6)') ' OLD  2', prin(:,2,iref )
!write(*,'(a,4f10.6)') ' OLD  3', prin(:,3,iref )
do j=1, 3  ! loop over three principal vectors
   v = matmul(emat, prin_in(1:3,j))    ! Transform principal vector to crystal space
   w = matmul(symm_mat,v)            ! Apply symmetry matrix
!write(*,'(a,i1,2(3f10.6,4x))') ' crys ',j, v, w
   v = matmul(eimat, w)              ! Transform back to cartesian space
   prin_new(1:3,j) = v
   prin_new(4  ,j) = prin_in(4,j)
enddo
! write(*,'(a,4f10.6)') ' pRIN 1', prin(:,1,nanis+1)
! write(*,'(a,4f10.6)') ' pRIN 2', prin(:,2,nanis+1)
! write(*,'(a,4f10.6)') ' pRIN 3', prin(:,3,nanis+1)
!
 xx_new = matmul(symm_mat, matmul(xx_in, transpose(symm_mat)))
call uij_to_xx(xx_new, ar_inv, uij_new)
!xx = xx_new
!anis_full_new(1,nanis+1) = uij_new(1,1)
!anis_full_new(2,nanis+1) = uij_new(2,2)
!anis_full_new(3,nanis+1) = uij_new(3,3)
!anis_full_new(4,nanis+1) = uij_new(2,3)
!anis_full_new(5,nanis+1) = uij_new(1,3)
!anis_full_new(6,nanis+1) = uij_new(1,2)
!
end subroutine anis_symm
!
!*******************************************************************************
!
subroutine old_anis_symm(ncatoms, nanis, iref, symm_mat, idim1, anis_full, &
   prin, xx, uij, emat, eimat, ar_inv)
!-
!  Apply the symmetry matrix symm_mat to the current displacement tensor XX and the 
!  principal vectors
!+
!
use matrix_mod
use precision_mod
!
implicit none
!
integer                                    , intent(in)    :: ncatoms   ! Number of atoms in unit cell, 
                                                                        ! dimension of cr_anis_full, cr_prin
integer                                    , intent(inout) :: nanis     ! Number of different anisotropic ADPs
integer                                    , intent(in)    :: iref      ! Reference atom
real(kind=PREC_DP), dimension(3,3)         , intent(in)    :: symm_mat  ! Symmetry matrix
integer                                    , intent(in)    :: idim1     ! Reference atom
real(kind=PREC_DP), dimension(6,   idim1  ), intent(inout) :: anis_full ! List of UIJ
real(kind=PREC_DP), dimension(4,3, idim1  ), intent(inout) :: prin      ! Principal vectors cartesian
real(kind=PREC_DP), dimension(3,3)         , intent(inout) :: uij       ! Full UIJ matrix
real(kind=PREC_DP), dimension(3,3)         , intent(inout) :: xx        ! Displacement tensor crystal space
real(kind=PREC_DP), dimension(3,3)         , intent(in)    :: emat      ! Transformation basis to Cartesian
real(kind=PREC_DP), dimension(3,3)         , intent(in)    :: eimat     ! Transformation basis to crystal
real(kind=PREC_DP), dimension(3)           , intent(in)    :: ar_inv    ! Reciprocal lattice parameters
!
integer                          :: j
real(kind=PREC_DP), dimension(3) :: v,w  ! Dummy vector
real(kind=PREC_DP), dimension(3,3) :: xx_new       ! Displacement tensor after symmetry operation
!real(kind=PREC_DP), dimension(3,3) :: uij          ! UIJ in (a*.a, b*.b, c*.c) space
real(kind=PREC_DP), dimension(3,3) :: symm_imat    ! Inverse Symmetry matrix
!
symm_imat = 0.0_PREC_DP
call matinv(symm_mat, symm_imat)
if(ier_num/=0) then
   ier_msg(1) = 'Could not invert symmetry operation'
   return
endif
!
!write(*,*) ' SYMMETRY ', ncatoms, nanis, iref
!write(*,'(a,4f10.6)') ' OLD  1', prin(:,1,iref )
!write(*,'(a,4f10.6)') ' OLD  2', prin(:,2,iref )
!write(*,'(a,4f10.6)') ' OLD  3', prin(:,3,iref )
do j=1, 3  ! loop over three principal vectors
   v = matmul(emat, prin(1:3, j, iref))    ! Transform principal vector to crystal space
!  w = matmul(symm_imat,v)                 ! Apply symmetry matrix
   w = matmul(symm_mat,v)                 ! Apply symmetry matrix
!write(*,'(a,i1,2(3f10.6,4x))') ' crys ',j, v, w
   v = matmul(eimat, w)                    ! Transform back to cartesian space
   prin(1:3, j, nanis+1) = v
   prin(4, j, nanis+1) = prin(4, j, iref)
enddo
! write(*,'(a,4f10.6)') ' pRIN 1', prin(:,1,nanis+1)
! write(*,'(a,4f10.6)') ' pRIN 2', prin(:,2,nanis+1)
! write(*,'(a,4f10.6)') ' pRIN 3', prin(:,3,nanis+1)
!
 xx_new = matmul(symm_mat, matmul(xx, transpose(symm_mat)))
!xx_new = matmul(transpose(symm_mat), matmul(xx,          (symm_mat)))
call uij_to_xx(xx_new, ar_inv, uij)
xx = xx_new
anis_full(1,nanis+1) = uij(1,1)
anis_full(2,nanis+1) = uij(2,2)
anis_full(3,nanis+1) = uij(3,3)
anis_full(4,nanis+1) = uij(2,3)
anis_full(5,nanis+1) = uij(1,3)
anis_full(6,nanis+1) = uij(1,2)
!
end subroutine old_anis_symm
!
!*******************************************************************************
!
subroutine test_symm(vec, xx, is_invar)
!-
!  Test if xx is invariant under local Wycoff symmetry at position vec
!+
!
use errlist_mod
use wyckoff_mod
use spcgr_apply, only:get_wyckoff
!
use precision_mod
!implicit none
!
real(kind=PREC_DP), dimension(3)  , intent(inout) :: vec   ! Atom position
real(kind=PREC_DP), dimension(3,3), intent(in)    :: xx    ! Displacement tensor <DELTAx^i . DELTAx^j> refers to standard basis
logical                           , intent(out)   :: is_invar
!
real(kind=PREC_DP), parameter :: TOL=0.00002_PREC_DP
!
integer :: i
real(kind=PREC_DP), dimension(3,3) :: mat   ! Wyckoff symmetry without translation
real(kind=PREC_DP), dimension(3,3) :: new   ! Transformed xx
!
call get_wyckoff(vec, .FALSE., 0)
!
is_invar = .TRUE.
!write(*,*) ' POSIT    ', vec, is_invar
loop_wyc: do i=2, wyc_n                     ! Loop over all Wyckoff symmetries no. 2 and up (1 is identity)
   mat = wyc_mat(1:3,1:3, i)
   new = matmul(mat, matmul(xx, transpose(mat)))
!  new = matmul(transpose(mat), matmul(xx,          (mat)))
!write(*,'(2(a, 3(2x,f9.6)),l4)') ' XX ', xx(1,:), ' new ', new(1,:), all(abs(xx-new)<TOL)
!write(*,'(2(a, 3(2x,f9.6)))') ' XX ', xx(2,:), ' new ', new(2,:)
!write(*,'(2(a, 3(2x,f9.6)))') ' XX ', xx(3,:), ' new ', new(3,:)
   is_invar = is_invar .and. all(abs(xx-new)<TOL)
!  if(.not. is_invar) exit loop_wyc
enddo loop_wyc
!s_invar = .true.
if(.not. is_invar) then
!write(*,*) ' XX ', xx
   ier_num = -195
   ier_typ = ER_APPL
   write(ier_msg(1),'(a,3f10.6)') 'at atom pos: ',vec
   ier_msg(2) = ' Uij do not comply with Wyckoff symmetry'
!write(*,*) ' POSIT    ', vec, is_invar
!write(*,'(2(a, 3(2x,f9.6)))') ' XX ', xx(1,:), ' new ', new(1,:)
!write(*,'(2(a, 3(2x,f9.6)))') ' XX ', xx(2,:), ' new ', new(2,:)
!write(*,'(2(a, 3(2x,f9.6)))') ' XX ', xx(3,:), ' new ', new(3,:)
endif
!
end subroutine test_symm
!
!*******************************************************************************
!
subroutine do_anis(line, length)
!-
!  Set anisotropic parameters from command line
!+
!
use crystal_mod
!
use errlist_mod
use get_params_mod
use precision_mod
use take_param_mod
!
character(len=*), intent(inout) :: line
integer         , intent(inout) :: length
!
integer, parameter :: MAXW  = 6
integer, parameter :: MAXWW = 6
real(kind=PREC_DP), parameter :: TOL=0.00002_PREC_DP
!
character(len=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
integer                                  :: ianz
character(len=PREC_STRING)               :: ccpara
integer                                  :: llpara
real(kind=PREC_DP) , dimension(6)        :: wwerte   ! Calculated values
!
logical                                  :: loutput  ! Output to screen ?
!
integer :: ianis, j
real(kind=PREC_DP)                 :: rlen    ! Vector length
real(kind=PREC_DP), dimension(3)   :: ar_inv  ! ( 1/a*, 1/b*, 1/c*)
real(kind=PREC_DP), dimension(3)   :: vec     ! a vector
real(kind=PREC_DP), dimension(3,3) :: uij   ! U_ij as in SHELX, reference to ((a*).a, (b*).b, (c*))c
real(kind=PREC_DP), dimension(3,3) :: xx    ! Displacement tensor <DELTAx^i . DELTAx^j> refers to standard basis
real(kind=PREC_DP), dimension(3,3) :: ucij  ! Displacement tensor <DELTAx^i . DELTAx^j> at cartesian coordinates
!
integer, parameter :: NOPTIONAL = 5
integer, parameter :: O_TYPE    = 1
integer, parameter :: O_VALUES  = 2
integer, parameter :: O_WYCKOFF = 3
integer, parameter :: O_OUTPUT  = 4
integer, parameter :: O_PRIN    = 5
character(LEN=   7), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 1 ! Number of values to calculate 
!
data oname  / 'type', 'values', 'wyckoff', 'output' , 'prin'  /
data loname /  4,      6      ,  7       ,  6       ,  4      /
opara  =  (/ '1.0000       ', '[0.01]       ', '[0.0,0.0,0.0]', 'screen       ', '[0.0,0.0,0.0]' /)   ! Always provide fresh default values
lopara =  (/  6,                6            ,  13            ,  6             ,  13             /)
owerte =  (/  0.0           ,   0.0          ,  0.0           ,  0.0           ,  0.0            /)
!
call get_params (line, ianz, cpara, lpara, MAXW, length)
if(ier_num /= 0) return
!
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num /= 0) return
loutput = opara(O_OUTPUT) == 'screen'
!
if(lpresent(O_WYCKOFF)) then
   ccpara =  opara(O_WYCKOFF)
   llpara = lopara(O_WYCKOFF)
   call get_optional_multi(MAXWW, ccpara, llpara, wwerte, ianz)
   if(ier_num /= 0) return
!
   call do_anis_wyc(MAXWW, wwerte, ianz, loutput)
   return
elseif(lpresent(O_PRIN)) then
   ccpara =  opara(O_PRIN)
   llpara = lopara(O_PRIN)
   call get_optional_multi(MAXWW, ccpara, llpara, wwerte, ianz)
   if(ier_num /= 0) return
!
   j     = nint(wwerte(1))
   ianis = nint(wwerte(2))
   vec                   = matmul(cr_eimat, wwerte(3:5))  ! Direct space to cartesian
   rlen                  = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
   cr_prin(1:3,j, ianis) = vec/rlen       ! Normalize to length 1
   cr_prin(4  ,j, ianis) = wwerte(6)
!  cr_dw(itype) = (cr_anis_full(1,j) + cr_anis_full(2,j) + cr_anis_full(3,j))/3.0_PREC_DP*8.0_PREC_DP*pi*pi
   return
endif
!
ccpara =  opara(O_VALUES)
llpara = lopara(O_VALUES)
call get_optional_multi(MAXWW, ccpara, llpara, wwerte, ianz)
if(ier_num /= 0) return
!
uij  = 0.0_PREC_DP
ucij = 0.0_PREC_DP
!
if(ianz==6) then     ! All six values, ANISOTROPIC
   uij(1,1) = wwerte(1)
   uij(2,2) = wwerte(2)
   uij(3,3) = wwerte(3)
   uij(3,2) = wwerte(4)
   uij(2,3) = wwerte(4)
   uij(1,3) = wwerte(5)
   uij(3,1) = wwerte(5)
   uij(1,2) = wwerte(6)
   uij(2,1) = wwerte(6)
   call uij_to_xx(uij, cr_ar, xx)      ! Transform Uij to XX tensor
   call xx_to_cart(xx, cr_eimat, ucij) ! Transform XX to cartesian basis
elseif(ianz==1) then       ! One parameter ISOTROPIC
!
   ar_inv(1) = 1./cr_ar(1)
   ar_inv(2) = 1./cr_ar(2)
   ar_inv(3) = 1./cr_ar(3)
   ucij(1,1) = wwerte(1)                ! Set the cartesian diagonal values
   ucij(2,2) = wwerte(1)
   ucij(3,3) = wwerte(1)
   call xx_to_cart(ucij, cr_emat, xx)   ! Transform XX from cartesian basis
   call uij_to_xx(xx, ar_inv, uij)      ! Transform Uij from XX tensor
else
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = '''value:'' needs 1 or 6 parameters'
endif
!
ianis = nint(owerte(O_TYPE))
cr_anis(:,ianis) = wwerte(1:6)   ! WORK  
call calc_prin_3x3(ianis, cr_nanis, ucij, ubound(cr_prin,3), cr_prin)
cr_anis_full(1, ianis) = uij(1,1)
cr_anis_full(2, ianis) = uij(2,2)
cr_anis_full(3, ianis) = uij(3,3)
cr_anis_full(4, ianis) = uij(2,3)
cr_anis_full(5, ianis) = uij(1,3)
cr_anis_full(6, ianis) = uij(1,2)
   do j=1, 6
!write(*,*) abs(cr_anis_full(j,ianis)), abs(cr_anis_full(j,ianis))<TOL
      if(abs(cr_anis_full(j,ianis))<TOL) then
         cr_anis_full(j,ianis) = 0.0D0
      endif
   enddo
!
end subroutine do_anis
!
!*******************************************************************************
!
subroutine do_anis_wyc(MAXWW, wwerte, ianz, loutput)
!-
!  Construct UIJ that  comply with Wyckoff symmetry at position wwerte(1:3)
!+
!
use crystal_mod
use wyckoff_mod
use spcgr_apply, only:get_wyckoff
!
use precision_mod
use errlist_mod
use lib_metric_mod
use matrix_mod
use param_mod
use prompt_mod
!
integer                             , intent(in) :: MAXWW   ! Array dimension
real(kind=PREC_DP), dimension(MAXWW), intent(in) :: wwerte  ! Position vector
integer                             , intent(in) :: ianz    ! dimension should be 3
logical                             , intent(in) :: loutput ! Output to screen?
!
character(len=PREC_STRING), dimension(2)  :: line
character(len=PREC_STRING), dimension(10) :: rules
integer                   , dimension(10) :: irule
integer                          :: i,j     ! Dummy indices
real(kind=PREC_DP), dimension(3) :: vec     ! Position vector
real(kind=PREC_DP), dimension(3) :: w,u     ! Position vector
real(kind=PREC_DP)               :: vvv     ! Position vector, length
integer           , dimension(-6:6) :: nsym    ! Number or axes with this symmetry ( -6, _, -4, -3, -2, -1, 0, 1, 2, 3, 4, m, 6 )
real(kind=PREC_DP), dimension(4,3) :: prin   ! Principal vectors
real(kind=PREC_DP), dimension(3,3) :: mat    ! Symmetry     matrix cartesian
real(kind=PREC_DP), dimension(3,3) :: matt   ! Symmetry^T   matrix 
real(kind=PREC_DP), dimension(3,3) :: ucij   ! Displacement matrix cartesian
real(kind=PREC_DP), dimension(3,3) :: uij    ! Displacement matrix a.a* base
real(kind=PREC_DP), dimension(3,3) :: xx     ! Displacement matrix crystal space
real(kind=PREC_DP), dimension(3,3) :: xxp    ! Displacement matrix crystal space
real(kind=PREC_DP), dimension(3,3) :: temp   ! Displacement matrix crystal space
real(kind=PREC_DP), dimension(3)   :: ar_inv  ! ( 1/a*, 1/b*, 1/c*)
!
temp = 0.0D0
temp(1,1) = -1.0D0
temp(2,2) = -1.0D0
temp(3,3) = -1.0D0
if(ianz /= 3) then
   ier_num = -6
   ier_typ =ER_COMM
   ier_msg(1) = 'Wyckoff parameter requires 3 coordinates'
   return
endif 
!
ar_inv(1) = 1./cr_ar(1)
ar_inv(2) = 1./cr_ar(2)
ar_inv(3) = 1./cr_ar(3)
!
vec = wwerte(1:3)
prin = 0.0_PREC_DP
call get_wyckoff (vec, .false., 0)
!
! Build a "random" displacement matrix
xxp     = 0.0000000_PREC_DP
ucij = 0.0D0
vec = (/ 0.95D0, 0.10D0, -0.05d0 /)
vvv = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
vec = vec/vvv
ucij(1,:) = vec * 0.0062000D0
w   = (/ 0.10D0, 0.53D0,  0.77D0 /)
call lib_vector_product(vec, w, u)
vvv = sqrt(  u(1)**2 +   u(2)**2 +   u(3)**2)
u   = u/vvv
ucij(2,:) =  -u * 0.0050000D0
call lib_vector_product(vec, u, w)
vvv = sqrt(  w(1)**2 +   w(2)**2 +   w(3)**2)
w   = w/vvv
ucij(3,:) =  -w * 0.0070000D0
!
call xx_to_cart(ucij, cr_emat, xx)   ! Transform cartesian UCij to crystal base xx 
!
nsym = 0
!
do j=1, 2
   do i=1, wyc_n
      mat = wyc_mat(1:3,1:3,i)
      matt = transpose(mat)
      xxp = matmul(mat,matmul(xx,matt))
      xx = (xx + xxp)*0.5000_PREC_DP
      vvv = (xx(1,2) + xx(2,1))*0.5000_PREC_DP
      xx(1,2) = vvv
      xx(2,1) = vvv
      vvv = (xx(1,3) + xx(3,1))*0.5000_PREC_DP
      xx(1,3) = vvv
      xx(3,1) = vvv
      vvv = (xx(2,3) + xx(3,2))*0.5000_PREC_DP
      xx(2,3) = vvv
      xx(3,2) = vvv
   if(j==1) then
!      write( *,'(a,3f7.3,2x, ''>'',a3,''<'')') ' Wyckoff Prin ', wyc_axis(:, i), wyc_char(i)(1:3)
      if(index(wyc_char(i)(1:3),' 6P') >0) nsym( 6) = nsym( 6) + 1
      if(index(wyc_char(i)(1:3),'-6P') >0) nsym(-6) = nsym(-6) + 1
      if(index(wyc_char(i)(1:3),' 4P') >0) nsym( 4) = nsym( 4) + 1
      if(index(wyc_char(i)(1:3),'-4P') >0) nsym(-4) = nsym(-4) + 1
      if(index(wyc_char(i)(1:3),' 3P') >0) nsym( 3) = nsym( 3) + 1
      if(index(wyc_char(i)(1:3),'-3P') >0) nsym(-3) = nsym(-3) + 1
      if(index(wyc_char(i)(1:3),' 2 ') >0) nsym( 2) = nsym( 2) + 1
      if(index(wyc_char(i)(1:3),' m ') >0) nsym( 5) = nsym( 5) + 1
      if(index(wyc_char(i)(1:3),' 1 ') >0) nsym( 1) = nsym( 1) + 1
      if(index(wyc_char(i)(1:3),'-1 ') >0) nsym(-1) = nsym(-1) + 1
   endif
   enddo
enddo
!
line = ' '
rules = ' '
if(loutput) write(output_io,'(a)') ' '
   if(    nsym(6) >  0) then  ! Hexagonal  6..
      call do_anis_wyc_6 (nsym, loutput, xx, line)
   elseif(nsym(4) >  2 .or. nsym(-4) >  2) then  ! Cubic      4..
      call do_anis_wyc_4c(nsym, loutput, xx, line)
   elseif(nsym(4) == 1 .or. nsym(-4) == 1) then  ! tetragonal 4..
      call do_anis_wyc_4t(nsym, loutput, xx, line)
   elseif(nsym(3) >  2 .or. nsym(-3) >  2) then  ! Cubic      .3.
      call do_anis_wyc_3c(nsym, loutput, xx, line)
   elseif(nsym(3) == 1 .or. nsym(-3) == 1) then  ! trigonal   3. 
      call do_anis_wyc_3t(nsym, loutput, xx, line)
   elseif(nsym(2) == 1 .and. nsym(5) == 0) then  ! monoclinic 2.. 
      call do_anis_wyc_mmm(nsym, loutput, xx, line)
   elseif(nsym(2) == 0 .and. nsym(5) == 1) then  ! monoclinic m.. 
      call do_anis_wyc_mmm(nsym, loutput, xx, line)
   elseif(nsym(2) == 1 .and. nsym(5) == 1) then  ! monoclinic m.. 
      call do_anis_wyc_mmm(nsym, loutput, xx, line)
   elseif((nsym(2)+nsym(5)) >= 3) then           ! orthorhombic mmm or 2mm or 222
      call do_anis_wyc_mmm(nsym, loutput, xx, line)
   elseif(nsym(1) == 1 .or. nsym(-1) == 1) then  ! triclinic
      call do_anis_wyc_1(nsym, loutput, xx, line)
   else                                          ! Undetected ???
      write(output_io,'(a)') ' Unknown point group '
   endif
!endif
call xx_to_cart(xx, cr_eimat, ucij)  ! Transform XX to cartesian basis
call uij_to_xx(xx, ar_inv, uij)      ! Transform Uij from XX tensor
call calc_prin_3x3(1    , 1       , ucij, 1, prin)
call build_rules(uij, line, rules, irule)
!
write(output_io,'(  a)')      line(1)(1:len_trim(line(1)))
write(output_io,'(a,a)')  ' Constraints:  ',line(2)(1:len_trim(line(2)))
write(output_io,'(a)') ' Example:'
write(output_io,'(5x,a)') 'U11       U22       U33       U23       U13       U12'
write(output_io,'(a, 6f10.6)') '  ', uij(1,1), uij(2,2), uij(3,3), uij(2,3), uij(1,3), uij(1,2)
res_para(0) = 16
res_para(1) = uij(1,1)
res_para(2) = uij(2,2)
res_para(3) = uij(3,3)
res_para(4) = uij(2,3)
res_para(5) = uij(1,3)
res_para(6) = uij(1,2)
do i=1, 10
   res_para(6+i) = real(irule(i), kind=PREC_DP)
enddo
!
end subroutine do_anis_wyc
!
!*******************************************************************************
!
subroutine do_anis_wyc_6 (nsym, loutput, xx, line)
!-
!  Determine Eigenvectors and possible Eigenvalues for a 6.. point group
!-
!
use crystal_mod
!
use precision_mod
use prompt_mod
!
implicit none
!
integer           , dimension(-6:6), intent(in)  :: nsym   ! Number or axes with this symmetry ( 1, 2, 3, 4, m, 6 )
logical                            , intent(in) :: loutput ! Output to screen?
real(kind=PREC_DP), dimension(3,3) , intent(out) :: xx     ! Principal vectors
character(len=PREC_STRING), dimension(2) , intent(inout) :: line
!
real(kind=PREC_DP), dimension(3,3) :: mat    ! Principal vectors
real(kind=PREC_DP)                 :: xaver  ! 
!
mat = 0.0D0
mat(1,1) =  1.0D0
mat(2,2) =  1.0D0
mat(3,3) = -1.0D0
call do_anis_wyc_sym(mat, xx)
!
xaver = (abs(xx(1,1)) + abs(xx(2,2)) )/2.0_PREC_DP
xx(1,1) = sign(xaver,xx(1,1))
xx(2,2) = xx(1,1)
xx(1,2) = 0.5_PREC_DP*xx(1,1)
xx(2,1) = 0.5_PREC_DP*xx(1,1)
xx(1,3) = 0.0_PREC_DP
xx(3,1) = 0.0_PREC_DP
xx(3,2) = 0.0_PREC_DP
xx(2,3) = 0.0_PREC_DP
!
line(1) =               ' Hexagonal point group 6..'
!
end subroutine do_anis_wyc_6
!
!*******************************************************************************
!
subroutine do_anis_wyc_4c(nsym, loutput, xx, line)
!-
!  Determine Eigenvectors and possible Eigenvalues for a 4.. cubic point group
!-
!
use precision_mod
use prompt_mod
!
implicit none
!
integer           , dimension(-6:6), intent(in)  :: nsym   ! Number or axes with this symmetry ( 1, 2, 3, 4, m, 6 )
logical                            , intent(in) :: loutput ! Output to screen?
real(kind=PREC_DP), dimension(3,3) , intent(out) :: xx     ! Principal vectors
character(len=PREC_STRING), dimension(2) , intent(inout) :: line
!
real(kind=PREC_DP) :: xaver
!
xaver = (abs(xx(1,1)) + abs(xx(2,2)) + abs(xx(3,3)))/3.0_PREC_DP
xx(1,1) = sign(xaver,xx(1,1))
xx(2,2) = sign(xaver,xx(2,2))
xx(3,3) = sign(xaver,xx(3,3))
xx(1,2) = 0.0_PREC_DP
xx(2,1) = 0.0_PREC_DP
xx(1,3) = 0.0_PREC_DP
xx(3,1) = 0.0_PREC_DP
xx(3,2) = 0.0_PREC_DP
xx(2,3) = 0.0_PREC_DP
!
line(1) =               ' Cubic point group 4.. => isotropic'
!
end subroutine do_anis_wyc_4c
!
!*******************************************************************************
!
subroutine do_anis_wyc_3c(nsym, loutput, xx, line)
!-
!  Determine Eigenvectors and possible Eigenvalues for a .3. cubic point group
!  Fine tune displacement tensor to adhere to local symmetry
!-
!
use precision_mod
use prompt_mod
!
implicit none
!
integer           , dimension(-6:6), intent(in)  :: nsym   ! Number or axes with this symmetry ( 1, 2, 3, 4, m, 6 )
logical                            , intent(in)  :: loutput ! Output to screen?
real(kind=PREC_DP), dimension(3,3) , intent(out) :: xx     ! Displacement tensor 
character(len=PREC_STRING), dimension(2) , intent(inout) :: line
!
real(kind=PREC_DP) :: xaver
!
xaver = (abs(xx(1,1)) + abs(xx(2,2)) + abs(xx(3,3)))/3.0_PREC_DP
xx(1,1) = sign(xaver,xx(1,1))
xx(2,2) = sign(xaver,xx(2,2))
xx(3,3) = sign(xaver,xx(3,3))
!
xx(1,1) = 0.0_PREC_DP
xx(2,2) = 0.0_PREC_DP
xx(3,3) = 0.0_PREC_DP
!
line(1) =               ' Cubic point group .3. => isotropic'
!
end subroutine do_anis_wyc_3c
!
!*******************************************************************************
!
subroutine do_anis_wyc_4t(nsym, loutput, xx, line  )
!-
!  Determine Eigenvectors and possible Eigenvalues for a 4.. point group
!  Needs to consider tetragonal and cubic cases
!-
use crystal_mod
use wyckoff_mod
!
use lib_metric_mod
use precision_mod
use prompt_mod
!
implicit none
!
integer           , dimension(-6:6), intent(in)  :: nsym   ! Number or axes with this symmetry ( 1, 2, 3, 4, m, 6 )
logical                            , intent(in) :: loutput ! Output to screen?
real(kind=PREC_DP), dimension(3,3) , intent(out) :: xx     ! Principal vectors
character(len=PREC_STRING), dimension(2) , intent(inout) :: line
!
integer :: i  ! Dummy index
real(kind=PREC_DP)               :: xaver
!
loop_search: do i=1, wyc_n
   if(index(wyc_char(i)(1:3),'4') >0) exit loop_search
enddo loop_search
!
xx(2,3) = 0.0_PREC_DP
xx(3,2) = 0.0_PREC_DP
xx(1,3) = 0.0_PREC_DP
xx(3,1) = 0.0_PREC_DP
xx(1,2) = 0.0_PREC_DP
xx(2,1) = 0.0_PREC_DP
!
if(abs(lib_bang(cr_gten, wyc_axis(:,i), (/ 0.0_PREC_DP, 0.0_PREC_DP, 1.0_PREC_DP /))) <   5.0_PREC_DP .or. &
   abs(lib_bang(cr_gten, wyc_axis(:,i), (/ 0.0_PREC_DP, 0.0_PREC_DP, 1.0_PREC_DP /))) > 175.0_PREC_DP ) then    ! parallel [0 0 1] 
   xaver = (abs(xx(1,1)) + abs(xx(2,2)) )/3.0_PREC_DP
   xx(1,1) = xaver
   xx(2,2) = xaver
   line(1) =               ' Tetragonal point group 4.. parallel [ 0 0 1 ]'
elseif(abs(lib_bang(cr_gten, wyc_axis(:,i), (/ 0.0_PREC_DP, 1.0_PREC_DP, 0.0_PREC_DP /))) <   5.0_PREC_DP .or. &
       abs(lib_bang(cr_gten, wyc_axis(:,i), (/ 0.0_PREC_DP, 1.0_PREC_DP, 0.0_PREC_DP /))) > 175.0_PREC_DP ) then    ! parallel [0 0 1] 
   xaver = (abs(xx(1,1)) + abs(xx(3,3)) )/3.0_PREC_DP
   xx(1,1) = xaver
   xx(3,3) = xaver
   line(1) =               ' Tetragonal point group 4.. parallel [ 0 1 0 ]'
elseif(abs(lib_bang(cr_gten, wyc_axis(:,i), (/ 1.0_PREC_DP, 0.0_PREC_DP, 0.0_PREC_DP /))) <   5.0_PREC_DP .or. &
       abs(lib_bang(cr_gten, wyc_axis(:,i), (/ 1.0_PREC_DP, 0.0_PREC_DP, 0.0_PREC_DP /))) > 175.0_PREC_DP ) then    ! parallel [1 0 0] 
   xaver = (abs(xx(2,2)) + abs(xx(3,3)) )/3.0_PREC_DP
   xx(2,2) = xaver
   xx(3,3) = xaver
   line(1) =               ' Tetragonal point group 4.. parallel [ 1 0 0 ]'
endif
!
end subroutine do_anis_wyc_4t
!
!*******************************************************************************
!
subroutine do_anis_wyc_3t(nsym, loutput, xx, line)
!-
!  Fine tune xx for a 3.. point group
!  Needs to consider trigonal [001]   and cubic <111> cases
!-
use crystal_mod
use wyckoff_mod
!
use lib_metric_mod
use precision_mod
use prompt_mod
!
implicit none
!
integer           , dimension(-6:6), intent(in)  :: nsym   ! Number or axes with this symmetry ( 1, 2, 3, 4, m, 6 )
logical                            , intent(in) :: loutput ! Output to screen?
real(kind=PREC_DP), dimension(3,3) , intent(out) :: xx   ! Principal vectors
character(len=PREC_STRING), dimension(2) , intent(inout) :: line
!
integer :: i   ! Dummy index
real(kind=PREC_DP), dimension(3,3) :: mat    ! Principal vectors
real(kind=PREC_DP)               :: xaver
!
!
loop_search: do i=1, wyc_n
   if(index(wyc_char(i)(1:3),'3') >0) exit loop_search
enddo loop_search
!
! TRIGONAL case
!
if(abs(lib_bang(cr_gten, wyc_axis(:,i), (/ 0.0_PREC_DP, 0.0_PREC_DP, 1.0_PREC_DP /))) <   5.0_PREC_DP .or. &
   abs(lib_bang(cr_gten, wyc_axis(:,i), (/ 0.0_PREC_DP, 0.0_PREC_DP, 1.0_PREC_DP /))) > 175.0_PREC_DP ) then    ! parallel [0 0 1] TRIGONAL
!
   mat = 0.0D0
   mat(1,1) =  1.0D0
   mat(2,2) =  1.0D0
   mat(3,3) = -1.0D0     ! Mirror normal to c
   call do_anis_wyc_sym(mat, xx)
!
!
   mat = 0.0D0
   mat(1,1) = -1.0D0     ! 2-fold around c
   mat(2,2) = -1.0D0
   mat(3,3) =  1.0D0
   call do_anis_wyc_sym(mat, xx)
!
   line(1) =               ' Trigonal point group 3.'
!
!  CUBIC case
!
else                     ! Cubic case
   xaver = (abs(xx(1,1)) + abs(xx(2,2)) + abs(xx(3,3)))/3.0_PREC_DP
   xx(1,1) = xaver
   xx(2,2) = xaver
   xx(3,3) = xaver
   if(abs(lib_bang(cr_gten, wyc_axis(:,i), (/ 1.0_PREC_DP, 1.0_PREC_DP, 1.0_PREC_DP /))) <   5.0_PREC_DP .or. &
      abs(lib_bang(cr_gten, wyc_axis(:,i), (/ 1.0_PREC_DP, 1.0_PREC_DP, 1.0_PREC_DP /))) > 175.0_PREC_DP ) then    ! parallel [1 1 1] TRIGONAL
      xaver = (abs(xx(2,3)) + abs(xx(3,2)) +                &
               abs(xx(1,3)) + abs(xx(3,1)) +                &
               abs(xx(1,2)) + abs(xx(2,1))   )/3.0_PREC_DP
      xx(2,3) = xaver
      xx(3,2) = xaver
      xx(1,3) = xaver
      xx(3,1) = xaver
      xx(1,2) = xaver
      xx(2,1) = xaver
!
      line(1) =               ' Cubic point group .3. parallel [  1  1  1 ]'
   elseif(abs(lib_bang(cr_gten, wyc_axis(:,i), (/-1.0_PREC_DP, 1.0_PREC_DP, 1.0_PREC_DP /))) <   5.0_PREC_DP .or. &
          abs(lib_bang(cr_gten, wyc_axis(:,i), (/-1.0_PREC_DP, 1.0_PREC_DP, 1.0_PREC_DP /))) > 175.0_PREC_DP ) then ! parallel [-1 1 1] TRIGONAL
      xaver = (abs(xx(2,3)) + abs(xx(3,2)) +                &
               abs(xx(1,3)) + abs(xx(3,1)) +                &
               abs(xx(1,2)) + abs(xx(2,1))   )/3.0_PREC_DP
      xx(2,3) = sign(xaver, xx(2,3))
      xx(3,2) = sign(xaver, xx(2,3))
      xx(1,3) = sign(xaver, xx(1,3))
      xx(3,1) = sign(xaver, xx(1,3))
      xx(1,2) = sign(xaver, xx(1,2))
      xx(2,1) = sign(xaver, xx(1,2))
!
      line(1) =               ' Cubic point group .3. parallel [ -1  1  1 ]'
   elseif(abs(lib_bang(cr_gten, wyc_axis(:,i), (/-1.0_PREC_DP,-1.0_PREC_DP, 1.0_PREC_DP /))) <   5.0_PREC_DP .or. &
          abs(lib_bang(cr_gten, wyc_axis(:,i), (/-1.0_PREC_DP,-1.0_PREC_DP, 1.0_PREC_DP /))) > 175.0_PREC_DP ) then ! parallel [-1 -1 1] TRIGONAL
      xaver = (abs(xx(2,3)) + abs(xx(3,2)) +                &
               abs(xx(1,3)) + abs(xx(3,1)) +                &
               abs(xx(1,2)) + abs(xx(2,1))   )/3.0_PREC_DP
      xx(2,3) = sign(xaver, xx(2,3))
      xx(3,2) = sign(xaver, xx(2,3))
      xx(1,3) = sign(xaver, xx(1,3))
      xx(3,1) = sign(xaver, xx(1,3))
      xx(1,2) = sign(xaver, xx(1,2))
      xx(2,1) = sign(xaver, xx(1,2))
!
      line(1) =               ' Cubic point group .3. parallel [ -1 -1  1 ]'
   elseif(abs(lib_bang(cr_gten, wyc_axis(:,i), (/ 1.0_PREC_DP,-1.0_PREC_DP, 1.0_PREC_DP /))) <   5.0_PREC_DP .or. &
          abs(lib_bang(cr_gten, wyc_axis(:,i), (/ 1.0_PREC_DP,-1.0_PREC_DP, 1.0_PREC_DP /))) > 175.0_PREC_DP ) then ! parallel [-1 1 1] TRIGONAL
      xaver = (abs(xx(2,3)) + abs(xx(3,2)) +                &
               abs(xx(1,3)) + abs(xx(3,1)) +                &
               abs(xx(1,2)) + abs(xx(2,1))   )/3.0_PREC_DP
      xx(2,3) = sign(xaver, xx(2,3))
      xx(3,2) = sign(xaver, xx(2,3))
      xx(1,3) = sign(xaver, xx(1,3))
      xx(3,1) = sign(xaver, xx(1,3))
      xx(1,2) = sign(xaver, xx(1,2))
      xx(2,1) = sign(xaver, xx(1,2))
!
      line(1) =               ' Cubic point group .3. parallel [  1 -1  1 ]'
   endif
endif
!
end subroutine do_anis_wyc_3t
!
!*******************************************************************************
!
subroutine do_anis_wyc_mmm(nsym, loutput, xx, line)
!-
!  Determine Eigenvectors and possible Eigenvalues for a mmm point group
!  Needs to consider tetragonal and cubic cases
!-
use crystal_mod
use wyckoff_mod
!
use blanks_mod
use lib_metric_mod
use precision_mod
use prompt_mod
!
implicit none
!
integer           , dimension(-6:6), intent(in)  :: nsym   ! Number or axes with this symmetry ( 1, 2, 3, 4, m, 6 )
logical                            , intent(in) :: loutput ! Output to screen?
real(kind=PREC_DP), dimension(3,3) , intent(out) :: xx   ! Principal vectors
character(len=PREC_STRING), dimension(2) , intent(inout) :: line
!character(len=PREC_STRING), dimension(10), intent(inout) :: rules
!
character(len=12), dimension(9)  :: system_name
data system_name / 'Triclinic   ', 'Monoclinic B', 'Monoclinic C', 'Orthorhombic',  &
                   'Tetragonal  ', 'Rhombohedral', 'Trigonal    ', 'Hexagonal   ',  &
                   'Cubic       '/
!
!loop_search: do i=1, wyc_n
!   if(index(wyc_char(i)(1:3),'2') > 0 .or. index(wyc_char(i)(1:3),'m') > 0) exit loop_search
!enddo loop_search
!
if(cr_syst==cr_monoclinicB) then        ! Monoclinic, unique [ 0 1 0 ]
   xx(3,2) = 0.0_PREC_DP
   xx(2,3) = 0.0_PREC_DP
   xx(1,2) = 0.0_PREC_DP
   xx(2,1) = 0.0_PREC_DP
   line(1) = ' Monoclinic point group .2., .2/m.      '
elseif(cr_syst==cr_monoclinicC) then        ! Monoclinic, unique [ 0 0 1 ]
   xx(1,2) = 0.0_PREC_DP
   xx(2,1) = 0.0_PREC_DP
   xx(1,3) = 0.0_PREC_DP
   xx(3,1) = 0.0_PREC_DP
   line(1) = ' Monoclinic point group ..2, ..2/m      '
elseif(cr_syst==cr_ortho      ) then        ! orthorhombic
   if(sum(nsym(2:5))==1) then
      if(nint(abs(wyc_axis(1,1)))==1) then  !  [ 1 0 0 ]
         xx(1,2) = 0.0_PREC_DP
         xx(2,1) = 0.0_PREC_DP
         xx(1,3) = 0.0_PREC_DP
         xx(3,1) = 0.0_PREC_DP
         line(1) = ' Orthorhombic point group m.., 2.., or 2/m..'
      elseif(nint(abs(wyc_axis(2,1)))==1) then  !  [ 0 1 0 ]
         xx(1,2) = 0.0_PREC_DP
         xx(2,1) = 0.0_PREC_DP
         xx(2,3) = 0.0_PREC_DP
         xx(3,2) = 0.0_PREC_DP
         line(1) = ' Orthorhombic point group .m., .2., or .2/m.'
      elseif(nint(abs(wyc_axis(3,1)))==1) then  !  [ 0 0 1 ]
         xx(1,3) = 0.0_PREC_DP
         xx(3,1) = 0.0_PREC_DP
         xx(2,3) = 0.0_PREC_DP
         xx(3,2) = 0.0_PREC_DP
         line(1) = ' Orthorhombic point group ..m, ..2, or ..2/m'
      endif
   elseif(sum(nsym(2:5))>1) then
      xx(2,3) = 0.0_PREC_DP
      xx(3,2) = 0.0_PREC_DP
      xx(1,2) = 0.0_PREC_DP
      xx(2,1) = 0.0_PREC_DP
      xx(1,3) = 0.0_PREC_DP
      xx(3,1) = 0.0_PREC_DP
      line(1) = ' Orthorhombic point group mmm, 222, or 2mm'  
   endif
elseif(cr_syst==cr_tetragonal .or. cr_syst==cr_cubic) then        ! tetragonal / cubic
   line = ' '
   if(nsym(2)==3 .and. nsym(5)==3) then     !  mmm
      line(1) = system_name(cr_syst)(1:len_trim(system_name(cr_syst))) // ' ' // &
                'point group mmm, 222, or 2mm'
   elseif(nsym(2)==1 .and. nsym(5)==2) then     !  2mm
      line(1) = system_name(cr_syst)(1:len_trim(system_name(cr_syst))) // ' ' // &
                'point group 2mm'
   elseif(nsym(2)==1 .and. nsym(5)==1) then     !  2/m
      line(1) = system_name(cr_syst)(1:len_trim(system_name(cr_syst))) // ' ' // &
                'point group 2/m'
   elseif(nsym(2)==1 .and. nsym(5)==0) then     !  2..
      line(1) = system_name(cr_syst)(1:len_trim(system_name(cr_syst))) // ' ' // &
                'point group 2..'
   elseif(nsym(2)==0 .and. nsym(5)==2) then     !  m..
      line(1) = system_name(cr_syst)(1:len_trim(system_name(cr_syst))) // ' ' // &
                'point group m..'
   endif
endif
!
!
end subroutine do_anis_wyc_mmm
!
!*******************************************************************************
!
subroutine do_anis_wyc_1 (nsym, loutput, xx, line)
!-
!  Determine Eigenvectors and possible Eigenvalues for a 1.. point group
!-
!
use lib_metric_mod
use precision_mod
use prompt_mod
!
implicit none
!
integer           , dimension(-6:6), intent(in)  :: nsym   ! Number or axes with this symmetry ( 1, 2, 3, 4, m, 6 )
logical                            , intent(in) :: loutput ! Output to screen?
real(kind=PREC_DP), dimension(3,3) , intent(out) :: xx   ! Principal vectors
character(len=PREC_STRING), dimension(2) , intent(inout) :: line
!
!
line(1) = ' Triclinic point group 1..'
!
end subroutine do_anis_wyc_1
!
!*******************************************************************************
!
subroutine do_anis_wyc_sym(mat, xx)
!-
!   Apply a symmetry operation to the Displacement tensor xxx, 
!   ensure matrix is symmetric
!+
use precision_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3,3) , intent(out) :: xx   ! Principal vectors
!
real(kind=PREC_DP), dimension(3,3) :: mat    ! Principal vectors
real(kind=PREC_DP), dimension(3,3) :: xxp    ! Principal vectors
real(kind=PREC_DP)                 :: vvv
!
xxp = matmul(mat,matmul(xx,mat))
xx = (xx + xxp)*0.5000_PREC_DP
vvv = (xx(1,2) + xx(2,1))*0.5000_PREC_DP
xx(1,2) = vvv
xx(2,1) = vvv
vvv = (xx(1,3) + xx(3,1))*0.5000_PREC_DP
xx(1,3) = vvv
xx(3,1) = vvv
vvv = (xx(2,3) + xx(3,2))*0.5000_PREC_DP
xx(2,3) = vvv
xx(3,2) = vvv
!
end subroutine do_anis_wyc_sym
!
!*******************************************************************************
!
subroutine build_rules(uij, line, rules, irule)
!-
!  Build Rules for Uij
!+
use crystal_mod
!
use blanks_mod
use precision_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3,3), intent(inout) :: uij
character(len=*), dimension(2) , intent(inout) :: line
character(len=*), dimension(10), intent(out)   :: rules
integer         , dimension(10), intent(out)   :: irule
!
real(kind=PREC_DP), parameter :: TOL = 1.0D-7
real(kind=PREC_DP), parameter :: TOLH= 1.0D-2
integer :: i
!
rules = ' '
irule = 10
!  Rules   U11 == U22
if    (abs(abs(uij(1,1))-abs(uij(2,2)))<TOL) then   !  U11 == U22
   rules(1) = 'U11 = U22'
   irule(1) =  1
elseif(abs(abs(uij(1,1))+abs(uij(2,2)))<TOL) then   !  U11 == -U22
   rules(1) = 'U11 = -U22'
   irule(1) = -1
endif
!  Rules   U11 == U33
if    (abs(abs(uij(1,1))-abs(uij(3,3)))<TOL) then   !  U11 == U22
   rules(2) = 'U11 = U33'
   irule(2) =  1
elseif(abs(abs(uij(1,1))+abs(uij(3,3)))<TOL) then   !  U11 == -U22
   rules(2) = 'U11 = -U33'
   irule(2) = -1
endif
!  Rules   U22 == U33
if    (abs(abs(uij(2,2))-abs(uij(3,3)))<TOL) then   !  U11 == U22
   rules(3) = 'U22 = U33'
   irule(3) =  1
elseif(abs(abs(uij(1,1))+abs(uij(3,3)))<TOL) then   !  U11 == -U22
   rules(3) = 'U22 = -U33'
   irule(3) = -1
endif
!  Rule  U23 == 0
if    (abs(uij(2,3))<TOL) then                     ! U23 == 0
   rules(4) = 'U23 = 0'
   irule(4) =  0
endif
!  Rule  U13 == 0
if    (abs(uij(1,3))<TOL) then                     ! U13 == 0
   rules(5) = 'U13 = 0'
   irule(5) =  0
endif
!  Rule  U12 == 0
if    (abs(uij(1,2))<TOL) then                     ! U12 == 0
   rules(6) = 'U12 = 0'
   irule(6) =  0
endif
!  Rules   U23 == U13
if    (abs(abs(uij(2,3))-abs(uij(1,3)))<TOL) then   !  U23 == U13
   rules(7) = 'U23 = U13'
   irule(7) =  1
elseif(abs(abs(uij(2,3))+abs(uij(1,3)))<TOL) then   !  U23 == -U13
   rules(7) = 'U23 = -U13'
   irule(7) = -1
endif
!  Rules   U23 == U12
if    (abs(abs(uij(2,3))-abs(uij(1,2)))<TOL) then   !  U23 == U12
   rules(8) = 'U23 = U12'
   irule(8) =  1
elseif(abs(abs(uij(2,3))+abs(uij(1,2)))<TOL) then   !  U23 == -U12
   rules(8) = 'U23 = -U12'
   irule(8) = -1
endif
!  Rules   U13 == U12
if    (abs(abs(uij(1,3))-abs(uij(1,2)))<TOL) then   !  U13 == U12
   rules(9) = 'U13 = U12'
   irule(9) =  1
elseif(abs(abs(uij(1,3))+abs(uij(1,2)))<TOL) then   !  U13 == -U12
   rules(9) = 'U13 = -U12'
   irule(9) = -1
endif
!
if(cr_syst==cr_trigonal .or. cr_syst==cr_hexagonal) then
   if(abs(abs(uij(1,1))-2.0_PREC_DP*abs(uij(1,2)))<TOLH) then   ! U11 = 2.*U12
      uij(1,2) = sign(0.50_PREC_DP*uij(1,1), uij(1,2))
      uij(2,1) = uij(1,2)
      rules(10) = 'U11 = 2.*U12'
      irule(10) =  2
   endif
endif
!
   do i=1, 10
      line(2) = line(2)(1:len_trim(line(2))) //'    ' // rules(i)(1:len_trim(rules(i)))
   enddo
   i = len_trim(line(2))
call rem_leading_bl(line(2),i)
!
end subroutine build_rules
!
!*******************************************************************************
!
end module prep_anis_mod
