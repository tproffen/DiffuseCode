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
!public calc_prin
!
interface lookup_anis
   module procedure lookup_anis_6, lookup_anis_3x3, lookup_anis_6_prin,         &
                    lookup_anis_1
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
subroutine prep_anis(natom)
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
!
implicit none
!
integer, intent(in) :: natom  ! Number of atoms to cycle 
!
integer :: itype    ! Loop index atom types
integer :: iatom    ! Loop index atoms
integer :: iref     ! Reference atom, symmetrically equivalent atoms take this atom as model
integer :: j, l     ! Dummy loop indices
logical :: lanis       ! This atom has anisotropic ADP
logical :: lsuccess    ! Success in subroutines
real(kind=PREC_DP), dimension(3)   :: ar_inv  ! ( 1/a*, 1/b*, 1/c*)
real(kind=PREC_DP), dimension(3,3) :: imat  ! Unit matrix
real(kind=PREC_DP), dimension(3,3) :: uij   ! U_ij as in SHELX, reference to ((a*).a, (b*).b, (c*))c
real(kind=PREC_DP), dimension(3,3) :: xx    ! Displacement tensor <DELTAx^i . DELTAx^j> refers to standard basis
real(kind=PREC_DP), dimension(3,3) :: ucij  ! Displacement tensor <DELTAx^i . DELTAx^j> at cartesian coordinates
real(kind=PREC_DP), dimension(3,3) :: symm_mat    ! Symmetry matrix
!
data imat / 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0/
!
ar_inv(1) = 1./cr_ar(1)
ar_inv(2) = 1./cr_ar(2)
ar_inv(3) = 1./cr_ar(3)
!
if(allocated(cr_anis_full)) deallocate(cr_anis_full)
if(allocated(cr_prin     )) deallocate(cr_prin     )
allocate(cr_anis_full(6,    cr_ncatoms))
allocate(cr_prin     (4, 3, cr_ncatoms))
cr_anis_full = 0.0D0
cr_prin      = 0.0D0
uij          = 0.0D0
ucij         = 0.0D0
xx           = 0.0D0
!
! Step 0:  Determine number of different atoms for a given scattering type
!
j = 0
lsuccess= .false.
cr_ndiffer = 0
cr_nanis = 0
cr_is_anis = .false.
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
loop_atoms:do iatom=1,natom 
   j = 0
   lsuccess= .false.
   itype = cr_iscat(1,iatom)
   cr_ianis(iatom) = iatom
!  write(*,*) ' NANIS ', cr_nanis
!  write(*,'(a,i3,2i3, 3f7.3, i4)') ' ATOM ', iatom, itype, cr_iscat(3,iatom), cr_pos(:,iatom), cr_is_sym(iatom)
!  write(*,'(a,        6f7.3    )') ' UANI ', cr_anis(:,itype)
   if(maxval(abs(cr_anis(2:,itype)))> 0.0D0) then     ! Elements 1 to 6 are given Anisotropic atom type
      lanis =.true.
      uij(1,1) = cr_anis(1, itype)
      uij(2,2) = cr_anis(2, itype)
      uij(3,3) = cr_anis(3, itype)
      uij(2,3) = cr_anis(4, itype)
      uij(3,2) = cr_anis(4, itype)
      uij(1,3) = cr_anis(5, itype)
      uij(3,1) = cr_anis(5, itype)
      uij(1,2) = cr_anis(6, itype)
      uij(2,1) = cr_anis(6, itype)
!     cr_iscat(3,iatom) = cr_is_sym(iatom)
      cr_iscat(3,iatom) = cr_iscat(2,iatom)
      call uij_to_xx(uij, cr_ar, xx)      ! Transform Uij to XX tensor
      call xx_to_cart(xx, cr_eimat, ucij) ! Transform XX to cartesian basis
      cr_is_anis = .true.
   else                       ! Isotropic atom set cartesian UCij and transform back to UIJ
      lanis =.false.
      ucij      = 0.0D0
      if(cr_anis(1,itype)> 0.0D0) then                ! Element U11 is the only one take as Uiso
         ucij(1,1) = cr_anis(1,itype)
         ucij(2,2) = cr_anis(1,itype)
         ucij(3,3) = cr_anis(1,itype)
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
!
!  if(cr_is_sym(iatom)>1) then             ! Symmetrically equivalent atom
!if(cr_iscat(2,iatom)==1) then
!write(*,'(a,i3, a,6f10.5, 3i4, l2)') ' UIJ ', 1, ' | ', cr_anis(:,itype), cr_iscat(:, iatom), lanis
!endif
   if(cr_iscat(2,iatom)>1) then            ! Symmetrically equivalent atom
      if(lanis) then                       ! Anisotropic ADP, apply symmetry to obtain new ADP, principal vectors
!        symm_mat = spc_mat(1:3,1:3, cr_is_sym(iatom))
         symm_mat = spc_mat(1:3,1:3, cr_iscat(2,iatom))
!write(*,'(i3)') cr_is_sym(iatom)
!write(*,'(3f7.2)') spc_mat(1,1:3, cr_is_sym(iatom))
!write(*,'(3f7.2)') spc_mat(2,1:3, cr_is_sym(iatom))
!write(*,'(3f7.2)') spc_mat(3,1:3, cr_is_sym(iatom))
         call anis_symm(cr_ncatoms, cr_nanis, iref , symm_mat, ubound(cr_anis_full,2), &
              cr_anis_full, cr_prin, xx, uij, cr_emat, cr_eimat, ar_inv)
!write(*,'(a,i3, a,6f10.5)') ' uij ', cr_iscat(2,iatom), ' | ', uij(1,1), uij(2,2), uij(3,3), uij(2,3), uij(1,3), uij(1,2)
         call lookup_anis(cr_nanis, cr_anis_full, cr_prin, uij, ucij, .FALSE., j, lsuccess)   ! Check if old ADPs are reproduced
!write(*,'(a,i3, a,4f10.5,3i6)') ' prin', 1, ' | ', cr_prin(:,1,j), j, cr_nanis, ubound(cr_prin,3)
!write(*,'(a,i3, a,6f10.5)') ' prin', 2, ' | ', cr_prin(:,2,j)
!write(*,'(a,i3, a,6f10.5)') ' prin', 3, ' | ', cr_prin(:,3,j)
         cr_iscat(3,iatom) = j            ! Assign ADP type
         cycle loop_atoms
      endif 
!  else
!     iref = iatom    ! Reference atom for ADP values 
   endif 
!
   call lookup_anis(cr_nanis, cr_anis_full, cr_prin, uij, ucij, .FALSE., j, lsuccess) ! Check if previous ADPs exist with identical uij
   cr_iscat(3,iatom) = j            ! Assign ADP type
   iref = j           ! Reference atom for ADP values
   if(lsuccess) then
      cycle loop_atoms
   endif
!  This atom has new UIJ, we need to calculate the principal vectors and Eigenvalues
   call calc_prin(cr_iscat(3,iatom), cr_nanis, ucij, ubound(cr_prin,3), cr_prin)
!write(*,'(a,i3, a,4f10.5,i6)') ' PRIN', 1, ' | ', cr_prin(:,1,cr_nanis), iref
!write(*,'(a,i3, a,6f10.5)') ' PRIN', 2, ' | ', cr_prin(:,2,cr_nanis)
!write(*,'(a,i3, a,6f10.5)') ' PRIN', 3, ' | ', cr_prin(:,3,cr_nanis)
enddo loop_atoms
!do j=1, cr_nanis
!write(*,'(a, 6f9.6)') 'PR 2 ', cr_anis_full(1:6,j)
!enddo
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
subroutine lookup_anis_6_prin(cr_nanis, cr_anis_full, cr_prin, uij_6, prin, ientry, lsuccess)
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
call lookup_anis_3x3(cr_nanis, cr_anis_full, cr_prin, uij, ucij, lcopy, ientry, lsuccess)
cr_prin(:,:,ientry) = prin(:,:)
!
end subroutine lookup_anis_6_prin
!
!*******************************************************************************
!
subroutine lookup_anis_6(cr_nanis, cr_anis_full, cr_prin, uij_6, cr_eimat, cr_ar, &
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
call lookup_anis_3x3(cr_nanis, cr_anis_full, cr_prin, uij, ucij, lcopy, ientry, lsuccess)
!
end subroutine lookup_anis_6
!
!*******************************************************************************
!
subroutine lookup_anis_1(cr_nanis, cr_anis_full, cr_prin, cr_dw, cr_emat, cr_ar, &
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
call xx_to_cart(ucij, cr_emat, xx)   ! Transform to XX from cartesian basis
call uij_to_xx(xx, ar_inv, uij)      ! Transform to Uij from XX tensor
!
lcopy = .true.
!
call lookup_anis_3x3(cr_nanis, cr_anis_full,  cr_prin, uij, ucij, lcopy, ientry, lsuccess)
!
end subroutine lookup_anis_1
!
!*******************************************************************************
!
subroutine lookup_anis_3x3(cr_nanis, cr_anis_full,  cr_prin, uij, ucij, lcopy, ientry, lsuccess)
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
!
!write(*,*) ' SUCCESS ', lsuccess, cr_nanis, lcopy
if(ientry==-1) then   ! new entry
   cr_nanis = cr_nanis + 1
   ientry   = cr_nanis
   if(cr_nanis>=ubound(cr_anis_full,2)) then
      call alloc_arr(cr_anis_full, 1,6, 1, cr_nanis, all_status, 0.0D0      )
      call alloc_arr(cr_prin, 1, 4,1,3, 1, cr_nanis, all_status, 0.0D0      )
   endif
   cr_anis_full(1,ientry) = uij(1,1)
   cr_anis_full(2,ientry) = uij(2,2)
   cr_anis_full(3,ientry) = uij(3,3)
   cr_anis_full(4,ientry) = uij(3,2)
   cr_anis_full(5,ientry) = uij(3,1)
   cr_anis_full(6,ientry) = uij(2,1)
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
subroutine xx_to_cart(xx, emat, ucij) ! Transform XX to cartesian basis
!-
! A simple 1 line routine that transforms XX to cartesian space or back
!  ucij = matmul((eimat), matmul(xx,   transpose(eimat)))  ! With eimat into cartesian
!  xx   = matmul((emat) , matmul(ucij, transpose(emat)))   ! With emat into crystal space
!+
!
use precision_mod
!implicit none
!
real(kind=PREC_DP), dimension(3,3), intent(in)  :: xx   ! Input matrix
real(kind=PREC_DP), dimension(3,3), intent(in)  :: emat ! transformation matrix
real(kind=PREC_DP), dimension(3,3), intent(out) :: ucij ! Result matrix
!
ucij = matmul((emat), matmul(xx, transpose(emat)))
!
end subroutine xx_to_cart
!
!*******************************************************************************
!
subroutine anis_symm(ncatoms, nanis, iref, symm_mat, idim1, anis_full, &
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
real(kind=PREC_DP), dimension(3,3)         , intent(in)    :: xx        ! Displacement tensor crystal space
real(kind=PREC_DP), dimension(3,3)         , intent(in)    :: emat      ! Transformation basis to Cartesian
real(kind=PREC_DP), dimension(3,3)         , intent(in)    :: eimat     ! Transformation basis to crystal
real(kind=PREC_DP), dimension(3)           , intent(in)    :: ar_inv    ! Reciprocal lattice parameters
!
integer                          :: j
real(kind=PREC_DP), dimension(3) :: v,w  ! Dummy vector
real(kind=PREC_DP), dimension(3,3) :: xx_new       ! Displacemetn tensor after symmetry operation
!real(kind=PREC_DP), dimension(3,3) :: uij          ! UIJ in (a*.a, b*.b, c*.c) space
real(kind=PREC_DP), dimension(3,3) :: symm_imat    ! Inverse Symmetry matrix
!
call matinv(symm_mat, symm_imat)
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
call uij_to_xx(xx_new, ar_inv, uij)
anis_full(1,nanis+1) = uij(1,1)
anis_full(2,nanis+1) = uij(2,2)
anis_full(3,nanis+1) = uij(3,3)
anis_full(4,nanis+1) = uij(2,3)
anis_full(5,nanis+1) = uij(1,3)
anis_full(6,nanis+1) = uij(1,2)
!
end subroutine anis_symm
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
integer, parameter :: MAXW  = 2
integer, parameter :: MAXWW = 6
character(len=PREC_STRING), dimension(2) :: cpara
integer                   , dimension(2) :: lpara
integer                                  :: ianz
character(len=PREC_STRING)               :: ccpara
integer                                  :: llpara
real(kind=PREC_DP) , dimension(6)        :: wwerte   ! Calculated values
!
integer :: ianis
real(kind=PREC_DP), dimension(3)   :: ar_inv  ! ( 1/a*, 1/b*, 1/c*)
real(kind=PREC_DP), dimension(3,3) :: uij   ! U_ij as in SHELX, reference to ((a*).a, (b*).b, (c*))c
real(kind=PREC_DP), dimension(3,3) :: xx    ! Displacement tensor <DELTAx^i . DELTAx^j> refers to standard basis
real(kind=PREC_DP), dimension(3,3) :: ucij  ! Displacement tensor <DELTAx^i . DELTAx^j> at cartesian coordinates
!
integer, parameter :: NOPTIONAL = 2
integer, parameter :: O_TYPE    = 1
integer, parameter :: O_VALUES  = 2
character(LEN=   6), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 1 ! Number of values to calculate 
!
data oname  / 'type', 'values'   /
data loname /  4,      6         /
opara  =  (/ '1.0000', '[0.01]' /)   ! Always provide fresh default values
lopara =  (/  6,        6       /)
owerte =  (/  0.0,      0.0     /)
!
call get_params (line, ianz, cpara, lpara, MAXW, length)
if(ier_num /= 0) return
!
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num /= 0) return
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
   ucij(1,1) = wwerte(1)                ! Set the caresian diagonal values
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
!
end subroutine do_anis
!
!*******************************************************************************
!
end module prep_anis_mod
