module prep_anis_mod
!-
! Prepare the anisotropic ADPs
!   Determine eigenvectors and Eigenvalues from U_ij
!   Apply symmetry to Eigenvectors for all sites
!+
private
public prep_anis
public lookup_anis
public calc_prin
!
interface lookup_anis
   module procedure lookup_anis_6, lookup_anis_3x3
end interface lookup_anis
!
interface calc_prin
   module procedure calc_prin_6, calc_prin_3x3
end interface calc_prin
!
contains
!
!*******************************************************************************
!
subroutine prep_anis
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
real(kind=PREC_DP), parameter :: TOL = 1.0D-6
integer :: itype    ! Loop index atom types
integer :: iatom    ! Loop index atoms
integer :: j, l, m, n  ! Dummy loop indices
integer :: neigen      ! Number distinct eigenvectors
logical :: lsuccess    ! Success in subroutines
real(kind=PREC_DP), dimension(3,3) :: uij   ! U_ij in crystal dimensions
real(kind=PREC_DP), dimension(3,3) :: uijs  ! U_ij in crystal dimensions after symmetry
real(kind=PREC_DP), dimension(3,3) :: ucij   ! U_ij in cartesian coordinates
real(kind=PREC_DP), dimension(3,3) :: symm_mat    ! Symmetry matrix
real(kind=PREC_DP), dimension(3,3) :: symm_imat    ! Inverse Symmetry matrix
!
!write(*,*) ' PREP_ANIS cr_anis  ', ubound(cr_anis)
!write(*,*) ' PREP_ANIS cr_anis_f', ubound(cr_anis_full)
!write(*,*) ' PREP_ANIS cr_prin  ', ubound(cr_prin,3)
!write(*,*)
!
if(allocated(cr_anis_full)) deallocate(cr_anis_full)
if(allocated(cr_prin     )) deallocate(cr_prin     )
allocate(cr_anis_full(6,    cr_ncatoms))
allocate(cr_prin     (3, 4, cr_ncatoms))
cr_anis_full = 0.0D0
cr_prin      = 0.0D0
!
! Step 0:  Determine number of different atoms for a given scattering type
!
j = 0
lsuccess= .false.
cr_ndiffer = 0
cr_nanis = 0
!if(maxval(abs(cr_anis(:,:)))<TOL) then
!   cr_iscat(1:cr_ncatoms,3) = 1
!   return
!endif
!do iatom=1,cr_ncatoms
!write(*,'( i3, 6f10.6)') iatom, cr_anis(:,iatom)
!enddo
do iatom=1,cr_ncatoms
   l = cr_iscat(iatom,1)        ! Current atom type
   cr_ndiffer(l) = cr_ndiffer(l) + 1
enddo
!
! Step 1:  Determine Eigenvalues and Eigenvectors 
loop_atoms:do iatom=1,cr_ncatoms 
   j = 0
   lsuccess= .false.
   itype = cr_iscat(iatom,1)
   cr_ianis(iatom) = iatom
!  write(*,*)
!  write(*,'(a,i3,2i3, 3f7.3, i4)') ' ATOM ', iatom, itype, cr_iscat(iatom,3), cr_pos(:,iatom), cr_is_sym(iatom)
!  write(*,'(a,        6f7.3    )') ' UANI ', cr_anis(:,itype)
   if(maxval(abs(cr_anis(:,itype)))> 0.0D0) then
      uij(1,1) = cr_anis(1, itype)
      uij(2,2) = cr_anis(2, itype)
      uij(3,3) = cr_anis(3, itype)
      uij(2,3) = cr_anis(4, itype)
      uij(3,2) = cr_anis(4, itype)
      uij(1,3) = cr_anis(5, itype)
      uij(3,1) = cr_anis(5, itype)
      uij(1,2) = cr_anis(6, itype)
      uij(2,1) = cr_anis(6, itype)
      cr_iscat(iatom, 3) = cr_is_sym(iatom)
   else
      ucij      = 0.0D0
      if(cr_dw(itype)>0.0_PREC_DP) then
      ucij(1,1) = cr_dw(itype)/8.0D0/pi**2    ! WORK !?!??!
      ucij(2,2) = cr_dw(itype)/8.0D0/pi**2
      ucij(3,3) = cr_dw(itype)/8.0D0/pi**2
      else
      ucij(1,1) = 1.0/8.0D0/pi**2
      ucij(2,2) = 1.0/8.0D0/pi**2
      ucij(3,3) = 1.0/8.0D0/pi**2
      endif
      uij = matmul((cr_dimat), matmul(ucij, transpose(cr_dimat)))
!      uij       = 0.0_PREC_DP
!     uij (1,1) = cr_dw(itype)/8.0D0/pi**2
!     uij (2,2) = cr_dw(itype)/8.0D0/pi**2
!     uij (3,3) = cr_dw(itype)/8.0D0/pi**2
      cr_iscat(iatom, 3) = 1
!    write(*,*) 
!    write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(1,:), ucij(1,:)
!    write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(2,:), ucij(2,:)
!    write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(3,:), ucij(3,:)
   endif
!
   if(cr_is_sym(iatom)>1) then
      symm_mat = spc_mat(1:3,1:3, cr_is_sym(iatom))
      call matinv(symm_mat, symm_imat)
!write(*,'(a,i3, 3f7.2)') ' Symm Row', 1, symm_mat(1,:)
!write(*,'(a,i3, 3f7.2)') ' Symm Row', 2, symm_mat(2,:)
!write(*,'(a,i3, 3f7.2)') ' Symm Row', 3, symm_mat(3,:)
      uijs = matmul((symm_imat), matmul(uij, (symm_mat)))
      uij = uijs
   endif 
!
!write(*,'(a,i3, 3f10.6)') ' Uij  Row', 1, uij(1,:)
!write(*,'(a,i3, 3f10.6)') ' Uij  Row', 2, uij(2,:)
!write(*,'(a,i3, 3f10.6)') ' Uij  Row', 3, uij(3,:)
   call lookup_anis(cr_nanis, cr_anis_full, cr_prin, uij, j, lsuccess)
!write(*,'(i3, 6f10.6)') j, cr_anis_full(:,j)
   cr_iscat(iatom,3) = j            ! Assign ADP type
   if(lsuccess) then
      cycle loop_atoms
   endif
!
   call calc_prin(cr_iscat(iatom, 3), cr_nanis, uij, cr_dmat, cr_emat, cr_gten, cr_prin)
enddo loop_atoms
!
!cr_nanis = cr_ncatoms
!write(*,*) ' cr_anis_full ', ubound(cr_anis_full)
!write(*,*) ' cr_prin      ', ubound(cr_prin     )
!do j=1, cr_nanis
!write(*,'(i3, 6f10.6)') j, cr_anis_full(:,j)
!enddo
!write(*,*)
!do j=1, cr_nanis
!write(*,'(a,3f9.5,a,f9.5,a,f9.5)') 'E(1,0,0) => ', cr_prin(1,1:3,j),  ' VAL ', cr_prin(1,4,j)
!write(*,'(a,3f9.5,a,f9.5,a,f9.5)') 'E(0,1,0) => ', cr_prin(2,1:3,j),  ' VAL ', cr_prin(2,4,j)
!write(*,'(a,3f9.5,a,f9.5,a,f9.5)') 'E(0,0,1) => ', cr_prin(3,1:3,j),  ' VAL ', cr_prin(3,4,j)
!write(*,*)
!enddo
!
!
end subroutine prep_anis
!
!*******************************************************************************
!
subroutine lookup_anis_6(cr_nanis, cr_anis_full,  cr_prin, uij_6, ientry, lsuccess)
!-
!  Check if Uij are different from any cr_anis_full entry, if so set new entry
!  Version for Uij in SHELX style
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
integer                                          , intent(out)   :: ientry
logical                                          , intent(out)   :: lsuccess
!
real(kind=PREC_DP), dimension(3,3) :: uij
uij(1,1) = uij_6(1)
uij(2,2) = uij_6(2)
uij(3,3) = uij_6(3)
uij(3,2) = uij_6(4)
uij(2,3) = uij_6(4)
uij(3,1) = uij_6(5)
uij(1,3) = uij_6(5)
uij(2,1) = uij_6(6)
uij(1,2) = uij_6(6)
call lookup_anis_3x3(cr_nanis, cr_anis_full,  cr_prin, uij, ientry, lsuccess)
!
end subroutine lookup_anis_6
!
!*******************************************************************************
!
subroutine lookup_anis_3x3(cr_nanis, cr_anis_full,  cr_prin, uij, ientry, lsuccess)
!-
!  Check if Uij are different from any cr_anis_full entry, if so set new entry
!  Version for Uij in as 3x3 matrix
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
real(kind=PREC_DP), dimension(3,3)               , intent(in)    :: uij
integer                                          , intent(out)   :: ientry
logical                                          , intent(out)   :: lsuccess
!
real(kind=PREC_DP), parameter :: TOL = 1.0D-6
!
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
!WORK         cr_iscat(iatom,3) = j            ! Assign ADP type
!WORK         cycle loop_atoms
   endif
enddo loop_nanis
!
if(ientry==-1) then   ! new entry
   cr_nanis = cr_nanis + 1
   ientry   = cr_nanis
   if(cr_nanis>=ubound(cr_anis_full,2)) then
      call alloc_arr(cr_anis_full, 1,6, 1, cr_nanis, all_status, 0.0D0      )
      call alloc_arr(cr_prin, 1, 3,1,4, 1, cr_nanis, all_status, 0.0D0      )
   endif
   cr_anis_full(1,ientry) = uij(1,1)
   cr_anis_full(2,ientry) = uij(2,2)
   cr_anis_full(3,ientry) = uij(3,3)
   cr_anis_full(4,ientry) = uij(3,2)
   cr_anis_full(5,ientry) = uij(3,1)
   cr_anis_full(6,ientry) = uij(2,1)
endif
!
end subroutine lookup_anis_3x3
!
!*******************************************************************************
!
subroutine calc_prin_6(ientry, nanis, uij_6, dmat, emat, gten, prin)
!-
! Calculate principal axis for Uij
!-
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
real(kind=PREC_DP), dimension(6)      , intent(in)  :: uij_6
real(kind=PREC_DP), dimension(3,3)    , intent(in)  :: dmat
real(kind=PREC_DP), dimension(3,3)    , intent(in)  :: emat
real(kind=PREC_DP), dimension(3,3)    , intent(in)  :: gten
real(kind=PREC_DP), dimension(3,4,nanis), intent(out) :: prin
!
real(kind=PREC_DP), dimension(3,3) :: uij
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
call calc_prin_3x3(ientry, nanis, uij, dmat, emat, gten, prin)
!
end subroutine calc_prin_6
!
!*******************************************************************************
!
subroutine calc_prin_3x3(ientry, nanis, uij, dmat, emat, gten, prin)
!-
! Calculate principal axis for Uij
!-
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
real(kind=PREC_DP), dimension(3,3)    , intent(in)  :: uij
real(kind=PREC_DP), dimension(3,3)    , intent(in)  :: dmat
real(kind=PREC_DP), dimension(3,3)    , intent(in)  :: emat
real(kind=PREC_DP), dimension(3,3)    , intent(in)  :: gten
real(kind=PREC_DP), dimension(3,4,nanis), intent(out) :: prin
!
integer :: neigen
real(kind=PREC_DP), dimension(3)   :: eigen_val   ! U_ij in cartesian coordinates
real(kind=PREC_DP), dimension(3,3) :: eigen_car   ! Eigenvectors in cartesian coordinates
real(kind=PREC_DP), DIMENSION(3  ) ::      vector
real(kind=PREC_DP), DIMENSION(3  ) ::      u, v, w
real(kind=PREC_DP), dimension(3,3) :: ucij
real(kind=PREC_DP), DIMENSION(3,3)          ::      imat     = &
   reshape((/1.0D0,0.0D0,0.0D0, 0.0D0,1.0D0,0.0D0, 0.0D0,0.0D0,1.0D0/),SHAPE(imat))
!
ucij = matmul((dmat), matmul(uij, transpose(dmat)))
!     write(*,*) 
!     write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(1,:), ucij(1,:)
!     write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(2,:), ucij(2,:)
!     write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(3,:), ucij(3,:)
!    ucij = matmul((dimat), matmul(uij, transpose(dimat)))
!    write(*,*) 
!    write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(1,:), ucij(1,:)
!    write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(2,:), ucij(2,:)
!    write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(3,:), ucij(3,:)
call eigen_value(ucij, eigen_val, eigen_car, imat, neigen)
prin  (1,1:3, ientry) = eigen_car(:,1)
prin  (2,1:3, ientry) = eigen_car(:,2)
prin  (3,1:3, ientry) = eigen_car(:,3)
prin  (:,4  , ientry) = eigen_val
!write(*,*) ' GOT EIGENVALUES ', ier_num, ier_typ
!write(*,*) ' GOT EIGENVALUES ', eigen_val, neigen
!write(*,*) ' GOT EIGENVECTOR ', eigen_car(:,1)
!write(*,*) ' GOT EIGENVECTOR ', eigen_car(:,2)
!write(*,*) ' GOT EIGENVECTOR ', eigen_car(:,3)
!debug
!     call test_eigen_value(ucij, eigen_val, eigen_car, imat )
!vector = eigen_car(:,1)
!u = matmul((emat), vector)
!vector = eigen_car(:,2)
!v = matmul((emat), vector)
!vector = eigen_car(:,3)
!w = matmul((emat), vector)
!write(*,'(a,3f9.5,a,f9.5,a,f9.5)') 'E(1,0,0) => ', prin(1,1:3,ientry), ' |u|', lib_blen(gten, u), ' VAL ', prin(1,4,ientry)
!write(*,'(a,3f9.5,a,f9.5,a,f9.5)') 'E(0,1,0) => ', prin(2,1:3,ientry), ' |v|', lib_blen(gten, v), ' VAL ', prin(2,4,ientry)
!write(*,'(a,3f9.5,a,f9.5,a,f9.5)') 'E(0,0,1) => ', prin(3,1:3,ientry), ' |w|', lib_blen(gten, w), ' VAL ', prin(3,4,ientry)
!write(*,'(a, f9.5)') ' Ang u,v    ', lib_bang(gten, u,v)
!write(*,'(a, f9.5)') ' Ang u,w    ', lib_bang(gten, u,w)
!write(*,'(a, f9.5)') ' Ang v,w    ', lib_bang(gten, v,w)
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
end module prep_anis_mod
