module prep_anis_mod
!-
! Prepare the anisotropic ADPs
!   Determine eigenvectors and Eigenvalues from U_ij
!   Apply symmetry to Eigenvectors for all sites
!+
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
real(kind=PREC_DP), dimension(3,3) :: uij   ! U_ij in crystal dimensions
real(kind=PREC_DP), dimension(3,3) :: uijs  ! U_ij in crystal dimensions after symmetry
real(kind=PREC_DP), dimension(3,3) :: ucij   ! U_ij in cartesian coordinates
real(kind=PREC_DP), dimension(3)   :: eigen_val   ! U_ij in cartesian coordinates
real(kind=PREC_DP), dimension(3,3) :: eigen_car   ! Eigenvectors in cartesian coordinates
real(kind=PREC_DP), dimension(3,3) :: eigen_cry   ! Eigenvectors in crystal space 
real(kind=PREC_DP), dimension(3,3) :: symm_mat    ! Symmetry matrix
real(kind=PREC_DP), dimension(3,3) :: symm_imat    ! Inverse Symmetry matrix
real(kind=PREC_DP), DIMENSION(3  ) ::      vector
real(kind=PREC_DP), DIMENSION(3  ) ::      u, v, w
real(kind=PREC_DP), DIMENSION(3,3)          ::      imat     = &
   reshape((/1.0D0,0.0D0,0.0D0, 0.0D0,1.0D0,0.0D0, 0.0D0,0.0D0,1.0D0/),SHAPE(imat))
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
cr_ndiffer = 0
cr_nanis = 0
if(maxval(abs(cr_anis(:,:)))<TOL) return
do iatom=1,cr_ncatoms
   l = cr_iscat(iatom)        ! Current atom type
   cr_ndiffer(l) = cr_ndiffer(l) + 1
enddo
!
! Step 1:  Determine Eigenvalues and Eigenvectors 
do iatom=1,cr_ncatoms 
   itype = cr_iscat(iatom)
   cr_ianis(iatom) = iatom
!  write(*,*)
!  write(*,'(a,i3,2i3, 3f7.3, i4)') ' ATOM ', iatom, itype, cr_ianis(iatom), cr_pos(:,iatom), cr_is_sym(iatom)
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
   else
      ucij      = 0.0D0
      ucij(1,1) = cr_dw(itype)/8.0D0/pi**2
      ucij(2,2) = cr_dw(itype)/8.0D0/pi**2
      ucij(3,3) = cr_dw(itype)/8.0D0/pi**2
      uij = matmul((cr_dimat), matmul(ucij, transpose(cr_dimat)))
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
   cr_anis_full(1,iatom) = uij(1,1)
   cr_anis_full(2,iatom) = uij(2,2)
   cr_anis_full(3,iatom) = uij(3,3)
   cr_anis_full(4,iatom) = uij(3,2)
   cr_anis_full(5,iatom) = uij(3,1)
   cr_anis_full(6,iatom) = uij(2,1)
   ucij = 0.0D0
!     do j=1, 3
!       do l=1, 3
!          do m=1, 3
!             do n=1, 3
!                ucij(j,l) = ucij(j,l) + cr_dmat(j,m)*cr_dmat(l,n)*uij(m,n)
!             enddo
!          enddo
!       enddo
!    enddo
!    write(*,'(a,i4)') 'Atom type ', itype
!    write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(1,:), ucij(1,:)
!    write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(2,:), ucij(2,:)
!    write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(3,:), ucij(3,:)
!
     ucij = matmul((cr_dmat), matmul(uij, transpose(cr_dmat)))
!     write(*,*) 
!     write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(1,:), ucij(1,:)
!     write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(2,:), ucij(2,:)
!     write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(3,:), ucij(3,:)
!    ucij = matmul((cr_dimat), matmul(uij, transpose(cr_dimat)))
!    write(*,*) 
!    write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(1,:), ucij(1,:)
!    write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(2,:), ucij(2,:)
!    write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(3,:), ucij(3,:)
     call eigen_value(ucij, eigen_val, eigen_car, imat, neigen)
!write(*,*) ' GOT EIGENVALUES ', ier_num, ier_typ
!write(*,*) ' GOT EIGENVALUES ', eigen_val, neigen
!write(*,*) ' GOT EIGENVECTOR ', eigen_car(:,1)
!write(*,*) ' GOT EIGENVECTOR ', eigen_car(:,2)
!write(*,*) ' GOT EIGENVECTOR ', eigen_car(:,3)
!     call test_eigen_value(ucij, eigen_val, eigen_car, imat )
vector = eigen_car(:,1)
u = matmul((cr_emat), vector)
vector = eigen_car(:,2)
v = matmul((cr_emat), vector)
vector = eigen_car(:,3)
w = matmul((cr_emat), vector)
     cr_prin  (1,1:3, iatom) = eigen_car(:,1)
     cr_prin  (2,1:3, iatom) = eigen_car(:,2)
     cr_prin  (3,1:3, iatom) = eigen_car(:,3)
     cr_prin  (:,4  , iatom) = eigen_val
!write(*,'(a,3f9.5,a,f9.5,a,f9.5)') 'E(1,0,0) => ', cr_prin(1,1:3,iatom), ' |u|', lib_blen(cr_gten, u), ' VAL ', cr_prin(1,4,iatom)
!write(*,'(a,3f9.5,a,f9.5,a,f9.5)') 'E(0,1,0) => ', cr_prin(2,1:3,iatom), ' |v|', lib_blen(cr_gten, v), ' VAL ', cr_prin(2,4,iatom)
!write(*,'(a,3f9.5,a,f9.5,a,f9.5)') 'E(0,0,1) => ', cr_prin(3,1:3,iatom), ' |w|', lib_blen(cr_gten, w), ' VAL ', cr_prin(3,4,iatom)
!write(*,'(a, f9.5)') ' Ang u,v    ', lib_bang(cr_gten, u,v)
!write(*,'(a, f9.5)') ' Ang u,w    ', lib_bang(cr_gten, u,w)
!write(*,'(a, f9.5)') ' Ang v,w    ', lib_bang(cr_gten, v,w)
!     call test_eigen_value(uij, eigen_val, eigen_cry, cr_gten )
enddo
!
cr_nanis = cr_ncatoms
!
write(*,*) ' CR_NANIS ', cr_nanis
!
end subroutine prep_anis
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
