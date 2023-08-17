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
!
use errlist_mod
use math_sup, only:eigen_value, test_eigen_value
use precision_mod
use wink_mod
!
implicit none
!
integer :: itype    ! Loop index atom types
integer :: iatom    ! Loop index atoms
integer :: j, l, m, n  ! Dummy loop indices
integer :: neigen      ! Number distinct eigenvectors
real(kind=PREC_DP), dimension(3,3) :: uij   ! U_ij in crystal dimensions
real(kind=PREC_DP), dimension(3,3) :: ucij   ! U_ij in cartesian coordinates
real(kind=PREC_DP), dimension(3)   :: eigen_val   ! U_ij in cartesian coordinates
real(kind=PREC_DP), dimension(3,3) :: eigen_vec   ! U_ij in cartesian coordinates
real(kind=PREC_DP), DIMENSION(3,3)          ::      imat     = &
   reshape((/1.0D0,0.0D0,0.0D0, 0.0D0,1.0D0,0.0D0, 0.0D0,0.0D0,1.0D0/),SHAPE(imat))
!
write(*,*) ' PREP_ANIS cr_anis  ', ubound(cr_anis)
write(*,*) ' PREP_ANIS cr_anis_f', ubound(cr_anis_full)
write(*,*) ' PREP_ANIS cr_prin  ', ubound(cr_prin,3)
write(*,*)
!
if(allocated(cr_anis_full)) deallocate(cr_anis_full)
if(allocated(cr_prin     )) deallocate(cr_prin     )
allocate(cr_anis_full(6,    cr_ncatoms))
allocate(cr_prin     (3, 3, cr_ncatoms))
!
! Step 1:  Determine Eigenvalues and Eigenvectors 
do iatom=1,cr_ncatoms 
   itype = cr_iscat(iatom)
   write(*,*)
   write(*,'(a,i3, i3, 3f7.3, i4)') ' ATOM ', itype, cr_iscat(iatom), cr_pos(:,iatom), cr_is_sym(iatom)
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
      ucij(1,1) = cr_dw(itype)/8.0D0/pi
      ucij(2,2) = cr_dw(itype)/8.0D0/pi
      ucij(3,3) = cr_dw(itype)/8.0D0/pi
      uij = matmul((cr_dimat), matmul(ucij, transpose(cr_dimat)))
     write(*,*) 
     write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(1,:), ucij(1,:)
     write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(2,:), ucij(2,:)
     write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(3,:), ucij(3,:)
   endif
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
     write(*,*) 
     write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(1,:), ucij(1,:)
     write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(2,:), ucij(2,:)
     write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(3,:), ucij(3,:)
!    ucij = matmul((cr_dimat), matmul(uij, transpose(cr_dimat)))
!    write(*,*) 
!    write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(1,:), ucij(1,:)
!    write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(2,:), ucij(2,:)
!    write(*,'(a,2(2x,3(2x,f9.6)))') 'UC_ij ', uij(3,:), ucij(3,:)
     call eigen_value(ucij, eigen_val, eigen_vec, imat, neigen)
write(*,*) ' GOT EIGENVALUES ', ier_num, ier_typ
write(*,*) ' GOT EIGENVALUES ', eigen_val, neigen
write(*,*) ' GOT EIGENVECTOR ', eigen_vec(:,1)
write(*,*) ' GOT EIGENVECTOR ', eigen_vec(:,2)
write(*,*) ' GOT EIGENVECTOR ', eigen_vec(:,3)
     call test_eigen_value(ucij, eigen_val, eigen_vec, imat )
enddo
!
!
!
end subroutine prep_anis
!
!*******************************************************************************
!
end module prep_anis_mod
