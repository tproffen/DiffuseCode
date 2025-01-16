module fourier_form_generic
!
!  Multiply a complex scattering function with:
!     isotropic / anisotropic  Debye -Waller
!     tabulted form factors or DISCAMB form factors
!
!
contains
!
!*******************************************************************************
!
subroutine tcsf_form_generic(lnufft, ldiscamb, lform, is_anis, iscat, isym, ientry, fnum, tcsfp)
!-
!  Call appropriate special function
!+
!
use diffuse_mod
!
use precision_mod
!
implicit none
!
logical              , intent(in) :: lnufft   ! Fourier via NUFFT==T or TURBO==F
logical              , intent(in) :: ldiscamb ! Fourier via DISCAMB form factors==T or analytic form factors==F
logical              , intent(in) :: lform    ! use tabulated form factors
logical              , intent(in) :: is_anis  ! At least one atom hase anisotropic ADPs
integer              , intent(in) :: iscat    ! Use this scattering type
integer              , intent(in) :: isym     ! Use this symmetry type
integer              , intent(in) :: ientry   ! Use this ADP entry type
integer, dimension(3), intent(in) :: fnum     ! Grid sizes in reciprocal space
complex(kind=PREC_DP), dimension(num(1), num(2), num(3)), intent(inout) :: tcsfp   ! Partial structure factor
!
integer, dimension(2) :: ldim_c
integer, dimension(2) :: udim_c
integer, dimension(3) :: ldim_i
integer, dimension(3) :: udim_i
integer               :: nanis
!
ldim_c = lbound(cfact)
udim_c = ubound(cfact)
ldim_i = lbound(istl)
udim_i = ubound(istl)
nanis  = ubound(four_dbw,4)
!
if(ldiscamb) then
!  call tcsf_form_discamb(iscat, isym, ientry, fnum, num, tcsfp)
   call tcsf_form_discamb(iscat, isym, ientry, fnum, num, eck, vi,  &
           nanis, four_dbw, tcsfp)
elseif(lform) then
   if(is_anis) then
!      call tcsf_form_aniso(iscat, ientry, fnum, tcsfp)
      call tcsf_form_aniso(iscat, ientry, fnum, ldim_c, udim_c, cfact, &
           ldim_i, udim_i, istl, nanis, four_dbw, num, tcsfp)
   else
      call tcsf_form_iso(iscat, fnum, ldim_c, udim_c, cfact, ldim_i, udim_i, istl, num, tcsfp)
   endif
endif
!
end subroutine tcsf_form_generic
!
!*******************************************************************************
!
subroutine tcsf_form_discamb(iscat, isym, ientry, fnum, num, eck, vi,  &
           nanis, four_dbw, tcsfp)
!-
!  Multiply tcsf with DISCAMB atomic form factor, anisotropic Uij
!+
!
!use diffuse_mod
use discamb_mod
!
use precision_mod
!
implicit none
!
integer              , intent(in) :: iscat   ! Scattering type
integer              , intent(in) :: isym    ! Symmetry type
integer              , intent(in) :: ientry  ! ADP type
integer, dimension(3), intent(in) :: fnum
integer, dimension(3), intent(in) ::  num
real(kind=PREC_DP), dimension(3,4), intent(in) :: eck
real(kind=PREC_DP), dimension(3,3), intent(in) :: vi 
integer              , intent(in) :: nanis   ! Number of ADPs == Upper boundary dimension4 in four_dbw
real   (kind=PREC_DP), dimension(num(1), num(2), num(3), nanis), intent(in) :: four_dbw ! Debye Waller term
complex(kind=PREC_DP), dimension(num(1), num(2), num(3)), intent(inout) :: tcsfp   ! Partial structure factor
!
integer :: i,j,h
integer, dimension(3) :: ii
!
do h = 1, fnum(3)
   do j = 1, fnum(2)
      do i = 1, fnum(1)
         ii(1) = nint(eck(1,1) + (i-1)*vi(1,1) + (j-1)*vi(1,2) + (h-1)*vi(1,3))
         ii(2) = nint(eck(2,1) + (i-1)*vi(2,1) + (j-1)*vi(2,2) + (h-1)*vi(2,3))
         ii(3) = nint(eck(3,1) + (i-1)*vi(3,1) + (j-1)*vi(3,2) + (h-1)*vi(3,3))
         tcsfp(i,j,h) = (tcsfp(i,j,h) * disc_list(iscat)%four_form_tsc(ii(1),ii(2),ii(3),isym)) &
                                     * four_dbw(i,j,h,ientry)
      enddo
   enddo
enddo   
!
end subroutine tcsf_form_discamb
!
!*******************************************************************************
!
subroutine tcsf_form_aniso(iscat, ientry, fnum, ldim_c, udim_c, cfact, &
           ldim_i, udim_i, istl, nanis, four_dbw, num, tcsfp)
!-
!  Multiply tcsf with atomic form factor, anisotropic Uij
!+
!
use precision_mod
!
implicit none
!
integer              , intent(in) :: iscat   ! Scattering type
integer              , intent(in) :: ientry  ! ADP type
integer, dimension(3), intent(in) :: fnum    ! Range of indices to be calculated
integer, dimension(2), intent(in) :: ldim_c  ! Lower boundary cfact
integer, dimension(2), intent(in) :: udim_c  ! Upper boundary cfact
!
complex(kind=PREC_DP), dimension(ldim_c(1):udim_c(1), ldim_c(2):udim_c(2)), intent(in)    :: cfact
!
integer, dimension(3), intent(in) :: ldim_i  ! Lower boundary istl
integer, dimension(3), intent(in) :: udim_i  ! Upper boundary istl
integer, dimension(ldim_i(1):udim_i(1), ldim_i(2):udim_i(2), ldim_i(3):udim_i(3)), intent(in) :: istl ! sin(theta)/lambda lookup
integer              , intent(in) :: nanis   ! Number of ADPs == Upper boundary dimension4 in four_dbw
integer, dimension(3), intent(in) :: num     ! Full sizes Fourier range
!
real   (kind=PREC_DP), dimension(num(1), num(2), num(3), nanis), intent(in) :: four_dbw ! Debye Waller term
complex(kind=PREC_DP), dimension(num(1), num(2), num(3))       , intent(inout) :: tcsfp ! Partial structure factor
!
integer :: i,j,h
!
do h = 1, fnum(3)
   do j = 1, fnum(2)
      do i = 1, fnum(1)
         tcsfp(i,j,h) = tcsfp(i,j,h) * cfact(istl(i,j,h), iscat) * four_dbw(i,j,h, ientry)  ! Anisotropic ADP
      enddo
   enddo
enddo   
!
end subroutine tcsf_form_aniso
!
!*******************************************************************************
!
subroutine tcsf_form_iso(iscat, fnum, ldim_c, udim_c, cfact, ldim_i, udim_i, istl, num, tcsfp)
!-
!  Multiply tcsf with atomic form factor, isotropic Uij
!+
!
use precision_mod
!
implicit none
!
integer                                                           , intent(in) :: iscat
integer              , dimension(3)                               , intent(in) :: fnum
integer              , dimension(2)                               , intent(in) :: ldim_c
integer              , dimension(2)                               , intent(in) :: udim_c
complex(kind=PREC_DP), dimension(ldim_c(1):udim_c(1), ldim_c(2):udim_c(2)), intent(in)    :: cfact
integer              , dimension(3)                               , intent(in) :: ldim_i
integer              , dimension(3)                               , intent(in) :: udim_i
integer              , dimension(ldim_i(1):udim_i(1), ldim_i(2):udim_i(2), ldim_i(3):udim_i(3)), intent(in)    :: istl  ! sin(theta)/lambda lookup
integer              , dimension(3)                               , intent(in) :: num
complex(kind=PREC_DP), dimension( num(1),  num(2),  num(3))       , intent(inout) :: tcsfp ! Partial structure factor
!
integer :: i,j,h
!
do h = 1, fnum(3)
   do j = 1, fnum(2)
      do i = 1, fnum(1)
         tcsfp(i,j,h) = tcsfp(i,j,h) * cfact(istl(i,j,h), iscat)          ! Isotropic ADP
      enddo
   enddo
enddo   
!
end subroutine tcsf_form_iso
!
!*******************************************************************************
!
end module fourier_form_generic
