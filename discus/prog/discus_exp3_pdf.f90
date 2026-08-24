module exp3_pdf_menu
!
!  Convert 3D intensities to a 3D-PDF
!
contains
!
!*******************************************************************************
!
subroutine exp3_pdf
!-
! Main menu
!+
use exp3_pdf_load
!
use precision_mod
!
implicit none
character(len=PREC_STRING) :: infile
!
infile = 'DATA/gly_ala.hkl'
call load_merge_hklf4(infile)
!
end subroutine exp3_pdf
!
!*******************************************************************************
!
end module exp3_pdf_menu
