module exp3_pdf_data
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), dimension(:), allocatable :: exp3_x   ! Coordinates along axis
real(kind=PREC_DP), dimension(:), allocatable :: exp3_y
real(kind=PREC_DP), dimension(:), allocatable :: exp3_z
real(kind=PREC_DP), dimension(:,:,:), allocatable :: exp3_data
real(kind=PREC_DP), dimension(:,:,:), allocatable :: exp3_sigma
logical           , dimension(:,:,:), allocatable :: exp3_mask   ! is false if no data at this point
!
end module exp3_pdf_data
