module wyckoff_mod
!
!
!     This file contains definitions for the symmetry operations
!     and local site symmetry
!
!*****7****************************************************************
!
use precision_mod
!
save
!
integer, PARAMETER  ::   SPC_MAX  =  192
integer, PARAMETER  ::   WYC_MAX  =   48
!
character(LEN=65)   ::   spc_char(1:SPC_MAX)      ! Description of symmetry elements
character(LEN=87)   ::   spc_xyz (1:SPC_MAX)      ! Points on the symmetry elements
integer             ::   spc_n                    ! Number in group
real(kind=PREC_DP)  ::   spc_mat(4,4,1:SPC_MAX)   ! The matrices
real(kind=PREC_DP)  ::   spc_det (1:SPC_MAX)      ! Determinants of symmetry elements
real(kind=PREC_DP)  ::   spc_spur(1:SPC_MAX)      ! Trace of symmetry elements
real(kind=PREC_DP)  ::   spc_axis(3,1:SPC_MAX)    ! Rotation axis or normal to mirror plane
logical, dimension(1:SPC_MAX) :: spc_point        ! Rotation matrix is part of point grup symmetry
logical, dimension(1:SPC_MAX) :: spc_center       ! Rotation matrix is Centering vector
integer, dimension(1:SPC_MAX) :: spc_gen          ! Rotation matrix is  generator matrix
integer, dimension(1:SPC_MAX, 1:SPC_MAX) :: spc_table  ! Multiplication table
!
integer             ::   wyc_n
integer             ::   wyc_list(48)
character(LEN=65)   ::   wyc_char(1:WYC_MAX)      ! Description of symmetry elements
character(LEN=87)   ::   wyc_xyz (1:WYC_MAX)      ! Points on the symmetry elements
real(kind=PREC_DP)  ::   wyc_mat(4,4,1:WYC_MAX)   ! The matrices
real(kind=PREC_DP)  ::   wyc_det (1:WYC_MAX)      ! Determinants of symmetry elements
real(kind=PREC_DP)  ::   wyc_spur(1:WYC_MAX)      ! Trace of symmetry elements
real(kind=PREC_DP)  ::   wyc_axis(3,1:WYC_MAX)    ! Rotation axis or normal to mirror plane
!
end module wyckoff_mod
