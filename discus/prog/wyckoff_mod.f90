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
!
character(LEN=65)   ::   spc_char(1:SPC_MAX)
character(LEN=87)   ::   spc_xyz (1:SPC_MAX)
integer             ::   spc_n
real(kind=PREC_DP)  ::   spc_mat(4,4,1:SPC_MAX)
real(kind=PREC_DP)  ::   spc_det (1:SPC_MAX)
real(kind=PREC_DP)  ::   spc_spur(1:SPC_MAX)
logical, dimension(1:SPC_MAX) :: spc_point        ! Rotation matrix is part of point grup symmetry
!
integer             ::   wyc_n
integer             ::   wyc_list(48)
!
end module wyckoff_mod
