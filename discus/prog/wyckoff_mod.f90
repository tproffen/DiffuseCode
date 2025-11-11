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
character(len=8), dimension(0:32), parameter :: wyc_point_names =  (/ &
'UNKNOWN ', &
'1       ', &    !  1
'-1      ', &    !  2
'2       ', &    !  3
'm       ', &    !  4
'2/m     ', &    !  5
'2 2 2   ', &    !  6
'm m 2   ', &    !  7
'm m m   ', &    !  8
'4       ', &    !  9
'-4      ', &    ! 10
'4/m     ', &    ! 11
'4 2 2   ', &    ! 12
'4 m m   ', &    ! 13
'-4 2 m  ', &    ! 14
'4/m m m ', &    ! 15
'3       ', &    ! 16
'-3      ', &    ! 17
'3 2 1   ', &    ! 18
'3 m 1   ', &    ! 19
'-3 m 1  ', &    ! 20
'6       ', &    ! 21
'-6      ', &    ! 22
'6/m     ', &    ! 23
'6 2 2   ', &    ! 24
'6 m m   ', &    ! 25
'-6 m 2  ', &    ! 26
'6/m m m ', &    ! 27
'2 3     ', &    ! 28
'm -3    ', &    ! 29
'4 3 2   ', &    ! 30
'-4 3 m  ', &    ! 31
'm -3 m  '  &    ! 32
/)
!
character(LEN=65)   ::   spc_char(1:SPC_MAX)      ! Description of symmetry elements
character(LEN=87)   ::   spc_xyz (1:SPC_MAX)      ! Points on the symmetry elements
integer             ::   spc_n                    ! Number in group
real(kind=PREC_DP)  ::   spc_mat(4,4,1:SPC_MAX)   ! The matrices
real(kind=PREC_DP)  ::   spc_det (1:SPC_MAX)      ! Determinants of symmetry elements
real(kind=PREC_DP)  ::   spc_spur(1:SPC_MAX)      ! Trace of symmetry elements
real(kind=PREC_DP)  ::   spc_axis(3,1:SPC_MAX)    ! Rotation axis or normal to mirror plane
real(kind=PREC_DP)  ::   spc_1bar(3)              ! Position of -1 operation if centrosymmetric
logical, dimension(1:SPC_MAX) :: spc_point        ! Rotation matrix is part of point grup symmetry
logical, dimension(1:SPC_MAX) :: spc_center       ! Rotation matrix is Centering vector
integer, dimension(1:SPC_MAX) :: spc_gen          ! Rotation matrix is  generator matrix
integer, dimension(1:SPC_MAX, 1:SPC_MAX) :: spc_table  ! Multiplication table
!
integer             ::   wyc_n
integer             ::   wyc_list(48)
character(LEN=65)   ::   wyc_char(1:WYC_MAX)      ! Description of symmetry elements
character(LEN=87)   ::   wyc_xyz (1:WYC_MAX)      ! Points on the symmetry elements
integer             ::   wyc_point_no             ! Point group number
real(kind=PREC_DP)  ::   wyc_mat(4,4,1:WYC_MAX)   ! The matrices
real(kind=PREC_DP)  ::   wyc_det (1:WYC_MAX)      ! Determinants of symmetry elements
real(kind=PREC_DP)  ::   wyc_spur(1:WYC_MAX)      ! Trace of symmetry elements
real(kind=PREC_DP)  ::   wyc_axis(3,1:WYC_MAX)    ! Rotation axis or normal to mirror plane
logical           , dimension(3) :: wyc_fix          ! Coordinate is fixed
character(len=16) , dimension(3) :: wyc_fix_pos      ! Coordinate is fixed to these position
integer           , dimension(4,4) :: wyc_fix_mat   ! matrix for:  MAT*(x,y,z)
!
end module wyckoff_mod
