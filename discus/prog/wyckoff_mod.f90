MODULE wyckoff_mod
!
!
!     This file contains definitions for the symmetry operations
!     and local site symmetry
!
!*****7****************************************************************
!
      SAVE
!
      INTEGER, PARAMETER  ::   SPC_MAX  =  192
!
      CHARACTER(LEN=65)   ::   spc_char(1:SPC_MAX)
      CHARACTER(LEN=87)   ::   spc_xyz (1:SPC_MAX)
      INTEGER             ::   spc_n
      REAL                ::   spc_mat(4,4,1:SPC_MAX)
      REAL                ::   spc_det (1:SPC_MAX)
      REAL                ::   spc_spur(1:SPC_MAX)
!
      INTEGER             ::   wyc_n
      INTEGER             ::   wyc_list(48)
!
END MODULE wyckoff_mod
