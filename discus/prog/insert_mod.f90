MODULE insert_mod
!-
!     Variables needed to insert objects
!+
!
SAVE
!
INTEGER, PARAMETER  ::  INS_NEWTYPE  =  -1
!
CHARACTER(LEN=4  )    ::  ins_obj_atom
CHARACTER(LEN=200)    ::  ins_file
INTEGER               ::  ins_CHARACTER
INTEGER               ::  ins_type
REAL   , DIMENSION(3) ::  ins_origin  ! (3)
REAL   , DIMENSION(3) ::  ins_cent    ! (3)
REAL                  ::  ins_density
REAL                  ::  ins_fuzzy
REAL                  ::  ins_adp
REAL   , DIMENSION(3) ::  ins_xaxis   ! (3)
REAL   , DIMENSION(3) ::  ins_yaxis   ! (3)
REAL   , DIMENSION(3) ::  ins_zaxis   ! (3)
REAL   , DIMENSION(3) ::  ins_xdim    ! (3)
REAL   , DIMENSION(3) ::  ins_ydim    ! (3)
REAL   , DIMENSION(3) ::  ins_zdim    ! (3)
!
END MODULE insert_mod
