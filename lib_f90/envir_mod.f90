MODULE envir_mod
!+
!     Helpfile location etc.
!-
   IMPLICIT NONE
   PUBLIC
   SAVE
!
   CHARACTER(LEN= 256) :: home_dir
   CHARACTER(LEN= 256) :: appl_dir
   CHARACTER(LEN= 256) :: mac_dir
   CHARACTER(LEN= 256) :: umac_dir
   CHARACTER(LEN= 256) :: start_dir
   CHARACTER(LEN= 80 ) :: deffile
   CHARACTER(LEN= 80 ) :: hlpfile
   CHARACTER(LEN= 80 ) :: colorfile
   CHARACTER(LEN= 80 ) :: nullfile
   INTEGER             :: home_dir_l
   INTEGER             :: appl_dir_l
   INTEGER             :: deffile_l
   INTEGER             :: hlpfile_l
   INTEGER             :: colorfile_l
   INTEGER             :: mac_dir_l
   INTEGER             :: umac_dir_l
   INTEGER             :: start_dir_l
   INTEGER             :: lines
!
END MODULE envir_mod
