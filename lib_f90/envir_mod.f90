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
   CHARACTER(LEN= 256) :: current_dir
   CHARACTER(LEN= 256) :: operating
   CHARACTER(LEN= 256) :: user_profile
   CHARACTER(LEN= 256) :: term_scheme_file
   CHARACTER(LEN= 80 ) :: deffile
   CHARACTER(LEN= 80 ) :: hlpfile
   CHARACTER(LEN= 80 ) :: hlpdir
   CHARACTER(LEN= 80 ) :: colorfile
   CHARACTER(LEN= 80 ) :: nullfile
   INTEGER             :: home_dir_l
   INTEGER             :: appl_dir_l
   INTEGER             :: term_scheme_file_l
   INTEGER             :: deffile_l
   INTEGER             :: hlpfile_l
   INTEGER             :: hlp_dir_l
   INTEGER             :: colorfile_l
   INTEGER             :: mac_dir_l
   INTEGER             :: umac_dir_l
   INTEGER             :: start_dir_l
   INTEGER             :: current_dir_l
   INTEGER             :: lines
   INTEGER             :: PID
   LOGICAL             :: term_scheme_exists = .FALSE.
!
END MODULE envir_mod
