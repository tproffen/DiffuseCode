MODULE prompt_mod
!+
!     Variables for program prompt
!-
   IMPLICIT NONE
   PUBLIC
   SAVE
!
   INTEGER, PARAMETER :: PROMPT_ON       = 0
   INTEGER, PARAMETER :: PROMPT_OFF      = 1
   INTEGER, PARAMETER :: PROMPT_MACRO    = 2
   INTEGER, PARAMETER :: PROMPT_REDIRECT = 3
!
   INTEGER, PARAMETER :: OUTPUT_SCREEN   = 1
   INTEGER, PARAMETER :: OUTPUT_NONE     = 2
   INTEGER, PARAMETER :: OUTPUT_FILE     = 3
!
   CHARACTER(LEN= 1024) :: blank
   CHARACTER(LEN= 1024) :: input_gui
   CHARACTER(LEN= 80  ) :: s_ipallowed
   CHARACTER(LEN= 40  ) :: prompt
   CHARACTER(LEN= 40  ) :: oprompt
   CHARACTER(LEN= 10  ) :: version
   CHARACTER(LEN= 7   ) :: pname,pname_cap
   INTEGER              :: s_port
   INTEGER              :: s_sock
   INTEGER              :: s_conid
   INTEGER              :: s_remote
   INTEGER              :: prompt_status     = PROMPT_ON
   INTEGER              :: prompt_status_old = PROMPT_ON
   INTEGER              :: output_status     = OUTPUT_SCREEN
   INTEGER              :: output_status_old = OUTPUT_SCREEN
   INTEGER              :: output_io
   INTEGER              :: socket_status     = OUTPUT_SCREEN
   INTEGER              :: socket_status_old = OUTPUT_SCREEN
   LOGICAL              :: first_input
   LOGICAL              :: lsocket           = .false.
   LOGICAL              :: lconn             = .false.
   LOGICAL              :: lremote           = .false.
   LOGICAL              :: lsetup_done       = .false.
   LOGICAL              :: lstandalone       = .true.
   LOGICAL              :: linteractive      = .true.
!
END MODULE prompt_mod
