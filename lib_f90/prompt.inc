!+
!     Variables for program prompt
!-
      INTEGER         PROMPT_ON,PROMPT_OFF,PROMPT_MACRO
      INTEGER         PROMPT_REDIRECT
!
      INTEGER            OUTPUT_SCREEN,OUTPUT_NONE,OUTPUT_FILE
!
      PARAMETER      (PROMPT_ON       = 0)
      PARAMETER      (PROMPT_OFF      = 1)
      PARAMETER      (PROMPT_MACRO    = 2)
      PARAMETER      (PROMPT_REDIRECT = 3)
!
      PARAMETER      (OUTPUT_SCREEN   = 1)
      PARAMETER      (OUTPUT_NONE     = 2)
      PARAMETER      (OUTPUT_FILE     = 3)
!
      CHARACTER*1024  blank
      CHARACTER*80      s_ipallowed
      CHARACTER*40      prompt,oprompt
      CHARACTER*10      version
      CHARACTER*7     pname,pname_cap
      INTEGER            s_port,s_sock,s_conid,s_remote
      INTEGER         prompt_status
      INTEGER         prompt_status_old
      INTEGER         output_status
      INTEGER         output_status_old
      INTEGER         output_io
      LOGICAL            first_input
      LOGICAL            lsocket,lconn,lremote
!
      COMMON /anzeig/ blank,prompt,oprompt,pname,version,pname_cap
      COMMON /anzeig_s/ output_status,output_status_old,output_io,      &
     &                  prompt_status,prompt_status_old,first_input
      COMMON /remote/   s_ipallowed,s_port,s_sock,s_conid,              &
     &                  s_remote,lsocket,lconn,lremote
