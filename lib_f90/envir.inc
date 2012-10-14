!+
!     Helpfile location etc.
!-
      CHARACTER*256     home_dir,appl_dir,mac_dir,umac_dir
      CHARACTER*256     start_dir
      CHARACTER*80      deffile,hlpfile,colorfile,nullfile
      INTEGER           home_dir_l,appl_dir_l,deffile_l,hlpfile_l
      INTEGER           colorfile_l,mac_dir_l,umac_dir_l
      INTEGER           start_dir_l
      INTEGER           lines
!
      COMMON /d_envir/  deffile,home_dir,appl_dir,mac_dir,umac_dir,     &
     &                  start_dir,                                      &
     &                  hlpfile,colorfile,nullfile,deffile_l,           &
     &                  home_dir_l,appl_dir_l,hlpfile_l,colorfile_l,    &
     &                  mac_dir_l,umac_dir_l,lines,start_dir_l
