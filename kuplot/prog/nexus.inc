!+
!     Variables for NeXus file handling
!-
      CHARACTER*200   run_title,run_stime,run_etime
      CHARACTER*200   nxs_fname
      INTEGER*4       run_iexp,run_irun
      INTEGER*4       nxs_id(NXHANDLESIZE)
      LOGICAL         nxs_open
!
      COMMON /nxsi/   run_title,run_stime,run_etime,                    &
     &       run_iexp,run_irun
      COMMON /nxs/    nxs_fname,nxs_id,nxs_open

