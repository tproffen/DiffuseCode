!+
!     This file contains COMMON block for time related variables
!-
      INTEGER            int_time(3)
      INTEGER            int_date(3)
      INTEGER            millisec
      INTEGER            midnight
!
      CHARACTER*24      f_modt
      CHARACTER*24      f_date
!
!
      COMMON /times/ int_time,int_date,f_date,f_modt,                   &
     &               millisec,midnight
!
