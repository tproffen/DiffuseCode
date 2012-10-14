!
      INTEGER            VAR_MAX
      PARAMETER      (VAR_MAX=150)
!
      INTEGER            VAR_TYPE_INTE
      PARAMETER      (VAR_TYPE_INTE=0)
      INTEGER            VAR_TYPE_REAL
      PARAMETER      (VAR_TYPE_REAL=1)
      INTEGER            VAR_TYPE_CHAR
      PARAMETER      (VAR_TYPE_CHAR=2)
!
      INTEGER            var_num
      CHARACTER*200      var_name(VAR_MAX)
      CHARACTER*200      var_char(VAR_MAX)
      INTEGER            var_l   (VAR_MAX)
      INTEGER            var_type(VAR_MAX)
      REAL            var_val (VAR_MAX)
!
      COMMON /variables/ var_num,var_name,var_char,var_l,var_type,      &
     &                              var_val
