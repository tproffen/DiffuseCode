!+
!     Common block for variablves used during do and if execution
!-
      CHARACTER*1024 do_comm(0:MAXCOM,0:MAXLEV)

      INTEGER       do_leng(0:MAXCOM,0:MAXLEV)
      INTEGER       level,nlevel(0:MAXLEV),ilevel(0:MAXLEV)
      INTEGER       jump(0:MAXLEV)
      LOGICAL       ltest(0:MAXLEV)
!
      COMMON /doexec/ do_comm,do_leng,ilevel,jump,level,ltest,nlevel
!
