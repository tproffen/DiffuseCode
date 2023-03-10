MODULE jsu_readline
!-------------------------------------------------------------------------------
   USE ISO_C_BINDING
   IMPLICIT NONE
   PRIVATE
   PUBLIC iso_readline
   public iso_readline_add
!-------------------------------------------------------------------------------
! define the call to the C routine
! extern char     *Freadline(int ilen, char *buf, char prompt[]);
  PUBLIC ::  Freadline
   INTERFACE
      SUBROUTINE Freadline(ilen,buf,prompt) BIND(C,NAME='FCreadline')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(KIND=C_INT),INTENT(IN),VALUE      ::  ilen
         CHARACTER(KIND=C_CHAR),intent(out)  ::  buf(*)
         CHARACTER(KIND=C_CHAR),intent(in)   ::  prompt(*)
      END SUBROUTINE Freadline
   END INTERFACE
! extern char     *Freadline_add(int ilen, char *buf);
  PUBLIC ::  Freadline_add
   INTERFACE
      SUBROUTINE Freadline_add(ilen,buf) BIND(C,NAME='FCreadline_add')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(KIND=C_INT),INTENT(IN),VALUE      ::  ilen
         CHARACTER(KIND=C_CHAR),intent(in)  ::  buf(*)
      END SUBROUTINE Freadline_add
   END INTERFACE
!-------------------------------------------------------------------------------
contains
! the routine that calls the C routine
SUBROUTINE iso_readline(line,prompt)
   USE ISO_C_BINDING
   IMPLICIT NONE
   CHARACTER(KIND=C_CHAR,LEN=*),INTENT(OUT) :: line
   CHARACTER(KIND=C_CHAR,LEN=*),INTENT(IN)  :: prompt

   ! trim to last non-blank character and append null for C
   CALL Freadline(LEN(line),line,prompt(:LEN_TRIM(prompt))//ACHAR(0))

 END SUBROUTINE iso_readline
!-------------------------------------------------------------------------------
subroutine iso_readline_add(line)
!-
! Add line to the readline history
!+
use iso_c_binding
implicit none
character(kind=C_CHAR, len=*), intent(in) :: line
!
call Freadline_add(len(line),line(:LEN_TRIM(line))//achar(0))
!
end subroutine iso_readline_add
!-------------------------------------------------------------------------------
END MODULE jsu_readline
