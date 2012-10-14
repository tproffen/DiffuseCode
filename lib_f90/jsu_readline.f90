!-------------------------------------------------------------------------------
MODULE jsu_readline
   USE ISO_C_BINDING
   IMPLICIT NONE
   PRIVATE
   PUBLIC iso_readline
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
END MODULE jsu_readline
!-------------------------------------------------------------------------------
