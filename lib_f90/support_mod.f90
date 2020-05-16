MODULE support_mod
!
! Contains basic support routines 
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE remove_comment (line, ll) 
!                                                                       
!     removes trailing in line comments                                 
!                                                                       
USE charact_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: line 
INTEGER          , INTENT(INOUT) ::ll 
!                                                                       
INTEGER :: i 
LOGICAL :: quote 
LOGICAL :: search 
!                                                                       
search = .true. 
DO i = 1, ll 
   quote = line (i:i) .eq.'"'.or.line (i:i) .eq.'''' 
   IF (quote) THEN 
      search = .not.search 
   ENDIF 
   IF (search) THEN 
      IF (line (i:i) .eq.'#'.or.line (i:i) .eq.'!') THEN 
         line (i:ll) = ' ' 
         ll = i - 1 
      ENDIF 
   ENDIF 
ENDDO 
ll = LEN_TRIM(line)
IF(ll>0) THEN
   DO WHILE(line(ll:ll)==TAB)      !  Remove trailing TABs
      line(ll:ll) = ' '
      ll = LEN_TRIM(line)
   ENDDO
ENDIF
!                                                                       
END SUBROUTINE remove_comment                 
!
!*******************************************************************************
!
END MODULE support_mod
