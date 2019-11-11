MODULE blanks_mod
!
CONTAINS
!
!*****7**************************************************************** 
!
SUBROUTINE rem_bl (line, ll) 
!                                                                       
!     Removes all blanks from a string                                  
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: line 
INTEGER          , INTENT(INOUT) :: ll
!
CHARACTER(LEN=1024) :: zeile 
INTEGER :: ibl 
!                                                                       
ibl = INDEX (line (1:ll) , ' ') 
DO while (ibl.gt.0) 
   zeile = ' ' 
   IF (ibl.gt.1) zeile (1:ibl - 1) = line (1:ibl - 1) 
   IF (ibl.lt.ll) then 
      zeile (ibl:ll - 1) = line (ibl + 1:ll) 
   ENDIF 
   ll = ll - 1 
   line = zeile 
   ibl = INDEX (line (1:ll) , ' ') 
ENDDO 
!
END SUBROUTINE rem_bl                         
!
!*****7**************************************************************** 
!
SUBROUTINE rem_dbl_bl(line,ll)
!
IMPLICIT NONE
!
CHARACTER (LEN=*), INTENT(INOUT) :: line
INTEGER          , INTENT(INOUT) :: ll
!
CHARACTER(LEN=LEN(line)) :: string
INTEGER                  :: i,j
LOGICAL                  :: l_bl
!
string=' '
l_bl = .FALSE.
j    = 0
!
DO i=1,ll
  IF(l_bl) THEN   ! CURRENTLY IN MULTIPLE BLANK SECTION
     IF(line(i:i)/=' ') THEN    ! Back to copy
        j = j + 1
        string(j:j) = line(i:i)
        l_bl = .FALSE.
     ENDIF
  ELSE
     j = j + 1
     string(j:j) = line(i:i)
     l_bl = line(i:i)==' '
  ENDIF
ENDDO
line = ' '
line = string
ll   = LEN_TRIM(line)
!
END SUBROUTINE rem_dbl_bl
!
!*****7**************************************************************** 
!
SUBROUTINE rem_leading_bl (line, ll) 
!-
!     Removes all leading blanks from a string                                  
!+
USE charact_mod
!
IMPLICIT none 
!                                                                       
CHARACTER ( LEN=* ), INTENT(INOUT) :: line 
INTEGER            , INTENT(INOUT) :: ll
!
CHARACTER(LEN=1024)  :: zeile 
INTEGER              :: i,j
!                                                                       
zeile = ' '
j     = 1
main: DO i = 1, ll
   j = i
   IF( .NOT. (line(i:i)==' ' .or. line(i:i)==tab)) THEN
      EXIT main
   ENDIF
ENDDO main
zeile = line(j:ll)
ll    = ll - j + 1
line  = zeile
!
END SUBROUTINE rem_leading_bl                         
!
!*****7**************************************************************** 
!
SUBROUTINE rem_insig_bl (line, ll) 
!+
!     Removes all insignificant blanks from a string                                  
!-
USE charact_mod
!
IMPLICIT none 
!                                                                       
CHARACTER ( LEN=* ), INTENT(INOUT) :: line 
INTEGER            , INTENT(INOUT) :: ll
CHARACTER(LEN=1024)  :: zeile 
INTEGER              :: i
INTEGER :: skip
LOGICAL :: l_hyp
LOGICAL :: l_quo
!                                                                       
zeile = ' '
l_hyp = .false.
l_quo = .false.
skip  = 0
main: DO i = 1, ll
   IF(l_hyp) THEN                     ! Within a '' pair keep blanks
      IF(line(i:i)=='''') THEN
         l_hyp = .false.
         zeile(i-skip:i-skip) = line(i:i)
      ELSE
         zeile(i-skip:i-skip) = line(i:i)
      ENDIF
   ELSEIF(l_quo) THEN                 ! Within a "" pair keep blanks
      IF(line(i:i)=='"') THEN
         l_quo = .false.
         zeile(i-skip:i-skip) = line(i:i)
      ELSE
         zeile(i-skip:i-skip) = line(i:i)
      ENDIF
   ELSE                               ! Regular mode delete blanks
      IF(line(i:i)==' '.or. line(i:i)==tab) THEN
         skip = skip + 1
      ELSEIF(line(i:i)=='''') THEN
         l_hyp = .true.
         zeile(i-skip:i-skip) = line(i:i)
      ELSEIF(line(i:i)=='"') THEN
         l_quo = .true.
         zeile(i-skip:i-skip) = line(i:i)
      ELSE
         zeile(i-skip:i-skip) = line(i:i)
      ENDIF
   ENDIF
ENDDO main
ll    = ll - skip
line  = zeile
!
END SUBROUTINE rem_insig_bl                         
!
!*****7**************************************************************** 
!
SUBROUTINE tab2blank(line, ll)
!+
!     Replaces all tabs by a group of 1 blank
!-
USE charact_mod
IMPLICIT none
!                                                                       
CHARACTER ( LEN=* ), INTENT(INOUT) :: line
INTEGER            , INTENT(INOUT) :: ll
!
INTEGER              :: i
!
DO i=1, ll
   IF(line(i:i)== tab) line(i:i) = ' '
ENDDO
!
END SUBROUTINE tab2blank
!
!*****7**************************************************************** 
!
END MODULE blanks_mod
