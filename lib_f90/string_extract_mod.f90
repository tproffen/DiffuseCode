MODULE string_extract_mod
!
!
IMPLICIT NONE
PRIVATE
PUBLIC string_extract
!
CONTAINS
!
!###############################################################################
!
SUBROUTINE string_extract(string, istart, iend, substring, s1, s2, s3)
!
!  Extract the first part of the string that is not enclosed in
!  single or double quotation marks
!
!  Returns:
!  substring   First part that is not in quotes
!  s1          start point within string
!  s2          end   point within string
!  s3          Next point after paired quotes
!  
IMPLICIT NONE
!
CHARACTER(LEN=*) , INTENT(IN)  :: string
INTEGER          , INTENT(IN)  :: istart
INTEGER          , INTENT(IN)  :: iend
CHARACTER(LEN=*) , INTENT(OUT) :: substring
INTEGER          , INTENT(OUT) :: s1
INTEGER          , INTENT(OUT) :: s2
INTEGER          , INTENT(OUT) :: s3
!
CHARACTER(LEN=1) :: quote
INTEGER :: ising
INTEGER :: idbl
INTEGER :: i
!
ising = INDEX(string(istart:iend),'''')
idbl  = INDEX(string(istart:iend),'"' )
!
IF(ising==0) THEN
   IF(idbl==0) THEN                    ! ising==0  idbl==0
      s1 = istart
      s2 = iend
   ELSE                                ! ising==0  idbl/=0
      s1 = istart
      s2 = idbl-1 + istart - 1
   ENDIF
ELSE
   IF(idbl==0) THEN                    ! ising/=0  idbl==0
      s1 = istart
      s2 = ising-1 + istart - 1
   ELSE                                ! ising/=0  idbl/=0
      s1 = istart
      s2 = MIN(ising, idbl) + 1 + istart - 1
   ENDIF
ENDIF
substring = string(s1:s2)
!
s3 = iend
IF(s2<iend-2) THEN      ! Find end of section in quote
   quote = string(s2+1:s2+1)
   s3    = s2 + INDEX(string(s2+2:iend),quote) + 2
ENDIF
!
END SUBROUTINE string_extract
!
!###############################################################################
!

END MODULE string_extract_mod
