!*****7***************************************************************  
!
INTEGER FUNCTION len_str ( string )
!
! Returns the length of the string without trailing blanks, TABS, CR or LF
!
IMPLICIT NONE
!
CHARACTER (LEN=*), INTENT(IN) :: string
CHARACTER (LEN=3), PARAMETER  :: WHITE = achar(9)//achar(10)//achar(13)
INTEGER                       :: itab
INTEGER                       :: laenge
!
laenge = LEN_TRIM(string)               ! intrinsic LEN_TRIM, length without trailing blanks
itab   = SCAN(string(1:laenge),WHITE,.true.) ! to LEN_TRIM a TAB,CR, LF is non-white !@!%!?!
!
DO
  IF ( itab<laenge ) EXIT                ! TAB,CR, LF is not last character thats it
  IF ( laenge==0   ) EXIT                ! empty string
!
!  IF ( laenge==1 ) THEN                 ! length one AND tab is last character
!    laenge = 0                          ! No need to test this, as Fortran 2003 yields
!    EXIT                                ! a string of zero length for 1:laenge-1
!  ENDIF                                 ! with laenge == 1

  laenge = LEN_TRIM(string(1:laenge-1))
  itab   = SCAN(string(1:laenge),WHITE,.true.) ! .true. means scan backwards
ENDDO
!
len_str = laenge
!
END FUNCTION len_str
!
!*****7***************************************************************  
!
INTEGER FUNCTION length_com (string, ikl) 
!-                                                                      
!     Determines the length of the variable or intrinsic function       
!     by searching backwards from the bracket to the first non          
!     character and non '_' character.                                  
!+                                                                      
USE charact_mod
IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER (LEN=*), INTENT(IN) :: string 
INTEGER          , INTENT(IN) :: ikl
!                                                                       
INTEGER :: i, c 
LOGICAL :: lchar 
!                                                                       
i = ikl 
lchar = .true. 
!
DO while (lchar.and.i.gt.1) 
   i = i - 1 
   c = iachar (string (i:i) ) 
   lchar = a.le.c.and.c.le.z.or.c.eq.u .OR. aa.le.c.and.c.le.zz
ENDDO 
!                                                                       
IF (i.eq.1.and.lchar) then 
   length_com = ikl - 1 
ELSE 
   length_com = ikl - i - 1 
ENDIF 
!                                                                       
END FUNCTION length_com                       
!
!*****7***************************************************************  
!
