MODULE sys_compiler
!
CONTAINS
!                                                                       
!     This file contains subroutines for:                               
!     Compiler specific routines GNU gfortran version                   
!
!*****7**************************************************************** 
!
SUBROUTINE sys_holecwd (cwd, dummy) 
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(LEN=*), INTENT(OUT) :: cwd 
INTEGER         , INTENT(OUT) :: dummy 
INTEGER  :: getcwd
!                                                                       
dummy = getcwd (cwd )
!                                                                       
END SUBROUTINE sys_holecwd                        
!
!*****7*****************************************************************
!
SUBROUTINE sys_chdir(cwd, dummy)
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN)  :: cwd
INTEGER         , INTENT(OUT) :: dummy
!
dummy = chdir(cwd)
!
END SUBROUTINE sys_chdir
!
!*****7*****************************************************************
!
SUBROUTINE sys_file_info (ifile) 
!+                                                                      
!     Gets and stores file modification time                            
!-                                                                      
USE errlist_mod 
USE times_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: ifile 
!
INTEGER, DIMENSION(13) :: buff
CHARACTER :: ctime
!
CALL fstat (ifile, buff, ier_num) 
IF (ier_num.ne.0) return 
f_modt = ctime(buff(10)) 
!                                                                       
END SUBROUTINE sys_file_info                      
!
!*****7*****************************************************************
!
SUBROUTINE sys_inquire_directory(MAXW, cpara, lpara,lexist)
!+
!   Compiler specific inquiry for directories
!   GNU there is no difference in INQUIRE for files and directories
!-
USE precision_mod
!
INTEGER                                    , INTENT(IN)  :: MAXW
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW), INTENT(IN)  :: cpara
INTEGER                   , DIMENSION(MAXW), INTENT(IN)  :: lpara
LOGICAL                                    , INTENT(OUT) :: lexist
!
INQUIRE (FILE=cpara(1)(1:lpara(1)), EXIST=lexist)
!
END SUBROUTINE sys_inquire_directory
!
!*****7*****************************************************************
!
END MODULE sys_compiler
