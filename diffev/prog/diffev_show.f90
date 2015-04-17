MODULE diffev_show_mod
!
CONTAINS
!*****7*****************************************************************
!                                                                       
!     This routine shows all settings of DIFFEV parameters. Some        
!     parameter settings can also be displayed by entering the          
!     corresponding command without parameters.                         
!                                                                       
!*****7*****************************************************************
SUBROUTINE diffev_do_show (line, lp) 
!                                                                       
!     Main show menu                                                    
!                                                                       
USE diffev_allocate_appl
!
USE errlist_mod 
USE prompt_mod 
!
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER   :: maxw = 2
!                                                                       
!                                                                       
CHARACTER (LEN= *  ), INTENT(INOUT) :: line 
INTEGER             , INTENT(INOUT) :: lp 
!
CHARACTER (LEN=1024), DIMENSION(maxw)  :: cpara (maxw) 
INTEGER             , DIMENSION(maxw)  :: lpara (maxw)
INTEGER                                :: ianz
LOGICAL                                :: str_comp 
!                                                                       
CALL get_params (line, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
IF (ianz.eq.1) then 
   IF (str_comp (cpara (1) , 'config', 1, lpara (1) , 6) )  THEN
      CALL diffev_show_config
   ELSE
!                                                                       
!     -- try generic show commands                                      
!                                                                       
      CALL do_show_generic (cpara, lpara, maxw) 
   ENDIF
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
END SUBROUTINE diffev_do_show                        
END MODULE diffev_show_mod
