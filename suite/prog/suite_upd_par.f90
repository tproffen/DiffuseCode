SUBROUTINE suite_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz) 
!                                                                       
!-                                                                      
!       replaces a substring in an expression by the value of the       
!       appropriate Parameter.                                          
!     This is needed if the parameter is read.                          
!+                                                                      
USE blanks_mod
USE errlist_mod 
USE param_mod 
!
IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER (LEN= * )  , INTENT(INOUT) :: string 
INTEGER              , INTENT(IN   ) :: ikl
INTEGER              , INTENT(IN   ) :: iklz
INTEGER              , INTENT(INOUT) :: ll
INTEGER              , INTENT(IN   ) :: maxw
INTEGER              , INTENT(IN   ) :: ianz
REAL, DIMENSION(MAXW), INTENT(IN   ) :: ww
!
CHARACTER (LEN=1024)                 :: zeile 
!                                                                       
INTEGER                              :: laenge, ltyp, kpara, kpara2
INTEGER                              :: lcomm 
INTEGER                              :: length_com 
!                                                                       
laenge = ll 
ltyp = 1 
zeile = ' ' 
kpara = nint (ww (1) ) 
kpara2 = 1
IF (maxw.ge.2) then 
   kpara2 = nint (ww (2) ) 
ENDIF 
!                                                                 
lcomm = length_com (string, ikl) 
!                                                                 
IF (lcomm.eq.1) then 
!                                                                 
   IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string (1: ikl - lcomm - 1)                                               
   IF (string (ikl - 1:ikl - 1) .eq.'i') then 
      IF (ianz.eq.1) then 
         IF (0.le.kpara.and.kpara.le.MAXPAR) then 
            WRITE (zeile (ikl - 1:ikl + 13) , '(i15)') inpara ( kpara)                                                
         ELSE 
            ier_num = - 8 
            ier_typ = ER_FORT 
         ENDIF 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 1:ikl - 1) .eq.'r') then 
      IF (ianz.eq.1) then 
         IF (0.le.kpara.and.kpara.le.MAXPAR) then 
            WRITE (zeile (ikl - 1:ikl + 13) , '(e15.8e2)') rpara ( kpara)                                                
            zeile (ikl + 10:ikl + 10) = 'e' 
         ELSE 
            ier_num = - 8 
            ier_typ = ER_FORT 
         ENDIF 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ENDIF 
!                                                                 
ELSEIF (lcomm.eq.3) then 
!                                                                 
   IF (string (ikl - 3:ikl - 1) .eq.'res') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string (1:ikl - lcomm - 1)                                      
         IF (0.le.kpara.and.kpara.le.MAXPAR_RES) then 
            WRITE (zeile (ikl - 3:ikl + 13) , '(e15.8e2)') res_para (kpara)                                      
            zeile (ikl + 8:ikl + 8) = 'e' 
         ELSE 
            ier_num = - 8 
            ier_typ = ER_FORT 
         ENDIF 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSE 
      ier_num = - 2 
      ier_typ = ER_FORT 
   ENDIF 
!                                                                 
ELSE 
   ier_num = - 2 
   ier_typ = ER_FORT 
ENDIF 
IF (ier_num.eq.0) then 
   ll = laenge+15 - ltyp - (iklz - ikl + 1) 
   IF (iklz + 1.le.laenge) zeile (ikl + 14:ll) = string (iklz + 1: laenge)                                                        
   string = zeile 
!ELSE 
!   WRITE ( *, * ) string 
ENDIF 
CALL rem_bl (string, ll) 
END SUBROUTINE suite_ersetz_para                    
!*****7*****************************************************************
SUBROUTINE suite_upd_para (ctype, ww, maxw, wert, ianz) 
!-                                                                      
!       updates the parameter spezified by ctype, index ww  to the      
!       new value of wert                                               
!+                                                                      
!USE allocate_appl
USE errlist_mod 
USE param_mod 
!
IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER (LEN=* ), INTENT(IN   )    :: ctype 
INTEGER           , INTENT(IN   )    :: maxw
INTEGER           , INTENT(IN   )    :: ianz 
INTEGER           , INTENT(IN   )    :: ww (maxw)
REAL              , INTENT(IN   )    :: wert 
!
IF (ctype.eq.'i') then 
   IF (ianz.eq.1) then 
      IF (0.le.ww (1) .and.ww (1) .le.MAXPAR) then 
         inpara (ww (1) ) = int (wert) 
      ELSE 
         ier_num = - 8 
         ier_typ = ER_FORT 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'r') then 
   IF (ianz.eq.1) then 
      IF (0.le.ww (1) .and.ww (1) .le.MAXPAR) then 
         rpara (ww (1) ) = wert 
      ELSE 
         ier_num = - 8 
         ier_typ = ER_FORT 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'res') then 
   IF (ianz.eq.1) then 
      IF (0.le.ww (1) .and.ww (1) .le.MAXPAR_RES) then 
         res_para (ww (1) ) = wert 
      ELSE 
         ier_num = - 8 
         ier_typ = ER_FORT 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSE 
   ier_num = - 2 
   ier_typ = ER_FORT 
   WRITE ( *, * ) ctype 
ENDIF 
! 2000 FORMAT  (' Integer Parameter: ',I1,' : ',i15) 
! 2010 FORMAT  (' Real    Parameter: ',I1,' : ',e15.8e2) 
!
END SUBROUTINE suite_upd_para                       
!*****7***************************************************************  
SUBROUTINE suite_calc_intr_spec (string, line, ikl, iklz, ww, laenge, lp)                                                               
!-                                                                      
!     These are special intrinsic function for the DISCSU_SUITE. Any          
!     intrinsic function that references crystallographic values        
!     is found in this subroutine.                                      
!     Currently empty, needed for formal reasons.
!+                                                                      
!
USE calc_expr_mod
USE errlist_mod 
USE ersetz_mod
USE param_mod 
IMPLICIT none 
!                                                                       
!                                                                       
INTEGER, PARAMETER   :: maxw = 9
!                                                                       
CHARACTER (LEN= * ), INTENT(INOUT) :: string
CHARACTER (LEN= * ), INTENT(INOUT) :: line 
INTEGER            , INTENT(IN   ) :: ikl
INTEGER            , INTENT(IN   ) :: iklz
INTEGER            , INTENT(INOUT) :: laenge
INTEGER            , INTENT(INOUT) :: lp
REAL               , INTENT(INOUT) :: ww
!
INTEGER              :: i, lcomm
REAL                 :: werte (maxw)
!                                                                       
INTEGER              :: length_com 
!                                                                       
lcomm = length_com (string(1:lp), ikl) 
ier_num = - 1 
ier_typ = ER_FORT 
DO i = 1, maxw 
   werte (i) = 0.0 
ENDDO 
!                                                                 
IF (lcomm.eq.0) then 
   CALL ersetz2 (string, ikl, iklz, ww, 0, laenge) 
ELSE 
   ier_num = - 3 
   ier_typ = ER_FORT 
ENDIF 
!                                                                 
IF (ier_num.ne.0) then 
   WRITE ( *, * ) string 
   WRITE ( *, * ) line 
ENDIF 
!                                                                       
END SUBROUTINE suite_calc_intr_spec                 
!*****7**************************************************************** 
SUBROUTINE suite_validate_var_spec (zeile, lp) 
!-                                                                      
!       checks whether the variable name is legal, DIFFEV specific part 
!                                                                       
!       Version : 1.0                                                   
!       Date    : 22 Oct 03                                             
!                                                                       
!       Author  : R.B. Neder  (reinhard.neder@krist.uni-erlangen.de)    
!+                                                                      
!
USE errlist_mod 
IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER (LEN= * ), INTENT(IN   ) :: zeile 
INTEGER            , INTENT(IN   ) :: lp 
!                                                                       
INTEGER, PARAMETER                       :: reserved_n = 1
                                                                        
CHARACTER (LEN=12),DIMENSION(reserved_n) :: reserved = &
              (/'            '/)
INTEGER                                  :: i 
!                                                                       
!                                                                       
ier_num = 0 
ier_typ = ER_NONE 
!                                                                       
DO i = 1, reserved_n 
   IF (index (reserved (i), zeile (1:lp) ) .ne.0) then 
      ier_num = - 25 
      ier_typ = ER_FORT 
   ENDIF 
ENDDO 
!                                                                       
END SUBROUTINE suite_validate_var_spec              
