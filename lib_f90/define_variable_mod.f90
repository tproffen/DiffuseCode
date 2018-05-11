MODULE define_variable_mod
!
CONTAINS
!
!*****7**************************************************************** 
!
SUBROUTINE define_variable (zeile, lp, is_diffev) 
!-                                                                      
!       Allows the user to define variables                             
!                                                                       
!       Version : 1.0                                                   
!       Date    : 22 Oct 03                                             
!                                                                       
!       Author  : R.B. Neder  (reinhard.neder@mail.uni-wuerzburg.de)    
!+                                                                      
USE ber_params_mod
USE build_name_mod
USE do_variable_mod
USE errlist_mod 
USE get_params_mod
USE prompt_mod 
USE variable_mod
!
IMPLICIT none 
!
CHARACTER (LEN=* ), INTENT(INOUT) :: zeile 
INTEGER           , INTENT(INOUT) :: lp 
LOGICAL           , INTENT(IN)    :: is_diffev
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 3) 
      LOGICAL , PARAMETER :: no_diffev = .FALSE.
                                                                        
      INTEGER ianz 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      REAL werte (maxw) 
!                                                                       
      CHARACTER(1024) c_type, c_temp, c_init 
      INTEGER l_type, l_temp 
      INTEGER ccc_type 
      INTEGER i, j 
      LOGICAL l_init 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      ccc_type = 0
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) RETURN 
!                                                                       
      IF (str_comp (cpara (1) , 'real', 3, lpara (1) , 4) .or.str_comp (&
      cpara (1) , 'inte', 2, lpara (1) , 4) .or.str_comp (cpara (1) ,   &
      'char', 2, lpara (1) , 4) ) THEN                                  
!                                                                       
!     --A new variable is being defined                                 
!                                                                       
         IF (ianz.eq.2.or.ianz.eq.3) THEN 
            IF (var_num.lt.VAR_MAX) THEN 
!                                                                       
!     ----- If a free slot is available validate the name against       
!     ----- illegal names like "cos", "sin" etc.                        
!                                                                       
               CALL validate_variable (cpara (2), lpara (2) ) 
               IF (ier_num.ne.0) RETURN 
!                                                                       
!     ----- temporarily store the variable type and name and evaluate   
!     ----- the optional initialising parameter                         
!                                                                       
               c_type = cpara (1) 
               l_type = lpara (1) 
               c_temp (1:lpara (2) ) = cpara (2) (1:lpara (2) ) 
               l_temp = lpara (2) 
               werte (1) = 0.0 
               IF (str_comp (c_type, 'real', 3, l_type, 4) ) THEN 
                  ccc_type = VAR_TYPE_REAL 
               ELSEIF (str_comp (c_type, 'inte', 2, l_type, 4) ) THEN 
                  ccc_type = VAR_TYPE_INTE 
               ELSEIF (str_comp (c_type, 'char', 2, l_type, 4) ) THEN 
                  ccc_type = VAR_TYPE_CHAR 
               ENDIF 
               l_init = .false. 
               c_init = ' ' 
               IF (ianz.eq.3) THEN 
                  CALL del_params (2, ianz, cpara, lpara, maxw) 
                  IF (ier_num.ne.0) RETURN 
                  IF (ccc_type.eq.VAR_TYPE_CHAR) THEN 
                     CALL do_build_name (ianz, cpara, lpara, werte,     &
                     maxw, 1)                                           
                     c_init = cpara (1) (1:lpara (1) ) 
                  ELSE 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  ENDIF 
                  l_init = .true. 
                  IF (ier_num.ne.0) RETURN 
               ENDIF 
!                                                                       
!     ----- Make sure the variable name has not yet been defined as     
!           other variable type.                                        
!           And initialisation value is not used on old variables       
!                                                                       
               ier_num = 0 
               ier_typ = ER_NONE 
               DO i = 1, var_num 
               IF (c_temp (1:l_temp) .eq.var_name (i) ) THEN 
                  IF (ccc_type.ne.var_type (i) ) THEN 
                     ier_num = - 32 
                     ier_typ = ER_FORT 
                  ELSE 
                     IF (l_init) THEN 
                        ier_num = - 33 
                        ier_typ = ER_FORT 
                     ELSE 
                        RETURN 
                     ENDIF 
                  ENDIF 
               ENDIF 
               ENDDO 
               IF (ier_num.ne.0) THEN 
                  CALL errlist 
                  RETURN 
               ENDIF 
!                                                                       
!     ----- sort the new variable name in descending length and         
!           descending alphabetical order                               
!                                                                       
               i = var_sys + 1 
               DO while (l_temp.lt.var_l (i) .and.i.le.var_num) 
               i = i + 1 
               ENDDO 
               DO while (l_temp.eq.var_l (i) .and.llt (c_temp, var_name &
               (i) ) .and.i.le.var_num)                                 
               i = i + 1 
               ENDDO 
               DO j = var_num, i, - 1 
                  var_name (j + 1) = var_name (j) 
                  var_l    (j + 1) = var_l    (j) 
                  var_type (j + 1) = var_type (j) 
                  var_val  (j + 1) = var_val  (j) 
                  var_char (j + 1) = var_char (j) 
                  var_diff (j + 1) = var_diff (j)    ! true if refine paaram from diffev
               ENDDO 
!                                                                       
!     ----- found the proper slot, store value, name and type           
!                                                                       
               var_num = var_num + 1 
               var_name (i) (1:l_temp) = c_temp 
               var_l (i) = l_temp 
               IF (str_comp (c_type, 'real', 3, l_type, 4) ) THEN 
                  var_type (i) = VAR_TYPE_REAL 
                  var_val (i) = werte (1) 
                  var_diff (i) = is_diffev      ! true if refine param from diffev
               ELSEIF (str_comp (c_type, 'inte', 2, l_type, 4) ) THEN 
                  var_type (i) = VAR_TYPE_INTE 
                  var_val (i) = nint (werte (1) ) 
                  var_diff (i) = is_diffev      ! true if refine param from diffev
               ELSEIF (str_comp (c_type, 'char', 2, l_type, 4) ) THEN 
                  var_type (i) = VAR_TYPE_CHAR 
                  var_val (i) = 0.0 
                  var_char (i) = c_init (1:len(var_char))
                  var_diff (i) = is_diffev      ! true if refine param from diffev
               ENDIF 
            ELSE 
               ier_num = - 23 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSEIF (str_comp (cpara (1) , 'show', 2, lpara (1) , 4) ) THEN 
         CALL show_variables 
      ELSEIF (str_comp (cpara (1) , 'reset', 5, lpara (1) , 5) ) THEN 
         CALL rese_variables (no_diffev)
      ELSEIF (str_comp (cpara (1) , 'delete', 3, lpara (1) , 6) ) THEN 
         CALL del_variables(MAXW, ianz, cpara, lpara, no_diffev)
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
END SUBROUTINE define_variable                
!
!*****7**************************************************************** 
END MODULE define_variable_mod
