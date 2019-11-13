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
USE constants_mod
USE do_variable_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE prompt_mod 
USE take_param_mod
USE variable_mod
!
IMPLICIT none 
!
CHARACTER (LEN=* ), INTENT(INOUT) :: zeile 
INTEGER           , INTENT(INOUT) :: lp 
LOGICAL           , INTENT(IN)    :: is_diffev
!                                                                       
INTEGER, PARAMETER :: maxw = 4
LOGICAL, PARAMETER :: no_diffev = .FALSE.
                                                                        
INTEGER :: ianz 
!                                                                       
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte
!
INTEGER, PARAMETER :: MAXF=2
CHARACTER(LEN=1024), DIMENSION(MAXF) :: ccpara
INTEGER            , DIMENSION(MAXF) :: llpara
REAL(KIND=PREC_DP) , DIMENSION(MAXF) :: wwerte
!                                                                       
CHARACTER(LEN=1024) :: c_type, c_temp, c_init , string
INTEGER :: l_type, l_temp 
INTEGER :: ccc_type 
INTEGER :: i, j , length, iianz
INTEGER :: n1, n2, n_data    ! Dimensions of arrays, total size
INTEGER :: place             ! location of variable arrays
LOGICAL :: l_init 
!                                                                       
LOGICAL :: str_comp 
!
INTEGER, PARAMETER :: NOPTIONAL = 1
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
DATA oname  / 'dim' /
DATA loname /  3    /
!
opara  =  (/ 'scalar' /)   ! Always provide fresh default values
lopara =  (/  6       /)
owerte =  (/  0.0     /)
!                                                                       
ccc_type = 0
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) RETURN 
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                              oname, loname, opara, lopara, lpresent, owerte)
!
n1 = 0
n2 = 0
n_data = n1 * n2
IF(opara(1) == 'scalar') THEN
   n1 = 0
   n2 = 0
   n_data = 0
ELSEif(opara(1)(1:1) == '[' .AND. opara(1)(lopara(1):lopara(1)) == ']') THEN
   string = opara(1)(2:lopara(1)-1)
   length = lopara(1)-2
   ccpara(:) = ' '
   llpara(:) = 0
   wwerte(:) = 0.0
   CALL get_params (string, iianz, ccpara, llpara, MAXF, length)
   IF(ier_num /= 0) THEN
      ier_msg(1) = 'Incorrect ''dim:[]'' parameter'
      ier_msg(2) = 'Variables can only be arrays with'
      ier_msg(3) = 'one or two dimensions '
      RETURN
   ENDIF
   CALL ber_params (iianz, ccpara, llpara, wwerte, MAXF)
   IF(ier_num /= 0) THEN
      ier_msg(1) = 'Incorrect ''dim:[]'' parameter'
      ier_msg(2) = 'Variables can only be arrays with'
      ier_msg(3) = 'one or two dimensions '
      RETURN
   ENDIF
   IF(iianz==1) THEN
      n1 = NINT(wwerte(1))
      n2 = 1
   ELSEIF(iianz==2) THEN
      n1 = NINT(wwerte(1))
      n2 = NINT(wwerte(2))
   ELSE
      ier_num = -40
      ier_typ = ER_FORT
      ier_msg(1) = 'Incorrect ''dim:[]'' parameter'
      ier_msg(2) = 'Variables can only be arrays with'
      ier_msg(3) = 'one or two dimensions '
      RETURN
   ENDIF
   n_data = n1 * n2
ENDIF
!                                                                       
iftype: IF(str_comp(cpara (1) , 'real', 3, lpara (1) , 4) .or.  &
   str_comp(cpara (1) , 'inte', 2, lpara (1) , 4) .or.  &
   str_comp(cpara (1) , 'char', 2, lpara (1) , 4) ) THEN                                  
!                                                                       
!     --A new variable is being defined                                 
!                                                                       
   ifianz: IF (ianz.eq.2.or.ianz.eq.3) THEN 
      ifvarnum: IF (var_num.lt.VAR_MAX) THEN 
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
            ccc_type =       IS_REAL 
         ELSEIF (str_comp (c_type, 'inte', 2, l_type, 4) ) THEN 
            ccc_type = IS_INTE 
         ELSEIF (str_comp (c_type, 'char', 2, l_type, 4) ) THEN 
            ccc_type =       IS_CHAR 
         ENDIF 
         l_init = .false. 
         c_init = ' ' 
         IF (ianz.eq.3) THEN 
            CALL del_params (2, ianz, cpara, lpara, maxw) 
            IF (ier_num.ne.0) RETURN 
            IF (ccc_type.eq.      IS_CHAR) THEN 
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
                  ier_msg(1) = c_temp(1:l_temp)
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
         DO WHILE (l_temp.lt.var_l (i) .AND.i.le.var_num) 
            i = i + 1 
         ENDDO 
         DO  WHILE(l_temp.eq.var_l (i) .AND.LLT(c_temp, var_name(i))  &
                                       .AND.i.le.var_num)                                 
            i = i + 1 
         ENDDO 
         DO j = var_num, i, - 1 
            var_name (j+1) = var_name (j)    ! Shift fields by size of new variable
            var_l    (j+1) = var_l    (j) 
            var_entry(j+1) = var_entry(j)
            var_type (j+1) = var_type (j) 
            var_val  (j+1) = var_val  (j) 
            var_char (j+1) = var_char (j) 
            var_diff (j+1) = var_diff (j)    ! true if refine param from diffev
         ENDDO 
!                                                                       
!     ----- found the proper slot, store value, name and type           
!                                                                       
         var_num = var_num + 1 
         var_name (i) (1:l_temp) = c_temp 
         var_l (i) = l_temp 
!
         place = 1
         IF(n_data>0) THEN             ! We have an array
            search_entry: DO j=1,VAR_MAX
               IF(var_field(j)%var_shape(1)==0) THEN   !Found free entry
                  place = j
                  var_n_arr = MAX(var_n_arr, place)
                  EXIT search_entry
               ENDIF
            ENDDO search_entry
            var_entry(i) = place       ! Keep track in which entry the array is stored
            var_field(place)%var_shape(1) = n1
            var_field(place)%var_shape(2) = n2
            IF(ALLOCATED(var_field(place)%var_value)) DEALLOCATE(var_field(place)%var_value)
            ALLOCATE(var_field(place)%var_value(n1,n2))
            IF(ALLOCATED(var_field(place)%var_char)) DEALLOCATE(var_field(place)%var_char)
            ALLOCATE(var_field(place)%var_char(n1,n2))
         ELSE
            var_entry(i) = 0
         ENDIF
         IF (str_comp (c_type, 'real', 3, l_type, 4) ) THEN 
            var_type (i) =       IS_REAL 
            var_val  (i) = werte (1) 
            IF(n_data>0) THEN             ! We have an array
                var_field(place)%var_value(:,:) = werte(1)
            ENDIF
            var_diff (i) = is_diffev      ! true if refine param from diffev
         ELSEIF (str_comp (c_type, 'inte', 2, l_type, 4) ) THEN 
            var_type (i) =       IS_INTE 
            var_val  (i) = NINT(werte (1) ) 
            IF(n_data>0) THEN             ! We have an array
                var_field(place)%var_value(:,:) = NINT(werte(1))
            ENDIF
            var_diff (i) = is_diffev      ! true if refine param from diffev
         ELSEIF (str_comp (c_type, 'char', 2, l_type, 4) ) THEN 
            var_type (i) =       IS_CHAR 
            var_val  (i) = 0.0 
            var_char (i) = c_init (1:len(var_char))
            IF(n_data>0) THEN             ! We have an array
                var_field(place)%var_char(:,:) = ' '
            ENDIF
            var_diff (i) = is_diffev      ! true if refine param from diffev
         ENDIF 
      ELSE  ifvarnum
         ier_num = - 23 
         ier_typ = ER_FORT 
      ENDIF  ifvarnum
   ELSE  ifianz
      ier_num = - 6 
      ier_typ = ER_COMM 
   ENDIF  ifianz
ELSEIF (str_comp (cpara (1) , 'show', 2, lpara (1) , 4) ) THEN iftype
   CALL show_variables 
ELSEIF (str_comp (cpara (1) , 'reset', 5, lpara (1) , 5) ) THEN iftype
   CALL rese_variables (no_diffev)
ELSEIF (str_comp (cpara (1) , 'delete', 3, lpara (1) , 6) ) THEN iftype
   CALL del_variables(MAXW, ianz, cpara, lpara, no_diffev)
ELSE iftype
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF iftype
!                                                                       
END SUBROUTINE define_variable                
!
!*****7**************************************************************** 
!
SUBROUTINE def_set_variable(v_type, v_name, v_value, IS_DIFFEV)
!
USE calc_expr_mod
USE errlist_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: v_type
CHARACTER(LEN=*), INTENT(IN) :: v_name
REAL            , INTENT(IN) :: v_value
LOGICAL         , INTENT(IN) :: IS_DIFFEV
!
CHARACTER(LEN=1024) :: string
INTEGER             :: length
INTEGER             :: indxg
!
string = v_type(1:LEN_TRIM(v_type)) // ', ' // v_name(1:LEN_TRIM(v_name))
length = LEN_TRIM(string)
CALL define_variable(string, length, IS_DIFFEV)           ! Define as user variable
WRITE(string,'(a,a,G20.8E3)') v_name(1:LEN_TRIM(v_name)), ' = ', v_value
indxg  = INDEX(string, '=')
length = LEN_TRIM(string)
CALL do_math (string, indxg, length)
!
END SUBROUTINE def_set_variable
!
!*****7**************************************************************** 
!
END MODULE define_variable_mod
