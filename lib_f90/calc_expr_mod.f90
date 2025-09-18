MODULE calc_expr_mod
!
private
!
public do_math
!
CONTAINS
!
!****7***************************************************************** 
!
SUBROUTINE do_math (line, indxg, length) 
!-                                                                      
!     Calculates the value of an expression and stores the result       
!     in the proper variable                                            
!                                                                       
USE blanks_mod
USE berechne_mod
USE ber_params_mod
USE constants_mod
USE do_variable_mod
use do_replace_expr_mod
USE errlist_mod 
USE do_string_alloc_mod
USE get_params_mod
USE precision_mod
USE set_sub_generic_mod
USE variable_mod
!
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: maxw = 10
character(len=128)  :: CLOWER_LIMIT = '-2147483647  '  ! == -2**32
character(len=128)  :: CUPPER_LIMIT = ' 2147483647  '  ! ==  2**32
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: line 
INTEGER          , INTENT(IN   ) :: indxg
INTEGER          , INTENT(INOUT) :: length
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line)))   :: zeile, cpara (maxw) , cdummy
character(len=PREC_STRING)                  :: line_expression       ! If EXPR[] is used
INTEGER, DIMENSION(3) :: var_is_type
!                                                                       
INTEGER, DIMENSION(MAXW) :: lpara
INTEGER :: i, ikk,             ianz, lll , laenge
integer :: j
integer, dimension(MAXW) :: lower_limit     ! Lower limit for variable range 
integer, dimension(MAXW) :: upper_limit     ! Upper limit for variable range 
INTEGER :: indxb    ! Location of an "["
INTEGER :: indxp    ! Location of an "("
INTEGER, DIMENSION(2) :: substr = (/0,VAR_CLEN/)
!
logical :: lrange    ! User provided a range of indices
logical :: lexpr     ! User provided n expression on the right side
!                                                                       
REAL(KIND=PREC_DP) ::  wert
REAL(KIND=PREC_DP), DIMENSION(MAXW) ::  werte
!
write(CLOWER_LIMIT, '(i128)') -huge(lower_limit)
write(CUPPER_LIMIT, '(i128)')  huge(upper_limit)
lrange = .FALSE.
cdummy = ' '
line_expression = ' '
lexpr           = .false.
lll = indxg -1
indxb = INDEX(line(1:lll),'[')
indxp = INDEX(line(1:lll),'(')
IF(indxb >0) THEN
   lll = indxb - 1
ELSE
   IF(indxp > 0) lll = indxp - 1
ENDIF

CALL p_get_var_type(line, lll, var_is_type)
IF(var_is_type(1)== IS_UNKNOWN) THEN     ! unknown variable name on left side
   ier_num = -2
   ier_typ = ER_FORT
   ier_msg(1) = 'Offending name is '//line(1:lll)
   RETURN
elseif(var_is_type(1)== IS_EXPR) then
   call do_set_expression (line, indxg, length)
   return
ENDIF
IF(var_is_type(3)==IS_READ) THEN
   ier_num = -41
   ier_typ = ER_FORT
   ier_msg(1) = 'Offending name is '//line(1:lll)
   RETURN
ENDIF
!
!     String substitution???
!
zeile = line(indxg+1:length)
laenge = length - indxg
CALL rem_leading_bl(zeile,laenge)
IF(zeile(1:1)=='"' .OR. zeile(1:1)=='''') THEN    ! This is a character substitution
!IF (INDEX (line, '"') .gt.0.or.INDEX (line, '''') .gt.0) then 
   IF(var_is_type(1)==IS_CHAR) THEN
      CALL do_string_alloc (line, indxg, length) 
      RETURN 
   ELSE
      ier_num = -42
      ier_typ = ER_FORT
      ier_msg(1) = 'Offending name is '//line(1:lll)
      RETURN 
   ENDIF 
ENDIF 
!                                                                       
!     Get the expression                                                
!                                                                       
lll = length - (indxg + 1) + 1 
CALL get_params (line (indxg + 1:length), ianz, cpara, lpara, maxw, lll)
IF (ier_num.ne.0) then 
   RETURN 
ELSEIF (ianz.eq.0) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
zeile = '('//cpara(1)(1:lpara(1)) //')' 
i     = lpara(1) + 2 
!                                                                       
!     Calculate the expression                                          
!                                                                       
IF(var_is_type(1)/=IS_CHAR) THEN     ! Variable on left side is numeric
!write(*,*) ' ZEILE RECHTS >>', zeile(1:len_trim(zeile))
   if(index(zeile(1:len_trim(zeile)),'EXPR')>0) then    ! something like ...= EXPR[index]
      line_expression = zeile(1:len_trim(zeile))
      lexpr           = .true.
      call do_replace_expr(zeile, i)
   endif
      wert = berechne (zeile, i) 
      IF (ier_num.eq.0) then 
!                                                                       
!-----evaluate the index of the variable                                
!                                                                       
         lll = indxg - 1 
         CALL get_params (line (1:indxg - 1), ianz, cpara, lpara, maxw, lll)
         IF (ier_num.eq.0) then 
            line = cpara (1) 
            i = lpara (1) 
            ikk = INDEX (line, '[') 
            IF (ikk.lt.i.and.ikk.gt.0) then 
               IF (line (i:i) .eq.']') then 
                  IF (i.gt.ikk + 1) then 
                     zeile = ' ' 
                     zeile (1:i - ikk - 1) = line (ikk + 1:i - 1) 
                     lll = i - ikk - 1 
                     CALL get_params (zeile, ianz, cpara, lpara, maxw, lll)
!write(*,*) ' zeile ', zeile(1:len_trim(zeile)), ianz, cpara(1)(1:lpara(1)), ier_num, ier_typ
                     lrange = .false.
                     if(any(index(cpara(1:ianz),':')>0)) then    ! At least one ':'
                        lrange = .true.
                        do i = ianz, 1, -1          ! Loop over all possible dimensions
                           j = index(cpara(i), ':')
                           if(j>0) then
                              if(j<lpara(i)) then 
                                 cpara(2*i) = cpara(i)(j+1:)
                              else
                                 cpara(2*i) = CUPPER_LIMIT
                              endif
                              if(j>1       ) then 
                                 cpara(2*i-1) = cpara(i)(1:j-1)
                              else
                                 cpara(2*i-1) = CLOWER_LIMIT
                              endif
                           else
                              cpara(2*i  ) = cpara(i)
                              cpara(2*i-1) = cpara(i)
                           endif
                           lpara(2*i  ) = len_trim(cpara(2*i  ))
                           lpara(2*i-1) = len_trim(cpara(2*i-1))
                        enddo
                     endif
                     if(lrange) ianz = 2*ianz
!do i = 1, ianz                 ! Loop over all possible dimensions
!write(*,*) 'PARA ', i, cpara(i)(1:len_trim(cpara(i)))
!enddo
                     IF (ier_num.eq.0) then 
                        IF((lrange .and.(ianz >= 2 .and. ianz <= 6)) .or.  & ! then 
                                        (ianz.ge.1 .and. ianz.le.3)      ) then 
                           CALL ber_params (ianz, cpara, lpara, werte, maxw)
                           IF (ier_num.eq.0) then 
!                             DO i = 1, ianz 
!                             iii (i) = nint (werte (i) ) 
!                             ENDDO 
                              lower_limit = -huge(lower_limit)
                              upper_limit =  huge(upper_limit)
                              if(lrange) then
                                 lower_limit(1:ianz/2) = nint(werte(1:ianz:2))
                                 upper_limit(1:ianz/2) = nint(werte(2:ianz:2))
                                 ianz = ianz / 2           ! Back to number of dimensions
                              else
                                 lower_limit(1:ianz  ) = nint(werte(1:ianz  ))
                                 upper_limit(1:ianz  ) = nint(werte(1:ianz  ))
                              endif
!write(*,'(a, i3, 2i20)') 'HUGE ', 0, -huge(i), huge(i) 
!do i = 1, 3                 ! Loop over all possible dimensions
!write(*,'(a, i3, 2i20)') 'LIMI ', i, lower_limit(i), upper_limit(i)
!enddo
!write(*,*) ' wert ', wert
!                                                                       
!     ------------Store result in the variable                          
!                                                                       
                              CALL p_upd_para(line(1:ikk - 1), lower_limit, upper_limit, ianz, lrange, wert, ianz, cdummy, substr, lexpr, line_expression)
                           ENDIF 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_FORT 
                        ENDIF 
                     ENDIF 
                  ELSE 
                     ier_num = - 10 
                     ier_typ = ER_FORT 
                  ENDIF 
               ELSE 
                  ier_num = - 9 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSEIF (ikk.eq.0) then 
               CALL upd_variable (line (1:i), i, wert, cpara (1), lpara (1), substr )
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
            ENDIF 
         ENDIF 
      ENDIF 
ELSE             ! Variable on left side is CHARACTER
         zeile = line (1:indxg) //'"%c",'//line (indxg + 1:length) 
         length = length + 5 
         CALL do_string_alloc (zeile, indxg, length) 
ENDIF 
!                                                                       
END SUBROUTINE do_math                        
!
!*****7**************************************************************** 
!
END MODULE calc_expr_mod
