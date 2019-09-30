MODULE do_string_alloc_mod
!
CONTAINS
!
!****7***************************************************************** 
!
SUBROUTINE do_string_alloc (line, indxg, length) 
!-                                                                      
!     Evaluates the parameters and stores the result                    
!     in the proper string variable                                     
!                                                                       
USE ber_params_mod
USE do_variable_mod
      USE build_name_mod
!     USE calc_intr_mod
      USE errlist_mod 
      USE get_params_mod
      USE param_mod 
      USE set_sub_generic_mod
USE precision_mod
!
      IMPLICIT none 
!                                                                       
      INTEGER, PARAMETER :: maxw= 10 
!                                                                       
      CHARACTER (LEN= * ), INTENT(INOUT) :: line 
      INTEGER            , INTENT(IN)    :: indxg
      INTEGER            , INTENT(INOUT) :: length
!
      CHARACTER(LEN=1024)                  :: zeile
      CHARACTER(LEN=1024), DIMENSION(maxw) :: cpara
      CHARACTER(LEN=1024)                  :: string
!                                                                       
      INTEGER, DIMENSION(maxw)   :: lpara (maxw) 
      INTEGER                    :: i, ikk, ianz, lll 
      INTEGER, DIMENSION(1:maxw) :: iii = 0
      INTEGER                    :: ising , indx_ind, indx_len, indx_env, indx_cwd
      INTEGER                    :: indx_par, indx_isv, indx_ise, indx_fmd
      INTEGER                    :: l_string 
      INTEGER                    :: ikl, iklz, ll, laenge
      INTEGER                    :: ios  ! I/O error status
!                                                                       
REAL(KIND=PREC_DP)    :: wert
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
!                                                                       
!     for flexibility                                                   
!                                                                       
      ising = INDEX (line, '''') 
      DO while (ising.gt.0) 
      line (ising:ising) = '"' 
      ising = INDEX (line, '''') 
      ENDDO 
!                                                                       
!     Get the expression                                                
!                                                                       
      lll = - (length - (indxg + 1) + 1) 
      CALL get_params (line (indxg + 1:length), ianz, cpara, lpara, maxw, lll)
      IF (ier_num.ne.0) then 
         RETURN 
      ELSEIF (ianz.eq.0) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!     locate first ""
      ising=INDEX(line (indxg + 1:length),'"')
      indx_ind=INDEX(line (indxg + 1:length),'index')   ! locate index function
      indx_len=INDEX(line (indxg + 1:length),'length')  ! locate length function
      indx_env=INDEX(line (indxg + 1:length),'getenv')  ! locate length function
      indx_cwd=INDEX(line (indxg + 1:length),'getcwd')  ! locate length function
      indx_isv=INDEX(line (indxg + 1:length),'isvar')   ! locate isvar function
      indx_ise=INDEX(line (indxg + 1:length),'isexp')   ! locate isexp function
      indx_fmd=INDEX(line (indxg + 1:length),'fmodt')   ! locate fmodt function
      indx_par=INDEX(line (indxg + 1:length),'par_name')  ! locate par_name function
!
      IF((indx_ind>0 .AND. indx_ind<ising)  .OR. & ! We got a function of a string argument
         (indx_len>0 .AND. indx_len<ising)  .OR. & ! We got a function of a string argument
         (indx_env>0 .AND. indx_env<ising)  .OR. & ! We got a function of a string argument
         (indx_cwd>0 .AND. indx_cwd<ising)  .OR. & ! We got a function of a string argument
         (indx_isv>0 .AND. indx_isv<ising)  .OR. & ! We got a function of a string argument
         (indx_ise>0 .AND. indx_ise<ising)  .OR. & ! We got a function of a string argument
         (indx_fmd>0 .AND. indx_fmd<ising)  .OR. & ! We got a function of a string argument
         (indx_par>0 .AND. indx_par<ising)) THEN
         string = line (indxg + 1:length)
         laenge = length - indxg
         ikl = INDEX(string,'(')
         iklz= INDEX(string,')',.TRUE.)
         zeile = string (ikl + 1:iklz - 1) 
         ll = iklz - ikl - 1 
         CALL calc_intr (string, zeile, ikl, iklz, laenge, ll) 
         READ(string(1:LEN_TRIM(string)),*,IOSTAT=ios) wert
         IF(ios/=0) THEN
            ier_msg(1) = string(1:42)
            ier_num = -6
            ier_typ = ER_FORT
            RETURN
         ENDIF
      ELSE
!                                                                       
!     Construct the regular string, not a function
!                                                                       
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         string = cpara (1) (1:lpara (1) ) 
         l_string = lpara (1) 
         res_para (0) = 1 
         res_para (1) = lpara (1) 
      ENDIF
      IF(indxg==0) THEN
         line = string
         RETURN
      ENDIF
      IF (ier_num.eq.0) then 
!                                                                       
!-----evaluate the INDEX of the variable                                
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
                     IF (ier_num.eq.0) then 
                        IF (ianz.ge.1.or.ianz.le.3) then 
                           CALL ber_params (ianz, cpara, lpara, werte, maxw)
                           IF (ier_num.eq.0) then 
                              DO i = 1, ianz 
                              iii (i) = nint (werte (i) ) 
                              ENDDO 
!                                                                       
!     ------------Store result in the variable                          
!                                                                       
                              CALL p_upd_para (line (1:ikk - 1), iii, ianz, wert, ianz, string)
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
               CALL upd_variable (line (1:i), i, wert, string, l_string) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
END SUBROUTINE do_string_alloc                
!
!****7***************************************************************** 
!
END MODULE do_string_alloc_mod
