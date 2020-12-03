!MODULE calc_intr_mod
!
!CONTAINS
!
!*****7**************************************************************** 
!
RECURSIVE      SUBROUTINE calc_intr (string, line, ikl, iklz, lll, lp) 
!                                                                       
!     Evaluate all intrinsic functions                                  
!                                                                       
USE build_name_mod
USE ersetz_mod
USE do_read_number_mod
USE errlist_mod 
USE get_params_mod
USE lib_length
USE lib_random_func
USE precision_mod
USE random_mod
USE variable_mod
USE wink_mod
USE times_mod
USE trig_degree_mod
USE set_sub_generic_mod
USE support_mod
!
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: MAXW = 30
!                                                                       
CHARACTER(LEN=*), INTENT(INOUT) :: string
CHARACTER(LEN=*), INTENT(INOUT) :: line 
INTEGER         , INTENT(INOUT) :: ikl
INTEGER         , INTENT(INOUT) :: iklz
INTEGER         , INTENT(INOUT) :: lll
INTEGER         , INTENT(INOUT) :: lp 
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(string))), DIMENSION(MAXW) :: cpara
CHARACTER(LEN=MAX(PREC_STRING,LEN(string))) :: zeile 
CHARACTER(LEN=MAX(PREC_STRING,LEN(string))) :: answer , search
CHARACTER(LEN=24)   :: fmodt
INTEGER, DIMENSION(MAXW)  :: lpara
INTEGER  :: ikom, i, ianz 
INTEGER  :: lcom 
INTEGER  :: ihyp 
INTEGER  :: dummy 
LOGICAL :: BACK   ! FLAG for index intrinsic
REAL  :: fl1, fl2, fl3, gbox_k, gbox_w, gbox_x 
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
REAL(KIND=PREC_DP)  :: ww, a , skew
REAL(KIND=PREC_DP)  :: ww1, ww2
REAL(KIND=PREC_DP), DIMENSION(MAXW)  :: wwerte
!     REAL sind, cosd, tand, asind, acosd, atand 
!     REAL atan2, atan2d 
!                                                                       
!
!                                                                       
      ier_num = -1 
      ier_typ = ER_FORT 
      ikom = INDEX (line, ',') 
      ihyp = MAX(INDEX(line, '''') , INDEX(line, '"') ) 
      IF (ikom.eq.0.and.ihyp.eq.0) then 
         ww = do_read_number (line, lp) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
      ENDIF 
      ier_num = 0 
      ier_typ = ER_NONE 
      lcom = length_com (string, ikl) 
      IF (lcom.eq.0) then 
         CALL ersetz2 (string, ikl, iklz, ww, 0, lll) 
      ELSEIF (lcom.eq.6) then 
         IF (string (ikl - 6:ikl - 1) .eq.'getcwd') then 
            CALL holecwd (zeile, dummy) 
            i = len_str (zeile) 
            CALL ersetzc (string, ikl, iklz, zeile, i, 6, lll) 
         ELSEIF (string (ikl - 6:ikl - 1) .eq.'getenv') then 
            zeile = line (2:lp - 1) 
            IF(line /= ' ' .AND. lp>2) THEN
               CALL holeenv (zeile, answer) 
               i = len_str (answer) 
               CALL ersetzc (string, ikl, iklz, answer, i, 6, lll) 
            ELSE
               ier_num = -6
               ier_typ = ER_FORT
            ENDIF
         ELSEIF (string (ikl - 6:ikl - 1) .eq.'length') then 
            zeile = line (2:lp - 1)
            i = lp - 2
            IF(zeile(1:1)=='''' .and. zeile(i:i)=='''' .OR. & !) THEN
               zeile(1:1)=='"'  .and. zeile(i:i)=='"' ) THEN
               ww = REAL(i-2, PREC_DP) 
               CALL ersetz2 (string, ikl, iklz, ww, 6, lll) 
            ELSE
               CALL get_params (line, ianz, cpara, lpara, maxw, lp) 
               IF(ier_num==0) THEN
                  CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
                     IF(ier_num==0) THEN
                     ww = lpara(1)
                  CALL ersetz2 (string, ikl, iklz, ww, 6, lll) 
                  ENDIF
               ENDIF
            ENDIF
         ELSE 
            CALL p_calc_intr_spec (string, line, ikl, iklz, ww, lll, lp) 
         ENDIF 
      ELSEIF (lcom.eq.5) then 
         IF (string (ikl - 5:ikl - 1) .eq.'asind') then 
            IF (abs (ww) .le.1.0D0) then 
               ww = asind(ww) 
               CALL ersetz2 (string, ikl, iklz, ww, 5, lll) 
            ELSE 
               ier_num = - 7 
               ier_typ = ER_FORT 
            ENDIF 
         ELSEIF (string (ikl - 5:ikl - 1) .eq.'acosd') then 
            IF (abs (ww) .le.1.0D0) then 
               ww = acosd (ww) 
               CALL ersetz2 (string, ikl, iklz, ww, 5, lll) 
            ELSE 
               ier_num = - 7 
               ier_typ = ER_FORT 
            ENDIF 
         ELSEIF (string (ikl - 5:ikl - 1) .eq.'atand') then 
            CALL get_params (line, ianz, cpara, lpara, 2, lp) 
            IF (ianz.eq.1) then 
               ww = atand (ww) 
            ELSEIF (ianz.eq.2) then 
               CALL        eval           (cpara (1), lpara (1) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               CALL        eval           (cpara (2), lpara (2) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               ww1       = do_read_number (cpara (1), lpara (1) ) 
               ww2       = do_read_number (cpara (2), lpara (2) ) 
               ww = atan2d(ww1, ww2)
            ELSE 
               ier_num = - 27 
               ier_typ = ER_FORT 
            ENDIF 
            CALL ersetz2 (string, ikl, iklz, ww, 5, lll) 
         ELSEIF (string (ikl - 5:ikl - 1) .eq.'gskew') then 
            CALL get_params (line, ianz, cpara, lpara, 2, lp) 
            IF (ianz.ge.2) then 
               CALL eval (cpara (1), lpara (1) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               ww1       = do_read_number (cpara (1), lpara (1) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (ww1       .lt.0.0) then 
                  ier_num = - 35 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
!
               CALL eval (cpara (2), lpara (2) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               ww2       = do_read_number (cpara (2), lpara (2) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (ABS(ww2      ) .gt.1.0) then 
                  ier_num = - 36 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF
               skew = ww2      
!                                                                       
               IF (ianz.eq.2.or.cpara (3) .eq.'s') then 
                  a = ww1       
               ELSEIF (ianz.eq.3.and.cpara (3) .eq.'f') then 
                  a = ww1       / sqrt (8.D0 * log (2.D0) ) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
            ww = gasskew (a,skew) 
            CALL ersetz2 (string, ikl, iklz, ww, 5, lll) 
         ELSEIF (string (ikl - 5:ikl - 1) .eq.'fmodt') then 
            CALL get_params (line, ianz, cpara, lpara, maxw, lp) 
            IF(ier_num==0) THEN
               IF(ianz == 0) THEN
                  CALL ersetzc (string, ikl, iklz, f_modt, 24, 5, lll) 
               ELSEIF(ianz == 1 .AND. cpara(1)=='0') THEN
                  CALL ersetzc (string, ikl, iklz, f_modt, 24, 5, lll) 
               ELSE
                  CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
                  IF(ier_num==0) THEN
                     IF(cpara(1)(1:1)=='''' .AND.               &
                        cpara(1)(lpara(1):lpara(1))=='''') THEN
                        zeile = cpara(1)(2:lpara(1)-1)
                     ELSE
                        zeile = cpara(1)(1:lpara(1))
                     ENDIF
                     CALL file_info_disk(zeile,fmodt)
                     CALL ersetzc (string, ikl, iklz, fmodt, 24, 5, lll) 
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF (string (ikl - 5:ikl - 1) .eq.'fdate') then 
            CALL datum 
            CALL ersetzc (string, ikl, iklz, f_date, 24, 5, lll) 
         ELSEIF (string (ikl - 5:ikl - 1) .eq.'index') then 
            CALL get_params (line, ianz, cpara, lpara, maxw, lp) 
            IF(ier_num==0) THEN
               IF(ianz > 1) THEN
                  BACK = .FALSE.
                  IF(cpara(ianz)(1:lpara(ianz))=='BACK') THEN
                     BACK = .TRUE.
                     ianz = ianz -1
                  ENDIF
                  IF(ianz > 1) THEN
                     CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
                     IF(ier_num==0) THEN
                        IF(cpara(1)(1:1)=='''' .AND.               &
                           cpara(1)(lpara(1):lpara(1))=='''') THEN
                           zeile = cpara(1)(2:lpara(1)-1)
                        ELSE
                           zeile = cpara(1)(1:lpara(1))
                        ENDIF
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        IF(ier_num==0) THEN
                           IF(ianz >= 1) THEN
                              CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
                              IF(ier_num==0) THEN
                                 IF(cpara(1)(1:1)=='''' .AND.               &
                                    cpara(1)(lpara(1):lpara(1))=='''') THEN
                                    search= cpara(1)(2:lpara(1)-1)
                                 ELSE
                                    search= cpara(1)(1:lpara(1))
                                 ENDIF
                                 ww = REAL(INDEX(zeile (1:LEN_TRIM(ZEILE )), &
                                                 search(1:LEN_TRIM(search)), &
                                                  BACK ), PREC_DP)
                                 CALL ersetz2 (string, ikl, iklz, ww, 5, lll) 
                              ENDIF
                           ELSE 
                              ier_num = - 6 
                              ier_typ = ER_FORT 
                              RETURN 
                           ENDIF
                        ENDIF
                     ENDIF
                  ELSE
                     ier_num = - 6 
                     ier_typ = ER_FORT 
                     RETURN 
                  ENDIF
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
            ENDIF
   ELSEIF (string (ikl - 5:ikl - 1) .eq.'psvgt') then 
      CALL intr_psvgt(line, lp, ww)
      CALL ersetz2 (string, ikl, iklz, ww, 5, lll) 
   ELSE 
            CALL p_calc_intr_spec (string, line, ikl, iklz, ww, lll, lp) 
   ENDIF 
ELSEIF (lcom.eq.4) then 
         IF (string (ikl - 4:ikl - 1) .eq.'asin') then 
            IF (abs (ww) .le.1.0D0) then 
               ww = asin (ww) 
               CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
            ELSE 
               ier_num = - 7 
               ier_typ = ER_FORT 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'acos') then 
            IF (abs (ww) .le.1.0D0) then 
               ww = acos (ww) 
               CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
            ELSE 
               ier_num = - 7 
               ier_typ = ER_FORT 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'atan') then 
            CALL get_params (line, ianz, cpara, lpara, 2, lp) 
            IF (ianz.eq.1) then 
               ww = atan (ww) 
            ELSEIF (ianz.eq.2) then 
               ww1       = do_read_number (cpara (1), lpara (1) ) 
               ww2       = do_read_number (cpara (2), lpara (2) ) 
               ww = atan2(ww1, ww2)
            ELSE 
               ier_num = - 27 
               ier_typ = ER_FORT 
            ENDIF 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'sind') then 
            ww = sind (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'cosd') then 
            ww = cosd (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'tand') then 
            ww = tand (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'tanh') then 
            ww = tanh (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'sinh') then 
            ww = sinh (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'cosh') then 
            ww = cosh (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'sqrt') then 
            IF (ww.ge.0D0) then 
               ww = sqrt (ww) 
               CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
            ELSE 
               ier_num = - 5 
               ier_typ = ER_FORT 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'nint') then 
            ww = REAL(NINT (ww), PREC_DP ) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'frac') then 
            ww = ww - REAL(INT(ww), PREC_DP) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'gbox') then 
            CALL get_params (line, ianz, cpara, lpara, 3, lp) 
            IF (ianz.eq.3) then 
               DO i = 1, 3 
               CALL eval (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               wwerte(i) = do_read_number (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (wwerte(i) .lt.0.0) then 
                  ier_num = - 35 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
               ENDDO 
!           Determine relative size of exponentials and box             
               fl1 = REAL(0.5D0 * REAL(wwerte(1), KIND=KIND(0.0D0)) * sqrt (zpi), KIND=KIND(0.0E0)) 
               fl2 = wwerte(2) 
               fl3 = REAL(0.5D0 * REAL(wwerte(3), KIND=KIND(0.0D0)) * sqrt (zpi), KIND=KIND(0.0E0)) 
!           Normalize to 1                                              
               gbox_k = 1.D0 / (fl1 + fl2 + fl3) 
               gbox_w = 1.D0 - (fl1 + fl3) * gbox_k 
!           Get random number                                           
               gbox_x = ran1 (idum) 
!                                                                       
               IF (gbox_x.lt.fl1 * gbox_k) then 
                  ww = - wwerte(2) * 0.5D0 - ABS (gasdev (wwerte(1)) ) 
               ELSEIF (gbox_x.lt. (fl1 + fl2) * gbox_k) then 
                  ww = - wwerte(2) * 0.5D0 + (gbox_x - fl1 * gbox_k)      &
                       * wwerte(2) / gbox_w
               ELSE 
                  ww = wwerte(2) * 0.5D0 + ABS (gasdev (wwerte(3) ) ) 
               ENDIF 
               CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'gran') then 
            CALL get_params (line, ianz, cpara, lpara, 2, lp) 
            IF (ianz.ge.1) then 
               CALL eval (cpara (1), lpara (1) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               ww1       = do_read_number (cpara (1), lpara (1) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (ww1       .lt.0.0) then 
                  ier_num = - 35 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
!                                                                       
               IF (ianz.eq.1.or.cpara (2) .eq.'s') then 
                  a = ww1       
               ELSEIF (ianz.eq.2.and.cpara (2) .eq.'f') then 
                  a = ww1       / sqrt (8.D0 * log (2.D0) ) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
            ww = gasdev (a) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'logn') then 
            CALL get_params (line, ianz, cpara, lpara, 2, lp) 
            IF (ianz.ge.2) then 
               CALL eval (cpara (1), lpara (1) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               ww1       = do_read_number (cpara (1), lpara (1) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (ww1       .lt.0.0D0) then 
                  ier_num = - 37 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
!                                                                       
               CALL eval (cpara (2), lpara (2) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               ww2       = do_read_number (cpara (2), lpara (2) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (ww2       .lt.0.0D0) then 
                  ier_num = - 35
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
!                                                                       
               IF (ianz.eq.2.or.cpara (3) .eq.'s') then 
                  a = ww2       
               ELSEIF (ianz.eq.3.and.cpara (3) .eq.'f') then 
                  a = ww2       / sqrt (8.D0 * log (2.D0) ) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
            ww = exp (log (ww1       ) + gasdev (a) ) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'pois') then 
            ww = poidev (REAL(ww), idum) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'date') then 
            CALL datum_intrinsic 
            CALL ersetzc (string, ikl, iklz, f_date, 24, 4, lll) 
         ELSE 
            CALL p_calc_intr_spec (string, line, ikl, iklz, ww, lll, lp) 
         ENDIF 
      ELSEIF (lcom.eq.3) then 
         IF (string (ikl - 3:ikl - 1) .eq.'sin') then 
            ww = sin (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'cos') then 
            ww = cos (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'tan') then 
            ww = tan (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'exp') then 
            ww = exp (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'abs') then 
            ww = abs (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'int') then 
            ww = REAL(INT(ww), PREC_DP ) 
            CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'max') then 
            CALL get_params (line, ianz, cpara, lpara, MAXW, lp) 
            IF (ianz>=1 .AND. ianz <=MAXW) THEN 
               CALL eval(cpara(1), lpara(1) ) 
               IF(ier_num /= 0) THEN 
                  RETURN 
               ENDIF 
               werte(1) = do_read_number(cpara(1), lpara(1) ) 
               IF (ier_num /= 0) THEN 
                  RETURN 
               ENDIF 
               ww = werte(1)
               DO i = 2, ianz 
                  CALL eval(cpara (i), lpara (i) ) 
                  IF (ier_num /= 0) THEN 
                     RETURN 
                  ENDIF 
                  wwerte(i) = do_read_number (cpara (i), lpara (i) ) 
                  IF (ier_num /= 0) THEN 
                     RETURN 
                  ENDIF 
                  ww = MAX(ww, wwerte(i) ) 
               ENDDO 
               CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'min') then 
            CALL get_params (line, ianz, cpara, lpara, MAXW, lp) 
            IF (ianz>=1 .AND. ianz <=MAXW) THEN 
               CALL eval(cpara(1), lpara(1) ) 
               IF(ier_num /= 0) THEN 
                  RETURN 
               ENDIF 
               werte(1) = do_read_number(cpara(1), lpara(1) ) 
               IF (ier_num /= 0) THEN 
                  RETURN 
               ENDIF 
               ww = werte(1)
               DO i = 2, ianz 
                  CALL eval(cpara (i), lpara (i) ) 
                  IF (ier_num /= 0) THEN 
                     RETURN 
                  ENDIF 
                  wwerte(i) = do_read_number (cpara (i), lpara (i) ) 
                  IF (ier_num /= 0) THEN 
                     RETURN 
                  ENDIF 
                  ww = MIN(ww, wwerte(i) ) 
               ENDDO 
               CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'mod') then 
            CALL get_params (line, ianz, cpara, lpara, 2, lp) 
            IF (ianz.eq.2) then 
               DO i = 1, ianz 
               CALL eval (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               wwerte(i) = do_read_number (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               ENDDO 
               ww = MOD(wwerte(1), wwerte(2) ) 
               CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'ran') then 
            ww = ran1 (idum) 
            CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
         ELSE 
            CALL p_calc_intr_spec (string, line, ikl, iklz, ww, lll, lp) 
         ENDIF 
      ELSEIF (lcom.eq.2) then 
         IF (string (ikl - 2:ikl - 1) .eq.'ln') then 
            IF (ww.gt.0.0) then 
               ww = log (ww) 
               CALL ersetz2 (string, ikl, iklz, ww, 2, lll) 
            ELSE 
               ier_num = - 20 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            CALL p_calc_intr_spec (string, line, ikl, iklz, ww, lll, lp) 
         ENDIF 
      ELSE 
         CALL p_calc_intr_spec (string, line, ikl, iklz, ww, lll, lp) 
      ENDIF 
!                                                                       
      END SUBROUTINE calc_intr                      
!
!*****7************************************************************************* 
!
SUBROUTINE intr_psvgt(line, lp, ww)
!-
! Calculate a pseudovoigt function
!+
USE ber_params_mod
USE element_data_mod
USE errlist_mod
USE get_params_mod
USE lib_f90_profile
USE precision_mod
USE take_param_mod
USE trig_degree_mod
USE wink_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: lp
REAL(KIND=PREC_DP), INTENT(OUT)   :: ww
!
INTEGER, PARAMETER :: MAXW = 30
CHARACTER(LEN=MAX(PREC_STRING,LEN(LINE))), DIMENSION(MAXW) :: cpara
INTEGER                                  , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP)                       , DIMENSION(MAXW) :: werte
INTEGER :: ianz          ! number of parameters
REAL(KIND=PREC_DP)  :: xx
REAL(KIND=PREC_DP)  :: P_eta
REAL(KIND=PREC_DP)  :: P_inte
REAL(KIND=PREC_DP)  :: P_pos
REAL(KIND=PREC_DP)  :: P_fwhm
REAL(KIND=PREC_DP)  :: P_asym1
REAL(KIND=PREC_DP)  :: P_asym2
REAL(KIND=PREC_DP)  :: itwo           ! Ratio Intensity Ka2/Ka1
!REAL(KIND=PREC_DP)  :: xw             ! Deviation from position 
!REAL(KIND=PREC_DP)  :: zz             ! Value for asymmetry function
!REAL(KIND=PREC_DP)  :: fa             ! Value for asymmetry function
!REAL(KIND=PREC_DP)  :: fb             ! Value for asymmetry function
!REAL(KIND=PREC_DP)  :: asym           ! Value for asymmetry function
LOGICAL             :: axis           ! TRUE==TTH; FALSE=Q
INTEGER             :: nwave          ! Entry in per_wavel
CHARACTER  (LEN=4) :: symbol
REAL               :: lambda1
REAL               :: lambda2
REAL               :: lam_1_2         ! Ratio lam(Ka1)/lam(Ka2)
!
INTEGER, PARAMETER :: NOPTIONAL = 3
INTEGER, PARAMETER :: O_ITWO    = 1
INTEGER, PARAMETER :: O_WAVE    = 2
INTEGER, PARAMETER :: O_AXIS    = 3
CHARACTER(LEN=   7)       , DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 1 ! Number of values to calculate
!
!
DATA oname  / 'itwo  ' , 'wave  '  ,  'axis   ' /
DATA loname /  4       ,  4        ,   4      /
!
opara  =  (/ '0.500000', 'CU      ',  'tth     '/)   ! Always provide fresh default values
lopara =  (/  8        ,  2        ,   3        /)
owerte =  (/  0.500000 ,  1.541800 ,   0.000000 /)
!
CALL get_params (line, ianz, cpara, lpara, MAXW, lp) 
IF(ier_num/=0) RETURN
!
IF(ianz>=5) THEN
   CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                     oname, loname, opara, lopara, lpresent, owerte)
   IF(ier_num/=0) RETURN
   CALL ber_params(ianz, cpara, lpara, werte, MAXW)
   IF(ier_num/=0) RETURN
!
   xx      = werte(1)
   P_eta   = werte(2)
   P_inte  = werte(3)
   P_pos   = werte(4)
   P_fwhm  = werte(5)
   P_asym1 = 0.0D0
   P_asym2 = 0.0D0
!
   IF(ianz>5) P_asym1 = werte(6)
   IF(ianz>6) P_asym2 = werte(7)
   IF(opara(O_AXIS)=='tth') THEN
      axis = .TRUE.
   ELSEIF(opara(O_AXIS)=='Q') THEN
      axis = .FALSE.
   ELSE
      ier_num = -6
      ier_typ = ER_FORT
      RETURN
   ENDIF
   nwave = get_wave_number(opara(O_WAVE))
   IF(nwave>0) THEN                   ! Found wave length entry
      CALL get_sym_length(nwave, symbol, lambda1)
   ELSE
      cpara(1) = opara(O_WAVE)
      lpara(1) = lopara(O_WAVE)
      ianz = 1
      CALL ber_params(ianz, cpara, lpara, werte, MAXW)
      IF(ier_num/=0) RETURN
      lambda1 = werte(1)
   ENDIF
!
   ww = pseudo_voigt(xx, P_eta, P_inte, P_pos, P_fwhm, P_asym1, P_asym2, axis)
   IF(symbol(3:4)=='12') THEN            ! Kalpha 1,2 doublet
      IF(lpresent(O_ITWO)) THEN          ! User supplied intensity ratios
         itwo = owerte(O_ITWO)
      ELSE
         itwo = get_ka21_inte(nwave)  ! Look up default intensity ratio
      ENDIF
      P_inte = P_inte * itwo
      lam_1_2 = get_ka12_len(nwave)
      lambda2 = lambda1/lam_1_2
      IF(axis) THEN                      ! 2Theta axis
         P_pos = 2.0D0*ASIND(lambda2/2.0D0/(lambda1/2.0D0/SIND(0.5D0*P_pos)))
      ELSE
         P_pos = P_pos/lam_1_2
      ENDIF
!      P_inte = P_inte * get_ka21_inte(nwave)
      ww = ww + pseudo_voigt(xx, P_eta, P_inte, P_pos, P_fwhm, P_asym1, P_asym2, axis)
   ENDIF
ELSE
   ier_num = -6
   ier_typ = ER_FORT
ENDIF
!
END SUBROUTINE intr_psvgt
!
!*****7************************************************************************* 
!
!END MODULE calc_intr_mod
