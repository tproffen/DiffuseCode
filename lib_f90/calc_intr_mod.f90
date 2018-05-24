!MODULE calc_intr_mod
!
!CONTAINS
!
!*****7**************************************************************** 
!
      SUBROUTINE calc_intr (string, line, ikl, iklz, lll, lp) 
!                                                                       
!     Evaluate all intrinsic functions                                  
!                                                                       
      USE build_name_mod
      USE ersetz_mod
      USE do_read_number_mod
      USE errlist_mod 
      USE get_params_mod
      USE random_mod
      USE variable_mod
      USE wink_mod
      USE times_mod
      USE trig_degree_mod
      USE set_sub_generic_mod
      IMPLICIT none 
!                                                                       
      INTEGER, PARAMETER :: maxw = 9
!                                                                       
      CHARACTER(LEN=*), INTENT(INOUT) :: string
      CHARACTER(LEN=*), INTENT(INOUT) :: line 
      INTEGER         , INTENT(INOUT) :: ikl
      INTEGER         , INTENT(INOUT) :: iklz
      INTEGER         , INTENT(INOUT) :: lll
      INTEGER         , INTENT(INOUT) :: lp 
!
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(1024) zeile 
      CHARACTER(1024) answer , search
      INTEGER lpara (maxw) 
      INTEGER ikom, i, ianz 
      INTEGER lcom 
      INTEGER ihyp 
      INTEGER dummy 
      LOGICAL :: BACK   ! FLAG for index intrinsic
      REAL fl1, fl2, fl3, gbox_k, gbox_w, gbox_x 
      REAL werte (maxw), ww, a , skew
!     REAL sind, cosd, tand, asind, acosd, atand 
!     REAL atan2, atan2d 
!                                                                       
      INTEGER length_com 
      INTEGER len_str 
      REAL gasdev, ran1, poidev , gasskew
!
!                                                                       
      ier_num = - 1 
      ier_typ = ER_FORT 
      ikom = INDEX (line, ',') 
      ihyp = max (INDEX (line, '''') , INDEX (line, '"') ) 
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
               ww = float (i-2) 
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
            IF (abs (ww) .le.1.0) then 
               ww = asind (ww) 
               CALL ersetz2 (string, ikl, iklz, ww, 5, lll) 
            ELSE 
               ier_num = - 7 
               ier_typ = ER_FORT 
            ENDIF 
         ELSEIF (string (ikl - 5:ikl - 1) .eq.'acosd') then 
            IF (abs (ww) .le.1.0) then 
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
               werte (1) = do_read_number (cpara (1), lpara (1) ) 
               werte (2) = do_read_number (cpara (2), lpara (2) ) 
               ww = atan2d (werte (1), werte (2) ) 
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
               werte (1) = do_read_number (cpara (1), lpara (1) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (werte (1) .lt.0.0) then 
                  ier_num = - 35 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
!
               CALL eval (cpara (2), lpara (2) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               werte (2) = do_read_number (cpara (2), lpara (2) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (ABS(werte (2)) .gt.1.0) then 
                  ier_num = - 36 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
               skew = werte(2)
!                                                                       
               IF (ianz.eq.2.or.cpara (3) .eq.'s') then 
                  a = werte (1) 
               ELSEIF (ianz.eq.3.and.cpara (3) .eq.'f') then 
                  a = werte (1) / sqrt (8. * log (2.) ) 
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
            CALL ersetzc (string, ikl, iklz, f_modt, 24, 5, lll) 
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
                                 ww = FLOAT(INDEX(zeile (1:LEN_TRIM(ZEILE )), &
                                                  search(1:LEN_TRIM(search)), &
                                                  BACK ))
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
         ELSE 
            CALL p_calc_intr_spec (string, line, ikl, iklz, ww, lll, lp) 
         ENDIF 
      ELSEIF (lcom.eq.4) then 
         IF (string (ikl - 4:ikl - 1) .eq.'asin') then 
            IF (abs (ww) .le.1.0) then 
               ww = asin (ww) 
               CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
            ELSE 
               ier_num = - 7 
               ier_typ = ER_FORT 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'acos') then 
            IF (abs (ww) .le.1.0) then 
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
               werte (1) = do_read_number (cpara (1), lpara (1) ) 
               werte (2) = do_read_number (cpara (2), lpara (2) ) 
               ww = atan2 (werte (1), werte (2) ) 
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
            IF (ww.ge.0) then 
               ww = sqrt (ww) 
               CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
            ELSE 
               ier_num = - 5 
               ier_typ = ER_FORT 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'nint') then 
            ww = float (nint (ww) ) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'frac') then 
            ww = ww - float (int (ww) ) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'gbox') then 
            CALL get_params (line, ianz, cpara, lpara, 3, lp) 
            IF (ianz.eq.3) then 
               DO i = 1, 3 
               CALL eval (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               werte (i) = do_read_number (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (werte (i) .lt.0.0) then 
                  ier_num = - 35 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
               ENDDO 
!           Determine relative size of exponentials and box             
               fl1 = REAL(0.5D0 * REAL(werte (1), KIND=KIND(0.0D0)) * sqrt (zpi), KIND=KIND(0.0E0)) 
               fl2 = werte (2) 
               fl3 = REAL(0.5D0 * REAL(werte (3), KIND=KIND(0.0D0)) * sqrt (zpi), KIND=KIND(0.0E0)) 
!           Normalize to 1                                              
               gbox_k = 1. / (fl1 + fl2 + fl3) 
               gbox_w = 1. - (fl1 + fl3) * gbox_k 
!           Get random number                                           
               gbox_x = ran1 (idum) 
!                                                                       
               IF (gbox_x.lt.fl1 * gbox_k) then 
                  ww = - werte (2) * 0.5 - abs (gasdev (werte (1) ) ) 
               ELSEIF (gbox_x.lt. (fl1 + fl2) * gbox_k) then 
                  ww = - werte (2) * 0.5 + (gbox_x - fl1 * gbox_k)      &
                       * werte (2) / gbox_w
               ELSE 
                  ww = werte (2) * 0.5 + abs (gasdev (werte (3) ) ) 
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
               werte (1) = do_read_number (cpara (1), lpara (1) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (werte (1) .lt.0.0) then 
                  ier_num = - 35 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
!                                                                       
               IF (ianz.eq.1.or.cpara (2) .eq.'s') then 
                  a = werte (1) 
               ELSEIF (ianz.eq.2.and.cpara (2) .eq.'f') then 
                  a = werte (1) / sqrt (8. * log (2.) ) 
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
               werte (1) = do_read_number (cpara (1), lpara (1) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (werte (1) .lt.0.0) then 
                  ier_num = - 37 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
!                                                                       
               CALL eval (cpara (2), lpara (2) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               werte (2) = do_read_number (cpara (2), lpara (2) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (werte (2) .lt.0.0) then 
                  ier_num = - 35
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
!                                                                       
               IF (ianz.eq.2.or.cpara (3) .eq.'s') then 
                  a = werte (2) 
               ELSEIF (ianz.eq.3.and.cpara (3) .eq.'f') then 
                  a = werte (2) / sqrt (8. * log (2.) ) 
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
            ww = exp (log (werte (1) ) + gasdev (a) ) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'pois') then 
            ww = poidev (ww, idum) 
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
            ww = float (int (ww) ) 
            CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'max') then 
            CALL get_params (line, ianz, cpara, lpara, 2, lp) 
            IF (ianz.eq.2) then 
               DO i = 1, ianz 
               CALL eval (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               werte (i) = do_read_number (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               ENDDO 
               ww = max (werte (1), werte (2) ) 
               CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'min') then 
            CALL get_params (line, ianz, cpara, lpara, 2, lp) 
            IF (ianz.eq.2) then 
               DO i = 1, ianz 
               CALL eval (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               werte (i) = do_read_number (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               ENDDO 
               ww = min (werte (1), werte (2) ) 
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
               werte (i) = do_read_number (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               ENDDO 
               ww = amod (werte (1), werte (2) ) 
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
!*****7**************************************************************** 
!END MODULE calc_intr_mod
