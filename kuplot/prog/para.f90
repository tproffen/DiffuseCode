module kuplot_para_mod
!
contains
!*****7**************************************************************   
!       Here are all the parameter setting routines ..                  
!*****7*****************************************************************
SUBROUTINE para_seti (zeile, lp, iarray, nia, nie, bef, imi, ima, lnull)                                                            
!+                                                                      
!     Sets option value in integer arrays(nia:nie)                      
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE kuplot_config 
USE kuplot_mod 
USE kuplot_words_mod
USE lib_help
USE lib_length
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp
integer         , intent(in)    :: nia
integer         , intent(in)    :: nie
INTEGER         , intent(out)   ::iarray (maxwin, maxframe, nia:nie) 
CHARACTER(len=*), intent(in   ) :: bef 
integer         , intent(in)    :: imi
integer         , intent(in)    :: ima
LOGICAL         , intent(in)    :: lnull 
!
INTEGER ,  PARAMETER :: maxw = 3 
!                                                                       
CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
CHARACTER(len=20) :: cdummy 
REAL(KIND=PREC_DP) :: werte (maxw) 
INTEGER :: ianz, ik, iw 
INTEGER :: lpara (maxw)
!                                                                       
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
IF (ianz.eq.0) then 
   cdummy = 'kuplot '//bef 
   CALL do_hel (cdummy, - len_str (cdummy) ) 
ELSEIF (ianz.eq.2) then 
   CALL get_words  (ianz, cpara, lpara, maxw, 2, bef) 
   CALL ber_params (ianz, cpara, lpara, werte, maxw) 
   IF (ier_num.ne.0) return 
   ik = nint (werte (1) ) 
   iw = nint (werte (2) ) 
   IF ( (ik.le. (iz - 1) .and.ik.ge.1) .or. (ik.eq.0.and.lnull) ) then
      IF (iw.ge.imi.and.iw.le.ima) then 
         iarray (iwin, iframe, ik) = iw 
      ELSE 
         ier_num = - 7 
         ier_typ = ER_APPL 
      ENDIF 
   ELSE 
      ier_num = - 4 
      ier_typ = ER_APPL 
   ENDIF 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
END SUBROUTINE para_seti                      
!
!*****7*****************************************************************
      SUBROUTINE para_setii (zeile, lp, iarray, nia, nib, befehl, imi,  &
      ima)                                                              
!+                                                                      
!     Sets option value in integer arrays(nia,nib)                      
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
      USE kuplot_words_mod
USE lib_help
USE lib_length
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp 
integer         , intent(in)    :: nia
integer         , intent(in)    :: nib
INTEGER         , intent(out)   :: iarray (maxwin, maxframe, nia, nib)
CHARACTER(len=*), intent(in)    :: befehl
integer         , intent(in)    :: imi
integer         , intent(in)    :: ima
      INTEGER maxw 
      PARAMETER (maxw = 4) 
!                                                                       
!     CHARACTER ( * ) zeile, befehl 
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      CHARACTER(20) cdummy 
      INTEGER lpara (maxw)!, nie
      INTEGER ianz, ik, iw, ip 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         cdummy = 'kuplot '//befehl 
         CALL do_hel (cdummy, - len_str (cdummy) ) 
      ELSEIF (ianz.eq.3) then 
         CALL get_words  (ianz, cpara, lpara, maxw, 2, befehl) 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ik = nint (werte (1) ) 
         ip = nint (werte (2) ) 
         iw = nint (werte (3) ) 
         IF (ik.le. (iz - 1) .and.ik.ge.1) then 
            IF (ip.ge.1.and.ip.le.maxhl) then 
               IF (iw.ge.imi.and.iw.le.ima) then 
                  iarray (iwin, iframe, ik, ip) = iw 
               ELSE 
                  ier_num = - 7 
                  ier_typ = ER_APPL 
               ENDIF 
            ELSE 
               ier_num = - 14 
               ier_typ = ER_APPL 
            ENDIF 
         ELSE 
            ier_num = - 4 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE para_setii                     
!*****7*****************************************************************
      SUBROUTINE para_setr (zeile, lp, array, befehl, awert) 
!+                                                                      
!     Sets option value in integer arrays                               
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp 
REAL            , intent(out)   :: array
CHARACTER(len=*), intent(in)    :: befehl
REAL            , intent(in)    :: awert 
      INTEGER maxw 
      PARAMETER (maxw = 2) 
!                                                                       
!     CHARACTER ( * ) zeile, befehl 
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw)
      INTEGER ianz 
      REAL(KIND=PREC_DP) :: werte (maxw)
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         array = awert 
      ELSEIF (ianz.eq.1) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         array = werte (1) 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE para_setr                      
!*****7**************************************************************   
      SUBROUTINE para_set_title (zeile, lp, string) 
!+                                                                      
!     Sets strings for title                                            
!-                                                                      
      USE build_name_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp 
CHARACTER(len=*), intent(out)   :: string
!
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw), ianz 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      lp = -lp
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         string = ' ' 
      ELSE 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.ne.0) return 
         string = cpara (1) 
      ENDIF 
!                                                                       
      END SUBROUTINE para_set_title                 
!*****7**************************************************************   
      SUBROUTINE para_set_achse (zeile, lp, string, lflag) 
!+                                                                      
!     Sets strings for axes                                             
!-                                                                      
      USE build_name_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
USE str_comp_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp 
CHARACTER(len=*), intent(out)   :: string
LOGICAL         , intent(out)   :: lflag
!
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw), ianz 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      lp = -lp
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         string = ' ' 
      ELSE 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.ne.0) return 
         string = cpara (1) 
!                                                                       
         IF (ianz.eq.2) then 
            lflag = str_comp (cpara (2) , 'log', 2, lpara (2) , 3) 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE para_set_achse                 
!*****7*****************************************************************
      SUBROUTINE set_fill (zeile, lp) 
!+                                                                      
!     Sets options for filling                                          
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
      USE kuplot_words_mod
USE lib_help
USE lib_length
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp 
!
      INTEGER maxw 
      PARAMETER (maxw = 7) 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw), cdummy 
      INTEGER lpara (maxw)
      INTEGER ianz, ik, ic, it 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         cdummy = 'kuplot fill' 
         CALL do_hel (cdummy, - len_str (cdummy) ) 
!                                                                       
      ELSEIF (ianz.eq.3.or.ianz.eq.5.or.ianz.eq.7) then 
         CALL get_words  (ianz, cpara, lpara, maxw, 2, 'fcol') 
         CALL get_words  (ianz, cpara, lpara, maxw, 3, 'fill') 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         ik = nint (werte (1) ) 
         ic = nint (werte (2) ) 
         it = nint (werte (3) ) 
!                                                                       
         IF (ik.gt. (iz - 1) .or.ik.lt.1) then 
            ier_num = - 4 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
!                                                                       
         IF (ic.gt.15.or.ic.lt.1.or.it.gt.8.or.it.lt.1) then 
            ier_num = - 7 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
!                                                                       
         ifillcol (iwin, iframe, ik) = ic 
         ifilltyp (iwin, iframe, ik) = it 
         fillrange (iwin, iframe, ik, 1) = - 9999. 
         fillrange (iwin, iframe, ik, 2) = - 9999. 
         fillrange (iwin, iframe, ik, 3) = - 9999. 
         fillrange (iwin, iframe, ik, 4) = - 9999. 
!                                                                       
         IF (ianz.eq.5) then 
            fillrange (iwin, iframe, ik, 1) = werte (4) 
            fillrange (iwin, iframe, ik, 2) = werte (5) 
         ENDIF 
         IF (ianz.eq.7) then 
            fillrange (iwin, iframe, ik, 1) = werte (4) 
            fillrange (iwin, iframe, ik, 2) = werte (5) 
            fillrange (iwin, iframe, ik, 3) = werte (6) 
            fillrange (iwin, iframe, ik, 4) = werte (7) 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE set_fill                       
!*****7*****************************************************************
      SUBROUTINE set_bond (zeile, lp) 
!+                                                                      
!     Sets options for bond drawings                                    
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp 
!
      INTEGER maxw 
      PARAMETER (maxw = 6) 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw)
      INTEGER ib, ianz 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         WRITE (output_io, 1000) iframe 
         DO ib = 1, maxbond 
         IF (bond_rad (iwin, iframe, ib) .gt.0.0) then 
            WRITE (output_io, 1100) ib, bond_rad (iwin, iframe, ib),    &
            bond_sig (iwin, iframe, ib), bond_lcol (iwin, iframe, ib),  &
            bond_ltyp (iwin, iframe, ib), bond_lwid (iwin, iframe, ib)  
         ENDIF 
         ENDDO 
      ELSEIF (ianz.ge.2) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         ib = nint (werte (1) ) 
         IF (ib.ge.1.and.ib.le.maxbond) then 
            IF (ianz.eq.2) then 
               bond_rad (iwin, iframe, ib) = werte (2) 
            ELSEIF (ianz.eq.3) then 
               bond_rad (iwin, iframe, ib) = werte (2) 
               bond_sig (iwin, iframe, ib) = werte (3) 
            ELSEIF (ianz.eq.6) then 
               bond_rad (iwin, iframe, ib) = werte (2) 
               bond_sig (iwin, iframe, ib) = werte (3) 
               bond_ltyp (iwin, iframe, ib) = nint (werte (4) ) 
               bond_lcol (iwin, iframe, ib) = nint (werte (5) ) 
               bond_lwid (iwin, iframe, ib) = werte (6) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSE 
            ier_num = - 33 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1000 FORMAT     (' Bond drawing settings for frame ',I3,' ...',//,     &
     &                   '    #   distance  sigma  col typ width ',/,   &
     &                   '   ----------------------------------- ')     
 1100 FORMAT     (3X,I3,2X,F8.4,2X,F5.3,3X,I1,2X,I1,3X,F4.2) 
      END SUBROUTINE set_bond                       
!*****7*****************************************************************
      SUBROUTINE set_linewidth (zeile, lp) 
!+                                                                      
!     Sets linewidth for given data set                                 
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
USE lib_help
USE lib_length
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp 
!
      INTEGER maxw 
      PARAMETER (maxw = 3) 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      CHARACTER(20) cdummy 
      INTEGER lpara (maxw)
      INTEGER ianz, ik 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         cdummy = 'kuplot lwid' 
         CALL do_hel (cdummy, - len_str (cdummy) ) 
      ELSEIF (ianz.eq.2) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ik = nint (werte (1) ) 
         IF (ik.le. (iz - 1) .and.ik.ge.0) then 
            linewid (iwin, iframe, ik) = werte (2) 
         ELSE 
            ier_num = - 4 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE set_linewidth                  
!*****7*****************************************************************
      SUBROUTINE set_ibox (zeile, lp) 
!+                                                                      
!     Sets type of box & axis for active frame                          
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
USE lib_help
USE lib_length
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp 
!
      INTEGER maxw 
      PARAMETER (maxw = 2) 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      CHARACTER(20) cdummy 
      INTEGER lpara (maxw)
      INTEGER ianz, it 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         cdummy = 'kuplot fset' 
         CALL do_hel (cdummy, - len_str (cdummy) ) 
      ELSEIF (ianz.eq.1) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         it = nint (werte (1) ) 
         IF (abs (it) .ge.0.and.abs (it) .le.3) then 
            ibox (iwin, iframe) = it 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE set_ibox                       
!*****7*****************************************************************
      SUBROUTINE set_sizemark (zeile, lp) 
!+                                                                      
!     Sets size of markers                                              
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
USE str_comp_mod
      USE string_convert_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp 
!
      INTEGER maxw 
      PARAMETER (maxw = 3) 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw)
      INTEGER ianz, ik 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz == 3) THEN
         CALL do_cap(cpara(3))
         IF( str_comp(cpara(3), 'X', 1, lpara(3), 1)) THEN
            cpara(3) = '-1'
            lpara(3) = 2
         ELSEIF( str_comp(cpara(3), 'Y', 1, lpara(3), 1)) THEN
            cpara(3) = '-2'
            lpara(3) = 2
         ENDIF
      ENDIF
      IF (ianz.eq.2 .OR. ianz==3) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ik = nint (werte (1) ) 
         IF (ik.le. (iz - 1) .and.ik.ge.0) then 
            sizemark (iwin, iframe, ik) = werte (2) 
         ELSE 
            ier_num = - 4 
            ier_typ = ER_APPL 
         ENDIF 
         IF(ianz==3) THEN
            IF(NINT(werte(3)) < 0) THEN
               rel_mark(iwin, iframe, ik) = NINT(werte (3))
            ELSEIF(NINT(werte(3)) == 0) THEN
               rel_mark(iwin, iframe, ik) = 0
            ELSE
               IF(lenc(ik) == lenc(NINT(werte(3)))) THEN
                  rel_mark(iwin, iframe, ik) = NINT(werte (3))
               ELSE
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  ier_msg(1) = 'The reference data set for markers '
                  ier_msg(2) = ' and the current data set differ in length'
                  rel_mark(iwin, iframe, ik) = 0
               ENDIF
            ENDIF
         ELSE
            rel_mark(iwin, iframe, ik) = 0
         ENDIF
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE set_sizemark                   
!*****7*****************************************************************
      SUBROUTINE set_legend (zeile, lp) 
!+                                                                      
!     Sets possible caption for given data set                          
!-                                                                      
      USE ber_params_mod
      USE build_name_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
USE lib_help
USE lib_length
USE precision_mod
USE str_comp_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp 
!
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      CHARACTER(20) cdummy 
      INTEGER lpara (maxw)
      INTEGER ianz, ik 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         cdummy = 'kuplot sleg' 
         CALL do_hel (cdummy, - len_str (cdummy) ) 
      ELSEIF (ianz.ge.2) then 
         CALL ber_params (1, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ik = nint (werte (1) ) 
         IF (ik.gt.0.and.ik.le. (iz - 1) ) then 
            IF (str_comp (cpara (2) , 'off', 3, lpara (2) , 3) ) then 
               ilegend (iwin, iframe, ik) = 0 
               info_orig (iwin, iframe, 1) = - 9999. 
            ELSE 
               lp = -lp
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp)
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
               IF (ier_num.ne.0) return 
               ilegend (iwin, iframe, ik) = 1 
               clegend (iwin, iframe, ik) = cpara (1)(1:MIN(40,LEN_TRIM(cpara(1))))
               IF (ianz.eq.3) then 
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  info_orig (iwin, iframe, 1) = werte (1) 
                  info_orig (iwin, iframe, 2) = werte (2) 
               ELSE 
                  info_orig (iwin, iframe, 1) = - 9999. 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 4 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE set_legend                     
!*****7*****************************************************************
      SUBROUTINE set_annotation (zeile, lp) 
!+                                                                      
!     Sets possible annotations                                         
!-                                                                      
      USE ber_params_mod
      USE build_name_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_show_mod
USE precision_mod
USE str_comp_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp 
!
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw)
      INTEGER ianz, ik 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         CALL show_annotation 
      ELSEIF (ianz.ge.2) then 
         CALL ber_params (1, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ik = nint (werte (1) ) 
         IF (ik.gt.0.and.ik.le.maxan) then 
            IF (str_comp (cpara (2) , 'OFF', 3, lpara (2) , 3) ) then 
               antext (iwin, iframe, ik) = 'OFF' 
            ELSE 
               lp = -lp
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp)
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
               IF (ier_num.ne.0) return 
               antext (iwin, iframe, ik) = cpara (1) (1:MIN(40,LEN_TRIM(cpara(1))))
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ianz.eq.2.or.ianz.eq.3.or.ianz.eq.4.or.ianz.eq.6)    &
               then                                                     
                  CALL ber_params (2, cpara, lpara, werte, maxw) 
                  antx (iwin, iframe, ik) = werte (1) 
                  anty (iwin, iframe, ik) = werte (2) 
!                                                                       
                  IF (ianz.ge.3) then 
                     IF (str_comp (cpara (3) , 'left', 1, lpara (3) , 4)&
                     ) then                                             
                        anjust (iwin, iframe, ik) = if_left 
                     ELSEIF (str_comp (cpara (3) , 'righ', 1, lpara (3) &
                     , 4) ) then                                        
                        anjust (iwin, iframe, ik) = if_right 
                     ELSEIF (str_comp (cpara (3) , 'cent', 1, lpara (3) &
                     , 4) ) then                                        
                        anjust (iwin, iframe, ik) = if_centre 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     anjust (iwin, iframe, ik) = if_left 
                  ENDIF 
!                                                                       
                  IF (ianz.ge.4) then 
                     cpara (3) = '0.0' 
                     lpara (3) = 3 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     anangle (iwin, iframe, ik) = werte (4) 
                  ELSE 
                     anangle (iwin, iframe, ik) = 0.0 
                  ENDIF 
!                                                                       
                  IF (ianz.eq.6) then 
                     anx (iwin, iframe, ik) = werte (5) 
                     any (iwin, iframe, ik) = werte (6) 
                  ELSE 
                     anx (iwin, iframe, ik) = antx (iwin, iframe, ik) 
                     any (iwin, iframe, ik) = anty (iwin, iframe, ik) 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 28 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE set_annotation                 
!*****7*****************************************************************
      SUBROUTINE set_buff (zeile, lp) 
!+                                                                      
!     Command buff ..                                                   
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp 
!
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw)
      INTEGER ianz 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.1) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ibuf (iwin, iframe, 1) = werte (1) 
         ibuf (iwin, iframe, 2) = werte (1) 
         ibuf (iwin, iframe, 3) = werte (1) 
         ibuf (iwin, iframe, 4) = werte (1) 
      ELSEIF (ianz.eq.4) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ibuf (iwin, iframe, 1) = werte (1) 
         ibuf (iwin, iframe, 2) = werte (2) 
         ibuf (iwin, iframe, 3) = werte (3) 
         ibuf (iwin, iframe, 4) = werte (4) 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE set_buff                       
!*****7*****************************************************************
      SUBROUTINE set_skal (zeile, lp) 
!+                                                                      
!     Command skal ..                                                   
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp 
!
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw)
      INTEGER ianz 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      t (iwin, iframe, 1) = - 9999.0 
      t (iwin, iframe, 2) = - 9999.0 
!                                                                       
      IF (ianz.eq.0) then 
         ex (iwin, iframe, 1) = - 9999.0 
         ey (iwin, iframe, 1) = - 9999.0 
      ELSEIF (ianz.eq.2.or.ianz.eq.4) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         IF (werte (1) .lt.werte (2) ) then 
            ex (iwin, iframe, 1) = werte (1) 
            ex (iwin, iframe, 2) = werte (2) 
            IF (ianz.eq.4) then 
               IF (werte (3) .lt.werte (4) ) then 
                  ey (iwin, iframe, 1) = werte (3) 
                  ey (iwin, iframe, 2) = werte (4) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ELSE 
               ey (iwin, iframe, 1) = - 9999.0 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      IF (ier_num.ne.0) return 
      CALL skalieren 
!                                                                       
      END SUBROUTINE set_skal                       
!*****7*****************************************************************
      SUBROUTINE set_mark (zeile, lp) 
!+                                                                      
!     Command mark ..                                                   
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
USE lib_errlist_func
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp 
!
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw)
      INTEGER ianz 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         t (iwin, iframe, 1) = - 9999.0 
      ELSEIF (ianz.eq.1.or.ianz.eq.2.or.ianz.eq.4) then 
         IF (ianz.eq.4) write ( *, 1000) 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         CALL skalieren 
         CALL no_error 
         IF (werte (1) .gt.0.0) then 
            t (iwin, iframe, 1) = werte (1) 
            IF (ianz.eq.2) then 
               IF (werte (2) .gt.0.0) then 
                  t (iwin, iframe, 2) = werte (2) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1000 FORMAT    (1x,'****WARN**** Only first two parameters',           &
     &                     '  used (see help mark) !')                  
      END SUBROUTINE set_mark                       
!*****7*****************************************************************
      SUBROUTINE set_hpak (zeile, lp) 
!+                                                                      
!     Setting of contour line intervalls                                
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_show_mod
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp 
!
      INTEGER maxw 
      PARAMETER (maxw = 2) 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw)
      INTEGER ianz, ii 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
!------ get parameters                                                  
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ no parameters : show actual settings                            
!                                                                       
      IF (ianz.eq.0) then 
         CALL show_hlin 
!                                                                       
!------ set number of contour line sets                                 
!                                                                       
      ELSEIF (ianz.eq.1) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ii = nint (werte (1) ) 
         IF (ii.gt.0.and.ii.lt.maxhl) then 
            iho (iwin, iframe) = ii 
         ELSE 
            ier_num = - 14 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE set_hpak                       
!
!*****7*****************************************************************
!
SUBROUTINE set_hlin(zeile, lp) 
!+                                                                      
!     Setting of contour line intervalls                                
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE kuplot_config 
USE kuplot_mod 
use kuplot_extrema_mod
use kuplot_low_mod
use kuplot_show_mod
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: MAXW = 6
!                                                                       
CHARACTER(LEN=*), INTENT(INOUT) :: zeile    ! incoming command
INTEGER,          INTENT(INOUT) :: lp       ! length
!
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
INTEGER                   , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP)        , DIMENSION(MAXW) :: werte
REAL(kind=PREC_SP)    :: zzmin, zzmax, zhub     ! MIN/MAX/Delta values for date/contour lines
REAL(kind=PREC_SP)    :: zlow    ! Minimum value for contour lines
INTEGER :: ianz, ihl, i 
!LOGICAL :: k_in_f  ! Current data set is in frame
LOGICAL :: proz    ! User specified a "%" als last parameter
LOGICAL :: labs    ! Percent scale is Based on zmax only
!                                                                       
!------ get parameters                                                  
!                                                                       
CALL get_params(zeile, ianz, cpara, lpara, MAXW, lp) 
IF(ier_num.ne.0) RETURN 
!                                                                       
!------ no parameters : show actual settings                            
!                                                                       
IF (ianz.eq.0) THEN 
   CALL show_hlin 
!                                                                       
!------ set contour line base, intervall, number and optinally mode     
!                                                                       
ELSE
   IF(ianz>=3 .AND. ianz<=5) THEN
      labs = .FALSE.
      proz = .FALSE.
      IF(cpara(ianz)(1:1) == '%') THEN     ! Data in percent
         proz = .TRUE.
         IF(cpara(ianz)(2:2) == 'm') labs = .TRUE.
         ianz = ianz - 1
      ELSE
         proz = .FALSE.
      ENDIF
      IF(ianz == 4.OR.ianz == 3) THEN
!        proz = .false. 
!        IF (ianz.eq.5) then 
!           proz = (cpara (5) (1:1) .eq.'%') 
!           ianz = ianz - 1 
!        ENDIF 
!                                                                       
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF(ier_num /= 0) RETURN 
!                                                                       
         ihl = NINT(werte(1) ) 
         IF(ihl >= 1.AND.ihl <= maxhl) THEN 
            iho(iwin, iframe) = MAX(iho(iwin, iframe), ihl) 
            nz(iwin, iframe, ihl) = NINT(werte(ianz) ) 
!                                                                       
               CALL get_extrema 
               zzmax = -1.e19 
               zzmin =  1.e19 
               DO i = 1, iz - 1 
                  IF(k_in_f(i) ) THEN 
                     zzmax = max (zzmax, real(zmax (i)) ) 
                     zzmin = min (zzmin, real(zmin (i)) ) 
                  ENDIF 
               ENDDO 
!                                                                       
            IF(proz) THEN                  ! Contour lines in percent
               IF(ianz==4) THEN            ! User specified "hmin, delta, n, %"
                  IF(labs) THEN            ! User specified percentage of MAXIMUM "%m"
                     zhub = zzmax
                     zlow = 0.0
                  ELSE
                     zhub = zzmax - zzmin 
                     zlow = zzmin
                  ENDIF
                  IF (zhub.eq.0) zhub = 100.0 
                  z_min(iwin, iframe, ihl) = zlow + 0.01 * zhub * werte(2)
                  z_inc(iwin, iframe, ihl) = 0.01 * zhub * ABS(werte(3))
               ELSE                        ! User specified "delta, n, %"
                  IF(zzmin>= 0.0) THEN     ! positive values only
                     zhub = zzmax
                     IF (zhub.eq.0) zhub = 100.0 
                     z_min(iwin, iframe, ihl) = MAX(zzmin, 0.00)
                     z_inc(iwin, iframe, ihl) = 0.01 * zhub * ABS(werte(2) )
                  ELSEIF(zzmin<0.0 .AND. zzmax > 0.0) THEN    ! neg and positive
                     zhub = zzmax          ! The maximum intensity
                     IF (zhub.eq.0) zhub = 100.0 
                     z_min(iwin, iframe, ihl) = -0.50 * zhub * ABS(werte(2))
                     z_inc(iwin, iframe, ihl) =  0.01 * zhub * ABS(werte(2))
!write(*,*) ' User specified "delta, n, %" Positive and negative values ', zhub
                  ELSE                     ! Negative values only
                     zhub = -zzmin 
                     IF (zhub.eq.0) zhub = 100.0 
                     z_min(iwin, iframe, ihl) = MIN(zzmin, 0.00)
                     z_inc(iwin, iframe, ihl) = 0.01 * zhub * ABS(werte(2))
                  ENDIF
               ENDIF
            ELSE                           ! Contour lines in absolute values
               IF(ianz==4) THEN            ! User specified "hmin, delta, n"
                  z_min(iwin, iframe, ihl) = werte (2) 
                  z_inc(iwin, iframe, ihl) = werte (3) 
               ELSE                        ! User specified "delta, n"
                  IF(zzmin>= 0.0) THEN     ! positive values only
                     z_min(iwin, iframe, ihl) = MAX(zzmin, 0.00)
                     z_inc(iwin, iframe, ihl) = werte(2) 
                  ELSEIF(zzmin<0.0 .AND. zzmax > 0.0) THEN    ! neg and positive
                     z_min(iwin, iframe, ihl) = -ABS(werte(2))*nz(iwin, iframe, ihl)/2.
                     z_inc(iwin, iframe, ihl) = werte(2) 
                  ELSE                     ! Negative values only
                     z_min(iwin, iframe, ihl) = -ABS(werte(2))*nz(iwin, iframe, ihl)
                     z_inc(iwin, iframe, ihl) = werte(2) 
                  ENDIF
               ENDIF
            ENDIF 
         ELSE 
            ier_num = - 14 
            ier_typ = ER_APPL 
         ENDIF 
      ENDIF 
   ELSE 
      ier_num = - 6 
      ier_typ = ER_COMM 
   ENDIF 
ENDIF 
!write(*,*) ' HLIN ', z_min(iwin, iframe, ihl), z_inc(iwin, iframe, ihl), &
!  z_min(iwin, iframe, ihl)+z_inc(iwin, iframe, ihl)*nz(iwin, iframe, ihl), &
!  nz(iwin, iframe, ihl)
!                                                                       
END SUBROUTINE set_hlin                       
!
!*****7*****************************************************************
!
SUBROUTINE set_font (zeile, lp) 
!+                                                                      
!     Setting of font size, type and color                              
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE kuplot_config 
USE kuplot_mod 
USE kuplot_words_mod
use kuplot_show_mod
USE precision_mod
USE str_comp_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: maxw = 4 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp 
!
CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
CHARACTER(LEN=4)           :: bef
INTEGER            :: lpara (maxw)
INTEGER            :: ianz, icol, ifon, ifid 
REAL(KIND=PREC_DP) :: werte (maxw) 
REAL               :: fsiz 
!                                                                       
!------ get parameters                                                  
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
!------ no parameters : show actual settings                            
!                                                                       
IF (ianz.eq.0) then 
   CALL show_font 
   return
endif
!
!------ Actual parameters
!
IF (ianz.ge.2) then 
!                                                                       
!------ - set font size                                                 
!                                                                       
   IF (str_comp (cpara (1) , 'size', 2, lpara (1) , 4) ) then 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      IF (ier_num.ne.0) return 
      IF (ianz.eq.2) then 
         ifon = nint (werte (1) ) 
         fsiz = werte (2) 
      IF (ifon.gt.0.and.ifon.le.nfon.and.fsiz.ge.4.0.and.fsiz.le.128.0) then
            fonsize (iwin, iframe, ifon) = 10.0 * fsiz 
         ELSE 
            ier_num = - 7 
            ier_typ = ER_APPL 
         ENDIF 
      ELSEIF (ianz.eq.1) then 
         fonscal (iwin, iframe) = werte (1) 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
!------ - set font color                                                
!                                                                       
   ELSEIF(str_comp(cpara(1), 'color', 2, lpara(1), 5)) then
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      bef = 'fcol'
      CALL get_words  (ianz, cpara, lpara, maxw, 2, bef) 
      CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      IF (ier_num.ne.0) return 
      IF (ianz.eq.2) then 
         ifon = nint (werte (1) ) 
         icol = nint (werte (2) ) 
         IF(ifon.gt.0.and.ifon.le.nfon.and.icol.gt.0.and.icol.le.15) then
            foncol (iwin, iframe, ifon) = icol 
         ELSE 
            ier_num = - 7 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
!------ - set font typ                                                  
!                                                                       
   ELSEIF (str_comp (cpara (1) , 'typ', 2, lpara (1) , 3) ) then 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      IF (ier_num.ne.0) return 
      IF (ianz.eq.2) then 
         ifon = nint (werte (1) ) 
         ifid = nint (werte (2) ) 
         IF(ifon.gt.0.and.ifon.le.nfon.and.ifid.gt.0.and.ifid.le.4) then                                                     
            fon_id (iwin, iframe, ifon) = ifid 
         ELSE 
            ier_num = - 7 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
!------ - set text justification for title / textframes                 
!                                                                       
   ELSEIF (str_comp (cpara (1) , 'just', 2, lpara (1) , 4) ) then 
      IF (ier_num.ne.0) return 
      IF (ianz.eq.2) then 
         IF(str_comp (cpara (2) , 'cen', 1, lpara (2) , 3) ) then
                  frjust (iwin, iframe) = if_centre 
         ELSEIF(str_comp(cpara(2), 'left', 1, lpara(2), 4)) then
            frjust (iwin, iframe) = if_left 
         ELSE 
            ier_num = - 7 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
   ELSE 
      ier_num = - 6 
      ier_typ = ER_COMM 
   ENDIF 
!                                                                       
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
END SUBROUTINE set_font                       
!
!****7******************************************************************
!
SUBROUTINE skalieren 
!+                                                                      
!       setting up plotting window                                      
!-                                                                      
USE errlist_mod 
USE koordinate_mod
USE kuplot_config 
USE kuplot_mod 
use kuplot_extrema_mod
use kuplot_low_mod
!
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(len=8) :: ct 
INTEGER :: i 
REAL               :: dummy_x (2), dummy_y (2) 
REAL(kind=PREC_DP) :: ymi, yma 
REAL(kind=PREC_DP) :: delta
!     LOGICAL k_in_f 
!                                                                       
!------ get data set extrema                                            
!                                                                       
      CALL get_extrema 
!                                                                       
!------ set plotting window                                             
!                                                                       
      IF (ex (iwin, iframe, 1) .eq. - 9999.) then 
         ex (iwin, iframe, 1) = 1.0e19 
         ex (iwin, iframe, 2) = - 1.0e19 
         DO i = 1, iz - 1 
         IF (k_in_f (i) ) then 
            ex (iwin, iframe, 1) = min (real(ex (iwin, iframe, 1),kind=PREC_DP), xmin (i) ) 
            ex (iwin, iframe, 2) = max (real(ex (iwin, iframe, 2),kind=PREC_DP), xmax (i) ) 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      IF (ey (iwin, iframe, 1) .eq. - 9999.) then 
         ey (iwin, iframe, 1) = 1.0e19 
         ey (iwin, iframe, 2) = - 1.0e19 
         DO i = 1, iz - 1 
         IF (k_in_f (i) ) then 
            IF ( lni(i) ) THEN
               ymi = ymin(i)
               yma = ymax(i)
            ELSE
              CALL get_extrema_xy_local (i, ymi, yma) 
            ENDIF
            ey (iwin, iframe, 1) = min (real(ey (iwin, iframe, 1),kind=PREC_DP), ymi) 
            ey (iwin, iframe, 2) = max (real(ey (iwin, iframe, 2),kind=PREC_DP), yma) 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      IF (ey (iwin, iframe, 2) .eq.ey (iwin, iframe, 1) ) then 
         IF(ABS(ey(iwin, iframe, 1)) > 0.0) THEN
            delta = 0.2*ABS(ey(iwin, iframe, 1))
         ELSE
            delta = 1.0
         ENDIF
         ey (iwin, iframe, 1) = ey (iwin, iframe, 1) - delta
         ey (iwin, iframe, 2) = ey (iwin, iframe, 2) + delta
!        ey (iwin, iframe, 1) = ey (iwin, iframe, 1) - 0.2 * ey (iwin,  &
!        iframe, 1)                                                     
!        ey (iwin, iframe, 2) = ey (iwin, iframe, 2) + 0.2 * ey (iwin,  &
!        iframe, 2)                                                     
      ENDIF 
!                                                                       
      IF (ex (iwin, iframe, 2) .eq.ex (iwin, iframe, 1) ) then 
         IF(ABS(ex(iwin, iframe, 1)) > 0.0) THEN
            delta = 0.2*ABS(ex(iwin, iframe, 1))
         ELSE
            delta = 1.0
         ENDIF
         ex (iwin, iframe, 1) = ex (iwin, iframe, 1) - delta
         ex (iwin, iframe, 2) = ex (iwin, iframe, 2) + delta
!        ex (iwin, iframe, 1) = ex (iwin, iframe, 1) - 0.2 * ex (iwin,  &
!        iframe, 1)                                                     
!        ex (iwin, iframe, 2) = ex (iwin, iframe, 2) + 0.2 * ex (iwin,  &
!        iframe, 2)                                                     
      ENDIF 
!                                                                       
!------ set tick marks                                                  
!                                                                       
      IF (t (iwin, iframe, 1) .eq. - 9999.) then 
         t (iwin, iframe, 1) = (ex (iwin, iframe, 2) - ex (iwin, iframe,&
         1) ) * 0.125                                                   
         t (iwin, iframe, 2) = (ey (iwin, iframe, 2) - ey (iwin, iframe,&
         1) ) * 0.125                                                   
         WRITE (ct, '(g8.1)') t (iwin, iframe, 1) 
         READ (ct, * ) t (iwin, iframe, 1) 
         WRITE (ct, '(g8.1)') t (iwin, iframe, 2) 
         READ (ct, * ) t (iwin, iframe, 2) 
      ENDIF 
!                                                                       
!------ Now extract the values for plotting                             
!                                                                       
      dummy_x (1) = ex (iwin, iframe, 1) 
      dummy_x (2) = ex (iwin, iframe, 2) 
      dummy_y (1) = ey (iwin, iframe, 1) 
      dummy_y (2) = ey (iwin, iframe, 2) 
!                                                                       
      IF (lachse (iwin, iframe, 1) .and.dummy_x (1) .le.0.0) then 
         ier_num = - 43 
         ier_typ = ER_APPL 
         dummy_x (1) = 1.0 
      ENDIF 
!                                                                       
      IF (lachse (iwin, iframe, 2) .and.dummy_y (1) .le.0.0) then 
         ier_num = - 43 
         ier_typ = ER_APPL 
         dummy_y (1) = 1.0 
      ENDIF 
!                                                                       
      CALL koor_log (2, dummy_x, dummy_y) 
!                                                                       
      pex (iwin, iframe, 1) = dummy_x (1) 
      pex (iwin, iframe, 2) = dummy_x (2) 
      pey (iwin, iframe, 1) = dummy_y (1) 
      pey (iwin, iframe, 2) = dummy_y (2) 
!                                                                       
      dummy_x (1) = t (iwin, iframe, 1) 
      dummy_y (1) = t (iwin, iframe, 2) 
!                                                                       
      CALL koor_log (1, dummy_x, dummy_y) 
!                                                                       
      IF (lachse (iwin, iframe, 1) ) dummy_x (1) = REAL(nint (dummy_x &
      (1) ) )                                                           
      IF (lachse (iwin, iframe, 2) ) dummy_y (1) = REAL(nint (dummy_y &
      (1) ) )                                                           
      IF (dummy_x (1) .eq.0) then 
         dummy_x (1) = 1.0 
      ENDIF 
      IF (dummy_y (1) .eq.0) then 
         dummy_y (1) = 1.0 
      ENDIF 
      pt (iwin, iframe, 1) = dummy_x (1) 
      pt (iwin, iframe, 2) = dummy_y (1) 
!                                                                       
END SUBROUTINE skalieren                      
!
!****7******************************************************************
!
SUBROUTINE write_para 
!+                                                                      
!       write parameter file 'kupl.par'                                 
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE support_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER ifil 
      PARAMETER (ifil = 21) 
!                                                                       
      INTEGER ifr, ik, i 
      LOGICAL lho 
!                                                                       
      lho = .false. 
      CALL oeffne (ifil, 'kupl.par', 'unknown') 
      IF (ier_num.ne.0) return 
!                                                                       
      WRITE (ifil, 1000) 
      WRITE (ifil, 1002) yskal_u (iwin, iframe), shear (iwin, iframe) 
      DO ik = 1, iz - 1 
      IF (lni (ik) ) then 
         WRITE (ifil, 1005) fname (ik), xmin (ik), xmax (ik), ymin (ik),&
         ymax (ik), zmin (ik), zmax (ik)                                
         lho = .true. 
      ELSE 
         WRITE (ifil, 1010) fname (ik), xmin (ik), xmax (ik), ymin (ik),&
         ymax (ik)                                                      
      ENDIF 
      ENDDO 
!                                                                       
      IF (lho) then 
         DO ifr = 1, iaf (iwin) 
         WRITE (ifil, 1020) ifr 
         DO i = 1, iho (iwin, ifr) 
         WRITE (ifil, 1030) z_min (iwin, ifr, i), z_inc (iwin, ifr, i), &
         nz (iwin, ifr, i)                                              
         ENDDO 
         ENDDO 
      ENDIF 
      CLOSE (ifil) 
!                                                                       
 1000 FORMAT     ( ' Plot - parameter : '/) 
 1002 FORMAT     ( ' Ratio y/x-axis   : ',g11.4e1,/                      &
     &                    ' Angle x/y axis   : ',g11.4e1,/)              
 1005 FORMAT     ( ' File : ',a40,/                                     &
     &           '    x : ',g11.4e1,' to ',g11.4e1/                       &
     &           '    y : ',g11.4e1,' to ',g11.4e1/                       &
     &           '    z : ',g11.4e1,' to ',g11.4e1)                       
 1010 FORMAT     ( ' File : ',a40,/                                     &
     &           '    x : ',g11.4e1,' to ',g11.4e1/                       &
     &           '    y : ',g11.4e1,' to ',g11.4e1)                       
 1020 FORMAT        (/' Contours     [frame ',i2,'] : '/                &
     &           ' minimum     step        #')                          
 1030 FORMAT     ( 1x,g11.4,1x,g11.4,2x,i2) 
!                                                                       
      END SUBROUTINE write_para                     
!
!*******************************************************************************
!
subroutine get_aver_angle(idim, ik, werte)
!-
!  Calculate aver for global 3D data set 
!+
!
use kuplot_mod
!
use errlist_mod
use lib_data_struc_h5
use lib_metric_mod
use matrix_mod
use precision_mod
!
implicit none
!
integer                            , intent(in)    :: idim
integer                            , intent(in)    :: ik
real(kind=PREC_DP), dimension(idim), intent(inout) :: werte   ! IN: data number; OUT: aver
!
integer                            :: node_number  ! Global Data node number
logical                            :: is_direct    ! Direct or reciprocal space
integer           , dimension(3)   :: inc          ! Number data points along the three axes
logical                            :: ltop         ! top vector exists for 3D data
real(kind=PREC_DP), dimension(6)   :: lattice      ! Lattice parameters
real(kind=PREC_DP), dimension(3,3) :: gten         ! metric tensot     
real(kind=PREC_DP), dimension(3,3) :: rten         ! reciprocal metric tensor
real(kind=PREC_DP), dimension(3,3) :: vi           ! Vectors along the three axes
real(kind=PREC_DP), DIMENSION(3)   ::  length
real(kind=PREC_DP)            ::  angle_vh
real(kind=PREC_DP)            ::  ratio_vh
real(kind=PREC_DP)            ::   aver_vh
real(kind=PREC_DP)            ::  angle_ht
real(kind=PREC_DP)            ::  ratio_ht
real(kind=PREC_DP)            ::   aver_ht
real(kind=PREC_DP)            ::  angle_tv
real(kind=PREC_DP)            ::  ratio_tv
real(kind=PREC_DP)            ::   aver_tv
!
!
ltop = .false.
if(ik>0 .and. ik<iz) then
   call dgl5_set_pointer(ik , ier_num, ier_typ, node_number)
   if(ier_num /= 0) return
   if(dgl5_get_ndims()==3) ltop = .true.
   is_direct = dgl5_get_direct()
   call dgl5_get_lattice(ik, lattice)
   call lib_tensor(gten, lattice(1:3), lattice(4:6))
   call matinv(gten, rten)
   call dgl5_get_dims(ik, inc)
   call dgl5_get_steps(ik, vi)
   if(is_direct) then
      call lib_angles(ltop, length, &
               angle_vh, ratio_vh, aver_vh, &
               angle_ht, ratio_ht, aver_ht, &
               angle_tv, ratio_tv, aver_tv, &
               gten, inc, vi)
   else
      call lib_angles(ltop, length, &
               angle_vh, ratio_vh, aver_vh, &
               angle_ht, ratio_ht, aver_ht, &
               angle_tv, ratio_tv, aver_tv, &
               rten, inc, vi)
   endif
   werte(1) = aver_vh
   werte(2) = angle_vh
else
   ier_num = -4
   ier_typ = ER_APPL
endif
!
end subroutine get_aver_angle
!
!*******************************************************************************
!
end module kuplot_para_mod
