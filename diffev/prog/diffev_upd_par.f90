module diffev_update_mod
!
contains
!
SUBROUTINE diffev_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz) 
!                                                                       
!-                                                                      
!       replaces a substring in an expression by the value of the       
!       appropriate Parameter.                                          
!     This is needed if the parameter is read.                          
!+                                                                      
USE diff_evol
USE population
USE blanks_mod
USE errlist_mod 
USE lib_errlist_func
USE lib_upd_mod
USE lib_length
USE param_mod 
USE precision_mod
USE precision_command_mod
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
REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(IN   ) :: ww
!
CHARACTER (LEN=MAX(PREC_STRING,LEN(string))) :: zeile 
!                                                                       
INTEGER                              :: laenge, ltyp, kpara, kpara2
INTEGER                              :: lcomm 
!
CALL lib_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
IF(ier_num == 0) RETURN
CALL no_error
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
   IF(ikl.gt.lcomm + 1) zeile(1:ikl - lcomm - 1) = string(1: ikl - lcomm - 1)
!
   IF (string (ikl - 1:ikl - 1) .eq.'p') then 
      IF (ianz.eq.1) then 
         IF (0.le.kpara.and.kpara.le.MAXDIMX .and. kpara<=pop_dimx) THEN 
            WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_REAL) pop_para (kpara)
            zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
ELSEIF (lcomm.eq.5) then 
   IF (string (ikl - 5:ikl - 1) .eq.'pop_n') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.6) zeile (1:ikl - 6) = string (1:ikl - 6) 
         WRITE (zeile (ikl - 5:ikl + PREC_WIDTH-2) , PREC_F_INTE) pop_n
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 5:ikl - 1) .eq.'pop_c') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.6) zeile (1:ikl - 6) = string (1:ikl - 6) 
         WRITE (zeile (ikl - 5:ikl + PREC_WIDTH-2) , PREC_F_INTE) pop_c
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 5:ikl - 1) .eq.'bestm') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.6) zeile (1:ikl - 6) = string (1:ikl - 6) 
         WRITE (zeile (ikl - 5:ikl + PREC_WIDTH-2) , PREC_F_INTE) pop_best 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 5:ikl - 1) .eq.'bestr') then 
      IF (ianz.eq.1) then 
         IF (0<pop_best .and. pop_best<=pop_n .and. pop_best<= MAXPOP) THEN
            IF (ikl.gt.6) zeile (1:ikl - 6) = string (1:ikl - 6) 
            WRITE (zeile (ikl - 5:ikl + PREC_WIDTH-2) , PREC_F_REAL) child_val (pop_best,0)
            zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
         ELSE
            ier_num = -9
            ier_typ = ER_APPL
            ier_msg(1) = 'The population has not been defined yet' 
         ENDIF
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 5:ikl - 1) .eq.'pop_v') then 
      IF (ianz.eq.2) then 
         IF( 0<kpara .and. kpara<=pop_dimx .and. kpara<=MAXDIMX) THEN
            IF( 0<kpara2 .and. kpara2<=pop_n .and. kpara2<=MAXPOP) THEN
               IF (ikl.gt.6) zeile (1:ikl - 6) = string (1:ikl - 6) 
               WRITE (zeile (ikl - 5:ikl + PREC_WIDTH-2) , PREC_F_REAL) child ( kpara, kpara2)
               zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
            ELSE
               ier_num = -14
               ier_typ = ER_APPL
               ier_msg(1) = 'Possible reasons:'
               ier_msg(2) = 'The population has not been defined' 
               ier_msg(3) = 'The second index is <0 or > pop_n'
            ENDIF
         ELSE
            ier_num = -9
            ier_typ = ER_APPL
            ier_msg(1) = 'Possible reasons:'
            ier_msg(2) = 'The value of pop_dimx[1] is not defined' 
            ier_msg(3) = 'The first index is <0 or > pop_dimx'
         ENDIF
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 5:ikl - 1) .eq.'pop_t') then 
      IF (ianz.eq.2) then 
         IF( 0<kpara .and. kpara<=pop_dimx .and. kpara<=MAXDIMX) THEN
            IF( 0<kpara2 .and. kpara2<=pop_n .and. kpara2<=MAXPOP) THEN
               IF (ikl.gt.6) zeile (1:ikl - 6) = string (1:ikl - 6) 
               WRITE (zeile (ikl - 5:ikl + PREC_WIDTH-2) , PREC_F_REAL) pop_t ( kpara, kpara2)
               zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
            ELSE
               ier_num = -14
               ier_typ = ER_APPL
               ier_msg(1) = 'Possible reasons:'
               ier_msg(2) = 'The population has not been defined' 
               ier_msg(3) = 'The second index is <0 or > pop_n'
            ENDIF
         ELSE
            ier_num = -9
            ier_typ = ER_APPL
            ier_msg(1) = 'Possible reasons:'
            ier_msg(2) = 'The value of pop_dimx[1] is not defined' 
            ier_msg(3) = 'The first index is <0 or > pop_dimx'
         ENDIF
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSE 
      ier_num = - 2 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
!                                                                 
ELSEIF (lcomm.eq.6) then 
   IF (string (ikl - 6:ikl - 1) .eq.'diff_f') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.7) zeile (1:ikl - 7) = string (1:ikl - 7) 
         WRITE (zeile (ikl - 6:ikl + PREC_WIDTH-2) , PREC_F_REAL) diff_f
         zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 6:ikl - 1) .eq.'diff_k') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.7) zeile (1:ikl - 7) = string (1:ikl - 7) 
         WRITE (zeile (ikl - 6:ikl + PREC_WIDTH-2) , PREC_F_REAL) diff_k 
         zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 6:ikl - 1) .eq.'rvalue') then 
      IF (ianz.eq.1) then 
         IF (0<kpara .and. kpara<=pop_n .and. kpara<= MAXPOP) THEN
            IF (ikl.gt.7) zeile (1:ikl - 7) = string (1:ikl - 7) 
            WRITE (zeile (ikl - 6:ikl + PREC_WIDTH-2) , PREC_F_REAL) parent_val (kpara,0)
            zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
         ELSE
            ier_num = -9
            ier_typ = ER_APPL
            ier_msg(1) = 'The population has not been defined yet' 
         ENDIF
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 6:ikl - 1) .eq.'worstm') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.7) zeile (1:ikl - 7) = string (1:ikl - 7) 
         WRITE (zeile (ikl - 6:ikl + PREC_WIDTH-2) , PREC_F_INTE) pop_worst 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 6:ikl - 1) .eq.'worstr') then 
      IF (ianz.eq.1) then 
         IF (0<pop_worst .and. pop_worst<=pop_n .and. pop_worst<= MAXPOP) THEN
            IF (ikl.gt.7) zeile (1:ikl - 7) = string (1:ikl - 7) 
            WRITE (zeile (ikl - 6:ikl + PREC_WIDTH-2) , PREC_F_REAL) child_val (pop_worst,0)
            zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
         ELSE
            ier_num = -9
            ier_typ = ER_APPL
            ier_msg(1) = 'The population has not been defined yet' 
         ENDIF
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSE 
      ier_num = - 2 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
!
ELSEIF (lcomm.eq.7) then 
   IF (string (ikl - 7:ikl - 1) .eq.'diff_cr') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.8) zeile (1:ikl - 8) = string (1:ikl - 8) 
         WRITE (zeile (ikl - 7:ikl + PREC_WIDTH-2) , PREC_F_REAL) diff_cr 
         zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 7:ikl - 1) .eq.'diff_lo') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.8) zeile (1:ikl - 8) = string (1:ikl - 8) 
         WRITE (zeile (ikl - 7:ikl + PREC_WIDTH-2) , PREC_F_REAL) diff_local
         zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 7:ikl - 1) .eq.'pop_gen') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.8) zeile (1:ikl - 8) = string (1:ikl - 8) 
            WRITE (zeile (ikl - 7:ikl + PREC_WIDTH-2) , PREC_F_INTE) pop_gen 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 7:ikl - 1) .eq.'pop_sig') then 
      IF (ianz.eq.1) then 
         IF (0<kpara .and. kpara<=pop_dimx .and. kpara<= MAXDIMX) THEN
            IF (ikl.gt.8) zeile (1:ikl - 8) = string (1:ikl - 8) 
            WRITE (zeile (ikl - 7:ikl + PREC_WIDTH-2) , PREC_F_REAL) pop_sigma (kpara)
            zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
         ELSE
            ier_num = -14
            ier_typ = ER_APPL
            ier_msg(1) = 'Possible reasons:'
            ier_msg(2) = 'The value of pop_dimx[1] is not yet defined' 
            ier_msg(3) = 'The index is <0 or > pop_dimx'
         ENDIF
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
!
   ELSE 
      ier_num = - 2 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
!                                                                 
ELSEIF (lcomm.eq.8) then 
   IF (string (ikl - 8:ikl - 1) .eq.'diff_sel') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.8) zeile (1:ikl - 9) = string (1:ikl - 9) 
         WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_INTE) diff_sel_mode
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 8:ikl - 1) .eq.'pop_dimx') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.9) zeile (1:ikl - 9) = string (1:ikl - 9) 
            WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_INTE) pop_dimx 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 8:ikl - 1) .eq.'pop_lsig') then 
      IF (0<kpara .and. kpara<=pop_dimx .and. kpara<= MAXDIMX) THEN
         IF (ikl.gt.9) zeile (1:ikl - 9) = string (1:ikl - 9) 
         WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_REAL) pop_lsig ( kpara)
         zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
      ELSE
         ier_num = -14
         ier_typ = ER_APPL
         ier_msg(1) = 'Possible reasons:'
         ier_msg(2) = 'The value of pop_dimx[1] is not yet defined' 
         ier_msg(3) = 'The index is <0 or > pop_dimx'
      ENDIF
   ELSEIF (string (ikl - 8:ikl - 1) .eq.'pop_xmin') then 
      IF (0<kpara .and. kpara<=pop_dimx .and. kpara<= MAXDIMX) THEN
         IF (ikl.gt.9) zeile (1:ikl - 9) = string (1:ikl - 9) 
         WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_REAL) pop_xmin ( kpara)
         zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
      ELSE
         ier_num = -14
         ier_typ = ER_APPL
         ier_msg(1) = 'Possible reasons:'
         ier_msg(2) = 'The value of pop_dimx[1] is not yet defined' 
         ier_msg(3) = 'The index is <0 or > pop_dimx'
      ENDIF
   ELSEIF (string (ikl - 8:ikl - 1) .eq.'pop_xmax') then 
      IF (0<kpara .and. kpara<=pop_dimx .and. kpara<= MAXDIMX) THEN
         IF (ikl.gt.9) zeile (1:ikl - 9) = string (1:ikl - 9) 
         WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_REAL) pop_xmax ( kpara)
         zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
      ELSE
         ier_num = -14
         ier_typ = ER_APPL
         ier_msg(1) = 'Possible reasons:'
         ier_msg(2) = 'The value of pop_dimx[1] is not yet defined' 
         ier_msg(3) = 'The index is <0 or > pop_dimx'
      ENDIF
   ELSEIF (string (ikl - 8:ikl - 1) .eq.'pop_smin') then 
      IF (0<kpara .and. kpara<=pop_dimx .and. kpara<= MAXDIMX) THEN
         IF (ikl.gt.9) zeile (1:ikl - 9) = string (1:ikl - 9) 
         WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_REAL) pop_smin ( kpara)
         zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
      ELSE
         ier_num = -14
         ier_typ = ER_APPL
         ier_msg(1) = 'Possible reasons:'
         ier_msg(2) = 'The value of pop_dimx[1] is not yet defined' 
         ier_msg(3) = 'The index is <0 or > pop_dimx'
      ENDIF
   ELSEIF (string (ikl - 8:ikl - 1) .eq.'pop_smax') then 
      IF (0<kpara .and. kpara<=pop_dimx .and. kpara<= MAXDIMX) THEN
         IF (ikl.gt.9) zeile (1:ikl - 9) = string (1:ikl - 9) 
         WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_REAL) pop_smax ( kpara)
         zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
      ELSE
         ier_num = -14
         ier_typ = ER_APPL
         ier_msg(1) = 'Possible reasons:'
         ier_msg(2) = 'The value of pop_dimx[1] is not yet defined' 
         ier_msg(3) = 'The index is <0 or > pop_dimx'
      ENDIF
   ELSE 
      ier_num = - 2 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
!                                                                 
!                                                                 
ELSE 
   ier_num = - 2 
   ier_typ = ER_FORT 
ENDIF 
IF (ier_num.eq.0) then 
   ll = laenge+PREC_WIDTH - ltyp - (iklz - ikl + 1) 
   IF (iklz + 1.le.laenge) zeile (ikl + PREC_WIDTH-1:ll) = string (iklz + 1: laenge)                                                        
   string = zeile 
ELSE
ll = min (40, laenge)
    WRITE (ier_msg (1), '(a)') string (1:ll)
ENDIF 
ll = LEN_TRIM(string)
!
END SUBROUTINE diffev_ersetz_para                    
!*****7*****************************************************************
SUBROUTINE diffev_upd_para (ctype, ww, maxw, wert, ianz, cstring, substr) 
!-                                                                      
!       updates the parameter spezified by ctype, index ww  to the      
!       new value of wert                                               
!+                                                                      
USE diffev_allocate_appl
USE create_trial_mod
USE diff_evol
USE population
USE errlist_mod 
USE param_mod 
USE lib_errlist_func
USE lib_f90_allocate_mod
USE lib_upd_mod
USE variable_mod
USE precision_mod
!
IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER (LEN=* ), INTENT(IN   )    :: ctype 
INTEGER           , INTENT(IN   )    :: maxw
INTEGER           , INTENT(IN   )    :: ianz 
INTEGER           , INTENT(IN   )    :: ww (maxw)
REAL(KIND=PREC_DP), INTENT(IN   )    :: wert 
CHARACTER (LEN=* ), INTENT(IN   )    :: cstring
INTEGER, DIMENSION(2), INTENT(IN)    :: substr ! Indices of substring
!
INTEGER               :: i
INTEGER               :: pop_neu
REAL(kind=PREC_DP)    :: highest_r
!
CALL lib_upd_para (ctype, ww, maxw, wert, ianz, cstring, substr)
IF(ier_num==0 .OR. (ier_num==-40 .AND. ier_typ==ER_FORT)) RETURN
CALL no_error
!                                                                       
IF (ctype.eq.'pop_n') then 
   IF (ianz.eq.1) then 
      IF (nint (wert) .lt.4) then 
         ier_num = - 3
         ier_typ = ER_APPL 
!           ELSEIF (MAXPOP.lt.nint (wert) ) then 
!              ier_num = -2 
!              ier_typ = ER_APPL 
      ELSE 
         IF ( pop_gen == 0) THEN    ! We are still in generation zero
            pop_n = nint (wert) 
            var_val(var_ref+1) = pop_n     ! Update global user variable
            IF(NINT(wert) > MAXPOP) THEN
               CALL alloc_population( pop_n, MAXDIMX )
               CALL alloc_ref_para(MAXDIMX)
               IF(ier_num < 0) THEN
                  RETURN
               ENDIF
            ENDIF
         ELSE                       ! During refinement 
            IF ( pop_n > nint (wert) ) THEN ! New population is smaller
               pop_n = nint (wert) 
               var_val(var_ref+1) = pop_n     ! Update global user variable
               CALL write_genfile
            ELSE                    ! New population has increased, needs initialization
               IF(NINT(wert) > MAXPOP) THEN
                  CALL alloc_population( NINT(wert), MAXDIMX )
                  CALL alloc_ref_para(MAXDIMX)
                  IF(ier_num < 0) THEN
                     RETURN
                  ENDIF
               ENDIF
               pop_neu = nint (wert)
               highest_r = maxval ( parent_val(1:pop_n,0))
               parent_val(pop_n+1:pop_neu,0) = 10.*highest_r
               FORALL(i=pop_n+1:pop_neu) pop_x (:,i)     = pop_x (:,1)
               pop_n = pop_neu
               var_val(var_ref+1) = pop_n     ! Update global user variable
               CALL write_genfile
            ENDIF 
         ENDIF 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'pop_c') then 
   IF (ianz.eq.1) then 
      IF (nint (wert) .lt.4) then 
         ier_num = -3 
         ier_typ = ER_APPL 
!           ELSEIF (MAXPOP.lt.nint (wert) ) then 
!              ier_num = -2 
!              ier_typ = ER_APPL 
      ELSE 
         IF ( pop_gen == 0) THEN    ! We are still in generation zero
            pop_c = nint (wert) 
            var_val(var_ref+2) = pop_c     ! Update global user variable
            IF(NINT(wert) > MAXPOP) THEN
               CALL alloc_population( pop_c, MAXDIMX )
               CALL alloc_ref_para(MAXDIMX)
               IF(ier_num < 0) THEN
                  RETURN
               ENDIF
            ENDIF
         ELSE                       ! During refinement 
            IF ( pop_c > nint (wert) ) THEN ! New population is smaller
               pop_c = nint (wert) 
               var_val(var_ref+2) = pop_c     ! Update global user variable
               CALL create_trial                 ! Make a new set
               CALL write_genfile
            ELSEIF ( pop_c < nint (wert) ) THEN  ! New population has increased, needs initialization
               IF(NINT(wert) > MAXPOP) THEN
                  CALL alloc_population( NINT(wert), MAXDIMX )
                  CALL alloc_ref_para(MAXDIMX)
                  IF(ier_num < 0) THEN
                     RETURN
                  ENDIF
               ENDIF
               pop_c = nint (wert) 
               var_val(var_ref+2) = pop_c     ! Update global user variable
               CALL create_trial                 ! Make a new set
               CALL write_genfile                ! Write the "GENERATION" file
            ENDIF 
         ENDIF 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'pop_t') then 
   IF (ianz.eq.2) then 
      IF (0.lt.ww (1) .and.ww (1) .le.MAXPOP .and. ww(1)<=pop_c) then 
         IF (0.lt.ww (2) .and.ww (2) .le.MAXDIMX .and. ww(2)<=pop_dimx) then 
            pop_t (ww (1),ww(2) ) = wert 
         ELSE 
            ier_num = -14 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = -9 
         ier_typ = ER_APPL 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'pop_dimx') then 
   IF (ianz.eq.1) then 
      IF (nint (wert) .lt.1) then 
         ier_num = - 10 
         ier_typ = ER_APPL 
      ELSE 
         IF ( pop_gen > 0 ) THEN
            pop_dimx_new = .true.
         ENDIF
         pop_dimx = nint (wert) 
         IF(pop_dimx  > MAXDIMX) THEN
            CALL alloc_population( MAXPOP, pop_dimx )
            CALL alloc_ref_para(pop_dimx)
            IF(ier_num < 0) THEN
               RETURN
            ENDIF
         ENDIF
         DO i=NINT(wert), MAXDIMX
            IF(pop_name(i)==' ') THEN
               WRITE(pop_name(i),'(a4,i4.4)') 'PARA',i
            ENDIF
         ENDDO
         IF ( pop_gen == 0 ) THEN
            pop_dimx_new = .false.
         ENDIF
         pop_dimx_init = .true.      ! The dimension has been initialized in this run
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'pop_gen') then 
   IF (ianz.eq.1) then 
      IF (0.le.nint (wert) ) then 
         pop_gen = nint (wert) 
         var_val(var_ref+0) = pop_gen   ! Update global user variable
      ELSE 
         ier_num = 1 
         ier_typ = ER_APPL 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'pop_sig') then 
   IF (ianz.eq.1) then 
      IF (0.lt.ww (1) .and.ww (1) .le.MAXDIMX .and. ww(1)<=pop_dimx) then 
         pop_sigma (ww (1) ) = wert 
         pop_dimx_init = .TRUE.      ! The dimension hase been initialized in this run
      ELSE 
         ier_num = -14 
         ier_typ = ER_APPL 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'pop_lsig') then 
   IF (ianz.eq.1) then 
      IF (0.lt.ww (1) .and.ww (1) .le.MAXDIMX .and. ww(1)<=pop_dimx) then 
         pop_lsig (ww (1) ) = wert 
         pop_dimx_init = .TRUE.      ! The dimension hase been initialized in this run
      ELSE 
         ier_num = -14 
         ier_typ = ER_APPL 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'pop_xmin') then 
   IF (ianz.eq.1) then 
      IF (0.lt.ww (1) .and.ww (1) .le.MAXDIMX .and. ww(1)<=pop_dimx) then 
         pop_xmin (ww (1) ) = wert 
         pop_dimx_init = .TRUE.      ! The dimension hase been initialized in this run
      ELSE 
         ier_num = -14 
         ier_typ = ER_APPL 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'pop_xmax') then 
   IF (ianz.eq.1) then 
      IF (0.lt.ww (1) .and.ww (1) .le.MAXDIMX .and. ww(1)<=pop_dimx) then 
         pop_xmax (ww (1) ) = wert 
         pop_dimx_init = .TRUE.      ! The dimension hase been initialized in this run
      ELSE 
         ier_num = -14 
         ier_typ = ER_APPL 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'pop_smin') then 
   IF (ianz.eq.1) then 
      IF (0.lt.ww (1) .and.ww (1) .le.MAXDIMX .and. ww(1)<=pop_dimx) then 
         pop_smin (ww (1) ) = wert 
         pop_dimx_init = .TRUE.      ! The dimension hase been initialized in this run
      ELSE 
         ier_num = -14 
         ier_typ = ER_APPL 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'pop_smax') then 
   IF (ianz.eq.1) then 
      IF (0.lt.ww (1) .and.ww (1) .le.MAXDIMX .and. ww(1)<=pop_dimx) then 
         pop_smax (ww (1) ) = wert 
         pop_dimx_init = .TRUE.      ! The dimension hase been initialized in this run
      ELSE 
         ier_num = -14 
         ier_typ = ER_APPL 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'child_val') then 
   IF (ianz.eq.1) then 
      IF (0.lt.ww (1) .and.ww (1) .le.MAXPOP .and. ww(1)<=pop_c) then 
         child_val (ww (1),0 ) = wert 
      ELSE 
         ier_num = -14 
         ier_typ = ER_APPL 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'diff_cr') then 
   IF (ianz.eq.1) then 
      IF (0.le.wert.and.wert.le.1.0) then 
         diff_cr = wert 
      ELSE 
         ier_num = - 5 
         ier_typ = ER_APPL 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'diff_lo') then 
   IF (ianz.eq.1) then 
      IF (0.le.wert.and.wert.le.1.0) then 
         diff_local = wert 
      ELSE 
         ier_num = - 5 
         ier_typ = ER_APPL 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'diff_f') then 
   IF (ianz.eq.1) then 
      IF (0.lt.wert) then 
         diff_f = wert 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_APPL 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
ELSEIF (ctype.eq.'diff_k') then 
   IF (ianz.eq.1) then 
      IF (0.lt.wert) then 
         diff_k = wert 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_APPL 
      ENDIF 
   ELSE 
      ier_num = - 13 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
!
ELSE 
   ier_num = - 2 
   ier_typ = ER_FORT 
   WRITE (ier_msg (1), '(a)') ctype
ENDIF 
!
END SUBROUTINE diffev_upd_para                       
!*****7***************************************************************  
SUBROUTINE diffev_calc_intr_spec(string, line, ikl, iklz, ww, laenge, lp)
!-                                                                      
!     These are special intrinsic function for the DIFFEV. Any          
!     intrinsic function that references crystallographic values        
!     is found in this subroutine.                                      
!     Currently empty, needed for formal reasons.
!+                                                                      
!
USE berechne_mod
USE calc_expr_mod
USE do_read_number_mod
USE ersetz_mod
USE population
USE errlist_mod 
USE get_params_mod
USE lib_length
USE param_mod 
IMPLICIT none 
!                                                                       
!                                                                       
INTEGER, PARAMETER   :: MAXW = 9
!                                                                       
CHARACTER (LEN= * ), INTENT(INOUT) :: string
CHARACTER (LEN= * ), INTENT(INOUT) :: line 
INTEGER            , INTENT(IN   ) :: ikl
INTEGER            , INTENT(IN   ) :: iklz
INTEGER            , INTENT(INOUT) :: laenge
INTEGER            , INTENT(INOUT) :: lp
REAL(KIND=PREC_DP) , INTENT(INOUT) :: ww
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(string))), DIMENSION(1:MAXW) :: cpara
INTEGER            , DIMENSION(1:MAXW) :: lpara
REAL(kind=PREC_DP) , DIMENSION(1:MAXW) :: werte
CHARACTER(LEN=MAX(PREC_STRING,LEN(string))) :: parstring
INTEGER              :: ianz
INTEGER              :: i, j, lcomm
!                                                                       
!                                                                       
lcomm = length_com (string(1:ikl), ikl) 
ier_num = - 1 
ier_typ = ER_FORT 
DO i = 1, maxw 
   werte (i) = 0.0 
ENDDO 
!                                                                 
IF (lcomm.eq.10) then 
   IF(string (ikl - lcomm:ikl - 1) .eq.'par_number') THEN
      CALL get_params (line, ianz, cpara, lpara, 6, lp)
      IF (ier_num.eq.0) THEN
         IF (ianz==1) THEN
            IF(line(1:1) == '''' .AND. line(LEN_TRIM(line):LEN_TRIM(line)) == '''') THEN
               parstring = line(2:lp-1)
            ELSE
               CALL eval(line    , lp      )
               IF(ier_num == 0) THEN
                  parstring = line
               ENDIF
            ENDIF
            ww = -2 
            loop:DO i=1,pop_dimx
               IF(parstring==pop_name(i)) THEN
                  ww = i
                  EXIT loop
               ENDIF
            ENDDO loop
            IF(0<=ww .AND. ww <= pop_dimx) THEN
               CALL ersetz2 (string, ikl, iklz, ww, lcomm, laenge)
            ELSE
               ier_num = - 6
               ier_typ = ER_COMM
            ENDIF
         ELSE
            ier_num = - 6
            ier_typ = ER_COMM
         ENDIF
      ENDIF
   ELSE
      ier_num = - 3
      ier_typ = ER_FORT
   ENDIF
ELSEIF (lcomm.eq.8 ) then 
   IF(string (ikl - lcomm:ikl - 1) .eq.'par_name') THEN
      CALL get_params (line, ianz, cpara, lpara, MAXW, lp)
      IF (ier_num.eq.0) THEN
         IF (ianz==1) THEN
            werte(:) = 0
            CALL eval(cpara(1), lpara(1) )
            IF(ier_num==0) THEN
               j = NINT(do_read_number (cpara (1), lpara (1) ))
               IF(0<=j .AND. j<=pop_dimx) THEN
                  parstring = pop_name(j)
                  CALL ersetzc(string, ikl, iklz, parstring, LEN_TRIM(parstring), lcomm, laenge)
               ELSE
                  ier_num = - 6
                  ier_typ = ER_COMM
               ENDIF
            ELSE
               ier_num = - 6
               ier_typ = ER_COMM
            ENDIF
         ELSE
            ier_num = - 6
            ier_typ = ER_COMM
         ENDIF
      ENDIF
   ELSE
      ier_num = - 3
      ier_typ = ER_FORT
   ENDIF
ELSEIF (lcomm.eq.0) then 
   CALL ersetz2 (string, ikl, iklz, ww, 0, laenge) 
ELSE 
   ier_num = - 3 
   ier_typ = ER_FORT 
ENDIF 
!                                                                 
!IF (ier_num.ne.0) then 
!   WRITE ( *, * ) string (1:len_trim(string))
!   WRITE ( *, * ) line (1:len_trim(line))
!ENDIF 
!                                                                       
END SUBROUTINE diffev_calc_intr_spec                 
!
!*****7**************************************************************** 
!
SUBROUTINE diffev_calc_intr_log_spec(string, length)
!
IMPLICIT NONE
CHARACTER(LEN=*) , INTENT(INOUT) :: string
INTEGER          , INTENT(INOUT) :: length
!
END SUBROUTINE diffev_calc_intr_log_spec
!*****7**************************************************************** 
SUBROUTINE diffev_validate_var_spec (zeile, lp) 
!-                                                                      
!       checks whether the variable name is legal, DIFFEV specific part 
!                                                                       
!       Version : 1.0                                                   
!       Date    : 22 Oct 03                                             
!                                                                       
!       Author  : R.B. Neder  (reinhard.neder@fau.de)    
!+                                                                      
!
USE errlist_mod 
USE reserved_mod
IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER (LEN= * ), INTENT(IN   ) :: zeile 
INTEGER            , INTENT(IN   ) :: lp 
!                                                                       
INTEGER                            :: i , length
!                                                                       
!                                                                       
ier_num = 0 
ier_typ = ER_NONE 
!                                                                       
main: DO i = 1, diffev_reserved_n 
!  IF (index (diffev_reserved (i), zeile (1:lp) ) .ne.0) THEN 
   length = MAX(LEN_TRIM(diffev_reserved(i)), LEN_TRIM(zeile(1:lp)))
   length = MIN(length, LEN(diffev_reserved), LEN(zeile))
   IF(diffev_reserved (i)(1:length)== zeile(1:length) ) THEN           
      ier_num = - 25 
      ier_typ = ER_FORT 
      EXIT main
   ENDIF 
ENDDO  main
!                                                                       
END SUBROUTINE diffev_validate_var_spec              
!
!*******************************************************************************
!
SUBROUTINE diffev_get_var_type(line,length, var_is_type)
!
! Returns the variable type : INTEGER, REAL, CHARACTER, and Scalar versus field
!
USE constants_mod
USE lib_get_var
USE variable_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*)     , INTENT(IN)  :: line
INTEGER              , INTENT(IN)  :: length
INTEGER, DIMENSION(3), INTENT(OUT) :: var_is_type
!
!I
INTEGER, PARAMETER :: MAXPAR = 24
CHARACTER(LEN=16), DIMENSION(MAXPAR) :: diffev_names
INTEGER          , DIMENSION(MAXPAR) :: diffev_type
INTEGER          , DIMENSION(MAXPAR) :: diffev_dim
LOGICAL          , DIMENSION(MAXPAR) :: diffev_ro 
INTEGER :: i
!
DATA diffev_names  &
    /'child_val', 'pop_xmin ', 'pop_xmax ', 'pop_smin ', 'pop_smax ',  &
     'pop_lsig ', 'pop_dimx ', 'diff_sel ', 'pop_sig  ', 'pop_gen  ',  &
     'diff_lo  ', 'diff_cr  ', 'worstr   ', 'worstm   ', 'rvalue   ',  &
     'diff_k   ', 'diff_f   ', 'pop_v    ', 'pop_t    ', 'pop_n    ',  &
     'pop_c    ', 'bestr    ', 'bestm    ', 'p        '                &
    /
DATA diffev_type &
    /  IS_REAL ,   IS_REAL ,   IS_REAL ,   IS_REAL ,   IS_REAL , &
       IS_REAL ,   IS_REAL ,   IS_INTE ,   IS_REAL ,   IS_INTE , &
       IS_REAL ,   IS_REAL ,   IS_REAL ,   IS_INTE ,   IS_REAL , &
       IS_REAL ,   IS_REAL ,   IS_REAL ,   IS_REAL ,   IS_INTE , &
       IS_INTE ,   IS_REAL ,   IS_INTE ,   IS_REAL               &
    /
DATA diffev_dim  &
    /  IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  , &
       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  , &
       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  , &
       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_ARR  ,   IS_VEC  , &
       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC                &
    /
DATA diffev_ro  &
    /  .FALSE. ,   .FALSE. ,   .FALSE. ,   .FALSE. ,   .FALSE. , &
       .FALSE. ,   .FALSE. ,   .TRUE.  ,   .FALSE. ,   .FALSE. , &
       .FALSE. ,   .FALSE. ,   .TRUE.  ,   .TRUE.  ,   .TRUE.  , &
       .FALSE. ,   .FALSE. ,   .TRUE.  ,   .FALSE. ,   .FALSE. , &
       .FALSE. ,   .TRUE.  ,   .TRUE.  ,   .TRUE.                &
    /
!
var_is_type(:) = IS_UNKNOWN
!
main: DO i=1, MAXPAR
   IF(line(1:length) == diffev_names(i)(1:LEN_TRIM(diffev_names(i)))) THEN
      var_is_type(1) = diffev_type(i)
      var_is_type(2) = diffev_dim (i)
      IF(diffev_ro(i)) THEN
         var_is_type(3) = IS_READ
      ELSE
         var_is_type(3) = IS_WRITE
      ENDIF
      RETURN
   ENDIF
ENDDO main
!
CALL lib_get_var_type(line, length, var_is_type)
!
!
END SUBROUTINE diffev_get_var_typE
!
end module diffev_update_mod
