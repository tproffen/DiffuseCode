SUBROUTINE diffev_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz) 
!                                                                       
!-                                                                      
!       replaces a substring in an expression by the value of the       
!       appropriate Parameter.                                          
!     This is needed if the parameter is read.                          
!+                                                                      
USE diff_evol
USE population
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
   ELSEIF (string (ikl - 1:ikl - 1) .eq.'p') then 
      IF (ianz.eq.1) then 
         IF (0.le.kpara.and.kpara.le.MAXDIMX .and. kpara<=pop_dimx) THEN 
            WRITE (zeile (ikl - 1:ikl + 13) , '(e15.8e2)') pop_para (kpara)                                      
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
   ELSE 
      ier_num = - 2 
      ier_typ = ER_FORT 
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
ELSEIF (lcomm.eq.5) then 
   IF (string (ikl - 5:ikl - 1) .eq.'pop_n') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.6) zeile (1:ikl - 6) = string (1:ikl - 6) 
         WRITE (zeile (ikl - 5:ikl + 13) , '(i15    )') pop_n 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 5:ikl - 1) .eq.'pop_c') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.6) zeile (1:ikl - 6) = string (1:ikl - 6) 
         WRITE (zeile (ikl - 5:ikl + 13) , '(i15    )') pop_c 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 5:ikl - 1) .eq.'bestm') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.6) zeile (1:ikl - 6) = string (1:ikl - 6) 
         WRITE (zeile (ikl - 5:ikl + 13) , '(i15    )') pop_best 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 5:ikl - 1) .eq.'bestr') then 
      IF (ianz.eq.1) then 
         IF (0<pop_best .and. pop_best<=pop_n .and. pop_best<= MAXPOP) THEN
            IF (ikl.gt.6) zeile (1:ikl - 6) = string (1:ikl - 6) 
            WRITE (zeile (ikl - 5:ikl + 13) , '(e15.8e2)') child_val (pop_best)                                               
            zeile (ikl + 6:ikl + 6) = 'e' 
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
               WRITE (zeile (ikl - 5:ikl + 13) , '(e15.8e2)') child ( kpara, kpara2)                                           
               zeile (ikl + 6:ikl + 6) = 'e' 
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
               WRITE (zeile (ikl - 5:ikl + 13) , '(e15.8e2)') pop_t ( kpara, kpara2)                                           
               zeile (ikl + 6:ikl + 6) = 'e' 
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
         WRITE (zeile (ikl - 6:ikl + 13) , '(e15.8e2)') diff_f 
         zeile (ikl + 5:ikl + 5) = 'e' 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 6:ikl - 1) .eq.'diff_k') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.7) zeile (1:ikl - 7) = string (1:ikl - 7) 
         WRITE (zeile (ikl - 6:ikl + 13) , '(e15.8e2)') diff_k 
         zeile (ikl + 5:ikl + 5) = 'e' 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 6:ikl - 1) .eq.'rvalue') then 
      IF (ianz.eq.1) then 
         IF (0<kpara .and. kpara<=pop_n .and. kpara<= MAXPOP) THEN
            IF (ikl.gt.7) zeile (1:ikl - 7) = string (1:ikl - 7) 
            WRITE (zeile (ikl - 6:ikl + 13) , '(e15.8e2)') parent_val (kpara)                                                  
            zeile (ikl + 5:ikl + 5) = 'e' 
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
         WRITE (zeile (ikl - 6:ikl + 13) , '(i15    )') pop_worst 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 6:ikl - 1) .eq.'worstr') then 
      IF (ianz.eq.1) then 
         IF (0<pop_worst .and. pop_worst<=pop_n .and. pop_worst<= MAXPOP) THEN
            IF (ikl.gt.7) zeile (1:ikl - 7) = string (1:ikl - 7) 
            WRITE (zeile (ikl - 6:ikl + 13) , '(e15.8e2)') child_val (pop_worst)                                              
            zeile (ikl + 5:ikl + 5) = 'e' 
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
         WRITE (zeile (ikl - 7:ikl + 13) , '(e15.8e2)') diff_cr 
         zeile (ikl + 4:ikl + 4) = 'e' 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 7:ikl - 1) .eq.'diff_lo') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.8) zeile (1:ikl - 8) = string (1:ikl - 8) 
         WRITE (zeile (ikl - 7:ikl + 13) , '(e15.8e2)') diff_local                                               
         zeile (ikl + 4:ikl + 4) = 'e' 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 7:ikl - 1) .eq.'pop_gen') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.8) zeile (1:ikl - 8) = string (1:ikl - 8) 
            WRITE (zeile (ikl - 7:ikl + 13) , '(i15    )') pop_gen 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 7:ikl - 1) .eq.'pop_sig') then 
      IF (ianz.eq.1) then 
         IF (0<kpara .and. kpara<=pop_dimx .and. kpara<= MAXDIMX) THEN
            IF (ikl.gt.8) zeile (1:ikl - 8) = string (1:ikl - 8) 
            WRITE (zeile (ikl - 7:ikl + 13) , '(e15.8e2)') pop_sigma (kpara)                                                  
            zeile (ikl + 4:ikl + 4) = 'e' 
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
         WRITE (zeile (ikl - 8:ikl + 13) , '(i15    )') diff_sel_mode 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 8:ikl - 1) .eq.'pop_dimx') then 
      IF (ianz.eq.1) then 
         IF (ikl.gt.9) zeile (1:ikl - 9) = string (1:ikl - 9) 
            WRITE (zeile (ikl - 8:ikl + 13) , '(i15    )') pop_dimx 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_FORT 
         RETURN 
      ENDIF 
   ELSEIF (string (ikl - 8:ikl - 1) .eq.'pop_lsig') then 
      IF (0<kpara .and. kpara<=pop_dimx .and. kpara<= MAXDIMX) THEN
         IF (ikl.gt.9) zeile (1:ikl - 9) = string (1:ikl - 9) 
         WRITE (zeile (ikl - 8:ikl + 13) , '(e15.8e2)') pop_lsig ( kpara)                                                      
         zeile (ikl + 3:ikl + 3) = 'e' 
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
         WRITE (zeile (ikl - 8:ikl + 13) , '(e15.8e2)') pop_xmin ( kpara)                                                      
         zeile (ikl + 3:ikl + 3) = 'e' 
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
         WRITE (zeile (ikl - 8:ikl + 13) , '(e15.8e2)') pop_xmax ( kpara)                                                      
         zeile (ikl + 3:ikl + 3) = 'e' 
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
         WRITE (zeile (ikl - 8:ikl + 13) , '(e15.8e2)') pop_smin ( kpara)                                                      
         zeile (ikl + 3:ikl + 3) = 'e' 
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
         WRITE (zeile (ikl - 8:ikl + 13) , '(e15.8e2)') pop_smax ( kpara)                                                      
         zeile (ikl + 3:ikl + 3) = 'e' 
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
   ll = laenge+15 - ltyp - (iklz - ikl + 1) 
   IF (iklz + 1.le.laenge) zeile (ikl + 14:ll) = string (iklz + 1: laenge)                                                        
   string = zeile 
!ELSE 
!   WRITE ( *, * ) string 
ENDIF 
CALL rem_bl (string, ll) 
END SUBROUTINE diffev_ersetz_para                    
!*****7*****************************************************************
SUBROUTINE diffev_upd_para (ctype, ww, maxw, wert, ianz) 
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
USE lib_f90_allocate_mod
USE variable_mod
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
INTEGER               :: i
INTEGER               :: pop_neu
REAL                  :: highest_r
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
ELSEIF (ctype.eq.'pop_n') then 
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
               highest_r = maxval ( parent_val(1:pop_n))
               parent_val(pop_n+1:pop_neu) = 10.*highest_r
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
!           ELSEIF (MAXDIMX.lt.nint (wert) ) then 
!              ier_num = -4 
!              ier_typ = ER_APPL 
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
         IF ( pop_gen == 0 ) THEN
            pop_dimx_new = .false.
         ENDIF
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
         child_val (ww (1) ) = wert 
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
END SUBROUTINE diffev_upd_para                       
!*****7***************************************************************  
SUBROUTINE diffev_calc_intr_spec (string, line, ikl, iklz, ww, laenge, lp)                                                               
!-                                                                      
!     These are special intrinsic function for the DIFFEV. Any          
!     intrinsic function that references crystallographic values        
!     is found in this subroutine.                                      
!     Currently empty, needed for formal reasons.
!+                                                                      
!
USE errlist_mod 
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
END SUBROUTINE diffev_calc_intr_spec                 
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
IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER (LEN= * ), INTENT(IN   ) :: zeile 
INTEGER            , INTENT(IN   ) :: lp 
!                                                                       
INTEGER, PARAMETER                       :: reserved_n = 8
                                                                        
CHARACTER (LEN=12),DIMENSION(reserved_n) :: reserved = &
              (/'pop_gen     ', 'pop_n       ', 'pop_t       ', 'pop_dimx    ', &
                'pop_xmin    ', 'pop_xmax    ', 'diff_cr     ', 'diff_f      '/)                      
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
END SUBROUTINE diffev_validate_var_spec              
