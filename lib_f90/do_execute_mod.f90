MODULE do_execute_mod
!
CONTAINS
!
!*****7**************************************************************** 
!
SUBROUTINE do_execute (lreg, line, laenge) 
!-                                                                      
!     This subroutine increments the commands in the block structure.   
!     Block structure commands are executed in this subroutine,         
!     regular commands are returned to the calling subroutine.          
!+                                                                      
USE ber_params_mod
USE blanks_mod
USE doexec_mod 
USE doloop_mod 
USE errlist_mod 
USE get_params_mod
USE precision_mod
!
IMPLICIT none 
!
CHARACTER(LEN=*), INTENT(INOUT) :: line 
LOGICAL         , INTENT(OUT)   :: lreg 
INTEGER         , INTENT(INOUT) :: laenge
!                                                                       
INTEGER, PARAMETER :: maxw = 2 
!                                                                       
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
CHARACTER(LEN=1024) :: string
INTEGER, DIMENSION(MAXW) ::  lpara
INTEGER :: i, ianz, length
INTEGER :: ikl, ithen , idummy
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
!
IF(ier_ctrlc) THEN
   ier_num = -14
   ier_typ = ER_COMM
ENDIF
IF(ier_num/=0) RETURN
!                                                                       
!     LOGICAL if_test 
!                                                                       
lreg = .false. 
!                                                                       
ilevel (level) = ilevel (level) + 1 
line = do_comm (ilevel (level), level) 
laenge = do_leng (ilevel (level), level) 
IF(line(1:7)=='run_mpi') THEN
    level_mpi = level
   nlevel_mpi = ilevel(level)
ENDIF
IF (line (1:2) .eq.'&&') then 
   level = level + 1 
   READ (line (3:6), * ) jump (level) 
   ilevel (level) = jump (level) 
   line = do_comm (ilevel (level), level) 
   laenge = do_leng (ilevel (level), level) 
ENDIF 
!
CALL do_value(line,laenge)     ! Handle value() operation
!                                                                       
IF(line(1:4) ==  'else') THEN 
   IF (INDEX (line, 'if') >   0) THEN 
      CALL rem_insig_bl(line, laenge) 
!        ELSE
!           ier_num = - 31 
!           ier_typ = ER_FORT 
!           RETURN
   ENDIF 
ELSEIF(line(1:3) ==  'end') THEN 
   IF (INDEX (line, 'if') >   0) THEN 
      CALL rem_insig_bl(line, laenge) 
   ELSEIF (INDEX (line, 'do') >   0) THEN 
      CALL rem_insig_bl(line, laenge) 
   ENDIF
ENDIF 
!     if do-loop command, evaluate counter                              
!                                                                       
IF (line (1:3) .eq.'do ') then 
   IF (level.eq.0) then 
      jump (level) = ilevel (level) 
   ENDIF 
   CALL do_do (line, level, laenge) 
   IF (ier_num.ne.0) then 
      RETURN 
   ENDIF 
   IF (ldostart (level) ) then 
!     ----do loop is done                                               
      level = level - 1 
   ENDIF 
!                                                                       
!     enddo command, check for enddo until                              
!                                                                       
ELSEIF (line (1:5) .eq.'enddo') then 
   CALL do_end (line, level, laenge) 
   IF (ier_num.ne.0) then 
      RETURN 
   ENDIF 
   IF (ldostart (level) ) then 
      ilevel (level) = ilevel (level) + 1 
      level = level - 1 
   ELSE 
      DO while (do_comm (ilevel (level) , level) (1:3) .ne.'do ') 
         ilevel (level) = ilevel (level) - 1 
      ENDDO 
         ilevel (level) = jump (level) - 1 
   ENDIF 
!                                                                       
!     if or elseif command                                              
!                                                                       
ELSEIF (line (1:2) .eq.'if'.or.line (1:6) .eq.'elseif') then 
   IF (INDEX (line, 'then') .eq.0) then 
      ier_num = - 31 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
   CALL rem_insig_bl (line, laenge) 
   IF (ltest (level) ) then 
!........... previous block of the current if-elseif had been true      
      ltest (level) = .false. 
      level = level - 1 
   ELSE 
!...........This is the first if statement thats true                   
      ikl = INDEX (line, '(') 
      ithen = INDEX (line, 'then') - 1 
      string = line (ikl:ithen) 
      idummy = ithen - ikl + 1
      ltest (level) = if_test (string, idummy)
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      IF (ier_num.ne.0) goto 999 
      IF (.not.ltest (level) ) then 
!.............If statement is false, ignore the commands that follow    
         DO while(do_comm(ilevel(level)+1,level)(1:4).ne.'else'   &
             .and.do_comm(ilevel(level)+1,level)(1:5).ne.'endif')
         ilevel (level) = ilevel (level) + 1 
         ENDDO 
      ENDIF 
   ENDIF 
!                                                                       
!     else command                                                      
!                                                                       
ELSEIF (line (1:5) .eq.'else ') then 
   IF (ltest (level) ) then 
!...........A previous block of the currrent if-elseif had been true    
      ltest (level) = .false. 
      level = level - 1 
   ENDIF 
!                                                                       
!     elseif command                                                    
!                                                                       
ELSEIF (line.eq.'endif') then 
   ltest (level) = .false. 
   level = level - 1 
!                                                                       
!     break command                                                     
!                                                                       
ELSEIF (line (1:5) .eq.'break') then 
   length = laenge-5
   CALL get_params (line (6:laenge), ianz, cpara, lpara, maxw, length)
   IF (ier_num.eq.0) then 
      IF (ianz.eq.1) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.eq.0) then 
            IF (level.ge.nint (werte (1) ) - 1) then 
               DO i = level, level - nint (werte (1) ) + 1, - 1
                  ldostart (i) = .true. 
                  ltest (i) = .false. 
               ENDDO 
               level = level - nint (werte (1) ) 
            ELSE 
               ier_num = - 28 
               ier_typ = ER_FORT 
            ENDIF 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_FORT 
      ENDIF 
   ENDIF 
!                                                                       
!     regular command                                                   
!                                                                       
ELSE 
   lreg = .true. 
ENDIF 
!                                                                       
999 CONTINUE 
!                                                                       
END SUBROUTINE do_execute                     
!
!*****7**************************************************************** 
!
LOGICAL FUNCTION if_test (string, laenge) 
!-                                                                      
!     Tests the logical condition given in 'string'                     
!+                                                                      
USE berechne_mod
USE build_name_mod
USE calc_expr_mod
USE do_variable_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE search_string_mod
USE set_sub_generic_mod
!
IMPLICIT none 
!
INTEGER, PARAMETER :: MAXW = 20 
!
CHARACTER(LEN=*), INTENT(INOUT) ::  string 
INTEGER         , INTENT(INOUT) ::  laenge 
!
CHARACTER(LEN=2)  :: comp 
CHARACTER(LEN=1024) :: zeile, line, oldstr 
CHARACTER(LEN=1024) :: string1, string2 
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara  !(MAXW) 
INTEGER, DIMENSION(MAXW) :: lpara ! (MAXW) 
INTEGER :: ianz 
!     INTEGER suche_vor, suche_nach , suche_vor_hoch
INTEGER :: icom, iz1, iz2 , iz3
INTEGER :: ic1, ic2, ic3
INTEGER :: iper
INTEGER :: lcom    ! comparator length
INTEGER :: ikl, iklz, ikla, ikla1, ikla2 
INTEGER :: ll, i, inot, lll 
INTEGER :: istring1, istring2 
INTEGER :: istring1_len 
INTEGER :: istring2_len 
INTEGER :: omask ! , nmask might be needed later
LOGICAL :: lscr, lscr1 
LOGICAL :: lstring1, lstring2 
INTEGER :: ios
REAL(KIND=PREC_DP) :: werte (MAXW) 
REAL(KIND=PREC_DP) ::  w1, w2 
LOGICAL  , DIMENSION(1024,0:1) :: lmask
!                                                                       
      lmask = .TRUE.
      ier_num = 0 
      ier_typ = ER_NONE 
      oldstr = string 
      w1 = 0
      w2 = 0
      lstring1 = .FALSE.
      lstring2 = .FALSE.
      CALL ersetz_variable(string, laenge, lmask, omask)
!     CALL rem_bl (string, laenge) 
      if_test = .false. 
!     IF (laenge.eq.0.or.string.eq.' '.or.ier_num.ne.0) then 
      IF (laenge <= 2.or.string == ' '.or.ier_num /= 0) THEN 
         ier_num = -12
         ier_typ = ER_FORT
         WRITE (ier_msg (2), '(a41)') line (1:41) 
      ELSE 
         icom = MAX (INDEX (string, '.lt.', .TRUE.) , INDEX (string, '.le.', .TRUE.) ,  &
                     INDEX (string, '.gt.', .TRUE.) , INDEX (string, '.ge.', .TRUE.) ,  &
                     INDEX (string, '.eq.', .TRUE.) , INDEX (string, '.ne.', .TRUE.) ,  &
                     INDEX (string, '<',    .TRUE.) , INDEX (string, '<=',   .TRUE.) ,  &
                     INDEX (string, '>',    .TRUE.) , INDEX (string, '<=',   .TRUE.) ,  &
                     INDEX (string, '==',   .TRUE.) , INDEX (string, '/=',   .TRUE.)    &
                    )
icom_sig:         DO while (icom /= 0) 
!                                                                       
!     --Found an operator, search for numbers before and after          
!                                                                       
            lcom = 4
            IF(string(icom:icom)=='.') THEN
               comp = string (icom + 1:icom + 2) 
               lcom = 4
            ELSEIF(string(icom+1:icom+1) == '=') THEN
               comp = string(icom:icom+1)
               lcom = 2
            ELSEIF(string(icom:icom)=='<' .OR. string(icom:icom)=='>') THEN
               comp = string(icom:icom)
               lcom = 1
            ENDIF
         ic1 = 0
         ic2 = 0
         lll = icom - 1 
         iz1 = suche_vor (string (1:icom - 1), lll) 
         IF(iz1>2) THEN
!
!      -- Search for a comma, to detect expressions like "%c",variable
!
            IF(string(iz1-1:iz1-1)==',') THEN  !Found comma, get string expression
               lll = iz1-2
               ic1 = suche_vor_hoch (string (1:iz1  - 2), lll)   ! Search the "
               iz1 = ic1                       ! Needed for ersetzt_log
               istring1 = ic1
               string1  = string(ic1:icom-1)
               istring1_len = icom - ic1
               lstring1 = .true. 
            ENDIF
         ENDIF
         IF(ic1==0) THEN                       ! No " found, regular string or number?
         zeile = '('//string (iz1:icom - 1) //')' 
         ll = icom - 1 - iz1 + 3 
         istring1 = INDEX (zeile, '''') 
         lstring1 = .false. 
         ENDIF

         IF (istring1.gt.1) THEN
            IF(ic1==0) THEN    ! Normal single string
!                                                                       
!     ----found a string variable                                       
!                                                                       
            istring2 = INDEX (zeile (istring1 + 1:ll) , '''') + istring1                                                    
            IF (istring2.eq.istring1) then 
!                                                                       
!     ----Missing second ', first ' is there                            
!                                                                       
               ier_num = - 21 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
            string1 = zeile (istring1 + 1:istring2 - 1) 
            istring1_len = istring2 - istring1 - 1 
            lstring1 = .true. 
            ENDIF 
         ELSE 
            w1 = berechne (zeile, ll) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ENDIF 
!        lll = (laenge) - (icom + 4) + 1 
!        ic2 = INDEX(string(icom+4:laenge),'"')   ! Search for "  in post string
         lll = (laenge) - (icom + lcom) + 1 
         ic2 = INDEX(string(icom+lcom:laenge),'"')   ! Search for "  in post string
         IF(ic2 > 0 ) THEN                        ! Found a "
!          zeile = string(icom+4+ic2:laenge)
!          lll = laenge-(icom+4+ic2) + 1
           zeile = string(icom+lcom+ic2:laenge)
           lll = laenge-(icom+lcom+ic2) + 1
           ic3 = INDEX(zeile(1:lll),'"')          ! Search closing "
           iper = 0
           DO i=1, ic3-1                          ! Count %
              IF(zeile(i:i).eq.'%') THEN
                 iper = iper + 1
              ENDIF
           ENDDO
!          ic3 = icom + 4 + ic2-1 + ic3 +1
           ic3 = icom + lcom + ic2-1 + ic3 +1
           DO i=1,iper                            ! jump to last argument
              ic3 = ic3 + INDEX(string(ic3:laenge),',') 
           ENDDO
           lll = (laenge) - ic3 + 1
           iz3 = suche_nach (string (ic3:laenge), lll) 
!          iz2 = ic3 + iz3 -icom - 4              ! Needed for ersetz_log
!          string2 = string(icom+4+ic2-1:ic3+iz3-1)  ! Cut the string expression
!          istring2_len = (ic3+iz3-1)-(icom+4+ic2-1) + 1
           iz2 = ic3 + iz3 -icom - lcom           ! Needed for ersetz_log
           string2 = string(icom+lcom+ic2-1:ic3+iz3-1)  ! Cut the string expression
           istring2_len = (ic3+iz3-1)-(icom+lcom+ic2-1) + 1
           lstring2 = .true.
         ENDIF
         IF(ic2==0) THEN                          ! No " found, regular string or number
!           iz2 = suche_nach (string (icom + 4:laenge), lll)
!           zeile = '('//string (icom + 4:icom + 4 + iz2 - 1) //')' 
!           ll = icom + 4 + iz2 - 1 - (icom + 4) + 3 
            iz2 = suche_nach (string (icom + lcom:laenge), lll)
            zeile = '('//string (icom + lcom:icom + lcom + iz2 - 1) //')' 
            ll = icom + lcom + iz2 - 1 - (icom + lcom) + 3 
            istring1 = INDEX (zeile, '''') 
            lstring2 = .false. 
!        ENDIF
            IF (istring1.gt.1) then 
!                                                                       
!     ----found a string variable                                       
!                                                                       
               istring2 = INDEX (zeile (istring1 + 1:ll) , '''') + istring1                                                    
               IF (istring2.eq.istring1) then 
!                                                                       
!     ----Missing second ', first ' is there                            
!                                                                       
                  ier_num = - 21 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
               string2 = zeile (istring1 + 1:istring2 - 1) 
               istring2_len = istring2 - istring1 - 1 
               lstring2 = .true. 
            ELSE 
               w2 = berechne (zeile, ll) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
            ENDIF 
         ENDIF
         IF (lstring1.and.lstring2) then   ! Compare two strings
            CALL get_params (string1 (1:istring1_len), ianz, cpara,     &
            lpara, maxw, istring1_len)                                  
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
            string1      = cpara (1) (1:lpara (1) ) 
            istring1_len = lpara (1) 
            CALL get_params (string2 (1:istring2_len), ianz, cpara,     &
            lpara, maxw, istring2_len)                                  
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
            string2      = cpara (1) (1:lpara (1) ) 
            istring2_len = lpara (1) 
!                                                                       
            IF (comp.eq.'eq'.OR.comp=='==') then 
               lscr = string1 (1:istring1_len) .eq.string2 (1:          &
               istring2_len)                                            
            ELSEIF (comp.eq.'ne'.OR. comp == '/=' ) then 
               lscr = string1 (1:istring1_len) .ne.string2 (1:          &
               istring2_len)                                            
            ELSEIF (comp.eq.'lt'.OR. comp == '<'  ) then 
               lscr = string1 (1:istring1_len) .lt.string2 (1:          &
               istring2_len)                                            
            ELSEIF (comp.eq.'le'.OR. comp == '<=' ) then 
               lscr = string1 (1:istring1_len) .le.string2 (1:          &
               istring2_len)                                            
            ELSEIF (comp.eq.'gt'.OR. comp == '>'  ) then 
               lscr = string1 (1:istring1_len) .gt.string2 (1:          &
               istring2_len)                                            
            ELSEIF (comp.eq.'ge'.OR. comp == '>=' ) then 
               lscr = string1 (1:istring1_len) .ge.string2 (1:          &
               istring2_len)                                            
            ENDIF 
!
         ELSE 
            IF (    comp.eq.'lt' .OR. comp == '<'  ) then 
               lscr = w1.lt.w2 
            ELSEIF (comp.eq.'le' .OR. comp == '<=' ) then 
               lscr = w1.le.w2 
            ELSEIF (comp.eq.'gt' .OR. comp == '>'  ) then 
               lscr = w1.gt.w2 
            ELSEIF (comp.eq.'ge' .OR. comp == '>=' ) then 
               lscr = w1.ge.w2 
            ELSEIF (comp.eq.'eq' .OR. comp == '==') then 
               lscr = w1.eq.w2 
            ELSEIF (comp.eq.'ne' .OR. comp == '/=' ) then 
               lscr = w1.NE.w2 
            ENDIF 
         ENDIF 
!        lcom = 4 
         CALL ersetz_log (string, iz1, iz2, icom, laenge, lcom, lscr) 
         icom = MAX (INDEX (string, '.lt.') , INDEX (string, '.le.') ,  &
                     INDEX (string, '.gt.') , INDEX (string, '.ge.') ,  &
                     INDEX (string, '.eq.') , INDEX (string, '.ne.') ,  &
                     INDEX (string, '<')    , INDEX (string, '<=')   ,  &
                     INDEX (string, '>')    , INDEX (string, '<=')   ,  &
                     INDEX (string, '==')   , INDEX (string, '/=')      &
                    )
         ENDDO  icom_sig

         IF(laenge>3) THEN    ! String is long enough for a logical function
            CALL calc_intr_log(string,laenge)
            CALL p_calc_intr_log_spec(string,laenge)
!           ikla = 1 + INDEX (string(2:laenge), '(') 
!           iklz = ikla + INDEX(string(ikla+1:laenge), ')')
!           line = string(1:laenge-1)
!           CALL calc_intr_log (string, ikla, iklz, laenge)
!           IF(ier_num /= -1 .AND. ier_typ == ER_FORT .AND. &
!              (INDEX(string,'.not.')>0 .OR. INDEX(string,'.and.')>0 .OR. &
!               INDEX(string,'.or.' )>0 .OR. INDEX(string,'.xor.')>0 .OR. &
!               INDEX(string,'.eqv.')>0                                  )&
!              ) THEN   
!              ier_num = 0    ! Ignore errors from calc_intr
!              ier_typ = 0
!           ENDIF
         ENDIF
         ikla = INDEX (string, '(') 
         DO while (ikla.ne.0) 
         iklz = INDEX (string (ikla + 1:laenge) , ')') + ikla 
         IF (iklz.eq.ikla) then 
            ier_num = - 9 
            ier_typ = ER_FORT 
            if_test = .false. 
            RETURN 
         ENDIF 
         ikla2 = INDEX (string (ikla + 1:iklz) , '(') + ikla 
         ikla1 = ikla 
         DO while (ikla2.lt.iklz.and.ikla2.gt.ikla1) 
         ikla1 = ikla2 
         ikla2 = INDEX (string (ikla1 + 1:iklz) , '(') + ikla1 
         ENDDO 
         ikl = max (ikla1, ikla2) 
         IF (ikl.ne.0) then 
            line = string (ikl + 1:iklz - 1) 
            ll = iklz - ikl - 1 
!                                                                       
!     ----Found a set of corresponding brackets                         
!                                                                       
!         Evaluate any '.not.'                                          
!                                                                       
            inot = INDEX (line, '.not.') 
            DO while (inot.ne.0) 
            ios = 0
            READ (line (inot + 5:inot + 5) , '(l1)',IOSTAT=ios) lscr 
            IF(ios/=0) THEN
               ier_num = -18
               ier_typ = ER_FORT
               RETURN
            ENDIF
            lscr = .not.lscr 
            iz1 = max (inot - 1, 1) 
            iz2 = 1 
            lcom = 5 
            CALL ersetz_log (line, iz1, iz2, inot, ll, lcom, lscr) 
            inot = INDEX (line, '.not.') 
            ENDDO 
!                                                                       
!         Evaluate any '.and.'                                          
!                                                                       
            inot = INDEX (line, '.and.') 
            DO while (inot.ne.0) 
            ios = 0
            READ (line (inot - 1:inot - 1) , '(l1)') lscr1 
            IF(ios/=0) THEN
               ier_num = -18
               ier_typ = ER_FORT
               RETURN
            ENDIF
            READ (line (inot + 5:inot + 5) , '(l1)') lscr 
            IF(ios/=0) THEN
               ier_num = -18
               ier_typ = ER_FORT
               RETURN
            ENDIF
            lscr = lscr1.and.lscr 
            iz1 = inot - 1 
            iz2 = 1 
            lcom = 5 
            CALL ersetz_log (line, iz1, iz2, inot, ll, lcom, lscr) 
            inot = INDEX (line, '.and.') 
            ENDDO 
!                                                                       
!         Evaluate any '.eqv.'                                          
!                                                                       
            inot = INDEX (line, '.eqv.') 
            DO while (inot.ne.0) 
            ios = 0
            READ (line (inot - 1:inot - 1) , '(l1)',IOSTAT=ios) lscr1 
            IF(ios/=0) THEN
               ier_num = -18
               ier_typ = ER_FORT
               RETURN
            ENDIF
            READ (line (inot + 5:inot + 5) , '(l1)',IOSTAT=ios) lscr 
            IF(ios/=0) THEN
               ier_num = -18
               ier_typ = ER_FORT
               RETURN
            ENDIF
            lscr = (lscr1.and.lscr) .or. (.not.lscr1.and..not.lscr) 
            iz1 = inot - 1 
            iz2 = 1 
            lcom = 5 
            CALL ersetz_log (line, iz1, iz2, inot, ll, lcom, lscr) 
            inot = INDEX (line, '.eqv.') 
            ENDDO 
!                                                                       
!         Evaluate any '.xor.'                                          
!                                                                       
            inot = INDEX (line, '.xor.') 
            DO while (inot.ne.0) 
            ios = 0
            READ (line (inot - 1:inot - 1) , '(l1)',IOSTAT=ios) lscr1 
            IF(ios/=0) THEN
               ier_num = -18
               ier_typ = ER_FORT
               RETURN
            ENDIF
            READ (line (inot + 5:inot + 5) , '(l1)',IOSTAT=ios) lscr 
            IF(ios/=0) THEN
               ier_num = -18
               ier_typ = ER_FORT
               RETURN
            ENDIF
            lscr = (lscr1.and..not.lscr) .or. (.not.lscr1.and.lscr) 
            iz1 = inot - 1 
            iz2 = 1 
            lcom = 5 
            CALL ersetz_log (line, iz1, iz2, inot, ll, lcom, lscr) 
            inot = INDEX (line, '.xor.') 
            ENDDO 
!                                                                       
!         Evaluate any '.or.'                                           
!                                                                       
            inot = INDEX (line, '.or.') 
            DO while (inot.ne.0) 
            ios = 0
            READ (line (inot - 1:inot - 1) , '(l1)',IOSTAT=ios) lscr1 
            IF(ios/=0) THEN
               ier_num = -18
               ier_typ = ER_FORT
               RETURN
            ENDIF
            READ (line (inot + 4:inot + 4) , '(l1)',IOSTAT=ios) lscr 
            IF(ios/=0) THEN
               ier_num = -18
               ier_typ = ER_FORT
               RETURN
            ENDIF
            lscr = lscr1.or.lscr 
            iz1 = inot - 1 
            iz2 = 1 
            lcom = 4 
            CALL ersetz_log (line, iz1, iz2, inot, ll, lcom, lscr) 
            inot = INDEX (line, '.or.') 
            ENDDO 
            zeile = ' ' 
            IF (ikl.gt.1) zeile (1:ikl - 1) = string (1:ikl - 1) 
            zeile (ikl:ikl + ll - 1) = line (1:ll) 
            lll = ikl + ll - 1 
            IF (iklz + 1.le.laenge) then 
               zeile (ikl + ll:ikl + ll + laenge-iklz) =  &
                    string (iklz +  1:laenge)
               lll = ikl + ll + laenge-iklz 
            ENDIF 
            string = zeile 
            laenge = lll 
         ENDIF 
         ikla = INDEX (string, '(') 
         ENDDO 
         IF(string/=' ') THEN
         i = 1 
         DO while (string (i:i) .eq.' ') 
         i = i + 1 
         ENDDO 
         zeile = string (i:laenge) 
         string = zeile 
         ENDIF
         IF (string.eq.'T'.or.string.eq.'t') then 
            if_test = .true. 
         ELSEIF (string.eq.'F'.or.string.eq.'f') then 
            if_test = .false. 
         ELSE 
            ier_num = - 18 
            ier_typ = ER_FORT 
            if_test = .false. 
         ENDIF 
      ENDIF 
!                                                                       
END FUNCTION if_test                          
!
!*****7**************************************************************** 
!
SUBROUTINE ersetz_log (string, iz1, iz2, icom, laenge, lcom, lscr) 
!-                                                                      
!       Replaces the result of a logical operation within the           
!       string                                                          
!                                                                       
!       Version  : 0.0                                                  
!                                                                       
!       Date     : September 1992                                       
!                                                                       
!       Author   : r.n.                                                 
!                                                                       
!       modified :                                                      
!                                                                       
!+                                                                      
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) string 
      CHARACTER(1024) zeile 
      INTEGER iz1, iz2, icom, laenge, lcom, ll 
      LOGICAL lscr 
!                                                                       
      zeile = ' ' 
      IF (iz1.gt.1) zeile (1:iz1 - 1) = string (1:iz1 - 1) 
      WRITE (zeile (iz1:iz1) , '(L1)') lscr 
      ll = laenge- (lcom + iz2 + icom - iz1) + 1 
      IF (icom + lcom + iz2.le.laenge) then 
         zeile (iz1 + 1:ll) = string (icom + lcom + iz2:laenge) 
      ENDIF 
      string = zeile 
      laenge = ll 
!                                                                       
END SUBROUTINE ersetz_log                     
!
!*****7**************************************************************** 
!
SUBROUTINE do_do (line, level, laenge) 
!-                                                                      
!     reads the 'do' command and evaluates the corresponding counter    
!+                                                                      
USE ber_params_mod
USE doloop_mod 
USE do_variable_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE set_sub_generic_mod
USE precision_mod
USE variable_mod
!
IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 3) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) zeile, cpara (maxw) , cdummy
      INTEGER lpara (maxw) 
      INTEGER ipos, ikp, ianz, level, laenge, lll 
      INTEGER ianz_d, i 
INTEGER :: idummy
INTEGER, DIMENSION(2) :: substr = (/0,VAR_CLEN/)    ! Dummy substring indices
!     LOGICAL if_test 
      LOGICAL l_var 
      REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
REAL(KIND=PREC_DP) :: wert
!                                                                       
cdummy = ' '
      ier_num = - 6 
      ier_typ = ER_FORT 
!                                                                       
!     search for argument separator on the do-loop command line         
!                                                                       
      ipos = INDEX (line, '=') 
      ikp = INDEX (line, '[') 
!                                                                       
!     Do-loop of type: do counter = start,end[,increment]               
!                                                                       
      IF (ipos.gt.0) then 
!                                                                       
!       The counter is a user defined variable name ?                   
!                                                                       
         l_var = (ikp.eq.0.or.ipos.lt.ikp) 
         lll = ipos - ikp - 2 
         IF (.not.l_var) then 
            IF (lll.gt.0) then 
               zeile = ' ' 
               zeile (1:lll) = line (ikp + 1:ipos - 2) 
               CALL get_params (zeile, ianz_d, cpara, lpara, maxw, lll) 
               IF (ier_num.eq.0) then 
                  CALL ber_params (ianz_d, cpara, lpara, werte, maxw) 
                  IF (ier_num.ne.0) then 
                     RETURN 
                  ENDIF 
                  IF (ianz_d.ge.1.and.ianz_d.le.2) then 
                     DO i = 1, ianz_d 
                     do_kpara (i) = nint (werte (i) ) 
                     ENDDO 
                  ELSE 
                     ier_num = - 17 
                     ier_typ = ER_FORT 
                     RETURN 
                  ENDIF 
               ELSE 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 14 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ENDIF 
         ier_num = 0 
         ier_typ = ER_NONE 
         IF (ldostart (level) ) then 
!.........Do loop needs to be started                                   
            IF (ipos.ne.0) then 
               lll = laenge- (ipos + 1) + 1 
               CALL get_params (line (ipos + 1:laenge), ianz, cpara,    &
                                lpara, maxw, lll)                                        
               IF (ier_num.eq.0.and.ianz.ge.2) then 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.ne.0) then 
                     RETURN 
                  ENDIF 
                  glow (level) = werte (1) 
                  ghigh (level) = werte (2) 
                  IF (ianz.eq.2) then 
                     ginc (level) = 1.0 
                  ELSE 
                     ginc (level) = werte (3) 
                  ENDIF 
                  nloop (level) = max (int ( (ghigh (level) - glow (    &
                  level) ) / ginc (level) ) + 1, 0)                     
                  iloop (level) = 0 
                  ldostart (level) = .false. 
                  ier_num = 0 
                  ier_typ = ER_NONE 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_FORT 
                  GOTO 999 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               GOTO 999 
            ENDIF 
         ENDIF 
         iloop (level) = iloop (level) + 1 
         IF (iloop (level) .le.nloop (level) ) then 
!.........Do loop is running                                            
            wert = glow (level) + (iloop (level) - 1) * ginc (level) 
            IF (l_var) then 
               CALL upd_variable (line (4:ipos - 1), ipos - 4, wert,    &
               cpara (1), lpara (1), substr )
            ELSE 
               CALL p_upd_para (line (4:ikp - 1), do_kpara, 1, wert, ianz_d, cdummy, substr)
            ENDIF 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE 
!.........Do loop is done                                               
            ldostart (level) = .true. 
         ENDIF 
!                                                                       
!     do...until loop                                                   
!                                                                       
      ELSEIF (line (4:4) .eq.'    ') then 
         ier_num = 0 
         ier_typ = ER_NONE 
         ldostart (level) = .false. 
!                                                                       
!     DO WHILE loop                                                     
!                                                                       
      ELSEIF (laenge.gt.4.and.INDEX(line (4:laenge),'while') .ne.0) then
         ier_num = 0 
         ier_typ = ER_NONE 
         ipos = INDEX (line, '(') 
         zeile = line (ipos:laenge) 
         idummy = laenge-ipos + 1
         ldostart (level) = .not.if_test (zeile, idummy)
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
!                                                                       
!     error, wrong parameter on the do loop command line                
!                                                                       
      ELSE 
         ier_num = - 8 
         ier_typ = ER_COMM 
         ldostart (level) = .true. 
      ENDIF 
!                                                                       
  999 CONTINUE 
!                                                                       
      END SUBROUTINE do_do                          
!
!*****7**************************************************************** 
!
SUBROUTINE do_end (line, level, laenge) 
!
!-                                                                      
!     reads the 'enddo'. If the 'until' is found, the expression is     
!     evaluated                                                         
!+                                                                      
USE doloop_mod 
USE errlist_mod 
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(LEN=*), INTENT(INOUT) :: line 
INTEGER         , INTENT(INOUT) :: laenge
INTEGER         , INTENT(INOUT) :: level
!
CHARACTER(LEN=1024) :: zeile 
INTEGER :: ipos
INTEGER :: idummy
!     LOGICAL if_test 
!                                                                       
ipos = INDEX (line, 'until') 
IF (ipos.ne.0) then 
   zeile = line (ipos + 5:laenge) 
   idummy = laenge- (ipos + 5) + 1
   ldostart (level) = if_test (zeile, idummy)
   IF (ier_num.ne.0) then 
      RETURN 
   ENDIF 
ENDIF 
!
END SUBROUTINE do_end                         
!
!*****7**************************************************************** 
!
SUBROUTINE calc_intr_log (string, length)
!                                                                       
!     Evaluate all intrinsic logical functions within string
!                                                                       
      USE errlist_mod 
      USE ersetzl_mod
      USE variable_mod
!
      CHARACTER (LEN=* ), INTENT(INout) :: string
      INTEGER           , INTENT(INOUT) :: length
!
      INTEGER, PARAMETER :: N_LF = 2
!
      CHARACTER (LEN=5), DIMENSION(1:N_LF) :: f_names
      INTEGER          , DIMENSION(1:N_LF) :: l_names
      INTEGER :: j
      INTEGER :: ikl, iklz
      INTEGER :: ihyp, ihyp2
      INTEGER :: n_isexp
      LOGICAL :: lres
!
!     LOGICAL :: is_expression
!
      DATA f_names /'isexp', 'isvar' /
      DATA l_names / 5     ,  5      /
!
ier_num = 0
ier_typ = ER_NONE
!
DO j=1,N_LF                              ! Loop over all defined functions
   any_isexp: DO                         ! Search unitil no more found
      n_isexp = INDEX(string,f_names(j))
      IF(n_isexp > 0 ) THEN              ! We found a function
         ikl  = n_isexp + INDEX(string(n_isexp+1:length),'(')
         IF(ikl > n_isexp) THEN          ! We found an opening '('
            ihyp = ikl + INDEX(string(ikl+1:length),'''')
            IF(ihyp > ikl) THEN          ! We found an opening ''''
               ihyp2 = ihyp + INDEX(string(ihyp+1:length),'''')
               IF(ihyp2 > ihyp) THEN     ! We found an closing ''''
                  iklz = ihyp2 + INDEX(string(ihyp2+1:length-1),')')
                  IF(iklz > ihyp2) THEN  ! We found an closing ')'
                  SELECT CASE (j)
                  CASE (1)
                     lres = is_expression(string(ikl+1:iklz-1))
                  CASE (2)
                     lres = is_variable(string(ikl+1:iklz-1))
                  END SELECT
                  CALL ersetzl (string, ikl, iklz, lres, l_names(j), length)
                  ELSE
                     ier_num = -11
                     ier_typ = ER_FORT
                     RETURN
                  ENDIF
               ELSE
                  ier_num = -38
                  ier_typ = ER_FORT
                  RETURN
               ENDIF
            ELSE
               ier_num = -38
               ier_typ = ER_FORT
               RETURN
            ENDIF
         ELSE
            ier_num = -11
            ier_typ = ER_FORT
            RETURN
         ENDIF
      ELSE
         EXIT any_isexp
      ENDIF
   ENDDO any_isexp
ENDDO
!
END SUBROUTINE calc_intr_log
!
!*****7**************************************************************** 
!
LOGICAL FUNCTION is_expression(string)
!
      USE berechne_mod
      USE calc_expr_mod
      USE errlist_mod
!
      CHARACTER(LEN=*), INTENT(IN) :: string
!
      CHARACTER(LEN=1024) :: line
      INTEGER             :: ihyp, ihyp2, i1, i2, length
      REAL                :: ww
!
      i1 = 1
      i2 = LEN_TRIM(string)
      ihyp = MAX (INDEX (string, '''') , INDEX (string, '"') )
         IF(ihyp > 0) THEN
         i1 = ihyp+1
         ihyp2 = ihyp + MAX (INDEX (string(ihyp+1:i2), '''') ,  &
                             INDEX (string(ihyp+1:i2), '"')   )
         i2 = ihyp2-1
      ENDIF
      is_expression = .FALSE.
      IF(i2 >= i1) THEN
         line = '(' // string(i1:i2) // ')'
         length = LEN_TRIM(line)
         ww = berechne(line, length)
         IF(ier_num == 0) then
            is_expression = .TRUE.
         ELSE
            is_expression = .FALSE.
            CALL no_error
         ENDIF
      ENDIF
!
END FUNCTION is_expression
!
!*****7**************************************************************** 
!
SUBROUTINE do_value(line, laenge)
!
!  Replaces a string "value(expression)" by the value of the expression
!
USE do_eval_mod
USE errlist_mod
USE search_string_mod
!
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: laenge
!
CHARACTER(LEN=1024) :: zeile, string
INTEGER             :: iv
INTEGER             :: lll
INTEGER             :: ikl
!
main: DO
   iv = INDEX(line, 'value(', .TRUE.)  ! Search last occurence of "value("
   IF(iv<=0) EXIT main
   zeile  = line(iv+6:laenge)          ! String starting after 'value('
   lll    = laenge - iv - 5
   ikl = suche_nach(zeile, lll)        ! ikl is last character prior to ')'
   zeile=zeile(1:ikl  )
   lll  = ikl
   CALL do_eval(zeile, lll, .FALSE.)
   IF(ier_num/=0) RETURN
   string = ' '
   string = line(1:iv-1) // zeile(1:LEN_TRIM(zeile)) // line(iv+5+ikl+2:laenge)
   line = string
   laenge = LEN_TRIM(string)
ENDDO main
!
END SUBROUTINE do_value
!
!*****7**************************************************************** 
!
END MODULE do_execute_mod
