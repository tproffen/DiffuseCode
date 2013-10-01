      SUBROUTINE ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz) 
!-                                                                      
!       Replaces a substring in an expression by the value of the       
!       appropriate parameter.                                          
!+                                                                      
      USE errlist_mod 
      USE param_mod 
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER,                    INTENT(IN   ) :: ikl
      INTEGER,                    INTENT(IN   ) :: iklz
      CHARACTER (LEN=*),          INTENT(OUT  ) :: string 
      INTEGER,                    INTENT(OUT  ) :: ll
      INTEGER,                    INTENT(IN   ) :: maxw
      REAL   , DIMENSION(1:maxw), INTENT(IN   ) :: ww
      INTEGER,                    INTENT(IN   ) :: ianz
!
      CHARACTER(1024) zeile 
!                                                                       
      INTEGER laenge, ltyp, kpara, kpara2, kpara3 
      INTEGER lcomm 
      INTEGER length_com 
      INTEGER i1, i2 
!                                                                       
      laenge = ll 
      ltyp = 1 
      zeile = ' ' 
      kpara = nint (ww (1) ) 
!                                                                       
      IF (maxw.ge.2) then 
         kpara2 = nint (ww (2) ) 
      ENDIF 
!                                                                       
      IF (maxw.ge.3) then 
         kpara3 = nint (ww (3) ) 
      ENDIF 
!                                                                       
      lcomm = length_com (string, laenge, ikl) 
!                                                                       
!------ Variables of length 1                                           
!                                                                       
      IF (lcomm.eq.1) then 
         IF (ikl.gt.2) zeile (1:ikl - 1) = string (1:ikl - 2) 
!                                                                       
!------ - General integer var i[n]                                      
!                                                                       
         IF (string (ikl - 1:ikl - 1) .eq.'i') then 
            IF (ianz.eq.1) then 
               IF (0.le.kpara.and.kpara.le.MAXPAR) then 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(i15)') inpara (   &
                  kpara)                                                
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
!                                                                       
!------ - General real var r[n]                                         
!                                                                       
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'r') then 
            IF (ianz.eq.1) then 
               IF (0.le.kpara.and.kpara.le.MAXPAR) then 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(e15.8e2)') rpara (&
                  kpara)                                                
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
!                                                                       
!-------- Program information n[n]                                      
!                                                                       
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'n') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara.and.kpara.le.1) then 
                  IF (kpara.eq.1) then 
                     WRITE (zeile (ikl - 1:ikl + 13) , '(i15)') exp_nd 
                  ENDIF 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
         ENDIF 
!                                                                       
!------ Variables of length 3                                           
!                                                                       
      ELSEIF (lcomm.eq.3) then 
         IF (ikl.gt.4) zeile (1:ikl - 3) = string (1:ikl - 4) 
!                                                                       
!-------- result array res[n]                                           
!                                                                       
         IF (string (ikl - 3:ikl - 1) .eq.'res') then 
            IF (ikl.gt.4) zeile (1:ikl - 3) = string (1:ikl - 4) 
            IF (ianz.eq.1) then 
               IF (0.le.kpara.and.kpara.le.MAXPAR_RES) then 
                  WRITE (zeile (ikl - 3:ikl + 13) , '(e15.8e2)')        &
                  res_para (kpara)                                      
                  zeile (ikl + 8:ikl + 8) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 2 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 2 
         ier_typ = ER_FORT 
      ENDIF 
!                                                                       
      IF (ier_num.eq.0) then 
         ll = laenge+15 - ltyp - (iklz - ikl + 1) 
         IF (iklz + 1.le.laenge) zeile (ikl + 14:ll) = string (iklz + 1:&
         laenge)                                                        
         string = zeile 
      ELSE 
         ier_msg (1) = string 
      ENDIF 
      CALL rem_bl (string, ll) 
      END SUBROUTINE ersetz_para                    
!*****7*****************************************************************
      SUBROUTINE upd_para (ctype, ww, maxw, wert, ianz) 
!-                                                                      
!       updates the parameter spezified by ctype, index ww  to the      
!       new value of wert                                               
!+                                                                      
      USE errlist_mod 
      USE param_mod 
      USE config_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN=*),          INTENT(IN) :: ctype 
      INTEGER,                    INTENT(IN) :: maxw
      INTEGER,                    INTENT(IN) :: ianz 
      INTEGER, DIMENSION(1:MAXW), INTENT(IN) :: ww
      REAL   ,                    INTENT(IN) :: wert 
!                                                                       
!------ Setting i[n]                                                    
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
!                                                                       
!------ Setting r[n]                                                    
!                                                                       
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
!                                                                       
      ELSE 
         ier_num = - 2 
         ier_typ = ER_FORT 
         ier_msg (1) = 'Variable '//ctype 
      ENDIF 
!                                                                       
      END SUBROUTINE upd_para                       
!*****7***************************************************************  
      SUBROUTINE calc_intr_spec (string, line, ikl, iklz, ww, laenge,   &
      lp)                                                               
!-                                                                      
!     These are special intrinsic function for the PDFFIT. Any          
!     intrinsic function that references crystallographic values        
!     is found in this subroutine.                                      
!+                                                                      
      USE errlist_mod 
      USE param_mod 
      USE config_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN=*), INTENT(IN   ) :: string
      CHARACTER (LEN=*), INTENT(INOUT) :: line 
      INTEGER,           INTENT(IN)    :: ikl
      INTEGER,           INTENT(IN)    :: iklz
      INTEGER,           INTENT(IN)    :: laenge
      INTEGER,           INTENT(IN)    :: lp
      REAL   ,           INTENT(OUT)   :: ww
      END SUBROUTINE calc_intr_spec                 
!*****7**************************************************************** 
      SUBROUTINE validate_var_spec (zeile, lp) 
!-                                                                      
!       checks whether the variable name is legal, MIXSCAT specific     
!                                                                       
!+                                                                      
      USE errlist_mod 
      USE config_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN=*), INTENT(IN) :: zeile 
      INTEGER,           INTENT(IN) :: lp 
!                                                                       
      INTEGER, PARAMETER :: reserved_n = 52 
                                                                        
      CHARACTER(12) reserved (reserved_n) 
      INTEGER i 
!                                                                       
      LOGICAL str_comp 
      DATA reserved / 'bang', 'blen', 'dstar', 'rang', 'gran', 'gbox',  &
      'rval', 'x', 'y', 'z', 'm', 'u', 'o', 'lat', 'delt', 'rcut',      &
      'srat', 'csca', 'rhoz', 'qsig', 'dsca', 'bave', 'dx', 'dy', 'dz', &
      'du', 'do', 'dlat', 'ddelt', 'dsrat', 'dcsca', 'drhoz', 'ddsca',  &
      'dqsig', 'n', 'p', 'dp', 'pf', 'cl', 'rw', 'np', 'pc', 'po', 'pw',&
      'delr', 'qmax', 'rang', 'rmin', 'rmax', 'fa', 'fb', 'fc' /        
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
      END SUBROUTINE validate_var_spec              
