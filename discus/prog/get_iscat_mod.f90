MODULE get_iscat_mod
!
PRIVATE
PUBLIC get_iscat
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE get_iscat (ianz, cpara, lpara, werte, maxw, lnew) 
!-                                                                      
!     Determines the scattering type of the parameter                   
!+                                                                      
      USE discus_config_mod 
      USE charact_mod
      USE crystal_mod 
      USE berechne_mod
      USE do_variable_mod
      USE errlist_mod 
USE precision_mod
      USE string_convert_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER ( * ) cpara (maxw) 
      CHARACTER(1024) zeile 
      INTEGER lpara (maxw) 
      INTEGER i, j, l, ianz, jj, jp 
      LOGICAL lnew 
      REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
!                                                                       
!                                                                       
!     ----Select which atoms are included in the wave                   
!                                                                       
!
!     Attempt to replace a string from a variable
!
      DO i=1, ianz
         zeile = ' '
         zeile = cpara(i)(1:lpara(i))
         j     = lpara(1)
         CALL ersetz_variable (zeile, j)
         jj = LEN_TRIM(zeile)
         IF(ier_num == 0 .AND. jj>2 .AND. zeile(1:1)=='''' .AND. &
            zeile(jj:jj) == ''''                                ) THEN
            cpara(i) = ' '
            cpara(i)(1:jj-2) = zeile(2:jj-1)
         ENDIF
      ENDDO
!
      ier_num = 0 
      ier_typ = ER_NONE 
      DO i = 1, maxw 
      werte (i) = 0.0 
      ENDDO 
      jj = 1 
      jp = 0 
      j = 1 
      DO while (j.le.ianz.and.ier_num.eq.0) 
      i = ichar (cpara (j) (1:1) ) 
      IF (cpara (j) .eq.'all') then 
         werte (1) = - 1 
         RETURN 
      ELSEIF ( ( (a.le.i.and.i.le.z) .or. (aa.le.i.and.i.le.zz) ) .and. &
      (index (cpara (j) , '[') .eq.0) ) then                            
         CALL do_cap (cpara (j) ) 
         ier_num = - 27 
         ier_typ = ER_APPL 
         place: DO i = 0, cr_nscat 
            IF (cpara (j)(1:lpara(j)) .eq.cr_at_lis (i) ) then 
               DO l=1,jj-1
                  IF(NINT(werte(l))== i) THEN
                     ier_num = 0 
                     ier_typ = ER_NONE 
                     CYCLE place    ! avoid double placement
                  ENDIF
               ENDDO
               werte (jj) = i 
               jj = jj + 1 
               jp = jp + 1 
               ier_num = 0 
               ier_typ = ER_NONE 
            ENDIF 
         ENDDO place
         IF (lnew) then 
            IF (j.gt.1) then 
               IF (cr_nscat.lt.MAXSCAT) then 
                  cr_nscat = cr_nscat + 1 
                  cr_at_lis (cr_nscat) = cpara (j) 
                  cr_dw (cr_nscat) = 0.0 
                  werte (jj) = cr_nscat 
                  ier_num = 0 
                  ier_typ = ER_NONE 
               ELSE 
                  ier_num = - 26 
                  ier_typ = ER_APPL 
               ENDIF 
            ENDIF 
         ENDIF 
      ELSE 
         zeile = ' ' 
         l = lpara (j) 
         zeile (1:1) = '(' 
         zeile (2:l + 1) = cpara (j) (1:lpara(j))
         zeile (l + 2:l + 2) = ')' 
         l = l + 2 
         werte (jj) = berechne (zeile, l) 
         IF (ier_num.eq.0) then 
            IF (0.le.nint (werte (jj) ) .and.nint (werte (jj) )         &
            .le.cr_nscat) then                                          
               jj = jj + 1 
               ier_num = 0 
               ier_typ = ER_NONE 
            ELSE 
               ier_num = - 27 
               ier_typ = ER_APPL 
            ENDIF 
         ENDIF 
      ENDIF 
      j = j + 1 
      ENDDO 
      ianz = max (ianz, jp) 
!                                                                       
END SUBROUTINE get_iscat                      
!
!*******************************************************************************
!
END MODULE get_iscat_mod
