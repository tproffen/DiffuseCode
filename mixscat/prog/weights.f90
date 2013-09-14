!*****7*****************************************************************
      SUBROUTINE dlink (is) 
!-                                                                      
!     This routine dlinks between the MIXSCAT package and the routines  
!     that generate the atomic form factors. These routines have been   
!     adapted from the program LAZY.                                    
!+                                                                      
      USE errlist_mod 
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      COMMON / anomal / delfr, delfi 
      COMMON / neu / fneu 
      COMMON / scat / fa, fb, fc 
!                                                                       
      CHARACTER(4) line 
      INTEGER element, symwl, kodlp, i, j, is 
      LOGICAL lxray 
      REAL fa (4, 8), fb (4, 8), fc (8) 
      REAL delfr, delfi, fneu 
      REAL wave 
!                                                                       
      ier_num = 0 
      ier_typ = ER_NONE 
!                                                                       
      WRITE (line, 100) '    ' 
      READ (line, 100) symwl 
      IF (.not.exp_lxray (is) ) then 
         kodlp = 2 
      ELSE 
         kodlp = 1 
      ENDIF 
      IF (elem_n.gt.0) then 
         DO i = 1, elem_n 
         IF (pdf_scat_int (i, is) ) then 
            IF (.not.exp_lxray (is) ) then 
               WRITE (line, 100) elem_name (i) 
      line (3:4)  = '  ' 
               READ (line, 100) element 
               CALL anoma (symwl, element, kodlp) 
               pdf_scat (1, i, is) = fneu 
               DO j = 2, 9 
               pdf_scat (j, i, is) = 0.0 
               ENDDO 
            ELSE 
               WRITE (line, 100) elem_name (i) 
               READ (line, 100) element 
               CALL FSC (element, 1) 
               pdf_scat (1, i, is) = fc (1) 
               DO j = 1, 4 
               pdf_scat (2 * j, i, is) = fa (j, 1) 
               pdf_scat (2 * j + 1, i, is) = fb (j, 1) 
               ENDDO 
            ENDIF 
!                                                                       
            IF (element.eq.0) then 
               ier_num = - 20 
               ier_typ = ER_APPL 
               GOTO 999 
            ENDIF 
         ENDIF 
         ENDDO 
      ELSE 
         ier_num = - 21 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
  999 CONTINUE 
!                                                                       
  100 FORMAT      (a4) 
  101 FORMAT      (a2) 
      END SUBROUTINE dlink                          
!*****7*****************************************************************
      REAL function quad (h, k, rten) 
!+                                                                      
!           Calculates the scalar product of h and k.                   
!           1/d**2 = h(i)*k(j)*rten(i,j)                                
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER idim 
      PARAMETER (idim = 3) 
!                                                                       
      INTEGER i, j 
      REAL h (idim), k (idim), rten (idim, idim) 
!                                                                       
      quad = 0.0 
      DO i = 1, idim 
      DO j = 1, idim 
      quad = quad+h (i) * k (j) * rten (i, j) 
      ENDDO 
      ENDDO 
      END FUNCTION quad                             
!*****7*****************************************************************
      REAL function form (ll, h2, is) 
!+                                                                      
!       calculates the form factor                                      
!-                                                                      
      USE config_mod 
      USE mixscat_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER ll, i, is 
      REAL h2 
!                                                                       
      form = pdf_scat (1, ll, is) 
      IF (exp_lxray (is) ) then 
         DO i = 1, 4 
         form = form + pdf_scat (2 * i, ll, is) * exp ( - pdf_scat (2 * &
         i + 1, ll, is) * h2)                                           
         ENDDO 
      ENDIF 
!                                                                       
      END FUNCTION form                             
