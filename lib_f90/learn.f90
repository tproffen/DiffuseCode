!*****7**************************************************************** 
      SUBROUTINE start_learn (zeile, lcomm) 
!-                                                                      
!     These routines control if the commands typed in are               
!     within a learning sequence or not ..                              
!                                                                       
!-                                                                      
      USE errlist_mod 
      USE learn_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 12) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lcomm 
!                                                                       
      CHARACTER(1024) fname 
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ip, ianz, len_str 
      REAL werte (maxw) 
!                                                                       
      IF (llearn) then 
         ier_num = - 7 
         ier_typ = ER_IO 
      ELSE 
!                                                                       
         CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm) 
         IF (ier_num.eq.0) then 
            IF (ianz.eq.1) then 
               CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
               IF (ier_num.eq.0) then 
                  fname = cpara (1) 
               ENDIF 
            ELSEIF (ianz.eq.0) then 
               fname = pname (1:len_str (pname) ) //'.mac' 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
            IF (ier_num.eq.0) then 
               ip = index (fname, '.') 
               IF (ip.eq.0) then 
                  ip = index (fname, ' ') 
                  fname = fname (1:ip - 1) //'.mac' 
               ENDIF 
               CALL oeffne (33, fname, 'unknown') 
               IF (ier_num.ne.0) return 
               llearn = .true. 
               WRITE (output_io, 1000) fname (1:len_str (fname) ) 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
 1000 FORMAT    (' ------ > Learning started, makrofile: ',a) 
      END SUBROUTINE start_learn                    
!*****7**************************************************************** 
      SUBROUTINE ende_learn 
!+                                                                      
!     End of learning sequenze                                          
!-                                                                      
      USE errlist_mod 
      USE learn_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      IF (.not.llearn) then 
         ier_num = - 8 
         ier_typ = ER_IO 
      ELSE 
         llearn = .false. 
         CLOSE (33) 
      ENDIF 
!                                                                       
      END SUBROUTINE ende_learn                     
