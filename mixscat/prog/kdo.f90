!*****7*****************************************************************
      SUBROUTINE mixscat_mache_kdo (line, lend, length) 
!+                                                                      
!     This routine interprets the commands and executes the             
!     corresponding function.                                           
!-                                                                      
      USE blanks_mod
      USE calc_expr_mod
      USE charact_mod
      USE errlist_mod 
      USE kdo_all_mod 
      USE learn_mod 
!                                                                       
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN= * ), INTENT(INOUT) :: line 
      LOGICAL            , INTENT(OUT)   :: lend
      INTEGER            , INTENT(INOUT) :: length 
!                                                                       
      CHARACTER(1024) zei 
      CHARACTER(70) command 
      CHARACTER(4) bef 
      INTEGER i, lc, lbef
      INTEGER ianz, indxg, indxb, indxt
!                                                                       
      LOGICAL str_comp 
!                                                                       
      SAVE command 
!                                                                       
      CALL no_error 
!                                                                       
!------ If a commentary return immediately                              
!                                                                       
      IF (line (1:1) .eq.' '.or.line (1:1) .eq.'#'.or.line (1:1) .eq.'!'&
     &.or.length.eq.0) return                                           
!                                                                       
!     Only the first 4 characters are significant. The command consists 
!     of the four nonblank characters                                   
!                                                                       
      bef = '    ' 

      indxt = INDEX (line, tab)       ! find a tabulator
      IF(indxt==0) indxt = length + 1
      indxb = index (line, ' ')       ! find a blank
      IF(indxb==0) indxb = length + 1
      indxb = MIN(indxb,indxt)
      lbef = min (indxb - 1, 4) 
      bef = line (1:lbef) 
!                                                                       
! command parameters start at the first character following the blank   
!                                                                       
      zei = ' ' 
      lc = 0 
      IF (indxb + 1.le.length) then 
         zei = line (indxb + 1:length) 
         lc = length - indxb 
         call rem_leading_bl(zei,lc)
      ENDIF 
!                                                                       
!-------Suche nach einem "="                                            
!                                                                       
      indxg = index (line, '=') 
      IF (indxg.ne.0.and..not. (str_comp (bef, 'echo', 2, lbef, 4) ) .an&
     &d..not. (str_comp (bef, 'syst', 2, lbef, 4) ) .and..not. (str_comp&
     & (bef, 'par', 2, lbef, 3) ) .and..not. (str_comp (bef, 'help', 2, &
     &lbef, 4) .or.str_comp (bef, '?   ', 1, lbef, 4) ) ) then          
!                                                                       
!-------Zuweisung eines Funktionswertes                                 
!                                                                       
         CALL do_math (line, indxg, length) 
      ELSE 
!                                                                       
!     --execute a macro file                                            
!                                                                       
         IF (bef (1:1) .eq.'@') then 
            IF (length.ge.2) then 
               CALL file_kdo (line (2:length), length - 1) 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_MAC 
            ENDIF 
!                                                                       
!------ Terminate MIXSCAT 'exit'                                        
!                                                                       
         ELSEIF (str_comp (bef, 'exit', 2, lbef, 4) ) then 
            CALL mixscat_do_exit 
            lend = .true. 
!                                                                       
!-------  --------------------------------------------------------------
!-------  Here we have MIXSCAT specific commands                        
!-------  --------------------------------------------------------------
!                                                                       
!-------  Reading structure, data and PDFFIT files                      
!                                                                       
         ELSEIF (str_comp (bef, 'read', 3, lbef, 4) ) then 
            CALL do_read (zei, lc) 
!                                                                       
!-------  Setting sample composition                                    
!                                                                       
         ELSEIF (str_comp (bef, 'elements', 2, lbef, 8) ) then 
            CALL do_elements (zei, lc) 
!                                                                       
!-------  Setting Q for f(xray)                                         
!                                                                       
         ELSEIF (str_comp (bef, 'xray', 2, lbef, 4) ) then 
            CALL do_xray (zei, lc) 
!                                                                       
!-------  Reset (only data set counter)                                 
!                                                                       
         ELSEIF (str_comp (bef, 'reset', 3, lbef, 5) ) then 
            exp_nd = 0 
!                                                                       
!-------  Setting atom pair to be removed                               
!                                                                       
         ELSEIF (str_comp (bef, 'remove', 3, lbef, 6) ) then 
            CALL do_remove (zei, lc) 
!                                                                       
!-------  Performing the calculation                                    
!                                                                       
         ELSEIF (str_comp (bef, 'calc', 2, lbef, 4) ) then 
            CALL do_calc 
!                                                                       
!-------  Matching rho0 slopes for auto scaling                         
!                                                                       
         ELSEIF (str_comp (bef, 'match', 2, lbef, 5) ) then 
            CALL do_match (zei, lc) 
!                                                                       
!-------  Setting scale factor for data set                             
!                                                                       
         ELSEIF (str_comp (bef, 'scale', 4, lbef, 5) ) then 
            CALL do_scale (zei, lc) 
!                                                                       
!-------  Setting scattering lengths manually                           
!                                                                       
         ELSEIF (str_comp (bef, 'scat', 4, lbef, 4) ) then 
            CALL do_scat (zei, lc) 
!                                                                       
!-------  Display various settings                                      
!                                                                       
         ELSEIF (str_comp (bef, 'show', 2, lbef, 4) ) then 
            CALL mixscat_do_show (zei, lc) 
!                                                                       
!-------  Save various files                                            
!                                                                       
         ELSEIF (str_comp (bef, 'save', 2, lbef, 4) ) then 
            CALL do_save (zei, lc) 
!                                                                       
!-------  Check for general commands                                    
!                                                                       
         ELSE 
            CALL kdo_all (bef, lbef, zei, lc) 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE mixscat_mache_kdo                      
