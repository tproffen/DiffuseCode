!*****7***************************************************************  
      SUBROUTINE init_devices 
!-                                                                      
!       Platform specific defaults                                      
!       X11 version ..                                                  
!+                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(80) prn_name 
      INTEGER ll 
      INTEGER len_str 
!                                                                       
      CALL get_environment_variable ('KUPLOT_PRINTER', prn_name) 
      ll = len_str (prn_name) 
!                                                                       
!------ X11 device definitions for PGPLOT                               
!                                                                       
      dev_name (x11) = '/XSERV' 
      dev_prn (x11) = ' ' 
      dev_name (ps) = '/CPS' 
      dev_name (vps) = '/VCPS' 
      IF (ll.eq.0) then 
         dev_prn (ps) = 'lpr ' 
         dev_prn (vps) = 'lpr ' 
      ELSE 
         dev_prn (ps) = 'lpr -P'//prn_name (1:len_str (prn_name) ) 
         dev_prn (vps) = 'lpr -P'//prn_name (1:len_str (prn_name) ) 
      ENDIF 
      dev_name (pic) = '/GIF' 
      dev_prn (pic) = ' ' 
      dev_name (vpic) = '/VGIF' 
      dev_prn (vpic) = ' ' 
      dev_name (png) = '/PNG' 
      dev_prn (png) = ' ' 
!                                                                       
!------ PGPLOT returns different characters for                         
!------ mouse clicks ;-(                                                
!                                                                       
      butt_l = 'A' 
      butt_m = 'D' 
      butt_r = 'X' 
!                                                                       
      END SUBROUTINE init_devices                   
