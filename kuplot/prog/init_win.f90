!*****7***************************************************************  
      SUBROUTINE init_devices 
!-                                                                      
!       Platform specific defaults                                      
!       WINDOWS version ..                                              
!+                                                                      
      USE kuplot_config 
      USE kuplot_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER(80) prn_name 
      INTEGER id, irtn, len_str 
!                                                                       
!------ WINDOWS device definitions for PGPLOT                           
!                                                                       
      CALL getenv ('KUPLOT_PRINTER', prn_name) 
      IF (prn_name.eq.' ') prn_name = 'LPT1' 
!                                                                       
      dev_name (x11) = '/GW' 
      dev_prn (x11) = ' ' 
      dev_name (ps) = '/CPS' 
      dev_prn (ps) = 'print /d:'//prn_name (1:len_str (prn_name) ) 
      dev_name (vps) = '/VCPS' 
      dev_prn (vps) = 'print /d:'//prn_name (1:len_str (prn_name) ) 
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
      butt_m = 'X' 
      butt_r = 'D' 
      keyb_l = 'l' 
      keyb_m = 'm' 
      keyb_r = 'r'
!                                                                       
      END SUBROUTINE init_devices                   
