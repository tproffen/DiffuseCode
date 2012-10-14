!                                                                       
!     Application independent blockdat for command language and         
!     fortran interpreter                                               
!                                                                       
      SUBROUTINE init_sysarrays 
!                                                                       
      IMPLICIT none 
!                                                                       
      include'charact.inc' 
      include'debug.inc' 
      include'doloop.inc' 
      include'learn.inc' 
      include'macro.inc' 
      include'param.inc' 
      include'prompt.inc' 
      include'variable.inc' 
!                                                                       
      INTEGER i 
!                                                                       
!                                                                       
!     DATA statements for the various common blocks                     
!     Commented lines with the common block definition are repeated     
!     for clarity                                                       
!                                                                       
!     COMMON /DEBUG/   DBG                                              
!                                                                       
      dbg = .false. 
!                                                                       
!     COMMON /DOLOOP/  ILOOP,NLOOP,GLOW,GHIGH,GINC,KPARA                
!                                                                       
!                                                                       
!     COMMON /LEARN/   LLEARN                                           
!                                                                       
      llearn = .false. 
!                                                                       
!------ COMMON /MAKRO/     LMAKRO,MAC_LEVEL,MAC_LINE,MAC_NAME           
!                                                                       
      lmakro = .false. 
      mac_level = 0 
!                                                                       
      DO i = 1, MAC_MAX_LEVEL 
      mac_line (i) = 0 
      mac_name (i) = ' ' 
      lmacro_dbg (i) = .false. 
      ENDDO 
!                                                                       
      DO i = 1, MAC_MAX_IO 
      io_unit (i) = 86 + i 
      io_open (i) = .false. 
      io_eof (i) = .false. 
      io_get_sub (i, 1) = 1 
      io_get_sub (i, 2) = - 1 
      ENDDO 
!                                                                       
      DO i = 1, MAC_MAX_FORM 
      io_out_format (i) = '(*)' 
      ENDDO 
!                                                                       
      lsocket = .false. 
      lremote = .false. 
      lconn = .false. 
!                                                                       
!     COMMON  /PARAMS/ INPARA,RPARA,RES_PARA                            
!                                                                       
      DO i = 0, MAXPAR 
      inpara (i) = 0 
      rpara (i) = 0. 
      ENDDO 
!                                                                       
      DO i = 0, MAXPAR_RES 
      res_para (i) = 0. 
      ENDDO 
!                                                                       
!     charact.inc                                                       
!                                                                       
      a      = iachar ('a') 
      z      = iachar ('z') 
      aa     = iachar ('A') 
      zz     = iachar ('Z') 
!                                                                       
      zero   = iachar ('0') 
      nine   = iachar ('9') 
      period = iachar ('.') 
      u      = iachar ('_') 
      blank1 = iachar (' ') 
      TAB    = achar  (9) 
!                                                                       
!------ prompt.inc                                                      
!                                                                       
      output_status = OUTPUT_SCREEN 
      output_status_old = OUTPUT_SCREEN 
      output_io = 6 
      first_input = .TRUE. 
!                                                                       
!     variable.inc                                                      
!                                                                       
      var_num = 0 
      DO i = 1, VAR_MAX 
      var_name (i) = ' ' 
      var_l (i) = 0 
      var_val (i) = 0.0 
      ENDDO 
!                                                                       
!     Setup readline                                                    
!                                                                       
!      CALL cinit 
      END SUBROUTINE init_sysarrays                 
