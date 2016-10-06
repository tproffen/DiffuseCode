!                                                                       
!     Application independent blockdat for command language and         
!     fortran interpreter                                               
!                                                                       
      SUBROUTINE init_sysarrays 
!                                                                       
      USE charact_mod
      USE debug_mod 
      USE doloop_mod
      USE learn_mod 
      USE macro_mod 
      USE param_mod 
      USE prompt_mod 
      USE variable_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER i 
!                                                                       
!                                                                       
!     DATA statements for the various modules 
!     Commented lines with the common block definition are repeated     
!     for clarity                                                       
!                                                                       
!     /DEBUG/   DBG                                              
!                                                                       
      dbg = .false. 
!                                                                       
!     /DOLOOP/  ILOOP,NLOOP,GLOW,GHIGH,GINC,KPARA                
!                                                                       
!                                                                       
!     /LEARN/   LLEARN                                           
!                                                                       
      llearn = .false. 
!                                                                       
!------ /MAKRO/     LMAKRO,MAC_LEVEL,MAC_LINE,MAC_NAME           
!                                                                       
!     lmakro = .false. 
!     mac_level = 0 
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
!     /PARAMS/ INPARA,RPARA,RES_PARA                            
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
!------ prompt_mod                                                      
!                                                                       
      output_status = OUTPUT_SCREEN 
      output_status_old = OUTPUT_SCREEN 
      output_io = 6 
      first_input = .TRUE. 
!                                                                       
!     variable.inc                                                      
!                                                                       
      CALL variable_init
!     var_num = 0 
!     DO i = 1, VAR_MAX 
!     var_name = ' ' 
!     var_l    = 0 
!     var_val  = 0.0 
!     ENDDO 
!                                                                       
!     Setup readline                                                    
!                                                                       
!      CALL cinit 
      END SUBROUTINE init_sysarrays                 
