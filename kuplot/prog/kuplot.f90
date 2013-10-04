      PROGRAM kuplot 
!                                                                       
      USE setup_mod
      USE kuplot_loop_mod 
      IMPLICIT none 
!*****7*****************************************************************
!       This is the universal plot program KUPLOT. It sets up most      
!     variables and calls the loop interpreting the commands.           
!*****7*****************************************************************
!                                                                       
!                                                                       
!                                                                       
      CALL setup 
      CALL kuplot_loop 
!                                                                       
!                                                                       
      END PROGRAM kuplot                            
