PROGRAM diffev 
!                                                                       
USE setup_mod
USE diffev_loop_mod
!                                                                       
IMPLICIT none 
!
!*****7*****************************************************************
!                                                                       
!     Main program for DIFFEV                                           
!                                                                       
!     This is the main program for DIFFEV. It sets up most              
!     variables and calls the loop interpreting the commands.           
!                                                                       
!     Authors : R.B. Neder  (reinhard.neder@fau.de)      
!                                                                       
!*****7*****************************************************************
!
CALL setup
CALL diffev_set_sub
CALL diffev_loop
!                                                                       
!                                                                       
END PROGRAM diffev                            
