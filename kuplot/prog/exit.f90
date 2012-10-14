!*****7**************************************************************** 
      SUBROUTINE do_exit 
!                                                                       
!       Things to do when KUPLOT exits                                  
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'prompt.inc' 
      include'kuplot.inc' 
      include'errlist.inc' 
!                                                                       
!------ close PGPLOT devices                                            
!                                                                       
      CALL PGEND 
!                                                                       
!------ call system wide exit routine                                   
!                                                                       
      CALL exit_all 
!                                                                       
      END SUBROUTINE do_exit                        
