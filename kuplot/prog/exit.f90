!*****7**************************************************************** 
      SUBROUTINE do_exit 
!                                                                       
!       Things to do when KUPLOT exits                                  
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      USE config_mod 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
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
