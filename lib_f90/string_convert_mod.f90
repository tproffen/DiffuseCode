module string_convert_mod
!
contains
!
!*****7*****************************************************************
!
subroutine do_cap (str) 
!+                                                                      
!       Converts string 'str' to upper letters                          
!-                                                                      
implicit none 
!                                                                       
character (len=* ), intent(inout) :: str 
integer  :: i 
integer  :: len 
!                                                                       
do i = 1, len (str) 
   if (iachar (str (i:i) ) >=  iachar ('a') .and. &
       iachar (str (i:i) ) <=  iachar ('z')       ) then                                            
      str(i:i) = achar(iachar(str(i:i)) - iachar('a') + iachar('A'))
   endif 
enddo 
!
end subroutine do_cap                         
!
!*****7*****************************************************************
!
subroutine do_low (str) 
!+                                                                      
!       Converts string 'str' to lower case letters                          
!-                                                                      
implicit none 
!                                                                       
character (len=* ), intent(inout) :: str 
integer  :: i 
integer  :: len 
!                                                                       
do i = 1, len (str) 
   if (iachar (str (i:i) ) >=  iachar ('A') .and. &
       iachar (str (i:i) ) <=  iachar ('Z')       ) then                                            
      str(i:i) = achar(iachar(str(i:i)) - iachar('A') + iachar('a'))
   endif 
enddo 
!
end subroutine do_low
!
!*****7*****************************************************************
!
subroutine do_str(str) 
!+                                                                      
!       Converts string 'str' to pure letters, blanks are not removed
!-                                                                      
implicit none 
!                                                                       
character (len=* ), intent(inout) :: str 
integer  :: i 
integer  :: len 
!                                                                       
main: do i = 1, len (str) 
   if ((iachar (str (i:i) ) >=  iachar ('a') .and.                              &
        iachar (str (i:i) ) <=  iachar ('z')      ) .or.                        &
       (iachar (str (i:i) ) >=  iachar ('A') .and.                              &
        iachar (str (i:i) ) <=  iachar ('Z')      )       ) then                                            
      cycle main
   else
      str(i:i) = ' '
   endif 
enddo  main
!
end subroutine do_str                         
!
!*****7***********************************************************      
!
end module string_convert_mod
