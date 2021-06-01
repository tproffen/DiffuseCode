module do_replace_expr_mod
!-
!   Replace any "EXPR[]" string by the appropriate character string 
!+
!
contains
!
!*******************************************************************************
!
subroutine do_replace_expr(line, ll)
!
use berechne_mod
use errlist_mod
use precision_mod
use search_string_mod
use variable_mod
!
implicit none
!
character(len=*), intent(inout) :: line
integer         , intent(inout) :: ll       ! Length of line
!
character(len=PREC_STRING) :: string
character(len=PREC_STRING) :: zeile 
integer :: length   ! Length of string
integer :: indxe    ! Last occurence of "EXPR"
integer :: ibro     ! position of "["
integer :: ibrc     ! position of "]"
integer :: kpara    ! Current index in "EXPR"
!
ll = len_trim(line)
indxe = index(line, 'EXPR', .TRUE.)
!
loop_main: do
   if(indxe==0) exit loop_main           ! no (more) "EXPR"
   ibro = indxe + index(line(indxe:), '[') - 1
   ibrc = indxe + index(line(indxe:), ']', .TRUE.) - 1  !NEEDS WORK
   zeile = line(ibro:)                   ! Copy line "[...."
   length = len_trim(zeile)
   ibrc = suche_pair(zeile, length, '[', ']') 
   if(ibrc>0) then
      ibrc = ibrc + ibro -1
   endif
   string = '(' // line(ibro+1:ibrc-1) // ')'         ! Index at [...] should have no "EXPR" in it
   length = ibrc - ibro + 1
   kpara = berechne(string, length)       ! Calculate the current index
   zeile = ' '
   if(indxe>1) then
      zeile(1:indxe-1) = line(1:indxe-1)  ! Prepend start of line
   endif
   zeile(indxe:) = var_expr(kpara)        ! Add actual string
   if(ibrc<ll) then                       ! Add remainder of line
      zeile(len_trim(zeile)+1:) = line(ibrc+1:)
   endif
   ll = len_trim(zeile)
   line = zeile
   indxe = index(line, 'EXPR', .TRUE.)
enddo loop_main
!
end subroutine do_replace_expr
!
!*******************************************************************************
!
end module do_replace_expr_mod
