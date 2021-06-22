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
!-
!  Replace the EXPR[] string by the actual character string
!+
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
subroutine do_evaluate_expr
!-
!  Evaluates the current EXPR'essions
!  Values are stored in var_expr_val(index)
!  If EXPR reference each other cyclically an error is flagged
!+
!
use berechne_mod
use errlist_mod
use precision_mod
use prompt_mod
use search_string_mod
use variable_mod
!
implicit none
!
character(len=PREC_STRING) :: line     ! dummy string
character(len=PREC_STRING) :: string   ! dummy string
character(len=PREC_STRING) :: zeile    ! dummy string
integer :: i
integer :: length
integer :: n_expr   ! Number of EXPR[] currently defined
integer :: indxe    ! Last occurence of "EXPR"
integer :: ibro     ! position of "["
integer :: ibrc     ! position of "]"
integer :: kpara    ! Current index in "EXPR"
integer :: nchange  ! Number of EXPR changed in current loop
logical, dimension(ubound(var_expr,1)) :: done
real(kind=PREC_DP) :: wert
!
n_expr = ubound(var_expr,1)    ! Determine number of expressions
done   = .false.               ! No expresions have been done
nchange = 0
!
loop_done: do                  ! Grand loop to allow arbitrary sequence of nested EXPR
   nchange = 0
   loop_repl: do i=1, n_expr   ! Loop over all EXPR
      if(done(i)) cycle loop_repl              ! This one is done, skip
!
      line = var_expr(i)
      indxe = index(line, 'EXPR', .TRUE.)
!
!     Search for EXPR embedded in current EXPR string
!
      loop_expr: do
         if(indxe==0) exit loop_expr           ! no (more) "EXPR"
         ibro = indxe + index(line(indxe:), '[') - 1
         zeile = line(ibro:)                   ! Copy line "[...."
         length = len_trim(zeile)
         ibrc = suche_pair(zeile, length, '[', ']') 
         if(ibrc>0) then
            ibrc = ibrc + ibro -1
         endif
         string = '(' // line(ibro+1:ibrc-1) // ')'         ! Index at [...] should have no "EXPR" in it
         length = ibrc - ibro + 1
         kpara = berechne(string, length)       ! Calculate the current index
         if(done(kpara)) then                   ! EXPR has been done, OK
            zeile = ' '
            if(ibro>5) zeile = line(1:ibro-5)      
            length = len_trim(zeile)
            write(zeile(length+1:length+16),  '(E16.8E3)') var_expr_val(kpara)
            length = len_trim(zeile)
            if(ibrc<len_trim(line)) zeile(length+1:) = line(ibrc+1:)
            line = zeile
            length = len_trim(line)
         else                                   ! Not yet done, skip current EXPR[i] till later
            cycle loop_repl
         endif
         indxe = index(line, 'EXPR', .TRUE.)    ! Any further EXPR within current?
      enddo loop_expr
      length = len_trim(line)
      line = '(' // line(1:length) // ')'
      length = length + 2
      if(line=='()') then              ! Empty EXPR, ignore
         wert = 0.0D0
      else
         wert = berechne(line, length) ! Evaluate EXPR at this moment
      endif
      var_expr_val(i) = wert           ! Place evaluated number into storage
      if(ier_num/=0) then              ! Error in berechne
         ier_num = -54
         ier_num = ER_FORT
         write(ier_msg(1),'(a,i4)') 'EXPR ',i
         return
      endif
      done(i) = .true.                 ! This expression is done
      nchange = nchange + 1            ! Increment success flag
   enddo loop_repl
!
   if(all(done)) then                  ! All EXPR have been replaced by numerical value
      exit loop_done
   else                                ! Not yet finished
      if(nchange == 0) then            ! No changes; error
         ier_num = -53
         ier_typ = ER_FORT
         ier_msg(1) = 'See list of expressions listed above'
         exit loop_done
      endif
   endif
enddo loop_done
!
do i=1,n_expr
   if(.not.done(i)) then
      write(output_io,'(a,i3,a,a)') ' ***FORT*** EXPR[',i,'] = ',var_expr(i)(1:len_trim(var_expr(i)))
!  else
!     write(output_io,'(a,i3,a,a, a, f10.8)') ' ***FORT*** EXPR[',i,'] = ',var_expr(i)(1:len_trim(var_expr(i))), &
!     ' = ',var_expr_val(i)
   endif
enddo
!
end subroutine do_evaluate_expr
!
!*******************************************************************************
!
subroutine do_use_expr(ianz, cpara, lpara, maxpara)
!-
!   use the evaluated expressions and replace the EXPR[] by this value
!+
use errlist_mod
use berechne_mod
use precision_mod
use search_string_mod
use variable_mod
!
implicit none
!
integer, intent(in) :: ianz
integer, intent(in) :: maxpara
character(LEN=*   ), dimension(MAXPARA), intent(inout)  :: cpara
integer            , dimension(MAXPARA), intent(inout)  :: lpara
!
character(len=PREC_STRING) :: line     ! dummy string
character(len=PREC_STRING) :: string   ! dummy string
character(len=PREC_STRING) :: zeile    ! dummy string
integer :: i
integer :: length
integer :: indxe    ! Last occurence of "EXPR"
integer :: ibro     ! position of "["
integer :: ibrc     ! position of "]"
integer :: kpara    ! Current index in "EXPR"
!
call do_evaluate_expr                      ! Determine current EXPR value
!
loop_main: do i=1, ianz
  line   = cpara(i)
  length = lpara(i)
  indxe  = index(line, 'EXPR', .TRUE.)
!
!  Search for EXPR with lower index
!
   loop_expr: do
      if(indxe==0) exit loop_expr           ! no (more) "EXPR"
      ibro = indxe + index(line(indxe:), '[') - 1
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
      if(ibro>5) zeile = line(1:ibro-5)
      length = len_trim(zeile)
      write(zeile(length+1:length+16),  '(E16.8E3)') var_expr_val(kpara)
      length = len_trim(zeile)
      if(ibrc<len_trim(line)) zeile(length+1:) = line(ibrc+1:)
      line = zeile
      length = len_trim(line)
      indxe  = index(line, 'EXPR', .TRUE.)
   enddo loop_expr
   length = len_trim(line)
   cpara(i) = line
   lpara(i) = length
enddo loop_main
!
end subroutine do_use_expr
!
!*******************************************************************************
!
end module do_replace_expr_mod
