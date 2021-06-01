module do_string_alloc_mod
!
contains
!
!****7***************************************************************** 
!
subroutine do_string_alloc (line, indxg, length) 
!-                                                                      
!     Evaluates the parameters and stores the result                    
!     in the proper string variable                                     
!                                                                       
use ber_params_mod
use do_variable_mod
use build_name_mod
!     use calc_intr_mod
use errlist_mod 
use get_params_mod
use param_mod 
use set_sub_generic_mod
use precision_mod
use variable_mod
!
implicit none 
!                                                                       
integer, parameter :: maxw= 10 
!
character (len= * ), intent(INOUT) :: line 
integer            , intent(IN)    :: indxg
integer            , intent(INOUT) :: length
!
character(len=max(PREC_STRING,len(line))) :: zeile
character(len=max(PREC_STRING,len(line))), dimension(maxw) :: cpara
character(len=max(PREC_STRING,len(line)))                  :: string
!                                                                       
integer, dimension(maxw)   :: lpara (maxw) 
integer                    :: i, ikk, ianz, lll 
integer, dimension(1:maxw) :: iii = 0
integer                    :: ising , indx_ind, indx_len, indx_env, indx_cwd
integer                    :: indx_par, indx_isv, indx_ise, indx_fmd
integer                    :: l_string 
integer                    :: ikl, iklz, ll, laenge
integer                    :: ios  ! I/O error status
integer                    :: indxpl ! Location of a substring '('
integer                    :: indxpr ! Location of a substring ')'
integer                    :: indxc  ! Location of a substring ':'
integer, dimension(2)      :: substr ! Indices of substring
!                                                                       
real(kind=PREC_DP)    :: wert
real(kind=PREC_DP), dimension(MAXW) :: werte
!                                                                       
!     for flexibility                                                   
!                                                                       
ising = index (line, '''') 
do while (ising.gt.0) 
   line (ising:ising) = '"' 
   ising = index (line, '''') 
enddo 
!                                                                       
if(line(1:4)=='EXPR') then
  call do_set_expression (line, indxg, length) 
  return
endif
!     Get the expression                                                
!                                                                       
lll = - (length - (indxg + 1) + 1) 
call get_params (line (indxg + 1:length), ianz, cpara, lpara, maxw, lll)
if (ier_num.ne.0) then 
   return 
elseif (ianz.eq.0) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   return 
endif 
!     locate first ""
ising=index(line (indxg + 1:length),'"')
indx_ind=index(line (indxg + 1:length),'index')   ! locate index function
indx_len=index(line (indxg + 1:length),'length')  ! locate length function
indx_env=index(line (indxg + 1:length),'getenv')  ! locate length function
indx_cwd=index(line (indxg + 1:length),'getcwd')  ! locate length function
indx_isv=index(line (indxg + 1:length),'isvar')   ! locate isvar function
indx_ise=index(line (indxg + 1:length),'isexp')   ! locate isexp function
indx_fmd=index(line (indxg + 1:length),'fmodt')   ! locate fmodt function
indx_par=index(line (indxg + 1:length),'par_name')  ! locate par_name function
!
if((indx_ind>0 .AND. indx_ind<ising)  .OR. & ! We got a function of a string argument
   (indx_len>0 .AND. indx_len<ising)  .OR. & ! We got a function of a string argument
   (indx_env>0 .AND. indx_env<ising)  .OR. & ! We got a function of a string argument
   (indx_cwd>0 .AND. indx_cwd<ising)  .OR. & ! We got a function of a string argument
   (indx_isv>0 .AND. indx_isv<ising)  .OR. & ! We got a function of a string argument
   (indx_ise>0 .AND. indx_ise<ising)  .OR. & ! We got a function of a string argument
   (indx_fmd>0 .AND. indx_fmd<ising)  .OR. & ! We got a function of a string argument
   (indx_par>0 .AND. indx_par<ising)) then
   string = line (indxg + 1:length)
   laenge = length - indxg
   ikl = index(string,'(')
   iklz= index(string,')',.TRUE.)
   zeile = string (ikl + 1:iklz - 1) 
   ll = iklz - ikl - 1 
   call calc_intr (string, zeile, ikl, iklz, laenge, ll) 
   READ(string(1:len_TRIM(string)),*,IOSTAT=ios) wert
   if(ios/=0) then
      ier_msg(1) = string(1:42)
      ier_num = -6
      ier_typ = ER_FORT
      return
   endif
else
!                                                                       
!     Construct the regular string, not a function
!                                                                       
   call do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
   string = cpara (1) (1:lpara (1) ) 
   l_string = lpara (1) 
   res_para (0) = 1 
   res_para (1) = lpara (1) 
endif
if(indxg==0) then
   line = string
   return
endif
if (ier_num.eq.0) then 
!
!-----evaluate a substring range ( : )
!
   substr(1) = 1     ! Always clear the substring indices
   substr(2) = VAR_CLEN  ! Always clear the substring indices
   indxpl = index(line(1:indxg-1), '(')
   if(indxpl > 0) then
      indxpr = index(line(1:indxg-1), ')', .TRUE.)  ! Find last ')'
      indxc  = index(line(1:indxg-1), ':', .TRUE.)  ! Find last ':'
      if(indxpl <= indxpr-2) then ! Enough space for at least (:)
         line(indxc:indxc) = ','
         zeile =  ' '
         if(indxc==indxpl+1 .AND. indxc==indxpr-1) then    ! == (:)
            WRITE(zeile,'(a,i4)') '1,',len(line)
         elseif(indxc==indxpl+1 .and. indxc<indxpr-1) then ! == (:number)
            zeile = '1,' // line(indxc+1:indxpr-1)
         elseif(indxc> indxpl+1 .and. indxc==indxpr-1) then ! == (number:)
            WRITE(zeile,'(a,i4)') line(indxpl+1:indxc-1) , len(line)
         else
            zeile(1:indxpr-indxpl-1) = line(indxpl+1:indxpr-1)  ! ==(number:number)
         endif
         lll = len_trim(zeile)
         call get_params (zeile, ianz, cpara, lpara, maxw, lll)
         if(ier_num.eq.0) then 
            if(ianz==2) then      ! Need exactly 2 params
               call ber_params (ianz, cpara, lpara, werte, maxw)
               if(ier_num.eq.0) then 
                  substr(1) = NINT(werte(1))
                  substr(2) = NINT(werte(2))
                  line(indxpl:indxpr) = ' '   ! Clear substring range in input line
               else
                  return
               endif 
            else 
               ier_num = -6 
               ier_typ = ER_FORT 
               return
            endif 
         else
            return
         endif
      else
         ier_num = -6
         ier_typ = ER_FORT
         ier_msg(1) = 'Substring values appear too short'
         return
      endif
   endif
!                                                                       
!-----evaluate the index of the variable                                
!                                                                       
   lll = indxg - 1 
   call get_params (line (1:indxg - 1), ianz, cpara, lpara, maxw, lll)
   if (ier_num.eq.0) then 
      line = cpara (1) 
      i = lpara (1) 
      ikk = index (line, '[') 
      if (ikk.lt.i.and.ikk.gt.0) then 
         if (line (i:i) .eq.']') then 
            if (i.gt.ikk + 1) then 
               zeile = ' ' 
               zeile (1:i - ikk - 1) = line (ikk + 1:i - 1) 
               lll = i - ikk - 1 
               call get_params (zeile, ianz, cpara, lpara, maxw, lll)
               if (ier_num.eq.0) then 
                  if (ianz.ge.1.or.ianz.le.3) then 
                     call ber_params (ianz, cpara, lpara, werte, maxw)
                     if (ier_num.eq.0) then 
                        do i = 1, ianz 
                        iii (i) = nint (werte (i) ) 
                        enddo 
!                                                                       
!     ------------Store result in the variable                          
!                                                                       
                        call p_upd_para (line (1:ikk - 1), iii, ianz, wert, ianz, string, substr)
                     endif 
                  else 
                     ier_num = - 6 
                     ier_typ = ER_FORT 
                  endif 
               endif 
            else 
               ier_num = - 10 
               ier_typ = ER_FORT 
            endif 
         else 
            ier_num = - 9 
            ier_typ = ER_FORT 
         endif 
      elseif (ikk.eq.0) then 
         call upd_variable (line (1:i), i, wert, string, l_string, substr) 
      else 
         ier_num = - 6 
         ier_typ = ER_FORT 
      endif 
   endif 
endif 
!                                                                       
end subroutine do_string_alloc                
!
!****7***************************************************************** 
!
subroutine do_set_expression(line, indxg, length) 
!-
!  set the values in the "EXPR" variables
!+
!
use ber_params_mod
use errlist_mod 
use lib_f90_allocate_mod
use precision_mod
use search_string_mod
use variable_mod
!
implicit none
!
character (len= * ), intent(INOUT) :: line 
integer            , intent(IN)    :: indxg
integer            , intent(INOUT) :: length
!
integer, parameter :: MAXW = 1
character(len=max(PREC_STRING,len(line))), dimension(MAXW) :: cpara
character(len=PREC_STRING) :: dummy
!                                                                       
integer, dimension(MAXW)   :: lpara (maxw) 
integer :: ianz
integer :: ibro    ! brace open
integer :: ibrc    ! brace close
integer :: n_expr  ! Number of expressions
integer :: ll
real(kind=PREC_DP), dimension(MAXW) :: werte
!
ibro = index(line(1:indxg-1), '[')           ! Locate square braces []
!ibrc = index(line(1:indxg-1), ']')
dummy = line(ibro:)
ll = len_trim(dummy)
ibrc = suche_pair(dummy, ll, '[', ']') + ibro - 1
if(ibro==0 .or. ibrc==0 .or. ibrc<ibro) then ! Check for consistency 
   ier_num = -9
   ier_typ = ER_FORT
   ier_msg(1) = 'Braces [] do not match at EXPR[index]'
   return
endif
ianz     = 1                                 ! Prepare and evaluate index
cpara(1) = line(ibro+1:ibrc-1)
lpara(1) = ibrc-ibro
call ber_params (ianz, cpara, lpara, werte, MAXW)
if(ier_num/=0) return
!
ianz = nint(werte(1))                        ! Override temporary variable
if(ianz > ubound(var_expr,1)) then           ! If needed allocate space
   n_expr = ubound(var_expr,1) + 5
   call alloc_expr(n_expr)
   if(ier_num/=0) return
endif
!
var_expr(ianz) = line(indxg+1:length)        ! place expression into proper entry, keep ""
!
end subroutine do_set_expression
!
!****7***************************************************************** 
!
end module do_string_alloc_mod
