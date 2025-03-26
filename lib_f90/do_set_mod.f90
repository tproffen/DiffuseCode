module do_set_mod
!
private
!
public do_set
!
!
contains
!
!*****7***********************************************************      
!
subroutine do_set(zeile, lp) 
!-                                                                      
!     Sets the value of status variables                                
!+                                                                      
use debug_mod 
use envir_mod 
use errlist_mod 
use do_wait_mod
use get_params_mod
use lib_length
use precision_mod
use prompt_mod 
use str_comp_mod
use support_mod
!                                                                       
implicit none 
!                                                                       
integer, parameter :: MAXW = 5
!                                                                       
character (len=*), intent(inout) :: zeile 
integer          , intent(inout) :: lp
!
character(len=MAX(PREC_STRING,len(zeile))), dimension(MAXW) :: cpara
integer            , dimension(MAXW) :: lpara
character(len=PREC_STRING) :: message
character(len=80)   :: logfile 
integer :: ios
integer :: ianz 
logical :: llog 
!                                                                       
!                                                                       
if (zeile.ne.' ') then 
   call get_params (zeile, ianz, cpara, lpara, maxw, lp) 
   if (ier_num.eq.0) then 
!                                                                       
!----- ---- set error                                                   
!                                                                       
      if (str_comp (cpara (1) , 'error', 1, lpara (1) , 5) ) then 
         if(ianz.eq.3) then 
            if(str_comp(cpara(3), 'save', 2, lpara(3), 4)) then
               ier_sav = ier_sta
            else 
               ier_num = - 6 
               ier_typ = ER_COMM 
               return
            endif 
            ianz = ianz -1
         endif
         if(ianz.eq.2) then 
            if(str_comp(cpara(2), 'continue', 2, lpara(2), 8)) then
               ier_sta = ER_S_CONT 
            elseif(str_comp(cpara(2), 'exit', 2, lpara (2), 4)) then
               ier_sta = ER_S_EXIT 
            elseif(str_comp(cpara(2), 'live', 2, lpara (2), 4)) then
               ier_sta = ER_S_LIVE 
            elseif(str_comp(cpara(2), 'old' , 2, lpara (2), 3)) then
               ier_sta = ier_sav
            else 
               ier_num = - 6 
               ier_typ = ER_COMM 
            endif 
         else 
            ier_num = - 6 
            ier_typ = ER_COMM 
         endif 
!                                                                       
!----- ---- set authors                                                 
!                                                                       
      elseif(str_comp(cpara(1), 'author', 3, lpara(1), 6)) then
         call do_set_author(ianz, cpara, lpara, MAXW)
!                                                                       
!----- ---- set parallel                                                
!                                                                       
      elseif(str_comp(cpara(1), 'parallel', 3, lpara(1), 8)) then
         call do_set_parallel(ianz, cpara, lpara, MAXW)
!                                                                       
!----- ---- set prompt                                                  
!                                                                       
      elseif(str_comp(cpara(1), 'prompt', 3, lpara(1), 6)) then
         if (ianz.ge.2) then 
!                                                                       
!------ ------- Third+1  optional parameter: save previous setting      
!                                                                       
            if (ianz.eq.4) then 
               prompt_status_old = prompt_status 
               output_status_old = output_status 
            endif 
!                                                                       
!------ ------- First+1 parameter PROMPT                                
!                                                                       
            if(str_comp(cpara(2), 'on', 2, lpara(2) , 2)) then
               prompt_status = PROMPT_ON 
!               socket_status = PROMPT_ON 
!                                                                       
            elseif(str_comp(cpara(2), 'off', 2, lpara(2), 3)) then
               prompt_status = PROMPT_OFF 
!               socket_status = PROMPT_OFF 
!                                                                       
            elseif(str_comp(cpara(2), 'redirect', 1, lpara(2), 8)) then
               prompt_status = PROMPT_REDIRECT 
!               socket_status = PROMPT_REDIRECT 
!                                                                       
            elseif(str_comp(cpara(2), 'old', 1, lpara(2), 3)) then
               if (output_status.ne.OUTPUT_SCREEN) then 
                  close (output_io) 
               endif 
               output_status = output_status_old 
               prompt_status = prompt_status_old 
!               socket_status = socket_status_old 
               if (output_status.eq.OUTPUT_SCREEN) then 
                  output_io = 6 
               elseif (output_status.eq.OUTPUT_none) then 
                  output_io = 37 
                  open (unit = output_io, file = nullfile, status &
                        = 'unknown',IOSTAT=ios,IOMSG=message)                                    
                  if(ios/=0) then
                     ier_num = -2
                     ier_typ = ER_IO
                     ier_msg(1) ='Could not open the NULL file'
                     ier_msg(2) = message(1:MAX(1,MIN(80,len_TRIM(message))))
                     return
                  endif
               elseif (output_status.eq.OUTPUT_FILE) then 
                  output_io = 37 
                  logfile = pname (1:len_str (pname) ) //'.log' 
                  inquire (file = logfile, exist = llog) 
                  if (llog) then 
                     call oeffne_append (output_io, logfile, 'old')
                     if (ier_num.ne.0) return 
!DBG                    open(unit=output_io,file=logfile,               
!DBG     &                   status='old',access='append')              
!DBG95     &                 status='old',access='sequential')          
                  else 
                     open (unit = output_io, file = logfile,      &
                           status = 'new',IOSTAT=ios,IOMSG=message)
                     if(ios/=0) then
                        ier_num = -2
                        ier_typ = ER_IO
                        ier_msg(1) ='Could not open the logfile'
                        ier_msg(2) = logfile
                        ier_msg(3) = message(1:MAX(1,MIN(80,len_TRIM(message))))
                        return
                     endif
                  endif 
               endif 
            else 
               ier_num = - 6 
               ier_typ = ER_COMM 
            endif 
!                                                                       
!------ ------- Second+1 optional parameter general output              
!                                                                       
            if (ianz.ge.3) then 
               if(str_comp(cpara(3), 'on', 2, lpara(3), 2)) then
                  if (output_status.ne.OUTPUT_SCREEN) then 
                     close (output_io) 
                  endif 
                  output_status = OUTPUT_SCREEN 
                  output_io = 6 
!                                                                       
               elseif(str_comp(cpara(3), 'off', 2, lpara(3), 3)) then
                  if (output_status.ne.OUTPUT_SCREEN) then 
                     close (output_io) 
                  endif 
                  output_status = OUTPUT_none 
                  output_io = 37 
                  open (unit = output_io, file = nullfile, status &
                        = 'unknown',IOSTAT=ios,IOMSG=message)                                    
                  if(ios/=0) then
                     ier_num = -2
                     ier_typ = ER_IO
                     ier_msg(1) ='Could not open the logfile'
                     ier_msg(2) = logfile
                     ier_msg(3) = message(1:MAX(1,MIN(80,len_TRIM(message))))
                     return
                  endif
!                                                                       
               elseif(str_comp(cpara(3), 'file', 2, lpara(3), 4)) then
                  if (output_status.ne.OUTPUT_SCREEN) then 
                     close (output_io) 
                  endif 
                  output_status = OUTPUT_FILE 
                  output_io = 37 
                  logfile = pname (1:len_str (pname) ) //'.log' 
                  inquire (file = logfile, exist = llog) 
                  if (llog) then 
                     call oeffne_append (output_io, logfile, 'old')
                     if (ier_num.ne.0) return 
                  else 
                     open (unit = output_io, file = logfile,      &
                           status = 'new',IOSTAT=ios,IOMSG=message)
                     if(ios/=0) then
                        ier_num = -2
                        ier_typ = ER_IO
                        ier_msg(1) ='Could not open the logfile'
                        ier_msg(2) = logfile
                        ier_msg(3) = message(1:MAX(1,MIN(80,len_TRIM(message))))
                        return
                     endif
                  endif 
               else 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               endif 
            endif 
!                                                                       
            if (dbg) then 
               write ( *, 5000) prompt_status, prompt_status_old 
               write ( *, 5010) output_status, output_status_old 
            endif 
         else 
            ier_num = - 6 
            ier_typ = ER_COMM 
         endif 
!                                                                       
!----- ---- set debug                                                   
!                                                                       
      elseif (str_comp (cpara (1) , 'debug', 2, lpara (1) , 5) ) then
         if (ianz.eq.2) then 
            dbg = (cpara (2) (1:2) .eq.'on') 
         else 
            ier_num = - 6 
            ier_typ = ER_COMM 
         endif 
!                                                                       
!----- ---- set wait                                                   
!                                                                       
      elseif (str_comp (cpara (1) , 'wait', 2, lpara (1) , 4) ) then
         if (ianz.eq.2) then 
            wait_active = (cpara (2) (1:2) .eq.'on') 
         else 
            ier_num = - 6 
            ier_typ = ER_COMM 
         endif 
!                                                                       
      else 
         ier_num = - 6 
         ier_typ = ER_COMM 
      endif 
   endif 
else 
   ier_num = - 6 
   ier_typ = ER_COMM 
endif 
!                                                                       
 5000 format    (' debug  > Prompt (current, old) : ',2I3) 
 5010 format    (' debug  > Output (current, old) : ',2I3) 
!                                                                       
end subroutine do_set                         
!
!*****7*****************************************************************
!
subroutine do_set_author(ianz, cpara, lpara, MAXW)
!
use envir_mod
use precision_mod
use take_param_mod
!
implicit none
!
integer, intent(inout) :: ianz
integer, intent(IN)    :: MAXW
character(len=*)  , dimension(MAXW), intent(inout) :: cpara
integer           , dimension(MAXW), intent(inout) :: lpara
!
integer, parameter :: NOPTIONAL = 1
integer, parameter :: O_AUTHOR  = 1                  ! Maximum difference
character(len=   7), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(len=MAX(PREC_STRING,len(cpara))), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte    ! Calculated values
integer, parameter                        :: ncalc = 0 ! Number of values to calculate
!integer :: nthreads    ! local variable number of threads
!
data oname  / 'author' /   ! 
data loname /  6       /
!
opara  = (/'DISCUS' /)
lopara = (/ 6       /)
owerte = (/ 0.0     /)
!
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                              oname, loname, opara, lopara, lpresent, owerte)
if(opara(O_AUTHOR)(1:1) == '[' .and.                                           &
   opara(O_AUTHOR)(len_trim(opara(O_AUTHOR)):len_trim(opara(O_AUTHOR)))==']' ) then
   author = opara(O_AUTHOR)(2:len_trim(opara(O_AUTHOR))-1)
else
   author = opara(O_AUTHOR)(1:len_trim(opara(O_AUTHOR)))
endif
!
end subroutine do_set_author
!
!*****7*****************************************************************
!
subroutine do_set_parallel(ianz, cpara, lpara, MAXW)
!-
! Get parameters for 'set parallel, ...' command
!+
use ber_params_mod
use errlist_mod
use parallel_mod
use param_mod
use precision_mod
use take_param_mod
!
implicit none
!
integer, intent(inout) :: ianz
integer, intent(IN)    :: MAXW
character(len=*)  , dimension(MAXW), intent(inout) :: cpara
integer           , dimension(MAXW), intent(inout) :: lpara
!
integer, parameter :: NOPTIONAL = 2
integer, parameter :: O_NTHread = 1                  ! Maximum difference
integer, parameter :: O_USEOMP  = 2                  ! Number of feedbacks to average
character(len=   7), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(len=MAX(PREC_STRING,len(cpara))), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte    ! Calculated values
integer, parameter                        :: ncalc = 0 ! Number of values to calculate
integer :: nthreads    ! local variable number of threads
!
data oname  / 'nthread', 'useomp ' /   ! 
data loname /  7       ,  6        /
!
opara  = (/'all   ', 'use   ' /)
lopara = (/ 3      ,  3       /)
owerte = (/ -1.    ,  0.0000  /)
!
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                              oname, loname, opara, lopara, lpresent, owerte)
par_omp_use        = .TRUE. ! DEFAULT: User wants to use OMP
par_omp_maxthreads = -1     ! DEFAULT: Maximum number of threads to use, -1==use all available
!
if(ier_num==0) then
   if(opara(O_USEOMP) == 'use' .or. opara(O_USEOMP) == 'parallel') then
      par_omp_use = .TRUE.
   elseif(opara(O_USEOMP) == 'serial'.or. opara(O_USEOMP) == 'off') then
      par_omp_use = .FALSE.
   else
      ier_num = -6
      ier_typ = ER_COMM
      ier_msg(1) = 'Optional OMP must be useomp:use ; useomp:parallel'
      ier_msg(2) = '                  or useomp:off ; useomp:serial'
      return
   endif
   if(par_omp_use) then
      if(opara(O_NTHread) == 'all') then
         par_omp_maxthreads = -1     ! DEFAULT: Maximum number of threads to use, -1==use all available
      elseif(opara(O_NTHread) == 'physical') then
         par_omp_maxthreads = par_omp_phys
      elseif(opara(O_NTHread) == 'logical') then
         par_omp_maxthreads = par_omp_logi
      else
         opara (O_USEOMP) = '0.0'
         lopara(O_USEOMP) = 3
         call ber_params (ianz, opara, lopara, owerte, MAXW)
         if(ier_num==0) then
            par_omp_maxthreads = NINT(owerte(O_NTHread))
            if(par_omp_maxthreads<1) then
               ier_num = -18
               ier_typ =  ER_COMM
            endif
         endif
      endif
      if(ier_num==0) then
         if(par_omp_maxthreads==-1) then
            nthreads = MAX(1, par_omp_phys, par_omp_logi)
         else
            nthreads = par_omp_maxthreads
         endif
!$       call OMP_SET_NUM_THreadS(nthreads)
      endif
   endif
endif
res_para(0) = 4
res_para(1) = 1
res_para(2) = real(par_omp_maxthreads)
res_para(3) = real(par_omp_phys)
res_para(4) = real(par_omp_logi)
!
end subroutine do_set_parallel
!
!*****7*****************************************************************
!
end module do_set_mod
