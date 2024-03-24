subroutine experi_mache_kdo(line, lend, length)
!-
!  Main command handling routine to EXPERI
!+
!
use experi_reset
!
use exp2pdf_menu
!
use blanks_mod
use charact_mod
use calc_expr_mod
use errlist_mod
use exit_para_mod
use lib_macro_func
use kdo_all_mod
use precision_mod
use str_comp_mod
!
implicit none
!
character(len=*), intent(inout) :: line     ! User command line
logical         , intent(inout) :: lend     ! End flag
integer         , intent(inout) :: length   ! command line length
!
character(len=PREC_STRING) :: zeile
character(len=5) :: cmd      ! User command
integer          :: lcmd     ! Command length
integer          :: indxt    ! Position of a tabulator
integer          :: indxb    ! Position of a blank
integer          :: indxg    ! Position of a '='
integer          :: lcomm    ! length command parameter
!
! Comment, return immediately
!
if(line(1:1)==' ' .or. line(1:1)=='#' .or. line(1:1)=='!' .or. line==char(13)) return
!
cmd = ' '     ! User command
indxt = index(line, tab)       ! find a tabulator
if(indxt==0) indxt = length + 1
indxb = index (line, ' ')       ! find a blank
if(indxb==0) indxb = length + 1
indxb = MIN(indxb,indxt)
lcmd  = min(indxb - 1, len(cmd))
cmd   = line(1:lcmd)
!
!------ command parameters start at the first character following       
!       the blank                                                       
!                                                                       
zeile = ' '
lcomm = 0
if(indxb + 1<=length) then
   zeile = line(indxb + 1:length)
   lcomm = length - indxb
   call rem_leading_bl(zeile, lcomm)
endif
!
! look for '='
indxg = index(line,'=')
is_math: if(indxg/=0                                                 &
            .and. .not.(str_comp (cmd, 'echo',   2, lcmd, 4) )       &
            .and. .not.(str_comp (cmd, 'system', 2, lcmd, 6) )       &
            .and. .not.(str_comp (cmd, 'fput',   2, lcmd, 4) )       &
            .and. .not.(str_comp (cmd, 'help',   2, lcmd, 4) .or.    &
                        str_comp (cmd, '?   ',   2, lcmd, 4) )       &
            .and. index(line,'==') == 0                            ) then
   call do_math(line, indxg, length)
   return
endif is_math
!
!  cpara(:) = ' '
!  lpara(:) = 0
!  werte(:) = 0.0
!                                                                 
!     --execute a macro file                                      
!                                                                 
is_befehl: if(cmd(1:1) .eq.'@') then
   if(length >= 2) then
       line = line(2:length)
       length = length -1
       call file_kdo(line, length)
   else
      ier_num = -13
      ier_typ = ER_MAC
   endif
elseif(str_comp(cmd, 'exit', 2, lcmd, 4)) then  is_befehl
   lend = .TRUE.
                                                                       
!-------Transform experimental powder pattern to PDF 'exp2pdf'
!                                                                       
elseif(str_comp(cmd, 'exp2pdf', 5, lcmd, 7)) then
   call exp2pdf
!                                                                 
!------- reset to system start 'reset'
!
elseif(str_comp(cmd, 'reset', 3, lcmd, 5)) then
   call experi_do_reset
!
!--- Try general command
!
else is_befehl
   call kdo_all(cmd, lcmd, zeile, lcomm)
   if(zeile == 'EXIT') then ! kdo_all detected "continue suite"
      lend = .TRUE.
   endif
endif is_befehl
!                                                                       
if(ex_do_exit) lend = .true.   ! A global exit was flagged
!
end subroutine experi_mache_kdo
