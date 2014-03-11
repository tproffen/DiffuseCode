
#include "debug.h"

       SUBROUTINE errlist
!*****7****************************************************************
!
!       Distributes the error codes into the special error routines.
!       Contains other error handling routines.
!
!*****7****************************************************************
       USE errlist_mod
       USE param_mod 
       IMPLICIT      none
!
!
       if    (ier_typ.eq.ER_NONE) then
         call errlist_none
       elseif(ier_typ.eq.ER_COMM) then
         call errlist_comm
       elseif(ier_typ.eq.ER_FORT) then
         call errlist_fort
       elseif(ier_typ.eq.ER_IO  ) then
         call errlist_io
       elseif(ier_typ.eq.ER_MAC ) then
         call errlist_mac
       elseif(ier_typ.eq.ER_MATH) then
         call errlist_math
       else
         call errlist_appl
       endif
!
!------       Terminate program if an error occured and the 
!       error status is set to ER_S_EXIT
!
       if(ier_num.lt.0 .and. ier_sta.eq.ER_S_EXIT) then
         write(*,1000) char(7)
         stop
       elseif(ier_num.lt.0 .and. ier_sta.eq.ER_S_CONT .or.              &
     &        ier_sta.eq.ER_S_LIVE                  ) then
         res_para(0) = -1
         res_para(1) = ier_num
         res_para(2) = ier_typ
       endif
!
       call no_error
!
1000  format(' ****EXIT**** Program terminated by error status',        &
     &       '        ****',a1/)
       end
!*****7****************************************************************
       SUBROUTINE no_error
!+
!       Resets error variables
!-
       USE errlist_mod 
       IMPLICIT       none
!
!
       MSG("no_error()")
       ier_num = 0
       ier_typ = ER_NONE
!
       ier_msg(1) = ' '
       ier_msg(2) = ' '
       ier_msg(3) = ' '
!
       end
!*****7****************************************************************
       SUBROUTINE disp_error (typ,error,iu,io)
!-
!       Displays the error messages 
!+
       USE errlist_mod 
       USE prompt_mod
       IMPLICIT      none
!
!
       integer       i,iu,io
       CHARACTER*80  estr
       CHARACTER*45  error(iu:io)
       CHARACTER*4   typ
       integer             le
!
       integer       len_str
!
       if(iu.le.ier_num .and. ier_num.le.io) THEN
         IF(    error(ier_num).ne.' ') then
            write(*      ,1000) typ,error(ier_num),ier_num,char(7)
            write(ier_out,1500) typ,error(ier_num)
            do i=1,3
              if (ier_msg(i)     .ne.' ')               &
     &                write(*,1500) typ,ier_msg(i),ier_num
            ENDDO
         else
           write(*,2000) ier_num,typ,char(7)
         endif
       else
         write(*,2000) ier_num,typ,char(7)
       endif
!
       if (lconn .and. lsocket) then
         write(estr,1500) typ,error(ier_num),ier_num
         le=len_str(estr)
         call socket_send(s_conid,estr,le)
       endif
!
1000  format(' ***',a,'*** ',a45,' ***',i4,' ***',a1)
1500  format(' ***',a,'*** ',a45,' ***',i4,' ***')
2000  format(' !!!! No error message for error no.:',I8,' !!!!'/        &
     &       '      Error type:    ',a,'            '/                  &
     &       '      Please document and report to the author',a1)
       end
!*****7****************************************************************
       SUBROUTINE errlist_none
!-
!       Displays error Messages for the error type NONE
!+
       USE errlist_mod 
       IMPLICIT      none
!
!
       integer       iu,io
       parameter    (iu=0,io=3)
!
       CHARACTER*45  error(iu:io)
!
       DATA ERROR (  0:  3) /                         &
     &  'No error',                                 & !  0  ! command
     &  'Done',                                     & !  1  ! command
     &  ' ',                                        & !  2  ! command
     &  'Task aborted'                              & !  3  ! command
     &       /
!
       call disp_error ('NONE',error,iu,io)
       end
!*****7****************************************************************
       SUBROUTINE errlist_comm
!-
!       Displays error Messages for the error type COMMands
!+
       USE errlist_mod 
       USE prompt_mod
       IMPLICIT      none
!
!
       integer       iu,io
       parameter    (IU=-17,IO=0)
!
       CHARACTER*45  error(IU:IO)
!
       DATA ERROR (-17:-01) /                         &
     &  'Too many parameters',                      & !-17  ! command
     &  ' ',                                        & !-16  ! command
     &  ' ',                                        & !-15  ! command
     &  ' ',                                        & !-14  ! command
     &  ' ',                                        & !-13  ! command
     &  ' ',                                        & !-12  ! command
     &  'Error in subroutine',                      & !-11  ! command
     &  ' ',                                        & !-10  ! command
     &  ' ',                                        & ! -9  ! command
     &  'Unknown command',                          & ! -8  ! command
     &  ' ',                                        & ! -7  ! command
     &  'Missing or wrong parameters for command',  & ! -6  ! command,fortran
     &  'Error in operating system command',        & ! -5  ! command
     &  ' ',                                        & ! -4  ! command
     &  'Could not allocate arrays',                & ! -3  ! command
     &  'Command parameter has zero length ',       & ! -2  ! io
     &  '       directory not defined '             & ! -1  ! io
     &       /
       DATA ERROR (  0:  0) /                         &
     &  ' '                                         & !  0  ! command
     &       /
!
       ERROR(-1)(1:6) = pname_cap
!
       call disp_error ('COMM',error,iu,io)
       end
!*****7****************************************************************
       SUBROUTINE errlist_fort
!-
!       Displays error Messages for the error type FORTran
!+
       USE errlist_mod 
       IMPLICIT      none
!
!
       integer       iu,io
       parameter    (IU=-34,IO=1)
!
       CHARACTER*45  ERROR(IU:IO)
!
       DATA ERROR (-34:-21) /                       &
     &  'String has length zero',                   & !-34  ! fortran
     &  'Variable in use; cannot initialize value', & !-33  ! fortran
     &  'Variable name is already defined',         & !-32  ! fortran
     &  'Incomplete (do,if) statement          ',   & !-31  ! fortran
     &  'Right quotation mark missing in format',   & !-30  ! fortran
     &  'Character substring out of bounds',        & !-29  ! fortran
     &  'Too deeply leveled break command',         & !-28  ! fortran
     &  'Function with wrong number of arguments ', & !-27  ! fortran
     &  'Variable name contains illegal characters',& !-26  ! fortran
     &  'Variable name is illegal',                 & !-25  ! fortran
     &  'Variable is not defined',                  & !-24  ! fortran
     &  'Maximum number of int. variables defined', & !-23  ! fortran
     &  'Maximum number of real variables defined', & !-22  ! fortran
     &  'Missing '' while comparing stings'         & !-21  ! fortran
     &   /
       DATA ERROR (-20:-01) /                       &
     &  'Illegal argument for ln(x) function',      & !-20  ! fortran
     &  'Illegal nesting of controll commands',     & !-19  ! fortran
     &  'Unresolvable condition',                   & !-18  ! fortran
     &  'Too many indizes for DO-loop counter',     & !-17  ! fortran
     &  'Too deeply leveled (do,if) construction',  & !-16  ! fortran
     &  'Too many commands',                        & !-15  ! fortran
     &  'Index of DO-loop counter is missing',      & !-14  ! fortran
     &  'Wrong number of indizes for array',        & !-13  ! fortran
     &  'Expression between () is missing',         & !-12  ! fortan
     &  'Number of parentheses is not matching',    & !-11  ! fortan
     &  'Index for array element is missing',       & !-10  ! fortan
     &  'Number of brackets is not matching',       & ! -9  ! fortan
     &  'Index outside array limits',               & ! -8  ! fortran
     &  'Argument for asin,acos greater 1',         & ! -7  ! fortran
     &  'Missing or wrong Parameters for command',  & ! -6  ! command,fortran
     &  'Square root of negative number',           & ! -5  ! fortran
     &  'Division by zero',                         & ! -4  ! fortran
     &  'Unknown intrinsic function',               & ! -3  ! fortran
     &  'Unknown Variable',                         & ! -2  ! fortran
     &  'Nonnumerical Parameters in expression'     & ! -1  ! fortran
     &   /
       DATA ERROR (  0:  1) /                       &
     &  ' ',                                        & !  0  ! command
     &  'Variable name is already defined'          & !  1  ! fortran
     &       /
!
       call disp_error ('FORT',error,iu,io)
       end
!*****7****************************************************************
       SUBROUTINE errlist_io
!-
!       Displays error Messages for the error type IO
!+
       USE errlist_mod
       IMPLICIT      none
!
!
       integer       iu,io
       parameter    (IU= -26,IO=0)
!
       CHARACTER*45  ERROR(IU:IO)
!
       DATA ERROR (-26:-21) /                       &
     &  'Second parameter must be >= first Param.', & ! -26 ! io
     &  'Socket listen problem ',                   & ! -25 ! command
     &  'Socket bind problem ',                     & ! -24 ! command
     &  'Scocket: Rejected connection  ',           & ! -23 ! command
     &  'Socket accept problem ',                   & ! -22 ! command
     &  'Received null string from socket receive'  & ! -21 ! io
     &   /
       DATA ERROR (-20:-01) /                       &
     &  'Problem receiving from socket',            & ! -20 ! io
     &  'Problem sending to socket',                & ! -19 ! io
     &  'Could not open socket connection',         & ! -18 ! io
     &  'Could not grab socket',                    & ! -17 ! io
     &  'Could not resolve hostname for socket',    & ! -16 ! io
     &  'No socket connected',                      & ! -15 ! io
     &  'Filename has zero length',                 & ! -14 ! io
     &  'I/O stream number outside valid range',    & ! -13 ! io
     &  'Error writing to file',                    & ! -12 ! io
     &  'No IO stream open to close',               & ! -11 ! io
     &  'IO stream already open',                   & ! -10 ! io
     &  'Error reading user input',                 & ! -9  ! io
     &  'Nothing learned - no macro written',       & ! -8  ! io
     &  'Learning sequence already in progress',    & ! -7  ! io
     &  'Unexpected end of file',                   & ! -6  ! io
     &  'No such entry in online help',             & ! -5  ! io
     &  'File already exists',                      & ! -4  ! io
     &  'Error reading file',                       & ! -3  ! io
     &  'Error opening file',                       & ! -2  ! io
     &  'File does not exist'                       & ! -1  ! io
     &       /
       DATA ERROR (  0:  0) /                       &
     &  ' '                                         & !  0  ! command
     &       /
!
       call disp_error ('I/O ',error,iu,io)
       end
!*****7****************************************************************
       SUBROUTINE errlist_mac
!-
!       Displays error Messages for the error type MACros
!+
       USE errlist_mod
       IMPLICIT      none
!
!
       integer       iu,io
       parameter    (IU=-41,IO=3)
!
       CHARACTER*45  ERROR(IU:IO)
!
       DATA ERROR (-41:-21) /                       &
     &  'Not enough macro parameters given ',       & !-41  ! macro
     &  ' ',                                        & !-40  ! discus
     &  ' ',                                        & !-39  ! discus
     &  ' ',                                        & !-38  ! math
     &  ' ',                                        & !-37  ! discus
     &  'Unexpected EOF in macro file',             & !-36  ! macro
     &  'Too deeply leveled macros',                & !-35  ! macro
     &  ' ',                                        & !-34  ! discus
     &  ' ',                                        & !-33  ! discus
     &  ' ',                                        & !-32  ! discus
     &  ' ',                                        & !-31  ! discus
     &  ' ',                                        & !-30  ! discus
     &  ' ',                                        & !-29  ! discus
     &  ' ',                                        & !-28  ! fortran
     &  ' ',                                        & !-27  ! discus
     &  ' ',                                        & !-26  ! discus
     &  ' ',                                        & !-25  ! discus
     &  ' ',                                        & !-24  ! discus
     &  ' ',                                        & !-23  ! discus
     &  ' ',                                        & !-22  ! discus
     &  ' '                                         & !-21  ! discus
     &   /
       DATA ERROR (-20:-01) /                       &
     &  ' ',                                        & !-20  ! discus
     &  ' ',                                        & !-19  ! fortran
     &  ' ',                                        & !-18  ! fortran
     &  ' ',                                        & !-17  ! command
     &  ' ',                                        & !-16  ! fortran
     &  ' ',                                        & !-15  ! fortran
     &  ' ',                                        & !-14  ! discus
     &  'Macro name missing on command line',       & !-13  ! macro
     &  'Macro not found',                          & !-12  ! macro
     &  ' ',                                        & !-11  ! command
     &  ' ',                                        & !-10  ! discus
     &  ' ',                                        & ! -9  ! io
     &  ' ',                                        & ! -8  ! command
     &  ' ',                                        & ! -7  ! discus
     &  ' ',                                        & ! -6  ! command,fortran
     &  ' ',                                        & ! -5  ! command
     &  ' ',                                        & ! -4  ! io
     &  ' ',                                        & ! -3  ! io
     &  ' ',                                        & ! -2  ! io
     &  'Too many macro parameters given'           & ! -1  ! macro
     &   /
       DATA ERROR (  0:  3) /                       &
     &  ' ',                                        & !  0  ! command
     &  ' ',                                        & !  1  ! command
     &  ' ',                                        & !  2  ! command
     &  ' '                                         & !  3  ! command
     &       /
!
       call disp_error ('MAC ',error,iu,io)
       end
!*****7****************************************************************
       SUBROUTINE errlist_math
!-
!       Displays error Messages for the error type MATHematical errors
!+
       USE errlist_mod
       IMPLICIT      none
!
!
       integer       iu,io
       parameter    (IU= -1,IO=0)
!
       CHARACTER*45  ERROR(IU:IO)
!
       DATA ERROR/                                  &
     &  'Singular Matrix',                          & !- 1  ! math
     &  ' '                                         & !  0  ! command
     &       /
!
       call disp_error ('MATH',error,iu,io)
       end
!*****7****************************************************************
