MODULE lib_errlist_func
!
CONTAINS
SUBROUTINE errlist
!*****7****************************************************************
!
!       Distributes the error codes into the special error routines.
!       Contains other error handling routines.
!
!*****7****************************************************************
       USE errlist_mod
       USE mpi_slave_mod
       USE param_mod 
       USE prompt_mod 
       USE set_sub_generic_mod
USE class_macro_internal
USE macro_mod
USE terminal_mod
!
       IMPLICIT      NONE
!
       IF(mpi_is_slave .AND. ier_num /= 0 .AND. ier_sta /= ER_S_LIVE) THEN
          mpi_slave_error = ier_num   ! Transfer error for MPI to signal
       ELSE
          mpi_slave_error = 0
       ENDIF
!
       IF    (ier_typ == ER_NONE) THEN
         CALL errlist_none
       ELSEIF(ier_typ == ER_COMM) THEN
         CALL errlist_comm
       ELSEIF(ier_typ == ER_FORT) THEN
         CALL errlist_fort
       ELSEIF(ier_typ == ER_IO  ) THEN
         CALL errlist_io
       ELSEIF(ier_typ == ER_MAC ) THEN
         CALL errlist_mac
       ELSEIF(ier_typ == ER_MATH) THEN
         CALL errlist_math
       ELSE
         CALL p_errlist_appl
       ENDIF
!
IF(lmakro .AND. lmakro_disp) THEN
   WRITE(output_io, 2000) TRIM(color_err), mac_tree_active%current,                 &
     mac_tree_active%level, TRIM(color_fg),                                         &
     TRIM(color_err),                                                               &
     mac_tree_active%active%macrofile(1:LEN_TRIM(mac_tree_active%active%macrofile)),&
     TRIM(color_fg)
   lmakro_disp = .FALSE.   ! Turn macro display off to avoid multiple displays
ENDIF
!------       Terminate program if an error occured and the 
!       error status is set to ER_S_EXIT
!
       IF(ier_num.lt.0 .AND. ier_sta == ER_S_EXIT .AND. &
          .NOT. mpi_is_slave ) THEN
         WRITE( error_io,1000) CHAR(7)
         IF(error_io /=0 ) WRITE(*,1000) CHAR(7)
         STOP
       ELSEIF(ier_num.lt.0 .AND. ier_sta == ER_S_CONT .OR.              &
     &        ier_sta == ER_S_LIVE                  ) THEN
         res_para(0) = -1
         res_para(1) = ier_num
         res_para(2) = ier_typ
       ENDIF
!
       CALL no_error
!
1000  FORMAT(' ****EXIT**** Program terminated by error status',        &
     &       '        ****',a1/)
2000 FORMAT(a,' ***MAC *** Error occured in line:',i5,' Level',i3, &
            '          ***',a,/,                              &
            a,' ***MAC *** ',a,a                                   &
     )
       END
!*****7****************************************************************
!
SUBROUTINE no_error
!+
!       Resets error variables
!-
USE errlist_mod 
USE class_macro_internal
!
IMPLICIT       NONE
!
ier_num = 0
ier_typ = ER_NONE
!
ier_msg(:) = ' '
!
END
!
!*****7****************************************************************
       SUBROUTINE disp_error (typ,error,iu,io)
!-
!       Displays the error messages 
!+
       USE errlist_mod 
USE lib_length
       USE terminal_mod
       USE prompt_mod
       IMPLICIT      NONE
!
!
       INTEGER          , INTENT(IN) :: iu
       INTEGER          , INTENT(IN) :: io
       CHARACTER(LEN=45), INTENT(IN) :: error(iu:io)
       CHARACTER(LEN=4) , INTENT(IN) :: typ
!
       INTEGER            :: i
       CHARACTER(LEN=80)  :: estr
       INTEGER            :: le
!
!
       IF(ier_ctrlc) THEN     ! Avoid multiple message in case of a ctrl-c
          IF(ier_rep .AND. ier_num==-14) RETURN
          ier_rep = .TRUE.
       ENDIF
       IF(iu <= ier_num .AND. ier_num <= io) THEN
         IF(    error(ier_num).ne.' ') THEN
            WRITE(error_io,1000) TRIM(color_err),typ,error(ier_num),ier_num,TRIM(color_fg)!,CHAR(7)
            IF(error_io/=0) &
               WRITE(*    ,1000) TRIM(color_err),typ,error(ier_num),ier_num,TRIM(color_fg),CHAR(7)
            IF(ier_mpi) &
            WRITE(output_io,1100)                 typ,error(ier_num),ier_num
!           WRITE(ier_out,1500) TRIM(color_err),typ,error(ier_num),ier_num,TRIM(color_fg)
            DO i=1,UBOUND(ier_msg,1)
              IF (ier_msg(i) /= ' ')  &
     &                WRITE(*,1500) TRIM(color_err),typ,ier_msg(i),ier_num,TRIM(color_fg)
            ENDDO
         ELSE
           WRITE(error_io,2000) TRIM(color_err),ier_num,typ,TRIM(color_fg),CHAR(7)
           IF(error_io /= 0) WRITE(*,2000) TRIM(color_err),ier_num,typ,TRIM(color_fg),CHAR(7)
         ENDIF
       ELSE
         WRITE(error_io,2000) TRIM(color_err),ier_num,typ,TRIM(color_fg),CHAR(7)
         IF(error_io /=0 ) WRITE(*,2000) TRIM(color_err),ier_num,typ,TRIM(color_fg),CHAR(7)
       ENDIF
!
       IF (lconn .AND. lsocket) THEN
         WRITE(estr,1500) TRIM(color_err),typ,error(ier_num),ier_num,TRIM(color_bg)
         le=len_str(estr)
         CALL socket_send(s_conid,estr,le)
       ENDIF
!
!
1000  FORMAT(a,' ***',a,'*** ',a45,' ***',i4,' ***',a,a1)
1100  FORMAT(  ' ***',a,'*** ',a45,' ***',i4,' ***'     )
1500  FORMAT(a,' ***',a,'*** ',a45,' ***',i4,' ***',a)
2000  FORMAT(a,' !!!! No error message for error no.:',I8,' !!!!'/        &
     &         '      Error type:    ',a,'            '/                  &
     &         '      Please document and report to the author',a,a1)
       END
!*****7****************************************************************
       SUBROUTINE errlist_none
!-
!       Displays error Messages for the error type NONE
!+
       USE errlist_mod 
       IMPLICIT      NONE
!
!
       INTEGER, PARAMETER :: iu = 0
       INTEGER, PARAMETER :: io = 3
!
       CHARACTER(LEN=45)  :: error(iu:io)
!
       DATA ERROR (  0:  3) /                         &
     &  'No error',                                 & !  0  ! command
     &  'Done',                                     & !  1  ! command
     &  ' ',                                        & !  2  ! command
     &  'Task aborted'                              & !  3  ! command
     &       /
!
       CALL disp_error ('NONE',error,iu,io)
       END
!*****7****************************************************************
       SUBROUTINE errlist_comm
!-
!       Displays error Messages for the error type COMMands
!+
       USE errlist_mod 
       USE prompt_mod
       IMPLICIT      NONE
!
!
       INTEGER, PARAMETER :: iu = -17
       INTEGER, PARAMETER :: io =   0
!
       CHARACTER(LEN=45)  error(IU:IO)
!
       DATA ERROR (-17:-01) /                         &
     &  'Too many parameters',                      & !-17  ! command
     &  'Continued line is too long (> 512 charac.',& !-16  ! command
     &  'Illegal branch to section',                & !-15  ! command
     &  'Program was interrupted by CTRL-C ',       & !-14  ! command
     &  'Optional parameter value not recognized',  & !-13  ! command
     &  'Optional parameter name not recognized',   & !-12  ! command
     &  'Error in subroutine',                      & !-11  ! command
     &  'Sub menu returned with error',             & !-10  ! command
     &  'Program section returned with error',      & ! -9  ! command
     &  'Unknown command',                          & ! -8  ! command
     &  'Branch is active within suite only ',      & ! -7  ! command
     &  'Missing or wrong parameters for command',  & ! -6  ! command,fortran
     &  'Error in operating system command',        & ! -5  ! command
     &  'Right quotation mark missing in FORMAT',   & ! -4  ! command
     &  'Could not allocate arrays',                & ! -3  ! command
     &  'Command parameter has zero length ',       & ! -2  ! io
     &  '       directory not defined '             & ! -1  ! io
     &       /
       DATA ERROR (  0:  0) /                         &
     &  ' '                                         & !  0  ! command
     &       /
!
       ERROR(-1)(1:6) = pname_cap(1:6)
!
       CALL disp_error ('COMM',error,iu,io)
       END
!*****7****************************************************************
       SUBROUTINE errlist_fort
!-
!       Displays error Messages for the error type FORTran
!+
       USE errlist_mod 
       IMPLICIT      NONE
!
!
       INTEGER, PARAMETER :: iu = -50
       INTEGER, PARAMETER :: io =   1
!
       CHARACTER(LEN=45)  ERROR(IU:IO)
!
       DATA ERROR ( iu:-41) /                       &
     &  'Value is outside limits',                  & !-50  ! fortran
     &  'Value is Nan',                             & !-49  ! fortran
     &  'Parameters must be matrix variables',      & !-48  ! fortran
     &  'INVMAT only implemented for dimension <=4',& !-47  ! fortran
     &  'Input variable is not a matrix',           & !-46  ! fortran
     &  'Matrix cannot be inverted DET==0',         & !-45  ! fortran
     &  'Matrix dimensions do not match',           & !-44  ! fortran
     &  'Expression is not a Character string  ',   & !-43  ! fortran
     &  'Variable on left side is not CHARACTER',   & !-42  ! fortran
     &  'Variable on left side is READ only'        & !-41  ! fortran
     &   /
       DATA ERROR (-40:-21) /                       &
     &  'Erroneous dimensions for user variable',   & !-40  ! fortran
     &  'Refine variable, use reset in diffev',     & !-39  ! fortran
     &  'Number of hyphenations is not matching',   & !-38  ! fortran
     &  'Mean of lognormal is less than zero',      & !-37  ! fortran
     &  'Skew parameter outside [-1:1]',            & !-36  ! fortran
     &  'Sigma/FWHM is less than zero',             & !-35  ! fortran
     &  'String has length zero',                   & !-34  ! fortran
     &  'Variable in use; cannot initialize value', & !-33  ! fortran
     &  'Variable name is already defined',         & !-32  ! fortran
     &  'Incomplete (do,if) statement          ',   & !-31  ! fortran
     &  'Right quotation mark missing in FORMAT',   & !-30  ! fortran
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
     &  'Illegal nesting of control commands',      & !-19  ! fortran
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
       CALL disp_error ('FORT',error,iu,io)
       END
!*****7****************************************************************
       SUBROUTINE errlist_io
!-
!       Displays error Messages for the error type IO
!+
       USE errlist_mod
       IMPLICIT      NONE
!
!
       INTEGER, PARAMETER :: iu = -29
       INTEGER, PARAMETER :: io =   0
!
       CHARACTER(LEN=45)  ERROR(IU:IO)
!
       DATA ERROR (-29:-21) /                       &
     &  'Number of data points exceeds limit',      & ! -29 ! io
     &  'Format number outside allowed range',      & ! -28 ! io
     &  'Not enough parameter for filename format ',& ! -27 ! io
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
       CALL disp_error ('I/O ',error,iu,io)
       END
!*****7****************************************************************
       SUBROUTINE errlist_mac
!-
!       Displays error Messages for the error type MACros
!+
       USE errlist_mod
       IMPLICIT      NONE
!
!
       INTEGER, PARAMETER :: iu = -41
       INTEGER, PARAMETER :: io =   3
!
       CHARACTER(LEN=45)  ERROR(IU:IO)
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
       CALL disp_error ('MAC ',error,iu,io)
       END
!*****7****************************************************************
       SUBROUTINE errlist_math
!-
!       Displays error Messages for the error type MATHematical errors
!+
       USE errlist_mod
       IMPLICIT      NONE
!
!
       INTEGER, PARAMETER :: iu =  -1
       INTEGER, PARAMETER :: io =   0
!
       CHARACTER(LEN=45)  ERROR(IU:IO)
!
       DATA ERROR/                                  &
     &  'Singular Matrix',                          & !- 1  ! math
     &  ' '                                         & !  0  ! command
     &       /
!
       CALL disp_error ('MATH',error,iu,io)
       END
!*****7****************************************************************
!
SUBROUTINE errlist_save
!-
! Saves error status into temporary error list
!+
USE errlist_mod
!
IMPLICIT      NONE
!
ier_num_tmp    = ier_num
ier_typ_tmp    = ier_typ
ier_msg_tmp(:) = ier_msg(:)
!
END SUBROUTINE errlist_save
!
!*****7****************************************************************
!
SUBROUTINE errlist_restore
!-
! Saves error status into temporary error list
!+
USE errlist_mod
!
IMPLICIT      NONE
!
ier_num    = ier_num_tmp
ier_typ    = ier_typ_tmp
ier_msg(:) = ier_msg_tmp(:)
!
ier_num_tmp    = 0  
ier_typ_tmp    = 0  
ier_msg_tmp(:) = ' '
!
END SUBROUTINE errlist_restore
!
!*****7****************************************************************
!
END MODULE lib_errlist_func
