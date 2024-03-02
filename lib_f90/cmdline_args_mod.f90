MODULE cmdline_args_mod
!
private
!
public cmdline_args
!
CONTAINS
!                                                                       
!*****7*****************************************************************
SUBROUTINE cmdline_args (local_mpi_myid)
!                                                                       
!     This routine checks for command line arguments and                
!     executes given macros ..                                          
!                                                                       
USE prompt_mod 
USE debug_mod 
USE envir_mod
USE errlist_mod 
USE lib_length
USE lib_macro_func
USE precision_mod
!USE sockets_mod
USE support_mod
use terminal_mod
use jsu_readline
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: local_mpi_myid
!                                                                       
INTEGER, PARAMETER :: marg = 20 
!                                                                       
CHARACTER(LEN=PREC_STRING), DIMENSION(MARG) :: arg ! (marg) 
CHARACTER(LEN=2048)                  :: line = ' '
CHARACTER(LEN=2048)                  :: string = ' '
!CHARACTER(LEN=40)                    :: str 
INTEGER :: iarg, i, j,  ilen , ilena
LOGICAL :: lexist, lautorun, luser
!
!IF(start_dir(1:5)=='/mnt/') THEN 
!   line = 'autorun.mac'
!   ilen = LEN_TRIM(line)
!   INQUIRE(FILE=line(1:ilen), EXIST=lexist)
!   IF(lexist) THEN
!      CALL file_kdo(line(1:ilen), ilen) ! Execute macro and return to normal prompt
!   ENDIF
!   RETURN    ! Ignore all command line args at WINDOWS 
!ENDIF
!
lautorun = .TRUE.                ! Assume to run macro: "autorun.mac"
luser    = .FALSE.               ! Assume no user macro
!                                                                       
CALL do_getargs (iarg, arg, marg) 
!
IF(iarg>0) THEN                  ! First test for "-noautorun"
   is_noauto: DO i=1, iarg
      IF(INDEX(arg(i), '-noautorun') > 0) THEN
         lautorun = .FALSE.      ! Turn off execution of autorun macro
         arg(i) = ' '            ! Clear current argument
         DO j=i+1,iarg           ! Shift down remaining arguments
            arg(j-1) = arg(j)
         ENDDO
         EXIT is_noauto
      ENDIF
   ENDDO is_noauto
   IF(.NOT.lautorun) iarg = iarg - 1   ! We deleted one argument, reduce total number
ENDIF
!
IF (iarg.gt.0) THEN 
   IF (index (arg (1) , '-macro') .ne. 0) THEN ! execute a macro with optional parameters
      IF (iarg.gt.1) THEN 
         ilena = len_str(arg(2)) ! arg2 is the actual macro name
         inquire(file=arg(2)(1:len_trim(arg(2))), exist=lexist)
         if(lexist) then         ! Macro file was found
         IF(iarg > 2) THEN       ! There are macro parameter(s)
            line  = arg(2)(1:ilena) // ' '
            ilen  = ilena + 1
         ELSE                    ! No macro parameters follow
            line  = arg(2)(1:ilena)
            ilen  = ilena
         ENDIF
         loop_args: DO i = 3, iarg          ! all further args are macro parameters
           ilena = len_str(arg(i))
           if(arg(i)(ilena:ilena)==',') then   ! User added comma
              arg(i)(ilena:ilena) = ' '
              ilena = ilena -1  ! Drop user comma
           endif
           if(arg(i)(1:1)==',') then  !User did ' 1 ,2'
              arg(i)(1:1)= ' '
           endif
           if(arg(i) == ' ') then     ! User did ' 1 , 2 ', skip single comma
              cycle loop_args
           endif
           line  = line(1:ilen) // arg(i)(1:ilena)
           ilen  = ilen + ilena
            IF ( i.lt. iarg) THEN ! separate all but last by ','
               line = line(1:ilen) // ','
               ilen = ilen + 1
            ENDIF
         ENDDO loop_args
         luser = .TRUE.           ! Found user macro
         string = '@' // line(1:len_trim(line))
         call iso_readline_add(string)
         IF(local_mpi_myid==0) WRITE(output_io, 1000) trim(color_high), &
               trim(color_fg), line (1:ilen)
         CALL file_kdo(line(1:ilen), ilen) ! Execute macro and return to normal prompt
         else
            IF(local_mpi_myid==0) WRITE(output_io, 1100) trim(color_high), &
               trim(color_fg), arg(2)(1:len_trim(arg(2)))
         endif
      ENDIF
   ELSE ! all other command line arguments
      DO i = 1, iarg 
!         IF (index (arg (i) , '-remote') .ne.0) THEN 
!            lsocket = .true. 
!         ELSEIF (index (arg (i) , '-port') .ne.0) THEN 
!            str = arg (i) (index (arg (i) , '=') + 1:len_str (arg (i) ))
!            READ (str, * ) s_port 
!         ELSEIF (index (arg (i) , '-access') .ne.0) THEN 
!            s_ipallowed = arg (i) (index (arg(i),'=')+1:len_str(arg(i)))
!         ELSEIF (index (arg (i) , '-help') .ne.0) THEN 
         IF (index (arg (i) , '-help') .ne.0) THEN 
            IF(local_mpi_myid==0) WRITE(output_io, 2000) pname (1:len_str (pname) ) 
         ELSEIF (index (arg (i) , '-debug') .ne.0) THEN 
            IF(local_mpi_myid==0) WRITE(output_io, 1500) 
            dbg = .true. 
         ELSE         ! several macros WITHOUT parameters
            ilen = len_str (arg (i) ) 
            IF(local_mpi_myid==0) WRITE(output_io, 1000) arg (i) (1:ilen) 
            CALL file_kdo(arg(i)(1:ilen), ilen) 
         ENDIF 
      ENDDO 
   ENDIF
!   IF (lsocket) call remote_init 
ENDIF 
!
IF(lautorun) THEN
   line = 'autorun.mac'
   ilen = LEN_TRIM(line)
   INQUIRE(FILE=line(1:ilen), EXIST=lexist)
   IF(lexist) THEN
      if(local_mpi_myid==0) then
         write(output_io,'(5a)') trim(color_high),' ------ > Reading ', &
             current_dir(1:current_dir_l),'autorun.mac', trim(color_fg)
         if(luser) then
            write(output_io,'(4a)') trim(color_high), &
                ' ------ > autorun will be run prior to: ', &
                trim(color_fg),  arg(2)(1:len_trim(arg(2)))
         endif
      endif
      CALL file_kdo(line(1:ilen), ilen) ! Execute macro and return to normal prompt
   ENDIF
ENDIF
!                                                                       
 1000 FORMAT     (a,' ------ > Reading macro: ',a,a) 
 1100 FORMAT     (a, ' ------ > Did not find macro: ',a,a) 
 1500 FORMAT     (' ------ > Running in debug mode ..') 
 2000 FORMAT     (' Usage: ',a,' [-remote] [-debug] [-port=p]',         &
     &                   ' [-access=ip]')                               
!
END SUBROUTINE cmdline_args                   
!
!*****7*****************************************************************
!
END MODULE cmdline_args_mod
