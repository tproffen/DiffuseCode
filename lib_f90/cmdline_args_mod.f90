MODULE cmdline_args_mod
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
USE errlist_mod 
USE precision_mod
USE sockets_mod
USE sys_compiler
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: local_mpi_myid
!                                                                       
INTEGER, PARAMETER :: marg = 20 
!                                                                       
CHARACTER(LEN=PREC_STRING), DIMENSION(MARG) :: arg ! (marg) 
CHARACTER(LEN=2048)                  :: line = ' '
CHARACTER(LEN=40)                    :: str 
INTEGER :: iarg, i, ilen , ilena
INTEGER :: len_str 
!                                                                       
CALL do_getargs (iarg, arg, marg) 
IF (iarg.gt.0) THEN 
   IF (index (arg (1) , '-macro') .ne. 0) THEN ! execute a macro with optional parameters
      IF (iarg.gt.1) THEN 
         ilena = len_str(arg(2)) ! arg2 is the actual macro name
         IF(iarg > 2) THEN       ! There are macro parameter(s)
            line  = arg(2)(1:ilena) // ' '
            ilen  = ilena + 1
         ELSE                    ! No macro parameters follow
            line  = arg(2)(1:ilena)
            ilen  = ilena
         ENDIF
         DO i = 3, iarg          ! all further args are macro parameters
           ilena = len_str(arg(i))
           line  = line(1:ilen) // arg(i)(1:ilena)
           ilen  = ilen + ilena
            IF ( i.lt. iarg) THEN ! seperate all but last by ','
               line = line(1:ilen) // ','
               ilen = ilen + 1
            ENDIF
         ENDDO
         IF(local_mpi_myid==0) WRITE ( *, 1000) line (1:ilen)
         CALL file_kdo(line(1:ilen), ilen) ! Execute macro and return to normal prompt
      ENDIF
   ELSE ! all other command line arguments
      DO i = 1, iarg 
         IF (index (arg (i) , '-remote') .ne.0) THEN 
            lsocket = .true. 
         ELSEIF (index (arg (i) , '-port') .ne.0) THEN 
            str = arg (i) (index (arg (i) , '=') + 1:len_str (arg (i) ))
            READ (str, * ) s_port 
         ELSEIF (index (arg (i) , '-access') .ne.0) THEN 
            s_ipallowed = arg (i) (index (arg(i),'=')+1:len_str(arg(i)))
         ELSEIF (index (arg (i) , '-help') .ne.0) THEN 
            IF(local_mpi_myid==0) WRITE ( *, 2000) pname (1:len_str (pname) ) 
         ELSEIF (index (arg (i) , '-debug') .ne.0) THEN 
            IF(local_mpi_myid==0) WRITE ( *, 1500) 
            dbg = .true. 
         ELSE         ! several macros WITHOUT parameters
            ilen = len_str (arg (i) ) 
            IF(local_mpi_myid==0) WRITE ( *, 1000) arg (i) (1:ilen) 
            CALL file_kdo (arg (i) (1:ilen), ilen) 
         ENDIF 
      ENDDO 
   ENDIF
   IF (lsocket) call remote_init 
ENDIF 
!                                                                       
 1000 FORMAT     (' ------ > Reading macro ',a) 
 1500 FORMAT     (' ------ > Running in debug mode ..') 
 2000 FORMAT     (' Usage: ',a,' [-remote] [-debug] [-port=p]',         &
     &                   ' [-access=ip]')                               
!
END SUBROUTINE cmdline_args                   
!
!*****7*****************************************************************
!
END MODULE cmdline_args_mod
