MODULE do_set_mod
!
CONTAINS
!*****7***********************************************************      
SUBROUTINE do_set (zeile, lp) 
!-                                                                      
!     Sets the value of status variables                                
!+                                                                      
USE debug_mod 
USE envir_mod 
USE errlist_mod 
USE do_wait_mod
USE get_params_mod
USE lib_length
USE precision_mod
USE prompt_mod 
USE str_comp_mod
USE support_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: MAXW = 5
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
INTEGER          , INTENT(INOUT) :: lp
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
CHARACTER(LEN=PREC_STRING) :: message
CHARACTER(LEN=80)   :: logfile 
INTEGER :: ios
INTEGER :: ianz 
LOGICAL :: llog 
!                                                                       
!                                                                       
IF (zeile.ne.' ') THEN 
   CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
   IF (ier_num.eq.0) THEN 
!                                                                       
!----- ---- set error                                                   
!                                                                       
      IF (str_comp (cpara (1) , 'error', 1, lpara (1) , 5) ) THEN 
         IF(ianz.eq.3) THEN 
            IF(str_comp(cpara(3), 'save', 2, lpara(3), 4)) THEN
               ier_sav = ier_sta
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
               RETURN
            ENDIF 
            ianz = ianz -1
         ENDIF
         IF(ianz.eq.2) THEN 
            IF(str_comp(cpara(2), 'cont', 2, lpara(2), 4)) THEN
               ier_sta = ER_S_CONT 
            ELSEIF(str_comp(cpara(2), 'exit', 2, lpara (2), 4)) THEN
               ier_sta = ER_S_EXIT 
            ELSEIF(str_comp(cpara(2), 'live', 2, lpara (2), 4)) THEN
               ier_sta = ER_S_LIVE 
            ELSEIF(str_comp(cpara(2), 'old' , 2, lpara (2), 3)) THEN
               ier_sta = ier_sav
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
!----- ---- set parallel                                                
!                                                                       
      ELSEIF(str_comp(cpara(1), 'parallel', 3, lpara(1), 8)) THEN
         CALL do_set_parallel(ianz, cpara, lpara, MAXW)
!                                                                       
!----- ---- set prompt                                                  
!                                                                       
      ELSEIF(str_comp(cpara(1), 'prompt', 3, lpara(1), 6)) THEN
         IF (ianz.ge.2) THEN 
!                                                                       
!------ ------- Third+1  optional parameter: save previous setting      
!                                                                       
            IF (ianz.eq.4) THEN 
               prompt_status_old = prompt_status 
               output_status_old = output_status 
            ENDIF 
!                                                                       
!------ ------- First+1 parameter PROMPT                                
!                                                                       
            IF(str_comp(cpara(2), 'on', 2, lpara(2) , 2)) THEN
               prompt_status = PROMPT_ON 
               socket_status = PROMPT_ON 
!                                                                       
            ELSEIF(str_comp(cpara(2), 'off', 2, lpara(2), 2)) THEN
               prompt_status = PROMPT_OFF 
               socket_status = PROMPT_OFF 
!                                                                       
            ELSEIF(str_comp(cpara(2), 'redirect', 1, lpara(2), 8)) THEN
               prompt_status = PROMPT_REDIRECT 
               socket_status = PROMPT_REDIRECT 
!                                                                       
            ELSEIF(str_comp(cpara(2), 'old', 1, lpara(2), 3)) THEN
               IF (output_status.ne.OUTPUT_SCREEN) THEN 
                  CLOSE (output_io) 
               ENDIF 
               output_status = output_status_old 
               prompt_status = prompt_status_old 
               socket_status = socket_status_old 
               IF (output_status.eq.OUTPUT_SCREEN) THEN 
                  output_io = 6 
               ELSEIF (output_status.eq.OUTPUT_NONE) THEN 
                  output_io = 37 
                  OPEN (unit = output_io, file = nullfile, status &
                        = 'unknown',IOSTAT=ios,IOMSG=message)                                    
                  IF(ios/=0) THEN
                     ier_num = -2
                     ier_typ = ER_IO
                     ier_msg(1) ='Could not open the NULL file'
                     ier_msg(2) = message(1:MAX(1,MIN(80,LEN_TRIM(message))))
                     RETURN
                  ENDIF
               ELSEIF (output_status.eq.OUTPUT_FILE) THEN 
                  output_io = 37 
                  logfile = pname (1:len_str (pname) ) //'.log' 
                  INQUIRE (file = logfile, exist = llog) 
                  IF (llog) THEN 
                     CALL oeffne_append (output_io, logfile, 'old')
                     IF (ier_num.ne.0) RETURN 
!DBG                    open(unit=output_io,file=logfile,               
!DBG     &                   status='old',access='append')              
!DBG95     &                 status='old',access='sequential')          
                  ELSE 
                     OPEN (unit = output_io, file = logfile,      &
                           status = 'new',IOSTAT=ios,IOMSG=message)
                     IF(ios/=0) THEN
                        ier_num = -2
                        ier_typ = ER_IO
                        ier_msg(1) ='Could not open the logfile'
                        ier_msg(2) = logfile
                        ier_msg(3) = message(1:MAX(1,MIN(80,LEN_TRIM(message))))
                        RETURN
                     ENDIF
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ ------- Second+1 optional parameter general output              
!                                                                       
            IF (ianz.ge.3) THEN 
               IF(str_comp(cpara(3), 'on', 2, lpara(3), 2)) THEN
                  IF (output_status.ne.OUTPUT_SCREEN) THEN 
                     CLOSE (output_io) 
                  ENDIF 
                  output_status = OUTPUT_SCREEN 
                  output_io = 6 
!                                                                       
               ELSEIF(str_comp(cpara(3), 'off', 2, lpara(3), 3)) THEN
                  IF (output_status.ne.OUTPUT_SCREEN) THEN 
                     CLOSE (output_io) 
                  ENDIF 
                  output_status = OUTPUT_NONE 
                  output_io = 37 
                  OPEN (unit = output_io, file = nullfile, status &
                        = 'unknown',IOSTAT=ios,IOMSG=message)                                    
                  IF(ios/=0) THEN
                     ier_num = -2
                     ier_typ = ER_IO
                     ier_msg(1) ='Could not open the logfile'
                     ier_msg(2) = logfile
                     ier_msg(3) = message(1:MAX(1,MIN(80,LEN_TRIM(message))))
                     RETURN
                  ENDIF
!                                                                       
               ELSEIF(str_comp(cpara(3), 'file', 2, lpara(3), 4)) THEN
                  IF (output_status.ne.OUTPUT_SCREEN) THEN 
                     CLOSE (output_io) 
                  ENDIF 
                  output_status = OUTPUT_FILE 
                  output_io = 37 
                  logfile = pname (1:len_str (pname) ) //'.log' 
                  INQUIRE (file = logfile, exist = llog) 
                  IF (llog) THEN 
                     CALL oeffne_append (output_io, logfile, 'old')
                     IF (ier_num.ne.0) RETURN 
                  ELSE 
                     OPEN (unit = output_io, file = logfile,      &
                           status = 'new',IOSTAT=ios,IOMSG=message)
                     IF(ios/=0) THEN
                        ier_num = -2
                        ier_typ = ER_IO
                        ier_msg(1) ='Could not open the logfile'
                        ier_msg(2) = logfile
                        ier_msg(3) = message(1:MAX(1,MIN(80,LEN_TRIM(message))))
                        RETURN
                     ENDIF
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
!                                                                       
            IF (dbg) THEN 
               WRITE ( *, 5000) prompt_status, prompt_status_old 
               WRITE ( *, 5010) output_status, output_status_old 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
!----- ---- set debug                                                   
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'debug', 2, lpara (1) , 5) ) THEN
         IF (ianz.eq.2) THEN 
            dbg = (cpara (2) (1:2) .eq.'on') 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
!----- ---- set wait                                                   
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'wait', 2, lpara (1) , 4) ) THEN
         IF (ianz.eq.2) THEN 
            wait_active = (cpara (2) (1:2) .eq.'on') 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
   ENDIF 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
 5000 FORMAT    (' debug  > Prompt (current, old) : ',2I3) 
 5010 FORMAT    (' debug  > Output (current, old) : ',2I3) 
!                                                                       
END SUBROUTINE do_set                         
!
!*****7*****************************************************************
!
SUBROUTINE do_set_parallel(ianz, cpara, lpara, MAXW)
!-
! Get parameters for 'set parallel, ...' command
!+
USE ber_params_mod
USE errlist_mod
USE parallel_mod
USE precision_mod
USE take_param_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(INOUT) :: ianz
INTEGER, INTENT(IN)    :: MAXW
CHARACTER(LEN=*)  , DIMENSION(MAXW), INTENT(INOUT) :: cpara
INTEGER           , DIMENSION(MAXW), INTENT(INOUT) :: lpara
!
INTEGER, PARAMETER :: NOPTIONAL = 2
INTEGER, PARAMETER :: O_NTHREAD = 1                  ! Maximum difference
INTEGER, PARAMETER :: O_USEOMP  = 2                  ! Number of feedbacks to average
CHARACTER(LEN=   7), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(cpara))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte    ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate
!
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
!
DATA oname  / 'nthread', 'useomp ' /   ! 
DATA loname /  7       ,  6        /
!
opara  = (/'all   ', 'use   ' /)
lopara = (/ 3      ,  3       /)
owerte = (/ -1.    ,  0.0000  /)
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                              oname, loname, opara, lopara, lpresent, owerte)
par_omp_use        = .TRUE. ! DEFAULT: User wants to use OMP
par_omp_maxthreads = -1     ! DEFAULT: Maximum number of threads to use, -1==use all available
!
IF(ier_num==0) THEN
   IF(opara(O_USEOMP) == 'use') THEN
      par_omp_use = .TRUE.
   ELSEIF(opara(O_USEOMP) == 'serial') THEN
      par_omp_use = .FALSE.
   ELSE
      ier_num = -6
      ier_typ = ER_COMM
      ier_msg(1) = 'Optional OMP must be omp:use or omp:serial'
      RETURN
   ENDIF
   IF(par_omp_use) THEN
      IF(opara(O_NTHREAD) == 'all') THEN
         par_omp_maxthreads = -1     ! DEFAULT: Maximum number of threads to use, -1==use all available
      ELSE
         opara (O_USEOMP) = '0.0'
         lopara(O_USEOMP) = 3
         CALL ber_params (ianz, cpara, lpara, werte, MAXW)
         IF(ier_num==0) THEN
            par_omp_maxthreads = NINT(werte(O_NTHREAD))
            IF(par_omp_maxthreads<1) THEN
               ier_num = -18
               ier_typ =  ER_COMM
            ENDIF
         ENDIF
      ENDIF
   ENDIF
ENDIF
!
END SUBROUTINE do_set_parallel
!
!*****7*****************************************************************
!
END MODULE do_set_mod
