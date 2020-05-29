MODULE sockets_mod
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE remote_init 
!                                                                       
!     This initializes the remote communication via SOCKETS             
!                                                                       
      USE debug_mod 
      USE envir_mod 
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      LOGICAL llocal 
      INTEGER socket_init
!                                                                       
      llocal = (s_ipallowed.eq.'localhost') 
      llocal = llocal.or. (s_ipallowed.eq.'127.0.0.1') 
!                                                                       
      WRITE (output_io, 2000) 
      ier_num = socket_init (s_sock, s_port, llocal) 
!     CALL socket_init (s_sock, s_port, llocal) 
      lconn = .false. 
!                                                                       
 2000 FORMAT     (' ------ > Running in SERVER mode') 
      END SUBROUTINE remote_init                    
!*****7*****************************************************************
      SUBROUTINE do_socket (zeile, lcomm) 
!                                                                       
!     This part contains all the program independent commands.          
!                                                                       
      USE berechne_mod
      USE blanks_mod
      USE calc_expr_mod
      USE debug_mod 
      USE envir_mod 
      USE errlist_mod 
      USE get_params_mod
      USE learn_mod 
      USE precision_mod
      USE prompt_mod 
USE str_comp_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lcomm 
!                                                                       
      INTEGER MAXW 
      PARAMETER (MAXW = 20) 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: string 
      CHARACTER(LEN=PREC_STRING) :: line 
      CHARACTER(LEN=PREC_STRING) :: cpara (MAXW) 
      INTEGER lpara (MAXW) 
      INTEGER ianz 
      INTEGER il, i, j 
      INTEGER port 
!                                                                       
      REAL wert 
!                                                                       
      INTEGER socket_connect
      INTEGER socket_get
      INTEGER socket_send
!                                                                       
!     Get parameters                                                    
!                                                                       
      lcomm = -lcomm
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm) 
      IF (ier_num.ne.0) THEN 
         RETURN 
      ENDIF 
!                                                                       
      IF (.not.str_comp (cpara (1) , 'open', 2, lpara (1) , 4)          &
          .and..not.lremote) THEN                                           
         ier_num = - 15 
         ier_typ = ER_IO 
         RETURN 
      ENDIF 
!                                                                       
!     This initializes the remote communication via SOCKETS             
!                                                                       
      IF (str_comp (cpara (1) , 'open', 2, lpara (1) , 4) ) THEN 
         IF (ianz.eq.3) THEN 
            port = nint (berechne (cpara (3), lpara (3) ) ) 
            CALL rem_bl (cpara (2), lpara (2) ) 
            ier_num = socket_connect (s_remote, cpara (2), lpara (2), port)
            IF(ier_num /=  0 ) THEN
               ier_typ = ER_IO
               lremote = .false. 
               RETURN
            ENDIF
!           CALL socket_connect (s_remote, cpara (2), lpara (2),        &
!           port)                                                       
            ier_num = socket_get (s_remote, line, il) 
            IF(socket_status == PROMPT_ON) THEN
            IF(line/='ready') THEN
              WRITE (output_io, 3000) line (1:il)
            ENDIF 
            ENDIF 
            IF(ier_num /=  0 ) THEN
               ier_typ = ER_IO
               lremote = .false. 
               RETURN
            ENDIF
            lremote = .true. 
            WRITE (output_io, 2100) 
         ELSE 
            ier_typ = ER_COMM 
            ier_num = - 6 
         ENDIF 
!                                                                       
!       Close the socket connection. The server keeps running           
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'close', 2, lpara (1) , 5) ) THEN 
         ier_num = socket_send (s_remote, 'bye', 3) 
         IF(ier_num < 0) THEN
            ier_num = -19
            RETURN
         ELSE
            ier_num = 0
         ENDIF
         CALL socket_close (s_remote) 
         lremote = .false. 
         WRITE (output_io, 2120) 
!                                                                       
!       Shut down the server                                            
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'exit', 2, lpara (1) , 4) ) THEN 
         ier_num = socket_send (s_remote, cpara (1), lpara (1) ) 
         IF(ier_num < 0) THEN
            ier_num = -19
            RETURN
         ELSE
            ier_num = 0
         ENDIF
         CALL socket_close (s_remote) 
         lremote = .false. 
         WRITE (output_io, 2110) 
!                                                                       
!       Transfer the value of a client variable to the server variable  
!       the input parameter 3 is evaluated and a line with              
!       <cpara(2)> = <value> is send to the server                      
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'transfer', 2, lpara (1) , 8) )     &
      THEN                                                              
         IF (ianz.eq.3) THEN 
            i = lpara (3) 
            zeile = '('//cpara (3) (1:i) //')' 
            i = i + 2 
            wert = berechne (zeile, i) 
!                                                                       
            IF (ier_num.eq.0) THEN 
               zeile (1:lpara (2) ) = cpara (2) (1:lpara (2) ) 
               WRITE (zeile (lpara (2) + 1:lpara (2) + 16), 4000) wert 
               lcomm = lpara (2) + 16 
               ier_num = socket_send (s_remote, zeile, lcomm) 
               IF(ier_num < 0) THEN
                  ier_num = -19
                  RETURN
               ELSE
                  ier_num = 0
               ENDIF
               ier_num = socket_get (s_remote, line, il) 
               IF(socket_status == PROMPT_ON  ) THEN
               IF(line/='ready') THEN
                 WRITE (output_io, 3000) line (1:il)
               ENDIF 
               ENDIF 
               IF(ier_num /=  0 ) THEN
                  ier_typ = ER_IO
                  lremote = .false. 
                  RETURN
               ENDIF
            ENDIF 
!                                                                       
         ELSEIF (ianz.gt.3) THEN 
            IF (index (cpara (3) , '"') .gt.0) THEN 
               zeile = cpara (4) (1:lpara (4) ) //',' 
               i = lpara (4) + 1 
               DO j = 5, ianz 
               zeile = zeile (1:i) //cpara (j) (1:lpara (j) ) //',' 
               i = i + lpara (j) + 1 
               ENDDO 
               zeile (i:i) = ' ' 
               i = i - 1 
               CALL berechne_char (zeile, i) 
               string = cpara (2) (1:lpara (2) ) //'='//cpara (3)       &
               (1:lpara (3) ) //", "//zeile (1:i)                       
               lcomm = lpara (2) + 1 + lpara (3) + 1 + i 
               ier_num = socket_send (s_remote, string, lcomm) 
               IF(ier_num < 0) THEN
                  ier_num = -19
                  RETURN
               ELSE
                  ier_num = 0
               ENDIF
               ier_num = socket_get (s_remote, line, il) 
               IF(socket_status == PROMPT_ON) THEN
               IF(line/='ready') THEN
                 WRITE (output_io, 3000) line (1:il)
               ENDIF 
               ENDIF 
               IF(ier_num /=  0 ) THEN
                  ier_typ = ER_IO
                  lremote = .false. 
                  RETURN
               ENDIF
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
               RETURN 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
            RETURN 
         ENDIF 
!                                                                       
!     socket send subcommand, send the rest of the line as is           
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'send', 2, lpara (1) , 4) ) THEN 
         ier_num = socket_send (s_remote, cpara (2), lpara (2) ) 
         IF(ier_num < 0) THEN
            ier_num = -19
            lremote = .false. 
            RETURN
         ELSE
            ier_num = 0
         ENDIF
         CALL socket_wait 
         IF(ier_num < 0 ) THEN
            lremote = .false.
            RETURN
         ENDIF
!                                                                       
!     No socket command so send the line as is                          
!                                                                       
      ELSE 
         ier_num = socket_send (s_remote, zeile, lcomm) 
         IF(ier_num < 0) THEN
            ier_num = -19
            ier_typ = ER_IO
            lremote = .false. 
            RETURN
         ELSE
            ier_num = 0
         ENDIF
         CALL socket_wait 
         IF(ier_num < 0 ) THEN
            lremote = .false.
            RETURN
         ENDIF
      ENDIF 
!                                                                       
 2100 FORMAT    (1x,'Connected ..') 
 2110 FORMAT    (1x,'Terminated by client ..') 
 2120 FORMAT    (1x,'Connection closed by client ..') 
 3000 FORMAT    (1x,'Server : ',a) 
 4000 FORMAT    ('=',e15.8e2) 
      END SUBROUTINE do_socket                      
!*****7*****************************************************************
      SUBROUTINE socket_wait 
!                                                                       
!     Waits for 'ready' from server                                     
!                                                                       
      USE errlist_mod 
USE lib_length
      USE prompt_mod 
USE str_comp_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER(256) line, cstr 
      INTEGER i, il, icr 
!                                                                       
      INTEGER socket_get 
!                                                                       
      line = '' 
      il = 0 
      icr = 0 
!                                                                       
      DO WHILE (.not.str_comp (line, 'ready', 5, il, 5) ) 
         ier_num = socket_get (s_remote, line, il) 
         IF( ier_num /= 0 ) THEN
            ier_typ = ER_IO
            RETURN
         ENDIF
         DO i = 1, il 
            IF (iachar (line (i:i) ) .eq.10.or.iachar (line (i:i) ) .eq.13) icr=i
         ENDDO 
         IF (icr.ne.0) THEN 
            IF(socket_status == PROMPT_ON) THEN
               IF(line/='ready') THEN
                 WRITE (output_io, 3000) line (1:icr - 1)
               ENDIF 
            ENDIF 
            cstr = line (icr + 1:il) 
            line = '' 
            il = len_str (cstr) 
            line = cstr (1:il) 
         ELSE 
            IF(socket_status == PROMPT_ON) THEN
               IF(line/='ready') THEN
                 WRITE (output_io, 3000) line (1:il)
               ENDIF 
            ENDIF 
         ENDIF 
      ENDDO 
!                                                                       
 3000 FORMAT    (1x,'Server : ',a) 
      END SUBROUTINE socket_wait                    
!*****7*****************************************************************
END MODULE sockets_mod
