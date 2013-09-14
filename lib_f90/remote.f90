      PROGRAM client 
!                                                                       
      INTEGER cid 
      INTEGER port 
      INTEGER ilen 
      CHARACTER(80) host, str, resp 
!                                                                       
      str = '' 
!                                                                       
      WRITE ( * , * ) '** DISCUS SOCKET CLIENT **' 
      WRITE ( *, * ) 
      WRITE ( *, 1) 
      READ ( * , '(a)') host 
      WRITE ( *, 2) 
      READ ( *, * ) port 
      WRITE ( *, * ) 
!                                                                       
    1 FORMAT (1x,'Host name of server to connect to : ',$) 
    2 FORMAT (1x,'Port number (> 8000)              : ',$) 
!                                                                       
      ilen = len_str (host) 
      CALL socket_connect (sock, host, ilen, port) 
!                                                                       
      WRITE ( * , * ) 'Connected to server ..' 
!                                                                       
      DO while (str.ne.'exit'.and.str.ne.'bye') 
      WRITE ( *, 3) 
    3 FORMAT   (1x,'Command (exit or bye to end) > ',$) 
      READ ( * , '(a)') str 
      ilen = len_str (str) 
      IF (ilen.gt.0) then 
         CALL socket_send (sock, str (1:ilen), ilen) 
         ilen = 80 
         resp = '' 
         DO while (resp (1:ilen)                                        &
         .ne.'ready'.and.str.ne.'exit'.and.str.ne.'bye')                
         CALL socket_get (sock, resp, ilen) 
         WRITE ( *, 4) resp (1:ilen), ilen 
         ENDDO 
    4 FORMAT     (1x,'Server responded: ',a,' (',i3,' chars)') 
      ENDIF 
      ENDDO 
!                                                                       
      CALL socket_close (sock) 
      END PROGRAM client                            
!*****7***************************************************************  
      INTEGER function len_str (string) 
!-                                                                      
      CHARACTER ( * ) string 
      INTEGER laenge, i 
!                                                                       
      laenge = len (string) 
      i = laenge 
      DO while (i.gt.0.and.string (i:i) .eq.' ') 
      i = i - 1 
      ENDDO 
      IF (i.ge.1) then 
         IF (iachar (string (i:i) ) .eq.13) i = i - 1 
      ENDIF 
      len_str = i 
      END FUNCTION len_str                          
