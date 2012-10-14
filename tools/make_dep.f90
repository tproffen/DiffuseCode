      PROGRAM make_dep
!
      IMPLICIT none
!
      CHARACTER(LEN= 80)  :: infile
      CHARACTER(LEN= 80)  :: line
      CHARACTER(LEN= 80), DIMENSION(500)  :: inc_files
      CHARACTER(LEN= 80), DIMENSION(500)  :: use_files
      CHARACTER(LEN= 80)  :: dummy
      CHARACTER(LEN= 800) :: string
      INTEGER           , DIMENSION(500)  :: l_inc_files
      INTEGER           , DIMENSION(500)  :: l_use_files
      INTEGER             :: io_stat
      INTEGER             :: i
      INTEGER             :: idot
      INTEGER             :: iapp1
      INTEGER             :: iapp2
      INTEGER             :: n_inc
      INTEGER             :: n_use
      LOGICAL             :: l_old
!
      OPEN(3,file='dependencies.list')
      OPEN(1,file='00.DEP')
!
      READ(1,1000,IOSTAT=io_stat) infile
main: DO WHILE ( io_stat == 0) 
!write(*,*) 'INFILE ', infile
         OPEN(2,file=infile)
         idot = index(infile,'.')
         infile(idot+1:80) = 'o'
         n_inc = 0
         n_use = 0
         inc_files = ' ' 
         use_files = ' ' 
         l_inc_files = 0
         l_use_files = 0
         READ(2,1000,IOSTAT=io_stat) line
         DO WHILE ( io_stat == 0) 
            CALL lower_case ( line )
            l_old = .false.
            iapp1 = index(line,'include')
            IF(line(1:iapp1-1) == ' ') THEN
               IF(iapp1 > 0 ) THEN
                  iapp1 = index(line,'''')
                  iapp2 = index(line(iapp1+1:80),'''')
                  dummy = line(iapp1+1:iapp1+iapp2-1)
l_inc:            DO i=1,n_inc
                     IF(inc_files(i) == dummy) THEN
                       l_old = .true.
                       EXIT l_inc
                     ENDIF
                  ENDDO l_inc
                  IF(.not. l_old) THEN
                    n_inc = n_inc + 1
                    inc_files(n_inc)   = dummy
                    l_inc_files(n_inc) = iapp2-1
                  ENDIF
               ENDIF
            ENDIF
            iapp1 = index(line,'use')
            IF(line(1:iapp1-1) == ' ') THEN
               IF(iapp1 > 0 ) THEN
                  iapp1 = index(line,'use')
                  iapp2 = index(line(iapp1+4:80),' ')
                  dummy = line(iapp1+4:iapp1+4+iapp2-2)//'.mod'
                  IF( dummy /= 'mpi.mod') THEN
l_use:               DO i=1,n_use
                        IF(use_files(i) == dummy) THEN
                          l_old = .true.
                          EXIT l_use
                        ENDIF
                     ENDDO l_use
                     IF(.not. l_old) THEN
                       n_use = n_use + 1
                       use_files(n_use)   = dummy(1:iapp2-1+4)
                       l_use_files(n_use) = iapp2-1+4
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            READ(2,1000,IOSTAT=io_stat) line
         ENDDO
         CLOSE(2)
         IF(n_inc > 0) THEN
            WRITE(string,3000) infile(1:idot+1),(inc_files(i)(1:l_inc_files(i)),i=1,n_inc)
            DO i=800,1,-1
               if(string(i:i).ne.' ') EXIT
            ENDDO
            IF(n_use == 0) THEN
              WRITE(3,1000) string(1:i)
            ELSE
              WRITE(3,2000) string(1:i),' \'
              WRITE(3,3100) achar(9)       ,(use_files(i)(1:l_use_files(i)),i=1,n_use)
            ENDIF
         ELSE
            IF(n_use > 0) THEN
               WRITE(3,4100) infile(1:idot+1),(use_files(i)(1:l_use_files(i)),i=1,n_use)
            ENDIF
         ENDIF
!exit main
         READ(1,1000,IOSTAT=io_stat) infile
      ENDDO main
!
      CLOSE(1)
      CLOSE(3)
!
1000  FORMAT(a)
2000  FORMAT(a,a)
3000  FORMAT(a,': ',50(a:,' '),'\')
3100  FORMAT(a,' ',50(a:,' '))
4100  FORMAT(a,': ',50(a:,' '))
      END PROGRAM make_dep
!
      SUBROUTINE lower_case (line )
!
      IMPLICIT NONE
!
      CHARACTER (LEN=*), INTENT(INOUT) :: line
!
      INTEGER   :: laenge
      INTEGER   :: buch
      INTEGER   :: i
      INTEGER, PARAMETER :: aaa = IACHAR('a')
      INTEGER, PARAMETER :: zzz = IACHAR('z')
      INTEGER, PARAMETER :: AA  = IACHAR('A')
      INTEGER, PARAMETER :: ZZ  = IACHAR('Z')
      INTEGER, PARAMETER :: shift = aaa - AA
!
      laenge = LEN(line)
!
      DO i=1,laenge
         buch = IACHAR(line(i:i))
         IF( AA <= buch .and. buch <= ZZ ) THEN
            line(i:i) = ACHAR(buch + shift)
         ENDIF
      ENDDO
!
      END SUBROUTINE lower_case
