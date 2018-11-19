MODULE diffev_set_gen_mod
!
!*******************************************************************************
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE set_gen(gen)
!
USE compare
USE population
USE errlist_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: gen
!
INTEGER, PARAMETER :: IRD = 17
INTEGER, PARAMETER :: IWR = 18
INTEGER, PARAMETER :: lpar=LEN(pop_name)
!
CHARACTER(LEN=1024)  :: line
CHARACTER(LEN=1024)  :: logfiles
CHARACTER(LEN=1024)  :: sumfiles
CHARACTER(LEN=1024)  :: curfiles
CHARACTER(LEN=1024)  :: outfile
CHARACTER(LEN=LEN(pop_name)), DIMENSION(:), ALLOCATABLE :: par_name
INTEGER :: i
INTEGER :: generation
INTEGER :: member
INTEGER :: children
INTEGER :: dimensions
INTEGER :: nlines
INTEGER :: lines_gen
INTEGER :: ios
INTEGER :: length
LOGICAL :: lexist
!
IF(gen==0) RETURN         ! No need to do anything
!
INQUIRE(FILE='GENERATION',  EXIST=lexist)
IF(.NOT.lexist) THEN
   ier_num = -35          ! GENERATION DOES NOT EXIST
   ier_typ = ER_APPL
   RETURN
ENDIF
OPEN(IRD, FILE='GENERATION', STATUS='OLD', ACCESS='SEQUENTIAL')
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) GOTO 9999   ! Error reading GENERATION
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) GOTO 9999   ! Error reading GENERATION
READ(line,*, IOSTAT=ios) generation, member, children, dimensions
IF(ios/=0) GOTO 9999   ! Error reading GENERATION
IF(gen==generation) THEN
   CLOSE(IRD)             ! No need to do anything
   RETURN
ELSEIF(gen>GENERATION) THEN
   ier_num = - 36
   ier_typ = ER_APPL
   CLOSE(IRD)
   RETURN
ENDIF
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) GOTO 9999   ! Error reading GENERATION
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) GOTO 9999   ! Error reading GENERATION
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) GOTO 9999   ! Error reading GENERATION
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) GOTO 9999   ! Error reading GENERATION
!
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) GOTO 9999   ! Error reading GENERATION
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) GOTO 9999   ! Error reading GENERATION
logfiles = line(1:LEN_TRIM(line))
!
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) GOTO 9999   ! Error reading GENERATION
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) GOTO 9999   ! Error reading GENERATION
sumfiles = line(1:LEN_TRIM(line))
!
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) GOTO 9999   ! Error reading GENERATION
READ(IRD,'(a)', IOSTAT=ios) line 
   IF(ios/=0) GOTO 9999   ! Error reading GENERATION
curfiles = line(1:LEN_TRIM(line))
!
ALLOCATE(par_name(-1:dimensions))
lines_gen = 12
par_loop: DO
   READ(IRD,'(a)', IOSTAT=ios) line
   IF(ios/=0) GOTO 9999   ! Error reading GENERATION
   lines_gen = lines_gen + 1
   IF(line(1:11)=='# Parameter') THEN           ! Got the list of parameter names
      DO i=-1, dimensions
         READ(IRD,'(a)', IOSTAT=ios) line
         IF(ios/=0) GOTO 9999   ! Error reading GENERATION
         lines_gen = lines_gen + 1
         par_name(i) = line(1:MIN(lpar, LEN_TRIM(line)))
      ENDDO
      EXIT par_loop
   ENDIF
ENDDO par_loop
!
CLOSE(IRD)
!
! Create the new (shortened) summary files
!
IF(sumfiles/=' ') THEN
   nlines = 3 + gen
   sum_loop: DO i=0, dimensions
      WRITE(line,'(a,a,a,a,a,i5,a,a,a,a,a)') 'cat ',                   &
         sumfiles(1:LEN_TRIM(sumfiles)), '.',                          &
         par_name(i)(1:LEN_TRIM(par_name(i))),                         &
         ' | head -n ',nlines, ' > ', sumfiles(1:LEN_TRIM(sumfiles)),  &
         '.',                                                          &
         par_name(i)(1:LEN_TRIM(par_name(i))), '.NEW'
      length = LEN_TRIM(line)
      CALL do_operating(line, length)
      WRITE(line,'(a,a,a,a,a,a,a,a)') 'mv ',                           &
         sumfiles(1:LEN_TRIM(sumfiles)), '.',                          &
         par_name(i)(1:LEN_TRIM(par_name(i))), '.NEW ',                &
         sumfiles(1:LEN_TRIM(sumfiles)), '.',                          &
         par_name(i)(1:LEN_TRIM(par_name(i)))
      length = LEN_TRIM(line)
      CALL do_operating(line, length)
   ENDDO sum_loop
ENDIF
!
! Create the new (shortened) parameter log files
!
IF(logfiles/=' ') THEN
   nlines = 1 + gen*(member+2)
   log_loop: DO i=0, dimensions
      WRITE(line,'(a,a,a,a,a,i5,a,a,a,a,a)') 'cat ',                   &
         logfiles(1:LEN_TRIM(logfiles)), '.',                          &
         par_name(i)(1:LEN_TRIM(par_name(i))),                         &
         ' | head -n ',nlines, ' > ', logfiles(1:LEN_TRIM(logfiles)),  &
         '.',                                                          &
         par_name(i)(1:LEN_TRIM(par_name(i))), '.NEW'
      length = LEN_TRIM(line)
      CALL do_operating(line, length)
      WRITE(line,'(a,a,a,a,a,a,a,a)') 'mv ',                           &
         logfiles(1:LEN_TRIM(logfiles)), '.',                          &
         par_name(i)(1:LEN_TRIM(par_name(i))), '.NEW ',                &
         logfiles(1:LEN_TRIM(logfiles)), '.',                          &
         par_name(i)(1:LEN_TRIM(par_name(i)))
      length = LEN_TRIM(line)
      CALL do_operating(line, length)
   ENDDO log_loop
ENDIF
!
! Create the new (shortened) current log files
!
IF(curfiles/=' ') THEN
   nlines = (member+2)
   cur_loop: DO i=0, dimensions
      WRITE(outfile,'(a,a,a,a)') curfiles(1:LEN_TRIM(curfiles)),'.',   &
         par_name(i)(1:LEN_TRIM(par_name(i))), '.NEW'
      OPEN(IWR, FILE=outfile, STATUS='unknown')
      WRITE(IWR, '(a)') '#C Current file by DIFFEV'
      CLOSE(IWR)
      WRITE(line,'(a,a,a,a,a,a,i5,a,a,a,a,a)') 'cat ',                 &
         logfiles(1:LEN_TRIM(logfiles)), '.',                          &
         par_name(i)(1:LEN_TRIM(par_name(i))),'    ',                  &
         ' | tail -n ',nlines, ' >> ', curfiles(1:LEN_TRIM(curfiles)), &
         '.',                                                          &
         par_name(i)(1:LEN_TRIM(par_name(i))), '.NEW'
      length = LEN_TRIM(line)
      CALL do_operating(line, length)
      WRITE(line,'(a,a,a,a,a,a,a,a)') 'mv ',                           &
         curfiles(1:LEN_TRIM(curfiles)), '.',                          &
         par_name(i)(1:LEN_TRIM(par_name(i))), '.NEW ',                &
         curfiles(1:LEN_TRIM(curfiles)), '.',                          &
         par_name(i)(1:LEN_TRIM(par_name(i)))
      length = LEN_TRIM(line)
      CALL do_operating(line, length)
   ENDDO cur_loop
ENDIF
!
! Create the new GENERATION file
!
nlines = lines_gen - 2
outfile = 'GENERATION.NEW'
OPEN(IWR, FILE=outfile, STATUS='unknown')
WRITE(IWR, '(a)') '# generation members children parameters'
WRITE(IWR, '(I8,3I10)') gen, member, children, dimensions
CLOSE(IWR)
WRITE(line,'(a,i5,a)') 'cat GENERATION | tail -n ',nlines,' >> GENERATION.NEW'
length = LEN_TRIM(line)
CALL do_operating(line, length)
line = 'mv GENERATION.NEW GENERATION'
length = LEN_TRIM(line)
CALL do_operating(line, length)
!
DEALLOCATE(par_name)
CALL do_read_values(.TRUE.)   ! Read the new refinement state
RETURN                    ! SUCCESS
!
9999 CONTINUE             ! Error reading GENERATION
CLOSE(IRD)
IF(ALLOCATED(par_name)) DEALLOCATE(par_name)
ier_num = -19
ier_typ = ER_APPL
ier_msg(1) = 'Last line read successfully from GENERATION'
ier_msg(2) = line(1:MIN(len_trim(line), LEN(ier_msg)))
RETURN
!
!
END SUBROUTINE set_gen
!
!*******************************************************************************
!
END MODULE diffev_set_gen_mod
