MODULE diffev_set_gen_mod
!
!*******************************************************************************
!
CONTAINS
!
!*******************************************************************************
!
subroutine purge_files(line, length)
!-
!  Main interface to purge old generation lines from backup
!+
!
!
use ber_params_mod
use errlist_mod
use get_params_mod
use precision_mod
!
implicit none
!
character(len=*), intent(inout) :: line
integer         , intent(inout) :: length
!
integer, parameter :: MAXW =  2
!
CHARACTER (LEN=PREC_STRING), DIMENSION(MAXW) :: cpara   = ' '
integer                    , DIMENSION(MAXW) :: lpara   = 0
real(kind=PREC_DP)         , DIMENSION(MAXW) :: werte   = 0
integer :: gen
integer :: ianz
!
call get_params(line, ianz, cpara, lpara, maxw, length)
if(ier_num/=0) return
call ber_params(ianz, cpara, lpara, werte, maxw)
if(ier_num/=0) return
gen = nint(werte(1))
call omit_gen(gen)
!
end subroutine purge_files
!
!*******************************************************************************
!
SUBROUTINE set_gen(gen)
!
USE compare
USE population
USE errlist_mod
USE lib_do_operating_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: gen
!
INTEGER, PARAMETER :: IRD = 17
INTEGER, PARAMETER :: IWR = 18
INTEGER, PARAMETER :: lpar=LEN(pop_name)
!
CHARACTER(LEN=PREC_STRING)  :: line
CHARACTER(LEN=PREC_STRING)  :: logfiles
CHARACTER(LEN=PREC_STRING)  :: sumfiles
CHARACTER(LEN=PREC_STRING)  :: curfiles
CHARACTER(LEN=PREC_STRING)  :: outfile
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
!LOGICAL :: lexist
!
IF(gen==0) RETURN         ! No need to do anything
!
ios = 0
ier_num = 0
call check_gen_file(gen, lines_gen, generation, member, children, dimensions, lpar, &
                          logfiles, sumfiles, curfiles,  par_name, &
                          ubound(ier_msg, 1), ier_num, ier_typ, ier_msg, ER_APPL)
if(ier_num/=0 .or. ios/=0) return
!
!QQ INQUIRE(FILE='GENERATION',  EXIST=lexist)
!QQ cond_exist: IF(.NOT.lexist) THEN
!QQ    ier_num = -35          ! GENERATION DOES NOT EXIST
!QQ    ier_typ = ER_APPL
!QQ    RETURN
!QQ else cond_exist
!QQ OPEN(IRD, FILE='GENERATION', STATUS='OLD', ACCESS='SEQUENTIAL')
!QQ READ(IRD,'(a)', IOSTAT=ios) line 
!QQ IF(ios/=0) exit cond_exist !! Error reading GENERATION
!QQ READ(IRD,'(a)', IOSTAT=ios) line 
!QQ IF(ios/=0) exit cond_exist !! Error reading GENERATION
!QQ READ(line,*, IOSTAT=ios) generation, member, children, dimensions
!QQ IF(ios/=0) exit cond_exist !! Error reading GENERATION
!QQ IF(gen==generation) THEN
!QQ    CLOSE(IRD)             ! No need to do anything
!QQ    RETURN
!QQ ELSEIF(gen>GENERATION) THEN
!QQ    ier_num = - 36
!QQ    ier_typ = ER_APPL
!QQ    CLOSE(IRD)
!QQ    RETURN
!QQ ENDIF
!QQ READ(IRD,'(a)', IOSTAT=ios) line 
!QQ IF(ios/=0) exit cond_exist !! Error reading GENERATION
!QQ READ(IRD,'(a)', IOSTAT=ios) line 
!QQ IF(ios/=0) exit cond_exist !! Error reading GENERATION
!QQ READ(IRD,'(a)', IOSTAT=ios) line 
!QQ IF(ios/=0) exit cond_exist !! Error reading GENERATION
!QQ READ(IRD,'(a)', IOSTAT=ios) line 
!QQ IF(ios/=0) exit cond_exist !! Error reading GENERATION
!QQ !
!QQ READ(IRD,'(a)', IOSTAT=ios) line 
!QQ IF(ios/=0) exit cond_exist !! Error reading GENERATION
!QQ READ(IRD,'(a)', IOSTAT=ios) line 
!QQ IF(ios/=0) exit cond_exist !! Error reading GENERATION
!QQ logfiles = line(1:LEN_TRIM(line))
!QQ !
!QQ READ(IRD,'(a)', IOSTAT=ios) line 
!QQ IF(ios/=0) exit cond_exist !! Error reading GENERATION
!QQ READ(IRD,'(a)', IOSTAT=ios) line 
!QQ IF(ios/=0) exit cond_exist !! Error reading GENERATION
!QQ sumfiles = line(1:LEN_TRIM(line))
!QQ !
!QQ READ(IRD,'(a)', IOSTAT=ios) line 
!QQ IF(ios/=0) exit cond_exist !! Error reading GENERATION
!QQ READ(IRD,'(a)', IOSTAT=ios) line 
!QQ    IF(ios/=0) exit cond_exist !! Error reading GENERATION
!QQ curfiles = line(1:LEN_TRIM(line))
!QQ !
!QQ ALLOCATE(par_name(-1:dimensions))
!QQ lines_gen = 12
!QQ par_loop: DO
!QQ    READ(IRD,'(a)', IOSTAT=ios) line
!QQ    IF(ios/=0) exit cond_exist !! Error reading GENERATION
!QQ    lines_gen = lines_gen + 1
!QQ    IF(line(1:11)=='# Parameter') THEN           ! Got the list of parameter names
!QQ       DO i=-1, dimensions
!QQ          READ(IRD,'(a)', IOSTAT=ios) line
!QQ          IF(ios/=0) exit cond_exist ! Error reading GENERATION
!QQ          lines_gen = lines_gen + 1
!QQ          par_name(i) = line(1:MIN(lpar, LEN_TRIM(line)))
!QQ       ENDDO
!QQ       EXIT par_loop
!QQ    ENDIF
!QQ ENDDO par_loop
!QQ !
!QQ CLOSE(IRD)
!
! Create the new (shortened) summary files
!
IF(sumfiles/=' ') THEN
   nlines = 3 + gen
   sum_loop: DO i=0, dimensions
      WRITE(line,'(a,a,a,a,a,i10,a,a,a,a,a)') 'cat ',                   &
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
      WRITE(line,'(a,a,a,a,a,i10,a,a,a,a,a)') 'cat ',                   &
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
      WRITE(line,'(a,a,a,a,a,a,i10,a,a,a,a,a)') 'cat ',                 &
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
!QQ endif cond_exist
!9999 CONTINUE             ! Error reading GENERATION
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
subroutine check_gen_file(gen, lines_gen,                                 &
                          generation, member, children, dimensions, lpar, &
                          logfiles, sumfiles, curfiles,  par_name, &
                          NMSG, ier_num, ier_typ, ier_msg, ER_APPL)
!-
! read and check the generation file
!+
use precision_mod
!
implicit none
!
integer, intent(in) :: gen
integer, intent(out) :: lines_gen
integer, intent(out) :: generation
integer, intent(out) :: member
integer, intent(out) :: children
integer, intent(out) :: dimensions
integer, intent(in ) :: lpar
character(len=*), intent(out) :: logfiles
character(len=*), intent(out) :: sumfiles
character(len=*), intent(out) :: curfiles
character(len=*), dimension(:), allocatable, intent(out) :: par_name
integer, intent(in ) :: NMSG
integer, intent(out) :: ier_num
integer, intent(out) :: ier_typ
character(len=*), dimension(NMSG), intent(out) :: ier_msg
integer, intent(in ) :: ER_APPL
!
INTEGER, PARAMETER :: IRD = 17
!
character(len=PREC_STRING) :: line
integer :: ios
integer :: i
logical              :: lexist
!
lexist = .false.
ier_num = 0
ier_typ = 0
ios = 0
!
INQUIRE(FILE='GENERATION',  EXIST=lexist)
cond_exist: IF(.NOT.lexist) THEN
   ier_num = -35          ! GENERATION DOES NOT EXIST
   ier_typ = ER_APPL
   RETURN
else cond_exist
OPEN(IRD, FILE='GENERATION', STATUS='OLD', ACCESS='SEQUENTIAL')
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) exit cond_exist !! Error reading GENERATION
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) exit cond_exist !! Error reading GENERATION
READ(line,*, IOSTAT=ios) generation, member, children, dimensions
IF(ios/=0) exit cond_exist !! Error reading GENERATION
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
IF(ios/=0) exit cond_exist !! Error reading GENERATION
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) exit cond_exist !! Error reading GENERATION
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) exit cond_exist !! Error reading GENERATION
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) exit cond_exist !! Error reading GENERATION
!
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) exit cond_exist !! Error reading GENERATION
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) exit cond_exist !! Error reading GENERATION
logfiles = line(1:LEN_TRIM(line))
!
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) exit cond_exist !! Error reading GENERATION
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) exit cond_exist !! Error reading GENERATION
sumfiles = line(1:LEN_TRIM(line))
!
READ(IRD,'(a)', IOSTAT=ios) line 
IF(ios/=0) exit cond_exist !! Error reading GENERATION
READ(IRD,'(a)', IOSTAT=ios) line 
   IF(ios/=0) exit cond_exist !! Error reading GENERATION
curfiles = line(1:LEN_TRIM(line))
!
ALLOCATE(par_name(-1:dimensions))
lines_gen = 12
par_loop: DO
   READ(IRD,'(a)', IOSTAT=ios) line
   IF(ios/=0) exit cond_exist !! Error reading GENERATION
   lines_gen = lines_gen + 1
   IF(line(1:11)=='# Parameter') THEN           ! Got the list of parameter names
      DO i=-1, dimensions
         READ(IRD,'(a)', IOSTAT=ios) line
         IF(ios/=0) exit cond_exist ! Error reading GENERATION
         lines_gen = lines_gen + 1
         par_name(i) = line(1:MIN(lpar, LEN_TRIM(line)))
      ENDDO
      EXIT par_loop
   ENDIF
ENDDO par_loop
endif cond_exist
!
CLOSE(IRD)
!
if(ios/=0) then
   IF(ALLOCATED(par_name)) DEALLOCATE(par_name)
   ier_num = -19
   ier_typ = ER_APPL
   ier_msg(1) = 'Last line read successfully from GENERATION'
   ier_msg(2) = line(1:MIN(len_trim(line), LEN(ier_msg)))
endif
!
end subroutine check_gen_file
!
!*******************************************************************************
!
subroutine omit_gen(gen)
!-
!  Omit all generations previous to gen
!+
!
use population
!
use errlist_mod
use precision_mod
!
implicit none
!
integer, intent(in) :: gen
!
INTEGER, PARAMETER :: lpar=LEN(pop_name)
!
CHARACTER(LEN=PREC_STRING)  :: logfiles
CHARACTER(LEN=PREC_STRING)  :: sumfiles
CHARACTER(LEN=PREC_STRING)  :: curfiles
!CHARACTER(LEN=PREC_STRING)  :: outfile
CHARACTER(LEN=LEN(pop_name)), DIMENSION(:), ALLOCATABLE :: par_name
INTEGER :: generation
INTEGER :: member
INTEGER :: children
INTEGER :: dimensions
integer :: lines_gen
!
character(len=PREC_STRING) :: line
character(len=PREC_STRING) :: message
integer             :: exit_msg
integer             :: ier_cmd
integer :: ios
!
IF(gen==0) RETURN         ! No need to do anything
!
ier_num = 0
ios     = 0
call check_gen_file(gen, lines_gen, generation, member, children, dimensions, lpar, &
                          logfiles, sumfiles, curfiles,  par_name, &
                          ubound(ier_msg, 1), ier_num, ier_typ, ier_msg, ER_APPL)
if(ier_num/=0 .or. ios/=0) return
!
line = 'mkdir -p DIFFEV_NEW'
call execute_command_line(line, wait=.true., cmdstat=ier_cmd, cmdmsg=message, exitstat=exit_msg)
!
if(sumfiles/= ' ') then
   call omit_summary(gen, generation, dimensions, par_name)
endif
call omit_current  (gen, generation, member, dimensions, par_name)
call omit_parameter(gen, generation, member, dimensions, par_name)
call omit_generation(gen, generation, dimensions, par_name)
!
line = 'mv       DIFFEV_NEW/GENERATION .'
call execute_command_line(line, wait=.true., cmdstat=ier_cmd, cmdmsg=message, exitstat=exit_msg)
!
line = 'mv       DIFFEV_NEW/*.* DIFFEV/' 
call execute_command_line(line, wait=.true., cmdstat=ier_cmd, cmdmsg=message, exitstat=exit_msg)
!
line = 'rm -rf   DIFFEV_NEW' 
call execute_command_line(line, wait=.true., cmdstat=ier_cmd, cmdmsg=message, exitstat=exit_msg)
!
end subroutine omit_gen
!
!*******************************************************************************
!
subroutine omit_summary(gen, generation, dimensions, par_name)
!
use precision_mod
!
implicit none
!
integer, intent(in) :: gen
integer, intent(in) :: generation
integer, intent(in) :: dimensions
character(len=*), dimension(-1:dimensions), intent(in) :: par_name
!
INTEGER, PARAMETER :: IRD = 17
INTEGER, PARAMETER :: IWR = 18
!
character(len=PREC_STRING) ::  line
character(len=PREC_STRING) ::  infile
character(len=PREC_STRING) :: outfile
integer :: i, j
integer :: o_gen
real(kind=PREC_DP) :: p_ave
real(kind=PREC_DP) :: p_min
real(kind=PREC_DP) :: p_max
real(kind=PREC_DP) :: p_sig
!
do i=0, dimensions
   infile  = 'DIFFEV/Summary.' // par_name(i)(1:len_trim(par_name(i)))
   outfile = 'DIFFEV_NEW/Summary.' // par_name(i)(1:len_trim(par_name(i)))
   open(unit=IRD,file=infile, status='old')
   open(unit=IWR,file=outfile, status='unknown')
   do j= 1, 3                  ! Copy header lines
      read(IRD,'(a)') line
      write(IWR,'(a)') line(1:len_trim(line))
   enddo
   do j=0,gen-1                ! Omit lines 0 to gen-1
      read(IRD,'(a)') line
   enddo
   do j=gen, generation - 1    ! read remaining lines
      read(IRD,*) o_gen, p_ave, p_min, p_max, p_sig
      write(IWR, '(i7, 4G18.10E3)') o_gen-gen, p_ave, p_min, p_max, p_sig
   enddo
   close(IRD)
   close(IWR)
enddo
!
end subroutine omit_summary
!
!*******************************************************************************
!
subroutine omit_current(gen, generation, member, dimensions, par_name)
!
use precision_mod
!
implicit none
!
integer, intent(in) :: gen
integer, intent(in) :: generation
integer, intent(in) :: member
integer, intent(in) :: dimensions
character(len=*), dimension(-1:dimensions), intent(in) :: par_name
!
INTEGER, PARAMETER :: IRD = 17
INTEGER, PARAMETER :: IWR = 18
!
character(len=PREC_STRING) ::  line
character(len=PREC_STRING) ::  infile
character(len=PREC_STRING) :: outfile
integer :: i, j
!integer :: o_gen
!real(kind=PREC_DP) :: p_ave
!real(kind=PREC_DP) :: p_min
!real(kind=PREC_DP) :: p_max
!real(kind=PREC_DP) :: p_sig
!
do i=0, dimensions
   infile  = 'DIFFEV/Current.' // par_name(i)(1:len_trim(par_name(i)))
   outfile = 'DIFFEV_NEW/Current.' // par_name(i)(1:len_trim(par_name(i)))
   open(unit=IRD,file=infile, status='old')
   open(unit=IWR,file=outfile, status='unknown')
!  Copy three header lines
   read(IRD,'(a)') line
   write(IWR,'(a)') line(1:len_trim(line))
   read(IRD,'(a)') line
   write(line(4:10), '(i6)') generation - 1 - gen
   write(IWR,'(a)') line(1:len_trim(line))
   read(IRD,'(a)') line
   write(IWR,'(a)') line(1:len_trim(line))
   do j=1,member               ! Copy remaining lines
      read(IRD,'(a)') line
      write(IWR,'(a)') line(1:len_trim(line))
   enddo
   close(IRD)
   close(IWR)
enddo
!
end subroutine omit_current
!
!*******************************************************************************
!
subroutine omit_parameter(gen, generation, member, dimensions, par_name)
!
use precision_mod
!
implicit none
!
integer, intent(in) :: gen
integer, intent(in) :: generation
integer, intent(in) :: member
integer, intent(in) :: dimensions
character(len=*), dimension(-1:dimensions), intent(in) :: par_name
!
INTEGER, PARAMETER :: IRD = 17
INTEGER, PARAMETER :: IWR = 18
!
character(len=PREC_STRING) ::  line
character(len=PREC_STRING) ::  infile
character(len=PREC_STRING) :: outfile
integer :: i, j, k, igen
!integer :: o_gen
!real(kind=PREC_DP) :: p_ave
!real(kind=PREC_DP) :: p_min
!real(kind=PREC_DP) :: p_max
!real(kind=PREC_DP) :: p_sig
!
do i=0, dimensions
   infile  = 'DIFFEV/Parameter.' // par_name(i)(1:len_trim(par_name(i)))
   outfile = 'DIFFEV_NEW/Parameter.' // par_name(i)(1:len_trim(par_name(i)))
   open(unit=IRD,file=infile, status='old')
   open(unit=IWR,file=outfile, status='unknown')
!
   read(IRD,'(a)') line
   write(IWR,'(a)') line(1:len_trim(line))
   do j=1, gen*(member+2)     ! Skip "old" generations
     read(IRD,'(a)') line
   enddo
   igen = 0
   do j=gen, generation - 1   ! Copy remaining generations
      igen = igen + 1
      read(IRD,'(a)') line
      write(line(4:10), '(i6)') igen !generation - 1 - gen
      write(IWR,'(a)') line(1:len_trim(line))
      do k=0,member               ! Copy remaining lines
         read(IRD,'(a)') line
         write(IWR,'(a)') line(1:len_trim(line))
      enddo
   enddo
   close(IRD)
   close(IWR)
enddo
!
end subroutine omit_parameter
!
!*******************************************************************************
!
subroutine omit_generation(gen, generation, dimensions, par_name)
!
use precision_mod
!
implicit none
!
integer, intent(in) :: gen
integer, intent(in) :: generation
integer, intent(in) :: dimensions
character(len=*), dimension(-1:dimensions), intent(in) :: par_name
!
INTEGER, PARAMETER :: IRD = 17
INTEGER, PARAMETER :: IWR = 18
!
character(len=PREC_STRING) ::  line
character(len=PREC_STRING) ::  infile
character(len=PREC_STRING) :: outfile
integer :: i, j
!integer :: o_gen
!real(kind=PREC_DP) :: p_ave
!real(kind=PREC_DP) :: p_min
!real(kind=PREC_DP) :: p_max
!real(kind=PREC_DP) :: p_sig
!
do i=0, dimensions
   infile  = 'GENERATION'
   outfile = 'DIFFEV_NEW/GENERATION'
   open(unit=IRD,file=infile, status='old')
   open(unit=IWR,file=outfile, status='unknown')
!
   read(IRD,'(a)') line
   write(IWR,'(a)') line(1:len_trim(line))
   read(IRD,'(a)') line
   write(line(4:10), '(i6)') generation - gen
   write(IWR,'(a)') line(1:len_trim(line))
   do j=1, 29 + dimensions+2  ! Copy remaining lines
      read(IRD,'(a)') line
      write(IWR,'(a)') line(1:len_trim(line))
   enddo
   close(IRD)
   close(IWR)
enddo
!
end subroutine omit_generation
!
!*******************************************************************************
!
END MODULE diffev_set_gen_mod
