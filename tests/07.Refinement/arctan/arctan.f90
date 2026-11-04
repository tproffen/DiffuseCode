program arctan
!
implicit  none
!
INTEGER :: i,ii
!
real      :: r(3)
real      :: exper(2,2001)
real      :: x,y
real      :: rval
real      :: sumrz,sumrn
INTEGER :: kid
INTEGER :: iarg
INTEGER :: ex_stat
!
CHARACTER(LEN=18)  :: result
CHARACTER(LEN=128), DIMENSION(4) ::  arg
!
iarg = COMMAND_ARGUMENT_COUNT()
IF(iarg == 4) THEN
   DO i=1,4
      CALL GET_COMMAND_ARGUMENT(i,arg(i))
   ENDDO
   READ(arg(1),*) kid
   READ(arg(2),*) r(1)
   READ(arg(3),*) r(2)
   READ(arg(4),*) r(3)
ENDIF
!
open(7,file='DATA/data.noisy',status='old')
read(7,*)
read(7,*)
do i=1,2001
   read(7,*) exper(1,i),exper(2,i)
enddo
close(7)
!
sumrz = 0.0
sumrn = 0.0
do ii=1,2001
  x = exper(1,ii)
  y = r(1)*atan(abs(x-r(2))/r(3))
  sumrz = sumrz + (exper(2,ii)-y)**2
  sumrn = sumrn + (exper(2,ii)  )**2
enddo
rval = sqrt(sumrz/sumrn)
WRITE(result,1100) kid
4000 FORMAT('TEMP/result.',i4.4)
!
open(8,file=result,status='unknown')
write(8,2000) kid,rval
close(8)
ex_stat = NINT(1000*rval)
write(*,*) ' ex_stat ', ex_stat
CALL exit(ex_stat)
!
1100  format('DIFFEV/Result.',i4.4)
2000  format(i4,f12.8)
end
