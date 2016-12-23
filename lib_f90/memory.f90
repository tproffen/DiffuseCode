!
subroutine memory_message(cpara)
!
use mpi_slave_mod
IMPLICIT NONE
character(Len=*), intent(IN) :: cpara
character(len=1024)          :: zei
!DBG_MEM
!DBG_MEM
character(len=1024):: pid_status
character(len=8   ):: pid_char
integer::my_pid
integer::ios
integer::vmpeak = 0
integer::vmsize = 0
integer::vmlck = 0
integer::vmpin = 0
integer::vmhwm = 0
integer::vmrss = 0
integer::vmdata = 0
integer::vmstk = 0
integer::vmexe = 0
integer::vmlib = 0
integer::vmpte = 0
integer::vmswap = 0
my_pid=getpid()
write(pid_char,'(I8)') my_pid
pid_status = '/proc/'//TRIM(ADJUSTL(pid_char))//'/status'
open(unit=100,file=pid_status,action='read')
memory: do
   read(100,'(a)',IOSTAT=ios) zei
   if(ios/=0) exit memory
   if(zei(1:7)=='VmPeak:') read(zei(8:),*) vmpeak
   if(zei(1:7)=='VmSize:') read(zei(8:),*) vmsize
   if(zei(1:6)=='VmLck:' ) read(zei(7:),*) vmlck 
   if(zei(1:6)=='VmPin:' ) read(zei(7:),*) vmpin 
   if(zei(1:6)=='VmHWM:' ) read(zei(7:),*) vmhwm 
   if(zei(1:6)=='VmRSS:' ) read(zei(7:),*) vmrss 
   if(zei(1:7)=='VmDATA:') read(zei(8:),*) vmdata
   if(zei(1:6)=='VmStk:' ) read(zei(7:),*) vmstk 
   if(zei(1:6)=='VmExe:' ) read(zei(7:),*) vmexe 
   if(zei(1:6)=='VmLib:' ) read(zei(7:),*) vmlib 
   if(zei(1:6)=='VmPTE:' ) read(zei(7:),*) vmpte 
   if(zei(1:7)=='VmSwap:') read(zei(8:),*) vmswap
enddo memory
close(100)
!write(101,'(a20,12(1x,I10))') cpara(1:20),vmpeak,vmsize,vmlck, vmpin,vmhwm,vmrss,&
!                             vmdata,vmstk,vmexe,vmlib,vmpte,vmswap
write(101,'(a20, 5(1x,I10))') cpara(1:20),vmpeak,vmsize, vmhwm,vmrss, vmpte
end subroutine memory_message
