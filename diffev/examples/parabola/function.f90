      program function
!
      implicit      none
!
      integer            generation
      integer            member
      integer            children
      integer            npar
      integer            i,ii,j
      integer            npoints
!
      real            r(50)
      real            exper(2,2001)
      real            calc (2,2001)
      real            x,y
      real            rval
      real            sumrz,sumrn
!
      character*19      trials
      character*19      result
      character*19      berech
!
      open(7,file='GENERATION',status='old')
      read(7,*)
      read(7,*) generation,member,children,npar
      close(7)
!
      open(7,file='DATA/function.data',status='old')
      read(7,*)
      read(7,*)
      i = 0
      do 
        read(7,*,end=7777) exper(1,i+1),exper(2,i+1)
        i = i+1
      enddo
7777  continue
      npoints = i
      close(7)
!
!
      do i=1,children
        write(trials,1000) i
        write(result,1100) i
        write(berech,1200) i
        do ii=15,18
          if(trials(ii:ii).eq.' ') trials(ii:ii) = '0'
          if(result(ii:ii).eq.' ') result(ii:ii) = '0'
        enddo
        open(7,file=trials,status='old')
        read(7,*)
        read(7,*)
        read(7,*)
        read(7,*)
        read(7,*)
        do ii=1,npar
          read(7,*) r(ii)
        enddo
        close(7)
        sumrz = 0.0
        sumrn = 0.0
        open(8,file=berech,status='unknown')
        write(8,3000)
        do ii=1, npoints
          x = exper(1,ii)
          y = r(1)
          do j=2,npar
             y = y + r(j)*x**(j-1)
          enddo
          WRITE(8,3100) x,y
          sumrz = sumrz + (exper(2,ii)-y)**2
          sumrn = sumrn + (exper(2,ii)  )**2
        enddo
        close(8)
        rval = sqrt(sumrz/sumrn)
        open(8,file=result,status='unknown')
        write(8,2000) i,rval
        close(8)
      enddo
!
1000      format('DIFFEV/Trials.',i4)
1100      format('DIFFEV/Results.',i4)
1200      format('TEMP/calc.',i4.4)
2000      format(i4,f12.8)
3000      format('#',/,'#')
3100      format(2(G20.8,2x))
      end
