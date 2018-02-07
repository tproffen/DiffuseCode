	program arctan
c
	implicit	none
c
	integer		generation
	integer		member
	integer		children
	integer		npar
	integer		i,ii
c
	real		r(3)
	real		exper(2,2001)
	real		x,y
	real		rval
	real		sumrz,sumrn
c
	character*18	trials
	character*18	result
c
	open(7,file='GENERATION',status='old')
	read(7,*)
	read(7,*) generation,member,children,npar
	close(7)
c
	open(7,file='DATA/data.noisy',status='old')
	read(7,*)
	read(7,*)
	do i=1,2001
	  read(7,*) exper(1,i),exper(2,i)
	enddo
	close(7)
c
c
	do i=1,children
	  write(trials,1000) i
	  write(result,1100) i
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
	  do ii=1,2001
	    x = exper(1,ii)
	    y = r(1)*atan(abs(x-r(2))/r(3))
	    sumrz = sumrz + (exper(2,ii)-y)**2
	    sumrn = sumrn + (exper(2,ii)  )**2
	  enddo
	  rval = sqrt(sumrz/sumrn)
	  open(8,file=result,status='unknown')
	  write(8,2000) i,rval
	  close(8)
	enddo
c
1000	format('DIFFEV/Trials.',i4)
1100	format('DIFFEV/Result.',i4)
2000	format(i4,f12.8)
	end
