module discamb_mod
!-
!  structure for the DISCAMB generated atomic form factors
!
!  discamb_alloc_list        ! Allocate the basic array
!  discamb_find_nsymm        ! Determine number symmetry operations
!  discamb_read              ! Read the DISCAMB '*.tsc' file
!  discamb_set_form_tsc      ! Copy inital read into global storage
!+
!
use precision_mod
!
private
!
type :: disc_form
!  private
   integer                                                :: iscat    ! The scattering type
   integer              , dimension(:), allocatable       :: isymm    ! List of symmetry operations
!  integer              , dimension(:), allocatable       :: use_symm      ! Use this symmetry
   complex(kind=PREC_DP), dimension(:,:,:,:), allocatable :: four_form_tsc ! Form factors from DISCAMB
end type disc_form
!
type(disc_form), dimension(:), pointer :: disc_list => NULL()
integer, dimension(:,:), allocatable   :: disc_nsymm          ! Number of symmetry operations for each scattering type
!
complex(kind=PREC_DP), dimension(:,:,:,:), allocatable :: four_form_tsc ! Form factors from DISCAMB as read
integer, dimension(3), save :: maxhkl   ! Maximum abs(hkl)
!
public discamb_read
public disc_list
!
contains
!
!*******************************************************************************
!
subroutine discamb_alloc_list(nscat)
!-
! Allocate the basic array 
!+
implicit none
!
integer, intent(in) :: nscat
integer :: i
!
if(associated(disc_list)) then
   do i=1, ubound(disc_list,1)
     if(allocated(disc_list(i)%isymm))         deallocate(disc_list(i)%isymm)
     if(allocated(disc_list(i)%four_form_tsc)) deallocate(disc_list(i)%four_form_tsc)
   enddo
   deallocate(disc_list)
endif
allocate(disc_list(nscat))
!
end subroutine discamb_alloc_list
!
!*******************************************************************************
!
subroutine discamb_find_nsymm
!-
!  Determine the number of symmetry operations for each atom type
!
!  At the moment, all symmetry operations for a given atom type are treated
!  individually. A reduction could be achieved, if the local Wyckoff site 
!  symmetry is determined and a dummy tensor for this point group symmetry is
!  transformed by all symmetry operations.
!+
!
use crystal_mod
use wyckoff_mod
!
implicit none
!
!real(kind=PREC_DP), parameter :: TOL = 1.0D-3
!
integer :: i,j!,k   ! Dummy loop index
integer :: nsymm
!real(kind=PREC_DP), dimension(3,3), parameter :: base = reshape((/.1,.2,.3,  .2, .4, .5, .3, .5, .6/), shape(base))
!real(kind=PREC_DP), dimension(:,:,:), allocatable :: tensors
!real(kind=PREC_DP), dimension(3,3)            :: temp 
!real(kind=PREC_DP), dimension(3,3)            :: symm_mat
!
allocate(disc_nsymm(0:cr_nscat, 0:192))
disc_nsymm = 0
!
loop_atoms: do i=1, cr_natoms
   loop_symm: do j=1, disc_nsymm(cr_iscat(1,i),0)
      if(disc_nsymm(cr_iscat(1,i),j)==cr_iscat(2,i)) cycle loop_atoms
   enddo loop_symm
   disc_nsymm(cr_iscat(1,i),0) = disc_nsymm(cr_iscat(1,i),0) + 1
   disc_nsymm(cr_iscat(1,i),disc_nsymm(cr_iscat(1,i),0)) = cr_iscat(2,i)
enddo loop_atoms
!do i=0,cr_nscat
!  write(*,*) i,' | ', disc_nsymm(i,0), ' : ', disc_nsymm(i,1:disc_nsymm(i,0))
!enddo
!
do i=1, cr_nscat
   if(allocated(disc_list(i)%isymm)) deallocate(disc_list(i)%isymm)
!  if(allocated(disc_list(i)%use_symm)) deallocate(disc_list(i)%use_symm)
   nsymm = disc_nsymm(i,0)
   allocate(disc_list(i)%isymm(0:nsymm))
!  allocate(disc_list(i)%use_symm(0:nsymm))
   disc_list(i)%isymm(0:nsymm) = disc_nsymm(i,0:nsymm)
   disc_list(i)%iscat = i
!  allocate(tensors(3,3,nsymm))
!  tensors = 0.0_PREC_DP
!  tensors(:,:,1) = base
!  disc_list(i)%use_symm(1) = 1    ! Identity is always reproduced by identity
!  loop_symmt: do j=2, nsymm       ! Test if symmetry operation changes the base
!     symm_mat = spc_mat(1:3,1:3, disc_list(i)%isymm(j))
!     temp = matmul(symm_mat, matmul(base, transpose(symm_mat)))
!if(i==3) then
!write(*,'(4(2x,3f7.3))') temp(1,:), symm_mat(1,:), base(1,:)!, transpose(sym_mat)(:,1)
!write(*,'(4(2x,3f7.3))') temp(2,:), symm_mat(2,:), base(2,:)!, transpose(sym_mat)(:,2)
!write(*,'(4(2x,3f7.3))') temp(3,:), symm_mat(3,:), base(3,:)!, transpose(sym_mat)(:,3)
!endif
!     do k=1, j-1                 ! Compare to previous tensors
!        if(all( abs(temp-tensors(:,:,k))<TOL)) then
!           disc_list(i)%use_symm(j) = k
!           cycle loop_symmt
!        endif
!     enddo
!     tensors(:,:,j) = temp
!     disc_list(i)%use_symm(j) = j    ! This is a new symmetry result
!  enddo loop_symmt
!  deallocate(tensors)
!  write(*,'(a,i3,a,i3,a,10i3)') 'Iscat ', i, ' SYMM ', disc_list(i)%isymm(0), ' | ', disc_list(i)%isymm(1:disc_list(i)%isymm(0))
!  write(*,'(a,i3,a,i3,a,10i3)') 'Iscat ', i, ' use  ', disc_list(i)%isymm(0), ' | ', disc_list(i)%use_symm(1:disc_list(i)%isymm(0))
enddo
!
deallocate(disc_nsymm)
!
end subroutine discamb_find_nsymm
!
!**********************************************************************
!
subroutine discamb_read(tscfile, trust)
!-
!  Read the user provided file with atomic form factors from DISCAMB
!+
!
use crystal_mod
use diffuse_mod
!
use blanks_mod
use errlist_mod
use get_params_mod
use matrix_mod
use precision_mod
use support_mod
!
implicit none
!
character(len=*), intent(in) :: tscfile    ! Input file
character(len=*), intent(in) :: trust      ! Input file
!character(len=PREC_STRING)   :: tscfile    ! Input file
!
integer, parameter :: IRD=33
!real(kind=PREC_DP), parameter :: TOL = 0.01_PREC_DP
!
character(len=PREC_STRING) :: line
character(len=PREC_STRING) :: string
integer :: ianz
character(len=PREC_STRING), dimension(:), allocatable :: cpara
integer                   , dimension(:), allocatable :: lpara
integer :: DIM_CPARA        ! Array size
integer :: ios              ! IO signal
integer :: i, ii,jj,kk,l    ! Dummy index
integer :: length           ! string lengths
integer :: nr               ! Line number
integer :: nhdr             ! header Line number
integer :: ntotal           ! total lines in input file
integer :: num_scatterers   ! Number of different atoms types in tscfile
integer, dimension(3) ::   ihkl   !            (hkl)
logical :: ldata   ! got DATA: line

real(kind=PREC_DP), dimension(3)   :: hkl
real(kind=PREC_DP), dimension(3,3) :: vi_user! User copy  vi matrix
real(kind=PREC_DP), dimension(3,3) :: viinv  ! inverse of vi matrix
real(kind=PREC_DP)               :: a,b    ! real, imag part
!complex(KIND=PREC_DP), dimension(:,:,:,:), allocatable :: four_form_tsc
!
if(associated(disc_list) .and. trust=='use') then 
   return                  ! We can trust the old tsc file
endif
call discamb_alloc_list(cr_nscat)
call discamb_find_nsymm
!
call oeffne(IRD, tscfile, 'old')
if(ier_num /=0) then
   ier_msg(1) = 'DISCAMB file does not exist'
   ier_msg(2) = 'FILE: '//tscfile(1:len_trim(tscfile))
   return
endif
!
nr = 0
maxhkl = 0
ldata  = .false.
loop_init: do
   read(IRD, '(a)', iostat=ios) line
   if(is_iostat_end(ios)) exit loop_init
   nr = nr + 1
   if(line(1:10)=='SCATTERERS') then  ! We got the list of atom names
      length=len_trim(line)
      jj = 0
      ii = 13                         ! First non character after :
      loop_types: do i = 13, length+1   ! Count all SCATTERERS, compare to atom names
         if(line(i:i)==' ') then
            jj = jj + 1
            if(line(ii:i-1)/= cr_at_lis(jj)) then
               ier_num = -201
               ier_typ = ER_APPL
               string = 'Crystal:'
               do l=1, cr_nscat
                  string(len_trim(string)+1:len_trim(string)+6) = ' ' // cr_at_lis(l)
               enddo
               l = len_trim(string)
               ier_msg(1) = string(1:min(len(ier_msg),l))
               write(ier_msg(2), '(2a)') 'TSC    : ', line(13:length)
               close(IRD)
               return
            else
               ii = i+1
               cycle loop_types
            endif
         endif
      enddo loop_types
      if(jj /= cr_nscat) then          ! Number of scatterers differs from crystal structure
         ier_num = -202
         ier_typ = ER_APPL
         write(ier_msg(1), '(a, i5)') 'Crystal structure ', cr_nscat
         write(ier_msg(2), '(a, i5)') 'DISCAMB TSC file  ', jj
         close(IRD)
         return
      endif
   endif
   if(line(1:5)=='DATA:') then
      ldata = .true.
      cycle loop_init
   endif
   if(ldata) then        ! We got the data line, read hkls
      read(line,*) ihkl
      maxhkl(1) = max(maxhkl(1),abs(ihkl(1)))
      maxhkl(2) = max(maxhkl(2),abs(ihkl(2)))
      maxhkl(3) = max(maxhkl(3),abs(ihkl(3)))
   endif
enddo loop_init
ntotal = nr
close(IRD) 
!
if(any(maxhkl<int(diff_maxhkl))) then   ! User requested hkl limis outside tcs file!
   ier_num = -200
   ier_typ = ER_APPL
   write(ier_msg(1),'(a,3i4)') 'Fourier  limits ', int(diff_maxhkl)
   write(ier_msg(2),'(a,3i4)') 'TSC file limits ', maxhkl
   return
endif
!
call oeffne(IRD, tscfile, 'old')
if(ier_num /=0) then
   ier_msg(1) = 'DISCAMB file does not exist'
   ier_msg(2) = 'FILE: '//tscfile(1:len_trim(tscfile))
   return
endif
!
vi_user = vi    ! Create backup
if(inc(3)==1) vi(3,3) = 1.0
if(inc(2)==1) vi(2,2) = 1.0
if(inc(1)==1) vi(1,1) = 1.0
!
!
call matinv(vi, viinv)
if(ier_num/=0) then
   ier_msg(1) = 'Matrix of increment vectors, all parallel?'
   vi = vi_user     ! use  backup
   return
endif
nr=0
loop_header: do
   nr = nr + 1
   read(IRD, '(a)', iostat=ios) line
   if(ios/=0) then
      write(ier_msg(1),'(a,i9)') 'Error in DISCAMB file, line: ', nr
      ier_num = -3
      ier_typ = ER_IO
      close(IRD)
   vi = vi_user     ! use  backup
      return
   endif
   if(index(line, 'SCATTERERS:') > 0) then   ! Found "SCATTERERS" line
      string = ' '//line(12:len_trim(line))      ! Extract substring with atom names
      length = len_trim(string)
      call rem_dbl_bl(string,length)        ! Ensure that we have no multiple blanks
      num_scatterers = 0
      do i=1,length
         if(string(i:i) == ' ') num_scatterers = num_scatterers + 1
      enddo
      DIM_CPARA = num_scatterers+3
      allocate(cpara(DIM_CPARA)) 
      allocate(lpara(DIM_CPARA)) 
      call get_params_blank(string, ianz, cpara, lpara, DIM_CPARA, length)
      exit loop_header
   endif
enddo loop_header
! Analyse that TSC names match 'cr_at_lis' WORK
if(allocated(four_form_tsc)) deallocate(four_form_tsc)
allocate(four_form_tsc(-maxhkl(1):maxhkl(1), -maxhkl(2):maxhkl(2), -maxhkl(3):maxhkl(3), num_scatterers))
!
!  Read remaining header
!
loop_data: do                         ! Read remaining header until 'DATA:'
   nr = nr + 1
   read(IRD, '(a)', iostat=ios) line
   if(ios/=0) then
      write(ier_msg(1),'(a,i9)') 'Error in DISCAMB file, line: ', nr
      ier_num = -3
      ier_typ = ER_IO
      close(IRD)
   vi = vi_user     ! use  backup
      return
   endif
   if(index(line, 'DATA:') > 0) then   ! Found "DATA" line
      nhdr = nr                        ! Store header line number
      exit loop_data
   endif
enddo loop_data
!
l=0
loop_values: do                        ! Read actual values
   nr = nr + 1
   read(IRD, '(a)', iostat=ios) line
   if(is_iostat_end(ios)) then
      exit loop_values
   elseif(ios/=0) then
      write(ier_msg(1),'(a,i9)') 'Error in DISCAMB file, line: ', nr
      ier_num = -3
      ier_typ = ER_IO
      close(IRD)
      return
   endif
   length= len_trim(line)
   call get_params_blank(line, ianz, cpara, lpara, DIM_CPARA, length)
   do i=1, 3
      read(cpara(i),*) hkl(i)
   enddo
   ii = nint(hkl(1))
   jj = nint(hkl(2))
   kk = nint(hkl(3))
   do i=1, num_scatterers
      read(cpara(i+3), *) a,b
      four_form_tsc(ii,jj,kk, i) = cmplx(a,b, kind=PREC_DP)
   enddo
!if(abs(ii)==1 .and. abs(jj)==1 .and. abs(kk)==1) then
!write(*,'(3i4, 2f10.6)') ii,jj,kk, real(four_form_tsc(ii,jj,kk, 1)), imag(four_form_tsc(ii,jj,kk, 1))
!endif
!  hkl = hkl - eck(:,1)            ! Subtract lower left bottom corner
!  uvw = matmul(viinv, hkl)
!  if(all(uvw-nint(uvw)<TOL)) then   ! All indices are integer
!     ii = nint(uvw(1)) + 1
!     jj = nint(uvw(2)) + 1
!     kk = nint(uvw(3)) + 1
!     if(ii>0 .and. ii<=num(1) .and. jj>0 .and. jj<=num(2) .and. kk>0 .and. kk<=num(3)) then ! Within range
!        l = l+1
!        do i=1, num_scatterers
!           read(cpara(i+3), *) a,b
!           four_form_tsc(ii,jj,kk, i) = cmplx(a,b, kind=PREC_DP)
!        enddo
!     endif
!  endif
enddo loop_values
four_form_tsc(0,0,0, :) = cmplx(1., 1.)
!
close(unit=IRD)
!do kk=-maxhkl(3), maxhkl(3)
!do jj=-maxhkl(2), maxhkl(2)
!do ii=-maxhkl(1), maxhkl(1)
!do i=1, num_scatterers
!if( abs(real(four_form_tsc(ii,jj,kk,i)))<0.000001 .or. &
!    abs(imag(four_form_tsc(ii,jj,kk,i)))<0.000001     ) then
!write(*,'(3i4, 2f9.6)') ii,jj,kk, four_form_tsc(ii,jj,kk,i)
!endif
!enddo
!enddo
!enddo
!enddo
!write(*,*) ' ZEROS ', any(abs(real(four_form_tsc(:,:,:,:)))<0.000001)
!write(*,*) ' ZEROS ', any(abs(imag(four_form_tsc(:,:,:,:)))<0.000001)
!
call discamb_set_form_tsc(cr_nscat)
!write(*,*) ' FOUND SCATTERERS ', num_scatterers
!write(*,*) ' FOUND POINTS     ', nr -1 - nhdr, l
!write(*,'(a,8(''('',f12.6,'','',f12.6,'')''))') ' LAST ', four_form_tsc(num(1),num(2),num(3), 1:min(8,num_scatterers))

i= 1  ! ATOM TYPE 1
!do j=1, disc_list(i)%isymm(0)
!  write(*,'(i4, 4(2x,2f10.6))') j, &
!real(disc_list(i)%four_form_tsc( 1, 1,1   , j)), &
!imag(disc_list(i)%four_form_tsc( 1, 1,1   , j)), &
!real(disc_list(i)%four_form_tsc(-1, 1,1   , j)), &
!imag(disc_list(i)%four_form_tsc(-1, 1,1   , j)), &
!real(disc_list(i)%four_form_tsc( 1,-1,1   , j)), &
!imag(disc_list(i)%four_form_tsc( 1,-1,1   , j)), &
!real(disc_list(i)%four_form_tsc(-1,-1,1   , j)), &
!imag(disc_list(i)%four_form_tsc(-1,-1,1   , j))
!enddo
!write(*,*)

!do j=1, disc_list(i)%isymm(0)
!  write(*,'(i4, 4(2x,2f10.6))') j, &
!real(disc_list(i)%four_form_tsc( 1, 1,-1  , j)), &
!imag(disc_list(i)%four_form_tsc( 1, 1,-1  , j)), &
!real(disc_list(i)%four_form_tsc(-1, 1,-1  , j)), &
!imag(disc_list(i)%four_form_tsc(-1, 1,-1  , j)), &
!real(disc_list(i)%four_form_tsc( 1,-1,-1  , j)), &
!imag(disc_list(i)%four_form_tsc( 1,-1,-1  , j)), &
!real(disc_list(i)%four_form_tsc(-1,-1,-1  , j)), &
!imag(disc_list(i)%four_form_tsc(-1,-1,-1  , j))
!enddo
!write(*,*) ' SEARCH FOR ZEROS '
!do kk=-maxhkl(3), maxhkl(3)
!do jj=-maxhkl(2), maxhkl(2)
!do ii=-maxhkl(1), maxhkl(1)
!do i=1, num_scatterers
!do j=1,4
!if( abs(real(disc_list(i)%four_form_tsc(ii,jj,kk,j)))<0.000001 .or. &
!    abs(imag(disc_list(i)%four_form_tsc(ii,jj,kk,j)))<0.000001     ) then
!!write(*,'(2i6,3i4, 2f9.6)') i, j, ii,jj,kk, disc_list(i)%four_form_tsc(ii,jj,kk,j)
!endif
!enddo
!enddo
!enddo
!enddo
!enddo
!write(*,*) ' ZEROS ',  any(abs(real(disc_list(1)%four_form_tsc(:,:,:,1:4)))<0.000001)
!write(*,*) ' ZEROS ',  any(abs(imag(disc_list(1)%four_form_tsc(:,:,:,1:4)))<0.000001)
deallocate(cpara)
deallocate(lpara)
vi = vi_user     ! use  backup
!deallocate(four_form_tsc)
cr_is_anis = .true.
!read(*,*) nr
!
end subroutine discamb_read
!
!*******************************************************************************
!
!
subroutine discamb_set_form_tsc(nscat)
!-
!  Copy the initial atomic form factors into the full table, apply symmetry
!+
use crystal_mod
use diffuse_mod
use wyckoff_mod
!
use matrix_mod
!
integer, intent(in) :: nscat    ! Number of scattering types
!
real(kind=PREC_DP), parameter :: TOL = 0.01_PREC_DP!
!
integer :: i,j   ! Dummy loop index
integer :: h,k,l   ! Dummy loop index
integer :: hl,kl,ll ! Lower bound
integer :: hh,kh,lh ! upper bound
integer :: ii,jj,kk ! Dummy indices
real(kind=PREC_DP), dimension(3)   :: veca, uvw    ! vectors
real(kind=PREC_DP), dimension(3,3) :: sym   ! Symmetry matrix
!real(kind=PREC_DP), dimension(3,3) :: viinv ! Inverse matrix to vsteps reciprocal space
!
!call matinv(vi, viinv)
!
hl = lbound(four_form_tsc,1)
kl = lbound(four_form_tsc,2)
ll = lbound(four_form_tsc,3)
hh = ubound(four_form_tsc,1)
kh = ubound(four_form_tsc,2)
lh = ubound(four_form_tsc,3)
loop_scat: do i=1, nscat
   allocate(disc_list(i)%four_form_tsc(hl:hh, kl:kh, ll:lh, disc_list(i)%isymm(0)))
   disc_list(i)%four_form_tsc(:,:,:,1) = four_form_tsc(:,:,:,i)
!write(*,*)
!write(*,'(a,2i3, 8f10.6)') ' iscat',i, 1, disc_list(i)%four_form_tsc(6,6,6,1), &
!                                          disc_list(i)%four_form_tsc(8,6,6,1), &
!                                          disc_list(i)%four_form_tsc(6,8,6,1), &
!                                          disc_list(i)%four_form_tsc(8,8,6,1)
!write(*,'(a,2i3, 8f10.6)') ' iscat',i, 1, disc_list(i)%four_form_tsc(6,6,8,1), &
!                                          disc_list(i)%four_form_tsc(8,6,8,1), &
!                                          disc_list(i)%four_form_tsc(6,8,8,1), &
!                                          disc_list(i)%four_form_tsc(8,8,8,1)
!
   loop_sym: do j=2, disc_list(i)%isymm(0)     ! Loop over all symmetry operations 
      sym = spc_mat(1:3,1:3, disc_list(i)%isymm(j))
!if(i==1) then   ! Atom type 1
!write(*,'(a, 3f5.2)') ' symm(1,:) ', sym(1,:)
!write(*,'(a, 3f5.2)') ' symm(2,:) ', sym(2,:)
!write(*,'(a, 3f5.2)') ' symm(3,:) ', sym(3,:)
!endif
      do l=ll, lh                    ! Loop over all points in reciprocal space
         veca(3) = l
         do k=kl, kh    
            veca(2) = k
            do h=hl, hh    
               veca(1) = h
               uvw = matmul(veca, sym)
!if(i==1) then   ! Atom type 1
!if(abs(h)==11 .and. abs(k)==13 .and. abs(l)==14) then
!write(*,'(3i4, 3i4,l3)') h,k,l, nint(uvw), all(abs(uvw-nint(uvw))<TOL)
!endif
!endif
               if(all(abs(uvw-nint(uvw))<TOL)) then   ! All indices are integer
                  ii = nint(uvw(1))
                  jj = nint(uvw(2))
                  kk = nint(uvw(3))
                  if(ii>=hl .and. ii<=hh .and. jj>=kl.and. jj<=kh .and. kk>=ll.and. kk<=lh) then ! Within range
!if(h==6 .and. k==6 .and. l==6) then
!write(*,'(a,3i3,3i3,3(2x,3f6.2))') 'VALID ', h,k,l, ii,jj,kk, veca, vecb, uvw
!endif

                     disc_list(i)%four_form_tsc(h,k,l   , j) = disc_list(i)%four_form_tsc(ii,jj,kk,1)
                  endif
               endif
            enddo
         enddo
      enddo
!write(*,'(a,2i3, 4f10.6, i6)') ' iscat',i, 1, disc_list(i)%four_form_tsc(6,6,6,1), disc_list(i)%four_form_tsc(num(1)-5,num(2)-5,num(3)-7,1), m
!write(*,'(a,2i3, 8f10.6)') ' iscat',i, j, disc_list(i)%four_form_tsc(6,6,6,j), &
!                                          disc_list(i)%four_form_tsc(8,6,6,j), &
!                                          disc_list(i)%four_form_tsc(6,8,6,j), &
!                                          disc_list(i)%four_form_tsc(8,8,6,j)
   enddo loop_sym
!disc_list(i)%four_form_tsc(0,0,0,:) = cmplx(1.,1.)
!write(*,*)
!write(*,*) ' ZEROS ', i, any(abs(real(disc_list(i)%four_form_tsc(:,:,:,1  )))<0.000001),  &
!                         any(abs(imag(disc_list(i)%four_form_tsc(:,:,:,1  )))<0.000001),  &
!' ',&
!                         any(abs(real(disc_list(i)%four_form_tsc(:,:,:,2  )))<0.000001),  &
!                         any(abs(imag(disc_list(i)%four_form_tsc(:,:,:,2  )))<0.000001),  &
!' ',&
!                         any(abs(real(disc_list(i)%four_form_tsc(:,:,:,3  )))<0.000001),  &
!                         any(abs(imag(disc_list(i)%four_form_tsc(:,:,:,3  )))<0.000001),  &
!' ',&
!                         any(abs(real(disc_list(i)%four_form_tsc(:,:,:,4  )))<0.000001),  &
!                         any(abs(imag(disc_list(i)%four_form_tsc(:,:,:,4  )))<0.000001)
enddo loop_scat
!
end subroutine discamb_set_form_tsc
!
!*******************************************************************************
!
!
subroutine old_discamb_set_form_tsc(nscat)
!-
!  Copy the initial atomic form factors into the full table, apply symmetry
!+
use crystal_mod
use diffuse_mod
use wyckoff_mod
!
use matrix_mod
!
integer, intent(in) :: nscat    ! Number of scattering types
!
real(kind=PREC_DP), parameter :: TOL = 0.01_PREC_DP!
!
integer :: i,j   ! Dummy loop index
integer :: h,k,l,m ! Dummy loop index
integer :: il,jl,kl ! Dummy indices
integer :: ii,jj,kk ! Dummy indices
real(kind=PREC_DP), dimension(3)   :: veca, vecb, uvw    ! vectors
real(kind=PREC_DP), dimension(3,3) :: sym   ! Symmetry matrix
real(kind=PREC_DP), dimension(3,3) :: rsym  ! Symmetry matrix reciprocal space
real(kind=PREC_DP), dimension(3,3) :: viinv ! Inverse matrix to vsteps reciprocal space
!
call matinv(vi, viinv)
!
il = lbound(four_form_tsc,1)
jl = lbound(four_form_tsc,2)
kl = lbound(four_form_tsc,3)
ii = ubound(four_form_tsc,1)
jj = ubound(four_form_tsc,2)
kk = ubound(four_form_tsc,3)
loop_scat: do i=1, nscat
!   allocate(disc_list(i)%four_form_tsc(num(1), num(2), num(3), disc_list(i)%isymm(0)))
   allocate(disc_list(i)%four_form_tsc(il:ii, jl:jj, kl:kk, disc_list(i)%isymm(0)))
   disc_list(i)%four_form_tsc(:,:,:,1) = four_form_tsc(:,:,:,i)
!write(*,*)
!write(*,'(a,2i3, 8f10.6)') ' iscat',i, 1, disc_list(i)%four_form_tsc(6,6,6,1), &
!                                          disc_list(i)%four_form_tsc(8,6,6,1), &
!                                          disc_list(i)%four_form_tsc(6,8,6,1), &
!                                          disc_list(i)%four_form_tsc(8,8,6,1)
!write(*,'(a,2i3, 8f10.6)') ' iscat',i, 1, disc_list(i)%four_form_tsc(6,6,8,1), &
!                                          disc_list(i)%four_form_tsc(8,6,8,1), &
!                                          disc_list(i)%four_form_tsc(6,8,8,1), &
!                                          disc_list(i)%four_form_tsc(8,8,8,1)
!
   loop_sym: do j=2, disc_list(i)%isymm(0)     ! Loop over all symmetry operations 
      sym = spc_mat(1:3,1:3, disc_list(i)%isymm(j))
      rsym = matmul(cr_gten, matmul(sym, cr_rten))
 m = 0
      do l=1, num(1)                 ! Loop over all points in reciprocal space
         do k=1, num(2)
            do h=1, num(3)
               veca(1) = eck(1,1) + (h-1)*vi(1,1) + (k-1)*vi(1,2) + (l-1)*vi(1,3)
               veca(2) = eck(2,1) + (h-1)*vi(2,1) + (k-1)*vi(2,2) + (l-1)*vi(2,3)
               veca(3) = eck(3,1) + (h-1)*vi(3,1) + (k-1)*vi(3,2) + (l-1)*vi(3,3)
               vecb = matmul(rsym, veca)                  ! Apply reciprocal space symmetry operation
               uvw = matmul(viinv, vecb  - eck(:,1) )! Determine indices in form factor array
               if(all(uvw-nint(uvw)<TOL)) then   ! All indices are integer
                  ii = nint(uvw(1)) + 1
                  jj = nint(uvw(2)) + 1
                  kk = nint(uvw(3)) + 1
                  if(ii>0 .and. ii<=num(1) .and. jj>0 .and. jj<=num(2) .and. kk>0 .and. kk<=num(3)) then ! Within range
 m = m+1
!if(h==6 .and. k==6 .and. l==6) then
!write(*,'(a,3i3,3i3,3(2x,3f6.2))') 'VALID ', h,k,l, ii,jj,kk, veca, vecb, uvw
!endif

                        disc_list(i)%four_form_tsc(ii,jj,kk, j) = disc_list(i)%four_form_tsc(h,k,l,1)
!else
!write(*,*) 'RANGE ', h,k,l, ii,jj,kk, uvw
                  endif
!else
!write(*,*) 'INTEG ', h,k,l, uvw
               endif
            enddo
         enddo
      enddo
!write(*,'(a,2i3, 4f10.6, i6)') ' iscat',i, 1, disc_list(i)%four_form_tsc(6,6,6,1), disc_list(i)%four_form_tsc(num(1)-5,num(2)-5,num(3)-7,1), m
!write(*,'(a,2i3, 8f10.6)') ' iscat',i, j, disc_list(i)%four_form_tsc(6,6,6,j), &
!                                          disc_list(i)%four_form_tsc(8,6,6,j), &
!                                          disc_list(i)%four_form_tsc(6,8,6,j), &
!                                          disc_list(i)%four_form_tsc(8,8,6,j)
   enddo loop_sym
!write(*,*)
enddo loop_scat
!
end subroutine old_discamb_set_form_tsc
!
!*******************************************************************************
!
!subroutine discamb_check_hkl( hkl, ier_num )
!!-
!! Check if the user intended hkl max values are covered by discamb file
!!+
!!
!use precision_mod
!!
!implicit none
!!
!real(kind=PREC_DP), dimension(3), intent(in) :: hkl  ! User maximum hkl
!integer , intent(out) :: ier_num
!!
!ier_num = 0
!if(any(maxhkl< int(hkl))) then
!   ier_num = -1
!endif
!!
!end subroutine discamb_check_hkl
!
!*******************************************************************************
!
end module discamb_mod
