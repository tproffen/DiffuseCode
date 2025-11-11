module lib_h5fortran_mod
!
! functions to read/write arrays into from an HDF5 file via h5fortran library
! Since the hflibrary stops upon errors, these intermediate functions use the
! checks available in the library to gracefully return to the calling program.
!
use precision_mod
use hdf5, only : HSIZE_T
!
private
public h5f_open
public h5f_read
public h5f_write
!
interface h5f_read
   module procedure read_char,         & !  Read a character string
                    read_int,          & !  Read a single integer number 
                    read_int_1D,       & !  Read directly into 1D array
                    read_int_2D,       & !  Read directly into 2D array
                    read_real_1D,      & !  Real valued 1D array
                    read_real_2D,      & !  Real valued 2D array
                    read_real_3D         !  Real valued 3D array
end interface h5f_read
!
interface h5f_write
   module procedure write_char,         & !  Write a character string
                    write_int,          & !  Write a single integer number 
                    write_int_1D,       & !  Write directly into 1D array
                    write_int_2D,       & !  Write directly into 2D array
                    write_real_1D,      & !  Write valued 1D array
                    write_real_2D,      & !  Write valued 2D array
                    write_real_3D         !  Write valued 3D array
end interface h5f_write
!
! error numbers:
! -1  File does not exist
! -2  Data set does not exist
! -3  Wrong rank
! -4  Wrong dimensions
! -5  Not an HDF5 file
!
!*******************************************************************************
!
contains
!
!*******************************************************************************
!
subroutine h5f_open(h5f, infile, read_write, NMSG, ier_num, ier_msg)
!-
!  Check for existence and open infile
!  read_write flag is either 'r' or 'w'
!+
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)                      , intent(inout)  :: h5f
character(len=*)                     , intent(in)    :: infile
character(len=*)                     , intent(in)    :: read_write
integer                              , intent(in)    :: NMSG
integer                              , intent(out)   :: ier_num
character(len=*)  , dimension(NMSG)  , intent(out)   :: ier_msg
!
integer, parameter :: IRD =65
character(len=16) :: line
logical :: lexist
!
if(read_write=='r') then
   inquire(file=infile, exist=lexist)
   if(lexist) then
      open(unit=IRD, file=infile, status='old')
      read(IRD,'(a)') line
      close(IRD)
      if(line(2:5)/='HDF') then
         ier_num = -5
         ier_msg(1) = 'File is not an HDF5 file'
         ier_msg(2) = infile
         return
      endif
   endif
endif
if(lexist.or. read_write=='w') then         ! Input file exists or write
   call h5f%open(infile, action=read_write)
else
   ier_num = -1
   ier_msg(1) = infile
endif
! 
end subroutine h5f_open
!
!*******************************************************************************
!
subroutine read_char(h5f, data_set, res_string, NMSG, ier_num, ier_msg)
!-
!   Read a character string
!+
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)           , intent(in)  :: h5f
character(len=*)          , intent(in)  :: data_set
character(len=*)          , intent(out) :: res_string
integer                              , intent(in)  :: NMSG
integer                              , intent(out) :: ier_num
character(len=*)  , dimension(NMSG)  , intent(out) :: ier_msg
!
character(len=PREC_STRING) :: string
!
ier_num = 0
ier_msg = ' '
!
string = '/entry/data/' // data_set(1:len_trim(data_set))
if(h5f%exist(string(1:len_trim(string)))) then
   call h5f%read(string(1:len_trim(string)), res_string)
else
   ier_num    = -2             ! Does not exist
   ier_msg(1) = data_set
endif
!
end subroutine read_char
!
!*******************************************************************************
!
subroutine read_int(h5f, data_set, iresult, NMSG, ier_num, ier_msg)
!-
!  Read a single integer number
!+
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)           , intent(in)  :: h5f
character(len=*)          , intent(in)  :: data_set
integer                   , intent(out) :: iresult
integer                              , intent(in)  :: NMSG
integer                              , intent(out) :: ier_num
character(len=*)  , dimension(NMSG)  , intent(out) :: ier_msg
!
character(len=PREC_STRING) :: string
!
string = '/entry/data/' // data_set(1:len_trim(data_set))
if(h5f%exist(string(1:len_trim(string)))) then
  call h5f%read(string(1:len_trim(string)), iresult)
else
   ier_num    = -2             ! Does not exist
   ier_msg(1) = data_set
endif
!
end subroutine read_int
!
!*******************************************************************************
!
subroutine read_int_1D(h5f, data_set, MAXDIM, iarray, NMSG, ier_num, ier_msg)
!-
!  Read an integer 1D array
!+
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)           , intent(in)  :: h5f
character(len=*)          , intent(in)  :: data_set
integer                   , intent(in)  :: MAXDIM
integer, dimension(MAXDIM), intent(out) :: iarray
integer                              , intent(in)  :: NMSG
integer                              , intent(out) :: ier_num
character(len=*)  , dimension(NMSG)  , intent(out) :: ier_msg
!
character(len=PREC_STRING) :: string
!
string = '/entry/data/' // data_set(1:len_trim(data_set))
if(h5f%exist(string(1:len_trim(string)))) then
  call h5f%read(string(1:len_trim(string)), iarray)
else
   ier_num    = -2             ! Does not exist
   ier_msg(1) = data_set
endif
!
end subroutine read_int_1D
!
!*******************************************************************************
!
subroutine read_int_2D(h5f, data_set, MAXDIM, iarray, NMSG, ier_num, ier_msg)
!-
!   Read an integer valued 2D array, HDF5 array needs to be transposed
!+
use precision_mod
!
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)                      , intent(in)  :: h5f
character(len=*)                     , intent(in)  :: data_set
integer           , dimension(2)     , intent(in)  :: MAXDIM
integer           , dimension(MAXDIM(1), MAXDIM(2)), intent(out) :: iarray
integer                              , intent(in)  :: NMSG
integer                              , intent(out) :: ier_num
character(len=*)  , dimension(NMSG)  , intent(out) :: ier_msg
!
character(len=PREC_STRING) :: string
integer                    :: ndim     ! Rank of actual array
integer(kind=HSIZE_T), dimension(:)  , allocatable :: dims
integer              , dimension(:,:), allocatable ::  int_array_2D
!
ier_num = 0
ier_msg = ' '
!
string = '/entry/data/' // data_set(1:len_trim(data_set))
if(h5f%exist(string(1:len_trim(string)))) then
   ndim = h5f%ndim(string(1:len_trim(string)))
   if(ndim==2) then                  ! Correct rank
      call h5f%shape(string(1:len_trim(string)),dims)
      if(dims(1) == MAXDIM(2) .and. dims(2) == MAXDIM(1)) then     ! Correct dimensions
         allocate( int_array_2D(dims(1), dims(2)))
         call h5f%read(string(1:len_trim(string)),  int_array_2D)
         iarray = transpose( int_array_2D)
         deallocate( int_array_2D)
      else                           ! Wrong dimensions
         ier_num = -4
         ier_msg(1) = data_set
         write(ier_msg(2),'(a,i3,a,i3)') 'Found ', dims(1), ' Expected ', MAXDIM
      endif
   else
      ier_num = -3                   ! Wrong rank
      ier_msg(1) = data_set
      write(ier_msg(2),'(a,2i3)') 'Found    ', dims
      write(ier_msg(3),'(a,2i3)') 'Expected ', MAXDIM
   endif
else
   ier_num    = -2                  ! Does not exist
   ier_msg(1) = data_set
endif
!
end subroutine read_int_2D
!
!*******************************************************************************
!
subroutine read_real_1D(h5f, data_set, MAXDIM, rarray, NMSG, ier_num, ier_msg)
!-
!   Read a real valued 1D array
!+
use precision_mod
!
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)                      , intent(in)  :: h5f
character(len=*)                     , intent(in)  :: data_set
integer                              , intent(in)  :: MAXDIM
real(KIND=PREC_DP), dimension(MAXDIM), intent(out) :: rarray
integer                              , intent(in)  :: NMSG
integer                              , intent(out) :: ier_num
character(len=*)  , dimension(NMSG)  , intent(out) :: ier_msg
!
character(len=PREC_STRING) :: string
integer                    :: ndim     ! Rank of actual array
integer(kind=HSIZE_T), dimension(:), allocatable :: dims
!
ier_num = 0
ier_msg = ' '
!
string = '/entry/data/' // data_set(1:len_trim(data_set))
if(h5f%exist(string(1:len_trim(string)))) then
   ndim = h5f%ndim(string(1:len_trim(string)))
   if(ndim==1) then                  ! Correct rank
      call h5f%shape(string(1:len_trim(string)),dims)
      if(dims(1) == MAXDIM) then     ! Correct dimensions
         call h5f%read(string(1:len_trim(string)), rarray)
      else                           ! Wrong dimensions
         ier_num = -4
         ier_msg(1) = data_set
         write(ier_msg(2),'(a,i3,a,i3)') 'Found ', dims(1), ' Expected ', MAXDIM
      endif
   elseif(ndim/=1) then  
      ier_num = -3                   ! Wrong rank
      ier_msg(1) = data_set
      write(ier_msg(2),'(a,i3,a,i3)') 'Found ', ndim, ' Expected ', 1
   endif
else
   ier_num    = -2                  ! Does not exist
   ier_msg(1) = data_set
endif
!
end subroutine read_real_1D
!
!*******************************************************************************
!
subroutine read_real_2D(h5f, data_set, MAXDIM, rarray, NMSG, ier_num, ier_msg)
!-
!   Read a real valued 2D array, HDF5 array needs to be transposed
!+
use precision_mod
!
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)                      , intent(in)  :: h5f
character(len=*)                     , intent(in)  :: data_set
integer           , dimension(2)     , intent(in)  :: MAXDIM
real(KIND=PREC_DP), dimension(MAXDIM(1), MAXDIM(2)), intent(out) :: rarray
integer                              , intent(in)  :: NMSG
integer                              , intent(out) :: ier_num
character(len=*)  , dimension(NMSG)  , intent(out) :: ier_msg
!
character(len=PREC_STRING) :: string
integer                    :: ndim     ! Rank of actual array
integer(kind=HSIZE_T), dimension(:)  , allocatable :: dims
real(kind=PREC_DP)   , dimension(:,:), allocatable :: real_array_2D
!
ier_num = 0
ier_msg = ' '
!
string = '/entry/data/' // data_set(1:len_trim(data_set))
if(h5f%exist(string(1:len_trim(string)))) then
   ndim = h5f%ndim(string(1:len_trim(string)))
   if(ndim==2) then                  ! Correct rank
      call h5f%shape(string(1:len_trim(string)),dims)
      if(dims(1) == MAXDIM(2) .and. dims(2) == MAXDIM(1)) then     ! Correct dimensions
         allocate(real_array_2D(dims(1), dims(2)))
         call h5f%read(string(1:len_trim(string)), real_array_2D)
         rarray = transpose(real_array_2D)
         deallocate(real_array_2D)
      else                           ! Wrong dimensions
         ier_num = -4
         ier_msg(1) = data_set
         write(ier_msg(2),'(a,i3,a,i3)') 'Found ', dims(1), ' Expected ', MAXDIM
      endif
   else
      ier_num = -3                   ! Wrong rank
      ier_msg(1) = data_set
      write(ier_msg(2),'(a,2i3)') 'Found    ', dims
      write(ier_msg(3),'(a,2i3)') 'Expected ', MAXDIM
   endif
else
   ier_num    = -2                  ! Does not exist
   ier_msg(1) = data_set
endif
!
end subroutine read_real_2D
!
!*******************************************************************************
!
subroutine read_real_3D(h5f, data_set, MAXDIM, rarray, NMSG, ier_num, ier_msg)
!-
!   Read a real valued 3D array, HDF5 array needs to be transposed
!+
use precision_mod
!
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)                      , intent(in)  :: h5f
character(len=*)                     , intent(in)  :: data_set
integer           , dimension(3)     , intent(in)  :: MAXDIM
real(KIND=PREC_DP), dimension(MAXDIM(1), MAXDIM(2), MAXDIM(3)), intent(out) :: rarray
integer                              , intent(in)  :: NMSG
integer                              , intent(out) :: ier_num
character(len=*)  , dimension(NMSG)  , intent(out) :: ier_msg
!
character(len=PREC_STRING) :: string
integer                    :: ndim     ! Rank of actual array
integer                    :: i,j,k
integer(kind=HSIZE_T), dimension(:)  , allocatable :: dims
real(kind=PREC_DP)   , dimension(:,:,:), allocatable :: real_array_3D
!
ier_num = 0
ier_msg = ' '
!
string = '/entry/data/' // data_set(1:len_trim(data_set))
if(h5f%exist(string(1:len_trim(string)))) then
   ndim = h5f%ndim(string(1:len_trim(string)))
   if(ndim==3) then                  ! Correct rank
      call h5f%shape(string(1:len_trim(string)),dims)
      if(dims(1) == MAXDIM(3) .and. dims(2) == MAXDIM(2) .and. dims(3) == MAXDIM(1)) then     ! Correct dimensions
         allocate(real_array_3D(dims(1), dims(2), dims(3)))
         call h5f%read(string(1:len_trim(string)), real_array_3D)
         do i=1,maxdim(1)
            do j=1,maxdim(2)
               do k=1,maxdim(3)
                  rarray(i,j,k) = real_array_3D(k,j,i)
               enddo
            enddo
         enddo
         deallocate(real_array_3D)
      else                           ! Wrong dimensions
         ier_num = -4
         ier_msg(1) = data_set
         write(ier_msg(2),'(a,3i3)') 'Found    ', dims
         write(ier_msg(3),'(a,3i3)') 'Expected ', MAXDIM
      endif
   else
      ier_num = -3                   ! Wrong rank
      ier_msg(1) = data_set
      write(ier_msg(2),'(a,i3,a,i3)') 'Found ', ndim, ' Expected ', 3
   endif
else
   ier_num    = -2                  ! Does not exist
   ier_msg(1) = data_set
endif
!
end subroutine read_real_3D
!
!*******************************************************************************
!*******************************************************************************
!
!*******************************************************************************
!
subroutine write_char(h5f, data_set,  in_string, NMSG, ier_num, ier_msg)
!-
!   Write a character string
!+
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)           , intent(in)  :: h5f
character(len=*)          , intent(in)  :: data_set
character(len=*)          , intent(in ) ::  in_string
integer                              , intent(in)  :: NMSG
integer                              , intent(out) :: ier_num
character(len=*)  , dimension(NMSG)  , intent(out) :: ier_msg
!
character(len=PREC_STRING) :: string
!
ier_num = 0
ier_msg = ' '
!
string = '/entry/data/' // data_set(1:len_trim(data_set))
!if(h5f%exist(string(1:len_trim(string)))) then
   call h5f%write(string(1:len_trim(string)),  in_string(1:len_trim( in_string)))
!else
!   ier_num    = -2             ! Does not exist
!   ier_msg(1) = data_set
!endif
!
end subroutine write_char
!
!*******************************************************************************
!
subroutine write_int(h5f, data_set, iresult, NMSG, ier_num, ier_msg)
!-
!  Write a single integer number as array of length 1
!+
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)           , intent(in)  :: h5f
character(len=*)          , intent(in)  :: data_set
integer                   , intent(in)  :: iresult
integer                              , intent(in)  :: NMSG
integer                              , intent(out) :: ier_num
character(len=*)  , dimension(NMSG)  , intent(out) :: ier_msg
!
character(len=PREC_STRING) :: string
integer, dimension(1) :: iarray
!
string = '/entry/data/' // data_set(1:len_trim(data_set))
!if(h5f%exist(string(1:len_trim(string)))) then
  iarray(1) = iresult 
  call h5f%write(string(1:len_trim(string)), iarray)
!else
!   ier_num    = -2             ! Does not exist
!   ier_msg(1) = data_set
!endif
!
end subroutine write_int
!
!*******************************************************************************
!
subroutine write_int_1D(h5f, data_set, MAXDIM, iarray, NMSG, ier_num, ier_msg)
!-
!  Write an integer 1D array
!+
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)           , intent(in)  :: h5f
character(len=*)          , intent(in)  :: data_set
integer                   , intent(in)  :: MAXDIM
integer, dimension(MAXDIM), intent(in)  :: iarray
integer                              , intent(in)  :: NMSG
integer                              , intent(out) :: ier_num
character(len=*)  , dimension(NMSG)  , intent(out) :: ier_msg
!
character(len=PREC_STRING) :: string
!
string = '/entry/data/' // data_set(1:len_trim(data_set))
!if(h5f%exist(string)) then
  call h5f%write(string(1:len_trim(string)), iarray)
!else
!   ier_num    = -2             ! Does not exist
!   ier_msg(1) = data_set
!endif
!
end subroutine write_int_1D
!
!*******************************************************************************
!
subroutine write_int_2D(h5f, data_set, MAXDIM, iarray, NMSG, ier_num, ier_msg)
!-
!   Write a integer valued 2D array, HDF5 array needs to be transposed
!+
use precision_mod
!
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)                      , intent(in)  :: h5f
character(len=*)                     , intent(in)  :: data_set
integer           , dimension(2)     , intent(in)  :: MAXDIM
integer           , dimension(MAXDIM(1), MAXDIM(2)), intent(in) :: iarray
integer                              , intent(in)  :: NMSG
integer                              , intent(out) :: ier_num
character(len=*)  , dimension(NMSG)  , intent(out) :: ier_msg
!
character(len=PREC_STRING) :: string
!integer                    :: ndim     ! Rank of actual array
!integer                    :: i,j,k
!integer(kind=HSIZE_T), dimension(:)  , allocatable :: dims
integer              , dimension(:,:), allocatable ::  int_array_2D
!
ier_num = 0
ier_msg = ' '
!
string = '/entry/data/' // data_set(1:len_trim(data_set))
!if(h5f%exist(string(1:len_trim(string)))) then
!   ndim = h5f%ndim(string(1:len_trim(string)))
!   if(ndim==3) then                  ! Correct rank
!      call h5f%shape(string(1:len_ltrim(string)),dims)
!      if(dims(1) == MAXDIM(3) .and. dims(2) == MAXDIM(2) .and. dims(3) == MAXDIM(1)) then     ! Correct dimensions
         allocate( int_array_2D(MAXDIM(2), MAXDIM(1)))
          int_array_2D = transpose(iarray)
         call h5f%write(string(1:len_trim(string)),  int_array_2D)
         deallocate( int_array_2D)
!      else                           ! Wrong dimensions
!         ier_num = -4
!         ier_msg(1) = data_set
!         write(ier_msg(2),'(a,3i3)') 'Found    ', dims
!         write(ier_msg(3),'(a,3i3)') 'Expected ', MAXDIM
!      endif
!   else
!      ier_num = -3                   ! Wrong rank
!      ier_msg(1) = data_set
!      write(ier_msg(2),'(a,i3,a,i3)') 'Found ', ndim, ' Expected ', 3
!   endif
!else
!   ier_num    = -2                  ! Does not exist
!   ier_msg(1) = data_set
!endif
!
end subroutine write_int_2D
!
!*******************************************************************************
!
subroutine write_real_1D(h5f, data_set, MAXDIM, rarray, NMSG, ier_num, ier_msg)
!-
!  Write an real valued 1D array
!+
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)           , intent(in)  :: h5f
character(len=*)          , intent(in)  :: data_set
integer                   , intent(in)  :: MAXDIM
real(kind=PREC_DP), dimension(MAXDIM), intent(in)  :: rarray
integer                              , intent(in)  :: NMSG
integer                              , intent(out) :: ier_num
character(len=*)  , dimension(NMSG)  , intent(out) :: ier_msg
!
character(len=PREC_STRING) :: string
!
!if(h5f%exist(string(1:len_trim(string)))) then
   string = '/entry/data/' // data_set(1:len_trim(data_set))
   call h5f%write(string(1:len_trim(string)), rarray)
!else
!   ier_num    = -5             ! Does not exist
!   ier_msg(1) = '/entry/data'
!endif
!
end subroutine write_real_1D
!
!*******************************************************************************
!
subroutine write_real_2D(h5f, data_set, MAXDIM, rarray, NMSG, ier_num, ier_msg)
!-
!   Write a real valued 2D array, HDF5 array needs to be transposed
!+
use precision_mod
!
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)                      , intent(in)  :: h5f
character(len=*)                     , intent(in)  :: data_set
integer           , dimension(2)     , intent(in)  :: MAXDIM
real(KIND=PREC_DP), dimension(MAXDIM(1), MAXDIM(2)), intent(in) :: rarray
integer                              , intent(in)  :: NMSG
integer                              , intent(out) :: ier_num
character(len=*)  , dimension(NMSG)  , intent(out) :: ier_msg
!
character(len=PREC_STRING) :: string
!integer                    :: ndim     ! Rank of actual array
!integer                    :: i,j,k
!integer(kind=HSIZE_T), dimension(:)  , allocatable :: dims
real(kind=PREC_DP)   , dimension(:,:), allocatable :: real_array_2D
!
ier_num = 0
ier_msg = ' '
!
string = '/entry/data/' // data_set(1:len_trim(data_set))
!if(h5f%exist(string(1:len_trim(string)))) then
!   ndim = h5f%ndim(string(1:len_trim(string)))
!   if(ndim==3) then                  ! Correct rank
!      call h5f%shape(string(1:len_trim(string)),dims)
!      if(dims(1) == MAXDIM(3) .and. dims(2) == MAXDIM(2) .and. dims(3) == MAXDIM(1)) then     ! Correct dimensions
         allocate(real_array_2D(MAXDIM(2), MAXDIM(1)))
         real_array_2D = transpose(rarray)
         call h5f%write(string(1:len_trim(string)), real_array_2D)
         deallocate(real_array_2D)
!      else                           ! Wrong dimensions
!         ier_num = -4
!         ier_msg(1) = data_set
!         write(ier_msg(2),'(a,3i3)') 'Found    ', dims
!         write(ier_msg(3),'(a,3i3)') 'Expected ', MAXDIM
!      endif
!   else
!      ier_num = -3                   ! Wrong rank
!      ier_msg(1) = data_set
!      write(ier_msg(2),'(a,i3,a,i3)') 'Found ', ndim, ' Expected ', 3
!   endif
!else
!   ier_num    = -2                  ! Does not exist
!   ier_msg(1) = data_set
!endif
!
end subroutine write_real_2D
!
!*******************************************************************************
!
subroutine write_real_3D(h5f, data_set, MAXDIM, rarray, NMSG, ier_num, ier_msg)
!-
!   Write a real valued 3D array, HDF5 array needs to be transposed
!+
use precision_mod
!
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)                      , intent(in)  :: h5f
character(len=*)                     , intent(in)  :: data_set
integer           , dimension(3)     , intent(in)  :: MAXDIM
real(KIND=PREC_DP), dimension(MAXDIM(1), MAXDIM(2), MAXDIM(3)), intent(in) :: rarray
integer                              , intent(in)  :: NMSG
integer                              , intent(out) :: ier_num
character(len=*)  , dimension(NMSG)  , intent(out) :: ier_msg
!
character(len=PREC_STRING) :: string
!integer                    :: ndim     ! Rank of actual array
integer                    :: i,j,k
!integer(kind=HSIZE_T), dimension(:)  , allocatable :: dims
real(kind=PREC_DP)   , dimension(:,:,:), allocatable :: real_array_3D
!
ier_num = 0
ier_msg = ' '
!
string = '/entry/data/' // data_set(1:len_trim(data_set))
!if(h5f%exist(string(1:len_trim(string)))) then
!   ndim = h5f%ndim(string(1:len_trim(string)))
!   if(ndim==3) then                  ! Correct rank
!      call h5f%shape(string(1:len_trim(string)),dims)
!      if(dims(1) == MAXDIM(3) .and. dims(2) == MAXDIM(2) .and. dims(3) == MAXDIM(1)) then     ! Correct dimensions
         allocate(real_array_3D(MAXDIM(3), MAXDIM(2), MAXDIM(1)))
         do i=1,maxdim(1)
            do j=1,maxdim(2)
               do k=1,maxdim(3)
                  real_array_3D(k,j,i) = rarray(i,j,k)
               enddo
            enddo
         enddo
!do i=1,MAXDIM(3)
!write(*,'(a,i8    )') ' Symmetry matrix No. ', i
!write(*,'(a,8f12.6)') ' Symmetry matrix     ', rarray(1,:,i), real_array_3D(i,:,1)
!write(*,'(a,8f12.6)') ' Symmetry matrix     ', rarray(2,:,i), real_array_3D(i,:,2)
!write(*,'(a,8f12.6)') ' Symmetry matrix     ', rarray(3,:,i), real_array_3D(i,:,3)
!enddo
!         call h5f%write(string(1:len_trim(string)), rarray)
         call h5f%write(string(1:len_trim(string)), real_array_3D)
         deallocate(real_array_3D)
!      else                           ! Wrong dimensions
!         ier_num = -4
!         ier_msg(1) = data_set
!         write(ier_msg(2),'(a,3i3)') 'Found    ', dims
!         write(ier_msg(3),'(a,3i3)') 'Expected ', MAXDIM
!      endif
!   else
!      ier_num = -3                   ! Wrong rank
!      ier_msg(1) = data_set
!      write(ier_msg(2),'(a,i3,a,i3)') 'Found ', ndim, ' Expected ', 3
!   endif
!else
!   ier_num    = -2                  ! Does not exist
!   ier_msg(1) = data_set
!endif
!
end subroutine write_real_3D
!
!*******************************************************************************
!
end module lib_h5fortran_mod
