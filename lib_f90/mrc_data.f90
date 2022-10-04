module lib_mrc_mod
!-
!  Routines to write and read MRC file format
!+
!
use errlist_mod
!
implicit none
!
private
public mrc_write
public mrc_read
!
!
contains
!
!*******************************************************************************
!
subroutine mrc_write(outfile,                                                   &
                     dims, mode, nxyzstart, cr_a, cr_win, extr_abs, extr_ord,   &
                     extr_top, odata                                            &
                    )
!
use envir_mod
use errlist_mod
use lib_errlist_func
use lib_length
use precision_mod
!
implicit none
!
character(len=*)                , intent(in) :: outfile
integer           , dimension(3), intent(in) :: dims
integer                         , intent(in) :: mode
integer           , dimension(3), intent(in) :: nxyzstart
real(kind=PREC_DP), dimension(3), intent(in) :: cr_a
real(kind=PREC_DP), dimension(3), intent(in) :: cr_win
integer                         , intent(in) :: extr_abs
integer                         , intent(in) :: extr_ord
integer                         , intent(in) :: extr_top
real(kind=PREC_DP), dimension(dims(1), dims(2), dims(3)), intent(in) :: odata
!
integer, parameter :: IMRC = 97
!
!integer, parameter :: mode = 2
!
!integer, parameter :: nxstart = 0
!integer, parameter :: nystart = 0
!integer, parameter :: nzstart = 0
!
integer, parameter :: record   =    4   ! Record is four Byte long
integer, parameter :: base_hdr =  256   ! Standard header size in Words=1024 Byte
integer, parameter :: extend_b = 3072   ! = (1024 - 256)*4 extended in Bytes
!                                          ! extend_b must be multiple of 4
integer, parameter :: extend_w = extend_b/record ! =(1024 - 256)   extended in Words
!
character(LEN=PREC_STRING) :: filename
character(LEN=PREC_STRING) :: line
character(LEN=PREC_STRING) :: message
integer            :: i, j,l, irec, ios
integer            :: l_datei
!real(kind=PREC_DP)               :: sqq
!
!
call no_error
!
l_datei = len_str(outfile)
if(outfile(1:1) .eq.'~') then 
   line = ' '
   line = home_dir(1:home_dir_l) // outfile(2:l_datei)
   filename = line(1:200)
else
   filename = outfile
endif
 
open(unit=IMRC,file=filename,status='unknown',        &
     form='unformatted', access='direct',recl=record, &
     iostat=ios,iomsg=message)
if(ios/=0) then
   ier_num = -2
   ier_typ = ER_IO
   ier_msg(3) = message(1:80)
   return
endif
!
!  Write header is identical to all 
!
write(IMRC,rec= 1) dims(1)              ! Dimension column
write(IMRC,rec= 2) dims(2)              ! Dimension row
write(IMRC,rec= 3) dims(3)              ! Dimension slices
write(IMRC,rec= 4) 2                    ! Data type = 2 32-bit real
write(IMRC,rec= 5) nxyzstart(1)         ! Number of first column in map = 0
write(IMRC,rec= 6) nxyzstart(2)         ! Number of first row    in map = 0
write(IMRC,rec= 7) nxyzstart(3)         ! Number of first slice  in map = 0
write(IMRC,rec= 8) dims(1)-1            ! Number of intervals along x
write(IMRC,rec= 9) dims(2)-1            ! Number of intervals along y
write(IMRC,rec=10) dims(3)-1            ! Number of intervals along z
write(IMRC,rec=11) real(cr_a(1))        ! Reciprocal lattice parameter a
write(IMRC,rec=12) real(cr_a(2))        ! Reciprocal lattice parameter b
write(IMRC,rec=13) real(cr_a(3))        ! Reciprocal lattice parameter c
write(IMRC,rec=14) real(cr_win(1))      ! Reciprocal lattice parameter alpha
write(IMRC,rec=15) real(cr_win(2))      ! Reciprocal lattice parameter beta
write(IMRC,rec=16) real(cr_win(3))      ! Reciprocal lattice parameter gamma
write(IMRC,rec=17) real(extr_abs)       ! Reciprocal axis along columns 
write(IMRC,rec=18) real(extr_ord)       ! Reciprocal axis along rows 
write(IMRC,rec=19) real(extr_top)       ! Reciprocal axis along sections 
write(IMRC,rec=20) real(minval(odata))  ! Minimum density
write(IMRC,rec=21) real(maxval(odata))  ! Maximum density
write(IMRC,rec=22) real(sum(odata)/(dims(1)*dims(2)*dims(3)))        ! Average density
write(IMRC,rec=23) 0                    ! 0 = image
write(IMRC,rec=24) extend_b             ! Extended header length in Bytes
do I = 25,49
   write(IMRC,rec= i) 0                 ! 0
enddo
write(IMRC,rec=50) 0.0                  ! Origin
write(IMRC,rec=51) 0.0                  ! Origin
write(IMRC,rec=52) 0.0                  ! Origin
write(IMRC,rec=53) 'MAP'                ! Identify file type as MAP
write(IMRC,rec=54) 4+4*16+68*256               ! Origin, 
do I = 55,base_hdr                      ! No info in the rest of base header
   write(IMRC,rec= i) 0                 ! 0
enddo
do I = base_hdr+1,base_hdr+extend_w     ! No info in the extended header
   write(IMRC,rec= i) 0                 ! 0
enddo
!
irec = base_hdr + extend_w              ! Image is after base + extended header
do l = 1, dims(3)
   do j = 1, dims(2)
      do i = 1, dims(1)
         irec = irec + 1
         write(IMRC,rec=irec) real(odata(i,j,l))
      enddo
   enddo
enddo
!
close(IMRC)
!
end subroutine mrc_write
!
!*******************************************************************************
!
subroutine mrc_read(outfile, node_number)
!                                       &
! ) 
!                    dims, mode, nxyzstart, cr_a, cr_win                       &
! ) 
!                   )
!
use envir_mod
use errlist_mod
use lib_errlist_func
use lib_data_struc_h5
use lib_length
use param_mod
use precision_mod
!
implicit none
!
character(len=*)                , intent(in)  :: outfile
integer                         , intent(out) :: node_number      ! Node in global storage
!
real(kind=PREC_DP), dimension(:,:,:), allocatable              :: odata
!character(len=3 )                             :: cmap
integer                                       :: ispg
integer           , dimension(3)              :: dims_step
!integer :: machst
integer :: extend_b        ! extended header in Bytes
!integer :: extend_w        ! extended header in word
integer(kind=PREC_INT_BYTE) :: ibyte
!
integer, parameter :: IMRC = 97
!
!integer, parameter :: mode = 2
!
!integer, parameter :: nxstart = 0
!integer, parameter :: nystart = 0
!integer, parameter :: nzstart = 0
!
!integer, parameter :: record   =    4   ! Record is four Byte long
!integer, parameter :: base_hdr =  256   ! Standard header size in Words=1024 Byte
!integer, parameter :: extend_b = 3072   ! = (1024 - 256)*4 extended in Bytes
!                                          ! extend_b must be multiple of 4
!integer, parameter :: extend_w = extend_b/record ! =(1024 - 256)   extended in Words
!
character(LEN=PREC_STRING) :: filename
character(LEN=PREC_STRING) :: line
character(LEN=PREC_STRING) :: message
character(len=80)  :: label
integer            :: i, ios
integer            :: l_datei
integer            :: ilabel      ! Number of text labels in header
integer           , dimension(3) :: dims
integer                          :: mode
integer           , dimension(3) :: nxyzstart
real(kind=PREC_DP), dimension(3) :: cr_a
real(kind=PREC_DP), dimension(3) :: cr_win
integer           , dimension(3) :: extr_abs
real(kind=PREC_SP), dimension(3) :: m_cr_a
real(kind=PREC_SP), dimension(3) :: m_cr_win
integer(kind=PREC_INT_SHORT), dimension(2) :: short       ! If MODE==1  16 Bit signed integer
integer(kind=PREC_INT_SHORT)               :: numint      ! If MODE==1  16 Bit signed integer
integer(kind=PREC_INT_SHORT)               :: numfloat    ! If MODE==1  16 Bit signed integer
!real(kind=PREC_SP)               :: sqq         ! IF MODE==2  32Byte floating
real(kind=PREC_SP)               :: maxdata, mindata
real(kind=PREC_SP)               :: avdata
integer(kind=PREC_INT_SHORT), dimension(:,:,:), allocatable :: short_map
!
! Variable for global storage
!
integer   :: nndims          ! Number of dimensions 
integer                 :: mrc_layer
logical                 :: mrc_direct
logical                 :: mrc_is_grid
logical                 :: mrc_has_dxyz
logical                 :: mrc_has_dval
real(kind=PREC_DP), dimension(3,4) :: mrc_corners
real(kind=PREC_DP), dimension(3,4) :: mrc_vectors
real(kind=PREC_DP), dimension(6  ) :: mrc_unit
real(kind=PREC_DP), dimension(:  ), allocatable :: mrc_x
real(kind=PREC_DP), dimension(:  ), allocatable :: mrc_y
real(kind=PREC_DP), dimension(:  ), allocatable :: mrc_z
real(kind=PREC_DP), dimension(:  ), allocatable :: mrc_dx
real(kind=PREC_DP), dimension(:  ), allocatable :: mrc_dy
real(kind=PREC_DP), dimension(:  ), allocatable :: mrc_dz
real(kind=PREC_DP), dimension(:,:,:), allocatable :: osigma
real(kind=PREC_DP), dimension(3  ) :: mrc_llims
real(kind=PREC_DP), dimension(3  ) :: mrc_steps
real(kind=PREC_DP), dimension(3,3) :: mrc_steps_full
real(kind=PREC_DP)                 :: mrc_pixel_size
real(kind=PREC_SP) :: pxsize
!
!
call no_error
!
l_datei = len_str(outfile)
if(outfile(1:1) .eq.'~') then 
   line = ' '
   line = home_dir(1:home_dir_l) // outfile(2:l_datei)
   filename = line(1:200)
else
   filename = outfile
endif
 
open(unit=IMRC,file=filename,status='unknown',        &
     form='unformatted', access='stream',             &
     iostat=ios,iomsg=message)
if(ios/=0) then
   ier_num = -2
   ier_typ = ER_IO
   ier_msg(3) = message(1:80)
   return
endif
!
!  Read header; is identical to all 
!
read(IMRC) dims                 ! Dimension column, row, slices
read(IMRC) mode                 ! Data type = 2 32-bit real
read(IMRC) nxyzstart            ! Number of first column/row/slice in map = 0
read(IMRC) dims_step            ! Number of intervals along x
read(IMRC) m_cr_a               ! Lattice parameters or Overall pixel numbers 
read(IMRC) m_cr_win          ! Cell angles or image angles
read(IMRC) extr_abs        ! Reciprocal axis along columns 
read(IMRC) mindata              ! Minimum density
read(IMRC) maxdata              ! Maximum density
read(IMRC) avdata                                            ! Average density
read(IMRC) short(1)             ! 0 = image
read(IMRC) short(2)             ! number of bytes for storing symmetry information
ispg = short(1)
!write(*,*) ' MRC DIMS     ', dims
!write(*,*) ' MRC MODE     ', mode
!write(*,*) ' MRC START    ', nxyzstart
!write(*,*) ' MRC STEPS    ', dims_step
!write(*,*) ' MRC A0       ', m_cr_a
!write(*,*) ' MRC win      ', m_cr_win
!write(*,*) ' MRC EXTR     ', extr_abs
!write(*,*) ' MRC min      ', mindata
!write(*,*) ' MRC max      ', maxdata
!write(*,*) ' MRC ave      ',  avdata
!write(*,*) ' MRC ispg     ', 88, ispg
!write(*,*) ' MRC symm     ', 90, short(2)
read(IMRC) extend_b                ! Extended header length in Bytes
!write(*,*) ' MRC Extended ', 92, extend_b
read(IMRC) short(1)             ! 0 = image
!write(*,*) ' MRC Creator  ', 96, short(1)
do i = 98,127
   read(IMRC) ibyte             ! 0
!  write(*,*) ' MRC BYTE     ',  i, ibyte
enddo
read(IMRC) numint, numfloat    ! 0 = image
!write(*,*) ' MRC NI NF    ', 128, numint, numfloat
do i = 132, 219
   read(IMRC) ibyte             ! 0
!  write(*,*) ' MRC BYTE     ',  i, ibyte
enddo
read(IMRC) ilabel
!write(*,*) ' MRC LABELS   ', 220, ilabel
do i=1, ilabel
   read(IMRC) label
!   write(*,*) ' MRC LABELS   ', 224 +(i-1)*80, label
enddo
do i=1, 10-ilabel
   read(IMRC) label
enddo
!write(*,*) 'Starting to read extended header'
!  End of standard header at 1024 Bytes
do i=1,44   ! Ignore initial 44 bytes
   read(IMRC) ibyte
enddo
   read(IMRC) pxsize
!write(*,*) ' pixel size ' , pxsize*1.0E-10 ,' a^-1/pixel'
do i=1,extend_b - 48
   read(IMRC) ibyte
enddo
!!
!write(*,*) 'Starting to read map'
allocate(odata(dims(1), dims(2), dims(3)))
odata = 0.0D0
!
if(mode==1) then                         ! Image is short integer
   allocate(short_map(dims(1), dims(2), dims(3)))
   short_map = 0
   read(IMRC) short_map
   odata = real(short_map, kind=PREC_DP)
elseif(mode==2) then
   read(IMRC) odata
endif
close(imrc)
!
!  Place into global storage
!
call dgl5_new_node
node_number = dgl5_get_number()
nndims = 0
if(dims(3)>1) nndims = nndims + 1
if(dims(2)>1) nndims = nndims + 1
if(dims(1)>1) nndims = nndims + 1
mrc_layer = dims(3)/2
mrc_direct = .false.
mrc_is_grid  = .true.
mrc_has_dxyz = .false.
mrc_has_dval = .false.
!
mrc_corners(1,1) = 1.0         ! Lower left corner
mrc_corners(2,1) = 1.0
mrc_corners(3,1) = 1.0
!
mrc_corners(1,2) = dims(1)     ! Lower right corner
mrc_corners(2,2) = 1.0
mrc_corners(3,2) = 1.0
!
mrc_corners(1,3) = 1.0
mrc_corners(2,3) = dims(2)     ! upper left  corner
mrc_corners(3,3) = 1.0
!
mrc_corners(1,4) = 1.0
mrc_corners(2,4) = 1.0
mrc_corners(3,4) = dims(3)     ! top   left  corner
!
mrc_vectors(1,1) = 1.000       ! Vector along abscissa
mrc_vectors(2,1) = 0.0
mrc_vectors(3,1) = 0.0
!
mrc_vectors(1,2) = 0.0         ! Vector along ordinate
mrc_vectors(2,2) = 1.0
mrc_vectors(3,2) = 0.0
!
mrc_vectors(1,3) = 0.0         ! Vector along top axis
mrc_vectors(2,3) = 0.0
mrc_vectors(3,3) = 1.0
!
mrc_unit(1:3)    = 1.0D0
mrc_unit(4:6)    = 90.0D0
mrc_llims        = 1.0D0
mrc_steps        = 1.0D0
!
mrc_steps_full(1,1) = 1.000       ! Vector along abscissa
mrc_steps_full(2,1) = 0.0
mrc_steps_full(3,1) = 0.0
!
mrc_steps_full(1,2) = 0.0         ! Vector along ordinate
mrc_steps_full(2,2) = 1.0
mrc_steps_full(3,2) = 0.0
!
mrc_steps_full(1,3) = 0.0         ! Vector along top axis
mrc_steps_full(2,3) = 0.0
mrc_steps_full(3,3) = 1.0
!
allocate(mrc_x(1:dims(1)))
allocate(mrc_y(1:dims(2)))
allocate(mrc_z(1:dims(3)))
allocate(mrc_dx(1:dims(1)))
allocate(mrc_dy(1:dims(2)))
allocate(mrc_dz(1:dims(3)))
mrc_dx = 0.0D0
mrc_dy = 0.0D0
mrc_dz = 0.0D0
do i=1, dims(1)
  mrc_x(i) = mrc_llims(1) + (i-1)*mrc_steps_full(1,1)
enddo
do i=1, dims(2)
  mrc_y(i) = mrc_llims(2) + (i-1)*mrc_steps_full(2,2)
enddo
do i=1, dims(3)
  mrc_z(i) = mrc_llims(3) + (i-1)*mrc_steps_full(3,3)
enddo
mrc_pixel_size = pxsize*1.0E-10
rpara(500) = mrc_pixel_size
!
!write(*,*) ' MAKE GLOBAL STORAGE '
call dgl5_set_node(filename , mrc_layer, mrc_direct, nndims,    dims ,         &
                   mrc_is_grid, mrc_has_dxyz, mrc_has_dval, mrc_corners, mrc_vectors,&
                   mrc_unit(1:3), mrc_unit(4:6), mrc_x, mrc_y, mrc_z, mrc_dx, mrc_dy,  &
                   mrc_dz,        odata               ,   osigma, mrc_llims,      &
                   mrc_steps, mrc_steps_full)
!write(*,*) ' DONE GLOBAL STORAGE '
!
deallocate(mrc_x)
deallocate(mrc_y)
deallocate(mrc_z)
!
deallocate(mrc_dx)
deallocate(mrc_dy)
deallocate(mrc_dz)
deallocate(short_map)
deallocate(odata)
!
end subroutine mrc_read
!
!*******************************************************************************
!
end module lib_mrc_mod
