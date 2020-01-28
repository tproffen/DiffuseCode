!------------------------------------------------------------------------------
! NeXus - Neutron & X-ray Common Data Format
!  
! Application Program Interface (Fortran 90)
!
! Copyright (C) 1999-2002, Ray Osborn
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!  For further information, see <http://www.nexusformat.org>
!
!$Id$
!------------------------------------------------------------------------------

MODULE NXmodule

   IMPLICIT NONE

   PUBLIC
! *** NeXus version parameter
   CHARACTER(len=*), PARAMETER, PUBLIC :: NeXus_version = "4.4.0"
! *** NeXus file access parameters
   INTEGER, PARAMETER, PUBLIC :: NXACC_READ = 1
   INTEGER, PARAMETER, PUBLIC :: NXACC_RDWR = 2
   INTEGER, PARAMETER, PUBLIC :: NXACC_CREATE = 3
   INTEGER, PARAMETER, PUBLIC :: NXACC_CREATE4 = 4
   INTEGER, PARAMETER, PUBLIC :: NXACC_CREATE5 = 5
! *** NeXus status parameters
   INTEGER, PARAMETER, PUBLIC :: NX_OK = 1
   INTEGER, PARAMETER, PUBLIC :: NX_ERROR = 0
   INTEGER, PARAMETER, PUBLIC :: NX_EOD = -1
! *** NeXus datatype parameters
   INTEGER, PARAMETER, PUBLIC :: NX_CHAR    = 4
   INTEGER, PARAMETER, PUBLIC :: NX_FLOAT32 = 5
   INTEGER, PARAMETER, PUBLIC :: NX_FLOAT64 = 6
   INTEGER, PARAMETER, PUBLIC :: NX_INT8    = 20
   INTEGER, PARAMETER, PUBLIC :: NX_UINT8   = 21
   INTEGER, PARAMETER, PUBLIC :: NX_INT16   = 22
   INTEGER, PARAMETER, PUBLIC :: NX_UINT16  = 23
   INTEGER, PARAMETER, PUBLIC :: NX_INT32   = 24
   INTEGER, PARAMETER, PUBLIC :: NX_UINT32  = 25
! *** NeXus compression parameters
   INTEGER, PARAMETER, PUBLIC :: NX_COMP_NONE = 100
   INTEGER, PARAMETER, PUBLIC :: NX_COMP_LZW  = 200
   INTEGER, PARAMETER, PUBLIC :: NX_COMP_RLE  = 300
   INTEGER, PARAMETER, PUBLIC :: NX_COMP_HUF  = 400
! *** NeXus Unlimited parameters
   INTEGER, PARAMETER, PUBLIC :: NX_UNLIMITED = -1
! *** NeXus limits
   INTEGER, PARAMETER, PUBLIC :: NX_MAXRANK = 32
   INTEGER, PARAMETER, PUBLIC :: NX_MAXNAMELEN = 64
   INTEGER, PARAMETER, PUBLIC :: NX_MAXSTACK = 20
! *** Kind parameters for different byte lengths (not guaranteed to work)
   INTEGER, PARAMETER, PUBLIC :: NXi1 = selected_int_kind(2)
   INTEGER, PARAMETER, PUBLIC :: NXi2 = selected_int_kind(4)
   INTEGER, PARAMETER, PUBLIC :: NXi4 = selected_int_kind(8)
   INTEGER, PARAMETER, PUBLIC :: NXr4 = kind(1.0)
   INTEGER, PARAMETER, PUBLIC :: NXr8 = kind(1.0D0)
! *** NeXus type definitions
   TYPE, PUBLIC :: NXlink
      INTEGER(kind=NXi4) :: dummy(1040) ! at least as large as in napi.h
   END TYPE
   TYPE, PUBLIC :: NXhandle
      INTEGER(kind=NXi4) :: dummy(9124) ! at least as large as in nxstack.c
   END TYPE
! *** Buffers for each type of parameter
   INTEGER(KIND=NXi1), ALLOCATABLE, PRIVATE :: buffer_i1(:)
   INTEGER(KIND=NXi2), ALLOCATABLE, PRIVATE :: buffer_i2(:)
   INTEGER(KIND=NXi4), ALLOCATABLE, PRIVATE :: buffer_i4(:)
   REAL(KIND=NXr4),    ALLOCATABLE, PRIVATE :: buffer_r4(:)
   REAL(KIND=NXr8),    ALLOCATABLE, PRIVATE :: buffer_r8(:)
   INTEGER, PRIVATE :: NXrank, NXdims(NX_MAXRANK), NXtype, NXsize
! *** NeXus core functions ***
   PUBLIC :: NXopen, NXclose, NXflush
   PUBLIC :: NXmakegroup, NXopengroup, NXclosegroup
   PUBLIC :: NXmakedata, NXopendata, NXcompress, NXclosedata
   PUBLIC :: NXgetdata, NXgetslab, NXgetattr, NXputdata, NXputslab, NXputattr
   PUBLIC :: NXgetinfo, NXgetnextentry, NXgetnextattr
   PUBLIC :: NXgetgroupID, NXgetdataID, NXsameID, NXmakelink 
   PUBLIC :: NXgetgroupinfo, NXinitgroupdir, NXgroupdir
   PUBLIC :: NXgetattrinfo, NXinitattrdir, NXattrdir
   PUBLIC :: NXreverse, NXCstring, NXFstring, NXdatatype, NXerror
! *** NeXus generic interfaces ***
   INTERFACE NXgetdata
      MODULE PROCEDURE NXgeti1, NXgeti2, NXgeti4, NXgetr4, NXgetr8, NXgetchar
   END INTERFACE
   INTERFACE NXgetslab
      MODULE PROCEDURE NXgeti1slab, NXgeti2slab, NXgeti4slab, &
                        NXgetr4slab, NXgetr8slab
   END INTERFACE
   INTERFACE NXgetattr
      MODULE PROCEDURE NXgeti1attr, NXgeti2attr, NXgeti4attr, NXgetr4attr, &
                        NXgetr8attr, NXgetcharattr
   END INTERFACE
   INTERFACE NXputdata
      MODULE PROCEDURE NXputi1, NXputi2, NXputi4, NXputr4, NXputr8, NXputchar
   END INTERFACE
   INTERFACE NXputslab
      MODULE PROCEDURE NXputi1slab, NXputi2slab, NXputi4slab, &
                        NXputr4slab, NXputr8slab
   END INTERFACE
   INTERFACE NXputattr
      MODULE PROCEDURE NXputi1attr, NXputi2attr, NXputi4attr,  &
                        NXputr4attr, NXputr8attr, NXputcharattr
   END INTERFACE

CONTAINS
!------------------------------------------------------------------------------
!NXopen opens a NeXus file and returns a file ID 
   FUNCTION NXopen (file_name, access_method, file_id) RESULT (status)

      CHARACTER(len=*), INTENT(in)  :: file_name
      INTEGER,          INTENT(in)  :: access_method
      TYPE(NXhandle),   INTENT(out) :: file_id
      INTEGER :: status, nxifopen
      EXTERNAL nxifopen

      status = nxifopen (NXCstring(file_name), access_method, file_id)

   END FUNCTION NXopen
!------------------------------------------------------------------------------
!NXclose closes a NeXus file defined by its file ID
   FUNCTION NXclose (file_id) RESULT (status)

      TYPE(NXhandle), INTENT(in) :: file_id
      INTEGER :: status, nxifclose
      EXTERNAL nxifclose

      status = nxifclose (file_id)

   END FUNCTION NXclose
!------------------------------------------------------------------------------
!NXflush flushes all pending data to disk
   FUNCTION NXflush (file_id) RESULT (status)

      TYPE(NXhandle), INTENT(inout) :: file_id
      INTEGER :: status, nxifflush
      EXTERNAL nxifflush

      status = nxifflush (file_id)

   END FUNCTION NXflush
!------------------------------------------------------------------------------
!NXmakegroup creates a NeXus group
   FUNCTION NXmakegroup (file_id, group_name, group_class) RESULT (status)

      TYPE(NXhandle),   INTENT(in) :: file_id
      CHARACTER(len=*), INTENT(in) :: group_name, group_class
      INTEGER :: status, nximakegroup
      EXTERNAL nximakegroup

      status = nximakegroup(file_id, NXCstring(group_name), &
                        NXCstring(group_class))

   END FUNCTION NXmakegroup
!------------------------------------------------------------------------------
!NXopengroup opens an existing NeXus group for input/output
   FUNCTION NXopengroup (file_id, group_name, group_class) RESULT (status)

      TYPE(NXhandle),   INTENT(in) :: file_id
      CHARACTER(len=*), INTENT(in) :: group_name, group_class
      INTEGER :: status, nxiopengroup
      EXTERNAL nxiopengroup

      status = nxiopengroup(file_id, NXCstring(group_name), &
                        NXCstring(group_class))

   END FUNCTION NXopengroup
!------------------------------------------------------------------------------
!NXclosegroup closes a NeXus group
   FUNCTION NXclosegroup (file_id) RESULT (status)

      TYPE(NXhandle), INTENT(in) :: file_id
      INTEGER :: status, nxiclosegroup
      EXTERNAL nxiclosegroup

      status = nxiclosegroup(file_id)

   END FUNCTION NXclosegroup
!------------------------------------------------------------------------------
!NXmakedata creates a NeXus data set (optionally with compression)
   FUNCTION NXmakedata (file_id, data_name, data_type, data_rank, &
                        data_dimensions, compress_type, chunk_size) &
                        RESULT (status)

      TYPE(NXhandle),   INTENT(in) :: file_id
      CHARACTER(len=*), INTENT(in) :: data_name
      INTEGER,          INTENT(in) :: data_type,data_rank,data_dimensions(:)
      INTEGER, OPTIONAL,INTENT(in) :: compress_type, chunk_size(:)
      INTEGER, ALLOCATABLE :: NXchunk_size(:)
      INTEGER :: status, i, nxifmakedata, nxifcompmakedata
      EXTERNAL nxifmakedata, nxifcompmakedata

      IF (PRESENT(compress_type)) THEN
         IF (PRESENT(chunk_size)) THEN
            ALLOCATE (NXchunk_size(data_rank))
            NXchunk_size = chunk_size
         ELSE
            ALLOCATE (NXchunk_size(data_rank))
            NXchunk_size = (/(data_dimensions(i),i=1,data_rank)/)
         END IF
         status = nxifcompmakedata(file_id, NXCstring(data_name), data_type, &
                        data_rank, data_dimensions, compress_type, NXchunk_size)
         DEALLOCATE (NXchunk_size)
      ELSE
         status = nxifmakedata(file_id, NXCstring(data_name), data_type, &
                        data_rank, data_dimensions)
      END IF

   END FUNCTION NXmakedata
!------------------------------------------------------------------------------
!NXopendata opens an existing NeXus data set for input/output
   FUNCTION NXopendata (file_id, data_name) RESULT (status)

      TYPE(NXhandle),   INTENT(in) :: file_id
      CHARACTER(len=*), INTENT(in) :: data_name
      INTEGER :: status, nxiopendata
      EXTERNAL nxiopendata

      status = nxiopendata(file_id, NXCstring(data_name))

   END FUNCTION NXopendata
!------------------------------------------------------------------------------
!NXcompress sets the compression algorithm for the open NeXus data set
   FUNCTION NXcompress (file_id, compress_type) RESULT (status)

      TYPE(NXhandle),   INTENT(in) :: file_id
      INTEGER,          INTENT(in) :: compress_type
      INTEGER :: status, nxifcompress
      EXTERNAL nxifcompress

      status = nxifcompress(file_id, compress_type)

   END FUNCTION NXcompress
!------------------------------------------------------------------------------
!NXclosedata closes a NeXus data set
   FUNCTION NXclosedata (file_id) RESULT (status)

      TYPE(NXhandle), INTENT(in) :: file_id
      INTEGER :: status, nxiclosedata
      EXTERNAL nxiclosedata

      status = nxiclosedata(file_id)

   END FUNCTION NXclosedata
!------------------------------------------------------------------------------
!NXgetdata reads data from the open data set
!
!The following routines define the generic function NXgetdata
!------------------------------------------------------------------------------
!NXgeti1 reads an integer*1 array from the open data set
   FUNCTION NXgeti1 (file_id, data) RESULT (status)

      TYPE(NXhandle),     INTENT(in)  :: file_id
      INTEGER(KIND=NXi1), INTENT(out) :: data(:)
      INTEGER :: status, nxigetdata
      EXTERNAL nxigetdata

      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      IF (status /= NX_OK) RETURN
      NXsize = PRODUCT(NXdims(1:NXrank))
      IF (NXsize > size(data)) THEN
         CALL NXerror ("The supplied array is not large enough for the data")
         status = NX_ERROR
      ELSE IF (NXtype == NX_INT8 .OR. NXtype == NX_UINT8) THEN
         ALLOCATE (buffer_i1(NXsize))
         status = nxigetdata(file_id, buffer_i1)
         data = buffer_i1
         DEALLOCATE (buffer_i1)
      ELSE IF (NXtype == NX_INT16 .OR. NXtype == NX_UINT16) THEN
         ALLOCATE (buffer_i2(NXsize))
         status = nxigetdata(file_id, buffer_i2)
         IF (abs(maxval(buffer_i2)) <= HUGE(data)) THEN
            data = buffer_i2
         ELSE
            CALL NXerror ("Input values too large for data type")
            status = NX_ERROR
         END IF
         DEALLOCATE (buffer_i2)
      ELSE IF (NXtype == NX_INT32 .OR. NXtype == NX_UINT32) THEN
         ALLOCATE (buffer_i4(NXsize))
         status = nxigetdata(file_id, buffer_i4)
         IF (abs(maxval(buffer_i4)) <= HUGE(data)) THEN
            data = buffer_i4
         ELSE
            CALL NXerror ("Input values too large for data type")
            status = NX_ERROR
         END IF
         DEALLOCATE (buffer_i4)
      ELSE
         call NXerror &
              ("The datatype is incompatible with the supplied variable")
         status = NX_ERROR
      END IF

   END FUNCTION NXgeti1
!------------------------------------------------------------------------------
!NXgeti2 reads an integer*2 array from the open data set
   FUNCTION NXgeti2 (file_id, data) RESULT (status)

      TYPE(NXhandle),     INTENT(in)  :: file_id
      INTEGER(KIND=NXi2), INTENT(out) :: data(:)
      INTEGER :: status, nxigetdata
      EXTERNAL nxigetdata

      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      IF (status /= NX_OK) RETURN
      NXsize = PRODUCT(NXdims(1:NXrank))
      IF (NXsize > size(data)) THEN
         CALL NXerror ("The supplied array is not large enough for the data")
         status = NX_ERROR
      ELSE IF (NXtype == NX_INT8 .OR. NXtype == NX_UINT8) THEN
         ALLOCATE (buffer_i1(NXsize))
         status = nxigetdata(file_id, buffer_i1)
         data = buffer_i1
         DEALLOCATE (buffer_i1)
      ELSE IF (NXtype == NX_INT16 .OR. NXtype == NX_UINT16) THEN
         ALLOCATE (buffer_i2(NXsize))
         status = nxigetdata(file_id, buffer_i2)
         data = buffer_i2
         DEALLOCATE (buffer_i2)
      ELSE IF (NXtype == NX_INT32 .OR. NXtype == NX_UINT32) THEN
         ALLOCATE (buffer_i4(NXsize))
         status = nxigetdata(file_id, buffer_i4)
         IF (abs(maxval(buffer_i4)) <= HUGE(data)) THEN
            data = buffer_i4
         ELSE
            CALL NXerror ("Input values too large for data type")
            status = NX_ERROR
         END IF
         DEALLOCATE (buffer_i4)
      ELSE
         call NXerror &
              ("The datatype is incompatible with the supplied variable")
         status = NX_ERROR
      END IF

   END FUNCTION NXgeti2
!------------------------------------------------------------------------------
!NXgeti4 reads an integer*4 array from the open data set
   FUNCTION NXgeti4 (file_id, data) RESULT (status)

      TYPE(NXhandle),     INTENT(in)  :: file_id
      INTEGER(KIND=NXi4), INTENT(out) :: data(:)
      INTEGER :: status, nxigetdata
      EXTERNAL nxigetdata

      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      IF (status /= NX_OK) RETURN
      NXsize = PRODUCT(NXdims(1:NXrank))
      IF (NXsize > size(data)) THEN
         CALL NXerror ("The supplied array is not large enough for the data")
         status = NX_ERROR
      ELSE IF (NXtype == NX_INT8 .OR. NXtype == NX_UINT8) THEN
         ALLOCATE (buffer_i1(NXsize))
         status = nxigetdata(file_id, buffer_i1)
         data = buffer_i1
         DEALLOCATE (buffer_i1)
      ELSE IF (NXtype == NX_INT16 .OR. NXtype == NX_UINT16) THEN
         ALLOCATE (buffer_i2(NXsize))
         status = nxigetdata(file_id, buffer_i2)
         data = buffer_i2
         DEALLOCATE (buffer_i2)
      ELSE IF (NXtype == NX_INT32 .OR. NXtype == NX_UINT32) THEN
         ALLOCATE (buffer_i4(NXsize))
         status = nxigetdata(file_id, buffer_i4)
         data = buffer_i4
         DEALLOCATE (buffer_i4)
      ELSE
         call NXerror &
              ("The datatype is incompatible with the supplied variable")
         status = NX_ERROR
      END IF

   END FUNCTION NXgeti4
!------------------------------------------------------------------------------
!NXgetr4 reads a real*4 array from the open data set
   FUNCTION NXgetr4 (file_id, data) RESULT (status)

      TYPE(NXhandle),  INTENT(in)  :: file_id
      REAL(KIND=NXr4), INTENT(out) :: data(:)
      INTEGER :: status, nxigetdata
      EXTERNAL nxigetdata

      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      IF (status /= NX_OK) RETURN
      NXsize = PRODUCT(NXdims(1:NXrank))
      IF (NXsize > size(data)) THEN
         CALL NXerror ("The supplied array is not large enough for the data")
         status = NX_ERROR
      ELSE IF (NXtype == NX_FLOAT32) THEN
         ALLOCATE (buffer_r4(NXsize))
         status = nxigetdata(file_id, buffer_r4)
         data = buffer_r4
         DEALLOCATE (buffer_r4)
      ELSE IF (NXtype == NX_FLOAT64) THEN
         ALLOCATE (buffer_r8(NXsize))
         status = nxigetdata(file_id, buffer_r8)
         data = buffer_r8
         DEALLOCATE (buffer_r8)
      ELSE
         call NXerror &
              ("The datatype is incompatible with the supplied variable")
         status = NX_ERROR
      END IF

   END FUNCTION NXgetr4
!------------------------------------------------------------------------------
!NXgetr8 reads a real*8 array from the open data set
   FUNCTION NXgetr8 (file_id, data) RESULT (status)

      TYPE(NXhandle),  INTENT(in)  :: file_id
      REAL(KIND=NXr8), INTENT(out) :: data(:)
      INTEGER :: status, nxigetdata
      EXTERNAL nxigetdata

      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      IF (status /= NX_OK) RETURN
      NXsize = PRODUCT(NXdims(1:NXrank))
      IF (NXsize > size(data)) THEN
         CALL NXerror ("The supplied array is not large enough for the data")
         status = NX_ERROR
      ELSE IF (NXtype == NX_FLOAT32) THEN
         ALLOCATE (buffer_r4(NXsize))
         status = nxigetdata(file_id, buffer_r4)
         data = buffer_r4
         DEALLOCATE (buffer_r4)
      ELSE IF (NXtype == NX_FLOAT64) THEN
         ALLOCATE (buffer_r8(NXsize))
         status = nxigetdata(file_id, buffer_r8)
         IF (abs(maxval(buffer_r8)) <= HUGE(data)) THEN
            data = buffer_r8
         ELSE
            CALL NXerror ("Input values too large for data type")
            status = NX_ERROR
         END IF
         DEALLOCATE (buffer_r8)
      ELSE
         call NXerror &
              ("The datatype is incompatible with the supplied variable")
         status = NX_ERROR
      END IF

   END FUNCTION NXgetr8
!------------------------------------------------------------------------------
!NXgetchar reads a character string from the open data set
   FUNCTION NXgetchar (file_id, data) RESULT (status)

      TYPE(NXhandle),   INTENT(in)  :: file_id
      CHARACTER(len=*), INTENT(out) :: data
      INTEGER :: status, nxigetdata
      INTEGER(kind=NXi1) :: Cstring(255)
      EXTERNAL nxigetdata

      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      IF (status /= NX_OK) RETURN
      NXsize = PRODUCT(NXdims(1:NXrank))
      IF (NXsize > len(data)) THEN
         CALL NXerror ("The supplied string is not large enough for the data")
         status = NX_ERROR
      ELSE IF (NXtype == NX_CHAR) THEN
         Cstring = 0 !HDF does not add null termination so ensure it's there
         status = nxigetdata(file_id, Cstring)
         IF (status == NX_OK) data = trim(NXFstring(Cstring))
      ELSE
         call NXerror &
              ("The datatype is incompatible with the supplied variable")
         status = NX_ERROR
      END IF

   END FUNCTION NXgetchar
!------------------------------------------------------------------------------
!NXgetslab reads a slab of the open data set
!
!The following routines define the generic function NXgetslab
!------------------------------------------------------------------------------
!NXgeti1slab reads a slab of integer*1 data from the open data set 
   FUNCTION NXgeti1slab (file_id, data, data_start, data_size) RESULT (status)

      TYPE(NXhandle),     INTENT(in)  :: file_id
      INTEGER,            INTENT(in)  :: data_start(:), data_size(:)
      INTEGER(KIND=NXi1), INTENT(out) :: data(:)
      INTEGER :: status, nxigetslab
      EXTERNAL nxigetslab

      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      IF (status /= NX_OK) RETURN
      NXsize = PRODUCT(data_size(1:NXrank))
      IF (NXsize > size(data)) THEN
         CALL NXerror ("The supplied array is not large enough for the data")
         status = NX_ERROR
      ELSE IF (NXtype == NX_INT8 .OR. NXtype == NX_UINT8) THEN
         ALLOCATE (buffer_i1(NXsize))
         status = nxigetslab(file_id, buffer_i1, &
                   NXreverse(NXrank,data_start)-1, NXreverse(NXrank,data_size))
         data = buffer_i1
         DEALLOCATE (buffer_i1)
      ELSE IF (NXtype == NX_INT16 .OR. NXtype == NX_UINT16) THEN
         ALLOCATE (buffer_i2(NXsize))
         status = nxigetslab(file_id, buffer_i2, &
                   NXreverse(NXrank,data_start)-1, NXreverse(NXrank,data_size))
         IF (abs(maxval(buffer_i2)) <= HUGE(data)) THEN
            data = buffer_i2
         ELSE
            CALL NXerror ("Input values too large for data type")
            status = NX_ERROR
         END IF
         DEALLOCATE (buffer_i2)
      ELSE IF (NXtype == NX_INT32 .OR. NXtype == NX_UINT32) THEN
         ALLOCATE (buffer_i4(NXsize))
         status = nxigetslab(file_id, buffer_i4, &
                   NXreverse(NXrank,data_start)-1, NXreverse(NXrank,data_size))
         IF (abs(maxval(buffer_i4)) <= HUGE(data)) THEN
            data = buffer_i4
         ELSE
            CALL NXerror ("Input values too large for data type")
            status = NX_ERROR
         END IF
         DEALLOCATE (buffer_i4)
      ELSE
         call NXerror &
              ("The datatype is incompatible with the supplied variable")
         status = NX_ERROR
      END IF

   END FUNCTION NXgeti1slab
!------------------------------------------------------------------------------
!NXgeti2slab reads a slab of integer*2 data from the open data set 
   FUNCTION NXgeti2slab (file_id, data, data_start, data_size) RESULT (status)

      TYPE(NXhandle),     INTENT(in)  :: file_id
      INTEGER,            INTENT(in)  :: data_start(:), data_size(:)
      INTEGER(KIND=NXi2), INTENT(out) :: data(:)
      INTEGER :: status, nxigetslab
      EXTERNAL nxigetslab

      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      IF (status /= NX_OK) RETURN
      NXsize = PRODUCT(data_size(1:NXrank))
      IF (NXsize > size(data)) THEN
         CALL NXerror ("The supplied array is not large enough for the data")
         status = NX_ERROR
      ELSE IF (NXtype == NX_INT8 .OR. NXtype == NX_UINT8) THEN
         ALLOCATE (buffer_i1(NXsize))
         status = nxigetslab(file_id, buffer_i1, &
                   NXreverse(NXrank,data_start)-1, NXreverse(NXrank,data_size))
         data = buffer_i1
         DEALLOCATE (buffer_i1)
      ELSE IF (NXtype == NX_INT16 .OR. NXtype == NX_UINT16) THEN
         ALLOCATE (buffer_i2(NXsize))
         status = nxigetslab(file_id, buffer_i2, &
                   NXreverse(NXrank,data_start)-1, NXreverse(NXrank,data_size))
         data = buffer_i2
         DEALLOCATE (buffer_i2)
      ELSE IF (NXtype == NX_INT32 .OR. NXtype == NX_UINT32) THEN
         ALLOCATE (buffer_i4(NXsize))
         status = nxigetslab(file_id, buffer_i4, &
                   NXreverse(NXrank,data_start)-1, NXreverse(NXrank,data_size))
         IF (abs(maxval(buffer_i4)) <= HUGE(data)) THEN
            data = buffer_i4
         ELSE
            CALL NXerror ("Input values too large for data type")
            status = NX_ERROR
         END IF
         DEALLOCATE (buffer_i4)
      ELSE
         call NXerror &
              ("The datatype is incompatible with the supplied variable")
         status = NX_ERROR
      END IF

   END FUNCTION NXgeti2slab
!------------------------------------------------------------------------------
!NXgeti4slab reads a slab of integer*4 data from the open data set 
   FUNCTION NXgeti4slab (file_id, data, data_start, data_size) RESULT (status)

      TYPE(NXhandle),     INTENT(in)  :: file_id
      INTEGER,            INTENT(in)  :: data_start(:), data_size(:)
      INTEGER(KIND=NXi4), INTENT(out) :: data(:)
      INTEGER :: status, nxigetslab
      EXTERNAL nxigetslab

      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      IF (status /= NX_OK) RETURN
      NXsize = PRODUCT(data_size(1:NXrank))
      IF (NXsize > size(data)) THEN
         CALL NXerror ("The supplied array is not large enough for the data")
         status = NX_ERROR
      ELSE IF (NXtype == NX_INT8 .OR. NXtype == NX_UINT8) THEN
         ALLOCATE (buffer_i1(NXsize))
         status = nxigetslab(file_id, buffer_i1, &
                   NXreverse(NXrank,data_start)-1, NXreverse(NXrank,data_size))
         data = buffer_i1
         DEALLOCATE (buffer_i1)
      ELSE IF (NXtype == NX_INT16 .OR. NXtype == NX_UINT16) THEN
         ALLOCATE (buffer_i2(NXsize))
         status = nxigetslab(file_id, buffer_i2, &
                   NXreverse(NXrank,data_start)-1, NXreverse(NXrank,data_size))
         data = buffer_i2
         DEALLOCATE (buffer_i2)
      ELSE IF (NXtype == NX_INT32 .OR. NXtype == NX_UINT32) THEN
         ALLOCATE (buffer_i4(NXsize))
         status = nxigetslab(file_id, buffer_i4, &
                   NXreverse(NXrank,data_start)-1, NXreverse(NXrank,data_size))
         data = buffer_i4
         DEALLOCATE (buffer_i4)
      ELSE
         call NXerror &
              ("The datatype is incompatible with the supplied variable")
         status = NX_ERROR
      END IF

   END FUNCTION NXgeti4slab
!------------------------------------------------------------------------------
!NXgetr4slab reads a slab of real*4 data from the open data set 
   FUNCTION NXgetr4slab (file_id, data, data_start, data_size) RESULT (status)

      TYPE(NXhandle),  INTENT(in)  :: file_id
      INTEGER,         INTENT(in)  :: data_start(:), data_size(:)
      REAL(KIND=NXr4), INTENT(out) :: data(:)
      INTEGER :: status, nxigetslab
      EXTERNAL nxigetslab

      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      IF (status /= NX_OK) RETURN
      NXsize = PRODUCT(data_size(1:NXrank))
      IF (NXsize > size(data)) THEN
         CALL NXerror ("The supplied array is not large enough for the data")
         status = NX_ERROR
      ELSE IF (NXtype == NX_FLOAT32) THEN
         ALLOCATE (buffer_r4(NXsize))
         status = nxigetslab(file_id, buffer_r4, &
                   NXreverse(NXrank,data_start)-1, NXreverse(NXrank,data_size))
         data = buffer_r4
         DEALLOCATE (buffer_r4)
      ELSE IF (NXtype == NX_FLOAT64) THEN
         ALLOCATE (buffer_r8(NXsize))
         status = nxigetslab(file_id, buffer_r8, &
                   NXreverse(NXrank,data_start)-1, NXreverse(NXrank,data_size))
         IF (abs(maxval(buffer_r8)) <= HUGE(data)) THEN
            data = buffer_r8
         ELSE
            CALL NXerror ("Input values too large for data type")
            status = NX_ERROR
         END IF
         DEALLOCATE (buffer_r8)
      ELSE
         call NXerror &
              ("The datatype is incompatible with the supplied variable")
         status = NX_ERROR
      END IF

   END FUNCTION NXgetr4slab
!------------------------------------------------------------------------------
!NXgetr8slab reads a slab of real*8 data from the open data set 
   FUNCTION NXgetr8slab (file_id, data, data_start, data_size) RESULT (status)

      TYPE(NXhandle),  INTENT(in)  :: file_id
      INTEGER,         INTENT(in)  :: data_start(:), data_size(:)
      REAL(KIND=NXr8), INTENT(out) :: data(:)
      INTEGER :: status, nxigetslab
      EXTERNAL nxigetslab

      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      IF (status /= NX_OK) RETURN
      NXsize = PRODUCT(data_size(1:NXrank))
      IF (NXsize > size(data)) THEN
         CALL NXerror ("The supplied array is not large enough for the data")
         status = NX_ERROR
      ELSE IF (NXtype == NX_FLOAT32) THEN
         ALLOCATE (buffer_r4(NXsize))
         status = nxigetslab(file_id, buffer_r4, &
                   NXreverse(NXrank,data_start), NXreverse(NXrank,data_size))
         data = buffer_r4
         DEALLOCATE (buffer_r4)
      ELSE IF (NXtype == NX_FLOAT64) THEN
         ALLOCATE (buffer_r8(NXsize))
         status = nxigetslab(file_id, buffer_r8, &
                   NXreverse(NXrank,data_start), NXreverse(NXrank,data_size))
         data = buffer_r8
         DEALLOCATE (buffer_r8)
      ELSE
         call NXerror &
              ("The datatype is incompatible with the supplied variable")
         status = NX_ERROR
      END IF

   END FUNCTION NXgetr8slab
!------------------------------------------------------------------------------
!NXgetattr reads attributes from the open data set 
!
!The following routines define the generic function NXgetattr
!------------------------------------------------------------------------------
!NXgeti1attr reads an integer*1 attribute from the open data set 
   FUNCTION NXgeti1attr (file_id, attr_name, value, attr_length, attr_type) &
                        RESULT (status)

      TYPE(NXhandle),     INTENT(in)    :: file_id
      CHARACTER(len=*),   INTENT(in)    :: attr_name
      INTEGER(KIND=NXi1), INTENT(out)   :: value
      INTEGER, OPTIONAL,  INTENT(inout) :: attr_length
      INTEGER, OPTIONAL,  INTENT(in)    :: attr_type
      INTEGER :: status, nxigetattr, value_length, value_type
      EXTERNAL nxigetattr

      value_length = 1; value_type = NX_INT8
      status = nxigetattr(file_id, NXCstring(attr_name), value, value_length, &
                        value_type)

   END FUNCTION NXgeti1attr
!------------------------------------------------------------------------------
!NXgeti2attr reads an integer*2 attribute from the open data set 
   FUNCTION NXgeti2attr (file_id, attr_name, value, attr_length, attr_type) &
                        RESULT (status)

      TYPE(NXhandle),     INTENT(in)    :: file_id
      CHARACTER(len=*),   INTENT(in)    :: attr_name
      INTEGER(KIND=NXi2), INTENT(out)   :: value
      INTEGER, OPTIONAL,  INTENT(inout) :: attr_length
      INTEGER, OPTIONAL,  INTENT(in)    :: attr_type
      INTEGER :: status, nxigetattr, value_length, value_type
      EXTERNAL nxigetattr

      value_length = 1; value_type = NX_INT16
      status = nxigetattr(file_id, NXCstring(attr_name), value, value_length, &
                        value_type)

   END FUNCTION NXgeti2attr
!------------------------------------------------------------------------------
!NXgeti4attr reads an integer*4 attribute from the open data set 
   FUNCTION NXgeti4attr (file_id, attr_name, value, attr_length, attr_type) &
                        RESULT (status)

      TYPE(NXhandle),     INTENT(in)    :: file_id
      CHARACTER(len=*),   INTENT(in)    :: attr_name
      INTEGER(KIND=NXi4), INTENT(out)   :: value
      INTEGER, OPTIONAL,  INTENT(inout) :: attr_length
      INTEGER, OPTIONAL,  INTENT(in)    :: attr_type
      INTEGER :: status, nxigetattr, value_length, value_type
      EXTERNAL nxigetattr

      value_length = 1; value_type = NX_INT32
      status = nxigetattr(file_id, NXCstring(attr_name), value, value_length, &
                        value_type)

   END FUNCTION NXgeti4attr
!------------------------------------------------------------------------------
!NXgetr4attr reads a real*4 attribute from the open data set 
   FUNCTION NXgetr4attr (file_id, attr_name, value, attr_length, attr_type) &
                        RESULT (status)

      TYPE(NXhandle),     INTENT(in)    :: file_id
      CHARACTER(len=*),   INTENT(in)    :: attr_name
      REAL(KIND=NXr4),    INTENT(out)   :: value
      INTEGER, OPTIONAL,  INTENT(inout) :: attr_length
      INTEGER, OPTIONAL,  INTENT(in)    :: attr_type
      INTEGER :: status, nxigetattr, value_length, value_type
      EXTERNAL nxigetattr

      value_length = 1; value_type = NX_FLOAT32
      status = nxigetattr(file_id, NXCstring(attr_name), value, value_length, &
                        value_type)

   END FUNCTION NXgetr4attr
!------------------------------------------------------------------------------
!NXgetr8attr reads a real*8 attribute from the open data set 
   FUNCTION NXgetr8attr (file_id, attr_name, value, attr_length, attr_type) &
                        RESULT (status)

      TYPE(NXhandle),     INTENT(in)    :: file_id
      CHARACTER(len=*),   INTENT(in)    :: attr_name
      REAL(KIND=NXr8),    INTENT(out)   :: value
      INTEGER, OPTIONAL,  INTENT(inout) :: attr_length
      INTEGER, OPTIONAL,  INTENT(in)    :: attr_type
      INTEGER :: status, nxigetattr, value_length, value_type
      EXTERNAL nxigetattr

      value_length = 1; value_type = NX_FLOAT64
      status = nxigetattr(file_id, NXCstring(attr_name), value, value_length, &
                        value_type)

   END FUNCTION NXgetr8attr
!------------------------------------------------------------------------------
!NXgetcharattr reads a character attribute from the open data set 
   FUNCTION NXgetcharattr (file_id, attr_name, value, attr_length, attr_type) &
                        RESULT (status)

      TYPE(NXhandle),     INTENT(in)    :: file_id
      CHARACTER(len=*),   INTENT(in)    :: attr_name
      CHARACTER(len=*),   INTENT(out)   :: value
      INTEGER, OPTIONAL,  INTENT(inout) :: attr_length
      INTEGER, OPTIONAL,  INTENT(in)    :: attr_type
      INTEGER :: status, nxigetattr, value_length, value_type
      INTEGER(kind=NXi1) :: Cstring(255)
      EXTERNAL nxigetattr

      value_length = len(value); value_type = NX_CHAR
      Cstring = 0
      status = nxigetattr(file_id, NXCstring(attr_name), Cstring, &
                        value_length, value_type)
      value = trim(NXFstring(Cstring))

   END FUNCTION NXgetcharattr
!------------------------------------------------------------------------------
!NXputdata writes data into the open data set
!
!The following routines define the generic function NXputdata
!------------------------------------------------------------------------------
!NXputi1 writes an integer*1 array to the open data set
   FUNCTION NXputi1 (file_id, data) RESULT (status)

      TYPE(NXhandle),     INTENT(in) :: file_id
      INTEGER(KIND=NXi1), INTENT(in) :: data(:)
      INTEGER :: status, nxiputdata
      EXTERNAL nxiputdata

      status = nxiputdata(file_id, data)

   END FUNCTION NXputi1
!------------------------------------------------------------------------------
!NXputi2 writes an integer*2 array to the open data set
   FUNCTION NXputi2 (file_id, data) RESULT (status)

      TYPE(NXhandle),     INTENT(in) :: file_id
      INTEGER(KIND=NXi2), INTENT(in) :: data(:)
      INTEGER :: status, nxiputdata
      EXTERNAL nxiputdata

      status = nxiputdata(file_id, data)

   END FUNCTION NXputi2
!------------------------------------------------------------------------------
!NXputi1 writes an integer*4 array to the open data set
   FUNCTION NXputi4 (file_id, data) RESULT (status)

      TYPE(NXhandle),     INTENT(in) :: file_id
      INTEGER(KIND=NXi4), INTENT(in) :: data(:)
      INTEGER :: status, nxiputdata
      EXTERNAL nxiputdata

      status = nxiputdata(file_id, data)

   END FUNCTION NXputi4
!------------------------------------------------------------------------------
!NXputreal writes a real*4 array to the open data set
   FUNCTION NXputr4 (file_id, data) RESULT (status)

      TYPE(NXhandle),  INTENT(in) :: file_id
      REAL(KIND=NXr4), INTENT(in) :: data(:)
      INTEGER :: status, nxiputdata
      EXTERNAL nxiputdata

      status = nxiputdata(file_id, data)

   END FUNCTION NXputr4
!------------------------------------------------------------------------------
!NXputr8 writes a real*8 array to the open data set
   FUNCTION NXputr8 (file_id, data) RESULT (status)

      TYPE(NXhandle),  INTENT(in) :: file_id
      REAL(KIND=NXr8), INTENT(in) :: data(:)
      INTEGER :: status, nxiputdata
      EXTERNAL nxiputdata

      status = nxiputdata(file_id, data)

   END FUNCTION NXputr8
!------------------------------------------------------------------------------
!NXputchar writes a character string to the open data set
   FUNCTION NXputchar (file_id, data) RESULT (status)

      TYPE(NXhandle),   INTENT(in) :: file_id
      CHARACTER(len=*), INTENT(in) :: data
      INTEGER :: status, nxiputdata
      EXTERNAL nxiputdata

      status = nxiputdata(file_id, NXCstring(data))

   END FUNCTION NXputchar
!------------------------------------------------------------------------------
!NXputslab writes a slab of data into the open data set
!
!The following routines define the generic function NXputslab
!------------------------------------------------------------------------------
!NXputi1slab writes a slab of integer*1 data into the open data set 
   FUNCTION NXputi1slab (file_id, data, data_start, data_size) RESULT (status)

      TYPE(NXhandle),     INTENT(in) :: file_id
      INTEGER,            INTENT(in) :: data_start(:), data_size(:)
      INTEGER(KIND=NXi1), INTENT(in) :: data(:)
      INTEGER :: status, nxiputslab
      EXTERNAL nxiputslab

      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      status = nxiputslab(file_id, data, &
                   NXreverse(NXrank,data_start)-1, NXreverse(NXrank,data_size))

   END FUNCTION NXputi1slab
!------------------------------------------------------------------------------
!NXputi2slab writes a slab of integer*2 data into the open data set 
   FUNCTION NXputi2slab (file_id, data, data_start, data_size) RESULT (status)

      TYPE(NXhandle),     INTENT(in) :: file_id
      INTEGER,            INTENT(in) :: data_start(:), data_size(:)
      INTEGER(KIND=NXi2), INTENT(in) :: data(:)
      INTEGER :: status, nxiputslab
      EXTERNAL nxiputslab

      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      status = nxiputslab(file_id, data, &
                   NXreverse(NXrank,data_start)-1, NXreverse(NXrank,data_size))

   END FUNCTION NXputi2slab
!------------------------------------------------------------------------------
!NXputi4slab writes a slab of integer*4 data into the open data set 
   FUNCTION NXputi4slab (file_id, data, data_start, data_size) RESULT (status)

      TYPE(NXhandle),     INTENT(in) :: file_id
      INTEGER,            INTENT(in) :: data_start(:), data_size(:)
      INTEGER(KIND=NXi4), INTENT(in) :: data(:)
      INTEGER :: status, nxiputslab
      EXTERNAL nxiputslab

      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      status = nxiputslab(file_id, data, &
                   NXreverse(NXrank,data_start)-1, NXreverse(NXrank,data_size))

   END FUNCTION NXputi4slab
!------------------------------------------------------------------------------
!NXputr4slab writes a slab of real*4 data into the open data set 
   FUNCTION NXputr4slab (file_id, data, data_start, data_size) RESULT (status)

      TYPE(NXhandle),  INTENT(in) :: file_id
      INTEGER,         INTENT(in) :: data_start(:), data_size(:)
      REAL(KIND=NXr4), INTENT(in) :: data(:)
      INTEGER :: status, nxiputslab
      EXTERNAL nxiputslab

      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      status = nxiputslab(file_id, data, &
                   NXreverse(NXrank,data_start)-1, NXreverse(NXrank,data_size))

   END FUNCTION NXputr4slab
!------------------------------------------------------------------------------
!NXputr8slab writes a slab of real*8 data into the open data set 
   FUNCTION NXputr8slab (file_id, data, data_start, data_size) RESULT (status)

      TYPE(NXhandle),  INTENT(in) :: file_id
      INTEGER,         INTENT(in) :: data_start(:), data_size(:)
      REAL(KIND=NXr8), INTENT(in) :: data(:)
      INTEGER :: status, nxiputslab
      EXTERNAL nxiputslab

      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      status = nxiputslab(file_id, data, &
                   NXreverse(NXrank,data_start)-1, NXreverse(NXrank,data_size))

   END FUNCTION NXputr8slab
!------------------------------------------------------------------------------
!NXputattr writes an attribute of the open data set
!
!The following routines define the generic function NXputdata
!------------------------------------------------------------------------------
!NXputi1attr writes an integer*1 attribute of the open data set 
   FUNCTION NXputi1attr (file_id, name, value, value_length, value_type) &
                        RESULT (status)

      TYPE(NXhandle),     INTENT(in) :: file_id
      CHARACTER(len=*),   INTENT(in) :: name
      INTEGER(KIND=NXi1), INTENT(in) :: value
      INTEGER, OPTIONAL,  INTENT(in) :: value_length
      INTEGER, OPTIONAL,  INTENT(in) :: value_type
      INTEGER :: status, nxifputattr
      EXTERNAL nxifputattr

      status = nxifputattr(file_id, NXCstring(name), value, 1, NX_INT8)

   END FUNCTION NXputi1attr
!------------------------------------------------------------------------------
!NXputi2attr writes an integer*2 attribute of the open data set 
   FUNCTION NXputi2attr (file_id, name, value, value_length, value_type) &
                        RESULT (status)

      TYPE(NXhandle),     INTENT(in) :: file_id
      CHARACTER(len=*),   INTENT(in) :: name
      INTEGER(KIND=NXi2), INTENT(in) :: value
      INTEGER, OPTIONAL,  INTENT(in) :: value_length
      INTEGER, OPTIONAL,  INTENT(in) :: value_type
      INTEGER :: status, nxifputattr
      EXTERNAL nxifputattr

      status = nxifputattr(file_id, NXCstring(name), value, 1, NX_INT16)

   END FUNCTION NXputi2attr
!------------------------------------------------------------------------------
!NXputi4attr writes an integer*4 attribute of the open data set 
   FUNCTION NXputi4attr (file_id, name, value, value_length, value_type) &
                        RESULT (status)

      TYPE(NXhandle),     INTENT(in) :: file_id
      CHARACTER(len=*),   INTENT(in) :: name
      INTEGER(KIND=NXi4), INTENT(in) :: value
      INTEGER, OPTIONAL,  INTENT(in) :: value_length
      INTEGER, OPTIONAL,  INTENT(in) :: value_type
      INTEGER :: status, nxifputattr
      EXTERNAL nxifputattr

      status = nxifputattr(file_id, NXCstring(name), value, 1, NX_INT32)

   END FUNCTION NXputi4attr
!------------------------------------------------------------------------------
!NXputr4attr writes a real*4 attribute of the open data set 
   FUNCTION NXputr4attr (file_id, name, value, value_length, value_type) &
                        RESULT (status)

      TYPE(NXhandle),     INTENT(in) :: file_id
      CHARACTER(len=*),   INTENT(in) :: name
      REAL(KIND=NXr4),    INTENT(in) :: value
      INTEGER, OPTIONAL,  INTENT(in) :: value_length
      INTEGER, OPTIONAL,  INTENT(in) :: value_type
      INTEGER :: status, nxifputattr
      EXTERNAL nxifputattr

      status = nxifputattr(file_id, NXCstring(name), value, 1, NX_FLOAT32)

   END FUNCTION NXputr4attr
!------------------------------------------------------------------------------
!NXputr8attr writes a real*8 attribute of the open data set 
   FUNCTION NXputr8attr (file_id, name, value, value_length, value_type) &
                        RESULT (status)

      TYPE(NXhandle),     INTENT(in) :: file_id
      CHARACTER(len=*),   INTENT(in) :: name
      REAL(KIND=NXr8),    INTENT(in) :: value
      INTEGER, OPTIONAL,  INTENT(in) :: value_length
      INTEGER, OPTIONAL,  INTENT(in) :: value_type
      INTEGER :: status, nxifputattr
      EXTERNAL nxifputattr

      status = nxifputattr(file_id, NXCstring(name), value, 1, NX_FLOAT64)

   END FUNCTION NXputr8attr
!------------------------------------------------------------------------------
!NXputcharattr writes character attribute of the open data set 
   FUNCTION NXputcharattr (file_id, name, value, value_length, value_type) &
                        RESULT (status)

      TYPE(NXhandle),     INTENT(in) :: file_id
      CHARACTER(len=*),   INTENT(in) :: name
      CHARACTER(len=*),   INTENT(in) :: value
      INTEGER, OPTIONAL,  INTENT(in) :: value_length
      INTEGER, OPTIONAL,  INTENT(in) :: value_type
      INTEGER :: status, nxifputattr
      EXTERNAL nxifputattr

      status = nxifputattr(file_id, NXCstring(name), NXCstring(value), &
                        len_trim(value), NX_CHAR)

   END FUNCTION NXputcharattr
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!NXgetinfo gets the rank, dimensions and type of the open data set
   FUNCTION NXgetinfo (file_id, data_rank, data_dimensions, data_type) &
                        RESULT (status)

      TYPE(NXhandle), INTENT(in)  :: file_id
      INTEGER,        INTENT(out) :: data_rank, data_dimensions(:), data_type
      INTEGER :: status, nxigetinfo, i, j, dimensions(size(data_dimensions))
      EXTERNAL nxigetinfo

      status = nxigetinfo(file_id, data_rank, dimensions, data_type)
      IF (status == NX_OK) THEN
         data_dimensions = NXreverse (data_rank, dimensions)
      END IF

   END FUNCTION NXgetinfo
!------------------------------------------------------------------------------
!NXgetnextentry implements a directory search of the open group
   FUNCTION NXgetnextentry (file_id, name, class, data_type) RESULT (status)

      TYPE(NXhandle),   INTENT(in)  :: file_id
      CHARACTER(len=*), INTENT(out) :: name, class
      INTEGER,          INTENT(out) :: data_type
      INTEGER :: status, nxigetnextentry, i, j
      INTEGER(kind=NXi1) :: Cname(NX_MAXNAMELEN), Cclass(NX_MAXNAMELEN)
      EXTERNAL nxigetnextentry

      status = nxigetnextentry(file_id, Cname, Cclass, data_type)
      name = trim(NXFstring(Cname))
      class = trim(NXFstring(Cclass))

   END FUNCTION NXgetnextentry
!------------------------------------------------------------------------------
!NXgetnextattr implements a search of all the attributes of the open data set
   FUNCTION NXgetnextattr (file_id, attr_name, attr_length, attr_type) &
                        RESULT (status)

      TYPE(NXhandle),   INTENT(in)  :: file_id
      CHARACTER(len=*), INTENT(out) :: attr_name
      INTEGER,          INTENT(out) :: attr_length, attr_type
      INTEGER :: status, nxigetnextattr
      INTEGER(kind=NXi1) :: Cstring(NX_MAXNAMELEN)
      EXTERNAL nxigetnextattr

      status = nxigetnextattr(file_id, Cstring, attr_length, attr_type)
      attr_name = trim(NXFstring(Cstring))

   END FUNCTION NXgetnextattr
!------------------------------------------------------------------------------
!NXgetgroupID returns the identifier of the open group as an NXlink structure
   FUNCTION NXgetgroupID (file_id, group_id) RESULT (status)

      TYPE(NXhandle), INTENT(in)  :: file_id
      TYPE(NXlink),   INTENT(out) :: group_id
      TYPE(NXlink) :: current_id
      INTEGER :: status, nxigetgroupid
      EXTERNAL nxigetgroupid

      status = nxigetgroupid(file_id, current_id)
      group_id = current_id

   END FUNCTION NXgetgroupID
!------------------------------------------------------------------------------
!NXgetdataID returns the identifier of the open data set as an NXlink structure
   FUNCTION NXgetdataID (file_id, data_id) RESULT (status)

      TYPE(NXhandle), INTENT(in)  :: file_id
      TYPE(NXlink),   INTENT(out) :: data_id
      TYPE(NXlink) :: current_id
      INTEGER :: status, nxigetdataid
      EXTERNAL nxigetdataid

      status = nxigetdataid(file_id, current_id)
      data_id = current_id

   END FUNCTION NXgetdataID
!------------------------------------------------------------------------------
!NXsameID checks that two group or data ID's are the same
   FUNCTION NXsameID (file_id, first_id, second_id) RESULT (same)

      TYPE(NXhandle), INTENT(in) :: file_id
      TYPE(NXlink), INTENT(in)   :: first_id, second_id
      LOGICAL :: same
      INTEGER :: status, nxisameid
      EXTERNAL nxisameid

      status = nxisameid(file_id, first_id, second_id)
      IF (status == NX_OK) THEN
         same = .TRUE.
      ELSE
         same = .FALSE.
      ENDIF

   END FUNCTION NXsameID
!------------------------------------------------------------------------------
!NXmakelink links a data item (group or set) to another group 
   FUNCTION NXmakelink (file_id, link) RESULT (status)

      TYPE(NXhandle), INTENT(in) :: file_id
      TYPE(NXlink),   INTENT(in) :: link
      INTEGER :: status, nximakelink
      EXTERNAL nximakelink

      status = nximakelink(file_id, link)

   END FUNCTION NXmakelink
!------------------------------------------------------------------------------
!NXgetgroupinfo returns the number of entries, name and class of the open group
   FUNCTION NXgetgroupinfo (file_id, item_number, group_name, group_class) &
                        RESULT (status)

      TYPE(NXhandle),   INTENT(in)  :: file_id
      INTEGER,          INTENT(out) :: item_number
      CHARACTER(len=*), INTENT(out), OPTIONAL :: group_name, group_class
      TYPE(NXlink) :: group_id, new_id
      INTEGER :: status, nxigetgroupinfo
      INTEGER(kind=NXi1) :: Cname(NX_MAXNAMELEN), Cclass(NX_MAXNAMELEN)
      EXTERNAL nxigetgroupinfo

      status = nxigetgroupinfo (file_id, item_number, Cname, Cclass)
      IF (PRESENT(group_name)) group_name = trim(NXFstring(Cname))
      IF (PRESENT(group_class)) group_class = trim(NXFstring(Cclass))

   END FUNCTION NXgetgroupinfo
!------------------------------------------------------------------------------
!NXinitgroupdir initializes data searches using NXgetnextentry
   FUNCTION NXinitgroupdir (file_id) RESULT (status)

      TYPE(NXhandle), INTENT(inout) :: file_id
      INTEGER :: status, nxiinitgroupdir
      EXTERNAL nxiinitgroupdir

      status = nxiinitgroupdir (file_id)

  END FUNCTION NXinitgroupdir
!------------------------------------------------------------------------------
!NXgroupdir returns a list of items in the currently open group
   FUNCTION NXgroupdir (file_id, item_number, item_name, item_class) &
                        RESULT (status)

      TYPE(NXhandle),   INTENT(inout)  :: file_id
      INTEGER,          INTENT(out)    :: item_number
      CHARACTER(len=*)                 :: item_name(:), item_class(:)
      CHARACTER(len=len(item_name)) :: name
      CHARACTER(len=len(item_class)) :: class
      INTEGER :: status

      status = NXinitgroupdir (file_id)
      item_number = 0
      DO
         status = NXgetnextentry (file_id, name, class, NXtype)
         IF (status == NX_OK .AND. &
                        (class(1:2) == "NX" .OR. class(1:3) == "SDS")) THEN
            item_number = item_number + 1
            IF (item_number > size(item_name) .OR. &
                        item_number > size(item_class)) THEN
               CALL NXerror ("Number of items greater than array size")
               status = NX_ERROR
               RETURN
            END IF
            item_name(item_number) = trim(name)
            item_class(item_number) = trim(class)
         ELSE IF (status == NX_EOD) THEN
            EXIT
         ELSE IF (status == NX_ERROR) THEN
            RETURN
         END IF
      END DO
      status = NX_OK

   END FUNCTION NXgroupdir
!------------------------------------------------------------------------------
!NXgetattrinfo returns the number of attributes of the open data set
   FUNCTION NXgetattrinfo (file_id, attr_number) RESULT (status)

      TYPE(NXhandle),   INTENT(inout)  :: file_id
      INTEGER,          INTENT(out)    :: attr_number
      INTEGER :: status, nxigetattrinfo
      EXTERNAL nxigetattrinfo

      status = nxigetattrinfo (file_id, attr_number)

   END FUNCTION NXgetattrinfo
!------------------------------------------------------------------------------
!NXinitattrdir initializes attribute searches using NXgetnextattr
   FUNCTION NXinitattrdir (file_id) RESULT (status)

      TYPE(NXhandle), INTENT(inout) :: file_id
      INTEGER :: status, nxiinitattrdir
      EXTERNAL nxiinitattrdir

      status = nxiinitattrdir (file_id)

  END FUNCTION NXinitattrdir
!------------------------------------------------------------------------------
!NXattrdir returns a list of NeXus attributes of current data item
   FUNCTION NXattrdir (file_id, attr_number, attr_name) RESULT (status)

      TYPE(NXhandle),   INTENT(inout)  :: file_id
      INTEGER,          INTENT(out)    :: attr_number
      CHARACTER(len=*)    :: attr_name(:)
      CHARACTER(len=len(attr_name))    :: name
      INTEGER :: status

      status = NXinitattrdir (file_id)
      attr_number = 0
      DO
         status = NXgetnextattr (file_id, name, NXsize, NXtype)
         IF (status == NX_OK) THEN
            attr_number = attr_number + 1
            IF (attr_number > size(attr_name)) THEN
               CALL NXerror ("Number of attributes greater than array size")
               status = NX_ERROR
               RETURN
            ELSE
               attr_name(attr_number) = trim(name)
            END IF
         ELSE IF (status == NX_EOD) THEN
            EXIT
         ELSE IF (status == NX_ERROR) THEN
            RETURN
         END IF
      END DO
      status = NX_OK

   END FUNCTION NXattrdir
!------------------------------------------------------------------------------
!NXreverse reverses dimensions for transferring data from F90 to C
   FUNCTION NXreverse (rank, dimensions) RESULT (reversed_dimensions)

      INTEGER, INTENT(in) :: rank
      INTEGER, INTENT(in) :: dimensions(:)
      INTEGER :: reversed_dimensions(size(dimensions))
      INTEGER :: i

      DO i = 1,rank
         reversed_dimensions(i) = dimensions(rank-i+1)
      END DO

  END FUNCTION NXreverse
!------------------------------------------------------------------------------
!NXCstring converts a Fortran string into a C string
   FUNCTION NXCstring (string) RESULT (array)

      CHARACTER(len=*), INTENT(in) :: string
      INTEGER(kind=NXi1) :: array(255)
      INTEGER :: i

      DO i = 1,min(len_trim(string),(size(array)-1))
         array(i) = ichar(string(i:i))
      END DO
      array(len_trim(string)+1) = 0

  END FUNCTION NXCstring
!------------------------------------------------------------------------------
!NXFstring converts a C string into a Fortran string
   FUNCTION NXFstring (array) RESULT (string)

      INTEGER(kind=NXi1), INTENT(in) :: array(:)
      CHARACTER(len=255) :: string
      INTEGER :: i

      string = " "
      DO i = 1,size(array)
         IF (array(i) == 0) EXIT
         string(i:i) = char(array(i))
      END DO

  END FUNCTION NXFstring
!------------------------------------------------------------------------------
!NXdatatype converts a NeXus data type into a character string
   FUNCTION NXdatatype (int_type) RESULT (char_type)

      INTEGER, INTENT(in) :: int_type
      CHARACTER(len=10) :: char_type

      SELECT CASE (int_type)
         CASE(NX_CHAR); char_type = "NX_CHAR"
         CASE(NX_FLOAT32); char_type = "NX_FLOAT32"
         CASE(NX_FLOAT64); char_type = "NX_FLOAT64"
         CASE(NX_INT8); char_type = "NX_INT8"
         CASE(NX_INT16); char_type = "NX_INT16"
         CASE(NX_INT32); char_type = "NX_INT32"
         CASE(NX_UINT32); char_type = "NX_UINT32"
         CASE DEFAULT; char_type = "UNKNOWN"
      END SELECT

  END FUNCTION NXdatatype
!------------------------------------------------------------------------------
!NXerror prints out an error message to the default unit
   SUBROUTINE NXerror (message)

      CHARACTER(len=*), INTENT(in) :: message

      PRINT *, "NXerror : "//message

  END SUBROUTINE NXerror

END MODULE NXmodule
