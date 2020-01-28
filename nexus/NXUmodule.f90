!------------------------------------------------------------------------------
! NeXus - Neutron & X-ray Common Data Format
!  
! Fortran 90 Utilities
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

MODULE NXUmodule

   USE NXmodule
   PUBLIC
! *** NeXus utility functions ***
   PUBLIC :: NXUwriteglobals, NXUwritegroup, NXUwritedata, NXUreaddata
   PUBLIC :: NXUsetcompress
   PUBLIC :: NXUfindgroup, NXUfindclass, NXUfinddata, NXUfindattr
   PUBLIC :: NXUfindsignal, NXUfindaxis
   PUBLIC :: NXUfindlink, NXUresumelink
! *** NeXus utility internal functions
   PRIVATE :: NXUpreparedata, NXUconfirmdata, NXUsearchgroup
! *** NeXus utility global variables
   INTEGER, PRIVATE :: NXcompress_type = NX_COMP_NONE
   INTEGER, PRIVATE :: NXcompress_size = 1000
   INTEGER, PRIVATE :: group_level
   INTEGER, PRIVATE :: NXrank, NXdims(NX_MAXRANK), NXtype, NXsize
! *** NeXus generic interfaces ***
   INTERFACE NXUwritedata
       MODULE PROCEDURE NXUwritei4, NXUwriter4, NXUwriter8, NXUwritechar, &
                          NXUwritei4array, NXUwriter4array, &
                          NXUwriter8array, NXUwrite2Di4array, &
                          NXUwrite2Dr4array, NXUwrite2Dr8array, &
                          NXUwrite3Di4array, NXUwrite3Dr4array, &
                          NXUwrite3Dr8array
   END INTERFACE
   INTERFACE NXUreaddata
       MODULE PROCEDURE NXUreadi4, NXUreadr4, NXUreadr8, NXUreadchar, &
                          NXUreadi4array, NXUreadr4array, NXUreadr8array, &
                          NXUread2Di4array, NXUread2Dr4array, &
                          NXUread2Dr8array, NXUread3Di4array, &
                          NXUread3Dr4array, NXUread3Dr8array
   END INTERFACE

CONTAINS
!------------------------------------------------------------------------------
!NXUwriteglobals writes the global attributes to a file
   FUNCTION NXUwriteglobals (file_id, user, affiliation, address, phone, fax, &
                        email) RESULT (status)

      TYPE(NXhandle),   INTENT(in) :: file_id      
      CHARACTER(len=*), INTENT(in), OPTIONAL :: user, affiliation, address, &
                        phone, fax, email
      INTEGER :: status

      IF (PRESENT(user)) THEN
         status = NXputattr (file_id, "user", trim(user))
         IF (status /= NX_OK) RETURN
      END IF
      IF (PRESENT(affiliation)) THEN
         status = NXputattr (file_id, "affiliation", trim(affiliation))
         IF (status /= NX_OK) RETURN
      END IF
      IF (PRESENT(address)) THEN
         status = NXputattr (file_id, "address", trim(address))
         IF (status /= NX_OK) RETURN
      END IF
      IF (PRESENT(phone)) THEN
         status = NXputattr (file_id, "telephone_number", trim(phone))
         IF (status /= NX_OK) RETURN
      END IF
      IF (PRESENT(fax)) THEN
         status = NXputattr (file_id, "fax_number", trim(fax))
         IF (status /= NX_OK) RETURN
      END IF
      IF (PRESENT(email)) THEN
         status = NXputattr (file_id, "email", trim(email))
         IF (status /= NX_OK) RETURN
      END IF

   END FUNCTION NXUwriteglobals
!------------------------------------------------------------------------------
!NXUwritegroup creates and leaves open a group
   FUNCTION NXUwritegroup (file_id, group_name, group_class) RESULT (status)

      TYPE(NXhandle),   INTENT(in) :: file_id
      CHARACTER(len=*), INTENT(in) :: group_name, group_class
      INTEGER :: status

      status = NXmakegroup (file_id, group_name, group_class)
      IF (status == NX_OK) THEN
         status = NXopengroup (file_id, group_name, group_class)
      END IF

   END FUNCTION NXUwritegroup
!------------------------------------------------------------------------------
!NXUwritedata creates and writes a data set
!
!The following routines define the generic function NXUwritedata
!------------------------------------------------------------------------------
!NXUwritei4 writes a scalar integer*4 data item
   FUNCTION NXUwritei4 (file_id, data_name, data, units) RESULT (status)

      TYPE(NXhandle),     INTENT(inout) :: file_id
      CHARACTER(len=*),   INTENT(in)    :: data_name
      INTEGER(kind=NXi4), INTENT(in)    :: data
      CHARACTER(len=*),   INTENT(in), OPTIONAL :: units
      INTEGER :: status

      status = NXUpreparedata (file_id, data_name, NX_INT32, 1, (/1/))
      IF (status /= NX_OK) RETURN
      IF (PRESENT(units) .AND. NXUfindattr(file_id, "units") == NX_EOD) THEN
         status = NXputattr (file_id, "units", units)
         IF (status /= NX_OK) RETURN
      END IF
      status = NXputdata (file_id, (/ data /))

   END FUNCTION NXUwritei4
!------------------------------------------------------------------------------
!NXUwriter4 writes a scalar real*4 data item
   FUNCTION NXUwriter4 (file_id, data_name, data, units) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      REAL(kind=NXr4),  INTENT(in)    :: data
      CHARACTER(len=*), INTENT(in), OPTIONAL :: units
      INTEGER :: status

      status = NXUpreparedata (file_id, data_name, NX_FLOAT32, 1, (/1/))
      IF (status /= NX_OK) RETURN
      IF (PRESENT(units) .AND. NXUfindattr(file_id, "units") == NX_EOD) THEN
         status = NXputattr (file_id, "units", units)
         IF (status /= NX_OK) RETURN
      END IF
      status = NXputdata (file_id, (/ data /))

   END FUNCTION NXUwriter4
!------------------------------------------------------------------------------
!NXUwriter8 writes a scalar real*8 data item
   FUNCTION NXUwriter8 (file_id, data_name, data, units) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      REAL(kind=NXr8),  INTENT(in)    :: data
      CHARACTER(len=*), INTENT(in), OPTIONAL :: units
      INTEGER :: status

      status = NXUpreparedata (file_id, data_name, NX_FLOAT64, 1, (/1/))
      IF (status /= NX_OK) RETURN
      IF (PRESENT(units) .AND. NXUfindattr(file_id, "units") == NX_EOD) THEN
         status = NXputattr (file_id, "units", units)
         IF (status /= NX_OK) RETURN
      END IF
      status = NXputdata (file_id, (/ data /))

   END FUNCTION NXUwriter8
!------------------------------------------------------------------------------
!NXUwritechar writes a character data item
   FUNCTION NXUwritechar (file_id, data_name, data, units) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      CHARACTER(len=*), INTENT(in)    :: data
      CHARACTER(len=*), INTENT(in), OPTIONAL :: units
      INTEGER :: status

      status = NXUpreparedata (file_id, data_name, NX_CHAR, 1, &
                        (/len_trim(data)/))
      IF (status /= NX_OK) RETURN
      IF (PRESENT(units) .AND. NXUfindattr(file_id, "units") == NX_EOD) THEN
         status = NXputattr (file_id, "units", units, len_trim(units), NX_CHAR)
         IF (status /= NX_OK) RETURN
      END IF
      status = NXputdata (file_id, data)

   END FUNCTION NXUwritechar
!------------------------------------------------------------------------------
!NXUwritei4array writes 1D integer*4 array data
   FUNCTION NXUwritei4array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),     INTENT(inout) :: file_id
      CHARACTER(len=*),   INTENT(in)    :: data_name
      INTEGER(kind=NXi4), INTENT(in)    :: data(:)
      CHARACTER(len=*),   INTENT(in), OPTIONAL :: units
      INTEGER,            INTENT(in), OPTIONAL :: data_start(:), data_size(:)
      INTEGER :: status

      status = NXUpreparedata (file_id, data_name, NX_INT32, 1, (/size(data)/))
      IF (status /= NX_OK) RETURN
      IF (PRESENT(units) .AND. NXUfindattr(file_id, "units") == NX_EOD) THEN
         status = NXputattr (file_id, "units", units)
         IF (status /= NX_OK) RETURN
      END IF
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         status = NXputslab (file_id, data, data_start, data_size)
      ELSE
         status = NXputdata (file_id, data)
      END IF

   END FUNCTION NXUwritei4array
!------------------------------------------------------------------------------
!NXUwriter4array writes 1D real*4 array data
   FUNCTION NXUwriter4array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      REAL(kind=NXr4),  INTENT(in)    :: data(:)
      CHARACTER(len=*), INTENT(in), OPTIONAL :: units
      INTEGER,          INTENT(in), OPTIONAL :: data_start(:), data_size(:)
      INTEGER :: status

      status = NXUpreparedata (file_id, data_name, NX_FLOAT32, 1, &
                        (/size(data)/))
      IF (status /= NX_OK) RETURN
      IF (PRESENT(units) .AND. NXUfindattr(file_id, "units") == NX_EOD) THEN
         status = NXputattr (file_id, "units", units)
         IF (status /= NX_OK) RETURN
      END IF
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         status = NXputslab (file_id, data, data_start, data_size)
      ELSE
         status = NXputdata (file_id, data)
      END IF

   END FUNCTION NXUwriter4array
!------------------------------------------------------------------------------
!NXUwriter8array writes real*8 array data
   FUNCTION NXUwriter8array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      REAL(kind=NXr8),  INTENT(in)    :: data(:)
      CHARACTER(len=*), INTENT(in), OPTIONAL :: units
      INTEGER,          INTENT(in), OPTIONAL :: data_start(:), data_size(:)
      INTEGER :: status

      status = NXUpreparedata (file_id, data_name, NX_FLOAT64, 1, &
                        (/size(data)/))
      IF (status /= NX_OK) RETURN
      IF (PRESENT(units) .AND. NXUfindattr(file_id, "units") == NX_EOD) THEN
         status = NXputattr (file_id, "units", units)
         IF (status /= NX_OK) RETURN
      END IF
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         status = NXputslab (file_id, data, data_start, data_size)
      ELSE
         status = NXputdata (file_id, data)
      END IF

   END FUNCTION NXUwriter8array
!------------------------------------------------------------------------------
!NXUwrite2Di4array writes 2D integer*4 data
   FUNCTION NXUwrite2Di4array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),     INTENT(inout) :: file_id
      CHARACTER(len=*),   INTENT(in)    :: data_name
      INTEGER(kind=NXi4), INTENT(in)    :: data(:,:)
      CHARACTER(len=*),   INTENT(in), OPTIONAL :: units
      INTEGER,            INTENT(in), OPTIONAL :: data_start(:), data_size(:)
      INTEGER :: status
      INTEGER, ALLOCATABLE :: buffer(:)

      status = NXUpreparedata (file_id, data_name, NX_INT32, 2, shape(data))
      IF (status /= NX_OK) RETURN
      IF (PRESENT(units) .AND. NXUfindattr(file_id, "units") == NX_EOD) THEN
         status = NXputattr (file_id, "units", units)
         IF (status /= NX_OK) RETURN
      END IF
      ALLOCATE (buffer(size(data)))
      buffer = RESHAPE (data, (/ size(data) /))
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         status = NXputslab (file_id, buffer, data_start, data_size)
      ELSE
         status = NXputdata(file_id, buffer)
      END IF
      DEALLOCATE (buffer)

   END FUNCTION NXUwrite2Di4array
!------------------------------------------------------------------------------
!NXUwrite2Dr4array writes 2D real*4 data
   FUNCTION NXUwrite2Dr4array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      REAL(kind=NXr4),  INTENT(in)    :: data(:,:)
      CHARACTER(len=*), INTENT(in), OPTIONAL :: units
      INTEGER,          INTENT(in), OPTIONAL :: data_start(:), data_size(:)
      INTEGER :: status
      REAL(kind=NXr4), ALLOCATABLE :: buffer(:)

      status = NXUpreparedata (file_id, data_name, NX_FLOAT32, 2, shape(data))
      IF (status /= NX_OK) RETURN
      IF (PRESENT(units) .AND. NXUfindattr(file_id, "units") == NX_EOD) THEN
         status = NXputattr (file_id, "units", units)
         IF (status /= NX_OK) RETURN
      END IF
      ALLOCATE (buffer(size(data)))
      buffer = RESHAPE (data, (/ size(data) /))
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         status = NXputslab (file_id, buffer, data_start, data_size)
      ELSE
         status = NXputdata(file_id, buffer)
      END IF
      DEALLOCATE (buffer)

   END FUNCTION NXUwrite2Dr4array
!------------------------------------------------------------------------------
!NXUwrite2Dr8array writes 2D real*8 data
   FUNCTION NXUwrite2Dr8array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      REAL(kind=NXr8),  INTENT(in)    :: data(:,:)
      CHARACTER(len=*), INTENT(in), OPTIONAL :: units
      INTEGER,          INTENT(in), OPTIONAL :: data_start(:), data_size(:)
      INTEGER :: status
      REAL(kind=NXr8), ALLOCATABLE :: buffer(:)

      status = NXUpreparedata (file_id, data_name, NX_FLOAT64, 2, shape(data))
      IF (status /= NX_OK) RETURN
      IF (PRESENT(units) .AND. NXUfindattr(file_id, "units") == NX_EOD) THEN
         status = NXputattr (file_id, "units", units)
         IF (status /= NX_OK) RETURN
      END IF
      ALLOCATE (buffer(size(data)))
      buffer = RESHAPE (data, (/ size(data) /))
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         status = NXputslab (file_id, buffer, data_start, data_size)
      ELSE
         status = NXputdata(file_id, buffer)
      END IF
      DEALLOCATE (buffer)

   END FUNCTION NXUwrite2Dr8array
!------------------------------------------------------------------------------
!NXUwrite3Di4array writes 3D integer*4 data
   FUNCTION NXUwrite3Di4array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),     INTENT(inout) :: file_id
      CHARACTER(len=*),   INTENT(in)    :: data_name
      INTEGER(kind=NXi4), INTENT(in)    :: data(:,:,:)
      CHARACTER(len=*),   INTENT(in), OPTIONAL :: units
      INTEGER,            INTENT(in), OPTIONAL :: data_start(:), data_size(:)
      INTEGER :: status
      INTEGER, ALLOCATABLE :: buffer(:)

      status = NXUpreparedata (file_id, data_name, NX_INT32, 3, shape(data))
      IF (status /= NX_OK) RETURN
      IF (PRESENT(units) .AND. NXUfindattr(file_id, "units") == NX_EOD) THEN
         status = NXputattr (file_id, "units", units)
         IF (status /= NX_OK) RETURN
      END IF
      ALLOCATE (buffer(size(data)))
      buffer = RESHAPE (data, (/ size(data) /))
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         status = NXputslab (file_id, buffer, data_start, data_size)
      ELSE
         status = NXputdata(file_id, buffer)
      END IF
      DEALLOCATE (buffer)

   END FUNCTION NXUwrite3Di4array
!------------------------------------------------------------------------------
!NXUwrite3Dr4array writes 3D real*4 data
   FUNCTION NXUwrite3Dr4array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      REAL(kind=NXr4),  INTENT(in)    :: data(:,:,:)
      CHARACTER(len=*), INTENT(in), OPTIONAL :: units
      INTEGER,          INTENT(in), OPTIONAL :: data_start(:), data_size(:)
      INTEGER :: status
      REAL(kind=NXr4), ALLOCATABLE :: buffer(:)

      status = NXUpreparedata (file_id, data_name, NX_FLOAT32, 3, shape(data))
      IF (status /= NX_OK) RETURN
      IF (PRESENT(units) .AND. NXUfindattr(file_id, "units") == NX_EOD) THEN
         status = NXputattr (file_id, "units", units)
         IF (status /= NX_OK) RETURN
      END IF
      ALLOCATE (buffer(size(data)))
      buffer = RESHAPE (data, (/ size(data) /))
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         status = NXputslab (file_id, buffer, data_start, data_size)
      ELSE
         status = NXputdata(file_id, buffer)
      END IF
      DEALLOCATE (buffer)

   END FUNCTION NXUwrite3Dr4array
!------------------------------------------------------------------------------
!NXUwrite3Dr8array writes 3D real*8 data
   FUNCTION NXUwrite3Dr8array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      REAL(kind=NXr8),  INTENT(in)    :: data(:,:,:)
      CHARACTER(len=*), INTENT(in), OPTIONAL :: units
      INTEGER,          INTENT(in), OPTIONAL :: data_start(:), data_size(:)
      INTEGER :: status
      REAL(kind=NXr8), ALLOCATABLE :: buffer(:)

      status = NXUpreparedata (file_id, data_name, NX_FLOAT64, 3, shape(data))
      IF (status /= NX_OK) RETURN
      IF (PRESENT(units) .AND. NXUfindattr(file_id, "units") == NX_EOD) THEN
         status = NXputattr (file_id, "units", units)
         IF (status /= NX_OK) RETURN
      END IF
      ALLOCATE (buffer(size(data)))
      buffer = RESHAPE (data, (/ size(data) /))
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         status = NXputslab (file_id, buffer, data_start, data_size)
      ELSE
         status = NXputdata(file_id, buffer)
      END IF
      DEALLOCATE (buffer)

   END FUNCTION NXUwrite3Dr8array
!------------------------------------------------------------------------------
!NXUreaddata reads data
!
!The following routines define the generic function NXUreaddata
!------------------------------------------------------------------------------
!NXUreadi4 reads a scalar integer*4 data item
   FUNCTION NXUreadi4 (file_id, data_name, data, units) RESULT (status)

      TYPE(NXhandle),     INTENT(inout)  :: file_id
      CHARACTER(len=*),   INTENT(in)     :: data_name
      INTEGER(kind=NXi4), INTENT(out)    :: data
      CHARACTER(len=*),   INTENT(out), OPTIONAL :: units
      INTEGER :: status, dimensions(NX_MAXRANK)
      INTEGER(kind=NXi4) :: buffer(1)

      status = NXUconfirmdata (file_id, data_name, NX_INT32, 1, dimensions)
      IF (status /= NX_OK) RETURN
      IF (dimensions(1) /= 1) THEN
         status = NX_ERROR
         RETURN
      END IF
      status = NXgetdata(file_id, buffer)
      IF (status == NX_OK) THEN
         data = buffer(1)
         IF (PRESENT(units)) THEN
            status = NXgetattr (file_id, "units", units)
         END IF
      END IF

   END FUNCTION NXUreadi4
!------------------------------------------------------------------------------
!NXgetr4 reads a scalar real*4 data item
   FUNCTION NXUreadr4 (file_id, data_name, data, units) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      REAL(kind=NXr4),  INTENT(out)   :: data
      CHARACTER(len=*), INTENT(out), OPTIONAL :: units
      INTEGER :: status, dimensions(NX_MAXRANK)
      REAL(kind=NXr4) :: buffer(1)

      status = NXUconfirmdata (file_id, data_name, NX_FLOAT32, 1, dimensions)
      IF (status /= NX_OK) RETURN
      IF (dimensions(1) /= 1) THEN
         status = NX_ERROR
         RETURN
      END IF
      status = NXgetdata(file_id, buffer)
      IF (status == NX_OK) THEN
         data = buffer(1)
         IF (PRESENT(units)) THEN
            status = NXgetattr (file_id, "units", units)
         END IF
      END IF

   END FUNCTION NXUreadr4
!------------------------------------------------------------------------------
!NXgetr8 reads a scalar real*8 data item
   FUNCTION NXUreadr8 (file_id, data_name, data, units) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      REAL(kind=NXr8),  INTENT(out)   :: data
      CHARACTER(len=*), INTENT(out), OPTIONAL :: units
      INTEGER :: status, dimensions(NX_MAXRANK)
      REAL(kind=NXr8) :: buffer(1)

      status = NXUconfirmdata (file_id, data_name, NX_FLOAT64, 1, dimensions)
      IF (status /= NX_OK) RETURN
      IF (dimensions(1) /= 1) THEN
         status = NX_ERROR
         RETURN
      END IF
      status = NXgetdata(file_id, buffer)
      IF (status == NX_OK) THEN
         data = buffer(1)
         IF (PRESENT(units)) THEN
            status = NXgetattr (file_id, "units", units)
         END IF
      END IF

   END FUNCTION NXUreadr8
!------------------------------------------------------------------------------
!NXgetchar reads a character string
   FUNCTION NXUreadchar (file_id, data_name, data, units) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      CHARACTER(len=*), INTENT(out)   :: data
      CHARACTER(len=*), INTENT(out), OPTIONAL :: units
      INTEGER :: status, dimensions(NX_MAXRANK)

      status = NXUconfirmdata (file_id, data_name, NX_CHAR, 1, dimensions)
      IF (status /= NX_OK) RETURN
      IF (dimensions(1) > len(data)) THEN
         status = NX_ERROR
         RETURN
      END IF
      status = NXgetdata(file_id, data)
      IF (status == NX_OK .and. PRESENT(units)) THEN
            status = NXgetattr (file_id, "units", units)
      END IF

   END FUNCTION NXUreadchar
!------------------------------------------------------------------------------
!NXUreadi4array reads an integer*4 array
   FUNCTION NXUreadi4array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),     INTENT(inout) :: file_id
      CHARACTER(len=*),   INTENT(in)    :: data_name
      INTEGER(kind=NXi4), POINTER       :: data(:)
      CHARACTER(len=*),   INTENT(out), OPTIONAL :: units
      INTEGER,            INTENT(in),  OPTIONAL :: data_start(:), data_size(:)
      INTEGER :: status, dimensions(NX_MAXRANK) 

      status = NXUconfirmdata (file_id, data_name, NX_INT32, 1, dimensions)
      IF (status /= NX_OK) RETURN
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         ALLOCATE (data(data_size(1)))
         status = NXgetslab (file_id, data, data_start, data_size)
      ELSE
         ALLOCATE (data(dimensions(1)))
         status = NXgetdata (file_id, data)
      END IF
      IF (status == NX_OK .and. PRESENT(units)) THEN
         status = NXgetattr (file_id, "units", units)
      END IF

   END FUNCTION NXUreadi4array
!------------------------------------------------------------------------------
!NXUreadr4array reads a real*4 array
   FUNCTION NXUreadr4array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      REAL(kind=NXr4),  POINTER       :: data(:)
      CHARACTER(len=*), INTENT(out), OPTIONAL :: units
      INTEGER,          INTENT(in),  OPTIONAL :: data_start(:), data_size(:)
      INTEGER :: status, dimensions(NX_MAXRANK)

      status = NXUconfirmdata (file_id, data_name, NX_FLOAT32, 1, dimensions)
      IF (status /= NX_OK) RETURN
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         ALLOCATE (data(data_size(1)))
         status = NXgetslab (file_id, data, data_start, data_size)
      ELSE
         ALLOCATE (data(dimensions(1)))
         status = NXgetdata (file_id, data)
      END IF
      IF (status == NX_OK .and. PRESENT(units)) THEN
         status = NXgetattr (file_id, "units", units)
      END IF

   END FUNCTION NXUreadr4array
!------------------------------------------------------------------------------
!NXUreadr8array reads a real*8 array
   FUNCTION NXUreadr8array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      REAL(kind=NXr8),  POINTER       :: data(:)
      CHARACTER(len=*), INTENT(out), OPTIONAL :: units
      INTEGER,          INTENT(in),  OPTIONAL :: data_start(:), data_size(:)
      INTEGER :: status, dimensions(NX_MAXRANK)

      status = NXUconfirmdata (file_id, data_name, NX_FLOAT64, 1, dimensions)
      IF (status /= NX_OK) RETURN
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         ALLOCATE (data(data_size(1)))
         status = NXgetslab (file_id, data, data_start, data_size)
      ELSE
         ALLOCATE (data(dimensions(1)))
         status = NXgetdata (file_id, data)
      END IF
      IF (status == NX_OK .and. PRESENT(units)) THEN
         status = NXgetattr (file_id, "units", units)
      END IF

   END FUNCTION NXUreadr8array
!------------------------------------------------------------------------------
!NXUread2Di4array reads a 2D integer*4 array
   FUNCTION NXUread2Di4array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),     INTENT(inout) :: file_id
      CHARACTER(len=*),   INTENT(in)    :: data_name
      INTEGER(kind=NXi4), POINTER       :: data(:,:)
      CHARACTER(len=*),   INTENT(out), OPTIONAL :: units
      INTEGER,            INTENT(in),  OPTIONAL :: data_start(:), data_size(:)
      INTEGER :: status, dimensions(NX_MAXRANK), data_shape(2)
      INTEGER, ALLOCATABLE :: buffer(:)

      status = NXUconfirmdata (file_id, data_name, NX_INT32, 2, dimensions)
      IF (status /= NX_OK) RETURN
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         ALLOCATE (buffer(PRODUCT(data_size(1:2))))
         status = NXgetslab (file_id, buffer, data_start, data_size)
         IF (status == NX_OK) THEN
            ALLOCATE (data(data_size(1),data_size(2)))
            data_shape = data_size(1:2)
            data = RESHAPE (buffer, data_shape)
         END IF
      ELSE
         ALLOCATE (buffer(PRODUCT(dimensions(1:2))))
         status = NXgetdata(file_id, buffer)
         IF (status == NX_OK) THEN
            ALLOCATE (data(dimensions(1),dimensions(2)))
            data = RESHAPE (buffer, dimensions(1:2))
         END IF
      END IF
      IF (status == NX_OK .and. PRESENT(units)) THEN
         status = NXgetattr (file_id, "units", units)
      END IF
      DEALLOCATE (buffer)

   END FUNCTION NXUread2Di4array
!------------------------------------------------------------------------------
!NXUread2Dr4array reads a 2D real*4 array
   FUNCTION NXUread2Dr4array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      REAL(kind=NXr4),  POINTER       :: data(:,:)
      CHARACTER(len=*), INTENT(out), OPTIONAL :: units
      INTEGER,          INTENT(in),  OPTIONAL ::data_start(:), data_size(:)
      INTEGER :: status, dimensions(NX_MAXRANK), data_shape(2)
      REAL, ALLOCATABLE :: buffer(:)

      status = NXUconfirmdata (file_id, data_name, NX_FLOAT32, 2, dimensions)
      IF (status /= NX_OK) RETURN
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         ALLOCATE (buffer(PRODUCT(data_size(1:2))))
         status = NXgetslab (file_id, buffer, data_start, data_size)
         IF (status == NX_OK) THEN
            ALLOCATE (data(data_size(1),data_size(2)))
            data_shape = data_size(1:2)
            data = RESHAPE (buffer, data_shape)
         END IF
      ELSE
         ALLOCATE (buffer(PRODUCT(dimensions(1:2))))
         status = NXgetdata(file_id, buffer)
         IF (status == NX_OK) THEN
            ALLOCATE (data(dimensions(1),dimensions(2)))
            data = RESHAPE (buffer, dimensions(1:2))
         END IF
      END IF
      IF (status == NX_OK .and. PRESENT(units)) THEN
         status = NXgetattr (file_id, "units", units)
      END IF
      DEALLOCATE (buffer)

   END FUNCTION NXUread2Dr4array
!------------------------------------------------------------------------------
!NXUread2Dr8array reads a 2D real*8 precision array
   FUNCTION NXUread2Dr8array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      REAL(kind=NXr8),  POINTER       :: data(:,:)
      CHARACTER(len=*), INTENT(out), OPTIONAL :: units
      INTEGER,          INTENT(in),  OPTIONAL ::data_start(:), data_size(:)
      INTEGER :: status, dimensions(NX_MAXRANK), data_shape(2)
      REAL, ALLOCATABLE :: buffer(:)

      status = NXUconfirmdata (file_id, data_name, NX_FLOAT64, 2, dimensions)
      IF (status /= NX_OK) RETURN
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         ALLOCATE (buffer(PRODUCT(data_size(1:2))))
         status = NXgetslab (file_id, buffer, data_start, data_size)
         IF (status == NX_OK) THEN
            ALLOCATE (data(data_size(1),data_size(2)))
            data_shape = data_size(1:2)
            data = RESHAPE (buffer, data_shape)
         END IF
      ELSE
         ALLOCATE (buffer(PRODUCT(dimensions(1:2))))
         status = NXgetdata(file_id, buffer)
         IF (status == NX_OK) THEN
            ALLOCATE (data(dimensions(1),dimensions(2)))
            data = RESHAPE (buffer, dimensions(1:2))
         END IF
      END IF
      IF (status == NX_OK .and. PRESENT(units)) THEN
         status = NXgetattr (file_id, "units", units)
      END IF
      DEALLOCATE (buffer)

   END FUNCTION NXUread2Dr8array
!------------------------------------------------------------------------------
!NXUread3Di4array reads a 3D integer*4 array
   FUNCTION NXUread3Di4array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),     INTENT(inout) :: file_id
      CHARACTER(len=*),   INTENT(in)    :: data_name
      INTEGER(kind=NXi4), POINTER       :: data(:,:,:)
      CHARACTER(len=*),   INTENT(out), OPTIONAL :: units
      INTEGER,            INTENT(in),  OPTIONAL :: data_start(:), data_size(:)
      INTEGER :: status, dimensions(NX_MAXRANK), data_shape(3)
      INTEGER, ALLOCATABLE :: buffer(:)

      status = NXUconfirmdata (file_id, data_name, NX_INT32, 3, dimensions)
      IF (status /= NX_OK) RETURN
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         ALLOCATE (buffer(PRODUCT(data_size(1:3))))
         status = NXgetslab (file_id, buffer, data_start, data_size)
         IF (status == NX_OK) THEN
            ALLOCATE (data(data_size(1),data_size(2),data_size(3)))
            data_shape = data_size(1:3)
            data = RESHAPE (buffer, data_shape)
         END IF
      ELSE
         ALLOCATE (buffer(PRODUCT(dimensions(1:3))))
         status = NXgetdata(file_id, buffer)
         IF (status == NX_OK) THEN
            ALLOCATE (data(dimensions(1),dimensions(2),dimensions(3)))
            data = RESHAPE (buffer, dimensions(1:3))
         END IF
      END IF
      IF (status == NX_OK .and. PRESENT(units)) THEN
         status = NXgetattr (file_id, "units", units)
      END IF
      DEALLOCATE (buffer)

   END FUNCTION NXUread3Di4array
!------------------------------------------------------------------------------
!NXUread3Dr4array reads a 3D real*4 array
   FUNCTION NXUread3Dr4array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      REAL(kind=NXr4),  POINTER       :: data(:,:,:)
      CHARACTER(len=*), INTENT(out), OPTIONAL :: units
      INTEGER,          INTENT(in),  OPTIONAL :: data_start(:), data_size(:)
      INTEGER :: status, dimensions(NX_MAXRANK), data_shape(3)
      REAL, ALLOCATABLE :: buffer(:)

      status = NXUconfirmdata (file_id, data_name, NX_FLOAT32, 3, dimensions)
      IF (status /= NX_OK) RETURN
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         ALLOCATE (buffer(PRODUCT(data_size(1:3))))
         status = NXgetslab (file_id, buffer, data_start, data_size)
         IF (status == NX_OK) THEN
            ALLOCATE (data(data_size(1),data_size(2),data_size(3)))
            data_shape = data_size(1:3)
            data = RESHAPE (buffer, data_shape)
         END IF
      ELSE
         ALLOCATE (buffer(PRODUCT(dimensions(1:3))))
         status = NXgetdata(file_id, buffer)
         IF (status == NX_OK) THEN
            ALLOCATE (data(dimensions(1),dimensions(2),dimensions(3)))
            data = RESHAPE (buffer, dimensions(1:3))
         END IF
      END IF
      IF (status == NX_OK .and. PRESENT(units)) THEN
         status = NXgetattr (file_id, "units", units)
      END IF
      DEALLOCATE (buffer)

   END FUNCTION NXUread3Dr4array
!------------------------------------------------------------------------------
!NXUread3Dr8array reads a 3D real*8 array
   FUNCTION NXUread3Dr8array (file_id, data_name, data, units, data_start, &
                        data_size) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      REAL(kind=NXr8),  POINTER       :: data(:,:,:)
      CHARACTER(len=*), INTENT(out), OPTIONAL :: units
      INTEGER,          INTENT(in),  OPTIONAL :: data_start(:), data_size(:)
      INTEGER :: status, dimensions(NX_MAXRANK), data_shape(3)
      REAL, ALLOCATABLE :: buffer(:)

      status = NXUconfirmdata (file_id, data_name, NX_FLOAT64, 3, dimensions)
      IF (status /= NX_OK) RETURN
      IF (PRESENT(data_start) .AND. PRESENT(data_size)) THEN
         ALLOCATE (buffer(PRODUCT(data_size(1:3))))
         status = NXgetslab (file_id, buffer, data_start, data_size)
         IF (status == NX_OK) THEN
            ALLOCATE (data(data_size(1),data_size(2),data_size(3)))
            data_shape = data_size(1:3)
            data = RESHAPE (buffer, data_shape)
         END IF
      ELSE
         ALLOCATE (buffer(PRODUCT(dimensions(1:3))))
         status = NXgetdata(file_id, buffer)
         IF (status == NX_OK) THEN
            ALLOCATE (data(dimensions(1),dimensions(2),dimensions(3)))
            data = RESHAPE (buffer, dimensions(1:3))
         END IF
      END IF
      IF (status == NX_OK .and. PRESENT(units)) THEN
         status = NXgetattr (file_id, "units", units)
      END IF
      DEALLOCATE (buffer)

   END FUNCTION NXUread3Dr8array
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!NXUsetcompress sets the default compression type and minimum size
   FUNCTION NXUsetcompress (file_id, compress_type, compress_size) &
                        RESULT (status)

      TYPE(NXhandle), INTENT(inout) :: file_id
      INTEGER,        INTENT(in)    :: compress_type
      INTEGER,        INTENT(in), OPTIONAL :: compress_size
      INTEGER :: status

      IF (compress_type == NX_COMP_LZW .OR. compress_type == NX_COMP_HUF .OR. &
          compress_type == NX_COMP_RLE .OR. compress_type == NX_COMP_NONE) THEN
         NXcompress_type = compress_type
         IF (PRESENT(compress_size)) NXcompress_size = compress_size
         status = NX_OK
      ELSE
         call NXerror ("Invalid compression option")
         status = NX_ERROR
      END IF

   END FUNCTION NXUsetcompress
!------------------------------------------------------------------------------
!NXUfindgroup finds if a NeXus group of the specified name exists
   FUNCTION NXUfindgroup (file_id, group_name, group_class) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: group_name
      CHARACTER(len=*), INTENT(out), OPTIONAL :: group_class
      CHARACTER(len=NX_MAXNAMELEN), ALLOCATABLE :: name(:), class(:)
      INTEGER :: status, n, i

      status = NXgetgroupinfo (file_id, n)
      IF (status /= NX_OK) RETURN
      ALLOCATE (name(n), class(n), STAT=status)
      IF (status /= 0) THEN
         call NXerror ("Unable to allocate directory arrays")
         status = NX_ERROR
         RETURN
      END IF
      status = NXgroupdir (file_id, n, name, class)
      IF (status == NX_OK) THEN
         status = NX_EOD
         DO i = 1,n
            IF (trim(name(i)) == trim(group_name)) THEN
               group_class = trim(class(i))
               IF (class(i)(1:2) == "NX") THEN
                  status = NX_OK
               ELSE
                  CALL NXerror (trim(name(i))//" is not a group")
                  status = NX_ERROR
               END IF
               EXIT
            END IF
         END DO
      END IF
      DEALLOCATE (name, class)

   END FUNCTION NXUfindgroup
!------------------------------------------------------------------------------
!NXUfindclass finds if a NeXus group of the specified class exists
   FUNCTION NXUfindclass (file_id, group_class, group_name, find_index) &
                        RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: group_class
      CHARACTER(len=*), INTENT(out)   :: group_name
      INTEGER,          INTENT(in), OPTIONAL :: find_index
      CHARACTER(len=NX_MAXNAMELEN), ALLOCATABLE :: name(:), class(:)
      INTEGER :: status, n, i, j

      status = NXgetgroupinfo (file_id, n)
      IF (status /= NX_OK) RETURN
      ALLOCATE (name(n), class(n), STAT=status)
      IF (status /= 0) THEN
         CALL NXerror ("Unable to allocate directory arrays")
         status = NX_ERROR
         RETURN
      END IF
      status = NXgroupdir (file_id, n, name, class)
      IF (status == NX_OK) THEN
         j = 0
         status = NX_EOD
         DO i = 1,n
            IF (trim(class(i)) == trim(group_class)) THEN
               IF (PRESENT(find_index)) THEN
                  j = j + 1
                  IF (j < find_index) CYCLE
               END IF
               group_name = trim(name(i))
               status = NX_OK
               EXIT
            END IF
         END DO
      END IF
      DEALLOCATE (name, class)

   END FUNCTION NXUfindclass
!------------------------------------------------------------------------------
!NXUfinddata finds if a NeXus data item is in the current group
   FUNCTION NXUfinddata (file_id, data_name) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      CHARACTER(len=NX_MAXNAMELEN), ALLOCATABLE :: name(:), class(:)
      INTEGER :: status, n, i

      status = NXgetgroupinfo (file_id, n)
      IF (status /= NX_OK) RETURN
      ALLOCATE (name(n), class(n), STAT=status)
      IF (status /= 0) THEN
         call NXerror ("Unable to allocate directory arrays")
         status = NX_ERROR
         RETURN
      END IF
      status = NXgroupdir (file_id, n, name, class)
      IF (status == NX_OK) THEN
         status = NX_EOD
         DO i = 1,n
            IF (trim(name(i)) == trim(data_name)) THEN
               IF (class(i)(1:3) == "SDS") THEN
                  status = NX_OK
               ELSE
                  CALL NXerror (trim(name(i))//" is not a data item")
                  status = NX_ERROR
               END IF
               EXIT
            END IF
         END DO
      END IF
      DEALLOCATE (name, class)

   END FUNCTION NXUfinddata
!------------------------------------------------------------------------------
!NXUfindattr finds if a NeXus attribute exists
   FUNCTION NXUfindattr (file_id, attr_name) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: attr_name
      CHARACTER(len=NX_MAXNAMELEN), ALLOCATABLE :: name(:)
      INTEGER :: status, n, i

      status = NXgetattrinfo (file_id, n)
      IF (status /= NX_OK) RETURN
      ALLOCATE (name(n), STAT=status)
      IF (status /= 0) THEN
         call NXerror ("Unable to allocate directory arrays")
         status = NX_ERROR
         RETURN
      END IF
      status = NXattrdir (file_id, n, name)
      IF (status == NX_OK) THEN
         status = NX_EOD
         DO i = 1,n
            IF (trim(name(i)) == trim(attr_name)) status = NX_OK
         END DO
      END IF
      DEALLOCATE (name)

   END FUNCTION NXUfindattr
!------------------------------------------------------------------------------
!NXUfindsignal finds the NeXus data item containing the required signal
   FUNCTION NXUfindsignal (file_id, signal, data_name, data_rank, data_type, &
                        data_dimensions) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      INTEGER,          INTENT(in)    :: signal
      CHARACTER(len=*)                :: data_name
      INTEGER,          INTENT(out)   :: data_rank, data_type, data_dimensions(:)
      CHARACTER(len=len(data_name)) :: name
      CHARACTER(len=NX_MAXNAMELEN) :: class, attr_name
      INTEGER :: status, value

      status = NXinitgroupdir (file_id)
      IF (status /= NX_OK) RETURN
      DO
         status = NXgetnextentry (file_id, name, class, NXtype)
         IF (status == NX_OK .AND. class == "SDS") THEN
            status = NXopendata (file_id, name)
            IF (status /= NX_OK) RETURN
            status = NXUfindattr (file_id, "signal")
            IF (status == NX_OK) THEN
               status = NXgetattr (file_id, "signal", value)
               IF (status /= NX_OK) RETURN
               IF (value == signal) THEN
                  status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
                  IF (status == NX_OK) THEN
                     data_name = name
                     data_rank = NXrank
                     data_type = NXtype
                     data_dimensions = NXdims
                     RETURN
                  END IF
               END IF
            ELSE IF (status == NX_EOD) THEN
               CYCLE
            ELSE IF (status == NX_ERROR) THEN         
               RETURN
            END IF
         ELSE IF (status == NX_EOD) THEN
            CALL NXerror ("No data with the attribute ""signal"" found")
            status = NX_ERROR
            EXIT
         ELSE IF (status == NX_ERROR) THEN
            RETURN
         END IF
      END DO

   END FUNCTION NXUfindsignal
!------------------------------------------------------------------------------
!NXUfindaxis finds the NeXus data item containing the required axis
   FUNCTION NXUfindaxis (file_id, axis, primary, data_name, data_type, &
                        data_dimensions) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      INTEGER,          INTENT(in)    :: axis, primary
      CHARACTER(len=*)                :: data_name
      INTEGER,          INTENT(out)   :: data_type, data_dimensions(NX_MAXRANK)
      CHARACTER(len=len(data_name)) :: name
      CHARACTER(len=NX_MAXNAMELEN) :: class, attr_name
      CHARACTER(len=255) :: axis_list
      INTEGER :: status, signal=1, value, data_rank, C_axis, i, j, k

      !First find data with "signal" attribute to check for "axes" attribute
      status = NXUfindsignal (file_id, signal, data_name, data_rank, &
                        data_type, data_dimensions)
      IF (status /= NX_OK) RETURN
      !The axis number cannot be greater than the data rank
      IF (axis > data_rank) THEN
         CALL NXerror ("Axis number greater than the data rank")
         status = NX_ERROR
         RETURN
      END IF
      !Check for "axes" attribute
      status = NXopendata (file_id, data_name)
      IF (status /= NX_OK) RETURN
      status = NXUfindattr (file_id, "axes")
      IF (status == NX_ERROR) THEN
         RETURN
      ELSE IF (status == NX_OK) THEN !"axes" attribute found
         status = NXgetattr (file_id, "axes", axis_list)
         !Strip off brackets around axis list
         IF (index(axis_list,"[") > 0) THEN 
            axis_list = axis_list(index(axis_list,"[")+1:len(axis_list))
         END IF
         IF (index(axis_list,"]") > 0) THEN
            axis_list = axis_list(1:index(axis_list,"]")-1)
         END IF
         !"axes" lists the axes in C (row-major) order so the axis numbers are reversed
         C_axis = data_rank - axis + 1 
         !Find axis label by looking for the delimiting commas
         j = 1
         DO i = 1,C_axis
            k = scan(axis_list(j:),",:") - 1
            IF (k < 0) k = len(trim(axis_list)) - j + 1
            IF (k < 0) THEN !We've run out of delimiters
               CALL NXerror ("Data attribute ""axes"" is not correctly defined")
               status = NX_ERROR
               RETURN
            END IF
            name = adjustl(axis_list(j:j+k-1))
            j = j + k + 1
         END DO
         !Open data to retrieve information about the dimension scale
         status = NXopendata (file_id, name)
         IF (status /= NX_OK) RETURN
         status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
         IF (status == NX_OK) THEN
            data_name = name
            data_type = NXtype
            data_dimensions = NXdims(1)
            RETURN
         ELSE
            RETURN
         END IF
      END IF   
      !Otherwise, check for "axis" attribute in each NXdata item
      status = NXinitgroupdir (file_id)
      IF (status /= NX_OK) RETURN
      DO
         status = NXgetnextentry (file_id, name, class, NXtype)
         IF (status == NX_OK .AND. class == "SDS") THEN
            status = NXopendata (file_id, name)
            IF (status /= NX_OK) RETURN
            status = NXUfindattr (file_id, "axis")
            IF (status == NX_OK) THEN
               status = NXgetattr (file_id, "axis", value)
               IF (status /= NX_OK) RETURN
               IF (value == axis) THEN
                  status = NXUfindattr (file_id, "primary")
                  IF (status == NX_OK) THEN
                     status = NXgetattr (file_id, "primary", value)
                  ELSE IF (status == NX_EOD) THEN
                     value = 1
                  ELSE
                     RETURN
                  END IF
                  IF (value == primary) THEN
                     status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
                     IF (status == NX_OK) THEN
                        data_name = name
                        data_type = NXtype
                        data_dimensions = NXdims(1)
                        RETURN
                     ELSE
                        RETURN
                     END IF
                  END IF
                END IF
             END IF
         ELSE IF (status == NX_EOD) THEN
            CALL NXerror ("Requested axis not found")
            status = NX_ERROR
            EXIT
         ELSE IF (status == NX_ERROR) THEN
            RETURN
         END IF
      END DO

   END FUNCTION NXUfindaxis
!------------------------------------------------------------------------------
!NXUfindlink finds another link to a NeXus data item and opens the group
   FUNCTION NXUfindlink (file_id, group_id, group_class) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      TYPE(NXlink),     INTENT(out)   :: group_id
      CHARACTER(len=*), INTENT(in), OPTIONAL :: group_class
      TYPE(NXlink) :: data_id
      INTEGER :: status

      !Get current group and data IDs
      status = NXgetgroupID (file_id, group_id)
      IF (status /= NX_OK) RETURN
      status = NXgetdataID (file_id, data_id)
      IF (status /= NX_OK) RETURN
      !Start the search in the group one level up
      status = NXclosegroup (file_id)
      IF (status /= NX_OK) RETURN
      !Start recursive searches for this data ID within this group 
      group_level = 0
      status = NXUsearchgroup (file_id, group_id, data_id, group_class)

   END FUNCTION NXUfindlink
!------------------------------------------------------------------------------
!NXUresumelink reopens the original group from which NXUfindlink was called
   FUNCTION NXUresumelink (file_id, group_id) RESULT (status)

      TYPE(NXhandle), INTENT(inout) :: file_id
      TYPE(NXlink),   INTENT(in)    :: group_id
      TYPE(NXlink) :: new_id
      CHARACTER(len=NX_MAXNAMELEN), ALLOCATABLE :: name(:), class(:)
      INTEGER :: status, n, i

      !Return to group level from which the link search was started
      DO i = 1, group_level
         status = NXclosegroup (file_id)
         IF (status /= NX_OK) RETURN
      END DO
      !Obtain list of groups at this level
      status = NXgetgroupinfo (file_id, n)
      IF (status /= NX_OK) RETURN
      ALLOCATE (name(n), class(n), STAT=status)
      IF (status /= 0) THEN
         CALL NXerror ("Unable to allocate space for group info")
         status = NX_ERROR
         RETURN
      END IF
      status = NXgroupdir (file_id, n, name, class)
      IF (status == NX_OK) THEN
         DO i = 1,n
            IF (class(i)(1:2) == "NX") THEN
               status = NXopengroup (file_id, name(i), class(i))
               IF (status /= NX_OK) EXIT
               status = NXgetgroupID (file_id, new_id)
               IF (status /= NX_OK) EXIT
               IF (NXsameID (file_id, new_id, group_id)) EXIT !Original group found
               status = NXclosegroup (file_id)
               IF (status /= NX_OK) EXIT
            END IF       
            status = NX_EOD
         END DO
      END IF
      !None of the groups was the correct one
      DEALLOCATE (name, class)

   END FUNCTION NXUresumelink
!------------------------------------------------------------------------------
!NXUsearchgroup searches a group for the required data
   RECURSIVE FUNCTION NXUsearchgroup (file_id, group_id, data_id, &
                        group_class) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      TYPE(NXlink),     INTENT(in)    :: group_id, data_id
      CHARACTER(len=*), INTENT(in), OPTIONAL :: group_class
      TYPE(NXlink) :: new_id
      CHARACTER(len=NX_MAXNAMELEN), ALLOCATABLE :: name(:), class(:)
      CHARACTER(len=NX_MAXNAMELEN) :: current_group, current_class
      INTEGER :: status, n, i

      !Obtain list of groups contained within this group
      status = NXgetgroupinfo (file_id, n, current_group, current_class)
      IF (status /= NX_OK) RETURN
      ALLOCATE (name(n), class(n), STAT=status)
      IF (status /= 0) THEN
         CALL NXerror ("Unable to allocate space for group info")
         status = NX_ERROR
         RETURN
      END IF      
      status = NXgroupdir (file_id, n, name, class)
      IF (status == NX_OK) THEN
         DO i = 1,n
            IF (class(i)(1:3) == "SDS") THEN
               IF (PRESENT(group_class) .AND. &
                        trim(group_class) /= trim(current_class)) THEN
                  status = NX_EOD
                  CYCLE
               END IF
               status = NXopendata (file_id, name(i))
               IF (status /= NX_OK) EXIT       
               status = NXgetdataID (file_id, new_id)
               IF (status /= NX_OK) EXIT
               IF (NXsameID (file_id, new_id, data_id)) THEN !Linked item found
                  status = NX_OK
                  EXIT
               END IF
            ELSE IF (class(i)(1:2) == "NX") THEN
               status = NXopengroup (file_id, name(i), class(i))
               IF (status /= NX_OK) EXIT
               status = NXgetgroupID (file_id, new_id)
               IF (status /= NX_OK) EXIT
               !Skip this group if it's where we started
               IF (NXsameID (file_id, new_id, group_id)) THEN
                  status = NXclosegroup (file_id)
                  IF (status /= NX_OK) EXIT
                  CYCLE
               END IF
               group_level = group_level + 1
               status = NXUsearchgroup(file_id, group_id, data_id, group_class)
               IF (status == NX_OK) EXIT !The item must have been found
               status = NXclosegroup (file_id)
               group_level = group_level - 1
               IF (status /= NX_OK) EXIT
            END IF
            status = NX_EOD
         END DO
      END IF
      !Return an error status because nothing has been found in this group
      DEALLOCATE (name, class)

   END FUNCTION NXUsearchgroup
!------------------------------------------------------------------------------
!NXUpreparedata creates and opens a data set
   FUNCTION NXUpreparedata (file_id, data_name, data_type, data_rank, &
                        data_dimensions) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      INTEGER,          INTENT(in)    :: data_type, data_rank
      INTEGER,          INTENT(in)    :: data_dimensions(:)
      INTEGER :: status, i

      status = NXUfinddata (file_id, data_name)
      IF (status == NX_EOD) THEN       !Data item needs to be created
         IF (NXcompress_type /= NX_COMP_NONE .AND. &
             PRODUCT(data_dimensions(1:data_rank)) > NXcompress_size) THEN
            status = NXmakedata (file_id, data_name, data_type, data_rank, &
                        data_dimensions, NXcompress_type)
         ELSE
            status = NXmakedata (file_id, data_name, data_type, data_rank, &
                        data_dimensions)
         END IF
         IF (status == NX_OK) status = NXopendata (file_id, data_name)
      ELSE if (status == NX_OK) THEN   !Data item already exists
         status = NXopendata (file_id, data_name)
         IF (status /= NX_OK) RETURN
         status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
         IF (NXtype /= data_type) THEN
            CALL NXerror ("Type of existing data item does not match new data")
            status = NX_ERROR
         ELSE IF (NXrank /= data_rank) THEN
            CALL NXerror ("Rank of existing data item does not match new data")
            status = NX_ERROR
         ELSE
            DO i = 1,NXrank
               IF (data_dimensions(i) > NXdims(i)) THEN
                  call NXerror ("Size of new data too large for existing item")
                  status = NX_ERROR
                  EXIT
               END IF
            END DO
         END IF
      END IF
      
   END FUNCTION NXUpreparedata
!------------------------------------------------------------------------------
!NXUconfirmdata checks that a dataset has the expected type, rank & dimensions
   FUNCTION NXUconfirmdata (file_id, data_name, data_type, data_rank, &
                        data_dimensions) RESULT (status)

      TYPE(NXhandle),   INTENT(inout) :: file_id
      CHARACTER(len=*), INTENT(in)    :: data_name
      INTEGER,          INTENT(in)    :: data_type, data_rank
      INTEGER,          INTENT(out)   :: data_dimensions(:)
      INTEGER :: status

      status = NXopendata (file_id, data_name)
      IF (status /= NX_OK) RETURN
      status = NXgetinfo (file_id, NXrank, NXdims, NXtype)
      IF (status /= NX_OK) RETURN
      IF (NXrank == data_rank) THEN
         !Check that the types match, or that they are both integer or real
         IF (NXtype /= data_type .AND. (NXtype/10) /= (data_type/10)) THEN
            CALL NXerror ("Type of data does not match supplied array")
         ELSE
            data_dimensions(1:NXrank) = NXdims(1:NXrank)
            status = NX_OK
            RETURN
         END IF
      ELSE
         CALL NXerror ("Rank of data does not match supplied array")
      END IF
      status = NXclosedata(file_id)
      status = NX_ERROR

   END FUNCTION NXUconfirmdata
         
END MODULE NXUmodule
