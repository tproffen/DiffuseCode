set (H5SOURCES attr.f90 attr_read.f90 attr_write.f90 datatype.f90 deflate.f90
               interface.f90 read_ascii.f90 reader.f90 reader_lt.f90 read.f90
               read_scalar.f90 utils.f90 write.f90 writer.f90 writer_lt.f90
               write_scalar.f90 )

add_library ( h5fortran_lib ${H5SOURCES} )
