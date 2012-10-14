!+
!     include file for macro and file related variables
!-
      INTEGER            MAC_MAX_LEVEL
      PARAMETER         (MAC_MAX_LEVEL = 10)

      INTEGER            MAC_MAX_PARA
      PARAMETER         (MAC_MAX_PARA = 20)

      INTEGER            MAC_MAX_IO
      PARAMETER         (MAC_MAX_IO = 10)
!
      INTEGER            MAC_MAX_FORM
      PARAMETER         (MAC_MAX_FORM = 50)
!
      CHARACTER*1024     mac_name(MAC_MAX_LEVEL)
      CHARACTER*80       mac_para(0:MAC_MAX_PARA,MAC_MAX_LEVEL)
      INTEGER            mac_leng(0:MAC_MAX_PARA,MAC_MAX_LEVEL)
      INTEGER            mac_level
      INTEGER            mac_line(MAC_MAX_LEVEL)
      INTEGER            mac_n_par(MAC_MAX_LEVEL)
      LOGICAL            lmakro
      LOGICAL            lmacro_dbg(MAC_MAX_LEVEL)
!
      CHARACTER*1024     io_file(MAC_MAX_IO)
      CHARACTER*20       io_out_format(MAC_MAX_FORM)
      LOGICAL            io_open(MAC_MAX_IO)
      INTEGER            io_unit(MAC_MAX_IO)
      LOGICAL            io_eof(MAC_MAX_IO)
      INTEGER                 io_get_sub(MAC_MAX_IO,2)
!
      COMMON /makro/ lmakro,mac_level,mac_line,mac_name,                &
     &               mac_leng,mac_n_par,mac_para,lmacro_dbg
      COMMON /io/    io_out_format,io_file,io_open,io_unit,             &
     &               io_eof,io_get_sub
