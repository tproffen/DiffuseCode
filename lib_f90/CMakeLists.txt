# Create the lib_90 library for the command language

set (FSOURCES allocate_generic.f90 appl_unix.f90 blockdat.f90 dkdo.f90 
              errlist.f90 help.f90 jsu_readline.f90 learn.f90 macro.f90 
              random.f90 sup.f90 sys_gfortran.f90)

set (CSOURCES FCreadline.c csocket.c csup.c)

IF (WIN32 OR CYGWIN )
  set (CSOURCES ${CSOURCES} win32-glob.c)
ENDIF (WIN32 OR CYGWIN)

add_library (lib_f90  ${FSOURCES})
add_library (lib_f90c ${CSOURCES})
