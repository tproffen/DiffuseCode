from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules = [Extension("discus_cython", ["discus_cython.pyx"], 
                         extra_objects=["CMakeFiles/discus_all.dir/*.o",
                                        "../../lib_f90/liblib_f90.a",
                                        "../../lib_f90/liblib_f90c.a"],
                         libraries=["gfortran","readline"])]

setup(name = 'DiscusPy',
      cmdclass = {'build_ext': build_ext},
      ext_modules = ext_modules)
