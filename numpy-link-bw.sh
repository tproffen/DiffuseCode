#!/bin/bash

/opt/gcc/4.8.2/bin/gfortran -Wall -Wall -shared build/temp.linux-x86_64-2.7/numpy/linalg/lapack_litemodule.o build/temp.linux-x86_64-2.7/numpy/linalg/python_xerbla.o -L/u/sciteam/wozniak/downloads/lapack-3.5.0/SRC -L/u/sciteam/wozniak/downloads/lapack-3.5.0/BLAS/SRC -L/u/sciteam/wozniak/sfw/python-2.7.6/lib -Lbuild/temp.linux-x86_64-2.7 -llapack -lblas -lpython2.7 -lgfortran -o build/lib.linux-x86_64-2.7/numpy/linalg/lapack_lite.so -Wl,-rpath -Wl,${HOME}/downloads/lapack-3.5.0/BLAS/SRC -Wl,-rpath -Wl,${HOME}/downloads/lapack-3.5.0/SRC
