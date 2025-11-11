#!/bin/bash
#
#
rm -rf h5fortran-main
unzip -e h5fortran-main.zip
cd h5fortran-main
cmake -B build
cmake --build build
#
cd ..
mkdir -p ../DiffuseBuild/h5fortran/lib
mkdir -p ../DiffuseBuild/h5fortran/include
cp h5fortran-main/build/libh5fortran.a        ../DiffuseBuild/h5fortran/lib
cp h5fortran-main/build/include/h5fortran.mod ../DiffuseBuild/h5fortran/include
