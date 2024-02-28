#
# CMakeLists.txt file for finufft within DISCUS
# R. Neder 2023_05_25
# Uses the archive from "https://github.com/flatironinstitute/finufft/releases/tag/v*.*.*"
#  as subdirectory "finufft" within the project
# 
#
# The parent: "../ CMakeLists.txt" sets (among other things):
# find_package(FFTW REQUIRED)
# set (PROJECT_LIBS ${PROJECT_LIBS} ${FFTW_DOUBLE_LIB} ${FFTW_DOUBLE_OPENMP_LIB} )
# set (CMAKE_CXX_FLAGS "-O3 -fPIC -funroll-loops -march=native -fcx-limited-range -std=c++14 -fopenmp")
# set (CMAKE_C_FLAGS "-O3 -fPIC")
# add_subdirectory(finufft)
# Here PROJECT_LIBS is a placeholder for other libraries that might be required by the project
#
# On MacOS the gnu compiler suite is installed into /usr/local/bin/
#    and is enforced prior to creating the project with :
#    "export CC=/usr/local/bin/gcc-12" or whatever version 
#    "export CXX=/usr/local/bin/g++-12" or whatever version 
# It seems to be important to do this prior to the first creation of the project or else
# the cache will still have clang ...
#
# To ease installation of precompiled versions I prefer a static library

set (CPPSOURCES src/finufft.cpp
                src/simpleinterfaces.cpp 
		src/spreadinterp.cpp 
                src/utils.cpp 
		src/utils_precindep.cpp
	        fortran/finufftfort.cpp
                contrib/legendre_rule_fast.cpp	)

include_directories(include contrib ${FFTW_INCLUDE_DIRS})

if (DIFFUSE_SHARED_LIB)
  add_library (finufft_lib SHARED ${CPPSOURCES})
else (DIFFUSE_SHARED_LIB)
  add_library (finufft_lib ${CPPSOURCES})
endif (DIFFUSE_SHARED_LIB)
#add_library (finufft_lib_s SHARED ${CPPSOURCES})

if (DIFFUSE_SHARED_LIB)
  install (TARGETS finufft_lib DESTINATION lib)
endif (DIFFUSE_SHARED_LIB)
