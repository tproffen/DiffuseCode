#!/bin/bash
gfortran -O3 -march=native -ffast-math -fopenmp \
  -o benchmark_math benchmark.f90
gfortran -O3 -march=native -fopenmp \
  -o benchmark_nomath benchmark.f90
