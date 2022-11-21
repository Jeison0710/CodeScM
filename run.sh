#!/bin/bash
pwd
cd code
export OMP_NUM_THREADS=20
gfortran runer.f90 -o parallel -fopenmp 
time ./parallel
