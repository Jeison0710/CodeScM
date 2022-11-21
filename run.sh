#!/bin/bash
pwd
cd code
#export OMP_NUM_THREADS=42


	gfortran -c nrtype.f90
	gfortran -c modules.f90
	gfortran -c gaussint.f
	gfortran -c evaluate_bsp.f90 
	gfortran -c evaluate_bsp_v2.f90 
	gfortran -c bspline_functions.f90
	gfortran -c matrices_pbc.f90
	gfortran -c initial_wavepacketpbc.f90
	gfortran -c backwardsubstitution-delande.f90
	gfortran -c llt-delande.f90
	gfortran -c subr_scalefactor.f90
	gfortran -c subr2_scalefactor.f90
#	gfortran -o tt schrodinger_poisson_bsplines.f90 modules.o nrtype.o gaussint.o evaluate_bsp_v2.o bspline_functions.o initial_wavepacketpbc.o matrices_pbc.o llt-delande.o backwardsubstitution-delande.o subr_scalefactor.o -lblas -llapack
	gfortran -o tt schrodinger_poisson_bsplines.f90 modules.o nrtype.o gaussint.o evaluate_bsp_v2.o bspline_functions.o initial_wavepacketpbc.o matrices_pbc.o llt-delande.o backwardsubstitution-delande.o subr2_scalefactor.o -lblas -llapack

	./tt < inschrodinger_poisson_bsplines.dat
	rm *.o *.mod
