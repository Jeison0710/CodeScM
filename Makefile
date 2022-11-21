FF90=gfortran
#FF90=ifort
USR=/usr
# mkl library for ifort: lapack package
# LIBS = -L$MKLPATH -I$MKLINCLUDE -lmkl_lapack -lmkl -lguide -lpthread
# mkl library for version 11.3
#LIBS = -L$MKLPATH -I$MKLINCLUDE -lmkl_lapack -lmkl_intel -lguide -lpthread
#LIBS = -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
LIBS = -lblas -llapack
# lapack at PKS
# LIBS = -llapack -lblas -lg2c
#LIBS = -llapack -lblas 

nrtype.o:
	$(FF90) -c nrtype.f90
modules.o: nrtype.o
	$(FF90) -c modules.f90
gaussint.o: 
	$(FF90) -c gaussint.f
evaluate_bsp.o: nrtype.o modules.o
	$(FF90) -c evaluate_bsp.f90 
evaluate_bsp_v2.o: nrtype.o modules.o
	$(FF90) -c evaluate_bsp_v2.f90 
bspline_functions.o:
	$(FF90) -c bspline_functions.f90
matrices_pbc.o: modules.o nrtype.o
	$(FF90) -c matrices_pbc.f90
initial_wavepacketpbc.o: nrtype.o modules.o
	$(FF90) -c initial_wavepacketpbc.f90
backwardsubstitution-delande.o: nrtype.o modules.o
	$(FF90) -c backwardsubstitution-delande.f90
llt-delande.o: nrtype.o modules.o
	$(FF90) -c llt-delande.f90
subr_scalefactor.o: nrtype.o modules.o
	$(FF90) -c subr_scalefactor.f90
subr2_scalefactor.o: nrtype.o modules.o
	$(FF90) -c subr2_scalefactor.f90
test: cleano modules.o nrtype.o gaussint.o evaluate_bsp_v2.o bspline_functions.o initial_wavepacketpbc.o matrices_pbc.o llt-delande.o backwardsubstitution-delande.o subr_scalefactor.o
	$(FF90) -o tt schrodinger_poisson_bsplines.f90 modules.o nrtype.o gaussint.o evaluate_bsp_v2.o bspline_functions.o initial_wavepacketpbc.o matrices_pbc.o llt-delande.o backwardsubstitution-delande.o subr_scalefactor.o $(LIBS)

test2: cleano modules.o nrtype.o gaussint.o evaluate_bsp_v2.o bspline_functions.o initial_wavepacketpbc.o matrices_pbc.o llt-delande.o backwardsubstitution-delande.o subr2_scalefactor.o
	$(FF90) -o tt schrodinger_poisson_bsplines.f90 modules.o nrtype.o gaussint.o evaluate_bsp_v2.o bspline_functions.o initial_wavepacketpbc.o matrices_pbc.o llt-delande.o backwardsubstitution-delande.o subr2_scalefactor.o $(LIBS)
run:
	./tt < inschrodinger_poisson_bsplines.dat
cleano: 
	rm *.o *.mod
clean: 
	rm *.o *.mod *~

#gfortran -o tt schrodinger_poisson_bsplines.f90 modules.o nrtype.o gaussint.o evaluate_bsp_v2.o bspline_functions.o initial_wavepacketpbc.o matrices_pbc.o llt-delande.o backwardsubstitution-delande.o subr2_scalefactor.o -lblas -llapack

#gfortran -c subr_scalefactor.f90
#	gfortran -o sf subr_scalefactor.o -lblas -llapack
