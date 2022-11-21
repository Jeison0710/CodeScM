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
test: cleano modules.o nrtype.o gaussint.o evaluate_bsp_v2.o bspline_functions.o 
	$(FF90) -o tt writter.f90 modules.o nrtype.o gaussint.o evaluate_bsp_v2.o bspline_functions.o $(LIBS)
run:
	./tt < input
cleano: 
	rm *.o *.mod
clean: 
	rm *.o *.mod *~
