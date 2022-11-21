
!c*****************************************************************c!
!c                                                                 c!
!c              pladcbrssv.f                                       c!
!c                                                                 c!
!c   version 4 du 5-10-93                                          c!
!c                                                                 c!
!c   version 4.ss1 B - double real stored in sparse-symmetric way  c!
!c                                                                 c!
!c   compared to pladcbrimv - new multiplication by a sparse B     c!
!c   matrix in routine bdrssrc.f                                   c!
!AK2/2/2000   now bcdrssc for complex b-matrix!!!!!!!! ! ! ! !! !  c!
!c                                 original bxdrc.f removed and    c!
!c   playing a dummy role for compatibility                        c!
!c                                                                 c!
!c   package for solving a generalized eigenvalue problem          c!
!c   (A-lambda*B)x=0                                               c!
!c   where A is a double complex symmetric banded matrix           c!
!c   and B a double precision symmetric banded matrix.             c!
!c   This package computes the eigenvalues close to zero.          c!
!c   It computes the overlap of the eigenvectors with some fixed   c!
!c   initial vectors.                                              c!
!c                                                                 c!
!c   To obtain the eigenvalues close to some initial number,       c!
!c   first shift the A matrix, then run this package and add       c!
!c   the shift to all the eigenvalues.                             c!
!c   (This is the responsability of the user and not provided in   c!
!c   this package).                                                c!
!c                                                                 c!
!c*****************************************************************c!
!c                                                                 c!
!c   Author : Dominique Delande                                    c!
!c            Laboratoire de Spectroscopie Hertzienne              c!
!c            de l'Ecole Normale Superieure                        c!
!c            Tour 12 - Etage 1                                    c!
!c            Universite Pierre et Marie Curie                     c!
!c            4, Place Jussieu                                     C!
!c            F-75252 PARIS CEDEX 05                               C!
!c            FRANCE                                               c!
!c                                                                 c!
!c            Phone: (33)144274404                                 c!
!c            Fax:   (33)144273845                                 c!
!c            E-mail : delande@spectro.jussieu.fr                  c!
!c                     or dod@frunip62.bitnet                      c!
!c                                                                 c!
!c*****************************************************************c! 
!c                                                                 c!
!c   Version 4.                                                    c!
!c                                                                 c!
!c   October 1992                                                  c!
!c   Copyright Dominique Delande                                   c!
!c                                                                 c!
!c*****************************************************************c!
!c                                                                  !
!   This package is almost entirely written with double precision 
!   variables and
!   arrays, rather than using double complex numbers. It has be found to
!   be more efficient on all computers.
!   Consequently, the double complex variables and arrays are divided in
!   real and imaginary part, which can be recognized with their
!   suffixes r and i (in the names of the variables).
!   Example : ar(i,j) and ai(i,j) are respectively the real and
!   imaginary parts of the (i,j) element of the double complex array a.
!   In this package, such a pair will be called an element of
!   a double complex array.
!
!c********************************************************************
!c
!   To use this package, the user needs to do the
!   4 following things:
!
!
!   1) Build his own include file and add it at the beginning of
!      each routine of the package,
!      i.e. replace the lines "include 'inrm98l.h'" by
!      "include 'myfilesizeparn1n2.inc'", where myfile.h is the include file.
!      Example of include file:
! 
!	integer ngr,nlmax,ndmax,npas,ndz
!       parameter (ngr=1000,nlmax=100,ndmax=100,npas=100,ndz=5)
!	double precision epsloc
!       parameter(epsloc=4.d-15)
!
!      These parametrers have the following meaning:
!		ngr : maximum size of the matrix A and B
!               nlmax : maximum bandwidth of matrix A
!               ndmax : maximum bandwidth of matrix B
!               npas : maximum number of Lanczos steps
!               ndz : maximum number of vectors whose overlaps with
! 		      converged eigenvectors are wanted
!		epsloc : a small number depending on the machine
!                        (4.d-15 is ok for Cray, it can be smaller on
!                         usual workstations) 
!
! 
!   2) Add at the beginning of the calling routine the include file
!      and the needed common blocks:
! 
!	double precision ar,ai
!       common/detax/ar(nlmax,ngr),ai(nlmax,ngr)
!	double precision d
!	common/detbx/d(ndmax,ngr)
!	double precision zr,zi
!	common/vp/zr(ngr,ndz),zi(ngr,ndz)
!	integer nconv
!	double precision tconvr,tconvi,tzconvr,tzconvi
!	common/conv/tconvr(npas),tconvi(npas),
!     &tzconvr(npas,ndz),tzconvi(npas,ndz)
!      
!
!   3) Fill the A matrix
!           the B matrix 
!           (these two matrices have to be filled using the symmetric
!            banded matrices storage convention)
!           the vectors z whose overlaps with the eigenvectors are wanted
!	    (in the first ndze columns of array z)
!      in the calling routine,
!      then make the two following calls:
!
!	  call lddc(ntot,nlarg)
!         call ladcbrimv(ntot,nlarg,nlard,np,nort,ndze,nconv,iscal)
!
!
!   INPUT VARIABLES
!
!   INTEGER
!   ntot : size of the banded matrix (output: not changed)
!   nlarg : bandwidth of the matrix A (output: not changed)
!   nlard : bandwidth of the matrix B (output: not changed)
!   np : number of Lanczos steps to perform (output : not changed)
!   nort : number of reorthogonalization steps to perform
!          Usually equal to np.
!          (output : not changed)
!   ndze : number of vectors whose overlap with the converged eigenvectors
!          is wanted (output : not changed)
!   iscal : determines wether the overlap is computed directly (iscal=1)
!           or with the scalar product defined by the B matrix (iscal=0).
!
!   OUTPUT
!
!   INTEGER
!   nconv : number of converged eigenvalues
!
!   double precision ARRAY
!
!   The results (eigenvalues and overlaps) are output 
!   in the common : 
!	common/conv/tconvr(npas),tconvi(npas),
!     &tzconvr(npas,ndz),tzconvi(npas,ndz)
!
!   tconvr(nconv),tconvi(nconv) : converged eigenvalues
!   (WARNING : starting from version 2.0, the eigenvalues themselves
!   are the output. This is incompatible with version 1.xx, where
!   the inverses of the eigenvalues were the output).
!
!   tzconvr(npas,ndz),tzconvi(npas,ndz) : components of the converged 
!                                   eigenstates on the input vectors in z
!				    Actual size nconv*ndze 
!
!   INPUT VARIABLES TRANSMITTED BY COMMONS
!
!c   /detax/ar(nlmax,ngr),ai(nlmax,ngr) : double complex banded A matrix
!c                                        actual size : nlarg*ntot
!c                                        (ouput : changed)
!c
!c   /detbx/d(ndmax,ngr)		       : double precision banded B matrix
!c					 actual size : nlard*ntot
!c                                        (output : changed)
!c
!c In banded form, the lines of ar,ai and dr,di contain the subdiagonals
!c of the A and B matrices with the upper left corner not used.
!c The diagonal is in the last line of the matrix, 
!c e.g. (d(nlard,i),i=1,ntot)).
!c The element A(i,j) (i>=j) is stored in ar(nlarg+j-i,i) (real part) 
!c and ai(nlarg+j-i,i) (imaginary part). Of course, A(j,i) is at the 
!c same place.
!c B(i,j) (i>=j) is stored in d(nlard+j-i,i).
!c
!c THE A AND B MATRICES HAVE TO BE FILLED BY THE USER IN THIS
!c BANDED FORM, BEFORE ANY CALL TO THE ROUTINES OF THE PACKAGE:
!c
!c   /vp/zr(ngr,ndz),zi(ngr,ndz) : double complex array
!c                         actual size : ntot*ndze
!c	                  Used for computing the overlap of the 
!c                         eigenvectors with some initially fixed 
!c                         vectors. On input, it has to be filled by the
!c                         user with the fixed vectors into the first
!c			  ndze columns.
!c                         (output : changed if iscal=0)
!c				    not changed if iscal=1)
!c
!c
!c   4) Compile and link the calling routine(s) with the package.
!c      
!c   The program is now ready to run.
!c
!c   Enjoy it!
!c   Prenez beaucoup de plaisir!
!c   Viel Spass!
!c
!c*********************************************************************
!c
!c   Routines used : lddc.f
!c                   ladcbrimv.f calling axdc.f
!c                                       bxdrc.f
!c                                       qrdcim.f
!c					cvdc.f
!c                                       tridcimv.f
!c
!c*********************************************************************
!c   
!c   Possible changes : 
!c   
!c   1) If no overlaps are needed, just set ndz and ndze to zero.
!c
!c   2) If only one overlap is needed, there isa special version pladcbri1.f
!c
!c   3) If the eigenvectors are needed, use pladcbrimv.f
!c
!c   4) If B has some special form (tridiagonal for exemple),
!c      a different form of storage can be used. The only modification
!c      is in the bxdrc.f routine which has to be adapted, and in the
!c      calls of this routine. It is called at three different places,
!c      in the ladcbrim.f routine. The /detbx/ common has also to be 
!c      changed.
!c
!c   5) More generally, if the user is able to provide an efficient
!c      routine making the product Bx, where x is an arbitrary input
!c      vector (for exemple, using the sparsity of B), he can replace
!c      the present bxdrc.f routine by his own routine. The B matrix 
!c      is not needed outside this bxdrc routine. As in 4), the calls
!c      in ladcbrim.f and the /detbx/ common has to be modified.
!c
!c*********************************************************************
!c
!c   subroutine lddc.f
!c
!c   version 4.0
!c   11-10-92
!c   by dominique delande
!c
!c********************************************************************
!c
!c   used in the lanczos algorithm version 4
!c
!c   performs the L*Lt decomposition of a double complex banded matrix A.
!c
!c*********************************************************************
!c
!c   INPUT
!c
!c   INTEGER
!c   ntot : size of the banded matrix (output: not changed)
!c   nlarg : bandwidth of the matrix A (output: not changed)
!c   
!c   OUTPUT
!c
!c   none
!c
!c********************************************************************
!c
!c   COMMONS
!c
!c   /detax/ar(nlmax,ngr),ai(nlmax,ngr) : double complex banded matrix
!c                                        actual size : nlarg*ntot
!c                                        (ouput : changed)
!c
!c*******************************************************************
!c
!c   Routines called : none
!c
!c******************************************************************** 
!c   
!c   Called by : user, generally before calling a lanczos routine
!c		whose name starts with ladc
!c               This routine has to be called only once before
!c               any call to axdc.f
!c
!c********************************************************************
!c 
!c This subroutine is used
!c by the Lanczos algorithm for a generalized eigenvalue problem
!c with banded A matrix. A is double complex.
!c The A matrix is input in the usual banded form,
!c In the present routine, A is decomposed as L2*L2t stored in 
!c banded form in a (bandwidth nlarg). It is later used by the 
!c routine axdc.f
!c to solve a generalized linear equation by backward substitution.
!c
!c In banded form, the lines of a contain the subdiagonals
!c of the A matrix with the upper left corner not used.
!c The diagonal is in the last line of the matrix, 
!c e.g. (a(nlarg,i),i=1,ntot)).
!c The element A(i,j) (i>=j) is stored in ar(nlarg+j-i,i) (real part)
!c and ai(nlarg+j-i) (imaginary part).
!c Of course, A(j,i) is at the 
!c same place.
!c
!c This routine is vectorisable, although the gain is not
!c the best one.
!c
!c For huge matrices, where only few eigenvalues are obtained
!c (few Lanczos steps << bandwidth), most of the time is spent
!c in this decomposition.
!c
!c Starting from version 4, the (possible) variable bandwidth
!c of matrix A is used. If (as it is often the case), the bandwidth
!c is not constant over the whole matrix A, much time should be lost
!c when computing the LDLt decomposition outside the effective
!c bandwidth, where matrix elements are all equal to zero.
!c In this new version, the local bandwidth of the A matrix is 
!c recorded in arrays ntl and ntl2. ntl(i) contains the index
!c of the first non-zero element on a given line i, ntl2(i)
!c the index of the first non-zero element on column i.
!c These two indices are recorded in banded matrix storage form.
!c ntl2 is not needed in the following and is consequently
!c a local variable. ntl is needed for the backward substitution
!c and is recorded in a common named varindex.
!c The variable bandwidth is determined in this routine
!c from the zeros found in the A matrix,
!c the user does not need to provide it.
!c
!c********************************************************************
!c********************************************************************
!c********************************************************************


	SUBROUTINE lddc(ntot,nlarg)
        USE nrtype
	use detax
	use varindex

	IMPLICIT none

!c arguments
	INTEGER(ITT) :: ntot,nlarg 
!c local variables
	INTEGER(ITT) :: i,j,k,i0,k0
	REAL(RT) :: xr,xi,x2r,x2i,x,x2
	COMPLEX(CT) :: xc

	INTEGER(ITT):: ntl2(0:nnggrr)
!c
!c********************************************************************
!c
!c routine starts here
!c
!c *******************************************************************
!c performs some checks
!c on the dimensions
	write(6,1)
1	format(' Now starts LLt decomposition of matrix A...')
!#####################################################################!
	ALLOCATE(ntl(0:ntot))
!#####################################################################!
!c perform decomposition of A first
!c determination of the effective bandwidth
	do i=1,ntot
	  do j=MAX(1,nlarg+1-i),nlarg
!c looking for the zeros in A
	    if ((ar(j,i).ne.0.0_RT).or.(ai(j,i).ne.0.0_RT)) go to 79
	  enddo
!c first non-zero element found
79	  continue
	  ntl(i)=j
	enddo
	do i=1,ntot
	  do j=MAX(1,i+nlarg-ntot),nlarg
	    if (ntl(i+nlarg-j).le.j) go to 89
	  enddo
89	  continue
	  ntl2(i)=j
	enddo
!c on va calculer la decomposition colonne par colonne
!c dans chaque colonne, on descend a partir de la diagonale
	do 5000 i=1,ntot
!c*********************************************
!c start computation of diagonal element A(i,i)
!c determination de l'element diagonal i,i
!c on somme les elements convenables
!c 	  j0=i-nlarg               
	  xr=0.0_RT
	  xi=0.0_RT
	  do 1000 j=ntl(i),nlarg-1
	    xr=xr+ar(j,i)**2-ai(j,i)**2
	    xi=xi+2.0_RT*ar(j,i)*ai(j,i)
!c	write(6,999)xr,xi
!c999	format(2(1x,g22.15))
!c make x=x+a(j,i)**2
1000	  continue     
!c if some element is too small, send a warning    
	  xr=ar(nlarg,i)-xr
	  xi=ai(nlarg,i)-xi
!c	write(6,999)xr,xi

!         xc=cdsqrt(CMPLX(xr,xi,CT))
!c  has been changed to
	  xc=SQRT(CMPLX(xr,xi,CT))
!c  to support kind variables

!c	write(6,999)xc
	  ar(nlarg,i)=REAL(xc,RT)
	  ai(nlarg,i)=AIMAG(xc)
	  xc=1.0_RT/xc
	  x2r=REAL(xc,RT)
	  x2i=AIMAG(xc)
!c	write(6,999)x2r,x2i
!c make a(nlarg,i)=sqrt(a(nlarg,i)-x)
!c and  x2=1/a(nlarg,i)
!c end computation of diagonal element A(i,i)
!c*********************************************
!c start computation of non diagonal elements A(i,j)
!c on va calculer l'element de la ligne i+nlarg-j et 
!c de la colonne i
!c on soustrait les elements deja calcules
	  do 2000 j=ntl2(i),nlarg-1
	    i0=i+nlarg-j
	    k0=nlarg-j
	    xr=0.0_RT
	    xi=0.0_RT
!c compute the sum of previously obtained terms
	    do 1500 k=MAX(ntl(i0),ntl(i)-k0),j-1
	      xr=xr+ar(k,i0)*ar(k+k0,i)-ai(k,i0)*ai(k+k0,i)
	      xi=xi+ar(k,i0)*ai(k+k0,i)+ai(k,i0)*ar(k+k0,i)
!c make x=x+a(k,i0)*a(k+k0,i)
1500	    continue
!c on en deduit l'element cherche
!c and substract it
	    xr=ar(j,i0)-xr
	    xi=ai(j,i0)-xi
	    ar(j,i0)=xr*x2r-xi*x2i
	    ai(j,i0)=xr*x2i+xi*x2r
!c	    write(6,1501)j,i0,ar(j,i0),ai(j,i0)
!c1501	    format(1x,i4,1x,i4,1x,g22.15,1x,g22.15)
!c make a(j,i0)=(a(j,i0)-x)*x2
2000	  continue
!c end computation of non diagonal elements A(i,j)
!c ***********************************************
5000	continue    
!c
!c end decomposition of A
!c end routine
!c 
	return

	END SUBROUTINE lddc


!c*********************************************************************
!c*********************************************************************
!c*********************************************************************
!c*********************************************************************
!c
!c   subroutine ladcbrimv.f
!c
!c   version 4.0
!c   27-11-92
!c   by dominique delande
!c
!c********************************************************************
!c
!c   used in the lanczos algorithm version 4
!c
!c   performs the Lanczos algorithm itself for a generalized
!c   eigenvalue problem (A-lambda*B)x=0 for banded symmetric matrices
!c   A is double complex and B double precision.
!c   with computation of some fixed component of the eigenvectors.
!c   The eigenvalues are sorted by
!c   ascending real part.
!c   It works on the inverse of the original problem (previously
!c   shifted) in order to get the eigenvalues close to some 
!c   initial shift.
!c   The A and B matrices are not used in explicit form.
!c   It requires a user-supplied routine axdc.f
!c   solving the linear problem Ax=b.
!c   Usually, this is done by backward substitution in axdc.f
!c   using the L*Lt decompositions of A computed in lddc.f
!c   Starting from version 3, the B matrix itself is not needed,
!c   but only a user-supplied subroutine which computes
!c   the action of B on some arbitrary input vector.
!c   This routine is here named bxdrc.f.
!c   bxdrc can be of any form the user wants.
!c   The currently tested forms are banded (with variable bandwidth),
!c   tridiagonal and unit matrices.
!c   In this package, a version of bxdrc for double precision banded
!c   matrix is provided.
!c   
!c
!c*********************************************************************
!c
!c   INPUT
!c
!c   INTEGER
!c   ntot : size of the banded matrix (output: not changed)
!c   nlarg : bandwidth of the matrix A (output: not changed)
!c   nlard : bandwidth of the matrix B (output: not changed)
!c   np : number of Lanczos steps to perform (output : not changed)
!c   ndze : number of initial vectors whose overlap with the converged
!c          eigenstates will be computed (output : not changed)
!c   nort : number of reorthogonalization steps to perform
!c          (output : not changed)
!c   iscal : integer specifying if the overlap of the eigenstates
!c           with the fixed ndze vectors has to be evaluated
!c           with the B scalar product (iscal=0) or with
!c           the standard scalar product (iscal=1).
!c	    If iscal=0, the input vectors are just multiplied
!c           by B at the beginning.
!c
!c   OUTPUT
!c
!c   INTEGER
!c   nconv : number of converged eigenvalues
!c
!c   double precision ARRAY
!c 
!c   tconvr(nconv),tconvi(nconv) : converged eigenvalues
!c
!c   tzconvr(npas,ndz),tzconvi(npas,ndz) : components of the converged 
!c                        eigenstates on some input vectors. 
!c                        Only the first nconv lines and ndze columns
!c			 are filled.
!c
!c********************************************************************
!c
!c   COMMONS
!c
!c   /bidq/qr(ngr,npas),qi(ngr,npas) : double complex array
!c                                     actual size : ntot*np
!c                                     each column contains a Lanczos
!c                                     vector used during the 
!c                                     calculation 
!c                                     (ouput : changed)
!c
!!c   /vp/zr(ngr,ndz),zi(ngr,ndz) : double complex array
!c                         actual size : ntot*ndze
!c			  On iput, it contains in its first ndze
!c			  columns the ndze vectors whose overlap
!c			  with the eigenstates is wanted
!c                         (output : not changed if iscal=1
!c				        changed if iscal=0)
!c
!c*******************************************************************
!c
!c   Routines called : axdc.f to solve the linear problem
!c                     bxdrc.f to compute Bx 
!c                     qrdcim.f to diagonalize the tridiagonal
!c                               Lanczos algorithm 
!c                     tridcimv.f to sort the eigenvalues
!c		      cvdc.f
!c
!c******************************************************************** 
!c   
!c   Called by : user, generally after call to lddc.f
!c                     to perform L*Lt decomposition
!c
!c********************************************************************
!c 
!c This subroutine is used
!c by the Lanczos algorithm for a generalized eigenvalue problem
!c with banded A and B matrices. A is double complex and B double complex.
!c The A matrix has been decomposed as L2*L2t stored in banded form
!c in a (bandwidth nlarg).
!c This decomposition is usually performed in
!c a routine lddc...
!c
!c In banded form, the lines of ar,ai contain the subdiagonals
!c of the A matrix with the upper left corner not used.
!c The diagonal is in the last line of the matrix, 
!c e.g. (ar(nlarg,i),ai(nlarg,i),i=1,ntot)).
!c The element A(i,j) (i>=j) is stored in ar(nlarg+j-i,i) (real part) 
!c and ai(nlarg+j-i,i) (imaginary part). Of course, A(j,i) is at the 
!c! same place.
!c
!c This routine performs np steps of the Lanczos algorithm,
!c using some initial guess for the first vector,
!c (here, a uniform decomposition on the whole basis)
!c with full reorthogonalization of the Lanczos vectors
!c during the first nort steps. It is usually good to set nort=np.
!c This routines makes np calls to the axdc and bxdrc routine which solves
!c a linear problem by backward substitution. For this reason,
!c the A and B matrices are not required in the present routine.
!c After np steps, the resulting np*np tridiagonal matrix is 
!c diagonalized using a QR algorithm (calls qrdcim.f). To monitor
!c the accuracy of the computed eigenvalues, QR computes the
!c component of all the eigenvectors of the tridiagonal matrix
!c (NOT the eigenvectors of the initial problem) on the last vector
!c of the np*np matrix. As shown in Parlett and Scott, this makes it
!c possible to discriminate between converged and non converged
!c eigenvalues. In addition, the components on some other
!c input vectors are computed, which makes it possible to compute
!c later the overlap of the computed eigenvectors with some
!c initially fixed vectors.
!c The converged eigenvalues are then sorted by ascending real part,
!c 
!c The most consuming part is the reorthogonlization. For a large
!c number of steps (larger than the bandwidth), most of the CPU
!c time is spent in this stupid operation. Very probably, this
!c can be improved by calls to optimized library functions to 
!c perform scalar products.
!c
!c********************************************************************
!c********************************************************************
!c********************************************************************

        SUBROUTINE ladcbrv(ntot,nlarg,nlard,np,nort,nconv,nnggrr)
        USE nrtype
	use conv
	use bidq
	use work
        use varindexb

	IMPLICIT none

	INTEGER(ITT) :: nnggrr
!c arguments
        INTEGER(ITT) :: ntot,nlarg,nlard,np,nort,nconv
!c local variables
	INTEGER(ITT) :: i,j,k,iconv     
	REAL(RT) :: betjr,betji,betjm1r,betjm1i,alpjr,alpji
	REAL(RT) :: qjaqjr,qjaqji,rjrjr,rjrji
	REAL(RT) :: scalr,scali,xr,xi
	REAL(RT) :: qnorm
	REAL(RT) :: eps,epscon
	COMPLEX(CT) :: xc
!c arrays td and ts (double precision and imaginary parts)
!c will contain the diagonal and subdiagonal of the Lanczos
!c tridiagonal matrix (ts(1) is set to 0 and never used)
!c tz is used to compute the overlap with the initially fixed
!c vector contained in z.
!c aux is an array used to compute the overlap of the
!c eigenvectors of the tridiagonal Lanczos matrix with
!c the last vector of the basis, and is used to monitor the
!c accuracy of the computed eigenvalues.
!c rj, v1 and v2 are working arrays of actual length ntot.


	REAL(RT) :: tdr(np),tsr(np),auxr(np)
	REAL(RT) :: tdi(np),tsi(np),auxi(np)

	REAL(RT) :: rjr(nnggrr),rji(nnggrr) 
	REAL(RT) :: v1r(nnggrr),v1i(nnggrr),v2r(nnggrr)
	REAL(RT) :: v2i(nnggrr)

	REAL(RT) :: scalpr,scalpi
        
        INTEGER(ITT), PARAMETER :: preci=PRECISION(eps)
        REAL(RT), parameter :: epsloc=4.0_RT*10.0_RT**(-preci)
!        REAL(RT), parameter :: epsloc=4.e-15


	REAL(RT) :: tscalr(0:np),tscali(0:np)
!	double precision:: tscalr(0:npas),tscali(0:npas)
!c these eps control the accuracy of the computed eigenvalues
!	data eps/1.d-15/
!	data epscon/1.d-6/
	eps=10.0_RT**(-preci)
	epscon=10.0_RT**(-(preci-3)/2)

!c
!c*********************************************************************
!c
!c routine starts here
!c
!c initialize some variables
!#####################################################################!
	ALLOCATE(qr(ntot,np))
	ALLOCATE(qi(ntot,np))
!#####################################################################!

	write(6,1)
1	format(' Now starts Lanczos algorithm...')

	  tdr=0.0_RT
	  tdi=0.0_RT
	  tsr=0.0_RT
	  tsi=0.0_RT
	  tscalr=0.0_RT
	  tscali=0.0_RT


	tscalr(0)=0.0_RT
	tscali(0)=0.0_RT
        betjr=0.0_RT
	betji=0.0_RT            
!c initialise norm of the Lanczos matrix to 0.
	qnorm=0.0_RT
!c end initialisation
!c set first Lanczos vector to 1 for all the components

	v1r=1.0_RT
	v1i=0.0_RT

!        do 100 i=1,ntot
!          v1r(i)=1.d0
!	  v1i(i)=0.d0
!100     continue
!c normalize this vector with respect to B

! ====================================================
!              ************     CHANGES BY J. MADRONERO    **********
! 24.5.02       NTLB(i) (i=1,2,...,NTOT) contains the index of the first 
!                                 non-zero element on a given line i of the matrix B.
! 24.5.02       If NTLB(0)=-1 the array NTLB hast not been defined yet.
                 ALLOCATE(ntlb(0:ntot))
                 ntlb(0)=-1
!AK2/2/2000	CALL bdrssrc(ntot,nlard,v1r,v1i,v2r,v2i)
! 24.5.02	      CALL bcdrssrc(ntot,nlard,rjr,rji,v2r,v2i)
! 24.5.02       Last routine was changed by the next one. 
!                    product of the matrix B on the vector V1R+I*V1I
                 CALL bxdrcak(ntot,nlard,v1r,v1i,v2r,v2i)
!
!
!
                 xr=0.0_RT
                 xi=0.0_RT
!	xr=DOT_PRODUCT(v1r,v2r)-DOT_PRODUCT(v1i,v2i)
!	xi=DOT_PRODUCT(v1r,v2i)+DOT_PRODUCT(v1i,v2r)
!ak9/6/99   xr,xi numbers

!ak9/6/99   MAYBE THE FOLLOWING LOOPS ARE FASTER USING 
!AK9/6/99   THE INTRINSIC F90 FCT DOT_PRODUCT(VECTOR_A,VECTOR_B)
!AK9/6/99 (---------------->MANCHESTER p.57)
!
	do 110 i=1,ntot
	  xr=xr+v1r(i)*v2r(i)-v1i(i)*v2i(i)
	  xi=xi+v1r(i)*v2i(i)+v1i(i)*v2r(i)
110	continue
!

	xc=1.0_RT/SQRT(CMPLX(xr,xi,CT))
	xr=REAL(xc,RT)
	xi=AIMAG(xc)
!c q (,1) will now contain the properly normalized first Lanczos vector
!c and v1 B times this vector
	qr(:,1)=v1r*xr-v1i*xi
	qi(:,1)=v1r*xi+v1i*xr
	v1r=v2r*xr-v2i*xi
	v1i=v2r*xi+v2i*xr
!ak9/6/99 v1r,v1i vectors


!	do 120 i=1,ntot
!	  qr(i,1)=v1r(i)*xr-v1i(i)*xi
!	  qi(i,1)=v1r(i)*xi+v1i(i)*xr
!	  v1r(i)=v2r(i)*xr-v2i(i)*xi
!	  v1i(i)=v2r(i)*xi+v2i(i)*xr
!120	continue
!c 
!c***********************************************
!c starts loop on lanczos steps     
!c iteration sur le nombre de pas de l'algorithme

        do 1000 j=1,np
	   write(6,1001)j
 1001	   format(' Starting Lanczos step ',i4)
!c at the beginning of each step, v1
!c contains the action of B on the previous 
!c Lanczos vector, but not reorthogonalized
!c against the previous ones
!c on decale les differents coefficients
	   betjm1r=betjr
	   betjm1i=betji
!c calcul de l'action de A**(-1) sur le vecteur v1
!c resultat dans rj
!c v2 sert de stockage intermediaire

	   CALL axdc(ntot,nlarg,v1r,v1i,rjr,rji,v2r,v2i)

!c on purge pour orthogonaliser
!c here only rj is reorthogonalized
!c all this is done in such a way
!c that the routine btxdr is called only once
!c per Lanczos step
	   do 7300 k=1,j-1
!ak                    j is the actual number of Lanczos step
!c tscal contains the residual scalar products
!c of the last Lanczos vector with the previous ones
!c it has been computed at the end of the previous step
!c and should be zero in exact arithmetic
!c cas generique
	      scalpr=tscalr(k)*tdr(k)-tscali(k)*tdi(k)+&
                   tscalr(k-1)*tsr(k)-tscali(k-1)*tsi(k)+&
                   tscalr(k+1)*tsr(k+1)-tscali(k+1)*tsi(k+1)
	      scalpi=tscalr(k)*tdi(k)+tscali(k)*tdr(k)+&
                   tscalr(k-1)*tsi(k)+tscali(k-1)*tsr(k)+&
                   tscalr(k+1)*tsi(k+1)+tscali(k+1)*tsr(k+1)
!c xr,xi contains the inverse of the previous betaj
	      scalr=scalpr*xr-scalpi*xi
	      scali=scalpr*xi+scalpi*xr

		 rjr=rjr-scalr*qr(:,k)+scali*qi(:,k)
		 rji=rji-scalr*qi(:,k)-scali*qr(:,k)
!ak9/6/99         (rjr, rji vectors)

!	      do 7200 i=1,ntot
!		 rjr(i)=rjr(i)-scalr*qr(i,k)+scali*qi(i,k)
!		 rji(i)=rji(i)-scalr*qi(i,k)-scali*qr(i,k)
! 7200	      continue
 7300	   continue	
!c determination du produit scalaire <qj|A|qj>
!ak9/6/99   MAYBE THE FOLLOWING LOOP IS FASTER USING 
!AK9/6/99   THE INTRINSIC F90 FCT DOT_PRODUCT(VECTOR_A,VECTOR_B)
!AK9/6/99 (---------------->MANCHESTER p.57)
	   qjaqjr=0.0_RT
	   qjaqji=0.0_RT
!	      qjaqjr=DOT_PRODUCT(v1r,rjr)-DOT_PRODUCT(v1i,rji)
!	      qjaqji=DOT_PRODUCT(v1r,rji)+DOT_PRODUCT(v1i,rjr)
!ak9/6/99      (qjajr and qjaji numbers)
!
!     USING THIS FCT. ONE GETS SLIGHTLY DIFFERENT RESULTS...
!        ........................
!
!
	   do 1100 i=1,ntot
	      qjaqjr=qjaqjr + v1r(i)*rjr(i) -v1i(i)*rji(i)
	      qjaqji=qjaqji + v1r(i)*rji(i) +v1i(i)*rjr(i)
!c make qjaqj=qjaqj + v1(i)*rj(i)
 1100	   continue  
!c <qj|A|qj> is now in variable alpj                              
	   alpjr=qjaqjr
	   alpji=qjaqji
	   if (j.eq.1) then
	      qnorm=ABS(alpjr)+ABS(alpji)
	   endif
!c affichage de alpj et betj
!c	  write(6,210)j,alpjr,alpji,betjr,betji
!c210	  format(' j = ',i3,4(5x,g22.15))
!c on remplit l'element courant de la matrice T
!c td est l'element diagonal et ts l'element non-diagonal
!c attention, il existe un decalage de une unite entre l'indice
!c betj et l'element de ts
	   tdr(j)=alpjr
	   tdi(j)=alpji
!c si on est au dernier pas, on peut s'arreter la
!c sinon, on continue
	   if (j.ne.np) then
!c on forme alors rj=Aqj-bet(j-1)*q(j-1)-alpj*qj
!c to compute q(j+1)
	      if (j.ne.1) then   
!c cas generique

		    rjr=rjr-betjm1r*qr(:,j-1)+betjm1i*qi(:,j-1)&
                         -alpjr*qr(:,j)+alpji*qi(:,j)
		    rji=rji-betjm1r*qi(:,j-1)-betjm1i*qr(:,j-1)&
                         -alpjr*qi(:,j)-alpji*qr(:,j)
!
!                         rjr,rji vectors
!                 
!		 do 1200 i=1,ntot
!		    rjr(i)=rjr(i)-betjm1r*qr(i,j-1)+betjm1i*qi(i,j-1)
!     &-alpjr*qr(i,j)+alpji*qi(i,j)
!		    rji(i)=rji(i)-betjm1r*qi(i,j-1)-betjm1i*qr(i,j-1)
!     &-alpjr*qi(i,j)-alpji*qr(i,j)
!c make rj(i)=rj(i)-betjm1*q(i,j-1)-alpj*q(i,j)
 1200		 continue
	      else            
!c cas j=1 ou on a un effet de bord

		    rjr=rjr-alpjr*qr(:,1)+alpji*qi(:,1)
		    rji=rji-alpjr*qi(:,1)-alpji*qr(:,1)

!                         rjr,rji vectors

!		 do 1201 i=1,ntot
!		    rjr(i)=rjr(i)-alpjr*qr(i,1)+alpji*qi(i,1)
!		    rji(i)=rji(i)-alpjr*qi(i,1)-alpji*qr(i,1)
!c make rj(i)=rj(i)-alpj*q(i,1)
! 1201		 continue
	      endif
!c on calcule l'action de B sur rj

!AK2/2/2000	      CALL bdrssrc(ntot,nlard,rjr,rji,v2r,v2i)
! 24.5.02	      CALL bcdrssrc(ntot,nlard,rjr,rji,v2r,v2i)
! 24.5.02       Last routine was changed by the next one. 
! 24.5.02       Aqui no es necesario definir ntlb(0)=-1.        
!                     Calculus of the action of the matrix B on the vector V1R+I*V1I
                 CALL bxdrcak(ntot,nlard,rjr,rji,v2r,v2i)

!c pour cela, on calcule les produits scalaires
	      do 4900 k=1,MIN(j-1,nort)
		 scalr=0.0_RT
		 scali=0.0_RT
!ak9/6/99   MAYBE THE FOLLOWING LOOPS ARE FASTER USING 
!AK9/6/99   THE INTRINSIC F90 FCT DOT_PRODUCT(VECTOR_A,VECTOR_B)
!AK9/6/99 (---------------->MANCHESTER p.57)

!		    scalr=DOT_PRODUCT(qr(:,k),v2r)-
!     &                    DOT_PRODUCT(qi(:,k),v2i)
!		    scali=DOT_PRODUCT(qr(:,k),v2i)+
!     &                    DOT_PRODUCT(qi(:,k),v2r)


		 do 4200 i=1,ntot
		    scalr=scalr+qr(i,k)*v2r(i)-qi(i,k)*v2i(i)
		    scali=scali+qr(i,k)*v2i(i)+qi(i,k)*v2r(i)
!c make scal=scal+q(i,k)*v2(i)
 4200		 continue
!c the scalar products are stored in tscal
!c and reused at the beginning of the next step
		 tscalr(k)=scalr
		 tscali(k)=scali
 4900	      continue	
!c on calcule la norme de rj 
	      rjrjr=0.0_RT
	      rjrji=0.0_RT
!9/6/99
!		 rjrjr = DOT_PRODUCT(rjr,v2r)-DOT_PRODUCT(rji,v2i)
!		 rjrji = DOT_PRODUCT(rjr,v2i)+DOT_PRODUCT(rji,v2r)
!9/6/99	      

	      do 1300 i=1,ntot 
		 rjrjr = rjrjr + rjr(i)*v2r(i) - rji(i)*v2i(i)
		 rjrji = rjrji + rjr(i)*v2i(i) + rji(i)*v2r(i)
 1300	      continue 
!c test if the norm of the vector is zero
!c this usually indicates too many Lanczos steps
!c or that A**n(q1) only spans some subspace
!c It produces a warning           
	      if ((ABS(rjrjr)+ABS(rjrji)).lt.(qnorm*epsloc)) then
		 rjrjr=qnorm*epsloc
		 rjrji=0.0_RT
		 write(6,1301)j
 1301		 format(' Warning : small number at Lanczos step #',i4,/,&
     ' This may cause overflow or inaccurate results')
	      endif
	      xc=SQRT(CMPLX(rjrjr,rjrji,CT))
	      betjr=REAL(xc,RT)
	      betji=AIMAG(xc)
	      xc=1.0_RT/xc
	      xr=REAL(xc,RT)
	      xi=AIMAG(xc)
!c make betj=sqrt(rjrj)
!c on divise rj par betj pour former q(j+1)
!c en fait, on multiplie par 1/betj=x
!c on multiplie aussi v2 par x pour avoir
!c l'action de B sur le vecteur de Lanczos dans v1
!c qui est reutilise au pas suivant

		 qr(:,j+1)=rjr*xr-rji*xi
		 qi(:,j+1)=rjr*xi+rji*xr
		 v1r=v2r*xr-v2i*xi
		 v1i=v2r*xi+v2i*xr
	      

!	      do 1400 i=1,ntot  
!		 qr(i,j+1)=rjr(i)*xr-rji(i)*xi
!		 qi(i,j+1)=rjr(i)*xi+rji(i)*xr
!		 v1r(i)=v2r(i)*xr-v2i(i)*xi
!		 v1i(i)=v2r(i)*xi+v2i(i)*xr
! 1400	      continue	
!c end computation of q(j+1)	
!c fill the subdiagonal element
	      tsr(j+1)=betjr
	      tsi(j+1)=betji
	   endif
	  qnorm=MAX(qnorm,ABS(betjr)+ABS(betji)+ABS(betjm1r)+&
               ABS(betjm1i)+ABS(alpjr)+ABS(alpji))
1000	continue
!c end of the Lanczos steps
!c********************************
!c start diagonalisation of the tridiagonal Lanczos matrix
!c 
!c fill the auxiliary array to compute overlap with the last
!c vector of the basis
!	do 2000 i=1,np
!  	  auxr(i)=0.d0
!	  auxi(i)=0.d0
!2000	continue   
	auxr=0.0_RT
	auxi=0.0_RT
	auxr(np)=1.0_RT
!#####################################################################!
	    ALLOCATE(workvr(np,np))
	    ALLOCATE(workvi(np,np))
!#####################################################################!
	write(6,2010)
2010	format(' Now starts QR diagonalization...')
!c diagonalisation de la matrice triangulaire
!c eps is the accuracy of the computed eigenvalues

	CALL qrdcv(np,eps,tdr,tdi,tsr,tsi,auxr,auxi)

!cab 201098
	write(*,*) 'qrdcv finished'
!c selection des valeurs propres convergees et elimination 
!c des multiplets
!c end QR diagonalization
!c *****************************
!#####################################################################!
!c start deleting non converged eigenvalues
!c the converged ones are put in the first elements
!c of the vectors
	iconv=1
	do 5000 i=1,np
!c test if eigenvalue is "converged" with a variant 
!c of the Parlett and Scott criterion
!c in fact, it is simple to consider "converged"
!c all the eigenstates of the Lanczos matrix
!c having an overlap with the last vector smaller
!c than some fixed small number epscon
!c usually, the error on the eigenvalue
!c will be of the order of epscon**2*(mean spacing)
!c this may cause a problem 
!c only if two eigenvalues are nearly degenerate

	  if ((ABS(auxr(i))+ABS(auxi(i))).lt.epscon) then

!c cas ou une vp est consideree comme convergee
!c on detecte si la valeur propre existe deja
!c try to detect possible doublets
!c as this does not happen with full reorthogonalization
!c these are very dangerous lines which may delete real
!c converged eigenvalues in case of quasidegenracy.
!c DO NOT USE.
!c	    if (iconv.ne.1) then
!c  	  if ((dabs(tdr(i)-tconvr(iconv-1))+dabs(tdi(i)-tconvi(iconv-1)))
!c      &.gt.(epsdef*(dabs(tdr(i))+dabs(tdi(i))))) then        
!c eigenvalue is Ok
!c invert it and store in tconv
!c store the intensity in tzconv

	    xc=1.0_RT/CMPLX(tdr(i),tdi(i),CT)
	    tconvr(iconv)=REAL(xc,RT)
	    tconvi(iconv)=AIMAG(xc)

!c on reordonne les lignes de la matrice workv
!c qui contient les rotations qr accumulees

	    do 4920 j=1,np
	      workvr(iconv,j)=workvr(i,j)
	      workvi(iconv,j)=workvi(i,j)
4920	    continue
	    iconv=iconv+1
	  endif
5000	continue
	nconv=iconv-1
!c end deleting non converged eigenvalues
!c***************************************
!c start sorting the converged eigenvalues in tconv,

	CALL tridcv(nconv,np)

!c end sorting the eigenvalues
!c***************************************
!c start computing eigenvectors
!c now apply matrix workv only on converged
!c Lanczos vectors and get the eigenvectors in ascending order

	CALL cvdc(ntot,np,nconv,auxr,auxi)

!c end computing eigenvectors
!c***************************************
!#####################################################################!
	DEALLOCATE(workvr)
	DEALLOCATE(workvi)
!#####################################################################!

!c	write(6,5010)iconv
!c5010	format(' nombre de valeurs propres convergees = ',i4)	
!c	do 3000 i=1,iconv
!c	write(6,2010)i,tconv(i)
!c2010	format(1x,i3,3x,g22.15)
!c3000	continue	
!c	write(6,7001)
!c7001	format(' End of Lanczos routine')
!c
!c********************************************************************
!c
!c end of routine
!c
        return

        END SUBROUTINE ladcbrv


!c*********************************************************************
!c*********************************************************************
!c
!c   subroutine axdc.f
!c
!c   version 4.0
!c   11-10-92
!c   by dominique delande
!!c
!c********************************************************************
!c
!c   used in the lanczos algorithm version 4
!c
!c   solves the problem Ax=b by backward substitution
!c   for an eigenvalue problem with A complex.
!c
!c   A is a banded symmetric matrix, b an arbitrary input
!c   vector and x the output solution
!c
!c*********************************************************************
!c
!c   INPUT
!c
!c   INTEGER
!c   ntot : size of the banded matrix (output: not changed)
!c   nlarg : bandwidth of the matrix A (output: not changed)
!c   
!c   double precision ARRAY
!c   xer,xei : double complex input vector (length ntot)(output: not changed)
!c   v1r,v1i : double complex working vector (length ntot) for intermediate
!c             storage (output: changed)
!c
!c   OUTPUT
!c
!c   double precision ARRAY
!c   xsr,xsi : double complex output vector (length ntot)
!c
!c********************************************************************
!c
!c   COMMONS
!c
!c   /detax/ar(nlmax,ngr),ai(nlmax,ngr) : double complex banded matrix
!c                                        actual size : nlarg*ntot
!c                                        (ouput : not changed)
!c
!c*******************************************************************
!c
!c   Routines called : none
!c
!c******************************************************************** 
!c   
!c   Called by : all the Lanczos routines whose name
!c		starts with ladc...
!c
!c********************************************************************
!c 
!c This subroutine solves a linear problem by backward substitution
!c for the Lanczos algorithm for a banded symmetric matrix A
!c A is double complex.
!c The A matrix has been decomposed as L2*L2t, stored in banded form
!c in ar,ai (bandwidth
!c nlarg). This decomposition is usually performed in
!c a routine lddc...
!c
!c In banded form, the lines of ar,ai contain the subdiagonals
!c of the A matrix with the upper left corner not used.
!c The diagonal is in the last line of the matrix, 
!c e.g. (ar(nlarg,i),ai(nlarg,i),i=1,ntot)).
!c The element A(i,j) (i>=j) is stored in ar(nlarg+j-i,i) (real part) 
!c and ai(nlarg+j-i,i) (imaginary part). Of course, A(j,i) is at the 
!c same place.
!c
!c The routine performs successively the actions of L2**(-1) and
!c L2t**(-1) on the initial vector.
!c the two steps are vectorizable using either accumulation
!c or direct sum.
!c This is extremely important, as most of the time is spent
!c in this routine.
!c Starting from version 4, the (possible) variable bandwidth
!c of matrix A is used. If (as it is often the case), the bandwidth
!c is not constant over the whole matrix A, much time should be lost
!c when computing the backward substitution outside the effective
!c bandwidth, where matrix elements are all equal to zero.
!c In this new version, the local bandwidth of the A matrix is 
!c recorded in array ntl. ntl(i) contains the index
!c of the first non-zero element on a given line i.
!c This index is recorded in banded matrix storage form.
!c The variable bandwidth is determined in the routine axdc
!c from the zeros found in the A matrix,
!c the user does not need to provide it.
!c It is passed to the present routine through the varindex common
!!c
!c********************************************************************
!c*********************************************************************
!c*********************************************************************

        SUBROUTINE axdc(ntot,nlarg,xer,xei,xsr,xsi,v1r,v1i)
        USE nrtype
	use detax
	use varindex

	IMPLICIT none

!c arguments
	INTEGER(ITT) :: ntot,nlarg

        REAL(RT) :: xer(ntot),xsr(ntot),v1r(ntot)
        REAL(RT) :: xei(ntot),xsi(ntot),v1i(ntot)
!c local variables
	INTEGER(ITT) :: i,j,j0
	REAL(RT) :: xr,xi
	COMPLEX(CT) :: xc
!c
!c********************************************************************
!c 
!c routine starts here
!c start first backward substitution
!c formation du produit           
!c on forme tout d'abord le produit L2**(-1)
!c que l'on place dans xs          
	do 2000 i=1,ntot
	  xr=0.0_RT
	  xi=0.0_RT
	  j0=i-nlarg
!c direct loop with first index increasing linearly
!c ok for vectorisation
	  do 1000 j=ntl(i),nlarg-1
	    xr=xr+ar(j,i)*v1r(j+j0)-ai(j,i)*v1i(j+j0)
	    xi=xi+ar(j,i)*v1i(j+j0)+ai(j,i)*v1r(j+j0)
!c makes x=x*a(j,i)*v1(j+j0)
1000	  continue
	  xc=CMPLX(xer(i)-xr,xei(i)-xi,CT)/CMPLX(ar(nlarg,i),ai(nlarg,i),CT)
	  v1r(i)=REAL(xc,RT)
	  v1i(i)=AIMAG(xc)
!c makes v1(i)=(xe(i)-x)/a(nlarg,i)
2000	continue
!c end first backward substitution
!!c **********************************
!c start second backward substitution
!c on fait ensuite le produit par L2**(-t)
!c et on place le resultat dans xs
!c nouvelle version rusee du 9-11-90
!c the computed contributions are substracted
!c from the upper components as soons as they are
!c obtained
!c this accumulation makes it possible efficient vectorization 
	do 4100 i=ntot,1,-1
	  xc=CMPLX(v1r(i),v1i(i),CT)/CMPLX(ar(nlarg,i),ai(nlarg,i),CT)
	  xr=REAL(xc,RT)
	  xi=AIMAG(xc)
	  xsr(i)=xr
	  xsi(i)=xi
	  j0=nlarg-i
	  do 3100 j=ntl(i)-j0,i-1
	    v1r(j)=v1r(j)-xr*ar(j+j0,i)+xi*ai(j+j0,i)
	    v1i(j)=v1i(j)-xr*ai(j+j0,i)-xi*ar(j+j0,i)
3100	  continue
4100	continue
!c end second backward substitution
!c
!c********************************************************************
!c 
!c end routine
!c 
!	DEALLOCATE(ntl)
       return

      END SUBROUTINE axdc

!c*********************************************************************
!c*********************************************************************
!c
!c   subroutine qrdcv.f
!c
!c   version 4.0
!c   11-10-92
!c   by dominique delande
!c
!c********************************************************************
!c
!c   used in the lanczos algorithm version 4
!c
!c   Performs the diagonalization of the symmetric complex
!c   tridiagonal Lanczos matrix using the QR algorithm.
!c   In this version, the components of the eigenstates on some
!c   fixed vectors are computed.
!c   Usually called by the lanczos routines starting with ladc,
!c   and ending by im.f.
!!c   
!c*********************************************************************
!c
!c   INPUT
!c
!c   INTEGER
!c   n : dimension of the matrix (=number of Lanczos steps)
!c       (output : not changed)
!c   ndze : number of input vectors whose overlap with the eigenstates
!c	   is wanted (output : not changed)
!c
!c   REAL
!c   eps : accuracy of the computed eigenvalues
!c         (ouput : not changed)
!c
!c   double precision ARRAY
!c   ar(n),ai(n) : double complex array containing the diagonal of the input
!c                 matrix (output : contains the eigenvalues not sorted)
!c   br(n),bi(n) : double complex array containing the subdiagonal of the
!c                 input matrix starting at b(2)
!c                 (b(1) is arbitrary and not used)
!!c                 (output : destroyed)
!c!   v1r(n),v1i(n) : double complex array containing some initial vector
!c                   (output : contains the overlap of the eigenvectors
!c                   with the initial vector)
!c
!c********************************************************************
!c
!c   COMMONS
!c
!c   /qr/tzr(npas,ndz),tzi(npas,ndz) : double complex array
!c				containing some initial vectors
!c				in its first ndze columns
!c				(output : contains the overlap
!c				of these vectors with the eigenstates)
!c  
!c*******************************************************************
!c
!c   Routines called : none
!c
!c******************************************************************** 
!c   
!c   Called by : lanczos routine whose name starts 
!c               with ladc and ends with im
!c
!c
!c********************************************************************
!c 
!c The QR algorithm performs successive orthogonal (eventually complex)
!c elementary rotations (i.e., involving two consecutive lines and 
!c columns) to transform an initial tridiagonal matrix in a diagonal
!c one. At each step, the non diagonal elements are decreased.
!c The process stops when all non diagonal elements are smaller than
!c some epsilon, ensuring this accuracy to the eigenvalues.
!c At each step, the overlap with fixed
!c initial vectors is computed, accumulating
!c the changes step after step.
!c To increase the speed, the whole matrix is shifted at each
!c step, but of course restored at the end.
!c
!c********************************************************************
!c*********************************************************************
!c*********************************************************************


        SUBROUTINE qrdcv(n,eps,ar,ai,br,bi,v1r,v1i)
        USE nrtype
	use work

	IMPLICIT none

!c arguments
        INTEGER(ITT) :: n
        REAL(RT) ::eps
        REAL(RT) :: ar(n),br(n),v1r(n)
        REAL(RT) :: ai(n),bi(n),v1i(n)
!c local variables
	REAL(RT) :: cr,sr,pr,alr,blr
	REAL(RT) :: ci,si,pi,ali,bli
        REAL(RT) :: difr,tr,dr,pambdar,gr,hr,shiftr
        REAL(RT) :: difi,ti,di,pambdai,gi,hi,shifti
	REAL(RT) :: rr,ri
	REAL(RT) :: small,anorm
	COMPLEX(CT) :: dc,tc
        INTEGER(ITT) :: i,j,k,id,if,iff,ifprec,idd   
!c
!c********************************************************************
!c
!c routine starts here
!c
!c       integer indic
!c       do 1 i=1,n
!c       write(4,2)a(i),b(i)
!c2      format(4(3x,g15.8))
!c1      continue
!c       indic=0
!c cas trivial
        if (n.le.1) return
!c initialisation de la matrice de transformation
!c and some other variables
	do 10 i=1,n
	  do 9 j=1,n
!c	    workvr(j,i)=0.
!c	    workvi(j,i)=0.
	    workvr(j,i)=0.0_RT
	    workvi(j,i)=0.0_RT
9	  continue
          workvr(i,i)=1.0_RT
!c          workvr(i,i)=1.0
10	continue
!cab	write(*,*) 'loop 10'
        ifprec=n
	if=n
        id=0
        shiftr=0.0_RT
        shifti=0.0_RT
	br(1)=0.0_RT
	bi(1)=0.0_RT
	anorm=0.0_RT
!c end initialisation
!c**************************************
!c start computation of the norm 
	do 80 i=1,n-1
	anorm=MAX(anorm,ABS(ar(i))+ABS(ai(i))+ABS(br(i))+&
             ABS(bi(i))+ABS(br(i+1))+ABS(bi(i+1)))
80	continue
!cab	write(*,*) 'loop 80'
	anorm=MAX(anorm,ABS(ar(n))+ABS(ai(n))+ABS(br(n))+ABS(bi(n)))
!c end computation of the norm 
	small=eps*anorm
!cab	write(*,*) 'small,anorm=',small,anorm
!c**************************************
!c starts QR itself
!c on procede par deflation a partir du bas et du haut
!c des qu'un element est plus petit que small, on deflate la matrice
101     continue
        id=id+1
!cab	write(*,*)'id=',id
110     continue
!c if the subdiagonal element is less than small, deflate the matrix
        if ((ABS(br(id+1))+ABS(bi(id+1))).le.small) go to 190
        do 120 if=ifprec,id+1,-1
          if ((ABS(br(if))+ABS(bi(if))).gt.small) go to 130
!c restore original values by adding shift
          ar(if)=ar(if)+shiftr
          ai(if)=ai(if)+shifti
120     continue
!cab	write(*,*) 'loop 120'
130     continue
!cab	write(*,*) 'if=',if
!c calcul d'une valeur approchee
!c on utilise ainsi un accelerateur de convergence
!c       do 11 i=1,n
!c       write(4,2)a(i),b(i)
!c11      continue
        ifprec=if
        iff=if-1
        idd=id+1
!c start estimation of best possible shift
!cab	write(*,*) 'if,iff,ar(if)=',if,iff,ar(if)
        difr=ar(iff)-ar(if)
        difi=ai(iff)-ai(if)
!cab	write(*,*) 'difr'
	rr=4.0_RT*(br(if)**2-bi(if)**2)+difr**2-difi**2
	ri=8.0_RT*br(if)*bi(if)+2.0_RT*difr*difi
!cab	write(*,*) 'rr'
	dc=SQRT(CMPLX(rr,ri,CT))
	dr=REAL(dc,RT)
	di=AIMAG(dc)
!cab	write(*,*) 'best possible shift'
!c make d=sqrt(4.0*b(if)*b(if)+dif*dif)
        IF (difr.GT.0.0_RT) then
	  Dr=-Dr
	  di=-di
	endif
        pambdar=0.50_RT*(ar(if)+ar(iff)+dr)
        pambdai=0.50_RT*(ai(if)+ai(iff)+di)
        shiftr=shiftr+pambdar
        shifti=shifti+pambdai
!cab	write(*,*) 'made dsqrt'
!c end estimation of shift, added to the previous one
!c substract shift to the "active" part of the matrix 
        do 104 i=id,if
          ar(i)=ar(i)-pambdar
          ai(i)=ai(i)-pambdai
104     continue
!cab	write(*,*) 'loop 104'
        pr=Ar(ID)
        pi=Ai(ID)
!c initialise the elementary rotation
!c at zero angle c=1, s=0.
        Cr=1.0_RT
	ci=0.0_RT
        sr=0.0_RT
        si=0.0_RT
        do 105 i=idd,if
!c start elementary QR step
!c calcul des cbl generes par la methode qr
!c voir par exemple, exercice 60 de Ralston
!c        indic =indic+1
          gr = cr * br(i) - ci * bi(i)
          gi = cr * bi(i) + ci * br(i)
          hr = cr * pr - ci * pi
          hi = cr * pi + ci * pr
	  rr = br(i)**2-bi(i)**2+pr**2-pi**2
	  ri = 2.0_RT*(br(i)*bi(i)+pr*pi)
	  tc=SQRT(CMPLX(rr,ri,CT))
	  tr=REAL(tc,RT)
	  ti=AIMAG(tc)
!c make t=sqrt(b(i)**2+p**2)
          br(i-1)=tr*sr-ti*si
          bi(i-1)=ti*sr+tr*si
	  tc=1.0_RT/tc
	  tr=REAL(tc,RT)
	  ti=AIMAG(tc)
	  cr=pr*tr-pi*ti
	  ci=pr*ti+pi*tr
	  sr=br(i)*tr-bi(i)*ti
	  si=br(i)*ti+bi(i)*tr
!c make c=p/t
!c make s=b(i)/t
          pr = cr * ar(i) - ci * ai(i) - sr * gr + si * gi
          pi = cr * ai(i) + ci * ar(i) - sr * gi - si * gr
 	  rr=cr*gr-ci*gi+sr*ar(i)-si*ai(i)
 	  ri=cr*gi+ci*gr+sr*ai(i)+si*ar(i)
	  ar(i-1)=hr+sr*rr-si*ri
	  ai(i-1)=hi+sr*ri+si*rr
!c make a(i-1) = h + s * ( c * g + s * a(i) )
!c perform same elementary rotation on input vectors
!c cbl sur le vecteur propre
          alr=v1r(i-1)
	  ali=v1i(i-1)
	  blr=v1r(i)
	  bli=v1i(i)
	  v1r(i-1)=alr*cr-ali*ci+blr*sr-bli*si
	  v1i(i-1)=alr*ci+ali*cr+blr*si+bli*sr
	  v1r(i)=blr*cr-bli*ci-alr*sr+ali*si 
	  v1i(i)=blr*ci+bli*cr-alr*si-ali*sr
!c now perform the same operation in the accumulated
!c transform matrix 
!c le 8-11-91 : nouvelle methode
!c ou l'on accumule les rotations dans la matrice
!c workv de taille n*n avant d'appliquer
!c le tout sur la matrice des vecteurs de Lanczos
          do 193 k=1,n
	    alr=workvr(i-1,k)
            blr=workvr(i,k)
	    ali=workvi(i-1,k)
            bli=workvi(i,k)
            workvr(i-1,k)=alr*cr-ali*ci+blr*sr-bli*si
            workvi(i-1,k)=alr*ci+ali*cr+blr*si+bli*sr
            workvr(i,k)=blr*cr-bli*ci-alr*sr+ali*si  
            workvi(i,k)=blr*ci+bli*cr-alr*si-ali*sr
193	  continue
!cab	  write(*,*) 'loop 193'
!c end of elementary QR step
105     continue
!cab	write(*,*) 'loop 105'
!c modif the last elements
        Ar(IF)=pr*Cr-pi*ci
        Ai(IF)=pr*Ci+pi*cr
        br(IF)=pr*Sr-pi*si
        bi(IF)=pr*Si+pi*sr
        go to 110
190     continue
!cab	write(*,*) 'exit 190'
!c if element is deflated, restore by adding the
!c shift
        ar(id)=ar(id)+shiftr
        ai(id)=ai(id)+shifti
        if (if.ne.(id+1)) go to 101
        ar(if)=ar(if)+shiftr
        ai(if)=ai(if)+shifti
!c fin du calcul
!c       write(6,201)indic
!c201    format(' nombre de cbl = ',i6)
!c end of QR algorithm
!c*********************
!c
!c end of routine
!c
!c********************************************************************
        RETURN
 
	END SUBROUTINE qrdcv



!c*********************************************************************
!c*********************************************************************
!c
!c   subroutine tridcv.f
!c
!c   version 4.0
!c   10-11-92
!c   by dominique delande
!c
!c********************************************************************
!c
!c   used in the lanczos algorithm version 4
!c
!c   Performs the sorting of the converged eigenvalues for
!c   a double complex eigenvalue problem, by ascending real parts.
!c   Perform the same sort on some matrix elements
!c
!c*********************************************************************
!c
!c   INPUT
!c
!c   INTEGER
!c   nconv : number of converged eigenvalues (output: not changed)
!c
!c   ndze : number of matrix elements sorted simultaneously 
!c				(output : not changed)
!c
!c********************************************************************
!c
!c   COMMONS
!c
!c   /conv/tconvr(npas),tconvi(npas),
!c     &tzconvr(npas,ndz),tzconvi(npas,ndz)
!c	tconvr,tconvi : double complex array containing
!c			the eigenvalues to be sorted (output : changed)
!c	tzconvr,tzconvi : double complex array containing
!c			  the matrix elements simultaneously sorted
!c			  actual size : nconv*ndze			
!c			  (output : changed)
!c  
!c*******************************************************************
!c
!c   Routines called : none
!c
!c******************************************************************** 
!c   
!c   Called by :  lanczos routines whose name starts 
!c                with ladc and ends with im.
!c
!c********************************************************************
!c 
!c tri des valeurs propres par valeur croissante
!c on trie simultanement les composantes contenues dans tzconv
!c 
!c*********************************************************************
!c*********************************************************************


        SUBROUTINE tridcv(nconv,np)
        USE nrtype
	use conv
	use work

	IMPLICIT none

!c arguments  
	INTEGER(ITT) :: nconv,np
!c local variables
        INTEGER(ITT) :: i,j,k
	REAL(RT) :: a,b


	REAL(RT) :: aux2r(np),aux2i(np)


!	double precision:: aux2r(npas),aux2i(npas)

!c
!c********************************************************************
!c
!c routine starts here 
!c
!c the eigenvalues are sorted by ascending real part
!c debut classement des vp

        do 300 j=2,nconv
    	  a=tconvr(j)
	  b=tconvi(j)
	  do 290 k=1,np
            aux2r(k)=workvr(j,k)
            aux2i(k)=workvi(j,k)
290	  continue
	  do 100 i=j-1,1,-1
	    if (tconvr(i).le.a) go to 200
	    tconvr(i+1)=tconvr(i)
	    tconvi(i+1)=tconvi(i)
            do 90 k=1,np
	      workvr(i+1,k)=workvr(i,k)
	      workvi(i+1,k)=workvi(i,k)
90	    continue
100	  continue
	  i=0  
200	  continue
	  tconvr(i+1)=a
	  tconvi(i+1)=b
	  do 190 k=1,np
	    workvr(i+1,k)=aux2r(k)
	    workvi(i+1,k)=aux2i(k)
190	  continue
300	continue
!c
!c*********************************************************************
!c
!c end routine
!c
        return

        END SUBROUTINE tridcv

!  ======================================================
!  ======================================================
!                     SUBROUTINE bxdrcak(ntot,nlard,xer,xei,xsr,xsi)
!  
!  Realiza el producto de la matriz "real" B, guardada en la variable COMUN 
!  D(NLARD,NTOT) de acuerdo con la convencion para matrices simetricas,
!  y el vector XER+I*XEI. El resultado es el vector XSR+I*XSI. 
!                                                VARIABLES
!  NTOT: dimension de B
!  NLARD: ancho de B
!  XER: parte real del vector de entrada (de dimension NTOT)
!  XER: parte imaginaria del vector de entrada (de dimension NTOT)
!  XSR: parte real del vector de salida (de dimension NTOT)
!  XSR: parte imaginaria del vector de salida (de dimension NTOT)
!             informacion sobre el ancho efectivo de la columna I de D que se 
!             que se calcula en esta rutina siempre y cuando NTLB(0)=-1
!
!  ======================================================
!  ======================================================

 SUBROUTINE bxdrcak(ntot,nlard,xer,xei,xsr,xsi)
   USE nrtype
   use detbx
   use varindexb

   INTEGER(ITT) :: ntot,nlard
   REAL(RT), DIMENSION(ntot) :: xer,xei,xsr,xsi
!   dimension xer(ntot),xsr(ntot)
!   dimension xei(ntot),xsi(ntot)

   INTEGER(ITT) :: i,j,j0
   REAL(RT) :: xr,xi

!  ======================================================
!  at the first call, ntlb(0)=-1 and the effective bandwidth of B is computed
! Then ntlb(0).ne.-1 and the computed bandwidth is used for later calls
!  ======================================================
   if(ntlb(0).eq.-1) then    
      write(6,1)
1     format(' Computes effective bandwidth of matrix B...')
      ntlb(0)=nlard-1
      do  i=1,ntot
         do  j=MAX(1,nlard+1-i),nlard
            if (d(j,i).ne.0.0_RT) exit
         end do
         ntlb(i)=MIN(j,ntlb(i-1)+1)
      end do
   endif
   do  i=1,ntot
      xr=0.0_RT
      xi=0.0_RT
      j0=i-nlard
      do j=ntlb(i),nlard
         xr=xr+xer(j+j0)*d(j,i)
         xi=xi+xei(j+j0)*d(j,i)
      end do
      xsr(i)=xr
      xsi(i)=xi	
   end do

   do  i=1,ntot
      xr=xer(i)
      xi=xei(i)
      j0=nlard-i
      do  j=ntlb(i)-j0,i-1
         xsr(j)=xsr(j)+xr*d(j+j0,i)
         xsi(j)=xsi(j)+xi*d(j+j0,i)
      end do
   end do
   

   return
 END SUBROUTINE bxdrcak



!c**********************************************************************
!c**********************************************************************
!c
!c   subroutine cvdc.f
!c
!c   version 2.0
!c   20-1-92
!c   by dominique delande
!c
!c********************************************************************
!c
!c   used in the lanczos algorithm version 2
!c
!c   Computes the eigenvectors of a generalized eigenvalue problem
!c   with banded matrices. A is complex, B is real tridiagonal.
!c   Used after the Lanczos tridiagonal matrix has been
!c   diagonalised by qrsci2v.f
!c
!c*********************************************************************
!c
!c   INPUT
!c
!c   INTEGER
!c   n : number of Lanczos steps performed, i.e.
!c       dimension of the tridiagonal matrix, which
!c       has been diagonalized previously.
!c       (ouput : not changed)
!c   ntot : dimension of the initial matrix
!c          equal to the length of the Lanczos vectors.
!!c          (output : not changed)
!c!   nconv : number of converged eigenvalues and 
!c           eigenvectors. Only the first nconv lines
!c           of matrix workv are used.
!c           (output : not changed)
!c
!c   REAL ARRAY
!c   br(n),bi(n) : working complex vector of length >= nconv,
!c                 used to vectorize the computation.
!c                 (output : changed)
!!c 
!c ********************************************************************
!c
!c   COMMONS
!c
!c   /bidq/qr(ngr,npas),qi(ngr,npas) : complex array
!c                                     actual size : ntot*np
!c                                     each column contains
!c                                     an eigenvector as output
!c                                     only the first nconv columns
!c                                     are significant (converged
!c                                     eigenvectors), the other ones
!c                                     being used during the computation
!c                                     (ouput : changed)
!c
!c   /workvr(npas,npas),workvi(npas,npas) :
!c                         complex array; actual size nconv*n
!c                         working array used for the computation
!c                         of eigenvectors.
!c                         It is set by the QR algorithm qrsci2v.f
!c                         to the transformation matrix that
!!c                         diagonalize the tridiagonal Lanczos matrix.
!c!                         It is then used here to compute the
!c                         converged eigenvectors.
!c                         (output : not changed)
!c 
!c*******************************************************************
!c
!c   Routines called : none
!c
!c******************************************************************** 
!c   
!c   Called by : lascbtri1v.f
!c               lascbtrv.f
!c
!c********************************************************************
!c 
!c This subroutine is used
!c by the Lanczos algorithm for a generalized eigenvalue problem
!c with banded A and B matrices. A is complex and B real tridiagonal.
!c These two matrices are 
!c decomposed as L2*L2t and L1*L1t
!c respectively, stored in banded form in ar,ai (bandwidth nlarg)
!c for A and in d (diagonal) and e (subdiagonal) for B.
!c
!c In banded form, the lines of ar,ai contain the subdiagonals
!c of the A matrix with the upper left corner not used.
!c The diagonal is in the last line of the matrix, 
!c e.g. (ar(nlarg,i),i=1,ntot)).
!c The element A(i,j) (i>=j) is stored in ar(nlarg+j-i,i) (real part) 
!c and ai(nlarg+j-i,i) (imaginary part). Of course, A(j,i) is at the 
!c same place.
!c
!c This routine computes the eigenvectors after
!c the Lanczos tridiagonal matrix has been diagonalized
!c by the subroutine qrsci2v.f
!c and sorted using trisci1v.f
!c it applies the accumulated rotations in array workv
!c on the Lanczos vectors contained in q
!c in order to get only converged eigenvectors.
!c As the eigenvalues have previously been sorted,
!c and the non converged ones deleted, the array 
!c workv has been modified and only the first nconv
!c lines are used to generate the converged eigenvectors.
!c this requires NSTEP*NCONV*N operations.
!c
!c******************************************************************** 
!c**********************************************************************
!c**********************************************************************

	SUBROUTINE cvdc(ntot,n,nconv,br,bi)
        USE nrtype
	use bidq
	use work

	IMPLICIT none

!c arguments
	INTEGER(ITT) :: n,ntot,nconv
	REAL(RT) :: br(nconv),bi(nconv)
!c local variables
	INTEGER(ITT) :: i,j,k
	REAL(RT):: xr,xi
!c*********************************************************************
!c
!c routine starts here
!c
!c starts apply workv on the converged Lanczos vectors
!c in order to get the eigenvetors
!c on applique tout d'abord la transformation sur les vecteurs
!c de Lanczos
!c on utilise le vecteur  b
!c comme espace de stockage temporaire
!c pour vectoriser correctement
!	do 4100 k=1,nconv
!	  br(k)=0.
!	  bi(k)=0.
!4100	continue
	br=0.0_RT
	bi=0.0_RT
	do 4800 i=1,ntot
	  do 4700 j=1,n
	    xr=qr(i,j)
	    xi=qi(i,j)
!c accumulate transformation in b
!	      br=br+workvr(1:nconv,j)*xr-workvi(1:nconv,j)*xi
!	      bi=bi+workvr(1:nconv,j)*xi+workvi(1:nconv,j)*xr
	    do 4600 k=1,nconv
               br(k)=br(k)+workvr(k,j)*xr-workvi(k,j)*xi
               bi(k)=bi(k)+workvr(k,j)*xi+workvi(k,j)*xr
 4600       continue
4700	  continue
!c copy b in the ith line of q
	  do 4750 k=1,nconv
	    qr(i,k)=br(k)
	    qi(i,k)=bi(k)
!	    br(k)=0.
!	    bi(k)=0.
	    br(k)=0.0_RT
	    bi(k)=0.0_RT


4750	  continue
4800 	continue
!c end apply workv on the converged Lanczos vectors
!c*************************************************
!c
!c end routine
!c 
	return

	END SUBROUTINE cvdc
