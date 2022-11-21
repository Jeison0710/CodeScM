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
!c   xer,xei : double complex input vector (length ntot)(output: not
!changed)
!c   v1r,v1i : double complex working vector (length ntot) for
!intermediate
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
!c              starts with ladc...
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
1000      continue
          xc=CMPLX(xer(i)-xr,xei(i)-xi,CT)/CMPLX(ar(nlarg,i),ai(nlarg,i),CT)
          v1r(i)=REAL(xc,RT)
          v1i(i)=AIMAG(xc)
!c makes v1(i)=(xe(i)-x)/a(nlarg,i)
2000    continue
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
3100      continue
4100    continue
!c end second backward substitution
!c
!c********************************************************************
!c 
!c end routine
!c 
!       DEALLOCATE(ntl)
       return

      END SUBROUTINE axdc


