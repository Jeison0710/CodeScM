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
!c   performs the L*Lt decomposition of a double complex banded matrix
!A.
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
!c              whose name starts with ladc
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
        !write(6,1)
1       format(' Now starts LLt decomposition of matrix A...')
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
79        continue
          ntl(i)=j
        enddo
        do i=1,ntot
          do j=MAX(1,i+nlarg-ntot),nlarg
            if (ntl(i+nlarg-j).le.j) go to 89
          enddo
89        continue
          ntl2(i)=j
        enddo
!c on va calculer la decomposition colonne par colonne
!c dans chaque colonne, on descend a partir de la diagonale
        do 5000 i=1,ntot
!c*********************************************
!c start computation of diagonal element A(i,i)
!c determination de l'element diagonal i,i
!c on somme les elements convenables
!c        j0=i-nlarg               
          xr=0.0_RT
          xi=0.0_RT
          do 1000 j=ntl(i),nlarg-1
            xr=xr+ar(j,i)**2-ai(j,i)**2
            xi=xi+2.0_RT*ar(j,i)*ai(j,i)
!c      write(6,999)xr,xi
!c999   format(2(1x,g22.15))
!c make x=x+a(j,i)**2
1000      continue     
!c if some element is too small, send a warning    
          xr=ar(nlarg,i)-xr
          xi=ai(nlarg,i)-xi
!c      write(6,999)xr,xi

!         xc=cdsqrt(CMPLX(xr,xi,CT))
!c  has been changed to
          xc=SQRT(CMPLX(xr,xi,CT))
!c  to support kind variables

!c      write(6,999)xc
          ar(nlarg,i)=REAL(xc,RT)
          ai(nlarg,i)=AIMAG(xc)
          xc=1.0_RT/xc
          x2r=REAL(xc,RT)
          x2i=AIMAG(xc)
!c      write(6,999)x2r,x2i
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
1500        continue
!c on en deduit l'element cherche
!c and substract it
            xr=ar(j,i0)-xr
            xi=ai(j,i0)-xi
            ar(j,i0)=xr*x2r-xi*x2i
            ai(j,i0)=xr*x2i+xi*x2r
!c          write(6,1501)j,i0,ar(j,i0),ai(j,i0)
!c1501      format(1x,i4,1x,i4,1x,g22.15,1x,g22.15)
!c make a(j,i0)=(a(j,i0)-x)*x2
2000      continue
!c end computation of non diagonal elements A(i,j)
!c ***********************************************
5000    continue    
!c
!c end decomposition of A
!c end routine
!c 
        return

        END SUBROUTINE lddc

