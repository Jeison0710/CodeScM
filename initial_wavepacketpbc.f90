!*******************************************************************
! subroutine initial_wavepacketpbc
! psi0=1/sqrt(sigma*sqrt(pi))*exp(-(x-x0)**2/(2*sigma**2))*
!      exp(i*p0*x)
! nmax: number of basis functions. These are bsplines for periodic
!   boundary conditions
! psi0(i)=psi0(xd(i)): initial function evaluated at the quadrature 
!   point xd(i),  i=1,..,nq 
! psi0r(1:nmax),psi0i(1:nmax): real and imaginary parts of the 
!   overlaps <psi0|B_j^kp>, con B_j^kp(x) the basis functions
!   for j=1,...,nmax
! cpsi0r(1:nmax),cpsi0i(1:nmax): real and imaginary parts of the 
!   expansion coefficients of psi0 in the bspline basis
! d: overlap matrix is also constructed here
! s. Notes on bsplines, pag. 4 and latex report
!*******************************************************************

SUBROUTINE initial_wavepacketpbc

  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module
  USE detax
  USE varindexb
  USE detbx 
  use detbxcx
  USE varindex
  use mtx_pbc

  IMPLICIT NONE
  
  integer :: indexb
  INTEGER(ITT) :: i,j,iq,s,i1,ib,ipvt(kp-1),info!,Nran
  REAL(RT) :: aux,suma,a0,maxr,maxi,vsol(kp-1,2),wmu,step,ispline,aux1,dump
  !REAL(RT), DIMENSION(:,:), ALLOCATABLE :: w1,w2,wauxr,wauxi,maux0,maux1
  !REAL(RT),ALLOCATABLE,DIMENSION(:) :: vdeltar,wdeltar,vdeltai,wdeltai
  REAL(RT),ALLOCATABLE :: v1r(:),v1i(:),xran(:),yran(:),bran(:),cran(:),dran(:)
  COMPLEX(CT) :: sumac,suma2
  real(RT) :: aaa(4,2),bbb(3,2),ccc(5,3)
  integer, parameter :: out_unit=20
  real(RT), allocatable :: repsi(:),impsi(:),xs(:)!,dump(:)!read initial condition variables
  real(RT), allocatable :: xo(:),xaux(:),reaux(:),imaux(:)
  integer :: m,mtot,morder,k,k2!,ktot
  REAL(RT) :: aux2,dy,dymax,txmax,txmin,tdx
  ALLOCATE(psi0(nq),psi0r(nmax),psi0i(nmax),cpsi0r(nmax),cpsi0i(nmax),v1r(nmax),v1i(nmax))
  ALLOCATE(xsr(nmax),xsi(nmax))
  psi0r=0.0_RT;psi0i=0.0_RT;v1r=0.0_RT;v1i=0.0_RT
  xsr=zero;xsi=zero

  ntot=nmax-kp+1  !for periodic boundary conditions
  nlarg=kp
  nlard=kp
  wmu=sqrt(abs(g))
  wmu=2.0d0*pi_RT
  wmu=wmu/(xmax-xmin)
  WRITE(nfileoutput,'(a)')
  WRITE(nfileoutput,'(a)')         '******  SUB. initial_wavepacketpbc  ******'
  WRITE(nfileoutput,'(a24,1x,i6)') 'Banded matrices    ntot=',ntot
  WRITE(nfileoutput,'(a24,1x,i6)') 'Banded matrices   nlarg=',nlarg
  WRITE(nfileoutput,'(a24,1x,i6)') 'Banded matrices   nlard=',nlard
  
  !### create random phase function 

  GO TO (91, 62,66), ran !select between  hardcode cosine, random phase or textfile for initial condition

 
  62 CONTINUE
  !**************************************************************** 
  !!random phase with constant density!!
  !takes Nran random points and interpolate them on [-Pi,Pi]
  !then uses them as initial phase 
  !****************************************************************
  step=(xmax-xmin)/(Nran-1)
  allocate (xran(Nran),yran(Nran),bran(Nran), cran(Nran), dran(Nran))
  xran=zero
  yran=zero
  !making sample of randompoints
  xran(1)=xmin
  yran(1)= PI_RT*(2*rand()-1.0d0)
  Do i=1,Nran-1
     xran(i+1)=xran(i)+step !((xmax-xmin)*rand()+xmin)!xran(i)+step
     yran(i+1)=PI_RT*(2*rand()-1.0d0)
  End do



  yran(Nran)=yran(1) ! periodic boundary condition
  !calculates coefficients for interpolation
  call splineperiodic (xran, yran, bran, cran, dran,Nran)
  write(*,*) 'Constant density and random phase initial condition selected' 

  !open (unit=out_unit,file="initphase.txt",action="write",status="replace")
  
  !interpolates and finds the values for the initial wave function
  DO i=1,nq

    aux=ispline(xd(i),xran, yran, bran, cran, dran,Nran)

    psi0(i)=CMPLX(COS(aux),SIN(aux),CT)
    write(101,*) xd(i),dble(psi0(i)),dble(-ione*psi0(i))

   ! write(out_unit,*) xd(i),aux
  END DO
  close (out_unit)
  !stop
  !write(*,*) 'done' 
  !93 CONTINUE
  GO TO 108
  !stop
  !****************************************************************
  !!textfile!!
  ! read input text asuming it is on a regular grid
  !txmin and txmax are the edges of the grid and ktot is the number of points
  !here we asume periodic boundary conditions 
  !
  ! The real and imaginary part of the function in the grid are saved on repsi and impsi respectively
  ! and the spacial grid is saved on xo.
!****************************************************************
  66 CONTINUE
  Write(*,*) 'after open 66'
  txmin=xmin !-2437.2399999999998_RT     
  txmax=xmax!2437.2399999999998_RT     
  
  tdx=(txmax-txmin)/ktot
  allocate(repsi(ktot+1))
  allocate(impsi(ktot+1))
  allocate(xo(ktot+1))
  open (unit=90, file=initialtext, status='old', action='read')
 

  write(*,*) 'Reading initial condition' 
  xo(1)=txmin
  do i=1,ktot
     read(90,*)repsi(i),impsi(i)
     xo(i+1)=xo(i)+tdx
  enddo

  repsi(ktot+1)=repsi(1) !assuming periodic where re(xini)=re(xmax) and re(xini) is repsi(1)
  impsi(ktot+1)=impsi(1)

  close(90)

  morder=6
  allocate(xaux(morder)) !auxiliary list needed for the interpolation
  allocate(reaux(morder))
  allocate(imaux(morder))
  dymax=zero
  k=0
  do i=1,nquad*l0
     aux1=zero!xs(i)
     aux2=zero!xs(i)
     !locates the indices k, k2 of the interpolation region
     call locate(xo,ktot+1,xd(i+nq0),j)
     k = min(max(j-(morder-1)/2,1),ktot+1+1-morder)
     k2=min(k+morder,ktot+1)
     !list with the interpolation points 
     xaux=xo(k:k2)
     reaux=repsi(k:k2)
     imaux=impsi(k:k2)
     !interpolates the real and im part at xd(i+nq0) and save it at aux1 and aux2
     call polint(xaux,reaux,morder,xd(i+nq0),aux1,dy)
     call polint(xaux,imaux,morder,xd(i+nq0),aux2,dy)
      !call ratint(xaux,reaux,morder,xs(i),aux1,dy)
      !call ratint(xaux,imaux,morder,xs(i),aux2,dy)

     if(Abs(dy).gt.dymax)then
          dymax=abs(dy)
     end if
     !write(80,*) xd(i+nq0), aux1,aux2
     psi0(i+nq0)=CMPLX(aux1,aux2,CT)

     write(101,*) xd(i+nq0),dble(psi0(i+nq0)),dble(-ione*psi0(i+nq0))
  enddo
   !write(*,*) xs(1),xs(2)
  deallocate(repsi,impsi,reaux,imaux,xaux,xo)
  Write(*,*) 'after ending cycle 66'
  

  write(*,*) 'textfile initial condition selected' 
  !stop
!****************************************************************
  GO TO 108


  91 CONTINUE
!****************************************************************
  !###calculate the initial condition 
  DO i=1,nq
    !aux=-(ONE/100_RT)  +1.0_RT/SQRT(sigma*SQRT(pi_RT))*EXP(-(xd(i)-x0)**2/(2.0_RT*sigma**2))
    !aux=1.0_RT/(COSH(wmu*xd(i)))
    !aux=SQRT(100.0_RT/(sigma*SQRT(pi_RT))*EXP(-(xd(i)-x0)**2/(sigma**2)))
    !aux=SQRT(1.0d0-(1.0d0/COSH(xd(i)))*(1.0d0/COSH(xd(i)))*(1.0d0/COSH(xd(i)))&!  ) -Sech[x]^3 + Sech[x] Tanh[x]^2
    !+(1.0d0/COSH(xd(i)))*TANH(xd(i))*TANH(xd(i)) +1.5708d0)
    aux=SQRT(1.0d0+(0.0999999926465715d0)*(COS(wmu*xd(i))) + 0.3834951875677391d-4*SIN(wmu*xd(i))  ) 
    psi0(i)=CMPLX(aux*COS(p0*xd(i)),aux*SIN(p0*xd(i)),CT)
    !psi0(i)=CMPLX(aux,0.0d0,CT)
    write(101,*) xd(i),dble(psi0(i)),dble(-ione*psi0(i))
    !write(107,*) xd(i),abs(psi0(i))**2
  END DO
  
  write(*,*) 'gaussian initial condition selected' 
 !****************************************************************
  108 CONTINUE !###
  write(*,*) 'computing initial wave function'
  !write(*,*) 'end'
  !stop
  !norm of the initial wave function
  suma=0.0_RT
  DO j=1,l
     DO iq=1,nquad
        suma=suma+wd(iq)*ABS(psi0((j-1)*nquad+iq))**2*dx/2.0_RT
     END DO
  END DO
  WRITE(*,*)'norm of the initial wave function must be 1:',suma
 
  
  a0=dx/2.0_RT
  DO j=1,nmax
     sumac=czero
     DO s=0,kp-1
        DO iq=1,nquad
           i1=nq0+(j+s-1)*nquad+iq
           i1=index_xd(i1)
           sumac=sumac+wd(iq)*psi0(i1)*bba(indexb(j,i1))*a0
        ENDDO
     ENDDO
     psi0r(j)=REAL(sumac,RT)
     psi0i(j)=REAL(-ione*sumac,RT)
     write(102,*)psi0r(j),psi0i(j)
!     write(*,*)psi0r(j),psi0i(j)
  end DO
  
  




  !coefficients for the expansion in bsplines (s. Notes on bsplines (4.1) and pdf report)
  !construction of the overlap matrix B: stored in d, bp and b0
  ALLOCATE(d(nlarg,ntot)) 
  write(*,*)nlarg,nlard,ntot
  d=zero
  CALL over_pbc     !Array overlap(1:kp) containing the upper diagonal elements
  !write(*,*)overlap
  write(*,*)'ntot,nlarg,nmax',ntot,nlarg,nmax
  DO i=0,kp-1
     DO j=i+1,ntot
        ib=nlard-i                  !position for the array d
        d(ib,j)=d(ib,j)+overlap(i+1)
     END DO
  END DO
  ALLOCATE(b0(kp-1,kp-1))
  ALLOCATE(b1(kp-1,kp-1),b2(kp-1,kp-1))
  b0=zero;
  b1=zero;b2=zero
  !Matrix b0
  do i=1,kp-1
     do j=i,kp-1
        b0(i,j)=overlap(j-i+1)
        b0(j,i)=b0(i,j)
     end do
  end do
  !Matrix b1
  do i=0,kp-2
     do j=1,kp-1-i
        b1(j,j+i)=overlap(kp-i)
     end do
     do j=1,i+1
        b2(i+1,j)=overlap(kp-1+j-i)
     enddo
  end do
  WRITE(*,*)'b1 in initial_wavepacketpbc'
  do i=1,kp-1
     write(*,*)b1(i,:)
  end do
  !matrices ar and ai containing B for LLT decomposition and backwards substitution
  ALLOCATE(ar(nlarg,ntot),ai(nlarg,ntot))

  call solver(psi0r,psi0i,d,b0,b1,b2) !solves Ax=y with A=B 
  !where d has the diagonal part and the inputs b0,b1,b2 are the ofdiagonal parts of B
  !the result x is saven in psi0r
  cpsi0r=psi0r
  cpsi0i=psi0i
  deallocate(d,ar,ai)
END SUBROUTINE initial_wavepacketpbc


!****************************************************************
!subroutines needed for random phase initial condition
!****************************************************************

  subroutine splineperiodic (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
! 
!----------------------------------------------------------------------
!  input:
!  x = the arrays of data abscissas (in strictly increasing order) dim(n) (real)
!  y = the arrays of data ordinates dim(n) (real)
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients dim(n) (real)
!======================================================================
USE nrtype
implicit none
INTEGER(ITT) :: n, ipvt(n-1)
REAL(RT) :: x(n), y(n), b(n), c(n), d(n),a(n+1)
INTEGER(ITT) :: i, j, gap,info
REAL(RT) :: h(n),xaux(n-1),yaux(n-1)
REAL(RT) :: m(n-1,n-1)
gap = n-1
if(n.gt.1000000) then
	write(*,*) 'the number of points to interpolate is to large' , n
	stop
end if
! preparing matrix
do i = 1, n
  a(i) = y(i) 
end do
a(n+1)=y(2) !<-periodic part
do i = 1, n-1
  h(i) = x(i+1)-x(i) 
end do
h(n)=h(2)   !<_periodic part
m=zero
xaux=zero
yaux=zero
!filling matrix
do i = 1, n-1
  if(i == 1) then
    m(i,1) = 2.0d0*(h(1)+h(2))
    m(i,2) = h(2)
    else
    m(i, i - 1) = h(i)
    m(i, i) = 2.0d0 * (h(i) + h(i + 1))
    if (i < n - 1) then
       m(i, i + 1) = h(i + 1)
    end if
  end if
  if ((h(i) .ne. 0.0d0) .AND. (h(i + 1) .ne. 0.0d0)) then
     yaux(i) = ((a(i + 2) - a(i + 1)) / h(i + 1) - (a(i + 1) - a(i)) / h(i)) * 3.0d0
     else
     yaux(i) = 0.0d0
  end if
end do
m(1,n-1)=h(1)
m(n-1,1)=h(1)
!solving matrix 
call  dgesv (n-1, 1, m, n-1, ipvt, yaux, n-1, info )
!saving coefficients
do i = 2, n
  c(i) = yaux(i-1) 
end do

c(1)= c(n) !boundary for periodicity

do i = 1, n
  if (h(i) .ne. 0.0d0) then
     d(i)=1.0d0 / 3.0d0 / h(i) * (c(i + 1) - c(i))
     b(i) = (1.0d0 / h(i))* (a(i + 1) - a(i)) -( h(i) / 3.0d0 )* (c(i + 1) + 2 * c(i));
  end if 
end do


!write(*,*) x
!write(*,*) 'finished'
!stop

end subroutine splineperiodic


  function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated (real)
! x, y    = the arrays of given data points dim(n) (real)
! b, c, d = arrays of spline coefficients computed by spline dim(n) (real)
! n       = the number of data points (integer)
! output:
! ispline = interpolated value at point u (real)
!=======================================================================
!USE nrtype
implicit none
double precision ispline
integer n
double precision  u, x(n), y(n), b(n), c(n), d(n)
integer i, j, k
double precision dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
  ispline = y(1)
  return
end if
if(u >= x(n)) then
  ispline = y(n)
  return
end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
    j=k
    else
    i=k
   end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline


!****************************************************************
!subroutines needed for interpolation
!****************************************************************

SUBROUTINE polint(xa,ya,n,x,y,dy)
!======================================================================
!
!Polynomial Interpolation and Extrapolation
!
!Given arrays xa and ya , each of length n , and given a value x , this routine returns a
!value y , and an error estimate dy . If P (x) is the polynomial of degree N − 1 such that
!P ( xa i ) = ya i , i = 1, . . . , n , then the returned value y = P ( x ).
!
! This routine was adapted from from "Numerical recipes in Fortran 77 : 
!the art of scientific computing / William H. Press vol 1" section $3.1
!----------------------------------------------------------------------
! input..
! xa,ya     = the arrays of given data points of the spacial grid and the function dim(n) (real)
! n         = the number of data points (int)
! x         = point to be evaluated (real)
! output:
! y 	    = result of the interpolation at x (real)
! dy	    = aproximate error (real)
!=======================================================================
use nrtype
INTEGER :: n,NMAX
DOUBLE PRECISION  ::  dy,x,y,xa(n),ya(n)
PARAMETER (NMAX=10)
!Largest anticipated value of n.
!Given arrays xa and ya , each of length n , and given a value x , this routine returns a
!value y , and an error estimate dy . If P (x) is the polynomial of degree N − 1 such that
!P ( xa i ) = ya i , i = 1, . . . , n , then the returned value y = P ( x ).
INTEGER :: i,m,ns
DOUBLE PRECISION  ::  den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
ns=1
dif=abs(x-xa(1))

do i=1,n
	dift=abs(x-xa(i))
	if (dift.lt.dif) then
		ns=i
		dif=dift
	end if
	c(i)=ya(i)
	d(i)=ya(i)
end do 
y=ya(ns)
ns=ns-1
do m=1,n-1
	do i=1,n-m
		ho=xa(i)-x
		hp=xa(i+m)-x
		w=c(i+1)-d(i)
		den=ho-hp
		if(den.eq.0) write(*,*)'failure in point'
		den=w/den
		d(i)=hp*den
		c(i)=ho*den
	end do 
	if (2*ns.lt.n-m)then
		dy=c(ns+1)
	else
		dy=d(ns)
		ns=ns-1
	endif
	y=y+dy
end do 
return
end SUBROUTINE polint

SUBROUTINE locate(xx,n,x,j)
!======================================================================
!locate and hunt return an index j such that your
!desired value lies between table entries xx(j) and xx(j+1), where xx(1:n) is the
!full length of the table.
! xx(j) <x< xx(j+1)
! This routine was adapted from from "Numerical recipes in Fortran 77 : 
!the art of scientific computing / William H. Press vol 1" section $3.4
!----------------------------------------------------------------------
! input..
! xx      = the arrays of given data points of the spacial grid dim(n) (real)
! x       = desired value (real)
! n       = the number of data points (integer)
! output:
! j 	  = index j (integer)
!=======================================================================
use nrtype
INTEGER :: j,n
DOUBLE PRECISION  ::  x,xx(n)
!Given an array xx(1:n) , and given a value x , returns a value j such that x is between
!xx(j) and xx(j+1) . xx(1:n) must be monotonic, either increasing or decreasing. j=0
!or j=n is returned to indicate that x is out of range.
INTEGER :: jl,jm,ju
jl=0
ju=n+1
10 if(ju-jl.gt.1)then
	jm=(ju+jl)/2
	if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
		jl=jm
	else
		ju=jm
	end if
   goto 10
   end if
if(x.eq.xx(1))then
	j=1
else if(x.eq.xx(n))then
	j=n-1
else
	j=jl
end if
return
end SUBROUTINE locate
!****************************************************************
!subroutines needed for initial condition and psipotential
!****************************************************************
!****************************************************************
!# Descending sort of sets of 3 numbers 
subroutine sort_ijk(v)
!======================================================================
!locate and hunt return an index j such that your
!desired value lies between table entries xx(j) and xx(j+1), where xx(1:n) is the
!full length of the table.
! xx(j) <x< xx(j+1)
! This routine was adapted from from "Numerical recipes in Fortran 77 : 
!the art of scientific computing / William H. Press vol 1" section $3.4
!----------------------------------------------------------------------
! input..
! v      = the arrays of points to be sorted in descending order dim(3)
! output:
! v 	 = sorted array dim(3)
!=======================================================================
  USE nrtype
  use bspline_module
  use prob_param
  
  implicit none

  integer(ITT) :: v(3),i,j,maximo,maximo0,vaux(3),i0,a
  !sort v(i)
  do j=2, 3
    a=v(j)
    do i=j-1,1,-1
      if (v(i)<=a) goto 10
      v(i+1)=v(i)
    end do
	i=0
10  v(i+1)=a
  end do

 
  do j=1,3
     vaux(j)=-v(1)+v(j)+1
  end do
  v=vaux
 

end subroutine sort_ijk
!#
!****************************************************************
! Solver of Ax=y
!****************************************************************


SUBROUTINE solver(u0r,u0i,dd,dd0,dd1,dd2)
!======================================================================
! Solves the problem Ax=y for x  where A must be real
! u0r,u0i are the real and imaginary parts of y and dd,dd0,dd1,dd2 are the 
! nonzero matrix blocks of A
! the solution x is saved on u0r and u0i
! for more details of A see figure 3 of report and section 1.5
!----------------------------------------------------------------------
! input..
! dd     	   = diagonal block of A dim(nlarg,ntot) (real)
! dd0,dd1,dd2      = of diagonal blocks of A dim(kp-1,kp-1)(real)
! u0r,u0i          = real and imaginary part of y dim(nmax)(real)
! output:
! u0r,u0i	   = real and imaginary parts of x dim(nmax)(real)
!=======================================================================
  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module
  USE detax
  USE varindexb
  USE detbx 
  use detbxcx
  USE varindex
  use mtx_pbc

  IMPLICIT NONE
  
  
  INTEGER(ITT) :: i,j,iq,s,i1,ib,ipvt(kp-1),info
  REAL(RT) :: vsol(kp-1,2),u0r(nmax),u0i(nmax), dd1(kp-1,kp-1), dd2(kp-1,kp-1), dd0(kp-1,kp-1),dd(nlarg,ntot)
  REAL(RT), DIMENSION(:,:), ALLOCATABLE :: w1,w2,wauxr,wauxi,maux0,maux1
  REAL(RT),ALLOCATABLE,DIMENSION(:) :: vdeltar,wdeltar,vdeltai,wdeltai
  REAL(RT),ALLOCATABLE :: v1r(:),v1i(:)
  !real(RT) :: aaa(4,2),bbb(3,2),ccc(5,3)

  ALLOCATE(v1r(nmax),v1i(nmax))
  v1r=0.0_RT;v1i=0.0_RT
  xsr=zero;xsi=zero

  IF(ALLOCATED(ntl))DEALLOCATE(ntl)!@!
  ALLOCATE(w1(kp-1,kp-1),w2(kp-1,kp-1),maux0(kp-1,kp-1))
  w1=zero
  w2=zero

  !matrices ar and ai containing B for LLT decomposition and backwards substitution
  
  ar=dd
  ai=zero
  IF(ALLOCATED(ntlb))DEALLOCATE(ntlb)
  ALLOCATE(ntlb(0:ntot))
  ntlb=0
  ntlb(0)=-1
  !L*LT DECOMPOSITION
  nnggrr=ntot
  
  CALL lddc(ntot,nlard)
  
  !backward substitution for the solution of the systems of equations (see pdf report):
  ! M_delta v=y_delta    y    M_delta w_j=m_j, j=1,...,kp-1
  ! M_delta is stored in ar,ai.  m_j are the column vectors of bp. y==psi0r,psi0i
  ALLOCATE(vdeltar(ntot),wdeltar(ntot),vdeltai(ntot),wdeltai(ntot))
  vdeltar=0.0_RT;wdeltar=0.0_RT;vdeltai=0.0_RT;wdeltai=0.0_RT;
  ALLOCATE(wauxr(ntot,kp-1),wauxi(ntot,kp-1))
  wauxr=zero
  wauxi=zero
  
  ! we start with  M_delta v=y_delta
  !b=xsr+I*xsi; y=psi0r+I*psi0i; v1r,v1i are auxiliar arrays
  CALL axdc(ntot,nlard,u0r,u0i,vdeltar,vdeltai,v1r,v1i)
  xsi=zero
  xsr=zero
  do j=1,kp-1
     xsr(1:kp-1)=dd1(:,j)
     xsr(ntot-kp+2:ntot)=dd2(:,j)
     CALL axdc(ntot,nlard,xsr,xsi,wauxr(:,j),wauxi(:,j),v1r,v1i)
  end do
  !d and bp are real, xsi is zero, therefore wauxi is zero
  !contsruction of the matrix M0-(M1^+)**T*W (pdf report)
  maux0=zero
  maux0=dd0
  w1=wauxr(1:kp-1,1:kp-1)
  w2=wauxr(ntot-kp+2:ntot,1:kp-1)
  call DGEMM('T','N',kp-1,kp-1,kp-1,mone,dd1,kp-1,w1,kp-1,one,maux0,kp-1)
  call DGEMM('T','N',kp-1,kp-1,kp-1,mone,dd2,kp-1,w2,kp-1,one,maux0,kp-1)
  !maux0 is overwritten: maux0=mone*(bp)**T*wauxr+one*maux0
  xsr(1:kp-1)=u0r(ntot+1:nmax)
  xsi(1:kp-1)=u0i(ntot+1:nmax)
  !matrix-vector product (M1^+)**T*vdelta
  call dgemv('T', kp-1, kp-1, mone, dd1,kp-1, vdeltar(1:kp-1), 1, one, xsr, 1)
  call dgemv('T', kp-1, kp-1, mone, dd2,kp-1, vdeltar(ntot-kp+2:ntot), 1, one, xsr, 1)
  
  call dgemv('T', kp-1, kp-1, mone, dd1,kp-1, vdeltai(1:kp-1), 1, one, xsi, 1)
  call dgemv('T', kp-1, kp-1, mone, dd2,kp-1, vdeltai(ntot-kp+2:ntot), 1, one, xsi, 1)
  !output: xs(1:kp-1)=xs-bp**T*vdelta
  !                  =x0-bp**T*vdelta
  !solving maux0*x0=xs; xs is now stored in vsol
  vsol(:,1)=xsr(1:kp-1)
  vsol(:,2)=xsi(1:kp-1)
  call  dgesv (kp-1, 2, maux0, kp-1, ipvt, vsol, kp-1, info )
  !write(*,*)'x',ipvt
  !vsolr and vsoli contain now the solutions of maux0*x0=xs
  !matrix vector product vdelta-wauxr*vsol; Note that wauxi=0
  call dgemv('N', ntot, kp-1, mone, wauxr,ntot, vsol(:,1), 1, one, vdeltar, 1)
  call dgemv('N', ntot, kp-1, mone, wauxr,ntot, vsol(:,2), 1, one, vdeltai, 1)
  u0r(1:ntot)=vdeltar
  u0i(1:ntot)=vdeltai
  u0r(ntot+1:ntot+kp-1)=vsol(:,1)
  u0i(ntot+1:ntot+kp-1)=vsol(:,2)
  

  ar=zero
  ai=zero
  IF(ALLOCATED(ntl))DEALLOCATE(ntl)
  DEALLOCATE(v1r,v1i)
  deallocate(wauxr,wauxi)
  deallocate(vdeltar,wdeltar,vdeltai,wdeltai,w1,w2,maux0)
  RETURN
END SUBROUTINE solver
!****************************************************************
