program FILEREADER
  use nrtype
  ! USE read_param
  USE prob_param
  USE bspline_module
  USE quad_module
!!
!this program takes the coeficients and give the wave function evaluated in a grid
!!
   real(RT), allocatable :: power(:),impsi(:),aux(:,:),dump(:),xo(:),poweraux(:),xaux(:)
   integer :: i,morder,k,k2,j  !n,m,i,mtot
  ! INTEGER :: hdferr,long
  !! INTEGER :: layout
  ! INTEGER(HID_T)  :: file, space, dset, dcpl ! Handles
   CHARACTER(LEN=37) :: filename1 ,filename2,auxstring
   !CHARACTER(LEN=3) , PARAMETER :: dataset  = "DS1"
   REAL(RT) :: timeunitary1,timeunitary2,timeunitary3,aux1,aux2,dxstep,dy,dymax

   READ(*,*)filename1
   READ(*,*)lmax   !number of points output grid
   READ(*,*)l0     !number of points in the grid
   READ(*,*)xmin   !start of the box
   READ(*,*)xmax   !end of the box (usually rmin=-rmax)


  

  dxstep=(xmax-xmin)/(1d0*lmax)
  filename2='GRID'//TRIM(filename1)
  l0=l0  
  
  L=100_RT
  P=L/2
  dxstep=(PI_RT)/P
  lmax=int((xmax-xmin)/dxstep)
  write(auxstring,*)L
  filename2=TRIM(auxstring)//filename2
  WRITE(*,*) 'lmax', lmax
!!!!!!!!
  allocate(power(l0+1))
  !allocate(impsi(l0+1))
  allocate(xo(l0+1))
  open (unit=90, file=filename1, status='old', action='read')
 

  write(*,*) 'Reading initial condition' 
  !#xo(1)=txmin
  read(90,*)
  do i=1,l0+1
     read(90,*) k2,xo(i),power(i)
     !#xo(i+1)=xo(i)+tdx
  enddo
  write(*,*) 'Readed',xo(l0),power(l0)
  !repsi(ktot+1)=repsi(1) !assuming periodic where re(xini)=re(xmax) and re(xini) is repsi(1)
  !impsi(ktot+1)=impsi(1)
  


  !xmin=xo(1)
  !dxstep=(xmax-xmin)/(1d0*lmax)

  close(90)
  open (unit=101, file=filename2, action='write')
  morder=6
  allocate(xaux(morder)) !auxiliary list needed for the interpolation , status='new'
  allocate(poweraux(morder))
  !allocate(imaux(morder))
  allocate(xd(lmax+1))
  dymax=zero
  k=0


  xd(1)=xmin
  do i=1,lmax
     xd(i+1)=xd(i)+dxstep
  enddo



  write(*,*) 'before'
  do i=1,lmax
     aux1=zero!xs(i)
     aux2=zero!xs(i)
     !xd(i+1)=xd(i)+dxstep
     !locates the indices k, k2 of the interpolation region
    ! write(*,*) 'before  locate',i
     call locate(xo,l0+1,xd(i),j)
     !write(*,*) 'after  locate',i
     k = min(max(j-(morder-1)/2,1),l0+1+1-morder)
     k2=min(k+morder,l0+1)
     !list with the interpolation points 
     xaux=xo(k:k2)
     poweraux=power(k:k2)
     !imaux=impsi(k:k2)
    ! write(*,*) 'before  polint',i
     !interpolates the real and im part at xd(i+nq0) and save it at aux1 and aux2
     call polint(xaux,poweraux,morder,xd(i),aux1,dy)
     !call polint(xaux,imaux,morder,xd(i+nq0),aux2,dy)
      !call ratint(xaux,reaux,morder,xs(i),aux1,dy)
      !call ratint(xaux,imaux,morder,xs(i),aux2,dy)
     !write(*,*) 'after  polint',i
     if(Abs(dy).gt.dymax)then
          dymax=abs(dy)
     end if
     !write(80,*) xd(i+nq0), aux1,aux2
     !psi0(i+nq0)=CMPLX(aux1,aux2,CT)

     write(101,*) xd(i),dble(aux1)
  enddo
   !write(*,*) xs(1),xs(2)
  write(*,*) 'after'
  deallocate(power,poweraux,xaux,xo)
  Write(*,*) 'after ending cycle 66'
  close(101)

  write(*,*) 'textfile initial condition selected' 






end program


!**************************************************************************
SUBROUTINE knots(v,ll)

  USE nrtype
  USE bspline_module
  USE prob_param
  
  IMPLICIT NONE
  !----------------------------------------------------------------------------
  ! Gives the knots {v}, i=1,...,n+kp with n=l+kp-1 in the interval [xmin,xmax]
  !  with l subdivisions from the breakpoints calculated with bps
  !  Endpoints with multiplicity kp
  !  Calculates the breakpoints sequence between [0,xmax] and stores
  !  in temporary vector xt. Then saves the total breakpoints in x, which will 
  !  run from (-xmax) to (xmax). If change to [0,xmax] only, change to store
  !  xt=x and lmidp=ll+1.
  !  INPUT
  !    ll:  number of breakpoints
  !  OUTPUT
  !    v: knots
  !----------------------------------------------------------------------------
  INTEGER(ITT) :: lmidp,m,ll,j
  REAL(RT) :: v(1:ll+2*kp-1)
  REAL(RT) :: x(1:ll+1)

  xmin0=xmin-(kp-1)*dx

  !WRITE(*,*)xmin0,xmin,(kp-1)*dx,dx,ll
  
  DO j=1,ll+1
     x(j)=(j-1)*dx+xmin0 !breakpoints
  ENDDO

  DO m=1,kp
     v(m)=x(1)
  ENDDO
  DO m=1,nb-kp
     v(m+kp)=x(m+1)
  ENDDO
  DO m=1,kp
     v(nb+m)=x(ll+1)
  ENDDO

END SUBROUTINE knots
!**************************************************************************
SUBROUTINE basis_construction

  USE nrtype
  USE bspline_module
  USE quad_module
  USE prob_param

  IMPLICIT NONE

  INTEGER(ITT) :: i,j,i0,maxi,cont
  REAL(RT),ALLOCATABLE :: bb_temp(:,:),dbb_temp(:,:)
  CHARACTER(25) :: auxfile


  WRITE(*,*)'matrix size nmax',nmax
  !Bspline j=2*kp-1 and its derivative is stored in bba and dbba
  !This function is non-zero on the interval (xmin,xmin+k*dx) which 
  ! contains nqk=nquad*kp quadrature points
  ALLOCATE(bba(0:nqk),dbba(0:nqk))
  DO i=1,nqk
     bba(i)=bb(2*kp-1,i)
     dbba(i)=dbb(2*kp-1,i)
  END DO
  bba(0)=zero;dbba(0)=zero
  DEALLOCATE(bb,dbb)
  !The rest of the basis functions are cyclic translations of this function
  !The information is stored in the index array indexb(:,:)
  !ALLOCATE(indexb(nmax,nq))
  !indexb=0
  !DO j=1,nb-3*kp+3
  !   i0=nq0+(j-1)*nquad
  !   cont=0
  !   DO i=1,nqk
  !      IF(i0+i.GT.nq-(kp-1)*nquad)THEN
  !         cont=cont+1
  !         indexb(j,nq0+cont)=i
  !      ELSE
  !         indexb(j,i0+i)=i
  !      END IF
        !IF(j==nb-2*kp+2)WRITE(106,*)i,indexb(j,i),'nb-2*kp+2=',nb-2*kp+2,i0,nq-(kp-1)*nquad
  !   END DO
  !END DO

  !array index_xd(1:nq)
  write(*,*) 'making index_xd'
  ALLOCATE(index_xd(nq))
  index_xd=0
  DO i=1,nq0+nquad*nmax
     index_xd(i)=i
  END DO
  DO i=nq0+nquad*nmax+1,nq0+nquad*(nmax+kp-1)
     index_xd(i)=i-nquad*nmax
     !write(*,*)i,index_xd(i),xd(i),xd(index_xd(i))
  END DO
  write(*,*) 'finished index_xd'
!!$  !(****  output: bspline basis for periodic boundary conditions
!!$  DO j=1,nmax
!!$     WRITE(auxfile,200)'basis_',j,'.dat'
!!$     !WRITE(auxfile,201)'bspline',j,'.table' !pgf format
!!$     WRITE(*,*)auxfile
!!$     OPEN(4,file=auxfile)
!!$200  FORMAT(a6,i3.3,a4)
!!$201  FORMAT(a7,i3.3,a6)  !pgf-format
!!$202  FORMAT(2f24.16,a3)
!!$     DO i=1,nq
!!$        WRITE(4,*) xd(i),bba(indexb(j,i)),dx 
!!$     END DO
!!$     !stop
!!$     CLOSE(4)
!!$  END DO
!!$  !****)

END SUBROUTINE basis_construction

!**************************************************************************

SUBROUTINE tiempo(timereal)

  IMPLICIT NONE

  CHARACTER(10) :: date,time,charx
  INTEGER ::  hh,mm,ss,sss,cc,yy,month,dd,timesss
  DOUBLE PRECISION  :: timereal
    
  CALL DATE_AND_TIME(date,time)
  READ(time,700)hh,mm,ss,charx,sss
  READ(date,701)cc,yy,month,dd
  timesss=sss+ss*1000+60000*mm+hh*3600000+dd*86400000
  timereal=DBLE(timesss)/1000
  !write(*,*) timereal
700 FORMAT(i2.2,i2.2,i2.2,a1,i3.3)                                      !time
701 FORMAT(i2.2,i2.2,i2.2,i2.2,i2.2)                                    !date

END SUBROUTINE tiempo
!**************************************************************************
SUBROUTINE indexebb(j,ip,iout)

  USE nrtype
  USE bspline_module
  USE quad_module
  USE prob_param

  IMPLICIT NONE

  INTEGER(ITT) :: i,j,i0,a,b,iout,ip
  !REAL(RT),ALLOCATABLE :: bb_temp(:,:),dbb_temp(:,:)
  CHARACTER(25) :: auxfile

  i0=nq0+(j-1)*nquad
  If (ip.lt.i0)then
     iout=ip+nquad*nmax
  else  
     iout=ip
  end if 

  IF(iout-i0.GT.nqk)THEN
     
     iout=0
  ELSE
     iout=iout-i0
  END IF

END SUBROUTINE indexebb
!**************************************************************************

SUBROUTINE store_distr(cont,infilename)

  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module
  USE filenames

  IMPLICIT NONE

  INTEGER(ITT) :: cont,i,j,iq,auxint
  COMPLEX(CT) :: sumac
  CHARACTER(100) :: filename
  CHARACTER(37) :: infilename
  

  WRITE(filename,100)cont,infilename,'.dat'
100 FORMAT(i6.6,a7,a4)
 ! IF(cont==-1)filename=filedistr    !storing the final distribution
  OPEN(14,file=filename)
  print_all='y'
  IF(print_all=='y') THEN
     DO i=nq0+1,nq0+nmax*nquad
        sumac=czero
        DO j=1,nmax
           call indexebb(j,i,auxint)

           sumac=sumac+CMPLX(cpsi0r(j),cpsi0i(j),CT)*bba(auxint)
        END DO
        WRITE(14,*) xd(i),ABS(sumac)**2,REAL(sumac,RT), REAL(-ione*sumac,RT)
     END DO
  ELSE IF (print_all.EQ.'n'.AND.cont.NE.-1) THEN
     !in this case the output file contains information of the wave function psi(x)
     !at the first quadrature point of each breakpoint intervals
     DO i=kp,l-kp+1  !l=number of breakpoints
        sumac=czero
        iq=(i-1)*nquad+1
        DO j=1,nmax
           call indexebb(j,iq,auxint)
           sumac=sumac+CMPLX(cpsi0r(j),cpsi0i(j),CT)*bba(auxint)
!           sumac=sumac+CMPLX(cpsi0r(j-1),cpsi0i(j-1),CT)*bb(j,(i-1)*nquad+1)
        END DO
        WRITE(14,*) xd(iq),ABS(sumac)**2,REAL(sumac,RT), REAL(-ione*sumac,RT)
     END DO
  ELSE IF (print_all.EQ.'n'.AND.cont.EQ.-1) THEN
     !the final output file 'd_s...' contains the full wave function psi(x)
     DO i=nq0+1,nq0+nmax*nquad
        sumac=czero
        DO j=1,nmax
           call indexebb(j,i,auxint)
           sumac=sumac+CMPLX(cpsi0r(j),cpsi0i(j),CT)*bba(auxint)
        END DO
        WRITE(14,*) xd(i),ABS(sumac)**2,REAL(sumac,RT), REAL(-ione*sumac,RT)
     END DO
  END IF
  CLOSE(14)
  write(*,*) 'saving done'
END SUBROUTINE store_distr

!****************************************************************
! Quadrature points and weights at these points are generated.
! Gauss-Legendre quadratre is used:       
! subroutine gaussq(kind,....)    with kind=1
! nquad: number of quadrature points
! xd(1:nquad): quadrature points: zeros of legendre polynomials
!              in (-1,1)
! wd(1:nquad): weights at the quadrature points
SUBROUTINE evaluate_quad

  USE nrtype
  USE quad_module
  USE prob_param

  IMPLICIT NONE

  INTEGER(ITT) :: i
  REAL(RT):: endpts(2)
  REAL(RT), ALLOCATABLE :: work1(:),xd1(:),wd1(:)
  WRITE(*,*)'in evaluate_quad nquad=',nquad

  ALLOCATE(xd1(1000),wd1(1000),work1(1000))
  CALL gaussq(1,nquad,0.0d0,0.0d0,0,endpts,work1,xd1,wd1)
  DEALLOCATE(work1)
  IF(ALLOCATED(xd)) DEALLOCATE(xd)
  IF(ALLOCATED(wd)) DEALLOCATE(wd)
  ALLOCATE(xd(nquad),wd(nquad))
  xd(1:nquad)=xd1(1:nquad)
  wd(1:nquad)=wd1(1:nquad)	
  DEALLOCATE(xd1,wd1)

!   WRITE(102,*)'quadrature points in (-1,1) and weights: i,xd(i),wd(i)'
!   DO i=1,nquad
!      WRITE(102,*)i,xd(i),wd(i)
!   END DO

END SUBROUTINE evaluate_quad


!****************************************************************
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
