program FILEREADER
  use nrtype
  ! USE read_param
  USE prob_param
  USE bspline_module
  USE quad_module
!!
!this program takes the coeficients and give the wave function evaluated in a grid
!!
   real(RT), allocatable :: repsi(:),impsi(:),aux(:,:),xs(:),dump(:)
   integer :: i!n,m,i,mtot
  ! INTEGER :: hdferr
  !! INTEGER :: layout
  ! INTEGER(HID_T)  :: file, space, dset, dcpl ! Handles
   CHARACTER(LEN=37) :: filename1 ,filename2
   !CHARACTER(LEN=3) , PARAMETER :: dataset  = "DS1"
   REAL(RT) :: timeunitary1,timeunitary2,timeunitary3,aux1,aux2,dxstep

   READ(*,*)filename1
   READ(*,*)kp     !order of the bsplines
   READ(*,*)nb     !number of bsplines for the initial knot sequence (assumed an uniform knot sequence).
   READ(*,*)xmin   !start of the box
   READ(*,*)xmax   !end of the box (usually rmin=-rmax)
   READ(*,*)nquad  !number of points in each interval nmax
   l=nb-kp+1
   l0=l-2*(kp-1)
   IF (l0.LT.1) THEN
      WRITE(*,*)'in print_output: bad choice of parameters. l0=',l0
      STOP
   END IF
   dx=(xmax-xmin)/(1d0*l0)
   nq=l*nquad      !total number of quadrature points
   nmax=l0
   nq0=nquad*(kp-1)
   nqk=nquad*kp
   dxstep=dx/nquad
   !ktot=4194304
   !xmin=-2437.24_RT
   !xmax=2437.24_RT
   ALLOCATE(tx(1:l+2*kp-1))   !knots will be generated in the sub. knots
   CALL knots(tx,l)
   WRITE(100,*)'l=',l,'tx:'
   DO i=1,l+2*kp-1
      WRITE(100,*)i,tx(i)
   END DO
   allocate(xd(nq))
   xd(1)=tx(1)
   !creates the grid to be evaluated
   DO i=1,nq-1
      xd(i+1)=xd(i)+dxstep
   ENDDO
   !write(*,*) 'evaluating quad'
   !CALL evaluate_quad            !Gauss quadrature weights and points
   !write(*,*) 'evaluating bsppbc'
   CALL evaluate_bsppbc
   CALL basis_construction
   write(*,*) ''
   ALLOCATE(cpsi0r(l0),cpsi0i(l0))
    write(*,*) 'reading coeficients cpsi'
   open (unit=32,file=filename1,action="read")!,status="old")
   DO i=1,l0
     !write(*,*) 'before read 32'
     Read(32,*) cpsi0r(i),cpsi0i(i)
    ! write(*,*) 'before read 32'
    !! psi0(i+nq0)=CMPLX(aux,aux1,CT)

   !  write(101,*) xd(i+nq0),dble(psi0(i+nq0)),dble(-ione*psi0(i+nq0))
    !write(107,*) xd(i),abs(psi0(i))**2
   END DO
   close (32)
   write(*,*) 'storing coef'
   call store_distr(9999,filename1)
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
