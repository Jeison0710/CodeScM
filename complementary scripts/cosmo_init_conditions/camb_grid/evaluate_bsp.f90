! ================================================================
! SUBROUTINE evaluate_bsp
! generation of the full grid and bsplines
!   l=nb-kp+1:	number of break points (over (xmin,xmax)
!   nquad:  	number of quadrature points in (-1,1)
!   nq=l*nquad: total number of quadrature points
!   xd(1:nq): 	quadrature points
!   bb(j,i):  	j-th bspline of order kp B_j^kp(xd(i))
!   dbb(j,i):	derivative of the j-th bspline of order kp:
!		d/dx(B_j^kp(x))_{x=xd(i)}
! ================================================================
SUBROUTINE evaluate_bsp

  USE nrtype
  USE bspline_module
  USE quad_module
  USE prob_param

  IMPLICIT NONE

  INTEGER(ITT) :: i,j
  INTEGER(ITT) :: left,mflag,leftmk,leftmk1,nderiv
  REAL(RT) :: xb,xa,xbpa,xbma
  REAL(RT),ALLOCATABLE :: fn(:,:),a(:,:)
  REAL(RT),ALLOCATABLE :: tmpx(:)
  CHARACTER(25) :: auxfile

  !First calculates the total grid in x for the quadrature and stores in tmpx
  ALLOCATE(tmpx(nq))
  DO i=1,l
     xb = tx(kp+i)
     xa = tx(kp+i-1)
     xbpa=(xb+xa)*0.5d0
     xbma=(xb-xa)*0.5d0
     DO j=1,nquad
	tmpx((i-1)*nquad+j) = xbpa+xbma*xd(j)
     ENDDO
  ENDDO
  DEALLOCATE(xd)
  !Saves the new quadrature in xd
  ALLOCATE(xd(nq))
  xd=tmpx
  DEALLOCATE(tmpx)

  !Evaluate the Bsplines in the quadrature points
  IF(ALLOCATED(bb)) DEALLOCATE(bb)
  IF(ALLOCATED(dbb)) DEALLOCATE(dbb)
  ALLOCATE(bb(n+2,nq),dbb(n+2,nq))
  bb=0d0
  dbb=0d0
  nderiv=2
  ALLOCATE(fn(kp,nderiv),a(kp,kp))

  DO i=1,nq
     CALL interv(tx,l+2*kp-1,xd(i),left,mflag)
     leftmk=left-kp
     CALL bsplvd(tx,kp,xd(i),left,a,fn,nderiv)
     leftmk1=leftmk+1
     !     write(103,*)'i,leftmk1,left,tx(i)=',i,leftmk1,left,tx(i)
     DO j=leftmk1,left
        bb(j,i)=fn(j-leftmk,1)
        dbb(j,i)=fn(j-leftmk,2)
     END DO
  END DO

  !bsplines and derivatives at the boundaries
  !bb_boundary(j,1): j-th bspline at xmin 
  !bb_boundary(j,2): j-th bspline at xmax 
  !dbb_boundary(j,1): derivative of the j-th bspline at xmin 
  !dbb_boundary(j,2): derivative of the j-th bspline at xmax 
  ALLOCATE(bb_boundary(n+2,2),dbb_boundary(n+2,2))
  bb_boundary=0.0_RT;dbb_boundary=0.0_RT
  CALL interv(tx,l+2*kp-1,tx(1),left,mflag)
  leftmk=left-kp
  CALL bsplvd(tx,kp,tx(1),left,a,fn,nderiv)
  leftmk1=leftmk+1
  DO j=leftmk1,left
     !write(*,*)tx(1),xmin,fn(j-leftmk,1),fn(j-leftmk,2)
     bb_boundary(j,1)=fn(j-leftmk,1)
     dbb_boundary(j,1)=fn(j-leftmk,2)
  END DO
  CALL interv(tx,l+2*kp-1,tx(l+2*kp-1),left,mflag)
  leftmk=left-kp
  CALL bsplvd(tx,kp,tx(l+2*kp-1),left,a,fn,nderiv)
  leftmk1=leftmk+1
  DO j=leftmk1,left
     !write(*,*)tx(l+2*kp-1),xmax,fn(j-leftmk,1),fn(j-leftmk,2)
     bb_boundary(j,2)=fn(j-leftmk,1)
     dbb_boundary(j,2)=fn(j-leftmk,2)
  END DO
  !tonto
  !bb_boundary=0.0_RT;dbb_boundary=0.0_RT


  WRITE(*,*)'in evaluate_bsp n,n+2,nmax=',n,n+2,nmax
!!$  !output: all bsplines
!!$  DO j=1,n+2
!!$     !write(auxfile,200)'b-temp',j,'.dat'
!!$     WRITE(auxfile,201)'bspline',j,'.table' !pgf format
!!$     write(*,*)auxfile
!!$     OPEN(4,file=auxfile)
!!$200  FORMAT(a6,i3.3,a4)
!!$201  FORMAT(a7,i3.3,a6)  !pgf-format
!!$202  FORMAT(2f24.16,a3)
!!$     DO i=1,nq
!!$        !        write(4,*) xd(i),bb(j,i),dx 
!!$        WRITE(4,202) xd(i),bb(j,i),'  i'
!!$        !        write(4,*) xd(i)/dx,bb(j,i),dx
!!$        !         write(4,*) xd(i),dbb(j,i)
!!$     END DO
!!$     CLOSE(4)
!!$  END DO

  DEALLOCATE(fn,a)

  RETURN
END SUBROUTINE evaluate_bsp

!***********************************************************************
SUBROUTINE evaluate_bsp1(x,nn,bb1)

  USE nrtype
  USE bspline_module
  USE quad_module
  USE prob_param

  IMPLICIT NONE

  INTEGER(ITT) :: i,j,nn
  INTEGER(ITT) :: left,mflag,leftmk,leftmk1,nderiv
  REAL(RT) :: x,bb1(nn+2)
  REAL(RT),ALLOCATABLE :: fn(:,:),a(:,:)

  !Evaluate the Bsplines in x
  bb1=0d0
  nderiv=1
  ALLOCATE(fn(kp,nderiv),a(kp,kp))

  CALL interv(tx,l+2*kp-1,x,left,mflag)
  leftmk=left-kp
  CALL bsplvd(tx,kp,x,left,a,fn,nderiv)
  leftmk1=leftmk+1
  DO j=leftmk1,left
     bb1(j)=fn(j-leftmk,1)
  END DO
  DEALLOCATE(fn,a)

  RETURN
END SUBROUTINE evaluate_bsp1
