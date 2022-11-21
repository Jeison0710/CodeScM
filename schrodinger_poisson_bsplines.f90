!===================================================================
!14.12.2018; Javier Madronero
! Code in devolopment
! Program to compute the evolution of a free Gaussian wave packet 
!   with the Gross-Pitaevskii equation. The space coordinate is 
!   discretized in a regular grid for the representation in b-splines. 
! Time propagation is achieved with the help of the Crank-Nicolson 
!   method with a predictor-corrector.
! The parametrized TDSE is
!   i*df/dt=(-1/2*d^2/dx^2+g|f|^2)f
! The initial Gaussian wave packet in these coordinates reads
!  psi(z)=1/sqrt(pi*sqrt(sigma))*exp(-(z-z0)^2/(2*sigma^2)+i*p0*z)
!===================================================================
!-----------------------------------------------------

PROGRAM schrodinger_poisson_bsplines 

  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module
  USE filenames
  
  IMPLICIT NONE

  INTEGER(ITT) :: i

  !****************** INPUT PARAMETERS *******************!
  ! All parameters are given in SI units unless stated otherwise !
  !****************** System parameters ******************!
  READ(*,*)       !#      --Data for the initial gaussian wave packet
  READ(*,*)sigma  !Gaussian dispersion             !
  READ(*,*)x0     !initial position of the Gaussian!
  READ(*,*)p0     !initial velocity
  READ(*,*)g      !Gross-Pitaevskii parameter
  READ(*,*)ran    !ran phase init
  READ(*,*)initialtext  !name text file containing initial condition
  READ(*,*)ktot
  READ(*,*)Nran      !Nran
  !*************** Numerical parameters ***************!
  READ(*,*)       !# !     --Data for the Bspline basis and time evolution
  READ(*,*)kp     !order of the bsplines
  READ(*,*)nb     !number of bsplines for the initial knot sequence (assumed an uniform knot sequence).
  READ(*,*)xmin   !start of the box
  READ(*,*)xmax   !end of the box (usually rmin=-rmax)
  READ(*,*)nquad  !number of quadrature points (Legendre-Gauss quadrature) in each interval
  READ(*,*)tfin   !final time in scaled units
  READ(*,*)dt0    !time step in scaled units
  READ(*,*)aini
  READ(*,*)aend   
  !******************* Other parameters ******************!
  READ(*,*)       !# !     --Data for the Bspline basis and time evolution
  READ(*,*)nshow  !Information (norm) is displayed every nshow steps
  READ(*,*)tprinti !from where start to save states on nshow steps
  READ(*,*)tprinte !where to end save states shown on nshow steps
  READ(*,*)print_all  !'y': output files contain psi(x_i) for all quadrature points
  !                   !'n': output files contain psi(x_i) for one quadrature point per breakpoint
  READ(*,*)
  READ(*,*)V0_switch  !# !V0_switch: 0/1 => switch off/on the potential
  READ(*,*)discoef  !save 0 to save dist or 1 to save coefficients
  V0=1.0_RT
  IF (V0_switch==0) V0=0.0_RT

  !**************** END INPUT PARAMETERS *****************!

  !all parameters are given in scaled units;
!   !##################      INPUT DATA      ###################
!   !##################     PROGRAM START     ###################
  CALL filescreate
  CALL print_output
  ALLOCATE(tx(1:l+2*kp-1))   !knots will be generated in the sub. knots
  CALL knots(tx,l)
  WRITE(100,*)'l=',l,'tx:'
  DO i=1,l+2*kp-1
     WRITE(100,*)i,tx(i)
  END DO
  CALL evaluate_quad            !Gauss quadrature weights and points
  CALL evaluate_bsppbc
  WRITE(*,'(a)')'Bspline functions are ready'
  CALL basis_construction
  WRITE(*,'(a)')'Bspline basis for periodic boundary conditions is ready'
  WRITE(*,*)'matrix size nmax',nmax
  CALL initial_wavepacketpbc
  WRITE(*,'(a)')'Starting initialization of the matrix elements arrays'
  CALL ini_mtxs
  WRITE(*,'(a)')'Matrices are ready'
  !***** time propagation *****
  CALL crank_nicolson_propagator_pbc

  WRITE(*,*)'end'
  CLOSE(nfileoutput)

END PROGRAM schrodinger_poisson_bsplines

!*******************************************************************
! Construction of the basis for periodic boundary conditions
! The bspline j=2k-1 and its derivative is stored in the array
!  bba(:) and dbba(:) 
!*******************************************************************
SUBROUTINE basis_construction

  USE nrtype
  USE bspline_module
  USE quad_module
  USE prob_param

  IMPLICIT NONE
  integer :: indexb
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
    ! i0=nq0+(j-1)*nquad
    ! cont=0
    ! DO i=1,nqk
    !    IF(i0+i.GT.nq-(kp-1)*nquad)THEN
     !      cont=cont+1
     !      indexb(j,nq0+cont)=i
     !   ELSE
     !      indexb(j,i0+i)=i
    !    END IF
        !IF(j==nb-2*kp+2)WRITE(106,*)i,indexb(j,i),'nb-2*kp+2=',nb-2*kp+2,i0,nq-(kp-1)*nquad
   !  END DO
  !END DO

  !array index_xd(1:nq)
  ALLOCATE(index_xd(nq))
  index_xd=0
  DO i=1,nq0+nquad*nmax
     index_xd(i)=i
  END DO
  DO i=nq0+nquad*nmax+1,nq0+nquad*(nmax+kp-1)
     index_xd(i)=i-nquad*nmax
     !write(*,*)i,index_xd(i),xd(i),xd(index_xd(i))
  END DO
  
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

!*******************************************************************
! Creation of the names of the output files
!*******************************************************************
SUBROUTINE filescreate
  USE prob_param
  USE filenames
  IMPLICIT NONE

  CHARACTER(3) :: sgnp,sgnx

  sgnp='+'
  IF(p0<0.0_RT)sgnp='-'
  sgnx='+'
  IF(xmin<0.0_RT)sgnx='-'
  WRITE(fileoutput,102)'o','_','s',sigma,'x',x0,'p',sgnp,ABS(p0),&
       'n',nb,'k',kp,'dt',dt0,'tf',tfin,'xm',sgnx,ABS(xmin),'xM',xmax&
       ,'g',g
  WRITE(filenorm,102)'n','_','s',sigma,'x',x0,'p',sgnp,ABS(p0),&
       'n',nb,'k',kp,'dt',dt0,'tf',tfin,'xm',sgnx,ABS(xmin),'xM',xmax&
       ,'g',g
  WRITE(filedistr,102)'d','_','s',sigma,'x',x0,'p',sgnp,ABS(p0),&
       'n',nb,'k',kp,'dt',dt0,'tf',tfin,'xm',sgnx,ABS(xmin),'xM',xmax&
       ,'g',g
102 FORMAT(a1,a1,2(a1,e8.3),a1,a1,e8.3,&
         a1,i6.6,a1,i3.3,a2,e8.3,a2,e8.3,a2,a1,e8.3,a2,e8.3,a1,e8.3)
  WRITE(*,*)fileoutput
  WRITE(*,*)filenorm
  WRITE(*,*)filedistr

END SUBROUTINE filescreate


!*******************************************************************
! Storing the distribution at a time given by cont*nshow*dt in a 
!   file created here
!*******************************************************************
SUBROUTINE store_distr(cont)

  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module
  USE filenames

  IMPLICIT NONE
  integer :: indexb
  INTEGER(ITT) :: cont,i,j,iq
  COMPLEX(CT) :: sumac
  CHARACTER(100) :: filename


  WRITE(filename,100)cont,'.dat'
100 FORMAT(i6.6,a4)
  IF(cont==-1)filename=filedistr    !storing the final distribution
  OPEN(14,file=filename)
  IF(cont==-1) THEN
     WRITE(14,'(a)')'# Step potential: |psi(x)|^2' 
     WRITE(14,'(3(a8,e14.6))')'# sigma=',sigma,'     x0=',x0,'      p0=',p0 
     WRITE(14,'(3(a8,i9))')'# ntot= ',ntot,'  nlarg=',nlarg,'  nlard=',nlard 
     WRITE(14,'(a5,i5,a6,i5,a9,i5,a15,a1)')'# nb=',nb,'   kp=',kp,'   nquad=',nquad,'   print_all=  ',print_all
     WRITE(14,'(a7,e14.6,a8,e14.6,a7,e14.6,a8,e14.6)')&
          '# xmin=',xmin,'   xmax=',xmax,'   dt0=',dt0,'   tfin=',tfin 
     WRITE(14,'(a)')&
          '#        x                        |psi(x)|^2                  Re psi(x)                Im psi(x)'  
  END IF
  IF(print_all=='y') THEN
     DO i=nq0+1,nq0+nmax*nquad
        sumac=czero
        DO j=1,nmax
           sumac=sumac+CMPLX(cpsi0r(j),cpsi0i(j),CT)*bba(indexb(j,i))
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
           sumac=sumac+CMPLX(cpsi0r(j),cpsi0i(j),CT)*bba(indexb(j,iq))
!           sumac=sumac+CMPLX(cpsi0r(j-1),cpsi0i(j-1),CT)*bb(j,(i-1)*nquad+1)
        END DO
        WRITE(14,*) xd(iq),ABS(sumac)**2,REAL(sumac,RT), REAL(-ione*sumac,RT)
     END DO
  ELSE IF (print_all.EQ.'n'.AND.cont.EQ.-1) THEN
     !the final output file 'd_s...' contains the full wave function psi(x)
     DO i=nq0+1,nq0+nmax*nquad
        sumac=czero
        DO j=1,nmax
           sumac=sumac+CMPLX(cpsi0r(j),cpsi0i(j),CT)*bba(indexb(j,i))
        END DO
        WRITE(14,*) xd(i),ABS(sumac)**2,REAL(sumac,RT), REAL(-ione*sumac,RT)
     END DO
  END IF
  CLOSE(14)

END SUBROUTINE store_distr




!*******************************************************************
!*******************************************************************
! Calculating the norm of the wave funtion
!*******************************************************************

SUBROUTINE calculate_norm(suma,cont)

  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module
  USE filenames

  IMPLICIT NONE
  integer :: indexb
  INTEGER(ITT) :: cont,i,j,iq
  REAL(RT) :: suma,suma1
  COMPLEX(CT) :: sumac(nq)

  !calculating the wave function at each quadrature point
  sumac=czero
  DO i=nq0+1,nq0+nmax*nquad
     !sumac(i)=czero
     DO j=1,nmax
        sumac(i)=sumac(i)+CMPLX(cpsi0r(j),cpsi0i(j),CT)*bba(indexb(j,i))
     END DO
  END DO
  !norm of wave function
  suma=0.0_RT
  suma1=0.0_RT
  DO j=1,l
     DO iq=1,nquad
        suma=suma+wd(iq)*ABS(sumac((j-1)*nquad+iq))**2*dx/2.0_RT
        IF(cont.EQ.0)&
             suma1=suma1+wd(iq)*ABS(psi0((j-1)*nquad+iq))**2*dx/2.0_RT
     END DO
  END DO
  IF(cont.EQ.0)&
       WRITE(*,*)'norm^2 of the wave function must be 1:',suma,suma1

END SUBROUTINE calculate_norm

!*******************************************************************
! SUBROUTINE crank_nicolson_propagator_pbc
! Time propagation of the wave packet using Caley form of the infinitesimal 
! time evolution operator. This is a unitary  method:
!  (1+i*dt/(2*hbar)H)psi(x,t+dt)=(1-i*dt/(2*hbar)H)psi(x,t)
! in scaled units:
!  (1+i*dt/2*H)psi(x,t+dt)=(1-i*dt/2*hbarH)psi(x,t)
! where the scaled Hamiltonian is, as described in the main program,
!  H'=-1/2*d^2/dx^2+V(x)
! The representation is achieved with bsplines. The identity is  
! the overlap matrix
! The solution of the problem requires backward substitution:
!   The subroutine lddc performs the L*LT decomposition of the complex
!   symmetric matrix (1+i*dt/(2*hbar)H) which is stored in the arrays
!   ar and ai;
!   The subroutine axdc solves the problem Ax=b by backward,
!   substitution where A is the L*LT decomposition of (1+i*dt/(2*hbar)H)
!   which is stored in ar and ai
! This subroutine uses Lapack subroutines for matrix-vector and matrix-matrix
!   products and for the solution of linear systems of equations
! REVISAR: tonto jma 30.10.2018 
!*******************************************************************
SUBROUTINE crank_nicolson_propagator_pbc

  USE nrtype
  USE prob_param
  USE filenames
  USE detbx
  USE detbxcx
  USE detax
  USE varindexb
  USE bspline_module
  USE mtx_pbc
  USE quad_module
  IMPLICIT NONE

  INTEGER(ITT) :: j,cont,it,i,ipvt(kp-1),info,time
  REAL(RT) :: dt,anorm1,t,sumac
  REAL(RT) :: timeunitary1,timeunitary2,timeunitary3
  REAL(RT),ALLOCATABLE :: v1r(:),v1i(:),wr(:),wi(:),mr(:),mi(:)
  COMPLEX(CT),ALLOCATABLE :: w1(:,:),w2(:,:),w(:,:),maux0(:,:),yc(:)
  COMPLEX(CT) :: z0(kp-1),z1(kp-1),z2(kp-1),z3(kp-1)
  COMPLEX(CT) :: xx0(kp-1),x1(kp-1),x2(kp-1),x3(kp-1)
  integer, parameter :: out_unit=20
  WRITE(*,*)'*********************************************************'
  WRITE(*,*)'******     STARTING CRANK-NICOLSON PROPAGATOR      ******' 
  WRITE(*,*)'*********************************************************'
  !scalefactor calculates the scale factor taking into acount the time step change
  CALL scalefactor
  dt=dt0
  !nt=int(tfin/dt)
  WRITE(nfileoutput,'(a)')         '******   CRANK-NICOLSON PROPAGATOR  ******'
  WRITE(nfileoutput,'(a24,e22.14)')'Time step    dt        =',dt
  WRITE(nfileoutput,'(a24,1x,i9)') 'Number of time steps nt=',nt   
  !matrices constructed in fill_mtxs
  WRITE(*,*)'nlarg,nlard,ntot',nlarg,nlard,ntot
  !ARRAYS: **************************************!
  IF(ALLOCATED(ntlb))DEALLOCATE(ntlb)
  ALLOCATE(ntlb(0:ntot))
  ntlb(0)=-1
  ALLOCATE(v1r(ntot),v1i(ntot),mr(ntot),mi(ntot),wr(ntot),wi(ntot))
  v1r=zero;v1i=zero;mr=zero;mi=zero;wr=zero;wi=zero;
  ALLOCATE(w(ntot,kp-1),w1(kp-1,kp-1),w2(kp-1,kp-1))
  w=czero;w1=czero;w2=czero;
  ALLOCATE(yc(ntot))
  yc=czero
  !END ARRAYS: **************************************!
  !matrices' construction
  t=0.0_RT
  !call fill_mtxs(t)
!stop  

  !L*LT decomposition of ar+i*ai is performed once and stored in the same arrays
  nnggrr=ntot
  !write(*,*)'ar',ar(:,4)
  !write(*,*)'ai',ai(:,4)
  !CALL lddc(ntot,nlard)
  !WRITE(*,*)'after LLT decomposition'
!!$  !***************************************************************
!!$  !checking backwards substitution
!!$  call axdc(nx1,nlard,psi0r,psi0i,xsr,xsi,v1r,v1i)
!!$  !end do
!!$  do i=1,nx1
!!$     write(303,101)i,xsr(i),xsi(i)
!!$  end do
!!$  !if everything is correct then x=xsr+I*xsi and b=psi0r+I*psi0i
!!$  !must satisfy the equation Ax=b
!!$  !xsr=1.0d0;xsi=0.0d0;v1r=0.d0;v1i=0.d0
!!$  call bxdrcakcx(nx1,nlard,xsr,xsi,v1r,v1i) !for this test dr=ar and di=ai
!!$  !we just performed the product Ax=v; v=v1 must be equal to b=psi0
!!$  do i=1,nx1
!!$     write(302,*)xmin+i*dx,v1r(i)**2+v1i(i)**2
!!$     write(300,101)i,v1r(i),v1i(i)
!!$     !write(300,101)i,v1r(i),psi0r(i),v1i(i),psi0i(i)
!!$     !write(300,101)i,(v1r(i)-psi0r(i))**2+(v1i(i)-psi0i(i))**2
!!$     !write(300,*)i,psi0r(i),psi0i(i)
!!$101 format(i5,4e20.12)
!!$  end do
!!$  stop
!!$  !***************************************************************
  !CALL filescreate
  !CALL tiempo(timeunitary1)
  OPEN(4,file=filenorm)
  cont=0
  it=0
  CALL calculate_norm(anorm1,cont)
  WRITE(*,*)'norm^2= ', it*dt, anorm1,nt
  WRITE(4,*)it,it*dt,anorm1
  CALL store_distr(cont)
  tfin=dble(dt0*(nt1)) +dble(dtada*(nt-nt1))
  CALL tiempo(timeunitary1)
!!$  WRITE(*,*)'ACHTUNG: nt=1'
!!$  nt=1
  DO it = 1 , nt 
     !# saves current state and reset potential for predictor corrector (eqs 53 and 54 of the report)
     psi2=zero
     cpsi0raux=cpsi0r
     cpsi0iaux=cpsi0i
     !#
     IF(it.EQ.11)THEN
        CALL tiempo(timeunitary2)
        timeunitary2=(timeunitary2-timeunitary1)/10.0d0
        !timeunitary2=timeunitary2*1.869_RT !take into account change in presicion
        WRITE(*,99)'Estimated time for iteration (h,m,s)=',&
             timeunitary2*nt/3.6d3,timeunitary2*nt/6.0d1,timeunitary2*nt
     END IF
     !predictor corrector steps
     DO pc=1,2
        
 
        !## calculates matrices for the predictor pc=1, and corrector pc=2 steps.
        call fill_mtxs(it)
        !## pc=1 does nothing, pc=2 uses the Phi for the time evolution 54 OF report
        cpsi0r=cpsi0raux
        cpsi0i=cpsi0iaux
        !###########################################################################3
        !The following  Block solves the schrodinger equation for the present predictor
        !corrector step
        !###########################################################################3
        CALL lddc(ntot,nlard)

        !product (1-i*dt/2*H)*psi0==(B-i*dt/2*H)*cpsi0=(dr+i*di)*cpsi0     
        !y_Delta=B_Delta*x_Delta, x_Delta=cpsi0
        CALL bxdrcakcx(ntot,nlard,cpsi0r,cpsi0i,xsr,xsi)
        !xs=(dr+i*di)*cpsi0
        !matrix-vector product (M1^+)**T*vdelta
        z0=cpsi0r(ntot+1:nmax)+ione*cpsi0i(ntot+1:nmax)
        z1=czero; z2=czero; xx0=czero
        CALL zgemv('N', kp-1, kp-1, cone, m1m,kp-1, z0, 1, czero, z1, 1)
        CALL zgemv('N', kp-1, kp-1, cone, m2m,kp-1, z0, 1, czero, z2, 1)
        x1=cpsi0r(1:kp-1)+ione*cpsi0i(1:kp-1)
        x2=cpsi0r(ntot-kp+2:ntot)+ione*cpsi0i(ntot-kp+2:ntot)
        CALL zgemv('T', kp-1, kp-1, cone, m1m,kp-1, x1, 1, czero, xx0, 1)
        CALL zgemv('T', kp-1, kp-1, cone, m2m,kp-1, x2, 1, cone, xx0, 1)
        CALL zgemv('N', kp-1, kp-1, cone, m0m,kp-1, z0, 1, cone, xx0, 1)
        !here: xx0=m0m*z0+m1m*x1+mm2*x2
        xsr(1:kp-1)=xsr(1:kp-1)+REAL(z1,RT)
        xsi(1:kp-1)=xsi(1:kp-1)+REAL(-ione*z1,RT)
        xsr(ntot-kp+2:ntot)=xsr(ntot-kp+2:ntot)+REAL(z2,RT)
        xsi(ntot-kp+2:ntot)=xsi(ntot-kp+2:ntot)+REAL(-ione*z2,RT)
        !here: xs(1:ntot) joined xx0 --> contains M^-*cpsi0
        !backward substitution: solution of the problem Ax=b, with A=L*LT
        !b=xsr+I*xsi; x=cpsi0r+I*cpsi0i; v1r,v1i are auxiliar arrays
        CALL axdc(ntot,nlard,xsr,xsi,cpsi0r,cpsi0i,v1r,v1i)
        DO j=1,kp-1
           mr(1:kp-1)=REAL(m1p(:,j),RT)
           mr(ntot-kp+2:ntot)=REAL(m2p(:,j),RT)
           mi(1:kp-1)=REAL(-ione*m1p(:,j),RT)
           mi(ntot-kp+2:ntot)=REAL(-ione*m2p(:,j),RT)
           CALL axdc(ntot,nlard,mr,mi,wr,wi,v1r,v1i)
           w(:,j)=CMPLX(wr,wi,CT)
           w1(:,j)=w(1:kp-1,j)
           w2(:,j)=w(ntot-kp+2:ntot,j)
        END DO
        !contsruction of the matrix M0^+-(M1^+)**T*W (pdf report)
        maux0=m0p
        CALL ZGEMM('T','N',kp-1,kp-1,kp-1,cmone,m1p,kp-1,w1,kp-1,cone,maux0,kp-1)
        CALL ZGEMM('T','N',kp-1,kp-1,kp-1,cmone,m2p,kp-1,w2,kp-1,cone,maux0,kp-1)
        !maux0 is overwritten: maux0=cmone*(m1p)**T*wauxr+cone*m0p (pdf report)
        z0=xx0
        z1=cpsi0r(1:kp-1)+ione*cpsi0i(1:kp-1)
        z2=cpsi0r(ntot-kp+2:ntot)+ione*cpsi0i(ntot-kp+2:ntot)
        !matrix-vector product (M1^+)**T*vdelta
        CALL zgemv('T', kp-1, kp-1, cmone, m1p,kp-1, z1, 1, cone, z0, 1)
        CALL zgemv('T', kp-1, kp-1, cmone, m2p,kp-1, z2, 1, cone, z0, 1)
        !output: z0=z0-(M1^+)**T*vdelta         
        !solving maux0*z0=y->z0; 
        CALL  zgesv (kp-1, 1, maux0, kp-1, ipvt, z0, kp-1, info )
        !write(*,*)'x',ipvt
        !z0 = cpsi0(ntot+1,nmax)
        !matrix vector product vdelta-W*z0=cpsi0-W*z0
        CALL zgemv('N', ntot, kp-1, cmone, w,ntot, z0, 1, czero, yc, 1)
        cpsi0r(1:ntot)=cpsi0r(1:ntot)+REAL(yc(1:ntot),RT)
        cpsi0i(1:ntot)=cpsi0i(1:ntot)+REAL(-ione*yc(1:ntot),RT)
        cpsi0r(ntot+1:nmax)=REAL(z0,RT)
        cpsi0i(ntot+1:nmax)=REAL(-ione*z0,RT)
        !###########################################################################3
     !cpsi0 contains the new expansion coefficients for pc=2
     !for pc=1 contains the intermediate step Phi eq. 53 of report
     END DO
     !save states in nshow before time step change
   !  IF((MOD(it,nshow)==0).and.(it.lt.nt1)) THEN 
   !     cont=cont+1
   !     time=cont*10
   !     CALL calculate_norm(anorm1,cont)
   !     WRITE(*,*)'norm^2= ', it*dt,time,anorm1
   !     time=cont*10
   !     WRITE(4,*)it,it*dt,anorm1
   !     If (discoef.eq.1)then
   !        CALL store_coef('evol',time)
   ! 	else  
    ! 	   CALL store_distr(time)
     !	end if 
       ! CALL store_coef('evol',time)
    ! END IF
     !save states in nshow after time step change
    ! IF((MOD(it,nshow)==0).and.(it.gt.nt1)) THEN
    !    cont=cont+1
    !    time=int(16.9/(dt0*nshow))
    !    time=int((10*nt1)/nshow)
    !    time=time+(cont-int((nt1)/nshow))
    !    CALL calculate_norm(anorm1,cont)
    !    WRITE(*,*)'norm^2= ', dble(dt0*(nt1)) +dble(dtada*(it-nt1)),time,anorm1
   !     WRITE(4,*)it,it*dt,anorm1
        
      !  If (discoef.eq.1)then
     !      CALL store_coef('evol',time)
    !	else  
     !	   CALL store_distr(time)
     !	end if 
     !  ! CALL store_coef('evol',time)
    ! END IF
    !##############################################
    !MODIFICACION
        !##############################################
      IF(((MOD(aat(it),0.02_RT).lt.aat(it+1)-aat(it)).and.(aat(it)-0.1_RT.le.aat(it+1)-aat(it)))) THEN 
        cont=cont+1
        time=2*cont!*10
        CALL calculate_norm(anorm1,cont)
        WRITE(*,*)'norm^2= ', it*dt,aat(it),time,anorm1
        !time=2*cont!*10
        WRITE(4,*)it,it*dt,anorm1
        If (discoef.eq.1)then
           CALL store_coef('evol',time)
    	else  
     	   CALL store_distr(time)
     	end if 
       ! CALL store_coef('evol',time)
     END IF
    ! CALL tiempo(taux6) 
     !taux6=taux6-taux
     !write(*,*)"writing time is ", taux6
     !stop
     !save states in nshow after time step change
    ! IF((MOD(it,nshow)==0).and.(it.gt.nt1)) THEN
     !if ((pc.eq.1).and.((MOD(aat(t),0.1_RT).lt.aat(t+1)-aat(t)).and.(aat(t)-0.1_RT.le.aat(t+1)-aat(t))) )then
  !if ((pc.eq.1).and.((MOD(aat(t),0.01_RT).lt.aat(t+1)-aat(t)).and.(aat(t)-0.1_RT.gt.aat(t+1)-aat(t))).or.(t.eq.nt) )then
     IF(((MOD(aat(it),0.1_RT).lt.aat(it+1)-aat(it)).and.(aat(it)-0.1_RT.gt.aat(it+1)-aat(it))).or.(it.eq.nt)) THEN
        !cont=cont+1
        !time=int(16.9/(dt0*nshow))
        !time=int((10*nt1)/nshow)
        !time=time+(cont-int((nt1)/nshow))
        cont=cont+10
        time=cont!*10
        CALL calculate_norm(anorm1,cont)
        WRITE(*,*)'norm^2= ', dble(dt0*(nt1)) +dble(dtada*(it-nt1)),aat(it),time,anorm1
        WRITE(4,*)it,it*dt,anorm1
        
        If (discoef.eq.1)then
           CALL store_coef('evol',time)
    	else  
     	   CALL store_distr(time)
     	end if 
       ! CALL store_coef('evol',time)
     END IF
     
        !##############################################
         !END MODIFICACION
        !##############################################
  END DO
  CALL calculate_norm(anorm1,cont)
  WRITE(*,*)'norm^2= ', (nt+1)*dt, anorm1
  WRITE(4,*)nt+1,(nt+1)*dt,anorm1
  CALL store_coef('final',cont)
   



  !storing the final distribution
  CALL store_distr(-1)
  CALL tiempo(timeunitary3)
  timeunitary3=(timeunitary3-timeunitary1)
  WRITE(*,99)'Time for series exp. propag. (h,m,s)=',&
       timeunitary3/3.6d3,timeunitary3/6.0d1,timeunitary3
  DEALLOCATE(ntlb)
  DEALLOCATE(v1r,v1i,mr,mi,wr,wi)
  DEALLOCATE(w,w1,w2,yc)

99 FORMAT(a37,f10.5,f13.5,f20.10)

END SUBROUTINE crank_nicolson_propagator_pbc


!****************************************************************
! SUBROUTINE ini_mtxs
! matrices ar(nlarg,ntot), ai(nlarg,ntot) and d(nlard,ntot)
!****************************************************************
SUBROUTINE ini_mtxs
!======================================================================
! Initializes different matrices used in subroutines fill_matrix and psi2potencial(al)
! and calculates the matrix representation of the laplace operator and the overlap matrix
!----------------------------------------------------------------------

! output:
! dr		= diagonal block of representation matrix of overlap matrix B dim(nlarg,ntot) (real) 
! dxx		= diagonal block of representation matrix of the laplace operator dim(nlarg,ntot) (real)
! dxx0,dxx1,dxx2= of diagonal block representations matrix of the laplace operator dim(kp-1,kp-1) (real)
! 27.05.2021; victor loaiza
!=======================================================================
  !construction of the matrices for diagonalization with Lanczos
  !and for time propagation using backwards substitution 
  USE nrtype
  USE bspline_module
  USE prob_param
  USE detax
  USE detbx
  USE detbxcx
  USE mtx_pbc

  IMPLICIT NONE

  INTEGER(ITT) :: i,j,ia,ib
  REAL(RT),ALLOCATABLE :: aaux(:,:)
!  COMPLEX(CT),ALLOCATABLE :: aauxc(:,:)

  ntot=nmax-kp+1  !for periodic boundary conditions
  nlarg=kp
  nlard=kp
  WRITE(*,*)'ntot=',ntot
  WRITE(nfileoutput,'(a)')         '******        SUB. ini_mtxs         ******'
  WRITE(nfileoutput,'(a24,1x,i6)') 'Banded matrices    ntot=',ntot
  WRITE(nfileoutput,'(a24,1x,i6)') 'Banded matrices   nlarg=',nlarg
  WRITE(nfileoutput,'(a24,1x,i6)') 'Banded matrices   nlard=',nlard

  IF (ntot.LE.kp)THEN
     WRITE(*,*)'stop: ntot.le.kp',ntot,kp
     STOP
  END IF
  ALLOCATE(ar(nlarg,ntot),ai(nlarg,ntot),dr(nlard,ntot),di(nlard,ntot),dxx(nlard,ntot))
  ALLOCATE(dxx0(kp-1,kp-1),dxx1(kp-1,kp-1),dxx2(kp-1,kp-1))
  ar=zero;ai=zero;dr=zero;di=zero
  !ar,ai: matrix M_Delta^+;  dr,di: matrix M_Delta^-
  !mauxb=block similar to $M_1^+$ corresponding to the overlap matrix B (report Bsplines) 
  !mauxh=block similar to $M_1^+$ corresponding to the hamiltonian matrix H (report Bsplines) 
  !m1p=$M_1^+$, m1m=$M_1^-$, m0p=$M_0^+$ and m0m=$M_0^-$ (report Bsplines)Matrix H

  ALLOCATE(m0m(kp-1,kp-1),m1m(kp-1,kp-1),m2m(kp-1,kp-1))
  ALLOCATE(m0p(kp-1,kp-1),m1p(kp-1,kp-1),m2p(kp-1,kp-1))
  ALLOCATE(m0aux(kp-1,kp-1),m1aux(kp-1,kp-1),m2aux(kp-1,kp-1))
  ALLOCATE(aa0(kp-1,kp-1),aa1(kp-1,kp-1),aa2(kp-1,kp-1))
  ALLOCATE(psi2(2*kp-1,nmax))  !!##?????
  psi2=zero
  m0p=czero;m1p=czero;m2p=czero;m0m=czero;m1m=czero;m2m=czero
  m0aux=czero;m1aux=czero;m2aux=czero
  aa0=zero;aa1=zero;aa2=zero

  !ALLOCATE(pot(2*kp-1,nmax))
  !pot=zero
  !call osc(pot)
  
  !construction of the overlap matrix initially stored in dr,di,m1m and m0m
  !Array overlap containing <B(jp,kp)|B(j,kp)> has been constructed in the 
  !   subroutine initial_wavepacketpbc
  !kinetic energy -0.5*<B(jp,kp)|d/dx^2|B(j,kp)>
  CALL cin_pbc   !Array kinetic(1:kp) containing the upper diagonal elements
  !!call integral_bibjbk !valid only for i
  call integral_bibjbk 
  
  !Overlap matrix ${\cal B}$ is composed of B_\Delta, B_1 y B_0 (see report) 
  DO i=0,kp-1
     DO j=i+1,ntot
        ib=nlard-i                  !position for the array d
        dr(ib,j)=dr(ib,j)+overlap(i+1)
     END DO
  END DO

  dxx=zero
  dxx1=zero;dxx2=zero;dxx0=zero
  
  DO i=0,kp-1
     DO j=i+1,ntot
        ib=nlard-i                  !position for the array d
        !dxx(ib,j)=dxx(ib,j)+pot(kp-i,j)!+kinetic(i+1)+dxx(ib,j)
        dxx(ib,j)=dxx(ib,j)+kinetic(i+1)
     END DO
  END DO
  
  !Matrix Dxx0
   DO i=1,kp-1
      DO j=i,kp-1
         !dxx0(i,j)=pot(kp+j-i,ntot+i) 
         dxx0(i,j)=kinetic(j-i+1) 
         dxx0(j,i)=dxx0(i,j)
      END DO
   END DO

!   !Matrices Dxx1 and Dxx2
   DO i=0,kp-2
      DO j=1,kp-1-i
         !dxx1(j,j+i)=pot(i+1,j)
         dxx1(j,j+i)=kinetic(kp-i)
      END DO
      DO j=1,i+1
         !dxx2(i+1,j)=pot(2*kp-2+j-i,ntot-kp+2+i)
         dxx2(i+1,j)=kinetic(kp-1+j-i)
      ENDDO
   END DO



END SUBROUTINE ini_mtxs

!****************************************************************
! SUBROUTINE fill_mtxs
! matrices ar(nlarg,ntot), ai(nlarg,ntot) and d(nlard,ntot)
!****************************************************************
SUBROUTINE fill_mtxs(t)
!======================================================================
! Calculation of the matrix representation of the hamiltonian
! more details at section 4 of report
!----------------------------------------------------------------------
! input..
! t		=current time (integer)

! hidden inputs:
! pc		= indicates step in predictor corrector scheme(integer)
! cpsi0r,cpsi0i	= real and imaginary parts of the vector of coeficients cpsi dim(nmax) (real)

! output:
! ar,ai		= diagonal block of representation matrix of the real and imaginary part of H1 dim(nlarg,ntot) (real)
!m0m,m1m,m2m,m0p,m1p,m2p = of diagonal block representations 
! 27.05.2021; victor loaiza
!=======================================================================
  !construction of the matrices for diagonalization with Lanczos
  !and for time propagation using backwards substitution 
  USE nrtype
  USE bspline_module
  USE prob_param
  USE detax
  USE detbx
  USE detbxcx
  USE mtx_pbc
  USE varindex
  IMPLICIT NONE

  INTEGER(ITT) :: i,j,ia,ib,t
  REAL(RT) :: At,dtaux
  REAL(RT),ALLOCATABLE :: aaux(:,:)

  IF (ntot.LE.kp)THEN
     WRITE(*,*)'stop: ntot.le.kp',ntot,kp
     STOP
  END IF
  !Overlap matrix ${\cal B}$ is composed of B_\Delta, B_1 y B_0 (see report) 
  !B_\Delta --> dr; B_1-->b1 and b2; B_0-->b0
  !matrices b0,b1,b2 constructed in subr. initial_wavepacketpbc
  !# calculates potential matrix elements
  
  !At=g*1.0d0!aat(t+(pc-1))
  At=aat(t+(pc-1)) !scale factor a(t), where a(t) has build in it the time step change (see more at subroutine scalefactor)
  !WRITE(*,*) At
  !stop
  call psi2potencial(At) !calculate matrix elements u_nm 
  

  !matrix elements of T+V with T=p^2/2
  !   ==> T+V=kin-0.5*aaux
  !A_\Delta --> ai; A_1-->aa1 and aa2; A_0-->aa0 (s. report)
  DO i=0,kp-1
     DO j=i+1,ntot
        ia=nlarg-i                  !position for the array ab
        ai(ia,j)=kinetic(i+1)+psi2(kp-i,j)
!        write (*,*) ai(ia,j),kinetic(i+1),psi2(kp-i,j)
     END DO
     
  END DO

  !Matrix aa0
  DO i=1,kp-1
     DO j=i,kp-1
        aa0(i,j)=kinetic(j-i+1)+psi2(kp+j-i,ntot+i)
        aa0(j,i)=aa0(i,j)
        !write (*,*) aa0(i,j),kinetic(j-i-1),psi2(kp-i+j,ntot+i)
     END DO
  END DO

  !Matrices aa1 and aa2
  DO i=0,kp-2
     DO j=1,kp-1-i
        aa1(j,j+i)=kinetic(kp-i)+psi2(i+1,j) 
       ! write (*,*) aa1(j,j+i),kinetic(kp-i),psi2(i+1,j)
     END DO
     DO j=1,i+1
        aa2(i+1,j)=kinetic(kp-1+j-i)+psi2(2*kp-2+j-i,ntot-kp+2+i)!
     END DO
  END DO
  dtaux=dt0
  IF ((t.GT.nt1) )THEN !time step change at nt1 see more at report section 2 or subroutine scalefactor
     dtaux=dtada
    
  END IF
   IF ((t.eq.nt1).AND.(pc.EQ.1) )THEN
     dtaux=dtada
    
  END IF
  ai=ai*dtaux*half
  ar=dr
  di=-ai
  aa0=aa0*dtaux*half
  aa1=aa1*dtaux*half
  aa2=aa2*dtaux*half
  m0m=b0-ione*aa0
  m1m=b1-ione*aa1
  m2m=b2-ione*aa2
  m0p=b0+ione*aa0
  m1p=b1+ione*aa1
  m2p=b2+ione*aa2
  !WRITE(*,*)'WARNING: STEP POTENTIAL IS NOT FULLY IMPLEMENTED'
  IF(ALLOCATED(ntl))DEALLOCATE(ntl)


END SUBROUTINE fill_mtxs



!****************************************************************
!         psi2potencial
!****************************************************************


subroutine psi2potencial(al)
!======================================================================
! Calculation of the matrix representation of the poisson potential
! by solving the poisson equation
! more details at section 5.1 of report
!----------------------------------------------------------------------
! input..
! al		=value of scale factor (real)

! hidden inputs: Are variables that heavely affect the subroutine.
! pc		= indicates step in predictor corrector scheme(integer)
! cpsi0r,cpsi0i	= real and imaginary parts of the vector of coeficients cpsi dim(nmax) (real)

! output:
! psi2		= representation matrix of the potential u_nm dim(nlarg,ntot) (real)
! 27.05.2021; victor loaiza
!=======================================================================
  use nrtype
  use prob_param
  use mtx_pbc
  use bspline_module
  use detbx !@!
  use detbxcx
  use quad_module
  implicit none 
  integer :: indexb
  integer(ITT) :: j,s,j1,j2,j3,sn,sm,ijk(3),po,i1
  complex(CT) :: auxc
  real(RT) :: auxb, ptom1, vp0r(nmax),vp0i(nmax),corr,al,vraux(nmax)
  integer, parameter :: out_unit=20

  CHARACTER(10) :: date,time,charx
  INTEGER ::  hh,mm,ss,sss,cc,yy,month,dd,timesss
  DOUBLE PRECISION  :: timereal1 ,timereal2
700 FORMAT(i2.2,i2.2,i2.2,a1,i3.3)                                      !time
701 FORMAT(i2.2,i2.2,i2.2,i2.2,i2.2) 
  ptom1=one/pc
  auxb=zero
  vp0r=zero
  vp0i=zero
  !Finds Abs(psi)^2 - 1 



 ! call tiempo(timereal1)
  DO j=1,nmax
     auxc=czero
     DO s=-(kp-1),kp-1
        
        j1=j+s
        IF(j1.GT.nmax)j1=j1-nmax
        IF(j1.LT.1)j1=nmax+j1
        do sm=MAX0(-(kp-1)+s,-(kp-1)),MIN0(kp-1,(kp-1)+s)
           j2=j+sm
           IF(j2.GT.nmax)j2=j2-nmax
           IF(j2.LT.1)j2=nmax+j2
           !ijk=(/0,s,sm/)
           !call sort_ijk(ijk)
           IF (sm.LT.0 .OR. s.LT.0) THEN
              !auxb=bijk(1,1+s-sm,1-sm )
              IF (s.LE.sm) THEN
                 auxb=bijk(1,1-s+sm,1-s )
              ELSE
                 auxb=bijk(1,1+s-sm,1-sm )
                 
              !ijk=(/1,1+s,1+sm /)
              END IF

           ELSE
              auxb=bijk(1,1+s,1+sm )
              !ijk=(/1,1+s,1+sm /)
           END IF
           !auxb=bijk(ijk(1),ijk(2),ijk(3))
           auxc=auxc+CMPLX(cpsi0r(j2),-cpsi0i(j2),CT)&
                   *CMPLX(cpsi0r(j1),cpsi0i(j1),CT)*auxb
        end do       
     end DO
     
     vp0r(j)=-bint + auxc
     
  end DO 


  !Solves Poisson with -(Abs(psi)^2 - 1) as source to find u(x)
  d=MONE*dxx
  call solver(vp0r,vp0i,d,MONE*dxx0,MONE*dxx1,MONE*dxx2)
  vp0r=HALF*vp0r!multiplies by 0.5 since the matrix dxx correspond to the operator K which is -0.5 (d/dx)^2

  !impose condition u(L)=0 on u(x)  
  i1=nq0+nquad*nmax
  !!!i1=index_xd(i1)
  auxb=zero
  !auxb=vp0r(1)*bba(indexb(1,i1))
  DO s=ntot,nmax
     auxb=auxb+vp0r(s)*bba(indexb(s,i1))
     !auxb=auxb+vp0r(s)*bba(indexb(s,nq0+))
  END DO

  !impose condition u(k=0)=0 on u(x)  
  auxb=zero
  DO j=1,nmax
     auxb=auxb+vp0r(j)
  END DO
  auxb=auxb/nmax
  DO j=1,nmax
     vp0r(j)=-auxb+vp0r(j)
  END DO





!calculate the matrix elements u_nm from u_i


  auxb=zero
  DO j=1,nmax
     auxc=czero
     DO s=0,kp-1
        auxc=czero
        j1=j+s
        IF(j1.GT.nmax)j1=j1-nmax
        do sm=-(kp-1)+s,kp-1
           j2=j+sm
           IF(j2.GT.nmax)j2=j2-nmax
           IF(j2.LT.1)j2=nmax+j2

           !ijk=(/j,j1,j2/)

           !call sort_ijk(ijk)
          
           !ijk=(/1,1+1,2+1 /)
           !ijk=(/1,s+1,sm+1 /)
           IF (sm.LT.0) THEN
              auxb=bijk(1,1+s-sm,1-sm )
           ELSE
              auxb=bijk(1,1+s,1+sm )
              !ijk=(/1,1+s,1+sm /)
           END IF
           !IF(s.GT.sm)ijk=(/1,sm+1,s+1 /)
           !IF(sm.LT.0)ijk=(/1,-sm+1,s+1-sm /)
           !ijk=(/1,1+1,2+1 /)
           !auxb=bijk(ijk(1),ijk(2),ijk(3))
           auxc=auxc+vp0r(j2)*auxb

        end do
        psi2(kp+s,j)=ptom1*(al*auxc+psi2(kp+s,j))
      
        psi2(kp-s,j1)=psi2(kp+s,j)!#j1
        
     end DO
     !stop
  end DO
 
!  stop
 

end subroutine psi2potencial

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
  END DO

  DO m=1,kp
     v(m)=x(1)
  END DO
  DO m=1,nb-kp
     v(m+kp)=x(m+1)
  END DO
  DO m=1,kp
     v(nb+m)=x(ll+1)
  END DO

END SUBROUTINE knots

!=================================================================!
! Calculation of important parameters
! Generates screen output of the input and other parameters
!=================================================================!
SUBROUTINE print_output

  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module
  USE filenames

  IMPLICIT NONE

  
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

  OPEN(nfileoutput,file=fileoutput, status='old', action='write')
  WRITE(nfileoutput,'(a)')         '=========================================='
  WRITE(nfileoutput,'(a)')         '******       SUB. print_output      ******'
  WRITE(nfileoutput,'(a)')         '******  DATA FOR THE INITIAL STATE  ******'
  WRITE(nfileoutput,'(a24,e22.14)')         'sigma                  =',sigma
  WRITE(nfileoutput,'(a24,e22.14)')         'x0                     =',x0
  WRITE(nfileoutput,'(a24,e22.14)')         'p0                     =',p0
  WRITE(nfileoutput,'(a)')         '******  DATA FOR THE BSPLINE BASIS  ******'
  WRITE(nfileoutput,'(a24,1x,i6)') 'Order kp               =',kp
  WRITE(nfileoutput,'(a24,1x,i6)') 'Number of bsplines nb  =',nb
  WRITE(nfileoutput,'(a24,1x,i6)') 'Number of breakpoints l=',l
  WRITE(nfileoutput,'(a24,1x,i6)') 'Inner breakpoints l0   =',l0
  WRITE(nfileoutput,'(a24,1x,i6)') 'Matrix size nmax       =',nmax
  WRITE(nfileoutput,'(a)') 'Uniform breakpoints'
  WRITE(nfileoutput,'(a26,f7.2,a1,f7.2,a1)')'Endpoints (xmin,xmax)  = (',xmin,',',xmax,')'
  WRITE(nfileoutput,'(a24,e22.14)')'Grid spacing dx        =',dx
  WRITE(nfileoutput,'(a)')         '******   DATA FOR THE QUADRATURES   ******'
  WRITE(nfileoutput,'(a24,1x,i6)') 'Num. of quad pts nquad =',nquad
  WRITE(nfileoutput,'(a24,1x,i6)') 'Total quadr. points nq =',nq
  WRITE(nfileoutput,'(a24,1x,i6)') 'Initial quad point nq0 =',nq0
  WRITE(nfileoutput,'(a24,1x,i6)') '# quad points on k nqk =',nqk
  WRITE(nfileoutput,'(a)')         '******       TIME PROPAGATION       ******'
  WRITE(nfileoutput,'(a24,e22.14)')'Time step      dt0     =',dt0
  WRITE(nfileoutput,'(a24,e22.14)')'Final time     tfin    =',tfin
  WRITE(nfileoutput,'(a)')         '******       OTHER PARAMETERS       ******'
  WRITE(nfileoutput,'(a24,1x,i6)') 'nshow                  =',nshow
  WRITE(nfileoutput,'(a24,1x,i6)') 'V0_switch              =',V0_switch
  WRITE(nfileoutput,'(a24,e22.14)')'V0                     =',V0

END SUBROUTINE print_output

!===================================================================
! multiplication of a band matrix and a vector.
! the real and imaginary parts of the banded matrix are stored in 
! the arrays dr and di, respectively. 
!INPUT
!   ntot: integer -> matrix dimension
!   nlard: integer -> band width of the matrix
!   xer,xei: arrays containing the real and imaginary parts of the 
!            input vector
!OUTPUT:
!   xsr,xsi: arrays containing the real and imaginary parts of the 
!            product of the matrix with the vector
!===================================================================
SUBROUTINE bxdrcakcx(ntot,nlard,xer,xei,xsr,xsi)

  USE nrtype
  USE detbxcx
  USE varindexb
  
  INTEGER(ITT) :: ntot,nlard
  REAL(RT) :: xer(ntot),xei(ntot),xsr(ntot),xsi(ntot)
  INTEGER i,j,j0
  DOUBLE PRECISION xr,xi
  
  !  ======================================================
  !  at the first call, ntlb(0)=-1 and the effective bandwidth of B is computed
  ! Then ntlb(0).ne.-1 and the computed bandwidth is used for later calls
  !  ======================================================
  IF(ntlb(0).EQ.-1) THEN    
     !WRITE(6,1)
1    FORMAT(' Computes effective bandwidth of matrix B...')
     ntlb(0)=nlard-1
     DO  i=1,ntot
        DO  j=max0(1,nlard+1-i),nlard
           IF (dr(j,i).NE.0.d0.OR.di(j,i).NE.0.d0) EXIT
        END DO
        ntlb(i)=min0(j,ntlb(i-1)+1)
     END DO
  END IF
  DO  i=1,ntot
     xr=0.d0
     xi=0.d0
     j0=i-nlard
     DO j=ntlb(i),nlard
        xr=xr+xer(j+j0)*dr(j,i)-xei(j+j0)*di(j,i)
        xi=xi+xei(j+j0)*dr(j,i)+xer(j+j0)*di(j,i)
     END DO
     xsr(i)=xr
     xsi(i)=xi	
  END DO
  
  DO  i=1,ntot
     xr=xer(i)
     xi=xei(i)
     j0=nlard-i
     DO  j=ntlb(i)-j0,i-1
        xsr(j)=xsr(j)+xr*dr(j+j0,i)-xi*di(j+j0,i)
        xsi(j)=xsi(j)+xi*dr(j+j0,i)+xr*di(j+j0,i)
     END DO
  END DO

  RETURN

END SUBROUTINE bxdrcakcx

! ****************************************************************
!  08.04.2003
!  Subroutine tiempo
!  OUTPUT: timereal=tiempo en segundos del mes correspondiente
! ****************************************************************
SUBROUTINE tiempo(timereal)

  IMPLICIT NONE

  CHARACTER(10) :: date,time,charx
  INTEGER ::  hh,mm,ss,sss,cc,yy,month,dd,timesss,monthdays(12)
  DOUBLE PRECISION  :: timereal
    
  CALL DATE_AND_TIME(date,time)
  READ(time,700)hh,mm,ss,charx,sss
  READ(date,701)cc,yy,month,dd
  !timesss=sss+ss*1000+60000*mm+hh*3600000+dd*86400000
  monthdays=(/0,31,59,90,120,151,181,212,243,273,304,334/)
  timesss=sss+ss*1000+60000*mm+hh*3600000+dd*86400000+monthdays(month)*86400000
  timereal=DBLE(timesss)/1000
  !timereal=DBLE(timesss)

700 FORMAT(i2.2,i2.2,i2.2,a1,i3.3)                                      !time
701 FORMAT(i2.2,i2.2,i2.2,i2.2,i2.2)                                    !date

END SUBROUTINE tiempo
!*******************************************************************
! Storing the distribution at a time given by cont in a 
!   file created here
!*******************************************************************
SUBROUTINE store_coef(texting,cont)

  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module
  USE filenames

  IMPLICIT NONE

  INTEGER(ITT) :: cont,i,j,iq
  COMPLEX(CT) :: sumac
  CHARACTER(100) :: filename
  CHARACTER(4) ::texting

  WRITE(filename,100) initialtext,'Y',texting,'Y',cont,'.dat'
100 FORMAT(a10,a1,a4,a1,i6.6,a4)!100 FORMAT(a4,i6.6,a4)
  !IF(cont==-1)filename=filedistr    !storing the final distribution
  OPEN(14,file=filename)
  Do j=1,nmax
     WRITE(14,*) cpsi0r(j),cpsi0i(j)
  End do

  CLOSE(14)

END SUBROUTINE store_coef
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! indexb
!======================================================================
! takes the integer j of the bspline and ip the point to be evaluated
! and return the integer indexb  which used in bba(indexb) gives 
!the correct value of the Bspline evaluated at ip
! Bj(x(ip))=bba(index(j,ip))
!----------------------------------------------------------------------
! input..
! j,ip           = Bspline and point to be evaluated (integer)
! output:
! indexb 	 = corresponding integer to be evaluated at bba (integer)
!=======================================================================
!!!!!!!!!!!!!
INTEGER FUNCTION indexb(j,ip)

     USE nrtype
     USE bspline_module
     USE quad_module
     USE prob_param

     IMPLICIT NONE
  !INTEGER(ITT) :: indexb
     INTEGER(ITT) :: i,j,i0,a,b,iout,ip
     i0=nq0+(j-1)*nquad
     If (ip.lt.i0) then
        iout=ip+nquad*nmax
     else  
        iout=ip
     end if 

     IF(iout-i0.GT.nqk) THEN
     
        iout=0
     ELSE
        iout=iout-i0
     END IF
    


     indexb = iout
     RETURN  
END FUNCTION indexb

