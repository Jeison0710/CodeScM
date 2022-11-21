!jma: ORIGINAL SUBROUTINES BY ANA FRAPICCINI
!    SUBROUTINES FOR THE CALCULATION OF THE MATRIX ELEMENTS OF SOME OPERATORS 
!    IN A BSPLINE BASIS
!    The matrix elements are stored in an array du(2*kp-1,nmax)
!    kp: order of the Bsplines
!    nmax: dimension of the matrices which is given by nb-2, with
!          nb the number of bsplines
!    du(kp,:)   -> diagonal elements
!    du(kp+i,:) -> elements of the ith lower diagonal
!    du(kp-i,:) -> elements of the ith upper diagonal
!
!    FOR THE UPPER BAND THE TRANSFORMATION FROM THE MATRIX (a_{ij})_{nxn}
!    TO THE ARRAY du(:,:) READS
!        du(k-i,j)=a_{j-i,j} for i=0,..,kp-1 and j>=i   
!***************************************************
!Kinetic matrix for periodic boundary conditions
!Since the bspline functions are translationally invariant, the overlaps is calculated
! only for the first basis function $\tilde{B}_1$
!***************************************************
SUBROUTINE cin_pbc
  !*****Kinetic energy:
  !     -0.5*<B(jp,kp)|d²/dx²|B(j,kp)>=0.5*<d/dx(B(jp,kp))|d/dx(B(j,kp))>
  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module

  IMPLICIT NONE
  integer :: indexb
  INTEGER(ITT) :: i,j,iq,s,i0,i0p,k
  REAL(RT) :: a0

  !Matrix elements of one row are stored in kinetic(:)
  !Diagonal elements:   kinetic(1)= <d/dx B(j,kp)|d/dx B(j,kp)>
  ALLOCATE(kinetic(kp))
  !derivative of the basis function is stored in dbba(1:nqk)
  a0=dx*HALF
  kinetic=zero
  DO k=1,kp
     DO s=0,kp-1
        i0=s*nquad
        i0p=i0-(k-1)*nquad
        DO iq=1,nquad
           j=i0p+iq
           IF(j.LT.0)j=0
           kinetic(k)=kinetic(k)+wd(iq)*dbba(i0+iq)*dbba(j)*a0
        ENDDO
     ENDDO
     !write(*,*)k,kinetic(k)
  END DO
  !stop
  kinetic=kinetic*HALF
  !write(*,*)kinetic
  RETURN
END SUBROUTINE cin_pbc

!================================================================'
!overlap matrix for periodic boundary conditions
!Since the bspline functions are translationally invariant, the overlaps is calculated
! only for the first basis function $\tilde{B}_1$
SUBROUTINE over_pbc
  !*****Overlap matrix <B(jp,kp)|B(j,kp)>
  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module

  IMPLICIT NONE
  integer :: indexb
  INTEGER(ITT) :: i,j,iq,s,i0,i0p,k
  REAL(RT) :: a0

  !Matrix elements of one row are stored in overlap(:)
  !Diagonal elements:   overlap(1)= <B(j,kp)|B(j,kp)>
  ALLOCATE(overlap(kp))
  !basis function is stored in bba(1:nqk)
  a0=dx/2.0_RT
  overlap=zero


  DO k=1,kp
     DO s=0,kp-1
        i0=s*nquad
        i0p=i0-(k-1)*nquad
        DO iq=1,nquad
           j=i0p+iq
           IF(j.LT.0)j=0
           overlap(k)=overlap(k)+wd(iq)*bba(i0+iq)*bba(j)*a0
        ENDDO

     ENDDO
     !WRITE(*,*) k, overlap(k)
  END DO

  RETURN 
END SUBROUTINE over_pbc

!================================================================'
SUBROUTINE steppot_pbc(du)
  !***** <B(jp,kp)|V|B(j,kp)>, with V=0, x>0 and V=-1, x<0
  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module

  IMPLICIT NONE
  integer :: indexb
  INTEGER(ITT) :: i,j,iq,s,i1,i2,j1,j0
  REAL(RT) :: du(2*kp-1,nmax),a0
  REAL(RT), ALLOCATABLE :: dut(:,:)
  CHARACTER(30)::filetonto

  !only valid for periodic boundary conditions on a regular grid
  a0=dx/2.0_RT
  !Diagonal elements
  du=00_RT
  !WRITE(200,'(a)')&
  !     '   j  j1   s  iq   tx(j+s)     tx(j+s+1)    i1  indexb    xd_pbc         bba(j)      bba(j1)'

  !only diagonal and upper diagonal elements are calculated
  j0=2*kp-2
  DO i=0,kp-1
     !write(200,*)'i=',i
     DO j=1,nmax
        j1=j+i
        IF(j1.GT.nmax)j1=j1-nmax
        !WRITE(201,*)'# j,j1=',j,j1
        DO s=i,kp-1
           DO iq=1,nquad
              i1=nq0+(j+s-1)*nquad+iq
              i1=index_xd(i1)
              IF(xd(i1).LE.0.0_RT) THEN
                 du(kp+i,j)=du(kp+i,j)+&
                      wd(iq)*bba(indexb(j,i1))*bba(indexb(j1,i1))*a0
              END IF
           ENDDO
        ENDDO
        du(kp-i,j1)=du(kp+i,j)
     ENDDO
  END DO

  du=-du  !V=-1, x<0

  RETURN 
END SUBROUTINE steppot_pbc


!================================================================'
!================================================================'
SUBROUTINE osc(du)
  !*****Harmonic oscillator 0.5*<B(jp,kp)|x²|B(j,kp)>
  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module

  IMPLICIT NONE
  integer :: indexb
  INTEGER(ITT) ::  i,j,iq,s,i1,i2,j1,j0
  REAL(RT) :: du(2*kp-1,nmax),a0
  REAL(RT), ALLOCATABLE :: dut(:,:)
  a0=dx/2.0_RT
  !Diagonal elements
  du=00_RT
  !WRITE(200,'(a)')&
  !     '   j  j1   s  iq   tx(j+s)     tx(j+s+1)    i1  indexb    xd_pbc         bba(j)      bba(j1)'

  !only diagonal and upper diagonal elements are calculated
  j0=2*kp-2
  DO i=0,kp-1
     !write(200,*)'i=',i
     DO j=1,nmax
        j1=j+i
        IF(j1.GT.nmax)j1=j1-nmax
        !WRITE(201,*)'# j,j1=',j,j1
        DO s=i,kp-1
           DO iq=1,nquad
              i1=nq0+(j+s-1)*nquad+iq
              i1=index_xd(i1)
             
              du(kp+i,j)=du(kp+i,j)+&
                   wd(iq)*bba(indexb(j,i1))*bba(indexb(j1,i1))*a0*xd(i1)*xd(i1)
              
           ENDDO
        ENDDO
        du(kp-i,j1)=du(kp+i,j)
     ENDDO
  END DO
  
  du=du*0.5d0

  RETURN 
END SUBROUTINE osc

!================================================================'
SUBROUTINE cos2(du)
  !*****Cos**2 <B(jp,kp)|cos**2(x/2)|B(j,kp)>
  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module

  IMPLICIT NONE
  integer :: indexb
  INTEGER(ITT) ::  i,j,iq,s,i1,i2,j1,j0
  REAL(RT) :: du(2*kp-1,nmax),a0
  REAL(RT), ALLOCATABLE :: dut(:,:)
  a0=dx/2.0_RT
  !Diagonal elements
  du=00_RT
  !WRITE(200,'(a)')&
  !     '   j  j1   s  iq   tx(j+s)     tx(j+s+1)    i1  indexb    xd_pbc         bba(j)      bba(j1)'

  !only diagonal and upper diagonal elements are calculated
  j0=2*kp-2
  DO i=0,kp-1
     !write(200,*)'i=',i
     DO j=1,nmax
        j1=j+i
        IF(j1.GT.nmax)j1=j1-nmax
        !WRITE(201,*)'# j,j1=',j,j1
        DO s=i,kp-1
           DO iq=1,nquad
              i1=nq0+(j+s-1)*nquad+iq
              i1=index_xd(i1)
             
              du(kp+i,j)=du(kp+i,j)+&
                   wd(iq)*bba(indexb(j,i1))*bba(indexb(j1,i1))*a0*COS(0.5d0*xd(i1))*COS(0.5d0*xd(i1))
              
           ENDDO
        ENDDO
        du(kp-i,j1)=du(kp+i,j)
     ENDDO
  END DO
  
  

  RETURN 
END SUBROUTINE cos2


!================================================================'
!================================================================'


!================================================================'
SUBROUTINE mat_inter(du)
  !*****Dipole interaction (vel. gauge),nonsymmetric!!! mat_i(j,i)=-mat_i(i,j)
  !*****<B(jp,kp)|d/dx|B(j,kp)>
  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module

  IMPLICIT NONE
  integer :: indexb
  INTEGER(ITT) :: i,j,iq,s
  REAL(RT) :: du(2*kp-1,nmax),a0
  REAL(RT), ALLOCATABLE :: dut(:,:)

  !Diagonal elements
  du=0d0
  DO j=2,n+1
     DO s=0,kp-1
        IF(tx(j+s).NE.tx(j+s+1))THEN
           a0=(tx(j+s+1)-tx(j+s))/2d0
           DO iq=1,nquad
              du(kp,j-1)=du(kp,j-1)+&
                   wd(iq)*bb(j,(j+s-kp)*nquad+iq)*dbb(j,(j+s-kp)*nquad+iq)*a0
           ENDDO
        ENDIF
     ENDDO
  ENDDO

  !Upper diagonal elements
  DO i=1,kp-1
     DO j=i+2,n+1
        DO s=0,kp-i-1
           IF(tx(j+s).NE.tx(j+s+1))THEN
              a0=(tx(j+s+1)-tx(j+s))/2d0
              DO iq=1,nquad
                 du(kp-i,j-1)=du(kp-i,j-1)+wd(iq)*&
                      bb(j-i,(j+s-kp)*nquad+iq)*dbb(j,(j+s-kp)*nquad+iq)*a0
              ENDDO
           ENDIF
        ENDDO
        du(kp+i,j-i-1)=-du(kp-i,j-1)
     ENDDO
  ENDDO

  RETURN 
END SUBROUTINE mat_inter


!================================================================'
SUBROUTINE one_over_r2(du)
  !*****Harmonic oscillator <B(jp,kp)|1/x²|B(j,kp)>
  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module

  IMPLICIT NONE
  integer :: indexb
  INTEGER(ITT) :: i,j,iq,s
  REAL(RT) :: du(2*kp-1,nmax),a0
  REAL(RT), ALLOCATABLE :: dut(:,:)

  !Diagonal elements
  du=0d0
  DO j=2,n+1
     DO s=0,kp-1
        IF(tx(j+s).NE.tx(j+s+1))THEN
           a0=(tx(j+s+1)-tx(j+s))/2d0
           DO iq=1,nquad
              du(kp,j-1)=du(kp,j-1)+&
                   wd(iq)*bb(j,(j+s-kp)*nquad+iq)*bb(j,(j+s-kp)*nquad+iq)*&
                   1.0d0/(xd((j+s-kp)*nquad+iq))**2*a0
           ENDDO
        ENDIF
     ENDDO
  ENDDO

  !Upper diagonal elements
  DO i=1,kp-1
     DO j=i+2,n+1
        DO s=0,kp-i-1
           IF(tx(j+s).NE.tx(j+s+1))THEN
              a0=(tx(j+s+1)-tx(j+s))/2d0
              DO iq=1,nquad
                 du(kp-i,j-1)=du(kp-i,j-1)+&
                      wd(iq)*bb(j-i,(j+s-kp)*nquad+iq)*bb(j,(j+s-kp)*nquad+iq)*&
                      1.0d0/(xd((j+s-kp)*nquad+iq)**2)*a0
              ENDDO
           ENDIF
        ENDDO
        du(kp+i,j-i-1)=du(kp-i,j-1)
     ENDDO
  ENDDO

  RETURN 
END SUBROUTINE one_over_r2

!================================================================'
SUBROUTINE cin_surface(du)
  !*****Kinetic energy including the surface term:
  !     -0.5*<B(jp,kp)|d²/dx²|B(j,kp)>
  !     <B(jp,kp)|d²/dx²|B(j,kp)>=-<d/dx(B(jp,kp))|d/dx(B(j,kp))>
  !                               +B(jp,kp))*d/dx(B(j,kp))|_xmin^xmax
  ! The matrix is non-symmetric and therefore non hermitian
  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module

  IMPLICIT NONE
  integer :: indexb
  INTEGER(ITT) :: i,j,iq,s
  REAL(RT) :: du(2*kp-1,ntot),a0
  REAL(RT), ALLOCATABLE :: dut(:,:)

  !Diagonal elements
  du=0.d0
  DO j=1,n+2
     DO s=0,kp-1
        IF(tx(j+s).NE.tx(j+s+1))THEN
           a0=(tx(j+s+1)-tx(j+s))/2d0
           DO iq=1,nquad
              du(kp,j)=du(kp,j)+&
                   wd(iq)*dbb(j,(j+s-kp)*nquad+iq)*dbb(j,(j+s-kp)*nquad+iq)*a0
           ENDDO
        ENDIF
     ENDDO
     !now the surface term
     du(kp,j)=-du(kp,j)+bb_boundary(j,2)*dbb_boundary(j,2)-bb_boundary(j,1)*dbb_boundary(j,1)
  ENDDO

  !Upper diagonal elements
  DO i=1,kp-1
     DO j=i+2,n+2
        DO s=0,kp-i-1
           IF(tx(j+s).NE.tx(j+s+1))THEN
              a0=(tx(j+s+1)-tx(j+s))/2d0
              DO iq=1,nquad
                 du(kp-i,j)=du(kp-i,j)+wd(iq)*&
                      dbb(j-i,(j+s-kp)*nquad+iq)*dbb(j,(j+s-kp)*nquad+iq)*a0
              ENDDO
           ENDIF
        ENDDO
        du(kp-i,j)=-du(kp-i,j)&
             +bb_boundary(j-i,2)*dbb_boundary(j,2)-bb_boundary(j-i,1)*dbb_boundary(j,1)
        du(kp+i,j-i)=du(kp-i,j)
     ENDDO
  ENDDO

  du=du*(-0.5d0)

  RETURN
END SUBROUTINE cin_surface

!================================================================'
SUBROUTINE integral_bibjbkbl
  !*****integral:
  !     \int B_i(x)*B_j(x)*B_k(x)*B_l(x) dx
  !  stored in the array bijkl:
  !     bijkl(i,j,k,l) calculated ONLY FOR 1.LE.i.LE.j.LE.k.LE.l.LE.kp 
  !     otherwise it is zero
  !     see pdf report on bsplines
  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module

  IMPLICIT NONE
  integer :: indexb
  INTEGER(ITT) :: i,j,k,l1,i0,iq,s,cont
  REAL(RT) :: a0

  !Integrals are stored in b_ijkl(:)
  ALLOCATE(bijkl(kp,kp,kp,kp))
  !derivative of the basis function is stored in dbba(1:nqk)
  a0=dx/2.0_RT
  bijkl=zero
  cont=0
  DO i=1,kp
     DO j=i,kp
        DO k=j,kp
           DO l1=k,kp
              cont=cont+1
              !WRITE(*,*)cont,'|',i,j,k,l1
              DO s=0,kp-l1+i-1
                 i0=nq0+(l1-1+s)*nquad
                 !bba(indexb(i,ii))
                 DO iq=1,nquad
                    bijkl(i,j,k,l1)=bijkl(i,j,k,l1)+wd(iq)*bba(indexb(i,i0+iq))&
                         *bba(indexb(j,i0+iq))*bba(indexb(k,i0+iq))*bba(indexb(l1,i0+iq))*a0
                 END DO
               
              END DO
           END DO
        END DO
     END DO
  END DO


  RETURN

END SUBROUTINE integral_bibjbkbl

!================================================================'
SUBROUTINE integral_bibjbk
  !*****integral:
  !     \int B_i(x)*B_j(x)*B_k(x) dx
  !  stored in the array bijkl:
  !     bijk(i,j,k) calculated ONLY FOR 1.LE.i.LE.j.LE.k.LE.kp 
  !     otherwise it is zero
  !     see pdf report on bsplines
  USE nrtype
  USE prob_param
  USE bspline_module
  USE quad_module

  IMPLICIT NONE
  integer :: indexb
  INTEGER(ITT) :: i,j,k,i0,iq,s
  REAL(RT) :: a0

  !Integrals are stored in b_ijk(:)
  ALLOCATE(bijk(kp,kp,kp))
  !derivative of the basis function is stored in dbba(1:nqk)
  a0=dx/2.0_RT
  bijk=zero
  
  DO i=1,kp
     DO j=i,kp
        DO k=j,kp
           
           
           !WRITE(*,*)cont,'|',i,j,k
           bijk(i,j,k)=zero
           DO s=0,kp-k+i-1
              i0=nq0+(k-1+s)*nquad
              !bba(indexb(i,ii))
              DO iq=1,nquad
                 bijk(i,j,k)=bijk(i,j,k)+wd(iq)*bba(indexb(i,i0+iq))&
                      *bba(indexb(j,i0+iq))*bba(indexb(k,i0+iq))*a0
              END DO
               
           END DO
           bijk(j,i,k)=bijk(i,j,k)
           bijk(k,j,i)=bijk(i,j,k)
           bijk(i,k,j)=bijk(i,j,k)
        END DO
     END DO
  END DO
 ! write(*,*) bijk(3,2,1)
  !stop

  bint=zero
  DO s=0,kp-1
     i0=s*nquad
     DO iq=1,nquad
       j=i0+iq
        IF(j.LT.0)j=0
         bint=bint+wd(iq)*bba(j)*a0
     ENDDO

  ENDDO


  RETURN

END SUBROUTINE integral_bibjbk



