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
SUBROUTINE cin(du)
  !*****Kinetic energy:
  !     -0.5*<B(jp,kp)|d²/dx²|B(j,kp)>=0.5*<d/dx(B(jp,kp))|d/dx(B(j,kp))>
  use nrtype
  USE prob_param
  USE bspline_module
  USE quad_module

  IMPLICIT NONE

  INTEGER(ITT) :: i,j,iq,s
  REAL(RT) :: du(2*kp-1,nmax),a0
  REAL(RT), ALLOCATABLE :: dut(:,:)

  !Diagonal elements
  du=0.d0
  DO j=2,n+1
     DO s=0,kp-1
        IF(tx(j+s).NE.tx(j+s+1))THEN
           a0=(tx(j+s+1)-tx(j+s))/2d0
           DO iq=1,nquad
              du(kp,j-1)=du(kp,j-1)+&
                   wd(iq)*dbb(j,(j+s-kp)*nquad+iq)*dbb(j,(j+s-kp)*nquad+iq)*a0
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
                      dbb(j-i,(j+s-kp)*nquad+iq)*dbb(j,(j+s-kp)*nquad+iq)*a0
              ENDDO
           ENDIF
        ENDDO
        du(kp+i,j-i-1)=du(kp-i,j-1)
     ENDDO
  ENDDO

  du=du*0.5d0

  RETURN
END SUBROUTINE cin

!================================================================'
SUBROUTINE over(du)
  !*****Overlap matrix <B(jp,kp)|B(j,kp)>
  use nrtype
  USE prob_param
  USE bspline_module
  USE quad_module

  IMPLICIT NONE

  INTEGER(ITT) :: i,j,iq,s
  REAL(RT) :: du(2*kp-1,nmax),a0
  REAL(RT), ALLOCATABLE :: dut(:,:)

  !Diagonal elements
  du=0.d0
  DO j=2,n+1
     DO s=0,kp-1
        IF(tx(j+s).NE.tx(j+s+1))THEN
           a0=(tx(j+s+1)-tx(j+s))/2d0
           DO iq=1,nquad
              du(kp,j-1)=du(kp,j-1)+&
                   wd(iq)*bb(j,(j+s-kp)*nquad+iq)*bb(j,(j+s-kp)*nquad+iq)*a0
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
                 du(kp-i,j-1)=du(kp-i,j-1)+wd(iq)*bb(j-i,(j+s-kp)*nquad+iq)*&
                      bb(j,(j+s-kp)*nquad+iq)*a0
              ENDDO
           ENDIF
        ENDDO
        du(kp+i,j-i-1)=du(kp-i,j-1)
     ENDDO
  ENDDO

  RETURN 
END SUBROUTINE over

!================================================================'
SUBROUTINE steppot(du)
  !***** <B(jp,kp)|V|B(j,kp)>, with V=0, x>0 and V=-1, x<0
  use nrtype
  USE prob_param
  USE bspline_module
  USE quad_module

  IMPLICIT NONE

  INTEGER(ITT) :: i,j,iq,s
  REAL(RT) :: du(2*kp-1,nmax),a0
  REAL(RT), ALLOCATABLE :: dut(:,:)

  !Diagonal elements
  du=0.d0
  DO j=2,n+1
     DO s=0,kp-1
        IF(tx(j+s).NE.tx(j+s+1))THEN
           a0=(tx(j+s+1)-tx(j+s))/2d0
           DO iq=1,nquad
              IF(xd((j+s-kp)*nquad+iq).le.0.0_RT) then
		du(kp,j-1)=du(kp,j-1)+&
		   wd(iq)*bb(j,(j+s-kp)*nquad+iq)*bb(j,(j+s-kp)*nquad+iq)*a0
              end IF
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
                 IF(xd((j+s-kp)*nquad+iq).le.0.0_RT) then
                    du(kp-i,j-1)=du(kp-i,j-1)+wd(iq)*bb(j-i,(j+s-kp)*nquad+iq)*&
                         bb(j,(j+s-kp)*nquad+iq)*a0
		    if(j==11.and.i.eq.1)write(112,*)iq,xd((j+s-kp)*nquad+iq),du(kp-i,j-1)
                 end IF
              ENDDO
           ENDIF
        ENDDO
        du(kp+i,j-i-1)=du(kp-i,j-1)
     ENDDO
  ENDDO
  du=-du  !V=-1, x<0

  RETURN 
END SUBROUTINE steppot


!================================================================'
SUBROUTINE osc(du)
  !*****Harmonic oscillator 0.5*<B(jp,kp)|x²|B(j,kp)>
  use nrtype
  USE prob_param
  USE bspline_module
  USE quad_module

  IMPLICIT NONE

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
                   xd((j+s-kp)*nquad+iq)*xd((j+s-kp)*nquad+iq)*a0
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
                      xd((j+s-kp)*nquad+iq)*xd((j+s-kp)*nquad+iq)*a0
              ENDDO
           ENDIF
        ENDDO
        du(kp+i,j-i-1)=du(kp-i,j-1)
     ENDDO
  ENDDO
  du=du*0.5d0

  RETURN 
END SUBROUTINE osc


!================================================================'
SUBROUTINE mat_inter(du)
  !*****Dipole interaction (vel. gauge),nonsymmetric!!! mat_i(j,i)=-mat_i(i,j)
  !*****<B(jp,kp)|d/dx|B(j,kp)>
  use nrtype
  USE prob_param
  USE bspline_module
  USE quad_module

  IMPLICIT NONE

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
  use nrtype
  USE prob_param
  USE bspline_module
  USE quad_module

  IMPLICIT NONE

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

