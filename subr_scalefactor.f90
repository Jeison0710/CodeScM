!*******************************************************************
! Time and scale factor
!*******************************************************************

SUBROUTINE scalefactor
!======================================================================
!Find the scale factor and saves it in aat
!it takes into account the time step change which happend at time
!tc where a(tc)=0.1  for an initial condition 
!aini=9.90099009900990111d-3 and aend=1.0d0
!the step change happents at nt1
!----------------------------------------------------------------------
!=======================================================================
  USE nrtype
  USE prob_param

  implicit none


  REAL(RT) :: sigm0,aux,a3,f,t,tc,aa !tc is the time in which precision is increased
  INTEGER(ITT) :: i,flag,nm,nt2
  REAL(RT), ALLOCATABLE :: auxat(:)
  integer, parameter :: out_unit=20
  
  !external f
  dtada=dt0*1.0d-1  
  nm=100000!0000000 
  da=dble((aend-aini)/(nm)) 
  sigm0=0.3_RT!43382943_RT
  !Saves t(a) in array aat for grid a={aini,..aend} of length nm
  ALLOCATE (aat(nm+2))
  aat=ZERO
  WRITE(*,*) da
  flag=0
  Do i=1,nm+1
     aux=3.0_RT*sigm0
     a3=((aini+dble(i-1)*da)+ HALF*da)!!!!!!!!!!i-1  +1/2
     a3=a3*a3*a3
     aux=aux/(2.0_RT*a3*(sigm0+(ONE-sigm0)*a3))
     aux=DSQRT(aux)
     
     aat(i+1)=aat(i)+da*aux
     
     IF ((da*i+aini .GE. 0.1d0).and.(flag==0)) THEN
         tc=aat(i+1)
         write (*,*) 'step change time is',tc, i*da+aini
         flag=1
      !   stop
     END IF
  END DO
  flag=0
  WRITE(*,*) aat(nm+1), aini+(nm+1-1)*da
  WRITE(*,*) aat(nm+2), aini+(nm+1)*da

  tc=16.869303654267441_RT      ! 0.10000009009909999    


  !Saves a(t) in aat for t=0,... t(aend) with dt0 step size
   
  nt=int((aat(nm))/dt0)
  nt1=int(tc/dt0)
  nt2=int((aat(nm)-tc)/dtada)
  nt=nt1+nt2

  write(*,*) '--------------------------------------'
  write(*,*) nt1*dt0+nt2*dtada,aat(nm)

  WRITE(*,*) nt,(aat(nm+1)),(aini+(nm-1)*da)
  open (unit=out_unit,file="Scalefactor.txt",action="write",status="replace")
  !DO i=1,nm+1
  !   write (out_unit,*) aat(i)
  !end DO 
  !stop
  ALLOCATE (auxat(nt+1))
  auxat=zero
  auxat(1)=aini
  Do i= 1 , nt1 !calculate a(t) until tc with dt=dt0 and save it in auxat
     
     t=dble(dt0*(i))
     call bisection(t,aini,aend,1.0d-9,aux,flag)
     auxat(i+1)=aux
     aux=zero
     write (out_unit,*) t,auxat(i+1)
  END DO
  Do i= nt1+1 , nt !calculate a(t) from tc to the end with dt=dt0/10 and save it in auxat
     
     t=dble(dt0*(nt1)) +dble(dtada*(i-nt1))
     call bisection(t,aini,aend,1.0d-9,aux,flag)
     auxat(i+1)=aux
     aux=zero
     write (out_unit,*) t,auxat(i+1)
  END DO
  !stop
  WRITE(*,*)  nt,nt1,nt2,nt1*dt0+nt2*dtada,auxat(nt+1)
 
  close (out_unit)
  !stop
  !saves auxat in aat
  DEALLOCATE (aat)
  ALLOCATE (aat(nt+1))
  aat=auxat 
  DEALLOCATE (auxat)
 ! WRITE(*,*) nt*
  !*******************************
   

  !********************
 

  

END SUBROUTINE scalefactor

!*******************************************************************
! linear interpolation function for scale factor
!*******************************************************************
  subroutine f(x)
!======================================================================
! linear interpolation function for scale factor
! it has aat() as fixed points
!aat() was defined in SUBROUTINE scalefactor
!----------------------------------------------------------------------
! input..
! x     	   = input to evaluate (real)

! Hidden input:	   
! aat              = array of values t(ai) dim(nmax)(real) 
! aini,aend,da	   = initial,final value and distance of the grid (real)

! output:
! f	  	   = value of the function(real)
!=======================================================================
  USE nrtype
  USE prob_param
  
  IMPLICIT NONE

  REAL(RT) :: x
  INTEGER(ITT) :: k
 

  IF (x .LT. aini) THEN
     x=aat(1)
  ELSE IF (x .GT. aend) THEN
     x=aat(1000000+1) !!MUST BE EQUAL TO nm
  ELSE
    k=FLOOR((x-aini)/da)+ 1
    x=((x-(aini+(k-1)*da))*aat(k+1)+(-x+(aini+(k)*da))*aat(k))/da
  END IF
  return
  !f=((x-(aini+(k-1)*da))*aat(k+1)+(-x+(aini+(k)*da))*aat(k))/da
  end subroutine f

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Subroutine bisection(fo,x1,x2,eps,Root,flag)
!============================================================
! Solutions of equation f(x)=0 on [x1,x2] interval
! Method: Bisectional (closed domain) (a single root)
! Alex G. January 2010
!------------------------------------------------------------
! input ...
! f   - function - evaluates f(x) for any x in [x1,x2]
! x1  - left endpoint of initial interval
! x2  - right endpoint of initial interval
! eps - desired uncertainity of the root as |b-a|<eps
! output ...
! Root  - root of the equation f(x)=0
! flag  - indicator of success
!         >0 - a single root found, flag=number of iterations
!          0 - no solutions for the bisectional method
!         <0 - not a root but singularity, flag=number of iterations
!
! Comments: Function f(x) has to change sign between x1 and x2
!           Max number of iterations - 200 (accuracy (b-a)/2**200)
!====================================================================
use nrtype
implicit none

REAL(RT) ::  x1, x2, eps, Root,fo,auxa,auxc
REAL(RT) :: a, b, c
INTEGER(ITT) i, flag
INTEGER(ITT), parameter:: iter=20000
auxa=x1
auxc=x2
call f(auxa)
call f(auxc)
!* check the bisection condition
if((auxa-fo)*(auxc-fo).gt.zero) then
  flag = 0
  return
end if

!* initialize calculations
a=x1
b=x2

!* Iterative refining the solution 
do i=1,iter
  c=(b+a)*HALF
  auxa=a
  auxc=c
  call f(auxa)
  call f(auxc)
  if((auxa-fo)*(auxc-fo).le.zero) then
      b = c
    else
      a=c
  end if
! condition(s) to stop iterations)
  if(ABS(b-a).le. eps) exit  
end do
Root=(b+a)*HALF

!* check if it is a root or singularity
auxa=root
call f(auxa)
if (ABS(auxa-fo) .le. ONE) then
  flag=i
  else
  flag = -i
end if
if(flag == 0) then 
  write(*,104)
  stop
end if
if(flag < 0) then
    write(*,103) Root
end if

104 format(' no roots for Bisectional method')
103 format(' Singularity = ',1pe12.5)

end subroutine bisection
