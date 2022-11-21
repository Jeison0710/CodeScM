MODULE prob_param
  USE nrtype
  INTEGER(ITT) ::  kp,nb,nquad,ang_mom,npoints,pc,Nran,ran
  INTEGER(ITT) ::  ntot,nlarg,nlard,nt,nt1,ktot
  INTEGER(ITT) ::  nshow,V0_switch,discoef
  REAL(RT) :: x0,sigma,p0,V0,g,aini,aend,da,tprinti,tprinte
  REAL(RT) :: xmin0,xmin,xmax
  REAL(RT) :: dt0,tfin,dtada
  complex(CT), allocatable :: psi0(:)
  real(RT), allocatable :: psi0r(:),psi0i(:),cpsi0r(:),cpsi0i(:),cpsi0raux(:),cpsi0iaux(:)
  real(RT), allocatable :: xsr(:),xsi(:),aat(:)
  character(1) :: print_all
  character(20) :: initialtext
  !FILE PARAMETERS
  INTEGER(ITT), PARAMETER :: nfileoutput=20 
END MODULE prob_param

MODULE bspline_module
  USE nrtype
  !external indexb
  !INTEGER(ITT) ::indexb
  INTEGER(ITT) :: nmax,n,l,l0!,indexb
  INTEGER(ITT), ALLOCATABLE :: index_xd(:)
  REAL(RT) :: dx,bint
  REAL(RT), ALLOCATABLE :: tx(:),bb(:,:),dbb(:,:),bba(:),dbba(:)
  REAL(RT), ALLOCATABLE :: kinetic(:),overlap(:),bb_boundary(:,:),dbb_boundary(:,:)
  REAL(RT), ALLOCATABLE :: bijkl(:,:,:,:),bijk(:,:,:)
  


END MODULE bspline_module



MODULE quad_module
  USE nrtype
  INTEGER(ITT) :: nq,nqk,nq0
  REAL(RT),ALLOCATABLE :: xd(:),wd(:)
END MODULE quad_module

module filenames
  USE nrtype
  CHARACTER(100) :: filenorm
  CHARACTER(100) :: filedistr,fileoutput
end module filenames

!************  modules for backward substitution ************
 MODULE detax
   USE nrtype
   IMPLICIT NONE
   REAL(RT), DIMENSION(:,:), ALLOCATABLE :: ar,ai
   INTEGER(ITT), DIMENSION(:,:), ALLOCATABLE :: bitmtx
   INTEGER(ITT) :: nnllmmaaxx, nnggrr
 END MODULE detax
 
 MODULE varindex
   USE nrtype
   IMPLICIT NONE
   INTEGER(ITT), DIMENSION(:), ALLOCATAbLE :: ntl
 END MODULE varindex

 MODULE varindexb
   USE nrtype
   IMPLICIT NONE
   INTEGER(ITT), DIMENSION(:), ALLOCATABLE::ntlb
 end module varindexb

 module detbx
   USE nrtype
   IMPLICIT NONE
   REAL(RT), DIMENSION(:,:), ALLOCATABLE :: d,dxx,dxx1,dxx2,dxx0
   REAL(RT):: theta
 end module detbx

 module detbxcx
   use nrtype
   real(RT),dimension(:,:),allocatable:: dr, di
 end module detbxcx

!************  modules for periodic boundary conditions ************
 module mtx_pbc
   use nrtype
   complex(CT),dimension(:,:),allocatable:: m0p,m0m,m1p,m1m,m2p,m2m
   complex(CT),dimension(:,:),allocatable:: m0aux,m1aux,m2aux
   real(RT),dimension(:,:),allocatable:: mauxb,mauxh
   real(RT),dimension(:,:),allocatable:: b0,b1,b2,aa0,aa1,aa2,psi2,pot
 end module mtx_pbc



