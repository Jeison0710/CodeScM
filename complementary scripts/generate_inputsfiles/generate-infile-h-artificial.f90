program generate_input

  implicit none

  double precision :: xmin,xmax
  character(33) :: outfile,format_string,intext,printall
  double precision :: sigma,x0,p0,g,dt0,ainit,aend,nshowini,nshowend,tfin
  integer :: selecter,Nran,kp,nb,nshow,v0,discoef,i,nquad
  write(*,*)'aqui'
  READ(*,*)          !'y'->stops before filling matrices
  READ(*,*)sigma           !1 ->save matrices 
  READ(*,*)x0
  READ(*,*)p0
  READ(*,*)g
  READ(*,*)selecter
  READ(*,*)intext 
  READ(*,*)Nran
  READ(*,*) 
  READ(*,*)kp 
  READ(*,*)nb 
  READ(*,*)xmin 
  READ(*,*)xmax
  READ(*,*)nquad 
  READ(*,*)tfin                
  READ(*,*)dt0
  READ(*,*)ainit
  READ(*,*)aend 
  READ(*,*) 
  READ(*,*)nshow 
  READ(*,*)nshowini 
  READ(*,*)nshowend
  READ(*,*)printall
  READ(*,*)
  READ(*,*)v0  
  READ(*,*)discoef

  do i=1,20
    ! calp=0.163d0+i*0.002d0
     write(outfile,100) 'initial',i-1,'.dat'
100 format(a7,i2.2,a4)
     open(4,file=outfile)
99   format(a7,i2.2,a25)
     WRITE(4,*)          !'y'->stops before filling matrices
     WRITE(4,*)sigma           !1 ->save matrices 
     WRITE(4,*)x0
     WRITE(4,*)p0
     WRITE(4,*)g
     WRITE(4,*)selecter
     WRITE(4,99)'seed_00',i-1,'_N_1048576_L_500prein.txt'
     WRITE(4,*)Nran
     WRITE(4,*) 
     WRITE(4,*)kp 
     WRITE(4,*)nb 
     WRITE(4,*)xmin 
     WRITE(4,*)xmax
     WRITE(4,*)nquad 
     WRITE(4,*)tfin                
     WRITE(4,*)dt0
     WRITE(4,*)ainit
     WRITE(4,*)aend 
     WRITE(4,*) 
     WRITE(4,*)nshow 
     WRITE(4,*)nshowini 
     WRITE(4,*)nshowend
     WRITE(4,*)printall
     WRITE(4,*)
     WRITE(4,*)v0
     WRITE(4,*)discoef          
     close(4)
  end do
end program generate_input
