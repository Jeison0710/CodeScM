program runer
   USE OMP_LIB
   integer :: n,m,i,mtot,thread_id
  ! INTEGER :: hdferr
  !! INTEGER :: layout
  ! INTEGER(HID_T)  :: file, space, dset, dcpl ! Handles
   CHARACTER(LEN=7) :: prefile
   CHARACTER(LEN=4) :: postfile
   CHARACTER(LEN=11) :: filename
   CHARACTER(LEN=2) :: x1
   CHARACTER(LEN=11) :: formatto
   DOUBLE PRECISION  :: timeunitary1,timeunitary2
   !CHARACTER(LEN=3) , PARAMETER :: dataset  = "DS1
   prefile='initial'
   postfile= '.dat'
   !REAL(RT) :: timeunitary1,timeunitary2,timeunitary3,aux1,aux2
   mtot=30
   90 FORMAT(i2.2)
   CALL tiempo(timeunitary1)
   !$OMP PARALLEL DO PRIVATE(thread_id,x1,filename) 
   do i=0,mtot-1
      write(x1,90) i
      thread_id= OMP_GET_THREAD_NUM()
      filename=prefile//trim(x1)//postfile
      write(*,*) './tt<'//filename
      write(*,*) thread_id
      call execute_command_line('./tt<'//filename)
   end do
   !$OMP END PARALLEL DO
   99 FORMAT(a37,f10.5,f13.5,f20.10)

   CALL tiempo(timeunitary2)
   timeunitary2=(timeunitary2-timeunitary1)/100.0d0
   WRITE(*,99)'time taken (h,m,s)=',&
        timeunitary2/3.6d3,timeunitary2/6.0d1,timeunitary2
end program

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
