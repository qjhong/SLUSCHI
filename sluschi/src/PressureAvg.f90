program DetermineSize
implicit none

  integer i
  integer nwork
  double precision,allocatable,dimension(:,:)::pressure3
  double precision,dimension(6)::pressure3_avg !Pxx Pyy Pzz Pxy Pyz Pzx
  double precision pressure_kinetic,pressure_Pulay,pressure_target
  double precision,dimension(6)::pressure3_total !Pxx Pyy Pzz Pxy Pyz Pzx

  ! read pressure3.out; average
  open(7,file='pressure3.out')
  nwork=0
  do while (.true.)
    read(7,*,err=147,end=147)
    nwork = nwork+1
  enddo
  147 continue
  rewind 7
  allocate(pressure3(nwork,6))
  do i=1,nwork
    read(7,*) pressure3(i,1:6)
  enddo
  close(7)
  do i=2,nwork
    pressure3(i,1:6)=pressure3(i,1:6)+pressure3(i-1,1:6)
  enddo
  pressure3_avg = pressure3(nwork,1:6)/nwork

  ! read kinetic pressure, Pulay stress, and target pressure; calculate total pressure
  open(8,file='pressure_kinetic.out')
  read(8,*) pressure_kinetic
  close(8)
  open(9,file='pressure_Pulay.out')
  read(9,*) pressure_Pulay
  close(9)
  open(10,file='pressure_target.out')
  read(10,*) pressure_target
  close(10)
  pressure3_total = pressure3_avg
  pressure3_total(1:3) = pressure3_total(1:3) + pressure_kinetic + pressure_Pulay - pressure_target

  ! output pressure3_total
  open(11,file='pressure3_total.out')
  write(11,"(6F15.6)") pressure3_total
  close(11)

end
