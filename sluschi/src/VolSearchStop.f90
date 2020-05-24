program VolSearchStop
implicit none

  integer i
  integer nwork
  double precision,dimension(6)::pressure3_total
  double precision,dimension(15)::tmp
  integer stopflag
  double precision cnvg

  ! read pressure3.out; average
  open(7,file='lattice_pressure_history.out')
  open(8,file='thmexp.in')
  read(8,*) cnvg
  close(8)
  nwork=0
  do while (.true.)
    read(7,*,err=147,end=147)
    nwork = nwork+1
  enddo
  147 continue
  rewind 7
  do i=1,nwork
    read(7,*) tmp
  enddo
  close(7)
  pressure3_total(:) = tmp(10:15)

  stopflag = 1
  do i=1,6
    if (dabs(pressure3_total(i)) > cnvg) stopflag=0
  enddo

  open(9,file='VolSearchStop.out')
  write(9,"(I4)") stopflag
  close(9)

end

