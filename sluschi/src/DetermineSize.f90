program DetermineSize
implicit none

  integer i,j,k
  integer nwork
  double precision,allocatable,dimension(:,:)::pressure3_total
  double precision,dimension(15)::tmp
  double precision,dimension(3,3)::latt,latt_new
  double precision,dimension(3,3)::scal

  double precision factor_divide
  integer n_avg


  ! read pressure3.out; average
  open(7,file='lattice_pressure_history.out')
  nwork=0
  do while (.true.)
    read(7,*,err=147,end=147)
    nwork = nwork+1
  enddo
  147 continue
  rewind 7
  allocate(pressure3_total(nwork,6))
  do i=1,nwork
    read(7,*) tmp
    pressure3_total(i,:) = tmp(10:15)
  enddo
  latt(1,:)= tmp(1:3)
  latt(2,:)= tmp(4:6)
  latt(3,:)= tmp(7:9)
  close(7)

  open(8,file='param.in')
  read(8,*) n_avg
  read(8,*) factor_divide
  close(8)

  scal(1,1) = avg(pressure3_total(:,1),nwork,n_avg)
  scal(2,2) = avg(pressure3_total(:,2),nwork,n_avg)
  scal(3,3) = avg(pressure3_total(:,3),nwork,n_avg)
  scal(1,2) = avg(pressure3_total(:,4),nwork,n_avg)
  scal(2,1) = scal(1,2)
  scal(2,3) = avg(pressure3_total(:,5),nwork,n_avg)
  scal(3,2) = scal(2,3)
  scal(1,3) = avg(pressure3_total(:,6),nwork,n_avg)
  scal(3,1) = scal(1,3)

  scal = scal/factor_divide
  scal(1,1) = scal(1,1) + 1.d0
  scal(2,2) = scal(2,2) + 1.d0
  scal(3,3) = scal(3,3) + 1.d0

  latt_new = 0.d0
  do i=1,3
  do j=1,3
    do k=1,3
      latt_new(i,j) = latt_new(i,j) + latt(i,k)*scal(k,j)
    enddo
  enddo
  enddo

  open(9,file='lattice_predict.out')
  do i=1,3
    write(9,"(3F15.6)") latt_new(i,:)
  enddo
  close(9)

contains
double precision function avg(x,ndim,navg)
implicit none

  integer ndim
  double precision,dimension(ndim)::x
  integer i
  integer navg

  do i=1,navg
    if (ndim-navg+i>0) then
      avg = avg + x(ndim-navg+i)
    endif
  enddo

  avg = avg/min(navg,ndim)

end function

end

