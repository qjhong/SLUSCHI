program avg_pos
implicit none

  integer niter
  double precision,dimension(3,3)::latt
  double precision,dimension(3)::pos
  double precision,allocatable,dimension(:,:)::pos_history,pos_last
  integer natom
  double precision factor

  integer i,j,k
  integer s
  double precision r
  integer id,atom
  character str*40

  double precision,dimension(3)::pos_tmp,del_pos
  integer imin,jmin,kmin
  double precision min_norm,norm

  open(7,file='OUTCAR')
  open(8,file='POSCAR')
  open(9,file='strct')
  open(10,file='natom')

  read(8,*) str
  read(8,*) factor
  read(8,*) pos
  latt(:,1) = pos*factor
  read(8,*) pos
  latt(:,2) = pos*factor
  read(8,*) pos
  latt(:,3) = pos*factor
  read(10,*) natom
  
  allocate(pos_history(3,natom))
  allocate(pos_last(3,natom))

  niter = 0
  pos_history = 0.d0
  do while (.true.)
    read(7,*,err=147,end=147) str
    if (str(1:8)=="POSITION") then
      niter=niter+1
      read(7,*) str
      do s=1,natom
        read(7,*) pos
        if (niter == 1) then
          pos_last(:,s) = pos
        else
          min_norm = 100.d0
          do i=-1,1
          do j=-1,1
          do k=-1,1
            pos_tmp = pos + latt(:,1)*dble(i) + latt(:,2)*dble(j) + latt(:,3)*dble(k)
            del_pos = pos_tmp - pos_last(:,s)
            norm = dsqrt(del_pos(1)**2.d0 + del_pos(2)**2.d0 + del_pos(3)**2.d0)
            if ( norm < min_norm ) then
              min_norm = norm
              imin = i
              jmin = j
              kmin = k
            endif
          enddo
          enddo
          enddo
          pos = pos + latt(:,1)*dble(imin) + latt(:,2)*dble(jmin) + latt(:,3)*dble(kmin)
          pos_last(:,s) = pos
        endif
        pos_history(:,s) = pos_history(:,s) + pos(1:3)
      enddo
    endif
  enddo
  147 continue

  close(7)
  close(8)
  close(10)

  pos_history = pos_history / niter

  write(9,99) latt(:,1)
  write(9,99) latt(:,2)
  write(9,99) latt(:,3)
  do i=1,natom
    write(9,99) pos_history(:,i)
    99 format(3F12.5)
  enddo
  close(9)

end
