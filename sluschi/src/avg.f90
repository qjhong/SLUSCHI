program avg
implicit none

  integer i
  integer n
  double precision navg
  double precision,allocatable,dimension(:,:)::sum_pos
  double precision,dimension(3)::pos
  double precision,dimension(3,3)::latt
  double precision,dimension(3)::pos_tmp,del_pos
  double precision norm_min,norm
  integer k1,k2,k3,k1_min,k2_min,k3_min

  open(7,file='struct1')
  n = 0
  do while (.true.)
    read(7,*,end=147,err=147)
    n = n + 1
  enddo
  147 continue
  close(7)
  n = n - 3

  allocate(sum_pos(3,n))
  open(7,file='struct1')
  do i=1,3
    read(7,*) latt(:,i)
  enddo
  do i=1,n
    read(7,*) pos
    sum_pos(1:3,i) = pos
  enddo
  close(7)

  navg = 1.d0
  open(7,file='struct2')
  do i=1,3
    read(7,*) latt(:,i)
  enddo
  do i=1,n
    read(7,*) pos
    norm_min = 100.d0
    do k1=-1,1
    do k2=-1,1
    do k3=-1,1
      pos_tmp = pos + latt(:,1)*dble(k1) + latt(:,2)*dble(k2) + latt(:,3)*dble(k3)
      del_pos = pos_tmp - sum_pos(:,i)/navg
      norm = dsqrt(del_pos(1)**2.d0 + del_pos(2)**2.d0 + del_pos(3)**2.d0)
      if (norm < norm_min) then
        norm_min = norm
        k1_min = k1
        k2_min = k2
        k3_min = k3
      endif
    enddo
    enddo
    enddo
    pos = pos + latt(:,1)*dble(k1_min) + latt(:,2)*dble(k2_min) + latt(:,3)*dble(k3_min)
    sum_pos(:,i) = sum_pos(:,i)+pos
  enddo
  close(7)

  navg = navg + 1
  open(7,file='struct3')
  do i=1,3
    read(7,*) latt(:,i)
  enddo
  do i=1,n
    read(7,*) pos
    norm_min = 100.d0
    do k1=-1,1
    do k2=-1,1
    do k3=-1,1
      pos_tmp = pos + latt(:,1)*dble(k1) + latt(:,2)*dble(k2) + latt(:,3)*dble(k3)
      del_pos = pos_tmp - sum_pos(:,i)/navg
      norm = dsqrt(del_pos(1)**2.d0 + del_pos(2)**2.d0 + del_pos(3)**2.d0)
      if (norm < norm_min) then
        norm_min = norm
        k1_min = k1
        k2_min = k2
        k3_min = k3
      endif
    enddo
    enddo
    enddo
    pos = pos + latt(:,1)*dble(k1_min) + latt(:,2)*dble(k2_min) + latt(:,3)*dble(k3_min)
    sum_pos(:,i) = sum_pos(:,i)+pos
  enddo
  close(7)

  navg = navg + 1
  open(7,file='struct4')
  do i=1,3
    read(7,*) latt(:,i)
  enddo
  do i=1,n
    read(7,*) pos
    norm_min = 100.d0
    do k1=-1,1
    do k2=-1,1
    do k3=-1,1
      pos_tmp = pos + latt(:,1)*dble(k1) + latt(:,2)*dble(k2) + latt(:,3)*dble(k3)
      del_pos = pos_tmp - sum_pos(:,i)/navg
      norm = dsqrt(del_pos(1)**2.d0 + del_pos(2)**2.d0 + del_pos(3)**2.d0)
      if (norm < norm_min) then
        norm_min = norm
        k1_min = k1
        k2_min = k2
        k3_min = k3
      endif
    enddo
    enddo
    enddo
    pos = pos + latt(:,1)*dble(k1_min) + latt(:,2)*dble(k2_min) + latt(:,3)*dble(k3_min)
    sum_pos(:,i) = sum_pos(:,i)+pos
  enddo
  close(7)

  navg = navg + 1
  open(7,file='struct5')
  do i=1,3
    read(7,*) latt(:,i)
  enddo
  do i=1,n
    read(7,*) pos
    norm_min = 100.d0
    do k1=-1,1
    do k2=-1,1
    do k3=-1,1
      pos_tmp = pos + latt(:,1)*dble(k1) + latt(:,2)*dble(k2) + latt(:,3)*dble(k3)
      del_pos = pos_tmp - sum_pos(:,i)/navg
      norm = dsqrt(del_pos(1)**2.d0 + del_pos(2)**2.d0 + del_pos(3)**2.d0)
      if (norm < norm_min) then
        norm_min = norm
        k1_min = k1
        k2_min = k2
        k3_min = k3
      endif
    enddo
    enddo
    enddo
    pos = pos + latt(:,1)*dble(k1_min) + latt(:,2)*dble(k2_min) + latt(:,3)*dble(k3_min)
    sum_pos(:,i) = sum_pos(:,i)+pos
  enddo
  close(7)

  navg = navg + 1
  open(7,file='struct6')
  do i=1,3
    read(7,*) latt(:,i)
  enddo
  do i=1,n
    read(7,*) pos
    norm_min = 100.d0
    do k1=-1,1
    do k2=-1,1
    do k3=-1,1
      pos_tmp = pos + latt(:,1)*dble(k1) + latt(:,2)*dble(k2) + latt(:,3)*dble(k3)
      del_pos = pos_tmp - sum_pos(:,i)/navg
      norm = dsqrt(del_pos(1)**2.d0 + del_pos(2)**2.d0 + del_pos(3)**2.d0)
      if (norm < norm_min) then
        norm_min = norm
        k1_min = k1
        k2_min = k2
        k3_min = k3
      endif
    enddo
    enddo
    enddo
    pos = pos + latt(:,1)*dble(k1_min) + latt(:,2)*dble(k2_min) + latt(:,3)*dble(k3_min)
    sum_pos(:,i) = sum_pos(:,i)+pos
  enddo
  close(7)

  navg = navg + 1
  open(7,file='struct7')
  do i=1,3
    read(7,*) latt(:,i)
  enddo
  do i=1,n
    read(7,*) pos
    norm_min = 100.d0
    do k1=-1,1
    do k2=-1,1
    do k3=-1,1
      pos_tmp = pos + latt(:,1)*dble(k1) + latt(:,2)*dble(k2) + latt(:,3)*dble(k3)
      del_pos = pos_tmp - sum_pos(:,i)/navg
      norm = dsqrt(del_pos(1)**2.d0 + del_pos(2)**2.d0 + del_pos(3)**2.d0)
      if (norm < norm_min) then
        norm_min = norm
        k1_min = k1
        k2_min = k2
        k3_min = k3
      endif
    enddo
    enddo
    enddo
    pos = pos + latt(:,1)*dble(k1_min) + latt(:,2)*dble(k2_min) + latt(:,3)*dble(k3_min)
    sum_pos(:,i) = sum_pos(:,i)+pos
  enddo
  close(7)

  navg = navg + 1
  open(7,file='struct8')
  do i=1,3
    read(7,*) latt(:,i)
  enddo
  do i=1,n
    read(7,*) pos
    norm_min = 100.d0
    do k1=-1,1
    do k2=-1,1
    do k3=-1,1
      pos_tmp = pos + latt(:,1)*dble(k1) + latt(:,2)*dble(k2) + latt(:,3)*dble(k3)
      del_pos = pos_tmp - sum_pos(:,i)/navg
      norm = dsqrt(del_pos(1)**2.d0 + del_pos(2)**2.d0 + del_pos(3)**2.d0)
      if (norm < norm_min) then
        norm_min = norm
        k1_min = k1
        k2_min = k2
        k3_min = k3
      endif
    enddo
    enddo
    enddo
    pos = pos + latt(:,1)*dble(k1_min) + latt(:,2)*dble(k2_min) + latt(:,3)*dble(k3_min)
    sum_pos(:,i) = sum_pos(:,i)+pos
  enddo
  close(7)

  navg = navg + 1
  open(7,file='struct9')
  do i=1,3
    read(7,*) latt(:,i)
  enddo
  do i=1,n
    read(7,*) pos
    norm_min = 100.d0
    do k1=-1,1
    do k2=-1,1
    do k3=-1,1
      pos_tmp = pos + latt(:,1)*dble(k1) + latt(:,2)*dble(k2) + latt(:,3)*dble(k3)
      del_pos = pos_tmp - sum_pos(:,i)/navg
      norm = dsqrt(del_pos(1)**2.d0 + del_pos(2)**2.d0 + del_pos(3)**2.d0)
      if (norm < norm_min) then
        norm_min = norm
        k1_min = k1
        k2_min = k2
        k3_min = k3
      endif
    enddo
    enddo
    enddo
    pos = pos + latt(:,1)*dble(k1_min) + latt(:,2)*dble(k2_min) + latt(:,3)*dble(k3_min)
    sum_pos(:,i) = sum_pos(:,i)+pos
  enddo
  close(7)

  navg = navg + 1
  open(7,file='struct10')
  do i=1,3
    read(7,*) latt(:,i)
  enddo
  do i=1,n
    read(7,*) pos
    norm_min = 100.d0
    do k1=-1,1
    do k2=-1,1
    do k3=-1,1
      pos_tmp = pos + latt(:,1)*dble(k1) + latt(:,2)*dble(k2) + latt(:,3)*dble(k3)
      del_pos = pos_tmp - sum_pos(:,i)/navg
      norm = dsqrt(del_pos(1)**2.d0 + del_pos(2)**2.d0 + del_pos(3)**2.d0)
      if (norm < norm_min) then
        norm_min = norm
        k1_min = k1
        k2_min = k2
        k3_min = k3
      endif
    enddo
    enddo
    enddo
    pos = pos + latt(:,1)*dble(k1_min) + latt(:,2)*dble(k2_min) + latt(:,3)*dble(k3_min)
    sum_pos(:,i) = sum_pos(:,i)+pos
  enddo
  close(7)

  sum_pos = sum_pos/10
  open(11,file='avg')
  write(11,99) latt(:,1)
  write(11,99) latt(:,2)
  write(11,99) latt(:,3)
  do i=1,n
    write(11,99) sum_pos(:,i)
  enddo
  99 format(3F12.5)
  close(11)

end
