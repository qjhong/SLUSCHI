program avg_std
implicit none

  double precision,allocatable,dimension(:)::array,Abar
  double precision avg,avgsqr,std,ste,avgAbar,avgsqrAbar

  integer i,j,n,m
  integer s,s_end,nB,t
  double precision P
  double precision,allocatable,dimension(:)::PHis,AvgHis,StdHis

  integer navg
  integer natom
  double precision AvgPHis,StdPHis

  open(7,file='for_avg')
  open(8,file='avg_std_detail.out')
  open(10,file='avg_std.out')

  n = 0
  do while (.true.)
    read(7,*,end=147,err=147)
    n = n + 1
  enddo
  147 continue
  rewind 7

  allocate(array(n),Abar(n))

  do i=1,n
    read(7,*) array(i)
  enddo
  close(7)

  open(9,file='natom')
  read(9,*) natom
  close(9)

  m = n
  s=1
  s_end=min(dble(200),dble(m)/50.d0)
  allocate(PHis(s_end),AvgHis(s_end),StdHis(s_end))
  do while (s<=s_end)
    avg = 0.d0
    avgsqr = 0.d0
    nB = floor(dble(m)/dble(s))
    Abar=0.d0
    t=0
    j=1
    do i=n-nB*s+1,n
      t=t+1
      if ( t==s ) then
        Abar(j) = Abar(j)/dble(s)
        t=0
        j=j+1
      endif
      Abar(j) = Abar(j) + array(i)
      avg = avg + array(i)
      avgsqr = avgsqr + array(i)**2.d0
    enddo
    avg = avg/dble(nB*s)
    avgsqr = avgsqr/dble(nB*s)
!    std = dsqrt(avgsqr-avg**2.d0)*dsqrt(dble(nB*s-1)/dble(nB*s))
!    ste = std/dsqrt(dble(nB*s))
    avgAbar=0.d0
    avgsqrAbar=0.d0
    do j=1,nB
      avgAbar = avgAbar + Abar(j)
      avgsqrAbar = avgsqrAbar + Abar(j)**2.d0
    enddo
    avgAbar = avgAbar/dble(nB)
    avgsqrAbar = avgsqrAbar/dble(nB)
    P = dble(s)*(avgsqrAbar-avgAbar**2.d0)/(avgsqr-avg**2.d0)
    write(8,99) m,s,nB,P,avg,sqrt( 2.d0*s/m * (avgsqr-avg**2.d0) )
    PHis(s)=P
    AvgHis(s)=avg
    StdHis(s)=sqrt( 2.d0*s/m * (avgsqr-avg**2.d0) )
    s = s+1
  enddo
  navg=max(s_end/10,10)
  do i=1,s_end-3 ! Average PHis, reduce oscillation
    PHis(s_end+1-i)=(PHis(s_end+1-i)+PHis(s_end-i)+PHis(s_end-1-i)+PHis(s_end-2-i))/4.d0
  enddo
  AvgPHis=0.d0
  StdPHis=0.d0
  do i=s_end-navg+1,s_end
    AvgPHis=AvgPHis+PHis(i)
    StdPHis=StdPHis+PHis(i)**2.d0
  enddo
  AvgPHis=AvgPHis/dble(navg)
  StdPHis=StdPHis/dble(navg)
  StdPHis=dsqrt(StdPHis-AvgPHis**2.d0)
  do i=1,s_end-10
    if (dabs(PHis(s_end+1-i)-AvgPHis)>StdPHis*2.d0 .and. dabs(PHis(s_end-i)-AvgPHis)>StdPHis*2.d0) exit
  enddo
  if (i/=1) i=i-1
  write(10,98) m,AvgPHis,StdPHis,i,PHis(s_end+1-i),AvgHis(s_end+1-i)/dble(natom),StdHis(s_end+1-i)/dble(natom)
  98 format(I5,2F12.5,I7,3F12.5)
  deallocate(PHis,AvgHis,StdHis)
  99 format(3I7,3F12.5)

  m = floor(dble(n)/2.d0)
  s=1
  s_end=min(dble(200),dble(m)/50.d0)
  allocate(PHis(s_end),AvgHis(s_end),StdHis(s_end))
  do while (s<=s_end)
    avg = 0.d0
    avgsqr = 0.d0
    nB = floor(dble(m)/dble(s))
    Abar=0.d0
    t=0
    j=1
    do i=n-nB*s+1,n
      t=t+1
      if ( t==s ) then
        Abar(j) = Abar(j)/dble(s)
        t=0
        j=j+1
      endif
      Abar(j) = Abar(j) + array(i)
      avg = avg + array(i)
      avgsqr = avgsqr + array(i)**2.d0
    enddo
    avg = avg/dble(nB*s)
    avgsqr = avgsqr/dble(nB*s)
!    std = dsqrt(avgsqr-avg**2.d0)*dsqrt(dble(nB*s-1)/dble(nB*s))
!    ste = std/dsqrt(dble(nB*s))
    avgAbar=0.d0
    avgsqrAbar=0.d0
    do j=1,nB
      avgAbar = avgAbar + Abar(j)
      avgsqrAbar = avgsqrAbar + Abar(j)**2.d0
    enddo
    avgAbar = avgAbar/dble(nB)
    avgsqrAbar = avgsqrAbar/dble(nB)
    P = dble(s)*(avgsqrAbar-avgAbar**2.d0)/(avgsqr-avg**2.d0)
    write(8,99) m,s,nB,P,avg,sqrt( 2.d0*s/m * (avgsqr-avg**2.d0) )
    PHis(s)=P
    AvgHis(s)=avg
    StdHis(s)=sqrt( 2.d0*s/m * (avgsqr-avg**2.d0) )
    s = s+1
  enddo
  navg=max(s_end/10,10)
  do i=1,s_end-3 ! Average PHis, reduce oscillation
    PHis(s_end+1-i)=(PHis(s_end+1-i)+PHis(s_end-i)+PHis(s_end-1-i)+PHis(s_end-2-i))/4.d0
  enddo
  AvgPHis=0.d0
  StdPHis=0.d0
  do i=s_end-navg+1,s_end
    AvgPHis=AvgPHis+PHis(i)
    StdPHis=StdPHis+PHis(i)**2.d0
  enddo
  AvgPHis=AvgPHis/dble(navg)
  StdPHis=StdPHis/dble(navg)
  StdPHis=dsqrt(StdPHis-AvgPHis**2.d0)
  do i=1,s_end-10
    if (dabs(PHis(s_end+1-i)-AvgPHis)>StdPHis*2.d0 .and. dabs(PHis(s_end-i)-AvgPHis)>StdPHis*2.d0) exit
  enddo
  if (i/=1) i=i-1
  write(10,98) m,AvgPHis,StdPHis,i,PHis(s_end+1-i),AvgHis(s_end+1-i)/dble(natom),StdHis(s_end+1-i)/dble(natom)
  deallocate(PHis,AvgHis,StdHis)

  m = floor(dble(n)/4.d0)
  s=1
  s_end=min(dble(200),dble(m)/50.d0)
  allocate(PHis(s_end),AvgHis(s_end),StdHis(s_end))
  do while (s<=s_end)
    avg = 0.d0
    avgsqr = 0.d0
    nB = floor(dble(m)/dble(s))
    Abar=0.d0
    t=0
    j=1
    do i=n-nB*s+1,n
      t=t+1
      if ( t==s ) then
        Abar(j) = Abar(j)/dble(s)
        t=0
        j=j+1
      endif
      Abar(j) = Abar(j) + array(i)
      avg = avg + array(i)
      avgsqr = avgsqr + array(i)**2.d0
    enddo
    avg = avg/dble(nB*s)
    avgsqr = avgsqr/dble(nB*s)
!    std = dsqrt(avgsqr-avg**2.d0)*dsqrt(dble(nB*s-1)/dble(nB*s))
!    ste = std/dsqrt(dble(nB*s))
    avgAbar=0.d0
    avgsqrAbar=0.d0
    do j=1,nB
      avgAbar = avgAbar + Abar(j)
      avgsqrAbar = avgsqrAbar + Abar(j)**2.d0
    enddo
    avgAbar = avgAbar/dble(nB)
    avgsqrAbar = avgsqrAbar/dble(nB)
    P = dble(s)*(avgsqrAbar-avgAbar**2.d0)/(avgsqr-avg**2.d0)
    write(8,99) m,s,nB,P,avg,sqrt( 2.d0*s/m * (avgsqr-avg**2.d0) )
    PHis(s)=P
    AvgHis(s)=avg
    StdHis(s)=sqrt( 2.d0*s/m * (avgsqr-avg**2.d0) )
    s = s+1
  enddo
  navg=max(s_end/10,10)
  do i=1,s_end-3 ! Average PHis, reduce oscillation
    PHis(s_end+1-i)=(PHis(s_end+1-i)+PHis(s_end-i)+PHis(s_end-1-i)+PHis(s_end-2-i))/4.d0
  enddo
  AvgPHis=0.d0
  StdPHis=0.d0
  do i=s_end-navg+1,s_end
    AvgPHis=AvgPHis+PHis(i)
    StdPHis=StdPHis+PHis(i)**2.d0
  enddo
  AvgPHis=AvgPHis/dble(navg)
  StdPHis=StdPHis/dble(navg)
  StdPHis=dsqrt(StdPHis-AvgPHis**2.d0)
  do i=1,s_end-10
    if (dabs(PHis(s_end+1-i)-AvgPHis)>StdPHis*2.d0 .and. dabs(PHis(s_end-i)-AvgPHis)>StdPHis*2.d0) exit
  enddo
  if (i/=1) i=i-1
  write(10,98) m,AvgPHis,StdPHis,i,PHis(s_end+1-i),AvgHis(s_end+1-i)/dble(natom),StdHis(s_end+1-i)/dble(natom)
  deallocate(PHis,AvgHis,StdHis)

  close(8)
end
