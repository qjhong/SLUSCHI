program status2out
implicit none

  double precision,allocatable,dimension(:)::mtemp
  integer,allocatable,dimension(:)::msure
  integer,allocatable,dimension(:)::msolid
  integer m

  integer iter,nmd

  double precision,allocatable,dimension(:)::temp
  integer,allocatable,dimension(:)::nsolid,nliquid,ntotal
  double precision,allocatable,dimension(:)::nsolid2,nliquid2
  integer ntemp

  double precision r,s
  integer i

  open(7,file='jobs_status.out')
  open(8,file='jobs_results.out2')

  read(7,*) 
  m = 0
  do while ( .true. )
    read(7,*,end=149,err=149) r
    m = m + 1
  enddo
  149 continue
  allocate(mtemp(m),msure(m),msolid(m))
  rewind 7
  read(7,*)

  i = 0
  do while ( .true. )
    i = i + 1
    read(7,*,end=150,err=150) mtemp(i),iter,nmd,msure(i),msolid(i)
  enddo
  150 continue
  close(7)

  ntemp = 1
  do i = 1,m-1
    if ( dabs(mtemp(i)-mtemp(i+1)) > .1d0 ) then
      ntemp = ntemp + 1
    endif
  enddo
  allocate(temp(ntemp),nsolid(ntemp),nliquid(ntemp),ntotal(ntemp))
  ntemp = 1
  nsolid = 0
  nliquid = 0
  ntotal = 0
  do i = 1,m
    if ( i>1 .and. dabs(mtemp(i)-mtemp(i-1)) > .1d0 ) then
      ntemp = ntemp + 1
    endif
    temp(ntemp) = mtemp(i)
    ntotal(ntemp) = ntotal(ntemp) + 1
    if (msolid(i) == 1 .or. msolid(i) == 11 ) nsolid(ntemp)=nsolid(ntemp)+1
    if (msolid(i) == 0 .or. msolid(i) == 10 ) nliquid(ntemp)=nliquid(ntemp)+1
  enddo

  allocate(nsolid2(ntemp),nliquid2(ntemp))
  do i= 1,ntemp
    nsolid2(i)=dble(nsolid(i)) + 1.d-4
    nliquid2(i)=dble(nliquid(i))
  enddo
  do i= 1,ntemp-1
    nsolid2(i)=nsolid2(i) + dble(nsolid(i+1))*(50.d0/(temp(i+1)-temp(i)))
    nliquid2(i)=nliquid2(i) + dble(nliquid(i+1))*(50.d0/(temp(i+1)-temp(i)))
  enddo
  do i= 2,ntemp
    nsolid2(i)=nsolid2(i) + dble(nsolid(i-1))*(50.d0/(temp(i)-temp(i-1)))
    nliquid2(i)=nliquid2(i) + dble(nliquid(i-1))*(50.d0/(temp(i)-temp(i-1)))
  enddo
  
  do i=1,ntemp
    r = nsolid2(i)/(nsolid2(i)+nliquid2(i)) * dble(ntotal(i))
    s = nliquid2(i)/(nsolid2(i)+nliquid2(i)) * dble(ntotal(i))
    write(8,*) nint(temp(i)),nint(r),ntotal(i)-nint(r),ntotal(i)
  enddo
  close(8)
end
