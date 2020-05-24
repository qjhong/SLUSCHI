program MPFit
implicit none

  integer,parameter::n_interval=1001,n=n_interval-1
  integer,parameter::n_randomnumber=100000
  double precision,parameter:: eV=1.602d-19,kB=1.3806d-23,F_converge=160.d0

  integer,allocatable,dimension(:)::nsolid,nliquid,nsum
  double precision,allocatable,dimension(:)::temp
  double precision,allocatable,dimension(:,:)::F
  double precision,dimension(n_interval)::Tm,gam
  double precision,dimension(n_interval,n_interval)::P
  double precision Tm_max,Tm_min,gam_max,gam_min
  integer nwork
  double precision s
  integer i,j,k,l
  double precision T,delG,fexp,pr,factor
  double precision Pmax
  integer imax,jmax
  integer m
  integer iup,ilow,jup,jlow
  double precision,dimension(n_randomnumber)::T_rn
  double precision sum_P,sum_Q
  integer ncount

  open(7,file='jobs_results.out')
  nwork=0
  do while (.true.)
    read(7,*,err=147,end=147)
    nwork = nwork+1
  enddo
  147 continue
  rewind 7
  allocate(temp(nwork),nsolid(nwork),nliquid(nwork),nsum(nwork))
  do i=1,nwork
    read(7,*) temp(i),nsolid(i),nliquid(i)
    nsum(i) = nsolid(i)+nliquid(i)
  enddo
  close(7)
  ! To perform a fit, jobs_results.out has to satisfy
  ! 1. there are at least two T' with >4 duplictates. T': with fractional distribution (not 0 or 1)
  ! 2. each neighbor of T' with >4 duplicates
  open(8,file='FlagFitSuccess.out')
  open(9,file='MPFit.log')
  open(10,file='MPFit.out')
  ncount = 0
  if (nliquid(1)==0 .and. nsolid(nwork)==0) then
    do i=2,nwork-1
      if (nsolid(i)>0 .and. nliquid(i)>0 .and. nsum(i)>=4 .and. nsum(i-1)>=4 .and. nsum(i+1)>=4 ) ncount = ncount+1
    enddo
  endif
  if (ncount < 2) then
    write(8,*) '0'
  else
    write(8,*) '1'
  Tm_min = 100000.d0
  Tm_max = 0.d0
  do i=1,nwork
    if (temp(i)<Tm_min) Tm_min=temp(i)
    if (temp(i)>Tm_max) Tm_max=temp(i)
  enddo
  gam_min = 1.d0
  gam_max = 100.d0
  
  allocate(F(n_interval,nwork))
  do i=1,nwork
    F(:,i) = prob(nsolid(i),nliquid(i),n_interval)
  enddo

  do m=1,5

    P = 1.d0
    s = (Tm_max-Tm_min)/dble(n_interval-1)
    do i = 1,n_interval
      Tm(i) = Tm_min + s * dble(i-1)
    enddo
    s = (gam_max/gam_min)**(1.d0/dble(n_interval-1))
    do i = 1,n_interval
      gam(i) = gam_min * s**(i-1)
    enddo

    do j=1,n_interval
    do i=1,n_interval
      do k=1,nwork
        T = temp(k)
        delG = (Tm(i)-T)/Tm(i)
        if ( - delG * eV / (2.d0*kB*T) * gam(j) < 100.d0 ) then
          fexp = dexp( - delG * eV / (2.d0*kB*T) * gam(j) )
          pr = fexp / (1.d0 + fexp)
        else 
          fexp = dexp( delG * eV / (2.d0*kB*T) * gam(j) )
          pr = 1.d0 / (1.d0 + fexp)
        endif
        l = floor(pr*n) + 1
        if (l==n+1) l=n
        factor = F(l,k)*(dble(l)-pr*dble(n)) + F(l+1,k)*(pr*dble(n)-dble(l)+1.d0)
        P(i,j) = P(i,j) * factor;
      enddo
    enddo
    enddo
   
    ! find max P
    Pmax = 0.d0
    do i=1,n_interval
    do j=1,n_interval
      if (P(i,j)>Pmax) then
        Pmax = P(i,j)
        imax = i
        jmax = j
      endif
    enddo
    enddo
    ! locate boundary
    iup=0
    ilow=n_interval
    jup=0
    jlow=n_interval
    do i=1,n_interval
    do j=1,n_interval
      if (P(i,j)>PMax/F_converge) then
        if (i>iup) iup=i
        if (i<ilow) ilow=i
        if (j>jup) jup=j
        if (j<jlow) jlow=j
      endif
    enddo
    enddo
    if (iup==n_interval) then
      Tm_max = 2.d0*Tm(n_interval) - Tm(1)
    else
      Tm_max = Tm(iup)+10.d0
    endif
    if (ilow==1) then
      Tm_min = max(2.d0*Tm(1)-Tm(n_interval),0.d0)
    else
      Tm_min = Tm(ilow)-10.d0
    endif
    if (jup==n_interval) then
      gam_max = min(gam(n_interval)*gam(n_interval)/gam(1),1000.d0)
    else
      gam_max = gam(jup)*1.1
    endif
    if (jlow==1) then
      gam_min = gam(1)*gam(1)/gam(n_interval)
    else
      gam_min = gam(jlow)*0.9
    endif
    write(9,'(4F20.5)') Tm_min,Tm_max,gam_min,gam_max

  enddo

  sum_P=0.d0
  do i=1,n_interval
  do j=1,n_interval
    sum_P = sum_P + P(i,j)
  enddo
  enddo
  P=P/sum_P
  do i=1,n_interval
  do j=1,n_interval-1
    P(i,n_interval) = P(i,n_interval) + P(i,j)
  enddo
  enddo
  
  do k=1,n_randomnumber
    call random_number(s)
    sum_P = 0.d0
    do i=1,n_interval
      sum_P = sum_P + P(i,n_interval)
      if (sum_P > s) goto 148
    enddo
    148 continue
    T_rn(k) = Tm(i)
  enddo
  sum_P=0.d0
  sum_Q=0.d0
  do k=1,n_randomnumber
    sum_P = sum_P + T_rn(k)
    sum_Q = sum_Q + T_rn(k)**2.d0
  enddo
  sum_P = sum_P/n_randomnumber
  sum_Q = sum_Q/n_randomnumber
  if ( dsqrt(sum_Q-sum_P**2.d0) < 500.d0 ) then
    write(10,*) sum_P,dsqrt(sum_Q-sum_P**2.d0)
  else
    write(10,*) Tm(imax),dsqrt(sum_Q-sum_P**2.d0)
  endif
  endif
  
contains
function prob(n_s,n_l,n_interval)
implicit none
  double precision,dimension(n_interval):: p
  double precision,dimension(n_interval):: f
  double precision s 
  double precision sum_f 
  integer i
  integer n_s,n_l
  integer n_interval
  double precision,dimension(n_interval)::prob

  s = dble(1)/dble(n_interval-1)
  sum_f = 0.d0
  do i = 1,n_interval
    p(i) = s * dble(i-1)
    f(i) = p(i)**n_l * (1.d0-p(i))**n_s
    sum_f = sum_f + f(i)
  enddo
  f = f / sum_f
  prob = f

end function

end program

