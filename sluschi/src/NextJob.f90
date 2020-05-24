program NextJob
implicit none
! This code reads
! 1. jobs_results.out: temp, nsolid, nliquid.
! 2. confident.in: 0 or 1. Whether you are confident (1) or not (0) about MP you guess.
! 3. MPFit.out (optional): if MPFit was done previously, read melting temperature (MP) and error (STD_MP).
! 4. ErrorTarget.in: the fitting error specified in job.in.
! This code outputs
! 1. jobs_to_run.out: temp, nsum
! 2. FlagStop.out: 1 STD_MP < ErrorTarget, STOP!
  integer nwork
  double precision,allocatable,dimension(:)::temp,temp2
  integer,allocatable,dimension(:)::nsolid,nliquid,nsum,nsolid2,nliquid2,nsum2
  integer max_nsum
  integer sum_solid,sum_liquid

  integer i,j
  integer flag_confident
  double precision delt,delt0
  integer ncount,i_transit
  double precision MP,STD_MP,STD_TARGET
  double precision dT
  integer inext
  double precision maxgap

  integer return_fit
  double precision Tm_return,Tm_err_return
  integer imax,jmax
  double precision efficiency

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
    nsum(i) = nsolid(i) + nliquid(i)
  enddo
  close(7)
  max_nsum = 0
  do i=1,nwork
    if (nsum(i)>max_nsum) max_nsum = nsum(i)
  enddo

  open(8,file='jobs_to_run.out')
  delt0 = dble(ceiling(temp(1)/2000.d0))*200.d0
  do i=1,nwork
    if (nliquid(i)>0) then
      delt0 = dble(ceiling(temp(i)/2000.d0))*200.d0
      exit
    endif
  enddo
!  if (max_nsum < 4) then
    ! Search Mode
    ! determine step size: if not confident, dT = 800 K; else, dT = 1/2/3/400K
    open(9,file='confident.in')
    read(9,*) flag_confident
    close(9)
    if (flag_confident==1) then 
      delt = delt0
    else
      delt = 800.d0
    endif
    ! extend temperature region if
    ! (0) all solid / liquid
    ! (1) there is liquid in lowest temperature  or
    ! (2) there is solid in highest temperature 
    sum_solid = 0
    sum_liquid = 0
    do i=1,nwork
      sum_solid = sum_solid + nsolid(i)
      sum_liquid = sum_liquid + nliquid(i)
    enddo
    if ( nliquid(1) > 0 .and. sum_solid == 0 ) then
      write(8,'(F12.1,A)') temp(1)-delt,' 1'
    endif
    if ( nsolid(nwork) > 0 .and. sum_liquid == 0) then
      write(8,'(F12.1,A)') temp(nwork)+delt,' 1'
    endif
    ! T(i-1) is solid, T(i) is liquid, check T(i)-T(i-1). 
    ! If DelT = 100 K, search is done. Else, add point in between
    do i=2,nwork
      if (nsolid(i-1)>0 .and. (nliquid(i)>0 .or. nliquid(i-1)>0) ) then
        if (temp(i)-temp(i-1)<delt0+1.d0) then
          write(8,*) temp(i-1),' 4'
          write(8,*) temp(i),' 4'
        elseif ( mod(nint(temp(i)-temp(i-1)),nint(delt0)*4)==0 ) then
          write(8,*) (temp(i-1)+temp(i))/2.d0,' 1'
        else
          !write(8,*) temp(i-1)+delt0,' 1'
          write(8,*) (temp(i-1)+temp(i))/2.d0,' 1'
        endif
      endif
    enddo
!  else
  if (max_nsum >= 4) then
    delt=delt0
    ! extend temperature region if
    ! (1) there is liquid in lowest temperature  or
    ! (2) there is solid in highest temperature 
    if ( nliquid(1) > 0 ) then
      write(8,*) temp(1)-delt,' 1'
      goto 999
    endif
    if ( nsolid(nwork) > 0 ) then
      write(8,*) temp(nwork)+delt,' 1'
      goto 999
    endif
    ! if there is a solid-liquid transition
    do i=2,nwork
      if (nsolid(i-1)>0 .and. nliquid(i)>0 .and. (nsum(i-1)<4 .or. nsum(i)<4) ) then
        if (temp(i)-temp(i-1)<delt0+1.d0) then
          write(8,*) temp(i-1),' 4'
          write(8,*) temp(i),' 4'
        else
          write(8,*) (temp(i-1)+temp(i))/2.d0,' 1'
          !if (nsum(i-1)<4) then
          !  write(8,*) temp(i)-delt,' 4'
          !else
          !  write(8,*) temp(i-1)+delt,' 4'
          !endif
        endif
        goto 999
      endif
    enddo
    ! Strategy: T': T with fractional distribution (not 0 or 1)
    ! 1. for all T', run at least 4 duplicates
    ! 2. if number of T' <= 1, run one more T in between
    ! 2.5 if number of T' >= 2, fit
    ! 3. T' closest to MP, run more
    ! Case 1. for all T', run at least 4 duplicates
    do i=2,nwork-1
      if (nsolid(i)>0 .and. nliquid(i)>0 .and. (nsum(i-1)<4 .or. nsum(i)<4 .or. nsum(i+1)<4) ) then
        write(8,*) temp(i-1),' 4'
        write(8,*) temp(i),' 4'
        write(8,*) temp(i+1),' 4'
        goto 999
      endif
    enddo
    ! Case 2. if number of T' <= 1, run one more T in between
    ncount = 0
    do i=2,nwork
      if (nsolid(i)>0 .and. nliquid(i)>0) then
        ncount = ncount + 1
      endif
      if (nliquid(i-1)==0 .and. nliquid(i)>0) then
        i_transit = i
      endif
    enddo
    if (ncount == 0) then
      write(8,*) (temp(i_transit-1)+temp(i_transit))/2.d0,' 4'
      goto 999
    elseif (ncount == 1) then
      if (nliquid(i_transit)<=2) then
        write(8,*) (temp(i_transit)+temp(i_transit+1))/2.d0,' 4'
      elseif (nliquid(i_transit)==3) then
        write(8,*) (temp(i_transit-1)+temp(i_transit))/2.d0,' 4'
      endif
      goto 999
    endif
    
    ! Now Run a MP Fit
    call MPFit_NextJob(temp,nsolid,nliquid,nsum,nwork,return_fit,Tm_return,Tm_err_return)
    write(*,*) "current fit"
    write(*,*) return_fit,Tm_return,Tm_err_return
    write(*,*) "possibilities:"
    STD_MP = 10000.d0
    if ( return_fit == 1 ) STD_MP = Tm_err_return
    open(11,file='ErrorTarget.in')
    read(11,*) STD_TARGET
    close(11)
    ! 3. T' closest to MP, run more
    if (STD_MP > STD_TARGET) then
      ! If MPFit works, rely on it to predict the next job to run
      if (return_fit == 1) then 
        imax = 0
        jmax = 0
        ! 1. more samples at existing temperatures
        allocate(temp2(nwork),nsolid2(nwork),nliquid2(nwork),nsum2(nwork))
        efficiency = 0.d0
        do i=1,nwork
          if (nsum(i) < 16) then
            temp2 = temp
            nsolid2 = nsolid
            nliquid2 = nliquid
            nsum2 = nsum
            nsolid2(i) = nsolid2(i)*2
            nliquid2(i) = nliquid2(i)*2
            nsum2(i) = nsum2(i)*2
            write(*,*) temp2(i),nsolid2(i),nliquid2(i),nsum2(i)
            call MPFit_NextJob(temp2,nsolid2,nliquid2,nsum2,nwork,return_fit,Tm_return,Tm_err_return)
            write(*,*) return_fit,Tm_return,Tm_err_return
            if ( return_fit == 1 .and. (STD_MP-Tm_err_return)/dble(nsum2(i)/2) > efficiency) then
                efficiency = (STD_MP-Tm_err_return)/dble(nsum2(i)/2)
                imax = i
              endif
              if (i>1 .and. i<nwork) then
              nsolid2(i) = floor(nsum2(i) * ( dble(nsolid2(i-1))/dble(nsum2(i-1)) + dble(nsolid2(i+1))/dble(nsum2(i+1)) ) / 2.d0)
              nliquid2(i) = nsum2(i)-nsolid2(i)
              if (nsolid2(i) < nsolid(i)) then
                nsolid2(i) = nsolid(i)
                nliquid2(i) = nsum2(i)-nsolid2(i)
              endif
              if (nliquid2(i) < nliquid(i)) then
                nliquid2(i) = nliquid(i)
                nsolid2(i) = nsum2(i)-nliquid2(i)
              endif
              write(*,*) temp2(i),nsolid2(i),nliquid2(i),nsum2(i)
              call MPFit_NextJob(temp2,nsolid2,nliquid2,nsum2,nwork,return_fit,Tm_return,Tm_err_return)
              write(*,*) return_fit,Tm_return,Tm_err_return
              if ( return_fit == 1 .and. (STD_MP-Tm_err_return)/dble(nsum2(i)/2) > efficiency) then
                efficiency = (STD_MP-Tm_err_return)/dble(nsum2(i)/2)
                imax = i
              endif
            endif
          endif
        enddo
        deallocate(temp2,nsolid2,nliquid2,nsum2)
        ! 2. new temperature
        allocate(temp2(nwork+1),nsolid2(nwork+1),nliquid2(nwork+1),nsum2(nwork+1))
        do i=1,nwork-1
          temp2(1:i) = temp(1:i)
          temp2(i+2:nwork+1) = temp(i+1:nwork)
          nsolid2(1:i) = nsolid(1:i)
          nsolid2(i+2:nwork+1) = nsolid(i+1:nwork)
          nliquid2(1:i) = nliquid(1:i)
          nliquid2(i+2:nwork+1) = nliquid(i+1:nwork)
          nsum2(1:i) = nsum(1:i)
          nsum2(i+2:nwork+1) = nsum(i+1:nwork)
          temp2(i+1) = (temp2(i)+temp2(i+2))/2.d0
          nsum2(i+1) = 4
          nsolid2(i+1) = floor(4.d0 * ( dble(nsolid2(i))/dble(nsum2(i)) + dble(nsolid2(i+2))/dble(nsum2(i+2)) ) / 2.d0)
          if ( nsolid2(i)>0 .and. nsolid2(i+2)>0) then
            nsolid2(i+1) = max(nsolid2(i+1),1)
          endif
          nliquid2(i+1) = 4 - nsolid2(i+1)
          if ( nliquid2(i)>0 .and. nliquid2(i+2)>0) then
            nliquid2(i+1) = max(nliquid2(i+1),1)
            nsolid2(i+1) = 4 - nliquid2(i+1)
          endif
          write(*,*) temp2(i+1),nsolid2(i+1),nliquid2(i+1),nsum2(i+1)
          call MPFit_NextJob(temp2,nsolid2,nliquid2,nsum2,nwork+1,return_fit,Tm_return,Tm_err_return)
          write(*,*) return_fit,Tm_return,Tm_err_return
          if ( return_fit == 1 .and. (STD_MP-Tm_err_return)/4.d0 > efficiency) then
            efficiency = (STD_MP-Tm_err_return)/4.d0
            jmax = i
          endif
        enddo
        if (jmax == 0) then
          if (nsum(imax)==1 .or. nsum(imax)==2 .or. nsum(imax)==4 .or. nsum(imax)==8) then
            write(8,*) temp(imax),nsum(imax)*2
          endif
        else
          write(8,*) (temp(jmax)+temp(jmax+1))/2.d0,4
        endif
      else
        ! 3.1 If dabs(T1-MP)<dabs(T2-MP), nsum(at T1)>=nsum(at T2). This is reasonable.
        do i=1,nwork
          do j=1,nwork
            if (dabs(temp(j)-MP)>dabs(temp(i)-MP) .and. nsum(j) > nsum(i)) then
              write(8,*) temp(i),nsum(j)
              goto 999
            endif
          enddo
        enddo
        ! 3.2 run max_nsum duplicates for closest T'
        dT = 10000.d0
        inext = 100
        do i=1,nwork
          if (nsolid(i)>0 .and. nliquid(i)>0 .and. nsum(i)<max_nsum .and. dabs(temp(i)-MP) < dT) then
            dT = dabs(temp(i)-MP)
            inext = i
          endif
        enddo 
        if (inext < 100) then ! there exists T', such that nsum(inext)<max_nsum
          write(8,*) temp(inext),max_nsum
          goto 999
        endif
        ! 3.3 increase max_nsum
        if (max_nsum < 16) then ! these doesn't exist T', that nsum(i) < max_nsum
          dT = 10000.d0
          inext = 100
          do i=1,nwork
            if (dabs(temp(i)-MP) < dT) then
              dT = dabs(temp(i)-MP)
              inext = i
            endif
          enddo 
          if (max_nsum==1 .or. max_nsum==2 .or. max_nsum==4 .or. max_nsum==8) then
            write(8,*) temp(inext),max_nsum*2 ! number of duplicates: 4->8->16
            goto 999
          endif
        endif
        ! 3.4 add temp in between
        maxgap = 0.d0
        do i=2,nwork ! For all T', # of duplicates == 16
          if ( (nsolid(i-1)>0 .and. nliquid(i-1)>0 .or. nsolid(i)>0 .and. nliquid(i)>0) .and. temp(i)-temp(i-1)>maxgap) then
            maxgap = temp(i)-temp(i-1)
            inext = i
          endif
        enddo 
        write(8,*) (temp(inext-1)+temp(inext))/2.d0,' 4' ! Add more T in between
      endif
    else
      open(12,file='FlagStop.out')
      write(12,*) '1'
      close(12)
    endif
    
  endif 

  999 continue
  close(8)
    
contains

subroutine MPFit_NextJob(temp,nsolid,nliquid,nsum,nwork,return_fit,Tm_return,Tm_err_return)

implicit none

  integer,parameter::n_interval=1001,n=n_interval-1
  integer,parameter::n_randomnumber=100000
  double precision,parameter:: eV=1.602d-19,kB=1.3806d-23,F_converge=160.d0

  integer,dimension(nwork)::nsolid,nliquid,nsum
  double precision,dimension(nwork)::temp
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
  integer return_fit
  double precision Tm_return,Tm_err_return

  ! To perform a fit, jobs_results.out has to satisfy
  ! 1. there are at least two T' with >4 duplictates. T': with fractional distribution (not 0 or 1)
  ! 2. each neighbor of T' with >4 duplicates
  ncount = 0
  if (nliquid(1)==0 .and. nsolid(nwork)==0) then
    do i=2,nwork-1
      if (nsolid(i)>0 .and. nliquid(i)>0 .and. nsum(i)>=4 .and. nsum(i-1)>=4 .and. nsum(i+1)>=4 ) ncount = ncount+1
    enddo
  endif
  if (ncount < 2) then
    return_fit=0
  else
    return_fit=1
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
!    write(9,'(4F20.5)') Tm_min,Tm_max,gam_min,gam_max

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
    Tm_return = sum_P
    Tm_err_return = dsqrt(sum_Q-sum_P**2.d0)
  else
    Tm_return = Tm(imax)
    Tm_err_return = dsqrt(sum_Q-sum_P**2.d0)
  endif
  endif
  
end subroutine

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

