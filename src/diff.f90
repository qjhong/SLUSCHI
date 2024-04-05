program diff
implicit none

  integer i,j
  double precision tmp
  integer natom
  double precision,allocatable,dimension(:,:)::pos0,pos1
  double precision,allocatable,dimension(:,:)::diff_pos
  double precision,allocatable,dimension(:)::norm_diff
  double precision avg_left,avg_right

  integer n_diffhis
  double precision,allocatable,dimension(:,:)::norm_diffhis
  double precision,allocatable,dimension(:)::max_diff
  integer n_int,max_int

  integer k1,k2,k3,k1_min,k2_min,k3_min
  double precision norm,norm_min
  double precision,dimension(3)::diff_tmp
  double precision,dimension(3,3)::latt0,latt1
  double precision diff_solid_cri,diff_liquid_cri,diff_liquid_cri_count
  double precision energy
  double precision energy_solid_cri,energy_liquid_cri

  double precision,allocatable,dimension(:,:)::diff_info,diff_solid,diff_liquid
  double precision,allocatable,dimension(:)::diff_solid_v,diff_liquid_v
  integer,allocatable,dimension(:)::count_info,count_solid,count_liquid
  double precision,allocatable,dimension(:)::energy_info,energy_solid,energy_liquid
  integer,allocatable,dimension(:)::flag_info
  double precision max_energy,min_energy
  integer nsolid,nliquid,ninfo
  double precision avg_diff_solid,avg_diff_liquid
  integer avg_count_liquid

  integer flag_select
  integer,allocatable,dimension(:)::nselect
  integer nselect0,nselect1

  double precision potim

  avg_count_liquid = 9

  open(7,file='struct0')
  open(8,file='struct1')
  open(9,file='diff_his')
  open(11,file='diff')
  open(13,file='avghalf')
  open(12,file='flag_homo')
  open(14,file='flag_solid')
  open(15,file='diff.in')
  open(16,file='stop_info.in')
  open(17,file='natoms')
  open(18,file='energy.out')
  open(19,file='potim')

  read(19,*) potim 
  close(19)

  read(15,*) diff_solid_cri
  read(15,*) diff_liquid_cri_count
  flag_select=-1
  read(15,*,err=150,end=150) flag_select
  150 continue
  diff_liquid_cri=diff_liquid_cri_count*0.75
!  diff_liquid_cri=0.d0
  close(15)
  if (flag_select>0) then
    allocate(nselect(flag_select))
    read(17,*) nselect
    nselect0=0
    do i=1,flag_select-1
      nselect0=nselect0+nselect(i)
    enddo
    nselect1=nselect0+nselect(flag_select)
  endif

  read(18,*) energy
  close(18)

  ! Read stop_info.in
  ninfo = 0
  do while (.true.)
    read(16,*,err=149,end=149)
    ninfo = ninfo + 1
  enddo
  149 continue
  allocate(diff_info(2,ninfo),count_info(ninfo),energy_info(ninfo),flag_info(ninfo))
  allocate(diff_solid(2,ninfo),count_solid(ninfo),energy_solid(ninfo))
  allocate(diff_liquid(2,ninfo),count_liquid(ninfo),energy_liquid(ninfo))
  rewind 16
  max_energy = -1.d10
  min_energy = 1.d10
  nsolid=0
  nliquid=0
  ! Seperate solid and liquid info
  do i=1,ninfo
    read(16,*) diff_info(:,i),count_info(i),energy_info(i),flag_info(i)

    if (energy_info(i)>max_energy) max_energy=energy_info(i)
    if (energy_info(i)<min_energy) min_energy=energy_info(i)
    if ( flag_info(i) == 1 ) then
      nsolid=nsolid+1
      diff_solid(:,nsolid)=diff_info(:,i)
      count_solid(nsolid)=count_info(i)
      energy_solid(nsolid)=energy_info(i)
    elseif ( flag_info(i) == 0 ) then
      nliquid=nliquid+1
      diff_liquid(:,nliquid)=diff_info(:,i)
      count_liquid(nliquid)=count_info(i)
      energy_liquid(nliquid)=energy_info(i)
    endif
  enddo
  close(16)
  ! diff_solid_cri = medium(diff_solid_v)
  ! energy_solid_cri = medium(energy_solid)
  if ( nsolid > 2 ) then
    avg_diff_solid = sum(diff_solid(1,1:nsolid))+sum(diff_solid(2,1:nsolid))
    avg_diff_solid = avg_diff_solid/dble(2*nsolid)
    allocate(diff_solid_v(1000))
    diff_solid_v(1:nsolid)=diff_solid(1,1:nsolid)
    diff_solid_v(nsolid+1:nsolid*2)=diff_solid(2,1:nsolid)
    do i=1,nsolid*2-1
    do j=2,nsolid*2-i+1
      if (diff_solid_v(j-1)<diff_solid_v(j)) then
        tmp = diff_solid_v(j)
        diff_solid_v(j) = diff_solid_v(j-1)
        diff_solid_v(j-1) = tmp
      endif
    enddo
    enddo
    diff_solid_cri = (diff_solid_v(nsolid)+diff_solid_v(nsolid+1))/2.d0
    do i=1,nsolid-1
    do j=2,nsolid-i+1
      if (energy_solid(j-1)<energy_solid(j)) then
        tmp = energy_solid(j)
        energy_solid(j) = energy_solid(j-1)
        energy_solid(j-1) = tmp
      endif
    enddo
    enddo
    if (mod(nsolid,2)==1) then
      energy_solid_cri = energy_solid(nsolid/2+1)
    else
      energy_solid_cri = energy_solid(nsolid/2)
    endif
  endif
  ! diff_liquid_cri = medium(diff_liquid_v)
  ! energy_liquid_cri = medium(energy_liquid)
  if ( nliquid > 2 ) then
    avg_diff_liquid = sum(diff_liquid(1,1:nliquid))+sum(diff_liquid(2,1:nliquid))
    avg_diff_liquid = avg_diff_liquid/dble(2*nliquid)
    allocate(diff_liquid_v(1000))
    diff_liquid_v(1:nliquid)=diff_liquid(1,1:nliquid)
    diff_liquid_v(nliquid+1:nliquid*2)=diff_liquid(2,1:nliquid)
    do i=1,nliquid*2-1
    do j=2,nliquid*2-i+1
      if (diff_liquid_v(j-1)<diff_liquid_v(j)) then
        tmp = diff_liquid_v(j)
        diff_liquid_v(j) = diff_liquid_v(j-1)
        diff_liquid_v(j-1) = tmp
      endif
    enddo
    enddo
    diff_liquid_cri = (diff_liquid_v(nliquid)+diff_liquid_v(nliquid+1))/2.d0
    avg_count_liquid = sum(count_liquid(1:nliquid))
    avg_count_liquid = ceiling(dble(avg_count_liquid)/dble(nliquid))
    do i=1,nliquid-1
    do j=2,nliquid-i+1
      if (energy_liquid(j-1)<energy_liquid(j)) then
        tmp = energy_liquid(j)
        energy_liquid(j) = energy_liquid(j-1)
        energy_liquid(j-1) = tmp
      endif
    enddo
    enddo
    if (mod(nliquid,2)==1) then
      energy_liquid_cri = energy_liquid(nliquid/2+1)
    else
      energy_liquid_cri = energy_liquid(nliquid/2)
    endif
  endif

  natom = 0
  do while (.true.)
    read(7,"(3F12.5)",err=147,end=147) 
    natom = natom + 1
  enddo
  147 continue
  natom = natom - 3

  allocate(pos0(3,natom),pos1(3,natom))
  allocate(diff_pos(3,natom))
  allocate(norm_diff(natom))
  rewind 7
  read(7,"(3F12.5)") latt0(:,1)
  read(7,"(3F12.5)") latt0(:,2)
  read(7,"(3F12.5)") latt0(:,3)
  do i=1,natom
    read(7,"(3F12.5)") pos0(:,i)
  enddo
  close(7)
  read(8,"(3F12.5)") latt1(:,1)
  read(8,"(3F12.5)") latt1(:,2)
  read(8,"(3F12.5)") latt1(:,3)
  do i=1,natom
    read(8,"(3F12.5)") pos1(:,i)
  enddo
  close(8)

  do i=1,natom
    norm_min = 100.d0
    do k1=-2,2
    do k2=-2,2
    do k3=-2,2
      diff_tmp = pos0(:,i)-pos1(:,i) + latt1(:,1)*dble(k1) + latt1(:,2)*dble(k2) + latt1(:,3)*dble(k3)
      norm = dsqrt(diff_tmp(1)**2.d0 + diff_tmp(2)**2.d0 + diff_tmp(3)**2.d0)
      norm = norm*dsqrt(1.5d0/potim)
      if (norm < norm_min) then
        norm_min = norm
        k1_min = k1
        k2_min = k2
        k3_min = k3
      endif
    enddo
    enddo
    enddo
    !if (i==natom) then
    !    write(*,*) pos0(:,i),pos1(:,i)
    !endif
    diff_pos(:,i) = pos0(:,i)-pos1(:,i) + latt1(:,1)*dble(k1_min) + latt1(:,2)*dble(k2_min) + latt1(:,3)*dble(k3_min)
    norm_diff(i) = norm_min
  enddo

  do i=1,natom
    write(11,99) norm_diff(i)
  enddo
  99 format(F10.5)
  close(11)

  avg_left = 0.d0
  if ( flag_select > 0 ) then
    do i=nselect0+1,nselect1
      avg_left = avg_left + norm_diff(i)
    enddo
    avg_right = 0.d0
    do i=natom/2+nselect0+1,natom/2+nselect1
      avg_right = avg_right + norm_diff(i)
    enddo
    avg_left = avg_left/(nselect1-nselect0)
    avg_right = avg_right/(nselect1-nselect0)
  else
    do i=1,natom/2
      avg_left = avg_left + norm_diff(i)
    enddo
    avg_right = 0.d0
    do i=natom/2+1,natom
      avg_right = avg_right + norm_diff(i)
    enddo
    avg_left = avg_left/natom*2
    avg_right = avg_right/natom*2
  endif
  write(13,*) avg_left,avg_right

  ! for solid: left < 0.001 .and. right < 0.001
  ! for liquid: max diff history; count > 0.01; max interval < 5
  n_diffhis = 0
  do while (.true.)
    read(9,*,end=148,err=148)
    n_diffhis = n_diffhis + 1
  enddo
  148 continue
  n_diffhis = n_diffhis/natom + 1
  allocate(norm_diffhis(natom,n_diffhis))
  allocate(max_diff(natom))
  rewind 9
  do i=1,n_diffhis-1
  do j=1,natom
    read(9,*) norm_diffhis(j,i)
  enddo
  enddo
  norm_diffhis(1:natom,n_diffhis) = norm_diff
  max_diff = 0.d0
  do j=1,natom
  do i=1,n_diffhis
    if (norm_diffhis(j,i)>max_diff(j)) max_diff(j) = norm_diffhis(j,i)
  enddo
  enddo
  !write(*,*) max_diff
  n_int = 0
  max_int = 0
  do i=1,natom
    if ( max_diff(i) > diff_liquid_cri_count ) then
      if ( n_int > max_int ) max_int = n_int
      n_int = 0
    else
      n_int = n_int + 1
    endif
  enddo
  if ( n_int > max_int ) max_int = n_int
  !write(*,*) diff_solid_cri,diff_liquid_cri,diff_solid_cri*0.9d0+diff_liquid_cri*0.1d0
  !write(*,*) energy_solid_cri,energy_liquid_cri,energy_solid_cri*0.9d0+energy_liquid_cri*0.1d0
  if ( ( (avg_left < diff_solid_cri*1.1d0 .and. avg_right < diff_solid_cri*1.1d0) &
  & .and. avg_left/avg_right > 0.5d0 .and. avg_left/avg_right < 2.d0 ) &
  & .or. &
  & ( (nsolid>=3 .and. nliquid>=3 &
  & .and. avg_left < (diff_solid_cri*0.9d0+diff_liquid_cri*0.1d0) &
  & .and. avg_right < (diff_solid_cri*0.9d0+diff_liquid_cri*0.1d0) & 
  & .and. energy < (energy_solid_cri*0.9d0+energy_liquid_cri*0.1d0)) &
  & .and. avg_left/avg_right > 0.33d0 .and. avg_left/avg_right < 3.d0 ) ) then ! solid
    write(12,*) '1'
    write(14,*) '1'
  elseif ( ( (avg_left>diff_liquid_cri*0.9d0 .and. avg_right>diff_liquid_cri*0.9d0) &
  & .or. (nsolid>=3 .and. nliquid>=3 &
  & .and. avg_left > (diff_solid_cri*0.1d0+diff_liquid_cri*0.9d0) &
  & .and. avg_right > (diff_solid_cri*0.1d0+diff_liquid_cri*0.9d0) &
  & .and. energy > (energy_solid_cri*0.1d0+energy_liquid_cri*0.9d0)) ) &
  & .and. max_int < floor(dble(avg_count_liquid)*1.1d0) .and. &
  & avg_left/avg_right > 0.8d0 .and. avg_left/avg_right < 1.25d0 ) then ! liquid
    write(12,*) '1'
    write(14,*) '0'
  elseif (avg_left < (diff_solid_cri*0.8d0+diff_liquid_cri*0.2d0) .and. &
  & avg_right < (diff_solid_cri*0.8d0+diff_liquid_cri*0.2d0) ) then !maybe solid
    write(12,*) '0'
    write(14,*) '11'
  elseif (avg_left > (diff_solid_cri*0.2d0+diff_liquid_cri*0.8d0) .and. &
  & avg_right > (diff_solid_cri*0.2d0+diff_liquid_cri*0.8d0) ) then !maybe liquid
    write(12,*) '0'
    write(14,*) '10'
  else
    write(12,*) '0'
    write(14,*) '-1'
  endif
  write(13,*) max_int

end
