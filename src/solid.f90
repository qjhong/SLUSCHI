program solid
implicit none

  integer,parameter:: n_rot = 32, n_spcl = 10
  double precision,parameter:: PI=3.14159d0
  double precision factor
  double precision,dimension(3,3)::latt,latt0,latt_tmp,latt_new,latt_new_inv
  double precision,dimension(3,3)::Rx,Ry,Rz
  integer,dimension(3,3)::index_min,index_min_all

  double precision,dimension(3)::pos_tmp,rad_min,norm,angle,pos_cart,pos_frac
  double precision rad2,rad_min_all,norm_10

  character str*40
  integer ix,iy,iz,jx,jy,jz
  double precision thetax,thetay,thetaz
  integer i,j

  integer ntype
  integer,dimension(10)::natom,natom_type
  logical flag_vasp5
  double precision,allocatable,dimension(:,:,:)::pos
  double precision det
  integer max_natom
  double precision a

  open(8,file='POSCAR') 
  open(10,file='size.in')
  
  read(8,*) str
  read(8,*) factor
  read(8,*) latt(:,1)
  latt(:,1) = latt(:,1)*factor
  read(8,*) latt(:,2)
  latt(:,2) = latt(:,2)*factor
  read(8,*) latt(:,3)
  latt(:,3) = latt(:,3)*factor

  latt0(:,1) = (/ 1.0,0.0,0.0 /)
  latt0(:,2) = (/ 0.0,1.0,0.0 /)
  latt0(:,3) = (/ 0.0,0.5,0.866 /)
  !latt0(:,3) = (/ 0.0,0.0,1.0 /)
  read(10,*) a
  close(10)
  latt0 = latt0 * a

  rad_min_all = 1.d10
  do ix = 1,n_rot
    thetax = 2.d0*PI*dble(ix)/dble(n_rot)
    Rx(:,1) = (/ 1.d0,0.d0,0.d0/)
    Rx(:,2) = (/ 0.d0,dcos(thetax),dsin(thetax)/)
    Rx(:,3) = (/ 0.d0,-dsin(thetax),dcos(thetax)/)
  do iy = 1,n_rot
    thetay = 2.d0*PI*dble(iy)/dble(n_rot)
    Ry(:,1) = (/ dcos(thetay),0.d0,-dsin(thetay)/)
    Ry(:,2) = (/ 0.d0,1.d0,0.d0/)
    Ry(:,3) = (/ dsin(thetay),0.d0,dcos(thetay)/)
  do iz = 1,n_rot
    thetaz = 2.d0*PI*dble(iz)/dble(n_rot)
    Rz(:,1) = (/ dcos(thetaz),dsin(thetaz),0.d0/)
    Rz(:,2) = (/ -dsin(thetaz),dcos(thetaz),0.d0/)
    Rz(:,3) = (/ 0.d0,0.d0,1.d0/)

    latt_tmp = latt0
    ! generate ideal cell, after rotation
    do j=1,3
      rad_min(j) = 1.d10
      latt_tmp(:,j) = latt_tmp(1,j)*Rx(:,1)+latt_tmp(2,j)*Rx(:,2)+latt_tmp(3,j)*Rx(:,3)
      latt_tmp(:,j) = latt_tmp(1,j)*Ry(:,1)+latt_tmp(2,j)*Ry(:,2)+latt_tmp(3,j)*Ry(:,3)
      latt_tmp(:,j) = latt_tmp(1,j)*Rz(:,1)+latt_tmp(2,j)*Rz(:,2)+latt_tmp(3,j)*Rz(:,3)
      do jx=-n_spcl,n_spcl
      do jy=-n_spcl,n_spcl
      do jz=-n_spcl,n_spcl
        pos_tmp = latt(:,1)*dble(jx) + latt(:,2)*dble(jy) + latt(:,3)*dble(jz)
        norm_10 = dsqrt( pos_tmp(1)**2.d0 + pos_tmp(2)**2.d0 + pos_tmp(3)**2.d0 )
        pos_tmp = pos_tmp - latt_tmp(:,j)
        rad2 = pos_tmp(1)**2.d0 + pos_tmp(2)**2.d0 + pos_tmp(3)**2.d0
        if ( norm_10 < a ) rad2 = rad2*2.0
        if ( rad2 < rad_min(j) ) then
          rad_min(j) = rad2
          index_min(:,j) = (/jx,jy,jz/)
        endif
      enddo
      enddo
      enddo
    enddo

    if ( rad_min(1)+rad_min(2)+rad_min(3) < rad_min_all ) then
      rad_min_all = rad_min(1)+rad_min(2)+rad_min(3)
      index_min_all = index_min
    endif
    
  enddo
  enddo
  enddo

  do j=1,3
    latt_new(:,j) = latt(:,1)*dble(index_min_all(1,j)) + latt(:,2)*dble(index_min_all(2,j)) + latt(:,3)*dble(index_min_all(3,j))
  enddo
  do j=1,3
    norm(j) = dsqrt(latt_new(1,j)**2.d0+latt_new(2,j)**2.d0+latt_new(3,j)**2.d0)
  enddo
  pos_tmp = latt_new(:,2) - latt_new(:,3)
  angle(1) = dacos( ( norm(2)**2.d0 + norm(3)**2.d0 - pos_tmp(1)**2.d0 - &
  & pos_tmp(2)**2.d0 - pos_tmp(3)**2.d0 ) / (norm(2)*norm(3)*2.d0) ) / PI * 180.d0
  pos_tmp = latt_new(:,1) - latt_new(:,3)
  angle(2) = dacos( ( norm(1)**2.d0 + norm(3)**2.d0 - pos_tmp(1)**2.d0 - &
  & pos_tmp(2)**2.d0 - pos_tmp(3)**2.d0 ) / (norm(1)*norm(3)*2.d0) ) / PI * 180.d0
  pos_tmp = latt_new(:,1) - latt_new(:,2)
  angle(3) = dacos( ( norm(1)**2.d0 + norm(2)**2.d0 - pos_tmp(1)**2.d0 - &
  & pos_tmp(2)**2.d0 - pos_tmp(3)**2.d0 ) / (norm(1)*norm(2)*2.d0) ) / PI * 180.d0

  write(*,'(A)') '*** Generate a supercell from the current unitcell ***'
  write(*,'(A)') 'The supercell is:'
  write(*,*) latt_new(:,1)
  write(*,*) latt_new(:,2)
  write(*,*) latt_new(:,3)
  write(*,'(A)') '|a|,|b|,|c|,theta(bc),theta(ac),theta(ab):'
  write(*,'(6F10.3)') norm,angle

  ! one more line in vasp5
  flag_vasp5 = .true.
  rewind 8
  do i=1,5
    read(8,*)
  enddo
  read(8,*,err=148,end=148) natom(1)
  flag_vasp5 = .false.
  148 continue
  
  ntype = 0
  natom = 0
  do while (.true.)
    rewind 8
    do i=1,5
      read(8,*)
    enddo
    if (flag_vasp5) read(8,*)
    read(8,*,err=147,end=147) natom(1:ntype)
    ntype=ntype+1
  enddo
  147 continue
  ntype = ntype-1
  rewind 8
  do i=1,5
    read(8,*)
  enddo
  if (flag_vasp5) read(8,*)
  read(8,*) natom(1:ntype)
  max_natom = 0
  write(*,'(A)',advance='no') "In UNIT-cell, number of atoms:"
  do i=1,ntype
    if (natom(i)>max_natom) max_natom = natom(i)
    write(*,'(I5)',advance='no') natom(i)
  enddo
  write(*,'(A)',advance='no') " total: "
  write(*,'(I5)') sum(natom)
  read(8,*) str
  allocate(pos(3,max_natom,ntype))
  do i=1,ntype
  do j=1,natom(i)
    read(8,*) pos(:,j,i)
  enddo
  enddo
  if ( str == 'D' .or. str == ' D' .or. str == 'Direct' ) then
    do i=1,ntype
    do j=1,natom(i)
      pos(:,j,i) = pos(1,j,i)*latt(:,1) + pos(2,j,i)*latt(:,2) + pos(3,j,i)*latt(:,3)
!      pos(2,j,i) = pos(1,j,i)*latt(2,1) + pos(2,j,i)*latt(2,2) + pos(3,j,i)*latt(2,3)
!      pos(3,j,i) = pos(1,j,i)*latt(3,1) + pos(2,j,i)*latt(3,2) + pos(3,j,i)*latt(3,3)
    enddo
    enddo
  elseif ( str == 'C' .or. str == ' C' ) then
    pos = pos*factor
  endif

  ! inverse matrix of latt_new
  det = latt_new(1,1)*latt_new(2,2)*latt_new(3,3) + latt_new(1,2)*latt_new(2,3)*latt_new(3,1) &
  & + latt_new(1,3)*latt_new(2,1)*latt_new(3,2) - latt_new(1,3)*latt_new(2,2)*latt_new(3,1) &
  & - latt_new(1,2)*latt_new(2,1)*latt_new(3,3) - latt_new(1,1)*latt_new(2,3)*latt_new(3,2)
  latt_new_inv(1,1) = latt_new(2,2)*latt_new(3,3) - latt_new(2,3)*latt_new(3,2)
  latt_new_inv(1,2) = latt_new(1,3)*latt_new(3,2) - latt_new(3,3)*latt_new(1,2)
  latt_new_inv(1,3) = latt_new(1,2)*latt_new(2,3) - latt_new(2,2)*latt_new(1,3)
  latt_new_inv(2,1) = latt_new(2,3)*latt_new(3,1) - latt_new(3,3)*latt_new(2,1)
  latt_new_inv(2,2) = latt_new(1,1)*latt_new(3,3) - latt_new(3,1)*latt_new(1,3)
  latt_new_inv(2,3) = latt_new(1,3)*latt_new(2,1) - latt_new(2,3)*latt_new(1,1)
  latt_new_inv(3,1) = latt_new(2,1)*latt_new(3,2) - latt_new(2,2)*latt_new(3,1)
  latt_new_inv(3,2) = latt_new(1,2)*latt_new(3,1) - latt_new(3,2)*latt_new(1,1)
  latt_new_inv(3,3) = latt_new(1,1)*latt_new(2,2) - latt_new(2,1)*latt_new(1,2)
  latt_new_inv = latt_new_inv/det
  write(*,'(A)') 'Inverse Matrix is:'
  write(*,*) latt_new_inv(:,1)
  write(*,*) latt_new_inv(:,2)
  write(*,*) latt_new_inv(:,3)

  !do i=1,3
  !  write(*,*) latt_new_inv(:,1)*latt_new(1,i) + latt_new_inv(:,2)*latt_new(2,i) + latt_new_inv(:,3)*latt_new(3,i)
  !enddo

  open(9,file='POSCAR_STRCT')
  open(10,file='POSCAR_HEADER')
  write(10,'(A)') 'Structure by SLUSCHI'
  write(10,'(A)') '1.0'
  do i=1,3
  write(10,*) latt_new(:,i)
  enddo
  natom_type = 0
  do i = 1,ntype
  do ix=-n_spcl,n_spcl
  do iy=-n_spcl,n_spcl
  do iz=-n_spcl,n_spcl
  do j = 1,natom(i)
    pos_cart(1) = latt(1,1)*dble(ix) + latt(1,2)*dble(iy) + latt(1,3)*dble(iz) + pos(1,j,i)! Cartisian
    pos_cart(2) = latt(2,1)*dble(ix) + latt(2,2)*dble(iy) + latt(2,3)*dble(iz) + pos(2,j,i)! Cartisian
    pos_cart(3) = latt(3,1)*dble(ix) + latt(3,2)*dble(iy) + latt(3,3)*dble(iz) + pos(3,j,i)! Cartisian
    pos_frac = latt_new_inv(:,1)*pos_cart(1) + latt_new_inv(:,2)*pos_cart(2) + latt_new_inv(:,3)*pos_cart(3) ! fractional
    if ( pos_frac(1)>0.d0-1.2345d-6 .and. pos_frac(1)<1.d0-1.2345d-6 & 
    & .and.  pos_frac(2)>0.d0-1.2345d-6 .and. pos_frac(2)<1.d0-1.2345d-6 &
    & .and.  pos_frac(3)>0.d0-1.2345d-6 .and. pos_frac(3)<1.d0-1.2345d-6) then
      write(9,*) pos_cart
      write(180,*) pos_frac
      natom_type(i) = natom_type(i) + 1
    endif
  enddo
  enddo
  enddo
  enddo !j
  enddo !i
  close(9)
  write(*,'(A)',advance='no') "In SUPER-cell, number of atoms:"
  do i=1,ntype
    write(*,'(I5)',advance='no') natom_type(i)
  enddo
  write(*,'(A)',advance='no') " total:"
  write(*,'(I5)',advance='no') sum(natom_type)
  write(*,*)
  do i=1,ntype
    write(10,'(I5)',advance='no') natom_type(i)
  enddo
  write(10,*) 
  write(10,'(A)') 'C'
  close(10)
  close(8)

end program solid
