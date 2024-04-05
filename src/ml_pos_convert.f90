program solid
implicit none

  integer,parameter:: n_rot = 32, n_spcl = 10
  double precision,parameter:: PI=3.14159d0
  double precision factor
  double precision,dimension(3,3)::latt,latt0,latt_tmp,latt_new,latt_inv
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
  character,allocatable,dimension(:,:,:)::flag
  double precision det
  integer max_natom
  double precision a

  open(8,file='POSCAR') 
  
  read(8,*) str
  read(8,*) factor
  read(8,*) latt(:,1)
  latt(:,1) = latt(:,1)*factor
  read(8,*) latt(:,2)
  latt(:,2) = latt(:,2)*factor
  read(8,*) latt(:,3)
  latt(:,3) = latt(:,3)*factor

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
  do i=1,ntype                                                                  
    if (natom(i)>max_natom) max_natom = natom(i)                                
  enddo                                                                         
  max_natom = max_natom * 2
  read(8,*) 
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

  ! inverse matrix of latt                                                  
  det = latt(1,1)*latt(2,2)*latt(3,3) + latt(1,2)*latt(2,3)*latt(3,1) &
  & + latt(1,3)*latt(2,1)*latt(3,2) - latt(1,3)*latt(2,2)*latt(3,1) &
  & - latt(1,2)*latt(2,1)*latt(3,3) - latt(1,1)*latt(2,3)*latt(3,2)
  latt_inv(1,1) = latt(2,2)*latt(3,3) - latt(2,3)*latt(3,2) 
  latt_inv(1,2) = latt(1,3)*latt(3,2) - latt(3,3)*latt(1,2) 
  latt_inv(1,3) = latt(1,2)*latt(2,3) - latt(2,2)*latt(1,3) 
  latt_inv(2,1) = latt(2,3)*latt(3,1) - latt(3,3)*latt(2,1) 
  latt_inv(2,2) = latt(1,1)*latt(3,3) - latt(3,1)*latt(1,3) 
  latt_inv(2,3) = latt(1,3)*latt(2,1) - latt(2,3)*latt(1,1) 
  latt_inv(3,1) = latt(2,1)*latt(3,2) - latt(2,2)*latt(3,1) 
  latt_inv(3,2) = latt(1,2)*latt(3,1) - latt(3,2)*latt(1,1) 
  latt_inv(3,3) = latt(1,1)*latt(2,2) - latt(2,1)*latt(1,2) 
  latt_inv = latt_inv/det                                               
  write(*,'(A)') 'Inverse Matrix is:'                                           
  write(*,*) latt_inv(:,1)                                                  
  write(*,*) latt_inv(:,2)                                                  
  write(*,*) latt_inv(:,3)                                                  


  open(9,file='POSCAR_ml')
  write(9,'(A)') 'Structure Converted by SLUSCHI'
  write(9,'(A)') '1.0'
  do i=1,3
  write(9,*) latt(:,i)
  enddo
  do i=1,ntype/2
    write(9,'(I5)',advance='no') natom(i)*2
  enddo
  write(9,*) 
  write(9,*) 'Selective dynamics'
  write(9,'(A)') 'C'
  natom_type = 0
  do i = 1,ntype/2
  do ix=-n_spcl,n_spcl
  do iy=-n_spcl,n_spcl
  do iz=-n_spcl,n_spcl
  do j = 1,natom(i)
    pos_cart(1) = latt(1,1)*dble(ix) + latt(1,2)*dble(iy) + latt(1,3)*dble(iz) + pos(1,j,i)! Cartisian
    pos_cart(2) = latt(2,1)*dble(ix) + latt(2,2)*dble(iy) + latt(2,3)*dble(iz) + pos(2,j,i)! Cartisian
    pos_cart(3) = latt(3,1)*dble(ix) + latt(3,2)*dble(iy) + latt(3,3)*dble(iz) + pos(3,j,i)! Cartisian
    pos_frac = latt_inv(:,1)*pos_cart(1) + latt_inv(:,2)*pos_cart(2) + latt_inv(:,3)*pos_cart(3) ! fractional
    if ( pos_frac(1)>0.d0-1.2345d-6 .and. pos_frac(1)<1.d0-1.2345d-6 & 
    & .and.  pos_frac(2)>0.d0-1.2345d-6 .and. pos_frac(2)<1.d0-1.2345d-6 &
    & .and.  pos_frac(3)>0.d0-1.2345d-6 .and. pos_frac(3)<1.d0-1.2345d-6) then
      write(9,'(3F12.5)',advance='no') pos_cart
      write(9,*) 'T T T'
      natom_type(i) = natom_type(i) + 1
    endif
  enddo
  do j = 1,natom(i)
    pos_cart(1) = latt(1,1)*dble(ix) + latt(1,2)*dble(iy) + latt(1,3)*dble(iz) + pos(1,j,i+ntype/2)! Cartisian
    pos_cart(2) = latt(2,1)*dble(ix) + latt(2,2)*dble(iy) + latt(2,3)*dble(iz) + pos(2,j,i+ntype/2)! Cartisian
    pos_cart(3) = latt(3,1)*dble(ix) + latt(3,2)*dble(iy) + latt(3,3)*dble(iz) + pos(3,j,i+ntype/2)! Cartisian
    pos_frac = latt_inv(:,1)*pos_cart(1) + latt_inv(:,2)*pos_cart(2) + latt_inv(:,3)*pos_cart(3) ! fractional
    if ( pos_frac(1)>0.d0-1.2345d-6 .and. pos_frac(1)<1.d0-1.2345d-6 & 
    & .and.  pos_frac(2)>0.d0-1.2345d-6 .and. pos_frac(2)<1.d0-1.2345d-6 &
    & .and.  pos_frac(3)>0.d0-1.2345d-6 .and. pos_frac(3)<1.d0-1.2345d-6) then
      write(9,'(3F12.5)',advance='no') pos_cart
      write(9,*) 'F F F'
      natom_type(i) = natom_type(i) + 1
    endif
  enddo
  enddo
  enddo
  enddo !j
  enddo !i
  close(9)

end program solid
