program GetPos
implicit none

  integer,parameter::j=1
  integer n
  integer n_iter,n_count
  double precision,allocatable,dimension(:,:)::pos
  double precision,dimension(6)::temp

  integer i
  logical success_flag
  character str*40

  open(7,file='iteration')
  open(10,file='OUTCAR')
  open(11,file='natom')
  open(8,file='position')
  open(9,file='get_pos_success')

  read(7,*) n_iter

  read(11,*) n
  allocate(pos(n,3))

  n_count=0
  success_flag=.false.
  do while (success_flag .eqv. .false.)
    read(10,*,err=147,end=147) str
    if (str(1:8)=="POSITION") n_count=n_count+1
    if (n_count==n_iter) then
      read(10,*) str
      do i=1,n
        read(10,*,err=147,end=147) temp
        pos(i,1:3)=temp(1:3)
      enddo
      success_flag=.true.
    endif
  enddo

  147 continue
  if (success_flag) then
    do i=1,n
      write(8,99) i,j,pos(i,1:3)
    enddo
    99 format(I5,I5,3F15.8)
    write(8,*) 
    write(9,*) '1'
  else
    write(9,*) '0'
  endif

end
