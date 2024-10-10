program StopJob
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

      open(12,file='FlagStop.out')
      write(12,*) '1'
      close(12)
    
end program

