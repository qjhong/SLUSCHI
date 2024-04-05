#!/bin/csh

if ( -e anal_time ) rm anal_time

cd ../
@ temp = 1
@ iter = 1
@ step = 1
while ( -d $step ) 
  cd $step
  if ( -e OUTCAR.gz ) gunzip OUTCAR.gz 
  if ( -e OUTCAR ) then
    set time = `tail OUTCAR | grep Ela | awk '{print $4}'`
    head -4 OUTCAR | tail -1 | grep mpi > /dev/null
    if ($? == 1) then
      @ nodes = `head -4 OUTCAR | tail -1 | awk '{print $3}'`
    else
      @ nodes = `head -4 OUTCAR | tail -1 | awk '{print $2}'`
    endif
    echo $temp $iter $step $nodes $time >> ../timecollect/anal_time
  endif
  cd ..
  @ step = $step + 1
end

if ( -e OUTCAR_collect.gz ) gunzip OUTCAR_collect.gz 
set md_length = `grep POTIM OUTCAR_collect  | grep ion | awk '{sum += $3*80} END {print sum}'`
set n_md = `grep POTIM OUTCAR_collect  | grep ion | wc -l`
echo $md_length $n_md > timecollect/time
