#!/bin/csh

if ( -e anal_time ) rm anal_time

@ temp = 0
@ temp_end = 10000

while ( $temp < $temp_end )
  if ( -d ../$temp ) then
    cd ../$temp/$temp
    @ iter = 0
    @ iter_end = 10000
    while ( $iter < $iter_end )
      if ( -d $iter ) then
        cd $iter
        @ step = 1
        while ( -d $step ) 
          cd $step
          if ( -e OUTCAR.gz ) gunzip OUTCAR.gz 
          if ( -e OUTCAR ) then
            set time = `tail OUTCAR | grep Ela | awk '{print $4}'`
            @ nodes = `head -4 OUTCAR | tail -1 | awk '{print $3}'`
            echo $temp $iter $step $nodes $time >> ../../../../timecollect/anal_time
            gzip OUTCAR
          endif
          cd ..
          sleep 1
          @ step = $step + 1
        end
        cd ..
      endif
      @ iter = $iter + 1
    end
    cd ../../timecollect
  endif
  @ temp = $temp + 1
end

