#!/bin/csh

cp anal_time tmp

@ l = `cat tmp | wc -l`
@ m = 1
set sum = 0
while ( $m < $l )
  set ncore = `head -$m tmp | tail -1 | awk '{print $4}'`
  set time = `head -$m tmp | tail -1 | awk '{print $5}'`
  set sum = `awk 'BEGIN{print '$sum'+'$ncore'*'$time'}'`
  echo $sum
  @ m = $m + 1
end

#set sum = `echo "$sum * 80.0 / 3600.0" | awk '{printf "%.2f\n", $0}'`
set sum = `printf "%.10f / 3600.0\n" $sum | bc -l`
echo $sum

