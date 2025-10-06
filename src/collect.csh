#!/bin/csh

if ( $#argv >= 2 ) then
  @ iter = $2 + 1
else
  @ iter = 1
endif

if ( -e $1_collect ) rm $1_collect

while ( -d $iter )
  cat $iter/$1 >> $1_collect
  @ iter = $iter + 1
end
