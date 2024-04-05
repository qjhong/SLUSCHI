#!/bin/csh

@ iter = 1

if ( -e $1_collect ) rm $1_collect

while ( -d $iter )
  cat $iter/$1 >> $1_collect
  @ iter = $iter + 1
end
