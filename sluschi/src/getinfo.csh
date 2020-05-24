#!/bin/csh

cat $2 | cut -d'#' -f1 | grep "^ *$1 *=" > /dev/null
if ( $? == 0 ) then
  cat $2 | cut -d'#' -f1 | grep "^ *$1 *=" | cut -d= -f2
  exit 0
endif

exit 1
