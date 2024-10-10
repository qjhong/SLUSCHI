#!/bin/csh

find . -name log.out | grep mds > filenames

@ l = `cat filenames | wc -l`

while ( $l > 0 )

  set dir = `head -1 filenames | sed 's/log.out//'`
  cd $dir
  set dir = `pwd`
  echo 'working on ' $dir
  echo 'Svib entropy/1/2:'
  grep 'Constrained by ideal gas entropy. Use this value, not the line above.' $dir/log.out | awk '{print $2}' | sort -n 
  set Svib = `grep 'Constrained by ideal gas entropy. Use this value, not the line above.' $dir/log.out | awk '{print $2}' | sort -n | tail -2 | head -1`
  echo 'Sconf entropy/1/2:'
  grep 'The pair between element' $dir/log.out | awk '{print $NF}' | sort -n 
  set Sconf = `grep 'The pair between element' $dir/log.out | awk '{print $NF}' | sort -n | tail -2 | head -1`
  set temp = `cat ../phase_temp | cut -d'_' -f2`
  echo 'Summary: '$temp, $Svib, $Sconf
  cd -

  sed -i '1d' filenames
  @ l = `cat filenames | wc -l`

end
