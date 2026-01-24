#!/bin/csh
set sluschipath = `grep sluschipath ~/.sluschi.rc | cut -d'=' -f2`
set path_src = $sluschipath/mds_src/

find . -name 'entropy*' | grep -v entropy.out > foldernames

@ l = `cat foldernames | wc -l`

while ( $l > 0 )
  set folder = `head -1 foldernames`
  cd $folder
  pwd
  $path_src/summary.csh
  cd ..
  sed -i '1d' foldernames
  @ l = `cat foldernames | wc -l`
end
