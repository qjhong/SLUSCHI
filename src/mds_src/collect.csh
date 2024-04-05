#!/bin/csh
find . -name 'entropy*' | grep -v entropy.out > foldernames

@ l = `cat foldernames | wc -l`

while ( $l > 0 )
  set folder = `head -1 foldernames`
  cd $folder
  pwd
  /data/qhong7/qhong7/entropy/src/summary.csh
  cd -
  sed -i '1d' foldernames
  @ l = `cat foldernames | wc -l`
end
