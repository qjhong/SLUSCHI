#!/bin/csh

find . -name Dir_VolSearch > files

@ l = `cat files | wc -l`

while ( $l > 0 )
  set dir = `head -1 files`
  cd $dir
  diffusion.csh > diffusion.out 
  cd -
  sed -i '1d' files
  @ l = `cat files | wc -l`
end
