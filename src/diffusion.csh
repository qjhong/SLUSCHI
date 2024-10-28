#!/bin/csh

set path_cur = `pwd`
echo $path_cur

set sluschipath = `grep sluschipath ~/.sluschi.rc | cut -d'=' -f2`
set path_src = $sluschipath/mds_src/
 
# n files
@ n = 1
while ( -d $n )
  @ n = $n + 1
end
@ n = $n - 1
# run script
if ( -e n_iter_exclude ) then
  @ n_exclude = `cat n_iter_exclude`
else
  @ n_exclude = 50
endif
@ n = $n - $n_exclude
echo $n folders will be used to run analysis
$path_src/script_v4.csh $n
python $sluschipath/diffusion.py
