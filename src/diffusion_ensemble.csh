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

@ int = `echo "( ($n/4 + 0.5 ) / 1 )" | bc`
@ n_start = $n_exclude
@ n_end = $n_start + $int
mv $n_end folder_tmp
collect.csh OUTCAR $n_start
python $sluschipath/diffusion.py
mv diffusion.pdf diffusion_1.pdf
mv diffusion_avg.pdf diffusion_avg_1.pdf 
mv diffusion_coef.png diffusion_coef_1.png

mv folder_tmp $n_end
@ n_start = $n_end
@ n_end = $n_start + $int
mv $n_end folder_tmp
collect.csh OUTCAR $n_start
python $sluschipath/diffusion.py
mv diffusion.pdf diffusion_2.pdf
mv diffusion_avg.pdf diffusion_avg_2.pdf 
mv diffusion_coef.png diffusion_coef_2.png

mv folder_tmp $n_end
@ n_start = $n_end
@ n_end = $n_start + $int
mv $n_end folder_tmp
collect.csh OUTCAR $n_start
python $sluschipath/diffusion.py
mv diffusion.pdf diffusion_3.pdf
mv diffusion_avg.pdf diffusion_avg_3.pdf 
mv diffusion_coef.png diffusion_coef_3.png

mv folder_tmp $n_end
@ n_start = $n_end
@ n_end = $n_start + $int
mv $n_end folder_tmp
collect.csh OUTCAR $n_start
python $sluschipath/diffusion.py
mv diffusion.pdf diffusion_4.pdf
mv diffusion_avg.pdf diffusion_avg_4.pdf 
mv diffusion_coef.png diffusion_coef_4.png

mv folder_tmp $n_end
