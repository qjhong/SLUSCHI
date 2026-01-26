#!/bin/csh

set path_cur = `pwd`
echo $path_cur

set sluschipath = `grep sluschipath ~/.sluschi.rc | cut -d'=' -f2`              
set path_src = $sluschipath/mds_src/
echo $path_src

echo ""
echo "========================================"
echo "        SLUSCHI by Qi-Jun Hong"
echo "========================================"
echo ""
set TIMESTAMP = `date -u +%Y%m%dT%H%M%SZ`
set ORIG_SCRIPT = "${path_src}/mds_master.csh"
echo "=== SLUSCHI MDS ENTROPY ANALYSIS START: ${TIMESTAMP} ==="
echo "Running: ${ORIG_SCRIPT} ${ARG}"

# 1. python getfile.py <arg1>
set noglob
/home/qhong7/anaconda3/bin/python $sluschipath/api/getfile.py $argv[1]
unset noglob

  $path_src/script_v4.csh 0
  # mkdir and cp
  if ( -d entropy ) rm -r entropy
  if ( ! -d entropy ) then
    mkdir entropy
  endif
  cp $path_src/entropy/* entropy
  mv pos param latt step entropy
  #set temp = `grep temp job.in | head -1 | cut -d'=' -f2 | sed 's/ //g' | cut -d'.' -f1`
  set temp = `grep TEBEG OUTCAR_collect | tail -1 | cut -d'.' -f1 | awk '{print $3}'`
  #echo $temp
  #set liquid = `grep thmexp_liq job.in | cut -d'=' -f2 | sed 's/ //'`
  cd entropy
  mv main_api.m main.m
  set solid = s
  #if ( $liquid == 1 ) then
  #  set solid = l
  #endif
  sed -i "s/replace_here/$solid\_$temp/" jobsub_master
  sed -i "s/replace_here/$solid\_$temp/" main.m
  sed -i "s|replace_folder_here|$path_src|" main.m
  # run jobsub
  # rm sl*out
  # sbatch jobsub_master

#module load matlab

set filename = $solid\_$temp
#echo $filename

cp latt latt_$filename
cp step step_$filename
cp param param_$filename
cp pos pos_$filename

/home/qhong7/bin/matlab -r "run('main.m'); exit;" > entropy.out

cd ..

$sluschipath/mds_src/collect.csh

# Final marker line
set EXITCODE=0
echo ""
echo "=== SLUSCHI MDS ENTROPY ANALYSIS DONE: status=OK exit=${EXITCODE} ==="

