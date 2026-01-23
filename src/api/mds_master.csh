#!/bin/csh

set path_cur = `pwd`
echo $path_cur

set sluschipath = `grep sluschipath ~/.sluschi.rc | cut -d'=' -f2`              
set path_src = $sluschipath/mds_src/
echo $path_src

# 1. python getfile.py <arg1>
set noglob
/home/qhong7/anaconda3/bin/python getfile.py $argv[1]
unset noglob

  $path_src/script_v4.csh 0
  # mkdir and cp
  rm -r entropy
  if ( ! -d entropy ) then
    mkdir entropy
  endif
  date
  cp $path_src/entropy/* entropy
  mv pos param latt step entropy
  date
  #set temp = `grep temp job.in | head -1 | cut -d'=' -f2 | sed 's/ //g' | cut -d'.' -f1`
  set temp = `grep TEBEG OUTCAR_collect | tail -1 | cut -d'.' -f1 | awk '{print $3}'`
  echo $temp
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
  rm sl*out
  # sbatch jobsub_master

module load matlab

set filename = $solid\_$temp
echo $filename

cp latt latt_$filename
cp step step_$filename
cp param param_$filename
cp pos pos_$filename

matlab -r "run('main.m'); exit;" > entropy.out

cd ..

$sluschipath/mds_src/collect.csh
