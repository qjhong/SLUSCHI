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
set ORIG_SCRIPT = "${path_src}/master.csh"
echo "=== SLUSCHI MDS ENTROPY ANALYSIS START: ${TIMESTAMP} ==="
echo "Running: ${ORIG_SCRIPT}"

find . -name Dir_VolSearch > foldernames                                              
                                                                                
@ l = `cat foldernames | wc -l`                                                 
                                                                                
while ( $l > 0 )                                                                
  set folder = `head -1 foldernames`                                            
  #set foldername = `echo $folder | cut -d'/' -f1`
  cd $folder                                                                    
  pwd                                                                           
  # n files
  @ n = 1
  while ( -d $n )
    @ n = $n + 1
  end
  @ n = $n - 1
# run 1
  # run script
  if ( -e n_iter_exclude ) then
    @ n_exclude = `cat n_iter_exclude`
  else
    @ n_exclude = 50
  endif
  @ n = $n - $n_exclude
  echo $n folders will be used to run analysis
  $path_src/script_v4.csh $n
  # mkdir and cp
  rm -r entropy
  if ( ! -d entropy ) then
    mkdir entropy
  endif
  cp $path_src/entropy/* entropy
  mv pos param latt step entropy
  set temp = `grep temp job.in | head -1 | cut -d'=' -f2 | sed 's/ //g' | cut -d'.' -f1`
  #echo $temp
  set liquid = `grep thmexp_liq job.in | cut -d'=' -f2 | sed 's/ //'`
  cd entropy
  set solid = s
  if ( $liquid == 1 ) then
    set solid = l
  endif
  sed -i "s/replace_here/$solid\_$temp/" jobsub_master
  sed -i "s/replace_here/$solid\_$temp/" main.m
  sed -i "s|replace_folder_here|$path_src|" main.m
  # run jobsub
  rm sl*out
  # sbatch jobsub_master

module load matlab

set filename = $solid\_$temp
#echo $filename

cp latt latt_$filename
cp step step_$filename
cp param param_$filename
cp pos pos_$filename

matlab -r "run('main.m'); exit;" > entropy.out

  # go back
  cd -
# run 2
  # run script
  if ( $n_exclude < 50 ) then
    @ n = $n - $n_exclude
  else
    @ n = $n - 50
  endif
  echo $n folders will be used to run analysis
  $path_src/script_v4.csh $n
  # mkdir and cp
  rm -r entropy1
  if ( ! -d entropy1 ) then
    mkdir entropy1
  endif
  cp $path_src/entropy/* entropy1
  mv pos param latt step entropy1
  set temp = `grep temp job.in | head -1 | cut -d'=' -f2 | sed 's/ //g' | cut -d'.' -f1`
  echo $temp
  set liquid = `grep thmexp_liq job.in | cut -d'=' -f2 | sed 's/ //'`
  cd entropy1
  set solid = s
  if ( $liquid == 1 ) then
    set solid = l
  endif
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

  # go back
  cd -
# run 3
  # run script
  if ( $n_exclude < 50 ) then
    @ n = $n - $n_exclude
  else
    @ n = $n - 50
  endif
  echo $n folders will be used to run analysis
  $path_src/script_v4.csh $n
  # mkdir and cp
  rm -r entropy2
  if ( ! -d entropy2 ) then
    mkdir entropy2
  endif
  cp $path_src/entropy/* entropy2
  mv pos param latt step entropy2
  set temp = `grep temp job.in | head -1 | cut -d'=' -f2 | sed 's/ //g' | cut -d'.' -f1`
  echo $temp
  set liquid = `grep thmexp_liq job.in | cut -d'=' -f2 | sed 's/ //'`
  cd entropy2
  set solid = s
  if ( $liquid == 1 ) then
    set solid = l
  endif
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

  # go back
  cd $path_cur
  # /home/qhong7/data/qhong7/entropy/src/summary.csh                              
  sed -i '1d' foldernames                                                       
  @ l = `cat foldernames | wc -l`                                               
end                                                                             

# Final marker line
set EXITCODE=0
echo ""
echo "=== SLUSCHI MDS ENTROPY ANALYSIS DONE: status=OK exit=${EXITCODE} ==="

