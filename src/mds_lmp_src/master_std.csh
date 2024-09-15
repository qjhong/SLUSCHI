#!/bin/csh

set path_cur = `pwd`
echo $path_cur

set sluschipath = `grep sluschipath ~/.sluschi.rc | cut -d'=' -f2`              
set path_src = $sluschipath/mds_src/                                            

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
  if ( ! -d entropy ) then
    mkdir entropy
  endif
  cp $path_src/entropy/* entropy
  cp pos param latt step entropy
  set temp = `grep temp job.in | cut -d'=' -f2 | sed 's/ //g' | cut -d'.' -f1`
  echo $temp
  set liquid = `grep thmexp_liq job.in | cut -d'=' -f2 | sed 's/ //'`
  cd entropy
  set solid = s
  if ( $liquid == 1 ) then
    set solid = l
  endif
  sed -i "s/replace_here/$solid\_$temp/" jobsub_master
  sed -i "s/replace_here/$solid\_$temp/" main.m
  # run jobsub
  rm sl*out
  sbatch jobsub_master

  cd ..
  @ n = $n - 50
  echo $n folders will be used to run analysis
  $path_src/script_v4.csh $n
  # mkdir and cp
  if ( ! -d entropy1 ) then
    mkdir entropy1
  endif
  cp $path_src/entropy/* entropy1
  cp pos param latt step entropy1
  set temp = `grep temp job.in | cut -d'=' -f2 | sed 's/ //g' | cut -d'.' -f1`
  echo $temp
  set liquid = `grep thmexp_liq job.in | cut -d'=' -f2 | sed 's/ //'`
  cd entropy1
  set solid = s
  if ( $liquid == 1 ) then
    set solid = l
  endif
  sed -i "s/replace_here/$solid\_$temp/" jobsub_master
  sed -i "s/replace_here/$solid\_$temp/" main.m
  # run jobsub
  rm sl*out
  sbatch jobsub_master

  cd ..
  @ n = $n - 50
  echo $n folders will be used to run analysis
  $path_src/script_v4.csh $n
  # mkdir and cp
  if ( ! -d entropy2 ) then
    mkdir entropy2
  endif
  cp $path_src/entropy/* entropy2
  cp pos param latt step entropy2
  set temp = `grep temp job.in | cut -d'=' -f2 | sed 's/ //g' | cut -d'.' -f1`
  echo $temp
  set liquid = `grep thmexp_liq job.in | cut -d'=' -f2 | sed 's/ //'`
  cd entropy2
  set solid = s
  if ( $liquid == 1 ) then
    set solid = l
  endif
  sed -i "s/replace_here/$solid\_$temp/" jobsub_master
  sed -i "s/replace_here/$solid\_$temp/" main.m
  # run jobsub
  rm sl*out
  sbatch jobsub_master

  # go back
  cd $path_cur
  # /home/qhong7/data/qhong7/entropy/src/summary.csh                              
  sed -i '1d' foldernames                                                       
  @ l = `cat foldernames | wc -l`                                               
end                                                                             

