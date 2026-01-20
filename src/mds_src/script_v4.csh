#!/bin/csh
set sluschipath = `grep sluschipath ~/.sluschi.rc | cut -d'=' -f2`
set path_src = $sluschipath/mds_src/

@ l = `find . -name OUTCAR.gz | wc -l`
if ( $l > 0 ) then
  find . -name OUTCAR.gz | xargs gunzip
endif
@ l = `find . -name OUTCAR | wc -l`
if ( $l > 0 ) then
  $sluschipath/collect.csh OUTCAR
endif
if ( $1 == 0 ) then
  set n1 = `grep POSI OUTCAR_collect | wc -l`
  @ n =  $n1 / 80
else
  set n = $1
endif
set nlines = `awk 'BEGIN{print '$n'*3}'`
grep 'Primitive cell' OUTCAR_collect -A  7 | grep -v cell | grep -v latt | grep -v '^$' | grep -v '\-\-' | sed 's/\-/ \-/g' | tail -$nlines > latt
@ n_latt = `cat latt | wc -l`
#echo $n_latt
if ( $n_latt == 0 ) then 
  echo "trying another way to generate latt"
  grep "k-points in units of" OUTCAR_collect -B9 | grep lattice -A3 | grep -v lattice | grep -v '^$' | grep -v '\-\-' | sed 's/\-/ \-/g' | tail -$nlines > latt; 
endif


grep NIONS OUTCAR_collect | head -2 | tail -1 | awk '{print $12}' > natom
set natom = `cat natom`
set param1 = `awk 'BEGIN{print '$natom'+1}'`
set param2 = `awk 'BEGIN{print '$natom'*80*'$n'}'`
grep POSI OUTCAR_collect -A$param1 | grep -v POSI | grep -v '\-\-' | tail -$param2 > pos

if ( -e POTCAR ) then
  set l = `cat POTCAR | wc -l`
else
  @ l = 0
endif
if ( $l == 0 ) then
  head -10000 OUTCAR_collect > POTCAR
endif
set n_elm = `grep POMASS POTCAR | grep ZVAL | wc -l`
echo $n_elm > param
#head -7 POSCAR | tail -1 >> param
grep 'ions per type' OUTCAR_collect | head -1 | cut -d"=" -f2  >> param
grep POMASS POTCAR | grep ZVAL | awk '{print $3}' | cut -d';' -f1 >> param
grep POTIM OUTCAR_collect | head -2 | grep -v use | head -1 | awk '{print $3}'a >> param
cat natom >> param
grep PAW POTCAR | grep -v TIT | grep -v LPAW | grep -v radia | head -$n_elm | sed 's/POTCAR://' | awk '{print $2}' >> param

set nlines = `awk 'BEGIN{print '$n'*80}'`
grep 'energy  without entropy=' OUTCAR_collect | awk '{print $4}' | tail -$nlines > for_avg
$sluschipath/avg_std.x 
echo "MD steps |   C1   |   C2   | Block size | C3 | Average Potential Energy | Standard Error"
cat avg_std.out
set E = `head -1 avg_std.out | awk '{print $6}'`
set E2 = `head -2 avg_std.out | tail -1 | awk '{print $6}'`
set E3 = `head -3 avg_std.out | tail -1 | awk '{print $6}'`
#set E = `grep 'energy  with' OUTCAR_collect -B2 | tail -1 | awk '{print $4}'`
set nlines = `awk 'BEGIN{print '$n'*80}'`
grep 'free  energy   TOTEN' OUTCAR_collect | awk '{print $5}' | tail -$nlines > for_avg
$sluschipath/avg_std.x 
echo "MD steps |   C1   |   C2   | Block size | C3 | Average Potential Energy | Standard Error"
cat avg_std.out
set F = `head -1 avg_std.out | awk '{print $6}'`
set T = `grep TEBEG OUTCAR_collect | tail -1 | awk '{print $3}' | cut -d';' -f1`
set S = `awk -v var1="$T" -v var2="$natom" 'BEGIN{print ('$E'-('$F'))/var1*96485}'`
echo "Potential Energy with Electronic Entropy | Potential Energy | Temperature | Electronic Entropy"
echo $F, $E, $T, $S
echo $S >> param
echo $E >> param
echo $E2 >> param
echo $E3 >> param
echo $F >> param

grep vol OUTCAR_collect | grep ion | awk '{print $5}' | tail -$n > for_avg
$sluschipath/avg_std.x 
echo "MD steps | C1 | C2 | Block size | C3 | Average Volume | Standard Error"
cat avg_std.out | awk '{print $1,$2,$3,$4,$5,$6*'$natom',$7*'$natom'}'
set V = `head -1 avg_std.out | awk '{print $6*'$natom'}'`
set V2 = `head -2 avg_std.out | tail -1 | awk '{print $6*'$natom'}'`
set V3 = `head -3 avg_std.out | tail -1 | awk '{print $6*'$natom'}'`
echo $V >> param
echo $V2 >> param
echo $V3 >> param

#grep Ela OUTCAR_collect | awk '{sum+=$4;} END{print sum;}' | awk '{print Total CPU Hours Spent: $1*48/3600}'
grep Ela OUTCAR_collect | awk '{ sum += $4 } END { printf "Total Physical Hours Spent: %.3f\n", sum / 3600 }'

set n = $1
grep POTIM OUTCAR_collect | grep step| tail -$n | awk '{print $3}' > step
