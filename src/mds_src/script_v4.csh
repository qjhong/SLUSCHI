#!/bin/csh
gunzip */OUTCAR.gz
~/sluschi/sluschi_latest/src/collect.csh OUTCAR
set n = $1
set nlines = `awk 'BEGIN{print '$n'*3}'`
grep 'Primitive cell' OUTCAR_collect -A  7 | grep -v cell | grep -v latt | grep -v '^$' | grep -v '\-\-' | sed 's/\-/ \-/g' | tail -$nlines > latt
@ n_latt = `cat latt | wc -l`
echo $n_latt
if ( $n_latt == 0 ) then 
  echo "try another way to generate latt"
  grep "k-points in units of" OUTCAR_collect -B9 | grep lattice -A3 | grep -v lattice | grep -v '^$' | grep -v '\-\-' | sed 's/\-/ \-/g' | tail -$nlines > latt; 
endif


grep NIONS OUTCAR_collect | head -2 | tail -1 | awk '{print $12}' > natom
set natom = `cat natom`
set param1 = `awk 'BEGIN{print '$natom'+1}'`
set param2 = `awk 'BEGIN{print '$natom'*80*'$n'}'`
grep POSI OUTCAR_collect -A$param1 | grep -v POSI | grep -v '\-\-' | tail -$param2 > pos

grep POMASS POTCAR | wc -l > param
head -7 POSCAR | tail -1 >> param
grep POMASS POTCAR | awk '{print $3}' | cut -d';' -f1 >> param
grep POTIM OUTCAR_collect | head -2 | grep -v use | head -1 | awk '{print $3}'a >> param
cat natom >> param
grep PAW POTCAR | grep -v TIT | grep -v LPAW | grep -v radia | awk '{print $2}' >> param

set nlines = `awk 'BEGIN{print '$n'*80}'`
grep 'energy  without entropy=' OUTCAR_collect | awk '{print $4}' | tail -$nlines > for_avg
~/sluschi/sluschi_latest/src/avg_std.x 
cat avg_std.out
set E = `head -1 avg_std.out | awk '{print $6}'`
set E2 = `head -2 avg_std.out | tail -1 | awk '{print $6}'`
set E3 = `head -3 avg_std.out | tail -1 | awk '{print $6}'`
#set E = `grep 'energy  with' OUTCAR_collect -B2 | tail -1 | awk '{print $4}'`
set nlines = `awk 'BEGIN{print '$n'*80}'`
grep 'free  energy   TOTEN' OUTCAR_collect | awk '{print $5}' | tail -$nlines > for_avg
~/sluschi/sluschi_latest/src/avg_std.x 
cat avg_std.out
set F = `head -1 avg_std.out | awk '{print $6}'`
set T = `grep TEBEG OUTCAR_collect | tail -1 | awk '{print $3}' | cut -d';' -f1`
set S = `awk -v var1="$T" -v var2="$natom" 'BEGIN{print ('$E'-('$F'))/var1*96485}'`
echo $F, $E, $T, $S
echo $S >> param
echo $E >> param
echo $E2 >> param
echo $E3 >> param
echo $F >> param

grep vol OUTCAR_collect | grep ion | awk '{print $5}' | tail -$n > for_avg
~/sluschi/sluschi_latest/src/avg_std.x 
cat avg_std.out | awk '{print $1,$2,$3,$4,$5,$6*'$natom'}'
set V = `head -1 avg_std.out | awk '{print $6*'$natom'}'`
set V2 = `head -2 avg_std.out | tail -1 | awk '{print $6*'$natom'}'`
set V3 = `head -3 avg_std.out | tail -1 | awk '{print $6*'$natom'}'`
echo $V >> param
echo $V2 >> param
echo $V3 >> param

grep Ela OUTCAR_collect | awk '{sum+=$4;} END{print sum;}' | awk '{print $1*48/3600}'

set n = $1
grep POTIM OUTCAR_collect | grep step| tail -$n | awk '{print $3}' > step
