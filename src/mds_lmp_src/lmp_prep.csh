#!/bin/csh

set sluschipath = `grep sluschipath ~/.sluschi.rc | cut -d'=' -f2`

if ( ! -e POSCAR_0 ) then
  python lmp_pos.py
endif

#@ natoms = `tail -1 lmp.dump | awk '{print $1}'`
set datafile = `grep '^read_data' ../in* | awk '{print $2}'`
@ natoms = `head -3 ../$datafile | tail -1 | awk '{print $1}'`
set step = `grep timestep ../in* | awk '{print $2}'`
@ nelms = `grep mass ../in* | wc -l`
rm pos latt step param

echo $nelms > param
head -7 POSCAR_0 | tail -1 >> param
grep mass ../in* | awk '{print $3}' >> param
echo $step >> param
echo $natoms >> param
@ iter = 0
while ( $iter < $nelms )
  @ iter = $iter + 1
  echo elm_$iter >> param
end
echo 0.0 >> param
echo 0.0 >> param
echo 0.0 >> param
echo 0.0 >> param
echo 0.0 >> param
echo 0.0 >> param
echo 0.0 >> param
echo 0.0 >> param

@ iter = $1
while ( -e POSCAR_$iter )

  head -5 POSCAR_$iter | tail -3 | awk '{print $1,$2,$3,0,0,0}' >> latt
  tail -$natoms POSCAR_$iter | awk '{print $1,$2,$3,0,0,0}' >> pos
  echo $step >> step
  @ iter = $iter + 1

end

#rm POSCAR*
