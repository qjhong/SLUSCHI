#!/bin/csh

rm param step latt pos
set sluschipath = `grep sluschipath ~/.sluschi.rc | cut -d'=' -f2`
set path_src = $sluschipath/mds_lmp_src/

set n_elms = `head -1 param*`
@ n = $n_elms + 1
@ n_iter = `cat latt* | wc -l`
@ n_iter = $n_iter / 3
set T = `ls pos_* | cut -d'_' -f3`
set Svib1 = `grep S_p_1 entropy.out -A$n | tail -$n_elms | xargs`
set Svib2 = `grep S_p_2 entropy.out -A$n | tail -$n_elms | xargs`
set E = `tail -7 param* | head -1`
set E1 = `tail -6 param* | head -1`
set E2 = `tail -5 param* | head -1`
set V = `tail -3 param* | head -1`
set V1 = `tail -2 param* | head -1`
set V2 = `tail -1 param* | head -1`
set Selec = `tail -8 param* | head -1`
set Sconf = `grep med entropy.out -A5 | grep -v med | grep -v '\-\-' | grep -v '^$' | xargs `
set Sconf_min = `grep min_S_vec entropy.out -A3 | grep -v min | grep -v '\-\-' | grep -v '^$' | xargs `

echo $Sconf > Sconf.txt
echo $Sconf_min > Sconf_min.txt

echo '===== ANALYSIS STARTED ====='
echo 'Curent folder is ...'
pwd
echo 'inline: ' $T $n_iter $V $V1 $V2 $E $E1 $E2 $Svib1 $Svib2 $Selec $Sconf $Sconf_min
echo 'temperature: ' $T ' K'
echo 'iteration: ' $n_iter ' MD runs, each 80 ionic steps.'
echo 'volume: ' $V $V1 $V2 ' Ang^3 per atom. 0 for lammps results. please calculate from lammps output.'
echo 'energy: ' $E $E1 $E2 ' eV per atom. 0 for lammps results. please calculate from lammps output.'
echo 'Svib: ' $Svib1 ' J/K/mol atom. 1,2,...,n. Do NOT use this value.'
echo 'Svib: ' $Svib2 ' Constrained by ideal gas entropy. Use this value, not the line above.'
echo 'Selec: ' $Selec ' J/K/mol atom. 0 for lammps results.'
echo 'Sconf: ' $Sconf ' J/K/mol atom, 1-1, 1-2, ..., 2-1, 2-2,..., n-1, n-2,... n-n'
echo 'Sconf_min: ' $Sconf_min ' J/K/mol atom, 1-1, 1-2, ..., 2-1, 2-2,..., n-1, n-2,... n-n'
echo '---- Sconf ----'
echo 'I suggest that you use these values...'
python $path_src/Sconf.py
echo '---- Sconf ----'
echo '===== ANALYSIS COMPLETED ====='
