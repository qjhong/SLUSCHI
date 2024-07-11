#!/bin/csh

set n_elms = `head -1 param`
@ n = $n_elms + 1
@ n_iter = `cat latt | wc -l`
@ n_iter = $n_iter / 3
set T = `ls pos_* | cut -d'_' -f3`
set Svib1 = `grep S_p_1 entropy.out -A$n | tail -$n_elms | xargs`
set Svib2 = `grep S_p_2 entropy.out -A$n | tail -$n_elms | xargs`
set E = `tail -7 param | head -1`
set E1 = `tail -6 param | head -1`
set E2 = `tail -5 param | head -1`
set V = `tail -3 param | head -1`
set V1 = `tail -2 param | head -1`
set V2 = `tail -1 param | head -1`
set Selec = `tail -8 param | head -1`
set Sconf = `grep med entropy.out -A5 | grep -v med | grep -v '\-\-' | grep -v '^$' | xargs `
set Sconf_min = `grep min_S_vec entropy.out -A3 | grep -v min | grep -v '\-\-' | grep -v '^$' | xargs `
echo 'inline: ' $T $n_iter $V $V1 $V2 $E $E1 $E2 $Svib1 $Svib2 $Selec $Sconf $Sconf_min
echo 'temperature: ' $T ' K'
echo 'iteration: ' $n_iter ' MD runs, each 80 ionic steps'
echo 'volume: ' $V $V1 $V2 ' Ang^3 per atom. take the median.'
echo 'energy: ' $E $E1 $E2 ' eV per atom. take the median.'
echo 'Svib: ' $Svib1 ' J/K/mol atom. 1,2,...,n'
echo 'Svib: ' $Svib2 ' Constrained by ideal gas entropy'
echo 'Selec: ' $Selec ' J/K/mol atom'
echo 'Sconf: ' $Sconf ' J/K/mol atom, 1-1, 1-2, ..., 2-1, 2-2,..., n-1, n-2,... n-n'
echo 'Sconf_min: ' $Sconf_min ' J/K/mol atom, 1-1, 1-2, ..., 2-1, 2-2,..., n-1, n-2,... n-n'
