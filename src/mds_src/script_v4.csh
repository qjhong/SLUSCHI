#!/bin/csh
set sluschipath = `grep sluschipath ~/.sluschi.rc | cut -d'=' -f2`
set path_src = $sluschipath/mds_src/

echo ""
echo "========================================"
echo "        SLUSCHI by Qi-Jun Hong"
echo "========================================"
echo ""
set TIMESTAMP = `date -u +%Y%m%dT%H%M%SZ`
set ORIG_SCRIPT = "${path_src}/script_v4.csh"
# forward argument (n) if provided
# forward argument (n) if provided
if ( $#argv >= 1 ) then
    set ARG = "$argv[1]"
else
    set ARG = 0
endif
echo "=== SLUSCHI ENERGY and VOLUME ANALYSIS START: ${TIMESTAMP} ==="
echo "Running: ${ORIG_SCRIPT} ${ARG}"

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
#echo "MD steps |   C1   |   C2   | Block size | C3 | Average Potential Energy | Standard Error"
echo "MD steps |   C1   |   C2   | Block size | C3 | \033[1;31mAverage Potential Energy [eV/atom]\033[0m | Standard Error [eV/atom]"
cat avg_std.out
set E = `head -1 avg_std.out | awk '{print $6}'`
set E2 = `head -2 avg_std.out | tail -1 | awk '{print $6}'`
set E3 = `head -3 avg_std.out | tail -1 | awk '{print $6}'`
#set E = `grep 'energy  with' OUTCAR_collect -B2 | tail -1 | awk '{print $4}'`
set nlines = `awk 'BEGIN{print '$n'*80}'`
grep 'free  energy   TOTEN' OUTCAR_collect | awk '{print $5}' | tail -$nlines > for_avg
$sluschipath/avg_std.x 
#echo "MD steps |   C1   |   C2   | Block size | C3 | Average Potential Energy | Standard Error"
echo "MD steps |   C1   |   C2   | Block size | C3 | \033[1;31mAverage Potential Energy with Electronic Entropy [eV/atom]\033[0m | Standard Error [eV/atom]"
cat avg_std.out
set F = `head -1 avg_std.out | awk '{print $6}'`
set T = `grep TEBEG OUTCAR_collect | tail -1 | awk '{print $3}' | cut -d';' -f1`
set S = `awk -v var1="$T" -v var2="$natom" 'BEGIN{print ('$E'-('$F'))/var1*96485}'`
echo "=== Summary of Potential Energy ==="
echo "Potential Energy with Electronic Entropy [eV/atom] | Potential Energy [eV/atom] | Temperature [K] | Electronic Entropy [J/K/(mol atom)]"
echo $F, $E, $T, $S
echo $S >> param
echo $E >> param
echo $E2 >> param
echo $E3 >> param
echo $F >> param

grep vol OUTCAR_collect | grep ion | awk '{print $5}' | tail -$n > for_avg
$sluschipath/avg_std.x 
echo "MD jobs (80 MD steps each) | C1 | C2 | Block size | C3 | \033[1;31mAverage Volume [Ang^3/atom]\033[0m | Standard Error [Ang^3/atom]"
cat avg_std.out | awk '{print $1,$2,$3,$4,$5,$6*'$natom',$7*'$natom'}'
set V = `head -1 avg_std.out | awk '{print $6*'$natom'}'`
set V2 = `head -2 avg_std.out | tail -1 | awk '{print $6*'$natom'}'`
set V3 = `head -3 avg_std.out | tail -1 | awk '{print $6*'$natom'}'`
echo $V >> param
echo $V2 >> param
echo $V3 >> param

#grep Ela OUTCAR_collect | awk '{sum+=$4;} END{print sum;}' | awk '{print Total CPU Hours Spent: $1*48/3600}'
grep Ela OUTCAR_collect | awk '{ sum += $4 } END { printf "Total Physical Hours Spent: %.3f\n", sum / 3600 }'

#set n = $1
grep POTIM OUTCAR_collect | grep step| tail -$n | awk '{print $3}' > step

# ---- robust extraction of natoms and element names from param ----
set PARAM_FILE = param

if ( ! -e $PARAM_FILE ) then
    echo "ERROR: param file not found: $PARAM_FILE"
    exit 1
endif

# 1) read the 2nd non-empty line and normalize whitespace -> NATOMS_LINE
set NATOMS_LINE = `awk 'NF {count++; if(count==2){ for(i=1;i<=NF;i++) printf "%s ", $i; exit } }' $PARAM_FILE`
set NATOMS_LINE = `echo $NATOMS_LINE | tr -s ' ' | sed 's/^ //; s/ $//'`

# fallback if empty for some reason
if ( "$NATOMS_LINE" == "" ) then
    set NATOMS_LINE = `awk 'NF{c++; if(c==2){print; exit}}' $PARAM_FILE`
    set NATOMS_LINE = `echo $NATOMS_LINE | tr -s ' ' | sed 's/^ //; s/ $//'`
endif

echo "Number of Atoms: [ $NATOMS_LINE ]"

# 2) count how many numbers are on that line -> NTYPE
set NTYPE = `echo $NATOMS_LINE | wc -w | awk '{print $1}'`

if ( "$NTYPE" == "0" || "$NTYPE" == "" ) then
    echo "ERROR: could not parse natoms from second line of $PARAM_FILE"
    exit 1
endif

# 3) find a block of NTYPE consecutive lines whose first token starts with a letter.
#    This matches blocks like:
#      La
#      Li
#      O
#      Zr_sv
# --- find first block of NTYPE consecutive element-like lines using perl ---
#set TYPE_ATOMS = `perl -e ' $n = shift @ARGV; @L = <>; for ($i=0; $i < @L; $i++) { $ok = 1; for ($j=0; $j < $n; $j++) { $idx = $i + $j; if ($idx >= @L) { $ok = 0; last } ($tok) = split(/\s+/, $L[$idx]); if (!defined $tok || $tok !~ /^[A-Za-z]/) { $ok = 0; last } } if ($ok) { for ($j=0; $j < $n; $j++) { ($tok) = split(/\s+/, $L[$i+$j]); print $tok, "\n" } last } }' $NTYPE $PARAM_FILE`

# Fallback: if perl didn't find anything, try simple token scan (first NTYPE tokens that start with letters)
#if ( "$TYPE_ATOMS" == "" ) then
#    set TYPE_ATOMS = `awk -v n=$NTYPE '{ for(i=1;i<=NF;i++){ if ($i ~ /^[A-Za-z]/) { found[++c] = $i } if (c==n) break } } END { for(i=1;i<=n && i<=c;i++) print found[i] }' $PARAM_FILE`
#endif
#set TYPE_ATOMS = `awk -v n=$NTYPE ' NF { lines[++c] = $0 } END { for (i = 1; i <= c; i++) { ok = 1 for (j = 0; j < n; j++) { idx = i + j if (idx > c) { ok = 0; break } split(lines[idx], a, /[[:space:],;]+/) token = a[1] if (token !~ /^[A-Za-z]/) { ok = 0; break } } if (ok) { for (j = 0; j < n; j++) { split(lines[i+j], b, /[[:space:],;]+/) if (b[1] != "") print b[1] } exit } } } ' $PARAM_FILE`

# 4) fallback: try to find any tokens that look like element names and take the first NTYPE
#if ( "$TYPE_ATOMS" == "" ) then
    set TYPE_ATOMS = `awk -v n=$NTYPE '{ for(i=1;i<=NF;i++){ if ($i ~ /^[A-Za-z]/) { found[++c] = $i } if (c==n) break } } END { for(i=1;i<=n && i<=c;i++) print found[i] }' $PARAM_FILE`
#endif

# 5) Print elements nicely
echo -n "Type of Atoms: [ "
foreach t ($TYPE_ATOMS)
    echo -n "$t "
end
echo "]"
# ------------------------------------------------------------------
#
# Final marker line
set EXITCODE=0
echo ""
echo "=== SLUSCHI ENERGY and VOLUME ANALYSIS DONE: status=OK exit=${EXITCODE} ==="
