#!/bin/csh

# usage: script.csh <pattern> [ELEMENT] [PHASE] [TEMP_K] [STEPS]
# example: script.csh mypattern Al liquid 1000 16000

set sluschipath = `grep sluschipath ~/.sluschi.rc | cut -d'=' -f2`

find . -name 'latt*press*' | grep $1 > filenames

@ l = `cat filenames | wc -l`

while ( $l > 0 ) 
  set dir = `head -1 filenames | sed 's/lattice_pressure_history.out//' `
  cd $dir
  pwd
  gunzip */OUTCAR.gz
  collect.csh OUTCAR
  gzip */OUTCAR

# ---------- element detection from POTCAR (step 1) ----------
# ---------- element = chemical formula from POTCAR + OUTCAR_collect ----------
set element = "UnknownElement"

# get counts from OUTCAR_collect (ions per type = ...). trim spaces.
set counts_line = `grep 'ions per type' OUTCAR_collect | head -1 | cut -d"=" -f2 | awk '{print $1}'`
#echo $counts_line
if ( "$counts_line" != "" ) then
  # normalize whitespace and build array of counts
  set counts_line = `echo "$counts_line" | tr -s ' ' | sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//'`
  set counts = ( $counts_line )
else
  set counts = ( )
endif

# get element symbols from POTCAR using your pipeline
# adjust head -$n_elm later by using number of counts or number of symbols
set syms_line = `grep PAW POTCAR | grep -v TIT | grep -v LPAW | grep -v radia | sed 's/POTCAR://' | awk '{print $2}'`
#echo $syms_line
if ( "$syms_line" != "" ) then
  # put into array (one per token/line)
  set syms = ( $syms_line )
else
  set syms = ( )
endif

# determine how many pairs we can make
set nc = $#counts
set ns = $#syms

if ( $nc == 0 && $ns == 0 ) then
  set element = "UnknownElement"
else
  if ( $nc == 0 ) then
    set n_pairs = $ns
  else if ( $ns == 0 ) then
    set n_pairs = $nc
  else
    # use smaller of the two to avoid index mismatch
    if ( $nc < $ns ) then
      set n_pairs = $nc
    else
      set n_pairs = $ns
    endif
  endif

  # build formula
  set formula = ""
  @ i = 1
  while ( $i <= $n_pairs )
    # default count = 1 if missing or zero-length
    set sym = ""
    if ( $i <= $#syms ) then
      set sym = $syms[$i]
    endif

    set cnt = ""
    if ( $i <= $#counts ) then
      set cnt = $counts[$i]
    endif

    # sanitize cnt (remove non-digits)
    if ( "$cnt" != "" ) then
      set cnt = `echo $cnt | sed -E 's/[^0-9].*//'`
    endif
    if ( "$cnt" == "" || "$cnt" == "0" ) then
      set cnt = "1"
    endif

    # only append if symbol is non-empty
    if ( "$sym" != "" ) then
      if ( "$cnt" == "1" ) then
        set formula = "${formula}${sym}"
      else
        set formula = "${formula}${sym}${cnt}"
      endif
    endif

    @ i = $i + 1
  end

  if ( "$formula" != "" ) then
    set element = $formula
  else
    set element = "UnknownElement"
  endif
endif
#echo $element
# ---------- end chemical-formula detection ----------

# ---------- phase detection using getinfo.csh ----------
set phase = "UnknownPhase"

if ( -e job.in ) then
  @ flag_liq = `getinfo.csh thmexp_liq job.in`

  if ( $flag_liq > 0 ) then
    set phase = "liquid"
  else
    set phase = "solid"
  endif
endif
# ---------- end phase detection ----------

# ---------- temperature (K) detection from OUTCAR (step 3) ----------
set tempK = ""
set temp = ""

# 1) Lines containing the word "temperature" and a number
set temp = `grep TEBEG OUTCAR_collect | tail -1 | awk '{print $3}' | cut -d';' -f1`

# 2) Pattern like "1000 K" (number immediately followed by 'K' or ' K')
if ( "$temp" == "" ) then
  set temp = `grep -m1 -Eo '[0-9]+(\\.[0-9]+)?[[:space:]]*K' OUTCAR | sed -E 's/[[:space:]]*[Kk]//g'`
endif

# 3) Pattern like "T=1000" or "T = 1000"
if ( "$temp" == "" ) then
  set temp = `grep -m1 -Eo 'T[[:space:]]*=[[:space:]]*[0-9]+(\\.[0-9]+)?' OUTCAR | sed -E 's/[^0-9\\.]//g'`
endif

# 4) VASP-specific tags like TEBEG / TEEND: "TEBEG = 300.0"
if ( "$temp" == "" ) then
  set temp = `grep -m1 -E 'TEBEG|TEEND' OUTCAR | sed -E 's/.*=[[:space:]]*([0-9]+(\\.[0-9]+)?).*/\\1/'`
endif

# 5) fallback: try to find any number that looks like a temperature near "kinetic energy" lines
if ( "$temp" == "" ) then
  set temp = `grep -m1 -i "kinetic energy" OUTCAR | sed -E 's/.*([0-9]+\\.[0-9]+|[0-9]+).*/\\1/'`
endif

# normalize/round to integer K if we extracted a numeric value
if ( "$temp" != "" ) then
  # remove possible leading plus signs
  set temp = `echo $temp | sed -E 's/^\\+//g'`
  # round to nearest integer using awk
  set tempK = `echo $temp | awk '{printf("%d\n",$1+0.5)}'`
else
  set tempK = "unknownTemp"
endif

# sanity: if tempK is empty or non-numeric, fallback
if ( "$tempK" == "" || "$tempK" == "0" ) then
  set tempK = "unknownTemp"
endif
#echo $tempK
# ---------- end temperature detection ----------

# ---------- steps detection from OUTCAR (step 4) ----------
set steps = ""

# 1) Try to find NSW = <number> (common INCAR/OUTCAR appearance)
set steps = `grep POSI OUTCAR_collect | wc -l`

# 2) Some OUTCARs print "Number of ionic steps:" or "number of steps"
if ( "$steps" == "" ) then
  set steps = `grep -m1 -i -E 'number of (ionic )?steps|number of steps|Number of ionic steps' OUTCAR | sed -E 's/.*([0-9]{1,6}).*/\1/'`
endif

# 3) Try patterns like "Steps =" or "Step =" or lines with "STEP" followed by a number
if ( "$steps" == "" ) then
  set steps = `grep -m1 -E -o '[sS]teps?[[:space:]]*[:=][[:space:]]*[0-9]+' OUTCAR | sed -E 's/[^0-9]*([0-9]+).*/\1/'`
endif

# 4) Fallback: count the number of ionic-step energy prints in OUTCAR.
#    "free  energy   TOTEN" is printed once per ionic step in many VASP OUTCARs.
#    We'll count a few common energy-summary tokens and take the maximum count found.
if ( "$steps" == "" ) then
  set c1 = `grep -c -i "free .*energy" OUTCAR `
  set c2 = `grep -c -i "FREE ENERGIE TOTEN" OUTCAR `
  set c3 = `grep -c -i "energy  without entropy" OUTCAR `
  # choose the largest non-zero count (safer if some patterns appear in other contexts)
  set maxc = 0
  if ( "$c1" != "" && $c1 > $maxc ) set maxc = $c1
  if ( "$c2" != "" && $c2 > $maxc ) set maxc = $c2
  if ( "$c3" != "" && $c3 > $maxc ) set maxc = $c3

  if ( $maxc > 0 ) then
    set steps = $maxc
  endif
endif

# 5) If still not found, try reading INCAR (often has NSW) as a last resort
if ( "$steps" == "" && -e INCAR ) then
  set steps = `grep -m1 -E -o 'NSW[[:space:]]*=[[:space:]]*[0-9]+' INCAR  | sed -E 's/[^0-9]*([0-9]+).*/\1/'`
endif

# final fallback
if ( "$steps" == "" ) then
  set steps = "unknownSteps"
endif

# sanity: remove any leading zeros and ensure it's not empty
if ( "$steps" != "unknownSteps" ) then
  set steps = `echo $steps | sed -E 's/^0+([0-9])/\1/; s/[^0-9].*//'`
  if ( "$steps" == "" ) then
    set steps = "unknownSteps"
  endif
endif
#echo $steps
# ---------- end steps detection ----------
#
#
#
  # try to parse number of atoms from OUTCAR (commonly NIONS or number of ions)
  # fallback to "unknownAtoms" if not found
  set nat = `grep NIONS OUTCAR_collect | head -2 | tail -1 | awk '{print $12}'`
  if ( "$nat" == "" ) then
    # try another pattern (some OUTCARs use "number of ions" or "ions per type" lines)
    set nat = `grep -m1 -i "number of ions" OUTCAR | awk '{print $4}'`
  endif
  if ( "$nat" == "" ) then
    set nat = "unknownAtoms"
  endif

  # construct final name
  set destdir = "/data/qhong7/qhong7/OUTCARs"
  set newname = "OUTCAR_collect_${element}_${phase}_${tempK}K_${steps}steps_${nat}atoms"

  # move/rename. Add suffix if file already exists to avoid overwrite
  if ( -e ${destdir}/${newname} ) then
    # append timestamp to avoid clobbering
    set stamp = `date +%Y%m%d%H%M%S`
    mv OUTCAR_collect ${destdir}/${newname}_${stamp}
  else
    mv OUTCAR_collect ${destdir}/${newname}
  endif
  # ---------- end added code ----------

  cd -

  sed -i '1d' filenames
  @ l = `cat filenames | wc -l`
end
