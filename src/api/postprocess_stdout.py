#!/usr/bin/env python3
# postprocess_stdout.py
# Parse the raw stdout.log from the csh pipeline and produce a pretty human-readable summary
# and a machine-readable results.json.
#
# Usage (called by wrapper):
# python3 postprocess_stdout.py --stdout stdout.log --stderr stderr.log --pretty stdout_pretty.txt --json results.json

import argparse
import json
import os
import re
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument("--stdout", required=True)
parser.add_argument("--stderr", required=True)
parser.add_argument("--pretty", required=True)
parser.add_argument("--json", required=True)
parser.add_argument("--timestamp", required=False, default=None)
args = parser.parse_args()

with open(args.stdout, "r", errors="replace") as f:
    raw = f.read()

with open(args.stderr, "r", errors="replace") as f:
    rawerr = f.read()

# Helpers to extract values
def find_first(pattern, text, flags=0):
    m = re.search(pattern, text, flags)
    return m.group(1) if m else None

def find_all(pattern, text, flags=0):
    return re.findall(pattern, text, flags)

# 1) Extract Number of Atoms per type line
# look for 'Number of Atoms:' or 'Number of atoms' lines
natom_list = None
m = re.search(r'Number of Atoms:?\s*\[?([0-9,\s]+)\]?', raw)
if m:
    nums = m.group(1)
    try:
        natom_list = [int(x.strip()) for x in nums.split(",") if x.strip()]
    except Exception:
        natom_list = None

# 2) natoms single value: 'natoms:\s+(\d+)'
natoms = find_first(r'natoms:\s*([0-9]+)', raw)
if natoms:
    natoms = int(natoms)
else:
    # try other forms
    n1 = find_first(r'set natom = `cat natom`\s*', raw)
    # fallback: sum natom_list
    if natom_list:
        natoms = sum(natom_list)
    else:
        natoms = None

# 3) Total MD length: look for 'Total MD length: ... fs'
total_md_len = find_first(r'Total MD length:\s*[\t ]*([\d\.]+)\s*fs', raw)
if total_md_len:
    try:
        total_md_len = float(total_md_len)
    except:
        total_md_len = None

# 4) Extract per-element blocks
# The script prints blocks like:
# Element # 1
# ...
# Total MD length: \t\t8890.80 fs
# Diffusion coefficient is: \t1.63e-04 Ang^2/fs or 1.63e-05 cm^2/s
# R2 of the linear fitting is: \t0.91
# The uncertainty of the slope is: \t2.17e-07

elements = []
# find all occurrences of "Element # <n>" and capture following text up to next 'Element #' or end
blocks = re.split(r'\n\s*Element\s*#\s*\d+\s*\n', raw)
# The first split part may be header; subsequent parts correspond to element blocks only if 'Element #' present
# Instead use finditer to get start indices and then slice
it = list(re.finditer(r'Element\s*#\s*(\d+)', raw))
for idx, m in enumerate(it):
    el_index = int(m.group(1))
    start = m.end()
    end = it[idx+1].start() if idx+1 < len(it) else len(raw)
    block = raw[start:end]

    # extract diffusion coefficient (Ang^2/fs)
    d_val = find_first(r'Diffusion coefficient is:\s*([-\d.eE\+]+)\s*Ang\^2/fs', block)
    if not d_val:
        # sometimes printed with tabs and wording, try alternative
        d_val = find_first(r'Diffusion coefficient is:\s*([-\d.eE\+]+)', block)
    try:
        d_val = float(d_val) if d_val else None
    except:
        d_val = None

    # R2
    r2 = find_first(r'R2 of the linear fitting is:\s*([-\d.eE\+]+)', block)
    try:
        r2 = float(r2) if r2 else None
    except:
        r2 = None

    # uncertainty
    unc = find_first(r'The uncertainty of the slope is:\s*([-\d.eE\+]+)', block)
    try:
        unc = float(unc) if unc else None
    except:
        unc = None

    elements.append({
        "element_index": el_index,
        "D_A2_per_fs": d_val,
        "D_cm2_s": (d_val/10.0) if d_val is not None else None,  # conversion: Å^2/fs -> cm^2/s approximate scaling used earlier
        "R2": r2,
        "uncertainty": unc
    })

# If no element blocks found, try to find single-line diffusion coefficients across file
if not elements:
    # pattern: "Diffusion coefficient is:    1.63e-04 Ang^2/fs or 1.63e-05 cm^2/s"
    all_d = find_all(r'Diffusion coefficient is:\s*([-\d.eE\+]+)\s*Ang\^2/fs', raw)
    all_r2 = find_all(r'R2 of the linear fitting is:\s*([-\d.eE\+]+)', raw)
    all_unc = find_all(r'The uncertainty of the slope is:\s*([-\d.eE\+]+)', raw)
    for i, dv in enumerate(all_d):
        try:
            dvf = float(dv)
        except:
            dvf = None
        r2v = float(all_r2[i]) if i < len(all_r2) else None
        uncv = float(all_unc[i]) if i < len(all_unc) else None
        elements.append({
            "element_index": i+1,
            "D_A2_per_fs": dvf,
            "D_cm2_s": (dvf/10.0) if dvf is not None else None,
            "R2": r2v,
            "uncertainty": uncv
        })

# 5) Extract block-average tables (avg_std.out style) if present
# We look for repeated table lines like: 14400    69.69089     2.31040 ...
# We'll capture lines that look numeric and have 6-7 columns
block_averages = []
for line in raw.splitlines():
    # match lines with 6-7 numeric-ish columns
    parts = re.findall(r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?', line)
    if len(parts) >= 6:
        # Heuristically pick lines that look like the avg_std rows: first field integer MD steps
        try:
            mdsteps = int(float(parts[0]))
            # take relevant columns
            C1 = float(parts[1])
            C2 = float(parts[2])
            block_size = parts[3]
            avg_pot = float(parts[5])
            stderr = float(parts[6]) if len(parts) > 6 else None
            block_averages.append({
                "mdsteps": mdsteps,
                "C1": C1,
                "C2": C2,
                "block_size": block_size,
                "avg_potential": avg_pot,
                "stderr": stderr
            })
        except Exception:
            continue

# 6) Potential energy summary extraction
pot_E = find_first(r'Potential Energy with Electronic Entropy.*\n([-\d\.,\s]+)', raw, flags=re.IGNORECASE)
# fallback: look for a comma-separated line like '-6.96937, -6.96936, 1400.0, 0.000689179'
potline = find_first(r'(-?\d+\.\d+,\s*-?\d+\.\d+,\s*\d+\.?\d*,\s*\d+\.\d+)', raw)
if potline and not pot_E:
    pot_E = potline

# 7) File artifacts (if lines mention OUTCAR_collect or diffusion_coef.png)
artifacts = []
if "OUTCAR_collect" in raw or os.path.exists("OUTCAR_collect"):
    artifacts.append("OUTCAR_collect")
# find png/pdf file mentions in stdout (diffusion_coef.png etc)
for ext in [".png", ".pdf"]:
    found = re.findall(r'[\w\-/]*\w+' + re.escape(ext), raw)
    for f in found:
        if f not in artifacts:
            artifacts.append(os.path.basename(f))

# Build summary dict
summary = {
    "created_at": args.timestamp or datetime.utcnow().isoformat() + "Z",
    "natoms_per_type": natom_list,
    "natoms": natoms,
    "total_md_length_fs": total_md_len,
    "elements": sorted(elements, key=lambda x: x.get("element_index") or 0),
    "block_averages": block_averages,
    "potential_energy_line": pot_E,
    "artifacts": artifacts,
    "notes": "Generated by postprocess_stdout.py"
}

# Write results.json
with open(args.json, "w") as jf:
    json.dump(summary, jf, indent=2, sort_keys=False)

# Create a pretty human-readable stdout summary
lines = []
lines.append(f"=== SLUSCHI RUN SUMMARY ({summary['created_at']}) ===")
if natom_list:
    lines.append(f"Number of atoms per type: {', '.join(str(x) for x in natom_list)}")
if natoms:
    lines.append(f"Total number of atoms: {natoms}")
if total_md_len:
    lines.append(f"Total MD length (fs): {total_md_len}")
if summary["block_averages"]:
    lines.append("")
    lines.append("-- Block averages (potential energy) --")
    lines.append("{:8s} {:10s} {:10s} {:8s} {:12s} {:10s}".format("MDsteps","C1","C2","Block","AvgPot","StdErr"))
    for b in summary["block_averages"]:
        lines.append("{:8d} {:10.5f} {:10.5f} {:8s} {:12.5f} {:10s}".format(
            b["mdsteps"], b["C1"], b["C2"], str(b["block_size"]),
            b["avg_potential"], str(b["stderr"])
        ))
if summary["potential_energy_line"]:
    lines.append("")
    lines.append("-- Potential energy (with electronic entropy) --")
    lines.append(str(summary["potential_energy_line"]))

if summary["elements"]:
    lines.append("")
    lines.append("-- Per-element diffusion results --")
    lines.append("{:3s} {:15s} {:15s} {:8s} {:12s}".format("#", "D (Å^2/fs)", "D (cm^2/s)", "R2", "uncertainty"))
    for e in summary["elements"]:
        lines.append("{:3d} {:15s} {:15s} {:8s} {:12s}".format(
            e.get("element_index") or 0,
            ("{:.3e}".format(e["D_A2_per_fs"]) if e.get("D_A2_per_fs") is not None else "NA"),
            ("{:.3e}".format(e["D_cm2_s"]) if e.get("D_cm2_s") is not None else "NA"),
            (("{:.2f}".format(e["R2"])) if e.get("R2") is not None else "NA"),
            (("{:.3e}".format(e["uncertainty"])) if e.get("uncertainty") is not None else "NA")
        ))

if artifacts:
    lines.append("")
    lines.append("Artifacts detected: " + ", ".join(artifacts))

# include a short note about raw logs
lines.append("")
lines.append("Raw logs saved to stdout.log and stderr.log in the run directory.")
lines.append("Full raw output is available for advanced inspection.")

pretty_text = "\n".join(lines) + "\n"

with open(args.pretty, "w") as pf:
    pf.write(pretty_text)

# done
