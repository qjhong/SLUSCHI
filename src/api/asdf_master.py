#!/usr/bin/env python3

import subprocess
import sys
import os

def run(cmd, shell=False):
    print(f"Running: {cmd}")
    subprocess.check_call(cmd, shell=shell)

def main():
    # 1. python getfile.py
    run([sys.executable, "/home/qhong7/github/SLUSCHI/src/mds_src/getfile.py", sys.argv[1]])

    # 2. script_v4.csh 0
    env = os.environ.copy()
    env["PATH"] = "/home/qhong7/github/asdf/src/" + ":" + env.get("PATH", "")
    run(["cp", "/home/qhong7/github/asdf/src/run_scripts/driver.sh", "."])
    run(["cp", "/home/qhong7/github/asdf/src/run_scripts/jobsub_master.sh", "."])

    # 2. script_v4.csh 0
    n = 51
    run(["/usr/bin/bash", "./driver.sh", str(n), "OUTCAR_collect"])

    run(["grep", "Entropy per mol", "all_jobs.out"])

    import math
    
    # how many jobs you ran (must be odd)
    if n % 2 == 0:
        raise ValueError("n must be an odd number")
    
    values = []
    
    with open("all_jobs.out") as f:
        for line in f:
            if "Total:" in line:
                parts = line.split()
                correction = float(parts[1])
            if "Entropy per mol" in line:
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        values.append(float(parts[5]))
                    except ValueError:
                        pass
    
    if len(values) < n:
        raise RuntimeError(f"Fewer than {n} valid entropy values found")
    
    values.sort()
    
    median_index = n // 2
    result = values[median_index] + correction
    print(f"{'Median entropy (raw):':32s} {values[median_index]:12.6g}  J/K/mol")
    print(f"{'Correction:':32s} {correction:12.6g}  J/K/mol")
    print(f"{'Final entropy:':32s} {result:12.6g}  J/K/mol")

    import matplotlib
    matplotlib.use("Agg")  # important for servers / no GUI
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Convert to numpy array
    vals = np.array(values)
    
    fig, ax = plt.subplots(figsize=(6, 4))
    
    # Violin plot
    vp = ax.violinplot(
        vals,
        showmeans=True,
        showmedians=True,
        showextrema=True
    )
    
    # Styling (portable, no fancy backend tricks)
    for body in vp['bodies']:
        body.set_alpha(0.6)
    
    ax.set_ylabel("Entropy (J/K/mol)")
    ax.set_title("Entropy Distribution")
    
    # Annotate median value
    median_val = np.median(vals)
    ax.axhline(median_val, linestyle="--", linewidth=1)
    ax.text(
        1.02,
        median_val,
        f"median = {median_val:.4g}",
        va="center",
        fontsize=9
    )
    
    # Remove x-axis ticks (single distribution)
    ax.set_xticks([])
    
    plt.tight_layout()
    plt.savefig("entropy_violin.png", dpi=300)
    plt.close()
    
    print("Saved entropy_violin.png")

if __name__ == "__main__":
    main()


