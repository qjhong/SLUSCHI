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
    result = values[median_index]
    
    print('The median value of entropy is ' + str(result) + 'J/K/mol.')

if __name__ == "__main__":
    main()


