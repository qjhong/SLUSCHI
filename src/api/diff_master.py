#!/usr/bin/env python3

import subprocess
import sys

def run(cmd, shell=False):
    print(f"Running: {cmd}")
    subprocess.check_call(cmd, shell=shell)

def main():
    # 1. python getfile.py
    run([sys.executable, "/home/qhong7/github/SLUSCHI/src/mds_src/getfile.py", sys.argv[1]])

    # 2. script_v4.csh 0
    run(["/home/qhong7/github/SLUSCHI/src/mds_src/script_v4.csh", "0"])

    # 3. python diffusion.py
    run([sys.executable, "/home/qhong7/github/SLUSCHI/src/diffusion.py"])

if __name__ == "__main__":
    main()

