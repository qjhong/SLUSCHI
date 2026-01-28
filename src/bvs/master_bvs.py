#!/usr/bin/env python3
"""
master_bvs.py — simple runner for bvs module (API-friendly)

Assumptions:
 - Called from inside the job directory (app.py creates the job dir and copies the CIF there).
 - Argument is the CIF filename (basename or path relative to cwd).
 - bvs.py and plot_bvs.py are located at ../bvs/ or an absolute path. Adjust BVS_PY / PLOT_PY as needed.

Behavior:
 - Runs bvs.py; writes its stdout to bvs.txt (so callers can view it).
 - If bvs.py succeeds, runs plot_bvs.py.
 - Streams stdout/stderr to the process stdout/stderr so the API can capture partial logs.
 - Exits with non-zero on failure.
"""

from __future__ import annotations
import sys
import os
import subprocess
import traceback

# adjust these paths if your layout differs; they can be absolute or relative
BVS_PY = os.path.join(os.path.dirname(__file__), "bvs.py")          # src/bvs/bvs.py
PLOT_BVS_PY = os.path.join(os.path.dirname(__file__), "plot_bvs.py")  # src/bvs/plot_bvs.py

def run_cmd_streaming(cmd, cwd=None, env=None, stdout_capture_path=None):
    """
    Run command list `cmd`, stream stdout/stderr to the current process' stdout/stderr,
    and also capture stdout lines (returned). If stdout_capture_path is given, write the
    captured stdout to that path at the end.
    Returns (returncode, stdout_lines, stderr_lines)
    """
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=cwd, env=env)
    stdout_lines = []
    stderr_lines = []
    try:
        # read lines until process exits
        while True:
            out = proc.stdout.readline()
            err = proc.stderr.readline()
            if out:
                # stream to parent stdout and capture
                sys.stdout.write(out)
                sys.stdout.flush()
                stdout_lines.append(out)
            if err:
                sys.stderr.write(err)
                sys.stderr.flush()
                stderr_lines.append(err)
            if proc.poll() is not None and not out and not err:
                break
        rc = proc.wait()
        # write captured stdout to file if requested
        if stdout_capture_path:
            try:
                with open(stdout_capture_path, "w") as f:
                    f.write("".join(stdout_lines))
            except Exception as e:
                sys.stderr.write(f"Warning: failed to write {stdout_capture_path}: {e}\n")
        return rc, stdout_lines, stderr_lines
    except Exception:
        proc.kill()
        raise

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: master_bvs.py <str.cif>\n")
        return 2

    cif_arg = argv[1]
    # if user passed an absolute or relative path, consider basename
    cif_basename = os.path.basename(cif_arg)
    cif_path = cif_basename if os.path.exists(cif_basename) else cif_arg

    if not os.path.exists(cif_path):
        sys.stderr.write(f"Error: CIF file not found in cwd: {cif_path}\n")
        return 3

    # environment: inherit current env; caller (app.py) should set HOME/PATH if needed
    env = os.environ.copy()
    # write bvs stdout to bvs.txt in current dir
    bvs_txt = os.path.join(os.getcwd(), "bvs.txt")

    try:
        # 1) run bvs.py
        cmd1 = [sys.executable, BVS_PY, cif_path]
        sys.stdout.write(f"Running: {' '.join(cmd1)}\n")
        sys.stdout.flush()
        rc1, out1, err1 = run_cmd_streaming(cmd1, cwd=os.getcwd(), env=env, stdout_capture_path=bvs_txt)
        if rc1 != 0:
            sys.stderr.write(f"bvs.py failed (rc={rc1}); see stdout/stderr above\n")
            return rc1

        # 2) run plot_bvs.py (assumes it reads outputs in cwd)
        cmd2 = [sys.executable, PLOT_BVS_PY]
        sys.stdout.write(f"Running: {' '.join(cmd2)}\n")
        sys.stdout.flush()
        rc2, out2, err2 = run_cmd_streaming(cmd2, cwd=os.getcwd(), env=env)
        if rc2 != 0:
            sys.stderr.write(f"plot_bvs.py failed (rc={rc2}); see stdout/stderr above\n")
            return rc2

        # success
        sys.stdout.write("master_bvs: completed successfully\n")
        sys.stdout.flush()
        return 0

    except Exception:
        traceback.print_exc(file=sys.stderr)
        return 10

if __name__ == "__main__":
    rc = main(sys.argv)
    sys.exit(rc)
