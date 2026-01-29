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

# paste/replace this function into master_bvs.py
import subprocess
import threading
import sys
from typing import List, Tuple, Optional
import os

# Replace/insert this implementation for run_cmd_streaming
import subprocess
import os

def _tail_of_file(path: str, max_bytes: int = 8192) -> str:
    """Return the last max_bytes of file as text (utf-8 with errors replaced)."""
    try:
        size = os.path.getsize(path)
        start = max(0, size - max_bytes)
        with open(path, "rb") as fh:
            fh.seek(start)
            data = fh.read()
        # decode with replacement to avoid decode errors
        return data.decode("utf-8", errors="replace")
    except Exception:
        return ""

def run_cmd_streaming(cmd,
                      cwd: str = None,
                      env: dict = None,
                      stdout_capture_path: str = None,
                      stderr_capture_path: str = None,
                      tail_bytes: int = 8192,
                      check: bool = False,
                      timeout: float = None):
    """
    Run command while streaming stdout/stderr to files to avoid pipe-blocking.

    Parameters
    ----------
    cmd : list or str
        Command to run (list is preferred).
    cwd : str
        Working directory.
    env : dict
        Environment variables (if None defaults to os.environ).
    stdout_capture_path : str
        Path to file where stdout will be written. If None, uses ./cmd_stdout.log.
    stderr_capture_path : str
        Path to file where stderr will be written. If None, uses ./cmd_stderr.log.
    tail_bytes : int
        How many trailing bytes to read back and return for stdout/stderr.
    check : bool
        If True, raise CalledProcessError on non-zero exit (like subprocess.run).
    timeout : float
        Optional timeout in seconds to kill the process.

    Returns
    -------
    (rc, stdout_tail, stderr_tail)
        rc : int return code
        stdout_tail, stderr_tail : str (decoded tail of the file)
    """
    if env is None:
        env = os.environ.copy()
    else:
        # make a shallow copy so we don't mutate caller's env
        tmp_env = os.environ.copy()
        tmp_env.update(env)
        env = tmp_env

    if stdout_capture_path is None:
        stdout_capture_path = os.path.join(cwd or os.getcwd(), "cmd_stdout.log")
    if stderr_capture_path is None:
        stderr_capture_path = os.path.join(cwd or os.getcwd(), "cmd_stderr.log")

    # Ensure parent directory exists
    os.makedirs(os.path.dirname(os.path.abspath(stdout_capture_path)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(stderr_capture_path)), exist_ok=True)

    # Open files in append mode so previous content is preserved (or use 'wb' to overwrite)
    # Use binary mode to avoid issues with encoding while writing raw bytes.
    with open(stdout_capture_path, "ab") as fh_out, open(stderr_capture_path, "ab") as fh_err:
        try:
            # If cmd is a string, let shell=False call it directly; prefer list input when possible.
            # Here we use subprocess.run with stdout/stderr redirected to files -> no pipes to block.
            completed = subprocess.run(
                cmd,
                cwd=cwd,
                env=env,
                stdout=fh_out,
                stderr=fh_err,
                check=False,
                timeout=timeout,
                shell=isinstance(cmd, str)
            )
            rc = completed.returncode
        except subprocess.TimeoutExpired as te:
            rc = -9
            # write a timeout message into stderr file
            fh_err.write(f"\n*** TIMEOUT after {timeout} seconds ***\n".encode("utf-8", errors="replace"))
            fh_err.flush()
        except Exception as e:
            rc = -1
            fh_err.write(f"\n*** EXCEPTION running command: {e} ***\n".encode("utf-8", errors="replace"))
            fh_err.flush()

    # Return small tails so caller can log or inspect messages without huge memory usage
    stdout_tail = _tail_of_file(stdout_capture_path, max_bytes=tail_bytes)
    stderr_tail = _tail_of_file(stderr_capture_path, max_bytes=tail_bytes)

    if check and rc != 0:
        # mimic subprocess.run(check=True) behavior
        raise subprocess.CalledProcessError(rc, cmd, output=stdout_tail, stderr=stderr_tail)

    return rc, stdout_tail, stderr_tail

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
