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

def run_cmd_streaming(cmd: List[str],
                      cwd: Optional[str] = None,
                      env: Optional[dict] = None,
                      stdout_capture_path: Optional[str] = None,
                      timeout: Optional[float] = None
                      ) -> Tuple[int, str, str]:
    """
    Run cmd (list) and stream stdout/stderr live to the console while also
    capturing them. Writes stdout to stdout_capture_path if provided.

    Returns: (returncode, full_stdout, full_stderr)
    """
    # Ensure we call the right interpreter if cmd is a python script invocation.
    # If callers already include sys.executable, this is a no-op.
    if isinstance(cmd, (list, tuple)) and len(cmd) >= 1 and cmd[0].endswith("python"):
        cmd[0] = sys.executable

    # Make a safe copy of env
    env_copy = os.environ.copy() if env is None else dict(env)

    # Ensure child does not inherit an open stdin (prevents interactive hangs)
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        stdin=subprocess.DEVNULL,
        cwd=cwd,
        env=env_copy,
        bufsize=1,  # line-buffered
        universal_newlines=True  # text mode
    )

    stdout_lines = []
    stderr_lines = []

    # Helper to read stderr in a thread (so we can stream both without deadlock)
    def _read_stderr():
        try:
            for line in proc.stderr:
                if line is None:
                    break
                stderr_lines.append(line)
                # print stderr lines to stderr immediately
                sys.stderr.write(line)
                sys.stderr.flush()
        except Exception:
            pass

    stderr_thread = threading.Thread(target=_read_stderr, daemon=True)
    stderr_thread.start()

    # Stream stdout line-by-line; write to capture file if requested
    if stdout_capture_path:
        try:
            fh = open(stdout_capture_path, "w", encoding="utf-8")
        except Exception:
            fh = None
    else:
        fh = None

    try:
        if proc.stdout is None:
            # No stdout to read; wait and collect stderr instead
            proc.wait(timeout=timeout)
        else:
            for line in proc.stdout:
                if line is None:
                    break
                stdout_lines.append(line)
                # print to stdout immediately
                sys.stdout.write(line)
                sys.stdout.flush()
                if fh:
                    fh.write(line)
            # ensure stdout iterator drained
            proc.stdout.close()
            # wait for process to finish
            proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()
        stderr_thread.join(timeout=1.0)
        if fh:
            fh.close()
        full_out = "".join(stdout_lines)
        full_err = "".join(stderr_lines)
        return proc.returncode, full_out, full_err
    finally:
        if fh:
            fh.close()

    # Ensure stderr thread has finished reading
    stderr_thread.join(timeout=1.0)

    full_out = "".join(stdout_lines)
    full_err = "".join(stderr_lines)
    return proc.returncode, full_out, full_err
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
