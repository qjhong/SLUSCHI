#!/usr/bin/env python3
"""
master_bvs.py

Usage:
    python3 master_bvs.py path/to/structure.cif

This script runs:
    python bvs.py <cif>
then (if successful)
    python plot_bvs.py

It creates a job directory under /tmp/sluschi_jobs/<job_id>/ and writes:
 - stdout.log, stderr.log
 - summary.json
 - any artifacts produced (CSV, PNG, txt, etc.)

Finally it atomically updates /tmp/sluschi_last_result.json with the summary.
"""

import os
import sys
import json
import time
import uuid
import shutil
import subprocess
import glob
import traceback

JOB_BASE = "/tmp/sluschi_jobs"
LAST_RESULT_PATH = "/tmp/sluschi_last_result.json"

# Names of files the bvs/plot scripts commonly produce (adjust as needed)
COMMON_ARTIFACT_PATTERNS = [
    "bvs_results.csv",
    "bvs.txt",
    "bvs_summary.txt",
    "bvs_plot.png",
    "plot_bvs.png",
    "*.png",
    "*.csv",
    "*.txt",
]


def atomic_write(path, data):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2)
    os.replace(tmp, path)


def run_cmd_streaming(cmd, cwd=None, env=None, stdout_path=None, stderr_path=None):
    """
    Run cmd (list) streaming stdout/stderr lines into buffers and writing
    partial logs to stdout_path/stderr_path periodically. Returns (retcode, stdout_buf, stderr_buf).
    """
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=cwd, env=env)
    stdout_buf = []
    stderr_buf = []
    try:
        while True:
            out = proc.stdout.readline()
            err = proc.stderr.readline()
            if out:
                stdout_buf.append(out)
            if err:
                stderr_buf.append(err)

            # write partial logs (best-effort)
            if stdout_path:
                try:
                    with open(stdout_path, "w") as f:
                        f.write("".join(stdout_buf))
                except Exception:
                    pass
            if stderr_path:
                try:
                    with open(stderr_path, "w") as f:
                        f.write("".join(stderr_buf))
                except Exception:
                    pass

            if proc.poll() is not None and not out and not err:
                break
        ret = proc.wait()
        return ret, stdout_buf, stderr_buf
    except Exception:
        proc.kill()
        raise


def collect_artifacts(job_dir, src_dirs):
    """
    Collects files matching COMMON_ARTIFACT_PATTERNS from src_dirs and job_dir
    into job_dir (no-op for files already in job_dir). Returns list of basenames.
    """
    artifacts = []
    for d in src_dirs:
        if not d:
            continue
        for pat in COMMON_ARTIFACT_PATTERNS:
            # glob relative to d
            try:
                matches = glob.glob(os.path.join(d, pat))
            except Exception:
                matches = []
            for m in matches:
                if os.path.isfile(m):
                    try:
                        dst = os.path.join(job_dir, os.path.basename(m))
                        # copy only if source not same as dest or dest doesn't exist
                        if os.path.abspath(m) != os.path.abspath(dst):
                            shutil.copy2(m, dst)
                        artifacts.append(os.path.basename(dst))
                    except Exception:
                        pass

    # also include any files already in the job directory (except summary.json)
    for f in os.listdir(job_dir):
        if f not in artifacts and f != "summary.json":
            artifacts.append(f)

    # de-dupe while preserving order
    seen = []
    out = []
    for a in artifacts:
        if a not in seen:
            seen.append(a)
            out.append(a)
    return out


def main(argv):
    if len(argv) < 2:
        print("Usage: python3 master_bvs.py path/to/structure.cif", file=sys.stderr)
        return 2

    cif_path = argv[1]
    if not os.path.exists(cif_path):
        print("Error: cif file not found:", cif_path, file=sys.stderr)
        return 3

    os.makedirs(JOB_BASE, exist_ok=True)
    job_id = str(uuid.uuid4())
    job_dir = os.path.join(JOB_BASE, job_id)
    os.makedirs(job_dir, exist_ok=True)

    started_at = time.time()
    summary = {
        "job_id": job_id,
        "input": os.path.abspath(cif_path),
        "status": "running",
        "started_at": started_at,
        "finished_at": None,
        "duration_s": None,
        "returncode": None,
        "artifacts": [],
        "stdout": "stdout.log",
        "stderr": "stderr.log",
    }

    # copy the CIF into job dir so bvs runs locally
    try:
        cif_basename = os.path.basename(cif_path)
        cif_dest = os.path.join(job_dir, cif_basename)
        shutil.copy2(cif_path, cif_dest)
    except Exception:
        pass

    stdout_path = os.path.join(job_dir, "stdout.log")
    stderr_path = os.path.join(job_dir, "stderr.log")

    # write initial summary (shared)
    try:
        with open(os.path.join(job_dir, "summary.json"), "w") as f:
            json.dump(summary, f, indent=2)
        tmp_shared = LAST_RESULT_PATH + ".tmp"
        shutil.copy2(os.path.join(job_dir, "summary.json"), tmp_shared)
        os.replace(tmp_shared, LAST_RESULT_PATH)
    except Exception:
        pass

    env = os.environ.copy()
    env["HOME"] = env.get("HOME", "/root")

    retcode = 1
    stdout_buf_total = []
    stderr_buf_total = []
    error_trace = None

    try:
        # 1) run bvs.py
        # --- run bvs.py and tee stdout to bvs.txt ---
        bvs_txt_path = os.path.join(job_dir, "bvs.txt")
        
        cmd1 = [sys.executable, "/home/qhong7/github/SLUSCHI/src/bvs/bvs.py", cif_basename]
        
        rc1, out1, err1 = run_cmd_streaming(
            cmd1,
            cwd=job_dir,
            env=env,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
        )
        
        # write bvs stdout explicitly to bvs.txt
        try:
            with open(bvs_txt_path, "w") as f:
                f.write("".join(out1))
        except Exception:
            pass
        
        stdout_buf_total.extend(out1)
        stderr_buf_total.extend(err1)
        
        if rc1 != 0:
            retcode = rc1
            raise RuntimeError(f"bvs.py failed with return code {rc1}")

        # 2) run plot_bvs.py (assumes it reads bvs outputs in cwd)
        cmd2 = [sys.executable, "/home/qhong7/github/SLUSCHI/src/bvs/plot_bvs.py"]
        rc2, out2, err2 = run_cmd_streaming(cmd2, cwd=job_dir, env=env, stdout_path=stdout_path, stderr_path=stderr_path)
        stdout_buf_total.extend(out2)
        stderr_buf_total.extend(err2)

        if rc2 != 0:
            retcode = rc2
            raise RuntimeError(f"plot_bvs.py failed with return code {rc2}")

        retcode = 0

    except Exception as exc:
        error_trace = traceback.format_exc()
        stderr_buf_total.append("\n*** MASTER_BVS EXCEPTION ***\n")
        stderr_buf_total.append(error_trace)
        summary["status"] = "failed"

    finally:
        # final write logs
        try:
            with open(stdout_path, "w") as f:
                f.write("".join(stdout_buf_total))
        except Exception:
            pass
        try:
            with open(stderr_path, "w") as f:
                f.write("".join(stderr_buf_total))
        except Exception:
            pass

        finished_at = time.time()
        summary["finished_at"] = finished_at
        summary["duration_s"] = round(finished_at - started_at, 3)
        summary["returncode"] = retcode
        summary["status"] = "ok" if retcode == 0 else "failed"

        # collect artifacts from job dir and also current working dir (in case scripts write elsewhere)
        src_dirs = [job_dir, os.getcwd()]
        artifacts = collect_artifacts(job_dir, src_dirs)
        summary["artifacts"] = artifacts

        # write final summary and atomically update shared status
        try:
            with open(os.path.join(job_dir, "summary.json"), "w") as f:
                json.dump(summary, f, indent=2)
            tmp_shared = LAST_RESULT_PATH + ".tmp"
            shutil.copy2(os.path.join(job_dir, "summary.json"), tmp_shared)
            os.replace(tmp_shared, LAST_RESULT_PATH)
        except Exception:
            pass

    # print final summary JSON to stdout so caller can capture it
    print(json.dumps(summary, indent=2))
    return 0 if retcode == 0 else 1


if __name__ == "__main__":
    rc = main(sys.argv)
    sys.exit(rc)
