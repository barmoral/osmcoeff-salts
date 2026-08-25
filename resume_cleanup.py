#!/usr/bin/env python3
"""
resume_cleanup.py - make a chained SLURM job safe to resume.

    python resume_cleanup.py <case_dir> <expected_frames> [--report-only]

Why this exists
---------------
osmotic_sim_dispatch_membar_cont.py decides whether to redo a replicate with

    simulation_exists(): equil_sim_NVT_*, equil_sim_NPT_* and prod_sim_* exist?

A replicate killed by the wall clock has all three folders on disk, so that
check calls it finished and the next job in the chain skips it - leaving a
truncated trajectory silently in the dataset.

This script reads the frame count from each production DCD header and deletes
the folders of any replicate that did not reach `expected_frames`, so the next
job restarts it from the beginning.

Run it BEFORE the dispatch in every chained job. Safe on a fresh directory.

Only needs the standard library, so it works in either conda env.
"""
import argparse
import glob
import os
import shutil
import struct


def dcd_nframes(path):
    """Frame count from a DCD header.

    Layout: 4-byte record marker, 'CORD', then NSET (int32) at byte offset 8.
    OpenMM rewrites NSET on every frame, so a killed run reports what it wrote.
    Returns -1 if the file is unreadable.
    """
    try:
        with open(path, "rb") as fh:
            fh.seek(8)
            return struct.unpack("<i", fh.read(4))[0]
    except Exception:
        return -1


def main():
    p = argparse.ArgumentParser()
    p.add_argument("case")
    p.add_argument("expected", type=int)
    p.add_argument("--report-only", action="store_true",
                   help="list replicate status without deleting anything")
    a = p.parse_args()

    if not os.path.isdir(a.case):
        print(f"  {a.case}/ does not exist yet - fresh start, nothing to clean")
        return

    prods = sorted(glob.glob(os.path.join(a.case, "prod_sim_*")))
    if not prods:
        print("  no production folders yet - fresh start")
        return

    complete, removed = 0, 0
    for prod in prods:
        tag = os.path.basename(prod).replace("prod_sim_", "")
        dcds = glob.glob(os.path.join(prod, "*trajectory.dcd"))
        n = dcd_nframes(dcds[0]) if dcds else 0

        if n >= a.expected:
            print(f"  complete  {tag}: {n} frames")
            complete += 1
        elif a.report_only:
            print(f"  PARTIAL   {tag}: {n}/{a.expected} frames")
        else:
            print(f"  REDO      {tag}: {n}/{a.expected} frames - removing folders")
            for stage in ("equil_sim_NVT", "equil_sim_NPT", "prod_sim"):
                d = os.path.join(a.case, f"{stage}_{tag}")
                if os.path.isdir(d):
                    shutil.rmtree(d)
                    removed += 1

    if a.report_only:
        print(f"  {complete} replicate(s) complete at {a.expected} frames")
    else:
        print(f"  {complete} complete, removed {removed} folder(s) for restart")


if __name__ == "__main__":
    main()
