#!/usr/bin/env python3
"""Benchmark cadical-xor (Gauss-Jordan) against CryptoMiniSat on dense-XOR
(ApproxMC-style) formulas.  Reports per-instance wall-clock solve time and
totals, and flags any verdict disagreement (correctness guard)."""
import argparse
import glob
import os
import subprocess
import time

import oracle

CAD = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "build", "cadical"))


def run(cmd, timeout):
    t0 = time.perf_counter()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        dt = time.perf_counter() - t0
        v = {10: "SAT", 20: "UNSAT"}.get(r.returncode, f"E{r.returncode}")
        return v, dt
    except subprocess.TimeoutExpired:
        return "TIMEOUT", timeout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default=os.path.join(os.path.dirname(__file__), "perf"))
    ap.add_argument("--cms", default="/home/arijit/bins/trillium-bins/cryptominisat5")
    ap.add_argument("--timeout", type=float, default=60.0)
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.suite, "*.cnf")))
    tot_g = tot_c = 0.0
    print(f"{'instance':26s} {'verdict':7s} {'gauss(s)':>10s} {'cms(s)':>10s} {'speedup':>8s}")
    print("-" * 66)
    disagree = 0
    for p in files:
        gv, gt = run([CAD, "-q", "--gauss=1", "--xorblast=0", p], args.timeout)
        cv, ct = run([args.cms, "--verb=0", p], args.timeout)
        ref = None
        flag = ""
        if gv in ("SAT", "UNSAT") and cv in ("SAT", "UNSAT") and gv != cv:
            flag = "  <<< DISAGREE"
            disagree += 1
        sp = (ct / gt) if (gt > 0 and gv not in ("TIMEOUT",) and cv not in ("TIMEOUT",)) else 0
        if gv != "TIMEOUT":
            tot_g += gt
        if cv != "TIMEOUT":
            tot_c += ct
        print(f"{os.path.basename(p):26s} {gv:7s} {gt:10.3f} {ct:10.3f} {sp:7.2f}x{flag}")
    print("-" * 66)
    print(f"{'TOTAL':26s} {'':7s} {tot_g:10.3f} {tot_c:10.3f} {(tot_c/tot_g if tot_g else 0):7.2f}x")
    print(f"\n{len(files)} instances, {disagree} verdict disagreements")


if __name__ == "__main__":
    main()
