#!/usr/bin/env python3
"""Cross-check cadical-xor against the reference oracle (and CMS if available).

For every CNF-XOR file in the suite:
  * compute the ground-truth status with the blast-to-CNF oracle (plain cadical)
  * run cadical-xor (this repo's build/cadical, which understands 'x' lines)
  * optionally run CryptoMiniSat
and report any disagreement.  Exit code is non-zero if any mismatch is found.
"""
import argparse
import glob
import os
import subprocess
import sys

import oracle

REPO = oracle.REPO
CADICAL_XOR = os.path.join(REPO, "build", "cadical")
CMS = os.environ.get("CMS_BIN", "")


def run_solver(binary, path, extra=None):
    cmd = [binary, "-q"] + (extra or []) + [path]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return "TIMEOUT"
    except FileNotFoundError:
        return "NOBIN"
    if r.returncode == 10:
        return "SAT"
    if r.returncode == 20:
        return "UNSAT"
    return f"ERR(rc={r.returncode})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default=os.path.join(os.path.dirname(__file__), "suite"))
    ap.add_argument("--cadical-xor", default=CADICAL_XOR)
    ap.add_argument("--cms", default=CMS)
    ap.add_argument("--mode", choices=["gauss", "blast"], default="gauss",
                    help="cadical-xor XOR engine: Gauss-Jordan or CNF blasting")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    cx_flags = (["--gauss=1", "--xorblast=0"] if args.mode == "gauss"
                else ["--gauss=0", "--xorblast=1"])

    files = sorted(glob.glob(os.path.join(args.suite, "*.cnf")))
    mism = 0
    for p in files:
        _, _, _, decl = oracle.parse_cnfxor(p)
        ref = oracle.oracle_status(p)
        cx = run_solver(args.cadical_xor, p, cx_flags)
        row = f"{os.path.basename(p):28s} oracle={ref:6s} cadical-xor[{args.mode}]={cx:10s}"
        bad = cx not in ("SAT", "UNSAT") or cx != ref
        if args.cms:
            cms = run_solver(args.cms, p)
            row += f" cms={cms:6s}"
            bad = bad or (cms in ("SAT", "UNSAT") and cms != ref)
        if decl and decl != ref:
            row += f" !!DECL={decl}"
        if bad:
            row += "  <<< MISMATCH"
            mism += 1
        if args.verbose or bad:
            print(row)
    print(f"\n{len(files)} files, {mism} mismatches")
    sys.exit(1 if mism else 0)


if __name__ == "__main__":
    main()
