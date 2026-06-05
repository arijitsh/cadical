#!/usr/bin/env python3
"""Verify that cadical-xor's reported model satisfies the ORIGINAL CNF-XOR.

This is stronger than comparing SAT/UNSAT verdicts: it confirms the returned
assignment actually satisfies every ordinary clause AND every XOR constraint
(XOR(literals) == true) of the input file.
"""
import os
import subprocess
import sys

import oracle


EXTRA = os.environ.get("CADICAL_FLAGS", "").split()


def get_model(path, solver=oracle.DEFAULT_SOLVER):
    r = subprocess.run([solver] + EXTRA + [path], capture_output=True, text=True)
    if r.returncode == 20:
        return None, "UNSAT"
    if r.returncode != 10:
        return None, f"ERR(rc={r.returncode})"
    val = {}
    for line in r.stdout.splitlines():
        if line.startswith("v "):
            for tok in line[2:].split():
                lit = int(tok)
                if lit == 0:
                    continue
                val[abs(lit)] = 1 if lit > 0 else 0
    return val, "SAT"


def check(path):
    nvars, clauses, xors, decl = oracle.parse_cnfxor(path)
    val, status = get_model(path)
    if status != "SAT":
        return status, None  # nothing to check for UNSAT/ERR
    # ordinary clauses
    for c in clauses:
        if not any(val.get(abs(l), 0) == (1 if l > 0 else 0) for l in c):
            return "SAT", f"clause violated: {c}"
    # xor constraints: XOR(literals)==true
    for lits, rhs in xors:
        p = 0
        for l in lits:
            v = val.get(abs(l), 0)
            if l < 0:
                v ^= 1
            p ^= v
        if p != rhs:
            return "SAT", f"xor violated: {lits} parity={p} expected={rhs}"
    return "SAT", None


def main():
    files = sys.argv[1:]
    bad = 0
    sat = 0
    for p in files:
        st, err = check(p)
        if st == "SAT":
            sat += 1
            if err:
                bad += 1
                print(f"{os.path.basename(p)}: BAD MODEL -> {err}")
    print(f"\nchecked {len(files)} files, {sat} SAT, {bad} bad models")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
