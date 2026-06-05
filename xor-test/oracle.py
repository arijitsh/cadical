#!/usr/bin/env python3
"""Reference oracle for CNF-XOR formulas.

Parses a CMS-style CNF-XOR file, blasts every XOR constraint into ordinary
CNF (introducing Tseitin aux vars for long XORs), and solves the resulting
pure-CNF formula with a plain SAT solver (default: the repo's build/cadical).

The oracle's SAT/UNSAT verdict is independent of cadical-xor's XOR engine and
of CryptoMiniSat, so it is a trustworthy ground truth for cross-checking.
"""
import itertools
import subprocess
import sys
import tempfile
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_SOLVER = os.path.join(REPO, "build", "cadical")


def parse_cnfxor(path):
    nvars = 0
    clauses = []   # list of list[int]
    xors = []      # list of (list[int] literals, rhs) ; rhs always 1 here
    declared_status = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("c"):
                parts = line.split()
                if len(parts) >= 3 and parts[1] == "STATUS":
                    declared_status = parts[2]
                continue
            if line.startswith("p"):
                parts = line.split()
                nvars = int(parts[2])
                continue
            if line.startswith("x"):
                toks = line[1:].split()
                lits = [int(t) for t in toks]
                assert lits[-1] == 0, f"xor line not 0-terminated: {line}"
                lits = lits[:-1]
                # CMS semantics: XOR(literals) = true.
                xors.append((lits, 1))
                continue
            toks = line.split()
            lits = [int(t) for t in toks]
            assert lits[-1] == 0, f"clause not 0-terminated: {line}"
            clauses.append(lits[:-1])
    return nvars, clauses, xors, declared_status


def blast_xor(lits, rhs, fresh):
    """Return CNF clauses encoding XOR(lits) == rhs.

    `fresh` is a callable returning a new aux variable id.
    Long XORs are split with Tseitin aux vars (cut length 6).
    """
    CUT = 6
    out = []
    # Reduce long xor into chained sub-xors of length <= CUT sharing aux vars.
    work = list(lits)
    cur_rhs = rhs
    while len(work) > CUT:
        head = work[:CUT - 1]
        a = fresh()
        # head ... XOR == a  (define a as XOR of head)
        out += _blast_small(head + [a], 0)   # XOR(head, a) == 0  => a == XOR(head)
        work = [a] + work[CUT - 1:]
    out += _blast_small(work, cur_rhs)
    return out


def _blast_small(lits, rhs):
    """Direct CNF for XOR(lits) == rhs with len(lits) small."""
    k = len(lits)
    out = []
    if k == 0:
        # empty XOR == rhs: rhs==0 -> trivially true; rhs==1 -> empty clause
        if rhs == 1:
            out.append([])  # empty clause -> UNSAT
        return out
    # forbid every assignment whose parity != rhs
    for bits in itertools.product([0, 1], repeat=k):
        parity = 0
        for b, l in zip(bits, lits):
            # literal value: if l>0 then var=b means literal true when b==1
            lit_true = b if l > 0 else (1 - b)
            parity ^= lit_true
        if parity != rhs:
            # clause that rules out this assignment: for each var, the literal
            # that is FALSE under `bits`
            clause = []
            for b, l in zip(bits, lits):
                v = abs(l)
                # under bits, var v assigned b; the satisfying literal to add is
                # the one that is true when v != b, i.e. negation of current.
                clause.append(-v if b == 1 else v)
            out.append(clause)
    return out


def to_cnf(nvars, clauses, xors):
    counter = [nvars]

    def fresh():
        counter[0] += 1
        return counter[0]

    out = list(clauses)
    for lits, rhs in xors:
        out += blast_xor(lits, rhs, fresh)
    return counter[0], out


def write_dimacs(path, nvars, clauses):
    with open(path, "w") as f:
        f.write(f"p cnf {nvars} {len(clauses)}\n")
        for c in clauses:
            f.write(" ".join(str(x) for x in c) + " 0\n")


def solve(path, solver=DEFAULT_SOLVER):
    r = subprocess.run([solver, "-q", path], capture_output=True, text=True)
    if r.returncode == 10:
        return "SAT"
    if r.returncode == 20:
        return "UNSAT"
    return f"ERR(rc={r.returncode})"


def oracle_status(cnfxor_path, solver=DEFAULT_SOLVER):
    nvars, clauses, xors, _ = parse_cnfxor(cnfxor_path)
    nv2, cnf = to_cnf(nvars, clauses, xors)
    if any(len(c) == 0 for c in cnf):
        return "UNSAT"
    fd, tmp = tempfile.mkstemp(suffix=".cnf")
    os.close(fd)
    try:
        write_dimacs(tmp, nv2, cnf)
        return solve(tmp, solver)
    finally:
        os.unlink(tmp)


if __name__ == "__main__":
    for p in sys.argv[1:]:
        _, _, _, decl = parse_cnfxor(p)
        st = oracle_status(p)
        ok = "OK" if (decl is None or decl == st) else "MISMATCH"
        print(f"{p}: oracle={st} declared={decl} {ok}")
