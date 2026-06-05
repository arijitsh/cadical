#!/usr/bin/env python3
"""Generate ApproxMC-style CNF + dense-XOR formulas.

Mimics the instances ApproxMC feeds to the SAT solver: a base CNF plus 'm'
random parity (XOR) hash constraints, where each variable is included in each
XOR independently with probability 'density' (0.5 by default => ~n/2 vars per
XOR).  Such XORs are far too long to blast to CNF (2^(n/2) clauses), so they
exercise the Gauss-Jordan engine specifically.

By default the formula is SAT by construction: a random solution is planted,
the base CNF is satisfied by it, and every XOR's right-hand side is set to the
parity the plant induces.  Use --random-rhs for random (mostly UNSAT once
m ' #vars) right-hand sides.
"""
import argparse
import random
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vars", type=int, default=60)
    ap.add_argument("--cnf", type=int, default=0,
                    help="number of random base 3-CNF clauses")
    ap.add_argument("--xors", type=int, required=True,
                    help="number of dense XOR hash constraints")
    ap.add_argument("--density", type=float, default=0.5)
    ap.add_argument("--random-rhs", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("-o", "--out", default="-")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    n = args.vars
    sol = {v: rng.randint(0, 1) for v in range(1, n + 1)}

    clauses = []
    for _ in range(args.cnf):
        vs = rng.sample(range(1, n + 1), min(3, n))
        lits = [v if rng.random() < 0.5 else -v for v in vs]
        if not any((sol[abs(l)] == 1) == (l > 0) for l in lits):
            i = rng.randrange(len(lits))
            v = abs(lits[i])
            lits[i] = v if sol[v] == 1 else -v
        clauses.append(lits)

    xors = []
    status = "SAT"
    for _ in range(args.xors):
        vs = [v for v in range(1, n + 1) if rng.random() < args.density]
        if not vs:
            vs = [rng.randrange(1, n + 1)]
        if args.random_rhs:
            rhs = rng.randint(0, 1)
            status = "UNKNOWN"  # likely UNSAT for large m, but not guaranteed
            lits = list(vs)
            # encode rhs by flipping one literal's sign: XOR(literals)=true means
            # parity of literals is 1; to get var-parity == rhs we add a sign.
            parity_vars = 0  # we want XOR(vars)=rhs ; XOR(literals)=1 always
            # number of negative literals needed: (1 ^ rhs)
            if (1 ^ rhs) & 1:
                lits[0] = -lits[0]
        else:
            # SAT: set rhs to the plant's parity so XOR(literals)==true holds
            lits = list(vs)
            p = 0
            for v in vs:
                p ^= sol[v]
            # XOR(literals) must be 1; currently all-positive literal parity = p
            if p != 1:
                lits[0] = -lits[0]
        xors.append(lits)

    out = sys.stdout if args.out == "-" else open(args.out, "w")
    out.write(f"c STATUS {status}\n")
    out.write(f"c dense XORs: {args.xors} over {n} vars, density {args.density}\n")
    out.write(f"p cnf {n} {len(clauses) + len(xors)}\n")
    for c in clauses:
        out.write(" ".join(map(str, c)) + " 0\n")
    for x in xors:
        out.write("x " + " ".join(map(str, x)) + " 0\n")
    if out is not sys.stdout:
        out.close()


if __name__ == "__main__":
    main()
