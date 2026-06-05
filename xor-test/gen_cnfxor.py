#!/usr/bin/env python3
"""Generate random CNF-XOR formulas in CryptoMiniSat-compatible DIMACS.

XOR semantics (CMS-compatible):
    A line "x l1 l2 ... lk 0" asserts  XOR(l1,...,lk) = true,
    where a negative literal -v contributes (1 XOR v); i.e. each negated
    literal flips the right-hand-side parity.  Equivalently, over the
    *variables* v1..vk:  v1 XOR ... XOR vk = (1 XOR #negatives mod 2).

The header "p cnf V C" counts V = #variables and C = #(regular clauses)
+ #(xor lines), so both CMS and cadical-xor agree on the clause count.

Each generated file is SAT or UNSAT *by construction* and records the
intended status in a leading comment line "c STATUS SAT|UNSAT".

Generation strategies
---------------------
sat        : plant a random full assignment; emit XOR + CNF constraints all
             satisfied by it.  Guaranteed SAT.
unsat_xor  : build a linearly *inconsistent* XOR system (pure Gauss-Jordan
             refutation) optionally padded with satisfiable CNF.  Guaranteed
             UNSAT and specifically exercises GJ.
unsat_mix  : plant an assignment, emit consistent XORs, then add CNF clauses
             that forbid that assignment's projection while XORs pin it down,
             producing CNF+XOR UNSAT.
"""
import argparse
import random
import sys


def parity_rhs_for_lits(lits, assign):
    """Return XOR over given literals under assignment (dict var->0/1)."""
    p = 0
    for l in lits:
        v = abs(l)
        val = assign[v]
        if l < 0:
            val ^= 1
        p ^= val
    return p


def lits_from_vars(vars_, negate_count_to_set_rhs, rhs, rng):
    """Build literals over `vars_` whose XOR(literals)==rhs.

    Start all-positive (literal XOR == XOR of var-values placeholder); we only
    control signs here, the caller guarantees var-value parity. We flip an
    even/odd number of signs to realize the requested literal-level rhs offset.
    """
    lits = list(vars_)
    # randomly negate some literals, then fix parity with one extra flip.
    neg = 0
    for i in range(len(lits)):
        if rng.random() < 0.5:
            lits[i] = -lits[i]
            neg += 1
    # current literal-rhs contributed by signs is `neg % 2` relative to vars.
    # caller handles the var-value parity; here ensure (#neg %2) == desired.
    if (neg % 2) != (negate_count_to_set_rhs % 2):
        lits[0] = -lits[0]
    return lits


def emit(out, header_vars, clauses, xors, status):
    n = header_vars
    c = len(clauses) + len(xors)
    out.write(f"c STATUS {status}\n")
    out.write(f"p cnf {n} {c}\n")
    for cl in clauses:
        out.write(" ".join(str(x) for x in cl) + " 0\n")
    for xr in xors:
        out.write("x " + " ".join(str(x) for x in xr) + " 0\n")


def gen_sat(nv, n_xor, n_cnf, xlen, clen, rng):
    assign = {v: rng.randint(0, 1) for v in range(1, nv + 1)}
    xors = []
    for _ in range(n_xor):
        k = max(2, min(xlen, nv))
        vars_ = rng.sample(range(1, nv + 1), k)
        lits = [v if rng.random() < 0.5 else -v for v in vars_]
        # force XOR(lits)=true by flipping one literal sign if needed
        if parity_rhs_for_lits(lits, assign) != 1:
            lits[0] = -lits[0]
        xors.append(lits)
    clauses = []
    for _ in range(n_cnf):
        k = max(1, min(clen, nv))
        vars_ = rng.sample(range(1, nv + 1), k)
        lits = [v if rng.random() < 0.5 else -v for v in vars_]
        # ensure at least one literal is satisfied by assignment
        if not any((assign[abs(l)] == 1) == (l > 0) for l in lits):
            i = rng.randrange(len(lits))
            v = abs(lits[i])
            lits[i] = v if assign[v] == 1 else -v
        clauses.append(lits)
    return clauses, xors, "SAT"


def gen_unsat_xor(nv, n_xor, n_cnf, xlen, clen, rng):
    # Build a consistent random linear system, then add one row equal to the
    # XOR-sum of a random subset but with flipped RHS -> contradiction.
    assign = {v: rng.randint(0, 1) for v in range(1, nv + 1)}
    xors = []
    rows = []  # store (lits)
    for _ in range(n_xor):
        k = max(2, min(xlen, nv))
        vars_ = rng.sample(range(1, nv + 1), k)
        lits = [v if rng.random() < 0.5 else -v for v in vars_]
        if parity_rhs_for_lits(lits, assign) != 1:
            lits[0] = -lits[0]
        xors.append(lits)
        rows.append(lits)
    # contradiction: pick a subset, XOR them (variable-wise), produce a row
    # over the symmetric-difference of variables with the OPPOSITE rhs.
    subset = rng.sample(rows, max(1, min(len(rows), rng.randint(1, 3))))
    var_count = {}
    base_rhs = 0
    for r in subset:
        # rhs of each row (as literals) is 1; accumulate
        base_rhs ^= 1
        for l in r:
            v = abs(l)
            sign_neg = 1 if l < 0 else 0
            var_count[v] = var_count.get(v, [0, 0])
            var_count[v][0] ^= 1          # appearance parity
            var_count[v][1] ^= sign_neg   # accumulated negation parity
    contradiction_lits = []
    neg_parity = 0
    for v, (appear, negp) in var_count.items():
        if appear:
            contradiction_lits.append(-v if negp else v)
            neg_parity ^= negp
    if not contradiction_lits:
        # subset cancelled to empty; force a trivial contradiction x1 != x1
        v = rng.randint(1, nv)
        xors.append([v])           # v = true
        xors.append([-v])          # v = false  -> unsat
        clauses = _rand_sat_cnf(nv, n_cnf, clen, assign, rng)
        return clauses, xors, "UNSAT"
    # The honest XOR-sum of the subset literals equals base_rhs(=len(subset)%2)
    # Flip it to create an inconsistent constraint.
    if parity_rhs_for_lits(contradiction_lits, assign) == 1:
        # currently consistent with assign -> flip a sign to make rhs wrong
        contradiction_lits[0] = -contradiction_lits[0]
    xors.append(contradiction_lits)
    clauses = _rand_sat_cnf(nv, n_cnf, clen, assign, rng)
    return clauses, xors, "UNSAT"


def _rand_sat_cnf(nv, n_cnf, clen, assign, rng):
    clauses = []
    for _ in range(n_cnf):
        k = max(1, min(clen, nv))
        vars_ = rng.sample(range(1, nv + 1), k)
        lits = [v if rng.random() < 0.5 else -v for v in vars_]
        if not any((assign[abs(l)] == 1) == (l > 0) for l in lits):
            i = rng.randrange(len(lits))
            v = abs(lits[i])
            lits[i] = v if assign[v] == 1 else -v
        clauses.append(lits)
    return clauses


def gen_unsat_mix(nv, n_xor, n_cnf, xlen, clen, rng):
    # Plant assignment, pin a small set of vars with unit XORs, then add a CNF
    # clause that is falsified by exactly that pinned assignment.
    assign = {v: rng.randint(0, 1) for v in range(1, nv + 1)}
    xors = []
    for _ in range(n_xor):
        k = max(2, min(xlen, nv))
        vars_ = rng.sample(range(1, nv + 1), k)
        lits = [v if rng.random() < 0.5 else -v for v in vars_]
        if parity_rhs_for_lits(lits, assign) != 1:
            lits[0] = -lits[0]
        xors.append(lits)
    pin = rng.sample(range(1, nv + 1), max(1, min(4, nv)))
    for v in pin:
        # unit xor pinning v to assign[v]: "x v 0" means v=true; "x -v 0" v=false
        xors.append([v] if assign[v] == 1 else [-v])
    # blocking clause over pinned vars: falsified exactly when all == assign
    block = [(-v if assign[v] == 1 else v) for v in pin]
    clauses = _rand_sat_cnf(nv, n_cnf, clen, assign, rng)
    clauses.append(block)
    return clauses, xors, "UNSAT"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["sat", "unsat_xor", "unsat_mix"], required=True)
    ap.add_argument("--vars", type=int, default=30)
    ap.add_argument("--xors", type=int, default=20)
    ap.add_argument("--cnf", type=int, default=40)
    ap.add_argument("--xlen", type=int, default=4)
    ap.add_argument("--clen", type=int, default=3)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("-o", "--out", default="-")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    if args.mode == "sat":
        clauses, xors, status = gen_sat(args.vars, args.xors, args.cnf, args.xlen, args.clen, rng)
    elif args.mode == "unsat_xor":
        clauses, xors, status = gen_unsat_xor(args.vars, args.xors, args.cnf, args.xlen, args.clen, rng)
    else:
        clauses, xors, status = gen_unsat_mix(args.vars, args.xors, args.cnf, args.xlen, args.clen, rng)

    out = sys.stdout if args.out == "-" else open(args.out, "w")
    emit(out, args.vars, clauses, xors, status)
    if out is not sys.stdout:
        out.close()


if __name__ == "__main__":
    main()
