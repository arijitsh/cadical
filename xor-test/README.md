# cadical-xor test & validation harness

Validates CaDiCaL's CNF-XOR support (CMS-compatible `x` DIMACS lines and the
`add_xor_clause` API) against independent ground truth.

## XOR format (CryptoMiniSat-compatible)
A line `x l1 l2 ... lk 0` asserts `XOR(l1,...,lk) == true`, where a negative
literal `-v` contributes `(1 XOR v)` (each negation flips the rhs parity).
The `p cnf V C` header counts XOR lines toward `C`.

## Files
- `gen_cnfxor.py` — generate SAT/UNSAT CNF-XOR formulas by construction
  (`--mode sat|unsat_xor|unsat_mix`). Each file records `c STATUS SAT|UNSAT`.
- `oracle.py` — independent ground truth: blasts every XOR to CNF and solves
  with plain `build/cadical`. Verdict is independent of the XOR engine.
- `run_tests.py` — runs cadical-xor (and CMS if `--cms`/`$CMS_BIN`) on the
  suite and diffs against the oracle. Non-zero exit on any mismatch.
- `check_model.py` — verifies a returned model satisfies the *original* CNF
  and XOR constraints (stronger than a SAT/UNSAT match).
- `test_matrix.cpp` — standalone unit test for the ported GF(2) matrix core
  (`src/gauss/packed_matrix.hpp`) vs an independent reference eliminator.
- `suite/` — generated formulas.

## XOR engines in cadical-xor
Two options control how XOR constraints are handled:
- `--gauss=1` enables the in-solver Gauss-Jordan engine (deep BCP integration;
  `src/gauss.cpp`, matrix core in `src/gauss/packed_matrix.hpp`).
- `--xorblast=1` (default) blasts each XOR to CNF (correct fallback, also used
  in proof mode where the GJ engine is currently disabled).

To exercise the GJ engine alone, run with `--gauss=1 --xorblast=0`.

## Reproduce
```sh
# (re)build cadical-xor
( cd .. && ./configure && make cadical )

# cross-check the Gauss-Jordan engine against the oracle
python3 run_tests.py --mode gauss          # verdicts
CADICAL_FLAGS="--gauss=1 --xorblast=0" python3 check_model.py suite/*.cnf
# and the blasting fallback
python3 run_tests.py --mode blast

# regenerate suite (optional; suite/ is committed)
python3 - <<'PY'
import subprocess
for mode in ("sat","unsat_xor","unsat_mix"):
    for seed in range(6):
        for nv,nx,nc in [(15,10,15),(30,20,30),(50,35,60)]:
            subprocess.run(["python3","gen_cnfxor.py","--mode",mode,"--vars",str(nv),
                "--xors",str(nx),"--cnf",str(nc),"--seed",str(seed),
                "-o",f"suite/{mode}_v{nv}_s{seed}.cnf"],check=True)
PY

python3 run_tests.py          # verdict cross-check vs oracle
python3 check_model.py suite/*.cnf   # model-level check

# GF(2) matrix core unit test
g++ -O2 -std=c++17 -I../src test_matrix.cpp -o test_matrix && ./test_matrix
```

## Performance benchmark (ApproxMC-style dense XORs)
`gen_dense.py` builds CNF + dense (~50%-of-variables) XOR-hash formulas, the
kind ApproxMC feeds to the solver.  Such XORs are too long to blast to CNF, so
only the Gauss-Jordan engine handles them.  `bench.py` times cadical-xor
(`--gauss=1 --xorblast=0`) against CryptoMiniSat and guards verdict agreement:
```sh
for n in 30 40 50 60; do for s in 0 1; do
  python3 gen_dense.py --vars $n --cnf $((2*n)) --xors $((n/2)) --density 0.5 \
    --random-rhs --seed $s -o perf/b_v${n}_s${s}.cnf
done; done
python3 bench.py --suite perf --timeout 30 \
  --cms /home/arijit/bins/trillium-bins/cryptominisat5
```
Status: results agree with CMS (0 disagreements across the suites), and the
engine performs *full* Gaussian elimination of the rows over the unassigned
columns every round (`Internal::gauss_round`), giving the full propagation
strength of Gaussian reasoning.  On small instances it is at parity with CMS
(constant factors either way); on the harder/larger dense instances it is
substantially **faster** than CMS (e.g. n=110-130 with ~n/2 dense XORs on a
phase-transition base: ~0.3-1.5s vs CMS ~5-18s).  The elimination is rebuilt
each round and only redone when a matrix variable's assignment changed; the
remaining headroom (vs CMS's incremental watched-column scheme) is to maintain
the reduced matrix incrementally instead of rebuilding it.

## CryptoMiniSat
CMS's bundled `cadiback` references CaDiCaL symbols (`get_eqiv_lits`,
`traverse_red_clauses`) absent from this fork, so the CMS binary does not link
against this repo's `libcadical`. The blast-oracle is used as ground truth
instead. To include CMS, build it against a compatible CaDiCaL and pass
`--cms /path/to/cryptominisat5` to `run_tests.py`.
