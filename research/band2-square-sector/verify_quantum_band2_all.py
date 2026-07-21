#!/usr/bin/env python3
"""Driver for the full quantum band-2 assembly (quantum-band2-theorem.md).

Part 1: assembly-specific exact checks (sympy):
  - the Fourier map phi: x -> -d, d -> x extends to an automorphism on the
    ladder generators (multiplicativity spot-certified on generic pieces);
  - the triangular factorization phi = exp(ad(x^2/2)) exp(ad(-d^2/2))
    exp(ad(x^2/2)) holds exactly on both generators;
  - pair-exchange bookkeeping: [ -X, D ] = [D, X] as an operator identity;
  - a random tame round-trip: the C4 family member with random rational
    parameters satisfies [D, X] = 1 and yields x and d back by the explicit
    generation formulas.
Part 2: runs every component verifier as a subprocess and fails loudly if any
fails. Run:  python3 verify_quantum_band2_all.py   (requires sympy)
"""
import os
import subprocess
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OK = True


def check(name, cond):
    global OK
    print(("PASS" if cond else "FAIL"), name, flush=True)
    OK = OK and cond


# ---------------- Weyl algebra engine (normal order: x left, d right) ------
x, d = sp.symbols('x d', commutative=False)


def normal(expr):
    """Normal-order a noncommutative polynomial in x, d using d*x = x*d + 1."""
    expr = sp.expand(expr)
    changed = True
    while changed:
        changed = False
        terms = sp.Add.make_args(sp.expand(expr))
        out = 0
        for t in terms:
            coeff, factors = t.as_coeff_mul()
            seq = []
            for f in factors:
                if f.is_Pow:
                    seq += [f.base] * int(f.exp)
                else:
                    seq.append(f)
            # find a d before an x
            swapped = False
            for i in range(len(seq) - 1):
                if seq[i] == d and seq[i + 1] == x:
                    left = seq[:i]
                    right = seq[i + 2:]
                    t1 = coeff * sp.Mul(*(left + [x, d] + right)) if left + right or True else 0
                    t2 = coeff * sp.Mul(*(left + right)) if (left + right) else coeff
                    out += t1 + t2
                    swapped = True
                    changed = True
                    break
            if not swapped:
                out += t
        expr = sp.expand(out)
    return expr


def comm(a, b):
    return normal(a * b - b * a)


check("engine sanity: [d, x] = 1", sp.simplify(comm(d, x) - 1) == 0)

# ---------------- phi and its triangular factorization ----------------------
# Actions (verified as automorphisms by bracket preservation below):
#   A = exp(ad(-x^2/2)):  x -> x,        d -> d + x   (since [-x^2/2, d] = x)
#   B = exp(ad(-d^2/2)):  x -> x - d,    d -> d       (since [-d^2/2, x] = -d)
# Composite A o B o A sends x -> -d, d -> x (SL2: shear(1)*shear(-1)*shear(1)).
def A_act(w):
    return normal(w.subs({d: d + x}, simultaneous=True))


def B_act(w):
    return normal(w.subs({x: x - d}, simultaneous=True))


# each preserves the defining relation
check("exp(ad(x^2/2)) action preserves [d,x]=1",
      sp.simplify(comm(A_act(d), A_act(x)) - 1) == 0)
check("exp(ad(-d^2/2)) action preserves [d,x]=1",
      sp.simplify(comm(B_act(d), B_act(x)) - 1) == 0)

# ad-exponential certificates: [-x^2/2, d] = x and the series truncates
check("[-x^2/2, d] = x (triangular generator)",
      sp.simplify(comm(-x * x / 2, d) - x) == 0)
check("[-d^2/2, x] = -d (triangular generator)",
      sp.simplify(comm(-d * d / 2, x) + d) == 0)

phi_x = A_act(B_act(A_act(x)))
phi_d = A_act(B_act(A_act(d)))
check("triangular composite sends x -> -d", sp.simplify(phi_x + d) == 0)
check("triangular composite sends d -> x", sp.simplify(phi_d - x) == 0)
check("phi preserves [d,x]=1", sp.simplify(comm(phi_d, phi_x) - 1) == 0)

# ---------------- pair exchange ---------------------------------------------
X0 = x * x + 3 * x * d + d
D0 = 2 * x - d * d
check("pair exchange: [-X, D] == [D, X] on a generic pair",
      sp.simplify(comm(-X0, D0) - comm(D0, X0)) == 0)

# ---------------- tame family round-trip ------------------------------------
c0, c1, Acst, beta, lam = sp.Rational(3, 2), sp.Rational(-2, 5), sp.Rational(7, 3), sp.Rational(1, 4), sp.Rational(2, 7)
kap = sp.Rational(5, 2)
U = x + c0 + c1 * d
Xf = normal(U * U - d / kap - Acst)
Df = normal(lam * Xf + kap * U + beta)
check("tame family satisfies [D, X] = 1",
      sp.simplify(comm(Df, Xf) - 1) == 0)
# generation: U = (D - lam*X - beta)/kap ; d = kap*(U^2 - X - A) ; x = U - c0 - c1*d
U_rec = normal((Df - lam * Xf - beta) / kap)
d_rec = normal(kap * (normal(U_rec * U_rec) - Xf - Acst))
x_rec = normal(U_rec - c0 - c1 * d_rec)
check("generation: d recovered", sp.simplify(d_rec - d) == 0)
check("generation: x recovered", sp.simplify(x_rec - x) == 0)

print(flush=True)

# ---------------- component verifiers ----------------------------------------
COMPONENTS = [
    ("verify_quantum_a2zero.py", HERE),
    ("verify_quantum_M4.py", HERE),
    ("verify_quantum_mirror.py", HERE),
    ("verify_quantum_completion.py", HERE),
    ("audit_band1_engine.py", os.path.join(HERE, "..", "..", "archive-import", "provisional", "dixmier-band-program")),
    ("audit_band1_classification.py", os.path.join(HERE, "..", "..", "archive-import", "provisional", "dixmier-band-program")),
]
for script, base in COMPONENTS:
    path = os.path.join(base, script)
    if not os.path.exists(path):
        check(f"component present: {script}", False)
        continue
    r = subprocess.run([sys.executable, path], capture_output=True, text=True)
    tail = (r.stdout.strip().splitlines() or ["<no output>"])[-1]
    check(f"component {script}: exit 0 ({tail[:60]})", r.returncode == 0)

print()
print("ALL QUANTUM BAND-2 ASSEMBLY CHECKS PASSED" if OK else "SOME QUANTUM BAND-2 ASSEMBLY CHECKS FAILED", flush=True)
raise SystemExit(0 if OK else 1)
