#!/usr/bin/env python3
"""Exact verifier for astar-band3.md (DC1 program): the A*-band3 NEGATIVE-TAIL
closure of the CONSTANT-TOP (constant-h) band-3 sector on the Weyl/Dixmier face.

Target statement (shifted-power-residuals.md sec.3-4, RESIDUAL 3; the (kappa2-closure)):
    a3 = 1 (constant top), gauge b3 = 0, wall Q5 => b2 = kappa2 (constant).
    TARGET: exclude genuine Weyl pairs [D,X]=1 with kappa2 != 0. The displayed tame
    family lies in the kappa2 = 0 slice; no converse classification of that slice is claimed.

Conventions (frozen, DC1 house engine):
    A_1[x^{-1}] = (+)_k x^k C[E], E = x d, (x^a f)(x^b g) = x^{a+b} f(E+b) g(E),
    f^[r](E) = f(E+r), T f = f^[1],
    Q_m = sum_{k+l=m} ( b_l^[k] a_k - a_k^[l] b_l ),  [D,X]=1 <=> Q_m = delta_{m0},
    membership E(E-1)...(E-r+1) | a_{-r}, b_{-r}; gauge b_3 = 0;
    G = sum_{i=1}^{3} sum_{r=1}^{i} ( a_i^[-r] b_{-i}^[i-r] - a_{-i}^[i-r] b_i^[-r] ), Q_0 = (T-1)G.

EVIDENCE TIERS
  - committed (default run, exact SymPy over QQ, reproduced in-repo):
      * engine + telescoping potential identities;
      * constant-top positive-cascade + Lemma-P slope + potential-triviality identities;
      * negative-tail mu3-source decompositions + b-bilinearity + generic-X
        non-solvability (the moment/covector obstruction);
      * BOUNDED CERTIFICATE at cap d=1: the kappa2 != 0 sector is the UNIT IDEAL
        (EMPTY) over QQ; explicit kappa2 = 0 tame witness [D,X]=1.
  - HEAVY (behind HEAVY=1): cap d=2, d=3 emptiness (sympy QQ and/or msolve).
  - optional msolve corroboration (SKIPs cleanly if msolve absent): unit ideal over
    several primes at cap d=1 (integer coefficients only; msolve MISPARSES rational
    monomials, so denominators are always cleared -- a real trap documented here).

'Proved (arbitrary degree)' is claimed for the symbolic identities ONLY. The
emptiness of kappa2 != 0 is BOUNDED (cap d), never arbitrary-degree: the
arbitrary-degree closure is OPEN, exactly mirroring the classical e != 0 sibling
(classical-e-nonzero-closure.md), and the memo says so.

Run:  uv run --with sympy python research/dc1-program/verify_astar_band3.py
      HEAVY=1 uv run --with sympy python research/dc1-program/verify_astar_band3.py
The final summary distinguishes no-skip success from executed-check success with skips.
"""
import os, sys, shutil, subprocess, random, tempfile
import sympy as sp

HEAVY = os.environ.get("HEAVY", "") not in ("", "0", "false", "False")
E = sp.symbols("E")
kap, mu3 = sp.symbols("kappa mu3")
PASS = 0
SKIP = 0


def check(cond, label):
    global PASS
    if cond is not True and cond != True:  # noqa: E712
        raise AssertionError("FAIL: " + label)
    PASS += 1
    print("PASS", label, flush=True)


def skip(label):
    global SKIP
    SKIP += 1
    print("SKIP", label, flush=True)


def sh(f, r):
    return sp.expand(sp.sympify(f).subs(E, E + r))


def fall(r):
    return sp.prod([E - i for i in range(r)]) if r > 0 else sp.Integer(1)


def Qm(a, b, m):
    return sp.expand(sum(sh(b[l], k) * a[k] - sh(a[k], l) * b[l]
                         for k in range(-3, 4) for l in range(-3, 4) if k + l == m))


def Gpot(a, b):
    G = 0
    for i in range(1, 4):
        for r in range(1, i + 1):
            G += sh(a[i], -r) * sh(b[-i], i - r) - sh(a[-i], i - r) * sh(b[i], -r)
    return sp.expand(G)


def gp(name, d):
    return sum(sp.Symbol(f"{name}_{i}") * E**i for i in range(d + 1))


# ============================================================================
print("=" * 72)
print("SECTION 0.  Engine re-derivation: Q_m = crossed-product commutator; Q0=(T-1)G")
print("=" * 72)
# direct crossed-product commutator of two band-3 ladder elements with generic
# degree-2 coefficients, compared to the Q_m assembly.
Aq = {k: sum(sp.Symbol(f"A{k+3}_{i}") * E**i for i in range(3)) for k in range(-3, 4)}
Bq = {k: sum(sp.Symbol(f"B{k+3}_{i}") * E**i for i in range(3)) for k in range(-3, 4)}


def ladder_mul(P, Q):
    R = {}
    for a1, fa in P.items():
        for b1, fb in Q.items():
            R[a1 + b1] = R.get(a1 + b1, 0) + sh(fa, b1) * fb
    return {k: sp.expand(v) for k, v in R.items()}


DX = ladder_mul(Bq, Aq)
XD = ladder_mul(Aq, Bq)
for m in [6, 5, 1, 0, -1, -5, -6]:
    direct = sp.expand(DX.get(m, 0) - XD.get(m, 0))
    check(sp.expand(Qm(Aq, Bq, m) - direct) == 0,
          f"Q_{m} equals the direct crossed-product commutator [D,X]_{m}")
# telescoping potential on a random membership-respecting instance
random.seed(3)
def rp(dd=2): return sum(sp.Integer(random.randint(-3, 3)) * E**i for i in range(dd + 1))
At = {k: rp() for k in range(-3, 4)}
Bt = {k: rp() for k in range(-3, 4)}
check(sp.expand(Qm(At, Bt, 0) - (sh(Gpot(At, Bt), 1) - Gpot(At, Bt))) == 0,
      "Q_0 = (T-1) G  with the band-agnostic telescoping potential G (generic)")
# membership => G(0) = 0 (each term carries a membership-vanishing bottom factor)
Am = {3: rp(), 2: rp(), 1: rp(), 0: rp(),
      -1: fall(1) * rp(), -2: fall(2) * rp(), -3: fall(3) * rp()}
Bm = {3: rp(), 2: rp(), 1: rp(), 0: rp(),
      -1: fall(1) * rp(), -2: fall(2) * rp(), -3: fall(3) * rp()}
check(sp.expand(Gpot(Am, Bm)).subs(E, 0) == 0,
      "membership => G(0)=0, so Q_0=1 forces G=E (slope G(1) is the moment)")

# ============================================================================
print()
print("=" * 72)
print("SECTION 1.  Constant top a3=1, gauge b3=0: positive cascade + Lemma-P slope")
print("=" * 72)
a = {k: sp.Function(f"a{k}")(E) for k in range(-3, 4)}
b = {k: sp.Function(f"b{k}")(E) for k in range(-3, 4)}
a[3] = sp.Integer(1)
b[3] = sp.Integer(0)
b[2] = kap
# Q5 wall: for a3=1 it is exactly (T^3-1) b2, so b2 3-periodic => b2 = kappa2 constant
check(sp.expand(Qm(a, b, 5) - (sh(kap, 3) - kap)) == 0,
      "Q5 = (T^3-1) b2 for a3=1  =>  b2 = kappa2 CONSTANT (the gauged wall constant)")
# Q4 identity: the constant-top positive rung
check(sp.expand(Qm(a, b, 4) - (kap * (a[2] - sh(a[2], 2)) + (sh(b[1], 3) - b[1]))) == 0,
      "Q4 = kappa2*(a2 - a2^[2]) + (b1^[3] - b1)   (constant-top positive rung)")
# potential triviality for constant h: the shifted-power factor h^[-1] = 1, so
# G = M with NO nonconstant-h forcing available -- the obstruction is NOT here.
#   (In the nonconstant-h descent, G = h^[-1] M and Q0=1 kills nonconstant h.
#    Here h=1 => h^[-1]=1 => G=M=E is CONSISTENT: no kill from the top.)
check(sp.simplify(sh(sp.Integer(1), -1) - 1) == 0,
      "constant h=1: h^[-1] = 1, so G = h^[-1] M = M -- top potential-factorization "
      "is VACUOUS (no nonconstant-h kill; obstruction sits in the negative tail)")
# Lemma-P slope for the constant top (branch a_-3 != 0 with b_-3 = mu3 a_-3):
Ai = {3: sp.Integer(1), 2: rp(), 1: rp(), 0: rp(),
      -1: fall(1) * rp(), -2: fall(2) * rp(), -3: fall(3) * rp()}
kv = sp.Integer(5)
Bi = {3: sp.Integer(0), 2: kv, 1: rp(), 0: rp(),
      -1: fall(1) * rp(), -2: fall(2) * rp()}
Bi[-3] = sp.expand(mu3 * Ai[-3])
slope_formula = (Ai[1].subs(E, 0) * Bi[-1].subs(E, 1) - Ai[-1].subs(E, 1) * Bi[1].subs(E, 0)
                 + Ai[2].subs(E, 0) * Bi[-2].subs(E, 2) - Ai[-2].subs(E, 2) * kv
                 + mu3 * Ai[-3].subs(E, 3))
check(sp.expand(Gpot(Ai, Bi).subs(E, 1) - slope_formula) == 0,
      "Lemma P slope G(1)=a1(0)b-1(1)-a-1(1)b1(0)+a2(0)b-2(2)-a-2(2)*kappa2+mu3*a-3(3)"
      "  (Q0=1 => this = 1)")

# ============================================================================
print()
print("=" * 72)
print("SECTION 2.  Negative tail: bottom proportionality, inhomogeneous mu3-source,")
print("            b-bilinearity and the moment/covector obstruction")
print("=" * 72)
# Q_-6 => b_-3 = mu3 a_-3 (bottom 3-periodicity), independent of the spent top gauge
bmu = dict(b)
bmu[-3] = mu3 * a[-3]
check(sp.expand(Qm(a, bmu, -6)) == 0, "Q_-6: b_-3 = mu3 a_-3 solves the bottom Wronskian")
# the genuinely-new band-3 feature: Q_-5, Q_-4 carry a mu3-proportional SOURCE with
# no band-2 shadow (the lambda3-mu3 cross-coupling; quantum-band3-cascade.md sec 5)
pw = sh(bmu[-2], -3) * a[-3] - sh(a[-3], -2) * bmu[-2]
src5 = sh(a[-3], -2) * a[-2] - sh(a[-2], -3) * a[-3]
check(sp.expand(Qm(a, bmu, -5) - (pw + mu3 * src5)) == 0,
      "Q_-5 = [bottom wall in b_-2] + mu3*[a_-3,a_-2 source]  (INHOMOGENEOUS)")
d1 = sh(bmu[-1], -3) * a[-3] - sh(a[-3], -1) * bmu[-1]
d2 = sh(a[-3], -1) * a[-1] - sh(a[-1], -3) * a[-3]
d3 = sh(bmu[-2], -2) * a[-2] - sh(a[-2], -2) * bmu[-2]
check(sp.expand(Qm(a, bmu, -4) - (d1 + mu3 * d2 + d3)) == 0,
      "Q_-4 = [b_-1,a_-3 stagger] + mu3*[source] + [b_-2,a_-2 square]  (INHOMOGENEOUS)")
# bilinearity: [D,X] is linear in every b_l and every a_k separately.  Hence for a
# FIXED constant-top X the whole system Q_m=delta is LINEAR in the D-coefficients.
r1, r2 = rp(), rp()
b_lin = {k: sp.Function(f"b{k}")(E) for k in range(-3, 4)}
lhs = Qm(At, {**b_lin, 2: r1 + r2}, 1)
rhs = (Qm(At, {**b_lin, 2: r1}, 1) + Qm(At, {**b_lin, 2: r2}, 1)
       - Qm(At, {**b_lin, 2: sp.Integer(0)}, 1))
check(sp.expand(lhs - rhs) == 0,
      "[D,X] is linear in each b_l (and each a_k): the system is bilinear")
# the moment/covector obstruction: for a GENERIC constant-top X, the b-linear system
# Q_m = delta_{m0} is NOT solvable -- delta_{m0}=1 is outside the D-image; solvability
# is a positive-codimension condition (the cokernel covector on the RHS).
random.seed(11)
def rpi(dd): return sum(sp.Integer(random.randint(-3, 3)) * E**i for i in range(dd + 1))
dX, dD = 1, 3
AX = {3: sp.Integer(1), 2: rpi(dX), 1: rpi(dX), 0: rpi(dX),
      -1: fall(1) * rpi(dX), -2: fall(2) * rpi(dX), -3: fall(3) * rpi(dX)}
bsyms = []
def bp(name, fac=1):
    cs = [sp.Symbol(f"{name}_{i}") for i in range(dD + 1)]
    bsyms.extend(cs)
    return fac * sum(cs[i] * E**i for i in range(dD + 1))
BX = {3: sp.Integer(0), 2: bp("b2"), 1: bp("b1"), 0: bp("b0"),
      -1: bp("bm1", fall(1)), -2: bp("bm2", fall(2)), -3: bp("bm3", fall(3))}
rows = []
for m in range(-6, 7):
    q = sp.expand(Qm(AX, BX, m) - (1 if m == 0 else 0))
    if q != 0:
        rows += sp.Poly(q, E).all_coeffs()
Mmat, rhsv = sp.linear_eq_to_matrix(rows, bsyms)
rM = Mmat.rank()
rAug = Mmat.row_join(rhsv).rank()
check(rAug == rM + 1,
      f"generic constant-top X: b-image rank {rM}, augmented {rAug} => NOT solvable "
      "(an unspecified cokernel covector obstructs; this does not identify the Lemma-P moment functional)")
check(rM == Mmat.cols - 1,
      "b-kernel is 1-dimensional (the constant centralizer D=const); rank = cols-1")

# ============================================================================
print()
print("=" * 72)
print("SECTION 3.  BOUNDED EMPTINESS CERTIFICATE  (both branches, cap d)")
print("=" * 72)


def build_full(d, kv):
    """Full constant-top system Q_m - delta_{m0} at free-poly degree cap d, wall
    constant b2 = kv (integer).  Returns integer-coefficient polynomial equations."""
    A = {3: sp.Integer(1)}
    B = {3: sp.Integer(0), 2: sp.Integer(kv)}
    A[2] = gp("a2", d); A[1] = gp("a1", d); A[0] = gp("a0", d)
    A[-1] = fall(1) * gp("am1", d); A[-2] = fall(2) * gp("am2", d); A[-3] = fall(3) * gp("am3", d)
    B[1] = gp("b1", d); B[0] = gp("b0", d)
    B[-1] = fall(1) * gp("bm1", d); B[-2] = fall(2) * gp("bm2", d); B[-3] = fall(3) * gp("bm3", d)
    eqs = []
    for m in range(-6, 7):
        q = sp.expand(Qm(A, B, m) - (1 if m == 0 else 0))
        if q != 0:
            eqs += [e for e in sp.Poly(q, E).all_coeffs() if e != 0]
    unk = sorted(set().union(*[e.free_symbols for e in eqs]), key=str)
    return eqs, unk


def is_unit_ideal_QQ(eqs, unk):
    G = sp.groebner(eqs, *unk, order="grevlex")
    return list(G) == [sp.Integer(1)]


def _run_msolve(eqs, unk, characteristic, args, timeout_s):
    """Run msolve securely; timeout is an optional-check SKIP, other failures are fatal."""
    vs = ",".join(str(u) for u in unk)
    body = ",\n".join(str(sp.expand(e)).replace("**", "^") for e in eqs)
    with tempfile.TemporaryDirectory(prefix="astar-band3-") as tmp:
        fn = os.path.join(tmp, "input.ms")
        out = os.path.join(tmp, "output.txt")
        with open(fn, "w", encoding="utf-8") as fh:
            fh.write(f"{vs}\n{characteristic}\n{body}\n")
        try:
            rr = subprocess.run(["msolve", *args, "-f", fn, "-o", out],
                                capture_output=True, text=True, timeout=timeout_s)
        except subprocess.TimeoutExpired:
            return None
        if rr.returncode != 0:
            raise RuntimeError(f"msolve failed with status {rr.returncode}: {rr.stderr.strip()}")
        if not os.path.exists(out):
            raise RuntimeError(f"msolve produced no output file; stderr: {rr.stderr.strip()}")
        with open(out, encoding="utf-8") as fh:
            return fh.read()


def msolve_unit(eqs, unk, prime):
    """msolve -g 2 over F_prime; None means timeout, malformed output is fatal."""
    txt = _run_msolve(eqs, unk, prime, ["-g", "2"], 300)
    if txt is None:
        return None
    lines = [ln for ln in txt.splitlines() if not ln.lstrip().startswith("#")]
    parsed = "".join(lines).replace(" ", "")
    if not parsed.startswith("["):
        raise RuntimeError(f"malformed msolve finite-field output: {parsed[:200]!r}")
    return parsed.startswith("[1]:") or parsed == "[1]"


def msolve_empty_QQ(eqs, unk, timeout_s):
    """msolve over QQ; None means timeout, malformed output is fatal."""
    txt = _run_msolve(eqs, unk, 0, [], timeout_s)
    parsed = txt.strip() if txt is not None else None
    if parsed is None:
        return None
    if parsed.startswith("[-1]"):
        return True
    if parsed.startswith("["):
        return False
    raise RuntimeError(f"malformed msolve characteristic-zero output: {parsed[:200]!r}")


# --- kappa2 = 0 branch: an explicit tame witness (nonempty control) ---
# positive control U = x + c1*d, X = U^3 - d/kap, D = kap*U  (kap=1): a3=1, b2=0.
c1 = sp.Integer(2)
U = {1: sp.Integer(1), -1: c1 * E}
U2 = ladder_mul(U, U); U3 = ladder_mul(U2, U)
Xw = {k: sp.expand(v) for k, v in U3.items()}
Xw[-1] = sp.expand(Xw.get(-1, 0) - E)          # - d/kap with kap=1, d = x^{-1}E
Dw = {k: sp.expand(v) for k, v in U.items()}   # kap=1
Aw = {k: Xw.get(k, sp.Integer(0)) for k in range(-3, 4)}
Bw = {k: Dw.get(k, sp.Integer(0)) for k in range(-3, 4)}
check(all(sp.expand(Qm(Aw, Bw, m) - (1 if m == 0 else 0)) == 0 for m in range(-6, 7)),
      "kappa2 = 0 tame witness: explicit [D,X]=1 (positive control), a3=1, b2=0")
check(sp.expand(Bw[2]) == 0 and sp.expand(Aw[3]) == 1,
      "  witness has b2 = kappa2 = 0 (the slice contains this tame family; no converse claimed)")

# --- kappa2 != 0 branch: EMPTY at cap d=1 (normalize kappa2 = 1) ---
print("  building cap d=1 kappa2=1 system ...", flush=True)
eqs1, unk1 = build_full(1, 1)
print(f"  cap d=1: {len(eqs1)} eqs, {len(unk1)} unknowns; computing QQ Groebner ...", flush=True)
check(is_unit_ideal_QQ(eqs1, unk1),
      "COMMITTED: kappa2!=0 constant-top sector at cap d=1 is the UNIT IDEAL (EMPTY over QQ)")
# sanity: the kappa2 = 0 system at cap d=1 is NOT unit (tame family present)
eqs1z, unk1z = build_full(1, 0)
check(not is_unit_ideal_QQ(eqs1z, unk1z),
      "control: kappa2 = 0 at cap d=1 is NOT the unit ideal (tame pairs exist)")

# --- msolve corroboration (multi-prime, integer coeffs); SKIP if absent ---
if shutil.which("msolve") is None:
    skip("msolve corroboration (msolve not in PATH)")
else:
    ok = True
    seen = False
    for p in (32003, 1000003, 2147483647):
        r = msolve_unit(eqs1, unk1, p)
        if r is None:
            skip(f"msolve F_{p} corroboration (timed out)")
        else:
            seen = True
            ok = ok and r
    if seen:
        check(ok, "msolve corroboration: kappa2!=0 cap d=1 is unit ideal over each tested prime")

# --- HEAVY: cap d=2, d=3 (msolve char-0 over QQ; rigorous; runtimes are minutes-hours) ---
# NOTE (runtime): cap d=2 msolve char-0 took ~35 min in development; cap d=3 may
# exceed an hour and is attempted with a hard cap, SKIPping cleanly on timeout.
if not HEAVY:
    skip("HEAVY cap d=2 emptiness (set HEAVY=1; msolve char-0 over QQ, ~35 min)")
    skip("HEAVY cap d=3 emptiness (set HEAVY=1; msolve char-0, may exceed 1 h)")
elif shutil.which("msolve") is None:
    skip("HEAVY cap d=2,3 emptiness (msolve not in PATH)")
else:
    for d, tmo in ((2, 3000), (3, 5400)):
        print(f"  [HEAVY] building cap d={d} kappa2=1 ...", flush=True)
        eqs, unk = build_full(d, 1)
        print(f"  [HEAVY] cap d={d}: {len(eqs)} eqs {len(unk)} unk; msolve char-0 QQ "
              f"(cap {tmo}s) ...", flush=True)
        r = msolve_empty_QQ(eqs, unk, tmo)
        if r is None:
            skip(f"HEAVY cap d={d} emptiness (msolve char-0 timed out at {tmo}s)")
        else:
            check(r, f"HEAVY: kappa2!=0 constant-top sector at cap d={d} is EMPTY over QQ "
                     "(msolve char-0 = [-1])")

# ============================================================================
print()
print(f"({PASS} checks passed, {SKIP} skipped)")
if SKIP:
    print("ALL EXECUTED ASTAR-BAND3 (DC1) CHECKS PASSED; OPTIONAL CHECKS SKIPPED")
else:
    print("ALL ASTAR-BAND3 (DC1) CHECKS PASSED; NO SKIPS")
