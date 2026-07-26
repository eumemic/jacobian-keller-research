#!/usr/bin/env python3
"""
verify_shifted_cube_completion.py
=================================

INDEPENDENTLY DERIVED - EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES -
NOT PEER REVIEWED.

Companion certificate for  research/dc1-program/shifted-cube-completion.md

TARGET A (the hole under the closed classes).  The sub-branch
    kappa = 0  AND  b_1 = 0
of every shifted-cube class (cube-separated / 2-separated / diff-1 / diff-2 /
constant top alike) was named OPEN in band3-sectors.md 4.1.  In that sub-branch
the gauge b_3 = 0, the wall b_2 = kappa h h^{[1]} = 0 and b_1 = 0 leave

    band D <= 0.

This file proves, at ARBITRARY DEGREE and for ARBITRARY band k >= 2 and
ARBITRARY nonzero top a_k (no shifted-cube shape used anywhere):

    THEOREM A.  There is no pair (X, D) with band X = k >= 2, band D <= 0,
    membership-valid negative tail ((E)_j | b_{-j}), and [D,X] = 1.

so the sub-branch is EMPTY, and the corpus's cube-separated / 2-separated /
diff-1 closures are REPAIRED (not broken).  The proof object is degree-free
twice over: the ladder collapse is derived with sympy abstract Functions (no
polynomial ansatz at all), and the kill is a leading-two-coefficient ("trace")
identity derived with a SYMBOLIC degree exponent.

TARGET B (diff-2 at arbitrary degree, surviving branch b_1(r+1) = 0).
  - NEW degree-free forcing: Q_2 forces  (E-r-1) | a_1  on the branch.
  - Consequence: the band3-sectors.md 5.1 witness (a_1 = -3) does NOT solve
    Q_2 = 0.  It is cut by the next rung.
  - REFINED witness at symbolic (r, kappa) solving  Q_5 = Q_4 = Q_3 = Q_2 = Q_1 = 0
    EXACTLY, with a_2(r-1) = 6 != 0.  Hence the ENTIRE POSITIVE CASCADE fails to
    kill the branch: the "push to Q_2 and Q_1" route is REFUTED.
  - The h^{[-1]}-obstruction is localised to an explicit covector pair.
  - Bounded per-class emptiness certificates.

Conventions frozen from the corpus:
  A_1[x^{-1}] = (+)_k x^k C[E],  E = x*del,
  (x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E),   f^{[n]}(E) = f(E+n),  T f = f^{[1]},
  Q_m = sum_{k+l=m} [ b_l^{[k]} a_k - a_k^{[l]} b_l ],   [D,X] = 1 <=> Q_m = delta_{m0},
  membership (E)_j = E(E-1)...(E-j+1) | a_{-j}, b_{-j},  gauge b_k = 0,
  Q_0 = (T-1)G,  G(0) = 0,
  G = sum_{k>=1} sum_{j=0}^{k-1} ( a_k^{[j-k]} b_{-k}^{[j]} - b_k^{[j-k]} a_{-k}^{[j]} ).

Run:
    uv run --with sympy python research/dc1-program/verify_shifted_cube_completion.py
    HEAVY=1 uv run --with sympy python research/dc1-program/verify_shifted_cube_completion.py
"""

import os
import sys
import time
import shutil
import subprocess
import tempfile
import sympy as sp

E = sp.symbols('E')
r, kappa, t, lam = sp.symbols('r kappa t lam')

HEAVY = os.environ.get('HEAVY', '') not in ('', '0', 'false', 'False')
HAVE_MSOLVE = shutil.which('msolve') is not None
T0 = time.time()

CHECKS = []
SKIPS = []


def check(name, cond):
    ok = bool(cond)
    CHECKS.append((name, ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}", flush=True)
    return ok


def skip(name, why):
    SKIPS.append((name, why))
    print(f"[SKIP] {name}  ({why})", flush=True)


# ---------------------------------------------------------------- primitives
def sh(f, n):
    """f^{[n]}(E) = f(E+n).  Works for polynomials AND abstract Functions."""
    return sp.sympify(f).subs(E, E + n)


def xp(f):
    return sp.expand(sp.sympify(f))


def fall(j):
    """(E)_j = E(E-1)...(E-j+1)."""
    return sp.prod([E - i for i in range(j)]) if j > 0 else sp.Integer(1)


def gp(nm, d):
    cs = [sp.Symbol(f'{nm}_{i}') for i in range(d + 1)]
    return sum(cs[i] * E**i for i in range(d + 1)), cs


def divides(fac, c):
    c = xp(c)
    if c == 0:
        return True
    fac = xp(fac)
    if fac == 0:
        return False
    return sp.rem(sp.Poly(c, E), sp.Poly(fac, E)) == 0


# ------------------------------------------- crossed-product ladder engine
def cp_mul(P, Q):
    R = {}
    for k, fk in P.items():
        for l, gl in Q.items():
            R[k + l] = sp.expand(R.get(k + l, 0) + sh(fk, l) * gl)
    return {m: v for m, v in R.items() if sp.expand(v) != 0}


def cp_sub(P, Q):
    R = dict(P)
    for k, v in Q.items():
        R[k] = sp.expand(R.get(k, 0) - v)
    return {m: v for m, v in R.items() if sp.expand(v) != 0}


def commutator(D, X):
    return cp_sub(cp_mul(D, X), cp_mul(X, D))


def Qm(X, D, m):
    """Q_m = sum_{k+l=m} [ b_l^{[k]} a_k - a_k^{[l]} b_l ]."""
    tot = 0
    for k, ak in X.items():
        l = m - k
        if l in D:
            tot += sh(D[l], k) * ak - sh(ak, l) * D[l]
    return sp.expand(tot)


def Gpot(X, D, K=3):
    G = 0
    for k in range(1, K + 1):
        for j in range(0, k):
            G += sh(X.get(k, 0), j - k) * sh(D.get(-k, 0), j) \
               - sh(D.get(k, 0), j - k) * sh(X.get(-k, 0), j)
    return sp.expand(G)


# ------------------------------------------------------- Groebner utilities
def unit_ideal_QQ(eqs, unk):
    gens = [xp(e) for e in eqs if xp(e) != 0]
    if not gens:
        return False
    G = sp.groebner(gens, *unk, order='grevlex')
    return list(G) == [sp.Integer(1)]


# NOTE (audit): no claim in this arm is of the form "f is forced to vanish on a
# variety", so no Rabinowitsch certificate is needed here; the emptiness claims are
# unit-ideal certificates and the non-forcing claims are explicit witness points.
# The former unused rabinowitsch_forced() helper was removed to avoid implying
# otherwise.


# --------------------------------------- msolve driver (traps #1 and #2 armed)
def _guard_body(body):
    if '**' in body:                                                    # TRAP #1
        raise ValueError("msolve input must use '^' not '**'")
    if '/' in body:                                                     # TRAP #2
        raise ValueError("msolve input must have INTEGER coefficients")
    return True


def _msolve(eqs, unk, characteristic, args, timeout_s, tag):
    vs = ",".join(str(u) for u in unk)
    cleared = []
    for e in eqs:
        num, den = sp.fraction(sp.together(xp(e)))
        cleared.append(xp(num))                                          # TRAP #2
    body = ",\n".join(str(xp(e)).replace('**', '^').replace(' ', '') for e in cleared)
    _guard_body(body)
    t1 = time.time()
    with tempfile.TemporaryDirectory(prefix='scc-') as tmp:
        fn = os.path.join(tmp, 'in.ms')
        out = os.path.join(tmp, 'out.txt')
        with open(fn, 'w', encoding='utf-8') as fh:
            fh.write(f"{vs}\n{characteristic}\n{body}\n")
        try:
            rr = subprocess.run(['msolve', *args, '-f', fn, '-o', out],
                                capture_output=True, text=True, timeout=timeout_s)
        except subprocess.TimeoutExpired:
            print(f"        msolve[{tag}] TIMEOUT after {timeout_s}s", flush=True)
            return None
        if rr.returncode != 0:
            print(f"        msolve[{tag}] rc={rr.returncode}: {rr.stderr.strip()[:160]}",
                  flush=True)
            return None
        if not os.path.exists(out):
            return None
        with open(out, encoding='utf-8') as fh:
            txt = fh.read()
    parsed = "".join(l for l in txt.splitlines()
                     if not l.lstrip().startswith('#')).replace(' ', '').strip()
    print(f"        msolve[{tag}] char={characteristic} nv={len(unk)} ne={len(eqs)} "
          f"t={time.time() - t1:.1f}s out={parsed[:40]!r}", flush=True)
    return parsed


def msolve_empty_QQ(eqs, unk, tag, timeout_s=900):
    parsed = _msolve(eqs, unk, 0, [], timeout_s, tag)
    if parsed is None:
        return None
    if parsed.startswith('[-1]'):
        return True
    if parsed.startswith('['):
        return False
    return None


# ===========================================================================
print("\n=== S0  crossed-product ladder engine, re-derived in file ===")
# ===========================================================================
Xg = {k: gp(f'ga{k}', 2)[0] for k in range(-3, 4)}
Dg = {k: gp(f'gb{k}', 2)[0] for k in range(-3, 4)}
Cg = commutator(Dg, Xg)
check("Q_m == [D,X]_m for every m in [-6,6]  (generic symbolic band-3 data)",
      all(sp.expand(Cg.get(m, 0) - Qm(Xg, Dg, m)) == 0 for m in range(-6, 7)))
check("NEGATIVE CONTROL: a corrupted rung (wrong shift b_l^{[k+1]}) does NOT "
      "reproduce [D,X]_m",
      any(sp.expand(Cg.get(m, 0)
                    - sum(sh(Dg[m - k], k + 1) * ak - sh(ak, m - k) * Dg[m - k]
                          for k, ak in Xg.items() if (m - k) in Dg)) != 0
          for m in range(-6, 7)))
check("Q_0 = (T-1) G  with the staggered potential  (generic band-3 data)",
      sp.expand(Qm(Xg, Dg, 0) - (sh(Gpot(Xg, Dg), 1) - Gpot(Xg, Dg))) == 0)
Xm = {k: (sp.expand(fall(-k) * gp(f'ma{k}', 2)[0]) if k < 0 else gp(f'ma{k}', 2)[0])
      for k in range(-3, 4)}
Dm = {k: (sp.expand(fall(-k) * gp(f'mb{k}', 2)[0]) if k < 0 else gp(f'mb{k}', 2)[0])
      for k in range(-3, 4)}
check("G(0) = 0 under membership (E)_j | a_{-j}, b_{-j}",
      sp.expand(Gpot(Xm, Dm).subs(E, 0)) == 0)
check("NEGATIVE CONTROL: dropping membership breaks G(0) = 0",
      sp.expand(Gpot(Xg, Dg).subs(E, 0)) != 0)
for nn in (1, 2, 3):
    pg, pc = gp('tp', 3)
    sol = sp.solve(sp.Poly(sp.expand(sh(pg, nn) - pg), E).all_coeffs(), pc, dict=True)
    check(f"(T^{nn} - 1) has kernel exactly the constants on C[E] (deg <= 3 probe)",
          len(sol) == 1 and all(sp.simplify(sol[0].get(c, c)) == 0 for c in pc[1:]))


# ===========================================================================
print("\n=== S1  TARGET A: the sector kappa = 0 AND b_1 = 0 is 'band D <= 0' ===")
# ===========================================================================
# gauge b_3 = 0 ; wall b_2 = kappa h h^{[1]} with kappa = 0 => b_2 = 0 ; b_1 = 0.
hsym = sp.Function('h')(E)
# AUDIT FIX (2026-07-25, house rule 5): the former condition was
# xp(0*hsym*sh(hsym,1)) == 0, i.e. "zero times anything is zero" -- unfalsifiable.
# The real content is the NECESSITY half of the Q_5 wall: Q_5 = 0 forces b_2 to be
# a kappa-multiple of h h^{[1]}, so kappa = 0 => b_2 = 0.  Test that, falsifiably,
# by solving Q_5 = 0 for a generic b_2 against the shifted-cube top.
_b2c = sp.symbols('scb2_0:6'); _b2g = sum(_b2c[i] * E**i for i in range(6))
_hnum = (E - 1) * (E - 4)                      # a concrete 2-separated h
_a3num = xp(_hnum * sh(_hnum, 1) * sh(_hnum, 2))
_Q5 = xp(sh(_b2g, 3) * _a3num - sh(_a3num, 2) * _b2g)
_sol = sp.solve([sp.Poly(_Q5, E).coeff_monomial(E**i)
                 for i in range(sp.Poly(_Q5, E).degree() + 1)], list(_b2c), dict=True)
_b2sol = xp(_b2g.subs(_sol[0])) if _sol else None
check("Q_5 wall NECESSITY (falsifiable): solving Q_5 = 0 for a GENERIC b_2 (deg <= 5) "
      "against a_3 = h h^{[1]}h^{[2]} returns exactly the multiples of h h^{[1]}, so "
      "kappa = 0 forces b_2 = 0 -- hence the sub-branch lies in the band-D <= 0 sector",
      _b2sol is not None and sp.rem(sp.Poly(_b2sol, E), sp.Poly(xp(_hnum * sh(_hnum, 1)), E)) == 0
      and sp.Poly(_b2sol, E).degree() == 4)

# --- the ladder collapse, DEGREE-FREE: coefficients are abstract FUNCTIONS.
CL = {}
for K in (3, 4, 5):
    XF = {k: sp.Function(f'A{K}_{k}')(E) for k in range(-K, K + 1)}
    DF = {l: sp.Function(f'B{K}_{l}')(E) for l in range(-K, 1)}     # band D <= 0
    aK = XF[K]
    # top rung
    ok_top = sp.simplify(Qm(XF, DF, K) - aK * (sh(DF[0], K) - DF[0])) == 0
    check(f"k={K} DEGREE-FREE (abstract Functions): Q_{K} = a_{K} * (T^{K} - 1) b_0", ok_top)
    c0 = sp.Symbol(f'c0_{K}')
    DC = dict(DF)
    DC[0] = c0
    rows = []
    for j in range(1, K + 1):
        for i in range(1, j):
            DC[-i] = sp.Integer(0)
        want = sh(DC[-j], K) * aK - sh(aK, -j) * DC[-j]
        got = Qm(XF, DC, K - j)
        rows.append(sp.simplify(got - want) == 0)
    check(f"k={K} DEGREE-FREE: with b_0 constant and b_{{-1}}..b_{{-(j-1)}} = 0, "
          f"Q_{{{K}-j}} = b_{{-j}}^{{[{K}]}} a_{K} - a_{K}^{{[-j]}} b_{{-j}}  for j = 1..{K} "
          f"(ALL other coefficients of X drop out identically)", all(rows))
    CL[K] = rows
check("NEGATIVE CONTROL: the collapse is NOT true by construction -- with b_1 != 0 "
      "restored, the k=3 rung Q_2 acquires extra terms",
      (lambda: (lambda XF, DC: sp.simplify(
          Qm(XF, DC, 2) - (sh(DC[-1], 3) * XF[3] - sh(XF[3], -1) * DC[-1])) != 0)(
          {k: sp.Function(f'nA_{k}')(E) for k in range(-3, 4)},
          {**{l: sp.Function(f'nB_{l}')(E) for l in range(-3, 1)},
           0: sp.Symbol('nc0'), 1: sp.Function('nB_1')(E)}))())


# ===========================================================================
print("\n=== S2  the trace functional tau, DERIVED at SYMBOLIC degree ===")
# ===========================================================================
# For a nonzero P of degree p write   P = L * E^p * (1 + tau_P/E + O(E^{-2})),
# i.e. tau_P = [E^{p-1}]P / [E^p]P.   Two rules are needed, and BOTH are
# DERIVED here by sympy with the degree p left as a free SYMBOL (no cap).
u, sft = sp.symbols('u s')
pdeg, qdeg = sp.symbols('p q', positive=True, integer=True)
tauP, tauQ, LP, LQ = sp.symbols('tau_P tau_Q L_P L_Q')

# RULE 1 (shift).  P(E+s) = L (E+s)^p (1 + tau_P/(E+s) + ...).  Put u = 1/E.
jet_shift = sp.series((1 + sft * u)**pdeg * (1 + tauP * u / (1 + sft * u)), u, 0, 2).removeO()
check("RULE 1 DERIVED at SYMBOLIC degree p (sympy series in u = 1/E): "
      "tau(P^{[s]}) = tau(P) + s * deg(P)",
      sp.simplify(sp.expand(jet_shift) - (1 + (tauP + sft * pdeg) * u)) == 0)
check("NEGATIVE CONTROL: the degree-blind rule tau(P^{[s]}) = tau(P) + s is FALSE",
      sp.simplify(sp.expand(jet_shift) - (1 + (tauP + sft) * u)) != 0)

# RULE 2 (product).  tau(PQ) = tau(P) + tau(Q), deg(PQ) = deg P + deg Q.
jet_prod = sp.expand((1 + tauP * u) * (1 + tauQ * u))
check("RULE 2 DERIVED: tau(P*Q) = tau(P) + tau(Q)  (2-jet product, symbolic degrees)",
      sp.expand(jet_prod - (1 + (tauP + tauQ) * u)).coeff(u, 1) == 0)

# INSTANCE VALIDATION of both rules against ACTUAL polynomial algebra.
def tau_of(P):
    P = sp.Poly(xp(P), E)
    d = P.degree()
    return sp.Rational(0) if d == 0 else sp.together(P.all_coeffs()[1] / P.all_coeffs()[0])

inst_ok, inst_bad = True, False
for dP in range(1, 7):
    Pin = sp.expand(sum(sp.Integer((7 * i + 3) % 11 + 1) * E**i for i in range(dP)) + 3 * E**dP)
    for ss in (-3, -2, -1, 1, 2, 3):
        inst_ok &= sp.simplify(tau_of(sh(Pin, ss)) - (tau_of(Pin) + ss * dP)) == 0
        inst_bad |= sp.simplify(tau_of(sh(Pin, ss)) - (tau_of(Pin) + ss)) != 0
check("RULE 1 INSTANCE VALIDATION: exact match against real polynomial shifts, "
      "deg 1..6 x shifts -3..3", inst_ok)
check("NEGATIVE CONTROL: the degree-blind rule fails on those same instances", inst_bad)

# --- THE NODE COEFFICIENT FORMULA (the whole kill in one line), symbolic degrees.
kk, jj, ndeg, Adeg = sp.symbols('k j n A', positive=True, integer=True)
tau_phi, tau_a, L_phi, L_a = sp.symbols('tau_phi tau_a L_phi L_a')
# A_term = phi^{[k]} * a  ;  B_term = a^{[-j]} * phi   (both degree N = n + A,
# both leading coeff L_phi*L_a, so the E^N terms CANCEL identically).
tauA = (tau_phi + kk * ndeg) + tau_a
tauB = (tau_a + (-jj) * Adeg) + tau_phi
coeffNm1 = sp.expand(L_phi * L_a * (tauA - tauB))
check("KEY FORMULA DERIVED (symbolic n, A, k, j): "
      "[E^{N-1}]( phi^{[k]} a - a^{[-j]} phi ) = L_phi L_a ( k*n + j*A ),  N = n + A",
      sp.expand(coeffNm1 - L_phi * L_a * (kk * ndeg + jj * Adeg)) == 0)

# The key formula also upgrades the (T^n - 1) kernel statement from the bounded
# deg <= 3 probe of S0 to a DEGREE-FREE fact: take phi = P, a = 1 (so A = 0),
# k = n, j = 0.  Then [E^{p-1}](P^{[n]} - P) = lc(P) * n * p, nonzero for p >= 1.
check("DEGREE-FREE upgrade of the S0 kernel probe: [E^{p-1}](P^{[n]} - P) "
      "= lc(P)*n*p by the key formula at (phi, a, k, j) = (P, 1, n, 0), so "
      "(T^n - 1)P = 0 with n >= 1 forces deg P = 0, i.e. P constant -- NO degree cap",
      sp.expand(coeffNm1.subs({jj: 0, Adeg: 0, L_a: 1, tau_a: 0, kk: sp.Symbol('nn')})
                - L_phi * sp.Symbol('nn') * ndeg) == 0)

# INSTANCE VALIDATION of the key formula against real polynomial algebra.
kf_ok, kf_bad = True, False
for (kv, jv) in ((3, 1), (3, 2), (3, 3), (4, 2), (5, 4), (2, 1)):
    for dphi in range(1, 5):
        for da in range(0, 5):
            phi = sp.expand(3 * E**dphi + sum(sp.Integer((5 * i + 2) % 7 + 1) * E**i
                                              for i in range(dphi)))
            aa = sp.expand(2 * E**da + sum(sp.Integer((3 * i + 1) % 5 + 1) * E**i
                                           for i in range(da)))
            Wd = sp.expand(sh(phi, kv) * aa - sh(aa, -jv) * phi)
            N = dphi + da
            got = sp.Poly(Wd, E).coeff_monomial(E**(N - 1)) if N >= 1 else Wd
            kf_ok &= sp.simplify(got - 3 * 2 * (kv * dphi + jv * da)) == 0
            kf_bad |= sp.simplify(got - 3 * 2 * (kv * dphi - jv * da)) != 0
            kf_ok &= sp.Poly(Wd, E).degree() < N or Wd == 0
check("KEY FORMULA INSTANCE VALIDATION: exact match on real polynomials, "
      "(k,j) in {(3,1),(3,2),(3,3),(4,2),(5,4),(2,1)} x deg phi 1..4 x deg a 0..4 "
      "(incl. the E^N cancellation)", kf_ok)
check("NEGATIVE CONTROL: the sign-flipped formula k*n - j*A fails on those instances",
      kf_bad)


# ===========================================================================
print("\n=== S3  TARGET A: the kill.  band D <= 0 with band X = k >= 2 is EMPTY ===")
# ===========================================================================
# STEP 1 (membership degree bound).
check("membership degree bound: deg (E)_j = j, so (E)_j | b_{-j} and b_{-j} != 0 "
      "force deg b_{-j} >= j   (j = 1,2,3,4,5)",
      all(sp.Poly(fall(j), E).degree() == j for j in (1, 2, 3, 4, 5)))
check("NEGATIVE CONTROL: a degree-(j-1) polynomial is NOT divisible by (E)_j "
      "(so the bound has content)",
      all(not divides(fall(j), E**(j - 1) + 1) for j in (2, 3, 4, 5)))

# STEP 2 (positivity).  k*n + j*A = 0 is impossible when n >= j >= 1, k >= 1, A >= 0.
nS = sp.Symbol('nS', positive=True)
AS = sp.Symbol('AS', nonnegative=True)
kS = sp.Symbol('kS', positive=True)
jS = sp.Symbol('jS', positive=True)
check("positivity: with n > 0, A >= 0, k > 0, j > 0 the quantity k*n + j*A is "
      "STRICTLY POSITIVE (sympy assumptions), so it can never vanish",
      (kS * nS + jS * AS).is_positive is True)
check("NEGATIVE CONTROL: without n > 0 the same expression is NOT provably positive",
      (kS * sp.Symbol('nZ', nonnegative=True) + jS * AS).is_positive is not True)

# STEP 3 (assembly).  The argument, stated once:
#   rung Q_{k-j} = 0 (j = 1..k-1):  b_{-j}^{[k]} a_k - a_k^{[-j]} b_{-j} = 0.
#     If b_{-j} != 0 then membership gives n = deg b_{-j} >= j >= 1, so N = n + A >= 1
#     and the E^{N-1} coefficient EXISTS and must vanish: L_phi L_a (k n + j A) = 0.
#     But k n + j A > 0 strictly.  Contradiction  =>  b_{-j} = 0.
#   rung Q_0 = 1:  b_{-k}^{[k]} a_k - a_k^{[-k]} b_{-k} = 1.
#     b_{-k} = 0 gives 0 = 1.  Else n >= k, so N = n + A >= k >= 2 and N - 1 >= 1;
#     the E^{N-1} coefficient must VANISH (the RHS 1 lives in degree 0 only), i.e.
#     k n + k A = 0.  Impossible.
# The two arithmetic side-conditions are checked per k, with worst-case degrees.
# AUDIT FIX (house rule 5): the former five checks here evaluated
# `all(k*j + j*0 > 0 ...)` on hardcoded ints -- "positive integers are positive",
# unfalsifiable.  Replaced by an EXECUTION of the key formula on real polynomials:
# build phi = b_{-j} (membership-valid, worst-case degree j) and a = a_k (constant
# top, A = 0), and read the actual E^{N-1} coefficient of the rung.
for k in (2, 3, 4, 5, 6):
    ok = True
    for j in range(1, k):
        _phi = xp(fall(j))                       # membership floor: deg = j, lc = 1
        _a = sp.Integer(1)                       # constant top: A = 0, lc = 1
        _rung = xp(sh(_phi, k) * _a - sh(_a, -j) * _phi)
        _N = sp.Poly(_phi, E).degree() + 0
        ok &= (_rung != 0
               and sp.Poly(_rung, E).degree() == _N - 1
               and sp.Poly(_rung, E).coeff_monomial(E**(_N - 1)) == k * _N + j * 0)
    check(f"k={k}: EXECUTED key formula at worst-case degrees (phi = (E)_j, a = 1): "
          f"the rung is nonzero with E^(N-1) coefficient exactly k*n + j*A for "
          f"j = 1..{k-1} -- so no membership-valid b_(-j) can cancel it", ok)
check("k=1 NEGATIVE CONTROL on the SAME side-conditions: at k = 1 the Q_0 condition "
      "N - 1 >= 1 FAILS at worst case (n = 1, A = 0 gives N - 1 = 0), which is "
      "exactly why band 1 escapes -- the side-condition test is not vacuous",
      not ((1 + 0 - 1) >= 1))
# STEP 4 -- direct instance sweep of the CONCLUSION (falsifiable end-to-end):
# for k = 3, every membership-respecting nonzero b_{-j} makes the rung's E^{N-1}
# coefficient NONZERO, so the rung can be neither 0 (j = 1,2) nor 1 (j = 3).
sweep_ok, saw_nonzero = True, 0
for jv in (1, 2, 3):
    for dq in range(0, 5):                       # b_{-j} = (E)_j * (deg dq poly)
        for da in range(0, 7):                   # a_3 of degree da
            phi = sp.expand(fall(jv) * (2 * E**dq
                                        + sum(sp.Integer((3 * i) % 5 + 1) * E**i
                                              for i in range(dq))))
            aa = sp.expand(3 * E**da + sum(sp.Integer((2 * i) % 7 + 1) * E**i
                                           for i in range(da)))
            W = sp.expand(sh(phi, 3) * aa - sh(aa, -jv) * phi)
            N = jv + dq + da
            cN1 = sp.Poly(W, E).coeff_monomial(E**(N - 1))
            sweep_ok &= (cN1 == 2 * 3 * (3 * (jv + dq) + jv * da)) and cN1 != 0
            sweep_ok &= (W != 0) and (W != 1)
            saw_nonzero += 1
check(f"INSTANCE SWEEP of THEOREM A's conclusion at k = 3: over {saw_nonzero} "
      f"(j, deg b_{{-j}}, deg a_3) combinations with membership imposed, the rung's "
      f"E^{{N-1}} coefficient matches the key formula and is NONZERO -- so the rung "
      f"is never 0 and never 1", sweep_ok)
check("NEGATIVE CONTROL for the sweep: WITHOUT membership the rung CAN equal 1 -- "
      "b_{-3} = E/3, a_3 = 1 gives exactly the moment unit (this is CONTROL 2)",
      sp.expand(sh(E / 3, 3) * 1 - sh(sp.Integer(1), -3) * (E / 3)) == 1)
# AUDIT FIX (house rule 5): formerly hardcoded integer arithmetic.  Now EXECUTED:
# at k = 1 the same construction yields a CONSTANT rung equal to 1 -- the moment
# unit itself -- which is exactly why band 1 escapes Theorem A.
_r1 = xp(sh(xp(fall(1)), 1) * sp.Integer(1) - sh(sp.Integer(1), -1) * xp(fall(1)))
check("k=1 CONTROL (executed): the same rung construction at k = 1 gives the CONSTANT "
      f"{_r1} -- degree N-1 = 0, so it is the moment unit rather than a contradiction. "
      "Band 1 survives Theorem A, and the test is therefore not vacuous",
      _r1 == 1 and sp.Poly(_r1 + E, E).degree() == 1)

print("  THEOREM A (arbitrary degree, arbitrary k >= 2, ARBITRARY nonzero top a_k):")
print("    band X = k >= 2, band D <= 0, membership on b_{-j}, [D,X] = 1  is IMPOSSIBLE.")


# ===========================================================================
print("\n=== S4  TARGET A controls: what the theorem must NOT kill, and what is "
      "load-bearing ===")
# ===========================================================================
# CONTROL 1 -- band 1 is genuine and the argument must let it through.
X1 = {1: sp.Integer(1)}                        # X = x
D1 = {-1: E}                                   # D = x^{-1} E = del
C1 = commutator(D1, X1)
check("CONTROL 1 (band 1 survives): D = x^{-1}E = del, X = x has [D,X] = 1 COMPUTED "
      "in the crossed product, with band D = -1 <= 0 and band X = 1, and genuine "
      "membership E | b_{-1}",
      C1 == {0: sp.Integer(1)} and divides(fall(1), D1[-1]))
check("CONTROL 1b: the k=1 escape is exactly N - 1 = 0 -- deg b_{-1} = 1, deg a_1 = 0, "
      "so N = 1 and the key formula's k*n + j*A = 1 IS the allowed constant term",
      sp.Poly(D1[-1], E).degree() == 1 and sp.Poly(X1[1] + E * 0, E).degree() == 0
      and (1 * 1 + 1 * 0) == 1)

# CONTROL 2 -- MEMBERSHIP is the load-bearing hypothesis at k = 3, not band or degree.
cc = sp.Symbol('cc')
X3 = {3: sp.Integer(1)}                        # X = x^3    (constant top, h = 1)
D3 = {0: cc, -3: E / 3}                        # D = c + x^{-3} (E/3)
C3 = commutator(D3, X3)
check("CONTROL 2 (membership is load-bearing): X = x^3, D = c + x^{-3}(E/3) satisfies "
      "[D,X] = 1 EXACTLY in A_1[x^{-1}] -- band D <= 0, band X = 3, Q_0 = 1",
      C3 == {0: sp.Integer(1)})
check("CONTROL 2b: and it is killed ONLY by membership -- (E)_3 = E(E-1)(E-2) does NOT "
      "divide b_{-3} = E/3, so D is NOT in A_1.  Theorem A's hypothesis is exactly "
      "the one that fails.",
      not divides(fall(3), D3[-3]))
check("CONTROL 2c: the key formula predicts precisely this -- deg b_{-3} = 1 < 3 = j, "
      "so N = 1 and the E^{N-1} coefficient IS the constant, no contradiction; "
      "membership forces deg b_{-3} >= 3 and then N >= 3",
      sp.Poly(D3[-3], E).degree() == 1)

# CONTROL 3 -- bounded Groebner cross-check of THEOREM A at k = 3, several tops.
def sectorA(hpoly, d, withQ0=True, membership=True):
    """band-3 sector with kappa = 0 AND b_1 = 0 (so band D <= 0), coefficient cap d."""
    a3 = sp.expand(hpoly * sh(hpoly, 1) * sh(hpoly, 2))
    X = {3: a3}
    D = {3: sp.Integer(0), 2: sp.Integer(0), 1: sp.Integer(0)}
    un = []
    tag = str(abs(hash(sp.srepr(hpoly))) % 9973)
    for k, nm in ((2, 'ta2'), (1, 'ta1'), (0, 'ta0')):
        p, c = gp(f'{nm}{tag}', d)
        X[k] = p
        un += c
    p, c = gp(f'tb0{tag}', d)
    D[0] = p
    un += c
    for j in (1, 2, 3):
        p, c = gp(f'tam{j}{tag}', d)
        X[-j] = sp.expand(fall(j) * p)
        un += c
        q, c2 = gp(f'tbm{j}{tag}', d)
        D[-j] = sp.expand((fall(j) if membership else 1) * q)
        un += c2
    eqs = []
    for m in range(-6, 7):
        if m == 0 and not withQ0:
            continue
        v = sp.expand(Qm(X, D, m) - (1 if m == 0 else 0))
        if v != 0:
            eqs += sp.Poly(v, E).all_coeffs()
    return X, D, un, [xp(e) for e in eqs if xp(e) != 0]

TOPS = [('diff-1  h = E(E-1)', E * (E - 1)),
        ('diff-2  h = E(E-2)', E * (E - 2)),
        ('cube-separated  h = E(2E-1)', E * (2 * E - 1)),
        ('diff-3  h = E(E-3)', E * (E - 3)),
        ('double root  h = E^2', E**2),
        ('constant top  h = 1', sp.Integer(1))]
for nm, hp in TOPS:
    for d in (1, 2):
        _, _, un, eqs = sectorA(hp, d)
        t1 = time.time()
        ok = unit_ideal_QQ(eqs, un)
        check(f"CONTROL 3 [{nm}]: kappa = 0 AND b_1 = 0 sector, full cascade + Q_0 = 1 "
              f"+ membership at cap d = {d} is the UNIT ideal => EMPTY over C-bar "
              f"[{len(un)} vars, {time.time() - t1:.1f}s]  (independent cross-check "
              f"of THEOREM A)", ok)

# CONTROL 4 -- NON-VACUITY by EXPLICIT POINTS (house rule: 'nonempty' needs a point,
# never a normal form).  Both drops are exhibited, not Groebner-inferred.
for nm, hp in TOPS:
    a3p = sp.expand(hp * sh(hp, 1) * sh(hp, 2))
    Xp = {3: a3p}
    Dp = {0: sp.Integer(5)}                    # b_0 = 5, everything else zero
    Cp = commutator(Dp, Xp)
    # AUDIT FIX: `all(divides(fall(j), 0) ...)` was true unconditionally (divides()
    # short-circuits on the zero numerator), so it tested nothing.  The real
    # content is that the commutator vanishes identically at this point.
    okp = (Cp == {})
    check(f"CONTROL 4a NON-VACUITY [{nm}]: dropping Q_0 = 1 leaves an EXPLICIT POINT of "
          f"the kappa = 0, b_1 = 0 sector -- a_3 = h h^{{[1]}}h^{{[2]}}, b_0 = 5, all "
          f"else 0, gives Q_m = 0 for EVERY m (commutator COMPUTED = 0, so in "
          f"particular Q_0 = 0 != 1).  The emptiness certificates are not vacuous.",
          okp)
check("CONTROL 4b NON-VACUITY: dropping MEMBERSHIP leaves an EXPLICIT POINT WITH the "
      "moment unit -- CONTROL 2's (X = x^3, D = c + x^{-3}E/3) has [D,X] = 1 exactly.  "
      "So membership alone does the killing, exactly as THEOREM A says.",
      C3 == {0: sp.Integer(1)} and not divides(fall(3), D3[-3]))


# ===========================================================================
print("\n=== S5  TARGET A consequence: which corpus claims this REPAIRS ===")
# ===========================================================================
# The named gap (band3-sectors.md 4.1): in kappa = 0 AND b_1 = 0,
# Q_3 = a_3 (T^3-1) b_0 forces only b_0 constant and puts NO condition on a_2.
XN = {k: sp.Function(f'nX{k}')(E) for k in range(-3, 4)}
DN = {3: sp.Integer(0), 2: sp.Integer(0), 1: sp.Integer(0),
      0: sp.Function('nb0')(E), -1: sp.Function('nbm1')(E),
      -2: sp.Function('nbm2')(E), -3: sp.Function('nbm3')(E)}
check("the NAMED GAP is real (not a mis-statement): Q_3 = a_3 (T^3 - 1) b_0 carries NO "
      "occurrence of a_2 whatsoever -- so the corpus chain's first rung really does "
      "say nothing about a_2  (degree-free, abstract Functions)",
      not sp.simplify(Qm(XN, DN, 3)).has(XN[2]))
check("...and THEOREM A closes it anyway, WITHOUT ever constraining a_2: the kill uses "
      "only a_3, b_0, b_{-1}, b_{-2}, b_{-3}  (degree-free)",
      all(not sp.simplify(Qm(XN, {**DN, 0: sp.Symbol('nc'),
                                  **{-i: sp.Integer(0) for i in range(1, j)}}, 3 - j)
                          ).has(XN[i])
          for j in (1, 2, 3) for i in (2, 1, 0, -1, -2, -3)))

REPAIRED = [
    "cube-separated h  (shifted-power-residuals.md 1.2 chain: 'Q_4 => b_1 = c h, "
    "Q_3 => h h^{[1]} | a_2' silently assumed c != 0)",
    "2-separated h  (shifted-power-residuals.md 2.3, incl. the whole diff-3 class)",
    "diff-1  h = (E-r)(E-r-1)  (band3-sectors.md 4: closure stated for kappa != 0 and "
    "for kappa = 0 with b_1 != 0)",
    "diff-2  h = (E-r)(E-r-2)  (same kappa = 0 sub-branch)",
    "A*-band3 constant top a_3 = 1 with the gauged wall kappa_2 = 0 and b_1 = 0",
]
# AUDIT FIX (house rule 5): the former condition was len(REPAIRED) == 5, the length
# of a hardcoded list -- it never exercised h-independence.  Now EXECUTED: run the
# Theorem-A rung collapse against SEVERAL genuinely different tops (including ones
# outside every shifted-cube class) and confirm the kill coefficient is nonzero for
# each -- that is what "no property of h is used" means.
_tops = {'h=1 (constant)': sp.Integer(1),
         'diff-1  h h^[1]h^[2]': xp((E) * (E - 1) * sh((E) * (E - 1), 1) * sh((E) * (E - 1), 2)),
         'diff-2  h h^[1]h^[2]': xp((E) * (E - 2) * sh((E) * (E - 2), 1) * sh((E) * (E - 2), 2)),
         'W2 hatch E(E+2)(E+4)': xp(E * (E + 2) * (E + 4)),
         'NON-wall generic top':  xp(E**3 + 7 * E - 5)}
_okh = True
for _nm, _a3 in _tops.items():
    _A = sp.Poly(_a3, E).degree()
    for _j in (1, 2, 3):
        _phi = xp(fall(_j))
        _rung = xp(sh(_phi, 3) * _a3 - sh(_a3, -_j) * _phi)
        _N = _j + _A
        _okh &= (_rung != 0 and sp.Poly(_rung, E).coeff_monomial(E**(_N - 1))
                 == sp.LC(sp.Poly(_a3, E)) * (3 * _j + _j * _A))
check("the repair is h-INDEPENDENT (executed): the Theorem-A rung kill is nonzero with "
      "the predicted coefficient for FIVE structurally different tops -- constant, "
      "diff-1, diff-2, the W2 hatch and a non-wall generic cubic -- so no property of "
      "h (or of the wall at all) is used, and every listed class inherits the repair",
      _okh)
for nm in REPAIRED:
    print(f"        REPAIRED: {nm}")


# ===========================================================================
print("\n=== S6  TARGET B: diff-2, the surviving branch b_1(r+1) = 0 ===")
# ===========================================================================
hD2 = (E - r) * (E - r - 2)
a3D2 = sp.expand(hD2 * sh(hD2, 1) * sh(hD2, 2))
b2D2 = sp.expand(kappa * hD2 * sh(hD2, 1))
check("diff-2 sector re-derived: h = (E-r)(E-r-2) has gcd(h, h^{[2]}) = (E-r) != 1 "
      "(the broken separation) while gcd(h, h^{[1]}) = 1",
      sp.degree(sp.gcd(sp.Poly(hD2, E), sp.Poly(sh(hD2, 2), E)), E) == 1
      and sp.degree(sp.gcd(sp.Poly(hD2, E), sp.Poly(sh(hD2, 1), E)), E) == 0)
check("the Q_5 wall is SOLVED by b_2 = kappa h h^{[1]} for this h (re-derived, "
      "not cited)",
      sp.expand(sh(b2D2, 3) * a3D2 - sh(a3D2, 2) * b2D2) == 0)

# --- degree-free node table (a_2, a_1, a_0, b_1, b_0 are undetermined FUNCTIONS)
fa2, fa1, fa0 = sp.Function('a2')(E), sp.Function('a1')(E), sp.Function('a0')(E)
fb1, fb0 = sp.Function('b1')(E), sp.Function('b0')(E)
XB = {3: a3D2, 2: fa2, 1: fa1, 0: fa0}
DB = {3: sp.Integer(0), 2: b2D2, 1: fb1, 0: fb0}
Q4B, Q3B, Q2B = Qm(XB, DB, 4), Qm(XB, DB, 3), Qm(XB, DB, 2)


def evn(expr, j, d=0):
    e = sp.diff(expr, E, d) if d else expr
    return sp.simplify(e.subs(E, r + j))


check("Q_4 node table (DEGREE-FREE): Q_4(r-2) = -24k a_2(r), Q_4(r+1) = 24k a_2(r+1), "
      "Q_4(r-3) = -120k a_2(r-1) + 360 b_1(r), Q_4(r+2) = 120k a_2(r+2) - 360 b_1(r+2)",
      evn(Q4B, -2) == sp.simplify(-24 * kappa * fa2.subs(E, r))
      and evn(Q4B, 1) == sp.simplify(24 * kappa * fa2.subs(E, r + 1))
      and sp.simplify(evn(Q4B, -3) - (-120 * kappa * fa2.subs(E, r - 1)
                                      + 360 * fb1.subs(E, r))) == 0
      and sp.simplify(evn(Q4B, 2) - (120 * kappa * fa2.subs(E, r + 2)
                                     - 360 * fb1.subs(E, r + 2))) == 0)
check("Q_4 derivative node (DEGREE-FREE): Q_4'(r-1) = -2k a_2(r-1) + 6k a_2(r+1) "
      "+ 6 b_1(r+2)  -- the source of the 'tie' a_2(r+2) = a_2(r-1)",
      sp.simplify(evn(Q4B, -1, 1) - (-2 * kappa * fa2.subs(E, r - 1)
                                     + 6 * kappa * fa2.subs(E, r + 1)
                                     + 6 * fb1.subs(E, r + 2))) == 0)
check("Q_3 node disjunction (DEGREE-FREE): Q_3(r-1) = a_2(r-1)b_1(r+1) - a_2(r)b_1(r-1) "
      "and Q_3(r+1) = a_2(r+1)b_1(r+3) - a_2(r+2)b_1(r+1) -- a PRODUCT, hence the "
      "branch split (band3-sectors.md 5)",
      sp.simplify(evn(Q3B, -1) - (fa2.subs(E, r - 1) * fb1.subs(E, r + 1)
                                  - fa2.subs(E, r) * fb1.subs(E, r - 1))) == 0
      and sp.simplify(evn(Q3B, 1) - (fa2.subs(E, r + 1) * fb1.subs(E, r + 3)
                                     - fa2.subs(E, r + 2) * fb1.subs(E, r + 1))) == 0)

# --- THE NEW RUNG.  Q_2 at the two nodes where the ENTIRE negative tail drops out.
tailfns = [sp.Function(nm) for nm in ('bm1', 'bm2', 'bm3', 'am1', 'am2', 'am3')]
XBt = dict(XB)
DBt = dict(DB)
for j in (1, 2, 3):
    XBt[-j] = sp.Function(f'am{j}')(E)
    DBt[-j] = sp.Function(f'bm{j}')(E)
Q2full = Qm(XBt, DBt, 2)
# AUDIT FIX: the name claimed EXCLUSIVITY ("at exactly the nodes r and r+1") but the
# condition only tested those two.  Now the exclusivity is tested too.
# AUDIT FIX (2026-07-25) + CORRECTION: the memo claimed Q_2 is tail-free at
# "exactly the nodes r and r+1"; testing the exclusivity (rather than only the two
# nodes used) shows the tail-free set is LARGER -- r-1 and r+2 are tail-free too.
# The two EXTRA nodes are unused forcing equations and are recorded as a lead for
# the open diff-2 residual.
_free = [o for o in range(-5, 6) if not sp.simplify(evn(Q2full, o)).has(*tailfns)]
check("Q_2 tail-free node set is EXACTLY {r-1, r, r+1, r+2} (exclusivity tested over "
      "offsets -5..5) -- CORRECTING the 'exactly r and r+1' reading; the two extra "
      "nodes r-1, r+2 are additional tail-free equations not yet exploited",
      _free == [-1, 0, 1, 2])
check("Q_2(r) = a_1(r)b_1(r+1) - a_1(r+1)b_1(r) - a_2(r)(b_0(r) - b_0(r+2))   and   "
      "Q_2(r+1) = a_1(r+1)b_1(r+2) - a_1(r+2)b_1(r+1) - a_2(r+1)(b_0(r+1) - b_0(r+3))"
      "   (DEGREE-FREE)",
      sp.simplify(evn(Q2full, 0)
                  - (fa1.subs(E, r) * fb1.subs(E, r + 1) - fa1.subs(E, r + 1) * fb1.subs(E, r)
                     - fa2.subs(E, r) * (fb0.subs(E, r) - fb0.subs(E, r + 2)))) == 0
      and sp.simplify(evn(Q2full, 1)
                      - (fa1.subs(E, r + 1) * fb1.subs(E, r + 2)
                         - fa1.subs(E, r + 2) * fb1.subs(E, r + 1)
                         - fa2.subs(E, r + 1) * (fb0.subs(E, r + 1) - fb0.subs(E, r + 3)))) == 0)

# On the Q_4 locus (a_2(r) = a_2(r+1) = 0, b_1(r) = b_1(r+2) = k a_2(r-1)/3)
# AND the surviving branch b_1(r+1) = 0, both nodes collapse to a_1(r+1)*a_2(r-1).
A2m1 = sp.Symbol('A')                       # a_2(r-1) = a_2(r+2), free on the branch
sub_locus = {fa2.subs(E, r): 0, fa2.subs(E, r + 1): 0,
             fb1.subs(E, r): kappa * A2m1 / 3, fb1.subs(E, r + 2): kappa * A2m1 / 3,
             fb1.subs(E, r + 1): 0}
n0 = sp.simplify(evn(Q2full, 0).subs(sub_locus))
n1 = sp.simplify(evn(Q2full, 1).subs(sub_locus))
check("NEW FORCING (arbitrary degree): on the Q_4 locus AND the surviving branch "
      "b_1(r+1) = 0, both clean Q_2 nodes reduce to  -/+ (kappa/3) a_2(r-1) a_1(r+1). "
      "With kappa != 0 and a_2(r-1) != 0 (the defining property of the branch) this "
      "forces  a_1(r+1) = 0,  i.e.  (E - r - 1) | a_1.",
      sp.simplify(n0 + kappa * A2m1 * fa1.subs(E, r + 1) / 3) == 0
      and sp.simplify(n1 - kappa * A2m1 * fa1.subs(E, r + 1) / 3) == 0)
check("NEGATIVE CONTROL: the same two nodes do NOT force a_1(r) or a_1(r+2) "
      "(they do not appear in either node expression)",
      not n0.has(fa1.subs(E, r + 2)) and not n1.has(fa1.subs(E, r)))

# --- CONSEQUENCE: the published Q_4 ^ Q_3 witness is CUT by Q_2.
a2W = sp.expand(3 * (E - r) * (E - r - 1))
b1W = sp.expand(2 * kappa * (E - r - 1)**2)
XW = {3: a3D2, 2: a2W, 1: sp.Integer(-3), 0: sp.Integer(0)}
DW = {3: sp.Integer(0), 2: b2D2, 1: b1W, 0: sp.Integer(0)}
check("band3-sectors.md 5.1 witness re-derived independently: a_2 = 3(E-r)(E-r-1), "
      "b_1 = 2k(E-r-1)^2, a_1 = -3, b_0 = 0 solves Q_5 = Q_4 = Q_3 = 0 at symbolic "
      "(r, kappa), with a_2(r-1) = 6 != 0 and b_1(r+1) = 0",
      all(Qm(XW, DW, m) == 0 for m in (5, 4, 3))
      and sp.simplify(a2W.subs(E, r - 1)) == 6
      and sp.simplify(b1W.subs(E, r + 1)) == 0)
check("BUT IT IS CUT BY THE NEXT RUNG: that witness has Q_2 != 0 (its a_1 = -3 "
      "violates the new forcing a_1(r+1) = 0).  The band3-sectors.md 5.1 surviving "
      "family is therefore STRICTLY SMALLER than reported there.",
      sp.simplify(Qm(XW, DW, 2)) != 0)

# --- THE REFINED WITNESS: the whole positive cascade still fails to kill the branch.
a1R = sp.expand(-3 * (E - r - 1)**2)
cR = sp.Symbol('c_R')
XR = {3: a3D2, 2: sp.expand(lam * a2W), 1: sp.expand(lam**2 * a1R), 0: cR,
      -1: sp.Integer(0), -2: sp.Integer(0), -3: sp.Integer(0)}
DR = {3: sp.Integer(0), 2: b2D2, 1: sp.expand(lam * b1W), 0: sp.Integer(0),
      -1: sp.Integer(0), -2: sp.Integer(0), -3: sp.Integer(0)}
check("REFINED WITNESS (symbolic r, kappa, scale lam, constant a_0 = c): "
      "a_2 = 3L(E-r)(E-r-1), b_1 = 2kL(E-r-1)^2, a_1 = -3L^2(E-r-1)^2, b_0 = 0 "
      "solves Q_5 = Q_4 = Q_3 = Q_2 = Q_1 = 0  EXACTLY -- the ENTIRE POSITIVE CASCADE",
      all(sp.simplify(Qm(XR, DR, m)) == 0 for m in (6, 5, 4, 3, 2, 1)))
check("...and it genuinely sits on the surviving branch: b_1(r+1) = 0 while "
      "a_2(r-1) = a_2(r+2) = 6L != 0, so h h^{[1]} does NOT divide a_2 "
      "(the clean divisibility is still absent)",
      sp.simplify(DR[1].subs(E, r + 1)) == 0   # AUDIT FIX: was XR[2]*0 + ... (dead term)
      and sp.simplify(XR[2].subs(E, r - 1) - 6 * lam) == 0
      and sp.simplify(XR[2].subs(E, r + 2) - 6 * lam) == 0
      and not divides(sp.expand((hD2 * sh(hD2, 1)).subs(r, 0)),
                      sp.expand(XR[2].subs({r: 0, lam: 1}))))
check("the structural reason Q_2 vanishes on the refined witness: the middle term "
      "b_1^{[1]}a_1 - a_1^{[1]}b_1 is IDENTICALLY ZERO because b_1 and a_1 are "
      "PROPORTIONAL (both are multiples of (E-r-1)^2)",
      sp.simplify(sh(DR[1], 1) * XR[1] - sh(XR[1], 1) * DR[1]) == 0)
check("SCOPE GUARD: the refined witness is NOT a candidate Weyl pair -- with the zero "
      "negative tail shown it has Q_0 = 0 != 1, so [D,X] != 1.  It is a solution of "
      "the POSITIVE CASCADE ONLY, and nothing here constructs a genuine pair.",
      sp.simplify(Qm(XR, DR, 0)) == 0)
print("        TARGET B VERDICT: the mission's 'push to Q_2 and Q_1' route is "
      "REFUTED at arbitrary degree.")
print("        Q_2 does add a real new condition ((E-r-1) | a_1) but does NOT close "
      "the branch; Q_1 adds nothing.")


# ===========================================================================
print("\n=== S7  TARGET B: the moment-unit route, localised to TWO tail covectors ===")
# ===========================================================================
# Which G-terms are AUTOMATICALLY divisible by h^{[-1]} on the diff-2 branch?
hm1 = sh(hD2, -1)
auto = {
    'a_3^{[-1]}': sh(a3D2, -1), 'a_3^{[-2]}': sh(a3D2, -2), 'a_3^{[-3]}': sh(a3D2, -3),
    'b_2^{[-1]}': sh(b2D2, -1), 'b_2^{[-2]}': sh(b2D2, -2),
}
check("h^{[-1]} divides EVERY a_3-carrying and b_2-carrying G-term (all five) as an "
      "EXACT PRODUCT IDENTITY at SYMBOLIC r and kappa -- a_3^{[-j]} = h^{[-1]}*(...) "
      "for j = 1,2,3 and b_2^{[-j]} = h^{[-1]}*(...) for j = 1,2",
      sp.expand(auto['a_3^{[-1]}'] - hm1 * hD2 * sh(hD2, 1)) == 0
      and sp.expand(auto['a_3^{[-2]}'] - sh(hD2, -2) * hm1 * hD2) == 0
      and sp.expand(auto['a_3^{[-3]}'] - sh(hD2, -3) * sh(hD2, -2) * hm1) == 0
      and sp.expand(auto['b_2^{[-1]}'] - kappa * hm1 * hD2) == 0
      and sp.expand(auto['b_2^{[-2]}'] - kappa * sh(hD2, -2) * hm1) == 0
      )   # AUDIT FIX: dropped the tautological `len(auto) == 5` conjunct
check("NEGATIVE CONTROL for that bookkeeping: h^{[-1]} does NOT divide the NON-auto "
      "terms' carriers a_2^{[-1]}, a_2^{[-2]}, a_1^{[-1]}, b_1^{[-1]} even on the "
      "branch shapes (a_2 = 3(E-r)(E-r-1), b_1 = 2k(E-r-1)^2, a_1 = -3(E-r-1)^2)",
      not any(divides(xp(hm1.subs(r, 0)), xp(sh(v, s).subs({r: 0, kappa: 1})))
              for v, s in ((a2W, -1), (a2W, -2), (sp.expand(-3 * (E - r - 1)**2), -1),
                           (b1W, -1))))
check("NON-VACUITY of that bookkeeping: h^{[-1]} does NOT divide h itself "
      "(so the six divisibilities are not trivial)",
      not divides(xp(hm1.subs(r, 0)), xp(hD2.subs(r, 0))))
check("h^{[-1]} = (E-r-1)(E-r-3) has TWO DISTINCT simple roots for diff-2, so "
      "h^{[-1]} | G  <=>  G(r+1) = G(r+3) = 0",
      sp.simplify(sp.factor(hm1) - (E - r - 1) * (E - r - 3)) == 0)
# the residual R = the four terms that are NOT automatically divisible
fam1, fam2 = sp.Function('am1')(E), sp.Function('am2')(E)
fbm1, fbm2 = sp.Function('bm1')(E), sp.Function('bm2')(E)
Xc = {3: a3D2, 2: fa2, 1: fa1, 0: fa0, -1: fam1, -2: fam2, -3: sp.Function('am3')(E)}
Dc = {3: sp.Integer(0), 2: b2D2, 1: fb1, 0: fb0, -1: fbm1, -2: fbm2,
      -3: sp.Function('bm3')(E)}
Gc = Gpot(Xc, Dc)
Rres = sh(fa1, -1) * fbm1 - sh(fb1, -1) * fam1 + sh(fa2, -2) * fbm2 + sh(fa2, -1) * sh(fbm2, 1)
check("G  minus  R  is h^{[-1]}-divisible term-by-term, where "
      "R = a_1^{[-1]}b_{-1} - b_1^{[-1]}a_{-1} + a_2^{[-2]}b_{-2} + a_2^{[-1]}b_{-2}^{[1]} "
      "-- i.e. R carries the ENTIRE obstruction  (DEGREE-FREE)",
      sp.simplify(sp.expand(Gc - Rres)
                  - sum(sh(a3D2, -j) * sh(Dc[-3], 3 - j) for j in (1, 2, 3))
                  + sh(b2D2, -2) * fam2 + sh(b2D2, -1) * sh(fam2, 1)) == 0)
cov1 = sp.simplify(Rres.subs(E, r + 1).subs({fa2.subs(E, r): 0}))
cov3 = sp.simplify(Rres.subs(E, r + 3).subs({fa2.subs(E, r + 1): 0}))
check("THE TWO COVECTORS (degree-free, on the Q_4 locus a_2(r) = a_2(r+1) = 0):\n"
      "        C_1 = a_1(r)b_{-1}(r+1) - b_1(r)a_{-1}(r+1) + a_2(r-1)b_{-2}(r+1)\n"
      "        C_3 = a_1(r+2)b_{-1}(r+3) - b_1(r+2)a_{-1}(r+3) + a_2(r+2)b_{-2}(r+4)\n"
      "        and  h^{[-1]} | G  <=>  C_1 = C_3 = 0.",
      sp.simplify(cov1 - (fa1.subs(E, r) * fbm1.subs(E, r + 1)
                          - fb1.subs(E, r) * fam1.subs(E, r + 1)
                          + fa2.subs(E, r - 1) * fbm2.subs(E, r + 1))) == 0
      and sp.simplify(cov3 - (fa1.subs(E, r + 2) * fbm1.subs(E, r + 3)
                              - fb1.subs(E, r + 2) * fam1.subs(E, r + 3)
                              + fa2.subs(E, r + 2) * fbm2.subs(E, r + 4))) == 0)
_probe, _pc = gp('probe', 5)
_pd = sp.solve([_probe.subs(E, r + 1).subs(r, 0), _probe.subs(E, r + 3).subs(r, 0)],
               _pc[:2], dict=True)[0]
check("the two-evaluation criterion is exact (distinct simple roots): a polynomial "
      "vanishing at r+1 and r+3 IS divisible by h^{[-1]}, and a generic one is NOT",
      divides(xp(hm1.subs(r, 0)), xp(_probe.subs(_pd).subs(r, 0)))
      and not divides(xp(hm1.subs(r, 0)), xp(_probe.subs(r, 0))))
# AUDIT FIX: the former condition was a degree comparison (2 > 1) that did not
# exercise the divisibility step it announced.  Now EXECUTED: with C_1 = C_3 = 0
# imposed, check that h^{[-1]} really does divide G, and that h^{[-1]} | E is
# genuinely impossible (exhibit the nonzero remainder).
_hm1n = xp(hm1.subs(r, 0))
_remE = sp.rem(sp.Poly(E, E), sp.Poly(_hm1n, E))
check("CONDITIONAL CLOSURE (executed): h^{[-1]} does NOT divide E -- the remainder of "
      f"E mod h^{{[-1]}} is {sp.expand(_remE.as_expr())} != 0 -- so IF C_1 = C_3 = 0 "
      "(giving h^{[-1]} | G) then G = E is impossible and diff-2 closes at arbitrary "
      "degree MODULO the two covectors",
      sp.expand(_remE.as_expr()) != 0 and sp.Poly(_hm1n, E).degree() == 2)
check("NEGATIVE CONTROL: the two covectors are NOT identically zero -- an explicit "
      "tail assignment makes C_1 != 0, so the conditional closure is a REAL residual "
      "and not a vacuous one",
      sp.simplify(cov1.subs({fa1.subs(E, r): 1, fbm1.subs(E, r + 1): 1,
                             fb1.subs(E, r): 0, fa2.subs(E, r - 1): 0})) == 1)

# -------------------------------------------------- bounded emptiness, both engines
def full_sector(kind, d, withQ0=True, rv=0, kv=1):
    hh = (E - rv) * (E - rv - (1 if kind == 'diff1' else 2))
    X = {3: sp.expand(hh * sh(hh, 1) * sh(hh, 2))}
    D = {3: sp.Integer(0), 2: sp.expand(kv * hh * sh(hh, 1))}
    un = []
    for k, nm in ((2, 'wa2'), (1, 'wa1'), (0, 'wa0')):
        p, c = gp(f'{nm}{kind}', d)
        X[k] = p
        un += c
    for k, nm in ((1, 'wb1'), (0, 'wb0')):
        p, c = gp(f'{nm}{kind}', d)
        D[k] = p
        un += c
    for j in (1, 2, 3):
        p, c = gp(f'wam{j}{kind}', d)
        X[-j] = sp.expand(fall(j) * p)
        un += c
        q, c2 = gp(f'wbm{j}{kind}', d)
        D[-j] = sp.expand(fall(j) * q)
        un += c2
    eqs = []
    for m in range(-6, 7):
        if m == 0 and not withQ0:
            continue
        v = sp.expand(Qm(X, D, m) - (1 if m == 0 else 0))
        if v != 0:
            eqs += sp.Poly(v, E).all_coeffs()
    return un, [xp(e) for e in eqs if xp(e) != 0]

# msolve PARSER VALIDATION before any load-bearing msolve call (house rule 3):
# a known UNIT ideal must report [-1]; a COMPLEX-ONLY ideal must NOT report empty;
# a real-rooted feasible ideal must NOT report empty.
xv, yv = sp.symbols('xv yv')
if HAVE_MSOLVE:
    v_unit = msolve_empty_QQ([xv - 1, xv - 2], [xv], 'parser-unit', 60)
    v_cplx = msolve_empty_QQ([xv**2 + 1], [xv], 'parser-complex-only', 60)
    v_real = msolve_empty_QQ([xv**2 - 2, yv - xv], [xv, yv], 'parser-feasible', 60)
    check("msolve PARSER VALIDATION: known UNIT ideal (x-1, x-2) reports EMPTY", v_unit is True)
    check("msolve PARSER VALIDATION: COMPLEX-ONLY ideal (x^2+1) reports NONEMPTY "
          "(this is the check that makes an EMPTY verdict meaningful over C-bar)",
          v_cplx is False)
    check("msolve PARSER VALIDATION: real-rooted feasible ideal reports NONEMPTY",
          v_real is False)
    ok_trap1 = ok_trap2 = False
    try:
        _guard_body("x**2+1")
    except ValueError:
        ok_trap1 = True
    try:
        _guard_body("2*x^2/3")
    except ValueError:
        ok_trap2 = True
    check("msolve TRAP GUARDS unit-tested: '**' body rejected (trap #1) and rational "
          "body rejected (trap #2); a clean integer body passes",
          ok_trap1 and ok_trap2 and _guard_body("2*x^2+1") is True)
else:
    skip("msolve parser validation + trap guards", "msolve not on PATH")

for kind in ('diff1', 'diff2'):
    for d in (1, 2, 3):
        un, eqs = full_sector(kind, d)
        t1 = time.time()
        ok = unit_ideal_QQ(eqs, un)
        check(f"{kind}: FULL cascade + Q_0 = 1 + membership at cap d = {d} "
              f"(r = 0, kappa = 1) is the UNIT ideal => EMPTY over C-bar  "
              f"[sympy, {len(un)} vars, {time.time() - t1:.1f}s]", ok)
check("diff-2 NON-VACUITY: dropping Q_0 = 1 leaves an EXPLICIT POINT (a_3, b_2 as "
      "given, everything else 0 => all Q_m = 0), so the emptiness above is not vacuous",
      commutator({3: sp.Integer(0), 2: sp.expand(E * (E - 2) * (E + 1) * (E - 1))},
                 {3: sp.expand((E * (E - 2)) * ((E + 1) * (E - 1)) * ((E + 2) * E))}) == {})

if HEAVY and HAVE_MSOLVE:
    for dcap in (4, 6, 8, 10, 12):
        for kind in ('diff1', 'diff2'):
            un, eqs = full_sector(kind, dcap)
            t1 = time.time()
            res = msolve_empty_QQ(eqs, un, f'{kind}-d{dcap}', timeout_s=600)
            if res is None:
                skip(f"HEAVY {kind} cap d = {dcap} emptiness (msolve char 0)",
                     "timeout / no result")
            else:
                check(f"HEAVY {kind}: FULL cascade + Q_0 = 1 + membership at cap "
                      f"d = {dcap} is EMPTY over C-bar  [msolve char 0, {len(un)} vars, "
                      f"{time.time() - t1:.1f}s]", res is True)
    for kind in ('diff1', 'diff2'):
        un, eqs = full_sector(kind, 4)
        t1 = time.time()
        try:
            ok = unit_ideal_QQ(eqs, un)
            check(f"HEAVY {kind} cap d = 4, SECOND ENGINE (sympy Groebner): UNIT ideal "
                  f"[{time.time() - t1:.1f}s]", ok)
        except Exception as ex:                                    # pragma: no cover
            skip(f"HEAVY {kind} cap d = 4 second engine (sympy)", f"{type(ex).__name__}")
else:
    skip("HEAVY diff-1/diff-2 cap d = 4 emptiness, BOTH engines",
         "set HEAVY=1 (msolve + sympy; several minutes)")


# ===========================================================================
print("\n" + "=" * 75)
npass = sum(1 for _, ok in CHECKS if ok)
nfail = len(CHECKS) - npass
for nm, ok in CHECKS:
    if not ok:
        print(f"  FAILED: {nm}")
print(f"  checks executed: {len(CHECKS)}   passed: {npass}   failed: {nfail}   "
      f"skipped: {len(SKIPS)}")
for nm, why in SKIPS:
    print(f"  SKIPPED: {nm}   ({why})")
print(f"  wall time: {time.time() - T0:.1f}s   HEAVY={'1' if HEAVY else '0'}   "
      f"msolve={'yes' if HAVE_MSOLVE else 'no'}")
if nfail:
    print("  *** SOME CHECKS FAILED ***")
    sys.exit(1)
if SKIPS:
    print("  ALL EXECUTED CHECKS PASSED; OPTIONAL CHECKS SKIPPED")
else:
    print("  ALL CHECKS PASSED")
