#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_band3_sectors.py  --  band-3 residual sectors, closed radical-correctly.

INDEPENDENTLY DERIVED - EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES - NOT PEER REVIEWED

Companion memo: band3-sectors.md

Conventions (frozen from the corpus, re-derived in S0):
    A_1[x^{-1}] = (+)_k x^k C[E],  E = x*del,
    (x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E),      f^{[n]}(E) = f(E+n),  T f = f^{[1]}
    Q_m = sum_{k+l=m} [ b_l^{[k]} a_k  -  a_k^{[l]} b_l ],      [D,X] = 1  <=>  Q_m = delta_{m0}
    membership  (E)_j = E(E-1)...(E-j+1)  |  a_{-j}, b_{-j}
    gauge b_3 = 0,   Q_0 = (T-1)G,   G(0) = 0 under membership.

Sector (band-3 shifted cube):  a_3 = h h^{[1]} h^{[2]},  b_3 = 0,  b_2 = kappa h h^{[1]}.
    diff-1:  h = (E-r)(E-r-1)      gcd(h,h^{[1]}) = (E-r)
    diff-2:  h = (E-r)(E-r-2)      gcd(h,h^{[2]}) = (E-r)

WHAT IS PROVED HERE, AND AT WHAT TIER
-------------------------------------
S4/S5  diff-1  RESTORATION IS FORCED, and the class CLOSES at ARBITRARY DEGREE
       (kappa != 0, and kappa = 0 with b_1 != 0).  The certificate is a genuine
       *radical* certificate with exponent 2:  on the Q_4 locus
              Q_3(r) = (kappa/2) * a_2(r-1)^2 ,
       so a_2(r-1)^2 lies in the cascade ideal.  A nonzero Groebner normal form
       for a_2(r-1) is therefore EXPECTED and proves nothing -- exactly the
       radical-vs-ideal trap.  Degree-free: every node identity is derived with
       a_2, b_1, a_1, b_0 as UNDETERMINED FUNCTIONS (sympy Function), so no
       degree cap enters anywhere in S3-S5.

S6/S7  diff-2  RESTORATION IS **NOT** FORCED.  The Q_3 node obstruction is the
       PRODUCT  a_2(r-1) * b_1(r+1)  (not a power of a_2(r-1)), so the variety
       genuinely splits.  An EXPLICIT EXACT SOLUTION at symbolic (r, kappa) with
       a_2(r-1) = 6 != 0 is exhibited and verified.

S8     TASK B  A*-band3 constant-top kappa_2 != 0 sector (integer systems only;
       msolve parser validated in-file on a known unit and a known feasible ideal).
S9     TASK C  slope-forcing probe for the constant-top sector (Rabinowitsch).

HOUSE RULES ENFORCED
    * msolve trap #1: '^' not '**'.
    * msolve trap #2: msolve MISPARSES rational coefficients -> integer input only,
      guarded, and the parser is validated in-file on a known unit ideal and a
      known feasible ideal before any load-bearing call.
    * RADICAL vs IDEAL: no geometric claim is ever inferred from a normal form.
      "forced"   <- an exponent-k membership identity (radical) or Rabinowitsch.
      "not forced" <- an EXPLICIT EXACT WITNESS POINT.
      "empty"    <- unit ideal (weak Nullstellensatz), stated with its exact cap.
      "nonempty" <- proper ideal (weak Nullstellensatz) or an explicit point.

Run:
    uv run --with sympy python research/dc1-program/verify_band3_sectors.py
    uv run --with sympy python research/dc1-program/verify_band3_sectors.py --require-msolve
    HEAVY=1 uv run --with sympy python research/dc1-program/verify_band3_sectors.py --require-msolve

Without --require-msolve, unavailable or unsuccessful invoked solver legs are explicit
optional SKIPs.  With it, every solver leg actually invoked must finish successfully and
return a recognized verdict; HEAVY-gated legs are invoked, and thus required, only when
HEAVY=1.
"""

import argparse
import os
import sys
import time
import shutil
import subprocess
import tempfile
import re
import sympy as sp

E = sp.symbols('E')
r, kappa, t = sp.symbols('r kappa t')

HEAVY = os.environ.get('HEAVY', '') not in ('', '0', 'false', 'False')
_parser = argparse.ArgumentParser(
    description="Verify the exact band-3 sector identities and bounded solver probes.",
    epilog=("Without --require-msolve, unavailable or unsuccessful invoked msolve legs "
            "are optional SKIPs. With --require-msolve, every invoked msolve leg must "
            "finish successfully and return a recognized verdict; HEAVY-gated legs are "
            "invoked, and therefore required, only when HEAVY=1."),
)
_parser.add_argument(
    '--require-msolve', action='store_true',
    help=("require msolve to be on PATH and fail closed for every solver leg actually "
          "invoked (timeout, nonzero exit, missing output, or malformed/unrecognized "
          "result); HEAVY-only legs are required only when HEAVY=1"),
)
_args = _parser.parse_args()
REQUIRE_MSOLVE = _args.require_msolve
MSOLVE_PATH = shutil.which('msolve')
HAVE_MSOLVE = MSOLVE_PATH is not None
if REQUIRE_MSOLVE and not HAVE_MSOLVE:
    _parser.error("--require-msolve requested, but msolve is not on PATH")
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


def solver_problem(tag, why):
    """Record an unavailable solver verdict, failing closed when required."""
    message = f"msolve[{tag}] {why}"
    if REQUIRE_MSOLVE:
        check(message + " (--require-msolve fails closed)", False)
    else:
        print(f"        {message}; optional leg will be skipped", flush=True)
    return None


def skip_solver_leg(name, why):
    """Skip a failed invoked leg only in ordinary optional-solver mode."""
    if not REQUIRE_MSOLVE:
        skip(name, why)


# ---------------------------------------------------------------- primitives
def sh(f, n):
    """f^{[n]}(E) = f(E+n)."""
    return sp.expand(sp.sympify(f).subs(E, E + n))


def fall(j):
    """(E)_j = E(E-1)...(E-j+1)."""
    return sp.prod([E - i for i in range(j)]) if j > 0 else sp.Integer(1)


def gp(nm, d):
    """generic polynomial of degree <= d, plus its coefficient symbols."""
    cs = [sp.Symbol(f'{nm}_{i}') for i in range(d + 1)]
    return sum(cs[i] * E**i for i in range(d + 1)), cs


def divides(fac, c):
    c = sp.expand(sp.sympify(c))
    if c == 0:
        return True
    fac = sp.expand(sp.sympify(fac))
    if fac == 0:
        return False
    return sp.rem(sp.Poly(c, E), sp.Poly(fac, E)) == 0


def gcd_deg(a, b):
    return sp.gcd(sp.Poly(sp.expand(a), E), sp.Poly(sp.expand(b), E)).degree()


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
    """G = sum_{k=1..K} sum_{j=0..k-1} ( a_k^{[j-k]} b_{-k}^{[j]} - b_k^{[j-k]} a_{-k}^{[j]} )."""
    G = 0
    for k in range(1, K + 1):
        for j in range(0, k):
            G += sh(X.get(k, 0), j - k) * sh(D.get(-k, 0), j) \
               - sh(D.get(k, 0), j - k) * sh(X.get(-k, 0), j)
    return sp.expand(G)


# ------------------------------------------------------- Groebner utilities
def unit_ideal_QQ(eqs, unk):
    """[1] <=> unit ideal <=> V = empty over C-bar (weak Nullstellensatz)."""
    gens = [sp.expand(e) for e in eqs if sp.expand(e) != 0]
    if not gens:
        return False
    G = sp.groebner(gens, *unk, order='grevlex')
    return list(G) == [sp.Integer(1)]


def rabinowitsch_forced(eqs, unk, f, tvar):
    """f vanishes on V(eqs)  <=>  (eqs, 1 - t f) is the unit ideal.  Radical-correct."""
    return unit_ideal_QQ(list(eqs) + [1 - tvar * f], list(unk) + [tvar])


# --------------------------------------- msolve driver (traps #1 and #2 armed)
def _guard_body(body):
    """The two msolve traps, as one unit-testable guard."""
    if '**' in body:                                                     # TRAP #1
        raise ValueError("msolve input must use '^' not '**'")
    if '/' in body:                                                      # TRAP #2
        raise ValueError("msolve input must have INTEGER coefficients")
    return True


def _serialize_msolve_polynomial(expr, unknowns):
    """Clear coefficient-only denominators and emit one safe msolve polynomial."""
    num, den = sp.fraction(sp.together(sp.sympify(expr)))
    if den.free_symbols & set(unknowns):
        raise ValueError("msolve input has a variable-dependent denominator")
    text = str(sp.expand(num)).replace('**', '^').replace(' ', '')
    if '/' in text or '**' in text:
        raise ValueError("msolve input contains unsafe '/' or '**' serialization")
    return text


def _serialization_rejected(expr, unknowns):
    try:
        _serialize_msolve_polynomial(expr, unknowns)
    except ValueError:
        return True
    return False


def _msolve(eqs, unk, characteristic, args, timeout_s, tag):
    vs = ",".join(str(u) for u in unk)
    body = ",\n".join(_serialize_msolve_polynomial(e, unk) for e in eqs)
    _guard_body(body)
    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix='band3-sectors-') as tmp:
        fn = os.path.join(tmp, 'in.ms')
        out = os.path.join(tmp, 'out.txt')
        with open(fn, 'w', encoding='utf-8') as fh:
            fh.write(f"{vs}\n{characteristic}\n{body}\n")
        try:
            rr = subprocess.run([MSOLVE_PATH, *args, '-f', fn, '-o', out],
                                capture_output=True, text=True, timeout=timeout_s)
        except subprocess.TimeoutExpired:
            return solver_problem(tag, f"timed out after {timeout_s}s")
        if rr.returncode != 0:
            return solver_problem(tag, f"exited with status {rr.returncode}: "
                                  f"{rr.stderr.strip()[:160]}")
        if not os.path.exists(out):
            return solver_problem(tag, "completed without creating its output file")
        with open(out, encoding='utf-8') as fh:
            txt = fh.read()
    parsed = "".join(l for l in txt.splitlines()
                     if not l.lstrip().startswith('#')).replace(' ', '').strip()
    print(f"        msolve[{tag}] char={characteristic} nv={len(unk)} ne={len(eqs)} "
          f"t={time.time() - t0:.1f}s out={parsed[:44]!r}", flush=True)
    return parsed


def _parse_empty_QQ_record(parsed):
    """Parse only proven canonical char-zero msolve verdict records.

    [-1] is the unit-ideal/EMPTY record.  The NONEMPTY form accepted here is
    the complete positive-dimensional rational-solver header emitted in the
    historical and parser-validation runs: [1,<positive decimal>,-1,[]].
    Any other output is unknown, never a mathematical verdict.
    """
    if parsed == '[-1]':
        return True
    if re.fullmatch(r'\[1,[1-9][0-9]*,-1,\[\]\]', parsed or ''):
        return False
    return None


def _parse_unit_Fp_record(parsed):
    """Accept only msolve's documented reduced-GB unit-ideal record."""
    return True if parsed == '[1]' else None


def msolve_empty_QQ(eqs, unk, tag, timeout_s=900):
    """char-0 msolve prints the complete record [-1] iff V is empty over C-bar."""
    parsed = _msolve(eqs, unk, 0, [], timeout_s, tag)
    if parsed is None:
        return None
    verdict = _parse_empty_QQ_record(parsed)
    if verdict is None:
        return solver_problem(tag, "returned an unrecognized or malformed char-0 result")
    return verdict


def msolve_unit_Fp(eqs, unk, p, tag, timeout_s=600):
    """Reduced GB [1] over F_p is corroboration only, parsed as one complete record."""
    parsed = _msolve(eqs, unk, p, ['-g', '2'], timeout_s, tag)
    if parsed is None:
        return None
    verdict = _parse_unit_Fp_record(parsed)
    if verdict is None:
        return solver_problem(tag, "returned an unrecognized or malformed F_p result")
    return verdict


# ===========================================================================
print("\n=== S0  crossed-product ladder engine, re-derived in file ===")
# ===========================================================================
Xg = {k: gp(f'ga{k}', 2)[0] for k in range(-3, 4)}
Dg = {k: gp(f'gb{k}', 2)[0] for k in range(-3, 4)}
Cg = commutator(Dg, Xg)
check("Q_m == [D,X]_m for every m in [-6,6]  (generic symbolic band-3 data)",
      all(sp.expand(Cg.get(m, 0) - Qm(Xg, Dg, m)) == 0 for m in range(-6, 7)))
check("Q_0 = (T-1)G with the staggered potential  (generic symbolic band-3 data)",
      sp.expand(sh(Gpot(Xg, Dg), 1) - Gpot(Xg, Dg) - Qm(Xg, Dg, 0)) == 0)

# negative control: the potential formula is NOT a tautology -- a wrong shift fails
check("control: G with a corrupted shift does NOT satisfy Q_0 = (T-1)G",
      sp.expand(sh(Gpot(Xg, Dg), 2) - Gpot(Xg, Dg) - Qm(Xg, Dg, 0)) != 0)

# membership => G(0) = 0, degree-free in the free quotients
Xm = {k: gp(f'ma{k}', 2)[0] for k in range(0, 4)}
Dm = {k: gp(f'mb{k}', 2)[0] for k in range(0, 4)}
for j in (1, 2, 3):
    Xm[-j] = sp.expand(fall(j) * gp(f'mu{j}', 2)[0])
    Dm[-j] = sp.expand(fall(j) * gp(f'mv{j}', 2)[0])
check("membership (E)_j | a_{-j}, b_{-j}  =>  G(0) = 0",
      sp.expand(Gpot(Xm, Dm).subs(E, 0)) == 0)
check("control: without membership G(0) is NOT identically 0",
      sp.expand(Gpot(Xg, Dg).subs(E, 0)) != 0)

# (T^n - 1) is surjective on C[E] with kernel the constants -- used repeatedly
_surj_ok, _ker_ok = True, True
for n in (1, 3):
    for m in range(0, 5):
        tgt = E**m
        psi, cps = gp('sj', m + 1)
        s = sp.solve(sp.Poly(sp.expand(sh(psi, n) - psi - tgt), E).all_coeffs(),
                     cps, dict=True)
        if not s:
            _surj_ok = False
    psi, cps = gp('sk', 4)
    s = sp.solve(sp.Poly(sp.expand(sh(psi, n) - psi), E).all_coeffs(), cps, dict=True)
    if s:
        _ker_ok = _ker_ok and all(sp.expand(psi.subs(s[0]).diff(E)) == 0 for _ in [0])
check("(T^n - 1): C[E] -> C[E] is surjective (n = 1, 3; targets E^0..E^4)", _surj_ok)
check("ker(T^n - 1) = constants (n = 1, 3)", _ker_ok)


# ===========================================================================
print("\n=== S1  the sector: Q_5 wall and the general Q_4 / Q_3 rung identities ===")
# ===========================================================================
hs = gp('hh', 2)[0]                         # fully generic quadratic top h
a3s = sp.expand(hs * sh(hs, 1) * sh(hs, 2))

# Q_5 wall: with b_3 = 0, Q_5 = (T^3 - 1) b_2 * (leading normalisation)
Xw = {3: a3s, 2: gp('wa2', 2)[0]}
Dw = {3: sp.Integer(0), 2: gp('wb2', 2)[0]}
Q5w = Qm(Xw, Dw, 5)
check("Q_5 = a_3 * (T^3 - 1) b_2 / (shape)   i.e.  Q_5 = b_2^{[3]} a_3 - a_3^{[2]} b_2",
      sp.expand(Q5w - (sh(Dw[2], 3) * a3s - sh(a3s, 2) * Dw[2])) == 0)
# the wall solution b_2 = kappa h h^{[1]} solves Q_5 = 0 for ANY h  (arbitrary degree)
Dw2 = {3: sp.Integer(0), 2: kappa * sp.expand(hs * sh(hs, 1))}
check("wall solution b_2 = kappa h h^{[1]} satisfies Q_5 = 0  (generic quadratic h)",
      sp.expand(Qm({3: a3s}, Dw2, 5)) == 0)
check("control: b_2 = kappa h (wrong wall) does NOT satisfy Q_5 = 0",
      sp.expand(Qm({3: a3s}, {3: sp.Integer(0), 2: kappa * hs}, 5)) != 0)

# general Q_4 / Q_3 identities, NO shape assumption on a_2, b_1, a_1, b_0
Xf = {3: a3s, 2: gp('fa2', 2)[0], 1: gp('fa1', 2)[0], 0: gp('fa0', 2)[0]}
Df = {3: sp.Integer(0), 2: kappa * sp.expand(hs * sh(hs, 1)),
      1: gp('fb1', 2)[0], 0: gp('fb0', 2)[0]}
Q4id = sp.expand(kappa * (sh(hs, 2) * sh(hs, 3) * Xf[2] - hs * sh(hs, 1) * sh(Xf[2], 2))
                 + sh(hs, 1) * sh(hs, 2) * (sh(Df[1], 3) * hs - sh(hs, 3) * Df[1]))
Q3id = sp.expand(a3s * (sh(Df[0], 3) - Df[0])
                 + (sh(Df[1], 2) * Xf[2] - sh(Xf[2], 1) * Df[1])
                 + kappa * sh(hs, 1) * (sh(hs, 2) * Xf[1] - hs * sh(Xf[1], 2)))
check("Q_4 = h^{[1]}h^{[2]}(b_1^{[3]}h - h^{[3]}b_1) + kappa(h^{[2]}h^{[3]}a_2 - h h^{[1]}a_2^{[2]})"
      "  (no shape assumption)", sp.expand(Qm(Xf, Df, 4) - Q4id) == 0)
check("Q_3 = a_3(b_0^{[3]}-b_0) + (b_1^{[2]}a_2 - a_2^{[1]}b_1) + kappa h^{[1]}(h^{[2]}a_1 - h a_1^{[2]})"
      "  (no shape assumption)", sp.expand(Qm(Xf, Df, 3) - Q3id) == 0)

# the clean-shape cancellation the broken classes are accused of losing
gg = gp('cg', 2)[0]
bb = gp('cb', 2)[0]
mid_clean = sp.expand(sh(hs * bb, 2) * (hs * sh(hs, 1) * gg)
                      - sh(hs * sh(hs, 1) * gg, 1) * (hs * bb))
check("clean shapes a_2 = h h^{[1]}g, b_1 = h beta  =>  a_3 | (b_1^{[2]}a_2 - a_2^{[1]}b_1)",
      divides(a3s, mid_clean))
check("clean-shape quotient equals beta^{[2]}g - g^{[1]}beta",
      sp.expand(sp.cancel(mid_clean / a3s) - (sh(bb, 2) * gg - sh(gg, 1) * bb)) == 0)


# ===========================================================================
print("\n=== S2  normalisation: the (r, kappa) reduction is an exact bijection ===")
# ===========================================================================
# The Q_m are built from shifts and products only, so E -> E + r is a ring
# automorphism commuting with T; and Q_4, Q_3 are (quasi-)homogeneous in kappa.
_rho = sp.Symbol('rho')


def sector_data(kind, rr=r, kk=kappa):
    hh = (E - rr) * (E - rr - (1 if kind == 'diff1' else 2))
    return (sp.expand(hh), sp.expand(sh(hh, 1)), sp.expand(sh(hh, 2)),
            sp.expand(sh(hh, 3)), sp.expand(hh * sh(hh, 1) * sh(hh, 2)),
            sp.expand(kk * hh * sh(hh, 1)))


for kind in ('diff1', 'diff2'):
    h0, h1_, h2_, h3_, a30, b20 = sector_data(kind, 0, 1)
    hr, hr1, hr2, hr3, a3r, b2r = sector_data(kind, _rho, 1)
    check(f"{kind}: shift covariance  h(E)|_{{r=rho}} = h(E-rho)|_{{r=0}}",
          sp.expand(hr - h0.subs(E, E - _rho)) == 0)
    check(f"{kind}: shift covariance of the whole sector (a_3, b_2)",
          sp.expand(a3r - a30.subs(E, E - _rho)) == 0
          and sp.expand(b2r - b20.subs(E, E - _rho)) == 0)

# kappa scaling: (a_2, b_1) solves Q_4 at kappa  <=>  (a_2, s b_1) solves Q_4 at s*kappa
s_ = sp.Symbol('s')
for kind in ('diff1', 'diff2'):
    hh, h1_, h2_, h3_, a3_, _ = sector_data(kind, r, kappa)
    A2f = sp.Function('a2')
    B1f = sp.Function('b1')
    Q4a = (kappa * (h2_ * h3_ * A2f(E) - hh * h1_ * A2f(E + 2))
           + a3_ * B1f(E + 3) - sh(a3_, 1) * B1f(E))
    Q4b = ((s_ * kappa) * (h2_ * h3_ * A2f(E) - hh * h1_ * A2f(E + 2))
           + a3_ * s_ * B1f(E + 3) - sh(a3_, 1) * s_ * B1f(E))
    check(f"{kind}: Q_4(s*kappa; a_2, s b_1) = s * Q_4(kappa; a_2, b_1)  (degree-free)",
          sp.expand(Q4b - s_ * Q4a) == 0)


# ===========================================================================
print("\n=== S3  DEGREE-FREE node table (a_2, b_1, a_1, b_0 = undetermined functions) ===")
# ===========================================================================
# No degree cap anywhere below: a_2, b_1, a_1, b_0 are sympy Functions, so every
# identity holds for polynomials of ARBITRARY degree (indeed for any functions).
A2 = sp.Function('a2')
B1 = sp.Function('b1')
A1 = sp.Function('a1')
B0 = sp.Function('b0')


def rungs_df(kind):
    """Q_4, Q_3 of the sector with UNDETERMINED a_2, b_1, a_1, b_0."""
    hh, h1_, h2_, h3_, a3_, b2_ = sector_data(kind)
    Q4 = sp.expand(kappa * (h2_ * h3_ * A2(E) - hh * h1_ * A2(E + 2))
                   + a3_ * B1(E + 3) - sp.expand(sh(a3_, 1)) * B1(E))
    Q3 = sp.expand(a3_ * (B0(E + 3) - B0(E))
                   + (B1(E + 2) * A2(E) - A2(E + 1) * B1(E))
                   + kappa * h1_ * (h2_ * A1(E) - hh * A1(E + 2)))
    return hh, h1_, h2_, h3_, a3_, b2_, Q4, Q3


def ev(expr, j, deriv=0):
    """value (or E-derivative) of expr at E = r + j, exactly."""
    e = sp.diff(expr, E, deriv) if deriv else expr
    return sp.simplify(sp.expand(e.subs(E, r + j)))


D1 = rungs_df('diff1')
D2 = rungs_df('diff2')

# --- the h-evaluation tables (pure arithmetic, degree-free, and falsifiable) ---
check("diff-1: h(r+j) for j=-3..2 equals [12, 6, 2, 0, 0, 2]",
      [sp.simplify(D1[0].subs(E, r + j)) for j in range(-3, 3)] == [12, 6, 2, 0, 0, 2])
check("diff-2: h(r+j) for j=-3..2 equals [15, 8, 3, 0, -1, 0]",
      [sp.simplify(D2[0].subs(E, r + j)) for j in range(-3, 3)] == [15, 8, 3, 0, -1, 0])
check("diff-1: gcd(h,h^{[1]}) = (E-r) (deg 1); gcd(h,h^{[2]}) = gcd(h,h^{[3]}) = 1",
      gcd_deg(D1[0].subs(r, 0), D1[1].subs(r, 0)) == 1
      and gcd_deg(D1[0].subs(r, 0), D1[2].subs(r, 0)) == 0
      and gcd_deg(D1[0].subs(r, 0), D1[3].subs(r, 0)) == 0)
check("diff-2: gcd(h,h^{[2]}) = (E-r) (deg 1); gcd(h,h^{[1]}) = gcd(h,h^{[3]}) = 1",
      gcd_deg(D2[0].subs(r, 0), D2[2].subs(r, 0)) == 1
      and gcd_deg(D2[0].subs(r, 0), D2[1].subs(r, 0)) == 0
      and gcd_deg(D2[0].subs(r, 0), D2[3].subs(r, 0)) == 0)

# multiple-root top (E-r)^2 is cube-separated, hence NOT a broken class
hmr = (E - r)**2
check("mult-root top h = (E-r)^2 is CUBE-separated: gcd(h,h^{[j]}) = 1 for j = 1,2,3",
      all(gcd_deg(hmr.subs(r, 0), sh(hmr, j).subs(r, 0)) == 0 for j in (1, 2, 3)))

# --- diff-1 degree-free node identities -------------------------------------
_, _, _, _, _, _, Q4_1, Q3_1 = D1
id_d1 = [
    ("Q_4(r)      = 12 kappa a_2(r)",
     ev(Q4_1, 0), 12 * kappa * A2(r)),
    ("Q_4'(r-1)   = 2 kappa ( a_2(r-1) + a_2(r+1) )",
     ev(Q4_1, -1, 1), 2 * kappa * (A2(r - 1) + A2(r + 1))),
    ("Q_4(r-3)    = 144 b_1(r) - 72 kappa a_2(r-1)",
     ev(Q4_1, -3), 144 * B1(r) - 72 * kappa * A2(r - 1)),
    ("Q_4(r+1)    = 72 kappa a_2(r+1) - 144 b_1(r+1)",
     ev(Q4_1, 1), 72 * kappa * A2(r + 1) - 144 * B1(r + 1)),
    ("Q_4'(r)     = 28 kappa a_2(r) + 12 kappa a_2'(r) - 12 b_1(r)",
     ev(Q4_1, 0, 1),
     28 * kappa * A2(r) + 12 * kappa * sp.Subs(sp.Derivative(A2(E), E), E, r) - 12 * B1(r)),
    ("Q_3(r)      = a_2(r) b_1(r+2) - a_2(r+1) b_1(r)",
     ev(Q3_1, 0), A2(r) * B1(r + 2) - A2(r + 1) * B1(r)),
    ("Q_3(r-1)    = a_2(r-1) b_1(r+1) - a_2(r) b_1(r-1)",
     ev(Q3_1, -1), A2(r - 1) * B1(r + 1) - A2(r) * B1(r - 1)),
]
for nm, lhs, rhs in id_d1:
    check("diff-1 [degree-free]  " + nm, sp.simplify(sp.expand(lhs - rhs)) == 0)

# --- diff-2 degree-free node identities --------------------------------------
_, _, _, _, _, _, Q4_2, Q3_2 = D2
id_d2 = [
    ("Q_4(r-2)    = -24 kappa a_2(r)", ev(Q4_2, -2), -24 * kappa * A2(r)),
    ("Q_4(r+1)    =  24 kappa a_2(r+1)", ev(Q4_2, 1), 24 * kappa * A2(r + 1)),
    ("Q_4(r-3)    = -120 kappa a_2(r-1) + 360 b_1(r)",
     ev(Q4_2, -3), -120 * kappa * A2(r - 1) + 360 * B1(r)),
    ("Q_4(r+2)    =  120 kappa a_2(r+2) - 360 b_1(r+2)",
     ev(Q4_2, 2), 120 * kappa * A2(r + 2) - 360 * B1(r + 2)),
    ("Q_4'(r-1)   = -2 kappa a_2(r-1) + 6 kappa a_2(r+1) + 6 b_1(r+2)",
     ev(Q4_2, -1, 1),
     -2 * kappa * A2(r - 1) + 6 * kappa * A2(r + 1) + 6 * B1(r + 2)),
    ("Q_3(r)      = a_2(r) b_1(r+2) - a_2(r+1) b_1(r)",
     ev(Q3_2, 0), A2(r) * B1(r + 2) - A2(r + 1) * B1(r)),
    ("Q_3(r-1)    = a_2(r-1) b_1(r+1) - a_2(r) b_1(r-1)",
     ev(Q3_2, -1), A2(r - 1) * B1(r + 1) - A2(r) * B1(r - 1)),
    ("Q_3(r+1)    = a_2(r+1) b_1(r+3) - a_2(r+2) b_1(r+1)",
     ev(Q3_2, 1), A2(r + 1) * B1(r + 3) - A2(r + 2) * B1(r + 1)),
]
for nm, lhs, rhs in id_d2:
    check("diff-2 [degree-free]  " + nm, sp.simplify(sp.expand(lhs - rhs)) == 0)

# control: a node identity that is FALSE must be reported false
check("control: the false identity Q_4(r) = 13 kappa a_2(r) is rejected",
      sp.simplify(sp.expand(ev(Q4_1, 0) - 13 * kappa * A2(r))) != 0)


# ===========================================================================
print("\n=== S4  TASK A / diff-1: the RADICAL certificate (exponent 2, degree-free) ===")
# ===========================================================================
# Abstract node values.  Because S3 proved the node identities with a_2, b_1,
# a_1, b_0 UNDETERMINED, the relations below are relations between the values of
# ARBITRARY polynomial data -- no degree cap is involved.
Ar0, Ar1, Arm1, Br0, Br1, Br2, Brm1, dAr0 = sp.symbols(
    'A_r A_rp1 A_rm1 B_r B_rp1 B_rp2 B_rm1 dA_r')

F1 = 12 * kappa * Ar0                          # = Q_4(r)
F2 = 2 * kappa * (Arm1 + Ar1)                  # = Q_4'(r-1)
F3 = 144 * Br0 - 72 * kappa * Arm1             # = Q_4(r-3)
F4 = 72 * kappa * Ar1 - 144 * Br1              # = Q_4(r+1)
F5 = 28 * kappa * Ar0 + 12 * kappa * dAr0 - 12 * Br0   # = Q_4'(r)
G1 = Ar0 * Br2 - Ar1 * Br0                     # = Q_3(r)
G2 = Arm1 * Br1 - Ar0 * Brm1                   # = Q_3(r-1)

# ---- THE CERTIFICATE ------------------------------------------------------
# 144 kappa^2 a_2(r-1)^2  =  288 kappa Q_3(r) - 24 b_1(r+2) Q_4(r)
#                            + Q_4'(r-1) Q_4(r-3) + 72 kappa a_2(r-1) Q_4'(r-1)
#                            - 2 kappa a_2(r-1) Q_4(r-3)
CERT_LHS = 144 * kappa**2 * Arm1**2
CERT_RHS = (288 * kappa * G1 - 24 * Br2 * F1 + F2 * F3
            + 72 * kappa * Arm1 * F2 - 2 * kappa * Arm1 * F3)
check("diff-1 CERTIFICATE: 144 kappa^2 a_2(r-1)^2 = 288 kappa Q_3(r) - 24 b_1(r+2) Q_4(r)"
      " + Q_4'(r-1)Q_4(r-3) + 72 kappa a_2(r-1)Q_4'(r-1) - 2 kappa a_2(r-1)Q_4(r-3)",
      sp.expand(CERT_LHS - CERT_RHS) == 0)
check("control: the certificate FAILS if the Q_3(r) cofactor is perturbed (289 kappa)",
      sp.expand(CERT_LHS - (CERT_RHS + kappa * G1)) != 0)

# A second, independent exponent-2 certificate from the node r-1
CERT2 = (-288 * kappa * G2 - 24 * Brm1 * F1
         + 72 * kappa * Arm1 * F2 - 2 * kappa * Arm1 * F4)
check("diff-1 SECOND certificate (node r-1): 144 kappa^2 a_2(r-1)^2 = -288 kappa Q_3(r-1)"
      " - 24 b_1(r-1)Q_4(r) + 72 kappa a_2(r-1)Q_4'(r-1) - 2 kappa a_2(r-1)Q_4(r+1)",
      sp.expand(144 * kappa**2 * Arm1**2 - CERT2) == 0)
check("control: the second certificate FAILS if Q_4(r+1) is swapped for Q_4(r-3)",
      sp.expand(144 * kappa**2 * Arm1**2
                - (CERT2 + 2 * kappa * Arm1 * F4 - 2 * kappa * Arm1 * F3)) != 0)

# The certificate is exactly a RADICAL statement: a_2(r-1)^2 in I, a_2(r-1) NOT in I.
# Machine-check the non-membership of the linear form in the ideal generated by the
# node functionals -- so that a naive normal-form test would MISLEAD.
_nodevars = [Ar0, Ar1, Arm1, Br0, Br1, Br2, Brm1, dAr0]
_Inode = [F1, F2, F3, F4, F5, G1, G2]
_kfix = {kappa: sp.Integer(1)}
_Ifix = [sp.expand(f.subs(_kfix)) for f in _Inode]
_GB = sp.groebner(_Ifix, *_nodevars, order='grevlex')
check("diff-1: a_2(r-1) has NONZERO normal form mod the node ideal"
      "  (so a normal form alone would falsely suggest 'not forced')",
      sp.expand(_GB.reduce(Arm1)[1]) != 0)
check("diff-1: a_2(r-1)^2 HAS zero normal form mod the node ideal"
      "  (radical membership, exponent 2)",
      sp.expand(_GB.reduce(Arm1**2)[1]) == 0)
# Rabinowitsch, the independent radical-correct route
check("diff-1: Rabinowitsch -- (node ideal, 1 - t a_2(r-1)) is the UNIT ideal",
      rabinowitsch_forced(_Ifix, _nodevars, Arm1, t))
check("control: Rabinowitsch on a genuinely free node value (b_1(r-1)) is NOT unit",
      not rabinowitsch_forced(_Ifix, _nodevars, Brm1, t))

# ---- consequence: the FULL clean divisibility is restored -------------------
_sol = sp.solve([F1, F2, F3, F4, F5, sp.Eq(Arm1, 0)],
                [Ar0, Ar1, Br0, Br1, dAr0], dict=True)
check("diff-1: on Q_4 = 0 with a_2(r-1) = 0  =>  a_2(r) = a_2(r+1) = a_2'(r) = 0"
      " and b_1(r) = b_1(r+1) = 0",
      bool(_sol) and all(sp.simplify(_sol[0].get(v, v).subs(Arm1, 0)) == 0
                         for v in (Ar0, Ar1, Br0, Br1, dAr0)))
# and, before Q_3, a_2'(r) is genuinely NOT forced -- the memo's claim, confirmed
_solQ4 = sp.solve([F1, F2, F3, F4, F5], [Ar0, Ar1, Br0, Br1, dAr0], dict=True)
check("diff-1: from Q_4 ALONE, a_2'(r) = a_2(r-1)/2 -- free, so h h^{[1]} | a_2 FAILS at Q_4",
      bool(_solQ4) and sp.simplify(_solQ4[0][dAr0] - Arm1 / 2) == 0)


# --- jet -> divisibility, degree-free via the confluent Vandermonde ----------
def jet_matrix(fac, funcs):
    """fac of degree n; funcs = n linear functionals on the remainder space
       C[E]_{<n}.  Invertible matrix  <=>  (fac | P  <=>  all functionals kill P),
       for polynomials P of ARBITRARY degree (P and P mod fac have the same
       functional values when each functional factors through C[E]/(fac))."""
    n = sp.Poly(sp.expand(fac), E).degree()
    basis = [E**i for i in range(n)]
    M = sp.Matrix([[f(b) for b in basis] for f in funcs])
    return M


h1s = D1[0].subs(r, 0)                    # diff-1 h at r = 0 (shift-covariant)
hh1 = sp.expand(h1s * sh(h1s, 1))         # h h^{[1]}
Mjet = jet_matrix(hh1, [lambda p: p.subs(E, -1), lambda p: p.subs(E, 0),
                        lambda p: sp.diff(p, E).subs(E, 0), lambda p: p.subs(E, 1)])
check("diff-1 [degree-free]: the jet (P(r-1), P(r), P'(r), P(r+1)) determines P mod h h^{[1]}"
      "  (confluent Vandermonde det != 0)", Mjet.det() != 0)
Mh = jet_matrix(h1s, [lambda p: p.subs(E, 0), lambda p: p.subs(E, 1)])
check("diff-1 [degree-free]: (P(r), P(r+1)) determines P mod h   (det != 0)",
      Mh.det() != 0)
# Radical-correct conclusion: EVERY jet functional of the clean divisibility is
# forced (Rabinowitsch per functional), while a genuinely free value is not.
_forced = {v: rabinowitsch_forced(_Ifix, _nodevars, v, t)
           for v in (Ar0, Ar1, Arm1, dAr0, Br0, Br1)}
check("diff-1 [degree-free] CONCLUSION: Rabinowitsch forces the WHOLE clean jet --"
      " a_2(r-1), a_2(r), a_2'(r), a_2(r+1), b_1(r), b_1(r+1) all vanish on V"
      "  =>  h h^{[1]} | a_2 and h | b_1  (arbitrary degree)",
      all(_forced.values()) and Mjet.det() != 0 and Mh.det() != 0)
check("control: b_1(r-1) and b_1(r+2) are NOT forced (the conclusion is not vacuous)",
      (not rabinowitsch_forced(_Ifix, _nodevars, Brm1, t))
      and (not rabinowitsch_forced(_Ifix, _nodevars, Br2, t)))


# ===========================================================================
print("\n=== S5  TASK A consequence / diff-1: ARBITRARY-DEGREE CLOSURE ===")
# ===========================================================================
hD1 = D1[0]                                     # h = (E-r)(E-r-1), symbolic r
hD1_1, hD1_2, hD1_3 = sh(hD1, 1), sh(hD1, 2), sh(hD1, 3)
hD1_m1 = sh(hD1, -1)
a3D1 = sp.expand(hD1 * hD1_1 * hD1_2)

# (i) with the restored shapes, Q_3 = 0 forces h | a_1  (degree-free)
Rsrc = sp.expand(kappa * hD1_1 * (hD1_2 * A1(E) - hD1 * A1(E + 2)))
check("diff-1 [degree-free]: R := kappa h^{[1]}(h^{[2]}a_1 - h a_1^{[2]}) has R(r) = 0"
      " (double node) ", sp.simplify(ev(Rsrc, 0)) == 0)
check("diff-1 [degree-free]: R'(r) = 2 kappa a_1(r)   =>  a_1(r) = 0",
      sp.simplify(sp.expand(ev(Rsrc, 0, 1) - 2 * kappa * A1(r))) == 0)
check("diff-1 [degree-free]: R(r+1) = 12 kappa a_1(r+1)   =>  a_1(r+1) = 0",
      sp.simplify(sp.expand(ev(Rsrc, 1) - 12 * kappa * A1(r + 1))) == 0)
check("diff-1 [degree-free]: a_1(r) = a_1(r+1) = 0  =>  h | a_1  (Vandermonde det != 0)",
      Mh.det() != 0)

# (ii) the potential factorisation G = h^{[-1]} M, degree-free term by term
FK = {1: hD1, 2: sp.expand(hD1 * hD1_1), 3: a3D1}
_gdiv = []
for k in (1, 2, 3):
    for j in range(0, k):
        _gdiv.append((k, j, divides(hD1_m1, sh(FK[k], j - k))))
check("diff-1 [degree-free]: EVERY G-term carries h^{[-1]}  -- h^{[-1]} | F_k^{[j-k]}"
      " for all 1<=k<=3, 0<=j<k  (F_1 = h, F_2 = h h^{[1]}, F_3 = a_3; b_3 = 0)",
      all(ok for _, _, ok in _gdiv))
check("control: h^{[-1]} does NOT divide F_1^{[0]} = h  (the bookkeeping is not vacuous)",
      not divides(hD1_m1, FK[1]))

# instance cross-check of the factorisation on generic quotients
_g2 = gp('zg', 2)[0]
_be = gp('zb', 2)[0]
_p1 = gp('zp', 2)[0]
Xc = {3: a3D1, 2: sp.expand(hD1 * hD1_1 * _g2), 1: sp.expand(hD1 * _p1),
      0: gp('za0', 2)[0]}
Dc = {3: sp.Integer(0), 2: kappa * sp.expand(hD1 * hD1_1), 1: sp.expand(hD1 * _be),
      0: gp('zb0', 2)[0]}
for k in (1, 2, 3):
    Xc[-k] = gp(f'zam{k}', 2)[0]
    Dc[-k] = gp(f'zbm{k}', 2)[0]
Gc = Gpot(Xc, Dc)
check("diff-1: G = h^{[-1]} M on the restored shapes (generic quotients, symbolic r)",
      divides(hD1_m1, Gc))

# (iii) the affine kill: Q_0 = 1 and G(0) = 0 give G = E; deg h^{[-1]} = 2 > 1
_Gsym = gp('gk', 3)[0]
_gk = sp.solve(sp.Poly(sp.expand(sh(_Gsym, 1) - _Gsym - 1), E).all_coeffs()
               + [_Gsym.subs(E, 0)], gp('gk', 3)[1], dict=True)
check("Q_0 = 1 and G(0) = 0  =>  G = E  (uses ker(T-1) = constants)",
      bool(_gk) and sp.expand(_Gsym.subs(_gk[0]) - E) == 0)
check("diff-1 CLOSURE: deg h^{[-1]} = 2 > 1 = deg E, so h^{[-1]} | E is IMPOSSIBLE"
      "  =>  the diff-1 sector is EMPTY at ARBITRARY DEGREE (kappa != 0)",
      sp.Poly(sp.expand(hD1_m1.subs(r, 0)), E).degree() == 2 and not divides(hD1_m1, E))

# (iv) the kappa = 0 branch
# kappa = 0: Q_4 = h^{[1]}h^{[2]} W with W = h b_1^{[3]} - h^{[3]} b_1, and
# h^{[1]}h^{[2]} != 0, so Q_4 = 0 <=> W = 0.
Q4k0 = sp.expand(hD1_1 * hD1_2 * (hD1 * B1(E + 3) - hD1_3 * B1(E)))
Wk0 = sp.expand(hD1 * B1(E + 3) - hD1_3 * B1(E))
check("diff-1 kappa=0 [degree-free]: Q_4 = h^{[1]}h^{[2]} W,  W = h b_1^{[3]} - h^{[3]}b_1",
      sp.expand(Q4k0 - hD1_1 * hD1_2 * Wk0) == 0)
check("diff-1 kappa=0 [degree-free]: W(r) = -6 b_1(r), W(r+1) = -12 b_1(r+1)  =>  h | b_1",
      sp.simplify(sp.expand(ev(Wk0, 0) + 6 * B1(r))) == 0
      and sp.simplify(sp.expand(ev(Wk0, 1) + 12 * B1(r + 1))) == 0)
_psi = gp('kp', 4)
_pk = sp.solve(sp.Poly(sp.expand(sh(_psi[0], 3) - _psi[0]), E).all_coeffs(),
               _psi[1], dict=True)
check("diff-1 kappa=0: b_1 = h psi with psi^{[3]} = psi  =>  psi = c constant",
      bool(_pk) and sp.expand(sp.diff(_psi[0].subs(_pk[0]), E)) == 0)
Fk0 = sp.expand(hD1_2 * A2(E) - hD1 * A2(E + 1))
check("diff-1 kappa=0, c != 0 [degree-free]: a_3 | F := h^{[2]}a_2 - h a_2^{[1]} gives"
      " F(r) = 2 a_2(r), F(r+1) = 6 a_2(r+1), F(r-2) = -6 a_2(r-1)",
      sp.simplify(sp.expand(ev(Fk0, 0) - 2 * A2(r))) == 0
      and sp.simplify(sp.expand(ev(Fk0, 1) - 6 * A2(r + 1))) == 0
      and sp.simplify(sp.expand(ev(Fk0, -2) + 6 * A2(r - 1))) == 0)
check("diff-1 kappa=0, c != 0 [degree-free]: F'(r) = 3 a_2(r) + a_2(r+1) + 2 a_2'(r),"
      " so the three vanishings above force a_2'(r) = 0  =>  h h^{[1]} | a_2",
      sp.simplify(sp.expand(ev(Fk0, 0, 1)
                            - (3 * A2(r) + A2(r + 1)
                               + 2 * sp.Derivative(A2(r), r)))) == 0)
Vk0 = sp.expand(hD1_1 * A1(E) - hD1 * A1(E + 1))
_ck0 = sp.symbols('c_k0', nonzero=True)
_gk0 = gp('gk0', 2)[0]
Xk0 = {3: a3D1, 2: sp.expand(hD1 * hD1_1 * _gk0), 1: A1(E)}
Dk0 = {3: sp.Integer(0), 2: sp.Integer(0), 1: sp.expand(_ck0 * hD1)}
Q2true = sp.expand(Qm(Xk0, Dk0, 2))
check("diff-1 kappa=0, c != 0 [degree-free]: direct crossed-product Q_2 satisfies "
      "Q_2-c V in (h h^{[1]}), V=h^{[1]}a_1-h a_1^{[1]}",
      divides(sp.expand(hD1 * hD1_1), sp.expand(Q2true - _ck0 * Vk0)))
check("diff-1 kappa=0, c != 0 [degree-free]: Q_2 gives h h^{[1]} | V := h^{[1]}a_1 - h a_1^{[1]},"
      " and V(r+1) = 2 a_1(r+1), V(r-1) = -2 a_1(r)  =>  h | a_1",
      sp.simplify(sp.expand(ev(Vk0, 1) - 2 * A1(r + 1))) == 0
      and sp.simplify(sp.expand(ev(Vk0, -1) + 2 * A1(r))) == 0)
print("      NOTE (audit 2026-07-26): the sub-branch kappa = 0 AND b_1 = 0 is NOT closed")
print("      by this direct residual/congruence route -- it stays silent there.  It is")
print("      closed INDEPENDENTLY by the Nonpositive-D Exclusion Theorem of")
print("      shifted-cube-completion.md, so the diff-1 sector is fully closed at")
print("      arbitrary degree in all three branches:")
print("        * kappa != 0                : this file (radical-correct certificate)")
print("        * kappa = 0, b_1 != 0       : this file (direct cascade)")
print("        * kappa = 0, b_1 = 0        : Nonpositive-D Exclusion Theorem")


# ===========================================================================
print("\n=== S6  TASK A / diff-2: restoration is NOT forced -- EXPLICIT EXACT WITNESS ===")
# ===========================================================================
# The diff-2 node obstruction is a PRODUCT, not a power: hence no radical
# membership, and the variety genuinely splits into two branches.
Ar2 = sp.Symbol('A_rp2')
d2F = [
    -24 * kappa * Ar0,                                   # Q_4(r-2)
    24 * kappa * Ar1,                                    # Q_4(r+1)
    -120 * kappa * Arm1 + 360 * Br0,                     # Q_4(r-3)
    120 * kappa * Ar2 - 360 * Br2,                       # Q_4(r+2)
    -2 * kappa * Arm1 + 6 * kappa * Ar1 + 6 * Br2,       # Q_4'(r-1)
]
d2G = [
    Ar0 * Br2 - Ar1 * Br0,                               # Q_3(r)
    Arm1 * Br1 - Ar0 * Brm1,                             # Q_3(r-1)
    Ar1 * sp.Symbol('B_rp3') - Ar2 * Br1,                # Q_3(r+1)
]
_d2vars = [Ar0, Ar1, Ar2, Arm1, Br0, Br1, Br2, Brm1, sp.Symbol('B_rp3')]
_d2sol = sp.solve(d2F, [Ar0, Ar1, Ar2, Br0, Br2], dict=True)
check("diff-2: Q_4 forces a_2(r) = a_2(r+1) = 0, b_1(r) = b_1(r+2) = kappa a_2(r-1)/3,"
      " and the TIE a_2(r+2) = a_2(r-1)  (so only (E-r)(E-r-1) | a_2 is forced)",
      bool(_d2sol) and sp.simplify(_d2sol[0][Ar0]) == 0
      and sp.simplify(_d2sol[0][Ar1]) == 0
      and sp.simplify(sp.expand(_d2sol[0][Ar2] - Arm1)) == 0
      and sp.simplify(sp.expand(_d2sol[0][Br0] - kappa * Arm1 / 3)) == 0
      and sp.simplify(sp.expand(_d2sol[0][Br2] - kappa * Arm1 / 3)) == 0)
check("control: a_2(r-1) itself is NOT solved away by Q_4 (it stays a free node value)",
      Arm1 not in _d2sol[0])
_d2on = {Ar0: 0, Ar1: 0, Ar2: Arm1, Br0: kappa * Arm1 / 3, Br2: kappa * Arm1 / 3}
check("diff-2: on the Q_4 locus, Q_3(r) VANISHES identically (no constraint)",
      sp.simplify(sp.expand(d2G[0].subs(_d2on))) == 0)
check("diff-2 OBSTRUCTION IS A PRODUCT: Q_3(r-1) = a_2(r-1) b_1(r+1) and"
      " Q_3(r+1) = -a_2(r-1) b_1(r+1) on the Q_4 locus",
      sp.simplify(sp.expand(d2G[1].subs(_d2on) - Arm1 * Br1)) == 0
      and sp.simplify(sp.expand(d2G[2].subs(_d2on) + Arm1 * Br1)) == 0)
# radical-correct: a_2(r-1) is NOT in the radical of the diff-2 node ideal
_d2I = [sp.expand(f.subs(kappa, 1)) for f in d2F] + \
       [sp.expand(g.subs(kappa, 1)) for g in d2G]
check("diff-2: Rabinowitsch says a_2(r-1) is NOT forced (ideal + 1 - t a_2(r-1)"
      " is NOT the unit ideal)",
      not rabinowitsch_forced(_d2I, _d2vars, Arm1, t))
check("control: on the SAME machinery, diff-1's a_2(r-1) IS forced -- the two"
      " classes are separated by the method, not by the setup",
      rabinowitsch_forced(_Ifix, _nodevars, Arm1, t))

# ---------------- THE EXPLICIT EXACT WITNESS, symbolic (r, kappa) ------------
hW, hW1, hW2, hW3, a3W, b2W = sector_data('diff2', r, kappa)
a2W = sp.expand(3 * (E - r) * (E - r - 1))
b1W = sp.expand(2 * kappa * (E - r - 1)**2)
a1W = sp.Integer(-3)
b0W = sp.Integer(0)
Q4W = sp.expand(sh(b2W, 2) * a2W - sh(a2W, 2) * b2W + sh(b1W, 3) * a3W - sh(a3W, 1) * b1W)
Q3W = sp.expand(a3W * (sh(b0W, 3) - b0W) + (sh(b1W, 2) * a2W - sh(a2W, 1) * b1W)
                + kappa * hW1 * (hW2 * a1W - hW * sh(a1W, 2)))
check("diff-2 WITNESS  a_2 = 3(E-r)(E-r-1),  b_1 = 2 kappa (E-r-1)^2,  a_1 = -3,  b_0 = 0"
      "  satisfies Q_4 = 0 EXACTLY at symbolic (r, kappa)", sp.expand(Q4W) == 0)
check("diff-2 WITNESS satisfies Q_3 = 0 EXACTLY at symbolic (r, kappa)",
      sp.expand(Q3W) == 0)
check("diff-2 WITNESS has a_2(r-1) = 6 != 0 and a_2(r+2) = 6 != 0"
      "  =>  the missing node values are NOT forced",
      sp.expand(a2W.subs(E, r - 1)) == 6 and sp.expand(a2W.subs(E, r + 2)) == 6)
check("diff-2 WITNESS: the degraded Q_4 forcing holds -- (E-r)(E-r-1) | a_2",
      divides(sp.expand((E - r) * (E - r - 1)), a2W))
check("diff-2 WITNESS: the CLEAN divisibility h h^{[1]} | a_2 FAILS"
      "  =>  restoration is REFUTED for diff-2",
      not divides(sp.expand(hW * hW1), a2W))
check("diff-2 WITNESS lies on the surviving branch b_1(r+1) = 0",
      sp.expand(b1W.subs(E, r + 1)) == 0)
check("control: the witness is a genuine point, not the zero solution",
      sp.expand(a2W) != 0 and sp.expand(b1W.subs(kappa, 1)) != 0)
# the witness is integral and its whole one-parameter scaling family works
_lam = sp.Symbol('lam')
Q4Wl = sp.expand(sh(b2W, 2) * (_lam * a2W) - sh(_lam * a2W, 2) * b2W
                 + sh(_lam * b1W, 3) * a3W - sh(a3W, 1) * (_lam * b1W))
Q3Wl = sp.expand(a3W * 0 + (sh(_lam * b1W, 2) * (_lam * a2W)
                            - sh(_lam * a2W, 1) * (_lam * b1W))
                 + kappa * hW1 * (hW2 * (_lam**2 * a1W) - hW * _lam**2 * sh(a1W, 2)))
check("diff-2 WITNESS scales: (lam a_2, lam b_1, lam^2 a_1, 0) solves Q_4 = Q_3 = 0"
      " for every lam  (a one-parameter family of counterexamples to restoration)",
      sp.expand(Q4Wl) == 0 and sp.expand(Q3Wl) == 0)


# ===========================================================================
print("\n=== S7  TASK A consequence: what survives, and bounded emptiness ===")
# ===========================================================================
def full_sector(kind, d, withQ0=True, rv=0, kv=1, extra=()):
    """Full band-3 shifted-cube sector at coefficient-degree cap d, INTEGER data.
       Q_m = delta_{m0} for m in [-6,6]; gauge b_3 = 0; wall b_2 = kappa h h^{[1]};
       genuine membership (E)_j | a_{-j}, b_{-j}."""
    hh = (E - rv) * (E - rv - (1 if kind == 'diff1' else 2))
    X = {3: sp.expand(hh * sh(hh, 1) * sh(hh, 2))}
    D = {3: sp.Integer(0), 2: sp.expand(kv * hh * sh(hh, 1))}
    un = []
    for k, nm in ((2, 'sa2'), (1, 'sa1'), (0, 'sa0')):
        p, c = gp(f'{nm}{kind}', d)
        X[k] = p
        un += c
    for k, nm in ((1, 'sb1'), (0, 'sb0')):
        p, c = gp(f'{nm}{kind}', d)
        D[k] = p
        un += c
    for j in (1, 2, 3):
        p, c = gp(f'sam{j}{kind}', d)
        X[-j] = sp.expand(fall(j) * p)
        un += c
        q, c2 = gp(f'sbm{j}{kind}', d)
        D[-j] = sp.expand(fall(j) * q)
        un += c2
    eqs = []
    for m in range(-6, 7):
        if m == 0 and not withQ0:
            continue
        v = sp.expand(Qm(X, D, m) - (1 if m == 0 else 0))
        if v != 0:
            eqs += sp.Poly(v, E).all_coeffs()
    eqs = [sp.expand(e) for e in eqs if sp.expand(e) != 0] + list(extra)
    return X, D, un, eqs


print("  [diff-1's arbitrary-degree closure (S5) predicts EMPTY at every cap --")
print("   the caps below are an independent cross-check of that theorem.]")
CAPS = [1, 2, 3]
for kind in ('diff1', 'diff2'):
    for d in CAPS:
        _, _, un, eqs = full_sector(kind, d)
        t1 = time.time()
        ok = unit_ideal_QQ(eqs, un)
        check(f"{kind}: full cascade + Q_0 = 1 + membership at cap d = {d}"
              f" (r = 0, kappa = 1) is the UNIT ideal  =>  EMPTY over C-bar"
              f"  [{len(un)} vars, {time.time() - t1:.1f}s]", ok)

# NON-VACUITY control, now in-file (the prior memo had to import it from an audit)
for kind in ('diff1', 'diff2'):
    _, _, un, eqs = full_sector(kind, 2, withQ0=False)
    ok = unit_ideal_QQ(eqs, un)
    check(f"{kind} CONTROL: dropping Q_0 = 1 leaves a PROPER ideal at cap d = 2"
          f"  =>  the sector is NONEMPTY without the moment unit (weak Nullstellensatz);"
          f" the emptiness certificates are not vacuous", not ok)

# the diff-2 SURVIVOR is a genuine point of the (Q_4, Q_3) locus at r = 0, kappa = 1,
# and it is exactly what the propagation theorem needs and cannot have.
check("diff-2 survivor: at r = 0, kappa = 1 the witness has integer coefficients"
      " a_2 = 3E(E-1), b_1 = 2(E-1)^2, a_1 = -3, b_0 = 0",
      sp.expand(a2W.subs({r: 0, kappa: 1}) - 3 * E * (E - 1)) == 0
      and sp.expand(b1W.subs({r: 0, kappa: 1}) - 2 * (E - 1)**2) == 0)
check("diff-2 survivor: the potential factorisation G = h^{[-1]} M is UNAVAILABLE --"
      " h^{[-1]} does NOT divide a_2^{[-1]} for the witness"
      "  (this is exactly why the diff-1 closure route does not transfer)",
      not divides(sh(hW.subs({r: 0, kappa: 1}), -1),
                  sh(a2W.subs({r: 0, kappa: 1}), -1)))

if HEAVY:
    for kind in ('diff1', 'diff2'):
        t1 = time.time()
        _, _, un, eqs = full_sector(kind, 4)
        ok = unit_ideal_QQ(eqs, un)
        check(f"HEAVY {kind}: full cascade + Q_0 = 1 + membership at cap d = 4"
              f" is the UNIT ideal  [{len(un)} vars, {time.time() - t1:.1f}s]", ok)
else:
    skip("HEAVY diff-1/diff-2 cap d = 4 unit-ideal emptiness", "set HEAVY=1 (~90 s)")


# ===========================================================================
print("\n=== S8  TASK B: A*-band3, the constant-top (a_3 = 1) kappa_2 != 0 sector ===")
# ===========================================================================
# Structural facts first (fast, exact, re-derived -- not cited).
Xct = {3: sp.Integer(1), 2: gp('ca2', 2)[0], 1: gp('ca1', 2)[0], 0: gp('ca0', 2)[0]}
Dct = {3: sp.Integer(0), 2: gp('cb2', 2)[0], 1: gp('cb1', 2)[0], 0: gp('cb0', 2)[0]}
check("constant top: Q_5 = (T^3 - 1) b_2   =>   b_2 = kappa_2 constant",
      sp.expand(Qm(Xct, Dct, 5) - (sh(Dct[2], 3) - Dct[2])) == 0)
Dct2 = dict(Dct)
Dct2[2] = sp.Symbol('k2')
check("constant top: Q_4 = (b_1^{[3]} - b_1) + kappa_2 (a_2 - a_2^{[2]})  (exact rung)",
      sp.expand(Qm(Xct, Dct2, 4)
                - ((sh(Dct2[1], 3) - Dct2[1])
                   + sp.Symbol('k2') * (Xct[2] - sh(Xct[2], 2)))) == 0)
# Weyl-coordinate scaling over C-bar: x->s x, del->s^-1 del, hence E=x del
# is fixed.  Ladder coefficients scale as a_k->s^k a_k, b_k->s^k b_k;
# the additional operator normalization X->s^-3 X, D->s^3 D gives net
# a_k->s^(k-3)a_k and b_k->s^(k+3)b_k.
_ss = sp.symbols('s', nonzero=True)
Xscl = {k: sp.expand(_ss**(k - 3) * v) for k, v in Xct.items()}
Dscl = {k: sp.expand(_ss**(k + 3) * v) for k, v in Dct2.items()}
check("A*-band3 scaling: [D,X]=1 is preserved rungwise under x->sx, del->s^-1del, "
      "X->s^-3X, D->s^3D",
      all(sp.expand(Qm(Xscl, Dscl, m) - _ss**m * Qm(Xct, Dct2, m)) == 0
          for m in range(-6, 7)))
check("A*-band3 scaling preserves a_3=1, E-degree caps, and (E)_j negative-band membership",
      Xscl[3] == 1
      and all(sp.degree(Xscl[k], E) == sp.degree(Xct[k], E) for k in Xct)
      and all(sp.degree(Dscl[k], E) == sp.degree(Dct2[k], E) for k in Dct2)
      and all(divides(fall(j), sp.expand(_ss**(-j - 3) * fall(j) * gp(f'sca{j}', 1)[0]))
              and divides(fall(j), sp.expand(_ss**(-j + 3) * fall(j) * gp(f'scb{j}', 1)[0]))
              for j in (1, 2, 3)))
check("A*-band3 scaling sends kappa_2 to s^5 kappa_2; choosing s^5=kappa_2^-1 "
      "normalizes nonzero kappa_2 to 1 over C-bar (not generally over Q)",
      sp.expand(Dscl[2] - _ss**5 * Dct2[2]) == 0)
# AUDIT FIX (2026-07-25): the former check here evaluated sh(1,-1) == 1, i.e.
# 1 == 1 -- a tautology that could not fail.  Replaced by the falsifiable
# content: for a NONCONSTANT h the shifted factor is a genuine nonunit (so the
# nonconstant-h route has something to extract), while for the constant top it
# is a unit -- which is exactly why that route is vacuous here.
_hnc = E * (E + 1)
check("constant top: h^{[-1]} is a UNIT (so G = M and the top-level h-forcing is "
      "VACUOUS), whereas a nonconstant h has deg h^{[-1]} > 0 -- falsifiable contrast",
      sp.expand(sh(sp.Integer(1), -1)) == 1
      and sp.Poly(sp.expand(sh(_hnc, -1)), E).degree() == 2)
# bottom proportionality and the mu_3 source
Xb = {3: sp.Integer(1)}
Db = {3: sp.Integer(0)}
for j in (1, 2, 3):
    Xb[-j] = sp.expand(fall(j) * gp(f'ta{j}', 2)[0])
    Db[-j] = sp.expand(fall(j) * gp(f'tb{j}', 2)[0])
check("constant top: Q_{-6} = b_{-3}^{[-3]} a_{-3} - a_{-3}^{[-3]} b_{-3}"
      "  =>  b_{-3} = mu_3 a_{-3} on a_{-3} != 0",
      sp.expand(Qm(Xb, Db, -6)
                - (sh(Db[-3], -3) * Xb[-3] - sh(Xb[-3], -3) * Db[-3])) == 0)
mu = sp.Symbol('mu3')
Db2 = dict(Db)
Db2[-3] = sp.expand(mu * Xb[-3])
check("constant top: with b_{-3} = mu_3 a_{-3}, Q_{-5} carries the inhomogeneous"
      " mu_3-source  mu_3 (a_{-3}^{[-2]}a_{-2} - a_{-2}^{[-3]}a_{-3})",
      sp.expand(Qm({**Xb, 2: sp.Integer(0)}, {**Db2, 2: sp.Integer(0)}, -5)
                - ((sh(Db2[-2], -3) * Xb[-3] - sh(Xb[-3], -2) * Db2[-2])
                   + mu * (sh(Xb[-3], -2) * Xb[-2] - sh(Xb[-2], -3) * Xb[-3]))) == 0)
# Lemma-P moment slope for the constant top
Xs = {3: sp.Integer(1), 2: gp('sa2', 2)[0], 1: gp('sa1', 2)[0]}
Ds = {3: sp.Integer(0), 2: sp.Symbol('k2'), 1: gp('sb1', 2)[0]}
for j in (1, 2, 3):
    Xs[-j] = sp.expand(fall(j) * gp(f'sma{j}', 2)[0])
    Ds[-j] = sp.expand(fall(j) * gp(f'smb{j}', 2)[0])
Gs = Gpot(Xs, Ds)
slope_formula = sum(Xs[i].subs(E, 0) * Ds[-i].subs(E, i)
                    - Xs[-i].subs(E, i) * Ds[i].subs(E, 0) for i in (1, 2, 3))
check("constant top: Lemma-P moment slope  G(1) = sum_i (a_i(0)b_{-i}(i) - a_{-i}(i)b_i(0))"
      "  and G(0) = 0", sp.expand(Gs.subs(E, 1) - slope_formula) == 0
      and sp.expand(Gs.subs(E, 0)) == 0)


def astar_sector(d, k2, withQ0=True, branch='all'):
    """Constant-top A*-band3 sector, cap d, INTEGER coefficients throughout.
       branch: 'all' | 'a-3=0' | ('sat', i) saturate at the i-th coeff of a_{-3}."""
    X = {3: sp.Integer(1)}
    D = {3: sp.Integer(0), 2: sp.Integer(k2)}
    un = []
    for k, nm in ((2, 'xa2'), (1, 'xa1'), (0, 'xa0')):
        p, c = gp(nm, d)
        X[k] = p
        un += c
    for k, nm in ((1, 'xb1'), (0, 'xb0')):
        p, c = gp(nm, d)
        D[k] = p
        un += c
    am3c = None
    for j in (1, 2, 3):
        p, c = gp(f'xam{j}', d)
        X[-j] = sp.expand(fall(j) * p)
        un += c
        if j == 3:
            am3c = c
        q, c2 = gp(f'xbm{j}', d)
        D[-j] = sp.expand(fall(j) * q)
        un += c2
    eqs = []
    for m in range(-6, 7):
        if m == 0 and not withQ0:
            continue
        v = sp.expand(Qm(X, D, m) - (1 if m == 0 else 0))
        if v != 0:
            eqs += sp.Poly(v, E).all_coeffs()
    eqs = [sp.expand(e) for e in eqs if sp.expand(e) != 0]
    if branch == 'a-3=0':
        sub = {c: 0 for c in am3c}
        eqs = [sp.expand(e.subs(sub)) for e in eqs]
        eqs = [e for e in eqs if e != 0]
        un = [u for u in un if u not in am3c]
    elif isinstance(branch, tuple):
        eqs = eqs + [1 - t * am3c[branch[1]]]
        un = un + [t]
    return X, D, un, eqs


# ---- msolve parser validation, IN FILE, before any load-bearing call --------
_malformed_msolve_records = (
    '[-1]garbage', '[-1]\nSECOND_RECORD', '[totally malformed', '[]garbage',
    '[garbage]', '[[nonsense]]', '[1,23,-1,[garbage]]',
    '[1,23,-1,[]]garbage', '[1,23,-1,[]][1]', '[1,0,-1,[]]',
    '[1,-23,-1,[]]', '[1,23,1,[]]', '[1,23,-1,[],0]',
    '[1]:', '[1]:123',
)
for malformed in _malformed_msolve_records:
    check(f"strict msolve char-0 parser rejects malformed output {malformed!r}",
          _parse_empty_QQ_record(malformed.replace(' ', '')) is None)
    check(f"strict msolve F_p parser rejects malformed output {malformed!r}",
          _parse_unit_Fp_record(malformed.replace(' ', '')) is None)
check("strict msolve char-0 parser accepts exact canonical EMPTY record [-1]",
      _parse_empty_QQ_record('[-1]') is True)
check("strict msolve char-0 parser accepts canonical NONEMPTY record [1,23,-1,[]]",
      _parse_empty_QQ_record('[1,23,-1,[]]') is False)
check("strict msolve char-0 parser accepts historical NONEMPTY record [1,22,-1,[]]",
      _parse_empty_QQ_record('[1,22,-1,[]]') is False)
check("strict msolve F_p parser accepts exact canonical unit-ideal record [1]",
      _parse_unit_Fp_record('[1]') is True)
_ser_x = sp.Symbol('serialization_x')
check("msolve serialization accepts rational constant coefficients by clearing"
      " their common denominator",
      _serialize_msolve_polynomial(sp.Rational(2, 3) * _ser_x**2 - sp.Rational(5, 7),
                                   [_ser_x]) == '14*serialization_x^2-15')
check("msolve serialization rejects a denominator involving an unknown before subprocess",
      _serialization_rejected(1 / (_ser_x + 1), [_ser_x]))
if HAVE_MSOLVE:
    _ver = 'version unavailable'
    for _args in (['--version'], ['-v']):
        try:
            _vr = subprocess.run([MSOLVE_PATH, *_args], capture_output=True, text=True, timeout=10)
            _txt = (_vr.stdout or _vr.stderr).strip().splitlines()
            if _txt:
                _ver = _txt[0][:160]
                break
        except (subprocess.TimeoutExpired, OSError):
            pass
    print(f"      msolve executable: {MSOLVE_PATH}; identity/version: {_ver}")
    xv, yv = sp.symbols('pv_x pv_y')
    v_unit = msolve_empty_QQ([xv - 1, xv - 2], [xv], 'validate-unit', 60)
    v_feas = msolve_empty_QQ([xv**2 - 2, yv - xv], [xv, yv], 'validate-feasible', 60)
    if v_unit is None:
        skip_solver_leg("msolve PARSER VALIDATION: known UNIT ideal",
                        "invoked msolve leg did not return a recognized verdict")
    else:
        check("msolve PARSER VALIDATION: known UNIT ideal (x-1, x-2) reports EMPTY [-1]",
              v_unit is True)
    if v_feas is None:
        skip_solver_leg("msolve PARSER VALIDATION: known FEASIBLE ideal",
                        "invoked msolve leg did not return a recognized verdict")
    else:
        check("msolve PARSER VALIDATION: known FEASIBLE ideal (x^2-2, y-x) reports NONEMPTY",
              v_feas is False)
    # AUDIT FIX (2026-07-25): (x^2-2, y-x) is REAL-rooted, so it cannot distinguish
    # "no REAL roots" from "empty over C-bar" -- precisely the failure mode that
    # would invalidate every "EMPTY over C-bar" verdict below.  Add a
    # complex-only ideal to guard what the validation purports to guard.
    v_cplx = msolve_empty_QQ([xv**2 + 1], [xv], 'validate-complex-only', 60)
    if v_cplx is None:
        skip_solver_leg("msolve PARSER VALIDATION: complex-only feasible ideal",
                        "invoked msolve leg did not return a recognized verdict")
    else:
        check("msolve PARSER VALIDATION (complex-only): (x^2+1) has NO REAL roots but is "
              "NONEMPTY over C-bar -- msolve must NOT report it as the unit ideal",
              v_cplx is False)
    _trap1 = False
    try:
        _guard_body('pv_x^2-2')
    except ValueError:
        _trap1 = True
    check("control: the trap guards accept a clean integer system (no false alarm)",
          not _trap1)
    def _guard_raises(body):
        try:
            _guard_body(body)
            return False
        except ValueError:
            return True
    check("msolve TRAP #2 GUARD (unit test): a body carrying a rational monomial"
          " '2*x^2/3' is REJECTED, never silently handed to msolve",
          _guard_raises('2*x^2/3-1'))
    check("msolve TRAP #1 GUARD (unit test): a body carrying '**' is REJECTED",
          _guard_raises('2*x**2-3'))
    check("control: the guard PASSES a clean integer body (it is not always-reject)",
          _guard_body('2*x^2-3') is True)
    _num = sp.fraction(sp.together(sp.expand(sp.Rational(2, 3) * xv**2 - 1)))[0]
    check("msolve TRAP #2: denominators are cleared to integers before hand-off"
          "  (2/3 x^2 - 1  ->  2 x^2 - 3)", sp.expand(_num - (2 * xv**2 - 3)) == 0)
    # AUDIT FIX: the former check here asserted the postcondition of a string
    # replace performed on the line immediately above it (tautological).  The
    # real trap-#1 guard is the unit test a few lines up, which rejects a body
    # still carrying '**'.
else:
    skip("msolve parser validation + all load-bearing msolve legs",
         "NOT RERUN/SKIPPED: msolve not on PATH; use --require-msolve for strict validation")

# ---- the kappa_2 = 0 control: an explicit genuine pair -----------------------
U = {1: sp.Integer(1), -1: 2 * E}                      # U = x + 2*del
U2 = cp_mul(U, U)
U3 = cp_mul(U2, U)
Xtame = cp_sub(U3, {-1: E})                            # X = U^3 - del
Dtame = dict(U)
comm = commutator(Dtame, Xtame)
check("kappa_2 = 0 CONTROL: U = x + 2 del, X = U^3 - del, D = U satisfies [D,X] = 1",
      sp.expand(comm.get(0, 0) - 1) == 0
      and all(sp.expand(v) == 0 for m, v in comm.items() if m != 0))
check("kappa_2 = 0 CONTROL: that pair has a_3 = 1, b_3 = 0, b_2 = 0"
      "  =>  the kappa_2 = 0 slice is NONEMPTY (explicit point)",
      sp.expand(Xtame.get(3, 0) - 1) == 0 and sp.expand(Dtame.get(3, 0)) == 0
      and sp.expand(Dtame.get(2, 0)) == 0)
check("kappa_2 = 0 CONTROL: the pair satisfies genuine membership (E)_j | a_{-j}, b_{-j}",
      all(divides(fall(j), Xtame.get(-j, 0)) and divides(fall(j), Dtame.get(-j, 0))
          for j in (1, 2, 3)))

# explicit, instant non-vacuity control: the POSITIVE cascade alone permits
# kappa_2 != 0 (so the kappa_2 != 0 exclusion is genuinely a negative-tail fact).
_Xp = {3: sp.Integer(1), 2: sp.Integer(5)}
_Dp = {3: sp.Integer(0), 2: sp.Integer(1), 1: sp.Integer(7)}
check("A*-band3 CONTROL (explicit point, instant): a_3 = 1, b_2 = kappa_2 = 1,"
      " a_2 = 5, b_1 = 7 constants satisfy Q_6 = Q_5 = Q_4 = 0  =>  the POSITIVE"
      " cascade permits kappa_2 != 0; the exclusion is a negative-tail fact",
      all(sp.expand(Qm(_Xp, _Dp, m)) == 0 for m in (4, 5, 6)))
check("control: a NON-constant b_2 breaks the wall Q_5 = 0 (the control is not vacuous)",
      sp.expand(Qm(_Xp, {**_Dp, 2: E}, 5)) != 0)

# ---- the emptiness legs -----------------------------------------------------
# ENGINE NOTE (recorded for the corpus): SymPy's Groebner engine does NOT finish
# the A*-band3 cap-d=1 system in 20 min, although the RIGID diff-1/diff-2 sectors
# at cap d=3 (44 vars) finish in ~4 s.  The constant top a_3 = 1 removes all
# rigidity.  msolve, by contrast, does cap d=1 in ~0.1 s.  So every A*-band3
# emptiness leg here is msolve-only; d=1 is committed, d>=2 is HEAVY-gated.
if HAVE_MSOLVE:
    _, _, un1, eqs1 = astar_sector(1, 1)
    r1 = msolve_empty_QQ(eqs1, un1, 'astar d=1 kappa2=1', 240)
    if r1 is None:
        skip_solver_leg("A*-band3 kappa_2 = 1, cap d = 1: msolve char-0 emptiness over QQ",
                        "invoked msolve leg did not return a recognized verdict")
    else:
        check("A*-band3 kappa_2 = 1, cap d = 1: msolve char-0 over QQ reports"
              " [-1] = EMPTY over C-bar   (COMMITTED, default run)", r1 is True)
    # control 1: the kappa_2 = 0 slice at the SAME cap must be NONEMPTY
    _, _, un0, eqs0 = astar_sector(1, 0)
    r0 = msolve_empty_QQ(eqs0, un0, 'astar d=1 kappa2=0', 240)
    if r0 is None:
        skip_solver_leg("A*-band3 CONTROL kappa_2 = 0 at cap d = 1 nonempty",
                        "invoked msolve leg did not return a recognized verdict")
    else:
        check("A*-band3 CONTROL: the kappa_2 = 0 slice at cap d = 1 is NONEMPTY"
              "  (so the kappa_2 = 1 emptiness is a real separation, not a broken"
              " encoding)", r0 is False)
    # F_p corroboration -- explicitly NOT a QQ proof
    rp = msolve_unit_Fp(eqs1, un1, 1073741827, 'astar d=1 F_p', 240)
    if rp is None:
        skip_solver_leg("A*-band3 cap d = 1: F_p corroboration",
                        "invoked msolve leg did not return a recognized verdict")
    else:
        check("A*-band3 kappa_2 = 1, cap d = 1: unit ideal mod p = 2^30+3."
              "  CORROBORATION ONLY -- unit mod p does NOT imply unit over QQ"
              " (counterexample: (p x - 1)).", rp is True)

if HAVE_MSOLVE and HEAVY:
    # the msolve-based non-vacuity control (redundant with the two committed
    # controls above, and slow: the no-Q_0 variety is positive-dimensional)
    _, _, unv, eqsv = astar_sector(1, 1, withQ0=False)
    rv = msolve_empty_QQ(eqsv, unv, 'astar d=1 no-Q_0', 1200)
    if rv is None:
        skip_solver_leg("HEAVY A*-band3 CONTROL: non-vacuity (drop Q_0 = 1) at cap d = 1",
                        "invoked msolve leg did not return a recognized verdict")
    else:
        check("HEAVY A*-band3 CONTROL: dropping Q_0 = 1 at cap d = 1 leaves a"
              " NONEMPTY variety  =>  the emptiness certificate is not vacuous",
              rv is False)
    for d, tmo in ((2, 3000), (3, 3000)):
        for br, brn in (('a-3=0', 'a_{-3} = 0'), ('all', 'full')):
            _, _, un, eqs = astar_sector(d, 1, branch=br)
            res = msolve_empty_QQ(eqs, un, f'astar d={d} {brn}', tmo)
            if res is None:
                skip_solver_leg(f"HEAVY A*-band3 kappa_2 = 1, cap d = {d}, branch {brn}:"
                                f" msolve char-0 emptiness over QQ",
                                "invoked msolve leg did not return a recognized verdict")
            else:
                check(f"HEAVY A*-band3 kappa_2 = 1, cap d = {d}, branch {brn}:"
                      f" msolve char-0 reports"
                      f" {'EMPTY' if res else 'NONEMPTY'} over C-bar", res is True)
elif HAVE_MSOLVE:
    skip("HEAVY A*-band3: caps d = 2, 3 emptiness (branch a_{-3} = 0 and full)"
         " + the msolve non-vacuity control",
         "set HEAVY=1; msolve char-0, hard timeout 3000 s per leg")


# ===========================================================================
print("\n=== S9  TASK C: slope forcing for the constant-top sector (Rabinowitsch) ===")
# ===========================================================================
# Is the Lemma-P moment slope G(1) forced to 0 by the TAIL ALONE (no Q_0)?
# If yes, then Q_0 = 1 (which needs G(1) = 1) is impossible -- a degree-free target.
def astar_slope_system(d, k2):
    X = {3: sp.Integer(1)}
    D = {3: sp.Integer(0), 2: sp.Integer(k2)}
    un = []
    for k, nm in ((2, 'ya2'), (1, 'ya1'), (0, 'ya0')):
        p, c = gp(nm, d)
        X[k] = p
        un += c
    for k, nm in ((1, 'yb1'), (0, 'yb0')):
        p, c = gp(nm, d)
        D[k] = p
        un += c
    for j in (1, 2, 3):
        p, c = gp(f'yam{j}', d)
        X[-j] = sp.expand(fall(j) * p)
        un += c
        q, c2 = gp(f'ybm{j}', d)
        D[-j] = sp.expand(fall(j) * q)
        un += c2
    eqs = []
    for m in range(-6, 7):
        if m == 0:
            continue
        v = sp.expand(Qm(X, D, m))
        if v != 0:
            eqs += sp.Poly(v, E).all_coeffs()
    eqs = [sp.expand(e) for e in eqs if sp.expand(e) != 0]
    slope = sp.expand(sum(X[i].subs(E, 0) * D[-i].subs(E, i)
                          - X[-i].subs(E, i) * D[i].subs(E, 0) for i in (1, 2, 3)))
    return un, eqs, slope


# TASK C VERDICT LOGIC.  Rabinowitsch: G(1) vanishes on V(tail)  <=>
# (tail, 1 - t G(1)) is the UNIT ideal (msolve char-0 prints [-1]).  A NONEMPTY
# result is an existence statement -- there is a genuine point of the tail with
# G(1) != 0 -- and therefore RULES THE SLOPE ROUTE OUT at that cap.  That is a
# radical-correct negative verdict, not a normal-form inference.
if HAVE_MSOLVE:
    un, eqs, slope = astar_slope_system(1, 1)
    res = msolve_empty_QQ(eqs + [1 - t * slope], un + [t], 'slope-forcing d=1', 400)
    if res is None:
        skip_solver_leg("TASK C, cap d = 1: Rabinowitsch slope-forcing probe",
                        "invoked msolve leg did not return a recognized verdict")
    else:
        check("TASK C, formal constant-top A*-band3 sector, kappa_2=1, cap d=1: "
              "the Rabinowitsch system (tail, 1-tG(1)) is NOT the unit ideal, so "
              "G(1) is not forced to 0 and this slope route fails in exactly that sector/cap",
              res is False)
        print("      TASK C VERDICT (formal constant-top A*-band3, kappa_2=1, cap d=1):")
        print("      G(1) is NOT forced to 0 by the tail alone. This bounded formal probe")
        print("      does not support a claim beyond this sector and cap.")

if HAVE_MSOLVE and HEAVY:
    un, eqs, slope = astar_slope_system(2, 1)
    res = msolve_empty_QQ(eqs + [1 - t * slope], un + [t], 'slope-forcing d=2', 3000)
    if res is None:
        skip_solver_leg("HEAVY TASK C, cap d = 2: Rabinowitsch slope-forcing probe",
                        "invoked msolve leg did not return a recognized verdict")
    else:
        check(f"HEAVY TASK C, cap d = 2: G(1) is"
              f" {'FORCED to 0' if res else 'NOT forced to 0'} by the tail alone",
              res is False)
elif HAVE_MSOLVE:
    skip("HEAVY TASK C: slope-forcing Rabinowitsch probe at cap d = 2", "set HEAVY=1")
if not HAVE_MSOLVE:
    skip("TASK C: formal constant-top A*-band3, kappa_2=1, cap d=1 slope probe",
         "NOT RERUN/SKIPPED: msolve not on PATH")


# ===========================================================================
print("\n=== SUMMARY ===")
# ===========================================================================
npass = sum(1 for _, ok in CHECKS if ok)
nfail = len(CHECKS) - npass
print(f"  checks executed : {len(CHECKS)}   passed: {npass}   failed: {nfail}")
print(f"  skipped         : {len(SKIPS)}")
for nm, why in SKIPS:
    print(f"     SKIP  {nm}   ({why})")
print(f"  wall time       : {time.time() - T0:.1f}s   HEAVY={'1' if HEAVY else '0'}"
      f"   require-msolve={'1' if REQUIRE_MSOLVE else '0'}"
      f"   msolve={MSOLVE_PATH if HAVE_MSOLVE else 'NOT RERUN/SKIPPED'}")
if nfail:
    print("\n  *** SOME CHECKS FAILED ***")
    for nm, ok in CHECKS:
        if not ok:
            print(f"     FAIL  {nm}")
    sys.exit(1)
if SKIPS:
    print("\n  ALL EXECUTED CHECKS PASSED; OPTIONAL CHECKS SKIPPED (listed above).")
else:
    print("\n  ALL CHECKS PASSED.")
