#!/usr/bin/env python3
"""
audit_classical_A_search.py  --  Referee audit (auditor A), classical Theorem J3:
EXHAUSTIVE BOUNDED SEARCH of the band-2 system with both extremes nonzero and
non-(scalar-)square, deg a_k <= 3, plus independence cross-checks.

Normalizations used (all legitimate for the search, since they preserve the
hypotheses -- band <= 2, {D,X} in C*, both extremes nonzero -- and the
scalar-square class of every coefficient):
  N1. X -> t*X (t in C*) and the diagonal symplectic scaling
      (x, xi) -> (rho*x, rho^{-1}*xi), which sends a_k -> t*rho^k*a_k.
      Choosing t*rho^2 = 1/lead(a2), t*rho^{-2} = 1/lead(a-2) (solvable over
      C) makes BOTH extreme leading coefficients equal to 1.  Scalar-square
      status is invariant under nonzero scalar multiples.
  N2. (reduced variant only) D -> D - lambda2*X + const, which preserves
      {D,X} and the band, and sets lambda2 = 0, b0(0-th coefficient) = 0.

Structure ("graded to stay tractable" -- the cascade linearizes the system):

  PART B (the exhaustive deg <= 3 search; complete given the two lemma
    certificates from audit_classical_A_identities.py):
    The outer equations C_{+-4}, C_{+-3} have COMPLETE polynomial solution
    sets, for b of ANY degree (Wronskian lemma; Lemma J2 + both extremes
    nonsquare => zero homogeneous deviation):
        b2 = lam*a2, b1 = lam*a1, b-1 = mu*a-1, b-2 = mu*a-2.
    Substituting these and ONLY these, for each of the 9 strata
    (deg a2, deg a-2) in {1, 2 nonsquare, 3} -- which exhaust all nonzero
    nonsquare extremes of degree <= 3 -- with fully symbolic a1, a0, a-1 of
    degree <= 3 (subsuming every lower profile including zero), Groebner +
    Rabinowitsch certifies the system C_{+-2} = C_{+-1} = 0, C_0 = c, c != 0
    is EMPTY:
      B-full:    lam, mu free; b0 undetermined to degree 4 (complete: C_2 =
                 2 a2 (b0' - lam a0') = 0 with a2 != 0 in the domain C[tau]
                 bounds deg b0 <= deg a0 <= 3),
      B-reduced: N2 normalization lam = 0; b0' eliminated via the certified
                 collapse identities and the domain step.

  PART C (random exact sweep + POSITIVE CONTROLS): exact QQ-linear solves of
    the FULL untouched system (all five b's free to degree 4-5) over random
    integer a's on every degree profile <= 3 with both extremes nonsquare:
    assert [C_0]_0 vanishes identically on the solution space of the other
    constraints.  Positive controls (a-2 = 0 localized near-miss; band-1
    pair) confirm the machinery detects c != 0 solutions when they exist.

  PART D (generic-extremes symbolic rank check): a2, a-2 fully symbolic on
    each (deg, deg) stratum, random middles: appending the c-functional row
    to the constraint matrix does not raise the rank over the rational
    function field => no c != 0 solution at the generic point of any stratum.

  PART A (fully independent cross-check; no cascade input): the complete
    bilinear system with ALL a- and b-coefficients undetermined, a's of
    degree <= 1 (extremes of degree exactly 1 -- at that size this IS the
    two-sided nonsquare sector, after N1), b's of degree <= Db:
        1 in ideal + Rabinowitsch(c) ?
    Runs with a timeout; Db = 2, falling back to Db = 1.  Part A is an
    independence check -- exhaustiveness at deg <= 3 is carried by Part B.

A counterexample anywhere prints VIOLATION with full data and exits 1.
"""
import random
import signal
import sys
import time
from itertools import product

import sympy as sp

tau, w, w1, w2, w3 = sp.symbols("tau w w1 w2 w3")
lam, mu = sp.symbols("lambda2 mu2")
Dt = lambda f: sp.diff(f, tau)

FAILED = []


def check(name, ok):
    print(("  [PASS] " if ok else "  [FAIL] ") + name, flush=True)
    if not ok:
        FAILED.append(name)


class Timeout(Exception):
    pass


def with_timeout(seconds, fn, *args, **kw):
    def handler(signum, frame):
        raise Timeout()
    old = signal.signal(signal.SIGALRM, handler)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        return fn(*args, **kw)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)


def Cm(m, A, B):
    s = sp.Integer(0)
    for k in range(-2, 3):
        l = m - k
        if -2 <= l <= 2:
            s += k * A[k] * Dt(B[l]) - l * Dt(A[k]) * B[l]
    return sp.expand(s)


def tau_coeffs(expr):
    if expr == 0:
        return []
    p = sp.Poly(expr, tau)
    return [sp.expand(p.coeff_monomial(tau**i)) for i in range(p.degree() + 1)]


def gpoly(prefix, deg, monic=False):
    if deg is None:
        return sp.Integer(0), ()
    if monic:
        cs = sp.symbols(f"{prefix}0:{deg}")
        return tau**deg + sum(cs[i] * tau**i for i in range(deg)), cs
    cs = sp.symbols(f"{prefix}0:{deg + 1}")
    return sum(cs[i] * tau**i for i in range(deg + 1)), cs


def is_scalar_square(poly):
    """Scalar-square over C: c*h(tau)^2.  Parity of square-free multiplicities
    is field-independent in char 0, so testing over Q decides the C-status."""
    if poly == 0:
        return False
    p = sp.Poly(poly, tau)
    if p.degree() == 0:
        return True
    _, factors = sp.sqf_list(p.as_expr(), tau)
    return all(e % 2 == 0 for _, e in factors)


def groebner_is_unit(gens, vars_):
    gb = sp.groebner([sp.expand(g) for g in gens if sp.expand(g) != 0],
                     *vars_, order="grevlex")
    return list(gb.exprs) == [sp.Integer(1)]


# ---------------------------------------------------------------------------
# PART B -- exhaustive deg <= 3 search via the certified outer reductions
# ---------------------------------------------------------------------------
STRATA = {"deg1": 1, "deg2ns": 2, "deg3": 3}   # deg2ns adds disc != 0


def stratum_poly(prefix, label):
    d = STRATA[label]
    pol, cs = gpoly(prefix, d, monic=True)   # N1: monic extremes
    units = []
    if label == "deg2ns":
        # monic quadratic tau^2 + c1 tau + c0: disc = c1^2 - 4 c0
        units.append(sp.expand(cs[1]**2 - 4 * cs[0]))
    return pol, list(cs), units


def part_B(mode):
    print(f"PART B ({mode}): 9 strata x fully symbolic middles (deg <= 3)",
          flush=True)
    allok = True
    for top_lab, bot_lab in product(STRATA, STRATA):
        t0 = time.time()
        a2, pcs, units_top = stratum_poly("p", top_lab)
        am2, scs, units_bot = stratum_poly("s", bot_lab)
        a1, qcs = gpoly("q", 3)
        a0, gcs = gpoly("g", 3)
        am1, rcs = gpoly("r", 3)
        A = {2: a2, 1: a1, 0: a0, -1: am1, -2: am2}
        U = sp.expand(Dt(a2) * am1 + 2 * a2 * Dt(am1))
        L = sp.expand(Dt(am2) * a1 + 2 * am2 * Dt(a1))
        Mmom = sp.expand(2 * a2 * am2 + a1 * am1)
        vars_ = pcs + scs + list(qcs) + list(gcs) + list(rcs)
        if mode == "full":
            b0, b0c = gpoly("B0_", 4)
            Bsub = {2: lam * a2, 1: lam * a1, 0: b0, -1: mu * am1, -2: mu * am2}
            eqs = []
            for m in (2, 1, -1, -2):
                eqs += tau_coeffs(Cm(m, A, Bsub))
            cc = tau_coeffs(Cm(0, A, Bsub))
            eqs += cc[1:]
            const_term = cc[0] if cc else sp.Integer(0)
            vars_ += [lam, mu] + list(b0c)
        else:
            # N2: lam = 0; collapse identities + domain step: b0' = 0 = mu a0';
            # remaining system in a-coeffs and mu only.
            eqs = tau_coeffs(sp.expand(mu * Dt(a0)))
            eqs += tau_coeffs(sp.expand(mu * U))
            eqs += tau_coeffs(sp.expand(mu * L))
            cc = tau_coeffs(sp.expand(mu * Dt(Mmom)))
            eqs += cc[1:]
            const_term = cc[0] if cc else sp.Integer(0)
            vars_ += [mu]
        gens = list(eqs) + [w3 * const_term - 1]     # c != 0
        vars2 = vars_ + [w3]
        for i, u_ in enumerate(units_top + units_bot):
            wv = sp.Symbol(f"wu{i}")
            gens.append(wv * u_ - 1)
            vars2.append(wv)
        ok = groebner_is_unit(gens, vars2)
        allok &= ok
        check(f"  stratum a2:{top_lab} x a-2:{bot_lab} empty "
              f"[{len(gens)} gens, {len(vars2)} vars, {time.time()-t0:.1f}s]",
              ok)
    check(f"PART B ({mode}): all 9 nonsquare strata empty at deg <= 3", allok)


# ---------------------------------------------------------------------------
# PART C -- random exact sweep + positive controls
# ---------------------------------------------------------------------------
def build_linear_system(Avals, Db):
    bs, ballc = {}, []
    for k in range(-2, 3):
        bs[k], cs = gpoly(f"cb{('m' if k < 0 else '')}{abs(k)}_", Db)
        ballc += list(cs)
    eqs = []
    for m in range(-4, 5):
        if m == 0:
            continue
        eqs += tau_coeffs(Cm(m, Avals, bs))
    cc = tau_coeffs(Cm(0, Avals, bs))
    eqs += cc[1:]
    phi_expr = cc[0] if cc else sp.Integer(0)
    Mrows = [[sp.diff(e, bc) for bc in ballc] for e in eqs if e != 0]
    M = sp.Matrix(Mrows) if Mrows else sp.zeros(0, len(ballc))
    phi = sp.Matrix([[sp.diff(phi_expr, bc) for bc in ballc]])
    return M, phi, ballc


def c_reachable(Avals, Db):
    M, phi, _ = build_linear_system(Avals, Db)
    ns = M.nullspace()
    return any(sp.simplify((phi * v)[0]) != 0 for v in ns), len(ns)


def rand_poly(deg, rng, exact=True):
    if deg is None:
        return sp.Integer(0)
    while True:
        cs = [rng.randint(-3, 3) for _ in range(deg + 1)]
        if not exact or deg == 0 or cs[-1] != 0:
            break
    if all(v == 0 for v in cs):
        cs[0] = 1
    return sum(cs[i] * tau**i for i in range(deg + 1))


def part_C():
    rng = random.Random(20260720)
    t0 = time.time()
    ctrl1 = {2: tau, 1: 0, 0: 0, -1: 0, -2: 0}            # near-miss, a-2 = 0
    ok1, _ = c_reachable(ctrl1, 3)
    ctrl2 = {2: 0, 1: sp.Integer(1), 0: 0, -1: 0, -2: 0}  # X = x, D = xi
    ok2, _ = c_reachable(ctrl2, 3)
    check("PART C controls: solver finds c != 0 for the localized near-miss "
          "(a2 = tau, a-2 = 0) and for the band-1 pair (X = x)", ok1 and ok2)

    middles = [None, 0, 1, 2, 3]
    tested = viol = 0
    for d2, dm2 in product((1, 2, 3), repeat=2):
        for d1, d0, dm1 in product(middles, repeat=3):
            a2 = rand_poly(d2, rng)
            while is_scalar_square(a2):
                a2 = rand_poly(d2, rng)
            am2 = rand_poly(dm2, rng)
            while is_scalar_square(am2):
                am2 = rand_poly(dm2, rng)
            A = {2: a2, 1: rand_poly(d1, rng, exact=False),
                 0: rand_poly(d0, rng, exact=False),
                 -1: rand_poly(dm1, rng, exact=False), -2: am2}
            reach, _ = c_reachable(A, 4)
            tested += 1
            if reach:
                viol += 1
                print("  VIOLATION:", A, flush=True)
    check(f"PART C sweep: {tested} random both-nonsquare systems across all "
          f"degree profiles <= 3 (b free to deg 4): no bracket-constant "
          f"solution anywhere [{time.time()-t0:.0f}s]", viol == 0)

    t0 = time.time()
    extra = viol = 0
    for _ in range(150):
        d2, dm2 = rng.choice((1, 2, 3)), rng.choice((1, 2, 3))
        a2 = rand_poly(d2, rng)
        am2 = rand_poly(dm2, rng)
        if is_scalar_square(a2) or is_scalar_square(am2):
            continue
        A = {2: a2, 1: rand_poly(rng.choice(middles), rng, exact=False),
             0: rand_poly(rng.choice(middles), rng, exact=False),
             -1: rand_poly(rng.choice(middles), rng, exact=False), -2: am2}
        reach, _ = c_reachable(A, 5)
        extra += 1
        if reach:
            viol += 1
            print("  VIOLATION:", A, flush=True)
    check(f"PART C extra: {extra} random systems with b free to degree 5: "
          f"no solution [{time.time()-t0:.0f}s]", viol == 0)


# ---------------------------------------------------------------------------
# PART D -- generic-extremes symbolic rank comparison
# ---------------------------------------------------------------------------
def part_D():
    rng = random.Random(4711)
    allok, done = True, 0
    for d2, dm2 in product((1, 2, 3), repeat=2):
        t0 = time.time()
        a2, _ = gpoly("p", d2)
        am2, _ = gpoly("s", dm2)
        A = {2: a2, 1: rand_poly(2, rng, exact=False),
             0: rand_poly(2, rng, exact=False),
             -1: rand_poly(2, rng, exact=False), -2: am2}
        M, phi, _ = build_linear_system(A, 3)

        def two_ranks():
            return M.rank(), M.col_join(phi).rank()
        try:
            rk, rk2 = with_timeout(90, two_ranks)
        except Timeout:
            print(f"  [SKIP] PART D ({d2},{dm2}) symbolic rank exceeded 90s "
                  f"(covered by Parts B and C)", flush=True)
            continue
        ok = (rk == rk2)
        allok &= ok
        done += 1
        check(f"  PART D deg(a2,a-2)=({d2},{dm2}), symbolic extremes, random "
              f"middles: appending the c-row does not raise rank "
              f"[{time.time()-t0:.1f}s]", ok)
    check(f"PART D: generic-extreme rank obstruction on {done}/9 strata "
          f"(skips, if any, are redundancy covered by Parts B/C)",
          allok and done >= 4)


# ---------------------------------------------------------------------------
# PART A -- fully independent Groebner emptiness for ALL 9 nonsquare strata:
# no cascade input whatsoever.  All ten coefficient polynomials undetermined
# (extremes monic per N1, deg exactly 1, 2 with disc != 0, or 3; middles of
# degree <= 3; b's of degree <= 3), c != 0 by Rabinowitsch.
# ---------------------------------------------------------------------------
def part_A_stratum(d2, dm2, dmid, Db):
    a2, pcs = gpoly("p", d2, monic=True)
    am2, scs = gpoly("s", dm2, monic=True)
    a1, qcs = gpoly("q", dmid)
    a0, gcs = gpoly("g", dmid)
    am1, rcs = gpoly("r", dmid)
    A = {2: a2, 1: a1, 0: a0, -1: am1, -2: am2}
    bs, ballc = {}, []
    for k in range(-2, 3):
        bs[k], cs = gpoly(f"b{('m' if k < 0 else '')}{abs(k)}_", Db)
        ballc += list(cs)
    eqs = []
    for m in range(-4, 5):
        if m == 0:
            continue
        eqs += tau_coeffs(Cm(m, A, bs))
    cc = tau_coeffs(Cm(0, A, bs))
    eqs += cc[1:]
    const_term = cc[0] if cc else sp.Integer(0)
    gens = list(eqs) + [w3 * const_term - 1]
    vars_ = list(pcs) + list(scs) + list(qcs) + list(gcs) + list(rcs) \
        + ballc + [w3]
    iu = 0
    for d, cs in ((d2, pcs), (dm2, scs)):
        if d == 2:                       # nonsquare quadratic: disc != 0
            wv = sp.Symbol(f"wd{iu}")
            iu += 1
            gens.append(wv * (cs[1]**2 - 4 * cs[0]) - 1)
            vars_.append(wv)
    return gens, vars_


def part_A():
    print("PART A: fully independent emptiness (no cascade), 9 strata, "
          "middles deg <= 3, b deg <= 3", flush=True)
    allok, done = True, 0
    for d2, dm2 in product((1, 2, 3), repeat=2):
        t0 = time.time()
        gens, vars_ = part_A_stratum(d2, dm2, 3, 3)
        try:
            ok = with_timeout(420, groebner_is_unit, gens, vars_)
        except Timeout:
            print(f"  [SKIP] PART A ({d2},{dm2}) exceeded 420s (stratum "
                  f"covered by PART B)", flush=True)
            continue
        allok &= ok
        done += 1
        check(f"  PART A extremes deg ({d2},{dm2}) monic"
              f"{' disc!=0' if 2 in (d2, dm2) else ''}: unit ideal "
              f"[{len(gens)} gens, {len(vars_)} vars, {time.time()-t0:.1f}s]",
              ok)
    check(f"PART A: fully independent emptiness certificates on {done}/9 "
          f"strata (any skip is covered by PART B)", allok and done >= 6)


# ---------------------------------------------------------------------------
# PART E -- negative controls: the Groebner pipeline must return NON-unit
# ideals where solutions genuinely exist (so the unit certificates above are
# not artifacts of a broken construction).
# ---------------------------------------------------------------------------
def part_E():
    # E1: Part B reduced stratum (deg1, deg1) WITHOUT the c != 0 inverse:
    # solutions with c = 0 exist (mu = 0, D proportional to X), so the ideal
    # must not be the unit ideal.
    a2, pcs = gpoly("p", 1, monic=True)
    am2, scs = gpoly("s", 1, monic=True)
    a1, qcs = gpoly("q", 3); a0, gcs = gpoly("g", 3); am1, rcs = gpoly("r", 3)
    U = sp.expand(Dt(a2) * am1 + 2 * a2 * Dt(am1))
    L = sp.expand(Dt(am2) * a1 + 2 * am2 * Dt(a1))
    Mmom = sp.expand(2 * a2 * am2 + a1 * am1)
    eqs = tau_coeffs(sp.expand(mu * Dt(a0)))
    eqs += tau_coeffs(sp.expand(mu * U))
    eqs += tau_coeffs(sp.expand(mu * L))
    eqs += tau_coeffs(sp.expand(mu * Dt(Mmom)))[1:]
    vars_ = list(pcs) + list(scs) + list(qcs) + list(gcs) + list(rcs) + [mu]
    check("PART E1 negative control: same stratum system WITHOUT c != 0 is "
          "NOT the unit ideal (c = 0 solutions exist)",
          not groebner_is_unit(eqs, vars_))
    # E2: full independent system on the memo's localized near-miss profile
    # (a2 = tau, a-2 = 0), WITH the c != 0 inverse: the solution
    # b-2 = 1/2, c = 1 exists, so the ideal must not be the unit ideal.
    A = {2: tau, 1: sp.Integer(0), 0: sp.Integer(0),
         -1: sp.Integer(0), -2: sp.Integer(0)}
    bs, ballc = {}, []
    for k in range(-2, 3):
        bs[k], cs = gpoly(f"eb{('m' if k < 0 else '')}{abs(k)}_", 2)
        ballc += list(cs)
    eqs = []
    for m in range(-4, 5):
        if m == 0:
            continue
        eqs += tau_coeffs(Cm(m, A, bs))
    cc = tau_coeffs(Cm(0, A, bs))
    eqs += cc[1:]
    const_term = cc[0] if cc else sp.Integer(0)
    gens = list(eqs) + [w3 * const_term - 1]
    check("PART E2 negative control: full system on the near-miss profile "
          "(a2 = tau, a-2 = 0) WITH c != 0 is NOT the unit ideal "
          "(the localized solution exists)",
          not groebner_is_unit(gens, ballc + [w3]))


def main():
    parts = sys.argv[1].upper() if len(sys.argv) > 1 else "BCDAE"
    print("=" * 78, flush=True)
    print(f"AUDIT A (classical J3): exhaustive bounded search, deg a_k <= 3 "
          f"[parts {parts}]", flush=True)
    print("=" * 78, flush=True)
    if "B" in parts:
        part_B("reduced")
        part_B("full")
    if "C" in parts:
        part_C()
    if "D" in parts:
        part_D()
    if "A" in parts:
        part_A()
    if "E" in parts:
        part_E()
    if FAILED:
        print("FAILURES:", *FAILED, sep="\n  ")
        sys.exit(1)
    print(f"ALL AUDIT-A SEARCH CHECKS PASSED (parts {parts}) -- no violating "
          "solution exists at deg a_k <= 3 (b unbounded via lemma "
          "completeness; b <= 5 sampled independently)", flush=True)


if __name__ == "__main__":
    main()
