#!/usr/bin/env python3
"""Independent verification of the affine-quotient obstruction for F3, F4, F5.

This script is written independently of run_obstruction.py and certifies the
same claim (Outcome B: no affine-linear quotient (pi, G) exists for F3, F4, F5)
by DIFFERENT routes:

  V1. Reconstruct F3, F4, F5 from paper/main.tex from scratch (own antiderivative
      code, exact multivariate division with remainder-0 assertion, own Jacobian
      determinant), verify the paper's polynomiality criterion (c1)-(c3), the
      component degrees and leading forms, and compare the EXPANDED POLYNOMIALS
      term-by-term against the builder's family() construction.

  V2. Re-derive the structural facts the builder's method rests on:
      (a) chart completeness: every rank-2 real/rational 2x3 matrix has RREF of
          exactly one of the three shapes C12/C13/C23 (own Fraction RREF on
          random + adversarial matrices, gauge factor M = A*R checked exactly,
          kernel formulas checked symbolically);
      (b) Lemma 1 spot checks on random affine pi (both directions, with the
          converse checked through an explicit affine coordinate inversion,
          not through the lemma);
      (c) "existence of G for fixed pi is a LINEAR system": own exact
          Fraction-arithmetic Gaussian elimination; on a randomly gauged copy
          of a control map the system is solved and (pi, G) verified exactly;
          on F3 with a random translated rank-2 pi the system is refuted
          exactly over Q.

  V3. MAIN independent certification (derivative-free substitution route).
      For each chart, complete (l1, l2) to an invertible affine coordinate
      system and substitute; membership l_i o F in C[l1, l2] holds iff the
      transformed polynomial has NO monomial containing the third coordinate.
      The offending coefficients generate an ideal in Q[params]:
        C23 (0 params): all generators are explicit nonzero rationals ->
             infeasible by inspection;
        C13 (1 param): own univariate gcd over Fraction; gcd = nonzero
             constant -> ideal = Q[p] -> no complex p (no Groebner engine);
        C12 (2 params): a generator equal to a nonzero rational constant
             (unit certificate) -> infeasible for all complex (p, q).
      All 9 (map, chart) cells are covered; unit-certificate positions are
      checked against the leading-form prediction.

  V4. Third flavour for the 2-parameter chart C12 (exact point evaluation):
      if G existed, l_i o F would be constant on every fiber P + s*v of pi.
      The differences l_i o F(P+v) - l_i o F(P) at explicit integer points are
      polynomials in Q[p, q]; their Groebner basis is [1] (different
      generators than the builder's coefficient ideal). For C23, explicit
      integer point pairs in one pi-fiber with different images are printed.

  V5. Positive controls CTRL1/CTRL2 pushed through THIS script's substitution
      route: the unique feasible chart/parameter is found, G is reconstructed,
      pi o F = G o pi verified exactly, and Jac(G) classified.

  V6. Audit of the builder's artifacts: the memo table (generator counts, unit
      counts, sample certificates, GB = [1]) is recomputed with an own
      implementation of the fiber-derivative route and compared against both
      obstruction-memo.md (hardcoded here) and certificates.json.

Everything is exact (sympy rationals / Fractions); no floating point anywhere.

Run:  uv run --with sympy python verify_obstruction_independent.py
"""

import importlib.util
import json
import os
import random
import sys
import time
from fractions import Fraction

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))

x, y, z = sp.symbols("x y z")
NU, T = sp.symbols("NU T")
p, q = sp.symbols("p q")
X, Y, Z = sp.symbols("X Y Z")
w1, w2 = sp.symbols("w1 w2")

t0 = time.time()
FAILURES = []


def check(label, ok):
    status = "ok" if ok else "FAIL"
    print("  [%s] %s" % (status, label))
    if not ok:
        FAILURES.append(label)
    return ok


# ===========================================================================
# V1. Reconstruction from paper/main.tex, from scratch.
# ===========================================================================
print("== V1: reconstruct F3, F4, F5 from paper/main.tex ==")

def antider_xE(E):
    """A(nu) = int_0^nu xi*E(xi) dxi, by hand from the coefficient list."""
    Pn = sp.Poly(E, NU)
    return sp.Add(*[c * NU ** (k + 2) / (k + 2) for (k,), c in Pn.terms()])

def antider_2E(E):
    """C(nu) = 2*int_0^nu E(xi) dxi, by hand."""
    Pn = sp.Poly(E, NU)
    return sp.Add(*[2 * c * NU ** (k + 1) / (k + 1) for (k,), c in Pn.terms()])

def build_member(E, h):
    """F_{E,h} per main.tex eq. (family): exact division, remainder asserted 0."""
    A = antider_xE(E)
    C = antider_2E(E)
    g = sp.expand(h.subs(T, x * y) - x**2 * z)
    u = 1 + x * y
    nuval = sp.expand(u * g)
    S = sp.expand(A.subs(NU, nuval) + u * g**2)
    Tt = sp.expand(C.subs(NU, nuval) + 2 * g)
    q1, r1 = sp.div(S, sp.expand(x**2 * g**2), x, y, z)
    q2, r2 = sp.div(Tt, sp.expand(x * g), x, y, z)
    assert r1 == 0 and r2 == 0, "division not exact: F not polynomial?"
    return [sp.expand(q1), sp.expand(q2), sp.expand(x * g)]

def det_jac(F):
    """Own 3x3 cofactor expansion."""
    d = [[sp.diff(f, v) for v in (x, y, z)] for f in F]
    return sp.expand(
        d[0][0] * (d[1][1] * d[2][2] - d[1][2] * d[2][1])
        - d[0][1] * (d[1][0] * d[2][2] - d[1][2] * d[2][0])
        + d[0][2] * (d[1][0] * d[2][1] - d[1][1] * d[2][0]))

def leading_form(f):
    P = sp.Poly(f, x, y, z)
    d = P.total_degree()
    return d, sorted([(m, c) for m, c in P.terms() if sum(m) == d])

# The three members exactly as defined in paper/main.tex.
MEMBERS = {
    "F3": (2 - 6 * NU, 1 - sp.Rational(3, 2) * T, 1),
    "F4": (3 - 12 * NU + 6 * NU**2, 1 - 2 * T, 2),
    "F5": (4 - 15 * NU + 10 * NU**3, sp.Integer(1), 3),
}

# main.tex says the n>=2 Legendre members are E = 2 - 6 nu + L_n(2 nu - 1):
E4_legendre = sp.expand(2 - 6 * NU + sp.legendre(2, 2 * NU - 1))
check("F4's E equals paper's Legendre recipe 2 - 6 nu + L_2(2 nu - 1)",
      sp.expand(E4_legendre - MEMBERS["F4"][0]) == 0)

# Expected sanity data (independent statement of what the memo claims).
EXPECTED = {
    "F3": {"degs": (7, 6, 4),
           "LF": [((3, 3, 1), sp.Integer(2)), ((3, 2, 1), sp.Integer(6)),
                  ((3, 0, 1), sp.Integer(-1))]},
    "F4": {"degs": (12, 11, 4),
           "LF": [((6, 4, 2), sp.Rational(3, 2)), ((6, 3, 2), sp.Integer(4)),
                  ((3, 0, 1), sp.Integer(-1))]},
    "F5": {"degs": (17, 16, 4),
           "LF": [((9, 5, 3), sp.Integer(-2)), ((9, 4, 3), sp.Integer(-5)),
                  ((3, 0, 1), sp.Integer(-1))]},
}

built = {}
for name, (E, h, n) in MEMBERS.items():
    # Paper Theorem (telescoping family): polynomial iff a != 0 and (c1)-(c3).
    a = h.subs(T, 0)
    A = antider_xE(E)
    C = antider_2E(E)
    c1 = sp.simplify(A.subs(NU, a) + a**2) == 0
    c2 = sp.simplify(C.subs(NU, a) + 2 * a) == 0
    Ea = E.subs(NU, a)
    c3 = sp.simplify(1 + Ea + sp.diff(h, T).subs(T, 0) * (Ea + 2) / a) == 0
    check("%s satisfies (c1),(c2),(c3), a = h(0) = %s != 0" % (name, a),
          c1 and c2 and c3 and a != 0)
    F = build_member(E, h)
    built[name] = F
    check("%s: det Jac = -2 (own cofactor expansion)" % name, det_jac(F) == -2)
    degs, lfs = [], []
    for f in F:
        d, lf = leading_form(f)
        degs.append(d)
        lfs.append(lf)
    check("%s: component degrees %s" % (name, tuple(degs)),
          tuple(degs) == EXPECTED[name]["degs"])
    ok_lf = all(len(lf) == 1 and lf[0][0] == EXPECTED[name]["LF"][i][0]
                and lf[0][1] == EXPECTED[name]["LF"][i][1]
                for i, lf in enumerate(lfs))
    check("%s: leading forms match memo table (monomials, coeffs)" % name, ok_lf)
    check("%s: LF(F1), LF(F2) single monomials divisible by x*y*z" % name,
          all(all(e >= 1 for e in lfs[i][0][0]) for i in (0, 1)))
    check("%s: strictly decreasing component degrees" % name,
          degs[0] > degs[1] > degs[2])
    # geometric degree claim of the paper: n + 2 (not re-derived; recorded)
    print("      (paper: deg_geom = n + 2 = %d; taken from paper, not re-derived)"
          % (n + 2))

# Compare expanded polynomials against the builder's construction.
spec = importlib.util.spec_from_file_location(
    "run_obstruction", os.path.join(HERE, "run_obstruction.py"))
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
for name in MEMBERS:
    Eb, hb, _ = mod.MEMBERS[name]
    Fb = mod.family(Eb, hb)
    same = all(sp.expand(built[name][i] - Fb[i]) == 0 for i in range(3))
    check("%s: my reconstruction == builder's family() (expanded, all 3 comps)"
          % name, same)
    # and that the builder's (E, h) is the paper's (E, h)
    E, h, _ = MEMBERS[name]
    check("%s: builder's (E, h) == paper's (E, h)" % name,
          sp.expand(Eb - E.subs(NU, mod.nu)) == 0
          and sp.expand(hb - h.subs(T, mod.t)) == 0)

# ===========================================================================
# V2. Structural facts, re-derived.
# ===========================================================================
print("== V2a: chart completeness (own RREF over Fractions) ==")

def rref_frac(M):
    """RREF of a list-of-lists of Fractions; returns (R, pivot_columns)."""
    R = [row[:] for row in M]
    nr, nc = len(R), len(R[0])
    piv_cols, r = [], 0
    for c in range(nc):
        piv = next((i for i in range(r, nr) if R[i][c] != 0), None)
        if piv is None:
            continue
        R[r], R[piv] = R[piv], R[r]
        R[r] = [v / R[r][c] for v in R[r]]
        for i in range(nr):
            if i != r and R[i][c] != 0:
                f = R[i][c]
                R[i] = [a - f * b for a, b in zip(R[i], R[r])]
        piv_cols.append(c)
        r += 1
        if r == nr:
            break
    return R, piv_cols

def classify_chart(R, piv):
    """Return chart name if R matches one of the three patterns, else None."""
    if piv == [0, 1] and R[0][0] == 1 and R[0][1] == 0 and R[1][0] == 0 and R[1][1] == 1:
        return "C12"   # [[1,0,p],[0,1,q]]
    if piv == [0, 2] and R[0][0] == 1 and R[0][2] == 0 and R[1] == [0, 0, 1]:
        return "C13"   # [[1,p,0],[0,0,1]]
    if piv == [1, 2] and R[0] == [0, 1, 0] and R[1] == [0, 0, 1]:
        return "C23"   # [[0,1,0],[0,0,1]]
    return None

rng = random.Random(987654321)
trials = [[[Fraction(rng.randint(-6, 6)) for _ in range(3)] for _ in range(2)]
          for _ in range(600)]
trials += [  # adversarial shapes
    [[Fraction(0), Fraction(1), Fraction(0)], [Fraction(0), Fraction(0), Fraction(1)]],
    [[Fraction(0), Fraction(2), Fraction(3)], [Fraction(0), Fraction(5), Fraction(7)]],
    [[Fraction(1), Fraction(2), Fraction(3)], [Fraction(2), Fraction(4), Fraction(7)]],
    [[Fraction(1), Fraction(2), Fraction(3)], [Fraction(2), Fraction(4), Fraction(6)]],
    [[Fraction(0), Fraction(0), Fraction(5)], [Fraction(0), Fraction(3), Fraction(1)]],
]
seen, nrank2, allok = set(), 0, True
for M in trials:
    R, piv = rref_frac(M)
    if len(piv) < 2:
        continue
    nrank2 += 1
    name = classify_chart(R, piv)
    if name is None:
        allok = False
        print("    UNCLASSIFIED RREF:", R)
        continue
    seen.add(name)
    # gauge factor: A = M restricted to pivot columns; check M == A*R and det != 0
    A = [[M[i][piv[0]], M[i][piv[1]]] for i in range(2)]
    detA = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    prod = [[sum(A[i][k] * R[k][j] for k in range(2)) for j in range(3)]
            for i in range(2)]
    if detA == 0 or prod != M:
        allok = False
check("all %d rank-2 samples have RREF in {C12, C13, C23}, with exact gauge "
      "factor M = A*R, det A != 0" % nrank2, allok and nrank2 > 500)
check("all three chart types realized in sample", seen == {"C12", "C13", "C23"})
# pivot-set argument: a rank-2 RREF has pivot columns {1,2}, {1,3} or {2,3};
# the three patterns above are the general forms -> the decomposition is COMPLETE.

# kernel formulas, symbolically
kerok = True
for rows, v in [(((1, 0, p), (0, 1, q)), (-p, -q, 1)),
                (((1, p, 0), (0, 0, 1)), (-p, 1, 0)),
                (((0, 1, 0), (0, 0, 1)), (1, 0, 0))]:
    for r in rows:
        kerok = kerok and sp.expand(r[0] * v[0] + r[1] * v[1] + r[2] * v[2]) == 0
check("chart kernel vectors annihilate the chart rows (symbolic, all charts)", kerok)

print("== V2b: Lemma 1 spot checks on random affine pi ==")
rng2 = random.Random(24681357)
lem_ok = True
for trial in range(5):
    while True:
        M = [[sp.Rational(rng2.randint(-4, 4)) for _ in range(3)] for _ in range(2)]
        v = (M[0][1] * M[1][2] - M[0][2] * M[1][1],
             M[0][2] * M[1][0] - M[0][0] * M[1][2],
             M[0][0] * M[1][1] - M[0][1] * M[1][0])  # row1 x row2 spans ker
        if any(c != 0 for c in v):
            break
    u1, u2 = rng2.randint(-3, 3), rng2.randint(-3, 3)
    l1 = M[0][0] * x + M[0][1] * y + M[0][2] * z + u1
    l2 = M[1][0] * x + M[1][1] * y + M[1][2] * z + u2
    G = sum(sp.Rational(rng2.randint(-3, 3)) * w1**a * w2**b
            for a in range(3) for b in range(3 - a))
    L = sp.expand(G.subs({w1: l1, w2: l2}, simultaneous=True))
    dv = lambda f: sp.expand(v[0] * sp.diff(f, x) + v[1] * sp.diff(f, y)
                             + v[2] * sp.diff(f, z))
    lem_ok = lem_ok and dv(L) == 0
    # converse direction checked WITHOUT the lemma: pick lam affine with
    # lam(v) = 1, invert the affine system, substitute, require no 3rd variable
    lam = None
    for cand in (x, y, z):
        cv = {x: v[0], y: v[1], z: v[2]}[cand]
        if cv != 0:
            lam = cand / cv
            break
    sol = sp.solve([sp.Eq(X, l1), sp.Eq(Y, l2), sp.Eq(Z, lam)], [x, y, z], dict=True)
    assert len(sol) == 1
    Lt = sp.expand(L.subs(sol[0], simultaneous=True))
    lem_ok = lem_ok and sp.diff(Lt, Z) == 0          # L in C[l1,l2] -> no Z
    Lbad = sp.expand(L + lam**3)
    Lbt = sp.expand(Lbad.subs(sol[0], simultaneous=True))
    lem_ok = lem_ok and dv(Lbad) != 0 and sp.diff(Lbt, Z) != 0  # both tests agree
check("Lemma 1 agrees with explicit affine inversion on 5 random pi "
      "(membership <-> no 3rd coordinate <-> d_v = 0)", lem_ok)

print("== V2c: existence of G for fixed pi is a linear system (exact QQ) ==")

def monoms_upto(D):
    return [(a, b) for s in range(D + 1) for a in range(s + 1) for b in [s - a]]

def gauss_solve_frac(Arows, brhs):
    """Solve A c = b over Fraction. Return (consistent, particular solution)."""
    m = len(Arows)
    n = len(Arows[0]) if m else 0
    Aug = [row[:] + [rv] for row, rv in zip(Arows, brhs)]
    r, piv_of_col = 0, {}
    for c in range(n):
        piv = next((i for i in range(r, m) if Aug[i][c] != 0), None)
        if piv is None:
            continue
        Aug[r], Aug[piv] = Aug[piv], Aug[r]
        Aug[r] = [vv / Aug[r][c] for vv in Aug[r]]
        for i in range(m):
            if i != r and Aug[i][c] != 0:
                f = Aug[i][c]
                Aug[i] = [aa - f * bb for aa, bb in zip(Aug[i], Aug[r])]
        piv_of_col[c] = r
        r += 1
    for i in range(r, m):
        if Aug[i][n] != 0:
            return False, None
    solvec = [Fraction(0)] * n
    for c, rr in piv_of_col.items():
        solvec[c] = Aug[rr][n]
    return True, solvec

def membership_system(L, l1, l2, D):
    """Rows/rhs over Fraction for L = sum_{a+b<=D} c_ab l1^a l2^b."""
    basis = monoms_upto(D)
    cols = []
    monset = set()
    for (a, b) in basis:
        Pe = sp.Poly(sp.expand(l1**a * l2**b), x, y, z)
        dcol = {m: c for m, c in Pe.terms()}
        monset |= set(dcol)
        cols.append(dcol)
    PL = sp.Poly(sp.expand(L), x, y, z)
    dL = {m: c for m, c in PL.terms()}
    monset |= set(dL)
    mons = sorted(monset)
    Arows = [[Fraction(int(sp.Rational(col.get(m, 0)).p),
                       int(sp.Rational(col.get(m, 0)).q)) for col in cols]
             for m in mons]
    brhs = [Fraction(int(sp.Rational(dL.get(m, 0)).p),
                     int(sp.Rational(dL.get(m, 0)).q)) for m in mons]
    return basis, Arows, brhs

# (i) randomly gauged control: pi' = A*(x, y+z) + b still admits a quotient.
CTRL1 = [x + (y + z)**2, y + x**3, z - x**3]
Ag, bg = [[2, 1], [1, 1]], [3, -2]
l1g = sp.expand(Ag[0][0] * x + Ag[0][1] * (y + z) + bg[0])
l2g = sp.expand(Ag[1][0] * x + Ag[1][1] * (y + z) + bg[1])
okg = True
Gg = []
for l in (l1g, l2g):
    LoF = sp.expand(l.subs({x: CTRL1[0], y: CTRL1[1], z: CTRL1[2]},
                           simultaneous=True))
    D = sp.Poly(LoF, x, y, z).total_degree()
    basis, Arows, brhs = membership_system(LoF, l1g, l2g, D)
    cons, sol = gauss_solve_frac(Arows, brhs)
    okg = okg and cons
    if cons:
        Gi = sum(sp.Rational(sol[k].numerator, sol[k].denominator)
                 * w1**a * w2**b for k, (a, b) in enumerate(basis))
        Gg.append(sp.expand(Gi))
        okg = okg and sp.expand(
            Gi.subs({w1: l1g, w2: l2g}, simultaneous=True) - LoF) == 0
jacGg = sp.expand(sp.diff(Gg[0], w1) * sp.diff(Gg[1], w2)
                  - sp.diff(Gg[0], w2) * sp.diff(Gg[1], w1)) if len(Gg) == 2 else None
check("randomly gauged CTRL1 pi: linear system solvable (own exact Gauss), "
      "pi o F = G o pi verified, Jac(G) = %s nonzero constant" % jacGg,
      okg and jacGg is not None and jacGg.is_constant() and jacGg != 0)

# (ii) F3 with a random translated rank-2 pi: system refuted exactly over QQ.
l1r = x + 2 * y + 3 * z + 1
l2r = 4 * x + 5 * y - 6 * z - 1
refuted = True
for l in (l1r, l2r):
    LoF = sp.expand(l.subs({x: built["F3"][0], y: built["F3"][1],
                            z: built["F3"][2]}, simultaneous=True))
    D = sp.Poly(LoF, x, y, z).total_degree()
    _, Arows, brhs = membership_system(LoF, l1r, l2r, D)
    cons, _ = gauss_solve_frac(Arows, brhs)
    if cons:
        refuted = False
        break  # first component membership must already fail for this pi
check("F3 with pi = (x+2y+3z+1, 4x+5y-6z-1): exact QQ linear system "
      "INCONSISTENT (no G of any degree <= deg(l o F))", refuted)

# ===========================================================================
# V3. Main independent certification: substitution route, all 9 cells.
# ===========================================================================
print("== V3: substitution-route certification (derivative-free), 9 cells ==")

def subst_generators(F, chart):
    """Generators of the membership-obstruction ideal, derivative-free.

    Complete (l1, l2) to affine coordinates, substitute the inverse, and
    collect every coefficient of a monomial containing the third coordinate.
    Returns (list of sympy exprs in Q[params], list of (which L, monomial)).
    """
    if chart == "C12":     # l1 = x + p z, l2 = y + q z; inverse x=X-pZ, y=Y-qZ, z=Z
        Ls = [sp.expand(F[0] + p * F[2]), sp.expand(F[1] + q * F[2])]
        rep = {x: X - p * Z, y: Y - q * Z, z: Z}
        bad = 2  # index of Z exponent
    elif chart == "C13":   # l1 = x + p y, l2 = z; inverse x=X-pY, y=Y, z=Z
        Ls = [sp.expand(F[0] + p * F[1]), F[2]]
        rep = {x: X - p * Y, y: Y, z: Z}
        bad = 1  # Y exponent
    else:                  # C23: l1 = y, l2 = z; inverse x=X (3rd coord), y=Y, z=Z
        Ls = [F[1], F[2]]
        rep = {x: X, y: Y, z: Z}
        bad = 0  # X exponent
    gens, where = [], []
    for i, L in enumerate(Ls, 1):
        Lt = sp.expand(L.xreplace(rep))
        for m, c in sp.Poly(Lt, X, Y, Z).terms():
            if m[bad] >= 1:
                ce = sp.expand(sp.sympify(c))
                if ce != 0:
                    gens.append(ce)
                    where.append((i, m))
    return gens, where

def frac_poly_from_expr(e):
    """sympy univariate poly in p -> Fraction coeff list (ascending)."""
    Pu = sp.Poly(e, p)
    n = Pu.degree()
    out = [Fraction(0)] * (n + 1)
    for (k,), c in Pu.terms():
        r = sp.Rational(c)
        out[k] = Fraction(int(r.p), int(r.q))
    return out

def upoly_gcd(a, b):
    """gcd of two Fraction coefficient lists (ascending), monic result."""
    def norm(c):
        while c and c[-1] == 0:
            c.pop()
        return c
    a, b = norm(a[:]), norm(b[:])
    while b:
        # a mod b
        while len(a) >= len(b) and a:
            f = a[-1] / b[-1]
            sh = len(a) - len(b)
            for i in range(len(b)):
                a[sh + i] -= f * b[i]
            a = norm(a)
        a, b = b, a
    if a:
        lead = a[-1]
        a = [c / lead for c in a]
    return a

table = {}
all_infeasible = True
for name, F in built.items():
    lfF1 = EXPECTED[name]["LF"][0]
    for chart in ("C12", "C13", "C23"):
        gens, where = subst_generators(F, chart)
        units = [(g, w) for g, w in zip(gens, where) if sp.sympify(g).is_Rational]
        rec = {"n_gen_occurrences": len(gens), "n_unit_occurrences": len(units)}
        if chart == "C23":
            infeas = any(g != 0 for g in gens)   # all gens are rational constants
            assert all(sp.sympify(g).is_Rational for g in gens)
            method = "direct inspection (all generators rational constants)"
        elif chart == "C13":
            gcd = None
            for g in gens:
                c = frac_poly_from_expr(g)
                gcd = c if gcd is None else upoly_gcd(gcd, c)
                if gcd == [Fraction(1)]:
                    break
            infeas = (gcd == [Fraction(1)])
            method = "own univariate gcd over Fraction == 1 (ideal = Q[p])"
        else:  # C12
            infeas = len(units) > 0
            method = "unit-constant generator (valid for all complex p, q)"
        rec.update(infeasible=infeas, method=method)
        table[(name, chart)] = rec
        all_infeasible = all_infeasible and infeas
        check("%s %s: INFEASIBLE via %s  (%d gen occurrences, %d unit occurrences)"
              % (name, chart, method, len(gens), len(units)), infeas)
        # leading-form prediction: for C12 the coefficient of X^a Y^b Z^c
        # (a,b,c) = LF(F1) exponents must be exactly the LF coefficient
        if chart == "C12":
            pred = [g for g, w in zip(gens, where)
                    if w == (1, lfF1[0])]
            check("%s C12: unit certificate at LF(F1) position %s equals %s "
                  "(from-scratch leading-form collapse)"
                  % (name, lfF1[0], lfF1[1]),
                  len(pred) == 1 and sp.expand(pred[0] - lfF1[1]) == 0)

check("ALL 9 (map, chart) cells infeasible by the substitution route",
      all_infeasible)

# near-miss recorded in the memo: F5, C13, p = 0: l2 o F = F5_3 in C[x, z]
nm_ok = (sp.diff(built["F5"][2], y) == 0
         and sp.diff(built["F3"][2], y) != 0
         and sp.diff(built["F4"][2], y) != 0
         and sp.diff(built["F5"][0], y) != 0)
gens5, _ = subst_generators(built["F5"], "C13")
# generators coming from L2 = F3 must all vanish at p=0; check via full recompute
g5_L2 = [g for g, w in zip(*subst_generators(built["F5"], "C13")) if w[0] == 2]
nm_ok = nm_ok and all(sp.expand(g.subs(p, 0)) == 0 for g in g5_L2)
check("near-miss confirmed: F5_3 = x - x^3 z is y-free (2nd eq satisfiable at "
      "p=0 in C13) but F5_1 involves y; F3_3, F4_3 are NOT y-free", nm_ok)

# ===========================================================================
# V4. Third flavour: exact point-evaluation (fiber-constancy) for C12; and
#     explicit integer-point contradictions for C23.
# ===========================================================================
print("== V4: point-evaluation route ==")
BASE_PTS = [(0, 0, 0), (1, 1, 1), (1, -1, 2), (2, 1, -1)]
for name, F in built.items():
    v = (-p, -q, 1)
    gens = []
    for (a0, b0, c0) in BASE_PTS:
        for L in (sp.expand(F[0] + p * F[2]), sp.expand(F[1] + q * F[2])):
            Dp = sp.expand(
                L.xreplace({x: a0 + v[0], y: b0 + v[1], z: c0 + v[2]})
                - L.xreplace({x: sp.Integer(a0), y: sp.Integer(b0),
                              z: sp.Integer(c0)}))
            if Dp != 0:
                gens.append(Dp)
    gb = sp.groebner(gens, p, q, order="grevlex")
    check("%s C12: Groebner basis of %d point-difference polynomials = [1] "
          "(no complex (p,q) makes l_i o F fiber-constant)" % (name, len(gens)),
          list(gb.exprs) == [sp.Integer(1)])

for name, F in built.items():
    # pi = (y, z): the points (0,0,-1) and (1,0,-1) share pi but split under F3
    v0 = F[2].xreplace({x: sp.Integer(0), y: sp.Integer(0), z: sp.Integer(-1)})
    v1 = F[2].xreplace({x: sp.Integer(1), y: sp.Integer(0), z: sp.Integer(-1)})
    check("%s C23: exact contradiction at integer points: F_3(0,0,-1) = %s != "
          "%s = F_3(1,0,-1) though pi(y,z) agrees" % (name, v0, v1), v0 != v1)

# ===========================================================================
# V5. Positive controls through THIS pipeline.
# ===========================================================================
print("== V5: positive controls through the substitution route ==")
CTRL2 = [x + (y + z)**2, (y + z)**3 - z + x**5, z - x**5]
for cname, Fc, expect_const in (("CTRL1", CTRL1, True), ("CTRL2", CTRL2, False)):
    found = []
    for chart in ("C12", "C13", "C23"):
        gens, where = subst_generators(Fc, chart)
        if chart == "C23":
            if all(g == 0 for g in gens):
                found.append((chart, {}))
            continue
        if chart == "C13":
            gcd = None
            for g in gens:
                c = frac_poly_from_expr(g)
                gcd = c if gcd is None else upoly_gcd(gcd, c)
            if gcd != [Fraction(1)]:
                # solve for roots of gcd in Q (here: none expected for controls)
                found.append((chart, {"gcd": gcd}))
            continue
        gb = sp.groebner(gens, p, q, order="grevlex")
        if list(gb.exprs) == [sp.Integer(1)]:
            continue
        sols = sp.solve(gens, [p, q], dict=True)
        for s in sols:
            found.append((chart, s))
    ok = len(found) == 1 and found[0][0] == "C12" and \
        found[0][1] == {p: 0, q: 1}
    check("%s: unique feasible chart/parameters = C12, (p,q) = (0,1)" % cname, ok)
    # reconstruct G by setting Z = 0 in the transformed polynomial
    s = found[0][1]
    Gs = []
    idok = True
    for L in (sp.expand(Fc[0] + s[p] * Fc[2]), sp.expand(Fc[1] + s[q] * Fc[2])):
        Lt = sp.expand(L.xreplace({x: X - s[p] * Z, y: Y - s[q] * Z, z: Z}))
        idok = idok and sp.diff(Lt, Z) == 0
        Gi = Lt.xreplace({X: w1, Y: w2, Z: sp.Integer(0)})
        Gs.append(sp.expand(Gi))
        l1c, l2c = x + s[p] * z, y + s[q] * z
        idok = idok and sp.expand(
            Gi.xreplace({w1: l1c, w2: l2c}) - L) == 0
    jacg = sp.expand(sp.diff(Gs[0], w1) * sp.diff(Gs[1], w2)
                     - sp.diff(Gs[0], w2) * sp.diff(Gs[1], w1))
    isconst = jacg.is_constant() and jacg != 0
    check("%s: G = %s reconstructed, pi o F = G o pi exact, Jac(G) = %s, "
          "nonzero-constant classification correct" % (cname, Gs, jacg),
          idok and isconst == expect_const)

# ===========================================================================
# V6. Audit of builder's artifacts (memo table vs certificates.json vs own
#     recomputation of the fiber-derivative route).
# ===========================================================================
print("== V6: audit builder's reported artifacts ==")
MEMO_TABLE = {  # (map, chart): (n_generators, n_unit_certificates)
    ("F3", "C12"): (23, 5), ("F3", "C13"): (18, 7), ("F3", "C23"): (6, 8),
    ("F4", "C12"): (61, 17), ("F4", "C13"): (48, 15), ("F4", "C23"): (16, 16),
    ("F5", "C12"): (71, 30), ("F5", "C13"): (59, 19), ("F5", "C23"): (18, 19),
}
MEMO_SAMPLES = {  # (map, chart): (equation index, (x,y,z) monomial, value)
    ("F3", "C12"): (1, (3, 3, 0), "2"), ("F3", "C13"): (1, (3, 2, 1), "6"),
    ("F3", "C23"): (2, (0, 0, 0), "1"),
    ("F4", "C12"): (1, (6, 4, 1), "3"), ("F4", "C13"): (1, (6, 3, 2), "6"),
    ("F4", "C23"): (2, (0, 0, 0), "1"),
    ("F5", "C12"): (1, (9, 5, 2), "-6"), ("F5", "C13"): (1, (9, 4, 3), "-10"),
    ("F5", "C23"): (2, (0, 0, 0), "1"),
}
CHART_DEFS = {  # rows, kernel v -- restated here, not imported
    "C12": ([[1, 0, p], [0, 1, q]], (-p, -q, 1)),
    "C13": ([[1, p, 0], [0, 0, 1]], (-p, 1, 0)),
    "C23": ([[0, 1, 0], [0, 0, 1]], (1, 0, 0)),
}
cert = json.load(open(os.path.join(HERE, "certificates.json")))
audit_ok = True
for name, F in built.items():
    for chart, (rows, v) in CHART_DEFS.items():
        gens_set, units = set(), []
        for i, r in enumerate(rows, 1):
            L = sp.expand(r[0] * F[0] + r[1] * F[1] + r[2] * F[2])
            Ei = sp.expand(v[0] * sp.diff(L, x) + v[1] * sp.diff(L, y)
                           + v[2] * sp.diff(L, z))
            if Ei == 0:
                continue
            for m, c in sp.Poly(Ei, x, y, z).terms():
                ce = sp.expand(sp.sympify(c))
                if ce != 0:
                    gens_set.add(ce)
                    if ce.is_Rational:
                        units.append((i, m, str(ce)))
        n_g, n_u = len(gens_set), len(units)
        # 1) memo table
        if (n_g, n_u) != MEMO_TABLE[(name, chart)]:
            audit_ok = False
            print("    MISMATCH memo table %s %s: got (%d,%d), memo %s"
                  % (name, chart, n_g, n_u, MEMO_TABLE[(name, chart)]))
        # 2) memo sample unit certificate present
        eqi, mono, val = MEMO_SAMPLES[(name, chart)]
        if (eqi, mono, val) not in units:
            audit_ok = False
            print("    MISSING memo sample certificate %s %s: %s"
                  % (name, chart, MEMO_SAMPLES[(name, chart)]))
        # 3) certificates.json agreement
        jc = cert["members"][name]["charts"][chart]
        if (jc["n_generators"] != n_g
                or len(jc["unit_certificates"]) != n_u
                or jc["groebner_basis"] != ["1"]
                or jc["infeasible"] is not True):
            audit_ok = False
            print("    MISMATCH certificates.json %s %s" % (name, chart))
        # 4) GB recomputed on MY generator set (grevlex AND lex)
        if chart != "C23":
            unk = (p, q) if chart == "C12" else (p,)
            for order in ("grevlex", "lex"):
                gb = sp.groebner(sorted(gens_set, key=sp.default_sort_key),
                                 *unk, order=order)
                if list(gb.exprs) != [sp.Integer(1)]:
                    audit_ok = False
                    print("    GB(%s) != [1] for %s %s" % (order, name, chart))
check("memo table, memo sample certificates, certificates.json and own "
      "recomputation of the fiber-derivative route all agree; GB = [1] "
      "in grevlex AND lex for every parametric cell", audit_ok)

# ===========================================================================
print()
el = time.time() - t0
if FAILURES:
    print("INDEPENDENT VERIFICATION FAILED (%d failures, %.1f s):" % (len(FAILURES), el))
    for f in FAILURES:
        print("  -", f)
    sys.exit(1)
print("INDEPENDENT VERIFICATION PASSED (%.1f s)" % el)
print("Outcome B re-certified by derivative-free substitution route (all 9 cells),")
print("point-evaluation route (C12 x 3 maps, plus integer-point contradictions for")
print("C23), exact-QQ linear-system refutation on a sampled pi, chart completeness")
print("re-derived, positive controls reproduced, and builder artifacts audited.")
