#!/usr/bin/env python3
"""
verify_band4_experiment.py
==========================
Exact SymPy verification backing `band4-moment-unit-experiment.md`: the decisive
QUANTUM band-4 moment-unit stress test at the composite cyclotomic wall.

CAMPAIGN CONTEXT.  DC1 (every endomorphism of A_1 is an automorphism) is attacked
via a band filtration; the candidate uniform mechanism is the MOMENT-UNIT
principle -- the central potential G with Q_0 = (T-1)G is band-agnostic (W4),
membership pins G(0)=0, so Q_0 = 1 <=> G = E is the slope statement, and in every
quantum-resistant branch so far the branch structure forces slope 0, making the
unit in [D,X]=1 unrealizable.  Band 2 (84978b9/b9f9cf3) and band 3 (99fe6ee..9fa9f74)
confirmed this.  This script runs the k=4 test, where the necklace factor
(S^4-1)/(S-1) = S_4 = 1+S+S^2+S^3 = (1+S)(1+S^2) = Phi_2 * Phi_4 is COMPOSITE
(vs k=3's irreducible Phi_3).

Conventions (frozen, identical to every sibling quantum memo):
    A_1[x^{-1}] = (+)_k x^k C[E],  (x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E),
    f^[r](E) = f(E+r),  E = x d,  ladder-m coeff of [D,X]:
    Q_m = sum_{k+l=m} ( b_l^[k] a_k - a_k^[l] b_l ),  [D,X]=1 <=> Q_m = delta_{m0},
      m in [-8,8].
    Genuine A_1 membership: E(E-1)...(E-r+1) | a_{-r}, b_{-r}, r=1..4.
    Gauge b_4 = 0 (spent on the top, Q_8).

SECTIONS
  0. engine: Q_m == direct crossed-product commutator (band 4, m in [-8,8]);
     gauge b_4=0 descent operators L_m; Q_0 = (T-1)G (both W4 and closure forms);
     G(0)=0 under membership; slope := const coeff of Q_0 = G(1).
  1. THE WALL Q_7 = b_3^[4] a_4 - a_4^[3] b_3 (the k=4 necklace lemma).  Necklace
     form S_4 B = S*S_3*A; S_4 = Phi_2*Phi_4 COMPOSITE; gcd(S_3,S_4)=1 forces the
     FULL S_4 | A (the 'either factor' hypothesis is FALSE); shifted-4th-power
     sufficiency; degree law 4 deg b_3 = 3 deg a_4; the exotic (non-4th-power)
     witnesses solve the wall.
  2. ENUMERATION of minimal (deg a_4 = 4, single-coset) exotic tops: complete scan
     -> exactly 5 realizable tops = 1 shifted-4th-power + 4 EXOTIC (3 reflection
     classes).  NO collapse to a single AP class (contrast k=3's unique AP).  The
     3-tetromino == 4-triomino tiling picture.  Cofactors incl. Phi_6, Phi_6*Phi_12.
  3. THE EXPERIMENT.  For each exotic top: positive cascade Q_6..Q_1 forward-solves
     with free lower data; {cascade} u {Q_0=1} is INFEASIBLE (Groebner=[1]) while
     {cascade} u {Q_0=0} is FEASIBLE -- the moment UNIT is the killer.  Verified at
     d=1 and d=2 for all four exotic tops, plus integer-r instances of the step-3
     AP family {r,r+3,r+6,r+9}.  (Each Groebner=[1] is an exact per-instance
     emptiness PROOF.  The uniform-in-r / uniform-in-d closure is a stated residual.)
  4. the b_3=0 sub-branch: the L_0 degree obstruction 4(4+q) (arbitrary degree),
     and the Groebner kill; honest note that unlike k=3 the collapse is incomplete.
  5. POSITIVE CONTROL: a genuine band-4 pair (U=x+d, X=U^4-d, D=U) is reproduced by
     the cascade and satisfies Q_0=1 -- the pipeline does NOT falsely kill genuine
     pairs; the exotic infeasibility is real.

Run:  uv run --with sympy python research/dc1-program/verify_band4_experiment.py
Ends: ALL BAND4 EXPERIMENT CHECKS PASSED
"""
import sympy as sp
from itertools import combinations

E, S, r = sp.symbols('E S r')
K = 4


# ---------------------------------------------------------------- primitives
def sh(f, s):
    return sp.expand(sp.sympify(f).subs(E, E + s))


def poly(name, deg):
    if deg < 0:
        return sp.Integer(0), []
    cs = list(sp.symbols(f'{name}_0:{deg+1}'))
    return sp.expand(sum(cs[i] * E**i for i in range(deg + 1))), cs


def falling(k):
    return sp.prod([E - i for i in range(k)]) if k > 0 else sp.Integer(1)


def Qm(a, b, m):
    return sp.expand(sum(
        sh(b[l], k) * a[k] - sh(a[k], l) * b[l]
        for k in range(-K, K + 1) for l in range(-K, K + 1) if k + l == m))


def gauged():
    return {k: sp.Integer(0) for k in range(-K, K + 1)}


def az(expr, label):
    e = sp.expand(sp.sympify(expr))
    if e != 0:
        raise AssertionError(label + "  :  residual = " + str(e))
    print("PASS", label)


def istrue(cond, label):
    if not cond:
        raise AssertionError(label + "  :  FALSE")
    print("PASS", label)


# -------------------------------------------------- cascade solver (K=4)
def clean_solve(A, B, m, lkey, name, memb, target, bdeg_raw, rsym=None):
    """Solve Q_m = target for B[lkey] = falling(memb)*generic(deg bdeg_raw).
    Left-nullspace of the operator matrix gives solvability conditions; an
    independent-row solve gives a particular solution; genuine kernels are
    re-added as named free params.  rsym (if given) is the single allowed
    symbolic parameter (the AP shift r)."""
    braw, bc = poly(f'{name}c_', bdeg_raw)
    bnew = sp.expand(falling(memb) * braw)
    Bt = dict(B); Bt[lkey] = bnew
    eq = sp.expand(Qm(A, Bt, m) - target)
    coeffs = sp.Poly(eq, E).all_coeffs() if eq != 0 else [sp.Integer(0)]
    M, rhs = sp.linear_eq_to_matrix(coeffs, bc)
    allowed = set() if rsym is None else {rsym}
    if not all((e.free_symbols <= allowed) for e in M):
        raise AssertionError("operator matrix not numeric/r-parametric (bilinearity leaked)")
    conds = [sp.expand(n.dot(rhs)) for n in M.T.nullspace()]
    conds = [c for c in conds if c != 0]
    Mred = sp.zeros(0, len(bc)); rhs_red = []
    for i in range(M.shape[0]):
        cand = Mred.col_join(M[i, :])
        if cand.rank() > Mred.rank():
            Mred = cand; rhs_red.append(rhs[i])
    if Mred.shape[0] == 0:
        bc_vals = [sp.Integer(0)] * len(bc)
    else:
        sol, _ = Mred.gauss_jordan_solve(sp.Matrix(rhs_red))
        tau = [s for s in sol.free_symbols if str(s).startswith('tau')]
        bc_vals = [x.subs({t: 0 for t in tau}) for x in sol]
    bsol = sp.expand(bnew.subs(dict(zip(bc, bc_vals))))
    ker = []
    for idx, kv in enumerate(M.nullspace()):
        g = sp.symbols(f'{name}K{idx}')
        kp = sp.expand(falling(memb) * sum(kv[i] * E**i for i in range(len(bc))))
        bsol = sp.expand(bsol + g * kp); ker.append(g)
    return bsol, ker, conds


def build_cascade(a4, b3, d, rsym=None):
    """Forward-solve the positive cascade Q_6,Q_5,Q_4,Q_3,Q_2,Q_1 for
    b_2,b_1,b_0,b_-1,b_-2,b_-3 with free lower data of degree d (membership on
    negatives) and b_-4 = mu4 * a_-4.  Returns (A,B,pos_conditions,vars,cs)."""
    a3, ca3 = poly('a3', d); a2, ca2 = poly('a2', d); a1, ca1 = poly('a1', d); a0, ca0 = poly('a0', d)
    am1t, cam1 = poly('am1', d); am2t, cam2 = poly('am2', d)
    am3t, cam3 = poly('am3', d); am4t, cam4 = poly('am4', d)
    am1 = sp.expand(falling(1) * am1t); am2 = sp.expand(falling(2) * am2t)
    am3 = sp.expand(falling(3) * am3t); am4 = sp.expand(falling(4) * am4t)
    mu4 = sp.symbols('mu4')
    A = {4: a4, 3: a3, 2: a2, 1: a1, 0: a0, -1: am1, -2: am2, -3: am3, -4: am4}
    B = gauged(); B[3] = b3; B[-4] = sp.expand(mu4 * am4)
    pos = []; ks = []
    for (m, l, mem) in [(6, 2, 0), (5, 1, 0), (4, 0, 0), (3, -1, 1), (2, -2, 2), (1, -3, 3)]:
        bb, kk, cc = clean_solve(A, B, m, l, f'b{l}', mem, 0, 3 * d + 6, rsym=rsym)
        B[l] = bb; pos += cc; ks += kk
    allvars = (list(ca3) + list(ca2) + list(ca1) + list(ca0)
               + list(cam1) + list(cam2) + list(cam3) + list(cam4) + [mu4] + ks)
    cs = dict(ca3=ca3, ca2=ca2, ca1=ca1, ca0=ca0, cam1=cam1, cam2=cam2, cam3=cam3, cam4=cam4, mu4=mu4)
    return A, B, pos, allvars, cs


def q0_conditions(A, B, target):
    ex = sp.expand(Qm(A, B, 0) - target)
    return [sp.expand(c) for c in sp.Poly(ex, E).all_coeffs() if sp.expand(c) != 0] if ex != 0 else []


def mk(ar, br):
    return sp.expand(sp.prod([E - x for x in ar])), sp.expand(sp.prod([E - x for x in br]))


# The four minimal exotic tops (deg a_4 = 4), and the shifted-4th-power control top.
EXOTIC = [('{0,2,3,5}', [0, 2, 3, 5], [1, 3, 5]),
          ('{0,1,3,6}', [0, 1, 3, 6], [1, 2, 6]),
          ('{0,3,5,6}', [0, 3, 5, 6], [1, 5, 6]),
          ('{0,3,6,9}', [0, 3, 6, 9], [1, 5, 9])]

S4 = 1 + S + S**2 + S**3
S3 = 1 + S + S**2
Phi2 = 1 + S
Phi4 = 1 + S**2


# =====================================================================
print("=" * 72)
print("0. engine: Q_m == commutator (band 4); Q_0=(T-1)G; G(0)=0; slope = G(1)")
print("=" * 72)
A0 = {k: poly(f'A{k+K}', 2)[0] for k in range(-K, K + 1)}
B0 = {k: poly(f'B{k+K}', 2)[0] for k in range(-K, K + 1)}


def direct_commutator(m):
    return sp.expand(sum(sh(B0[l], k) * A0[k] - sh(A0[k], l) * B0[l]
                         for k in range(-K, K + 1) for l in range(-K, K + 1) if k + l == m))


for m in range(-8, 9):
    az(direct_commutator(m) - Qm(A0, B0, m), f"Q_{m} = direct crossed-product commutator ladder coeff")

# gauge b_4=0 descent operator L_m[b] = b^[4] a_4 - a_4^[m-4] b (the (4,m-4) pair):
Bg = dict(B0); Bg[4] = sp.Integer(0)
for m in [7, 6, 5, 4, 3, 2, 1, 0]:
    op = sh(Bg[m - 4], 4) * A0[4] - sh(A0[4], m - 4) * Bg[m - 4]
    rest = sp.expand(sum(sh(Bg[l], k) * A0[k] - sh(A0[k], l) * Bg[l]
                         for k in range(-K, K + 1) for l in range(-K, K + 1)
                         if k + l == m and not (k == 4 and l == m - 4)))
    az(Qm(A0, Bg, m) - op - rest,
       f"Q_{m}|_(b4=0) = L_{m}[b_{m-4}] + rest,  L_{m}[b] = b^[4] a4 - a4^[{m-4}] b")

# Telescoping potential G, two equivalent closed forms (W4 and the closure form):
def potential_G(A, B):
    return sp.expand(sum(sh(A[k], j - k) * sh(B[-k], j) - sh(B[k], j - k) * sh(A[-k], j)
                         for k in range(1, K + 1) for j in range(0, k)))


def potential_G_W4(A, B):
    return sp.expand(sum(sh(A[i], -rr) * sh(B[-i], i - rr) - sh(A[-i], i - rr) * sh(B[i], -rr)
                         for i in range(1, K + 1) for rr in range(1, i + 1)))


az(Qm(A0, B0, 0) - (sh(potential_G(A0, B0), 1) - potential_G(A0, B0)),
   "Q_0 = (T-1) G   (telescoping lemma, generic band-4, closure form)")
az(potential_G(A0, B0) - potential_G_W4(A0, B0),
   "closure-form G == W4-form G = sum_i sum_{r=1}^i (a_i^[-r] b_-i^[i-r] - a_-i^[i-r] b_i^[-r])")

# G(0)=0 under membership; slope := const coeff of Q_0 = G(1):
Am = {k: poly(f'Am{k+K}', 2)[0] for k in range(0, K + 1)}
Bm = {k: poly(f'Bm{k+K}', 2)[0] for k in range(0, K + 1)}
for k in range(1, K + 1):
    Am[-k] = sp.expand(falling(k) * poly(f'amt{k}', 2)[0])
    Bm[-k] = sp.expand(falling(k) * poly(f'bmt{k}', 2)[0])
az(potential_G(Am, Bm).subs(E, 0), "G(0) = 0 identically under membership  =>  Q_0=1 <=> G=E")
Gm = potential_G(Am, Bm)
az(sp.Poly(sp.expand(sh(Gm, 1) - Gm), E).nth(0) - Gm.subs(E, 1),
   "slope := const coeff of Q_0 = G(1)  (the moment-unit slot)")
print("   => Q_0=1 forces G=E; the unit is realized iff the moment slope G(1) can be 1.")


# =====================================================================
print("\n" + "=" * 72)
print("1. THE WALL Q_7 (k=4 necklace lemma): S_4 B = S*S_3*A, S_4 = Phi_2*Phi_4")
print("=" * 72)
az(S4 - Phi2 * Phi4, "S_4 = 1+S+S^2+S^3 = (1+S)(1+S^2) = Phi_2 * Phi_4  (COMPOSITE)")
istrue(sp.factor_list(S4)[1] and len(sp.factor_list(S4)[1]) == 2,
       "S_4 factors into TWO distinct cyclotomics (composite, unlike k=3's irreducible Phi_3)")
az(sp.gcd(sp.Poly(S3, S), sp.Poly(S4, S)).as_expr() - 1,
   "gcd(S_3, S_4) = 1  =>  wall forces the FULL S_4 | A (not just one factor)")

# closed form of Q_7 in gauge b_4=0 (only (4,3) pair):
aq = {k: poly(f'aa{k+K}', 3)[0] for k in range(-K, K + 1)}
bq = gauged(); bq[3] = poly('bb3', 3)[0]
az(Qm(aq, bq, 7) - (sh(bq[3], 4) * aq[4] - sh(aq[4], 3) * bq[3]),
   "Q_7|_(b4=0) = b_3^[4] a_4 - a_4^[3] b_3   (the WALL)")

# shifted-4th-power SUFFICIENCY: a_4 = h h^[1] h^[2] h^[3]  =>  b_3 = h h^[1] h^[2] solves
a4c = sp.prod([sp.Function('h')(E + j) for j in range(4)])
b3c = sp.prod([sp.Function('h')(E + j) for j in range(3)])
az(sh(b3c, 4) * a4c - sh(a4c, 3) * b3c,
   "shifted-4th-power a_4=prod_0^3 h^[j], b_3=prod_0^2 h^[j]: solves the wall (symbolic h)")

# the four exotic witnesses solve the wall and are NOT shifted 4th powers:
for tag, ar, br in EXOTIC:
    a4, b3 = mk(ar, br)
    az(sh(b3, 4) * a4 - sh(a4, 3) * b3, f"{tag}: b_3 (roots {br}) solves the wall b_3^[4] a_4 = a_4^[3] b_3")
    # non-4th-power certificate: no linear h with a_4 = c h h^[1] h^[2] h^[3]
    hg, hgc = poly('hg', 1); c4 = sp.symbols('c4')
    sols = sp.solve(sp.Poly(sp.expand(c4 * hg * sh(hg, 1) * sh(hg, 2) * sh(hg, 3) - a4), E).all_coeffs(),
                    hgc + [c4], dict=True)
    istrue(len(sols) == 0, f"{tag}: NOT a shifted 4th power c h h^[1] h^[2] h^[3]  (EXOTIC)")

# degree law: staggered lead coeff(E^{p+q-1}) of Q_7 = (4 deg b_3 - 3 deg a_4) lc lc,
# so Q_7=0 forces 4 deg b_3 = 3 deg a_4  (== necklace sum of S_4 B = S S_3 A at S=1).
for (pa, qb) in [(4, 3), (8, 6), (4, 2), (5, 3)]:
    a4v, _ = poly('dA', pa); b3v, _ = poly('dB', qb)
    co = sp.Poly(sp.expand(sh(b3v, 4) * a4v - sh(a4v, 3) * b3v), E).coeff_monomial(E**(pa + qb - 1))
    la = sp.Poly(a4v, E).nth(pa); lb = sp.Poly(b3v, E).nth(qb)
    az(co - (4 * qb - 3 * pa) * la * lb,
       f"Q_7 staggered lead: coeff(E^{{p+q-1}}), p={pa},q={qb} = (4q-3p) lc lc = {4*qb-3*pa} lc lc")
print("   => Q_7=0 forces 4 deg b_3 = 3 deg a_4; minimal nonzero deg a_4 = 4, deg b_3 = 3.")


# =====================================================================
print("\n" + "=" * 72)
print("2. ENUMERATION: minimal (deg a_4=4, single-coset) exotic tops; tiling; no AP collapse")
print("=" * 72)


def wall_admissible(roots):
    """roots: integers (single mod-Z coset).  Returns (S4-divisible, b3-effective,
    non-4th-power, b3_roots) via the necklace formula B = S*S_3*A / S_4."""
    A = sum(S**rr for rr in roots)
    qA, rA = sp.div(sp.expand(A), S4, S)
    if sp.expand(rA) != 0:
        return (False, False, False, None)
    num = sp.expand(S * S3 * A)
    qB, rB = sp.div(num, S4, S)
    if sp.expand(rB) != 0:
        return (True, False, False, None)
    Bp = sp.Poly(qB, S)
    eff = all(Bp.nth(i) >= 0 for i in range(Bp.degree() + 1))
    Broots = [i for i in range(Bp.degree() + 1) for _ in range(int(Bp.nth(i)))] if eff else None
    Cp = sp.Poly(qA, S)
    noncube = any(Cp.nth(i) < 0 for i in range(Cp.degree() + 1))
    return (True, eff, noncube, Broots)


# complete scan of deg-4 single-coset tops {0,a,b,c}, span bounded (12-cell tiling => span<=9):
realizable = {}
for a, b, c in combinations(range(1, 16), 3):
    roots = [0, a, b, c]
    s4, eff, nc, Br = wall_admissible(roots)
    if s4 and eff:
        a4, b3 = mk(roots, Br)
        az(sh(b3, 4) * a4 - sh(a4, 3) * b3, f"scan {roots}: b_3 roots {Br} solves the wall")
        realizable[tuple(roots)] = ('EXOTIC' if nc else '4th-power', Br)
istrue(len(realizable) == 5,
       "deg-4 single-coset realizable tops (min root 0, span<=15): EXACTLY 5")
exo = [k for k in realizable if realizable[k][0] == 'EXOTIC']
istrue(sorted(exo) == [(0, 1, 3, 6), (0, 2, 3, 5), (0, 3, 5, 6), (0, 3, 6, 9)],
       "the 4 EXOTIC deg-4 tops (up to translation) = {0,1,3,6},{0,2,3,5},{0,3,5,6},{0,3,6,9}")
istrue(realizable[(0, 1, 2, 3)][0] == '4th-power',
       "the lone non-exotic realizable top {0,1,2,3} is the consecutive shifted 4th power")


def reflect(roots):
    M = max(roots); return tuple(sorted(M - x for x in roots))


classes = []
used = set()
for k in sorted(exo):
    if k in used:
        continue
    cls = {k, reflect(k)}; used |= cls; classes.append(sorted(cls))
istrue(len(classes) == 3,
       "the 4 exotic tops form 3 reflection classes: {0,2,3,5}(palindromic), "
       "{0,1,3,6}~{0,3,5,6}(pair), {0,3,6,9}(palindromic step-3 AP)")

# NO collapse to a single AP class (contrast k=3, where the exotic branch was exactly
# the step-2 AP {r,r+2,r+4}): here 3 of the 4 exotic tops are NON-arithmetic-progression.
non_ap = [k for k in exo if not all(sorted(k)[i + 1] - sorted(k)[i] == (sorted(k)[1] - sorted(k)[0])
                                    for i in range(len(k) - 1))]
istrue(sorted(non_ap) == [(0, 1, 3, 6), (0, 2, 3, 5), (0, 3, 5, 6)],
       "3 of 4 exotic tops are NON-AP  =>  the k=3 AP/tiling collapse FAILS at k=4")

# tiling picture: S_4*b3_roots (3 tetrominoes) == S*S_3*(a4_roots) (4 triominoes), as multisets:
for k in sorted(exo):
    Br = realizable[k][1]
    tet = sorted([pp + j for pp in Br for j in range(4)])          # S_4 * B
    tri = sorted([(rr + 1) + j for rr in k for j in range(3)])     # S * S_3 * A
    istrue(tet == tri, f"tiling {list(k)}: 3 tetrominoes (S_4 b_3) == 4 triominoes (S S_3 a_4) as multisets")

# the cofactors C = A / S_4 (exotic = C non-effective): note Phi_6 and Phi_6*Phi_12 appear:
for k, expected in [((0, 2, 3, 5), sp.cyclotomic_poly(6, S)),
                    ((0, 3, 6, 9), sp.expand(sp.cyclotomic_poly(6, S) * sp.cyclotomic_poly(12, S)))]:
    A = sum(S**rr for rr in k); C, _ = sp.div(sp.expand(A), S4, S)
    az(sp.expand(C - expected), f"cofactor of {list(k)} = {'Phi_6' if k==(0,2,3,5) else 'Phi_6*Phi_12'} "
                                f"(non-effective => exotic; = k=3 exotic cofactor Phi_6 reused)")
print("   => composite S_4=Phi_2 Phi_4 opens a WIDER effective exotic class than k=3's single AP family.")


# =====================================================================
print("\n" + "=" * 72)
print("3. THE EXPERIMENT: cascade Q_6..Q_1 + moment-unit test Q_0=1 vs Q_0=0")
print("=" * 72)
# 3a. explicit positive-cascade solvability witness for {0,2,3,5} (a concrete feasible point):
a4w, b3w = mk([0, 2, 3, 5], [1, 3, 5])
A, B, pos, allvars, cs = build_cascade(a4w, b3w, 1)
istrue(len(pos) >= 1, "{0,2,3,5} d=1: positive cascade Q_6..Q_1 forward-solves (with solvability conds)")
Gh = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(0)), *allvars, order='grevlex'))
istrue(Gh != [sp.Integer(1)], "{0,2,3,5} d=1: positive cascade is SOLVABLE (Q_0=0 branch feasible)")

# 3b. THE KILL: for every exotic top, {cascade} u {Q_0=1} infeasible, {Q_0=0} feasible.  d=1 and d=2.
for d in (1, 2):
    for tag, ar, br in EXOTIC:
        a4, b3 = mk(ar, br)
        A, B, pos, allvars, cs = build_cascade(a4, b3, d)
        Gu = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(1)), *allvars, order='grevlex'))
        Gh = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(0)), *allvars, order='grevlex'))
        istrue(Gu == [sp.Integer(1)] and Gh != [sp.Integer(1)],
               f"{tag} d={d}: {{cascade}}u{{Q_0=1}} INFEASIBLE, {{Q_0=0}} feasible => moment UNIT kills it")
print("   => EVERY minimal exotic top is KILLED at Q_0; the moment unit is unrealizable.")
print("      NO exotic top admits Q_0=1  =>  no DC1 counterexample candidate at band-4 minimal degree.")

# 3c. integer-r instances of the step-3 AP family {r,r+3,r+6,r+9} (translation robustness):
for rr in (0, 1, -1, 2):
    ar = [rr, rr + 3, rr + 6, rr + 9]; br = [rr + 1, rr + 5, rr + 9]
    a4, b3 = mk(ar, br)
    az(sh(b3, 4) * a4 - sh(a4, 3) * b3, f"step-3 AP r={rr}: b_3 solves the wall")
    A, B, pos, allvars, cs = build_cascade(a4, b3, 1)
    Gu = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(1)), *allvars, order='grevlex'))
    Gh = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(0)), *allvars, order='grevlex'))
    istrue(Gu == [sp.Integer(1)] and Gh != [sp.Integer(1)],
           f"step-3 AP r={rr} d=1: {{cascade}}u{{Q_0=1}} INFEASIBLE, {{Q_0=0}} feasible => killed")


# =====================================================================
print("\n" + "=" * 72)
print("4. the b_3=0 sub-branch: L_0 degree obstruction + Groebner kill")
print("=" * 72)
# L_0[b] = b^[4] a_4 - a_4^[-4] b, deg a_4 = 4.  The E^{q+4} terms cancel; the next
# coefficient E^{q+3} is 4(4+q) lc(a_4) lc(b) != 0 (char 0), so L_0[b] has degree q+3 >= 3
# (or is 0 when b=0) -- it can NEVER equal a nonzero constant.  [PROVED, arbitrary deg b.]
a4gen, _ = poly('LA', 4)
for q in range(0, 7):
    b, bc = poly('Lb', q)
    L0 = sp.expand(sh(b, 4) * a4gen - sh(a4gen, -4) * b)
    deg = sp.Poly(L0, E).degree(); lead = sp.Poly(L0, E).nth(q + 3)
    la = sp.Poly(a4gen, E).nth(4); lb = sp.Poly(b, E).nth(q)
    istrue(deg == q + 3 and sp.expand(lead - 4 * (4 + q) * la * lb) == 0,
           f"L_0 obstruction: deg b={q} => deg L_0 = {q+3} >= 3, lead = 4(4+q) lc(a4) lc(b) != 0")
print("   => on a fully-collapsed b_3=0 branch (b_2=..=b_-3=0, b_0=const), Q_0 = L_0[b_-4] = 1")
print("      is IMPOSSIBLE (deg >= 3 or 0, never a nonzero const).  [PROVED, arbitrary degree]")
# HONEST caveat: unlike k=3, the b_3=0 collapse is INCOMPLETE at k=4 -- L_6 (the (4,2) pair,
# shift +2) has a NONtrivial kernel for the exotic a_4 (e.g. {0,2,3,5} = (1+S^2)(1+S^3) has
# a genuine b_2 wall-solution), so b_2 is not forced to 0.  We therefore also verify the kill
# directly by Groebner on the FULL b_3=0 cascade for each exotic top:
az(sum(S**rr for rr in [0, 2, 3, 5]) - (1 + S**2) * (1 + S**3),
   "note: exotic a_4={0,2,3,5} root data A = (1+S^2)(1+S^3) => L_6 (shift +2) has a b_2 kernel")
# L_6[b]=0 necklace: S_4 B = S^2(1+S) A; for {0,2,3,5} this gives a genuine nonzero b_2 (roots {2,5}):
b2ker = sp.expand((E - 2) * (E - 5)); a4_235 = sp.expand(E * (E - 2) * (E - 3) * (E - 5))
az(sh(b2ker, 4) * a4_235 - sh(a4_235, 2) * b2ker,
   "L_6 kernel witness: b_2=(E-2)(E-5) solves b_2^[4]a_4 = a_4^[2]b_2 (b_3=0 collapse is INCOMPLETE)")
for tag, ar, br in EXOTIC:
    a4, _ = mk(ar, br)
    A, B, pos, allvars, cs = build_cascade(a4, sp.Integer(0), 1)
    Gu = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(1)), *allvars, order='grevlex'))
    Gh = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(0)), *allvars, order='grevlex'))
    istrue(Gu == [sp.Integer(1)] and Gh != [sp.Integer(1)],
           f"{tag} b_3=0 d=1: {{cascade}}u{{Q_0=1}} INFEASIBLE, {{Q_0=0}} feasible => killed at Q_0 too")


# =====================================================================
print("\n" + "=" * 72)
print("5. POSITIVE CONTROL: a genuine band-4 pair is reproduced; Q_0=1 HOLDS (no false kill)")
print("=" * 72)
def mul(P, Q):
    R = {}
    for k1, p1 in P.items():
        for k2, p2 in Q.items():
            R[k1 + k2] = R.get(k1 + k2, 0) + sh(p1, k2) * p2
    return {k: sp.expand(v) for k, v in R.items() if sp.expand(v) != 0}


# U = x + d,  X = U^4 - d,  D = U ;  d = x^{-1} E ;  [D,X] = -[U,d] = 1.  a_4 = 1 (trivial 4th pow).
U = {1: sp.Integer(1), -1: sp.expand(E)}
U4 = mul(mul(U, U), mul(U, U))
Xc = gauged()
for k, v in U4.items():
    Xc[k] = Xc.get(k, 0) + v
Xc[-1] = sp.expand(Xc[-1] - E)              # subtract d = x^{-1}E
Dc = gauged(); Dc[1] = sp.Integer(1); Dc[-1] = sp.expand(E)
istrue(all(sp.expand(Qm(Xc, Dc, m) - (1 if m == 0 else 0)) == 0 for m in range(-8, 9)),
       "control U=x+d, X=U^4-d, D=U : [D,X] = 1 (all Q_m = delta_{m0}), a_4 = 1 trivial 4th power")
istrue(sp.rem(sp.Poly(Xc[-4], E), sp.Poly(falling(4), E)) == 0,
       "control membership: E(E-1)(E-2)(E-3) | a_-4")
# feed the control's X into the SAME cascade solver; reconstruct D, confirm Q_0=1 holds on it.
Av = dict(Xc); Bv = gauged(); Bv[-4] = Dc[-4]
kern = []
for (m, l, mem, bd) in [(6, 2, 0, 6), (5, 1, 0, 6), (4, 0, 0, 6), (3, -1, 1, 7), (2, -2, 2, 8), (1, -3, 3, 9)]:
    bb, kk, cc = clean_solve(Av, Bv, m, l, f'vb{l}', mem, 0, bd)
    Bv[l] = bb; kern += kk
    # every emitted condition is a constraint on the b-KERNEL freedom, not a false a-obstruction:
    istrue(all(set(sp.sympify(c).free_symbols) <= set(kern) for c in cc),
           f"control Q_{m}: emitted conditions involve only b-kernel params (no false a-obstruction)")
choose = {s: (1 if 'vb1K' in str(s) else 0) for s in kern}
Brec = gauged(); Brec[-4] = Dc[-4]
for l in range(-4, 5):
    Brec[l] = sp.expand(sp.sympify(Bv[l]).subs(choose))
istrue(all(sp.expand(Brec[l] - Dc[l]) == 0 for l in range(-4, 5)),
       "control: cascade reconstructs D exactly (b_1=1, b_-1=E, rest 0) => feasibility detected")
istrue(all(sp.expand(Qm(Av, Brec, m) - (1 if m == 0 else 0)) == 0 for m in range(-8, 9)),
       "control: Q_0 = 1 HOLDS on the reconstructed genuine pair => pipeline does NOT falsely kill")


print("\nALL BAND4 EXPERIMENT CHECKS PASSED")
