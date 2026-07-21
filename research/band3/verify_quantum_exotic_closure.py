#!/usr/bin/env python3
"""
verify_quantum_exotic_closure.py
================================
Exact SymPy verification backing `quantum-exotic-closure.md`: the DEGREE-FREE
closure of the QUANTUM band-3 exotic (non-shifted-cube) `b_2 != 0` sub-case --
the one residual gap left by `quantum-exotic-branch.md` (commit ebfc64d), and the
load-bearing input the quantum band-3 floor routes through (`astar-band3.md` sends
quantum A* into the exotic non-cube question).

Conventions (frozen, identical to every sibling quantum memo, commit ebfc64d):
    A_1[x^{-1}] = (+)_k x^k C[E],  (x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E),
    f^[r](E) = f(E+r),  E = x d,  ladder-m coeff of [D,X]:
    Q_m = sum_{k+l=m} ( b_l^[k] a_k - a_k^[l] b_l ),  [D,X]=1 <=> Q_m = delta_{m0}.
    Genuine A_1 membership: E(E-1)...(E-r+1) | a_{-r}, b_{-r}.
    Gauge b_3 = 0 (spent on the top, Q_6).  Exotic top: a_3 solves the Q_5 wall
    b_2^[3] a_3 = a_3^[2] b_2 with b_2 != 0 and a_3 NOT a shifted cube.

WHAT THIS SCRIPT ESTABLISHES (see the memo for the argument):
  0. engine: Q_m == direct commutator; Q_0 = (T-1)G telescoping; G(0)=0 under
     membership; hence Q_0 = 1 <=> G = E, and the "moment unit" is the SLOPE of G,
     slope = G(1) = constant coefficient of Q_0.
  1. STRUCTURAL: at deg a_3 = 3 the realizable exotic tops are EXACTLY the step-2
     arithmetic-progression class {r,r+2,r+4} (every other Phi_3-divisible non-cube
     multiset {0,a,b} fails b_2-effectivity, i.e. is not in the b_2!=0 branch at
     all).  So the whole degree-3 exotic branch IS the AP class.  [DEGREE-6 is
     richer: non-AP realizable exotic tops exist -- flagged as residual.]
  2. MEMBERSHIP PROTECTION: the level-3 (bottom) contribution to the slope is
     mu_3 * a_3(0) * a_{-3}(3); vanishes identically when 0 is a root of a_3.
  3. TROPICAL SKELETON: the k=3 and k=2 block leading coefficients (Lemma-R
     staggered identities, symbolic in degree) force the top filler coefficients
     to vanish -- the top of the annihilation cascade.
  4. MAIN THEOREM (moment carries no unit): for the AP class {r,r+2,r+4} with r
     SYMBOLIC and d=1 free data, the slope of G is forced to 0 -- hence Q_0 = 1 is
     impossible.  This closes the ENTIRE degree-3 exotic branch at d=1.  The
     explicit W1 triangular filler-annihilation is walked through.
  5. VERIFICATION across the class: W1, W2, AP r in {1,-1,2,3} at d=1 (Q_0=1
     Groebner-infeasible, Q_0=0 feasible); W1,W2 at d=2 with the exact
     {slope forced 0} certificate; arbitrary top-degree AP g=1,2,3 (deg a_3=3,6,9).
  6. positive control: the pipeline reproduces a genuine b_2=0 pair, no spurious
     conditions -- the infeasibility above is a real kill, not an artifact.

Run:  uv run --with sympy python research/band3/verify_quantum_exotic_closure.py
Ends: ALL QUANTUM EXOTIC CLOSURE CHECKS PASSED
"""
import sympy as sp
from itertools import combinations

E, S = sp.symbols('E S')


# ---------------------------------------------------------------- primitives
def sh(f, r):
    return sp.expand(sp.sympify(f).subs(E, E + r))


def poly(name, deg):
    if deg < 0:
        return sp.Integer(0), []
    cs = list(sp.symbols(f'{name}_0:{deg+1}'))
    return sp.expand(sum(cs[i] * E**i for i in range(deg + 1))), cs


def falling(r):
    return sp.prod([E - i for i in range(r)]) if r > 0 else sp.Integer(1)


def Qm(a, b, m, K=3):
    return sp.expand(sum(
        sh(b[l], k) * a[k] - sh(a[k], l) * b[l]
        for k in range(-K, K + 1) for l in range(-K, K + 1) if k + l == m))


def gauged():
    return {k: sp.Integer(0) for k in range(-3, 4)}


def az(expr, label):
    e = sp.expand(sp.sympify(expr))
    if e != 0:
        raise AssertionError(label + "  :  residual = " + str(e))
    print("PASS", label)


def istrue(cond, label):
    if not cond:
        raise AssertionError(label + "  :  FALSE")
    print("PASS", label)


def clean_solve(A, B, m, lkey, name, memb, target, bdeg_raw, rsym=None):
    """Solve Q_m = target for B[lkey] = falling(memb)*generic(deg bdeg_raw).
    Left-nullspace of the (numeric-or-r-parametric) operator matrix gives the
    solvability conditions; an independent-row solve gives a particular solution;
    genuine kernels are re-added as named free params.  `rsym`, if given, is the
    (single) allowed symbolic parameter in the operator matrix (the AP shift r)."""
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


def build_cascade(a3, b2, d, rsym=None):
    """Forward-solve the positive cascade Q_4,Q_3,Q_2,Q_1 with free lower data of
    degree d and membership.  Returns (A, B, pos_conditions, all_free_vars)."""
    a2, ca2 = poly('a2', d); a1, ca1 = poly('a1', d); a0, ca0 = poly('a0', d)
    am1t, cam1 = poly('am1', d); am2t, cam2 = poly('am2', d); am3t, cam3 = poly('am3', d)
    am1 = sp.expand(falling(1) * am1t); am2 = sp.expand(falling(2) * am2t); am3 = sp.expand(falling(3) * am3t)
    mu3 = sp.symbols('mu3')
    A = {3: a3, 2: a2, 1: a1, 0: a0, -1: am1, -2: am2, -3: am3}
    B = gauged(); B[2] = b2; B[-3] = sp.expand(mu3 * am3)
    pos = []; ks = []
    for (m, l, nm, mem, bd) in [(4, 1, 'b1', 0, d + 3), (3, 0, 'b0', 0, 2 * d + 2),
                                (2, -1, 'bm1', 1, 2 * d + 3), (1, -2, 'bm2', 2, 2 * d + 4)]:
        bb, kk, cc = clean_solve(A, B, m, l, nm, mem, 0, bd, rsym=rsym)
        B[l] = bb; pos += cc; ks += kk
    allvars = list(ca2) + list(ca1) + list(ca0) + list(cam1) + list(cam2) + list(cam3) + [mu3] + ks
    return A, B, pos, allvars, dict(ca2=ca2, ca1=ca1, ca0=ca0, cam1=cam1, cam2=cam2, cam3=cam3, mu3=mu3)


def q0_conditions(A, B, target):
    ex = sp.expand(Qm(A, B, 0) - target)
    return [sp.expand(c) for c in sp.Poly(ex, E).all_coeffs() if sp.expand(c) != 0] if ex != 0 else []


def phi3_realizable(roots):
    """roots: integers (single mod-Z coset).  Returns (Phi3-divisible, b2-effective, non-cube)."""
    Amult = sum(S**r for r in roots)
    q, rem = sp.div(sp.expand(Amult), 1 + S + S**2, S)
    if sp.expand(rem) != 0:
        return (False, False, False)
    cof = sp.Poly(q, S); noncube = any(cof.nth(i) < 0 for i in range(cof.degree() + 1))
    qb, remb = sp.div(sp.expand(S * (1 + S) * Amult), 1 + S + S**2, S)
    Bp = sp.Poly(qb, S)
    eff = (sp.expand(remb) == 0) and all(Bp.nth(i) >= 0 for i in range(Bp.degree() + 1))
    return (True, eff, noncube)


def ap_top(roots):
    """Build (a_3, b_2) for a step-2 AP top from its integer root list."""
    a3 = sp.expand(sp.prod([E - r for r in roots]))
    Amult = sum(S**r for r in roots)
    qb, remb = sp.div(sp.expand(S * (1 + S) * Amult), 1 + S + S**2, S)
    assert sp.expand(remb) == 0
    Bp = sp.Poly(qb, S)
    broots = [i for i in range(Bp.degree() + 1) if Bp.nth(i) != 0]
    return a3, sp.expand(sp.prod([E - i for i in broots])), broots


# =====================================================================
# 0. Engine: Q_m == commutator; Q_0 = (T-1)G; G(0)=0; slope = G(1).
# =====================================================================
print("--- 0. engine: Q_m == commutator; Q_0=(T-1)G; G(0)=0; slope = G(1) ---")
A0 = {k: poly(f'A{k+3}', 2)[0] for k in range(-3, 4)}
B0 = {k: poly(f'B{k+3}', 2)[0] for k in range(-3, 4)}


def direct_commutator(m, K=3):
    return sp.expand(sum(sh(B0[l], k) * A0[k] - sh(A0[k], l) * B0[l]
                         for k in range(-K, K + 1) for l in range(-K, K + 1) if k + l == m))


for m in range(-6, 7):
    az(direct_commutator(m) - Qm(A0, B0, m), f"Q_{m} = direct commutator ladder coeff")

# Telescoping potential G (band-agnostic form, K=3), gauge b_3=0 not yet imposed:
def potential_G(A, B, K=3):
    return sp.expand(sum(
        sh(A[k], j - k) * sh(B[-k], j) - sh(B[k], j - k) * sh(A[-k], j)
        for k in range(1, K + 1) for j in range(0, k)))


az(Qm(A0, B0, 0) - (sh(potential_G(A0, B0), 1) - potential_G(A0, B0)),
   "Q_0 = (T-1)G  (telescoping lemma, generic band-3)")

# G(0) = 0 identically under membership (a_{-r},b_{-r} carry E(E-1)...(E-r+1)):
Am = {k: poly(f'Am{k+3}', 2)[0] for k in range(0, 4)}
Bm = {k: poly(f'Bm{k+3}', 2)[0] for k in range(0, 4)}
for k in range(1, 4):
    Am[-k] = sp.expand(falling(k) * poly(f'amt{k}', 2)[0])
    Bm[-k] = sp.expand(falling(k) * poly(f'bmt{k}', 2)[0])
az(potential_G(Am, Bm).subs(E, 0), "G(0) = 0 identically under membership  =>  Q_0=1 <=> G=E")
# slope of G (= G(1)-G(0) = G(1)) is the constant coefficient of Q_0:
Gm = potential_G(Am, Bm)
az(sp.Poly(sp.expand(sh(Gm, 1) - Gm), E).nth(0) - Gm.subs(E, 1),
   "slope := constant coeff of Q_0 = G(1)   (the 'moment unit' slot)")
print("   => Q_0=1 forces G=E; the unit is realized iff the moment slope G(1) can be 1.")


# =====================================================================
# 1. STRUCTURAL: degree-3 realizable exotic tops = step-2 AP class only.
# =====================================================================
print("\n--- 1. degree-3 realizable exotic tops = step-2 AP class {r,r+2,r+4} ---")
realizable3 = []
for a in range(1, 13):
    for b in range(a + 1, 13):
        phi, eff, nc = phi3_realizable([0, a, b])
        if phi and eff:
            realizable3.append(([0, a, b], 'exotic' if nc else 'cube'))
istrue(realizable3 == [([0, 1, 2], 'cube'), ([0, 2, 4], 'exotic')],
       "deg-3 tops {0,a,b}, a<b<=12: realizable = {0,1,2}(cube), {0,2,4}(exotic) ONLY")
istrue(sum(1 for r, k in realizable3 if k == 'exotic') == 1,
       "deg-3 realizable EXOTIC top is UNIQUE up to translation = {0,2,4} (step-2 AP)")
# Non-AP realizable exotic tops DO appear at degree 6 (honest residual flag):
nonap6 = []
for extra in combinations(range(1, 15), 5):
    roots = [0] + list(extra)
    phi, eff, nc = phi3_realizable(roots)
    if phi and eff and nc and not all(roots[i + 1] - roots[i] == 2 for i in range(len(roots) - 1)):
        nonap6.append(roots)
istrue(len(nonap6) > 0,
       f"deg-6 has {len(nonap6)} realizable NON-AP exotic tops (e.g. {nonap6[0]}) -- residual class")
print("   => at deg a_3 = 3 the exotic b_2!=0 branch IS the AP class; closing AP closes it.")


# =====================================================================
# 2. MEMBERSHIP PROTECTION of the bottom contribution to the slope.
# =====================================================================
print("\n--- 2. membership-protection: level-3 slope contribution = mu3*a3(0)*a_-3(3) ---")
a3g, _ = poly('cA', 3)
am3t, _ = poly('cB', 2)
am3 = sp.expand(falling(3) * am3t)
mu3 = sp.symbols('mu3')
bm3 = sp.expand(mu3 * am3)
# level-3 part of G in gauge b_3=0: P3 = sum_{j=0}^2 a3^[j-3] b_{-3}^[j]
P3 = sp.expand(sum(sh(a3g, j - 3) * sh(bm3, j) for j in range(3)))
az(P3.subs(E, 0), "P3(0) = 0   (all of a_-3(0),a_-3(1),a_-3(2) killed by membership)")
az(P3.subs(E, 1) - mu3 * a3g.subs(E, 0) * am3.subs(E, 3),
   "P3(1) = mu3 * a3(0) * a_-3(3)   (the ONLY surviving bottom term)")
# for an AP top with 0 a root, a3(0)=0, so the whole bottom drops out of the slope:
a3_r0 = sp.expand(E * (E - 2) * (E - 4))
az((sp.expand(sum(sh(a3_r0, j - 3) * sh(bm3, j) for j in range(3)))).subs(E, 1),
   "a3(0)=0 (roots ni 0): bottom contribution to slope vanishes identically")


# =====================================================================
# 3. TROPICAL SKELETON: block leading coefficients (Lemma-R), symbolic in degree.
# =====================================================================
print("\n--- 3. tropical skeleton: block leading coefficients force filler tops to vanish ---")
# k=3 block operator K3[c] = a3^[-3] c + a3^[-2] c^[1] + a3^[-1] c^[2]; deg = deg a3 + deg c,
# leading coeff = 3*lc(a3)*lc(c) != 0.  (T-1) keeps it top-order-1 nonzero => top of a_-3 forced 0.
for (pp, nn) in [(3, 5), (3, 6), (6, 8), (3, 7)]:
    a3v, _ = poly('kA', pp); cv, _ = poly('kB', nn)
    K3 = sp.expand(sh(a3v, -3) * cv + sh(a3v, -2) * sh(cv, 1) + sh(a3v, -1) * sh(cv, 2))
    la3 = sp.Poly(a3v, E).nth(pp); lc = sp.Poly(cv, E).nth(nn)
    istrue(sp.Poly(K3, E).degree() == pp + nn and sp.expand(sp.Poly(K3, E).all_coeffs()[0] - 3 * la3 * lc) == 0,
           f"k=3 block: deg a3={pp}, deg c={nn} -> K3 deg={pp+nn}, lead = 3*lc(a3)*lc(c) != 0")
# k=2 block filler term  -b2^[-2] a_-2 - b2^[-1] a_-2^[1] : leading coeff = -2*lc(b2)*lc(a_-2) != 0
for (qq, mm) in [(2, 4), (2, 5), (4, 6)]:
    b2v, _ = poly('kB2', qq); av, _ = poly('kAm2', mm)
    K2 = sp.expand(-sh(b2v, -2) * av - sh(b2v, -1) * sh(av, 1))
    lb = sp.Poly(b2v, E).nth(qq); la = sp.Poly(av, E).nth(mm)
    istrue(sp.Poly(K2, E).degree() == qq + mm and sp.expand(sp.Poly(K2, E).all_coeffs()[0] + 2 * lb * la) == 0,
           f"k=2 filler: deg b2={qq}, deg a_-2={mm} -> deg={qq+mm}, lead = -2*lc(b2)*lc(a_-2) != 0")
print("   => the top Q_0-coefficients triangularly annihilate the filler top coefficients.")


# =====================================================================
# 4. MAIN THEOREM: the moment carries no unit (AP class, r SYMBOLIC, d=1).
# =====================================================================
print("\n--- 4. MAIN: moment carries no unit -- AP class {r,r+2,r+4}, r SYMBOLIC, d=1 ---")
r = sp.symbols('r')


SEQ_PIVOTS = []


def seq_solve(eqs, prefer, allowed):
    subs = {}; remaining = list(eqs); changed = True
    SEQ_PIVOTS.clear()
    while changed:
        changed = False
        for eq in list(remaining):
            e = sp.expand(sp.sympify(eq).subs(subs))
            if e == 0:
                remaining.remove(eq); changed = True; break
            for v in prefer:
                if v in subs:
                    continue
                cc = e.coeff(v, 1); c0 = e.coeff(v, 0)
                if cc != 0 and (cc.free_symbols <= allowed) and sp.expand(e - (cc * v + c0)) == 0:
                    subs[v] = sp.together(-c0 / cc); SEQ_PIVOTS.append(cc)
                    remaining.remove(eq); changed = True; break
            if changed:
                break
    return subs, remaining


def slope_residuals(a3, b2, d, rsym=None):
    """Return the residual conditions of {positive cascade} ∪ {higher Q_0 coeffs = 0}
    ∪ {const coeff of Q_0 = c}, after eliminating the free NEGATIVE data (the fillers
    a_-2 and mu3*a_-3, relaxed to independent params) and solving the cascade.
    If the only residual carrying c forces c = 0, the moment carries no unit."""
    A, B, pos, allvars, cs = build_cascade(a3, b2, d, rsym=rsym)
    Q0 = sp.expand(Qm(A, B, 0))
    p = list(sp.symbols(f'p0:{d+1}'))
    subm = {}
    for j in range(d + 1):
        subm[cs['cam3'][j] * cs['mu3']] = p[j]; subm[cs['mu3'] * cs['cam3'][j]] = p[j]
    Q0 = Q0.subs(subm)
    pc = sp.Poly(Q0, E); deg = pc.degree()
    c = sp.symbols('c')
    const = sp.expand(pc.nth(0) - c)
    higher = [sp.expand(pc.nth(k)) for k in range(1, deg + 1)]
    higher = [h for h in higher if h != 0]
    allowed = set() if rsym is None else {rsym}
    prefer = list(cs['ca2']) + list(cs['ca1']) + list(cs['ca0']) + list(cs['cam1'])
    subs, rem = seq_solve(pos, prefer, allowed)
    const = sp.together(const.subs(subs))
    higher = [sp.together(h.subs(subs)) for h in higher] + [sp.together(x) for x in rem]
    higher = [h for h in higher if sp.expand(h) != 0]
    fillers = list(p) + list(cs['cam2'])
    rows = [sp.numer(sp.together(h)) for h in higher] + [sp.numer(const)]
    Mc, rhsc = sp.linear_eq_to_matrix(rows, fillers)
    resid = [sp.expand(n.dot(rhsc)) for n in Mc.T.nullspace()]
    return [sp.factor(x) for x in resid if x != 0], c


a3r = sp.expand((E - r) * (E - r - 2) * (E - r - 4))
b2r = sp.expand((E - r - 1) * (E - r - 4))
az(sh(b2r, 3) * a3r - sh(a3r, 2) * b2r, "AP top {r,r+2,r+4}: b2 solves the wall for ALL r (symbolic)")
resid_r, c = slope_residuals(a3r, b2r, 1, rsym=r)
istrue(all(r not in sp.sympify(cc).free_symbols for cc in SEQ_PIVOTS) and len(SEQ_PIVOTS) >= 6,
       "d=1, r SYMBOLIC: the positive cascade solves through r-INDEPENDENT pivots (uniform in r)")
istrue(resid_r == [c], "d=1, r SYMBOLIC: sole residual = c  =>  moment slope FORCED to 0 for generic r")
print("   => Q_0=1 (slope=1) impossible for the whole degree-3 exotic (=AP) branch at d=1.  [PROVED]")

# Explicit W1 triangular annihilation walk-through (the mechanism, spelled out):
print("\n--- 4b. explicit W1 (r=0) triangular filler-annihilation at d=1 ---")
W1 = (sp.expand(E * (E - 2) * (E - 4)), sp.expand((E - 1) * (E - 4)))
A, B, pos, av, cs = build_cascade(*W1, 1)
# after pos, Q_0 depends only on the fillers am2_0,am2_1, p0=mu3 am3_0, p1=mu3 am3_1:
possol = {cs['ca2'][1]: 0, cs['ca0'][1]: 0, cs['cam1'][1]: 0, cs['cam1'][0]: 0,
          cs['ca1'][0]: cs['ca2'][0]**2, cs['ca1'][1]: -cs['ca2'][0]**2 / 3}
istrue(all(sp.expand(p.subs(possol)) == 0 for p in pos),
       "W1 d=1: the positive cascade solves as a2_1=0,a0_1=0,a_-1=0,a1_0=a2_0^2,a1_1=-a2_0^2/3")
Q0w = sp.expand(Qm(A, B, 0).subs(possol))
p0, p1 = sp.symbols('p0 p1')
Q0w = Q0w.subs({cs['cam3'][0] * cs['mu3']: p0, cs['mu3'] * cs['cam3'][0]: p0,
                cs['cam3'][1] * cs['mu3']: p1, cs['mu3'] * cs['cam3'][1]: p1})
pw = sp.Poly(Q0w, E)
u0, u1 = cs['cam2'][0], cs['cam2'][1]
# top coeffs cascade p1 -> p0 -> u1 -> u0 -> 0, all forced to vanish:
az(pw.nth(6) - 21 * p1, "W1 d=1: [E^6] = 21*p1        => p1 = 0")
az((pw.nth(5) - 9 * (2 * p0 - 15 * p1)), "W1 d=1: [E^5] = 9*(2 p0 - 15 p1)  => p0 = 0")
az((pw.nth(4) + 5 * (2 * u1 + 27 * p0 - 57 * p1)), "W1 d=1: [E^4] = -5*(2 u1 + 27 p0 - 57 p1) => u1 = 0")
az((pw.nth(3) + 8 * u0 - 40 * u1 - 336 * p0 + 441 * p1), "W1 d=1: [E^3] = -8 u0 +40 u1 +336 p0 -441 p1 => u0 = 0")
az(pw.nth(0).subs({u0: 0, u1: 0, p0: 0, p1: 0}),
   "W1 d=1: with fillers = 0 the SLOPE [E^0] = 0 -- but the unit needs 1.  KILLED.")


# =====================================================================
# 5. VERIFICATION across the exotic (AP) class and top-degrees.
# =====================================================================
print("\n--- 5. verification: Q_0=1 infeasible (=[1]), Q_0=0 feasible; exact d=2 certificate ---")
W2 = (sp.expand(E * (E + 2) * (E + 4)), sp.expand(E * (E + 3)))


def ap(rr):
    return (sp.expand((E - rr) * (E - rr - 2) * (E - rr - 4)), sp.expand((E - rr - 1) * (E - rr - 4)))


for (tag, top) in [('W1 r=0', W1), ('W2 r=-4', W2), ('r=1', ap(1)), ('r=-1', ap(-1)),
                   ('r=2', ap(2)), ('r=3', ap(3))]:
    a3, b2 = top
    az(sh(b2, 3) * a3 - sh(a3, 2) * b2, f"{tag}: b2 solves the wall")
    A, B, pos, allvars, cs = build_cascade(a3, b2, 1)
    Gu = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(1)), *allvars, order='grevlex'))
    Gh = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(0)), *allvars, order='grevlex'))
    istrue(Gu == [sp.Integer(1)] and Gh != [sp.Integer(1)],
           f"{tag} d=1: {{cascade}}∪{{Q_0=1}} INFEASIBLE, {{Q_0=0}} feasible  => unit is the killer")

# W1, W2 at d=2 (Q_0=1 infeasible) -- the free-degree beyond the d=1 proof:
for (tag, top) in [('W1', W1), ('W2', W2)]:
    A, B, pos, allvars, cs = build_cascade(*top, 2)
    Gu = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(1)), *allvars, order='grevlex'))
    istrue(Gu == [sp.Integer(1)], f"{tag} d=2: {{cascade}}∪{{Q_0=1}} INFEASIBLE (verified)")

# W1 at d=3, d=4: the slope residual is a NONZERO multiple of c alone -- slope forced
# to 0 directly (no controlling combination w survives; d=2 is the lone exception).
for dd in [3, 4]:
    resd, cc = slope_residuals(*W1, dd)
    istrue(len(resd) == 1 and cc in resd[0].free_symbols
           and sp.simplify(sp.cancel(resd[0] / cc)).free_symbols == set() and resd[0] != 0,
           f"W1 d={dd}: sole slope residual is a nonzero multiple of c  =>  slope forced to 0 (Q_0=1 impossible)")

# EXACT d=2 W1 certificate: the moment slope is forced to 0 via a controlling
# combination w = a1_2^2 (a2_0 - 4 a2_2); the unit demands 7 w = 9 -- contradictory.
resid2, c = slope_residuals(*W1, 2)
# w = a1_2^2 (a2_0 - 4 a2_2), built from build_cascade's coefficient symbols:
a1_2, a2_0, a2_2 = sp.symbols('a1_2 a2_0 a2_2')
wval = a1_2**2 * (a2_0 - 4 * a2_2)
homog = [x for x in resid2 if c not in x.free_symbols]
withc = [x for x in resid2 if c in x.free_symbols]
istrue(len(homog) >= 1 and any(sp.simplify(sp.cancel(h / wval)).free_symbols == set() for h in homog),
       "W1 d=2: a homogeneous residual is a nonzero multiple of w = a1_2^2 (a2_0 - 4 a2_2)  => w = 0")
istrue(len(withc) == 1,
       "W1 d=2: exactly one residual carries the slope c (the unit) -- couples c to w")
# on {w=0} the c-residual forces c = 0 (moment carries no unit):
cres = withc[0]
csol = sp.solve(cres, c)[0]
istrue(sp.simplify(csol.subs({a2_0: 4 * a2_2})) == 0 or sp.simplify(sp.cancel(csol / wval)).free_symbols == set(),
       "W1 d=2: with w=0 the slope c is forced to 0  =>  Q_0=1 impossible (exact certificate)")

# Arbitrary TOP-degree AP: deg a_3 = 3, 6, 9 (g = 1,2,3), d=1:
print("\n--- 5b. arbitrary top-degree AP: deg a_3 = 3,6,9 ---")
for g in [1, 2, 3]:
    roots = list(range(0, 2 * (3 * g), 2))          # {0,2,...,2(3g-1)}, 3g roots, step 2
    a3, b2, br = ap_top(roots)
    az(sh(b2, 3) * a3 - sh(a3, 2) * b2, f"g={g}: step-2 AP top deg {3*g} (roots {roots}) solves the wall")
    A, B, pos, allvars, cs = build_cascade(a3, b2, 1)
    Gu = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(1)), *allvars, order='grevlex'))
    istrue(Gu == [sp.Integer(1)], f"g={g}: {{cascade}}∪{{Q_0=1}} INFEASIBLE  => AP killed at deg a_3={3*g}")


# =====================================================================
# 6. Positive control: the pipeline reproduces a genuine b_2=0 pair.
# =====================================================================
print("\n--- 6. positive control: genuine b_2=0 pair reproduced, no spurious conditions ---")
def mul(P, Q):
    R = {}
    for k1, p1 in P.items():
        for k2, p2 in Q.items():
            R[k1 + k2] = R.get(k1 + k2, 0) + sh(p1, k2) * p2
    return {k: sp.expand(v) for k, v in R.items() if sp.expand(v) != 0}


U = {1: sp.Integer(1), -1: sp.expand(E)}
U3 = mul(mul(U, U), U)
Xc = gauged()
for k, v in U3.items():
    Xc[k] = Xc.get(k, 0) + v
Xc[-1] = sp.expand(Xc[-1] - E)
Dc = gauged(); Dc[1] = sp.Integer(1); Dc[-1] = sp.expand(E)
istrue(all(sp.expand(Qm(Xc, Dc, m) - (1 if m == 0 else 0)) == 0 for m in range(-6, 7)),
       "positive control U=x+d, X=U^3-d, D=U : [D,X]=1 (all Q_m = delta_{m0})")
Av = dict(Xc); Bv = gauged(); Bv[-3] = Dc[-3]
rb1, kb1, cv1 = clean_solve(Av, Bv, 4, 1, 'vb1', 0, 0, 3); Bv[1] = rb1
rbm1, _, cv2 = clean_solve(Av, Bv, 2, -1, 'vbm1', 1, 0, 4); Bv[-1] = rbm1
istrue(len(cv1) == 0 and len(cv2) == 0, "control: clean_solve emits NO spurious conditions")
istrue(sp.expand(rb1.subs({kb1[0]: 1}) - Dc[1]) == 0 and sp.expand(rbm1.subs({kb1[0]: 1}) - Dc[-1]) == 0,
       "control: clean_solve reconstructs D (b1=1, b-1=E) exactly => real feasibility detected")


print("\nALL QUANTUM EXOTIC CLOSURE CHECKS PASSED")
