#!/usr/bin/env python3
"""
verify_moment_unit_k.py
=======================
Exact SymPy checks backing `moment-unit-general-k.md`.  The script verifies
finite symbolic identities and selected fixed-top, fixed-degree polynomial
systems in the general-band-k moment-unit program.  It does not prove uniform
band-k rigidity or an unrestricted DC1 statement.

Conventions (frozen, identical to every sibling quantum memo):
    A_1[x^{-1}] = (+)_k x^k C[E],  (x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E),
    f^[r](E) = f(E+r),  E = x d,
    Q_m = sum_{i+l=m} ( b_l^[i] a_i - a_i^[l] b_l ),  [D,X]=1 <=> Q_m = delta_{m0}.
    Genuine A_1 membership: E(E-1)...(E-r+1) | a_{-r}, b_{-r}.
    Gauge b_k = 0 (spent on the top, Q_{2k}).  Bottom prop. b_{-k}=mu_k a_{-k}.

WHAT THIS SCRIPT CHECKS (see the memo for the written arguments and scope):
  0. ENGINE: Q_0 = (T-1)G for generic coefficients at k=2,3,4,5.
  1. LEMMA P: finite symbolic corroboration at k=2,3,4,5 of the arbitrary-k
     written moment-slope theorem under membership, orientation, gauge, nonzero
     extreme, and normalization hypotheses.
  2. WALL NECKLACE: exact factorization and universal-cofactor examples, plus
     bounded normalized single-coset integer-root scans for k<=5.
  3. LEMMA-R BLOCKS: each of two filler blocks has a nonzero leading coefficient
     in isolation.  The checks do not rule out cross-cancellation between fillers
     or collision with solved blocks, so they do not prove uniform triangularity.
  4. FINITE SYSTEMS: selected fixed-top positive-cascade-plus-Q_0 systems at the
     displayed free degrees, with solved-b raw degree cap 2d+deg(a_k)+2.
  5. r-SYMBOLIC SLICE: one k=4 translation family at d=1 on the forward-solvable
     pivots encoded here.
  6. POSITIVE CONTROL: a genuine band-4 tame pair is admitted (no false kill).

Run:  uv run --with sympy python research/dc1-program/verify_moment_unit_k.py
Ends: ALL MOMENT UNIT K CHECKS PASSED
"""
import sympy as sp
from itertools import combinations

E, S = sp.symbols('E S')


# ---------------------------------------------------------------- primitives
def sh(f, r):
    return sp.expand(sp.sympify(f).subs(E, E + r))


def Sr(r):
    return sp.expand(sum(S**j for j in range(r)))


def poly(name, deg):
    if deg < 0:
        return sp.Integer(0), []
    cs = list(sp.symbols(f'{name}_0:{deg+1}'))
    return sp.expand(sum(cs[i] * E**i for i in range(deg + 1))), cs


def falling(r):
    return sp.prod([E - i for i in range(r)]) if r > 0 else sp.Integer(1)


def Qm(a, b, m, K):
    return sp.expand(sum(sh(b[l], i) * a[i] - sh(a[i], l) * b[l]
                         for i in range(-K, K + 1) for l in range(-K, K + 1) if i + l == m))


def potential_G(a, b, K):
    # W4q band-agnostic potential:
    #   G = sum_{i=1}^K sum_{r=1}^i ( a_i^[-r] b_{-i}^[i-r] - a_{-i}^[i-r] b_i^[-r] )
    return sp.expand(sum(
        sh(a[i], -r) * sh(b[-i], i - r) - sh(a[-i], i - r) * sh(b[i], -r)
        for i in range(1, K + 1) for r in range(1, i + 1)))


def gauged(K):
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


# ---------------------------------------------------------- wall / necklace
def wall_realizable(roots, k):
    """roots: integers, single mod-Z coset.  Wall (root form):
    S_k(S) delta(u) = S S_{k-1}(S) delta(a_k).  Returns (S_k-divisible, u-effective,
    non-cube-cofactor)."""
    A = sp.expand(sum(S**r for r in roots))
    q, rem = sp.div(A, Sr(k), S)
    Sk_div = (sp.expand(rem) == 0)
    num = sp.expand(S * Sr(k - 1) * A)
    qb, remb = sp.div(num, Sr(k), S)
    if sp.expand(remb) != 0:
        u_eff = False
    else:
        Bp = sp.Poly(qb, S)
        u_eff = all(Bp.nth(i) >= 0 for i in range(Bp.degree() + 1))
    noncube = False
    if Sk_div:
        cof = sp.Poly(q, S)
        noncube = any(cof.nth(i) < 0 for i in range(cof.degree() + 1))
    return Sk_div, u_eff, noncube


def top_and_wall(roots, k):
    """Build (a_k, u=b_{k-1}) in E from an integer root list of a_k."""
    a_k = sp.expand(sp.prod([E - r for r in roots]))
    A = sp.expand(sum(S**r for r in roots))
    qb, remb = sp.div(sp.expand(S * Sr(k - 1) * A), Sr(k), S)
    assert sp.expand(remb) == 0, "top is not S_k-admissible"
    Bp = sp.Poly(qb, S)
    uroots = [i for i in range(Bp.degree() + 1) for _ in range(int(Bp.nth(i)))]
    u = sp.expand(sp.prod([E - i for i in uroots]))
    return a_k, u


# ---------------------------------------------------------- cascade solver
def clean_solve(A, B, m, lkey, name, memb, bdeg_raw, K, rsym=None):
    """Solve Q_m = 0 for B[lkey] = falling(memb) * generic(deg bdeg_raw).  Returns
    (particular+kernel solution, kernel free-param list, solvability conditions on the
    OTHER data).  `rsym`, if given, is the single allowed symbolic parameter (the
    translation r) in the operator matrix."""
    braw, bc = poly(f'{name}c_', bdeg_raw)
    bnew = sp.expand(falling(memb) * braw)
    Bt = dict(B); Bt[lkey] = bnew
    eq = sp.expand(Qm(A, Bt, m, K))
    coeffs = sp.Poly(eq, E).all_coeffs() if eq != 0 else [sp.Integer(0)]
    M, rhs = sp.linear_eq_to_matrix(coeffs, bc)
    allowed = set() if rsym is None else {rsym}
    if not all((e.free_symbols <= allowed) for e in M):
        raise AssertionError("operator matrix not numeric/r-parametric (bilinearity leaked)")
    conds = [c for c in (sp.expand(n.dot(rhs)) for n in M.T.nullspace()) if c != 0]
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


def build_cascade(a_k, u, K, d, rsym=None):
    """Gauge b_K=0; wall gives b_{K-1}=u (fixed).  Forward-solve Q_{2K-2}..Q_1 for
    b_{K-2}..b_{-(K-1)}.  Free positive middle a_{K-1}..a_0 deg d; negatives deg d with
    membership; b_{-K}=mu*a_{-K}.  Returns (A,B,pos_conditions,all_free_vars,ca,mu)."""
    p = sp.Poly(a_k, E).degree()
    bdeg = 2 * d + p + 2
    A = {K: a_k}; ca = {}
    for i in range(0, K):
        A[i], ca[i] = poly(f'a{i}', d)
    for i in range(1, K + 1):
        At, ca[-i] = poly(f'am{i}', d)
        A[-i] = sp.expand(falling(i) * At)
    mu = sp.symbols('mu')
    B = gauged(K); B[K - 1] = u; B[-K] = sp.expand(mu * A[-K])
    pos = []; ks = []
    for m in range(2 * K - 2, 0, -1):
        l = m - K; mem = max(0, -l)
        bb, kk, cc = clean_solve(A, B, m, l, f'b{l}', mem, bdeg, K, rsym=rsym)
        B[l] = bb; pos += cc; ks += kk
    allvars = []
    for i in list(range(0, K)) + list(range(-1, -K - 1, -1)):
        allvars += list(ca[i])
    allvars += [mu] + ks
    return A, B, pos, allvars, ca, mu


def q0_conditions(A, B, target, K):
    ex = sp.expand(Qm(A, B, 0, K) - target)
    return [sp.expand(c) for c in sp.Poly(ex, E).all_coeffs() if sp.expand(c) != 0] if ex != 0 else []


# =====================================================================
# 0. ENGINE: finite symbolic checks of Q_0 = (T-1)G.
# =====================================================================
print("--- 0. engine: finite checks of Q_0=(T-1)G (band-agnostic W4q potential) ---")
for K in [2, 3, 4, 5]:
    a = {i: poly(f'A{K}_{i}', 2)[0] for i in range(-K, K + 1)}
    b = {i: poly(f'B{K}_{i}', 2)[0] for i in range(-K, K + 1)}
    az(Qm(a, b, 0, K) - (sh(potential_G(a, b, K), 1) - potential_G(a, b, K)),
       f"k={K}: Q_0 = (T-1)G   (finite symbolic corroboration)")


# =====================================================================
# 1. LEMMA P (moment-slope formula) -- the general protection lemma.
# =====================================================================
print("\n--- 1. LEMMA P: G(0)=0 under membership; G(1)=sum_i[a_i(0)b_-i(i)-a_-i(i)b_i(0)] ---")
for K in [2, 3, 4, 5]:
    a = {}; b = {}
    for i in range(0, K + 1):
        a[i] = poly(f'AA{K}_{i}', 2)[0]; b[i] = poly(f'BB{K}_{i}', 2)[0]
    for i in range(1, K + 1):
        a[-i] = sp.expand(falling(i) * poly(f'Am{K}_{i}', 2)[0])
        b[-i] = sp.expand(falling(i) * poly(f'Bm{K}_{i}', 2)[0])
    G = potential_G(a, b, K)
    az(G.subs(E, 0), f"k={K}: G(0) = 0 identically under membership  =>  Q_0=1 <=> G=E")
    az(sp.Poly(sp.expand(sh(G, 1) - G), E).nth(0) - G.subs(E, 1),
       f"k={K}: slope := const coeff of Q_0 = G(1)")
    LP = sp.expand(sum(a[i].subs(E, 0) * b[-i].subs(E, i) - a[-i].subs(E, i) * b[i].subs(E, 0)
                       for i in range(1, K + 1)))
    az(G.subs(E, 1) - LP,
       f"k={K}: LEMMA P: G(1) = sum_i [ a_i(0) b_-i(i) - a_-i(i) b_i(0) ]")
# gauged form: b_K=0, b_{-K}=mu a_{-K}; level-K term = mu a_K(0) a_{-K}(K)
for K in [2, 3, 4, 5]:
    a = {}; b = {}; mu = sp.symbols('muP')
    for i in range(0, K + 1):
        a[i] = poly(f'GA{K}_{i}', 2)[0]; b[i] = poly(f'GB{K}_{i}', 2)[0]
    b[K] = sp.Integer(0)
    for i in range(1, K + 1):
        a[-i] = sp.expand(falling(i) * poly(f'GAm{K}_{i}', 2)[0])
        b[-i] = sp.expand(falling(i) * poly(f'GBm{K}_{i}', 2)[0])
    b[-K] = sp.expand(mu * a[-K])
    pred = sp.expand(sum(a[i].subs(E, 0) * b[-i].subs(E, i) - a[-i].subs(E, i) * b[i].subs(E, 0)
                         for i in range(1, K)) + mu * a[K].subs(E, 0) * a[-K].subs(E, K))
    az(potential_G(a, b, K).subs(E, 1) - pred,
       f"k={K}: gauged Lemma P: level-K slope term = mu_k a_k(0) a_-k(k) (membership-protected)")


# =====================================================================
# 2. WALL NECKLACE: S_k structure, universal cofactor, minimal exotics.
# =====================================================================
print("\n--- 2. wall necklace: S_k factorization, universal cofactor, minimal exotics ---")
# prime k => S_k = cyclotomic Phi_k; composite k => product (the new slack structure)
for k, expect_irred in [(2, True), (3, True), (4, False), (5, True), (6, False), (7, True)]:
    fac = sp.factor(Sr(k))
    is_irred = not isinstance(fac, sp.Mul)
    istrue(is_irred == expect_irred,
           f"S_{k} = {fac}   ({'prime-cyclotomic' if expect_irred else 'COMPOSITE cyclotomic'})")
# universal cofactor g = 1 - S + S^2 solves the wall for every k>=3
g = 1 - S + S**2
for k in range(3, 8):
    a_k, u = top_and_wall([i for i in range(sp.Poly(sp.expand(Sr(k) * g), S).degree() + 1)
                           for _ in range(int(sp.Poly(sp.expand(Sr(k) * g), S).nth(i)))], k)
    az(sh(u, k) * a_k - sh(a_k, k - 1) * u,
       f"k={k}: universal cofactor g=1-S+S^2 -> a_k solves the wall u^[k]a_k = a_k^[k-1]u")
# Bounded normalized distinct single-coset integer-root scans at degree k.
# The window max(root)<=5k is a finite experimental choice, not a span bound or
# an unrestricted classification; agreement with larger off-verifier windows is
# only a robustness check.
counts = {}
for k in range(2, 6):
    exo = []
    for extra in combinations(range(1, 5 * k + 1), k - 1):
        roots = [0] + list(extra)
        sk, eff, nc = wall_realizable(roots, k)
        if sk and eff and nc:
            exo.append(roots)
    counts[k] = exo
istrue(counts[2] == [], "k=2 bounded scan: no exotic minimal top")
istrue(counts[3] == [[0, 2, 4]], "k=3 bounded scan: one exotic top {0,2,4}")
istrue(len(counts[4]) == 4, f"k=4 bounded scan: four exotic tops {counts[4]}")
istrue(len(counts[5]) == 13, "k=5 bounded scan: thirteen exotic tops (including {0,4,8,12,16})")
print("   => these are finite-window observations, not unrestricted top classifications.")


# =====================================================================
# 3. LEMMA-R BLOCKS: isolated leading coefficients (finite k checks).
# =====================================================================
print("\n--- 3. Lemma-R blocks: isolated leading coefficients (Q_0 fillers) ---")
for k in range(2, 7):
    a_k, _ = poly('la', k); c, _ = poly('lc', k + 2)
    Kk = sp.expand(sum(sh(a_k, -r) * sh(c, k - r) for r in range(1, k + 1)))
    la = sp.Poly(a_k, E).nth(k); lcc = sp.Poly(c, E).nth(k + 2)
    istrue(sp.Poly(Kk, E).degree() == 2 * k + 2
           and sp.expand(sp.Poly(Kk, E).all_coeffs()[0] - k * la * lcc) == 0,
           f"k={k}: level-k block K_k[c]=sum a_k^[-r]c^[k-r] -> lead = k*lc(a_k)*lc(c) != 0")
    bK1, _ = poly('lb', k - 1); am, _ = poly('lam', k + 1)
    Fk = sp.expand(-sum(sh(bK1, -r) * sh(am, k - 1 - r) for r in range(1, k)))
    lb = sp.Poly(bK1, E).nth(k - 1); lam = sp.Poly(am, E).nth(k + 1)
    istrue(sp.expand(sp.Poly(Fk, E).all_coeffs()[0] + (k - 1) * lb * lam) == 0,
           f"k={k}: level-(k-1) filler block -> lead = -(k-1)*lc(b_{{k-1}})*lc(a_-(k-1)) != 0")
print("   => each filler does not self-cancel in isolation.")
print("      Cross-cancellation between fillers and collision with solved blocks remain unchecked;")
print("      these computations do not prove uniform L4 triangular annihilation.")


# =====================================================================
# 4. SELECTED FINITE SYSTEMS: positive cascade plus Q_0.
# =====================================================================
print("\n--- 4. selected finite positive-cascade-plus-Q_0 systems ---")
print("    solved-b raw degree cap: 2d + deg(a_k) + 2")

def kill_check(roots, k, d, tag):
    a_k, u = top_and_wall(roots, k)
    az(sh(u, k) * a_k - sh(a_k, k - 1) * u, f"{tag}: b_{{k-1}} solves the wall")
    A, B, pos, allvars, ca, mu = build_cascade(a_k, u, k, d)
    Gu = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(1), k), *allvars, order='grevlex'))
    Gh = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(0), k), *allvars, order='grevlex'))
    istrue(Gu == [sp.Integer(1)] and Gh != [sp.Integer(1)],
           f"{tag}: encoded positive cascade + Q_0=1 has unit ideal; Q_0=0 ideal is proper")

# k=3: one selected exact AP slice at d=1
kill_check([0, 2, 4], 3, 1, "k=3 AP {0,2,4} d=1 selected slice")
# k=4: the four tops returned by the bounded scan, at d=1
for roots in [[0, 1, 3, 6], [0, 2, 3, 5], [0, 3, 5, 6], [0, 3, 6, 9]]:
    kill_check(roots, 4, 1, f"k=4 exotic {roots} d=1")
# k=4: integer translates of each family (r-independence corroboration)
print("   (integer translates t in {-1,1,2} of each k=4 family:)")
for roots in [[0, 1, 3, 6], [0, 2, 3, 5], [0, 3, 5, 6], [0, 3, 6, 9]]:
    for t in [-1, 1, 2]:
        rr = [x + t for x in roots]
        a_k, u = top_and_wall(rr, 4)
        A, B, pos, allvars, ca, mu = build_cascade(a_k, u, 4, 1)
        Gu = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(1), 4), *allvars, order='grevlex'))
        istrue(Gu == [sp.Integer(1)], f"k=4 exotic {rr} (translate) d=1: Q_0=1 INFEASIBLE")
# k=4 universal family at d=2 (free-degree beyond d=1)
kill_check([0, 2, 3, 5], 4, 2, "k=4 universal {0,2,3,5} d=2")
# k=5 universal family (PRIME cyclotomic Phi_5) at d=1
kill_check([0, 2, 3, 4, 6], 5, 1, "k=5 universal {0,2,3,4,6} d=1 (prime cyclotomic)")


# =====================================================================
# 5. r-SYMBOLIC certificate for the k=4 universal translation family.
# =====================================================================
print("\n--- 5. r-SYMBOLIC slope certificate: k=4 universal family {r,r+2,r+3,r+5}, d=1 ---")
r, cc = sp.symbols('r c')
SEQ_PIVOTS = []


def seq_solve(eqs, prefer, allowed):
    subs = {}; remaining = list(eqs); changed = True; SEQ_PIVOTS.clear()
    while changed:
        changed = False
        for eq in list(remaining):
            e = sp.expand(sp.sympify(eq).subs(subs))
            if e == 0:
                remaining.remove(eq); changed = True; break
            for v in prefer:
                if v in subs:
                    continue
                q1 = e.coeff(v, 1); q0 = e.coeff(v, 0)
                if q1 != 0 and (q1.free_symbols <= allowed) and sp.expand(e - (q1 * v + q0)) == 0:
                    subs[v] = sp.together(-q0 / q1); SEQ_PIVOTS.append(q1)
                    remaining.remove(eq); changed = True; break
            if changed:
                break
    return subs, remaining


def slope_residuals_r(a_k, u, K, d, rsym):
    A, B, pos, allvars, ca, mu = build_cascade(a_k, u, K, d, rsym=rsym)
    Q0 = sp.expand(Qm(A, B, 0, K))
    pj = list(sp.symbols(f'pp0:{d+1}')); subm = {}
    for j in range(d + 1):
        subm[ca[-K][j] * mu] = pj[j]; subm[mu * ca[-K][j]] = pj[j]
    Q0 = Q0.subs(subm)
    pc = sp.Poly(Q0, E)
    const = sp.expand(pc.nth(0) - cc)
    higher = [h for h in (sp.expand(pc.nth(t)) for t in range(1, pc.degree() + 1)) if h != 0]
    allowed = {rsym}
    prefer = sum([list(ca[i]) for i in range(K - 1, -1, -1)], []) \
        + sum([list(ca[-i]) for i in range(1, K - 1)], [])
    subs, rem = seq_solve(pos, prefer, allowed)
    const = sp.together(const.subs(subs))
    higher = [sp.together(h.subs(subs)) for h in higher] + [sp.together(x) for x in rem]
    higher = [h for h in higher if sp.expand(h) != 0]
    fillers = list(pj) + list(ca[-(K - 1)])
    rows = [sp.numer(sp.together(h)) for h in higher] + [sp.numer(const)]
    Mc, rhsc = sp.linear_eq_to_matrix(rows, fillers)
    resid = [sp.expand(n.dot(rhsc)) for n in Mc.T.nullspace()]
    return [sp.factor(x) for x in resid if x != 0]


a4r = sp.expand((E - r) * (E - r - 2) * (E - r - 3) * (E - r - 5))
u4r = sp.expand((E - r - 1) * (E - r - 3) * (E - r - 5))
az(sh(u4r, 4) * a4r - sh(a4r, 3) * u4r, "k=4 universal family: wall holds for ALL r (symbolic)")
resid = slope_residuals_r(a4r, u4r, 4, 1, r)
istrue(all(r not in sp.sympify(pv).free_symbols for pv in SEQ_PIVOTS) and len(SEQ_PIVOTS) >= 8,
       "k=4 d=1, r SYMBOLIC: positive cascade solves through r-INDEPENDENT pivots (uniform in r)")
istrue(resid == [cc],
       "k=4 d=1, r SYMBOLIC: sole residual = c  =>  slope FORCED to 0 for the whole family")
print("   => Q_0=1 (slope=1) impossible for the entire k=4 universal exotic family at once.")


# =====================================================================
# 6. POSITIVE CONTROL: a genuine band-4 tame pair is admitted (no false kill).
# =====================================================================
print("\n--- 6. positive control: genuine band-4 tame pair admitted (no false kill) ---")


def mul(P, Q):
    R = {}
    for k1, p1 in P.items():
        for k2, p2 in Q.items():
            R[k1 + k2] = R.get(k1 + k2, 0) + sh(p1, k2) * p2
    return {k: sp.expand(v) for k, v in R.items() if sp.expand(v) != 0}


K = 4
U = {1: sp.Integer(1), -1: sp.expand(E)}
Uk = U
for _ in range(K - 1):
    Uk = mul(Uk, U)
Xc = gauged(K)
for k, v in Uk.items():
    Xc[k] = Xc.get(k, 0) + v
Xc[-1] = sp.expand(Xc[-1] - E)                       # X = U^4 - d
Dc = gauged(K); Dc[1] = sp.Integer(1); Dc[-1] = sp.expand(E)   # D = U
istrue(all(sp.expand(Qm(Xc, Dc, m, K) - (1 if m == 0 else 0)) == 0 for m in range(-2 * K, 2 * K + 1)),
       "control U=x+d, X=U^4-d, D=U : [D,X]=1 (all Q_m = delta_{m0}), a_4=1 trivial cube")
# forward solver: the only condition forces a periodicity-kernel param to 0 (b_2=0 tame),
# and setting the b_1-kernel param to 1 recovers D exactly -- so no spurious infeasibility.
Av = dict(Xc); Bv = gauged(K); Bv[-K] = Dc[-K]; Bv[K - 1] = sp.Integer(0)
conds = []; b1ker = None
for m in range(2 * K - 2, 0, -1):
    l = m - K; mem = max(0, -l)
    bb, kk, ccnd = clean_solve(Av, Bv, m, l, f'vb{l}', mem, 6, K)
    Bv[l] = bb; conds += ccnd
    if l == 1:
        b1ker = kk
istrue(all(sp.expand(cnd.subs({s: 0 for s in cnd.free_symbols})) == 0 for cnd in conds),
       "control: every forward-solver condition forces a kernel param to 0 (consistent with b_2=0)")
istrue(b1ker and sp.expand(Bv[1].subs({b1ker[0]: 1, **{s: 0 for s in Bv[1].free_symbols
                                                        if s != b1ker[0]}}) - Dc[1]) == 0,
       "control: setting the b_1 periodicity-kernel param to 1 recovers D (b_1=1) => genuine pair admitted")


print("\nALL MOMENT UNIT K CHECKS PASSED")
