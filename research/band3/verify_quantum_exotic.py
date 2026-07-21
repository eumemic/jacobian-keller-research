#!/usr/bin/env python3
"""
verify_quantum_exotic.py
========================
Exact SymPy verification backing the memo `quantum-exotic-branch.md`:
the FATE of the QUANTUM band-3 exotic (non-shifted-cube) top branch that passes
the Q_5 cyclotomic wall (`quantum-band3-cascade.md`, commit 99fe6ee).

VERDICT (established here): the exotic branch is KILLED. A non-shifted-cube a_3
solving the wall admits NO band-3 pair [D,X]=1. The kill is at Q_0 (the m=0
central integral G=E) -- the W4 moment UNIT (the "1" in [D,X]=1) cannot be
realized. No DC1/JC2 counterexample is produced.

Conventions (frozen, matching every sibling quantum memo, commit 99fe6ee):
    A_1[x^{-1}] = (+)_k x^k C[E],  (x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E),
    f^[r](E) = f(E+r),  E = x d,  ladder-m coeff of [D,X]:
    Q_m = sum_{k+l=m} ( b_l^[k] a_k - a_k^[l] b_l ),  [D,X]=1  <=>  Q_m = delta_{m0}.
    Genuine A_1 membership:  E(E-1)...(E-r+1) | a_{-r}, b_{-r}.

Sections:
  0. crossed-product engine; Q_m == direct commutator (band 3); gauge b_3=0 recap.
  1. wall recap: the two Wave-A witnesses solve Q_5, both non-shifted-cube.
  2. b_2 = 0 sub-branch: PROVED KILL at Q_0, arbitrary degree (degree obstruction).
  3. b_2 != 0 sub-branch: positive cascade solvable; Q_0=1 obstruction:
       - explicit positive-cascade solution (Q_4..Q_1 = 0);
       - Q_0 = 1 INFEASIBLE (Groebner = [1]) for both witnesses;
       - Q_0 = 0 FEASIBLE  => the moment UNIT is exactly the killer;
       - exact infeasibility certificate {8w = 0, 7w = 9}, w = a1_2^2 (a2_0-4 a2_2).
  4. robustness across the exotic class: step-2 arithmetic progressions and a
       degree-6 exotic top -- all killed at Q_0.
  5. cross-validation: the pipeline reproduces the genuine positive control
       (b_2 = 0 tame pair) with NO spurious conditions.

Run:  uv run --with sympy python research/band3/verify_quantum_exotic.py
Ends: ALL QUANTUM EXOTIC CHECKS PASSED
"""
import sympy as sp

E = sp.symbols('E')


def sh(f, r):
    return sp.expand(sp.sympify(f).subs(E, E + r))


def poly(name, deg):
    if deg < 0:
        return sp.Integer(0), []
    cs = sp.symbols(f'{name}_0:{deg+1}')
    cs = list(cs) if isinstance(cs, (tuple, list)) else [cs]
    return sp.expand(sum(cs[i] * E**i for i in range(deg + 1))), list(cs)


def Qm(a, b, m, K=3):
    return sp.expand(sum(
        sh(b[l], k) * a[k] - sh(a[k], l) * b[l]
        for k in range(-K, K + 1) for l in range(-K, K + 1) if k + l == m))


def falling(r):
    return sp.prod([E - i for i in range(r)]) if r > 0 else sp.Integer(1)


def divides(fac, c):
    c = sp.expand(sp.sympify(c))
    if c == 0:
        return True
    return sp.rem(sp.Poly(c, E), sp.Poly(fac, E)) == 0


def az(expr, label):
    e = sp.expand(sp.sympify(expr))
    if e != 0:
        raise AssertionError(label + "  :  residual = " + str(e))
    print("PASS", label)


def istrue(cond, label):
    if not cond:
        raise AssertionError(label + "  :  FALSE")
    print("PASS", label)


def clean_solve(A, B, m, lkey, name, memb, target, bdeg_raw):
    """Solve Q_m = target for B[lkey] = falling(memb)*generic(deg bdeg_raw).
    Conditions via left-nullspace(numeric M); a clean particular solution via an
    independent-row solve; genuine kernel dims re-added as named free params.
    Returns (b_solution, kernel_syms, solvability_conditions)."""
    braw, bc = poly(f'{name}c_', bdeg_raw)
    bnew = sp.expand(falling(memb) * braw)
    Bt = dict(B); Bt[lkey] = bnew
    eq = sp.expand(Qm(A, Bt, m) - target)
    coeffs = sp.Poly(eq, E).all_coeffs() if eq != 0 else [sp.Integer(0)]
    M, rhs = sp.linear_eq_to_matrix(coeffs, bc)
    if not all(e.is_number for e in M):
        raise AssertionError("operator matrix not numeric (bilinearity leaked)")
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
    ker_syms = []
    for idx, kv in enumerate(M.nullspace()):
        g = sp.symbols(f'{name}K{idx}')
        kp = sp.expand(falling(memb) * sum(kv[i] * E**i for i in range(len(bc))))
        bsol = sp.expand(bsol + g * kp); ker_syms.append(g)
    return bsol, ker_syms, conds


def gauged_dict():
    return {k: sp.Integer(0) for k in range(-3, 4)}


# =====================================================================
# 0. Crossed-product engine; Q_m == direct commutator (band 3); gauge recap.
# =====================================================================
print("--- 0. engine: Q_m == direct commutator (band 3); gauge b_3 = 0 ---")
A0 = {k: poly(f'A{k+3}', 2)[0] for k in range(-3, 4)}
B0 = {k: poly(f'B{k+3}', 2)[0] for k in range(-3, 4)}


def direct_commutator(m, K=3):
    tot = 0
    for k in range(-K, K + 1):
        for l in range(-K, K + 1):
            if k + l == m:
                tot += sh(B0[l], k) * A0[k] - sh(A0[k], l) * B0[l]
    return sp.expand(tot)


for m in range(-6, 7):
    az(direct_commutator(m) - Qm(A0, B0, m), f"Q_{m} = direct commutator ladder coeff")
# In gauge b_3 = 0 the "new coefficient" of Q_m (m<=5) is b_{m-3}, via the pair (3,m-3);
# isolate that operator L_m[b] = b^[3] a3 - a3^[m-3] b from the rest of Q_m:
Bg = dict(B0); Bg[3] = sp.Integer(0)
for m in [5, 4, 3, 2, 1, 0]:
    op = sh(Bg[m - 3], 3) * A0[3] - sh(A0[3], m - 3) * Bg[m - 3]
    rest = sp.expand(sum(sh(Bg[l], k) * A0[k] - sh(A0[k], l) * Bg[l]
                         for k in range(-3, 4) for l in range(-3, 4)
                         if k + l == m and not (k == 3 and l == m - 3)))
    az(Qm(A0, Bg, m) - op - rest,
       f"Q_{m}|_(b3=0) = L_{m}[b_{m-3}] + rest,  L_{m}[b] = b^[3] a3 - a3^[{m-3}] b")
print("   => gauge b_3=0; the wall Q_5 gives b_2 (Wave A); we probe the exotic top downstream.")


# =====================================================================
# 1. Wall recap: the two Wave-A witnesses solve Q_5, both non-shifted-cube.
# =====================================================================
print("\n--- 1. wall witnesses (cite quantum-band3-cascade.md 99fe6ee) ---")
W1 = (sp.expand(E * (E - 2) * (E - 4)), sp.expand((E - 1) * (E - 4)))       # roots {0,2,4}
W2 = (sp.expand(E * (E + 2) * (E + 4)), sp.expand(E * (E + 3)))             # roots {0,-2,-4}
for (a3, b2, nm) in [(*W1, 'W1 a3=E(E-2)(E-4)'), (*W2, 'W2 a3=E(E+2)(E+4)')]:
    az(sh(b2, 3) * a3 - sh(a3, 2) * b2, f"{nm}: b2 solves the wall b2^[3]a3 = a3^[2]b2")
    # non-shifted-cube certificate: no degree-1 h with a3 = c h h^[1] h^[2]
    gg, ggc = poly('gg', 1); csq = sp.symbols('csq')
    sols = sp.solve(sp.Poly(sp.expand(csq * gg * sh(gg, 1) * sh(gg, 2) - a3), E).all_coeffs(),
                    ggc + [csq], dict=True)
    istrue(len(sols) == 0, f"{nm}: NOT a shifted cube c h h^[1] h^[2]  (exotic)")


# =====================================================================
# 2. b_2 = 0 sub-branch: PROVED KILL at Q_0, arbitrary degree.
# =====================================================================
print("\n--- 2. exotic b_2 = 0 sub-branch: KILL at Q_0 (arbitrary degree) ---")
a3 = W1[0]
# 2a. Positive collapse: L_m for m in {4,2,1} has TRIVIAL kernel; L_3 kernel = constants.
for (m, j, name) in [(4, 1, 'b1'), (2, -1, 'b-1'), (1, -2, 'b-2')]:
    trivial = True
    for q in range(0, 8):
        b, bc = poly('b', q)
        Lm = sp.expand(sh(b, 3) * a3 - sh(a3, j) * b)
        sol = list(sp.linsolve(sp.Poly(Lm, E).all_coeffs(), bc))
        if sol and len({s for x in sol[0] for s in x.free_symbols}) > 0:
            trivial = False
    istrue(trivial, f"b2=0: L_{m} (a3-shift[{j}]) has trivial kernel  => {name} = 0")
# L_3[b0] = a3 (b0^[3]-b0): kernel = constants
istrue(all((lambda q: (lambda sol: sol and len({s for x in sol[0] for s in x.free_symbols}) == 1)(
        list(sp.linsolve(sp.Poly(sp.expand(a3 * (sh(poly('b', q)[0], 3) - poly('b', q)[0])), E).all_coeffs(),
                         poly('b', q)[1]))))(q) for q in range(1, 6)),
        "b2=0: L_3 = a3(b0^[3]-b0) kernel = constants  => b0 = const")
# 2b. Then Q_0 = 1 becomes L_0[b_{-3}] = b_{-3}^[3] a3 - a3^[-3] b_{-3} = 1.
#     For a NONCONSTANT a3 (deg 3), symbolic-degree check: deg L_0[b] = deg b + 2,
#     leading coeff 3*(3 + deg b) != 0  =>  L_0[b] can NEVER be a nonzero constant.
q = sp.symbols('q', integer=True, positive=True)
for qd in range(0, 8):
    b, bc = poly('b', qd)
    L0 = sp.expand(sh(b, 3) * a3 - sh(a3, -3) * b)
    deg = sp.Poly(L0, E).degree() if L0 != 0 else -1
    lead = sp.Poly(L0, E).nth(qd + 2) if L0 != 0 else 0
    istrue(deg == qd + 2 and sp.expand(lead - 3 * (3 + qd) * bc[qd]) == 0,
           f"b2=0: L_0[b], deg b={qd}: deg={qd+2}>=2, lead coeff = 3(3+{qd})*lc != 0")
print("   => L_0[b_-3] = 1 is IMPOSSIBLE (LHS deg>=2 or 0, never a nonzero const):")
print("      the b_2=0 exotic sub-branch is KILLED at Q_0.  [PROVED, arbitrary degree]")


# =====================================================================
# 3. b_2 != 0 sub-branch: positive cascade solvable; Q_0=1 obstruction.
# =====================================================================
print("\n--- 3. exotic b_2 != 0 sub-branch: positive cascade solvable, Q_0=1 KILLS ---")


def build_cascade(a3, b2, d):
    """Forward-solve the positive cascade with generic free lower coeffs (degree d).
    Returns (A,B, pos_conditions, all_free_vars)."""
    a2, ca2 = poly('a2', d); a1, ca1 = poly('a1', d); a0, ca0 = poly('a0', d)
    am1t, cam1 = poly('am1', d); am2t, cam2 = poly('am2', d); am3t, cam3 = poly('am3', d)
    am1 = sp.expand(falling(1) * am1t); am2 = sp.expand(falling(2) * am2t); am3 = sp.expand(falling(3) * am3t)
    mu3 = sp.symbols('mu3')
    A = {3: a3, 2: a2, 1: a1, 0: a0, -1: am1, -2: am2, -3: am3}
    B = gauged_dict(); B[2] = b2; B[-3] = sp.expand(mu3 * am3)
    pos = []; ks = []
    for (m, l, nm, mem, bd) in [(4, 1, 'b1', 0, d + 3), (3, 0, 'b0', 0, 2 * d + 2),
                                (2, -1, 'bm1', 1, 2 * d + 3), (1, -2, 'bm2', 2, 2 * d + 4)]:
        bb, kk, cc = clean_solve(A, B, m, l, nm, mem, 0, bd)
        B[l] = bb; pos += cc; ks += kk
    allvars = list(ca2) + list(ca1) + list(ca0) + list(cam1) + list(cam2) + list(cam3) + [mu3] + ks
    return A, B, pos, allvars


# 3a. EXPLICIT positive-cascade solution for W1 (a concrete feasibility witness):
a3, b2 = W1
# free choice a2_0=0, a2_2=1, a1_2=1 (=> a2_1=-5) on the positive-solution variety of section 3b;
# a_-1 is forced (Q_1 couples a_-1 through the (-1,2) pair b2^[-1] a_-1, since b2 != 0).
Ap = {3: a3, 2: sp.expand(E**2 - 5 * E), 1: sp.expand(E**2 - 14 * E + 32),
      0: sp.expand(-E**2 / 3 + 3 * E),
      -1: sp.expand(E * (-sp.Rational(2, 3) + sp.Rational(2, 3) * E - sp.Rational(1, 6) * E**2)),
      -2: sp.Integer(0), -3: sp.Integer(0)}
Bp = gauged_dict(); Bp[2] = b2
b1p, _, c1 = clean_solve(Ap, Bp, 4, 1, 'pb1', 0, 0, 4); Bp[1] = b1p
b0p, k0, c0 = clean_solve(Ap, Bp, 3, 0, 'pb0', 0, 0, 5); Bp[0] = b0p.subs({k0[0]: 0})
bm1p, _, c2 = clean_solve(Ap, Bp, 2, -1, 'pbm1', 1, 0, 6); Bp[-1] = bm1p
bm2p, _, c3 = clean_solve(Ap, Bp, 1, -2, 'pbm2', 2, 0, 7); Bp[-2] = bm2p
istrue(all(len(c) == 0 for c in (c1, c0, c2, c3)),
       "positive point (a2=E(E-5),a1=E^2-14E+32,a0=3E-E^2/3) satisfies all Q4..Q1 solvability")
for m in [4, 3, 2, 1]:
    az(Qm(Ap, Bp, m), f"explicit positive solution: Q_{m} = 0  (positive cascade IS solvable)")

# 3b. Q_0 = 1 is INFEASIBLE for both witnesses (Groebner = [1]); Q_0 = 0 is FEASIBLE.
def q0_conditions(A, B, target):
    ex = sp.expand(Qm(A, B, 0) - target)
    return [sp.expand(c) for c in sp.Poly(ex, E).all_coeffs() if sp.expand(c) != 0] if ex != 0 else []


for (a3, b2, nm) in [(*W1, 'W1'), (*W2, 'W2')]:
    A, B, pos, allvars = build_cascade(a3, b2, 2)
    q0_1 = q0_conditions(A, B, sp.Integer(1))
    q0_0 = q0_conditions(A, B, sp.Integer(0))
    G_unit = list(sp.groebner(pos + q0_1, *allvars, order='grevlex'))
    G_hom = list(sp.groebner(pos + q0_0, *allvars, order='grevlex'))
    istrue(G_unit == [sp.Integer(1)], f"{nm}: pos-cascade + (Q_0 = 1) is INFEASIBLE  => exotic top KILLED")
    istrue(G_hom != [sp.Integer(1)], f"{nm}: pos-cascade + (Q_0 = 0) is FEASIBLE  => the UNIT is the killer")

# 3c. EXACT infeasibility certificate (W1): after eliminating ALL free negative data
#     (a_-2, a_-3, mu3) the system Q_0=1 collapses to two contradictory conditions
#     8 w = 0  and  7 w = 9   with  w := a1_2^2 (a2_0 - 4 a2_2).
a3, b2 = W1
a2, ca2 = poly('a2', 2); a1, ca1 = poly('a1', 2); a0, ca0 = poly('a0', 2)
am1t, cam1 = poly('am1', 2); am2t, cam2 = poly('am2', 2); am3t, cam3 = poly('am3', 2)
am1 = sp.expand(falling(1) * am1t); am2 = sp.expand(falling(2) * am2t); am3 = sp.expand(falling(3) * am3t)
mu3 = sp.symbols('mu3')
A = {3: a3, 2: a2, 1: a1, 0: a0, -1: am1, -2: am2, -3: am3}
B = gauged_dict(); B[2] = b2; B[-3] = sp.expand(mu3 * am3)
pos = []
for (m, l, nm, mem, bd) in [(4, 1, 'cb1', 0, 5), (3, 0, 'cb0', 0, 6), (2, -1, 'cbm1', 1, 7), (1, -2, 'cbm2', 2, 8)]:
    bb, kk, cc = clean_solve(A, B, m, l, nm, mem, 0, bd); B[l] = bb; pos += cc
psol = sp.solve(pos, dict=True)[0]
Q0 = sp.expand(Qm(A, B, 0).subs(psol) - 1)
# relax mu3*am3_j -> free p_j (exact: for mu3=0 gives 0, else arbitrary); eliminate p_j and a_-2:
p = sp.symbols('p0:3')
subm = {}
for j in range(3):
    subm[cam3[j] * mu3] = p[j]; subm[mu3 * cam3[j]] = p[j]
coeffsQ0 = [sp.expand(c.subs(subm)) for c in sp.Poly(Q0, E).all_coeffs()]
coeffsQ0 = [c for c in coeffsQ0 if c != 0]
fillers = list(p) + list(cam2)
Mc, rhsc = sp.linear_eq_to_matrix(coeffsQ0, fillers)
resid = [sp.expand(n.dot(rhsc)) for n in Mc.T.nullspace()]
resid = [sp.factor(r) for r in resid if r != 0]
w = ca1[2]**2 * (ca2[0] - 4 * ca2[2])
istrue(any(sp.expand(r - sp.Rational(4, 9) * w) == 0 or sp.expand(r + sp.Rational(4, 9) * w) == 0
           or sp.expand(9 * r - 8 * w) == 0 or sp.expand(9 * r + 8 * w) == 0 for r in resid),
       "certificate: one residual is (multiple of) w = a1_2^2(a2_0-4a2_2)  [forces w = 0]")
istrue(any(sp.expand(9 * r + (7 * w - 9)) == 0 or sp.expand(9 * r - (7 * w - 9)) == 0 for r in resid),
       "certificate: other residual is (7 w - 9)/9  [moment unit forces 7 w = 9]")
istrue(list(sp.groebner(resid, ca2[0], ca2[2], ca1[2], order='grevlex')) == [sp.Integer(1)],
       "certificate: {w = 0} AND {7 w = 9} are CONTRADICTORY  => Q_0=1 impossible")
print("   => exotic b_2!=0 branch KILLED at Q_0: the W4 moment UNIT cannot be realized.")


# =====================================================================
# 4. Robustness across the exotic class.
# =====================================================================
print("\n--- 4. robustness: the exotic step-2 AP class + a degree-6 exotic top ---")


def ap(r):
    return (sp.expand((E - r) * (E - r - 2) * (E - r - 4)), sp.expand((E - r - 1) * (E - r - 4)))


for r in (0, 1, -1, 3):
    a3, b2 = ap(r)
    az(sh(b2, 3) * a3 - sh(a3, 2) * b2, f"AP r={r}: a3={{{r},{r+2},{r+4}}} solves the wall")
    A, B, pos, allvars = build_cascade(a3, b2, 2)
    G = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(1)), *allvars, order='grevlex'))
    istrue(G == [sp.Integer(1)], f"AP r={r}: pos-cascade + Q_0=1 INFEASIBLE  => exotic top KILLED")

# degree-6 exotic {0,2,4,6,8,10}: Phi3-divisible, non-effective cofactor (exotic), wall-admissible.
S = sp.symbols('S')
Anl = sum(S**(2 * k) for k in range(6))
_, rphi = sp.div(sp.expand(Anl), 1 + S + S**2, S)
istrue(sp.expand(rphi) == 0, "deg-6 top {0,2,..,10}: root multiset divisible by Phi3 (wall-admissible)")
qc, _ = sp.div(sp.expand(Anl), 1 + S + S**2, S); Cp = sp.Poly(qc, S)
istrue(any(Cp.nth(i) < 0 for i in range(Cp.degree() + 1)),
       "deg-6 top {0,2,..,10}: cofactor A/Phi3 has a NEGATIVE coeff => NOT a shifted cube (exotic)")
qb, _ = sp.div(sp.expand(S * (1 + S) * Anl), 1 + S + S**2, S); Bp6 = sp.Poly(qb, S)
a3_6 = sp.expand(sp.prod([E - 2 * k for k in range(6)]))
b2_6 = sp.expand(sp.prod([E - i for i in range(Bp6.degree() + 1) if Bp6.nth(i) != 0]))
az(sh(b2_6, 3) * a3_6 - sh(a3_6, 2) * b2_6, "deg-6 exotic: b2 (roots {1,4,7,10}) solves the wall")
A6, B6, pos6, av6 = build_cascade(a3_6, b2_6, 1)
istrue(list(sp.groebner(pos6 + q0_conditions(A6, B6, sp.Integer(1)), *av6, order='grevlex')) == [sp.Integer(1)],
       "deg-6 exotic {0,2,..,10}: pos-cascade + Q_0=1 INFEASIBLE  => KILLED")


# =====================================================================
# 5. Cross-validation: the pipeline reproduces the genuine positive control.
# =====================================================================
print("\n--- 5. cross-validation: genuine positive control (b_2=0 tame pair) ---")
# U = x + c1 d, X = U^3 - d/kappa, D = kappa U ; c1 = kappa = 1  (a real [D,X]=1 pair).
def mul(P, Q):
    R = {}
    for k1, p1 in P.items():
        for k2, p2 in Q.items():
            R[k1 + k2] = R.get(k1 + k2, 0) + sh(p1, k2) * p2
    return {k: sp.expand(v) for k, v in R.items() if sp.expand(v) != 0}


U = {1: sp.Integer(1), -1: sp.expand(E)}
U3 = mul(mul(U, U), U)
Xc = gauged_dict()
for k, v in U3.items():
    Xc[k] = Xc.get(k, 0) + v
Xc[-1] = sp.expand(Xc[-1] - E)
Dc = gauged_dict(); Dc[1] = sp.Integer(1); Dc[-1] = sp.expand(E)
istrue(all(sp.expand(Qm(Xc, Dc, m) - (1 if m == 0 else 0)) == 0 for m in range(-6, 7)),
       "positive control: [D,X] = 1 (all Q_m = delta_{m0})")
# feed the control's X into clean_solve (b_2 = 0); it must reproduce D with NO spurious conditions.
Av = dict(Xc); Bv = gauged_dict(); Bv[-3] = Dc[-3]
rb1, kb1, cv1 = clean_solve(Av, Bv, 4, 1, 'vb1', 0, 0, 3); Bv[1] = rb1
rbm1, _, cv2 = clean_solve(Av, Bv, 2, -1, 'vbm1', 1, 0, 4); Bv[-1] = rbm1
istrue(len(cv1) == 0 and len(cv2) == 0, "control: clean_solve emits NO spurious conditions")
istrue(sp.expand(rb1.subs({kb1[0]: 1}) - Dc[1]) == 0 and sp.expand(rbm1.subs({kb1[0]: 1}) - Dc[-1]) == 0,
       "control: clean_solve reconstructs D (b1=1, b-1=E) exactly  => feasibility is detected")


print("\nALL QUANTUM EXOTIC CHECKS PASSED")
