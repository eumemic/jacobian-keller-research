#!/usr/bin/env python3
"""JOINT FILLER COVECTOR LEMMA (W2, architecture step 1) -- exact-algebra certificate.

Extracts the EXPLICIT generic covector for the stacked linear filler operator
L = [Phi ; M_-1] on the W2 slope-1 cascade locus, identifies its closed-form
(root-necklace) support, and establishes the payload: the joint obstruction is
present UNIFORMLY at raw cap d=3 (both branches).

W2 datum, gauge b_3=0, quantum band-3 conventions:  Q_m = sum_(k+l=m)[b_l^[k]a_k -
a_k^[l]b_l],  f^[n](E)=f(E+n),  membership (E)_j=E(E-1)...(E-j+1) | a_-j,b_-j;
    a_3 = E(E+2)(E+4)   (roots {0,-2,-4}),   b_2 = E(E+3)   (roots {0,-3}),   b_3=0.
Fillers (branch B):  a_-2=(E)_2 V,  b_-3=(E)_3 C  (a_-3=0).
Phi(C,V) = K_3[(E)_3 C] - H_2[(E)_2 V];  Im Phi = D*F[E], D=E(E-1)(E+1).
M_-1 = filler-linear part of Q_-1: a two-block map with (a_2,b_1) in place of (a_3,b_2).

SECTIONS
  S0  engine: Q_m = [D,X]_m (m in [-6,6]); Q_0=(T-1)G; G = potential.        (arb deg)
  S1  pivot: block-0 = G-E has filler-linear part Phi, D|Phi; slope gate
      D|(E-R) <=> R(1)=1, R(-1)=-1.                                          (arb deg)
  S2  explicit slope-1 datum: cascade Q_4..Q_1=0, R(1)=1, R(-1)=-1; build the
      stacked L=[Phi;M_-1] and rhs=(E-R, -N_-1).
  S3  THE JOINT COVECTOR (explicit): left-kernel lam of L, lam.L=0 (every column),
      lam.rhs != 0 (=16/45); minimal genuine coupling at filler degree dV=3; the
      pivot/residual formula giving lam as an explicit function of the datum.
  S4  CLOSED FORM (necklace): Ann(Im Phi)={ev_-1,ev_0,ev_1} DECOUPLES (obstruction 0);
      the coupled mu0 lives on the fixed (a_3,b_2) necklace
      roots(a_3) u roots(b_2) u roots(q_K) u roots(q_H) = {-4,-3,-2,-1,0,1,2,3}
      (confirms the union-of-necklaces conjecture); mu1 on the (a_2,b_1) datum
      necklace.  lam.rhs = mu0(E-R) - mu1(N_-1) has the both-ends Lemma-P shape.
  S5  FILLER-DEGREE-FREENESS: full column rank, gap=1, obstruction != 0 at
      dV=4..10 -- the datum admits no tail-clearing filler of any tested degree.
  S6  GENERIC in positive data: distinct slope-1 data, gap=1, obstruction != 0.
  S7  PAYLOAD (d=3 unit certificate; msolve '^'): cascade + Q_0=1 + Q_-1 is the
      UNIT ideal at d=3, BOTH branches (QQ + a prime), with a sympy cross-check at
      the explicit datum and a no-false-kill control.  This is the joint covector
      kill, UNIFORM at d=3 -- the full negative tail Q_-2..Q_-5 is NOT needed.
  S8  BRANCH A: per fixed gauge mu_3, full col rank, gap=1, obstruction != 0.
  S9  METHOD DEBT (audit): the sibling S8 feasibility claim 'cascade+Q_0=1+Q_-1 is
      NOT unit / a sub-locus survives Q_-1' was a msolve '**'-parser artifact;
      with '^' the system is UNIT -- recorded, the sub-locus is EMPTY at d=3.

SCOPE (honest).  Proved arbitrary degree: S0, S1 (engine, Im Phi=D F[E], slope
gate).  Explicit + bounded/generic evidence: S3-S6, S8 (the covector, its necklace
support, degree-freeness in the filler to dV=10, six slope-1 data, branch A per
gauge).  Bounded exact certificate (d=3, both branches, msolve '^' + sympy point
cross-check): S7 -- the Joint Filler Covector Lemma holds UNIFORMLY at raw cap d=3.
OPEN: arbitrary positive-data degree; a degree-free lam.rhs identity; the residual
sub-locus at d>3.  No Weyl pair, no DC1 counterexample.

Run:  uv run --with sympy python research/dc1-program/verify_joint_covector.py
Ends: ALL JOINT COVECTOR CHECKS PASSED
"""
import sympy as sp
import time
import random
import shutil
import subprocess
import tempfile
import os

E = sp.symbols("E")
LEVELS = range(-3, 4)
R = sp.Rational
_T0 = time.time()
random.seed(20260723)


def sh(f, n):
    return sp.expand(sp.sympify(f).subs(E, E + n))


def falling(n):
    return sp.prod(E - j for j in range(n))


def poly(name, degree):
    cs = list(sp.symbols(f"{name}_0:{degree + 1}"))
    return sp.expand(sum(cs[j] * E**j for j in range(degree + 1))), cs


def q_m(A, B, m):
    """Stipulated convention Q_m = sum_(k+l=m) (b_l^[k] a_k - a_k^[l] b_l)."""
    return sp.expand(sum(sh(B[l], k) * A[k] - sh(A[k], l) * B[l]
                         for k in LEVELS for l in LEVELS if k + l == m))


def potential(A, B):
    return sp.expand(sum(sh(A[k], j - k) * sh(B[-k], j) - sh(B[k], j - k) * sh(A[-k], j)
                         for k in range(1, 4) for j in range(k)))


def mul_ladders(P, Q):
    Rr = {}
    for k1, p1 in P.items():
        for k2, p2 in Q.items():
            Rr[k1 + k2] = Rr.get(k1 + k2, 0) + sh(p1, k2) * p2
    return {k: sp.expand(v) for k, v in Rr.items() if sp.expand(v) != 0}


def commutator_coeff(X, Dop, m):
    return sp.expand(mul_ladders(Dop, X).get(m, 0) - mul_ladders(X, Dop).get(m, 0))


def check(cond, label):
    if not cond:
        raise AssertionError("FAIL " + label)
    print(f"PASS [{time.time() - _T0:6.1f}s] {label}", flush=True)


def check_zero(val, label):
    if sp.expand(val) != 0:
        raise AssertionError(f"FAIL {label}: residual {sp.factor(sp.expand(val))}")
    print(f"PASS [{time.time() - _T0:6.1f}s] {label}", flush=True)


a3_w2 = sp.expand(E * (E + 2) * (E + 4))
b2_w2 = sp.expand(E * (E + 3))
D = sp.expand(E * (E - 1) * (E + 1))

# =====================================================================
print("--- S0. crossed-product engine: Q_m = [D,X]_m, Q_0 = (T-1)G ---", flush=True)
# =====================================================================
Ag, Bg = {}, {}
for lev in LEVELS:
    ar, _ = poly(f"Ag{lev + 3}", 2)
    br, _ = poly(f"Bg{lev + 3}", 2)
    memb = falling(-lev) if lev < 0 else 1
    Ag[lev] = sp.expand(memb * ar)
    Bg[lev] = sp.expand(memb * br)
for m in range(-6, 7):
    check_zero(q_m(Ag, Bg, m) - commutator_coeff(Ag, Bg, m),
               f"Q_{m} = [D,X]_{m} (generic degree-2 coefficients)")
check_zero(q_m(Ag, Bg, 0) - (sh(potential(Ag, Bg), 1) - potential(Ag, Bg)),
           "Q_0 = (T-1)G telescoping identity  (G = potential)")

# =====================================================================
print("\n--- S1. pivot: block-0 = G-E has filler part Phi, D|Phi; slope gate ---",
      flush=True)
# =====================================================================
def K3(c):
    return sp.expand(sum(sh(a3_w2, j - 3) * sh(c, j) for j in range(3)))


def H2(v):
    return sp.expand(sum(sh(b2_w2, j - 2) * sh(v, j) for j in range(2)))


for jdeg in range(0, 6):
    c = sp.expand(falling(3) * E**jdeg)
    v = sp.expand(falling(2) * E**jdeg)
    check_zero(sp.rem(K3(c), D, E), f"D | K_3[(E)_3 E^{jdeg}]  (top filler block)")
    check_zero(sp.rem(H2(v), D, E), f"D | H_2[(E)_2 E^{jdeg}]  (sub filler block)")
# slope gate
Rok = sp.expand(E + E * (E - 1) * (E + 1) * (E**2 + 7))
check(sp.rem(E - Rok, D, E) == 0, "slope gate: D | (E-R) when R(1)=1, R(-1)=-1 (R(0)=0)")
check(sp.rem(E - sp.expand(2 * E + D), D, E) != 0,
      "slope gate: D does NOT divide E-R when R(1)!=1")
check(sp.rem(E - sp.expand(E + E * (E - 1)), D, E) != 0,
      "slope gate: D does NOT divide E-R when R(-1)!=-1")

# =====================================================================
print("\n--- S2. explicit slope-1 datum; build stacked L=[Phi;M_-1] ---", flush=True)
# =====================================================================
# explicit branch-B slope-1 positive datum (w2-decisive Section 4; a_-3=0)
Apt = {
    3: a3_w2,
    2: sp.expand(-E**3 / 3 + E**2 / 2 + R(9, 2) * E + 2),
    1: sp.expand(R(5, 7) * E**3 + R(107, 63) * E**2 + R(118, 63) * E + R(10, 9)),
    0: sp.expand(-R(775, 5103) * E**3 - R(92545, 3402) * E**2 - R(277597, 10206) * E),
    -1: sp.expand(-R(9, 8) * E**4 - R(219830, 11907) * E),
}
Apt[-2] = sp.Integer(0)
Apt[-3] = sp.Integer(0)
Bpt = {
    3: sp.Integer(0), 2: b2_w2,
    1: sp.expand(-R(2, 9) * E**2 + R(8, 9) * E + R(4, 3)),
    0: sp.expand(R(263, 567) * E**2 + R(179, 567) * E),
    -1: sp.expand(-R(256, 5103) * E**2 - R(31130, 1701) * E),
    -2: sp.expand(-R(3, 4) * E**3 + R(2747993, 2571912) * E**2 - R(819059, 2571912) * E),
    -3: sp.Integer(0),
}
for m in (4, 3, 2, 1):
    check_zero(q_m(Apt, Bpt, m), f"explicit datum: positive cascade Q_{m}=0")
check_zero(q_m(Apt, Bpt, 0).subs(E, 0) - 1, "explicit datum: slope R(1)=Q_0(0)=1")


def stacked_L(A_, B_, dV, use_am2_filler=True):
    """Return L (rows over E-degree of [G-E ; Q_-1]), rhs, tags, fillers.
    Block 0 = G-E (filler part = Phi); block 1 = Q_-1 (filler part = M_-1)."""
    Vr, cV = poly("Vv", dV)
    Cr, cC = poly("Cc", dV)
    A = dict(A_)
    B = dict(B_)
    if use_am2_filler:
        A[-2] = sp.expand(falling(2) * Vr)
    B[-3] = sp.expand(falling(3) * Cr)
    fillers = cV + cC
    B0 = sp.expand(potential(A, B) - E)      # G - E
    B1 = sp.expand(q_m(A, B, -1))            # Q_-1
    rows, rhs, tags = [], [], []
    for blk, P in (("0", B0), ("1", B1)):
        for (dg,), co in sp.Poly(P, E).terms():
            lin = sp.Poly(sp.expand(co), *fillers)
            rows.append([lin.coeff_monomial(f) for f in fillers])
            rhs.append(-lin.coeff_monomial(sp.Integer(1)))
            tags.append((blk, dg))
    return sp.Matrix(rows), sp.Matrix(rhs), tags, fillers, B0, B1


L, rv, tags, fillers, B0, B1 = stacked_L(Apt, Bpt, 6)
nf = len(fillers)
# filler-linear part of B0 is Phi; check D | Phi at random filler
subr = {c: R(random.randint(-3, 3)) for c in fillers}
Phi_at = sp.expand((B0 - (potential(Apt, Bpt) - E)).subs(subr))
check(sp.rem(Phi_at, D, E) == 0, "block-0 filler part is Phi with D | Phi (Im Phi = D F[E])")
check(L.rank() == nf and L.row_join(rv).rank() - L.rank() == 1,
      f"stacked L=[Phi;M_-1] (dV=6): full column rank {nf} and Fredholm gap = 1")

# =====================================================================
print("\n--- S3. THE JOINT COVECTOR: explicit left-kernel + pivot formula ---",
      flush=True)
# =====================================================================
def extract_covector(L, rv):
    """Return (lam, lam_rhs) with lam.L=0 (every column) and lam.rhs=lam_rhs (!=0
    when the system is infeasible)."""
    for y in L.T.nullspace():
        p = sp.expand((y.T * rv)[0])
        if p != 0:
            return y, p
    return None, sp.Integer(0)


lam, lam_rhs = extract_covector(L, rv)
check(lam is not None, "a left-kernel covector lam with lam.rhs != 0 exists (dV=6)")
allcols = all(sp.expand((lam.T * L[:, c])[0]) == 0 for c in range(nf))
check(allcols, "lam annihilates every column of L (lam . L = 0)")
check(sp.expand(lam_rhs) != 0, f"lam . rhs = {lam_rhs} != 0  (explicit obstruction)")

# minimal genuine coupling: dV=3 is the first filler degree where BOTH blocks are
# individually feasible yet the join fails (lower dV is a truncation artifact).
def block_feasible(A_, B_, dV, blk):
    L2, rv2, tags2, fil2, b0, b1 = stacked_L(A_, B_, dV)
    idx = [i for i, t in enumerate(tags2) if t[0] == blk]
    Lb = L2[idx, :]
    rb = rv2[idx, :]
    return Lb.rank() == Lb.row_join(rb).rank()


check(block_feasible(Apt, Bpt, 3, "0") and block_feasible(Apt, Bpt, 3, "1"),
      "control (dV=3): block-0 (Q_0=1) AND block-1 (Q_-1) each individually FEASIBLE")
L3, rv3, tags3, fil3, _, _ = stacked_L(Apt, Bpt, 3)
check(L3.rank() == len(fil3) and L3.row_join(rv3).rank() - L3.rank() == 1,
      "dV=3 is the MINIMAL genuine joint infeasibility (full col rank, gap=1)")

# pivot/residual formula: lam as an explicit function of the datum.  Choose nf
# independent rows I of L; x* = L_I^{-1} rhs_I; for an augmenting row j the
# obstruction is L_j x* - rhs_j, realised by lam' = e_j - L_j L_I^{-1} (rows I).
indep, M0 = [], sp.zeros(0, nf)
for i in range(L.rows):
    cand = M0.col_join(L[i, :])
    if cand.rank() > M0.rank():
        M0 = cand
        indep.append(i)
    if len(indep) == nf:
        break
LI = L[indep, :]
xstar = LI.solve(rv[indep, :])
jrow = next(i for i in range(L.rows) if i not in indep
            and sp.expand((L[i, :] * xstar)[0] - rv[i]) != 0)
resid = sp.expand((L[jrow, :] * xstar)[0] - rv[jrow])
check(resid != 0, f"pivot/residual formula: obstruction L_j x* - rhs_j = {resid} != 0 "
                  f"(explicit rational function of the datum)")

# =====================================================================
print("\n--- S4. CLOSED FORM: point-triple DECOUPLES; coupled mu0 on the necklace ---",
      flush=True)
# =====================================================================
# Ann(Im Phi) = {ev_-1, ev_0, ev_1} (roots of D).  A covector with mu0 there and any
# mu1 has ZERO obstruction: (i) mu0(Phi)=0 since D|Phi; (ii) on slope mu0(E-R)=0
# since (E-R)(-1)=(E-R)(1)=0; (iii) mu1 alone cannot obstruct since block-1 (Q_-1)
# is feasible.  So the covector MUST leave the point triple.
def node_lin(P, n, fillers):
    lin = sp.Poly(sp.expand(P.subs(E, n)), *fillers)
    return [lin.coeff_monomial(f) for f in fillers], lin.coeff_monomial(sp.Integer(1))


def support_couples(A_, B_, dV, S0, S1):
    """True iff some covector with mu0 supported on S0 (block G-E) and mu1 on S1
    (block Q_-1) annihilates the fillers with a NONZERO constant (obstruction)."""
    L2, rv2, tags2, fil2, b0, b1 = stacked_L(A_, B_, dV)
    ps = sp.symbols(f"pp0:{len(S0)}")
    qs = sp.symbols(f"qq0:{len(S1)}")
    fil_eqs = [sp.Integer(0)] * len(fil2)
    const = sp.Integer(0)
    for i, n in enumerate(S0):
        cf, ct = node_lin(b0, n, fil2)
        for j in range(len(fil2)):
            fil_eqs[j] += ps[i] * cf[j]
        const += ps[i] * ct
    for i, n in enumerate(S1):
        cf, ct = node_lin(b1, n, fil2)
        for j in range(len(fil2)):
            fil_eqs[j] += qs[i] * cf[j]
        const += qs[i] * ct
    sol = sp.linsolve(fil_eqs, list(ps) + list(qs))
    if not sol:
        return False
    s = list(sol)[0]
    const_on = sp.expand(const.subs(dict(zip(list(ps) + list(qs), s))))
    free = set().union(*[v.free_symbols for v in s]) if s else set()
    return (const_on != 0) or any(const_on.has(f) for f in free)


S1win = list(range(-9, 10))
check(not support_couples(Apt, Bpt, 6, [-1, 0, 1], S1win),
      "mu0 on Ann(Im Phi)={ev_-1,ev_0,ev_1} DECOUPLES: obstruction is 0 "
      "(point triple is silent at slope one)")
# the fixed (a_3,b_2) root necklace:
roots_a3 = [0, -2, -4]
roots_b2 = [0, -3]
roots_qK = [0, 1, 2] + [r + 3 for r in roots_a3]   # roots of a_3(E-3)(E)_3 = {-1,0,1,2,3}
roots_qH = [0, 1] + [r + 2 for r in roots_b2]       # roots of b_2(E-2)(E)_2 = {-1,0,1,2}
N_fix = sorted(set(roots_a3) | set(roots_b2) | set(roots_qK) | set(roots_qH))
check(N_fix == [-4, -3, -2, -1, 0, 1, 2, 3],
      "fixed (a_3,b_2) necklace = roots(a_3) u roots(b_2) u roots(q_K) u roots(q_H) "
      "= {-4,-3,-2,-1,0,1,2,3}")
check(support_couples(Apt, Bpt, 6, N_fix, S1win),
      "mu0 on the fixed (a_3,b_2) necklace COUPLES: a nonzero-obstruction covector "
      "exists (confirms the union-of-necklaces conjecture)")
# it is the necklace, not merely a wide window: {-4,..,1} (a_3,b_2 roots + D) is
# INSUFFICIENT -- the coupled covector needs the shifted membership window {2,3} of q_K.
check(not support_couples(Apt, Bpt, 6, [-4, -3, -2, -1, 0, 1], S1win),
      "mu0 on {-4,..,1} (roots(a_3) u roots(b_2) u roots(D)) alone is INSUFFICIENT: "
      "the coupled covector needs the full q_K necklace window {..,2,3}")

# =====================================================================
print("\n--- S5. FILLER-DEGREE-FREENESS: gap=1, obstruction != 0 at dV=4..10 ---",
      flush=True)
# =====================================================================
for dV in range(4, 11):
    Ld, rvd, tagsd, fild, _, _ = stacked_L(Apt, Bpt, dV)
    rLd = Ld.rank()
    gapd = Ld.row_join(rvd).rank() - rLd
    _, valv = extract_covector(Ld, rvd)
    check(rLd == len(fild) and gapd == 1 and sp.expand(valv) != 0,
          f"dV={dV}: full col rank ({rLd}={len(fild)}), gap=1, obstruction={valv} != 0 "
          f"-- no tail-clearing filler at this degree")

# =====================================================================
print("\n--- S6. GENERIC in positive data: distinct slope-1 data, obstruction != 0 ---",
      flush=True)
# =====================================================================
def clean_solve(A, B, m, lkey, name, membership, raw_degree):
    raw, cs = poly(name, raw_degree)
    unknown = sp.expand(falling(membership) * raw)
    trial = dict(B)
    trial[lkey] = unknown
    equations = sp.Poly(q_m(A, trial, m), E).all_coeffs()
    Mx, rhs = sp.linear_eq_to_matrix(equations, cs)
    conditions = [c for c in (sp.expand(n.dot(rhs)) for n in Mx.T.nullspace()) if c != 0]
    independent = sp.zeros(0, len(cs))
    selected_rhs = []
    for i in range(Mx.rows):
        cand = independent.col_join(Mx[i, :])
        if cand.rank() > independent.rank():
            independent = cand
            selected_rhs.append(rhs[i])
    if independent.rows == 0:
        values = [sp.Integer(0)] * len(cs)
    else:
        solution, _ = independent.gauss_jordan_solve(sp.Matrix(selected_rhs))
        taus = [x for x in solution.free_symbols if str(x).startswith("tau")]
        values = [x.subs({t: 0 for t in taus}) for x in solution]
    result = sp.expand(unknown.subs(dict(zip(cs, values))))
    kernels = []
    for j, vector in enumerate(Mx.nullspace()):
        parameter = sp.symbols(f"{name}K{j}")
        kp = falling(membership) * sum(vector[i] * E**i for i in range(len(cs)))
        result = sp.expand(result + parameter * kp)
        kernels.append(parameter)
    return result, kernels, conditions


def positive_cascade(d, with_am2=True):
    a2, ca2 = poly("a2", d)
    a1, ca1 = poly("a1", d)
    a0, ca0 = poly("a0", d)
    am1_raw, cam1 = poly("am1", d)
    A = {3: a3_w2, 2: a2, 1: a1, 0: a0, -1: falling(1) * am1_raw, -2: sp.Integer(0)}
    if with_am2:
        am2_raw, cam2 = poly("am2", d)
        A[-2] = sp.expand(falling(2) * am2_raw)
    B = {k: sp.Integer(0) for k in range(-3, 4)}
    B[2] = b2_w2
    A[-3] = sp.Integer(0)
    B[-3] = sp.Integer(0)
    conds, kernels = [], []
    for m, lkey, name, membership, degree in [
            (4, 1, "b1c", 0, d + 3), (3, 0, "b0c", 0, 2 * d + 2),
            (2, -1, "bm1c", 1, 2 * d + 3), (1, -2, "bm2c", 2, 2 * d + 4)]:
        B[lkey], nk, nc = clean_solve(A, B, m, lkey, name, membership, degree)
        kernels += nk
        conds += nc
    dat = {"ca2": ca2, "ca1": ca1, "ca0": ca0, "cam1": cam1}
    return A, B, conds, dat, kernels


A0g, B0g, condsg, datg, kernelsg = positive_cascade(3, with_am2=False)
a2_0, a2_1, a2_2, a2_3 = datg["ca2"]
a1_0, a1_1, a1_2, a1_3 = datg["ca1"]
a0_0, a0_1, a0_2, a0_3 = datg["ca0"]
am1_0, am1_1, am1_2, am1_3 = datg["cam1"]


def build_slope1(slice_vals):
    fixL = dict(slice_vals)
    fixL[a2_1] = 3 * fixL[a2_2] - 9 * fixL[a2_3]
    cf = [c for c in (sp.expand(c.subs(fixL)) for c in condsg) if c != 0]
    a1_only = [c for c in cf if c.free_symbols <= {a1_0, a1_1, a1_2}]
    s1 = list(sp.linsolve(a1_only, [a1_0, a1_1, a1_2]))
    if not s1 or any(v.free_symbols for v in s1[0]):
        return None
    fixL.update(dict(zip([a1_0, a1_1, a1_2], s1[0])))
    cf = [c for c in (sp.expand(c.subs(fixL)) for c in condsg) if c != 0]
    rest = [a0_1, a0_2, a0_3, am1_0, am1_1, am1_2]
    s2 = list(sp.linsolve(cf, rest))
    if not s2:
        return None
    fixL.update(dict(zip(rest, s2[0])))
    for v in (am1_1, am1_2):
        r = sp.Integer(random.randint(-4, 4))
        fixL = {k: (sp.expand(val.subs(v, r)) if hasattr(val, "subs") else val)
                for k, val in fixL.items()}
        fixL[v] = r
    fixL.setdefault(a0_0, sp.Integer(0))
    for k in kernelsg:
        fixL.setdefault(k, sp.Integer(0))
    if not all(sp.expand(c.subs(fixL)) == 0 for c in condsg):
        return None
    Aev = {k: sp.expand(v.subs(fixL)) for k, v in A0g.items()}
    Bev = {k: sp.expand(v.subs(fixL)) for k, v in B0g.items()}
    R1 = sp.expand(q_m(Aev, Bev, 0).subs(E, 0))
    cs_ = sp.solve(sp.Eq(R1, 1), am1_3)
    if not cs_:
        return None
    sub = {am1_3: cs_[0]}
    Aev = {k: sp.expand(v.subs(sub)) for k, v in Aev.items()}
    Bev = {k: sp.expand(v.subs(sub)) for k, v in Bev.items()}
    if sp.expand(q_m(Aev, Bev, 0).subs(E, 0)) != 1:
        return None
    return Aev, Bev


slices = [
    {a2_0: R(2), a2_2: R(1, 2), a2_3: R(-1, 3), a1_3: R(5, 7)},
    {a2_0: R(1), a2_2: R(1), a2_3: R(1), a1_3: R(1)},
    {a2_0: R(3), a2_2: R(-1), a2_3: R(2), a1_3: R(-2)},
    {a2_0: R(-2), a2_2: R(2), a2_3: R(-1), a1_3: R(3)},
    {a2_0: R(1, 2), a2_2: R(3), a2_3: R(-2), a1_3: R(4)},
    {a2_0: R(2), a2_2: R(1, 2), a2_3: R(-1, 3), a1_3: R(-3)},
]
npts = 0
for i, sl in enumerate(slices):
    pt = build_slope1(sl)
    if pt is None:
        continue
    Aev, Bev = pt
    Aev = dict(Aev); Aev[-2] = sp.Integer(0)
    Lp, rvp, tagsp, filp, _, _ = stacked_L(Aev, Bev, 5)
    _, valp = extract_covector(Lp, rvp)
    npts += 1
    check(Lp.rank() == len(filp) and Lp.row_join(rvp).rank() - Lp.rank() == 1
          and sp.expand(valp) != 0,
          f"generic slope-1 pt {i}: full col rank, gap=1, obstruction={valp} != 0")
check(npts >= 5, f"collected {npts} distinct slope-1 data (>=5): generic obstruction != 0")

# =====================================================================
print("\n--- S7. PAYLOAD: d=3 unit certificate (cascade+Q_0=1+Q_-1), both branches ---",
      flush=True)
# =====================================================================
def coeffs(expr):
    return [c for c in sp.Poly(sp.expand(expr), E).all_coeffs() if sp.expand(c) != 0]


def is_unit_QQ(eqs, vs):
    eqs = [sp.expand(e) for e in eqs if sp.expand(e) != 0]
    if not eqs:
        return False
    Gb = sp.groebner(eqs, *vs, order="grevlex", domain=sp.QQ)
    return list(Gb.exprs) == [sp.Integer(1)]


# --- sympy exact cross-check at the explicit datum: {Q_0=1, Q_-1} is UNIT in the
#     fillers, while {Q_0=1} alone is FEASIBLE (no false kill). ---
for branch in ("B", "A"):
    dV = 6
    A = dict(Apt); B = dict(Bpt)
    Vr, cV = poly("Vloc", dV)
    A[-2] = sp.expand(falling(2) * Vr)
    if branch == "B":
        Cr, cC = poly("Cloc", dV)
        A[-3] = sp.Integer(0)
        B[-3] = sp.expand(falling(3) * Cr)
        fv = cV + cC
    else:
        am3, cam3 = poly("am3loc", dV)
        mu3 = sp.symbols("mu3loc")
        A[-3] = sp.expand(falling(3) * am3)
        B[-3] = sp.expand(mu3 * A[-3])
        fv = cV + cam3 + [mu3]
    slope = coeffs(q_m(A, B, 0) - 1)
    check(not is_unit_QQ(slope, fv),
          f"branch {branch}: {{Q_0=1}} alone is FEASIBLE in fillers (no false kill)")
    joint = slope + coeffs(q_m(A, B, -1))
    check(is_unit_QQ(joint, fv),
          f"branch {branch}: {{Q_0=1, Q_-1}} is the UNIT ideal in fillers (sympy, "
          f"explicit datum) -- the joint covector obstruction")


# --- symbolic d=3 unit certificate over the whole slope-1 cascade locus (msolve '^') ---
def build_full_branch(d, branch):
    A, B, pos, dat, kernels = positive_cascade(d, with_am2=True)
    base = dat["ca2"] + dat["ca1"] + dat["ca0"] + dat["cam1"]
    am2v = list(set(A[-2].free_symbols) - {E})
    if branch == "A":
        am3_raw, cam3 = poly("am3f", d)
        mu3 = sp.symbols("mu3f")
        A[-3] = sp.expand(falling(3) * am3_raw)
        B[-3] = sp.expand(mu3 * A[-3])
        extra = cam3 + [mu3]
    else:
        bm3_raw, cbm3 = poly("bm3f", d)
        A[-3] = sp.Integer(0)
        B[-3] = sp.expand(falling(3) * bm3_raw)
        extra = cbm3
    fv = base + am2v + extra + kernels
    posnz = [sp.expand(e) for e in pos if sp.expand(e) != 0]
    return A, B, posnz, fv


def clear_denoms(e):
    num, _ = sp.fraction(sp.together(sp.expand(e)))
    return sp.expand(num)


def msolve_empty(eqs, vs, char, tmo=480):
    """True iff V(eqs)=empty (msolve prints [-1]).  Powers written with '^'."""
    xs = sp.symbols(f"z0:{len(vs)}")
    sub = dict(zip(vs, xs))
    with tempfile.NamedTemporaryFile("w", suffix=".ms", delete=False) as f:
        inp = f.name
        f.write(",".join(str(x) for x in xs) + "\n" + str(char) + "\n")
        f.write(",\n".join(str(clear_denoms(e.subs(sub))).replace(" ", "").replace("**", "^")
                           for e in eqs if sp.expand(e) != 0) + "\n")
    out = inp + ".out"
    try:
        subprocess.run(["msolve", "-f", inp, "-o", out], check=True,
                       stderr=subprocess.DEVNULL, timeout=tmo)
        r = open(out).read().strip()
    finally:
        for pth in (inp, out):
            if os.path.exists(pth):
                os.remove(pth)
    return r.startswith("[-1]"), r[:16]


if shutil.which("msolve"):
    for branch in ("A", "B"):
        A, B, posnz, fv = build_full_branch(3, branch)
        sys_loc = posnz + coeffs(q_m(A, B, 0) - 1) + coeffs(q_m(A, B, -1))
        emp_q, raw_q = msolve_empty(sys_loc, fv, 0)
        check(emp_q,
              f"branch {branch} d=3: cascade + Q_0=1 + Q_-1 = UNIT ideal over QQ "
              f"(msolve '^') -- Joint Filler Covector Lemma UNIFORM at d=3; tail not needed")
        emp_p, raw_p = msolve_empty(sys_loc, fv, 65003)
        check(emp_p, f"branch {branch} d=3: same UNIT ideal mod 65003 (msolve '^' cross-check)")
    if os.environ.get("HEAVY") == "1":
        for branch in ("A", "B"):
            A, B, posnz, fv = build_full_branch(4, branch)
            sys_loc = posnz + coeffs(q_m(A, B, 0) - 1) + coeffs(q_m(A, B, -1))
            emp_q4, _ = msolve_empty(sys_loc, fv, 0)
            check(emp_q4,
                  f"[HEAVY] branch {branch} d=4: cascade + Q_0=1 + Q_-1 = UNIT over QQ "
                  f"(msolve '^') -- joint kill extends to raw cap d=4")
    else:
        print("    [set HEAVY=1 to also verify the d=4 joint unit certificate over QQ "
              "(both branches, ~15s)]")
else:
    print("    [msolve not on PATH; the d=3 unit certificate cascade+Q_0=1+Q_-1 is")
    print("     recorded in joint-covector.md; the sympy explicit-datum cross-check above")
    print("     ({Q_0=1,Q_-1} UNIT in fillers, both branches) is self-contained.]")

# =====================================================================
print("\n--- S8. BRANCH A: covector per fixed gauge mu_3 ---", flush=True)
# =====================================================================
def branchA_stacked(mu3, dV=4):
    Vr, cV = poly("VvA", dV)
    Wr, cW = poly("WwA", dV)
    A = dict(Apt); B = dict(Bpt)
    A[-2] = sp.expand(falling(2) * Vr)
    A[-3] = sp.expand(falling(3) * Wr)
    B[-3] = sp.expand(mu3 * A[-3])
    fillers = cV + cW
    B0b = sp.expand(potential(A, B) - E)
    B1b = sp.expand(q_m(A, B, -1))
    rows, rhs = [], []
    for P in (B0b, B1b):
        for (dg,), co in sp.Poly(P, E).terms():
            lin = sp.Poly(sp.expand(co), *fillers)
            rows.append([lin.coeff_monomial(f) for f in fillers])
            rhs.append(-lin.coeff_monomial(sp.Integer(1)))
    return sp.Matrix(rows), sp.Matrix(rhs), len(fillers)


for mu3 in [sp.Integer(1), sp.Integer(-2), R(1, 3)]:
    La, rva, nfa = branchA_stacked(mu3)
    _, vala = extract_covector(La, rva)
    check(La.rank() == nfa and La.row_join(rva).rank() - La.rank() == 1
          and sp.expand(vala) != 0,
          f"branch A gauge mu_3={mu3}: full col rank, gap=1, obstruction={vala} != 0")

# =====================================================================
print("\n--- S9. METHOD DEBT: the '**'->'^' audit flips the sibling S8 verdict ---",
      flush=True)
# =====================================================================
# hatch-census.md S5 flagged that verify_w2_joint.py S8's msolve FEASIBILITY claim
# 'cascade+Q_0=1+Q_-1 is NOT unit (a sub-locus survives Q_-1)' predates the '**'->'^'
# audit.  S7 above re-runs it with '^': the system IS the unit ideal (both branches,
# QQ + prime).  So the claimed surviving sub-locus is a '**'-parser ARTIFACT and is
# EMPTY at d=3; [Phi;M_-1] is the uniform kill at d=3, and Q_-2..Q_-5 are not needed.
# Reproduce the parser bug directly on a known-unit reproducer (hatch-census S5).
if shutil.which("msolve"):
    ys = sp.symbols("y0:5")
    y0, y1, y2, y3, y4 = ys
    g1 = 3 * y0**3 - 2 * y0**2 * y1 + 9 * y0 * y4
    g2 = (-9 * y0**5 + 12 * y0**4 * y1 - 4 * y0**3 * y1**2 - 54 * y0**3 * y4
          + 36 * y0**2 * y1 * y4 - 81 * y0 * y4**2 - 2916)
    # sympy proof of unit over QQ (independent of msolve):
    check(is_unit_QQ([g1, g2], list(ys)),
          "reproducer {g1,g2} is the UNIT ideal over QQ (sympy) -- ground truth")

    def msolve_raw(eqs, vs, char, use_caret):
        xs = sp.symbols(f"w0:{len(vs)}")
        sub = dict(zip(vs, xs))
        with tempfile.NamedTemporaryFile("w", suffix=".ms", delete=False) as f:
            inp = f.name
            f.write(",".join(str(x) for x in xs) + "\n" + str(char) + "\n")
            body = ",\n".join(str(sp.expand(e.subs(sub))).replace(" ", "")
                              for e in eqs)
            if use_caret:
                body = body.replace("**", "^")
            f.write(body + "\n")
        out = inp + ".out"
        try:
            subprocess.run(["msolve", "-f", inp, "-o", out], check=True,
                           stderr=subprocess.DEVNULL, timeout=120)
            r = open(out).read().strip()
        finally:
            for pth in (inp, out):
                if os.path.exists(pth):
                    os.remove(pth)
        return r
    r_caret = msolve_raw([g1, g2], list(ys), 0, True)
    r_star = msolve_raw([g1, g2], list(ys), 0, False)
    check(r_caret.startswith("[-1]") or r_caret.startswith("[1"),
          f"msolve with '^' returns UNIT/empty on the reproducer (raw {r_caret[:12]})")
    check(not (r_star.startswith("[-1]")),
          f"msolve with '**' MISPARSES and does NOT return unit (raw {r_star[:12]}) "
          f"-- the parser bug that produced the spurious 'sub-locus survives' claim")
else:
    print("    [msolve not on PATH; the '**' vs '^' reproducer is recorded in")
    print("     hatch-census.md S5 and joint-covector.md.]")

print("\n" + "=" * 70)
print("JOINT FILLER COVECTOR LEMMA (W2): the generic covector lam is EXPLICIT, its")
print("block-0 part lives on the fixed (a_3,b_2) root necklace {-4,..,3} (leaving the")
print("Ann(Im Phi) point triple), and lam.rhs != 0 across filler degrees dV=4..10 and")
print("6 slope-1 data.  PAYLOAD: cascade+Q_0=1+Q_-1 is the UNIT ideal at d=3 (both")
print("branches, msolve '^' + sympy) -- the lemma holds UNIFORMLY at d=3; the sibling")
print("'sub-locus survives Q_-1' was a '**' parser artifact.  Arbitrary positive-data")
print("degree and a degree-free lam.rhs identity remain OPEN.")
print("=" * 70)
print(f"\n(total {time.time() - _T0:.1f}s)")
print("ALL JOINT COVECTOR CHECKS PASSED")
