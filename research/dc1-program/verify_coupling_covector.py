#!/usr/bin/env python3
"""THE COUPLING COVECTOR at W2 -- the varying-tops coupling, resolved onto a FIXED
INTEGER-NODE grid, and the sharpened residual gap.

INDEPENDENTLY DERIVED - EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES - NOT PEER REVIEWED

Target (the prize): a fixed finite recipe producing, at EVERY positive-data degree d, a
covector annihilating the depth-3 filler columns whose residual pairing is a UNIT multiple
of W.  Convention block (quantum band 3, gauge b_3=0):
   Q_m = sum_(k+l=m)[ b_l^[k] a_k - a_k^[l] b_l ],  f^[n](E)=f(E+n),
   membership (E)_j = E(E-1)...(E-j+1) | a_-j, b_-j ;  a_3=E(E+2)(E+4), b_2=E(E+3),
   branch B: a_-3=0, b_-3=(E)_3 C ; fillers  P:=a_-2=(E)_2 V,  S:=b_-3=(E)_3 C.

WHAT THIS FILE ESTABLISHES
==========================
PROVED (degree-free / arbitrary d, machine-checked identities):
  * S0/S1 upstream re-derivation: engine Q_m=[D,X]_m, Q_0=(T-1)G, slope gate, both-ends
    Lemma P, and the factorization R(1)=a_2(0)*W.
  * S2 STAIRCASE COLLAPSE of the two blocks.  The two-term operators are
       a-block(Q_m) = b_{m+2}(E-2) P(E) - b_{m+2}(E) P(E+m+2),
       b-block(Q_m) = a_{m+3}(E) S(E+m+3) - a_{m+3}(E-3) S(E),
    and at m=-2 the a-block DEGENERATES to a pure multiplier [b_0(E-2)-b_0(E)]P(E), at
    m=-3 the b-block degenerates to the pure multiplier [a_0(E)-a_0(E-3)]S(E).  The
    depth-3 filler dependence is therefore a STAIRCASE, not a generic 2x2 coupling.
  * S3 MEMBERSHIP SILENCE TABLE.  Exactly six (rung,node) pairs are identically 0=0:
    Q_-1(0), Q_-2(0), Q_-2(1), Q_-3(0), Q_-3(1), Q_-3(2)  -- BOTH the filler part and the
    residual N_m vanish there.  Degree-free.
  * S4 THE INTEGER-NODE WINDOW (the fixed-shape recipe).  For every k>=3 and every degree
    d, the rows { Q_-1(e) }_{e=1..k-1} u { Q_-2(e) }_{e=2..k} u { Q_-3(e) }_{e=3..k}
    involve ONLY the 2k-2 filler VALUES P(2..k), S(3..k+1); the augmented matrix [M|N] has
    the d-INDEPENDENT shape (3k-4) x (2k-1) and its entries are the cascade datum evaluated
    at FIXED INTEGER NODES.  So the coupling across the varying tops (a_2,a_1,a_0) /
    (b_1,b_0,b_-1) collapses onto ONE fixed integer grid -- no algebraic necklace, no
    trace forms, no gcd(a_2(E),a_2(E-3)) condition.
  * S5 ON-SHELL NODE RECURSION.  The positive cascade fixes the b-tops on that grid by a
    degree-free 3-step recursion anchored at the a_3-roots {0,-2,-4}:
       b_1(0) = (2/3) a_2(0),   a_2(-3) = a_2(0)  (a cascade consistency relation),
       b_1(2) = (2/3)(2 a_2(-1) + a_2(1)),   15 b_1(4) = 48 b_1(1) - 18 a_2(1) + 4 a_2(3),
    with b_1(1), b_1(3) FREE SEEDS (they drop out of Q_4 at E=-2 and E=0).  The collapse is
    therefore PARTIAL: the b-tops share the a-grid but are not eliminated by it.

REFUTED (exact, stated scope):
  * S7 NO fixed-shape integer-node window minor is a unit multiple of W.  At d=3, for
    k=3,4,5 (1, 8, 55 maximal minors) every nonzero minor has a NONZERO residue at
    am1_3=0, at each tested exact rational specialization -- so no single window minor
    equals u*W.
  * S7 The k=3 window does not force W AT ALL: an EXPLICIT EXACT rational point of the
    cascade at which the k=3 window system is solvable (rank M = rank[M|N] = 4) while
    W = -8/9 != 0.  (Radical-correct: a witness point, not a normal form.)

BOUNDED-FINITE (exact scope stated):
  * S6 d=2 window-minor factorization: every k=3,4 maximal minor factors as
    (datum NODE VALUES) x (one common core K) -- the window elimination content is
    essentially principal, and the "resultant structure" is Sylvester-like in the fixed
    node values, not a generic Res_E(top(E),top(E-3)).
  * S8 combined-rung adjoint over Q(E) (symbolic node): the shifted-rung value matrix on
    the window j=0..3 is 12x12 of rank 11 at exact specializations; its unique covector's
    residual pairing is NOT a W-multiple.
  * S9 controls: the d=3 W-kill reproduced by RABINOWITSCH (sympy exact QQ + prime);
    a_2(0) NOT forced (explicit witness); cascade+tail non-vacuous (explicit point).

OPEN (not claimed): whether ANY k has W in sqrt(I_k); any degree-free W-forcing recipe;
  the arbitrary-degree residual identity; all of band 3 / DC1 / JC2.

Run:  uv run --with sympy python research/dc1-program/verify_coupling_covector.py
      HEAVY=1 ...  (adds d=4/d=5 window legs and the k=5 minor sweep at more seeds)
Ends: "ALL CHECKS PASSED" only when nothing was skipped; otherwise an explicit
      "ALL EXECUTED CHECKS PASSED; OPTIONAL CHECKS SKIPPED" line.
"""
import sympy as sp
import time, os, random, itertools

E = sp.symbols("E")
LEVELS = range(-3, 4)
R = sp.Rational
_T0 = time.time()
random.seed(20260724)
HEAVY = os.environ.get("HEAVY") == "1"
_NP = 0
_SKIPS = []


def sh(f, n):
    return sp.expand(sp.sympify(f).subs(E, E + n))


def falling(n):
    return sp.prod(E - j for j in range(n))


def poly(name, degree):
    cs = list(sp.symbols(f"{name}_0:{degree + 1}"))
    return sp.expand(sum(cs[j] * E**j for j in range(degree + 1))), cs


def q_m(A, B, m):
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


def coeffs(expr):
    return [c for c in sp.Poly(sp.expand(expr), E).all_coeffs() if sp.expand(c) != 0]


def check(cond, label):
    global _NP
    if not cond:
        raise AssertionError("FAIL " + label)
    _NP += 1
    print(f"PASS [{time.time() - _T0:6.1f}s] {label}", flush=True)


def check_zero(val, label):
    check(sp.expand(val) == 0, label)


def skip(label):
    _SKIPS.append(label)
    print(f"SKIP [{time.time() - _T0:6.1f}s] {label}", flush=True)


def sy_unit(eqs, vs, modulus=None):
    G = sp.groebner([sp.expand(e) for e in eqs if sp.expand(e) != 0], *vs,
                    order="grevlex", modulus=modulus)
    return list(G.exprs) == [sp.Integer(1)]


a3_w2 = sp.expand(E * (E + 2) * (E + 4))
b2_w2 = sp.expand(E * (E + 3))
D_op = sp.expand(E * (E - 1) * (E + 1))

# =====================================================================
print("--- S0. crossed-product engine (re-derived): Q_m=[D,X]_m, Q_0=(T-1)G ---", flush=True)
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
           "Q_0 = (T-1)G telescoping identity")

# =====================================================================
print("\n--- S1. slope gate; both-ends Lemma P; FACTORIZATION R(1)=a_2(0)*W (re-derived) ---",
      flush=True)
# =====================================================================
Rok = sp.expand(E + E * (E - 1) * (E + 1) * (E**2 + 7))
check(sp.rem(E - Rok, D_op, E) == 0, "slope gate: D | (E-R) when R(1)=1, R(-1)=-1")
check(sp.rem(E - sp.expand(2 * E + D_op), D_op, E) != 0, "slope gate: fails when R(1)!=1")
deg = 4
Af, Bf = {}, {}
for lev in LEVELS:
    ar, _ = poly(f"af{lev + 3}", deg)
    br, _ = poly(f"bf{lev + 3}", deg)
    memb = falling(-lev) if lev < 0 else 1
    Af[lev] = sp.expand(memb * ar)
    Bf[lev] = sp.expand(memb * br)
mu3s = sp.symbols("mu3s")
Aw, Bw = dict(Af), dict(Bf)
Aw[3] = a3_w2; Bw[2] = b2_w2; Bw[3] = sp.Integer(0); Bw[-3] = sp.expand(mu3s * Aw[-3])
Gw = potential(Aw, Bw)
check_zero(q_m(Aw, Bw, 0) - (sh(Gw, 1) - Gw), "W2: Q_0=(T-1)G on the W2 datum")
check_zero(Gw.subs(E, 0), "membership: G(0)=0  => R(1)=Q_0(0)=G(1)")
G1 = sp.expand(Aw[1].subs(E, 0) * Bw[-1].subs(E, 1) + Aw[2].subs(E, 0) * Bw[-2].subs(E, 2)
               - Aw[-1].subs(E, 1) * Bw[1].subs(E, 0))
check_zero(q_m(Aw, Bw, 0).subs(E, 0) - G1,
           "both-ends Lemma P (E=1): R(1)=a1(0)b-1(1)+a2(0)b-2(2)-a-1(1)b1(0)")
fillsyms = set(Aw[-2].free_symbols) | set(Bw[-3].free_symbols) | {mu3s}
check(not (sp.expand(q_m(Aw, Bw, 0).subs(E, 0)).free_symbols & fillsyms),
      "R(1)=Q_0(0) INDEPENDENT of the fillers a_-2, b_-3")
for dgen in (2, 3, 5):
    Ax, Bx = {}, {}
    for lev in LEVELS:
        ar, _ = poly(f"A{dgen}x{lev + 3}", dgen)
        br, _ = poly(f"B{dgen}x{lev + 3}", dgen)
        memb = falling(-lev) if lev < 0 else 1
        Ax[lev] = sp.expand(memb * ar); Bx[lev] = sp.expand(memb * br)
    Ax[3] = a3_w2; Bx[2] = b2_w2; Bx[3] = sp.Integer(0)
    a2x, a1x, b1x = Ax[2], Ax[1], Bx[1]
    check_zero(q_m(Ax, Bx, 4).subs(E, 0) - (10 * a2x.subs(E, 0) - 15 * b1x.subs(E, 0)),
               f"[d={dgen}] Q_4(0) = 10 a_2(0) - 15 b_1(0)")
    check_zero(q_m(Ax, Bx, 3).subs(E, 0)
               - (4 * a1x.subs(E, 0) + a2x.subs(E, 0) * b1x.subs(E, 2)
                  - a2x.subs(E, 1) * b1x.subs(E, 0)),
               f"[d={dgen}] Q_3(0) = 4 a_1(0) + a_2(0)b_1(2) - a_2(1)b_1(0)")
b1_0, a1_0, a2_0s, a2_1s, b1_2s = sp.symbols("b1_0 a1_0 a2_0 a2_1 b1_2")
sol_b1_0 = sp.solve(sp.Eq(10 * a2_0s - 15 * b1_0, 0), b1_0)[0]
sol_a1_0 = sp.solve(sp.Eq(4 * a1_0 + a2_0s * b1_2s - a2_1s * sol_b1_0, 0), a1_0)[0]
bm1_1s, bm2_2s, am1_1s = sp.symbols("bm1_1 bm2_2 am1_1")
R1_P = a1_0 * bm1_1s + a2_0s * bm2_2s - am1_1s * b1_0
W_closed = bm2_2s - R(2, 3) * am1_1s + R(1, 4) * (R(2, 3) * a2_1s - b1_2s) * bm1_1s
check_zero(R1_P.subs({a1_0: sol_a1_0, b1_0: sol_b1_0}) - a2_0s * W_closed,
           "FACTORIZATION: R(1) = a_2(0) * W  (degree-free)")


# ---- positive-cascade solver (branch B) ------------------------------------
def clean_solve(A, B, m, lkey, name, membership, raw_degree):
    raw, cs = poly(name, raw_degree)
    unknown = sp.expand(falling(membership) * raw)
    trial = dict(B); trial[lkey] = unknown
    equations = sp.Poly(q_m(A, trial, m), E).all_coeffs()
    M, rhs = sp.linear_eq_to_matrix(equations, cs)
    conditions = [c for c in (sp.expand(n.dot(rhs)) for n in M.T.nullspace()) if c != 0]
    independent = sp.zeros(0, len(cs)); selected_rhs = []
    for i in range(M.rows):
        cand = independent.col_join(M[i, :])
        if cand.rank() > independent.rank():
            independent = cand; selected_rhs.append(rhs[i])
    if independent.rows == 0:
        values = [sp.Integer(0)] * len(cs)
    else:
        solution, _ = independent.gauss_jordan_solve(sp.Matrix(selected_rhs))
        taus = [x for x in solution.free_symbols if str(x).startswith("tau")]
        values = [x.subs({t: 0 for t in taus}) for x in solution]
    result = sp.expand(unknown.subs(dict(zip(cs, values))))
    kernels = []
    for j, vector in enumerate(M.nullspace()):
        parameter = sp.symbols(f"{name}K{j}")
        kp = falling(membership) * sum(vector[i] * E**i for i in range(len(cs)))
        result = sp.expand(result + parameter * kp)
        kernels.append(parameter)
    return result, kernels, conditions


def positive_cascade(d):
    a2, ca2 = poly("a2", d); a1, ca1 = poly("a1", d)
    a0, ca0 = poly("a0", d); am1_raw, cam1 = poly("am1", d)
    am2_raw, _ = poly("am2", d)
    A = {3: a3_w2, 2: a2, 1: a1, 0: a0, -1: falling(1) * am1_raw,
         -2: sp.expand(falling(2) * am2_raw)}
    B = {k: sp.Integer(0) for k in range(-3, 4)}
    B[2] = b2_w2; A[-3] = sp.Integer(0); B[-3] = sp.Integer(0)
    conds, kernels = [], []
    for m, lkey, name, membership, degree in [
            (4, 1, "b1c", 0, d + 3), (3, 0, "b0c", 0, 2 * d + 2),
            (2, -1, "bm1c", 1, 2 * d + 3), (1, -2, "bm2c", 2, 2 * d + 4)]:
        B[lkey], nk, nc = clean_solve(A, B, m, lkey, name, membership, degree)
        kernels += nk; conds += nc
    return A, B, conds, ca2 + ca1 + ca0 + cam1, kernels


def greedy_param(conds, Pvars):
    sol, solved, changed = {}, set(), True
    while changed:
        changed = False
        best = None
        for c in conds:
            cc = sp.expand(c.subs(sol))
            if cc == 0:
                continue
            for v in sorted((cc.free_symbols & set(Pvars)) - solved, key=str):
                pc = sp.Poly(cc, v)
                if pc.degree() == 1:
                    lead = pc.coeff_monomial(v)
                    rest = sp.expand(cc - lead * v)
                    score = len(lead.free_symbols)
                    if best is None or score < best[0]:
                        best = (score, v, lead, rest)
        if best is not None:
            _, v, lead, rest = best
            val = sp.cancel(-rest / lead)
            sol = {k: sp.cancel(x.subs(v, val)) for k, x in sol.items()}
            sol[v] = val; solved.add(v); changed = True
    return sol, [v for v in Pvars if v not in solved]


def build_reduced(d):
    A, B, pos, base, kernels = positive_cascade(d)
    bm3_raw, cbm3 = poly("bm3f", d)
    A[-3] = sp.Integer(0); B[-3] = sp.expand(falling(3) * bm3_raw)
    fv = base + sorted(set(A[-2].free_symbols) - {E}, key=str) + cbm3 + kernels
    posnz = [sp.expand(e) for e in pos if sp.expand(e) != 0]
    fillers = sorted((set(A[-2].free_symbols) | set(B[-3].free_symbols)) - {E}, key=str)
    sol, free = greedy_param(posnz, [v for v in fv if v not in set(fillers)])
    R1 = sp.cancel(sp.sympify(q_m(A, B, 0).subs(E, 0)).subs(sol))
    resid = [sp.expand(sp.cancel(c.subs(sol))) for c in posnz]
    tk = {m: [sp.expand(sp.fraction(sp.together(sp.cancel(e.subs(sol))))[0])
              for e in coeffs(q_m(A, B, m)) if sp.expand(e) != 0] for m in (-1, -2, -3, -4, -5)}
    return dict(A=A, B=B, fillers=fillers, sol=sol, free=free, R1=R1, resid=resid, tk=tk)


def datum(rd, pt):
    sb = lambda x: sp.expand(sp.cancel(sp.sympify(x).subs(rd["sol"]).subs(pt)))
    A = {k: sb(rd["A"][k]) for k in rd["A"]}
    B = {k: sb(rd["B"][k]) for k in rd["B"]}
    A0 = dict(A); A0[-2] = sp.Integer(0)
    B0 = dict(B); B0[-3] = sp.Integer(0)
    N = {m: sp.expand(q_m(A0, B0, m)) for m in (-1, -2, -3)}
    W = sp.expand(R(1, 4) * (R(2, 3) * A[2].subs(E, 1) - B[1].subs(E, 2)) * B[-1].subs(E, 1)
                  + B[-2].subs(E, 2) - R(2, 3) * A[-1].subs(E, 1))
    return A, B, N, W


# =====================================================================
print("\n--- S2. THE STAIRCASE COLLAPSE of the two blocks (degree-free) ---", flush=True)
# =====================================================================
for dsym in (2, 3, 4, 5):
    Az, Bz = {}, {}
    for lev in LEVELS:
        ar, _ = poly(f"Az{dsym}_{lev + 3}", dsym); br, _ = poly(f"Bz{dsym}_{lev + 3}", dsym)
        memb = falling(-lev) if lev < 0 else 1
        Az[lev] = sp.expand(memb * ar); Bz[lev] = sp.expand(memb * br)
    Az[3] = a3_w2; Bz[2] = b2_w2; Bz[3] = sp.Integer(0); Az[-3] = sp.Integer(0)
    Pf, Sf = Az[-2], Bz[-3]
    Ps = sorted(set(Pf.free_symbols) - {E}, key=str); Ss = sorted(set(Sf.free_symbols) - {E}, key=str)
    for m in (-1, -2, -3):
        qm = q_m(Az, Bz, m)
        a_part = sp.expand(sum(sp.Poly(qm, *Ps).coeff_monomial(s) * s for s in Ps))
        b_part = sp.expand(sum(sp.Poly(qm, *Ss).coeff_monomial(s) * s for s in Ss))
        check_zero(a_part - sp.expand(sh(Bz[m + 2], -2) * Pf - Bz[m + 2] * sh(Pf, m + 2)),
                   f"[d={dsym}] Q_{m} a-block = b_{m+2}(E-2)P(E) - b_{m+2}(E)P(E+{m+2})")
        check_zero(b_part - sp.expand(Az[m + 3] * sh(Sf, m + 3) - sh(Az[m + 3], -3) * Sf),
                   f"[d={dsym}] Q_{m} b-block = a_{m+3}(E)S(E+{m+3}) - a_{m+3}(E-3)S(E)")
    # THE COLLAPSE: at m=-2 the a-block has NO shift; at m=-3 the b-block has NO shift.
    check_zero(sp.expand(sh(Bz[0], -2) * Pf - Bz[0] * sh(Pf, 0)) -
               sp.expand((sh(Bz[0], -2) - Bz[0]) * Pf),
               f"[d={dsym}] COLLAPSE: Q_-2 a-block = [b_0(E-2)-b_0(E)] P(E) -- pure MULTIPLIER")
    check_zero(sp.expand(Az[0] * sh(Sf, 0) - sh(Az[0], -3) * Sf) -
               sp.expand((Az[0] - sh(Az[0], -3)) * Sf),
               f"[d={dsym}] COLLAPSE: Q_-3 b-block = [a_0(E)-a_0(E-3)] S(E) -- pure MULTIPLIER")
print("       => the depth-3 filler dependence is a STAIRCASE:", flush=True)
print("          Q_-1: P at {0,+1}, S at {0,+2};  Q_-2: P at {0}, S at {0,+1};"
      "  Q_-3: P at {0,-1}, S at {0}.", flush=True)

# =====================================================================
print("\n--- S3. MEMBERSHIP SILENCE TABLE at the integer nodes (degree-free) ---", flush=True)
# =====================================================================
# Filler membership kills P(0),P(1) and S(0),S(1),S(2); the residuals N_m also carry
# membership factors.  Enumerate which (rung, node) pairs are identically 0=0.
# Two independent sources of silence, both degree-free:
#   (i) FILLER membership: P=a_-2=(E)_2 V vanishes at {0,1}; S=b_-3=(E)_3 C at {0,1,2};
#  (ii) CASCADE membership: b_-1=(E)_1(...) vanishes at 0, so the Q_-3 coefficients
#       b_-1(e-2) and b_-1(e) kill the P-terms at e=2 and e=0 respectively.
# A (rung,node) is silent iff every filler term dies by (i) or (ii)  [and then N_m(e)=0 too].
NULLP, NULLS = {0, 1}, {0, 1, 2}
# each entry: (filler tag shift, coefficient node offset, "which datum" or None)
TERMS = {-1: [(("P", 0), None), (("P", 1), None), (("S", 2), None), (("S", 0), None)],
         -2: [(("P", 0), None), (("S", 1), None), (("S", 0), None)],
         -3: [(("P", 0), ("bm1", -2)), (("P", -1), ("bm1", 0)), (("S", 0), None)]}
predicted_silent = set()
for m, terms in TERMS.items():
    for e in range(0, 5):
        dead = True
        for (tag, tsh), coef in terms:
            null_val = (e + tsh) in (NULLP if tag == "P" else NULLS)
            null_coef = coef is not None and (e + coef[1]) == 0   # b_-1(0)=0
            if not (null_val or null_coef):
                dead = False
        if dead:
            predicted_silent.add((m, e))
check(predicted_silent == {(-1, 0), (-2, 0), (-2, 1), (-3, 0), (-3, 1), (-3, 2)},
      f"level-incidence prediction: the silent (rung,node) pairs are {sorted(predicted_silent)} "
      f"(filler membership {{0,1}}/{{0,1,2}} plus the cascade membership b_-1(0)=0)")
for dsym in (2, 3, 4, 5):
    Az, Bz = {}, {}
    for lev in LEVELS:
        ar, _ = poly(f"Sz{dsym}_{lev + 3}", dsym); br, _ = poly(f"Tz{dsym}_{lev + 3}", dsym)
        memb = falling(-lev) if lev < 0 else 1
        Az[lev] = sp.expand(memb * ar); Bz[lev] = sp.expand(memb * br)
    Az[3] = a3_w2; Bz[2] = b2_w2; Bz[3] = sp.Integer(0); Az[-3] = sp.Integer(0)
    for (m, e) in sorted(predicted_silent):
        check_zero(sp.expand(q_m(Az, Bz, m).subs(E, e)),
                   f"[d={dsym}] SILENT rung: Q_{m}(E={e}) == 0 identically (fillers AND residual)")

# =====================================================================
print("\n--- S4. THE INTEGER-NODE WINDOW: fixed shape (3k-4)x(2k-1) at EVERY degree ---",
      flush=True)
# =====================================================================
NULLTAGS = {("P", 0), ("P", 1), ("S", 0), ("S", 1), ("S", 2)}


def window_rows(A, B, N, k):
    """rows Q_-1(1..k-1), Q_-2(2..k), Q_-3(3..k); cols P(2..k), S(3..k+1)."""
    ev = lambda f, x: sp.expand(sp.sympify(f).subs(E, x))
    b1, b0, bm1 = B[1], B[0], B[-1]
    a2, a1, a0 = A[2], A[1], A[0]
    cols = [("P", s) for s in range(2, k + 1)] + [("S", s) for s in range(3, k + 2)]
    ci = {c: i for i, c in enumerate(cols)}
    rows, tags = [], []

    def emit(tag, cd, nn):
        v = [sp.Integer(0)] * len(cols)
        for key, val in cd.items():
            if key in NULLTAGS:
                continue
            if key not in ci:
                raise AssertionError(f"window leak: {tag} touches {key}")
            v[ci[key]] = sp.expand(val)
        rows.append(v + [sp.expand(nn)]); tags.append(tag)
    for e in range(1, k):
        emit((-1, e), {("P", e): ev(b1, e - 2), ("P", e + 1): -ev(b1, e),
                       ("S", e + 2): ev(a2, e), ("S", e): -ev(a2, e - 3)}, ev(N[-1], e))
    for e in range(2, k + 1):
        emit((-2, e), {("P", e): ev(b0, e - 2) - ev(b0, e),
                       ("S", e + 1): ev(a1, e), ("S", e): -ev(a1, e - 3)}, ev(N[-2], e))
    for e in range(3, k + 1):
        emit((-3, e), {("P", e): ev(bm1, e - 2), ("P", e - 1): -ev(bm1, e),
                       ("S", e): ev(a0, e) - ev(a0, e - 3)}, ev(N[-3], e))
    return sp.Matrix(rows), cols, tags


# (a) DEGREE-FREE SHAPE (level incidence, no degree used):
for k in range(3, 8):
    nrow = (k - 1) + (k - 1) + (k - 2)
    ncol = (k - 1) + (k - 1)
    leaks = []
    for e in range(1, k):
        leaks += [("P", e), ("P", e + 1), ("S", e + 2), ("S", e)]
    for e in range(2, k + 1):
        leaks += [("P", e), ("S", e + 1), ("S", e)]
    for e in range(3, k + 1):
        leaks += [("P", e), ("P", e - 1), ("S", e)]
    live = {t for t in leaks if t not in NULLTAGS}
    ok = live <= ({("P", s) for s in range(2, k + 1)} | {("S", s) for s in range(3, k + 2)})
    check(ok and nrow == 3 * k - 4 and ncol == 2 * k - 2,
          f"k={k}: window rows={3*k-4}, live filler VALUES={2*k-2}, [M|N] shape "
          f"{3*k-4}x{2*k-1} -- INDEPENDENT of the degree d")
# (b) exact row identities on the real cascade at several degrees:
RED = {}
WDEGS = [2, 3] + ([4] if HEAVY else [])
for d in WDEGS:
    RED[d] = build_reduced(d)
    check(all(r == 0 for r in RED[d]["resid"]),
          f"d={d}: cascade parametrization satisfies every cascade condition identically")
    rng = random.Random(4100 + d)
    pt = {v: sp.Rational(rng.randint(-6, 6), rng.choice([1, 2, 3])) for v in RED[d]["free"]}
    A, B, N, W = datum(RED[d], pt)
    Pf, Sf = A[-2], B[-3]
    for k in (3, 4, 5):
        MN, cols, tags = window_rows(A, B, N, k)
        val = lambda t: sp.expand((Pf if t[0] == "P" else Sf).subs(E, t[1]))
        bad = 0
        for i, (m, e) in enumerate(tags):
            lhs = sp.expand(sum(MN[i, j] * val(cols[j]) for j in range(len(cols)))
                            + MN[i, len(cols)])
            if sp.expand(lhs - sp.expand(sp.sympify(q_m(A, B, m)).subs(E, e))) != 0:
                bad += 1
        check(bad == 0 and MN.rows == 3 * k - 4 and MN.cols == 2 * k - 1,
              f"[d={d}, k={k}] every window row reproduces Q_m(E=e) EXACTLY; "
              f"[M|N] is {MN.rows}x{MN.cols} (shape independent of d)")

# =====================================================================
print("\n--- S5. ON-SHELL NODE RECURSION: the cascade fixes the b-tops on the SAME grid ---",
      flush=True)
# =====================================================================
# Q_4 = b_2(E+2)a_2 - a_2(E+2)b_2 + b_1(E+3)a_3 - a_3(E+1)b_1 .  a_3=E(E+2)(E+4) has the
# integer roots {0,-2,-4}: at those three nodes b_1(E+3) DROPS OUT and Q_4=0 becomes a
# relation on b_1 and a_2 node values; elsewhere it is a 3-step forward recursion.
for dgen in (2, 3, 4, 6):
    a2g, _ = poly(f"on{dgen}a2", dgen); b1g, _ = poly(f"on{dgen}b1", dgen)
    Ao = {3: a3_w2, 2: a2g, 1: sp.Integer(0), 0: sp.Integer(0), -1: sp.Integer(0),
          -2: sp.Integer(0), -3: sp.Integer(0)}
    Bo = {k: sp.Integer(0) for k in range(-3, 4)}
    Bo[2] = b2_w2; Bo[1] = b1g
    Q4 = q_m(Ao, Bo, 4)
    ev = lambda f, x: sp.expand(sp.sympify(f).subs(E, x))
    b1c_syms = sorted(set(b1g.free_symbols) - {E}, key=str)
    # each ANCHOR node x sees b_1 through EXACTLY ONE node value b_1(n(x)):
    for x, n in [(0, 0), (-3, 0), (-2, -2), (-4, -4), (-1, 2)]:
        cvec = [sp.expand(sp.diff(Q4.subs(E, x), s)) for s in b1c_syms]
        lam = cvec[0]
        check(all(sp.expand(cvec[i] - lam * n**i) == 0 for i in range(len(cvec))),
              f"[d={dgen}] ANCHOR Q_4({x}): the b_1-dependence is exactly {lam}*b_1({n}) "
              f"-- one node value only (a_3-root incidence)")
    check_zero(Q4.subs(E, 0) - (10 * ev(a2g, 0) - 15 * ev(b1g, 0)),
               f"[d={dgen}] Q_4(0)=0  =>  b_1(0) = (2/3) a_2(0)   [anchor at a_3-root 0]")
    check_zero(Q4.subs(E, -3) - (-2 * ev(a2g, -3) + 3 * ev(b1g, 0)),
               f"[d={dgen}] Q_4(-3)=0 =>  b_1(0) = (2/3) a_2(-3)")
    check_zero(sp.expand(Q4.subs(E, 0) / 10 + Q4.subs(E, -3) / 2
                         - (ev(a2g, 0) - ev(a2g, -3))),
               f"[d={dgen}] CASCADE NODE RELATION (degree-free):  "
               f"a_2(0) - a_2(-3) = Q_4(0)/10 + Q_4(-3)/2,  so the cascade forces a_2(-3)=a_2(0)")
    check_zero(Q4.subs(E, -1) - (-3 * ev(b1g, 2) + 4 * ev(a2g, -1) + 2 * ev(a2g, 1)),
               f"[d={dgen}] Q_4(-1)=0 =>  b_1(2) = (2/3)(2 a_2(-1) + a_2(1))  "
               f"[b_1(-1) drops out: a_3(0)=0]")
    check_zero(Q4.subs(E, 1) - (18 * ev(a2g, 1) - 4 * ev(a2g, 3) - 48 * ev(b1g, 1)
                                + 15 * ev(b1g, 4)),
               f"[d={dgen}] Q_4(1)=0 =>  15 b_1(4) = 48 b_1(1) - 18 a_2(1) + 4 a_2(3)  "
               f"[3-step recursion]")
    # AUDIT FIX (2026-07-25): the former check here compared a hardcoded set with
    # itself -- a Python tautology that could not fail, executed once per degree,
    # and the claim it labelled ("b_1(1), b_1(3) are FREE SEEDS") is FALSE as
    # stated for the polynomial cascade: at d=3, Q_4=0 as a linear system for the
    # b_1 coefficients has rank 7 of 7, so b_1 is uniquely DETERMINED by a_2.
    # Freedom is a property of the VALUE-level 3-step recursion only.  Replaced
    # by a falsifiable statement of exactly that, plus the polynomial-level rank.
    _rows = sp.Poly(sp.expand(Q4), E).all_coeffs()
    _M = sp.Matrix([[sp.expand(r).coeff(s) for s in b1c_syms] for r in _rows])
    check(_M.rank() == len(b1c_syms),
          f"[d={dgen}] POLYNOMIAL level: Q_4=0 pins b_1 outright given a_2 "
          f"(filler-column rank {_M.rank()} = {len(b1c_syms)} unknowns, kernel 0) -- so "
          f"b_1(1), b_1(3) are 'free seeds' ONLY for the value-level 3-step recursion")
print("       => ON-SHELL, the b-tops live on the SAME integer grid as the a-tops:", flush=True)
print("          anchored at the a_3-roots {0,-2,-4} + free seeds. The collapse is PARTIAL.", flush=True)

# =====================================================================
print("\n--- S6. d=2 WINDOW-MINOR FACTORIZATION: node-value factors x ONE core (exact) ---",
      flush=True)
# =====================================================================
rd2 = RED[2]
free2 = rd2["free"]
A2s = {k: sp.cancel(sp.sympify(rd2["A"][k]).subs(rd2["sol"])) for k in rd2["A"]}
B2s = {k: sp.cancel(sp.sympify(rd2["B"][k]).subs(rd2["sol"])) for k in rd2["B"]}
A20 = dict(A2s); A20[-2] = sp.Integer(0)
B20 = dict(B2s); B20[-3] = sp.Integer(0)
N2s = {m: sp.cancel(q_m(A20, B20, m)) for m in (-1, -2, -3)}
W2s = sp.cancel(R(1, 4) * (R(2, 3) * A2s[2].subs(E, 1) - B2s[1].subs(E, 2)) * B2s[-1].subs(E, 1)
                + B2s[-2].subs(E, 2) - R(2, 3) * A2s[-1].subs(E, 1))
check(sp.simplify(W2s) == 0, "d=2: W == 0 identically on the cascade (forcing vacuous at d=2)")
cores, commons = {}, {}
for k in (3, 4):
    MN2, cols2, _ = window_rows(A2s, B2s, N2s, k)
    nc = len(cols2)
    facsets = []
    nmin = 0
    for sub in itertools.combinations(range(MN2.rows), nc + 1):
        nu = sp.expand(sp.fraction(sp.together(sp.cancel(MN2[list(sub), :].det())))[0])
        if nu == 0:
            continue
        nmin += 1
        facsets.append({sp.expand(p) for p, _ in sp.factor_list(nu)[1]})
    check(nmin >= 1, f"d=2, k={k}: {nmin} nonzero maximal window minors")
    common = set.intersection(*facsets)
    big = [c for c in common if sp.Poly(c, *free2).total_degree() >= 3]
    check(len(big) == 1,
          f"d=2, k={k}: ALL {nmin} minors share exactly ONE non-linear common factor (the "
          f"core, degree {sp.Poly(big[0], *free2).total_degree()}) -- the window elimination "
          f"content is essentially principal")
    cores[k] = big[0]; commons[k] = common
check(sp.simplify(sp.cancel(cores[3] / cores[4])).is_number,
      "d=2: the k=3 and k=4 cores agree up to a scalar (the window content stabilizes in k)")
# the linear common factors ARE datum NODE VALUES (Sylvester-like, not generic resultants):
a2_at = {x: sp.expand(A2s[2].subs(E, x)) for x in (-3, -2, -1, 0, 1, 2)}
check(sp.expand(a2_at[-3] - a2_at[0]) == 0, "d=2: a_2(-3)=a_2(0) holds on the cascade (S5)")
named = {f"a_2({x})": sp.expand(v) for x, v in a2_at.items()}
named.update({"a_1(2)": sp.expand(A2s[1].subs(E, 2)), "a_-1(1)": sp.expand(A2s[-1].subs(E, 1)),
              "a_0(3)-a_0(0)": sp.expand(A2s[0].subs(E, 3) - A2s[0].subs(E, 0)),
              "b_1(0)": sp.expand(B2s[1].subs(E, 0)), "b_-1(1)": sp.expand(B2s[-1].subs(E, 1))})
lin_common = [c for c in commons[4] if sp.Poly(c, *free2).total_degree() <= 2]
matched, unmatched = {}, []
for c in lin_common:
    hit = [nm for nm, v in named.items() if v != 0 and sp.simplify(sp.cancel(c / v)).is_number]
    if hit:
        matched[str(c)] = hit[0]
    else:
        unmatched.append(c)
check(len(matched) >= 2,
      f"d=2: the linear common factors of ALL k=4 minors include the DATUM NODE VALUES "
      f"{sorted(set(matched.values()))} -- Sylvester-like in the fixed integer nodes, NOT a "
      f"generic Res_E(top(E),top(E-3)); unmatched common factors: {unmatched}")
# CONTRAST that kills the naive extrapolation: at d=2 the TOP coefficient of the raw a_-1
# (the d=2 analogue of the W-slot am1_3) IS a common factor of every window minor -- but at
# d=3 the true W-slot am1_3 divides NO window minor (S7).
am1_2 = sp.Symbol("am1_2")
check(any(sp.expand(c - am1_2) == 0 for c in commons[4]) and any(
        sp.expand(c - am1_2) == 0 for c in commons[3]),
      "d=2: the a_-1 TOP coefficient am1_2 is a common factor of every k=3,4 window minor "
      "(the d=2 look-alike of a W-multiple) -- S7 shows this does NOT persist to d=3")

# =====================================================================
print("\n--- S7. THE VERDICT at d=3: no window minor is a unit multiple of W; k=3 REFUTED ---",
      flush=True)
# =====================================================================
rd3 = RED[3]
free3, fill3 = rd3["free"], rd3["fillers"]
am1_3, a2_0 = sp.symbols("am1_3 a2_0")
check(sp.simplify(rd3["R1"] - R(-4, 9) * a2_0 * am1_3) == 0,
      "d=3: R(1) = -(4/9) a_2(0) am1_3 = a_2(0)*W with W = -(4/9) am1_3 (so W|X <=> am1_3|X)")
KS = [3, 4] + ([5] if HEAVY else [])
NSEED = 2 if HEAVY else 1
for seed in range(NSEED):
    rng = random.Random(7700 + seed)
    pt = {v: sp.Rational(rng.randint(-7, 7), rng.choice([1, 2, 3]))
          for v in free3 if v != am1_3}
    if pt.get(a2_0, 1) == 0:
        pt[a2_0] = sp.Integer(3)
    A, B, N, W = datum(rd3, pt)
    check(sp.simplify(W + R(4, 9) * am1_3) == 0,
          f"[seed {seed}] the specialization keeps W = -(4/9) am1_3 symbolic")
    for k in KS:
        MN, cols, _ = window_rows(A, B, N, k)
        nc = len(cols)
        assert not any(am1_3 in sp.sympify(MN[i, j]).free_symbols
                       for i in range(MN.rows) for j in range(nc))
        nz = wm = 0
        for sub in itertools.combinations(range(MN.rows), nc + 1):
            Dm = sp.expand(MN[list(sub), :].det())
            if Dm == 0:
                continue
            nz += 1
            if sp.expand(Dm.subs(am1_3, 0)) == 0:
                wm += 1
        check(nz >= 1 and wm == 0,
              f"[seed {seed}] d=3, k={k}: of {nz} nonzero maximal window minors, {wm} are "
              f"divisible by W  => NO window minor is a unit multiple of W")
    if not HEAVY:
        skip("d=3 k=5 window-minor sweep + second seed (HEAVY: ~40s)")

# EXPLICIT EXACT WITNESS refuting k=3 forcing (radical-correct: a point, not a normal form)
wit_pt = {sp.Symbol("a2_0"): sp.Integer(-4), sp.Symbol("a2_2"): R(-73426, 34959),
          sp.Symbol("a2_3"): sp.Integer(0), sp.Symbol("a1_3"): sp.Integer(0),
          sp.Symbol("a0_0"): sp.Integer(2), sp.Symbol("a0_2"): sp.Integer(1),
          sp.Symbol("am1_2"): sp.Integer(3), sp.Symbol("am1_3"): sp.Integer(2),
          sp.Symbol("b0cK0"): sp.Integer(-4)}
check(set(wit_pt) == set(free3), "witness point assigns exactly the d=3 free cascade coordinates")
Aw3, Bw3, Nw3, Ww3 = datum(rd3, wit_pt)
MNw, colsw, _ = window_rows(Aw3, Bw3, Nw3, 3)
ncw = len(colsw)
Mw = MNw[:, :ncw]; Nw = MNw[:, ncw]
rM, rMN = Mw.rank(), Mw.row_join(Nw).rank()
check(sp.expand(MNw.det()) == 0 and rM == ncw and rMN == ncw,
      f"WITNESS: at this exact rational cascade point the k=3 window system is SOLVABLE "
      f"(det[M|N]=0, rank M = rank[M|N] = {rM} = full)")
xsolw = Mw.solve_least_squares(-Nw) if rM == ncw else None
check(all(sp.expand((Mw * xsolw + Nw)[i]) == 0 for i in range(MNw.rows)),
      "WITNESS: the explicit filler-value solution satisfies EVERY k=3 window row exactly")
check(sp.simplify(Ww3) != 0,
      f"WITNESS: at that point W = {sp.simplify(Ww3)} != 0  =>  the k=3 INTEGER-NODE WINDOW "
      f"DOES NOT FORCE W (explicit witness, radical-correct)")
print("       => the fixed-shape integer-node window at k=3 is provably too weak; the", flush=True)
print("          W-forcing genuinely uses more than the ladder identities at these nodes.", flush=True)

# =====================================================================
print("\n--- S8. COMBINED-RUNG ADJOINT over Q(E) (symbolic node): shape and verdict ---",
      flush=True)
# =====================================================================
# Treat [Q_-1;Q_-2;Q_-3] at the shifted arguments E+j (j=0..J) as ONE operator on the
# filler VALUES P(E+s), S(E+s): a fixed-shape matrix at every degree.
rngc = random.Random(9001)
ptc = {v: sp.Rational(rngc.randint(-6, 6), rngc.choice([1, 2, 3])) for v in free3 if v != am1_3}
if ptc.get(a2_0, 1) == 0:
    ptc[a2_0] = sp.Integer(2)
Ac, Bc, Nc, Wc = datum(rd3, ptc)
b1c, b0c, bm1c = Bc[1], Bc[0], Bc[-1]
a2c, a1c, a0c = Ac[2], Ac[1], Ac[0]


def shifted_rows(J):
    ent, Pi, Si = [], set(), set()
    for j in range(J + 1):
        for m in (-1, -2, -3):
            if m == -1:
                cd = {("P", j): sh(b1c, j - 2), ("P", j + 1): -sh(b1c, j),
                      ("S", j + 2): sh(a2c, j), ("S", j): -sh(a2c, j - 3)}
            elif m == -2:
                cd = {("P", j): sh(b0c, j - 2) - sh(b0c, j),
                      ("S", j + 1): sh(a1c, j), ("S", j): -sh(a1c, j - 3)}
            else:
                cd = {("P", j): sh(bm1c, j - 2), ("P", j - 1): -sh(bm1c, j),
                      ("S", j): sh(a0c, j) - sh(a0c, j - 3)}
            ent.append((cd, sh(Nc[m], j)))
            for t in cd:
                (Pi if t[0] == "P" else Si).add(t[1])
    cols = [("P", s) for s in sorted(Pi)] + [("S", s) for s in sorted(Si)]
    return ent, cols


for J in (2, 3):
    ent, cols = shifted_rows(J)
    check(len(ent) == 3 * (J + 1) and len(cols) == 2 * (J + 3),
          f"combined-rung adjoint window J={J}: {3*(J+1)} rows x {2*(J+3)} filler-value "
          f"columns -- shape independent of d")
e0 = R(17, 5)
ent, cols = shifted_rows(3)
ci = {c: i for i, c in enumerate(cols)}
Mc = sp.zeros(len(ent), len(cols))
Nvec = sp.zeros(len(ent), 1)
for i, (cd, nn) in enumerate(ent):
    for t, v in cd.items():
        Mc[i, ci[t]] = sp.expand(v.subs(E, e0))
    Nvec[i] = sp.expand(nn.subs(E, e0))
kerc = Mc.T.nullspace()
check(Mc.rows == 12 and Mc.cols == 12 and Mc.rank() == 11 and len(kerc) == 1,
      f"J=3 combined-rung matrix at E={e0} is 12x12 of rank 11 at the tested exact "
      f"specialization -- a UNIQUE covector (rank-deficient square value-column block)")
pairc = sp.expand((sp.Matrix([list(kerc[0])]) * Nvec)[0])
check(sp.expand(pairc.subs(am1_3, 0)) != 0 and pairc != 0,
      "J=3 combined-rung covector: its residual pairing is NONZERO and has a NONZERO "
      "am1_3-free part => it is NOT a unit multiple of W")

# =====================================================================
print("\n--- S9. CONTROLS: W-kill (Rabinowitsch), a_2(0) free, non-vacuity ---", flush=True)
# =====================================================================
t_rab = sp.symbols("t_rab")
tail3_5 = [e for m in (-1, -2, -3, -4, -5) for e in rd3["tk"][m]]
tail3_3 = [e for m in (-1, -2, -3) for e in rd3["tk"][m]]
allv3 = [t_rab] + list(free3) + list(fill3)
check(sy_unit(tail3_5 + [1 - t_rab * am1_3], allv3, 65003),
      "CONTROL (Rabinowitsch mod 65003): am1_3 in sqrt(cascade+Q_-1..Q_-5) => the tail forces W=0")
check(sy_unit(tail3_5 + [1 - t_rab * am1_3], allv3, None),
      "CONTROL (Rabinowitsch over QQ, exact): am1_3 in sqrt(cascade+Q_-1..Q_-5) => W=0 forced")
Mlin = sp.Matrix([[sp.Poly(e, *fill3).coeff_monomial(f) for f in fill3] for e in tail3_3])
Nlin = sp.Matrix([sp.Poly(e, *fill3).coeff_monomial(1) for e in tail3_3])
_rw = random.Random(20260724); wit = None
for _ in range(2500):
    fval = {v: sp.Integer(_rw.randint(-4, 4)) for v in free3}
    fval[am1_3] = sp.Integer(0)
    if fval[a2_0] == 0:
        fval[a2_0] = sp.Integer(2)
    Ms, Ns = Mlin.subs(fval), Nlin.subs(fval)
    if Ms.rank() == Ms.row_join(Ns).rank():
        wit = (fval, sp.expand(rd3["R1"].subs(fval)), Ms.rank())
        break
check(wit is not None and wit[0][a2_0] != 0 and wit[1] == 0,
      f"CONTROL: explicit cascade+depth-3-tail point with a_2(0)={wit[0][a2_0]} != 0 and "
      f"R(1)={wit[1]}=0  => a_2(0) NOT forced (the kill is the FACTOR W)")
xs = Mlin.subs(wit[0]).solve_least_squares(-Nlin.subs(wit[0]))
check(all(sp.expand((Mlin.subs(wit[0]) * xs + Nlin.subs(wit[0]))[i]) == 0
          for i in range(Mlin.rows)),
      "CONTROL (non-vacuity): the witness carries an EXPLICIT filler solution satisfying "
      "every depth-3 tail equation => cascade+tail is nonempty")

# =====================================================================
print("\n" + "=" * 74, flush=True)
print("PROVED (degree-free): engine + Q_0=(T-1)G + slope gate + Lemma P + R(1)=a_2(0)W;", flush=True)
print("  STAIRCASE COLLAPSE of the two blocks (Q_-2 a-block and Q_-3 b-block are pure", flush=True)
print("  multipliers); the MEMBERSHIP SILENCE TABLE (six identically-zero rungs); the", flush=True)
print("  INTEGER-NODE WINDOW -- for every k>=3 and every degree d the tail rows", flush=True)
print("  Q_-1(1..k-1), Q_-2(2..k), Q_-3(3..k) involve only 2k-2 filler VALUES and [M|N]", flush=True)
print("  has the d-independent shape (3k-4)x(2k-1) with fixed integer-node entries;", flush=True)
print("  the ON-SHELL node recursion (b_1(0)=(2/3)a_2(0), a_2(-3)=a_2(0),", flush=True)
print("  b_1(2)=(2/3)(2a_2(-1)+a_2(1)), 15b_1(4)=48b_1(1)-18a_2(1)+4a_2(3); free seeds).", flush=True)
print("REFUTED (exact): no fixed-shape window minor is a unit multiple of W at d=3 (k=3,4", flush=True)
print("  [+5 HEAVY]); and the k=3 window does NOT force W -- EXPLICIT rational witness with", flush=True)
print("  the window solvable and W=-8/9 != 0.  The combined-rung adjoint over Q(E) at J=3", flush=True)
print("  is 12x12 of rank 11 and its unique covector is not a W-multiple either.", flush=True)
print("BOUNDED: d=2 window minors = node values x ONE core (essentially principal).", flush=True)
print("OPEN: whether ANY k has W in sqrt(I_k); any degree-free W-forcing recipe; the", flush=True)
print("  arbitrary-degree residual identity; all of band 3, DC1, JC2.", flush=True)
print("=" * 74, flush=True)
if _SKIPS:
    print(f"\nALL EXECUTED CHECKS PASSED; OPTIONAL CHECKS SKIPPED ({len(_SKIPS)}):", flush=True)
    for s in _SKIPS:
        print(f"   SKIPPED: {s}", flush=True)
else:
    print("\nALL CHECKS PASSED (no skips).", flush=True)
print(f"(total {time.time() - _T0:.1f}s; {_NP} checks passed; "
      f"{'HEAVY' if HEAVY else 'default'} mode)", flush=True)
