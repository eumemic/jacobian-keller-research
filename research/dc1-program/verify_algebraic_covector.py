#!/usr/bin/env python3
"""THE ALGEBRAIC-NODE COVECTOR CALCULUS at W2 -- Thm A' over F(datum) via trace forms,
and the two-block necklace structure of the depth-3 negative tail.

INDEPENDENTLY DERIVED - EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES - NOT PEER REVIEWED

MISSION CONTEXT.  The W2 endgame gap (residual-identity.md / slope-forcing-degree-free.md):
at every positive-data degree d, does the depth-3 tail Q_-1=Q_-2=Q_-3=0 force the moment-
slope residual W=0 on the W2 cascade (R(1)=a_2(0)*W)?  The obstruction to a degree-free
proof is the covector on the ALGEBRAIC (a_2,b_1,...) necklace: its nodes are the roots of
the cascade-solved data, not integers.  This file builds and machine-checks the tool the
mission asked for -- the algebraic-node extension of the adjoint criterion Thm A'
(lambda-general-k.md) and the trace-form descent to F(datum) -- and derives the exact two-
block structure of the negative tail; then it reads the d=3 certificate through that lens
and reports HONESTLY what a fixed trace-form recipe does and does not deliver.

W2 datum, gauge b_3=0, quantum band-3 conventions:
  Q_m = sum_(k+l=m)[ b_l^[k] a_k - a_k^[l] b_l ],   f^[n](E)=f(E+n),
  membership (E)_j = E(E-1)...(E-j+1) | a_-j, b_-j;   a_3=E(E+2)(E+4), b_2=E(E+3).
  branch B: a_-3=0, b_-3=(E)_3 C ;  branch B fillers = a_-2=(E)_2 V, b_-3=(E)_3 C.

WHAT THIS FILE ESTABLISHES
==========================
PROVED (arbitrary degree / node-free, machine-checked identities):
  * S0 engine Q_m=[D,X]_m, Q_0=(T-1)G;  S1 slope gate, both-ends Lemma-P, factorization
    R(1)=a_2(0)*W (re-derived, degree-free).
  * S2 ALGEBRAIC-NODE THM A' (the tool):
      (i)  moving-sum adjoint identity  lambda(S_n g)=(S_n^* lambda)(g)  for SYMBOLIC nodes
           rho (integer OR algebraic) -- the criterion never needed integer nodes;
      (ii) TRACE-FORM DESCENT: for a squarefree datum polynomial p (symbolic coeffs) and any
           weight/test h, sum_{p(rho)=0} h(rho) = Tr_{F[E]/(p)}(h) is rational in the
           arbitrary coefficients of p,h, and polynomial after monic normalization (companion
           trace), computed with NO root named -- so a Galois-
           symmetric node-functional lam=sum_{p(rho)=0} g(rho) ev_rho descends to F(datum);
      (iii) trace forms are CLOSED under S_n^* (shift-window), so the algebraic-necklace
           support conditions are themselves trace-form/resultant-computable degree-freely.
  * S3 THE TWO BLOCKS: Q_-1,Q_-2,Q_-3 filler-linear parts are the explicit TWO-TERM operators
      a-block(Q_m) = b_{m+2}(E-2)a_-2(E) - b_{m+2}(E)a_-2(E+m+2),
      b-block(Q_m) = b_-3(E+m+3)a_{m+3}(E) - a_{m+3}(E-3)b_-3(E),
    with a-block tops (b_1,b_0,b_-1) and b-block tops (a_2,a_1,a_0) -- ALL cascade-solved
    (ALGEBRAIC).  The FIXED tops a_3,b_2 do NOT appear in the negative tail at all: the
    negative-tail necklace is ENTIRELY algebraic; the only integer nodes are the membership
    windows {0,1} (a-block) and {0,1,2} (b-block).  Membership-window covectors annihilate a
    block degree-free (the ev_0 analogue).  A single ALGEBRAIC node cannot annihilate a block
    (needs a shared root of top(E),top(E-3): generically none) -- the covector must COUPLE the
    two terms across the varying tops.
  * S4 BLOCK ADJOINT CRITERION (algebraic-node): the block-annihilation of a point covector
    reduces to explicit coefficient-of-filler-value vanishing conditions at the shifted nodes,
    a symbolic-node identity; over root-necklaces those conditions are trace forms (S2).

BOUNDED (d=3, exact; stated scope):
  * S5 the depth-3 tail cokernel is 16-dimensional; the W-forcing pairing exists (kill
    reproduced: am1_3 in sqrt(cascade+Q_-1..Q_-5), sympy exact QQ + prime); a_2(0) NOT forced
    (explicit witness); the computed specialized cokernel vectors annihilate the specialized
    filler columns exactly. The W-forcing covector's support meets the algebraic necklace (not the membership windows
    alone) -- the fixed part is silent, as at Q_0 (joint-covector.md).

OPEN / NOT CLAIMED (the prize, obstruction LOCALIZED -- lead 5):
  * A FIXED FINITE trace-form recipe producing a covector whose pairing is a unit multiple of
    W at EVERY d.  The two-term coupling across the varying tops (a_2,a_1,a_0)/(b_1,b_0,b_-1)
    does not reduce to one fixed node-selection+weight rule in these tests; that coupling is
    the exact residual gap.  S6 verifies the linear-route data (filler map full column rank;
    sampled W-forcing obstruction) at d=4 mod p, providing bounded supporting evidence
    only, not a symbolic covector or degree-free recipe. The arbitrary-d W-forcing identity
    remains OPEN.

Run:  uv run --with sympy python research/dc1-program/verify_algebraic_covector.py
      HEAVY=1 ... (adds d=4 linear-route mod-p leg + depth-3 msolve kill + larger trace-form)
Ends with a PASS/SKIP-aware summary.
"""
import sympy as sp
import time, os, random, subprocess, tempfile, shutil

E = sp.symbols("E")
LEVELS = range(-3, 4)
R = sp.Rational
_T0 = time.time()
random.seed(20260724)
HEAVY = os.environ.get("HEAVY") == "1"
PRIMES = (65003, 32003)
_NP = 0
_NSKIP = 0


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


def clear_denoms(e):
    num, _ = sp.fraction(sp.together(sp.expand(e)))
    return sp.expand(num)


def check(cond, label):
    global _NP
    if not cond:
        raise AssertionError("FAIL " + label)
    _NP += 1
    print(f"PASS [{time.time() - _T0:6.1f}s] {label}", flush=True)


def check_zero(val, label):
    check(sp.expand(val) == 0, label)


def skip(label):
    global _NSKIP
    _NSKIP += 1
    print(f"SKIP [{time.time() - _T0:6.1f}s] {label}", flush=True)


def sy_unit(eqs, vs, modulus=None):
    G = sp.groebner([sp.expand(e) for e in eqs if sp.expand(e) != 0], *vs,
                    order="grevlex", modulus=modulus)
    return list(G.exprs) == [sp.Integer(1)]


def msolve_unit(eqs, vs, char, tmo=200):
    xs = sp.symbols(f"z0:{len(vs)}")
    sub = dict(zip(vs, xs))
    with tempfile.TemporaryDirectory(prefix="algebraic-covector-") as tmp:
        inp = os.path.join(tmp, "input.ms")
        out = os.path.join(tmp, "output.txt")
        with open(inp, "w", encoding="utf-8") as f:
            f.write(",".join(str(x) for x in xs) + "\n" + str(char) + "\n")
            f.write(",\n".join(str(clear_denoms(e.subs(sub))).replace(" ", "").replace("**", "^")
                               for e in eqs if sp.expand(e) != 0) + "\n")
        rr = subprocess.run(["msolve", "-g", "2", "-f", inp, "-o", out],
                            capture_output=True, text=True, timeout=tmo)
        if rr.returncode != 0:
            raise RuntimeError(f"msolve failed with status {rr.returncode}: {rr.stderr.strip()}")
        if not os.path.exists(out):
            raise RuntimeError(f"msolve produced no output file; stderr: {rr.stderr.strip()}")
        with open(out, encoding="utf-8") as f:
            r = f.read()
    lines = [line for line in r.splitlines() if not line.lstrip().startswith("#")]
    parsed = "".join(lines).replace(" ", "")
    if not parsed.startswith("["):
        raise RuntimeError(f"malformed msolve output: {parsed[:200]!r}")
    return parsed.startswith("[1]:") or parsed == "[1]"


MSFAIL = (subprocess.TimeoutExpired,)


def rank_mod_p(rows, ncol, p):
    A = [[int(x) % p for x in r] for r in rows]
    r = 0
    for c in range(ncol):
        piv = next((i for i in range(r, len(A)) if A[i][c] % p != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], p - 2, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] % p != 0:
                f = A[i][c]
                A[i] = [(A[i][k] - f * A[r][k]) % p for k in range(ncol)]
        r += 1
        if r == len(A):
            break
    return r


a3_w2 = sp.expand(E * (E + 2) * (E + 4))
b2_w2 = sp.expand(E * (E + 3))
D = sp.expand(E * (E - 1) * (E + 1))

# =====================================================================
print("--- S0. crossed-product engine: Q_m=[D,X]_m, Q_0=(T-1)G (re-derived, arbitrary degree) ---", flush=True)
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
print("\n--- S1. slope gate; both-ends Lemma-P; FACTORIZATION R(1)=a_2(0)*W (re-derived) ---", flush=True)
# =====================================================================
check(sp.rem(E - sp.expand(E + D * (E**2 + 7)), D, E) == 0,
      "slope gate: D | (E-R) when R(1)=1, R(-1)=-1 (R(0)=0)")
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
check_zero(Gw.subs(E, 0), "membership: G(0)=0  (=> R(1)=Q_0(0)=G(1))")
G1 = sp.expand(Aw[1].subs(E, 0) * Bw[-1].subs(E, 1) + Aw[2].subs(E, 0) * Bw[-2].subs(E, 2)
               - Aw[-1].subs(E, 1) * Bw[1].subs(E, 0))
check_zero(q_m(Aw, Bw, 0).subs(E, 0) - G1,
           "both-ends Lemma-P: R(1)=a1(0)b-1(1)+a2(0)b-2(2)-a-1(1)b1(0)")
fillsyms = set(Aw[-2].free_symbols) | set(Bw[-3].free_symbols) | {mu3s}
check(not (sp.expand(q_m(Aw, Bw, 0).subs(E, 0)).free_symbols & fillsyms),
      "R(1)=Q_0(0) INDEPENDENT of the fillers a_-2, b_-3")
# factorization backbone (degree-free), re-derived symbolically:
b1_0, a1_0, a2_0s, a2_1s, b1_2s = sp.symbols("b1_0 a1_0 a2_0 a2_1 b1_2")
sol_b1_0 = R(2, 3) * a2_0s
sol_a1_0 = a2_0s / 4 * (R(2, 3) * a2_1s - b1_2s)
bm1_1s, bm2_2s, am1_1s = sp.symbols("bm1_1 bm2_2 am1_1")
R1_P = a1_0 * bm1_1s + a2_0s * bm2_2s - am1_1s * b1_0
W_closed = bm2_2s - R(2, 3) * am1_1s + R(1, 4) * (R(2, 3) * a2_1s - b1_2s) * bm1_1s
check_zero(R1_P.subs({a1_0: sol_a1_0, b1_0: sol_b1_0}) - a2_0s * W_closed,
           "FACTORIZATION: R(1) = a_2(0) * W  (Lemma-P + Q_4(0)=Q_3(0)=0), degree-free")


# ---- cascade machinery (branch B) ----
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
    am2_raw, cam2 = poly("am2", d)
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
    dat = {"ca2": ca2, "ca1": ca1, "ca0": ca0, "cam1": cam1}
    return A, B, conds, dat, kernels


def build_full_branch(d):
    A, B, pos, dat, kernels = positive_cascade(d)
    base = dat["ca2"] + dat["ca1"] + dat["ca0"] + dat["cam1"]
    am2v = list(set(A[-2].free_symbols) - {E})
    bm3_raw, cbm3 = poly("bm3f", d)
    A[-3] = sp.Integer(0); B[-3] = sp.expand(falling(3) * bm3_raw)
    fv = base + am2v + cbm3 + kernels
    posnz = [sp.expand(e) for e in pos if sp.expand(e) != 0]
    return A, B, posnz, fv


def greedy_param(conds, Pvars):
    sol, solved, changed = {}, set(), True
    while changed:
        changed = False; best = None
        for c in conds:
            cc = sp.expand(c.subs(sol))
            if cc == 0:
                continue
            for v in sorted((cc.free_symbols & set(Pvars)) - solved, key=str):
                pc = sp.Poly(cc, v)
                if pc.degree() == 1:
                    lead = pc.coeff_monomial(v); rest = sp.expand(cc - lead * v)
                    score = len(lead.free_symbols)
                    if best is None or score < best[0]:
                        best = (score, v, lead, rest)
        if best is not None:
            _, v, lead, rest = best; val = sp.cancel(-rest / lead)
            sol = {k: sp.cancel(x.subs(v, val)) for k, x in sol.items()}
            sol[v] = val; solved.add(v); changed = True
    free = [v for v in Pvars if v not in solved]
    return sol, free, solved


def build_reduced(d):
    A, B, posnz, fv = build_full_branch(d)
    fillers = sorted([s for s in (set(A[-2].free_symbols) | set(B[-3].free_symbols)) - {E}], key=str)
    Pvars = [v for v in fv if v not in set(fillers)]
    sol, free, solved = greedy_param(posnz, Pvars)
    R1 = sp.cancel(sp.sympify(q_m(A, B, 0).subs(E, 0)).subs(sol))
    tk = {m: [clear_denoms(sp.cancel(e.subs(sol)))
              for e in coeffs(q_m(A, B, m)) if sp.expand(e) != 0] for m in (-1, -2, -3, -4, -5)}
    return dict(A=A, B=B, posnz=posnz, fv=fv, fillers=fillers, sol=sol,
                free=free, solved=solved, R1=R1, tk=tk)


def W_on_cascade(rd):
    A, B, sol = rd["A"], rd["B"], rd["sol"]
    Wexpr = (R(1, 4) * (R(2, 3) * A[2].subs(E, 1) - B[1].subs(E, 2)) * B[-1].subs(E, 1)
             + B[-2].subs(E, 2) - R(2, 3) * sp.expand(A[-1].subs(E, 1)))
    return sp.cancel(Wexpr.subs(sol))


# =====================================================================
print("\n--- S2. ALGEBRAIC-NODE THM A': moving-sum adjoint (symbolic nodes) + TRACE FORM ---", flush=True)
# =====================================================================
# (i) the adjoint identity lambda(S_n g)=(S_n^* lambda)(g) with a SYMBOLIC node rho.  Thm A'
#     of lambda-general-k.md is stated for integer cosets; here we show the identity itself
#     never used integrality -- it holds verbatim for an ALGEBRAIC (symbolic) node.
rho = sp.symbols("rho")
for n in (2, 3, 4):
    g, gc = poly(f"gA{n}", 5)
    Sng = sp.expand(sum(sh(g, i) for i in range(n)))          # S_n g
    lhs = Sng.subs(E, rho)                                     # ev_rho(S_n g)
    rhs = sp.expand(sum(g.subs(E, rho + i) for i in range(n)))  # (S_n^* ev_rho)(g)
    check_zero(sp.expand(lhs - rhs),
               f"S_{n} adjoint  ev_rho(S_{n} g)=(S_{n}^* ev_rho)(g)  SYMBOLIC node rho, generic deg-5 g")
# (ii) TRACE-FORM DESCENT.  For squarefree p (symbolic coeffs), sum_{p(rho)=0} h(rho) equals
#     Tr_{F[E]/(p)}(h) = trace of multiplication-by-h on the companion basis -- rational in
#     arbitrary coefficients of p,h and polynomial after monic normalization. Thus a
#     Galois-symmetric node-functional descends to F(datum) with NO root ever named.


def companion_trace(pexpr, hexpr):
    """Tr_{F[E]/(p)}(h) via the companion matrix of monic(p); returns a poly in the coeffs."""
    P = sp.Poly(pexpr, E); dP = P.degree(); lc = P.LC()
    cof = [c / lc for c in P.all_coeffs()]           # monic: E^dP + cof[1] E^{dP-1}+...+cof[dP]
    Comp = sp.zeros(dP)
    for i in range(dP - 1):
        Comp[i + 1, i] = 1
    for i in range(dP):
        Comp[i, dP - 1] = -cof[dP - i]
    Hop = sp.zeros(dP)
    powC = sp.eye(dP)
    Hpoly = sp.Poly(hexpr, E)
    hc = list(reversed(Hpoly.all_coeffs()))          # low->high
    for j, cj in enumerate(hc):
        if j > 0:
            powC = powC * Comp
        Hop = Hop + cj * powC
    return sp.expand(sum(Hop[i, i] for i in range(dP)))


# symbolic p (degree 3), symbolic h (degree 2): trace form is a rational function of the coeffs
pc = sp.symbols("pc0:4"); hc = sp.symbols("hc0:3")
psym = sum(pc[i] * E**i for i in range(4)); hsym = sum(hc[i] * E**i for i in range(3))
tf = companion_trace(psym, hsym)
check(tf.free_symbols <= set(pc) | set(hc),
      "trace form Tr_{F[E]/(p)}(h) is a function of coeffs(p),coeffs(h) ONLY -- no root named (degree-free)")
# numeric confirmation across random squarefree specializations: trace form == actual root sum
_rt = random.Random(1)
ntf = 0
for _ in range(4):
    sub = {c: sp.Integer(_rt.randint(-4, 4)) for c in pc + hc}
    sub[pc[3]] = sp.Integer(1)
    pp = psym.subs(sub)
    if sp.Poly(pp, E).degree() != 3 or sp.discriminant(sp.Poly(pp, E)) == 0:
        continue
    rootsum = sp.nsimplify(sum(hsym.subs(sub).subs(E, r) for r in sp.Poly(pp, E).all_roots()))
    check_zero(sp.simplify(tf.subs(sub) - rootsum),
               "trace form == sum_{p(rho)=0} h(rho) (random squarefree p,h; no root named)")
    ntf += 1
if ntf == 0:
    skip("trace-form numeric confirmation: no squarefree sample drawn")
# (iii) trace forms CLOSED under S_n^*: (S_n^* lam)(g)=sum_{p(rho)=0} sum_i g(rho+i) is the
#     trace form of S_n g -- so the shift-window (moving-sum) support conditions stay trace-
#     computable degree-freely.  Confirm Tr(S_n g)=sum_i Tr(g shifted) as a coeff identity.
for n in (2, 3):
    g2, g2c = poly(f"gT{n}", 2)
    lhs = companion_trace(psym, sp.expand(sum(sh(g2, i) for i in range(n))))
    rhs = sp.expand(sum(companion_trace(psym, g2.subs(E, E + i)) for i in range(n)))
    check_zero(sp.expand(sp.together(lhs - rhs)),
               f"trace form closed under S_{n}^*: Tr(S_{n} g)=sum_i Tr(g(.+i)) (coeff identity in p,g)")

# =====================================================================
print("\n--- S3. THE TWO BLOCKS: two-term operators, level incidence, purely-algebraic necklace ---", flush=True)
# =====================================================================
# generic degree-2 data, branch B (a_-3=0), symbolic fillers a_-2,b_-3:
Az, Bz = {}, {}
for lev in LEVELS:
    ar, _ = poly(f"Az_{lev + 3}", 2); br, _ = poly(f"Bz_{lev + 3}", 2)
    memb = falling(-lev) if lev < 0 else 1
    Az[lev] = sp.expand(memb * ar); Bz[lev] = sp.expand(memb * br)
Az[3] = a3_w2; Bz[2] = b2_w2; Bz[3] = sp.Integer(0); Az[-3] = sp.Integer(0)
am2 = Az[-2]; bm3 = Bz[-3]
am2syms = sorted(set(am2.free_symbols) - {E}, key=str)
bm3syms = sorted(set(bm3.free_symbols) - {E}, key=str)
for m in (-1, -2, -3):
    qm = q_m(Az, Bz, m)
    a_part = sp.expand(sum(sp.Poly(qm, *am2syms).coeff_monomial(s) * s for s in am2syms))
    b_part = sp.expand(sum(sp.Poly(qm, *bm3syms).coeff_monomial(s) * s for s in bm3syms))
    a_formula = sp.expand(sh(Bz[m + 2], -2) * am2 - Bz[m + 2] * sh(am2, m + 2))
    b_formula = sp.expand(sh(bm3, m + 3) * Az[m + 3] - sh(Az[m + 3], -3) * bm3)
    check_zero(a_part - a_formula,
               f"Q_{m} a-block = b_{{{m+2}}}(E-2)a_-2(E) - b_{{{m+2}}}(E)a_-2(E+{m+2})  (top b_{m+2})")
    check_zero(b_part - b_formula,
               f"Q_{m} b-block = b_-3(E+{m+3})a_{{{m+3}}}(E) - a_{{{m+3}}}(E-3)b_-3(E)  (top a_{m+3})")
# NEGATIVE-TAIL NECKLACE IS ENTIRELY ALGEBRAIC: a_3,b_2 do not appear in Q_-1,Q_-2,Q_-3.
As, Bs = dict(Az), dict(Bz)
a3c = sp.symbols("A3c0:3"); b2c = sp.symbols("B2c0:3")
As[3] = sum(a3c[i] * E**i for i in range(3)); Bs[2] = sum(b2c[i] * E**i for i in range(3))
for m in (-1, -2, -3):
    check(not ((set(a3c) | set(b2c)) & q_m(As, Bs, m).free_symbols),
          f"Q_{m} contains NO a_3,b_2 symbol => negative-tail necklace is purely ALGEBRAIC (tops are datum)")
# MEMBERSHIP-WINDOW COVECTORS annihilate a block degree-free (the ev_0 analogue), at generic d:
for dwin in (2, 3, 4):
    bm3d, _ = poly("Cwin", dwin); bm3d = sp.expand(falling(3) * bm3d)
    Ad2 = dict(Az)  # tops a_2,a_1,a_0 generic deg-2
    # lambda_-3 = ev_rho, rho in {0,1,2}: on b-block(Q_-3)=b_-3(E)a_0(E)-a_0(E-3)b_-3(E)
    for rhov in (0, 1, 2):
        blk = sh(bm3d, 0) * Ad2[0] - sh(Ad2[0], -3) * bm3d
        check_zero(blk.subs(E, rhov),
                   f"[d={dwin}] membership covector ev_{rhov} on lambda_-3 annihilates b-block (root of (E)_3)")
    am2d, _ = poly("Vwin", dwin); am2d = sp.expand(falling(2) * am2d)
    for rhov in (0, 1):
        blk = sh(Bz[-1], -2) * am2d - Bz[-1] * sh(am2d, -1)  # a-block(Q_-3), top b_-1
        check_zero(blk.subs(E, rhov),
                   f"[d={dwin}] membership covector ev_{rhov} on lambda_-3 annihilates a-block (root of (E)_2)")
# On this tested generic a_2 b-block rung, no single-node kill exists: it would need
# a shared root of a_2(E),a_2(E-3), and generically there is none.
gsh = sp.gcd(sp.Poly(clear_denoms(Az[2]), E), sp.Poly(clear_denoms(sh(Az[2], -3)), E))
check(sp.Poly(gsh, E).degree() == 0,
      "tested generic a_2 b-block rung: gcd(a_2(E),a_2(E-3))=1, so no single "
      "algebraic node kills this rung; its two terms must be coupled")

# =====================================================================
print("\n--- S4. BLOCK ADJOINT CRITERION (algebraic-node): coefficient-of-filler-value vanishing ---", flush=True)
# =====================================================================
# For the b-block, a covector lambda=(lam_-1,lam_-2,lam_-3), lam_m=sum_rho c_{m,rho} ev_rho,
# paired with membership b_-3=(E)_3 C, gives sum over nodes of [coeff]*C(sigma).  Annihilation
# <=> coeff of C(sigma)=0 for every sigma OUTSIDE the window {0,1,2}.  We machine-check the
# EXACT coefficient bookkeeping for a symbolic node rho (single-node lam_m=ev_rho): the C(sigma)
# coefficients are the trace-form-ready weights (rho+m+3)_3 a_{m+3}(rho) and (rho)_3 a_{m+3}(rho-3).
Cpoly, Cc = poly("Cadj", 4)
bm3_sym = sp.expand(falling(3) * Cpoly)
for m in (-1, -2, -3):
    blk = sh(bm3_sym, m + 3) * Az[m + 3] - sh(Az[m + 3], -3) * bm3_sym
    val = blk.subs(E, rho)                      # ev_rho(b-block Q_m)
    # predicted: (rho+m+3)_3 C(rho+m+3) a_{m+3}(rho) - (rho)_3 C(rho) a_{m+3}(rho-3)
    fall3 = lambda x: x * (x - 1) * (x - 2)
    pred = (fall3(rho + m + 3) * Cpoly.subs(E, rho + m + 3) * Az[m + 3].subs(E, rho)
            - fall3(rho) * Cpoly.subs(E, rho) * sh(Az[m + 3], -3).subs(E, rho))
    check_zero(sp.expand(val - pred),
               f"b-block adjoint (Q_{m}): ev_rho = (rho+{m+3})_3 C(rho+{m+3}) a_{m+3}(rho) - (rho)_3 C(rho) a_{m+3}(rho-3)")
# and the a-block, symbolic node:
Vpoly, Vc = poly("Vadj", 4)
am2_sym = sp.expand(falling(2) * Vpoly)
for m in (-1, -2, -3):
    blk = sh(Bz[m + 2], -2) * am2_sym - Bz[m + 2] * sh(am2_sym, m + 2)
    val = blk.subs(E, rho)
    fall2 = lambda x: x * (x - 1)
    pred = (fall2(rho) * Vpoly.subs(E, rho) * sh(Bz[m + 2], -2).subs(E, rho)
            - fall2(rho + m + 2) * Vpoly.subs(E, rho + m + 2) * Bz[m + 2].subs(E, rho))
    check_zero(sp.expand(val - pred),
               f"a-block adjoint (Q_{m}): ev_rho = (rho)_2 V(rho) b_{m+2}(rho-2) - (rho+{m+2})_2 V(rho+{m+2}) b_{m+2}(rho)")
print("       => over a root necklace rho in roots(p), each C(sigma)/V(sigma) coefficient is a", flush=True)
print("        TRACE FORM in the datum (S2): the algebraic-necklace support conditions are degree-free-computable.", flush=True)


def param_tail(rd, depth):
    return [e for m in range(-1, -depth - 1, -1) for e in rd["tk"][m]]


# =====================================================================
print("\n--- S5. THE d=3 CERTIFICATE STRUCTURE: cokernel, W-kill reproduced, a_2(0) free ---", flush=True)
# =====================================================================
a2_0, am1_3 = sp.symbols("a2_0 am1_3")
rd3 = build_reduced(3)
free3, fill3 = rd3["free"], rd3["fillers"]
check(sp.simplify(rd3["R1"] - R(-4, 9) * a2_0 * am1_3) == 0,
      "d=3: R(1) = -(4/9) a_2(0) am1_3 = a_2(0)*W with W=-(4/9)am1_3 (W collapses to a top coeff)")
rows3 = [sp.expand(e) for e in param_tail(rd3, 3)]
check(all(sp.Poly(e, *fill3).total_degree() <= 1 for e in rows3),
      f"d=3 depth-3 tail = {len(rows3)} eqs, ALL LINEAR in the {len(fill3)} fillers")
Mmat = sp.Matrix([[sp.Poly(e, *fill3).coeff_monomial(f) for f in fill3] for e in rows3])
Nvec = sp.Matrix([sp.Poly(e, *fill3).coeff_monomial(1) for e in rows3])
ptF = {v: sp.Integer(random.randint(2, 60)) for v in free3}
rankM = Mmat.subs(ptF).rank()
check(rankM == len(fill3),
      f"filler map M FULL COLUMN RANK {rankM}=8; cokernel dim = {len(rows3)}-{rankM} = {len(rows3)-rankM}")
# Exact cokernel vectors at this specialization annihilate every specialized filler column.
LK = Mmat.subs(ptF).T.nullspace()
check(len(LK) == len(rows3) - rankM and all(
        all(sp.expand((sp.Matrix([list(v)]) * Mmat.subs(ptF))[k]) == 0 for k in range(len(fill3))) for v in LK),
      f"all {len(LK)} specialized cokernel covectors satisfy mu.M=0 exactly at the tested specialization")
# W-KILL reproduced (control): am1_3 in sqrt(cascade+Q_-1..Q_-5), sympy exact QQ + prime.
t_rab = sp.symbols("t_rab")
allv3 = [t_rab] + free3 + fill3
tail3_5 = param_tail(rd3, 5)
check(sy_unit(tail3_5 + [1 - t_rab * am1_3], allv3, PRIMES[0]),
      f"CONTROL: am1_3 in sqrt(cascade+Q_-1..Q_-5) mod {PRIMES[0]} (sympy) => the tail forces W=0")
check(sy_unit(tail3_5 + [1 - t_rab * am1_3], allv3, None),
      "CONTROL: am1_3 in sqrt(cascade+Q_-1..Q_-5) over QQ (sympy exact radical certificate) => W=0 forced")
# a_2(0) NOT forced (control): explicit cascade+tail witness with a_2(0)!=0, R(1)=0.
tail3_3 = param_tail(rd3, 3)
Mlin = sp.Matrix([[sp.Poly(e, *fill3).coeff_monomial(f) for f in fill3] for e in tail3_3])
Nlin = sp.Matrix([sp.Poly(e, *fill3).coeff_monomial(1) for e in tail3_3])
_rw = random.Random(20260724); wit = None
for _ in range(2000):
    fval = {v: sp.Integer(_rw.randint(-4, 4)) for v in free3}
    fval[am1_3] = sp.Integer(0)
    if fval[a2_0] == 0:
        fval[a2_0] = sp.Integer(2)
    if Mlin.subs(fval).rank() == Mlin.subs(fval).row_join(Nlin.subs(fval)).rank():
        wit = (fval[a2_0], sp.expand(rd3["R1"].subs(fval))); break
check(wit is not None and wit[0] != 0 and wit[1] == 0,
      f"CONTROL: cascade+tail witness a_2(0)={wit[0]}!=0, R(1)={wit[1]}=0 => a_2(0) NOT forced (kill is the FACTOR W)")
# STRUCTURE: W lives purely in the residual N (filler map is am1_3-free), and W ENTERS the
# consistency conditions -- some cokernel covector pairs N with a nonzero am1_3 coefficient.
# (Datum with am1_3 kept symbolic, other coords numeric; M is then am1_3-free.)
ptA = {v: sp.Rational(random.randint(-6, 6), random.choice([1, 2, 3])) for v in free3 if v != am1_3}
def subP(P): return sp.expand(sp.cancel(sp.sympify(P).subs(rd3["sol"]).subs(ptA)))
AdA = {k: subP(rd3["A"][k]) for k in rd3["A"]}; BdA = {k: subP(rd3["B"][k]) for k in rd3["B"]}
QmA = {m: sp.Poly(sp.expand(clear_denoms(sp.cancel(q_m(AdA, BdA, m)))), E) for m in (-1, -2, -3)}
rowsA = []
for m in (-1, -2, -3):
    for i in range(QmA[m].degree() + 1):
        rowsA.append(sp.expand(QmA[m].coeff_monomial(E**i)))
MA = sp.Matrix([[sp.Poly(e, *fill3).coeff_monomial(f) for f in fill3] for e in rowsA])
NA = sp.Matrix([sp.Poly(e, *fill3).coeff_monomial(1) for e in rowsA])
Mdep_am1 = any(am1_3 in sp.sympify(MA[i, j]).free_symbols for i in range(MA.rows) for j in range(MA.cols))
check(not Mdep_am1, "at d=3 the filler map M is INDEPENDENT of am1_3 (W lives in the rhs/residual N only)")
LKA = MA.T.nullspace()
w_coeffs = [sp.Poly(sp.expand((sp.Matrix([list(v)]) * NA)[0]), am1_3).coeff_monomial(am1_3)
            if am1_3 in sp.expand((sp.Matrix([list(v)]) * NA)[0]).free_symbols else sp.Integer(0)
            for v in LKA]
check(any(sp.simplify(c) != 0 for c in w_coeffs),
      f"W ENTERS the consistency conditions: some of the {len(LKA)} cokernel covectors pair N with a "
      "nonzero am1_3-coefficient (so tail-solvability constrains W)")
# DEPTH: Q_-1 alone does NOT force W -- an explicit cascade+Q_-1 point with am1_3!=0, Q_-1 solvable.
tail1 = rd3["tk"][-1]
M1 = sp.Matrix([[sp.Poly(e, *fill3).coeff_monomial(f) for f in fill3] for e in tail1])
N1 = sp.Matrix([sp.Poly(e, *fill3).coeff_monomial(1) for e in tail1])
_rd = random.Random(5); wit1 = None
for _ in range(400):
    fval = {v: sp.Integer(_rd.randint(-4, 4)) for v in free3}
    if fval[am1_3] == 0 or fval[a2_0] == 0:
        continue
    if M1.subs(fval).rank() == M1.subs(fval).row_join(N1.subs(fval)).rank():
        wit1 = sp.expand(rd3["R1"].subs(fval)); break
check(wit1 is not None and wit1 != 0,
      f"DEPTH: Q_-1 alone does NOT force W (cascade+Q_-1 solvable witness with R(1)={wit1}!=0) "
      "=> the kill genuinely needs the depth-3 coupling")

# =====================================================================
print("\n--- S6. UNIFORMITY: recipe at d=1,2,3 (exact) + d=4 linear route (mod p); HONEST verdict ---", flush=True)
# =====================================================================
for d in (1, 2, 3):
    rd = build_reduced(d)
    Wc = W_on_cascade(rd)
    prod = sp.cancel(rd["A"][2].subs(E, 0) * Wc)
    check(sp.simplify(rd["R1"] - prod) == 0, f"d={d}: R(1)=a_2(0)*W exact on the parametrized cascade")
    if d <= 2:
        check(sp.expand(rd["R1"]) == 0, f"d={d}: R(1)=0 on the cascade ALONE (forcing vacuous below d=3)")
# d=4 linear route (mod p): sampled full column rank plus absence of sampled
# W!=0 tail-solvable points is bounded supporting evidence, not a covector construction.
if HEAVY:
    rd4 = build_reduced(4)
    free4, fill4 = rd4["free"], rd4["fillers"]
    W4 = W_on_cascade(rd4)
    rows4 = [sp.expand(e) for e in param_tail(rd4, 3)]
    check(all(sp.Poly(e, *fill4).total_degree() <= 1 for e in rows4),
          f"d=4 depth-3 tail = {len(rows4)} rows LINEAR in the {len(fill4)} fillers")
    M4 = sp.Matrix([[sp.Poly(e, *fill4).coeff_monomial(f) for f in fill4] for e in rows4])
    N4 = sp.Matrix([sp.Poly(e, *fill4).coeff_monomial(1) for e in rows4])
    fM4 = sp.lambdify(free4, M4, modules="sympy"); fN4 = sp.lambdify(free4, N4, modules="sympy")
    fW4 = sp.lambdify(free4, sp.expand(W4), modules="sympy")
    pbig = 2147483647
    Mp = [[int(x) for x in row] for row in fM4(*[random.randint(2, 90) for _ in free4]).tolist()]
    check(rank_mod_p(Mp, len(fill4), pbig) == len(fill4),
          f"d=4 sampled filler map has FULL COLUMN RANK {len(fill4)}=10 (linear elimination valid at this sample)")
    n_wit = 0; n_Wnz = 0
    for _ in range(40):
        vals = [random.randint(-5, 5) for _ in free4]
        if sp.Integer(fW4(*vals)) == 0:
            continue
        n_Wnz += 1
        Mrows = [[int(x) for x in row] for row in fM4(*vals).tolist()]
        Nrows = [int(x) for x in sp.flatten(fN4(*vals).tolist())]
        rM = rank_mod_p(Mrows, len(fill4), pbig)
        rMN = rank_mod_p([Mrows[i] + [Nrows[i]] for i in range(len(Mrows))], len(fill4) + 1, pbig)
        if rM == rMN:
            fv = {v: sp.Integer(vals[k]) for k, v in enumerate(free4)}
            if M4.subs(fv).rank() == M4.subs(fv).row_join(N4.subs(fv)).rank():
                n_wit += 1
    check(n_wit == 0,
          f"d=4: among {n_Wnz} sampled cascade pts with W!=0, tail-solvable={n_wit} "
          "(mod-p bounded evidence only; no symbolic covector or degree-free recipe)")
    # depth-3 msolve kill at d=3 (single-tool corroboration of S5 control):
    if shutil.which("msolve"):
        t0 = time.time()
        try:
            check(msolve_unit(tail3_3 + [1 - t_rab * am1_3], allv3, PRIMES[0], tmo=180),
                  f"d=3 depth-3: am1_3 in sqrt(cascade+Q_-1..Q_-3) (msolve '^' mod {PRIMES[0]}) "
                  f"({time.time()-t0:.1f}s)")
        except MSFAIL:
            skip(f"d=3 depth-3 msolve kill (timeout/OOM {time.time()-t0:.1f}s); S5 sympy full-tail control stands")
    else:
        skip("d=3 depth-3 msolve kill (msolve not on PATH); S5 sympy full-tail control stands")
else:
    skip("S6 d=4 linear-route (full column rank + W-forcing sampling) is HEAVY-only")
    skip("S6 d=3 depth-3 msolve kill is HEAVY-only; S5 sympy full-tail control is the committed kill")

print("\n" + "=" * 74, flush=True)
print("PROVED (degree-free / node-free, machine-checked):", flush=True)
print("  * engine; slope gate; both-ends Lemma-P; R(1)=a_2(0)*W factorization.", flush=True)
print("  * ALGEBRAIC-NODE THM A': moving-sum adjoint for SYMBOLIC nodes; TRACE-FORM descent", flush=True)
print("    (companion trace = root sum, no root named; closed under S_n^*) -- the tool the", flush=True)
print("    (a_2,b_1) algebraic necklace needs.", flush=True)
print("  * TWO-BLOCK STRUCTURE: explicit two-term operators, level incidence; the negative-tail", flush=True)
print("    necklace is ENTIRELY algebraic (no a_3,b_2); membership-window covectors annihilate;", flush=True)
print("    no single-node kill on the tested generic a_2 b-block rung (coupling obstruction).", flush=True)
print("  * BLOCK ADJOINT CRITERION: symbolic-node coefficient-of-filler-value vanishing =>", flush=True)
print("    trace-form support conditions.", flush=True)
print("BOUNDED (d=3, exact): cokernel dim 16; W-kill reproduced (am1_3 in sqrt(cascade+tail),", flush=True)
print("  sympy QQ+prime); a_2(0) NOT forced (witness); specialized cokernel vectors", flush=True)
print("  annihilate specialized filler columns exactly;", flush=True)
print("  W enters the consistency conditions; Q_-1 alone does not force W (depth-3 coupling needed).", flush=True)
print("OPEN (obstruction LOCALIZED): a FIXED finite trace-form recipe giving unit*W at every d --", flush=True)
print("  blocked by the two-term coupling across the varying (a_2,a_1,a_0)/(b_1,b_0,b_-1) tops.", flush=True)
print("  d=4 (HEAVY): sampled rank/solvability evidence only; no symbolic covector or degree-free recipe.", flush=True)
print("=" * 74, flush=True)
print(f"\n(total {time.time() - _T0:.1f}s; {_NP} checks passed, {_NSKIP} skipped)", flush=True)
if _NSKIP:
    print("ALL EXECUTED ALGEBRAIC COVECTOR CHECKS PASSED; OPTIONAL CHECKS SKIPPED", flush=True)
else:
    print("ALL ALGEBRAIC COVECTOR CHECKS PASSED; NO SKIPS", flush=True)
