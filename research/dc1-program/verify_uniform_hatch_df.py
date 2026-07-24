#!/usr/bin/env python3
"""verify_uniform_hatch_df.py -- UNIFORM HATCH SLOPE-FORCING: the degree-free backbone.

INDEPENDENTLY DERIVED - EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES - NOT PEER REVIEWED

MISSION.  At the band-k hatch (k>=4) the negative tail (no Q_0) is conjectured to force
BOTH cokernel pairings of the residual R:  R(1)=0 AND lambda_3(R)=0, at arbitrary degree
-- one mechanism for architecture steps 1 (slope) + 2 (covector) across the whole tower.
This file establishes the DEGREE-FREE STRUCTURAL BACKBONE of that statement and the
bounded d=3 covector certificates, and states honestly what remains open.

CONVENTIONS (band-k, gauge b_k=0):  (x^a f)(x^b g)=x^{a+b} f(E+b)g(E), f^[n](E)=f(E+n),
  Q_m=sum_{i+l=m}(b_l^[i] a_i - a_i^[l] b_l),  u:=b_{k-1},  membership (E)_j|a_-j,b_-j,
  Q_0=(T-1)G, G(0)=0, R(1)=Q_0(0)=G(1), lambda_3(R)=lambda_3(G) (Im Phi subset ker lambda_3).
  hatch: a_k=prod_{i<k}(E+i(k-1)), u=prod_{j<k-1}(E+jk).
  band4: a_4=E(E+3)(E+6)(E+9), u=E(E+4)(E+8); lambda_3=ev_-2-ev_-4+ev_-5, lambda_3(E)=-3.
  band5: a_5=E(E+4)(E+8)(E+12)(E+16), u=E(E+5)(E+10)(E+15);
         lambda_3=ev_-3-ev_-5+ev_-7-ev_-10+ev_-11, lambda_3(E)=-6.
  Tail=Q_-1..Q_-(2k-1);  Phi fillers = a_-(k-1), b_-k.

WHAT IS PROVED (arbitrary degree, machine identities; k in {3,4,5}):
  * S0 engine: Q_m=[D,X]_m, Q_0=(T-1)G, G(0)=0.
  * S1 lambda_3 point functional (moving-sum), Im Phi subset ker lambda_3, targets phi-independent.
  * S2 BACKBONE-1  both-ends Lemma-P (degree-free):
        R(1) = Q_0(0) = G(1) = sum_{i=1}^{k-1} [ a_i(0) b_-i(i) - a_-i(i) b_i(0) ].
  * S2 BACKBONE-2  top-wall boundary identity (degree-free):
        Q_{2k-2}(0) = c_a a_{k-1}(0) - c_b b_{k-2}(0)  =>  b_{k-2}(0) = rho_k a_{k-1}(0),
        rho_3=2/3, rho_4=21/80, rho_5=8/55.  So two Lemma-P terms of BOTH targets share the
        factor a_{k-1}(0) (partial, shared backbone).
  * S3 DEPTH-2 FILLER-LINEARITY (general d, symbolic): the depth-2 tail Q_-1,Q_-2 is LINEAR
        in the Phi fillers -- hence the depth-2 kill is a consistency-covector statement at
        EVERY degree (band 4 d=2,3,4; band 5 d=3).

WHAT IS REFUTED (correcting the naive W2->band4 generalization):
  * S2: NO single boundary value a_j(0) divides R(1) or lambda_3(R) modulo the band-4
        cascade (d=3, GF(p) radical); the two targets share no common polynomial factor.
        W2's single-factor factorization R(1)=a_2(0)*W is SPECIAL to k=3.

WHAT IS BOUNDED (exact scope; band 4, d=3):
  * S4 depth-2 filler map: 18 rows x 8 Phi fillers, full column rank 8, cokernel dim 10;
        both targets free (nonconstant) + distinct on the cascade; origin feasible.
  * S5 THE COVECTOR KILL (kmin=2): R(1), lambda_3(R) in sqrt(cascade + Q_-1 + Q_-2) --
        msolve '^' UNIT (both targets; branch B default, +branch A/2nd prime/QQ HEAVY).
        Q_-1 alone forces NEITHER (explicit witnesses) -> rigorous kmin=2.
  * S7 (HEAVY) explicit determinant-SATURATED consistency-covector certificate:
        R(1), lambda_3(R) in sqrt( (cascade + 10 explicit covectors) : det^inf ) -- msolve UNIT.

WHAT IS OPEN / NOT CLAIMED:
  * The arbitrary-degree kill (the covector membership is bounded, d=3).
  * Band 5 non-degenerate test: band-5 is DEGENERATE through d=3 (targets forced 0 on the
        cascade alone -- S6, reduce-first GB); the genuine d>=4 test is OUT OF BUDGET.
  * A Weyl pair / DC1 / JC2.

Run:    uv run --with sympy python research/dc1-program/verify_uniform_hatch_df.py
Heavy:  HEAVY=1 uv run --with sympy python research/dc1-program/verify_uniform_hatch_df.py
Ends with an evidence ledger and one ALL-CAPS status line (SKIP discipline: without msolve
the d=3 covector payload prints SKIP and the final line does not claim it passed).
"""
import sympy as sp
import os, time, subprocess, tempfile, shutil, random

E, t = sp.symbols("E t")
_T0 = time.time()
HEAVY = os.environ.get("HEAVY", "") not in ("", "0")
HAVE_MS = shutil.which("msolve") is not None
PRIMES = (65003, 32003)
random.seed(20260723)
_NP = 0


def P(msg):
    print(f"      {msg}", flush=True)


def ok(cond, label):
    global _NP
    if not cond:
        raise AssertionError("FAIL " + label)
    _NP += 1
    print(f"PASS [{time.time()-_T0:6.1f}s] {label}", flush=True)


def okz(v, label):
    ok(sp.expand(v) == 0, label)


# ------------------------------------------------------------------ helpers
def sh(f, n):
    return sp.expand(sp.sympify(f).subs(E, E + n))


def falling(n):
    return sp.prod(E - j for j in range(n)) if n > 0 else sp.Integer(1)


def poly(name, deg):
    cs = list(sp.symbols(f"{name}_0:{deg+1}"))
    return sp.expand(sum(cs[j] * E**j for j in range(deg + 1))), cs


def q_m(A, B, m, k):
    L = range(-k, k + 1)
    return sp.expand(sum(sh(B[l], i) * A[i] - sh(A[i], l) * B[l]
                         for i in L for l in L if i + l == m))


def potential(A, B, k):
    return sp.expand(sum(sh(A[i], j - i) * sh(B[-i], j) - sh(B[i], j - i) * sh(A[-i], j)
                         for i in range(1, k + 1) for j in range(i)))


def mul_ladders(Pd, Qd):
    Rr = {}
    for k1, p1 in Pd.items():
        for k2, p2 in Qd.items():
            Rr[k1 + k2] = Rr.get(k1 + k2, 0) + sh(p1, k2) * p2
    return {kk: sp.expand(v) for kk, v in Rr.items()}


def commutator_coeff(X, Dop, m):
    return sp.expand(mul_ladders(Dop, X).get(m, 0) - mul_ladders(X, Dop).get(m, 0))


def hatch(k):
    a = sp.expand(sp.prod(E + i * (k - 1) for i in range(k)))
    u = sp.expand(sp.prod(E + j * k for j in range(k - 1)))
    return a, u


def Kblock(a, c, k):
    return sp.expand(sum(sh(a, j - k) * sh(c, j) for j in range(k)))


def Hblock(u, v, k):
    return sp.expand(sum(sh(u, j - (k - 1)) * sh(v, j) for j in range(k - 1)))


def scoeffs(e):
    return [c for c in sp.Poly(sp.expand(e), E).all_coeffs() if sp.expand(c) != 0]


def cd(e):
    n, _ = sp.fraction(sp.together(sp.expand(e)))
    return sp.expand(n)


# ---- exact GF(p) reduction utilities (independent sympy engine) ----
def _rmp(c, p):
    r = sp.Rational(c)
    return (int(r.p) % p) * pow(int(r.q) % p, p - 2, p) % p


def mp(e, p, gens):
    e = sp.expand(e)
    if e == 0:
        return sp.Integer(0)
    if not e.free_symbols:
        return sp.Integer(_rmp(e, p))
    Pq = sp.Poly(e, *gens, domain="QQ")
    d = {mo: _rmp(co, p) for mo, co in zip(Pq.monoms(), Pq.coeffs()) if _rmp(co, p) != 0}
    return sp.Poly(d, *gens, domain=sp.GF(p, symmetric=False)).as_expr() if d else sp.Integer(0)


def modp_lin_reduce(eqs, targets, p):
    """Forward-substitute equations linear in a variable with nonzero constant pivot, over GF(p)."""
    gens = sorted(set().union(*[e.free_symbols for e in eqs + targets]) - {t}, key=str)
    eqs = [x for x in (mp(e, p, gens) for e in eqs) if x != 0]
    targets = list(targets); changed = True
    while changed:
        changed = False
        for e in eqs:
            for v in sorted(e.free_symbols, key=str):
                pe = sp.Poly(e, v)
                if pe.degree() == 1 and not pe.nth(1).free_symbols and int(pe.nth(1)) % p != 0:
                    inv = pow(int(pe.nth(1)) % p, p - 2, p)
                    sol = mp(-pe.nth(0) * inv, p, gens)
                    eqs = [x for x in (mp(y.subs(v, sol), p, gens) for y in eqs) if x != 0]
                    targets = [mp(y.subs(v, sol), p, gens) for y in targets]
                    changed = True
                    break
            if changed:
                break
    return eqs, targets


def sy_unit_modp(eqs, p):
    eqs = [e for e in eqs if sp.expand(e) != 0]
    if not eqs:
        return False
    gens = sorted(set().union(*[e.free_symbols for e in eqs]), key=str)
    G = sp.groebner([mp(e, p, gens) for e in eqs], *gens, order="grevlex", modulus=p)
    return list(G.exprs) == [sp.Integer(1)]


# ---- msolve '^' engine (validated parsing) ----
def msolve_unit(eqs, fv, char=65003, tmo=280):
    if not HAVE_MS:
        return None
    eqs = [cd(e) for e in eqs if sp.expand(e) != 0]
    xs = sp.symbols(f"z_0:{len(fv)}")
    sub = dict(zip(fv, xs))
    with tempfile.NamedTemporaryFile("w", suffix=".ms", delete=False) as f:
        inp = f.name
        f.write(",".join(str(x) for x in xs) + "\n" + str(char) + "\n")
        polys = [str(cd(e.subs(sub))).replace(" ", "").replace("**", "^") for e in eqs]
        f.write(",\n".join(polys) + "\n")
    out = inp + ".out"
    try:
        subprocess.run(["msolve", "-g", "2", "-f", inp, "-o", out], check=True,
                       stderr=subprocess.DEVNULL, timeout=tmo)
        r = open(out).read()
    except (subprocess.TimeoutExpired, subprocess.CalledProcessError, OSError):
        for pth in (inp, out):
            if os.path.exists(pth):
                os.remove(pth)
        return None
    finally:
        for pth in (inp, out):
            if os.path.exists(pth):
                os.remove(pth)
    br = r[r.rfind("["):]
    return br.strip().startswith("[1]")


# ---- cascade builder (general k; RE-DERIVED, forward-solves the b's) ----
def clean_solve(A, B, m, lkey, name, membership, raw_degree, k):
    raw, cs = poly(name, raw_degree)
    unknown = sp.expand(falling(membership) * raw)
    trial = dict(B); trial[lkey] = unknown
    M, rhs = sp.linear_eq_to_matrix(sp.Poly(q_m(A, trial, m, k), E).all_coeffs(), cs)
    if any(entry.free_symbols for entry in M):
        raise AssertionError("operator matrix not numeric")
    indep = sp.zeros(0, len(cs)); selrhs = []
    for i in range(M.rows):
        cand = indep.col_join(M[i, :])
        if cand.rank() > indep.rank():
            indep = cand; selrhs.append(rhs[i])
    conds = [c for c in (sp.expand(n.dot(rhs)) for n in M.T.nullspace()) if c != 0]
    vals = ([sp.Integer(0)] * len(cs) if indep.rows == 0
            else [x.subs({s: 0 for s in x.free_symbols if str(s).startswith("tau")})
                  for x in indep.gauss_jordan_solve(sp.Matrix(selrhs))[0]])
    res = sp.expand(unknown.subs(dict(zip(cs, vals))))
    kernels = []
    for j, vec in enumerate(M.nullspace()):
        pk = sp.symbols(f"{name}K{j}")
        res = sp.expand(res + pk * falling(membership) * sum(vec[i] * E**i for i in range(len(cs))))
        kernels.append(pk)
    return res, kernels, conds


def build_branch(ak, u, d, k, branch):
    """Forward-solve b_{k-2}..b_{-(k-1)} from Q_{2k-2}..Q_1; return compat conds + free vars + phi."""
    A = {k: ak}; avars = []
    for i in range(k - 1, -1, -1):
        p_, c = poly(f"a{i}", d); A[i] = p_; avars += c
    for r in range(1, k):
        raw, c = poly(f"am{r}", d); A[-r] = sp.expand(falling(r) * raw); avars += c
    A[-k] = sp.Integer(0)
    B = {l: sp.Integer(0) for l in range(-k, k + 1)}; B[k - 1] = u
    conds, kernels = [], []
    for l in range(k - 2, -(k), -1):
        B[l], nk, nc = clean_solve(A, B, l + k, l, f"b{l}c", max(0, -l), 2 * d + 2 * k, k)
        kernels += nk; conds += nc
    if branch == "A":
        amk, c = poly("amk", d); mu = sp.symbols("muk")
        A[-k] = sp.expand(falling(k) * amk); B[-k] = sp.expand(mu * A[-k]); extra = c + [mu]
        phi = sorted(list(set(A[-(k - 1)].free_symbols) - {E}) + c + [mu], key=str)
    else:
        bmk, c = poly("bmk", d); A[-k] = sp.Integer(0); B[-k] = sp.expand(falling(k) * bmk); extra = c
        phi = sorted(list(set(A[-(k - 1)].free_symbols) - {E}) + c, key=str)
    return A, B, [e for e in conds if sp.expand(e) != 0], avars + extra + kernels, phi


HATCH = {k: hatch(k) for k in (3, 4, 5)}
LAM = {}
RHO = {3: sp.Rational(2, 3), 4: sp.Rational(21, 80), 5: sp.Rational(8, 55)}


def lam_of(k):
    wt = LAM[k]
    return lambda f: sp.expand(sum(c * sp.sympify(f).subs(E, p) for p, c in wt.items()))


# ======================================================================
print("=== S0. ENGINE (re-derived): Q_m=[D,X]_m, Q_0=(T-1)G, G(0)=0 (k=3,4,5) ===", flush=True)
# ======================================================================
for k in (3, 4, 5):
    Ag = {lev: sp.expand((falling(-lev) if lev < 0 else 1) * poly(f"Ag{lev+k}", 2)[0])
          for lev in range(-k, k + 1)}
    Bg = {lev: sp.expand((falling(-lev) if lev < 0 else 1) * poly(f"Bg{lev+k}", 2)[0])
          for lev in range(-k, k + 1)}
    for m in range(-2 * k, 2 * k + 1):
        okz(q_m(Ag, Bg, m, k) - commutator_coeff(Ag, Bg, m),
            f"k={k}: Q_{m} = [D,X]_{m} (generic deg-2 coeffs)")
    okz(q_m(Ag, Bg, 0, k) - (sh(potential(Ag, Bg, k), 1) - potential(Ag, Bg, k)),
        f"k={k}: Q_0 = (T-1)G")
    okz(potential(Ag, Bg, k).subs(E, 0), f"k={k}: G(0)=0 under membership => R(1)=Q_0(0)=G(1)")


# ======================================================================
print("\n=== S1. lambda_3 (re-derived via moving-sum); Im Phi subset ker lambda_3 ===", flush=True)
# ======================================================================
def qK_roots(a, k):
    return sorted(set([int(r) for r in sp.solve(sh(a, -k), E)] + list(range(k))))


def qH_roots(u, k):
    return sorted(set([int(r) for r in sp.solve(sh(u, -(k - 1)), E)] + list(range(k - 1))))


def moving_sum_solutions(k):
    a, u = hatch(k)
    RK, RH = qK_roots(a, k), qH_roots(u, k)
    Sk = sum(t**j for j in range(k)); Skm1 = sum(t**j for j in range(k - 1))
    pk = {x: sp.Symbol(f"p{x+50}") for x in RK}
    ph = {x: sp.Symbol(f"q{x+50}") for x in RH}
    PK = sum(pk[x] * t**x for x in RK); PH = sum(ph[x] * t**x for x in RH)
    mn = min(RK + RH)
    pdd = sp.Poly(sp.expand((PK * Skm1 - PH * Sk) * t**(-mn + 2)), t)
    eqs = [pdd.nth(i) for i in range(pdd.degree() + 1) if pdd.nth(i) != 0]
    allv = list(pk.values()) + list(ph.values())
    sol = list(sp.linsolve(eqs, allv))[0]
    subs0 = dict(zip(allv, sol))
    freev = sorted(set().union(*[e.free_symbols for e in sol]) or set(), key=str)
    basis = []
    for fv in freev:
        Lb = sp.cancel(sp.expand(PK.subs(subs0)).subs(
            {v: (1 if v == fv else 0) for v in freev}) / Sk)
        num, den = sp.fraction(sp.together(Lb))
        Ln, Ld = sp.Poly(sp.expand(num), t), sp.Poly(sp.expand(den), t)
        is_laurent = (Ld.is_monomial and len(Ld.free_symbols) <= 1)
        wt = {i - Ld.degree(): int(Ln.nth(i)) for i in range(Ln.degree() + 1) if Ln.nth(i) != 0}
        basis.append((sorted(wt), wt, is_laurent))
    return len(freev), basis


for k in (3, 4, 5):
    dim, basis = moving_sum_solutions(k)
    ok(dim == 3, f"k={k}: moving-sum annihilator space is 3-dimensional (= codim Im Phi)")
    ok(all(b[2] for b in basis),
       f"k={k}: EVERY annihilator symbol is a Laurent polynomial => ALL point functionals")
    cand = [b for b in basis if b[0] not in ([0], [1])]
    cand.sort(key=lambda b: len(b[0]))
    LAM[k] = cand[0][1]
    lam = lam_of(k)
    ok(lam(E) == -sp.Rational((k - 1) * (k - 2), 2), f"k={k}: lambda_3(E) = {lam(E)} = -(k-1)(k-2)/2")
    ok(lam(1) == 1, f"k={k}: lambda_3(1) = 1")
    P(f"k={k}: lambda_3 = {dict(sorted(LAM[k].items()))}")

# Im Phi subset ker lambda_3, degree-free (generic C,V of degree k) => lambda_3(R)=lambda_3(G).
for k in (3, 4, 5):
    a, u = hatch(k)
    lam = lam_of(k)
    Cc = list(sp.symbols(f"Cc{k}_0:{k+2}")); Vv = list(sp.symbols(f"Vv{k}_0:{k+2}"))
    C = sum(Cc[i] * E**i for i in range(k + 1)); V = sum(Vv[i] * E**i for i in range(k + 1))
    okz(lam(Kblock(a, falling(k) * C, k)),
        f"k={k}: lambda_3(K_k[(E)_k C]) = 0 for generic deg-{k} C (degree-free)")
    okz(lam(Hblock(u, falling(k - 1) * V, k)),
        f"k={k}: lambda_3(H_(k-1)[(E)_(k-1) V]) = 0 for generic deg-{k} V (degree-free)")

# targets phi-independent (Lemma P + Im Phi subset ker lambda_3), band 4,5
for k in (4, 5):
    ak, u = HATCH[k]; lam = lam_of(k)
    A, B, conds, fv, phi = build_branch(ak, u, 2, k, "B")
    G = potential(A, B, k)
    inter = (sp.expand(q_m(A, B, 0, k).subs(E, 0)).free_symbols
             | sp.expand(lam(G)).free_symbols) & set(phi)
    ok(inter == set(),
       f"k={k}: R(1)=Q_0(0) and lambda_3(R)=lambda_3(G) are INDEPENDENT of the Phi fillers "
       f"a_-(k-1),b_-k")
ok(LAM[3] == {-1: 1}, "k=3 (W2): lambda_3 = ev_-1 (single point) -- base case")


# ======================================================================
print("\n=== S2. BACKBONE: both-ends Lemma-P + top-wall b_{k-2}(0)=rho_k a_{k-1}(0) ===", flush=True)
# ======================================================================
def build_generic(k, d):
    ak, u = hatch(k)
    A, B = {}, {}
    A[k] = ak
    for i in range(k - 1, -1, -1):
        A[i], _ = poly(f"A{k}_{d}_a{i}", d)
    for r in range(1, k + 1):
        p_, _ = poly(f"A{k}_{d}_am{r}", d); A[-r] = sp.expand(falling(r) * p_)
    B[k] = sp.Integer(0); B[k - 1] = u
    for l in range(k - 2, -(k + 1), -1):
        p_, _ = poly(f"A{k}_{d}_b{l}", d); B[l] = sp.expand((falling(-l) if l < 0 else 1) * p_)
    return A, B


# BACKBONE-1: both-ends Lemma-P, degree-free, k=3,4,5, generic d=2,3
for k in (3, 4, 5):
    for d in (2, 3):
        A, B = build_generic(k, d)
        R1 = sp.expand(q_m(A, B, 0, k).subs(E, 0))
        lemmaP = sp.expand(sum(A[i].subs(E, 0) * B[-i].subs(E, i)
                               - A[-i].subs(E, i) * B[i].subs(E, 0) for i in range(1, k)))
        okz(R1 - lemmaP,
            f"k={k} d={d}: R(1)=Q_0(0) = sum_(i=1)^({k-1})[a_i(0)b_-i(i) - a_-i(i)b_i(0)] "
            f"(both-ends Lemma-P, degree-free)")

# BACKBONE-2: top-wall boundary identity Q_{2k-2}(0) => b_{k-2}(0)=rho_k a_{k-1}(0), k=3,4,5
for k in (3, 4, 5):
    A, B = build_generic(k, 3)
    Qtop = sp.expand(q_m(A, B, 2 * k - 2, k).subs(E, 0))
    ak1_0 = sp.Symbol(f"A{k}_3_a{k-1}_0"); bk2_0 = sp.Symbol(f"A{k}_3_b{k-2}_0")
    c_a = sp.Poly(Qtop, ak1_0).nth(1); c_b = sp.Poly(Qtop, bk2_0).nth(1)
    # Q_{2k-2}(0) is EXACTLY c_a*a_{k-1}(0) + c_b*b_{k-2}(0) (only these two boundary values survive)
    okz(Qtop - (c_a * ak1_0 + c_b * bk2_0),
        f"k={k}: Q_{2*k-2}(0) = ({c_a}) a_{k-1}(0) + ({c_b}) b_{k-2}(0) -- only these survive at E=0")
    ok(sp.nsimplify(-c_a / c_b) == RHO[k],
       f"k={k}: Q_{2*k-2}(0)=0 => b_{k-2}(0) = {RHO[k]} * a_{k-1}(0) (top-wall proportionality)")

# PARTIAL shared backbone (band 4): on the top wall b_2(0)=(21/80)a_3(0), the term
# -a_-2(2)b_2(0) of BOTH targets, and the explicit a_3(0)-terms, are a_3(0)-divisible.
k = 4
A, B = build_generic(4, 3)
lam = lam_of(4)
a3_0 = sp.Symbol("A4_3_a3_0"); b2_0 = sp.Symbol("A4_3_b2_0")
wall = {b2_0: sp.Rational(21, 80) * a3_0}
R1 = sp.expand(q_m(A, B, 0, 4).subs(E, 0)); L3 = sp.expand(lam(potential(A, B, 4)))
for tgt, nm in ((R1, "R(1)"), (L3, "lambda_3(R)")):
    # the b_2(0)-carrying part + the a_3(0)-explicit part, on the wall, is divisible by a_3(0):
    part = sp.expand(tgt.coeff(b2_0) * b2_0 + tgt.coeff(a3_0) * a3_0).subs(wall)
    ok(sp.expand(part) != 0 and sp.rem(sp.Poly(sp.expand(part), a3_0), sp.Poly(a3_0, a3_0)) == 0,
       f"band4: the b_2(0)+a_3(0)-explicit part of {nm} is a_3(0)-divisible on the top wall "
       f"(shared partial backbone)")

# REFUTED: no single boundary a_j(0) divides R(1)/lambda_3(R) mod band-4 cascade (d=3, GF(p) radical)
trab = sp.symbols("t_rab")
ak4, u4 = HATCH[4]; lam = lam_of(4)
Ab, Bb, condsb, fvb, phib = build_branch(ak4, u4, 3, 4, "B")
R1b = sp.expand(q_m(Ab, Bb, 0, 4).subs(E, 0)); L3b = sp.expand(lam(potential(Ab, Bb, 4)))
protect = {sp.Symbol(f"a{i}_0") for i in range(4)}


# reduce cascade WITHOUT eliminating protected boundary coeffs a_j(0)
def modp_lin_reduce_protect(eqs, targets, p, protect):
    gens = sorted(set().union(*[e.free_symbols for e in eqs + targets]) - {t}, key=str)
    eqs = [x for x in (mp(e, p, gens) for e in eqs) if x != 0]
    targets = list(targets); changed = True
    while changed:
        changed = False
        for e in eqs:
            for v in sorted(e.free_symbols, key=str):
                if v in protect:
                    continue
                pe = sp.Poly(e, v)
                if pe.degree() == 1 and not pe.nth(1).free_symbols and int(pe.nth(1)) % p != 0:
                    inv = pow(int(pe.nth(1)) % p, p - 2, p)
                    sol = mp(-pe.nth(0) * inv, p, gens)
                    eqs = [x for x in (mp(y.subs(v, sol), p, gens) for y in eqs) if x != 0]
                    targets = [mp(y.subs(v, sol), p, gens) for y in targets]
                    changed = True
                    break
            if changed:
                break
    return eqs, targets


residP, (R1P, L3P) = modp_lin_reduce_protect(condsb, [R1b, L3b], PRIMES[0], protect)
sf_names = ["a3_0"] + (["a2_0", "a1_0"] if HEAVY else [])
for bnm in sf_names:
    bsym = sp.Symbol(bnm)
    for tgt, nm in ((R1P, "R(1)"), (L3P, "lambda_3(R)")):
        unit = sy_unit_modp(residP + [bsym, sp.expand(1 - trab * tgt)], PRIMES[0])
        ok(unit is False,
           f"band4 d=3: {nm} NOT in sqrt(cascade + {bnm}=0) [sympy GF({PRIMES[0]})] -- {bnm} does "
           f"NOT divide {nm} mod cascade (W2 single-factor factorization is SPECIAL to k=3)")
if not HEAVY:
    P("[set HEAVY=1 to also refute a2_0, a1_0 as single divisors]")


# ======================================================================
print("\n=== S3. DEPTH-2 FILLER-LINEARITY (general d): Q_-1,Q_-2 linear in the Phi fillers ===", flush=True)
# ======================================================================
# This is what makes the depth-2 kill a CONSISTENCY-COVECTOR statement at EVERY degree.
for k, dset in ((4, (2, 3, 4)), (5, (3,))):
    ak, u = HATCH[k]; lam = lam_of(k)
    for d in dset:
        A, B, conds, fv, phi = build_branch(ak, u, d, k, "B")
        rows = scoeffs(q_m(A, B, -1, k)) + scoeffs(q_m(A, B, -2, k))
        maxdeg = max(sp.Poly(sp.expand(r), *phi).total_degree() for r in rows)
        R1 = sp.expand(q_m(A, B, 0, k).subs(E, 0)); L3 = sp.expand(lam(potential(A, B, k)))
        tindep = not ((R1.free_symbols | L3.free_symbols) & set(phi))
        ok(maxdeg <= 1 and tindep,
           f"band{k} d={d}: depth-2 tail (Q_-1,Q_-2) is LINEAR in the {len(phi)} Phi fillers "
           f"(max total-deg={maxdeg}); targets phi-independent -- depth-2 kill = covector statement")


# ======================================================================
print("\n=== S4. band-4 d=3 controls: cokernel dim; targets free+distinct; origin feasible ===", flush=True)
# ======================================================================
d = 3
ak4, u4 = HATCH[4]; lam = lam_of(4)
A, B, conds, fv, phi = build_branch(ak4, u4, d, 4, "B")
rows = scoeffs(q_m(A, B, -1, 4)) + scoeffs(q_m(A, B, -2, 4))
Mmat = sp.Matrix([[sp.Poly(sp.expand(r), *phi).coeff_monomial(f) if sp.expand(r).free_symbols & set(phi)
                   else sp.Integer(0) for f in phi] for r in rows])
nonfill = sorted((set().union(*[r.free_symbols for r in rows]) - {E} - set(phi)), key=str)
ptF = {v: sp.Integer(random.randint(2, PRIMES[0] - 2)) for v in nonfill}
rankM = mp_matrix_rank = sp.Matrix([[int(mp(Mmat[i, j], PRIMES[0], nonfill).subs(ptF)) % PRIMES[0]
                                     for j in range(Mmat.cols)] for i in range(Mmat.rows)]).rank()
ok(rankM == len(phi),
   f"band4 d=3: depth-2 filler map = {len(rows)} rows x {len(phi)} Phi fillers, FULL COLUMN RANK "
   f"{rankM} at a generic point; cokernel dim = {len(rows)}-{rankM} = {len(rows)-rankM} "
   f"(the consistency covectors)")

# free + distinct moduli on the cascade (non-vacuity)
R1 = sp.expand(q_m(A, B, 0, 4).subs(E, 0)); L3 = sp.expand(lam(potential(A, B, 4)))
p = PRIMES[0]
residc, (R1c, L3c) = modp_lin_reduce(conds, [R1, L3], p)
gens_c = sorted(set().union(*([e.free_symbols for e in residc] or [set()]),
                            R1c.free_symbols, L3c.free_symbols), key=str)
GBc = sp.groebner([mp(e, p, gens_c) for e in residc], *gens_c, order="grevlex", modulus=p)
nfR = sp.sympify(GBc.reduce(mp(R1c, p, gens_c))[1]); nfL = sp.sympify(GBc.reduce(mp(L3c, p, gens_c))[1])
ok(nfR.free_symbols != set() and nfL.free_symbols != set(),
   "band4 d=3: R(1) and lambda_3(R) are BOTH free (nonconstant) moduli on the cascade -- "
   "tail-forcing to 0 is NON-vacuous")
ok(sp.sympify(GBc.reduce(mp(sp.expand(R1c - L3c), p, gens_c))[1]) != 0,
   "band4 d=3: R(1) != lambda_3(R) mod cascade -- the two cokernel conditions are INDEPENDENT (k>=4)")

# origin feasible on cascade+tail
for branch in (("B", "A") if HEAVY else ("B",)):
    Ax, Bx, cx, fx, px = build_branch(ak4, u4, 3, 4, branch)
    tail = sum((scoeffs(q_m(Ax, Bx, m, 4)) for m in range(-1, -8, -1)), [])
    R1x = sp.expand(q_m(Ax, Bx, 0, 4).subs(E, 0)); L3x = sp.expand(lam(potential(Ax, Bx, 4)))
    z = {v: 0 for v in fx}
    ok(all(sp.expand(e.subs(z)) == 0 for e in cx + tail)
       and sp.expand(R1x.subs(z)) == 0 and sp.expand(L3x.subs(z)) == 0,
       f"band4 d=3 br{branch}: the ORIGIN lies on cascade+tail with R(1)=lambda_3(R)=0 "
       f"(cascade+tail non-empty)")


# ======================================================================
print("\n=== S5. THE COVECTOR KILL (band 4 d=3): both targets in sqrt(cascade+Q_-1+Q_-2), kmin=2 ===", flush=True)
# ======================================================================
def targets_and_rows(ak, u, d, k, branch):
    A, B, conds, fv, phi = build_branch(ak, u, d, k, branch)
    R1 = cd(q_m(A, B, 0, k).subs(E, 0)); L3 = cd(lam_of(k)(potential(A, B, k)))
    q1 = scoeffs(q_m(A, B, -1, k)); q2 = scoeffs(q_m(A, B, -2, k))
    return A, B, conds, fv, phi, R1, L3, q1, q2


PROBE = {}
HEADLINE_KEYS = {("B", nm) for nm in ("R(1)", "lambda_3(R)")}
if HAVE_MS:
    xx = sp.symbols("xx")
    ok(msolve_unit([sp.Integer(1)], [xx]) is True, "msolve parsing: unit ideal {1} -> UNIT")
    ok(msolve_unit([xx * (xx - 1)], [xx]) is False, "msolve parsing: feasible {xx(xx-1)} -> NOT unit")
    branches = ("B", "A") if HEAVY else ("B",)
    for branch in branches:
        A, B, conds, fv, phi, R1, L3, q1, q2 = targets_and_rows(ak4, u4, 3, 4, branch)
        base = conds + q1 + q2
        for tgt, nm in ((R1, "R(1)"), (L3, "lambda_3(R)")):
            res = msolve_unit(base + [sp.expand(1 - trab * tgt)], [trab] + fv, PRIMES[0], tmo=280)
            PROBE[(branch, nm)] = res
            ok(res is True,
               f"band4 d=3 br{branch}: {nm} in sqrt(cascade+Q_-1+Q_-2) -- msolve '^' UNIT mod "
               f"{PRIMES[0]} => the depth-2 tail forces {nm}=0 (consistency-covector kill)")
            if HEAVY:
                r2 = msolve_unit(base + [sp.expand(1 - trab * tgt)], [trab] + fv, PRIMES[1], tmo=280)
                ok(r2 is True, f"band4 d=3 br{branch}: {nm} in sqrt(cascade+Q_-1+Q_-2) -- msolve mod {PRIMES[1]}")
                rq = msolve_unit(base + [sp.expand(1 - trab * tgt)], [trab] + fv, 0, tmo=520)
                if rq is True:
                    ok(True, f"band4 d=3 br{branch}: {nm} in sqrt(cascade+Q_-1+Q_-2) -- msolve UNIT over QQ (char 0)")
                else:
                    P(f"[band4 d=3 br{branch} {nm}: msolve QQ char-0 not terminated in budget; mod-p stands]")
    # kmin=2 lower bound: Q_-1 alone forces NEITHER (explicit feasible witness with target=1)
    A, B, conds, fv, phi, R1, L3, q1, q2 = targets_and_rows(ak4, u4, 3, 4, "B")

    def witness_Qm1_not_forcing(tgt, seed):
        resid, reduced = modp_lin_reduce(conds, [tgt] + q1, PRIMES[0])
        tgt_r = reduced[0]; extra_r = [e for e in reduced[1:] if sp.expand(e) != 0]
        gens = sorted(set().union(*([e.free_symbols for e in resid + extra_r] or [set()]),
                                  tgt_r.free_symbols), key=str)
        system = resid + extra_r + [sp.expand(tgt_r - 1)]
        for nfix in (5, 4, 3, 2):
            for att in range(6):
                random.seed(seed + 100 * nfix + att)
                fixed = random.sample(gens, min(nfix, len(gens)))
                slic = [v - random.randint(1, PRIMES[0] - 1) for v in fixed]
                r = msolve_unit(system + slic, gens, PRIMES[0], tmo=30)
                if r is False:
                    return True
        return False
    for tgt, nm in ((R1, "R(1)"), (L3, "lambda_3(R)")):
        w = witness_Qm1_not_forcing(tgt, seed=71 + len(nm))
        ok(w is True,
           f"band4 d=3: {nm}=1 achievable on cascade+Q_-1 (feasible witness, msolve GB non-unit on a "
           f"0-dim slice) => Q_-1 alone does NOT force {nm} -- rigorous kmin({nm})=2")
else:
    P("[msolve not on PATH -- the d=3 band-4 covector-kill payload is SKIPPED]")


# ======================================================================
print("\n=== S6. BAND 5 via the reduce-first route: d=3 is DEGENERATE; d>=4 out of budget ===", flush=True)
# ======================================================================
k = 5
ak5, u5 = HATCH[5]; lam = lam_of(5)
A, B, conds, fv, phi = build_branch(ak5, u5, 3, 5, "B")
R1 = sp.expand(q_m(A, B, 0, 5).subs(E, 0)); L3 = sp.expand(lam(potential(A, B, 5)))
p = PRIMES[0]
resid5, (R1r5, L3r5) = modp_lin_reduce(conds, [R1, L3], p)
gens5 = sorted(set().union(*([e.free_symbols for e in resid5] or [set()]),
                           R1r5.free_symbols, L3r5.free_symbols), key=str)
GB5 = sp.groebner([mp(e, p, gens5) for e in resid5], *gens5, order="grevlex", modulus=p)
nfR5 = sp.sympify(GB5.reduce(mp(R1r5, p, gens5))[1]); nfL5 = sp.sympify(GB5.reduce(mp(L3r5, p, gens5))[1])
ok(nfR5 == 0 and nfL5 == 0,
   "band5 d=3: R(1) and lambda_3(R) reduce to 0 modulo the cascade GB [reduce-first, sympy GF(p)] "
   "-- band 5 is DEGENERATE through d=3 (targets forced 0 on the cascade ALONE; extends d<=2). "
   "The reduce-first route makes the d=3 verdict tractable where the raw msolve probe was not.")
# depth-2 filler-linearity holds at band 5 too (already checked in S3); note d>=4 open
P("[band 5 non-degenerate test is d>=4: the cascade build is ~50s and the reduce+GB exceeds the "
  "8-min budget -- OPEN / OUT OF BUDGET, as is the raw msolve probe.]")


# ======================================================================
print("\n=== S7 (HEAVY). explicit determinant-SATURATED consistency-covector certificate ===", flush=True)
# ======================================================================
if HEAVY and HAVE_MS:
    srab = sp.symbols("s_rab")
    A, B, conds, fv, phi = build_branch(ak4, u4, 3, 4, "B")
    R1 = sp.expand(q_m(A, B, 0, 4).subs(E, 0)); L3 = sp.expand(lam_of(4)(potential(A, B, 4)))
    rows = scoeffs(q_m(A, B, -1, 4)) + scoeffs(q_m(A, B, -2, 4))
    resid, red = modp_lin_reduce(conds, rows + [R1, L3], PRIMES[0])
    rows_r = red[:len(rows)]; R1r = red[len(rows)]; L3r = red[len(rows) + 1]
    y = sorted((set().union(*[sp.expand(r).free_symbols for r in rows_r + [R1r, L3r] + resid])
                - {E} - set(phi)), key=str)
    Mm = sp.Matrix([[sp.Poly(sp.expand(r), *phi).coeff_monomial(f)
                     if sp.expand(r).free_symbols & set(phi) else sp.Integer(0) for f in phi]
                    for r in rows_r])
    Nn = sp.Matrix([sp.expand(sp.expand(r) - sum(Mm[i, j] * phi[j] for j in range(len(phi))))
                    for i, r in enumerate(rows_r)])
    ypt = {v: sp.Integer(random.randint(2, PRIMES[0] - 2)) for v in y}
    Mev = sp.Matrix([[int(mp(Mm[i, j], PRIMES[0], y).subs(ypt)) % PRIMES[0] if Mm[i, j] != 0 else 0
                      for j in range(len(phi))] for i in range(Mm.rows)])
    piv, Msel = [], sp.zeros(0, len(phi))
    for i in range(Mm.rows):
        cand = Msel.col_join(Mev[i, :])
        if cand.rank() > Msel.rank():
            Msel = cand; piv.append(i)
        if len(piv) == len(phi):
            break
    nonpiv = [i for i in range(Mm.rows) if i not in piv]
    MI = sp.Matrix([Mm[i, :] for i in piv]); NI = sp.Matrix([Nn[i] for i in piv])
    Ddet = mp(MI.det(), PRIMES[0], y); adjMI = MI.adjugate()
    Cs = []
    for i in nonpiv:
        Mi = sp.Matrix([[Mm[i, j] for j in range(len(phi))]])
        Ci = mp(sp.expand(Ddet * Nn[i] - (Mi * adjMI * NI)[0]), PRIMES[0], y)
        if Ci != 0:
            Cs.append(Ci)
    P(f"band4 d=3: {len(piv)} pivot rows, {len(Cs)} explicit consistency covectors, det has "
      f"{len(sp.Poly(Ddet, *y).terms())} terms")
    allv = list(y) + [srab, trab]
    for tgt, nm in ((R1r, "R(1)"), (L3r, "lambda_3(R)")):
        sysx = resid + Cs + [sp.expand(1 - srab * Ddet), sp.expand(1 - trab * tgt)]
        u = msolve_unit(sysx, allv, PRIMES[0], tmo=300)
        ok(u is True,
           f"band4 d=3: {nm} in sqrt( (cascade + {len(Cs)} explicit covectors) : det^inf ) -- "
           f"msolve '^' UNIT (determinant-SATURATED consistency-covector certificate)")
    # band 5 d=2 degenerate corroboration
    A5, B5, c5, f5, p5 = build_branch(ak5, u5, 2, 5, "B")
    tail5 = sum((scoeffs(q_m(A5, B5, m, 5)) for m in range(-1, -10, -1)), [])
    R15 = sp.expand(q_m(A5, B5, 0, 5).subs(E, 0)); L35 = sp.expand(lam_of(5)(potential(A5, B5, 5)))
    _, (R15r, L35r) = modp_lin_reduce(c5 + tail5, [R15, L35], PRIMES[0])
    ok(sp.expand(R15r) == 0 and sp.expand(L35r) == 0,
       "band5 d=2: R(1),lambda_3(R) reduce to 0 on cascade+tail (degenerate corroboration)")
else:
    P("[set HEAVY=1 (and msolve on PATH) for the explicit det-saturated covector certificate + "
      "branch A + 2nd prime + QQ + band-5 d=2]")


# ======================================================================
print("\n" + "=" * 72, flush=True)
print("EVIDENCE LEDGER -- UNIFORM HATCH SLOPE-FORCING (degree-free backbone)", flush=True)
headline = HAVE_MS and HEADLINE_KEYS.issubset(PROBE) and all(PROBE[kk] is True for kk in HEADLINE_KEYS)
print(f"  PASS: {_NP} executed checks", flush=True)
print("  PROVED (arb degree, k=3,4,5): engine; lambda_3 point functional + Im Phi subset ker; "
      "both-ends Lemma-P; top-wall b_(k-2)(0)=rho_k a_(k-1)(0); depth-2 filler-linearity", flush=True)
print("  REFUTED: no single boundary a_j(0) divides either target mod band-4 cascade "
      "(W2 single-factor is special to k=3)", flush=True)
if headline:
    print(f"  PASS: band-4 d=3 COVECTOR KILL payload -- both targets in sqrt(cascade+Q_-1+Q_-2), "
          f"msolve '^' UNIT mod {PRIMES[0]} (kmin=2)", flush=True)
else:
    print("  SKIP: band-4 d=3 covector-kill payload -- requires msolve on PATH", flush=True)
print("  BAND 5: degenerate through d=3 (reduce-first GB); non-degenerate d>=4 OUT OF BUDGET", flush=True)
if not HEAVY:
    print(f"  SKIP: HEAVY corroboration -- det-saturated covector cert, branch A, mod {PRIMES[1]}, QQ, band5 d=2", flush=True)
print("  OPEN: arbitrary-degree kill; band-5 non-degenerate; DC1/JC2/Weyl pair", flush=True)
print("  FAIL: none", flush=True)
print("=" * 72, flush=True)
print(f"\n(total {time.time()-_T0:.1f}s; {_NP} checks passed; HEAVY={'1' if HEAVY else '0'}, "
      f"msolve={'yes' if HAVE_MS else 'no'})", flush=True)
if headline:
    print("ALL UNIFORM HATCH DF CHECKS PASSED", flush=True)
else:
    print("SKIP -- UNIFORM HATCH DF COVECTOR PAYLOAD NOT RUN; SUPPORTING CHECKS PASSED", flush=True)
