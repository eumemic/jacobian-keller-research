#!/usr/bin/env python3
"""verify_nonpoint_covector.py -- exact certificate for nonpoint-covector.md.

TARGET (architecture step 2): the "degree-free non-point covector" lambda_np for
the band-k hatches, k>=4.  The hatch (unique common root translated onto the
anchor 0):
    band 4:  a_4 = E(E+3)(E+6)(E+9),   u = b_3 = E(E+4)(E+8),   b_4=0.
    band 5:  a_5 = E(E+4)(E+8)(E+12)(E+16), u=b_4=E(E+5)(E+10)(E+15), b_5=0.
    general k: a_k = prod_{i<k}(E+i(k-1)),  u = prod_{j<k-1}(E+jk).

HEADLINE CORRECTION TO hatch-census.md / band45-lambda.md:
  The hatch's THIRD annihilator of Im Phi is NOT non-point.  It is a POINT
  functional lambda_3 (finite support).  The census's "one NON-POINT (infinite
  support) annihilator" is a computational ARTIFACT of its point-annihilator
  search window `range(-3, 3k+5)`, which starts at -3 and MISSES the genuine
  finite support at large-negative positions (down to -(k^2-3k+1)).  Widening the
  window recovers a third POINT annihilator; the honest codim is still 3, so
  {ev_0, ev_1, lambda_3} EXHAUST Ann(Im Phi) -- there is NO non-point covector.

PROOF (degree-free, arbitrary k): every annihilator's symbol L(t)=sum c_x t^x has
S_k*L = P_K (finite, supported on roots q_K) and S_{k-1}*L = P_H (finite, supported
on roots q_H); since gcd(S_k, S_{k-1}) = 1 in F[t,t^{-1}], S_k | P_K, so
L = P_K/S_k is a Laurent POLYNOMIAL => finite support => lambda is a POINT
functional.  There is no room in coker Phi (codim 3) for an infinite-support one.

Conventions frozen from the band-k corpus (b7e85e8..e223d21):
  A_1[x^-1]=(+)_k x^k C[E], (x^a f)(x^b g)=x^{a+b} f(E+b)g(E), f^[n](E)=f(E+n),
  Q_m=sum_{i+l=m}(b_l^[i] a_i - a_i^[l] b_i), gauge b_k=0, u:=b_{k-1},
  membership (E)_j=E(E-1)...(E-j+1) | a_-j,b_-j, Q_0=(T-1)G, G(0)=0.

SECTIONS
  0. Engine (RE-DERIVED): Q_m=[D,X]_m, Q_0=(T-1)G, G(0)=0 (k=4,5).
  1. lambda_3 is a POINT functional: moving-sum symbol; DIVISIBILITY theorem
     (all annihilators point, arbitrary k); explicit closed forms k=3..7;
     lambda_3(E)=-(k-1)(k-2)/2, lambda_3(1)=1; honest left-nullspace codim=3
     = dim span{ev_0,ev_1,lambda_3}; the census narrow-window artifact demo.
  2. Im Phi subset ker(lambda_3) DEGREE-FREE: the two support identities +
     symbolic block annihilation in generic C,V (arbitrary degree), k=3..7.
  3. Residual lambda_3(R) on the cascade (bounded): the reduction
     Q_0=1 <=> {R(0)=0(auto), R(1)=1, lambda_3(R)=lambda_3(E)}; lambda_3(R)=0
     forced at d<=2 (both bands, both branches) => lambda_3(E-R)=lambda_3(E)!=0
     (single-functional kill, d<=2).  [HEAVY] at d=3 (band 4) lambda_3(R) is a
     FREE modulus, so the identity FAILS at d=3 and the kill becomes the JOINT
     {slope, lambda_3} infeasibility (Q_0=1 UNIT, reproducing the census).
  4. General k closed form; degree-free annihilation k=6,7; residual conjecture.
  5. Sanity: lambda_3(E) != 0 (all k>=3).

Run:  uv run --with sympy python research/dc1-program/verify_nonpoint_covector.py
Heavy legs (d=3 msolve/Groebner):  HEAVY=1 uv run ... .
Ends: ALL NONPOINT COVECTOR CHECKS PASSED
"""
import sympy as sp
import os, time, subprocess, tempfile, re, shutil

E, t = sp.symbols("E t")
_T0 = time.time()
HEAVY = os.environ.get("HEAVY", "") not in ("", "0")


def sh(f, n):
    return sp.expand(sp.sympify(f).subs(E, E + n))


def falling(n):
    return sp.prod(E - j for j in range(n)) if n > 0 else sp.Integer(1)


def poly(name, deg):
    cs = list(sp.symbols(f"{name}_0:{deg+1}"))
    return sp.expand(sum(cs[j] * E**j for j in range(deg + 1))), cs


def ok(cond, label):
    if not cond:
        raise AssertionError("FAIL " + label)
    print(f"PASS [{time.time()-_T0:6.1f}s] {label}", flush=True)


def okz(v, label):
    if sp.expand(v) != 0:
        raise AssertionError(f"FAIL {label}: residual {sp.factor(sp.expand(v))}")
    print(f"PASS [{time.time()-_T0:6.1f}s] {label}", flush=True)


# ------- ladder / potential (RE-DERIVED, not imported) -------
def q_m(A, B, m, k):
    L = range(-k, k + 1)
    return sp.expand(sum(sh(B[l], i) * A[i] - sh(A[i], l) * B[l]
                         for i in L for l in L if i + l == m))


def potential(A, B, k):
    return sp.expand(sum(sh(A[i], j - i) * sh(B[-i], j) - sh(B[i], j - i) * sh(A[-i], j)
                         for i in range(1, k + 1) for j in range(i)))


def mul_ladders(P, Q):
    R = {}
    for k1, p1 in P.items():
        for k2, p2 in Q.items():
            R[k1 + k2] = R.get(k1 + k2, 0) + sh(p1, k2) * p2
    return {k: sp.expand(v) for k, v in R.items()}


def commutator_coeff(X, D, m):
    return sp.expand(mul_ladders(D, X).get(m, 0) - mul_ladders(X, D).get(m, 0))


def hatch(k):
    a = sp.expand(sp.prod(E + i * (k - 1) for i in range(k)))
    u = sp.expand(sp.prod(E + j * k for j in range(k - 1)))
    return a, u


def Kblock(a, c, k):
    return sp.expand(sum(sh(a, j - k) * sh(c, j) for j in range(k)))


def Hblock(u, v, k):
    return sp.expand(sum(sh(u, j - (k - 1)) * sh(v, j) for j in range(k - 1)))


# ======================================================================
print("=== 0. ENGINE (re-derived): Q_m=[D,X]_m, Q_0=(T-1)G, G(0)=0 (k=4,5) ===")
# ======================================================================
for k in (4, 5):
    Ag = {lev: sp.expand((falling(-lev) if lev < 0 else 1) * poly(f"Ag{lev+k}", 2)[0])
          for lev in range(-k, k + 1)}
    Bg = {lev: sp.expand((falling(-lev) if lev < 0 else 1) * poly(f"Bg{lev+k}", 2)[0])
          for lev in range(-k, k + 1)}
    for m in range(-2 * k, 2 * k + 1):
        okz(q_m(Ag, Bg, m, k) - commutator_coeff(Ag, Bg, m),
            f"k={k}: Q_{m} = [D,X]_{m} (generic deg-2 coeffs)")
    okz(q_m(Ag, Bg, 0, k) - (sh(potential(Ag, Bg, k), 1) - potential(Ag, Bg, k)),
        f"k={k}: Q_0 = (T-1)G")
    okz(potential(Ag, Bg, k).subs(E, 0), f"k={k}: G(0)=0 under membership")


# ======================================================================
print("\n=== 1. lambda_3 is a POINT functional (census 'non-point' REFUTED) ===")
# ======================================================================
def qK_roots(a, k):
    return sorted(set([int(r) for r in sp.solve(sh(a, -k), E)] + list(range(k))))


def qH_roots(u, k):
    return sorted(set([int(r) for r in sp.solve(sh(u, -(k - 1)), E)] + list(range(k - 1))))


def moving_sum_solutions(k):
    """Solve S_k L = P_K (supp roots qK), S_{k-1} L = P_H (supp roots qH).
    Returns (RK, RH, dim, [(support, weights, L)] basis)."""
    a, u = hatch(k)
    RK, RH = qK_roots(a, k), qH_roots(u, k)
    Sk = sum(t**j for j in range(k))
    Skm1 = sum(t**j for j in range(k - 1))
    pk = {x: sp.Symbol(f"p{x+50}") for x in RK}
    ph = {x: sp.Symbol(f"q{x+50}") for x in RH}
    PK = sum(pk[x] * t**x for x in RK)
    PH = sum(ph[x] * t**x for x in RH)
    mn = min(RK + RH)
    dd = sp.expand((PK * Skm1 - PH * Sk) * t**(-mn + 2))
    pdd = sp.Poly(dd, t)
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
        # den must be a monomial t^d  (so L is a Laurent polynomial = finite support)
        is_laurent = (Ld.is_monomial and len(Ld.free_symbols) <= 1)
        wt = {i - Ld.degree(): int(Ln.nth(i)) for i in range(Ln.degree() + 1) if Ln.nth(i) != 0}
        basis.append((sorted(wt), wt, Lb, is_laurent))
    return RK, RH, len(freev), basis


LAM3 = {}   # canonical shortest-support representative per k
for k in range(3, 8):
    RK, RH, dim, basis = moving_sum_solutions(k)
    ok(dim == 3, f"k={k}: moving-sum annihilator space is 3-dimensional (= codim Im Phi)")
    ok(all(b[3] for b in basis),
       f"k={k}: EVERY annihilator symbol L=P_K/S_k is a LAURENT POLYNOMIAL "
       f"(finite support) => ALL annihilators are POINT (gcd(S_k,S_(k-1))=1)")
    # the shortest-support representative (not ev_0=t^0 nor ev_1=t^1):
    cand = [b for b in basis if b[0] not in ([0], [1])]
    cand.sort(key=lambda b: len(b[0]))
    _, wt, Lb, _ = cand[0]
    LAM3[k] = wt
    lamE = sum(c * p for p, c in wt.items())
    lam1 = sum(wt.values())
    ok(lamE == -sp.Rational((k - 1) * (k - 2), 2),
       f"k={k}: lambda_3(E) = {lamE} = -(k-1)(k-2)/2  (closed form)")
    ok(lam1 == 1, f"k={k}: lambda_3(1) = {lam1}")
    print(f"        k={k}: lambda_3 = { {p: wt[p] for p in sorted(wt)} }  (symbol L={Lb})", flush=True)

# gcd(S_k,S_{k-1})=1 -- the engine of the divisibility theorem
for k in range(3, 9):
    Sk = sum(t**j for j in range(k))
    Skm1 = sum(t**j for j in range(k - 1))
    okz(sp.gcd(sp.Poly(Sk, t), sp.Poly(Skm1, t)).as_expr() - 1,
        f"k={k}: gcd(S_k,S_(k-1))=1  (forces L Laurent => point annihilators only)")


# honest left-nullspace codim, and the census narrow-window ARTIFACT
def images_deg_le(a, u, k, N):
    ims = []
    for j in range(N + 3):
        g = Kblock(a, falling(k) * E**j, k)
        if sp.Poly(g, E).degree() <= N:
            ims.append(g)
        h = Hblock(u, falling(k - 1) * E**j, k)
        if sp.Poly(h, E).degree() <= N:
            ims.append(h)
    return ims


def codim_leftnull(a, u, k, N):
    ims = images_deg_le(a, u, k, N)
    M = sp.Matrix([[sp.Poly(g, E).coeff_monomial(E**i) for i in range(N + 1)] for g in ims])
    return (N + 1) - M.rank()


def point_annih_dim(a, u, k, supp, maxj):
    ims = []
    for j in range(maxj):
        ims += [Kblock(a, falling(k) * E**j, k), Hblock(u, falling(k - 1) * E**j, k)]
    cs = sp.symbols(f"c0:{len(supp)}")
    eqs = [sp.expand(sum(cs[i] * g.subs(E, supp[i]) for i in range(len(supp)))) for g in ims]
    M, _ = sp.linear_eq_to_matrix(eqs, cs)
    return len(M.nullspace())


for k in (4, 5):
    a, u = hatch(k)
    cods = [codim_leftnull(a, u, k, N) for N in (16, 18, 20)]
    ok(cods == [3, 3, 3], f"k={k}: honest codim Im Phi = 3 (stable N=16,18,20)")
    # lambda_3, ev_0, ev_1 independent on the cokernel basis {1,E,E^2}
    wt = LAM3[k]

    def lam(f):
        return sp.expand(sum(c * sp.sympify(f).subs(E, p) for p, c in wt.items()))
    Mdual = sp.Matrix([[1, 1, 1],
                       [0, 1, lam(E)],
                       [0, 1, lam(E**2)]])
    ok(Mdual.det() != 0,
       f"k={k}: {{ev_0, ev_1, lambda_3}} independent on {{1,E,E^2}} (det={Mdual.det()}) "
       f"=> they SPAN coker Phi (dim 3) => EXHAUST Ann(Im Phi): no non-point covector")
    # the artifact: census window range(-3,3k+5) undercounts; wide window recovers 3
    narrow = point_annih_dim(a, u, k, list(range(-3, 3 * k + 5)), 3 * k + 12)
    lo = -(k * k - 3 * k + 1) - 2
    wide = point_annih_dim(a, u, k, list(range(lo, 3 * k + 5)), 3 * k + 12)
    ok(narrow == 2 and wide == 3,
       f"k={k}: census narrow window [-3,..] finds {narrow} POINT annihs (=> falsely infers "
       f"1 non-point); wide window [{lo},..] finds {wide} POINT annihs => the 3rd is POINT, "
       f"not non-point (census ARTIFACT)")


# ======================================================================
print("\n=== 2. Im Phi subset ker(lambda_3) DEGREE-FREE (adjoint support identities) ===")
# ======================================================================
for k in range(3, 8):
    a, u = hatch(k)
    wt = LAM3[k]
    Sk = sum(t**j for j in range(k))
    Skm1 = sum(t**j for j in range(k - 1))
    # symbol L (as Laurent poly) from the weights
    L = sum(c * t**p for p, c in wt.items())
    PK = sp.expand(Sk * L)
    PH = sp.expand(Skm1 * L)
    RK, RH = qK_roots(a, k), qH_roots(u, k)
    suppPK = [i for i in range(sp.Poly(sp.expand(PK * t**30), t).degree() + 1)
              if sp.Poly(sp.expand(PK * t**30), t).nth(i) != 0]
    suppPK = [i - 30 for i in suppPK]
    suppPH = [i for i in range(sp.Poly(sp.expand(PH * t**30), t).degree() + 1)
              if sp.Poly(sp.expand(PH * t**30), t).nth(i) != 0]
    suppPH = [i - 30 for i in suppPH]
    ok(set(suppPK) <= set(RK),
       f"k={k}: support(S_k * L) = {sorted(suppPK)} subset roots(q_K)  (K-block support identity)")
    ok(set(suppPH) <= set(RH),
       f"k={k}: support(S_(k-1) * L) = {sorted(suppPH)} subset roots(q_H)  (H-block support identity)")

    # degree-free block annihilation: symbolic in generic degree-(k+1) C, V
    def lam(f):
        return sp.expand(sum(c * sp.sympify(f).subs(E, p) for p, c in wt.items()))
    Cc = list(sp.symbols(f"Cc{k}_0:{k+2}"))
    Vv = list(sp.symbols(f"Vv{k}_0:{k+2}"))
    C = sum(Cc[i] * E**i for i in range(k + 1))
    V = sum(Vv[i] * E**i for i in range(k + 1))
    okz(lam(Kblock(a, falling(k) * C, k)),
        f"k={k}: lambda_3(K_k[(E)_k C]) = 0 for GENERIC deg-{k} C (degree-free)")
    okz(lam(Hblock(u, falling(k - 1) * V, k)),
        f"k={k}: lambda_3(H_(k-1)[(E)_(k-1) V]) = 0 for GENERIC deg-{k} V (degree-free)")
    # the structural identity K_k[(E)_k C] = sum_{j=0}^{k-1} (q_K C)(E+j)
    qK = sp.expand(sh(a, -k) * falling(k))
    okz(Kblock(a, falling(k) * C, k) - sp.expand(sum(sh(qK * C, j) for j in range(k))),
        f"k={k}: K_k[(E)_k C] = sum_j (q_K C)(E+j)  (adjoint identity, q_K=a(E-k)(E)_k)")


# ======================================================================
print("\n=== 3. Residual lambda_3(R): reduction + bounded d<=2 kill; d=3 boundary ===")
# ======================================================================
# --- authoritative cascade builder (census build_branch, RE-INCLUDED) ---
def clean_solve(A, B, m, lkey, name, membership, raw_degree, k):
    raw, cs = poly(name, raw_degree)
    unknown = sp.expand(falling(membership) * raw)
    trial = dict(B); trial[lkey] = unknown
    M, rhs = sp.linear_eq_to_matrix(sp.Poly(q_m(A, trial, m, k), E).all_coeffs(), cs)
    if any(entry.free_symbols for entry in M):
        raise AssertionError("operator matrix not numeric")
    conds = [c for c in (sp.expand(n.dot(rhs)) for n in M.T.nullspace()) if c != 0]
    indep = sp.zeros(0, len(cs)); selrhs = []
    for i in range(M.rows):
        cand = indep.col_join(M[i, :])
        if cand.rank() > indep.rank():
            indep = cand; selrhs.append(rhs[i])
    vals = ([sp.Integer(0)] * len(cs) if indep.rows == 0
            else [x.subs({s: 0 for s in x.free_symbols if str(s).startswith("tau")})
                  for x in indep.gauss_jordan_solve(sp.Matrix(selrhs))[0]])
    res = sp.expand(unknown.subs(dict(zip(cs, vals))))
    kernels = []
    for j, vec in enumerate(M.nullspace()):
        p = sp.symbols(f"{name}K{j}")
        res = sp.expand(res + p * falling(membership) * sum(vec[i] * E**i for i in range(len(cs))))
        kernels.append(p)
    return res, kernels, conds


def build_branch(ak, u, d, k, branch):
    A = {k: ak}; avars = []
    for i in range(k - 1, -1, -1):
        p, c = poly(f"a{i}", d); A[i] = p; avars += c
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
    else:
        bmk, c = poly("bmk", d); A[-k] = sp.Integer(0); B[-k] = sp.expand(falling(k) * bmk); extra = c
    return A, B, [e for e in conds if sp.expand(e) != 0], avars + extra + kernels


def scoeffs(e):
    return [c for c in sp.Poly(sp.expand(e), E).all_coeffs() if sp.expand(c) != 0]


def lin_elim(eqs, targets):
    eqs = [sp.expand(e) for e in eqs if sp.expand(e) != 0]; ch = True
    while ch:
        ch = False
        for e in eqs:
            for v in list(e.free_symbols):
                pe = sp.Poly(e, v)
                if pe.degree() == 1 and not pe.nth(1).free_symbols and pe.nth(1) != 0:
                    sol = sp.expand(-pe.nth(0) / pe.nth(1))
                    eqs = [sp.expand(x.subs(v, sol)) for x in eqs]; eqs = [x for x in eqs if x != 0]
                    targets = [sp.expand(x.subs(v, sol)) for x in targets]; ch = True; break
            if ch:
                break
    return eqs, targets


HATCH = {4: (sp.expand(E * (E + 3) * (E + 6) * (E + 9)), sp.expand(E * (E + 4) * (E + 8))),
         5: (sp.expand(E * (E + 4) * (E + 8) * (E + 12) * (E + 16)),
             sp.expand(E * (E + 5) * (E + 10) * (E + 15)))}

# 3a. the reduction (structural): Q_0=1 <=> G=E <=> E-R in Im Phi <=> the 3 functionals
#     vanish on E-R.  Justified by section-1 basis/exhaustion; here reproduce the
#     census UNIT at d=2 both via coeffs(Q_0-1) and via the 3-condition form.
for k in (4, 5):
    ak, u = HATCH[k]; wt = LAM3[k]

    def lam(f):
        return sp.expand(sum(c * sp.sympify(f).subs(E, p) for p, c in wt.items()))
    lamE = lam(E)
    for d in (1, 2):
        for br in ("A", "B"):
            A, B, pos, fv = build_branch(ak, u, d, k, br)
            G = potential(A, B, k)
            lamG = lam(G); R0 = sp.expand(G.subs(E, 0))
            _, (lamGr, R0r) = lin_elim(pos, [lamG, R0])
            ok(R0r == 0 and lamGr == 0,
               f"k={k} d={d} br{br}: R(0)=0 and lambda_3(R)=lambda_3(G)=0 mod cascade "
               f"=> lambda_3(E-R)=lambda_3(E)={lamE} != 0 (single-functional kill at d<=2)")

# reproduce the census moment-unit UNIT at d=2 (sympy Groebner), and show it equals
# the 3-condition form -- both are the UNIT ideal (slope R(1) is pinned to 0 at d<=2).
def is_unit(eqs, vs):
    eqs = [sp.expand(sp.numer(sp.together(e))) for e in eqs if sp.expand(e) != 0]
    return list(sp.groebner(eqs, *vs, order="grevlex").exprs) == [sp.Integer(1)]


for k in (4, 5):
    ak, u = HATCH[k]; wt = LAM3[k]

    def lam(f):
        return sp.expand(sum(c * sp.sympify(f).subs(E, p) for p, c in wt.items()))
    A, B, pos, fv = build_branch(ak, u, 2, k, "A")
    G = potential(A, B, k); R1 = sp.expand(G.subs(E, 1))
    ok(is_unit(pos + scoeffs(q_m(A, B, 0, k) - 1), fv),
       f"k={k} d=2 brA: cascade + Q_0=1 = UNIT (reproduces hatch-census moment-unit kill) [sympy]")
    ok(is_unit(pos + [R1 - 1], fv),
       f"k={k} d=2 brA: cascade + R(1)=1 = UNIT (slope pinned to 0 at d<=2) [sympy]")

# 3c. HEAVY: d=3 band 4 -- the residual identity FAILS and the kill goes joint.
if HEAVY:
    print("    [HEAVY] d=3 band-4 boundary (msolve + mod-p Groebner)")
    HAVE_MS = shutil.which("msolve") is not None

    def cd(e):
        n, _ = sp.fraction(sp.together(sp.expand(e))); return sp.expand(n)

    def clr(e):
        n = cd(e)
        if not n.free_symbols:
            return n
        p = sp.Poly(n, *sorted(n.free_symbols, key=str))
        L = sp.ilcm(*[int(sp.denom(sp.nsimplify(c))) for c in p.coeffs()]) if p.coeffs() else 1
        return sp.expand(n * L)

    def msolve_unit(eqs, fv, char=65003):
        eqs = [cd(e) for e in eqs if sp.expand(e) != 0]
        xs = sp.symbols(f"z_0:{len(fv)}"); sub = dict(zip(fv, xs))
        with tempfile.NamedTemporaryFile("w", suffix=".ms", delete=False) as f:
            inp = f.name; f.write(",".join(str(x) for x in xs) + "\n" + str(char) + "\n")
            polys = [str(cd(e.subs(sub))).replace(" ", "").replace("**", "^") for e in eqs]
            f.write(",\n".join(polys) + "\n")
        out = inp + ".out"
        try:
            subprocess.run(["msolve", "-g", "2", "-f", inp, "-o", out],
                           check=True, stderr=subprocess.DEVNULL, timeout=280)
            r = open(out).read()
        finally:
            for p_ in (inp, out):
                if os.path.exists(p_):
                    os.remove(p_)
        m = re.search(r"length of basis:\s*(\d+)", r)
        if m and m.group(1) == "1" and "[1]" in r.split("---")[-1]:
            return "UNIT"
        return "FEASIBLE" if m else ("UNIT" if r.strip().startswith("[-1]") else "FEASIBLE")

    ak, u = HATCH[4]; wt = LAM3[4]

    def lam(f):
        return sp.expand(sum(c * sp.sympify(f).subs(E, p) for p, c in wt.items()))
    lamE = lam(E)
    A, B, pos, fv = build_branch(ak, u, 3, 4, "A")
    G = potential(A, B, 4); lamG = sp.expand(lam(G)); R1 = sp.expand(G.subs(E, 1))
    R0 = sp.expand(G.subs(E, 0))
    # (i) census Q_0=1 UNIT reproduced two ways
    if HAVE_MS:
        ok(msolve_unit(pos + scoeffs(q_m(A, B, 0, 4) - 1), fv) == "UNIT",
           "k=4 d=3 brA: cascade + Q_0=1 = UNIT (census kill) [msolve ^]")
        ok(msolve_unit(pos + [R0, R1 - 1, sp.expand(lamG - lamE)], fv) == "UNIT",
           "k=4 d=3 brA: cascade + R(0)=0 + R(1)=1 + lambda_3(R)=lambda_3(E) = UNIT "
           "(SAME kill via the 3-covector reduction) [msolve ^]")
    # (ii) but lambda_3(R) is now a FREE modulus (residual identity FAILS at d=3)
    eqs2, (lamGr,) = lin_elim(pos, [lamG])
    rv = sorted(set().union(*([e.free_symbols for e in eqs2] + [lamGr.free_symbols])), key=str)
    gb = sp.groebner([clr(e) for e in eqs2 if sp.expand(e) != 0], *rv,
                     order="grevlex", modulus=32003)
    red = gb.reduce(clr(lamGr))[1]
    ok(sp.sympify(red).free_symbols != set(),
       "k=4 d=3: lambda_3(R) is a FREE (non-constant) modulus mod cascade "
       "=> the residual identity lambda_3(R)=0 FAILS at d=3; the kill is the JOINT "
       "{R(1)=1, lambda_3(R)=lambda_3(E)} infeasibility, not a single covector")
else:
    print("    [set HEAVY=1 for the d=3 band-4 boundary: Q_0=1 UNIT two ways + lambda_3(R) free]")


# ======================================================================
print("\n=== 4. General k: closed form + degree-free annihilation (k=6,7) ===")
# ======================================================================
# (annihilation for k=6,7 already checked in section 2 loop; restate closed form)
for k in (6, 7):
    wt = LAM3[k]
    lamE = sum(c * p for p, c in wt.items())
    ok(lamE == -sp.Rational((k - 1) * (k - 2), 2) and lamE != 0,
       f"k={k}: lambda_3 kills Im Phi degree-free (sec 2) and lambda_3(E)={lamE}=-(k-1)(k-2)/2 != 0")
print("        general-k RESIDUAL CONJECTURE: lambda_3(R)=0 forced on the cascade at low")
print("        filler degree (d<=2 verified, bands 4,5) but acquires the freed cascade")
print("        moduli at higher d, so the single-covector kill is BOUNDED and the")
print("        arbitrary-degree kill is the joint {slope, lambda_3} obstruction.")


# ======================================================================
print("\n=== 5. Sanity: lambda_3(E) != 0 for all bands ===")
# ======================================================================
for k in range(3, 8):
    wt = LAM3[k]
    lamE = sum(c * p for p, c in wt.items())
    ok(lamE != 0, f"k={k}: lambda_3(E) = {lamE} != 0 (obstruction is live; slope-achievability alone cannot rescue Q_0=1)")


print("\n" + "=" * 70)
print("SUMMARY")
print("  * The band-k hatch's 3rd annihilator lambda_3 is a POINT functional")
print("    (finite support), NOT non-point: census's 'non-point covector' is a")
print("    search-window ARTIFACT. gcd(S_k,S_(k-1))=1 => ALL annihilators point.")
print("  * Explicit closed form + lambda_3(E)=-(k-1)(k-2)/2, lambda_3(1)=1 (k=3..7).")
print("  * Im Phi subset ker(lambda_3) DEGREE-FREE (adjoint support identities).")
print("  * Residual: lambda_3(R)=0 => lambda_3(E-R)=lambda_3(E)!=0 at d<=2 (single-")
print("    functional kill); at d=3 lambda_3(R) is a free modulus, the identity FAILS,")
print("    and the kill is the JOINT {slope, lambda_3} infeasibility (Q_0=1 UNIT).")
print("=" * 70)
print("\nALL NONPOINT COVECTOR CHECKS PASSED")
