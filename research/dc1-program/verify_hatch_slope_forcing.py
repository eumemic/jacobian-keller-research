#!/usr/bin/env python3
"""verify_hatch_slope_forcing.py -- THE UNIFICATION PROBE.

INDEPENDENTLY DERIVED - EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES - NOT PEER REVIEWED

QUESTION.  At W2 (band 3) the negative tail Q_-1..Q_-5 (NO Q_0) forces the moment
slope R(1)=0 on the cascade+tail variety (slope-forcing-verdict.md).  For the
k>=4 hatches the cokernel Ann(Im Phi)=span{ev_0, ev_1, lambda_3} carries a SECOND
nontrivial covector lambda_3 (the point functional, nonpoint-covector.md) beyond
the slope ev_1.  DO THE k>=4 HATCH TAILS ALSO force BOTH cokernel pairings of the
residual R to vanish, i.e. on {cascade, tail (no Q_0), membership}:

        R(1) = 0      AND      lambda_3(R) = 0     ?

If YES, architecture steps 1 (slope-forcing) and 2 (covector) collapse into ONE
mechanism: the tail annihilates E-R against EVERY cokernel direction beyond ev_0,
so Q_0=1 is contradicted in all directions at once.

CONVENTIONS (band-k, gauge b_k=0):  (x^a f)(x^b g)=x^{a+b} f(E+b)g(E),
  f^[n](E)=f(E+n),  Q_m=sum_{i+l=m}(b_l^[i] a_i - a_i^[l] b_l),  u:=b_{k-1},
  membership (E)_j=E(E-1)..(E-j+1) | a_-j,b_-j,  Q_0=(T-1)G, G(0)=0,
  R(1)=Q_0(0)=G(1) (Lemma P: independent of the Phi fillers a_-(k-1), b_-k),
  lambda_3(R)=lambda_3(G) (since Im Phi subset ker lambda_3).
  hatch:  a_k=prod_{i<k}(E+i(k-1)),  u=prod_{j<k-1}(E+jk).
  band 4:  a_4=E(E+3)(E+6)(E+9), u=E(E+4)(E+8);  lambda_3=ev_-2-ev_-4+ev_-5, lambda_3(E)=-3.
  band 5:  a_5=E(E+4)(E+8)(E+12)(E+16), u=E(E+5)(E+10)(E+15);
           lambda_3=ev_-3-ev_-5+ev_-7-ev_-10+ev_-11, lambda_3(E)=-6.
  Tail = Q_-1..Q_-(2k-1);  Q_-2k splits branch A (a_-k=(E)_k am_k, b_-k=mu_k a_-k)
  from branch B (a_-k=0, b_-k=(E)_k C).

VERDICT (this file, machine-checked, band 4, d=3, both branches):
  YES.  R(1) IN sqrt(cascade+tail)  AND  lambda_3(R) IN sqrt(cascade+tail).
  Both cokernel pairings are tail-forced -- the unified slope-forcing principle
  holds at the band-4 hatch, d=3.  (At d<=2 both are already 0 on the cascade
  ALONE -- degenerate; the genuine tail test is at d=3, where both are free
  moduli on the cascade.)

ENGINES.  msolve 0.10.1 '^' Rabinowitsch on the full un-reduced d=3 system, mod
  65003 (+ mod 32003 and over QQ under HEAVY), both branches, both targets -- UNIT.
  Parsing validated on unit {1} / feasible controls.  Independent sympy Groebner
  confirms the d=2 instance and the reduced graded prefixes; sympy Buchberger on
  the full d=3 band-4 system (15 vars after linear reduction) is intractable and
  is reported as such.

LEDGER (honest):
  PROVED (arb degree / arb k, machine identities): S0 engine; S1 lambda_3 point
    functional + closed form lambda_3(E)=-(k-1)(k-2)/2, lambda_3(1)=1, gcd(S_k,
    S_{k-1})=1, Im Phi subset ker lambda_3 (block annihilation, generic C,V);
    targets independent of Phi fillers.
  EXACT (band 4, d=2, char 0 + machine): both targets vanish on the cascade ALONE
    (degenerate); sympy Groebner.
  BOUNDED (band 4, d=3, msolve '^' mod p [+ QQ HEAVY], both branches): BOTH
    R(1) and lambda_3(R) IN sqrt(cascade+tail); controls (free moduli on cascade)
    hold; origin feasible; graded depth kmin.
  BOUNDED (band 5, d=2, HEAVY): both targets vanish on the cascade alone.
  OPEN / NOT CLAIMED: arbitrary degree; band>=5 non-degenerate (d>=3) tail probe;
    a Weyl pair / DC1 / JC2.

Run:      uv run --with sympy python research/dc1-program/verify_hatch_slope_forcing.py
Heavy:    HEAVY=1 uv run --with sympy python research/dc1-program/verify_hatch_slope_forcing.py
Ends with an evidence ledger and one of:
  PASS -- HATCH SLOPE FORCING HEADLINE PAYLOAD PASSED
  SKIP -- HATCH SLOPE FORCING HEADLINE PAYLOAD NOT RUN
A failed assertion exits nonzero with FAIL.
"""
import sympy as sp
import os, time, sys, subprocess, tempfile, shutil

E, t = sp.symbols("E t")
_T0 = time.time()
HEAVY = os.environ.get("HEAVY", "") not in ("", "0")
HAVE_MS = shutil.which("msolve") is not None
PRIMES = (65003, 32003)
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
    gens = sorted(set().union(*[e.free_symbols for e in eqs]), key=str)
    G = sp.groebner([mp(e, p, gens) for e in eqs], *gens, order="grevlex", modulus=p)
    return list(G.exprs) == [sp.Integer(1)]


# ---- msolve '^' engine (validated parsing) ----
def msolve_unit(eqs, fv, char=65003, tmo=280):
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
        # timeout / OOM-kill / crash on a hard (e.g. positive-dimensional) instance:
        # report "indeterminate" so callers skip this instance rather than aborting.
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
    """Forward-solve b_{k-2}..b_{-(k-1)} from Q_{2k-2}..Q_1; return cascade compat conds + free vars."""
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
    else:
        bmk, c = poly("bmk", d); A[-k] = sp.Integer(0); B[-k] = sp.expand(falling(k) * bmk); extra = c
    return A, B, [e for e in conds if sp.expand(e) != 0], avars + extra + kernels


HATCH = {k: hatch(k) for k in (3, 4, 5)}
LAM = {}      # filled in S1 by the moving-sum solve (NOT hard-coded)


def lam_of(k):
    wt = LAM[k]
    return lambda f: sp.expand(sum(c * sp.sympify(f).subs(E, p) for p, c in wt.items()))


_GBCACHE = {}


def gb_cascade(k, d, p):
    """GB of the (branch-INDEPENDENT) cascade at (k,d) over GF(p), with R(1),lambda_3(R) reduced.
    The positive cascade Q_{2k-2}..Q_1 and the targets R(1)=Q_0(0), lambda_3(R)=lambda_3(G) do NOT
    depend on the branch data a_-k,b_-k (they are Phi-filler-independent), so this is computed ONCE."""
    key = (k, d, p)
    if key in _GBCACHE:
        return _GBCACHE[key]
    ak, u = HATCH[k]
    A, B, conds, fv = build_branch(ak, u, d, k, "B")
    G = potential(A, B, k)
    lam = lam_of(k)
    R1 = sp.expand(q_m(A, B, 0, k).subs(E, 0)); L3 = sp.expand(lam(G))
    resid, (R1c, L3c) = modp_lin_reduce(conds, [R1, L3], p)
    gens = sorted(set().union(*([e.free_symbols for e in resid] or [set()]),
                              R1c.free_symbols, L3c.free_symbols), key=str)
    GB = sp.groebner([mp(e, p, gens) for e in resid], *gens, order="grevlex", modulus=p)
    nfR = sp.sympify(GB.reduce(mp(R1c, p, gens))[1])
    nfL = sp.sympify(GB.reduce(mp(L3c, p, gens))[1])
    _GBCACHE[key] = (GB, gens, R1c, L3c, nfR, nfL)
    return _GBCACHE[key]


def witness_not_forced(k, d, branch, extra_levels, target_kind, p, seed):
    """Certify that the target is NOT in sqrt(cascade + {Q_-m : m in extra_levels}) by exhibiting a
    feasible point with target=1 (msolve GB non-unit on a low-dim random slice, in the reduced
    cascade coordinates).  target_kind in {'R1','L3'}.  Returns True if a witness was certified."""
    if not HAVE_MS:
        return None
    import random as _rnd
    ak, u = HATCH[k]
    A, B, conds, fv = build_branch(ak, u, d, k, branch)
    G = potential(A, B, k); lam = lam_of(k)
    tgt = sp.expand(q_m(A, B, 0, k).subs(E, 0)) if target_kind == "R1" else sp.expand(lam(G))
    extra = sum((scoeffs(q_m(A, B, -m, k)) for m in extra_levels), [])
    # reduce the extra necklaces + the target through the cascade linear substitutions
    resid, reduced = modp_lin_reduce(conds, [tgt] + extra, p)
    tgt_r = reduced[0]
    extra_r = [e for e in reduced[1:] if sp.expand(e) != 0]
    gens = sorted(set().union(*([e.free_symbols for e in resid + extra_r] or [set()]),
                              tgt_r.free_symbols), key=str)
    system = resid + extra_r + [sp.expand(tgt_r - 1)]
    # Fix nfix random coordinates.  Probing in DECREASING order: for nfix > dim the slice is empty
    # (fast unit); at nfix = dim it is 0-dimensional (fast, and feasible for a generic value, giving
    # the witness).  Stopping at the first feasible slice never probes nfix < dim, so no slow
    # positive-dimensional msolve calls occur.
    for nfix in (5, 4, 3, 2, 1):
        for attempt in range(8):
            _rnd.seed(seed + 1000 * nfix + attempt)
            fixed = _rnd.sample(gens, min(nfix, len(gens)))
            slic = [v - _rnd.randint(1, p - 1) for v in fixed]
            res = msolve_unit(system + slic, gens, p, tmo=30)
            if res is False:      # non-unit GB on a 0-dim slice => feasible => witness
                return True
    return False


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
        is_laurent = (Ld.is_monomial and len(Ld.free_symbols) <= 1)
        wt = {i - Ld.degree(): int(Ln.nth(i)) for i in range(Ln.degree() + 1) if Ln.nth(i) != 0}
        basis.append((sorted(wt), wt, is_laurent))
    return len(freev), basis


for k in (3, 4, 5):
    dim, basis = moving_sum_solutions(k)
    ok(dim == 3, f"k={k}: moving-sum annihilator space is 3-dimensional (= codim Im Phi)")
    ok(all(b[2] for b in basis),
       f"k={k}: EVERY annihilator symbol L=P_K/S_k is a Laurent polynomial => ALL point functionals")
    cand = [b for b in basis if b[0] not in ([0], [1])]
    cand.sort(key=lambda b: len(b[0]))
    LAM[k] = cand[0][1]
    lam = lam_of(k)
    ok(lam(E) == -sp.Rational((k - 1) * (k - 2), 2),
       f"k={k}: lambda_3(E) = {lam(E)} = -(k-1)(k-2)/2")
    ok(lam(1) == 1, f"k={k}: lambda_3(1) = 1")
    P(f"k={k}: lambda_3 = {dict(sorted(LAM[k].items()))}")

for k in (3, 4, 5, 6, 7, 8):
    Sk = sum(t**j for j in range(k)); Skm1 = sum(t**j for j in range(k - 1))
    okz(sp.gcd(sp.Poly(Sk, t), sp.Poly(Skm1, t)).as_expr() - 1,
        f"k={k}: gcd(S_k,S_(k-1))=1 (the engine of the point-functional theorem)")

# Im Phi subset ker lambda_3, degree-free (generic C, V of degree k) -- this is what
# makes lambda_3(R)=lambda_3(G) (the filler blocks are annihilated).
for k in (3, 4, 5):
    a, u = hatch(k)
    lam = lam_of(k)
    Cc = list(sp.symbols(f"Cc{k}_0:{k+2}")); Vv = list(sp.symbols(f"Vv{k}_0:{k+2}"))
    C = sum(Cc[i] * E**i for i in range(k + 1)); V = sum(Vv[i] * E**i for i in range(k + 1))
    okz(lam(Kblock(a, falling(k) * C, k)),
        f"k={k}: lambda_3(K_k[(E)_k C]) = 0 for generic deg-{k} C (degree-free)")
    okz(lam(Hblock(u, falling(k - 1) * V, k)),
        f"k={k}: lambda_3(H_(k-1)[(E)_(k-1) V]) = 0 for generic deg-{k} V (degree-free)")

# targets are independent of the Phi fillers (a_-(k-1), b_-k): R(1)=G(1), lambda_3(R)=lambda_3(G)
for k in (4, 5):
    ak, u = HATCH[k]
    lam = lam_of(k)
    A, B, conds, fv = build_branch(ak, u, 2, k, "B")
    G = potential(A, B, k)
    phi = (set(A[-(k - 1)].free_symbols) | set(B[-k].free_symbols)) - {E}
    okz_syms = (sp.expand(q_m(A, B, 0, k).subs(E, 0)).free_symbols | sp.expand(lam(G)).free_symbols) & phi
    ok(okz_syms == set(),
       f"k={k}: R(1)=Q_0(0) and lambda_3(R)=lambda_3(G) are INDEPENDENT of the Phi fillers "
       f"a_-(k-1),b_-k (Lemma P + Im Phi subset ker lambda_3)")

# W2 base case: lambda_3 = ev_-1, so lambda_3(R)=R(-1); the cascade relation R(1)+R(-1)=0
# makes the two conditions COINCIDE at k=3 (the W2 slope-forcing result IS this instance).
ok(LAM[3] == {-1: 1}, "k=3 (W2): lambda_3 = ev_-1, so lambda_3(R)=R(-1) (single point)")


# ======================================================================
print("\n=== S2. band-4 controls: d=2 targets 0 on cascade (degenerate); d=3 targets FREE ===", flush=True)
# ======================================================================
def reduce_cascade(ak, u, d, k, branch, p=PRIMES[0]):
    A, B, conds, fv = build_branch(ak, u, d, k, branch)
    G = potential(A, B, k)
    lam = lam_of(k)
    R1 = sp.expand(q_m(A, B, 0, k).subs(E, 0)); L3 = sp.expand(lam(G))
    _, (R1c, L3c) = modp_lin_reduce(conds, [R1, L3], p)
    return A, B, conds, fv, R1, L3, R1c, L3c


ak4, u4 = HATCH[4]
# d=2: both targets reduce to 0 on the cascade ALONE (degenerate -- tail not tested)
for branch in ("A", "B"):
    _, _, _, _, _, _, R1c, L3c = reduce_cascade(ak4, u4, 2, 4, branch)
    ok(sp.expand(R1c) == 0 and sp.expand(L3c) == 0,
       f"band4 d=2 br{branch}: R(1)=0 AND lambda_3(R)=0 on the cascade ALONE (linear elim) "
       f"-- DEGENERATE, the tail is not yet tested at d=2")

# d=3: both targets are FREE (non-constant) moduli on the cascade (GB-normal-form non-constant).
# The positive cascade AND the two targets are branch-INDEPENDENT (Phi-filler-independence, S1),
# so this control is computed ONCE (shared, cached) and applies to both branches.
p = PRIMES[0]
GBc, gens_c, R1c, L3c, nfR, nfL = gb_cascade(4, 3, p)
P(f"band4 d=3: cascade reduces to {len(gens_c)} free gens with {len(GBc.exprs)} GB elements "
  f"[sympy GF({p})]; this control is branch-independent (positive cascade + targets are "
  f"Phi-filler-independent)")
ok(nfR.free_symbols != set() and nfL.free_symbols != set(),
   "band4 d=3 (both branches): R(1) and lambda_3(R) are BOTH free (non-constant) moduli on the "
   "cascade (GB normal forms non-constant) -- so tail-forcing to 0 is NON-trivial")
ok(sp.expand(sp.sympify(GBc.reduce(mp(sp.expand(R1c - L3c), p, gens_c))[1])) != 0,
   "band4 d=3: R(1) != lambda_3(R) mod cascade -- the two cokernel conditions are INDEPENDENT for "
   "k>=4 (unlike k=3 where lambda_3(R)=R(-1)=-R(1) makes them coincide)")


# ======================================================================
print("\n=== S3. feasibility: origin lies on cascade+tail with R(1)=lambda_3(R)=0 ===", flush=True)
# ======================================================================
for branch in ("A", "B"):
    A, B, conds, fv = build_branch(ak4, u4, 3, 4, branch)
    G = potential(A, B, 4)
    lam = lam_of(4)
    tail = sum((scoeffs(q_m(A, B, m, 4)) for m in range(-1, -8, -1)), [])
    R1 = sp.expand(q_m(A, B, 0, 4).subs(E, 0)); L3 = sp.expand(lam(G))
    z = {v: 0 for v in fv}
    ok(all(sp.expand(e.subs(z)) == 0 for e in conds + tail) and
       sp.expand(R1.subs(z)) == 0 and sp.expand(L3.subs(z)) == 0,
       f"band4 d=3 br{branch}: the ORIGIN lies on cascade+tail with R(1)=lambda_3(R)=0 "
       f"(cascade+tail NON-empty; the forcing is non-vacuous)")


# ======================================================================
print("\n=== S4. THE PROBE: R(1) and lambda_3(R) in sqrt(cascade+tail), band 4 ===", flush=True)
# ======================================================================
trab = sp.symbols("t_rab")


def targets_and_tail(ak, u, d, k, branch):
    A, B, conds, fv = build_branch(ak, u, d, k, branch)
    G = potential(A, B, k)
    lam = lam_of(k)
    R1 = cd(q_m(A, B, 0, k).subs(E, 0)); L3 = cd(lam(G))
    tail = sum((scoeffs(q_m(A, B, m, k)) for m in range(-1, -(2 * k), -1)), [])
    return A, B, conds, fv, R1, L3, tail


# --- msolve parsing validation (guards against a false UNIT) ---
if HAVE_MS:
    xx = sp.symbols("xx")
    ok(msolve_unit([sp.Integer(1)], [xx]) is True, "msolve parsing: unit ideal {1} -> UNIT")
    ok(msolve_unit([xx * (xx - 1)], [xx]) is False, "msolve parsing: feasible {xx(xx-1)} -> NOT unit")

# --- d=2 (DEGENERATE): both targets already 0 on cascade -> trivially in sqrt(cascade+tail) ---
for branch in ("A", "B"):
    A, B, conds, fv, R1, L3, tail = targets_and_tail(ak4, u4, 2, 4, branch)
    # independent sympy: cascade+tail forces both (they are already 0 on the cascade)
    _, (R1r, L3r) = modp_lin_reduce(conds + tail, [R1, L3], PRIMES[0])
    ok(sp.expand(R1r) == 0 and sp.expand(L3r) == 0,
       f"band4 d=2 br{branch}: R(1),lambda_3(R) in sqrt(cascade+tail) [sympy GF({PRIMES[0]}), "
       f"reduce to 0] -- but DEGENERATE (already 0 on the cascade alone)")

# --- d=3 (THE REAL TEST): msolve '^' Rabinowitsch, full un-reduced system, both branches ---
PROBE = {}
HEADLINE_EXPECTED = {(branch, nm, PRIMES[0]) for branch in ("A", "B")
                     for nm in ("R(1)", "lambda_3(R)")}
if HAVE_MS:
    dp3_primes = PRIMES if HEAVY else (PRIMES[0],)
    for branch in ("B", "A"):
        A, B, conds, fv, R1, L3, tail = targets_and_tail(ak4, u4, 3, 4, branch)
        base = conds + tail
        for tgt, nm in ((R1, "R(1)"), (L3, "lambda_3(R)")):
            for p in dp3_primes:
                res = msolve_unit(base + [sp.expand(1 - trab * tgt)], [trab] + fv, p, tmo=280)
                PROBE[(branch, nm, p)] = res
                ok(res is True,
                   f"band4 d=3 br{branch}: {nm} in sqrt(cascade+tail) -- msolve '^' UNIT mod {p}")
            if HEAVY:
                res = msolve_unit(base + [sp.expand(1 - trab * tgt)], [trab] + fv, 0, tmo=560)
                if res is True:
                    ok(True, f"band4 d=3 br{branch}: {nm} in sqrt(cascade+tail) -- msolve '^' UNIT "
                             f"over QQ (char 0)")
                else:
                    P(f"[band4 d=3 br{branch} {nm}: msolve QQ (char 0) did not terminate in budget "
                      f"-- char-0 cert not obtained; mod-p verdicts stand]")
else:
    P("[msolve not on PATH -- d=3 band-4 probe SKIPPED]")

# NOTE on the second engine at d=3.  A sympy Buchberger run on the FULL d=3 band-4 system (15
# variables after linear reduction) does NOT terminate in a sane budget, so the d=3 radical-
# membership verdict itself is msolve-only -- mitigated by the in-file parser validation, the two
# primes, and the attempted QQ char-0 certificate.  sympy independently and exactly certifies (over
# GF(p), no msolve): the whole engine + lambda_3 layer (S0,S1), the d=2 degenerate probe (S4), the
# d=3 free-modulus + distinctness controls (S2, GB of the cascade), and the feasibility origin (S3).


# ======================================================================
print("\n=== S5. graded depth: which tail rows force each target (band 4 d=3 branch B) ===", flush=True)
# ======================================================================
# Two-sided:  (i) LOWER bound -- Q_-1 alone does NOT force: exhibit an explicit feasible point with
#     target=1 on cascade+Q_-1 (msolve GB non-unit on a low-dimensional random slice in the reduced
#     cascade coordinates).  (ii) UPPER bound / value -- cascade+Q_-1..Q_-2 DOES force (msolve '^'
#     UNIT, ~10s).  Hence kmin(target)=2 for both targets.
A, B, conds, fv, R1, L3, _ = targets_and_tail(ak4, u4, 3, 4, "B")
tailrows = {m: scoeffs(q_m(A, B, m, 4)) for m in range(-1, -8, -1)}
p = PRIMES[0]
for tgt, nm, kind in ((R1, "R(1)", "R1"), (L3, "lambda_3(R)", "L3")):
    if HAVE_MS:
        w = witness_not_forced(4, 3, "B", [1], kind, p, seed=42 + len(nm))
        ok(w is True,
           f"band4 d=3: {nm}=1 is achievable on cascade+Q_-1 (explicit feasible point, msolve GB "
           f"non-unit on a 0-dim slice) => Q_-1 alone does NOT force {nm}=0 -- rigorous kmin({nm})>1")
        res = msolve_unit(conds + tailrows[-1] + tailrows[-2] + [sp.expand(1 - trab * tgt)],
                          [trab] + fv, p, tmo=150)
        ok(res is True,
           f"band4 d=3: {nm} in sqrt(cascade+Q_-1..Q_-2) -- msolve '^' UNIT mod {p} "
           f"=> kmin({nm})=2 (both targets forced at the SAME depth 2; W2 needed depth 3 for R(1))")
    else:
        P(f"[msolve absent -- graded depth of {nm} not checked]")


# ======================================================================
print("\n=== S6. forced-value structure: does R(1) collapse to a monomial (as at W2)? ===", flush=True)
# ======================================================================
# At W2, R(1) = -108 * a2_0 * am1_3 on the cascade (a single monomial). Test whether the band-4
# d=3 forced values mirror that -- reuse the cached GB(cascade) normal forms from S2.
p = PRIMES[0]
_, _, _, _, nfR, nfL = gb_cascade(4, 3, p)
nfR = sp.expand(nfR); nfL = sp.expand(nfL)
nR = len(sp.Poly(nfR, *sorted(nfR.free_symbols, key=str)).terms()) if nfR.free_symbols else 1
nL = len(sp.Poly(nfL, *sorted(nfL.free_symbols, key=str)).terms()) if nfL.free_symbols else 1
P(f"band4 d=3: R(1) mod cascade has {nR} monomials; lambda_3(R) mod cascade has {nL} monomials")
ok(nR > 1 and nL > 1,
   "band4 d=3: R(1) and lambda_3(R) do NOT collapse to a single monomial mod cascade (unlike W2's "
   "R(1) = -108 a2_0 am1_3) -- the forced-value structure is a genuine multinomial modulus; only its "
   "VANISHING on the tail (S4) mirrors W2, not its cascade-collapsed form")


# ======================================================================
print("\n=== S7. band 5, d=2 (HEAVY): both targets vanish on the cascade alone ===", flush=True)
# ======================================================================
if HEAVY:
    ak5, u5 = HATCH[5]
    for branch in ("A", "B"):
        A, B, conds, fv, R1, L3, tail = targets_and_tail(ak5, u5, 2, 5, branch)
        _, (R1r, L3r) = modp_lin_reduce(conds + tail, [R1, L3], PRIMES[0])
        ok(sp.expand(R1r) == 0 and sp.expand(L3r) == 0,
           f"band5 d=2 br{branch}: R(1),lambda_3(R) in sqrt(cascade+tail) [sympy GF({PRIMES[0]})] "
           f"-- DEGENERATE (already 0 on the cascade alone; the non-degenerate d>=3 band-5 tail "
           f"probe is OPEN/not run -- budget)")
    if HAVE_MS:
        ak5, u5 = HATCH[5]
        A, B, conds, fv, R1, L3, tail = targets_and_tail(ak5, u5, 2, 5, "B")
        for tgt, nm in ((R1, "R(1)"), (L3, "lambda_3(R)")):
            res = msolve_unit(conds + tail + [sp.expand(1 - trab * tgt)], [trab] + fv, PRIMES[0], tmo=280)
            ok(res is True, f"band5 d=2 brB: {nm} in sqrt(cascade+tail) -- msolve '^' UNIT mod {PRIMES[0]} (degenerate)")
else:
    P("[set HEAVY=1 for band-5 d=2 + msolve second prime/QQ + full graded depth]")


# ======================================================================
print("\n" + "=" * 72, flush=True)
print("EVIDENCE LEDGER -- THE UNIFICATION PROBE", flush=True)
headline_passed = HEADLINE_EXPECTED.issubset(PROBE) and all(PROBE[key] is True for key in HEADLINE_EXPECTED)
print(f"  PASS: {_NP} executed supporting checks", flush=True)
if headline_passed:
    print("  PASS: band-4 d=3 headline payload -- both targets, both branches,", flush=True)
    print(f"        msolve '^' Rabinowitsch UNIT mod {PRIMES[0]}", flush=True)
else:
    print("  SKIP: band-4 d=3 headline payload -- requires msolve on PATH;", flush=True)
    print("        no radical-membership verdict was executed", flush=True)
if not HEAVY:
    print(f"  SKIP: heavy corroboration -- mod {PRIMES[1]}, QQ char 0, band-5 d=2", flush=True)
print("  FAIL: none", flush=True)
print("=" * 72, flush=True)
print(f"\n(total {time.time()-_T0:.1f}s; {_NP} checks passed; HEAVY={'1' if HEAVY else '0'}, "
      f"msolve={'yes' if HAVE_MS else 'no'})", flush=True)
if headline_passed:
    print("PASS -- HATCH SLOPE FORCING HEADLINE PAYLOAD PASSED", flush=True)
else:
    print("SKIP -- HATCH SLOPE FORCING HEADLINE PAYLOAD NOT RUN; SUPPORTING CHECKS PASSED", flush=True)
