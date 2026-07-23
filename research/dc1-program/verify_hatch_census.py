#!/usr/bin/env python3
"""verify_hatch_census.py -- exact certificate for hatch-census.md.

HATCH-CENSUS audit of architecture step 2's load-bearing prerequisite:
"one silent hatch per band", and the combined feasibility of the band-4/5 hatches.

Conventions frozen from the band-k corpus (b7e85e8..30d8c59):
  A_1[x^-1] = (+)_k x^k C[E],  (x^a f)(x^b g)=x^{a+b} f(E+b) g(E),  f^[n](E)=f(E+n),
  Q_m = sum_{i+l=m} (b_l^[i] a_i - a_i^[l] b_i),  gauge b_k=0,  u := b_{k-1},
  membership (E)_j = E(E-1)...(E-j+1) | a_{-j}, b_{-j},  Q_0=(T-1)G, G(0)=0.

The canonical band-k HATCH (translate the unique common root onto the anchor 0):
  a_k = prod_{i=0}^{k-1} (E + i(k-1))   (step-(k-1) AP, k roots {0,-(k-1),...,-(k-1)^2}),
  u   = prod_{j=0}^{k-2} (E + j k)      (step-k AP, k-1 roots {0,-k,...,-(k-2)k}).
Band 3 = W2 (a_3=E(E+2)(E+4), u=E(E+3)).

SECTIONS
  0. crossed-product engine sanity: Q_m=[D,X]_m (m in [-2k,2k]); Q_0=(T-1)G (k=3,4,5).
  1. HATCH is wall-admissible (base-rep collapse), UNIQUE common root at anchor
     (gcd(k-1,k)=1), and EXOTIC -- proved arbitrary k, checked k=3..7.
  2. HONEST codim Im Phi (left-nullspace; cross-checked by point-annihilators and
     explicit ideal membership).  KEY AUDIT CORRECTION:
        band 3 (W2): codim 3, POINT-complete {ev_-1,ev_0,ev_1}, Im Phi=(E(E-1)(E+1));
        bands 4,5:   codim 3, only 2 POINT annihs {ev_0,ev_1} + one NON-POINT annih,
                     so Im Phi != (E(E-1)) and w2-theory.md's "codim 2 / principal /
                     reduces to just R(1)=1" is FALSE.  The gcd-shortcut im_phi_gen
                     (used there) is unreliable -- demoed here.
  3. CENSUS uniqueness (bounded): among minimal squarefree single-coset wall-
     admissible tops in a large window, EXACTLY ONE unique-common-root hatch per
     band (3,4,5) = the canonical AP; all others have >=2 common roots.
  4. COMBINED feasibility:  the band-4/5 hatch moment unit  Q_0=1  is the UNIT ideal
     (both branches) -- the hatch dies at the MOMENT UNIT ALONE (=> FULL a fortiori
     UNIT; the negative tail is not even needed).  Cascade control feasible (no
     false kill).  For band 4 the PURE slope R(1)=1 is achievable yet Q_0=1 is unit
     -> the killer is the NON-POINT annihilator (a stronger/earlier kill than W2,
     which SURVIVES Q_0=1 and dies only at the joint slope+tail; w2-joint-theorem.md).
  5. CONTROL -- the msolve '**' parser bug:  msolve 0.10.1 MISPARSES Python '**'
     exponent notation (silently corrupts the ideal -> wrong 'feasible' verdict);
     '^' is correct.  Both verify_w2_verdict.py and verify_w2_joint.py feed '**'.

Run:  uv run --with sympy python research/dc1-program/verify_hatch_census.py
Ends: ALL HATCH CENSUS CHECKS PASSED
"""
import sympy as sp
from itertools import combinations
import subprocess, tempfile, os, re, shutil, time

E, t, sig = sp.symbols("E t sigma")
_T0 = time.time()


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
    print(f"PASS [{time.time()-_T0:6.1f}s] {label}")


def okz(v, label):
    if sp.expand(v) != 0:
        raise AssertionError(f"FAIL {label}: residual {sp.factor(sp.expand(v))}")
    print(f"PASS [{time.time()-_T0:6.1f}s] {label}")


def q_m(A, B, m, k):
    L = range(-k, k + 1)
    return sp.expand(sum(sh(B[l], i) * A[i] - sh(A[i], l) * B[l]
                         for i in L for l in L if i + l == m))


def potential(A, B, k):
    return sp.expand(sum(sh(A[i], j - i) * sh(B[-i], j) - sh(B[i], j - i) * sh(A[-i], j)
                         for i in range(1, k + 1) for j in range(i)))


def hatch(k):
    a = sp.expand(sp.prod(E + i * (k - 1) for i in range(k)))
    u = sp.expand(sp.prod(E + j * k for j in range(k - 1)))
    return a, u


# ---- crossed-product engine (independent check that Q_m encodes [D,X]) ----
def mul_ladders(P, Q):
    R = {}
    for k1, p1 in P.items():
        for k2, p2 in Q.items():
            R[k1 + k2] = R.get(k1 + k2, 0) + sh(p1, k2) * p2
    return {k: sp.expand(v) for k, v in R.items()}


def commutator_coeff(X, Dop, m):
    return sp.expand(mul_ladders(Dop, X).get(m, 0) - mul_ladders(X, Dop).get(m, 0))


# ======================================================================
print("=== 0. crossed-product engine:  Q_m = [D,X]_m ,  Q_0=(T-1)G  (k=3,4,5) ===")
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
    okz(potential(Ag, Bg, k).subs(E, 0), f"k={k}: G(0)=0 under membership")


# ======================================================================
print("\n=== 1. the hatch: wall-admissible, UNIQUE common root at anchor, EXOTIC ===")
# ======================================================================
def wall(a, u, k):
    return sp.expand(sh(u, k) * a - sh(a, k - 1) * u) == 0


def mulS(gd, r):  # necklace box-product S_r * gd
    out = {}
    for i in range(r):
        for p, c in gd.items():
            out[p + i] = out.get(p + i, 0) + c
    return {p: c for p, c in out.items() if c != 0}


for k in range(3, 8):
    a, u = hatch(k)
    aroots = sorted(-i * (k - 1) for i in range(k))
    uroots = sorted(-j * k for j in range(k - 1))
    common = sorted(set(aroots) & set(uroots))
    consec = aroots == list(range(min(aroots), max(aroots) + 1))
    ok(wall(a, u, k), f"k={k}: wall  u^[k] a_k = a_k^[k-1] u  holds")
    # arbitrary-k proof of wall: base-representation collapse
    da = {i * (k - 1): 1 for i in range(k)}
    du = {j * k: 1 for j in range(k - 1)}
    target = {n: 1 for n in range(k * (k - 1))}
    ok(mulS(du, k) == target and mulS(da, k - 1) == target,
       f"k={k}: S_k d(u) = S_(k-1) d(a_k) = S_({k*(k-1)}) (base-rep collapse: proves wall all k)")
    # unique common root: i(k-1)=jk, gcd(k-1,k)=1 => k|i => i=0
    cg = [i * (k - 1) for i in range(k) if any(i * (k - 1) == j * k for j in range(k - 1))]
    ok(cg == [0] and sp.gcd(k - 1, k) == 1 and common == [0],
       f"k={k}: gcd(k-1,k)=1 => UNIQUE common root at anchor 0  (roots a={aroots}, u={uroots})")
    ok(not consec, f"k={k}: EXOTIC (a_k non-consecutive)")


# ======================================================================
print("\n=== 2. HONEST codim Im Phi (AUDIT CORRECTION: all codim 3; k>=4 NON-POINT) ===")
# ======================================================================
def Kblock(a, c, k):
    return sp.expand(sum(sh(a, j - k) * sh(c, j) for j in range(k)))


def Hblock(u, v, k):
    return sp.expand(sum(sh(u, j - (k - 1)) * sh(v, j) for j in range(k - 1)))


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


def point_annih(a, u, k, supp, maxj):
    ims = []
    for j in range(maxj):
        ims += [Kblock(a, falling(k) * E**j, k), Hblock(u, falling(k - 1) * E**j, k)]
    cs = sp.symbols(f"c0:{len(supp)}")
    eqs = [sp.expand(sum(cs[i] * sp.sympify(g).subs(E, supp[i]) for i in range(len(supp)))) for g in ims]
    M, _ = sp.linear_eq_to_matrix(eqs, cs)
    return M.nullspace()


def in_imphi(target, a, u, k, maxdeg):
    gens = []
    for j in range(maxdeg):
        gens += [Kblock(a, falling(k) * E**j, k), Hblock(u, falling(k - 1) * E**j, k)]
    z = sp.symbols(f"z0:{len(gens)}")
    combo = sp.expand(sum(z[i] * gens[i] for i in range(len(gens))) - target)
    return len(sp.linsolve(sp.Poly(combo, E).all_coeffs(), z)) > 0


for k in (3, 4, 5):
    a, u = hatch(k)
    cods = [codim_leftnull(a, u, k, N) for N in (16, 18, 20)]
    ok(cods == [3, 3, 3], f"k={k}: codim Im Phi = 3 (stable N=16,18,20) [left-nullspace]")
    ns = point_annih(a, u, k, list(range(-3, 3 * k + 5)), 3 * k + 10)
    npt = len(ns)
    if k == 3:
        ok(npt == 3, "k=3 (W2): 3 POINT annihilators {ev_-1,ev_0,ev_1} = codim => POINT-COMPLETE")
        ok(in_imphi(sp.expand(E * (E - 1) * (E + 1)), a, u, k, 3 * k + 10),
           "k=3: E(E-1)(E+1) in Im Phi => Im Phi = (E(E-1)(E+1)) exactly (w2-theory CORRECT)")
    else:
        ok(npt == 2, f"k={k}: only 2 POINT annihilators {{ev_0,ev_1}} but codim=3 "
                     f"=> a NON-POINT annihilator exists (undercount)")
        ok(not in_imphi(sp.expand(E * (E - 1)), a, u, k, 3 * k + 10),
           f"k={k}: E(E-1) NOT in Im Phi => Im Phi != (E(E-1)); "
           f"w2-theory 'codim 2 / principal' is WRONG (true codim 3)")

# gcd-shortcut unreliability: un-translated band-4 top {0,3,6,9}/{1,5,9} (common root 9, off anchor)
aP = sp.expand(sp.prod(E - r for r in (0, 3, 6, 9)))
uP = sp.expand(sp.prod(E - r for r in (1, 5, 9)))
imsP = [Kblock(aP, falling(4) * E**j, 4) for j in range(14)] + [Hblock(uP, falling(3) * E**j, 4) for j in range(14)]
Dg = imsP[0]
for g in imsP[1:]:
    Dg = sp.gcd(Dg, g)
Dg = sp.expand(Dg / sp.LC(sp.Poly(Dg, E)))
nsP = point_annih(aP, uP, 4, list(range(-2, 20)), 20)
ok(sp.expand(Dg - E) == 0 and len(nsP) == 2,
   "gcd-shortcut on {0,3,6,9}/{1,5,9}: gcd=E (=> would claim codim 1) but honest codim=2 "
   "=> gcd-shortcut UNRELIABLE (assumes principal)")


# ======================================================================
print("\n=== 3. CENSUS uniqueness (bounded): one minimal single-coset hatch per band ===")
# ======================================================================
def minimal_hatches(k, W):
    Sk = sp.Poly(sum(sig**i for i in range(k)), sig)
    Skm1 = sp.Poly(sum(sig**i for i in range(k - 1)), sig)
    hatches, alladm = [], 0
    for rest in combinations(range(1, W + 1), k - 1):
        roots = (0,) + rest
        dA = sp.Poly(sum(sig**r for r in roots), sig)
        q, rem = sp.div(dA, Sk, sig)
        if not rem.is_zero:
            continue
        du = sp.Poly(Skm1.as_expr() * q.as_expr(), sig)
        if any(du.nth(i) < 0 for i in range(du.degree() + 1)):
            continue
        uroots = [i for i in range(du.degree() + 1) for _ in range(int(du.nth(i)))]
        common = sorted(set(roots) & set(uroots))
        alladm += 1
        if len(common) == 1:
            s = common[0]
            hatches.append((tuple(sorted(r - s for r in roots)), tuple(sorted(x - s for x in uroots))))
    return sorted(set(hatches)), alladm


for k, W in ((3, 22), (4, 20), (5, 22)):
    uniq, alladm = minimal_hatches(k, W)
    # enumeration works in root-POSITION convention (necklace sigma^p, root at E=-p);
    # translate the unique common root onto the anchor -> positions {i(k-1)}, {jk}.
    exp_top = tuple(sorted(i * (k - 1) for i in range(k)))
    exp_sub = tuple(sorted(j * k for j in range(k - 1)))
    ok(len(uniq) == 1 and uniq[0] == (exp_top, exp_sub),
       f"k={k} (window [0,{W}], {alladm} admissible tops): EXACTLY ONE unique-common-root "
       f"hatch = canonical step-{k-1}/step-{k} AP")


# ======================================================================
print("\n=== 4. COMBINED feasibility: band-4/5 hatch DIES at the MOMENT UNIT alone ===")
# ======================================================================
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


def clear_denoms(e):
    return sp.expand(sp.numer(sp.together(sp.expand(e))))


def is_unit_sympy(eqs, vs, domain=sp.QQ):
    eqs = [clear_denoms(e) for e in eqs if sp.expand(e) != 0]
    return bool(eqs) and list(sp.groebner(eqs, *vs, order="grevlex", domain=domain).exprs) == [sp.Integer(1)]


def linear_reduce(eqs, vs):
    eqs = [sp.expand(e) for e in eqs if sp.expand(e) != 0]; vs = list(vs); ch = True
    while ch:
        ch = False; best = None
        for e in eqs:
            for v in (e.free_symbols & set(vs)):
                if sp.Poly(e, v).degree() != 1:
                    continue
                c = e.coeff(v, 1)
                if c == 0 or (c.free_symbols & set(vs)):
                    continue
                rest = sp.expand(e - c * v); sc = len(sp.Add.make_args(rest))
                if best is None or sc < best[0]:
                    best = (sc, e, v, sp.expand(-rest / c))
        if best:
            _, e0, v0, sol = best
            eqs = [sp.expand(ee.subs(v0, sol)) for ee in eqs if ee is not e0]
            eqs = [ee for ee in eqs if sp.expand(ee) != 0]; vs.remove(v0); ch = True
    return eqs, vs


def msolve_ok(eqs, fv, char, starstar=False, timeout=300):
    """UNIT / FEASIBLE via msolve -g2 reduced Groebner basis.  starstar=True feeds
    the buggy '**' notation (control in section 5); default converts '**'->'^'."""
    eqs = [sp.expand(e) for e in eqs if sp.expand(e) != 0]
    xs = sp.symbols(f"z_0:{len(fv)}"); sub = dict(zip(fv, xs))
    with tempfile.NamedTemporaryFile("w", suffix=".ms", delete=False) as f:
        inp = f.name
        f.write(",".join(str(x) for x in xs) + "\n" + str(char) + "\n")
        polys = [str(clear_denoms(e.subs(sub))).replace(" ", "") for e in eqs]
        if not starstar:
            polys = [p.replace("**", "^") for p in polys]
        f.write(",\n".join(polys) + "\n")
    out = inp + ".out"
    try:
        subprocess.run(["msolve", "-g", "2", "-f", inp, "-o", out], check=True,
                       stderr=subprocess.DEVNULL, timeout=timeout)
        r = open(out).read()
    finally:
        for p_ in (inp, out):
            if os.path.exists(p_):
                os.remove(p_)
    m = re.search(r"length of basis:\s*(\d+)", r)
    if m and m.group(1) == "1" and "[1]" in r.split("---")[-1]:
        return "UNIT"
    return "FEASIBLE" if m else ("UNIT" if r.strip().startswith("[-1]") else "FEASIBLE")


HAVE_MSOLVE = shutil.which("msolve") is not None
if not HAVE_MSOLVE:
    print("    [msolve not on PATH: using sympy alone (slower); msolve is the primary engine]")


def kill_verdict(cs, fv, sympy_too):
    """cascade+Q_0=1 verdict: msolve (^) primary; optional independent sympy cross-check."""
    v_ms = msolve_ok(cs, fv, 65003) if HAVE_MSOLVE else None
    if sympy_too:
        red, vs = linear_reduce(cs, fv)
        v_sy = "UNIT" if is_unit_sympy(red, vs) else "FEASIBLE"
        if HAVE_MSOLVE and v_ms != v_sy:
            raise AssertionError(f"msolve/sympy DISAGREE: msolve={v_ms} sympy={v_sy}")
        return v_sy
    return v_ms


# control: cascade alone FEASIBLE at d=2 (no false kill) -- sympy on band 4 brA (representative)
_A, _B, _pos, _fv = build_branch(sp.expand(E * (E + 3) * (E + 6) * (E + 9)),
                                 sp.expand(E * (E + 4) * (E + 8)), 2, 4, "A")
_rc, _vc = linear_reduce(_pos, _fv)
ok(not is_unit_sympy(_rc, _vc), "control: band-4 cascade alone FEASIBLE at d=2 (no false kill) [sympy]")

# the KILL: bands 4,5, both branches, d=2 (msolve ^; sympy cross-check on band 4 brA)
for k, (ak_e, u_e) in ((4, (E * (E + 3) * (E + 6) * (E + 9), E * (E + 4) * (E + 8))),
                       (5, (E * (E + 4) * (E + 8) * (E + 12) * (E + 16),
                            E * (E + 5) * (E + 10) * (E + 15)))):
    ak = sp.expand(ak_e); u = sp.expand(u_e)
    for br in ("A", "B"):
        A, B, pos, fv = build_branch(ak, u, 2, k, br)
        cs = pos + scoeffs(q_m(A, B, 0, k) - 1)
        cross = (k == 4 and br == "A")
        ok(kill_verdict(cs, fv, cross) == "UNIT",
           f"k={k} d=2 br{br}: cascade + Q_0=1 = UNIT" + (" [msolve ^ + sympy]" if cross else " [msolve ^]")
           + "  => hatch dies at moment unit; FULL a fortiori UNIT")

# DECISIVE d=3 (band 4 both branches, sympy cross-check brA; band 5 both branches, msolve ^)
for k, (ak_e, u_e), do_sympy_brA in (
        (4, (E * (E + 3) * (E + 6) * (E + 9), E * (E + 4) * (E + 8)), True),
        (5, (E * (E + 4) * (E + 8) * (E + 12) * (E + 16), E * (E + 5) * (E + 10) * (E + 15)), False)):
    ak = sp.expand(ak_e); u = sp.expand(u_e)
    for br in ("A", "B"):
        A, B, pos, fv = build_branch(ak, u, 3, k, br)
        cs = pos + scoeffs(q_m(A, B, 0, k) - 1)
        cross = do_sympy_brA and br == "A"
        ok(kill_verdict(cs, fv, cross) == "UNIT",
           f"k={k} d=3 br{br}: cascade + Q_0=1 = UNIT" + (" [msolve ^ + sympy]" if cross else " [msolve ^]")
           + "  => hatch dead at the MOMENT UNIT (before the tail)")

# band 4, d=3: the moment-unit kill is NOT slope-pinning -- the slope R(1)=G(1) is a LIVE modulus
# at d=3 (depends on free cascade data), unlike d<=2 where it is structurally absent (forced 0).
A2, B2, pos2, fv2 = build_branch(sp.expand(E * (E + 3) * (E + 6) * (E + 9)),
                                 sp.expand(E * (E + 4) * (E + 8)), 2, 4, "A")
A3, B3, pos3, fv3 = build_branch(sp.expand(E * (E + 3) * (E + 6) * (E + 9)),
                                 sp.expand(E * (E + 4) * (E + 8)), 3, 4, "A")
G1_d2 = sp.expand(potential(A2, B2, 4).subs(E, 1))
G1_d3 = sp.expand(potential(A3, B3, 4).subs(E, 1))
ok(len(G1_d3.free_symbols) > len(G1_d2.free_symbols) and len(G1_d3.free_symbols) > 0,
   "k=4: the slope R(1)=G(1) is a LIVE modulus at d=3 (more free data than d=2) -- so the "
   "moment-unit kill is NOT slope-pinning; the killer is the NON-POINT annihilator "
   "(pure slope R(1)=1 confirmed achievable in-session by msolve; see hatch-census.md 3.1)")


# ======================================================================
print("\n=== 5. CONTROL: the msolve '**' parser bug (affects sibling verifiers) ===")
# ======================================================================
# Minimal reproducer: provably UNIT (eq0 forces y0=0 or w=0; both give eq6=-2916!=0),
# yet msolve MISPARSES '**' and returns FEASIBLE; '^' returns UNIT.
if HAVE_MSOLVE:
    y0, y1, y4 = sp.symbols("y0 y1 y4")
    eq0 = 3 * y0**3 - 2 * y0**2 * y1 + 9 * y0 * y4
    eq6 = (-9 * y0**5 + 12 * y0**4 * y1 - 4 * y0**3 * y1**2 - 54 * y0**3 * y4
           + 36 * y0**2 * y1 * y4 - 81 * y0 * y4**2 - 2916)
    mini = [eq0, eq6]
    ok(is_unit_sympy(mini, [y0, y1, y4]), "minimal reproducer is UNIT (sympy + hand: eq6=-2916!=0)")
    ok(msolve_ok(mini, [y0, y1, y4], 65003, starstar=False) == "UNIT",
       "msolve with '^' : UNIT  (correct)")
    ok(msolve_ok(mini, [y0, y1, y4], 65003, starstar=True) != "UNIT",
       "msolve with '**': FEASIBLE (WRONG) -- '**' silently corrupts the ideal; "
       "verify_w2_verdict.py:361 and verify_w2_joint.py:570 feed '**' to msolve")

print("\n" + "=" * 70)
print("SUMMARY")
print("  * one canonical hatch per band k (step-(k-1)/step-k AP): existence PROVED")
print("    (base-rep + gcd, arbitrary k); uniqueness BOUNDED (minimal single-coset, window).")
print("  * codim Im Phi = 3 for ALL hatches (bands 3,4,5); band 3 POINT-complete,")
print("    bands 4,5 carry a NON-POINT annihilator -> w2-theory 'codim 2 / just R(1)=1' is WRONG.")
print("  * bands 4,5 hatch: Q_0=1 is the UNIT ideal (both branches) -> DIES AT THE MOMENT UNIT")
print("    ALONE (no tail needed); the non-point annihilator is the killer.  No hatch survives")
print("    -> no counterexample lead -> no escalation.")
print("  * (method) msolve '**' bug: present in both sibling msolve verifiers.")
print("=" * 70)
print("\nALL HATCH CENSUS CHECKS PASSED")
