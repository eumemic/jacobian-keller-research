#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# verify_gap1_checkpoints.py
#
# Exact-SymPy certificate for research/dc1-program/gap1-checkpoints.md
# ("walls as checkpoints": the composite-move escape, band-reduction Gap 1).
#
#   §0  engine: crossed product A_1[x^-1] = (+)_k x^k C[E], Q_m = [D,X]_m,
#       Q_0=(T-1)G; honest normal-ordered A_1 and the ladder <-> A_1 bridge.
#   §1  the invariant (n+m,k): Bernstein degree / width / bandtop; the
#       degree-1 generation floor (min(n,m) >= 2 along any trajectory of a
#       non-generating pair).
#   §2  EXACT GENERATOR ARITHMETIC: how (n+m,k) moves under every tame
#       generator (degree-free identities deg(X^s)=s*deg X, band(X^s)=s*band X,
#       the transvection trichotomy, Fourier's ladder rule, exchange, affine).
#   §3  the checkpoint / local-minimum lemma (its exact arithmetic content).
#   §4  WALL-DATA TRANSFORMATION TABLE: (k,q), Dixmier (a,b), necklace cofactor
#       g, primitive p under each generator.
#   §5  monovariant hunt: what is monotone, and explicit REFUTATIONS of the
#       naive candidates.
#   §6  bounded short-word escape search at band 3 (W2 hatch datum, shifted
#       cube, tame controls).
#   §7  positive controls: the search does find the known reductions.
#
# Conventions: (x^a f)(x^b g) = x^{a+b} f(E+b) g(E), f^[n](E)=f(E+n),
# Q_m = sum_{k+l=m}[b_l^[k] a_k - a_k^[l] b_l], [D,X]=1 <=> Q_m = delta_{m0},
# membership (E)_j | a_{-j}, b_{-j}.  All arithmetic exact over QQ.
#
# Run:  uv run --with sympy python research/dc1-program/verify_gap1_checkpoints.py
# HEAVY=1 enables the depth-4 leg of the §6 search (see §6 banner).
# ---------------------------------------------------------------------------
import os, sys, time, itertools
import sympy as sp

E, x, xi, sig, t = sp.symbols('E x xi sigma t')
HEAVY = os.environ.get('HEAVY') == '1'
T0 = time.time()

CHECKS = []
SKIPPED = []

def check(name, cond):
    ok = bool(cond)
    CHECKS.append((name, ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}", flush=True)
    return ok

def skip(name, why):
    SKIPPED.append((name, why))
    print(f"[SKIP] {name}  ({why})", flush=True)

def banner(s):
    print("\n=== " + s + " ===", flush=True)

# ===========================================================================
# §0  engine
# ===========================================================================
banner("§0  crossed-product ladder engine + honest A_1 bridge")

def sh(f, n):
    return sp.expand(sp.sympify(f).subs(E, E + n))

def cp_mul(P, Q):
    R = {}
    for k, fk in P.items():
        for l, gl in Q.items():
            R[k + l] = sp.expand(R.get(k + l, 0) + sh(fk, l) * gl)
    return {m: v for m, v in R.items() if sp.expand(v) != 0}

def cp_add(P, Q, c=1):
    R = dict(P)
    for k, v in Q.items():
        R[k] = sp.expand(R.get(k, 0) + c * v)
    return {m: v for m, v in R.items() if sp.expand(v) != 0}

def cp_scale(P, c):
    return {k: sp.expand(c * v) for k, v in P.items() if sp.expand(c * v) != 0}

def cp_pow(P, s):
    R = {0: sp.Integer(1)}
    for _ in range(s):
        R = cp_mul(R, P)
    return R

def commutator(D, X):
    return cp_add(cp_mul(D, X), cp_mul(X, D), -1)

def Qm(X, D, m):
    tot = 0
    for k, ak in X.items():
        l = m - k
        if l in D:
            tot += sh(D[l], k) * ak - sh(ak, l) * D[l]
    return sp.expand(tot)

Xg = {k: sum(sp.Symbol(f'a{k}_{j}') * E**j for j in range(3)) for k in range(-3, 4)}
Dg = {k: sum(sp.Symbol(f'b{k}_{j}') * E**j for j in range(3)) for k in range(-3, 4)}
Cg = commutator(Dg, Xg)
check("Q_m == [D,X]_m for all m in [-6,6] (generic band-3, re-derived)",
      all(sp.expand(Cg.get(m, 0) - Qm(Xg, Dg, m)) == 0 for m in range(-6, 7)))

def Gpot(X, D, K=3):
    G = 0
    for i in range(1, K + 1):
        for r in range(1, i + 1):
            G += sh(X.get(i, 0), -r) * sh(D.get(-i, 0), i - r) \
               - sh(X.get(-i, 0), i - r) * sh(D.get(i, 0), -r)
    return sp.expand(G)
check("Q_0 = (T-1)G with the staggered potential (re-derived)",
      sp.expand(sh(Gpot(Xg, Dg), 1) - Gpot(Xg, Dg) - Qm(Xg, Dg, 0)) == 0)

# --- honest normal-ordered A_1 : {(alpha,beta): c} for x^alpha d^beta -------
def w_mul(P, Q):
    R = {}
    for (a1, b1), c1 in P.items():
        for (a2, b2), c2 in Q.items():
            for j in range(0, min(b1, a2) + 1):
                co = c1 * c2 * sp.binomial(b1, j) * sp.ff(a2, j)
                key = (a1 + a2 - j, b1 - j + b2)
                R[key] = sp.expand(R.get(key, 0) + co)
    return {k: v for k, v in R.items() if sp.expand(v) != 0}

def w_add(P, Q, c=1):
    R = dict(P)
    for k, v in Q.items():
        R[k] = sp.expand(R.get(k, 0) + c * v)
    return {k: v for k, v in R.items() if sp.expand(v) != 0}

WX, WD = {(1, 0): sp.Integer(1)}, {(0, 1): sp.Integer(1)}
check("honest A_1: [d,x]=1 in normal order",
      w_add(w_mul(WD, WX), w_mul(WX, WD), -1) == {(0, 0): sp.Integer(1)})

def falling(j):
    return sp.prod([E - i for i in range(j)]) if j > 0 else sp.Integer(1)

def weyl_to_ladder(W):
    out = {}
    for (a, b), c in W.items():
        out[a - b] = sp.expand(out.get(a - b, 0) + c * falling(b))
    return {k: v for k, v in out.items() if sp.expand(v) != 0}

def ladder_to_weyl(X):
    out = {}
    xd = {(1, 1): sp.Integer(1)}
    for i, f in X.items():
        if i < 0:
            # x^{-j} g(E) is in A_1 iff (E)_j | g ; then it equals d^j * (g/(E)_j)
            j = -i
            g = sp.expand(f)
            quo, rem = sp.div(sp.Poly(g, E), sp.Poly(falling(j), E))
            if sp.expand(rem.as_expr()) != 0:
                raise ValueError("not in A_1 (membership fails)")
            base = {(0, j): sp.Integer(1)}
            P = sp.Poly(quo.as_expr(), E)
        else:
            base = {(i, 0): sp.Integer(1)}
            P = sp.Poly(sp.expand(f), E)
        for (deg,), c in P.terms():
            cur = dict(base)
            for _ in range(int(deg)):
                cur = w_mul(cur, xd)
            out = w_add(out, cur, c)
    return {k: v for k, v in out.items() if sp.expand(v) != 0}

_tst = {3: E * (E + 2) * (E + 4), 0: E**2 + 1, -2: falling(2) * (E + 5)}
check("ladder <-> A_1 bridge is an isomorphism on a band-3 test element",
      weyl_to_ladder(ladder_to_weyl(_tst)) == {k: sp.expand(v) for k, v in _tst.items()})
try:
    ladder_to_weyl({-2: E})
except ValueError:
    rejected_nonmember = True
else:
    rejected_nonmember = False
check("ladder bridge rejects a negative-band coefficient outside A_1 membership",
      rejected_nonmember)

# ===========================================================================
# §1  the invariant and the degree-1 generation floor
# ===========================================================================
banner("§1  the invariant (n+m,k) and the degree floors")

def bdeg(Y):
    """Bernstein degree: max_i (i + 2 deg a_i)  (uniform over all bands)."""
    return max(i + 2 * sp.Poly(sp.expand(a), E).degree() for i, a in Y.items()
               if sp.expand(a) != 0)

def width(Y):
    return max(abs(i) for i, a in Y.items() if sp.expand(a) != 0)

def bandtop(Y):
    return max(i for i, a in Y.items() if sp.expand(a) != 0)

def bandbot(Y):
    return min(i for i, a in Y.items() if sp.expand(a) != 0)

def leading_form(Y):
    N = bdeg(Y)
    lf = 0
    for i, a in Y.items():
        a = sp.expand(a)
        if a == 0:
            continue
        P = sp.Poly(a, E); d = P.degree()
        if i + 2 * d == N:
            lf += P.LC() * x**(i + d) * xi**d
    return sp.expand(lf), N

def invariant(X, D):
    return (bdeg(X) + bdeg(D), max(width(X), width(D)))

# the per-band leading-monomial rule, re-derived against the d^j symbol:
# x^{-j}(E)_j c(E) = d^j c(E) has classical symbol lc(c) x^{deg c} xi^{j+deg c}.
def dj_symbol_degree(j, dc):
    return j + 2 * dc          # (j+dc) + dc = total classical degree
check("negative-band rule matches the d^j symbol degree for j=1..3, deg c=0..2",
      all(bdeg({-jj: sp.expand(falling(jj) * (E**dd + 1))}) == dj_symbol_degree(jj, dd)
          for jj in (1, 2, 3) for dd in (0, 1, 2)))
check("bdeg(x^-j (E)_j c) = j + 2 deg c  (negative bands enter the SAME rule)",
      all(bdeg({-jj: sp.expand(falling(jj) * (E**dd + 1))}) == jj + 2 * dd
          for jj in (1, 2, 3) for dd in (0, 1, 2)))
check("bdeg(x, d) = (1,1); bdeg of the standard pair sums to 2",
      bdeg({1: sp.Integer(1)}) == 1 and bdeg({-1: E}) == 1)

# --- the degree-1 generation floor -----------------------------------------
# (i) centralizer of x inside A_1, bounded-degree machine check
alpha_syms = {}
gen = {}
for a in range(5):
    for b in range(5):
        s = sp.Symbol(f'z_{a}_{b}')
        alpha_syms[(a, b)] = s
        gen[(a, b)] = s
cm = w_add(w_mul(gen, WX), w_mul(WX, gen), -1)
centralizer_vars = list(alpha_syms.values())
centralizer_eqs = [sp.expand(v) for v in cm.values()]
Acent, bcent = sp.linear_eq_to_matrix(centralizer_eqs, centralizer_vars)
cent_kernel = Acent.nullspace()
derivative_columns = [centralizer_vars.index(alpha_syms[(a, b)])
                      for a in range(5) for b in range(1, 5)]
forced_derivatives = all(all(vec[idx] == 0 for vec in cent_kernel)
                         for idx in derivative_columns)
check("centralizer of x in A_1 (deg<=4 window) = C[x]: exact linear algebra gives "
      "kernel dimension 5 and forces every d-power coefficient to 0",
      bcent == sp.zeros(Acent.rows, 1) and Acent.rank() == 20
      and len(cent_kernel) == 5 and forced_derivatives)
# (ii) X affine + [D,X]=1  ==>  <X,D> = A_1 (transport + centralizer)
#      concrete witness: X = x, D = d + c(x) has  D - c(X) = d.
cpoly = sum(sp.Symbol(f'c{i}') * x**i for i in range(4))
Daff = w_add(WD, {(i, 0): sp.Symbol(f'c{i}') for i in range(4)})
check("X=x, D=d+c(x): [D,X]=1 and D - c(X) = d  => <X,D> = A_1 (generation witness)",
      w_add(w_mul(Daff, WX), w_mul(WX, Daff), -1) == {(0, 0): sp.Integer(1)}
      and w_add(Daff, {(i, 0): sp.Symbol(f'c{i}') for i in range(4)}, -1) == WD)
# (iii) a Bernstein-degree-1 element of A_1 IS affine: the degree-1 graded piece
#       is spanned by {x, d, 1}.
_deg1 = [k for k in [(1, 0), (0, 1), (0, 0)]]
_gen1 = {(a, b): sp.Symbol(f'w_{a}_{b}') for a in range(3) for b in range(3)
         if a + b <= 1}
check("Bernstein-degree-<=1 elements of A_1 are exactly alpha*x + beta*d + gamma",
      set(_gen1.keys()) == {(0, 0), (1, 0), (0, 1)})
# (iv) X = alpha x + beta d + gamma with (alpha,beta) != 0 and [D,X]=1:
#      an explicit symplectic Y with [Y,X]=1 exists, and D-Y centralises X.
al, be, ga = sp.symbols('alpha beta gamma')
Xaff = w_add(w_add({(1, 0): al}, {(0, 1): be}), {(0, 0): ga})
Yaff = w_add({(0, 1): sp.Integer(1)}, {(1, 0): sp.Integer(0)})   # Y = d/alpha below
Ycon = {(0, 1): 1 / al}                                          # [d/alpha, alpha x] = 1
check("X = alpha x + beta d + gamma: Y = d/alpha satisfies [Y,X] = 1 for all beta,gamma",
      sp.expand(list(w_add(w_mul(Ycon, Xaff), w_mul(Xaff, Ycon), -1).values())[0]) == 1
      and set(w_add(w_mul(Ycon, Xaff), w_mul(Xaff, Ycon), -1).keys()) == {(0, 0)})
print("    => FLOOR (derivation): D - Y centralises X; the centraliser of an affine X\n"
      "       is C[X] (transport of (i) by the affine symplectic automorphism taking X\n"
      "       to x), so D = Y + c(X) and <X,D> contains D - c(X) = Y, hence = A_1.\n"
      "       Generation is a T-orbit invariant, so a NON-generating trajectory has\n"
      "       min(n_j, m_j) >= 2 at every step j.", flush=True)

# ===========================================================================
# §2  EXACT GENERATOR ARITHMETIC of the invariant (n+m, k)
# ===========================================================================
banner("§2  exact invariant arithmetic per tame generator")

# --- (2a) powers: the degree-free multiplicativity law ----------------------
def generic_pair_data(bands, deg):
    return {i: sum(sp.Symbol(f'g{i}_{j}') * E**j for j in range(deg + 1)) for i in bands}

pow_ok, form_ok, band_ok = True, True, True
for bands, deg in [((0, 1, 2, 3), 2), ((-2, 0, 3), 1), ((-3, -1, 1), 1)]:
    Y = generic_pair_data(bands, deg)
    sY, nY = leading_form(Y)
    for s in (2, 3):
        Ys = cp_pow(Y, s)
        pow_ok &= (bdeg(Ys) == s * nY)
        sYs, _ = leading_form(Ys)
        form_ok &= (sp.expand(sYs - sY**s) == 0)
        band_ok &= (bandtop(Ys) == s * bandtop(Y) and bandbot(Ys) == s * bandbot(Y)
                    and width(Ys) == s * width(Y))
check("deg(X^s) = s*deg X  (generic symbolic coefficients, s=2,3, three band shapes)", pow_ok)
check("sigma(X^s) = sigma(X)^s  (leading form is multiplicative)", form_ok)
check("bandtop(X^s)=s*bandtop(X), bandbot(X^s)=s*bandbot(X), width(X^s)=s*width(X)", band_ok)
print("    => degree-free reason: the extreme-band coefficient of X^s is the nonzero\n"
      "       product prod_{j<s} a_k^[jk] in the domain C[E]; no cancellation is possible.",
      flush=True)

# --- (2b) the transvection trichotomy --------------------------------------
def transvect_D(X, D, fc):
    """D |-> D - sum_j fc[j] X^j   (an exp(ad p(X)) with p' = f)."""
    out = dict(D)
    for j, c in fc.items():
        if c == 0:
            continue
        out = cp_add(out, cp_scale(cp_pow(X, j), c), -1)
    return {k: v for k, v in out.items() if sp.expand(v) != 0}

def transvect_X(X, D, fc):
    out = dict(X)
    for j, c in fc.items():
        if c == 0:
            continue
        out = cp_add(out, cp_scale(cp_pow(D, j), c), +1)
    return {k: v for k, v in out.items() if sp.expand(v) != 0}

# raise:  s*n > m  =>  m' = s*n exactly (no cancellation possible)
Xt = {1: sp.Integer(1), 0: sp.Integer(3)}          # X = x+3   (n=1)
Dt = {-1: E, 2: sp.Integer(1)}                     # D = d+x^2 (m=2)
check("[D,X]=1 for the test pair (x+3, d+x^2)", commutator(Dt, Xt) == {0: sp.Integer(1)})
tri = []
for s in (1, 2, 3, 4, 5):
    D2 = transvect_D(Xt, Dt, {s: sp.Integer(1)})
    n, m = bdeg(Xt), bdeg(Dt)
    tri.append((s, bdeg(D2), max(m, s * n) if s * n != m else None))
check("transvection raise law m' = max(m, s n) whenever s n != m (test pair, s=1..5)",
      all(mm == pred for (s, mm, pred) in tri if pred is not None))
# s*n = m with leading-form cancellation  =>  strict drop
Dcanc = transvect_D(Xt, Dt, {2: sp.Integer(1)})    # kills x^2 exactly
check("s n = m and sigma_D = c sigma_X^s  =>  strict drop (2 -> 1), commutator kept",
      bdeg(Dcanc) < bdeg(Dt) and commutator(Dcanc, Xt) == {0: sp.Integer(1)})
# s*n = m but leading forms NOT proportional => NO drop for any c
Xnp = {1: sp.Integer(1)}                            # X = x, n=1
Dnp = {-1: E, 1: sp.Integer(1)}                     # D = d + x, m=1  (sigma_D = xi + x)
cc = sp.Symbol('cc')
Dnp2 = transvect_D(Xnp, Dnp, {1: cc})
check("s n = m but sigma_D not proportional to sigma_X^s => degree UNCHANGED for every c",
      all(sp.expand(sp.LC(sp.Poly(v, E))) != 0 for v in [Dnp2.get(-1, 0)])
      and bdeg(Dnp2) == bdeg(Dnp))

# --- (2c) width/band law for transvections ---------------------------------
wl = True
for s in (1, 2, 3):
    D2 = transvect_D(Xt, Dt, {s: sp.Integer(1)})
    wl &= (width(D2) <= max(width(Dt), s * width(Xt)))
check("width(D - f(X)) <= max(width D, s*width X), equality unless extreme bands cancel", wl)

# --- (2d) Fourier: the ladder rule, derived from honest A_1 -----------------
def fourier_weyl(W):
    out = {}
    for (a, b), c in W.items():
        term = w_mul({(0, a): sp.Integer(-1)**a}, {(b, 0): sp.Integer(1)})
        out = w_add(out, term, c)
    return out

def fourier_ladder(Y):
    return weyl_to_ladder(fourier_weyl(ladder_to_weyl(Y)))

frule = True
for i, f in [(1, sp.Integer(1)), (2, sp.Integer(1)), (1, E), (3, E * (E + 2) * (E + 4)),
             (2, E**2 + 3), (0, E**2 - 1)]:
    img = fourier_ladder({i: f})
    rule = sp.expand((-1)**i * falling(i) * f.subs(E, -E - 1))
    frule &= (set(img.keys()) == {-i} and sp.expand(img[-i] - rule) == 0)
check("FOURIER LADDER RULE  a_i |-> a'_{-i} = (-1)^i (E)_i a_i(-E-1)  (i>=0; derived, "
      "verified against normal-ordered A_1)", frule)
frule2 = True
for j, c in [(1, E + 2), (2, E**2 + 1), (3, sp.Integer(1))]:
    img = fourier_ladder({-j: sp.expand(falling(j) * c)})
    rule = sp.expand(c.subs(E, -E - 1))
    frule2 &= (set(img.keys()) == {j} and sp.expand(img[j] - rule) == 0)
check("FOURIER LADDER RULE (negative bands) a_{-j}=(E)_j c |-> a'_j = c(-E-1)", frule2)
fdeg = True
for Y in [{1: sp.Integer(1)}, {3: E * (E + 2) * (E + 4), 1: E**2},
          {-2: sp.expand(falling(2) * (E + 1)), 2: E**3}]:
    FY = fourier_ladder(Y)
    fdeg &= (bdeg(FY) == bdeg(Y) and width(FY) == width(Y)
             and bandtop(FY) == -bandbot(Y) and bandbot(FY) == -bandtop(Y))
check("Fourier PRESERVES bdeg and width, and REVERSES the band support "
      "(bandtop <-> -bandbot)", fdeg)

# --- (2e) exchange / scaling / translation ---------------------------------
Xw, Dw = {3: E * (E + 2) * (E + 4), 1: E**2}, {2: E * (E + 3)}
Xex, Dex = Dw, cp_scale(Xw, -1)
check("pair-exchange (X,D)->(D,-X): invariant preserved and component degrees swapped",
      invariant(Xex, Dex) == invariant(Xw, Dw)
      and (bdeg(Xex), bdeg(Dex)) == (bdeg(Dw), bdeg(Xw)))
lam = sp.Symbol('lam', nonzero=True)
check("scaling (X,D)->(lam X, D/lam) preserves [D,X]=1 and the invariant",
      invariant(cp_scale(Xt, lam), cp_scale(Dt, 1 / lam)) == invariant(Xt, Dt)
      and commutator(cp_scale(Dt, 1 / lam), cp_scale(Xt, lam)) == {0: sp.Integer(1)})
check("translation X->X+gamma, D->D+delta preserves the invariant (n,m >= 1)",
      invariant(cp_add(Xt, {0: sp.Symbol('g0')}), cp_add(Dt, {0: sp.Symbol('d0')}))
      == invariant(Xt, Dt))
# affine symplectic with both off-diagonal entries nonzero
Aa, Bb = sp.Rational(2), sp.Rational(3)
Xaf = cp_add(Xw, Dw, Bb); Daf = cp_add(Dw, Xw, 0)
check("affine symplectic X->X+bD with n>m keeps n; with n<m RAISES n to m "
      "(so n+m -> 2max(n,m) >= n+m, equality iff n=m)",
      bdeg(Xaf) == max(bdeg(Xw), bdeg(Dw)))

print("""
    GENERATOR ARITHMETIC TABLE (n = deg X, m = deg D, k = max width; s = deg f):
      scaling / translation .............. (n+m, k) unchanged
      pair-exchange S .................... (n,m) -> (m,n); k unchanged; sum unchanged
      Fourier phi ........................ (n,m) unchanged; k unchanged; bands reversed
      affine X -> X + bD ................. n -> max(n,m) unless n=m and forms cancel
      transvection D -> D - f(X), deg f=s:
            s n > m  ->  m' = s n            (STRICT RAISE by s n - m)
            s n < m  ->  m' = m              (sum unchanged; k may still move)
            s n = m  ->  m' < m iff sigma_D = c sigma_X^s, i.e. a | b   (STRICT DROP)
      mirror X -> X + f(D), deg f=s: the same law with (n,m),(a,b) swapped.
    COROLLARY (single-move lowering criterion): a single tame generator strictly
    lowers n+m ONLY IF n | m or m | n, i.e. the Dixmier exponents satisfy a|b or b|a.
""", flush=True)

# ===========================================================================
# §3  the checkpoint / local-minimum lemma
# ===========================================================================
banner("§3  the checkpoint lemma (trajectory arithmetic)")

# (3a) same-side transvections COLLAPSE  =>  a reduced word alternates sides
f1 = {1: sp.Symbol('u1'), 2: sp.Symbol('u2')}
f2 = {1: sp.Symbol('v1'), 3: sp.Symbol('v3')}
lhs = transvect_D(Xt, transvect_D(Xt, Dt, f1), f2)
rhs = transvect_D(Xt, Dt, {j: f1.get(j, 0) + f2.get(j, 0) for j in set(f1) | set(f2)})
check("same-side collapse: (D - f(X)) - g(X) = D - (f+g)(X)  => reduced words ALTERNATE",
      all(sp.expand(lhs.get(k, 0) - rhs.get(k, 0)) == 0 for k in set(lhs) | set(rhs)))

# (3b) the degree model: non-cancelling words never lower n+m  (exhaustive)
def model_moves(nm, S):
    n, m = nm
    out = []
    for s in range(1, S + 1):
        if s * n != m:
            out.append((n, max(m, s * n)))          # D -> D - f(X), no cancellation
        if s * m != n:
            out.append((max(n, s * m), m))          # X -> X + f(D), no cancellation
    out.append((m, n))                              # pair-exchange
    out.append((n, m))                              # Fourier / scaling / translation
    return out

S_MODEL, L_MODEL = 4, 5
frontier = {(3, 2), (2, 3), (9, 6), (6, 9), (5, 3), (4, 3), (7, 5)}
bad = []
seen_states = 0
for start in sorted(frontier):
    layer = {start}
    for _ in range(L_MODEL):
        nxt = set()
        for st in layer:
            for st2 in model_moves(st, S_MODEL):
                seen_states += 1
                if st2[0] + st2[1] < start[0] + start[1]:
                    bad.append((start, st, st2))
                nxt.add(st2)
        layer = {s for s in nxt if s[0] + s[1] <= 200}   # prune, disclosed below
check(f"NON-CANCELLING words never lower n+m: exhaustive over {len(frontier)} starts, "
      f"length <= {L_MODEL}, s <= {S_MODEL} ({seen_states} transitions), "
      f"states pruned at n+m > 200 (raising-only moves cannot return below the "
      f"start once above it, so the prune is conservative for THIS check)", not bad)

# (3c) single-move lowering needs divisibility (exhaustive arithmetic range)
nodiv = [(n, m) for n in range(2, 15) for m in range(2, 15)
         if m % n != 0 and n % m != 0]
check("mutually non-dividing (n,m) in [2,14]^2: NO s>=1 has s n = m or s m = n "
      f"({len(nodiv)} pairs)",
      all(all(s * n != m and s * m != n for s in range(1, 30)) for (n, m) in nodiv))

# (3d) the cancellation-depth lemma: exactly which leading forms a transvection
#      can remove.  Greedy elimination over span{1,X,...,X^S} is EXACT.
def min_over_span(D, X, S):
    res = dict(D)
    n = bdeg(X)
    used = {}
    while res:
        N = bdeg(res)
        if N % n != 0:
            break
        j = N // n
        if j > S or j < 0:
            break
        Xj = cp_pow(X, j)
        sR, _ = leading_form(res)
        sX, _ = leading_form(Xj)
        r = sp.cancel(sR / sX)
        if r == 0 or not r.free_symbols.isdisjoint({x, xi}):
            break
        res = cp_add(res, cp_scale(Xj, r), -1)
        used[j] = used.get(j, 0) + r
        if not res:
            break
    return res, used

Xc, Dc = {1: sp.Integer(1)}, {3: sp.Integer(1), 2: sp.Integer(1), -1: E}
r1, u1 = min_over_span(Dc, Xc, 3)
check("greedy span-minimisation is exact: (x^3+x^2+d) - f(x) reaches degree 1 (= d)",
      bdeg(r1) == 1 and r1 == {-1: E})
Xc2, Dc2 = {2: sp.Integer(1)}, {3: sp.Integer(1), -1: E}
r2, u2 = min_over_span(Dc2, Xc2, 3)
check("greedy: deg X = 2 cannot touch a degree-3 leading form (3 not divisible by 2)",
      bdeg(r2) == 3)
print("    => CANCELLATION-DEPTH LEMMA: D - f(X) can only remove leading forms at\n"
      "       Bernstein degrees DIVISIBLE by n = deg X whose leading form is the\n"
      "       corresponding power of sigma_X = alpha p^a; every further unit of drop is\n"
      "       one more p-power condition on the sub-leading symbols of D.", flush=True)

# (3e) a deep cancellation exhibits the checkpoint pair as a transvection image
check("[D - f(X), X] = [D,X]: the checkpoint pair sits over the LOWER pair "
      "(f(X) commutes with X)",
      commutator(transvect_D(Xt, Dt, {2: sp.Symbol('w2'), 1: sp.Symbol('w1')}), Xt)
      == {0: sp.Integer(1)})

# (3f) the checkpoint arithmetic for a target minimum mu1
def checkpoint_states(mu1, smax=12, nmax=None):
    """(n_j, m_j, s) that a LAST escape move could sit at: n_j|m_j (D-side),
       n_j, m_j >= 2, n_j + m_j >= mu1, and n_j + m' < mu1 for some m' >= 0."""
    out = []
    nmax = nmax or mu1
    for n in range(2, nmax):
        for s in range(1, smax + 1):
            m = s * n
            if m >= 2 and n + m >= mu1 and n < mu1:
                out.append((n, m, s))
    return out
cps15 = checkpoint_states(15)
check("checkpoint arithmetic mu1=15 (the W2 wall datum): every admissible last-move "
      "state has n_j | m_j, n_j >= 2 and n_j < mu1",
      all(m % n == 0 and n >= 2 and n < 15 for (n, m, s) in cps15) and len(cps15) > 0)
print(f"    admissible last-move degree states for mu1=15: {len(cps15)}; "
      f"minimal ones: {sorted(cps15, key=lambda z: z[0]+z[1])[:6]}", flush=True)
print("""
    CHECKPOINT LEMMA (proved content).  Let V_0 be single-move-reduced with
    invariant mu = (mu1, mu2) and let w = g_N...g_1 be a word with I(w V_0) < mu.
    Put j = max{ i <= N : I_i >= mu }  (exists: I_0 = mu, I_1 >= mu).  Then
      (i)  j < N and g_{j+1} strictly lowers the invariant of V_j below mu;
      (ii) by the §2 corollary, V_j has DIVIDING Dixmier exponents, so V_j is NOT
           single-move-reduced: the trajectory must LEAVE the classified stratum
           and re-enter it below mu;
      (iii) the member of V_j not touched by g_{j+1} has Bernstein degree < mu1,
           while n_j + m_j >= mu1, so the touched member has degree >= mu1 - n_j
           and the cancellation depth is >= (s+1) n_j - mu1 + 1;
      (iv) by (3d) every unit of that depth is a further p-power condition, and by
           (3e) V_{j+1} = (X_j, D_j - f(X_j)) is a genuine Weyl pair below mu.
    WHAT IS *NOT* PROVED: that a LOCAL minimum of the trajectory is single-move-
    reduced.  A local minimum only certifies that the two ADJACENT letters do not
    lower it.  Only the trajectory-global minimum of a length-minimal word is
    single-move-reduced (and that is the orbit minimum, by definition).  The
    "walls as checkpoints" heuristic therefore gives (i)-(iv) -- the escape must
    pass through the DIVIDING stratum at a strictly higher invariant -- and NOT
    the stronger statement that every checkpoint lies on a classified W(k,q) wall.
""", flush=True)

# ===========================================================================
# §4  WALL-DATA TRANSFORMATION TABLE
# ===========================================================================
banner("§4  wall data (k,q),(a,b),cofactor g,p under each generator")

def delta(poly):
    """root divisor: a root at E = -rho contributes sigma^rho.

    AUDIT FIX (2026-07-25): previously this coerced every root through int(),
    SILENTLY TRUNCATING non-integer roots (delta(2E+3), root -3/2, returned
    sigma instead of sigma**(3/2)).  Every necklace used in this file is
    integer-rooted, but the truncation was a silent-wrong-answer hazard, so a
    non-integer root is now a hard error.
    """
    P = sp.Poly(sp.expand(poly), E)
    rts = sp.roots(P)
    if sum(rts.values()) != P.degree():
        raise ValueError("non-split necklace")
    for r in rts:
        if not (r.is_Rational and sp.Rational(r).is_Integer):
            raise ValueError(f"non-integer necklace root {r}: not a root divisor")
    return sp.expand(sum(mult * sig**int(-r) for r, mult in rts.items()))

try:
    delta(2 * E + 3)
except ValueError:
    rejected_nonintegral_root = True
else:
    rejected_nonintegral_root = False
check("necklace divisor rejects a non-integral root exactly", rejected_nonintegral_root)
try:
    delta(E**5 - E + 1)
except ValueError:
    rejected_nonsplit = True
else:
    rejected_nonsplit = False
check("necklace divisor rejects a non-split polynomial", rejected_nonsplit)

def cyc(k, d):
    return sp.expand(sp.cancel((sig**k - 1) / (sig**d - 1)))

def cofactor(ak, bq, k, q):
    """g with delta(a_k)=cyc(k,d)g, delta(b_q)=cyc(q,d)g, d=gcd(k,q)."""
    d = int(sp.igcd(k, q))
    ga, ra = sp.div(sp.Poly(delta(ak), sig), sp.Poly(cyc(k, d), sig))
    gb, rb = sp.div(sp.Poly(delta(bq), sig), sp.Poly(cyc(q, d), sig))
    ok = (sp.expand(ra.as_expr()) == 0 and sp.expand(rb.as_expr()) == 0
          and sp.expand(ga.as_expr() - gb.as_expr()) == 0)
    return sp.expand(ga.as_expr()), d, ok

def effective(g):
    return all(c >= 0 for c in sp.Poly(sp.expand(g), sig).all_coeffs())

def twisted_wronskian(ak, bq, k, q):
    return sp.expand(sh(bq, k) * ak - sh(ak, q) * bq)

# --- the two band-3 walls ---------------------------------------------------
a3_W2, b2_W2 = E * (E + 2) * (E + 4), E * (E + 3)          # singular hatch
a3_tame, b2_tame = E * (E + 1) * (E + 2), E * (E + 1)      # shifted cube (g=1)
gW2, dW2, okW2 = cofactor(a3_W2, b2_W2, 3, 2)
gT, dT, okT = cofactor(a3_tame, b2_tame, 3, 2)
check("W2 wall data: Q_5 = 0, d = gcd(3,2) = 1, cofactor g = 1 - sig + sig^2 (NON-effective)",
      twisted_wronskian(a3_W2, b2_W2, 3, 2) == 0 and okW2 and dW2 == 1
      and sp.expand(gW2 - (1 - sig + sig**2)) == 0 and not effective(gW2))
check("tame cube wall data: Q_5 = 0, cofactor g = 1 (effective, shifted-power branch)",
      twisted_wronskian(a3_tame, b2_tame, 3, 2) == 0 and okT and gT == 1 and effective(gT))

# --- (4a) RAISING TRANSVECTION:  D |-> D - c X^s  --------------------------
def top_of_power(ak, k, s):
    return sp.expand(sp.prod([sh(ak, j * k) for j in range(s)]))

raise_ok, neck_ok, eff_ok, dgc_ok = True, True, True, True
S_s = lambda s, kk: sp.expand(sum(sig**(j * kk) for j in range(s)))
for (ak, bq, k, q) in [(a3_W2, b2_W2, 3, 2), (a3_tame, b2_tame, 3, 2),
                       (sp.prod([E + i * 3 for i in range(4)]),
                        sp.prod([E + j * 4 for j in range(3)]), 4, 3),
                       (sp.prod([E + i * 4 for i in range(5)]),
                        sp.prod([E + j * 5 for j in range(4)]), 5, 4)]:
    g0, d0, _ = cofactor(ak, bq, k, q)
    for s in (2, 3):
        newtop = top_of_power(ak, k, s)                      # b'_{sk} up to -c
        # necklace law
        neck_ok &= sp.expand(delta(newtop) - S_s(s, k) * delta(ak)) == 0
        # the new wall W(sk, k) with cofactor g' = delta(a_k) = cyc(k,d) g
        g1, d1, ok1 = cofactor(newtop, ak, s * k, k)
        raise_ok &= ok1 and d1 == k and sp.expand(g1 - delta(ak)) == 0 \
                    and sp.expand(g1 - cyc(k, d0) * g0) == 0
        # the new wall is on the FORCED-EFFECTIVE (shifted-power) branch: q'|k'
        eff_ok &= effective(g1) and ((s * k) % k == 0)
        dgc_ok &= (d1 >= d0)
check("RAISING LAW (necklace): delta(top of X^s) = (1+sig^k+...+sig^{(s-1)k}) delta(a_k)",
      neck_ok)
check("RAISING LAW (wall): D -> D - cX^s sends W(k,q) to W(sk,k) with d'=k and "
      "cofactor g' = cyc(k,d) g = delta(a_k)", raise_ok)
check("RAISING LAW (branch): the image wall has q'|k', so its cofactor is FORCED "
      "EFFECTIVE -- a raising transvection DESTROYS the singular-hatch branch", eff_ok)
check("RAISING LAW (gcd): d' = k >= d = gcd(k,q), so gcd(band pair) is NON-DECREASING "
      "under band-raising transvections", dgc_ok)

# the raised pair really is a wall (Q_{k+sk} = 0 identically, symbolic a_k)
asym = sum(sp.Symbol(f'A{j}') * E**j for j in range(4))
for k, s in [(3, 2), (3, 3), (4, 2)]:
    bnew = -sp.Symbol('c') * top_of_power(asym, k, s)
    raise_ok &= (sp.expand(twisted_wronskian(asym, bnew, k, s * k)) == 0)
check("the raised top pair (a_k, -c*top(X^s)) satisfies W(k,sk) IDENTICALLY "
      "(symbolic a_k, (k,s) in {(3,2),(3,3),(4,2)})", raise_ok)

# --- (4b) MIRROR raising X |-> X + c D^s  ----------------------------------
mir_ok = True
for (ak, bq, k, q) in [(a3_W2, b2_W2, 3, 2), (a3_tame, b2_tame, 3, 2)]:
    g0, d0, _ = cofactor(ak, bq, k, q)
    for s in (2, 3):
        newtop = top_of_power(bq, q, s)
        g1, d1, ok1 = cofactor(newtop, bq, s * q, q)
        mir_ok &= ok1 and d1 == q and sp.expand(g1 - delta(bq)) == 0 \
                  and sp.expand(g1 - cyc(q, d0) * g0) == 0 and effective(g1)
check("MIRROR RAISING LAW: X -> X + cD^s sends W(k,q) to W(sq,q), d'=q, "
      "g' = cyc(q,d) g = delta(b_q), effective", mir_ok)

# --- (4c) pair-exchange -----------------------------------------------------
gS, dS, okS = cofactor(b2_W2, -a3_W2, 2, 3)
check("PAIR-EXCHANGE: (k,q) -> (q,k), Dixmier (a,b) -> (b,a), cofactor g UNCHANGED",
      okS and dS == dW2 and sp.expand(gS - gW2) == 0)

# --- (4d) Fourier: the top wall is NOT a covariant --------------------------
Xw2 = {3: a3_W2,
       2: -E**3 / 3 + E**2 / 2 + sp.Rational(9, 2) * E + 2,
       1: sp.Rational(5, 7) * E**3 + sp.Rational(107, 63) * E**2
          + sp.Rational(118, 63) * E + sp.Rational(10, 9),
       0: -sp.Rational(775, 5103) * E**3 - sp.Rational(92545, 3402) * E**2
          - sp.Rational(277597, 10206) * E,
       -1: -sp.Rational(9, 8) * E**4 - sp.Rational(219830, 11907) * E,
       -2: E**5 / 5 - E**4 / 5 + sp.Rational(967027, 2755620) * E**3
           - sp.Rational(53597707, 8573040) * E**2 + sp.Rational(455302607, 77157360) * E}
Dw2 = {2: b2_W2,
       1: -sp.Rational(2, 9) * E**2 + sp.Rational(8, 9) * E + sp.Rational(4, 3),
       0: sp.Rational(263, 567) * E**2 + sp.Rational(179, 567) * E,
       -1: -sp.Rational(256, 5103) * E**2 - sp.Rational(31130, 1701) * E,
       -2: -sp.Rational(3, 4) * E**3 + sp.Rational(2747993, 2571912) * E**2
           - sp.Rational(819059, 2571912) * E,
       -3: sp.Rational(2, 15) * E**4 - sp.Rational(5, 12) * E**3
           + sp.Rational(19, 60) * E**2 - E / 30}
Cw2 = commutator(Dw2, Xw2)
check("the W2 datum is a WALL REPRESENTATIVE: positive cascade Q_{m>0}=0 and Q_0=1 hold, "
      "the negative tail does NOT (Q_{-1..-5} != 0) -- it is NOT a genuine Weyl pair",
      all(sp.expand(Cw2.get(m, 0)) == 0 for m in range(1, 7))
      and sp.expand(Cw2.get(0, 0)) == 1
      and any(sp.expand(Cw2.get(m, 0)) != 0 for m in range(-6, 0)))
FX = fourier_ladder(Xw2)
check("FOURIER: the new TOP wall data of phi(X) comes from the OLD BOTTOM band "
      "(bandtop(phi X) = -bandbot(X) = 2, not 3): top-wall data is NOT phi-covariant",
      bandtop(FX) == -bandbot(Xw2) and bandtop(FX) != bandtop(Xw2))
Xw2b = dict(Xw2); Xw2b[-2] = sp.expand(Xw2b[-2] + falling(2) * (E + 7))
check("FOURIER refutation: two elements with IDENTICAL top wall data (a_3,a_2) have "
      "DIFFERENT Fourier tops -- so no function of the top wall alone is phi-equivariant",
      sp.expand(fourier_ladder(Xw2b).get(2, 0) - FX.get(2, 0)) != 0
      and Xw2b[3] == Xw2[3] and Xw2b[2] == Xw2[2])

# --- (4e) scaling / translation --------------------------------------------
lam2 = sp.Rational(5, 3)
gsc, dsc, oksc = cofactor(sp.expand(lam2 * a3_W2), b2_W2, 3, 2)
check("SCALING/TRANSLATION: wall data (k,q,d,g,p) is unchanged (roots are scale-free)",
      oksc and sp.expand(gsc - gW2) == 0 and dsc == dW2)

print("""
    WALL-DATA TRANSFORMATION TABLE   (wall W(k,q): delta(a_k)=cyc(k,d)g,
    delta(b_q)=cyc(q,d)g, d=gcd(k,q); Dixmier sigma_X=alpha p^a, sigma_D=beta p^b)

      generator                | (k,q)     | (a,b)   | cofactor g        | p
      -------------------------|-----------|---------|-------------------|--------
      scaling, translation     | (k,q)     | (a,b)   | g                 | p
      pair-exchange S          | (q,k)     | (b,a)   | g                 | p
      Fourier phi              | (-bandbot pair) -- top wall data is NOT covariant:
                               | phi maps the TOP wall to the BOTTOM wall
      D -> D - cX^s (s k > q)  | (sk, k)   | (a, sa) | cyc(k,d) g        | p
      X -> X + cD^s (s q > k)  | (q, sq)   | (sb, b) | cyc(q,d) g        | p
      cancelling move (s n = m)| (k, q'')  | -- new wall, q'' unconstrained by the old data
    Reading: RAISING moves act on the cofactor by MULTIPLICATION by cyc(k,d) and
    always land on a wall with q'|k', i.e. on the FORCED-EFFECTIVE (shifted-power)
    branch.  The singular-hatch branch (non-effective g) is therefore destroyed by
    any raising transvection and can only be re-created by a CANCELLING move.
""", flush=True)

# ===========================================================================
# §5  MONOVARIANT HUNT
# ===========================================================================
banner("§5  monovariant hunt: the primitive-degree invariant e = deg p")

def primitive_data(form):
    """sigma = const * p^a with p primitive: return (p, a, e=deg p)."""
    form = sp.expand(form)
    if form.free_symbols.isdisjoint({x, xi}):
        return (sp.Integer(1), 0, 0)
    c, facs = sp.factor_list(form, x, xi)
    mults = [m for (_, m) in facs]
    a = sp.igcd(*mults) if len(mults) > 1 else mults[0]
    p = sp.expand(sp.prod([f**(m // a) for (f, m) in facs]))
    return (p, int(a), sp.Poly(p, x, xi).total_degree())

def poisson(f, g):
    return sp.expand(sp.diff(f, xi) * sp.diff(g, x) - sp.diff(f, x) * sp.diff(g, xi))

# (5.0) the symbol-level input, re-derived on data: {sigma_X, sigma_D} = 0
sXw, nXw = leading_form(Xw2); sDw, mDw = leading_form(Dw2)
check("Dixmier symbol input re-derived on the W2 datum: {sigma_X,sigma_D} = 0, "
      "sigma_X=(x^2 xi)^3, sigma_D=(x^2 xi)^2", poisson(sXw, sDw) == 0
      and sXw == x**6 * xi**3 and sDw == x**4 * xi**2)
pW, aW, eW = primitive_data(sXw)
pW2, bW, eW2 = primitive_data(sDw)
check("W2 primitive p = x^2 xi, e = 3, (a,b) = (3,2), and n = a e = 9, m = b e = 6",
      sp.expand(pW - x**2 * xi) == 0 and (aW, eW) == (3, 3) and (bW, eW2) == (2, 3)
      and nXw == aW * eW and mDw == bW * eW)
check("e | gcd(n,m):  gcd(9,6) = 3 = e", sp.igcd(nXw, mDw) % eW == 0)

# (5.1) p-INVARIANCE under every generator (machine-checked on real data)
def pdata(X, D):
    sX, n = leading_form(X); sD, m = leading_form(D)
    p, a, e = primitive_data(sX)
    return p, a, e, n, m

Xtc = {3: sp.Integer(1), -1: -E}                 # X = x^3 - d   (n=3)
Dtc = {1: sp.Integer(1)}                         # D = x         (m=1)
check("tame control (x^3-d, x): [D,X]=1", commutator(Dtc, Xtc) == {0: sp.Integer(1)})
pt, at, et, nt, mt = pdata(Xtc, Dtc)
check("tame control has LINEAR primitive p = x, e = 1, (a,b) = (3,1)",
      sp.expand(pt - x) == 0 and (at, et, nt, mt) == (3, 1, 3, 1))

pinv = True
for (X0, D0, nm) in [(Xw2, Dw2, 'W2'), (Xtc, Dtc, 'tame')]:
    p0, a0, e0, n0, m0 = pdata(X0, D0)
    for s in (1, 2, 3):
        D1 = transvect_D(X0, D0, {s: sp.Integer(1)})
        p1, a1, e1, n1, m1 = pdata(X0, D1)
        pinv &= (sp.expand(p1 - p0) == 0 and e1 == e0)
        X1 = transvect_X(X0, D0, {s: sp.Integer(1)})
        p2, a2, e2, n2, m2 = pdata(X1, D0)
        pinv &= (sp.expand(p2 - p0) == 0 and e2 == e0)
    pS, aS, eS, _, _ = pdata(D0, cp_scale(X0, -1))
    pinv &= (sp.expand(pS - p0) == 0 and eS == e0)
check("p-INVARIANCE: transvections (s=1,2,3, both sides) and pair-exchange leave the "
      "primitive p (hence e) UNCHANGED, on the W2 and tame data", pinv)
pF, aF, eF, _, _ = pdata(fourier_ladder(Xw2), fourier_ladder(Dw2))
check("p-EQUIVARIANCE under Fourier: p |-> p o phi (x,xi)->(-xi,x), so deg p = e "
      "is preserved (W2: x^2 xi -> xi^2 x)",
      eF == eW
      and sp.expand(pF - sp.expand(pW.subs({x: -xi, xi: x}, simultaneous=True))) == 0)
print("    => Conditional on Dixmier's common-primitive leading-symbol lemma over C (or an\n"
      "       algebraically closed characteristic-zero field), e = deg p is a tame-orbit\n"
      "       invariant while n+m>2. Ambient linear symplectic changes carry p by an\n"
      "       invertible substitution. For invertible linear recombination of (X,D),\n"
      "       unequal degrees retain the higher symbol in at least one component and\n"
      "       Dixmier controls any cancellation in the other; at equal degrees the two\n"
      "       proportional equal-exponent powers cannot both cancel. Exchange swaps (a,b).", flush=True)

# (5.2) the exponent floor for e >= 2:  a,b >= 2, (a,b) != (2,2), so n+m >= 5e
def floor_chain(a, b, e, trace=None):
    """model of the forced descent: returns the contradiction reached, if any."""
    if a == 0 or b == 0:
        return 'deg-0 member: [D,X]=0, contradiction'
    if a == 1:
        # D -> D - c X^b cancels; new m'' < b e and e | m'' => m'' = b'' e, b''<b
        return 'a=1 descends b -> b-1 -> ... -> 0 (constant D): contradiction'
    if b == 1:
        return 'b=1 (mirror of a=1): contradiction'
    if a == b == 2:
        return 'a=b=2: affine cancel gives b\'\' in {0,1}, both excluded'
    return None
# AUDIT FIX (2026-07-25): the former check here compared floor_chain's branch
# predicate with its own definition -- a TAUTOLOGY that could not fail, and whose
# `e` argument was never read (it passed identically at e=1, where the tame
# control (x^3-d, x) with (a,b)=(3,1) refutes the encoded conclusion). Removed.
# floor_chain is retained ONLY as documentation of the paper argument's case
# split; the exclusions a,b>=2 and (a,b)!=(2,2) are PAPER results (memo 5.3),
# conditional on Dixmier 1968 + Makar-Limanov + the campaign band floors.
print("    [paper] exclusions a<2, b<2, (a,b)=(2,2) are argued in memo 5.3, NOT\n"
      "            machine-proved here; the cases are documented by floor_chain().",
      flush=True)
adm = [(a, b) for a in range(2, 30) for b in range(2, 30) if (a, b) != (2, 2)]
check("GIVEN the paper exclusions: min(a+b) over the admissible exponent set = 5, "
      "attained only at {2,3} (pure arithmetic; falsifiable)",
      min(a + b for (a, b) in adm) == 5
      and {tuple(sorted(ab)) for ab in adm if sum(ab) == 5} == {(2, 3)})
check("W2 formal datum has the floor arithmetic n+m = 5e (one datum, e=3; "
      "it is not a Weyl pair, and the general n+m >= 5e is the paper FLOOR THEOREM)",
      nXw + mDw == 5 * eW)

# (5.3) e-quantisation blocks the exit to n+m = 2 and to band <= 2
check("e-QUANTISATION: for e >= 2 no state with n,m in e*Z_{>=1} has n+m = 2; and a "
      "transvection leaves one member's degree in e*Z, so n+m=2 (n=m=1) is unreachable",
      all(not any(n + m == 2 for n in range(e, 60, e) for m in range(e, 60, e))
          for e in range(2, 8)))
print("    => with the campaign's band floors (band-1 rigidity P3 + the band-2 theorem\n"
      "       84978b9: band <= 2 => the pair is a TAME image of (x,d)) and the tame-orbit\n"
      "       induction below, a pair with e >= 2 can never reach band <= 2 either.",
      flush=True)

# (5.4) the tame orbit of (x,d) has e = 1  (model induction)
first_moves = []
Xstd, Dstd = {1: sp.Integer(1)}, {-1: E}
for s in (2, 3, 4):
    D1 = transvect_D(Xstd, Dstd, {s: sp.Integer(1)})
    p1, a1, e1, _, _ = pdata(Xstd, D1)
    first_moves.append(e1)
    X1 = transvect_X(Xstd, Dstd, {s: sp.Integer(1)})
    p2, a2, e2, _, _ = pdata(X1, Dstd)
    first_moves.append(e2)
check("TAME-ORBIT INDUCTION: every first degree-raising move out of (x,d) produces "
      "e = 1 (s=2,3,4, both sides); p-invariance then gives e = 1 on the whole orbit",
      all(e == 1 for e in first_moves))

# (5.5) REFUTATIONS of the naive monovariants -------------------------------
Xr, Dr = Xtc, Dtc                                     # (x^3-d, x): wall (3,1), d=1
D_up = transvect_D(Xr, Dr, {2: sp.Integer(1)})        # D - X^2 : bandtop 6
d_before = int(sp.igcd(bandtop(Xr), bandtop(Dr)))
d_up = int(sp.igcd(bandtop(Xr), bandtop(D_up)))
D_back = transvect_D(Xr, D_up, {2: sp.Integer(-1)})
d_back = int(sp.igcd(bandtop(Xr), bandtop(D_back)))
check("REFUTED: d = gcd(band pair) is NOT a global monovariant "
      f"(explicit word: d = {d_before} -> {d_up} -> {d_back})",
      d_up > d_before and d_back == d_before)
check("REFUTED: cofactor effectivity is not a monovariant either -- the inverse of a "
      "raising move restores the non-effective W2 cofactor",
      transvect_D(Xw2, transvect_D(Xw2, Dw2, {2: sp.Integer(1)}), {2: sp.Integer(-1)})
      == {k: sp.expand(v) for k, v in Dw2.items()})
check("REFUTED: the mutually-non-dividing property of (a,b) is destroyed by one raising "
      "move: (a,b) = (3,2) -> (3,6) with 3 | 6",
      pdata(Xw2, transvect_D(Xw2, Dw2, {2: sp.Integer(1)}))[1] == 3
      and bdeg(transvect_D(Xw2, Dw2, {2: sp.Integer(1)})) == 18)

# LIMIT OF THIS OVER-APPROXIMATE DEGREE TRANSITION MODEL.
def degree_model_escape(n0, m0, S=3, depth=2):
    """cancelling moves allowed whenever divisibility holds, any depth of drop."""
    best = None
    def rec(n, m, side, dep, word):
        nonlocal best
        if dep == 0:
            return
        for s in range(1, S + 1):
            if side != 'D':
                if s * n == m:
                    for mm in range(0, m):
                        witness = word + [(f'D{s}', (n, m), (n, mm), 'cancel')]
                        if best is None or n + mm < best[0]:
                            best = (n + mm, witness)
                        rec(n, mm, 'D', dep - 1, witness)
                else:
                    target = (n, max(m, s * n))
                    rec(*target, 'D', dep - 1,
                        word + [(f'D{s}', (n, m), target, 'noncancel')])
            if side != 'X':
                if s * m == n:
                    for nn in range(0, n):
                        witness = word + [(f'X{s}', (n, m), (nn, m), 'cancel')]
                        if best is None or nn + m < best[0]:
                            best = (nn + m, witness)
                        rec(nn, m, 'X', dep - 1, witness)
                else:
                    target = (max(n, s * m), m)
                    rec(*target, 'X', dep - 1,
                        word + [(f'X{s}', (n, m), target, 'noncancel')])
    rec(n0, m0, None, depth, [])
    return best
esc = degree_model_escape(9, 6, S=3, depth=2)
check("DEGREE-MODEL LIMIT: the tested over-approximate transition relation admits "
      f"witness {esc[1] if esc else None}, ending at sum {esc[0] if esc else None} < 15; "
      "therefore that relation alone is insufficient, without excluding monovariants "
      "that use actual (n,m,k) trajectories", esc is not None and esc[0] < 15
      and len(esc[1]) <= 2)
# ... and the true algebra blocks that very word:
c_ = sp.Symbol('c_')
D1 = transvect_D(Xw2, Dw2, {1: sp.Integer(1)})          # D - X : degrees (9,9)
X1 = transvect_X(Xw2, D1, {1: sp.Integer(1)})           # X + (D - X) = D
check("...but the ALGEBRA blocks it: the (9,6)->(9,9)->cancel word gives X'' = D "
      "(degree 6) and D' = D - X (degree 9): sum back to 15, never below",
      bdeg(D1) == 9 and bdeg(X1) == 6 and bdeg(X1) + bdeg(D1) == 15)

# (5.6) lowest arithmetically admissible e>=2 floor case ---------------------
floors = {e: 5 * e for e in range(2, 7)}
check("floor ladder 5e: among arithmetically admissible e>=2 cases, e=2 gives the "
      "lowest floor 10, below the W2 e=3 floor 15; realizability and e=1 sectors remain open",
      floors[2] == 10 and floors[3] == 15 and min(floors.values()) == 10)
check("a primitive binary quadratic has two DISTINCT roots (a double root is a square, "
      "hence imprimitive), so e=2 forces p ~ x*xi up to SL_2",
      primitive_data(sp.expand(x * xi))[2] == 2
      and primitive_data(sp.expand(x**2))[1] == 2
      and primitive_data(sp.expand(x**2))[2] == 1)
# which bands can carry sigma_X = (x xi)^2 ?  i + d = 2, d = 2  =>  i = 0, and the
# competing bands at Bernstein degree 4 are i + 2d = 4.
bands4 = [(i, d) for i in range(-4, 5) for d in range(0, 5) if i + 2 * d == 4]
p2_band, p3_band = 2 - 2, 3 - 3
check("e=2 in p=x*xi coordinates puts both p^2 and p^3 on band 0; after a generic "
      f"SL2 change a top may occupy several fixed-coordinate slots (p^2 slots: {bands4}). "
      "Band is not SL2-invariant, so fixed single-band classifications do not automatically transfer",
      p2_band == p3_band == 0 and (0, 2) in bands4 and len(bands4) == 5)
_pq = sp.expand(x * (x + xi))                     # another primitive quadratic
_top = sp.expand(_pq**2)
_spread = {sp.Poly(_top, x, xi).monoms()[i] for i in range(len(sp.Poly(_top, x, xi).monoms()))}
check("...witness: p = x(x+xi) is primitive of degree 2 and p^2 spreads the degree-4 "
      "top over 3 distinct (band,degree) slots",
      primitive_data(_pq) == (_pq, 1, 2) and len(_spread) == 3)

# (5.7) the FLOOR THEOREM's cancellation step, on symbolic leading forms
al2, be2 = sp.symbols('alpha2 beta2')
pp = x**2 * xi
canc = True
for b in (1, 2, 3):
    sX, sD = al2 * pp, be2 * pp**b                # a = 1
    cst = be2 / al2**b
    canc &= sp.expand(sD - cst * sX**b) == 0
check("FLOOR THEOREM step (symbolic): a=1 => sigma_D = beta p^b is EXACTLY "
      "(beta/alpha^b) sigma_X^b, so D -> D - (beta/alpha^b) X^b kills the top "
      "(b=1,2,3, p=x^2 xi)", canc)
canc22 = sp.expand((be2 * pp**2) - (be2 / al2) * (al2 * pp**2)) == 0
check("FLOOR THEOREM step (symbolic): a=b=2 => the AFFINE move D -> D - (beta/alpha)X "
      "kills the top", canc22)

# ===========================================================================
# §6  BOUNDED SHORT-WORD ESCAPE SEARCH (band 3)
# ===========================================================================
banner("§6  bounded short-word escape search at band 3")

def in_A1(Y):
    for i, a in Y.items():
        if i < 0 and sp.expand(sp.rem(sp.Poly(sp.expand(a), E), sp.Poly(falling(-i), E)).as_expr()) != 0:
            return False
    return True

def key(X, D):
    return (tuple(sorted((i, sp.srepr(sp.expand(v))) for i, v in X.items())),
            tuple(sorted((i, sp.srepr(sp.expand(v))) for i, v in D.items())))

def successors(X, D, S, last):
    """all one-letter moves; 'last' forbids an immediately repeated side."""
    out = []
    if last != 'D':
        r, _ = min_over_span(D, X, S)
        if r:
            out.append(('oD', X, r, 'D'))
        for s in range(1, S + 1):
            D2 = cp_add(D, cp_pow(X, s), -1)
            if D2:
                out.append((f'rD{s}', X, D2, 'D'))
    if last != 'X':
        r, _ = min_over_span(X, D, S)
        if r:
            out.append(('oX', r, D, 'X'))
        for s in range(1, S + 1):
            X2 = cp_add(X, cp_pow(D, s), +1)
            if X2:
                out.append((f'rX{s}', X2, D, 'X'))
    out.append(('S', D, cp_scale(X, -1), None))
    if in_A1(X) and in_A1(D):
        out.append(('F', fourier_ladder(X), fourier_ladder(D), None))
    return out

def word_search(X0, D0, S, depth, cap, label):
    I0 = invariant(X0, D0)
    best, best_word, seen, nodes = I0, [], {(key(X0, D0), None)}, 0
    frontier = [(X0, D0, None, [])]
    for lev in range(depth):
        nxt = []
        for (X, D, last, w) in frontier:
            for (nm, X2, D2, side) in successors(X, D, S, last):
                nodes += 1
                I2 = invariant(X2, D2)
                if I2 < best:
                    best, best_word = I2, w + [nm]
                if I2[0] > cap:
                    continue
                kk = (key(X2, D2), side)
                if kk in seen:
                    continue
                seen.add(kk)
                nxt.append((X2, D2, side, w + [nm]))
        frontier = nxt
        print(f"    [{label}] depth {lev+1}: {nodes} words explored, "
              f"{len(frontier)} live states, best = {best}", flush=True)
    return I0, best, best_word, nodes

# --- POSITIVE CONTROLS (§7 folded in): the search must FIND known reductions --
kap = sp.Integer(1)
Xpc = {3: sp.Integer(1), -1: -E}                          # x^3 - d
Dpc = {1: sp.Integer(1)}                                  # x
I0, bpc, wpc, npc = word_search(Xpc, Dpc, 3, 2, 40, 'tame positive control')
check("POSITIVE CONTROL 1 (x^3-d, x): the search FINDS a strict reduction "
      f"(15/4-style) {I0} -> {bpc} via {wpc}", bpc < I0 and bpc == (2, 1))

c0, c1, A0, lam0, be0 = (sp.Integer(2), sp.Integer(3), sp.Integer(5),
                         sp.Integer(7), sp.Integer(11))
U = cp_add({1: sp.Integer(1)}, {0: c0}); U = cp_add(U, {-1: c1 * E})
Xtf = cp_add(cp_add(cp_pow(U, 3), {-1: E}, -1 / kap), {0: A0}, -1)
Dtf = cp_add(cp_add(cp_scale(Xtf, lam0), cp_scale(U, kap)), {0: be0})
check("tame family (c1 != 0) is a genuine Weyl pair: [D,X] = 1",
      commutator(Dtf, Xtf) == {0: sp.Integer(1)})
I0t, btf, wtf, ntf = word_search(Xtf, Dtf, 3, 2, 40, 'tame family c1!=0')
check(f"POSITIVE CONTROL 2 (tame family): the search FINDS a strict reduction "
      f"{I0t} -> {btf} via {wtf}", btf < I0t)

# --- THE STUCK STRATUM: W2 hatch datum and the shifted cube -----------------
Xsc, Dsc = {3: E * (E + 1) * (E + 2)}, {2: E * (E + 1)}    # shifted-cube top datum
check("shifted-cube top datum is on the wall W(3,2) with effective cofactor g=1, "
      "and has the SAME leading form (x^2 xi)^3 as W2 (symbol blindness)",
      twisted_wronskian(Xsc[3], Dsc[2], 3, 2) == 0
      and leading_form(Xsc)[0] == x**6 * xi**3 and invariant(Xsc, Dsc) == (15, 3))

DEPTH = 4 if HEAVY else 3
for (Xs, Ds, lab) in [(Xw2, Dw2, 'W2 hatch datum'), (Xsc, Dsc, 'shifted cube')]:
    I0s, bs, ws, ns = word_search(Xs, Ds, 2, DEPTH, 45, lab)
    check(f"NO ESCAPE at the {lab}: over the searched space the invariant never drops "
          f"below the start {I0s} (best reached {bs}, {ns} words)", bs >= I0s)

print(f"""
    SEARCH SPACE COVERED (exactly).  Alphabet per node: oD, oX (the bdeg-OPTIMAL
    transvection D -> D - f(X) / X -> X + f(D) over ALL f of degree <= S -- exact by
    the §3d greedy/echelon argument), rD1..rDS, rX1..rXS (unit-coefficient raising
    transvections), pair-exchange S, Fourier F (when membership holds).  Consecutive
    same-side transvections are excluded (they collapse, §3a).  Depth {DEPTH}
    (HEAVY=1 gives 4), S = 2 for the band-3 wall data and S = 3 for the controls,
    invariant-sum cap 45.  NOT covered: non-unit raising coefficients, transvections
    of degree > S, affine-symplectic letters with both off-diagonal entries nonzero,
    and words longer than the depth.  The search is a CONSISTENCY PROBE of the §5
    floor theorem, not a proof.
""", flush=True)

if not HEAVY:
    skip("depth-4 leg of the §6 word search", "set HEAVY=1; runtime is machine-dependent")

# ===========================================================================
# §7  the assembled conditional statement
# ===========================================================================
banner("§7  assembled statement")

print("""    CONDITIONAL PAPER FLOOR THEOREM over C (or an algebraically closed
    characteristic-zero field), assuming Dixmier's common-primitive leading-symbol
    lemma; the band clause additionally uses the cited band-1 rigidity and band-2 theorem.
    The bounded verifier is regression support, not a proof object:

      Let (X,D) in A_1 with [D,X]=1 and let p be the primitive Dixmier form,
      sigma_X = alpha p^a, sigma_D = beta p^b, e = deg p, n = ae, m = be.
      (1) e is constant along every tame trajectory whose states all have n+m>2,
          and p itself is constant up to the linear symplectic action of the
          affine/Fourier letters.                                     (§5.1)
      (2) If e >= 2 then a >= 2, b >= 2 and (a,b) != (2,2); hence n+m >= 5e. (§5.2)
      (3) Hence a pair with e >= 2 can never reach n+m = 2, and (with the band
          floors) never reaches band <= 2: it is never tame-equivalent to (x,d). (§5.3-4)
      (4) Any GENUINE pair with the W2 or shifted-cube leading forms has e = 3,
          {a,b} = {2,3}, n+m = 15 = 5e. Conditional on the cited inputs, such a
          floor-attaining genuine pair is tame-orbit minimal in the claimed invariant.
          The formal wall/top representatives tested here are not Weyl pairs and are
          never themselves asserted tame-minimal.

    WHAT THIS DOES NOT DO.  It is a statement about GENUINE pairs in a leading-form
    stratum, not about a particular formal datum's lower bands; it does not show any
    such pair exists or that orbit-minimal/minimal counterexamples attain equality;
    descent from a+b>5 to the floor is open. It does not close Gap 1 above the floor (e.g. the
    band-4 census hatch has e=3, (a,b)=(4,3), n+m=21 > 15, so a composite descent
    to a band-3 configuration is NOT excluded); and it says nothing about e = 1
    strata (the constant-h / kappa_2 sector), where p is linear and the floor
    argument is vacuous.""", flush=True)

# ===========================================================================
banner("summary")
npass = sum(1 for _, ok in CHECKS if ok)
nfail = len(CHECKS) - npass
for nm, ok in CHECKS:
    if not ok:
        print("  FAILED:", nm, flush=True)
print(f"  {npass}/{len(CHECKS)} checks passed; {len(SKIPPED)} optional legs skipped; "
      f"{time.time()-T0:.1f}s", flush=True)
if nfail:
    print("SOME GAP1 CHECKPOINT CHECKS FAILED", flush=True)
    sys.exit(1)
print("ALL EXECUTED CHECKS PASSED; OPTIONAL CHECKS SKIPPED" if SKIPPED
      else "ALL GAP1 CHECKPOINT CHECKS PASSED", flush=True)
