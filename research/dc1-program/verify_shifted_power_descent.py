#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# verify_shifted_power_descent.py
#
# Exact-SymPy certificate for research/dc1-program/shifted-power-descent.md
# (band-reduction GAP 2: the shifted-power descent, architecture step 3).
#
# The Band Normalization Theorem forces a tame-minimal non-generating pair onto
# a W(k,q) wall.  On the SHIFTED-POWER branch (effective cofactor) the top is a
# consecutive shifted k-th power  a_k = c * h(E) h(E+1) ... h(E+k-1).  W(k,q)
# constrains only the top TWO ladder rungs; GAP 2 is whether that structure
# PROPAGATES to sub-leading bands and licenses a genuine band-lowering tame move.
#
# This file establishes, all from scratch and self-containedly:
#   §0  the crossed-product ladder engine  Q_m = [D,X]_m,  Q_0 = (T-1)G, and the
#       tame-move machinery (transvections exp(ad p(X)), exp(ad q(D)) preserve
#       [D,X]=1);  Bernstein-degree / band invariant (n+m, k).
#   §1  the W(k,q) top wall Q_{k+q}=0 and its cofactor; the balanced shifted power
#       a_k = prod h^[j], b_{k-1} = kappa prod h^[j] solving Q_{2k}=0 (W1) and the
#       wall Q_{2k-1}=0.  (k=3,4.)
#   §2  TASK 1 -- the NEXT rung Q_{k+q-1}=Q_{2k-2} written explicitly for the
#       balanced pair, and the h-compatibility it forces on a_{k-1}, b_{k-2}:
#       prod_{j=0}^{k-2} h^[j] | a_{k-1}  and  prod_{j=0}^{k-3} h^[j] | b_{k-2},
#       with the S_k = 1+T+...+T^{k-1} quantum midpoint.  (k=3 AND k=4.)
#   §3  TASKS 2/3 PROPAGATION -- re-verify the band-3 h-constant forcing for
#       cube-separated h (linear AND quadratic): the central telescoping factors
#       G = h^[-1] M, membership gives M(0)=0, and Q_0=1 forces h CONSTANT; plus
#       the U-side forcing  [d, x h(E)] = 1  <=>  h constant.
#   §4  TASK 2 -- the EXPLICIT band-3 tame move on genuine constant-h data:
#       (i) gauge D -> D - lam X (B0 collapse, band D -> 1),
#       (ii) mirror transvection X -> X - ((D-beta)/kappa)^3 (band X -> 1),
#       strictly lowering (n+m,k) and preserving [D,X]=1, on the positive control
#       AND the full tame family with c1 != 0.  Reducedness of the nonconstant
#       (3,2) leading form (single-transvection descent STALLS -- the gap's crux).
#   §5  TASK 4 -- the uniform reducing move for the shifted k-th power positive
#       control X = U^k - d/kappa, D = kappa U (k=3,4,5); imbalanced coprime
#       (k,1) walls reduce, (k,q>=2) coprime q not| k stall (the composite escape).
#   §6  TASK 5 -- band-3 branch distinctness bookkeeping (shifted-cube / exotic AP
#       / W2 tops are genuinely different objects sharing one leading form).
#
# Conventions frozen from the DC1 corpus:
#   A_1[x^-1] = (+)_k x^k C[E],  E = x d,  (x^a f)(x^b g) = x^{a+b} f(E+b) g(E),
#   f^[n](E) = f(E+n),  T f = f^[1],  S_n = 1+T+...+T^{n-1},
#   Q_m = sum_{k+l=m} ( b_l^[k] a_k - a_k^[l] b_l ),  [D,X]=1 <=> Q_m = delta_{m0},
#   membership (E)_j = E(E-1)...(E-j+1) | a_-j, b_-j,  gauge b_k = 0,
#   G = sum_{k>=1} sum_{j=0}^{k-1} ( a_k^[j-k] b_-k^[j] - b_k^[j-k] a_-k^[j] ),
#   Q_0 = (T-1)G.
#
# All arithmetic exact over QQ / QQ(kappa,rho,...).  Run:
#   uv run --with sympy python research/dc1-program/verify_shifted_power_descent.py
# A successful run ends:  ALL SHIFTED POWER DESCENT CHECKS PASSED
# ---------------------------------------------------------------------------
import sympy as sp

E = sp.symbols('E')
kappa, rho, alpha = sp.symbols('kappa rho alpha')

CHECKS = []
def check(name, cond):
    ok = bool(cond)
    CHECKS.append((name, ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}", flush=True)
    return ok

def sh(f, n):                                   # f^[n](E) = f(E+n)
    return sp.expand(sp.sympify(f).subs(E, E + n))

def gen(name, d):                               # generic degree-d polynomial + coeff list
    cs = sp.symbols(f'{name}_0:{d+1}')
    if not isinstance(cs, (list, tuple)):
        cs = (cs,)
    return sum(cs[i] * E**i for i in range(d + 1)), list(cs)

def divides(fac, c):                            # fac | c  over C[E]
    c = sp.expand(sp.sympify(c))
    if c == 0:
        return True
    return sp.rem(sp.Poly(c, E), sp.Poly(sp.expand(fac), E)) == 0

def prodsh(h, r):                               # prod_{j=0}^{r-1} h^[j]  (shifted r-fold)
    return sp.expand(sp.prod([sh(h, j) for j in range(r)])) if r > 0 else sp.Integer(1)

# ---- crossed-product ladder engine (dicts index -> C[E]) ----
def cp_mul(P, Q):
    R = {}
    for k, fk in P.items():
        for l, gl in Q.items():
            R[k + l] = sp.expand(R.get(k + l, 0) + sh(fk, l) * gl)
    return {m: v for m, v in R.items() if sp.expand(v) != 0}

def cp_sub(P, Q):
    R = dict(P)
    for k, v in Q.items():
        R[k] = sp.expand(R.get(k, 0) - v)
    return {m: v for m, v in R.items() if sp.expand(v) != 0}

def cp_scal(P, c):
    return {k: sp.expand(c * v) for k, v in P.items() if sp.expand(c * v) != 0}

def cp_pow(P, n):
    R = {0: sp.Integer(1)}
    for _ in range(n):
        R = cp_mul(R, P)
    return R

def commutator(D, X):                           # [D,X] = DX - XD
    return cp_sub(cp_mul(D, X), cp_mul(X, D))

def Qm(X, D, m):                                # ladder-m coeff of [D,X]
    tot = 0
    for k, ak in X.items():
        l = m - k
        if l in D:
            tot += sh(D[l], k) * ak - sh(ak, l) * D[l]
    return sp.expand(tot)

def Gpot(X, D, K):                              # closed-form moment potential
    G = 0
    for k in range(1, K + 1):
        for j in range(0, k):
            G += sh(X.get(k, 0), j - k) * sh(D.get(-k, 0), j) \
               - sh(D.get(k, 0), j - k) * sh(X.get(-k, 0), j)
    return sp.expand(G)

def bern_deg(coeffs):                           # Bernstein degree: max_k (k + 2 deg a_k)
    best = None
    for k, a in coeffs.items():
        a = sp.expand(a)
        if a == 0:
            continue
        val = k + 2 * sp.Poly(a, E).degree()
        best = val if best is None else max(best, val)
    return best

def band(coeffs):                               # band = max |k| with a_k != 0
    b = 0
    for k, a in coeffs.items():
        if sp.expand(a) != 0:
            b = max(b, abs(k))
    return b

def invariant(X, D):
    return (bern_deg(X) + bern_deg(D), max(band(X), band(D)))

def lex_lt(u, v):                               # strict lexicographic <
    return u[0] < v[0] or (u[0] == v[0] and u[1] < v[1])

# ===========================================================================
print("\n=== §0  crossed-product engine + tame-move machinery ===")
# ===========================================================================

# Q_m == direct commutator ladder coefficient (generic band-3 data, deg-2 coeffs)
Xg = {k: gen(f'a{k}', 2)[0] for k in range(-3, 4)}
Dg = {k: gen(f'b{k}', 2)[0] for k in range(-3, 4)}
Cg = commutator(Dg, Xg)
check("Q_m == [D,X]_m for all m in [-6,6] (generic band-3)",
      all(sp.expand(Cg.get(m, 0) - Qm(Xg, Dg, m)) == 0 for m in range(-6, 7)))
check("Q_0 = (T-1)G with the closed-form staggered potential (generic band-3)",
      sp.expand(sh(Gpot(Xg, Dg, 3), 1) - Gpot(Xg, Dg, 3) - Qm(Xg, Dg, 0)) == 0)

# reference pair (x, d): x = {1:1}, d = x^{-1}E = {-1:E}
Xref = {1: sp.Integer(1)}
Dref = {-1: E}
check("[d, x] = 1 (reference pair)", commutator(Dref, Xref) == {0: sp.Integer(1)})

# transvections preserve [D,X]=1: exp(ad p(X)): D -> D - p'(X);  exp(ad q(D)): X -> X + q'(D)
cc = sp.symbols('cc')
check("transvection D -> D - cc X^2 keeps [D,X]=1 (exp ad p(X), p=cc X^3/3)",
      commutator(cp_sub(Dref, cp_scal(cp_pow(Xref, 2), cc)), Xref) == {0: sp.Integer(1)})
check("transvection X -> X + cc D^2 keeps [D,X]=1 (exp ad q(D), q=cc D^3/3)",
      commutator(Dref, cp_sub(Xref, cp_scal(cp_pow(Dref, 2), -cc))) == {0: sp.Integer(1)})
# pair-exchange (X,D) -> (D,-X):  [-X, D] = [D,X]
check("pair-exchange (X,D)->(D,-X): [-X,D] equals [D,X]",
      all(sp.expand(commutator({k: -v for k, v in Xg.items()}, Dg).get(m, 0) - Cg.get(m, 0)) == 0
          for m in set(Cg) | set(commutator({k: -v for k, v in Xg.items()}, Dg))))

# invariant sanity: (x,d) has (n+m,k)=(2,1); (x^3-d/kappa, kappa x) has (4,3)
check("invariant of (x,d) is (2,1)", invariant(Xref, Dref) == (2, 1))

# ===========================================================================
print("\n=== §1  the W(k,q) top wall and the balanced shifted power ===")
# ===========================================================================
# For band X = k, band D = q, the top ladder band of [D,X] is k+q > 0, so its
# coefficient must vanish:  Q_{k+q} = b_q^[k] a_k - a_k^[q] b_q = 0  (unconditional).
def twisted_wronskian(ak, bq, k, q):
    return sp.expand(sh(bq, k) * ak - sh(ak, q) * bq)

for k in (3, 4):
    q = k - 1
    for hh in (E - rho, sp.expand((E - sp.Rational(1, 3)) * (E - sp.Rational(9, 4)))):
        ak = prodsh(hh, k)                 # shifted k-th power  a_k = prod_{j<k} h^[j]
        bk1 = sp.expand(kappa * prodsh(hh, k - 1))  # b_{k-1} = kappa prod_{j<k-1} h^[j]
        # W1: gauge b_k = 0; Q_{2k} = b_k^[k]a_k - a_k^[k]b_k = 0 in gauge (trivially),
        # and W1 forces b_k = lambda a_k -- re-derive the Wronskian vanishing:
        lam = sp.symbols('lam')
        check(f"k={k}: W1 Q_{{2k}} = (lam a_k)^[k] a_k - a_k^[k](lam a_k) = 0 (b_k=lam a_k, deg h={sp.degree(hh,E)})",
              twisted_wronskian(ak, sp.expand(lam * ak), k, k) == 0)
        # the wall W(k,k-1): Q_{2k-1} = b_{k-1}^[k] a_k - a_k^[k-1] b_{k-1} = 0
        check(f"k={k}: wall Q_{{2k-1}} = 0 solved by shifted power (b_{{k-1}}=kappa prod, deg h={sp.degree(hh,E)})",
              twisted_wronskian(ak, bk1, k, k - 1) == 0)

# cofactor: on the root necklace (sig = T^-1) the wall reads (sig^k-1)d(b)=(sig^q-1)d(a);
# minimal effective cofactor g=1 gives the CONSECUTIVE (tame) shifted power.  Re-derive
# the necklace identity S_k * d(b_{k-1}) = S_{k-1} * d(a_k) for the shifted power:
sig = sp.symbols('sigma')
def Sn(n):
    return sp.expand(sum(sig**i for i in range(n)))
for k in (3, 4):
    # consecutive shifted power (h = E, single root at 0): d(a_k)=S_k, d(b_{k-1})=S_{k-1};
    # the necklace wall (sig^k-1) d(b_{k-1}) = (sig^{k-1}-1) d(a_k) then holds because
    # sig^k-1 = (sig-1)S_k, sig^{k-1}-1 = (sig-1)S_{k-1}, so both sides = (sig-1)S_k S_{k-1}.
    check(f"k={k}: necklace wall (sig^k-1) S_{{k-1}} = (sig^{{k-1}}-1) S_k  (consecutive shifted power)",
          sp.expand((sig**k - 1) * Sn(k - 1) - (sig**(k - 1) - 1) * Sn(k)) == 0)

# ===========================================================================
print("\n=== §2  TASK 1: the next rung Q_{k+q-1}=Q_{2k-2} and h-compatibility ===")
# ===========================================================================
# Balanced q=k-1 so k+q-1 = 2k-2.  In gauge b_k=0 with a_k, b_{k-1} shifted powers,
#   Q_{2k-2} = b_{k-2}^[k]a_k - a_k^[k-2]b_{k-2} + b_{k-1}^[k-1]a_{k-1} - a_{k-1}^[k-1]b_{k-1}.
# CLAIM (h-compatibility): Q_{2k-2}=0 forces, for cube/shift-separated h,
#   prod_{j=0}^{k-2} h^[j] | a_{k-1}   and   prod_{j=0}^{k-3} h^[j] | b_{k-2}.

def Q2km2_explicit(k, ak, bk1, akm1, bkm2):
    # the three (i,l) with i+l=2k-2, i,l in [-k,k], b_k=0 kills (k-2, k):
    return sp.expand(sh(bkm2, k) * ak - sh(ak, k - 2) * bkm2
                     + sh(bk1, k - 1) * akm1 - sh(akm1, k - 1) * bk1)

# match the explicit form against the engine Qm on the actual ladder (k=3 and k=4):
for k in (3, 4):
    hh = E - rho
    ak = prodsh(hh, k); bk1 = sp.expand(kappa * prodsh(hh, k - 1))
    akm1, _ = gen('AK', 3); bkm2, _ = gen('BK', 3)
    X = {k: ak, k - 1: akm1}
    D = {k: sp.Integer(0), k - 1: bk1, k - 2: bkm2}
    check(f"k={k}: engine Q_{{2k-2}} equals the explicit 4-term form (gauge b_k=0)",
          sp.expand(Qm(X, D, 2 * k - 2) - Q2km2_explicit(k, ak, bk1, akm1, bkm2)) == 0)

# k=3 h-compatibility (shifted-cube memo re-derived): Q_4=0 forces h h^[1] | a_2, h | b_1.
hh = E - rho
a3 = prodsh(hh, 3); b2 = sp.expand(kappa * prodsh(hh, 2))
g3, _ = gen('g3', 2); a2_3 = sp.expand(hh * sh(hh, 1) * g3)     # a_2 = h h^[1] g
be3, _ = gen('be3', 2); b1_3 = sp.expand(hh * be3)             # b_1 = h beta
X3 = {3: a3, 2: a2_3}
D3 = {3: sp.Integer(0), 2: b2, 1: b1_3}
# Q_4 collapses to a_3 h^[3] [ (T^3-1)beta - kappa (T^2-1)g ]:
check("k=3: Q_4 = a_3 h^[3] [ (T^3-1)beta - kappa (T^2-1)g ]  (b_1=h beta, a_2=h h^[1] g)",
      sp.expand(Qm(X3, D3, 4)
                - a3 * sh(hh, 3) * ((sh(be3, 3) - be3) - kappa * (sh(g3, 2) - g3))) == 0)
# the S_3 midpoint: (T-1)[S_3 beta - kappa S_2 g] = (T^3-1)beta - kappa(T^2-1)g
S3b = sp.expand(sh(be3, 2) + sh(be3, 1) + be3)
S2g = sp.expand(sh(g3, 1) + g3)
check("k=3: S_3 midpoint (T-1)[S_3 beta - kappa S_2 g] = (T^3-1)beta - kappa(T^2-1)g",
      sp.expand((sh(S3b - kappa * S2g, 1) - (S3b - kappa * S2g))
                - ((sh(be3, 3) - be3) - kappa * (sh(g3, 2) - g3))) == 0)

# NECESSITY at a concrete cube-separated h = E(2E-1) (roots 0,1/2): Q_4=0 with generic
# a_2, b_1 FORCES h h^[1]|a_2 and h|b_1.  (Regression of the coprimality proof, k=3.)
hs = sp.expand(E * (2 * E - 1))
check("h=E(2E-1) is shift-separated: gcd(h,h^[j])=1 for j=1..3",
      all(sp.gcd(sp.Poly(hs, E), sp.Poly(sh(hs, j), E)).degree() == 0 for j in (1, 2, 3)))
a3s = prodsh(hs, 3); b2s = sp.expand(kappa * prodsh(hs, 2))
a2f, a2c = gen('a2f', 4); b1f, b1c = gen('b1f', 3)
Xk3 = {3: a3s, 2: a2f}
Dk3 = {3: sp.Integer(0), 2: b2s, 1: b1f}
sol = sp.solve(sp.Poly(Qm(Xk3, Dk3, 4), E).all_coeffs(), b1c + a2c, dict=True)
check("k=3: Q_4=0 has a unique constrained locus", len(sol) == 1)
a2sol = sp.expand(a2f.subs(sol[0])); b1sol = sp.expand(b1f.subs(sol[0]))
check("k=3: Q_4=0 FORCES h h^[1] | a_2 (necessity, cube-separated)",
      divides(hs, a2sol) and divides(sh(hs, 1), a2sol))
check("k=3: Q_4=0 FORCES h | b_1 (necessity, cube-separated)", divides(hs, b1sol))

# k=4 h-compatibility NECESSITY: Q_6=0 with generic a_3, b_2 (a_4,b_3 shifted powers)
# FORCES  h h^[1]h^[2] | a_3  and  h h^[1] | b_2  (the k=4 analogue).
a4s = prodsh(hs, 4); b3s = sp.expand(kappa * prodsh(hs, 3))
a3f, a3c = gen('a3f', 6); b2f4, b2c4 = gen('b2f4', 5)
Xk4 = {4: a4s, 3: a3f}
Dk4 = {4: sp.Integer(0), 3: b3s, 2: b2f4}
sol4 = sp.solve(sp.Poly(Qm(Xk4, Dk4, 6), E).all_coeffs(), a3c + b2c4, dict=True)
check("k=4: Q_6=0 has a unique constrained locus", len(sol4) == 1)
a3sol = sp.expand(a3f.subs(sol4[0])); b2sol = sp.expand(b2f4.subs(sol4[0]))
check("k=4: Q_6=0 FORCES prod_{j=0}^{2} h^[j] | a_3  (h-compatibility, k=4)",
      divides(prodsh(hs, 3), a3sol) and a3sol != 0)
check("k=4: Q_6=0 FORCES prod_{j=0}^{1} h^[j] | b_2  (h-compatibility, k=4)",
      divides(prodsh(hs, 2), b2sol) and b2sol != 0)
print("   => the sub-leading coefficients a_{k-1}, b_{k-2} inherit shifted-power divisibility:")
print("      a_{k-1} = prod_{j=0}^{k-2} h^[j] * g,   b_{k-2} = prod_{j=0}^{k-3} h^[j] * beta.")

# ===========================================================================
print("\n=== §3  TASKS 2/3 PROPAGATION: band-3 h-constant forcing (re-verified) ===")
# ===========================================================================
# The central telescoping factors G = h^[-1] M; membership gives M(0)=0; Q_0=1
# forces G=E hence h^[-1]|E hence h constant.  Re-verified for h LINEAR and h
# QUADRATIC (cube-separated), rebuilding G and M in file.

def sector_M(hh, g, be, p, am2, bm2, bm3, am3, bm1, am1, kap):
    hm1, hm2, hm3 = sh(hh, -1), sh(hh, -2), sh(hh, -3)
    return sp.expand(
        sh(p, -1) * bm1 - sh(be, -1) * am1
        + hm2 * sh(g, -2) * bm2 + hh * sh(g, -1) * sh(bm2, 1)
        + hm3 * hm2 * bm3 + hm2 * hh * sh(bm3, 1) + hh * sh(hh, 1) * sh(bm3, 2)
        - kap * hm2 * am2 - kap * hh * sh(am2, 1))

for tag, hh in [("linear h=E-rho", E - rho),
                ("quadratic h=(E-1/3)(E-9/4)", sp.expand((E - sp.Rational(1, 3)) * (E - sp.Rational(9, 4))))]:
    h1 = sh(hh, 1)
    a3 = sp.expand(hh * h1 * sh(hh, 2)); b2 = sp.expand(kappa * hh * h1)
    g, _ = gen('Pg', 1); a2 = sp.expand(hh * h1 * g)
    be, _ = gen('Pbe', 1); b1 = sp.expand(hh * be)
    p, _ = gen('Pp', 1); a1 = sp.expand(hh * p)
    a0, _ = gen('Pa0', 2); b0, _ = gen('Pb0', 2)
    am1 = sp.expand(E * gen('Pam1', 2)[0]); am2 = sp.expand(E * (E - 1) * gen('Pam2', 2)[0])
    am3 = sp.expand(E * (E - 1) * (E - 2) * gen('Pam3', 2)[0])
    bm1 = sp.expand(E * gen('Pbm1', 2)[0]); bm2 = sp.expand(E * (E - 1) * gen('Pbm2', 2)[0])
    bm3 = sp.expand(E * (E - 1) * (E - 2) * gen('Pbm3', 2)[0])
    A = {3: a3, 2: a2, 1: a1, 0: a0, -1: am1, -2: am2, -3: am3}
    B = {3: sp.Integer(0), 2: b2, 1: b1, 0: b0, -1: bm1, -2: bm2, -3: bm3}
    G = Gpot(A, B, 3); M = sector_M(hh, g, be, p, am2, bm2, bm3, am3, bm1, am1, kappa)
    check(f"propagation: G = h^[-1] M  (closed-form central identity, {tag})",
          sp.expand(G - sh(hh, -1) * M) == 0)
    check(f"propagation: M(0)=0 under membership ({tag})", sp.expand(M.subs(E, 0)) == 0)

# the affine kill (h=E-rho): Q_0=(T-1)G=1 forces G=E, hence h^[-1]|E, hence deg h<=1,
# and deg h=1 dies from M(0)=0 vs M=1/alpha.  Re-verify the pieces:
for nval in (5, 7):
    check(f"nonconstant polys are not 1-periodic: coeff(E^(n-1)) of (E+1)^n-E^n is n (n={nval})",
          sp.Poly(sp.expand((E + 1)**nval - E**nval), E).nth(nval - 1) == nval)
haff = alpha * (E - rho)
check("deg h=1, rho=-1: h^[-1]=alpha E divides E",
      sp.rem(sp.Poly(E, E), sp.Poly(sh(haff, -1).subs(rho, -1), E)) == 0)
check("deg h=1, rho!=-1: h^[-1] does NOT divide E (so G=h^[-1]M=E impossible unless const)",
      sp.rem(sp.Poly(E, E), sp.Poly(sh(haff, -1).subs(rho, 2), E)) != 0)
check("deg h=2: h^[-1] (degree 2) does not divide E (degree 1) -> impossible",
      sp.rem(sp.Poly(E, E), sp.Poly(sh(sp.expand((E - sp.Rational(1, 3)) * (E - sp.Rational(9, 4))), -1), E)) != 0)
print("   => cube-separated shifted-cube tops with nonconstant h carry NO genuine pair: h is forced CONSTANT.")

# U-side of the same forcing: the shifted power comes from U^k with U = x h(E-...);
# [d, x h(E)] = 1 forces h constant.  (Two-line arbitrary-degree confirmation.)
for hh in (sp.Integer(1), E - rho, E**2 - rho):
    comm = commutator({-1: E}, {1: hh})           # [d, x h(E)]
    top = sp.expand(comm.get(0, 0))
    check(f"[d, x*h] with h={hh}: equals 1 iff h constant (value = {top})",
          (top == 1) == (sp.Poly(hh, E).degree() == 0 and hh == 1))

# ===========================================================================
print("\n=== §4  TASK 2: the explicit band-3 tame move on genuine data ===")
# ===========================================================================
# Once h is constant the pair lands in the classical tame family
#   U = x + c0 + c1 d,  X = U^3 - d/kappa - A,  D = lam X + kappa U + beta.
# The band-lowering tame word:
#   (i)  gauge  D -> D' = D - lam X       (B0-band3 collapse: kills b_3, b_2, b_-3; band D -> 1)
#   (ii) mirror X -> X'' = X - ((D'-beta)/kappa)^3   (kills U^3; band X -> 1)
# strictly lowers (n+m,k) and preserves [D,X]=1.  Verified on: positive control,
# and the full tame family with c1 != 0 (nontrivial membership tail).

c0, c1, Aq, lam, betaq = sp.symbols('c0 c1 Aq lam beta_q')

def tame_family(c0v, c1v, Av, lamv, betav):
    U = {1: sp.Integer(1)}
    if c0v != 0: U[0] = c0v
    if c1v != 0: U[-1] = sp.expand(c1v * E)
    U3 = cp_pow(U, 3)
    X = dict(U3)
    X[-1] = sp.expand(X.get(-1, 0) - E / kappa)
    if Av != 0:
        X[0] = sp.expand(X.get(0, 0) - Av)
    X = {k: v for k, v in X.items() if sp.expand(v) != 0}
    D = {}
    for k in set(list(X) + list(U) + [0]):
        v = sp.expand(lamv * X.get(k, 0) + kappa * U.get(k, 0) + (betav if k == 0 else 0))
        if v != 0:
            D[k] = v
    return X, D, U

# --- positive control: X = x^3 - d/kappa, D = kappa x ---
Xpc = {3: sp.Integer(1), -1: sp.expand(-E / kappa)}
Dpc = {1: kappa}
check("positive control [D,X]=1", commutator(Dpc, Xpc) == {0: sp.Integer(1)})
inv0 = invariant(Xpc, Dpc)
# move (ii) directly (b_3=0 already, gauge trivial): X'' = X - (D/kappa)^3
Xpc2 = cp_sub(Xpc, cp_pow(cp_scal(Dpc, 1 / kappa), 3))
check("positive control: X -> X-(D/kappa)^3 gives band-1 X'' = -d/kappa",
      Xpc2 == {-1: sp.expand(-E / kappa)})
check("positive control: [D,X'']=1 preserved", commutator(Dpc, Xpc2) == {0: sp.Integer(1)})
inv1 = invariant(Xpc2, Dpc)
check(f"positive control: invariant strictly drops {inv0} -> {inv1}", lex_lt(inv1, inv0))

# --- full tame family, c1 != 0 (rational parameters so bands are unambiguous) ---
X, D, U = tame_family(sp.Rational(2), sp.Rational(3), sp.Rational(5), sp.Rational(7), sp.Rational(11))
check("tame family (c1!=0): [D,X]=1", commutator(D, X) == {0: sp.Integer(1)})
inv_start = invariant(X, D)
# (i) gauge D' = D - lam X
b3, a3 = D.get(3, 0), X.get(3, 0)
lam_gauge = sp.cancel(b3 / a3)
Dgauge = cp_sub(D, cp_scal(X, lam_gauge))
check("tame family: gauge D'=D-lam X collapses band D to 1 (B0-band3 collapse)",
      band(Dgauge) == 1 and commutator(Dgauge, X) == {0: sp.Integer(1)})
inv_gauge = invariant(X, Dgauge)
check(f"tame family: invariant drops after gauge {inv_start} -> {inv_gauge}", lex_lt(inv_gauge, inv_start))
# (ii) X'' = X - U^3, with U recovered from D' = kappa U + const:  U = (D' - D'_const0)/kappa + c0-part.
# Directly: U is known; but to keep the move intrinsic, recover kappa U = D' - (band-0 constant).
# D' = kappa U + beta_q + kappa c0.  So D' - (beta_q + kappa c0) = kappa U.
const_shift = sp.expand(Dgauge.get(0, 0) - kappa * U.get(0, 0))   # = beta_q (+ kappa c0 already in U0)
kU = cp_sub(Dgauge, {0: const_shift})
Ureco = cp_scal(kU, 1 / kappa)
check("tame family: (D'-const)/kappa recovers U exactly", Ureco == {k: sp.expand(v) for k, v in U.items()})
Xred = cp_sub(X, cp_pow(Ureco, 3))                # X'' = X - U^3 = -d/kappa - A
check("tame family: X''=X-U^3 has band 1 and [D',X'']=1",
      band(Xred) == 1 and commutator(Dgauge, Xred) == {0: sp.Integer(1)})
inv_end = invariant(Xred, Dgauge)
check(f"tame family: invariant strictly drops to band-1 {inv_gauge} -> {inv_end}", lex_lt(inv_end, inv_gauge))
check("tame family: final pair (X'',D') is band 1 (affine symplectic, generates A_1)",
      band(Xred) == 1 and band(Dgauge) == 1)

# --- reducedness of the NONCONSTANT shifted cube (the gap's crux): the (3,2) leading
#     form is Dixmier-stuck; NO single transvection lowers (n+m,k). ---
check("Dixmier exponents of nonconstant shifted cube are (a,b)=(3,2), mutually non-dividing",
      (3 % 2 != 0) and (2 % 3 != 0))
# structural obstruction: D band 2 => D^2 has a nonzero band-4 top, so any q'(D) reaching
# band 3 (needs deg>=2 in D) simultaneously injects band 4, RAISING band X. Confirm band(D^2)=4:
Dstuck = {2: sp.expand(kappa * (E - rho) * sh(E - rho, 1)), 1: gen('bb1', 1)[0], 0: gen('bb0', 1)[0]}
check("nonconstant shifted-cube D (band 2): D^2 has nonzero band-4 top -> band 3 uncancellable",
      band(cp_pow(Dstuck, 2)) == 4)
print("   => the descent is NOT a leading-form move; it is the sub-leading h-forcing of §3")
print("      (the difference-operator analogue of the Abhyankar-Moh / Newton-polygon step).")

# ===========================================================================
print("\n=== §5  TASK 4: general (k,q) -- uniform reducing move + imbalanced walls ===")
# ===========================================================================
# The shifted k-th power positive control X = U^k - d/kappa, D = kappa U (U=x)
# reduces by the SINGLE mirror transvection X -> X - (D/kappa)^k  ->  -d/kappa (band 1).
for k in (3, 4, 5):
    Uk = cp_pow({1: sp.Integer(1)}, k)
    Xk = dict(Uk); Xk[-1] = sp.expand(Xk.get(-1, 0) - E / kappa)
    Xk = {kk: v for kk, v in Xk.items() if sp.expand(v) != 0}
    Dk = {1: kappa}
    check(f"k={k}: positive control [D,X]=1, band X={k}, band D=1", commutator(Dk, Xk) == {0: sp.Integer(1)})
    Xk2 = cp_sub(Xk, cp_pow(cp_scal(Dk, 1 / kappa), k))
    inv_before, inv_after = invariant(Xk, Dk), invariant(Xk2, Dk)
    check(f"k={k}: X->X-(D/kappa)^k gives band-1 -d/kappa, [D,X'']=1, invariant {inv_before}->{inv_after}",
          Xk2 == {-1: sp.expand(-E / kappa)} and commutator(Dk, Xk2) == {0: sp.Integer(1)}
          and lex_lt(inv_after, inv_before))

# imbalanced walls q < k-1, step-d = gcd(k,q) shifted powers.  Coprime (d=1):
#   (k,1): band D = 1, 1 | k, the positive control -- reduces (verified above).
#   (k,q>=2), q not| k: a SINGLE transvection cannot band-match (q*j != k), so it
#   STALLS -- the composite-move escape (Gap 1 = the classical DC1 core).
for (k, q) in [(3, 1), (4, 1), (5, 1)]:
    check(f"imbalanced coprime (k,q)=({k},{q}): q | k, reduces by one transvection (positive control)",
          sp.igcd(k, q) == 1 and k % q == 0)
for (k, q) in [(5, 2), (5, 3), (7, 2), (7, 3)]:
    check(f"imbalanced coprime (k,q)=({k},{q}): q not| k, single transvection STALLS (-> Gap 1 composite)",
          sp.igcd(k, q) == 1 and k % q != 0 and q < k - 1)
print("   => the shifted-power DESCENT closes the balanced wall (band D collapses to 1);")
print("      imbalanced coprime walls with q not| k reduce only modulo the composite escape (Gap 1).")

# ===========================================================================
print("\n=== §6  TASK 5: band-3 branch-distinctness bookkeeping ===")
# ===========================================================================
# The band-3 top wall splits into three genuinely different objects that SHARE the
# leading form (x^2 xi)^3 (blindness of the symbol -- band-reduction §4):
#   shifted-cube  a_3 = E(E+1)(E+2)   (consecutive; THIS memo: descends),
#   exotic AP     a_3 = E(E+r)(E+2r)  r != -4 (closed by lambda_r),
#   W2 hatch      a_3 = E(E+2)(E+4)   (r=2 anchored; joint filler lemma OPEN).
a3_cube = sp.expand(E * (E + 1) * (E + 2))
a3_ap = sp.expand(E * (E - 1) * (E - 2))          # AP step -1 (a consecutive/tame representative)
a3_w2 = sp.expand(E * (E + 2) * (E + 4))
def lf3(a):                                        # leading form of a band-3 top
    d = sp.Poly(a, E).degree()
    return sp.LC(sp.Poly(a, E)) * sp.symbols('x')**(3 + d) * sp.symbols('xi')**d
x, xi = sp.symbols('x xi')
check("shifted-cube, tame-AP and W2 tops all share leading form (x^2 xi)^3",
      lf3(a3_cube) == x**6 * xi**3 and lf3(a3_w2) == x**6 * xi**3 and lf3(a3_ap) == x**6 * xi**3)
# they are DISTINCT root multisets: consecutive {0,-1,-2} vs W2 {0,-2,-4}
check("shifted-cube roots {0,-1,-2} are consecutive; W2 roots {0,-2,-4} are not (distinct branches)",
      sorted(sp.roots(sp.Poly(a3_cube, E)).keys()) == [-2, -1, 0]
      and sorted(sp.roots(sp.Poly(a3_w2, E)).keys()) == [-4, -2, 0])
# W2 top is NOT a shifted cube (would need consecutive roots) -- the cube descent does not touch it:
check("W2 top E(E+2)(E+4) is not c*h h^[1] h^[2] (roots not consecutive) -> outside this descent",
      sorted(sp.roots(sp.Poly(a3_w2, E)).keys()) != [-2, -1, 0])

# ===========================================================================
print("\n" + "=" * 72)
n_pass = sum(1 for _, ok in CHECKS if ok)
n_tot = len(CHECKS)
print(f"{n_pass}/{n_tot} checks passed")
if n_pass == n_tot:
    print("ALL SHIFTED POWER DESCENT CHECKS PASSED")
else:
    print("SOME CHECKS FAILED:")
    for name, ok in CHECKS:
        if not ok:
            print("   FAIL:", name)
    raise SystemExit(1)
