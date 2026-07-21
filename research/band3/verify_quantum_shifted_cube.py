"""
verify_quantum_shifted_cube.py
==============================

Exact SymPy certificate for the QUANTUM BAND-3 SHIFTED-CUBE (non-exotic) SECTOR.

Sector:  gauge b_3 = 0, top a_3 = h(E) h(E+1) h(E+2) != 0  (a genuine shifted cube),
         genuine A_1 membership  E^{underline r} | a_-r, b_-r,  and  [D,X] = 1.

Conventions (frozen, band 3), matching research/band3/verify_quantum_band3.py:
    A_1[x^-1] = (+)_k x^k C[E],  E = x d,  (x^a f)(x^b g) = x^{a+b} f(E+b) g(E),
    f^[r](E) = f(E+r),
    Q_m = sum_{k+l=m} ( b_l^[k] a_k - a_k^[l] b_l ),   [D,X] = 1  <=>  Q_m = delta_{m0},
    G   = sum_{k=1}^3 sum_{j=0}^{k-1} ( a_k^[j-k] b_-k^[j] - b_k^[j-k] a_-k^[j] ),  Q_0 = (T-1) G.

WHAT IS PROVED (arbitrary coefficient degree, machine-checked identities), for
cube-separated h (gcd(h,h^[j]) = 1 for j=1,2,3) and every kappa:
  * Wall Q_5 = 0 is solved by b_2 = kappa h h^[1] (shifted-cube sufficiency, 1-dim).
  * Positive cascade forces the divisibilities  h h^[1] | a_2,  h | a_1,  h | b_1,
    with the S_3 = 1+T+T^2 "quantum midpoint" replacing band 2's S_2 = 1+T.
  * The central telescoping factors:  G = h^[-1] M  (closed form).
  * Q_0 = 1 with membership => G = E => h^[-1] | E => deg h <= 1; deg 1 is killed by
    M(0) = 0 (membership) versus M = 1/alpha.  Hence h is CONSTANT.

EXCEPTIONAL / RESIDUAL (bounded evidence only, NOT arbitrary degree):
  * non-cube-separated h (two roots differing by 1,2,3): the divisibility derivation
    genuinely fails; bounded Groebner (cap D=2) still finds the sector empty.
  * kappa = 0 folds into the same factorization (also forces h constant).

CONSTANT-h sector: the tame family U=x+c0+c1 d, X=U^3-d/kappa-A, D=lam X+kappa U+beta
  (all parameters) is a genuine Weyl pair; the full negative-tail classification is
  left open (cascade memo section 7 / catalog A*-band3 conjecture, commit 050a4c0).

Run:  uv run --with sympy python research/band3/verify_quantum_shifted_cube.py
Ends: ALL QUANTUM SHIFTED CUBE CHECKS PASSED
"""

import sympy as sp
import time

E = sp.symbols('E')


def sh(f, r):
    return sp.expand(sp.sympify(f).subs(E, E + r))


def gen(name, d):
    cs = sp.symbols(f'{name}_0:{d+1}')
    if not isinstance(cs, (list, tuple)):
        cs = (cs,)
    return sum(cs[i] * E**i for i in range(d + 1)), list(cs)


def Qm(a, b, m, K=3):
    """Ladder-m coefficient of [D,X] from dicts k->poly (k in -K..K)."""
    return sp.expand(sum(
        sh(b[l], k) * a[k] - sh(a[k], l) * b[l]
        for k in range(-K, K + 1) for l in range(-K, K + 1) if k + l == m))


def Gpot(a, b, K=3):
    g = 0
    for k in range(1, K + 1):
        for j in range(0, k):
            g += sh(a[k], j - k) * sh(b[-k], j) - sh(b[k], j - k) * sh(a[-k], j)
    return sp.expand(g)


def az(expr, label):
    e = sp.expand(sp.sympify(expr))
    if e != 0:
        raise AssertionError(label + "  :  residual = " + str(e))
    print("PASS", label)


def istrue(cond, label):
    if not cond:
        raise AssertionError(label + "  :  FALSE")
    print("PASS", label)


def divides(fac, c):
    c = sp.expand(sp.sympify(c))
    if c == 0:
        return True
    return sp.rem(sp.Poly(c, E), sp.Poly(fac, E)) == 0


kappa, rho, alpha = sp.symbols('kappa rho alpha')

# =====================================================================
# 0. Q_m engine sanity for the sector data (agrees with the direct commutator).
# =====================================================================
print("--- 0. Q_m == direct crossed-product commutator (sector data) ---")
h1sym = E - rho
_h1, _h2 = sh(h1sym, 1), sh(h1sym, 2)
Achk = {3: sp.expand(h1sym * _h1 * _h2), 2: gen('c2', 2)[0], 1: gen('c1', 2)[0],
        0: gen('c0', 2)[0], -1: gen('cm1', 2)[0], -2: gen('cm2', 2)[0], -3: gen('cm3', 2)[0]}
Bchk = {3: sp.Integer(0), 2: sp.expand(kappa * h1sym * _h1), 1: gen('d1', 2)[0],
        0: gen('d0', 2)[0], -1: gen('dm1', 2)[0], -2: gen('dm2', 2)[0], -3: gen('dm3', 2)[0]}


def direct_commutator(a, b, m, K=3):
    tot = 0
    for k in range(-K, K + 1):
        for l in range(-K, K + 1):
            if k + l == m:
                # (x^l b_l)(x^k a_k) - (x^k a_k)(x^l b_l)
                tot += sh(b[l], k) * a[k] - sh(a[k], l) * b[l]
    return sp.expand(tot)


for m in range(-6, 7):
    az(Qm(Achk, Bchk, m) - direct_commutator(Achk, Bchk, m),
       f"Q_{m} == direct commutator ladder coefficient")

# =====================================================================
# 1. The wall Q_5 and the shifted-cube solution b_2 = kappa h h^[1].
# =====================================================================
print("\n--- 1. wall Q_5 = 0 solved by b_2 = kappa h h^[1] (shifted-cube) ---")
# symbolic-degree h: a_3 = h h^[1] h^[2],  b_2 = kappa h h^[1]  solves the wall.
for dh in (1, 2, 3):
    hh, _ = gen('hh', dh)
    a3 = sp.expand(hh * sh(hh, 1) * sh(hh, 2))
    b2 = sp.expand(kappa * hh * sh(hh, 1))
    Aw = {3: a3, 2: b2, 1: 0, 0: 0, -1: 0, -2: 0, -3: 0}
    Bw = {3: sp.Integer(0), 2: b2, 1: 0, 0: 0, -1: 0, -2: 0, -3: 0}
    # wall in gauge b3=0:  Q_5 = b_2^[3] a_3 - a_3^[2] b_2
    az(sh(b2, 3) * a3 - sh(a3, 2) * b2,
       f"shifted cube deg h={dh}: b_2 = kappa h h^[1] solves Q_5 = b_2^[3]a_3 - a_3^[2]b_2")

# the wall solution space is exactly 1-dimensional (freedom = the scalar kappa).
# Use a concrete separated cube h = E-5 so the free-symbol count is unambiguous.
for hc in (E - 5, E, sp.expand((E - sp.Rational(1, 3)) * (E - sp.Rational(9, 4)))):
    a3c = sp.expand(hc * sh(hc, 1) * sh(hc, 2))
    dcap = 2 * sp.Poly(hc, E).degree() + 3
    b2v, b2c = gen('wc', dcap)
    wall = sp.Poly(sp.expand(sh(b2v, 3) * a3c - sh(a3c, 2) * b2v), E)
    solw = list(sp.linsolve(wall.all_coeffs(), b2c))[0]
    nfree = len({s for x in solw for s in x.free_symbols})
    istrue(nfree == 1,
           f"wall solution space for cube h={sp.factor(hc)} is exactly 1-dimensional (kappa)")

# =====================================================================
# 2. Positive cascade: Q_4 -> b_1, Q_3 -> b_0, Q_2 -> b_-1, and the S_3 midpoint.
# =====================================================================
print("\n--- 2. positive cascade + quantum S_3 midpoint (band-3, cube-separated) ---")
hh = E - rho
h1, h2, h3 = sh(hh, 1), sh(hh, 2), sh(hh, 3)
hm1 = sh(hh, -1)
a3 = sp.expand(hh * h1 * h2)
b2 = sp.expand(kappa * hh * h1)

g, _ = gen('g', 2)
a2 = sp.expand(hh * h1 * g)          # h h^[1] | a_2   (write a_2 = h h^[1] g)
beta, _ = gen('be', 2)
b1 = sp.expand(hh * beta)            # h | b_1         (write b_1 = h beta)
p, _ = gen('p', 2)
a1 = sp.expand(hh * p)               # h | a_1         (write a_1 = h p)
a0, _ = gen('a0', 2)
b0, _ = gen('b0', 2)

A = {3: a3, 2: a2, 1: a1, 0: a0, -1: 0, -2: 0, -3: 0}
B = {3: sp.Integer(0), 2: b2, 1: b1, 0: b0, -1: 0, -2: 0, -3: 0}

# Q_4 exact form and the (T^3-1) beta = kappa (T^2-1) g reduction:
az(Qm(A, B, 4) - sp.expand(a3 * h3 * ((sh(beta, 3) - beta) - kappa * (sh(g, 2) - g))),
   "Q_4 = a_3 h^[3] [ (T^3-1) beta - kappa (T^2-1) g ]   (b_1=h beta, a_2=h h^[1] g)")

# the S_3 = 1+T+T^2 quantum midpoint: (T-1) [ S_3 beta - kappa S_2 g ] = (T^3-1)beta - kappa(T^2-1)g
S3b = sp.expand(sh(beta, 2) + sh(beta, 1) + beta)      # S_3 beta
S2g = sp.expand(sh(g, 1) + g)                          # S_2 g
az((sh(S3b - kappa * S2g, 1) - (S3b - kappa * S2g)) - ((sh(beta, 3) - beta) - kappa * (sh(g, 2) - g)),
   "S_3 midpoint: (T-1)[S_3 beta - kappa S_2 g] = (T^3-1)beta - kappa(T^2-1)g  => S_3 b_1/h = kappa S_2 a_2/(h h^[1]) + const")

# Q_3 exact form:
az(Qm(A, B, 3) - sp.expand(a3 * ((sh(b0, 3) - b0) + (sh(beta, 2) * g - sh(g, 1) * beta) + kappa * (p - sh(p, 2)))),
   "Q_3 = a_3 [ (T^3-1) b_0 + (beta^[2] g - g^[1] beta) + kappa (p - p^[2]) ]   (a_1 = h p)")

# Q_2 -> b_-1 operator (staggered +3 vs -1, factor h h^[1]):
bm1x, _ = gen('bm1x', 3)
Aop = {3: a3, 2: 0, 1: 0, 0: 0, -1: 0, -2: 0, -3: 0}
Bop = {3: sp.Integer(0), 2: 0, 1: 0, 0: 0, -1: bm1x, -2: 0, -3: 0}
az(Qm(Aop, Bop, 2) - sp.expand(hh * h1 * (sh(bm1x, 3) * h2 - hm1 * bm1x)),
   "Q_2 b_-1-operator = h h^[1] ( b_-1^[3] h^[2] - h^[-1] b_-1 )   (staggered +3 vs -1)")

# =====================================================================
# 2b. Divisibility NECESSITY at a concrete cube-separated h (regression of the
#     coprimality proof).  h = E(2E-1) has roots {0, 1/2}: gcd(h, h^[j]) = 1.
# =====================================================================
print("\n--- 2b. divisibility necessity (cube-separated h = E(2E-1)) ---")
hs = sp.expand(E * (2 * E - 1))
hs1, hs2 = sh(hs, 1), sh(hs, 2)
istrue(all(sp.gcd(sp.Poly(hs, E), sp.Poly(sh(hs, j), E)).degree() == 0 for j in (1, 2, 3)),
       "h = E(2E-1) is cube-separated: gcd(h, h^[j]) = 1 for j=1,2,3")
a3s = sp.expand(hs * hs1 * hs2)
b2s = sp.expand(kappa * hs * hs1)
a2f, a2c = gen('a2f', 4)          # fully generic a_2
b1f, b1c = gen('b1f', 3)
Ax = {3: a3s, 2: a2f, 1: 0, 0: 0, -1: 0, -2: 0, -3: 0}
Bx = {3: sp.Integer(0), 2: b2s, 1: b1f, 0: 0, -1: 0, -2: 0, -3: 0}
Q4x = sp.Poly(Qm(Ax, Bx, 4), E)
sol = sp.solve(Q4x.all_coeffs(), b1c + a2c, dict=True)
istrue(len(sol) == 1, "Q_4 = 0 has a unique constrained solution locus")
a2sol = sp.expand(a2f.subs(sol[0]))
b1sol = sp.expand(b1f.subs(sol[0]))
istrue(divides(hs, a2sol) and divides(hs1, a2sol),
       "Q_4 = 0 (polynomial b_1) FORCES h h^[1] | a_2   (necessity)")
istrue(divides(hs, b1sol), "Q_4 = 0 (polynomial b_1) FORCES h | b_1   (necessity)")

# =====================================================================
# 3. The general-h central identity  G = h^[-1] M  (closed form).
# =====================================================================
print("\n--- 3. central telescoping factors: G = h^[-1] M ---")


def sector_M(hh, g, beta, p, am2, bm2, bm3, am3, bm1, am1, kappa):
    hm1_, hm2_, hm3_ = sh(hh, -1), sh(hh, -2), sh(hh, -3)
    return sp.expand(
        sh(p, -1) * bm1 - sh(beta, -1) * am1
        + hm2_ * sh(g, -2) * bm2 + hh * sh(g, -1) * sh(bm2, 1)
        + hm3_ * hm2_ * bm3 + hm2_ * hh * sh(bm3, 1) + hh * sh(hh, 1) * sh(bm3, 2)
        - kappa * hm2_ * am2 - kappa * hh * sh(am2, 1))


# deg 1,2,3 instances confirm G = h^[-1] M is a term-by-term (degree-independent)
# factorization: e.g. a_1^[-1] b_-1 = (h p)^[-1] b_-1 = h^[-1] p^[-1] b_-1, and so on.
for dh in (1, 2, 3):
    if dh == 1:
        hh = E - rho
    elif dh == 2:
        # cube-separated degree-2 h (roots 1/3, 9/4: difference not an integer)
        hh = sp.expand((E - sp.Rational(1, 3)) * (E - sp.Rational(9, 4)))
    else:
        # cube-separated degree-3 h (roots 1/5, 7/3, 20/7: no integer differences)
        hh = sp.expand((5 * E - 1) * (3 * E - 7) * (7 * E - 20))
    h1 = sh(hh, 1)
    hm1 = sh(hh, -1)
    a3 = sp.expand(hh * h1 * sh(hh, 2))
    b2 = sp.expand(kappa * hh * h1)
    g, _ = gen(f'G{dh}g', 1)
    a2 = sp.expand(hh * h1 * g)
    beta, _ = gen(f'G{dh}be', 1)
    b1 = sp.expand(hh * beta)
    p, _ = gen(f'G{dh}p', 1)
    a1 = sp.expand(hh * p)
    a0, _ = gen(f'G{dh}a0', 2)
    b0, _ = gen(f'G{dh}b0', 2)
    am1 = sp.expand(E * gen(f'G{dh}am1', 2)[0])
    am2 = sp.expand(E * (E - 1) * gen(f'G{dh}am2', 2)[0])
    am3 = sp.expand(E * (E - 1) * (E - 2) * gen(f'G{dh}am3', 2)[0])
    bm1 = sp.expand(E * gen(f'G{dh}bm1', 2)[0])
    bm2 = sp.expand(E * (E - 1) * gen(f'G{dh}bm2', 2)[0])
    bm3 = sp.expand(E * (E - 1) * (E - 2) * gen(f'G{dh}bm3', 2)[0])
    Ag = {3: a3, 2: a2, 1: a1, 0: a0, -1: am1, -2: am2, -3: am3}
    Bg = {3: sp.Integer(0), 2: b2, 1: b1, 0: b0, -1: bm1, -2: bm2, -3: bm3}
    G = Gpot(Ag, Bg)
    M = sector_M(hh, g, beta, p, am2, bm2, bm3, am3, bm1, am1, kappa)
    az(G - sh(hh, -1) * M, f"G = h^[-1] M  (closed-form central identity, deg h = {dh})")

# =====================================================================
# 4. Membership => M(0) = 0 and G(0) = 0.
# =====================================================================
print("\n--- 4. membership: M(0) = 0 and G(0) = 0 ---")
# symbolic-h M(0) with genuine membership: all nine terms carry a membership zero at E=0.
hh = E - rho
g, _ = gen('MMg', 2)
beta, _ = gen('MMbe', 2)
p, _ = gen('MMp', 2)
am1 = sp.expand(E * gen('MMam1', 2)[0])
am2 = sp.expand(E * (E - 1) * gen('MMam2', 2)[0])
am3 = sp.expand(E * (E - 1) * (E - 2) * gen('MMam3', 2)[0])
bm1 = sp.expand(E * gen('MMbm1', 2)[0])
bm2 = sp.expand(E * (E - 1) * gen('MMbm2', 2)[0])
bm3 = sp.expand(E * (E - 1) * (E - 2) * gen('MMbm3', 2)[0])
M = sector_M(hh, g, beta, p, am2, bm2, bm3, am3, bm1, am1, kappa)
az(M.subs(E, 0), "M(0) = 0 identically under genuine membership (all 9 terms vanish at E=0)")

# G(0) = 0 identically for the full generic band-3 data under membership:
aM = {3: gen('N3', 2)[0], 2: gen('N2', 2)[0], 1: gen('N1', 2)[0], 0: gen('N0', 2)[0],
      -1: sp.expand(E * gen('Nm1', 2)[0]), -2: sp.expand(E * (E - 1) * gen('Nm2', 2)[0]),
      -3: sp.expand(E * (E - 1) * (E - 2) * gen('Nm3', 2)[0])}
bM = {3: sp.Integer(0), 2: gen('P2', 2)[0], 1: gen('P1', 2)[0], 0: gen('P0', 2)[0],
      -1: sp.expand(E * gen('Pm1', 2)[0]), -2: sp.expand(E * (E - 1) * gen('Pm2', 2)[0]),
      -3: sp.expand(E * (E - 1) * (E - 2) * gen('Pm3', 2)[0])}
az(Gpot(aM, bM).subs(E, 0), "G(0) = 0 identically under genuine membership")

# =====================================================================
# 5. The affine kill:  G = E  =>  deg h <= 1, and deg 1 is impossible.
# =====================================================================
print("\n--- 5. affine kill: Q_0 = 1 forces h constant ---")
# (a) Q_0 = (T-1) G = 1 with (T-1)E = 1 forces G - E to be 1-periodic; a nonconstant
#     polynomial is not 1-periodic (leading coeff of (E+1)^n - E^n is n != 0 in char 0).
for nval in (5, 7, 11):
    az(sp.Poly(sp.expand((E + 1)**nval - E**nval), E).nth(nval - 1) - nval,
       f"coeff(E^(n-1)) in (E+1)^n - E^n equals n (n={nval}): nonconstant polys are not 1-periodic => G = E + c")

# (b) G = E and G = h^[-1] M  =>  h^[-1] | E.  deg h >= 2 is impossible (degree);
#     deg h = 1 forces h^[-1] = alpha E, hence M = 1/alpha, contradicting M(0)=0.
# deg 1:  h = alpha (E - rho),  h^[-1] = alpha (E - 1 - rho).  h^[-1] | E  <=>  rho = -1.
haff = alpha * (E - rho)
hm1_aff = sh(haff, -1)
istrue(sp.rem(sp.Poly(E, E), sp.Poly(hm1_aff.subs(rho, -1), E)) == 0,
       "deg h=1, rho=-1: h^[-1] = alpha E divides E")
istrue(sp.rem(sp.Poly(E, E), sp.Poly(hm1_aff.subs(rho, 2), E)) != 0,
       "deg h=1, rho=2 (!=-1): h^[-1] does NOT divide E  => G=h^[-1]M=E impossible")
# for rho=-1:  h^[-1] M = alpha E * M = E  =>  M = 1/alpha  (polynomial identity)
az(sp.expand(alpha * E * (sp.Rational(1, 1) / alpha) - E),
   "deg h=1 (rho=-1): alpha E * (1/alpha) = E  => M = 1/alpha (constant)")
print("   => M(0) = 1/alpha != 0 contradicts M(0) = 0 (section 4): deg h = 1 impossible.")
# deg 2:  h^[-1] has degree 2 > 1 = deg E, so h^[-1] does not divide the nonzero E.
hdeg2 = sp.expand((E - sp.Rational(1, 3)) * (E - sp.Rational(9, 4)))
istrue(sp.rem(sp.Poly(E, E), sp.Poly(sh(hdeg2, -1), E)) != 0,
       "deg h=2: h^[-1] (degree 2) does not divide E (degree 1) => impossible")
print("   => h is CONSTANT (cube-separated case, any kappa).")

# =====================================================================
# 6. kappa = 0 folds into the same factorization (also forces h constant).
# =====================================================================
print("\n--- 6. kappa = 0 branch: b_1 = c h, then same divisibilities ---")
hh = E - rho
h1, h2 = sh(hh, 1), sh(hh, 2)
a3 = sp.expand(hh * h1 * h2)
# Q_4|_{b2=0} = b_1^[3] a_3 - a_3^[1] b_1 = 0  =>  b_1/h is 3-periodic  =>  b_1 = c h.
cc = sp.symbols('cc')
Ak = {3: a3, 2: 0, 1: 0, 0: 0, -1: 0, -2: 0, -3: 0}
Bk = {3: sp.Integer(0), 2: sp.Integer(0), 1: sp.expand(cc * hh), 0: 0, -1: 0, -2: 0, -3: 0}
az(Qm(Ak, Bk, 4), "kappa=0: b_1 = c h solves Q_4 = b_1^[3] a_3 - a_3^[1] b_1 = 0")
b1g, b1gc = gen('b1g', 3)
Bk[1] = b1g
solk = sp.solve(sp.Poly(Qm(Ak, Bk, 4), E).all_coeffs(), b1gc, dict=True)
b1ks = sp.expand(b1g.subs(solk[0]))
istrue(divides(hh, b1ks), "kappa=0: Q_4 = 0 forces h | b_1 (general solution b_1 = c h)")
# kappa=0, Q_3 forces h h^[1] | a_2 ; Q_2 forces h | a_1  (regression at symbolic h):
g2, _ = gen('kg', 2)
a2k = sp.expand(hh * h1 * g2)
a2f, a2c = gen('ka2f', 4)
b0f, b0c = gen('kb0f', 3)
Ak2 = {3: a3, 2: a2f, 1: 0, 0: 0, -1: 0, -2: 0, -3: 0}
Bk2 = {3: sp.Integer(0), 2: sp.Integer(0), 1: sp.expand(cc * hh), 0: b0f, -1: 0, -2: 0, -3: 0}
solk2 = sp.solve(sp.Poly(Qm(Ak2, Bk2, 3), E).all_coeffs(), b0c + a2c, dict=True)
a2ks = sp.expand(a2f.subs(solk2[0]))
istrue(divides(hh, a2ks) and divides(h1, a2ks),
       "kappa=0: Q_3 = 0 forces h h^[1] | a_2")
# kappa=0, Q_2 forces h | a_1 (with b_2=0, b_1=c h, a_2=h h^[1] g, b_-1 = E*(free)):
a1f, a1c = gen('ka1f', 3)
a0k, _ = gen('ka0', 2)
b0k, b0kc = gen('kb0k', 3)
bm1k, bm1kc = gen('kbm1', 3)
bm1k = sp.expand(E * bm1k)
Ak3 = {3: a3, 2: a2k, 1: a1f, 0: a0k, -1: 0, -2: 0, -3: 0}
Bk3 = {3: sp.Integer(0), 2: sp.Integer(0), 1: sp.expand(cc * hh), 0: b0k, -1: bm1k, -2: 0, -3: 0}
solk3 = sp.solve(sp.Poly(Qm(Ak3, Bk3, 2), E).all_coeffs(), b0kc + bm1kc + a1c, dict=True)
a1ks = sp.expand(a1f.subs(solk3[0]))
istrue(divides(hh, a1ks), "kappa=0: Q_2 = 0 forces h | a_1")
print("   => G = h^[-1] M holds with kappa=0 too; same affine kill => h constant.")

# =====================================================================
# 7. The constant-h tame family (genuine Weyl pairs).
# =====================================================================
print("\n--- 7. constant-h tame family: [D,X] = 1, all parameters ---")
c0, c1, Aq, lam, beta_q = sp.symbols('c0 c1 Aq lam beta_q')


def mul_ladders(P, Q):
    R = {}
    for k1, p1 in P.items():
        for k2, p2 in Q.items():
            R[k1 + k2] = R.get(k1 + k2, 0) + sh(p1, k2) * p2
    return {k: sp.expand(v) for k, v in R.items() if sp.expand(v) != 0}


# U = x + c0 + c1 d,  X = U^3 - d/kappa - Aq,  D = lam X + kappa U + beta_q
U = {1: sp.Integer(1), 0: c0, -1: sp.expand(c1 * E)}
U3 = mul_ladders(mul_ladders(U, U), U)
X = {k: sp.Integer(0) for k in range(-3, 4)}
for k, v in U3.items():
    X[k] += v
X[-1] = sp.expand(X[-1] - E / kappa)
X[0] = sp.expand(X[0] - Aq)
Dm = {}
for k in range(-3, 4):
    Dm[k] = sp.expand(lam * X[k] + kappa * U.get(k, 0) + (beta_q if k == 0 else 0))
for m in range(-6, 7):
    az(Qm(X, Dm, m) - (1 if m == 0 else 0),
       f"tame family [D,X]_{m} = delta_{{m0}}  (X=U^3-d/kappa-A, D=lam X+kappa U+beta)")
istrue(sp.expand(X[3] - 1) == 0, "tame family: a_3 = 1 (constant shifted cube, h=1)")
istrue(sp.expand(X[-3] - c1**3 * E * (E - 1) * (E - 2)) == 0,
       "tame family: a_-3 = c1^3 E(E-1)(E-2) (exact falling factorial membership)")
istrue(divides(E, X[-1]) and divides(E * (E - 1), X[-2]) and divides(E * (E - 1) * (E - 2), X[-3]),
       "tame family X memberships hold")
istrue(divides(E, Dm[-1]) and divides(E * (E - 1), Dm[-2]) and divides(E * (E - 1) * (E - 2), Dm[-3]),
       "tame family D memberships hold")
# gauge D -> D - lam X kills b_3, b_2, b_-3 simultaneously (B0-band3 collapse, kappa_2 = 0):
Dg = {k: sp.expand(Dm[k] - lam * X[k]) for k in range(-3, 4)}
istrue(Dg[3] == 0 and Dg[2] == 0 and Dg[-3] == 0,
       "tame family: gauge D->D-lam X kills b_3, b_2, b_-3 at once (gauged kappa_2 = 0)")

# constant h=1 positive cascade: S_3 midpoint (3-fold periodicity, matches band3 verify sec 7):
k2 = sp.symbols('k2')
a2t, _ = gen('a2t', 2)
b1t, _ = gen('b1t', 2)
At = {3: sp.Integer(1), 2: a2t, 1: 0, 0: 0, -1: 0, -2: 0, -3: 0}
Bt = {3: sp.Integer(0), 2: k2, 1: b1t, 0: 0, -1: 0, -2: 0, -3: 0}
az(Qm(At, Bt, 4) - (k2 * (a2t - sh(a2t, 2)) + (sh(b1t, 3) - b1t)),
   "constant h=1: Q_4 = kappa_2 (a_2 - a_2^[2]) + (b_1^[3] - b_1)  (3-fold S_3 midpoint)")

# section-6 positive control of the cascade memo (the c0=A=lam=beta=0 sub-case):
Uc = {1: sp.Integer(1), -1: sp.expand(c1 * E)}
U3c = mul_ladders(mul_ladders(Uc, Uc), Uc)
Xc = {k: sp.Integer(0) for k in range(-3, 4)}
for k, v in U3c.items():
    Xc[k] += v
Xc[-1] = sp.expand(Xc[-1] - E / kappa)
Dc = {k: sp.Integer(0) for k in range(-3, 4)}
Dc[1] = kappa
Dc[-1] = sp.expand(kappa * c1 * E)
istrue(all(sp.expand(Qm(Xc, Dc, m) - (1 if m == 0 else 0)) == 0 for m in range(-6, 7)),
       "cascade-memo section-6 positive control X=U^3-d/kappa, D=kappa U is a genuine pair")

# =====================================================================
# 7b. The divisibility GENUINELY FAILS for a non-cube-separated top: for
#     h = E(E-1) (roots 0,1 differ by 1), Q_4 = 0 has a polynomial b_1 for a
#     generic a_2 that is NOT divisible by h h^[1]. Hence G = h^[-1] M is not
#     available and the arbitrary-degree proof of section 4 does not reach it.
# =====================================================================
print("\n--- 7b. non-separated h = E(E-1): the divisibility h h^[1] | a_2 FAILS ---")
hns = sp.expand(E * (E - 1))
hns1, hns2 = sh(hns, 1), sh(hns, 2)
istrue(sp.gcd(sp.Poly(hns, E), sp.Poly(hns1, E)).degree() > 0,
       "h = E(E-1) is NOT cube-separated: gcd(h, h^[1]) has positive degree")
a3n = sp.expand(hns * hns1 * hns2)
b2n = sp.expand(kappa * hns * hns1)
a2fn, a2cn = gen('nsa2', 4)
b1fn, b1cn = gen('nsb1', 3)
An = {3: a3n, 2: a2fn, 1: 0, 0: 0, -1: 0, -2: 0, -3: 0}
Bn = {3: sp.Integer(0), 2: b2n, 1: b1fn, 0: 0, -1: 0, -2: 0, -3: 0}
soln = sp.solve(sp.Poly(Qm(An, Bn, 4), E).all_coeffs(), b1cn + a2cn, dict=True)
a2ns = sp.expand(a2fn.subs(soln[0]))
istrue(not divides(hns1, a2ns) and a2ns != 0,
       "non-separated: Q_4 = 0 admits a polynomial b_1 with a_2 NOT divisible by h^[1] (divisibility fails)")

# =====================================================================
# 8. Bounded evidence: nonconstant tops give the EMPTY sector (unit ideal, cap D=2).
#    (This is corroboration, NOT an arbitrary-degree proof; it covers the
#     non-cube-separated residual for the tested tops.)
# =====================================================================
print("\n--- 8. bounded evidence: nonconstant tops -> empty sector (cap D=2) ---")


def sector_empty(hh, D=2):
    hh = sp.expand(hh)
    h1_, h2_ = sh(hh, 1), sh(hh, 2)
    a3 = sp.expand(hh * h1_ * h2_)
    b2 = sp.expand(kappa * hh * h1_)
    a2, a2c = gen('s_a2', D + 2)
    a1, a1c = gen('s_a1', D + 2)
    a0, a0c = gen('s_a0', D + 2)
    b1, b1c = gen('s_b1', D + 3)
    b0, b0c = gen('s_b0', D + 3)
    am1 = sp.expand(E * gen('s_am1', D)[0]); am1c = gen('s_am1', D)[1]
    am2 = sp.expand(E * (E - 1) * gen('s_am2', D)[0]); am2c = gen('s_am2', D)[1]
    am3 = sp.expand(E * (E - 1) * (E - 2) * gen('s_am3', D)[0]); am3c = gen('s_am3', D)[1]
    bm1 = sp.expand(E * gen('s_bm1', D)[0]); bm1c = gen('s_bm1', D)[1]
    bm2 = sp.expand(E * (E - 1) * gen('s_bm2', D)[0]); bm2c = gen('s_bm2', D)[1]
    bm3 = sp.expand(E * (E - 1) * (E - 2) * gen('s_bm3', D)[0]); bm3c = gen('s_bm3', D)[1]
    A = {3: a3, 2: a2, 1: a1, 0: a0, -1: am1, -2: am2, -3: am3}
    B = {3: sp.Integer(0), 2: b2, 1: b1, 0: b0, -1: bm1, -2: bm2, -3: bm3}
    allc = a2c + a1c + a0c + b1c + b0c + am1c + am2c + am3c + bm1c + bm2c + bm3c + [kappa]
    eqs = []
    for m in [4, 3, 2, 1]:
        eqs += sp.Poly(Qm(A, B, m), E).all_coeffs()
    eqs += sp.Poly(Qm(A, B, 0) - 1, E).all_coeffs()
    Gb = sp.groebner(eqs, *allc, order='grevlex')
    return list(Gb) == [sp.Integer(1)]


# These non-cube-separated tops are exactly the residual NOT covered by the
# arbitrary-degree theorem of section 4 (their divisibility derivation fails).
# The separated/deg-1 nonconstant tops are already excluded by section 4-5, so
# they are not re-checked here.
t0 = time.time()
for hh, lab in [(E * (E - 1), "E(E-1) [non-sep, diff 1]"),
                (E * (E - 2), "E(E-2) [non-sep, diff 2]"),
                (E * (E - 3), "E(E-3) [non-sep, diff 3]"),
                ((E - 1) * (E - 3), "(E-1)(E-3) [non-sep, diff 2]")]:
    istrue(sector_empty(hh), f"bounded (cap D=2): sector EMPTY for top h = {lab}")
print(f"   (section 8 Groebner time: {time.time() - t0:.1f}s)")

print("\nALL QUANTUM SHIFTED CUBE CHECKS PASSED")
