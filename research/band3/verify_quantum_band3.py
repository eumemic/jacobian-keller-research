#!/usr/bin/env python3
"""
verify_quantum_band3.py
=======================
Exact SymPy verification backing the memo `quantum-band3-cascade.md`:
the QUANTUM band-3 cascade for A_1 = C<x,d> ([d,x]=1), E = x d, in the
E-ladder presentation with supports contained in [-3,3] (ladders m in [-6,6]).

Convention (frozen, matching every sibling quantum memo):
    A_1[x^{-1}] = (+)_{k in Z} x^k C[E],   (x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E),
    f^[r](E) = f(E+r),   E = x d,   d = x^{-1} E,
    ladder-m coefficient of [D,X]:  Q_m = sum_{k+l=m} ( b_l^[k] a_k - a_k^[l] b_l ),
    Keller/Weyl relation [D,X] = 1  <=>  Q_m = delta_{m0}.
    Genuine A_1 membership:  E(E-1)...(E-r+1) | a_{-r}, b_{-r}.  Band 3:
      E | a_-1,b_-1 ;  E(E-1) | a_-2,b_-2 ;  E(E-1)(E-2) | a_-3,b_-3.

Sections:
  0. crossed-product engine; Q_m == separately looped crossed-product expansion (m=-6..6); closed forms.
  1. Q_6 top: b3^[3]a3 - a3^[3]b3; period-3 rigidity; gauge legitimacy => b3=0.
  2. THE WALL Q_5 (gauge b3=0) = b2^[3]a3 - a3^[2]b2:
       - shifted-cube sufficiency at symbolic degree;
       - shift-operator reduction  S(1+S)A = (1+S+S^2)B  (root-multiset calculus);
       - COUNTEREXAMPLE a3=E(E-2)(E-4): solves wall, NOT a shifted cube (refutes the
         naive conjecture); 1-dim freedom; degree law 2 deg a3 = 3 deg b2;
       - nonnegativity is a genuine extra condition (Phi3-div but no b2).
  3. Q_0 telescoping: general closed-form local potential G with Q_0=(T-1)G;
       reduces to band-2 G; G(0)=0 by membership => Q_0=1 gives G=E.
  4. bottom mirror Q_-6, Q_-5; reflection E->-E-1 flips shifts; cross-coupling of
       proportionality constants (lambda3 vs mu3): the first genuinely new phenomenon.
  5. negative-tail staggered rigidity lemmas (Lemma-R family) with leading coeffs.
  6. positive control: a GENUINE band-3 pair X=U^3-d/kappa, D=kappa U (U=x+c1 d).
  7. bounded corroboration.

Run:  uv run --with sympy python research/band3/verify_quantum_band3.py
Ends: ALL QUANTUM BAND3 CHECKS PASSED
"""
import sympy as sp
import itertools

E = sp.symbols('E')
S = sp.symbols('S')
lam, mu3, kappa, c1 = sp.symbols('lambda mu3 kappa c1')


def sh(f, r):
    return sp.expand(sp.sympify(f).subs(E, E + r))


def poly(name, deg):
    cs = sp.symbols(f'{name}_0:{deg+1}')
    return sum(cs[i] * E**i for i in range(deg + 1)), list(cs)


def Qm(a, b, m, K=3):
    """Ladder-m coefficient of [D,X] from dicts k->poly (k in -K..K)."""
    return sp.expand(sum(
        sh(b[l], k) * a[k] - sh(a[k], l) * b[l]
        for k in range(-K, K + 1) for l in range(-K, K + 1) if k + l == m))


def az(expr, label):
    e = sp.expand(sp.sympify(expr))
    if e != 0:
        raise AssertionError(label + "  :  residual = " + str(e))
    print("PASS", label)


def istrue(cond, label):
    if not cond:
        raise AssertionError(label + "  :  FALSE")
    print("PASS", label)


def falling(r):
    return sp.prod([E - i for i in range(r)]) if r > 0 else sp.Integer(1)


def divides(fac, c):
    c = sp.expand(sp.sympify(c))
    if c == 0:
        return True
    return sp.rem(sp.Poly(c, E), sp.Poly(fac, E)) == 0


# =====================================================================
# 0. Crossed-product engine: Q_m equals the direct commutator ladder coeff.
# =====================================================================
print("--- 0. Q_m == separately looped crossed-product expansion (band 3) ---")
A = {k: poly(f'A{k+3}', 2)[0] for k in range(-3, 4)}
Bc = {k: poly(f'B{k+3}', 2)[0] for k in range(-3, 4)}


def expanded_commutator(m, K=3):
    # [D,X] = DX - XD;  (x^l b_l)(x^k a_k) = x^{l+k} b_l(E+k) a_k(E),
    #                   (x^k a_k)(x^l b_l) = x^{k+l} a_k(E+l) b_l(E).
    tot = 0
    for k in range(-K, K + 1):
        for l in range(-K, K + 1):
            if k + l == m:
                tot += sh(Bc[l], k) * A[k] - sh(A[k], l) * Bc[l]
    return sp.expand(tot)


for m in range(-6, 7):
    az(expanded_commutator(m) - Qm(A, Bc, m),
       f"Q_{m} agrees with the separately looped crossed-product expansion")

az(Qm(A, Bc, 6) - (sh(Bc[3], 3) * A[3] - sh(A[3], 3) * Bc[3]),
   "Q_6 = b3^[3] a3 - a3^[3] b3")
Bg0 = dict(Bc); Bg0[3] = sp.Integer(0)
az(Qm(A, Bg0, 5) - (sh(Bg0[2], 3) * A[3] - sh(A[3], 2) * Bg0[2]),
   "Q_5|_{b3=0} = b2^[3] a3 - a3^[2] b2")
az(Qm(A, Bc, -6) - (sh(Bc[-3], -3) * A[-3] - sh(A[-3], -3) * Bc[-3]),
   "Q_-6 = b-3^[-3] a-3 - a-3^[-3] b-3")
Bg0m = dict(Bc); Bg0m[-3] = sp.Integer(0)
az(Qm(A, Bg0m, -5) - (sh(Bg0m[-2], -3) * A[-3] - sh(A[-3], -2) * Bg0m[-2]),
   "Q_-5|_{b-3=0} = b-2^[-3] a-3 - a-3^[-2] b-2")


# =====================================================================
# 1. Q_6 top proportionality + gauge legitimacy.
# =====================================================================
print("\n--- 1. Q_6 top proportionality; gauge D -> D - lambda X ---")
# Q_6 = 0 <=> b3^[3] a3 = a3^[3] b3 <=> b3/a3 is 3-periodic (as a rational function)
# <=> constant (rational periodicity: a pole propagates along an infinite 3Z orbit).
# Certify the rational-periodicity engine exactly, as in J2q:  with r = b3/a3,
#   r^[3] - r = (b3^[3] a3 - a3^[3] b3) / (a3 a3^[3])   [ = -Q_6/(a3 a3^[3]) ].
a3g, _ = poly('a3g', 3)
b3g, _ = poly('b3g', 3)
rr = b3g / a3g
az(sp.simplify((sh(rr, 3) - rr) - (sh(b3g, 3) * a3g - sh(a3g, 3) * b3g) / (a3g * sh(a3g, 3))),
   "rational 3-periodicity engine: r^[3]-r = (b3^[3]a3 - a3^[3]b3)/(a3 a3^[3]), r=b3/a3")
# period-3 rigidity for the polynomial quotient (3-step diff leading coeff = 3 d lc):
for d in (1, 2, 3, 4):
    g, gc = poly(f'p3_{d}', d)
    az(sp.Poly(sp.expand(sh(g, 3) - g), E).nth(d - 1) - 3 * d * gc[d],
       f"3-step diff: coeff(E^{d-1}) of (g^[3]-g) = 3d*lc  (deg={d}) -> period-3 => const")
print("   => b3 = lambda3 * a3.")

# gauge D -> D - lambda X preserves [.,.]=1, band 3, membership; realizes b3 = 0.
Xg = {k: poly(f'X{k+3}', 3)[0] for k in range(-3, 4)}
Dg = {l: poly(f'D{l+3}', 3)[0] for l in range(-3, 4)}
Dgau = {l: sp.expand(Dg[l] - lam * Xg[l]) for l in range(-3, 4)}
for m in range(-6, 7):
    az(Qm(Dgau, Xg, m) - (Qm(Dg, Xg, m) - lam * Qm(Xg, Xg, m)),
       f"[D-lam X, X]_{m} = [D,X]_{m} - lam [X,X]_{m}")
    az(Qm(Xg, Xg, m), f"[X,X]_{m} = 0")
for r in (1, 2, 3):
    am = sp.expand(falling(r) * poly(f'am{r}', 2)[0])
    bm = sp.expand(falling(r) * poly(f'bm{r}', 2)[0])
    az(sp.rem(sp.Poly(sp.expand(bm - lam * am), E), sp.Poly(falling(r), E)).as_expr(),
       f"E^underline({r}) | (b_-{r} - lam a_-{r})  (gauge preserves membership)")
print("   => work in the gauge b_3 = 0.")


# =====================================================================
# 2. THE QUANTUM WALL Q_5.
# =====================================================================
print("\n--- 2. THE WALL Q_5 = b2^[3] a3 - a3^[2] b2 (gauge b3=0) ---")

# (2a) shifted-cube SUFFICIENCY at symbolic degree: a3 = c h h^[1] h^[2], b2 = k h h^[1].
for dh in (1, 2, 3):
    hh, _ = poly(f'h{dh}', dh)
    cc, kk = sp.symbols('cc_ kk_')
    a3c = sp.expand(cc * hh * sh(hh, 1) * sh(hh, 2))
    b2c = sp.expand(kk * hh * sh(hh, 1))
    az(sh(b2c, 3) * a3c - sh(a3c, 2) * b2c,
       f"shifted-cube (deg h={dh}): b2 = k h h^[1] solves the wall")
print("   => shifted cube a3 = c h h^[1] h^[2] is ALWAYS a wall solution.")


# (2b) root-multiset (shift-operator) calculus:  wall  <=>  S(1+S) A = (1+S+S^2) B.
def rootmult(pol):
    p = sp.Poly(sp.expand(pol), E)
    r = {}
    for root, mult in sp.roots(p).items():
        if not root.is_integer:
            raise ValueError("test data must have integer roots (single coset)")
        r[int(root)] = r.get(int(root), 0) + mult
    return r


def op_apply(coeffs, ashift):
    out = {}
    for s, sc in ashift.items():
        for n, cx in coeffs.items():
            out[n + s] = out.get(n + s, 0) + sc * cx
    return {n: v for n, v in out.items() if v != 0}


def op_div_phi3(coeffs):
    """Divide the Laurent multiset by 1+S+S^2; None if not divisible."""
    if not coeffs:
        return {}
    mn = min(coeffs)
    p = sum(cx * S**(n - mn) for n, cx in coeffs.items())
    q, r = sp.div(sp.expand(p), 1 + S + S**2, S)
    if sp.expand(r) != 0:
        return None
    qp = sp.Poly(q, S)
    return {i + mn: qp.nth(i) for i in range(qp.degree() + 1) if qp.nth(i) != 0}


def wall_B_from_A(a3):
    return op_div_phi3(op_apply(rootmult(a3), {1: 1, 2: 1}))


def poly_from_mult(mult, var=E):
    return sp.expand(sp.prod([(var - n)**m for n, m in mult.items()])) if mult else sp.Integer(1)


# verify: for a valid (Phi3-div & nonneg) A, the reconstructed b2 actually solves the wall.
for a3 in [sp.expand(E * (E - 5) * (E + 1) * (E - 4) * (E - 3) * (E - 8)),  # shifted cube h=E(E-5)(...)? generic
           sp.expand(E * (E - 1) * (E - 2)),          # consecutive cube (h=E)
           sp.expand(E * (E - 2) * (E - 4))]:         # the counterexample top
    B = wall_B_from_A(a3)
    if B is not None and all(v > 0 for v in B.values()):
        b2 = poly_from_mult(B)
        az(sh(b2, 3) * a3 - sh(a3, 2) * b2,
           f"multiset B=S(1+S)/(1+S^2+S) A reconstructs a wall solution b2 for a3 deg {sp.Poly(a3,E).degree()}")

# (2c) THE COUNTEREXAMPLE: a3 = E(E-2)(E-4) solves the wall but is NOT a shifted cube.
a3_ce = sp.expand(E * (E - 2) * (E - 4))
b2_ce = sp.expand((E - 1) * (E - 4))
az(sh(b2_ce, 3) * a3_ce - sh(a3_ce, 2) * b2_ce,
   "COUNTEREXAMPLE: a3=E(E-2)(E-4), b2=(E-1)(E-4) solves the wall b2^[3]a3=a3^[2]b2")
# a3_ce is NOT a shifted cube of ANY degree-1 h (only possible degree since deg a3=3):
gg, ggc = poly('gg', 1)
csq = sp.symbols('csq')
sols = sp.solve(sp.Poly(sp.expand(csq * gg * sh(gg, 1) * sh(gg, 2) - a3_ce), E).all_coeffs(),
                ggc + [csq], dict=True)
istrue(len(sols) == 0,
       "a3=E(E-2)(E-4) is NOT a shifted cube c*h h^[1] h^[2]  => shifted-cube conjecture REFUTED")

# (2d) 1-dimensionality of the wall solution space (freedom = one scalar kappa2).
h6 = E * (E - 5)
a3_cube6 = sp.expand(h6 * sh(h6, 1) * sh(h6, 2))
b2v, b2c = poly('bb', 6)
sol = list(sp.linsolve(sp.Poly(sp.expand(sh(b2v, 3) * a3_cube6 - sh(a3_cube6, 2) * b2v), E).all_coeffs(), b2c))[0]
nfree = len({s for x in sol for s in x.free_symbols})
istrue(nfree == 1, "wall solution space for cube a3 (within b2 deg<=6) is 1-dimensional")
b2v2, b2c2 = poly('cc', 4)
sol2 = list(sp.linsolve(sp.Poly(sp.expand(sh(b2v2, 3) * a3_ce - sh(a3_ce, 2) * b2v2), E).all_coeffs(), b2c2))[0]
nfree2 = len({s for x in sol2 for s in x.free_symbols})
istrue(nfree2 == 1, "wall solution space for exotic a3 (within b2 deg<=4) is 1-dimensional")

# (2e) degree law 2 deg a3 = 3 deg b2 : leading E^{p+q} cancels, E^{p+q-1} coeff = (3q-2p) lc.
p_, q_ = sp.symbols('p q', positive=True, integer=True)
af, afc = poly('af', 4)     # a3 generic deg 4 stand-in (p=4)
bf, bfc = poly('bf', 3)     # b2 generic deg 3 stand-in (q=3)
wall_gen = sp.expand(sh(bf, 3) * af - sh(af, 2) * bf)
pp = sp.Poly(af, E).degree(); qq = sp.Poly(bf, E).degree()
az(sp.Poly(wall_gen, E).nth(pp + qq), "wall: top coeff E^{p+q} cancels identically")
lc_a = afc[pp]; lc_b = bfc[qq]
az(sp.Poly(wall_gen, E).nth(pp + qq - 1) - (3 * qq - 2 * pp) * lc_a * lc_b,
   "wall: coeff(E^{p+q-1}) = (3 deg b2 - 2 deg a3) lc(a3) lc(b2)  => 2 deg a3 = 3 deg b2")

# (2f) nonnegativity is a GENUINE extra condition: a3 with roots {1,5,6} is Phi3-divisible,
#      but the forced B has a negative coefficient, so NO nonzero b2 exists.
a3_bad = sp.expand((E - 1) * (E - 5) * (E - 6))
Bbad = wall_B_from_A(a3_bad)
istrue(Bbad is not None, "a3 with roots {1,5,6}: root-multiset IS divisible by 1+S+S^2 (Phi3)")
istrue(any(v < 0 for v in Bbad.values()),
       "a3 roots {1,5,6}: forced B = S(1+S)/Phi3 * A has a NEGATIVE coefficient (not a multiset)")
for d in (4, 6, 8):
    bv, bcf = poly(f'bad{d}', d)
    s = list(sp.linsolve(sp.Poly(sp.expand(sh(bv, 3) * a3_bad - sh(a3_bad, 2) * bv), E).all_coeffs(), bcf))[0]
    istrue(all(x == 0 for x in s),
           f"a3 roots {{1,5,6}}: bounded solve has only b2=0 (b2 deg<={d}); Phi3-div is not sufficient")


# =====================================================================
# 3. Q_0 telescoping: closed-form local potential G, Q_0 = (T-1) G.
# =====================================================================
print("\n--- 3. Q_0 = (T-1) G  (closed-form local potential; band-k) ---")


def Gpot(a, b, K=3):
    g = 0
    for k in range(1, K + 1):
        for j in range(0, k):
            g += sh(a[k], j - k) * sh(b[-k], j) - sh(b[k], j - k) * sh(a[-k], j)
    return sp.expand(g)


Ggen = Gpot(A, Bc)
az((sh(Ggen, 1) - Ggen) - Qm(A, Bc, 0),
   "Q_0 = G^[1] - G  with  G = sum_{k=1}^{3} sum_{j=0}^{k-1} (a_k^[j-k] b_-k^[j] - b_k^[j-k] a_-k^[j])")

# reduces to the band-2 potential (K=2, gauge b2=0, a2=h h^[1], a1=h p, b1=kappa h):
hb, _ = poly('hb', 2); pb, _ = poly('pb', 2)
ub, _ = poly('ub', 2); vb, _ = poly('vb', 2); wb, _ = poly('wb', 2)
kap = sp.symbols('kap')
a2 = {2: sp.expand(hb * sh(hb, 1)), 1: sp.expand(hb * pb), 0: poly('a0b', 2)[0], -1: ub, -2: sp.expand(E * (E - 1) * poly('sb', 1)[0])}
b2 = {2: sp.Integer(0), 1: sp.expand(kap * hb), 0: poly('b0b', 2)[0], -1: vb, -2: wb}
G2 = Gpot(a2, b2, K=2)
G2_memo = sp.expand(sh(hb, -1) * (hb * sh(wb, 1) + sh(hb, -2) * wb + sh(pb, -1) * vb - kap * ub))
az(G2 - G2_memo,
   "band-2 specialization: G = h^[-1]( h w^[1] + h^[-2] w + p^[-1] v - kappa u )  (matches quantum-completion.md)")

# G(0) = 0 under genuine membership => Q_0=1 gives G = E (the constant is 0).
am1 = sp.expand(E * poly('u1', 2)[0]); am2 = sp.expand(E * (E - 1) * poly('u2', 2)[0])
am3 = sp.expand(E * (E - 1) * (E - 2) * poly('u3', 2)[0])
bm1 = sp.expand(E * poly('w1', 2)[0]); bm2 = sp.expand(E * (E - 1) * poly('w2', 2)[0])
bm3 = sp.expand(E * (E - 1) * (E - 2) * poly('w3', 2)[0])
aM = {3: poly('a3', 2)[0], 2: poly('a2', 2)[0], 1: poly('a1', 2)[0], 0: poly('a0', 2)[0], -1: am1, -2: am2, -3: am3}
bM = {3: sp.Integer(0), 2: poly('b2', 2)[0], 1: poly('b1', 2)[0], 0: poly('b0', 2)[0], -1: bm1, -2: bm2, -3: bm3}
az(Gpot(aM, bM).subs(E, 0),
   "G(0) = 0 identically under genuine membership  => Q_0=(T-1)G=1 forces G = E + c with c=0")


# =====================================================================
# 4. Bottom mirror + cross-coupling of proportionality constants.
# =====================================================================
print("\n--- 4. bottom mirror Q_-6/Q_-5; reflection; cross-coupling ---")
# Q_-6 => b_-3 = mu3 a_-3 (same 3-periodicity/rational argument, reflected).
az(sp.simplify((sh(rr, -3) - rr) - (sh(b3g, -3) * a3g - sh(a3g, -3) * b3g) / (a3g * sh(a3g, -3))),
   "rational (-3)-periodicity engine: r^[-3]-r = (b-3^[-3]a-3 - a-3^[-3]b-3)/(a-3 a-3^[-3])")
print("   => b_-3 = mu3 a_-3  (an INDEPENDENT scalar; gauge already spent on the top).")

# Reflection E -> -E-1 (quantum Fourier phi, cf. quantum-a2-zero.md) flips shift signs,
# carrying the bottom wall to the top-wall form (so a_-3 obeys the SAME Phi3 structure).
Pt, _ = poly('Pt', 3)
Pref = sp.expand(Pt.subs(E, -E - 1))
az(sh(Pt, -3).subs(E, -E - 1) - sh(Pref, 3),
   "reflection E->-E-1: (P^[-3])|refl = (P|refl)^[+3]  (bottom shift -3 <-> top shift +3)")
az(sh(Pt, -2).subs(E, -E - 1) - sh(Pref, 2),
   "reflection E->-E-1: (P^[-2])|refl = (P|refl)^[+2]")
print("   => bottom wall Q_-5 forces the Phi3 root-structure on a_-3, mirror of the top.")

# CROSS-COUPLING (the first genuinely new band-3 phenomenon): after the SINGLE top gauge
# b3=0, the bottom keeps b_-3 = mu3 a_-3 with mu3 free; lambda3 (spent) and mu3 are
# independent, and Q_-5, Q_-4 are INHOMOGENEOUS in b_-2,b_-1 with mu3-driven source terms.
Bcp = dict(Bc); Bcp[3] = sp.Integer(0); Bcp[-3] = sp.expand(mu3 * A[-3])
qm5_cp = Qm(A, Bcp, -5)
az(qm5_cp - ((sh(Bcp[-2], -3) * A[-3] - sh(A[-3], -2) * Bcp[-2])
             + mu3 * (sh(A[-3], -2) * A[-2] - sh(A[-2], -3) * A[-3])),
   "Q_-5 raw decomposition has the mu3 cross-coupling source")
phi_m2 = sp.expand(Bcp[-2] - mu3 * A[-2])
az(qm5_cp - (sh(phi_m2, -3) * A[-3] - sh(A[-3], -2) * phi_m2),
   "Q_-5 = phi^[-3] a_-3 - a_-3^[-2] phi for phi=b_-2-mu3*a_-2 (minus sign verified)")


# =====================================================================
# 5. Negative-tail staggered rigidity lemmas (Lemma-R family).
# =====================================================================
print("\n--- 5. Lemma-R staggered rigidities (negative tail) ---")
# General staggered leading-coefficient fact: for f^[a] g^[b] - f^[a'] g^[b'] with equal
# top terms (deg f = p, deg g = q), coeff(E^{p+q-1}) = ((a-a')p + (b-b')q) lc(f) lc(g).
# We verify it on the EXACT shift-pairs of the memo's Lemma-R table, and read off each
# forced degree relation from the obstruction coefficient ((a-a')p + (b-b')q) = 0.
ff, ffc = poly('ff', 4); gg2, gg2c = poly('gg2', 3)
pf = 4; qf = 3
table = [
    ("gauged bottom wall Q_-5 [a-3^[0] phi^[-3] - a-3^[-2] phi^[0]]", 0, -3, -2, 0, "2 deg a-3 = 3 deg phi"),
    ("Q_-4 term         [a-3^[0] b-1^[-3] - a-3^[-1] b-1^[0]]", 0, -3, -1, 0, "  deg a-3 = 3 deg b-1"),
    ("Q_-4 term (mu3)   [a-3^[-1] a-1^[0] - a-3^[0] a-1^[-3]]", -1, 0, 0, -3, "  deg a-3 = 3 deg a-1"),
]
for (name, a, bsh, ap, bp, rel) in table:
    diff = sp.expand(sh(ff, a) * sh(gg2, bsh) - sh(ff, ap) * sh(gg2, bp))
    obstruct = (a - ap) * pf + (bsh - bp) * qf
    az(sp.Poly(diff, E).nth(pf + qf), f"staggered top E^{{p+q}} cancels: {name}")
    az(sp.Poly(diff, E).nth(pf + qf - 1) - obstruct * ffc[pf] * gg2c[qf],
       f"staggered coeff(E^{{p+q-1}}) = (({a}-{ap})p+({bsh}-{bp})q) lc lc  =>  forces {rel.strip()}")
print("   Lemma-R degree-forcings (each fires when its term dominates):")
for (name, a, bsh, ap, bp, rel) in table:
    print(f"     {name}:  {rel.strip()}")
print("     top wall Q_5      [a3^[0] b2^[3] - a3^[2] b2^[0]]:  2 deg a3 = 3 deg b2  (verify sec 2)")

# Q_-4 decomposition (documents the middle 2-step 'square classification one rung down').
qm4_cp = Qm(A, Bcp, -4)
az(qm4_cp - ((sh(Bcp[-1], -3) * A[-3] - sh(A[-3], -1) * Bcp[-1])
             + mu3 * (sh(A[-3], -1) * A[-1] - sh(A[-1], -3) * A[-3])
             + (sh(Bcp[-2], -2) * A[-2] - sh(A[-2], -2) * Bcp[-2])),
   "Q_-4 = [b-1,a-3 staggered] + mu3[a-3,a-1 staggered] + [b-2^[-2]a-2 - a-2^[-2]b-2]")
# the middle piece isolated is the J2q 'one rung down': b-2/a-2 2-periodic => b-2 = nu a-2.
mid2, _ = poly('mid2', 2); a2mid, _ = poly('a2mid', 2)
az(sp.simplify((sh(mid2 / a2mid, 2) - mid2 / a2mid)
               - (sh(mid2, 2) * a2mid - sh(a2mid, 2) * mid2) / (a2mid * sh(a2mid, 2))),
   "middle 2-step engine: (b-2/a-2)^[2] - (b-2/a-2) = (b-2^[2]a-2 - a-2^[2]b-2)/(a-2 a-2^[2])")


# =====================================================================
# 6. Positive control: a GENUINE band-3 pair.
# =====================================================================
print("\n--- 6. positive control: genuine band-3 pair ---")
# U = x + c1 d,  X = U^3 - d/kappa,  D = kappa U.  [D,X] = kappa[U,X] = kappa*(1/kappa) = 1.
def mul_ladders(P, Q):
    R = {}
    for k1, p1 in P.items():
        for k2, p2 in Q.items():
            R[k1 + k2] = R.get(k1 + k2, 0) + sh(p1, k2) * p2
    return {k: sp.expand(v) for k, v in R.items() if sp.expand(v) != 0}


U = {1: sp.Integer(1), -1: sp.expand(c1 * E)}          # x + c1 d,  d = x^{-1} E
U3 = mul_ladders(mul_ladders(U, U), U)
Xw = {k: sp.Integer(0) for k in range(-3, 4)}
for k, v in U3.items():
    Xw[k] += v
Xw[-1] = sp.expand(Xw[-1] - E / kappa)                  # - d/kappa = - x^{-1}E/kappa
Dw = {k: sp.Integer(0) for k in range(-3, 4)}
Dw[1] = kappa; Dw[-1] = sp.expand(kappa * c1 * E)       # kappa U
for m in range(-6, 7):
    # [D,X] with D=Dw, X=Xw:
    az(Qm(Xw, Dw, m) - (1 if m == 0 else 0),
       f"witness: [D,X]_{m} = delta_{{m0}}  (X=U^3-d/kappa, D=kappa U)")
istrue(divides(falling(1), Xw[-1]) and divides(falling(2), Xw[-2]) and divides(falling(3), Xw[-3]),
       "witness X memberships: E|a-1, E(E-1)|a-2, E(E-1)(E-2)|a-3")
istrue(sp.expand(Xw[3] - 1) == 0, "witness a_3 = 1 (trivial shifted cube h=1)")
istrue(sp.expand(Xw[-3] - c1**3 * E * (E - 1) * (E - 2)) == 0, "witness a_-3 = c1^3 E(E-1)(E-2)")


# =====================================================================
# 7. Trivial-cube positive middle; bounded corroboration.
# =====================================================================
print("\n--- 7. trivial-cube positive middle; bounded corroboration ---")
# When a3 = 1 (normalized constant), the wall Q5 = a3(b2^[3]-b2) forces b2 constant,
# and Q4 collapses to the band-3 positive cascade with a 3-fold (not 2-fold) periodicity.
k2c = sp.symbols('kappa2')
Atc = {3: sp.Integer(1), 2: poly('t2', 2)[0], 1: poly('t1', 2)[0], 0: poly('t0', 2)[0],
       -1: poly('tm1', 2)[0], -2: poly('tm2', 2)[0], -3: poly('tm3', 2)[0]}
Btc = {3: sp.Integer(0), 2: k2c, 1: poly('s1', 2)[0], 0: poly('s0', 2)[0],
       -1: poly('sm1', 2)[0], -2: poly('sm2', 2)[0], -3: poly('sm3', 2)[0]}
az(Qm(Atc, Btc, 5), "trivial cube a3=1: Q5 = a3(b2^[3]-b2)=0 forces b2 const (b2=kappa2 gives Q5=0)")
az(Qm(Atc, Btc, 4) - (k2c * (Atc[2] - sh(Atc[2], 2)) + (sh(Btc[1], 3) - Btc[1])),
   "trivial cube a3=1, b2=kappa2: Q4 = kappa2(a2 - a2^[2]) + (b1^[3] - b1)  (3-fold periodicity)")

# All shifted cubes of degree-<=1 h give wall solutions; enumerate a small grid.
cnt = 0
for r0 in range(-2, 3):
    hlin = E - r0
    a3cb = sp.expand(hlin * sh(hlin, 1) * sh(hlin, 2))
    b2cb = sp.expand(hlin * sh(hlin, 1))
    if sp.expand(sh(b2cb, 3) * a3cb - sh(a3cb, 2) * b2cb) != 0:
        raise AssertionError("shifted cube failed wall (grid)")
    cnt += 1
print(f"PASS {cnt} degree-1 shifted cubes all solve the wall (grid corroboration)")

# A modest sweep: which 3-root integer tops {0,a,b} admit a b2 (Phi3-div + nonneg)?
admit = []
for aa in range(1, 6):
    for bb in range(aa + 1, 8):
        a3t = sp.expand(E * (E - aa) * (E - bb))
        try:
            Bt = wall_B_from_A(a3t)
        except ValueError:
            Bt = None
        if Bt is not None and all(v > 0 for v in Bt.values()):
            admit.append((0, aa, bb))
print("PASS 3-root tops {0,a,b} admitting a wall b2 (corroboration):", admit)
istrue((0, 2, 4) in admit, "the counterexample top {0,2,4} is among wall-admitting non-cubes")


print("\nALL QUANTUM BAND3 CHECKS PASSED")
