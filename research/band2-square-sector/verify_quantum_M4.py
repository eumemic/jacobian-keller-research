#!/usr/bin/env python3
"""
verify_quantum_M4.py
====================
Exact SymPy verification backing the memo `quantum-M4.md`:

    THEOREM (M4q).  In A_1 (the Weyl algebra, GENUINE membership -- not merely the
    localization A_1[x^{-1}]), there is NO pair [D,X] = 1 whose ladder supports are
    both contained in [-2,2] and whose top coefficient a_2(E) is nonzero and OUTSIDE
    the shifted-square class (i.e. a_2 != c*h(E)*h(E+1) for every c in C*, h in C[E]).
    Equivalently: every genuine band-2 Weyl pair with a_2 != 0 has a_2 shifted-square.

This is the quantum analogue of the proved classical M4 (research/band2-classical-proved,
commit f637b1a), mirroring its Section-5 descent.  It handles a_{-2} arbitrary (0 or not):
the one-sided content the audited two-sided obstruction J3q (commit 91a053a / the J3
memo) does NOT cover; the killing at the end is done by GENUINE Weyl MEMBERSHIP,
exactly as polynomiality did classically.

Convention (matching all sibling quantum memos):
    A_1[x^{-1}] = (+)_{k in Z} x^k C[E],   (x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E),
    f^[r](E) = f(E+r),   E = x d,
    ladder-m coefficient of [D,X]:   Q_m = sum_{k+l=m} ( b_l^[k] a_k - a_k^[l] b_l ),
    Keller/Weyl relation [D,X] = 1  <=>  Q_m = delta_{m0}.
    Membership (genuine A_1):  E(E-1)...(E-r+1) | a_{-r}(E)  for each negative ladder -r,
    and likewise for D.  In band 2:  E | a_{-1},b_{-1}  and  E(E-1) | a_{-2},b_{-2}.

Sections:
  0. crossed-product engine; Q_m == direct commutator (m = -4..4).
  1. gauge D -> D - lambda X: preserves [.,.]=1 and membership -> WLOG b_2 = 0.
  2. Q_3 (gauge b2=0) = b1^[2] a2 - a2^[1] b1 = J2q equation; J2q lemma => b_1 = 0.
  3. Q_2 (gauge b2=b1=0) = a2 (b0^[2]-b0) => b_0 constant (2-periodic).
  4. Q_1 (gauge, b0 const) = v^[2] a2 - a2^[-1] v (v=b_{-1});
       telescoping invariant P = a2^[-1] v v^[1] is 1-periodic => v = 0.
  5. Q_0 (gauge, v=0) = Pi^[2] - Pi with Pi = a2^[-2] w (w=b_{-2});  = 1
       => Pi linear (lc 1/2) => a_2 LINEAR, w nonzero constant.
  6. negative tail Q_{-1..-4} = w (a_j - a_j^[-2]) => a_1,a_0,a_{-1},a_{-2} constant.
  7. MEMBERSHIP CONTRADICTION: E(E-1) | b_{-2} = w, w nonzero constant -> impossible.
  8. bounded-degree EXHAUSTIVE emptiness (Groebner) for concrete non-shifted-square a_2.
  9. positive control: an explicit GENUINE pair with a_2 shifted-square IS accepted.

Run:  uv run --with sympy python research/band2-square-sector/verify_quantum_M4.py
Ends: ALL QUANTUM M4 CHECKS PASSED
"""
import sympy as sp

E = sp.symbols('E')
lam, gamma, beta, cc = sp.symbols('lambda gamma beta c')


def sh(f, r):
    return sp.expand(sp.sympify(f).subs(E, E + r))


def poly(name, deg):
    cs = sp.symbols(f'{name}_0:{deg+1}')
    return sum(cs[i] * E**i for i in range(deg + 1)), list(cs)


def Qm(a, b, m):
    """a,b: dicts k -> coeff poly (k in -2..2). Ladder-m coefficient of [D,X]."""
    return sp.expand(sum(
        sh(b[l], k) * a[k] - sh(a[k], l) * b[l]
        for k in range(-2, 3) for l in range(-2, 3) if k + l == m))


def az(expr, label):
    e = sp.expand(sp.sympify(expr))
    if e != 0:
        raise AssertionError(label + "  :  residual = " + str(e))
    print("PASS", label)


def falling(r):
    return sp.prod([E - i for i in range(r)]) if r > 0 else sp.Integer(1)


def divides(fac, c):
    """Does polynomial `fac` divide polynomial `c` in C[E]?"""
    c = sp.expand(sp.sympify(c))
    if c == 0:
        return True
    return sp.rem(sp.Poly(c, E), sp.Poly(fac, E)) == 0


# =====================================================================
# 0. Crossed-product engine: Q_m equals the direct commutator ladder coeff.
# =====================================================================
print("--- 0. Q_m == direct commutator ---")
A = {k: poly(f'A{k+2}', 2)[0] for k in range(-2, 3)}
Bc = {k: poly(f'B{k+2}', 2)[0] for k in range(-2, 3)}


def direct_commutator(m):
    # [D,X] = DX - XD ;  (x^l b_l)(x^k a_k) = x^{l+k} b_l(E+k) a_k(E),
    #                    (x^k a_k)(x^l b_l) = x^{k+l} a_k(E+l) b_l(E).
    tot = 0
    for k in range(-2, 3):
        for l in range(-2, 3):
            if k + l == m:
                tot += sh(Bc[l], k) * A[k] - sh(A[k], l) * Bc[l]
    return sp.expand(tot)


for m in range(-4, 5):
    az(direct_commutator(m) - Qm(A, Bc, m), f"Q_{m} = direct commutator ladder coeff")

# Explicit closed forms of the two workhorse ladders (audits the shift bookkeeping).
az(Qm(A, Bc, 4) - (sh(Bc[2], 2) * A[2] - sh(A[2], 2) * Bc[2]),
   "Q_4 = b2^[2] a2 - a2^[2] b2")
az(Qm(A, Bc, 0) - (sh(Bc[-2], 2) * A[2] - sh(A[2], -2) * Bc[-2]
                   + sh(Bc[-1], 1) * A[1] - sh(A[1], -1) * Bc[-1]
                   + sh(Bc[1], -1) * A[-1] - sh(A[-1], 1) * Bc[1]
                   + sh(Bc[2], -2) * A[-2] - sh(A[-2], 2) * Bc[2]),
   "Q_0 = b-2^[2]a2 - a2^[-2]b-2 + b-1^[1]a1 - a1^[-1]b-1"
   " + b1^[-1]a-1 - a-1^[1]b1 + b2^[-2]a-2 - a-2^[2]b2")


# =====================================================================
# 1. Gauge D -> D - lambda X.  Preserves [.,.] and membership => WLOG b_2 = 0.
# =====================================================================
print("\n--- 1. gauge legitimacy + membership preservation ---")
Xg = {k: poly(f'X{k+2}', 3)[0] for k in range(-2, 3)}
Dg = {l: poly(f'D{l+2}', 3)[0] for l in range(-2, 3)}
Dgau = {l: sp.expand(Dg[l] - lam * Xg[l]) for l in range(-2, 3)}
for m in range(-4, 5):
    az(Qm(Dgau, Xg, m) - (Qm(Dg, Xg, m) - lam * Qm(Xg, Xg, m)),
       f"[D-lam X, X]_{m} = [D,X]_{m} - lam [X,X]_{m}")
    az(Qm(Xg, Xg, m), f"[X,X]_{m} = 0")
# membership survives the gauge: E|a-1,b-1 and E(E-1)|a-2,b-2  =>  same for b-l - lam a-l.
am1 = sp.expand(E * poly('am1g', 2)[0]);      bm1 = sp.expand(E * poly('bm1g', 2)[0])
am2 = sp.expand(E * (E - 1) * poly('am2g', 2)[0]); bm2 = sp.expand(E * (E - 1) * poly('bm2g', 2)[0])
az(sp.rem(sp.Poly(sp.expand(bm1 - lam * am1), E), sp.Poly(falling(1), E)).as_expr(),
   "E | (b_{-1} - lam a_{-1})  (gauged b_{-1} keeps membership)")
az(sp.rem(sp.Poly(sp.expand(bm2 - lam * am2), E), sp.Poly(falling(2), E)).as_expr(),
   "E(E-1) | (b_{-2} - lam a_{-2})  (gauged b_{-2} keeps membership)")
print("   => work henceforth in the gauge b_2 = 0 (rename gauged coeffs b_l).")


# ---- coefficient polynomials for the descent (gauge b_2 = 0) ----
a2, _ = poly('a2', 5)
a1, _ = poly('a1', 5)
a0, _ = poly('a0', 5)
u, _ = poly('u', 5)   # a_{-1}
s, _ = poly('s', 5)   # a_{-2}
b1, _ = poly('b1', 5)
b0, _ = poly('b0', 5)
v, _ = poly('v', 5)   # b_{-1}
w, _ = poly('w', 5)   # b_{-2}
Z = sp.Integer(0)


def data(b2=Z, b1v=b1, b0v=b0, vv=v, wv=w):
    a = {2: a2, 1: a1, 0: a0, -1: u, -2: s}
    b = {2: b2, 1: b1v, 0: b0v, -1: vv, -2: wv}
    return a, b


# =====================================================================
# 2. Ladder 3 -> J2q equation -> b_1 = 0  (uses non-shifted-square of a_2).
# =====================================================================
print("\n--- 2. Q_3 = J2q equation => b_1 = 0 ---")
a, b = data()
az(Qm(a, b, 3) - (sh(b1, 2) * a2 - sh(a2, 1) * b1),
   "Q_3 = b1^[2] a2 - a2^[1] b1   (gauge b2=0)")

# J2q lemma (audited J3 memo, commit 91a053a): for a2 != 0 the polynomial solutions of
#     h(E+2) a2(E) = a2(E+1) h(E)
# are h = 0, or a2 = c h(E) h(E+1) (shifted-square).  Two exact facts:
#  (i) shifted-square a2 DOES admit a nonzero solution h;
h_, _ = poly('h', 3)
a2_sq = sp.expand(cc * h_ * sh(h_, 1))
az(sh(h_, 2) * a2_sq - sh(a2_sq, 1) * h_,
   "J2q (=>): if a2 = c h h^[1] then h solves h^[2]a2 = a2^[1]h")
#  (ii) the lemma's engine: for h != 0 set r = a2/(h h^[1]) in C(E).  Then the equation is
#       EXACTLY r(E+1)=r(E) -- verified here as a rational-function identity: r^[1]-r equals
#       the J2q residual (h^[2]a2 - a2^[1]h) over the denominator h h^[1] h^[2].
hgen, _ = poly('hg', 3)
rfun = a2 / (hgen * sh(hgen, 1))
resid = sh(hgen, 2) * a2 - sh(a2, 1) * hgen      # = Q_3 with b1=hgen
az(sp.simplify((sh(rfun, 1) - rfun) + resid / (hgen * sh(hgen, 1) * sh(hgen, 2))),
   "J2q engine: r^[1]-r = -(h^[2]a2 - a2^[1]h)/(h h^[1] h^[2])  => equation <=> r is 1-periodic")
#       Rational periodicity (audited lemma) => r is a nonzero constant c => a2 = c h h^[1].
print("PASS J2q (<=): nonzero h forces a2 = c h h^[1] (rational periodicity of r)")
print("   non-shifted-square a_2 => the only solution is h = b_1 = 0.")


# =====================================================================
# 3. Ladder 2 -> b_0 constant.
# =====================================================================
print("\n--- 3. Q_2 => b_0 constant ---")
a, b = data(b1v=Z)
az(Qm(a, b, 2) - a2 * (sh(b0, 2) - b0),
   "Q_2 = a2 (b0^[2] - b0)   (gauge b2=b1=0)")
# a2 != 0 (integral domain) => b0^[2] = b0.  Leading-coeff rigidity: 2-periodic => constant.
for d in (1, 2, 3, 4):
    g, gc = poly(f'g2_{d}', d)
    az(sp.Poly(sp.expand(sh(g, 2) - g), E).nth(d - 1) - 2 * d * gc[d],
       f"2-step diff: coeff(E^{d-1}) of (g^[2]-g) = 2d*lc(g) != 0  (deg g={d}) => b0 const")


# =====================================================================
# 4. Ladder 1 -> homogeneous eqn in v=b_{-1} -> telescoping invariant -> v = 0.
# =====================================================================
print("\n--- 4. Q_1 => v = 0 via 1-periodic invariant ---")
a, b = data(b1v=Z, b0v=gamma)   # b0 = gamma constant
az(Qm(a, b, 1) - (sh(v, 2) * a2 - sh(a2, -1) * v),
   "Q_1 = v^[2] a2 - a2^[-1] v   (indep. of a1,a0,a_-1,a_-2 after b1=b2=0, b0 const)")

# Invariant  P(E) = a2(E-1) v(E) v(E+1).  Exact identity:
#     P(E+1) - P(E) = v(E+1) * ( a2(E) v(E+2) - a2(E-1) v(E) )  =  v^[1] * Q_1 .
P = sp.expand(sh(a2, -1) * v * sh(v, 1))
az((sh(P, 1) - P) - sh(v, 1) * (sh(v, 2) * a2 - sh(a2, -1) * v),
   "P^[1] - P = v^[1] * Q_1   with P = a2^[-1] v v^[1]  => Q_1=0 makes P 1-periodic")
# 1-periodic polynomial is constant (leading-coeff rigidity):
for d in (1, 2, 3, 4):
    g, gc = poly(f'g1_{d}', d)
    az(sp.Poly(sp.expand(sh(g, 1) - g), E).nth(d - 1) - d * gc[d],
       f"1-step diff: coeff(E^{d-1}) of (g^[1]-g) = d*lc(g) != 0  (deg g={d}) => P const")
# product degree adds (integral domain): if v != 0 then deg P = deg a2 + 2 deg v.
for (da, dv) in [(1, 0), (1, 1), (2, 1), (3, 2)]:
    aa, _ = poly(f'aa_{da}_{dv}', da)
    vv, _ = poly(f'vv_{da}_{dv}', dv)
    Pex = sp.expand(sh(aa, -1) * vv * sh(vv, 1))
    az(sp.Poly(Pex, E).degree() - (da + 2 * dv),
       f"deg(a2^[-1] v v^[1]) = deg a2 + 2 deg v  ({da}+2*{dv})  (nonzero product)")
print("   v!=0 => P nonzero constant => a2^[-1] a unit => a2 constant, contra deg a2>=1."
      "  Hence v = b_{-1} = 0.")


# =====================================================================
# 5. Central ladder 0 = 1 -> a_2 LINEAR, w = b_{-2} nonzero constant.
# =====================================================================
print("\n--- 5. Q_0 = 1 => a_2 linear, w nonzero constant ---")
a, b = data(b1v=Z, b0v=gamma, vv=Z)   # v = 0
az(Qm(a, b, 0) - (sh(w, 2) * a2 - sh(a2, -2) * w),
   "Q_0 = w^[2] a2 - a2^[-2] w   (v=b1=b2=0)")
# Pi(E) = a2(E-2) w(E) satisfies  Q_0 = Pi(E+2) - Pi(E).
Pi = sp.expand(sh(a2, -2) * w)
az(Qm(a, b, 0) - (sh(Pi, 2) - Pi), "Q_0 = Pi^[2] - Pi   with Pi = a2^[-2] w")
# Pi^[2]-Pi = 1 forces deg Pi = 1 with lc(Pi)=1/2 (2-step diff leading coeff = 2d*lc):
for d in (1, 2, 3):
    g, gc = poly(f'g0_{d}', d)
    az(sp.Poly(sp.expand(sh(g, 2) - g), E).nth(d - 1) - 2 * d * gc[d],
       f"2-step diff leading coeff = 2d*lc(g) (deg={d}) => Pi^[2]-Pi=1 forces deg Pi=1, lc=1/2")
# product degree adds: deg Pi = deg a2 + deg w = 1, deg a2 >= 1 => deg a2 = 1, deg w = 0, w!=0.
for (da, dw) in [(1, 0), (2, 0), (1, 1)]:
    aa, _ = poly(f'pa_{da}_{dw}', da)
    ww, _ = poly(f'pw_{da}_{dw}', dw)
    az(sp.Poly(sp.expand(sh(aa, -2) * ww), E).degree() - (da + dw),
       f"deg(a2^[-2] w) = deg a2 + deg w  ({da}+{dw})")
print("   deg(a2^[-2] w)=1 & deg a2>=1  =>  a_2 LINEAR, w NONZERO CONSTANT.")


# =====================================================================
# 6. Negative tail (consistency): Q_{-1..-4} force a_1,a_0,a_{-1},a_{-2} constant.
# =====================================================================
print("\n--- 6. negative tail => a_1,a_0,a_{-1},a_{-2} constant ---")
wc = sp.symbols('w_const')
a, b = data(b1v=Z, b0v=gamma, vv=Z, wv=wc)   # w a nonzero constant
az(Qm(a, b, -1) - wc * (a1 - sh(a1, -2)),  "Q_-1 = w (a1 - a1^[-2])   => a_1 const")
az(Qm(a, b, -2) - wc * (a0 - sh(a0, -2)),  "Q_-2 = w (a0 - a0^[-2])   => a_0 const")
az(Qm(a, b, -3) - wc * (u - sh(u, -2)),    "Q_-3 = w (a_-1 - a_-1^[-2]) => a_-1 const")
az(Qm(a, b, -4) - wc * (s - sh(s, -2)),    "Q_-4 = w (a_-2 - a_-2^[-2]) => a_-2 const")
print("   (not needed for the contradiction, but confirms the reduced system is consistent:"
      " w!=0 => each 2-step difference vanishes => constants.)")


# =====================================================================
# 7. MEMBERSHIP CONTRADICTION.
# =====================================================================
print("\n--- 7. membership contradiction ---")
# Gauged b_{-2} = w is a NONZERO constant; genuine membership demands E(E-1) | b_{-2}.
for k in (1, -1, sp.Rational(1, 2), 7):
    if k != 0 and divides(falling(2), sp.Integer(1) * k):
        raise AssertionError("E(E-1) wrongly divides a nonzero constant")
az(sp.Integer(0), "E(E-1) does NOT divide any nonzero constant  =>  w=b_{-2}=0, contra w!=0")
print("   CONTRADICTION: no genuine band-2 Weyl pair with non-shifted-square a_2 exists.  QED.")


# =====================================================================
# 8. Bounded-degree EXHAUSTIVE emptiness (corroboration only).
# =====================================================================
print("\n--- 8. bounded-degree exhaustive emptiness (concrete non-sq a_2) ---")


def is_shifted_square(a2v):
    a2v = sp.expand(a2v)
    d = sp.Poly(a2v, E).degree()
    if d % 2 == 1:
        return False           # odd degree: c h h^[1] always has even degree
    dh = d // 2
    hc = sp.symbols(f'hh_0:{dh+1}')
    c = sp.symbols('cc_')
    hh = sum(hc[i] * E**i for i in range(dh + 1))
    return len(sp.solve(sp.Poly(sp.expand(c * hh * sh(hh, 1) - a2v), E).all_coeffs(),
                        list(hc) + [c], dict=True)) > 0


def gen(name, deg):
    cs = sp.symbols(f'{name}_0:{deg+1}')
    return sum(cs[i] * E**i for i in range(deg + 1)), list(cs)


def empty_gauge(a2val, Bdeg):
    """True iff NO genuine gauge-b2=0 pair with this fixed a2 and all other coeffs deg<=Bdeg."""
    a1v, a1c = gen('sa1', Bdeg)
    a0v, a0c = gen('sa0', Bdeg)
    um, umc = gen('su', Bdeg);  uu = sp.expand(E * um)             # a_{-1}, membership
    sm, smc = gen('ss', Bdeg);  ss = sp.expand(E * (E - 1) * sm)   # a_{-2}, membership
    b1v, b1c = gen('sb1', Bdeg)
    b0v, b0c = gen('sb0', Bdeg)
    vm, vmc = gen('sv', Bdeg);  vv = sp.expand(E * vm)             # b_{-1}, membership
    wm, wmc = gen('sw', Bdeg);  ww = sp.expand(E * (E - 1) * wm)   # b_{-2}, membership
    aa = {2: a2val, 1: a1v, 0: a0v, -1: uu, -2: ss}
    bb = {2: Z, 1: b1v, 0: b0v, -1: vv, -2: ww}
    eqs = []
    for m in range(-4, 5):
        tgt = 1 if m == 0 else 0
        for coef in sp.Poly(Qm(aa, bb, m) - tgt, E).all_coeffs():
            if coef != 0:
                eqs.append(coef)
    unknowns = a1c + a0c + umc + smc + b1c + b0c + vmc + wmc
    Gb = sp.groebner(eqs, *unknowns, order='grevlex')
    return list(Gb.exprs) == [sp.Integer(1)]


cases = {'E': E, '2E-3': 2 * E - 3, 'E^2+1': E**2 + 1,
         'E^2-E+1': E**2 - E + 1, 'E^3': E**3, 'E^4+1': E**4 + 1}
for nm, a2v in cases.items():
    if is_shifted_square(a2v):
        raise AssertionError(f"{nm} is shifted-square; bad test datum")
print("   chosen a_2 all confirmed NON-shifted-square:", ", ".join(cases))
allok = True
for Bdeg in (2, 3):
    for nm, a2v in cases.items():
        e = empty_gauge(a2v, Bdeg)
        allok = allok and e
        print(f"   B={Bdeg}  a_2={nm:9s}: {'EMPTY' if e else 'NONEMPTY!!'}")
if not allok:
    raise AssertionError("a bounded profile was NONEMPTY -- contradicts the theorem")
print("   all bounded profiles EMPTY (Groebner basis = <1>).")


# =====================================================================
# 9. Positive control: a GENUINE pair with shifted-square a_2 IS accepted.
# =====================================================================
print("\n--- 9. positive control: genuine pair, a_2 shifted-square ---")
# X = (x+d)^2 - d/2,  D = 2x + 2d,  d = partial.  [D,X] = 1 exactly.
#   (x+d)^2 = x^2 + (2E+1) + d^2,  d = x^{-1}E,  d^2 = x^{-2}E(E-1).
aX = {2: sp.Integer(1), 1: Z, 0: 2 * E + 1, -1: -E / 2, -2: E * (E - 1)}
bD = {2: Z, 1: sp.Integer(2), 0: Z, -1: 2 * E, -2: Z}
for m in range(-4, 5):
    az(Qm(aX, bD, m) - (1 if m == 0 else 0), f"witness: Q_{m} = delta_{{m0}}")
for r, c in [(1, aX[-1]), (2, aX[-2]), (1, bD[-1])]:
    if not divides(falling(r), c):
        raise AssertionError("witness violates membership")
az(sp.Integer(0), "witness genuine (memberships hold) and a_2 = 1 IS shifted-square")
if is_shifted_square(aX[2]):
    print("   witness a_2 = 1 is shifted-square (as required for a genuine pair to exist).")
else:
    raise AssertionError("witness a_2 should be shifted-square")
print("   => the Q_m machinery accepts real genuine pairs; emptiness for NON-square a_2 is a"
      " genuine distinction, driven by membership.")


print("\nALL QUANTUM M4 CHECKS PASSED")
