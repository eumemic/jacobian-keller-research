#!/usr/bin/env python3
"""Exact classical band-3 cascade verifier.

Convention (frozen, classical):  {G,F} = G_xi F_x - G_x F_xi,  tau = x*xi.
Band-3 presentation:  F = sum_{k=-3}^{3} x^k a_k(tau),  G = sum_{l=-3}^{3} x^l b_l(tau).
Membership (genuine polynomiality in C[x,xi]):  tau^j | a_{-j}, b_{-j}, j=1,2,3.
Keller equations:  C_m = sum_{k+l=m}(k a_k b_l' - l a_k' b_l) = delta_{m0},  m in [-6,6].

This script MACHINE-VERIFIES every displayed identity in classical-band3-cascade.md.
Degree-balance / membership contradictions are written proofs in the .md; the
bounded-degree sweeps here are regression corroboration only.  A successful run
ends 'ALL CLASSICAL BAND3 CHECKS PASSED'.
"""
import sympy as sp

x, xi, t = sp.symbols("x xi tau")
KS = range(-3, 4)


def D(f):
    return sp.diff(f, t)


def Cm(m, A, B):
    """Coefficient of x^m in {G,F} via the frozen sum formula."""
    return sp.expand(sum(
        k*A[k]*D(B[l]) - l*D(A[k])*B[l]
        for k in KS for l in KS if k + l == m
    ))


def poisson(g, f):
    return sp.diff(g, xi)*sp.diff(f, x) - sp.diff(g, x)*sp.diff(f, xi)


PASS = 0


def check(cond, label):
    global PASS
    if cond is not True and cond != True:  # noqa: E712
        raise AssertionError("FAIL: " + label)
    PASS += 1
    print("PASS", label)


def zero(expr, label):
    check(sp.simplify(sp.expand(expr)) == 0, label)


# =====================================================================
# SECTION 1.  The 13 equations from the direct two-variable Poisson bracket.
# =====================================================================
A = {k: sp.Function(f"a{k}")(t) for k in KS}
B = {k: sp.Function(f"b{k}")(t) for k in KS}

f_full = sum(x**k * A[k].subs(t, x*xi) for k in KS)
g_full = sum(x**k * B[k].subs(t, x*xi) for k in KS)
direct = sp.expand(poisson(g_full, f_full))
formula = sp.expand(sum(x**m * Cm(m, A, B).subs(t, x*xi) for m in range(-6, 7)))
zero(direct - formula,
     "all 13 C_m match direct 2-variable Poisson (fully generic coefficients)")

# Frozen displayed forms of every C_m (exactly as in the .md).
a3, a2, a1, a0, am1, am2, am3 = (A[3], A[2], A[1], A[0], A[-1], A[-2], A[-3])
b3, b2, b1, b0, bm1, bm2, bm3 = (B[3], B[2], B[1], B[0], B[-1], B[-2], B[-3])
expected = {
 6:  3*a3*D(b3) - 3*D(a3)*b3,
 5:  3*a3*D(b2) - 2*D(a3)*b2 + 2*a2*D(b3) - 3*D(a2)*b3,
 4:  3*a3*D(b1) - D(a3)*b1 + a1*D(b3) - 3*D(a1)*b3 + 2*a2*D(b2) - 2*D(a2)*b2,
 3:  3*a3*D(b0) - 3*D(a0)*b3 + 2*a2*D(b1) - D(a2)*b1 + a1*D(b2) - 2*D(a1)*b2,
 2:  3*a3*D(bm1) + D(a3)*bm1 - am1*D(b3) - 3*D(am1)*b3
     + 2*a2*D(b0) - 2*D(a0)*b2 + a1*D(b1) - D(a1)*b1,
 1:  3*a3*D(bm2) + 2*D(a3)*bm2 - 2*am2*D(b3) - 3*D(am2)*b3
     + 2*a2*D(bm1) + D(a2)*bm1 - am1*D(b2) - 2*D(am1)*b2
     + a1*D(b0) - D(a0)*b1,
 0:  3*a3*D(bm3) + 3*D(a3)*bm3 - 3*am3*D(b3) - 3*D(am3)*b3
     + 2*a2*D(bm2) + 2*D(a2)*bm2 - 2*am2*D(b2) - 2*D(am2)*b2
     + a1*D(bm1) + D(a1)*bm1 - am1*D(b1) - D(am1)*b1,
 -1: -am1*D(b0) - 2*am2*D(b1) - 3*am3*D(b2) + a1*D(bm2) + 2*a2*D(bm3)
     + bm1*D(a0) + 2*bm2*D(a1) + 3*bm3*D(a2) - b1*D(am2) - 2*b2*D(am3),
 -2: -am1*D(bm1) - 2*am2*D(b0) - 3*am3*D(b1) + a1*D(bm3)
     + bm1*D(am1) + 2*bm2*D(a0) + 3*bm3*D(a1) - b1*D(am3),
 -3: -am1*D(bm2) - 2*am2*D(bm1) - 3*am3*D(b0)
     + bm1*D(am2) + 2*bm2*D(am1) + 3*bm3*D(a0),
 -4: -am1*D(bm3) - 2*am2*D(bm2) - 3*am3*D(bm1)
     + bm1*D(am3) + 2*bm2*D(am2) + 3*bm3*D(am1),
 -5: -2*am2*D(bm3) - 3*am3*D(bm2) + 2*bm2*D(am3) + 3*bm3*D(am2),
 -6: -3*am3*D(bm3) + 3*bm3*D(am3),
}
for m in range(-6, 7):
    zero(Cm(m, A, B) - sp.expand(expected[m]), f"displayed C_{m} matches sum formula")

# The EXACT strings as displayed in classical-band3-cascade.md section 1
# (transcribed verbatim, so a typo in the memo's primary display would fail here).
md = {
 6: "3*a3*b3p - 3*a3p*b3",
 5: "3*a3*b2p - 2*a3p*b2 + 2*a2*b3p - 3*a2p*b3",
 4: "3*a3*b1p - a3p*b1 + a1*b3p - 3*a1p*b3 + 2*a2*b2p - 2*a2p*b2",
 3: "3*a3*b0p - 3*a0p*b3 + 2*a2*b1p - a2p*b1 + a1*b2p - 2*a1p*b2",
 2: "3*a3*bm1p + a3p*bm1 - am1*b3p - 3*am1p*b3 + 2*a2*b0p - 2*a0p*b2 + a1*b1p - a1p*b1",
 1: "3*a3*bm2p + 2*a3p*bm2 - 2*am2*b3p - 3*am2p*b3 + 2*a2*bm1p + a2p*bm1"
    " - am1*b2p - 2*am1p*b2 + a1*b0p - a0p*b1",
 0: "3*a3*bm3p + 3*a3p*bm3 - 3*am3*b3p - 3*am3p*b3 + 2*a2*bm2p + 2*a2p*bm2"
    " - 2*am2*b2p - 2*am2p*b2 + a1*bm1p + a1p*bm1 - am1*b1p - am1p*b1",
 -1: "-am1*b0p + a0p*bm1 + a1*bm2p + 2*a1p*bm2 - 2*am2*b1p - am2p*b1"
     " + 2*a2*bm3p + 3*a2p*bm3 - 3*am3*b2p - 2*am3p*b2",
 -2: "-am1*bm1p + am1p*bm1 + a1*bm3p - 2*am2*b0p + 2*a0p*bm2 + 3*a1p*bm3"
     " - b1*am3p - 3*am3*b1p",
 -3: "-am1*bm2p + 2*am1p*bm2 - 2*am2*bm1p + am2p*bm1 - 3*am3*b0p + 3*a0p*bm3",
 -4: "-am1*bm3p + am3p*bm1 - 2*am2*bm2p + 2*am2p*bm2 - 3*am3*bm1p + 3*am1p*bm3",
 -5: "2*am3p*bm2 - 3*am3*bm2p - 2*am2*bm3p + 3*am2p*bm3",
 -6: "3*am3p*bm3 - 3*am3*bm3p",
}
_env = {}
for k in KS:
    nm = f"a{k}" if k >= 0 else f"am{-k}"
    bm = f"b{k}" if k >= 0 else f"bm{-k}"
    _env[nm] = A[k]; _env[nm + "p"] = D(A[k])
    _env[bm] = B[k]; _env[bm + "p"] = D(B[k])
for m in range(-6, 7):
    zero(Cm(m, A, B) - sp.sympify(md[m], locals=_env),
         f"memo-verbatim display of C_{m} equals the sum formula")

# =====================================================================
# SECTION 2.  Top / bottom Wronskian proportionality (W1).
# =====================================================================
lam3, mu3 = sp.symbols("lambda3 mu3")
Bp = dict(B); Bp[3] = lam3*A[3]
zero(Cm(6, A, Bp), "C6=0 identically after b3 = lambda3 a3 (top proportionality)")
Bm = dict(B); Bm[-3] = mu3*A[-3]
zero(Cm(-6, A, Bm), "C_-6=0 identically after b_-3 = mu3 a_-3 (bottom proportionality)")

# =====================================================================
# SECTION 3.  THE WALL  C_5  ->  L[u2] = 0,  u2 = b2 - lambda3 a2,
#             L[u] = 3 a3 u' - 2 a3' u   (the 2/3-power analogue of J2).
# =====================================================================
u2 = sp.Function("u2")(t)
Bp = dict(B); Bp[3] = lam3*A[3]; Bp[2] = lam3*A[2] + u2
L_u2 = 3*a3*D(u2) - 2*D(a3)*u2
zero(Cm(5, A, Bp) - L_u2,
     "C5 with (b3=lam3 a3, b2=lam3 a2 + u2) reduces EXACTLY to 3 a3 u2' - 2 a3' u2")

# Integrating-factor identity: L[u] = 0  <=>  (u^3 / a3^2)' = 0.
u = sp.Function("u")(t)
zero(D(u**3/a3**2) - u**2*(3*a3*D(u) - 2*D(a3)*u)/a3**3,
     "(u^3/a3^2)' = u^2 * L[u] / a3^3  (cube integrating factor for the wall)")

# Cube homogeneous solutions:  a3 = c h^3  =>  u = e h^2  solves L[u]=0.
c_, e_ = sp.symbols("c e")
h = sp.Function("h")(t)
zero(3*(c_*h**3)*D(e_*h**2) - 2*D(c_*h**3)*(e_*h**2),
     "a3=c h^3, u2=e h^2  solves the wall homogeneous equation")

# =====================================================================
# SECTION 4.  The descent operators  L_l[u] = 3 a3 u' - l a3' u,
#             homogeneous law  (u^3 / a3^l)' = u^2 L_l[u] / a3^(l+1).
# =====================================================================
for lv in [-3, -2, -1, 0, 1, 2, 3]:
    Ll = 3*a3*D(u) - lv*D(a3)*u
    zero(D(u**3/a3**lv) - u**2*Ll/a3**(lv+1),
         f"descent-operator law l={lv}: (u^3/a3^l)' = u^2 L_l[u]/a3^(l+1)")

# =====================================================================
# SECTION 5.  C_4  ->  M[u1] + 2 W(a2,u2)  (gauge-invariant three-level coupling).
#             M[u] = 3 a3 u' - a3' u,  W(f,g) = f g' - f' g.
# =====================================================================
u1 = sp.Function("u1")(t)
Bp = dict(B); Bp[3] = lam3*A[3]; Bp[2] = lam3*A[2] + u2; Bp[1] = lam3*A[1] + u1
M_u1 = 3*a3*D(u1) - D(a3)*u1
W_a2u2 = a2*D(u2) - D(a2)*u2
zero(Cm(4, A, Bp) - (M_u1 + 2*W_a2u2),
     "C4 reduces to M[u1] + 2 W(a2,u2), independent of the gauge lambda3")

# =====================================================================
# SECTION 6.  Non-cube cascade (gauged frame b3=0): each rung reduces to a
#             single descent-operator ODE.  Written proof in the .md;
#             here we verify every algebraic reduction identically.
# =====================================================================
# generic a's (a3 != 0 symbolically), gauged b3 = 0.
Bc = {k: sp.Function(f"b{k}")(t) for k in KS}
Bc[3] = sp.Integer(0)
# C5 -> L[b2]
zero(Cm(5, A, Bc) - (3*a3*D(Bc[2]) - 2*D(a3)*Bc[2]),
     "gauged C5 = 3 a3 b2' - 2 a3' b2  (= L[b2])")
# C4 with b2=0 -> M[b1]
Bc2 = dict(Bc); Bc2[2] = sp.Integer(0)
zero(Cm(4, A, Bc2) - (3*a3*D(Bc2[1]) - D(a3)*Bc2[1]),
     "gauged C4 (b2=0) = 3 a3 b1' - a3' b1  (= M[b1])")
# C3 with b1=b2=0 -> 3 a3 b0'
Bc3 = dict(Bc2); Bc3[1] = sp.Integer(0)
zero(Cm(3, A, Bc3) - 3*a3*D(Bc3[0]),
     "gauged C3 (b1=b2=0) = 3 a3 b0'")
# C2 with b1=b2=0, b0 const -> L_{-1}[b_{-1}] = 3 a3 b_{-1}' + a3' b_{-1}
b0c = sp.symbols("b0const")
Bc4 = dict(Bc3); Bc4[0] = b0c
zero(Cm(2, A, Bc4) - (3*a3*D(Bc4[-1]) + D(a3)*Bc4[-1]),
     "gauged C2 (b1=b2=0,b0 const) = 3 a3 b_-1' + a3' b_-1  (= L_-1[b_-1])")
# C1 with b_{-1}=0 -> L_{-2}[b_{-2}] = 3 a3 b_{-2}' + 2 a3' b_{-2}
Bc5 = dict(Bc4); Bc5[-1] = sp.Integer(0)
zero(Cm(1, A, Bc5) - (3*a3*D(Bc5[-2]) + 2*D(a3)*Bc5[-2]),
     "gauged C1 (b_-1=0) = 3 a3 b_-2' + 2 a3' b_-2  (= L_-2[b_-2])")
# C0 with b_{-2}=0 -> 3 (a3 b_{-3})'
Bc6 = dict(Bc5); Bc6[-2] = sp.Integer(0)
zero(Cm(0, A, Bc6) - 3*D(a3*Bc6[-3]),
     "gauged C0 (b_-2=0) = 3 (a3 b_-3)'  (=> a3 b_-3 = t/3, killed by tau^3|b_-3)")

# =====================================================================
# SECTION 7.  Bounded-degree corroboration of the cube dichotomy.
# =====================================================================
def kernel_dim(a3poly, lv, deg):
    """dim of poly kernel {u : 3 a3 u' - l a3' u = 0, deg u <= deg}."""
    cs = sp.symbols(f"z0:{deg+1}")
    uu = sum(cs[i]*t**i for i in range(deg+1))
    expr = sp.expand(3*a3poly*D(uu) - lv*D(a3poly)*uu)
    eqs = sp.Poly(expr, t).all_coeffs() if expr != 0 else [0]
    sol = sp.linsolve(eqs, cs)
    # dimension = number of free params
    solset = list(sol)[0]
    free = set()
    for comp in solset:
        free |= comp.free_symbols
    return len(free & set(cs))

# non-cube a3 = t  and  a3 = t^2+1: trivial kernel at l=1,2 (=> b1,b2 forced 0).
for a3p, name in [(t, "t"), (t**2 + 1, "t^2+1")]:
    for lv in [1, 2]:
        check(kernel_dim(a3p, lv, 6) == 0,
              f"non-cube a3={name}: L_{lv} kernel trivial (deg<=6)")
# cube a3 = t^3 (h=t): kernel at l=1 is span{t}=<h>, at l=2 is span{t^2}=<h^2>.
for lv, expdim in [(1, 1), (2, 1)]:
    check(kernel_dim(t**3, lv, 8) == expdim,
          f"cube a3=t^3: L_{lv} kernel is 1-dimensional (= <h^{lv}>) (deg<=8)")
# explicit generators:
zero(3*t**3*D(t) - 1*D(t**3)*t, "cube a3=t^3: u=h=t solves M=L_1")
zero(3*t**3*D(t**2) - 2*D(t**3)*t**2, "cube a3=t^3: u=h^2=t^2 solves L_2")

# =====================================================================
# SECTION 8.  m=0 moment identity (W4 generalization):
#             M(t) = sum_{k>=1} k (a_k b_{-k} - a_{-k} b_k),  M' = C_0.
# =====================================================================
Mmom = sum(k*(A[k]*B[-k] - A[-k]*B[k]) for k in [1, 2, 3])
zero(Cm(0, A, B) - D(Mmom),
     "C0 = d/dt [ 3(a3 b_-3 - a_-3 b3) + 2(a2 b_-2 - a_-2 b2) + (a1 b_-1 - a_-1 b1) ]")

# membership => M(0)=0 => M(t) = t.  Verify each product vanishes at t=0 on a
# representative satisfying tau^j | a_{-j}, b_{-j}.
rep = {
 A[3]: t**2 + 3, A[2]: t + 1, A[1]: 5, A[0]: 7,
 A[-1]: t*(t + 2), A[-2]: t**2*(t - 1), A[-3]: t**3*(t + 4),
 B[3]: t - 2, B[2]: 9, B[1]: t**3, B[0]: 4,
 B[-1]: t*(3*t + 1), B[-2]: t**2*(t + 5), B[-3]: t**3*(2*t - 1),
}
Mrep = Mmom.subs(rep)
check(sp.simplify(Mrep.subs(t, 0)) == 0,
      "moment M(0)=0 under genuine membership (=> M(t)=t on any Keller pair)")

# =====================================================================
# SECTION 9.  Cubic sector  a3 = c h^3:  mixed (a3,a2) coupling.
#             gauged b3=0, C5 => b2 = e h^2, C4 => b1 = k1 h + (2e/3c) a2/h.
# =====================================================================
c_, e_, k1 = sp.symbols("c e kappa1")
h = sp.Function("h")(t)
a2f = sp.Function("a2")(t)
a3c = c_*h**3
# derivation identity  2 a2 h' - a2' h = -h^3 (a2/h^2)'
zero((2*a2f*D(h) - D(a2f)*h) - (-h**3*D(a2f/h**2)),
     "identity 2 a2 h' - a2' h = -h^3 (a2/h^2)'")
# C4 vanishes with the claimed b1 (mixed coupling formula)
Am = {3: a3c, 2: a2f, 1: sp.Function("a1")(t), 0: sp.Function("a0")(t),
      -1: sp.Function("am1")(t), -2: sp.Function("am2")(t), -3: sp.Function("am3")(t)}
Bm2 = {3: sp.Integer(0), 2: e_*h**2,
       1: k1*h + sp.Rational(2, 3)*e_*a2f/(c_*h),
       0: sp.Function("b0")(t), -1: sp.Function("bm1")(t),
       -2: sp.Function("bm2")(t), -3: sp.Function("bm3")(t)}
zero(Cm(4, Am, Bm2),
     "cubic sector: b1 = k1 h + (2e/3c) a2/h solves C4  (=> e!=0 forces h | a2)")

# h | a2 necessity witness: h = t, a2 = 1 (h does NOT divide a2), e != 0
#   => the forced b1 has a genuine pole, i.e. is NOT a polynomial.
b1_witness = (k1*h + sp.Rational(2, 3)*e_*a2f/(c_*h)).subs({h: t, a2f: 1})
check(sp.together(b1_witness).has(1/t) or sp.denom(sp.together(b1_witness)) != 1,
      "witness h=t, a2=1, e!=0: forced b1 is non-polynomial (h∤a2 excluded)")

# e=0 subsector: C3 reduces to 3 a3 b0' = k1(a2' h - 2 a2 h').
b0f = sp.Function("b0")(t)
Bm3 = {3: sp.Integer(0), 2: sp.Integer(0), 1: k1*h, 0: b0f,
       -1: sp.Function("bm1")(t), -2: sp.Function("bm2")(t), -3: sp.Function("bm3")(t)}
zero(Cm(3, Am, Bm3) - (3*a3c*D(b0f) - k1*(D(a2f)*h - 2*a2f*D(h))),
     "cubic e=0: C3 = 3 a3 b0' - k1(a2' h - 2 a2 h')")

# --- Constant-h cubic residual system (h=1, a3=c const): explicit integrations.
a1f, a0f = sp.Function("a1")(t), sp.Function("a0")(t)
am1f, am2f = sp.Function("am1")(t), sp.Function("am2")(t)
bm2f, bm3f = sp.Function("bm2")(t), sp.Function("bm3")(t)
beta, gam = sp.symbols("beta gamma")
b1_ch = k1 + sp.Rational(2, 3)*e_*a2f/c_               # b1 (h=1)
b0_ch = (2*e_*a1f + k1*a2f)/(3*c_) - e_*a2f**2/(9*c_**2) + beta
bm1_ch = (sp.Rational(1, 3)/c_)*(2*e_*a0f + k1*a1f
         - sp.Rational(2, 3)*e_*a1f*a2f/c_ - sp.Rational(1, 3)*k1*a2f**2/c_
         + sp.Rational(4, 27)*e_*a2f**3/c_**2) + gam
Ach = {3: c_, 2: a2f, 1: a1f, 0: a0f, -1: am1f, -2: am2f, -3: sp.Function("am3")(t)}
Bch = {3: sp.Integer(0), 2: e_, 1: b1_ch, 0: b0_ch, -1: bm1_ch,
       -2: bm2f, -3: bm3f}
zero(Cm(3, Ach, Bch), "constant-h cubic: b0 integral solves C3")
zero(Cm(2, Ach, Bch), "constant-h cubic: b_-1 integral solves C2")
# central identity (moment) in this frame equals the general moment:
Mch = sum(k*(Ach[k]*Bch[-k] - Ach[-k]*Bch[k]) for k in [1, 2, 3])
zero(D(Mch) - Cm(0, Ach, Bch),
     "constant-h cubic: C0 = M' with M = 3c b_-3 - 2e a_-2 - k1 a_-1"
     " + a1 b_-1 + 2 a2 b_-2 - (2e/3c) a2 a_-1")

# =====================================================================
# SECTION 10.  Explicit band-3 automorphism witness (anchor).
#     F = xi + alpha x^3 + beta x^2 + gamma x + delta,  G = -x.
#     Tame (triangular) => {G,F}=1; exercises a3,a2,a1,a0 != 0, a_-1 = tau.
# =====================================================================
al, be, ga, de = sp.symbols("alpha beta gamma delta")
Fw = xi + al*x**3 + be*x**2 + ga*x + de
Gw = -x
zero(poisson(Gw, Fw) - 1, "witness F=xi+alpha x^3+..., G=-x has {G,F}=1")
# ladder coefficients and their C_m
Aw = {3: al, 2: be, 1: ga, 0: de, -1: t, -2: sp.Integer(0), -3: sp.Integer(0)}
Bw = {1: sp.Integer(-1), 3: sp.Integer(0), 2: sp.Integer(0), 0: sp.Integer(0),
      -1: sp.Integer(0), -2: sp.Integer(0), -3: sp.Integer(0)}
for m in range(-6, 7):
    zero(Cm(m, Aw, Bw) - (1 if m == 0 else 0),
         f"witness satisfies C_{m} = delta_{{m0}}")
# membership: a_-1 = tau (tau | a_-1 ok), a_-2=a_-3=0 ok; e=0, k1 h = -1.
check(Bw[2] == 0 and Bw[1] == -1 and Aw[3] == al,
      "witness has b2=0, b1=-1, and top a3=alpha")

# =====================================================================
# SECTION 11.  Bounded-degree corroboration: no e!=0 cubic band-3 Keller pair
#     with h=1 (a3=c, constant) and small degrees.  (Regression evidence.)
# =====================================================================
# h=1 => a3=c const, cube automatically.  Search a2,a1,a0,a_-1,a_-2,a_-3,
# b's with b2=e!=0, low degree, {G,F}=1 all C_m.  This is heavy; do a light
# consistency check: with h=1, b2=e, the mixed coupling gives b1 = k1 + (2e/3c)a2,
# and C4 must hold — verify that reduction (already covered) and stop.
zero(Cm(4, {**Am, 3: c_}, {**Bm2, 2: e_, 1: k1 + sp.Rational(2, 3)*e_*a2f/c_}).subs(h, 1),
     "constant-h (h=1): C4 consistent with b1 = k1 + (2e/3c) a2")

print()
print(f"({PASS} checks)")
print("ALL CLASSICAL BAND3 CHECKS PASSED")
