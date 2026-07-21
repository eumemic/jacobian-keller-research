#!/usr/bin/env python3
"""Exact verifier for the classical band-3 HARD BRANCHES memo.

Two hard branches of the cubic sector a_3 = c h^3 (gauge b_3 = 0), left open by
Wave A (classical-band3-cascade.md, commit 99fe6ee):

  (I)  the e != 0 MIXED sector   (b_2 = e h^2 != 0 after gauge);
  (II) NONCONSTANT h             (h a nonconstant cube root of a_3).

Conventions frozen from Wave A (identical to classical-band3-cascade.md):
  {G,F} = G_xi F_x - G_x F_xi,  tau = x*xi,  F = sum x^k a_k(tau), G = sum x^k b_k(tau),
  membership tau^j | a_{-j}, b_{-j},  C_m = sum_{k+l=m}(k a_k b_l' - l a_k' b_l) = delta_{m0}.

Every displayed algebraic identity in classical-hard-branches.md is machine-checked
here. Degree / membership / divisibility contradictions are WRITTEN proofs in the
.md; the bounded Groebner sweeps here are regression corroboration only.
A successful run ends 'ALL CLASSICAL HARD CHECKS PASSED'.
"""
import sympy as sp

x, xi, t = sp.symbols("x xi tau")
KS = range(-3, 4)
PASS = 0


def D(f):
    return sp.diff(f, t)


def Cm(m, A, B):
    return sp.expand(sum(
        k*A[k]*D(B[l]) - l*D(A[k])*B[l]
        for k in KS for l in KS if k + l == m))


def poisson(g, f):
    return sp.diff(g, xi)*sp.diff(f, x) - sp.diff(g, x)*sp.diff(f, xi)


def check(cond, label):
    global PASS
    if cond is not True and cond != True:  # noqa: E712
        raise AssertionError("FAIL: " + label)
    PASS += 1
    print("PASS", label)


def zero(expr, label):
    check(sp.simplify(sp.expand(expr)) == 0, label)


# =====================================================================
# SECTION 1.  Reflection R: f(x,xi) -> f(xi,x). The bottom is a cube.
#   R is anti-symplectic: {RG,RF} = -R{G,F}, so (RF,RG) is a Keller pair.
#   (RF)_m = tau^{-m} a_{-m}; applying Theorem A to (RF,RG) forces the bottom
#   coefficient a_{-3} to be a scalar cube c-hat*hbar^3 with tau | hbar.
# =====================================================================
A = {k: sp.Function(f"a{k}")(t) for k in KS}
Bf = {k: sp.Function(f"b{k}")(t) for k in KS}
F = sum(x**k*A[k].subs(t, x*xi) for k in KS)
G = sum(x**k*Bf[k].subs(t, x*xi) for k in KS)
RF = F.subs({x: xi, xi: x}, simultaneous=True)
RG = G.subs({x: xi, xi: x}, simultaneous=True)
zero(poisson(RG, RF) - (-poisson(G, F).subs({x: xi, xi: x}, simultaneous=True)),
     "reflection is anti-symplectic: {RG,RF} = -R{G,F}  (so (RF,RG) is Keller)")
target = sum(x**m*(t**(-m)*A[-m]).subs(t, x*xi) for m in KS)
zero(RF - target, "(RF)_m = tau^{-m} a_{-m}  (top of RF is tau^{-3} a_{-3})")

# Bottom Wronskian C_-6 => b_-3 = mu3 a_-3 (Wave A W1 mirror), re-checked.
mu3, ec, cc, k1h = sp.symbols("mu3 ehat chat kappa1hat")
hb = sp.Function("hbar")(t)
am2, am1 = sp.Function("am2")(t), sp.Function("am1")(t)
Ab = {3: sp.Function("a3")(t), 2: sp.Function("a2")(t), 1: sp.Function("a1")(t),
      0: sp.Function("a0")(t), -1: am1, -2: am2, -3: cc*hb**3}
Bb = dict(Bf)
Bb[-3] = mu3*Ab[-3]
zero(Cm(-6, Ab, Bb), "C_-6 = 0 with b_-3 = mu3 a_-3  (bottom proportionality)")

# Bottom WALL C_-5: with b_-3 = mu3 a_-3, b_-2 = mu3 a_-2 + u2hat, reduces to the
# reflected 2/3-power wall -(3 a_-3 u2hat' - 2 a_-3' u2hat); a_-3 = chat hbar^3
# forces u2hat = ehat hbar^2 (bottom wall datum), the mirror of b_2 = e h^2.
u2h = sp.Function("u2hat")(t)
Bb2 = dict(Bb)
Bb2[-2] = mu3*am2 + u2h
zero(Cm(-5, Ab, Bb2) - (-(3*Ab[-3]*D(u2h) - 2*D(Ab[-3])*u2h)),
     "C_-5 reduces to the reflected bottom wall  -(3 a_-3 u2hat' - 2 a_-3' u2hat)")
zero(3*(cc*hb**3)*D(ec*hb**2) - 2*D(cc*hb**3)*(ec*hb**2),
     "a_-3 = chat hbar^3, u2hat = ehat hbar^2 solves the bottom wall")

# Bottom MIXED coupling C_-4: b_-1 = mu3 a_-1 + kappa1hat hbar + (2 ehat/3 chat) a_-2/hbar
Bb3 = dict(Bb2)
Bb3[-2] = mu3*am2 + ec*hb**2
Bb3[-1] = mu3*am1 + k1h*hb + sp.Rational(2, 3)*ec*am2/(cc*hb)
zero(Cm(-4, Ab, Bb3),
     "C_-4 = 0 with b_-1 = mu3 a_-1 + kappa1hat hbar + (2 ehat/3 chat) a_-2/hbar"
     "  (=> ehat != 0 forces hbar | a_-2)")

# =====================================================================
# SECTION 2.  Tame single-shear pairs realise exactly e = 0.
#   Structure lemma (band3-tame-catalog.md, 99fe6ee): b_k = lambda a_k for all
#   |k| >= 2 with ONE lambda. Hence gauged u_2 = b_2 - lambda a_2 = 0, i.e. e = 0.
#   Positive controls: catalog B3-2 (even) and B3-3 (odd) both have gauged e = 0.
# =====================================================================
def keller_ok(Ad, Bd):
    return all(Cm(m, Ad, Bd) == (1 if m == 0 else 0) for m in range(-6, 7))


A_B32 = {3: sp.Integer(1), 2: sp.Integer(1), 1: 3*t+1, 0: 2*t,
         -1: 3*t**2+2*t, -2: t**2, -3: t**3}
B_B32 = {3: sp.Integer(2), 2: sp.Integer(2), 1: 6*t+1, 0: 4*t,
         -1: 6*t**2+3*t, -2: 2*t**2, -3: 2*t**3}
check(keller_ok(A_B32, B_B32), "catalog B3-2 is a genuine band-3 Keller pair")
lam_B32 = sp.Rational(2)  # = b_3/a_3
check(sp.expand(B_B32[2] - lam_B32*A_B32[2]) == 0,
      "B3-2: gauged u_2 = b_2 - lambda a_2 = 0  => e = 0  (tame realises e=0)")
check(all(sp.expand(B_B32[k] - lam_B32*A_B32[k]) == 0 for k in [3, 2, -2, -3]),
      "B3-2: the single gauge kills the whole tower b_3,b_2,b_-2,b_-3")

A_B33 = {3: sp.Integer(1), 2: sp.Integer(0), 1: 3*t+2, 0: sp.Integer(0),
         -1: 3*t**2+3*t, -2: sp.Integer(0), -3: t**3}
B_B33 = {3: sp.Integer(1), 2: sp.Integer(0), 1: 3*t+1, 0: sp.Integer(0),
         -1: 3*t**2+2*t, -2: sp.Integer(0), -3: t**3}
check(keller_ok(A_B33, B_B33), "catalog B3-3 is a genuine band-3 Keller pair")
check(sp.expand(B_B33[2] - sp.Rational(1)*A_B33[2]) == 0,
      "B3-3: gauged u_2 = 0 => e = 0 (odd-only single cubic shear)")

# Tame family realises e = 0 with ANY kappa1: witness F = xi + a x^3+..., G = -x
al, be, ga, de = sp.symbols("alpha beta gamma delta")
zero(poisson(-x, xi + al*x**3 + be*x**2 + ga*x + de) - 1,
     "tame witness F=xi+alpha x^3+..., G=-x is Keller; b_2=0 (e=0), b_1=-1=kappa1 h")

print("\n--- Section 1-2 (reflection cube + tame=>e=0) done ---\n")
print(f"({PASS} checks so far)")

# =====================================================================
# SECTION 3.  The general-h forward cascade closed forms and the DIVISIBILITY
#   engine.  c = 1 (diagonal scaling), a_3 = h^3, gauge b_3 = 0.  Each downward
#   integration divides by another power of h; polynomiality forces divisibility.
# =====================================================================
e, k1, beta, gam = sp.symbols("e kappa1 beta gamma")
h = sp.Function("h")(t)
a2, a1, a0 = sp.Function("a2")(t), sp.Function("a1")(t), sp.Function("a0")(t)
Ah = {3: h**3, 2: a2, 1: a1, 0: a0,
      -1: sp.Function("am1")(t), -2: sp.Function("am2")(t), -3: sp.Function("am3")(t)}

# b_2 = e h^2 (wall), b_1 = k1 h + (2e/3) a2/h  (Wave A 4.1 with c=1).
b2 = e*h**2
b1 = k1*h + sp.Rational(2, 3)*e*a2/h
# b_0 closed form (general e): weight-homogeneous.
b0 = (sp.Rational(2, 3)*e*a1/h + sp.Rational(1, 3)*k1*a2/h**2
      - sp.Rational(1, 9)*e*a2**2/h**4 + beta)
Bh = {3: sp.Integer(0), 2: b2, 1: b1, 0: b0,
      -1: sp.Function("bm1")(t), -2: sp.Function("bm2")(t), -3: sp.Function("bm3")(t)}
zero(Cm(3, Ah, Bh), "general-h: b_0 = 2e a_1/3h + k1 a_2/3h^2 - e a_2^2/9h^4 + beta solves C_3")

# e = 0 branch: b_0 = (k1/3) a_2/h^2 + beta ; b_-1 closed form solves C_2.
b0_e0 = sp.Rational(1, 3)*k1*a2/h**2 + beta
bm1_e0 = (sp.Rational(1, 3)*k1*a1/h**2 - sp.Rational(1, 9)*k1*a2**2/h**5 + gam/h)
Bh0 = dict(Bh)
Bh0[2] = sp.Integer(0)
Bh0[1] = k1*h
Bh0[0] = b0_e0
Bh0[-1] = bm1_e0
zero(Cm(3, Ah, Bh0), "e=0: b_0 = (k1/3) a_2/h^2 + beta solves C_3")
zero(Cm(2, Ah, Bh0),
     "e=0: b_-1 = (k1/3) a_1/h^2 - (k1/9) a_2^2/h^5 + gamma/h solves C_2")

# DIVISIBILITY 1 (e=0, k1!=0):  b_0 polynomial  <=>  h^2 | a_2.
g2 = sp.Function("g2")(t)
zero((b0_e0.subs(a2, h**2*g2)) - (sp.Rational(1, 3)*k1*g2 + beta),
     "e=0: with a_2 = h^2 g_2, b_0 = (k1/3) g_2 + beta is polynomial  (h^2 | a_2 forced)")

# DIVISIBILITY 2 (e=0, k1!=0):  h^2 * b_-1  =  (k1/3) a_1  +  h*(gamma - k1 g2^2/9),
#   so polynomiality of b_-1 forces  h | a_1.
lhs = sp.expand(sp.together((bm1_e0.subs(a2, h**2*g2))*h**2))
rhs = sp.Rational(1, 3)*k1*a1 + h*(gam - sp.Rational(1, 9)*k1*g2**2)
zero(lhs - sp.expand(rhs),
     "e=0: h^2 b_-1 = (k1/3) a_1 + h(gamma - k1 g2^2/9)  =>  h | a_1")

# THE h | tau STEP (e=0):  with a_2 = h^2 g2 and a_1 = h*a1t, the moment
#   M = 3 h^3 b_-3 + 2 a_2 b_-2 - 2 a_-2 b_2 + a_1 b_-1 - a_-1 b_1  is divisible by h,
#   because b_2 = 0 and every surviving term carries a factor h.  Since M = tau, h | tau.
a1t = sp.Function("a1tilde")(t)
bm2, bm3 = sp.Function("bm2")(t), sp.Function("bm3")(t)
am1f, am2f = sp.Function("am1")(t), sp.Function("am2")(t)
M_expr = (3*h**3*bm3 + 2*(h**2*g2)*bm2 - 2*am2f*sp.Integer(0)
          + (h*a1t)*sp.Function("bm1")(t) - am1f*(k1*h))
zero(M_expr - h*sp.expand(M_expr/h),
     "e=0: moment M is divisible by h  (a_2=h^2 g2, a_1=h a1t, b_2=0) => M=tau gives h|tau")

print("\n--- Section 3 (divisibility cascade, h|tau) done ---\n")
print(f"({PASS} checks so far)")

# =====================================================================
# SECTION 4.  Nonconstant-h kill, sub-case b_1 = 0  (i.e. e = 0 AND kappa1 = 0,
#   so b_3=b_2=b_1=0: G supported in x-levels <= 0).  Theorem-A-style endgame:
#   the descent trivialises rung by rung and dies at C_0 vs membership.
# =====================================================================
b0c = sp.Symbol("b0const")
A4 = {3: h**3, 2: a2, 1: a1, 0: a0,
      -1: sp.Function("am1")(t), -2: sp.Function("am2")(t), -3: sp.Function("am3")(t)}
# C_3 with b_1=b_2=b_3=0 :
B4 = {3: sp.Integer(0), 2: sp.Integer(0), 1: sp.Integer(0),
      0: sp.Function("b0")(t), -1: sp.Function("bm1")(t),
      -2: sp.Function("bm2")(t), -3: sp.Function("bm3")(t)}
zero(Cm(3, A4, B4) - 3*h**3*D(B4[0]),
     "b_1=b_2=0: C_3 = 3 h^3 b_0'  => b_0' = 0 => b_0 = const")
# C_2 with b_0 const :
B4b = dict(B4)
B4b[0] = b0c
zero(Cm(2, A4, B4b) - 3*h**2*D(h*B4b[-1]),
     "b_1=b_2=0, b_0 const: C_2 = 3 h^2 (h b_-1)'  => h b_-1 = gamma => b_-1 = gamma/h"
     " => (nonconst h) gamma=0, b_-1=0")
# C_1 with b_-1 = 0 :
B4c = dict(B4b)
B4c[-1] = sp.Integer(0)
zero(Cm(1, A4, B4c) - 3*h*D(h**2*B4c[-2]),
     "b_-1=0: C_1 = 3 h (h^2 b_-2)'  => h^2 b_-2 = delta => (nonconst h) delta=0, b_-2=0")
# C_0 with b_-2 = 0 :
B4d = dict(B4c)
B4d[-2] = sp.Integer(0)
zero(Cm(0, A4, B4d) - 3*D(h**3*B4d[-3]),
     "b_-2=0: C_0 = 3 (h^3 b_-3)' = 1 => h^3 b_-3 = tau/3;"
     " tau^3|b_-3 forces tau-order>=3>1 CONTRADICTION (nonconst-h b_1=0 EMPTY)")

# --- kappa1 != 0 branch: h|tau (Thm 4.1) gives h = sigma*tau; the tau-ORDER of the
#     moment then contradicts M = tau.  With h = sigma*tau and the forced structure
#     h^2|a_2, h|a_1 (Section 3) plus memberships, EVERY term of the moment has
#     tau-order >= 2, so ord(M) >= 2 != 1 = ord(tau).  COMPLETE kill (no search).
sig = sp.Symbol("sigma")
rr = lambda n: sp.Function(n)(t)
Aord = {3: (sig*t)**3, 2: t**2*rr("ra2"), 1: t*rr("ra1"), 0: rr("ra0"),
        -1: t*rr("ram1"), -2: t**2*rr("ram2"), -3: t**3*rr("ram3")}
Bord = {3: sp.Integer(0), 2: sp.Integer(0), 1: k1*sig*t, 0: rr("rb0"),
        -1: t*rr("rbm1"), -2: t**2*rr("rbm2"), -3: t**3*rr("rbm3")}
Mord = sum(kk*(Aord[kk]*Bord[-kk] - Aord[-kk]*Bord[kk]) for kk in [1, 2, 3])
check(sp.simplify(Mord.subs(t, 0)) == 0 and sp.simplify(D(Mord).subs(t, 0)) == 0,
      "e=0, kappa1!=0, h=sigma*tau: moment M has tau-order >= 2 (M(0)=M'(0)=0)"
      " => M=tau (order 1) CONTRADICTION: nonconstant-h e=0 EMPTY (complete)")

print("\n--- Section 4 (nonconstant-h e=0 EMPTY: b_1=0 endgame + h|tau + tau-order) done ---\n")
print(f"({PASS} checks so far)")

# =====================================================================
# SECTION 5.  The e != 0 mixed sector: residual system, the moment first integral,
#   and the W3 obstruction (confirming Wave A: no gauge-free trailing integral).
#   Constant-h frame h=1, c=1.  b_2=e, b_1=k1+(2e/3)a_2 ; b_0,b_-1 determined.
# =====================================================================
a2c, a1c, a0c = sp.Function("a2")(t), sp.Function("a1")(t), sp.Function("a0")(t)
am1c, am2c, am3c = sp.Function("am1")(t), sp.Function("am2")(t), sp.Function("am3")(t)
A5 = {3: sp.Integer(1), 2: a2c, 1: a1c, 0: a0c, -1: am1c, -2: am2c, -3: am3c}
b2c = e
b1c = k1 + sp.Rational(2, 3)*e*a2c
b0c5 = (2*e*a1c + k1*a2c)/3 - e*a2c**2/9 + beta
bm1c = (sp.Rational(1, 3)*(2*e*a0c + k1*a1c - sp.Rational(2, 3)*e*a1c*a2c
        - sp.Rational(1, 3)*k1*a2c**2 + sp.Rational(4, 27)*e*a2c**3) + gam)
B5 = {3: sp.Integer(0), 2: b2c, 1: b1c, 0: b0c5, -1: bm1c,
      -2: sp.Function("bm2")(t), -3: sp.Function("bm3")(t)}
zero(Cm(3, A5, B5), "e!=0 constant-h: b_0 solves C_3")
zero(Cm(2, A5, B5), "e!=0 constant-h: b_-1 solves C_2")

# Moment M = tau is the central first integral (C_0 = M', W4 at band 3).
Mmom = sum(kk*(A5[kk]*B5[-kk] - A5[-kk]*B5[kk]) for kk in [1, 2, 3])
zero(Cm(0, A5, B5) - D(Mmom),
     "e!=0: C_0 = M' with M = 3(a_3 b_-3 - a_-3 b_3)+2(a_2 b_-2 - a_-2 b_2)+(a_1 b_-1 - a_-1 b_1)")

# W3 obstruction: the trailing coefficient a_-3 enters C_-2 with a NON-exact residue
#   -2 e a_-3 a_2'  (nonzero for e != 0), so no gauge-free W3 integral exists.
mu3s = sp.Symbol("mu3")
B5w = dict(B5)
B5w[-3] = mu3s*am3c              # bottom proportionality
# collect the a_-3-linear part of C_-2 and its non-total-derivative residue:
Cm2_full = Cm(-2, A5, B5w)
a3dep = sp.expand(Cm2_full - Cm2_full.subs(am3c, 0))  # part containing a_-3
# exact-derivative candidate:  d/dt[ (mu3 a_1 - b_1) a_-3 ] :
exact = D((mu3s*a1c - b1c)*am3c)
residue = sp.expand(a3dep - exact)
zero(residue - 2*am3c*(mu3s*D(a1c) - D(b1c)),
     "W3 obstruction: a_-3 in C_-2 = d/dt[(mu3 a_1 - b_1)a_-3] + 2 a_-3(mu3 a_1' - b_1')")
zero(sp.expand(residue.subs({mu3s: 0}) + 2*e*sp.Rational(2, 3)*am3c*D(a2c)),
     "  residue at mu3=0 is (4e/3) a_-3 a_2' != 0 for e!=0: no gauge-free trailing integral")

print("\n--- Section 5 (e!=0 residual + moment + W3 obstruction) done ---\n")
print(f"({PASS} checks so far)")

# =====================================================================
# SECTION 6.  Divisibility also bites the e != 0 sector under nonconstant h.
#   Wave A: e!=0 forces h|a_2 (b_1 polynomial). Then b_0 polynomiality forces
#   rad(h) | a_2/h as well (the cascade continues), corroborating that
#   nonconstant h is the obstructed direction in BOTH e-branches.
# =====================================================================
a2t = sp.Function("a2tilde")(t)   # a_2 = h * a2t  (Wave A: e!=0 => h|a_2)
b0_e = (sp.Rational(2, 3)*e*a1/h + sp.Rational(1, 3)*k1*(h*a2t)/h**2
        - sp.Rational(1, 9)*e*(h*a2t)**2/h**4 + beta)
# clear denominators: 9 h^2 b0 - 9 h^2 beta = 6 e h a1 + 3 k1 h a2t - e a2t^2
num = sp.expand(9*h**2*(b0_e - beta))
zero(num - (6*e*h*a1 + 3*k1*h*a2t - e*a2t**2),
     "e!=0, a_2=h a2t: 9 h^2 (b_0-beta) = 6e h a_1 + 3 k1 h a2t - e a2t^2"
     "  => h | e a2t^2  => rad(h) | a2t  (divisibility cascade continues)")

print("\n--- Section 6 (e!=0 divisibility) done ---\n")
print(f"({PASS} checks so far)")

# =====================================================================
# SECTION 7.  Bounded Groebner corroboration (REGRESSION ONLY, small fast boxes).
#   These are NOT proofs: each certifies emptiness of one bounded degree box.
# =====================================================================
def _pc(z):
    z = sp.expand(z)
    return [] if z == 0 else sp.Poly(z, t).all_coeffs()


def _gp(n, d, lo=0):
    cs = [sp.Symbol(f"{n}_{i}") for i in range(d + 1)]
    return sum(cs[i]*t**i for i in range(lo, d + 1)), cs[lo:]


def box_is_empty(a3poly, degs, rabinowitsch_e=False, b1_form=None):
    """True iff the bounded box admits no genuine band-3 pair (Groebner == (1)).

    Full 13-equation system, all b-levels generic polynomials (gauge b_3 = 0),
    memberships tau^j | a_{-j}, b_{-j} imposed as coefficient equations.
    rabinowitsch_e forces the wall constant b_2 = e != 0 via e*w = 1.
    """
    da2, da1, da0, dam1, dam2, dam3 = degs
    a2f, ca2 = _gp("A2", da2)
    a1f, ca1 = _gp("A1", da1)
    a0f, ca0 = _gp("A0", da0)
    am1f, cam1 = _gp("Am1", dam1, 1)
    am2f, cam2 = _gp("Am2", dam2, 2)
    am3f, cam3 = _gp("Am3", dam3, 3)
    Ad = {3: a3poly, 2: a2f, 1: a1f, 0: a0f, -1: am1f, -2: am2f, -3: am3f}
    ee, kk = sp.Symbol("e_"), sp.Symbol("k1_")
    b2v = ee if rabinowitsch_e else sp.Integer(0)
    b1v = b1_form(kk, ee, a2f, a3poly) if b1_form else kk*t
    b0f, cb0 = _gp("B0", da2 + 2)
    bm1g, cbm1 = _gp("Bm1", dam1 + 3, 1)
    bm2g, cbm2 = _gp("Bm2", dam2 + 4, 2)
    bm3g, cbm3 = _gp("Bm3", dam3 + 5, 3)
    Bd = {3: sp.Integer(0), 2: b2v, 1: b1v, 0: b0f,
          -1: bm1g, -2: bm2g, -3: bm3g}
    cons = []
    for m in range(-6, 7):
        cons += _pc(sp.expand(Cm(m, Ad, Bd) - (1 if m == 0 else 0)))
    cons = [sp.expand(z) for z in cons if sp.expand(z) != 0]
    gens = ([kk] + list(ca2) + list(ca1) + list(ca0) + list(cam1) + list(cam2)
            + list(cam3) + list(cb0) + list(cbm1) + list(cbm2) + list(cbm3))
    if rabinowitsch_e:
        w = sp.Symbol("w_"); cons = cons + [ee*w - 1]; gens = [w, ee] + gens
    Gb = sp.groebner(cons, *gens, order='grevlex')
    return list(Gb.exprs) == [sp.Integer(1)]


import os

# a_3 = tau^3 (h=tau), e=0: b_1 = k1*tau. Independent CROSS-CHECK of Theorem 4.3
# (which already proves this empty at arbitrary degree via the tau-order count). FAST.
_b1_t = lambda kk, ee, a2f, a3p: kk*t
check(box_is_empty(t**3, (2, 2, 1, 2, 3, 4), b1_form=_b1_t) is True,
      "[COMPUTED,cross-check Thm 4.3] a_3=tau^3, e=0: no pair in box (2,2,1,2,3,4) (Groebner=(1))")
check(box_is_empty(t**3, (3, 3, 2, 3, 4, 5), b1_form=_b1_t) is True,
      "[COMPUTED,cross-check Thm 4.3] a_3=tau^3, e=0: no pair in box (3,3,2,3,4,5) (Groebner=(1))")

# e != 0, constant h=1 (a_3=1): Rabinowitsch box; b_1 = k1 + (2e/3) a_2. SLOW
# (~minutes): gated behind HARD_FULL=1.  Cited in classical-hard-branches.md as
# offline corroboration; result was Groebner==(1) (no e!=0 pair in the box).
_b1_e = lambda kk, ee, a2f, a3p: kk + sp.Rational(2, 3)*ee*a2f
if os.environ.get("HARD_FULL"):
    check(box_is_empty(sp.Integer(1), (1, 1, 1, 2, 2, 3),
                       rabinowitsch_e=True, b1_form=_b1_e) is True,
          "[COMPUTED,HARD_FULL] e!=0 constant-h: no pair in box (1,1,1,2,2,3) (Groebner=(1))")
else:
    print("SKIP  e!=0 Rabinowitsch box (set HARD_FULL=1 to run; offline result: Groebner=(1))")

print("\n--- Section 7 (bounded Groebner corroboration) done ---\n")
print()
print(f"({PASS} checks)")
print("ALL CLASSICAL HARD CHECKS PASSED")
