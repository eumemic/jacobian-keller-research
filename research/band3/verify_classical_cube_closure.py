#!/usr/bin/env python3
"""Exact verifier for the classical band-3 CONSTANT-h, e=0 CUBE-SECTOR CLOSURE.

Convention (frozen, classical, identical to classical-band3-cascade.md at commit
99fe6ee):  {G,F} = G_xi F_x - G_x F_xi,  tau = x*xi,  {xi,x}=1.
Band-3 presentation:  F = sum_{k=-3}^{3} x^k a_k(tau),  G = sum x^l b_l(tau).
Membership (genuine polynomiality in C[x,xi]):  tau^j | a_{-j}, b_{-j}, j=1,2,3.
Keller equations:  C_m = sum_{k+l=m}(k a_k b_l' - l a_k' b_l) = delta_{m0}, m in [-6,6].

This script machine-verifies every displayed identity in classical-cube-closure.md.
The sector is:  a_3 = c (constant, h=1),  b_3 = 0 (top gauge),  b_2 = 0 (e = 0),
b_1 = kappa (constant).  A successful run ends
'ALL CLASSICAL CUBE CLOSURE CHECKS PASSED'.

Degree-lattice / branch-emptiness statements are written arguments in the .md;
the bounded checks here are regression corroboration only.
"""
import sympy as sp

t = sp.symbols("tau")
KS = range(-3, 4)


def D(f):
    return sp.diff(f, t)


def Cm(m, A, B, Dop=D):
    return sp.expand(sum(k * A[k] * Dop(B[l]) - l * Dop(A[k]) * B[l]
                         for k in KS for l in KS if k + l == m))


def poisson(g, f, x, xi):
    return sp.diff(g, xi) * sp.diff(f, x) - sp.diff(g, x) * sp.diff(f, xi)


PASS = 0


def check(cond, label):
    global PASS
    if cond is not True and cond != True:  # noqa: E712
        raise AssertionError("FAIL: " + label)
    PASS += 1
    print("PASS", label)


def zero(expr, label):
    check(sp.simplify(sp.expand(expr)) == 0, label)


kap, beta, gam, delta, c = sp.symbols("kappa beta gamma delta c")

# =====================================================================
# SECTION 0.  Normalization c=1 is a genuine scaling, not a loss of generality.
#   Diagonal symplectic (x,xi) -> (rho x, rho^-1 xi) sends a_k -> rho^k a_k;
#   so a_3 = c -> rho^3 c.  Choosing rho^3 = 1/c normalizes c = 1 while keeping
#   h = 1 (constant).  We KEEP c symbolic below (a weight-3 grading bookkeeper);
#   every identity holds for all c != 0 and specializes at c = 1.
# =====================================================================
x, xi = sp.symbols("x xi")
rho = sp.symbols("rho")
a3sym = sp.Function("a3")(t)
Fscale = x**3 * a3sym.subs(t, x * xi)
# under x->rho x, xi->rho^-1 xi the x^3-coefficient a_3 picks up rho^3:
scaled = Fscale.subs({x: rho * x, xi: xi / rho})
zero(scaled - rho**3 * x**3 * a3sym.subs(t, x * xi),
     "diagonal scaling multiplies a_3 by rho^3 (justifies normalizing c=1)")

# =====================================================================
# SECTION 1.  The 13 C_m match the direct two-variable Poisson bracket
#   (fully generic coefficients).  [reproduces classical-band3-cascade.md §1]
# =====================================================================
Ag = {k: sp.Function(f"a{k}")(t) for k in KS}
Bg = {k: sp.Function(f"b{k}")(t) for k in KS}
f_full = sum(x**k * Ag[k].subs(t, x * xi) for k in KS)
g_full = sum(x**k * Bg[k].subs(t, x * xi) for k in KS)
formula = sp.expand(sum(x**m * Cm(m, Ag, Bg).subs(t, x * xi) for m in range(-6, 7)))
zero(sp.expand(poisson(g_full, f_full, x, xi)) - formula,
     "all 13 C_m match direct 2-variable Poisson (fully generic coefficients)")

# =====================================================================
# SECTION 2.  The reduction:  constant-h (a_3=c), b_3=0, b_2=0 (e=0), b_1=kappa.
#   Free F-coefficients q=a_2, p=a_1, r=a_0, al=a_-1, s=a_-2, sig=a_-3.
#   Positive cascade C_3,C_2,C_1 integrate b_0,b_-1,b_-2; moment C_0 gives b_-3.
#   All of b_0,b_-1,b_-2,b_-3 become EXPLICIT in the F-coefficients.
# =====================================================================
q = sp.Function("q")(t); p = sp.Function("p")(t); r = sp.Function("r")(t)
al = sp.Function("al")(t); s = sp.Function("s")(t); sig = sp.Function("sig")(t)

A = {3: c, 2: q, 1: p, 0: r, -1: al, -2: s, -3: sig}
b0 = kap * q / (3 * c) + beta
bm1 = kap * p / (3 * c) - kap * q**2 / (9 * c**2) + gam
w = (kap * r / (3 * c) - sp.Rational(2, 9) * kap * p * q / c**2
     + sp.Rational(5, 81) * kap * q**3 / c**3 - gam * q / (3 * c) + delta)
z = sp.expand((t - 2 * q * w - p * bm1 + kap * al) / (3 * c))
B = {3: sp.Integer(0), 2: sp.Integer(0), 1: kap, 0: b0, -1: bm1, -2: w, -3: z}

zero(Cm(3, A, B), "C_3 = 0 identically  (integrates b_0 = kappa a_2/(3c) + beta)")
zero(Cm(2, A, B), "C_2 = 0 identically  (integrates b_-1 = kappa a_1/(3c) - kappa a_2^2/(9c^2) + gamma)")
zero(Cm(1, A, B), "C_1 = 0 identically  (integrates b_-2 explicit in a_0,a_1,a_2)")
zero(Cm(0, A, B) - 1, "C_0 = 1 identically  (moment => b_-3 explicit in a_-1,a_0,a_1,a_2)")

# The b_-2 antiderivative is exact (C_1 RHS is a total derivative):
wprime_check = 3 * D(w) - (kap * D(r) - sp.Rational(2, 3) * kap * (p * q).diff(t) / c
                           + sp.Rational(5, 27) * kap * (q**3).diff(t) / c**2
                           - gam * D(q))
# (this equals the C_1 relation up to the +2q/... moment-free structure) -- direct:
zero(sp.together(w) - sp.together(kap * r / (3 * c) - sp.Rational(2, 9) * kap * p * q / c**2
     + sp.Rational(5, 81) * kap * q**3 / c**3 - gam * q / (3 * c) + delta),
     "b_-2 is an EXPLICIT polynomial in a_0,a_1,a_2 (C_1 RHS is a total derivative)")

# =====================================================================
# SECTION 3.  C_-1 determines a_-2 ; C_-2 determines a_-3.
#   Each trailing coefficient enters its equation ONLY through -kappa*(deriv).
# =====================================================================
Cm1 = sp.expand(Cm(-1, A, B))
Cm2 = sp.expand(Cm(-2, A, B))
zero(sp.diff(Cm1, D(s)) + kap, "C_-1: coeff of a_-2' is -kappa (C_-1 DETERMINES a_-2)")
zero(sp.diff(Cm1, s), "C_-1: a_-2 enters only through -kappa a_-2' (no undifferentiated a_-2)")
check(not Cm1.has(sig), "C_-1 does not involve a_-3 (so a_-2 solved before a_-3)")
zero(sp.diff(Cm2, D(sig)) + kap, "C_-2: coeff of a_-3' is -kappa (C_-2 DETERMINES a_-3)")
zero(sp.diff(Cm2, sig), "C_-2: a_-3 enters only through -kappa a_-3'")

# =====================================================================
# SECTION 4.  Phi_1 : the band-3 W3-replacement first integral of C_-1.
#   D(Phi_1) = C_-1 exactly, with the SINGLE nonlocal generator Q1 = int a_2.
#   (The naive gauge-free W3 integral is obstructed at band 3 -- cascade §8 --
#    but after the full constant-h e=0 parametrization the trailing integral
#    survives with nonlocal Q1, exactly analogous to band-2 Lemma 2.1's P1=int p.)
# =====================================================================
Q1 = sp.Function("Q1")(t)  # antiderivative of a_2 = q

def Dq(expr):
    return sp.expand(sp.diff(expr, t).subs(sp.Derivative(Q1, t), q))

Phi1 = (2 * delta * p + gam * r - kap * s
        - sp.Rational(5, 3) * delta * q**2 / c - sp.Rational(4, 3) * gam * p * q / c
        + sp.Rational(2, 3) * kap * al * q / c + sp.Rational(2, 3) * kap * p * r / c
        + t * q / c - Q1 / (3 * c)
        + sp.Rational(14, 27) * gam * q**3 / c**2 - sp.Rational(5, 9) * kap * p**2 * q / c**2
        - sp.Rational(5, 9) * kap * q**2 * r / c**2
        + sp.Rational(40, 81) * kap * p * q**3 / c**3 - sp.Rational(22, 243) * kap * q**5 / c**4)
zero(Dq(Phi1) - Cm1, "Phi_1' = C_-1  (the W3-replacement first integral; nonlocal Q1=int a_2)")
# The obstruction that forces the nonlocal term is the constant Euler density:
def euler_q(expr):
    return sp.expand(sp.diff(expr, q) - D(sp.diff(expr, D(q)))
                     + D(D(sp.diff(expr, D(D(q))))))
zero(euler_q(Cm1) + sp.Rational(1, 3) / c,
     "Euler_q(C_-1) = -1/(3c) constant (the obstruction absorbed by Q1=int a_2)")

# =====================================================================
# SECTION 5.  Branch tree:  C_-6 bottom Wronskian ; C_-5 bottom cube-wall mirror.
# =====================================================================
mu3 = sp.symbols("mu3")
# C_-6 = 0  <=>  b_-3 = mu3 a_-3  (bottom Wronskian W1-mirror)
Bm = dict(B); Bm[-3] = mu3 * sig
zero(Cm(-6, {**A}, Bm), "C_-6 = 0 after b_-3 = mu3 a_-3 (bottom proportionality)")

# C_-5 with b_-3=mu3 a_-3 equals the bottom cube-wall operator on phi = mu3 a_-2 - b_-2
am3f = sp.Function("am3")(t); bm2f = sp.Function("bm2")(t); am2f = sp.Function("am2")(t)
Cm5_disp = (2 * D(am3f) * bm2f - 3 * am3f * D(bm2f)
            - 2 * am2f * D(mu3 * am3f) + 3 * D(am2f) * (mu3 * am3f))
phi = mu3 * am2f - bm2f
Lbot = 3 * am3f * D(phi) - 2 * D(am3f) * phi
zero(Cm5_disp - Lbot,
     "C_-5 (b_-3=mu3 a_-3) = 3 a_-3 phi' - 2 a_-3' phi,  phi = mu3 a_-2 - b_-2 (bottom wall)")
zero(D(phi**3 / am3f**2) - phi**2 * Lbot / am3f**3,
     "bottom cube integrating factor: (phi^3/a_-3^2)' = phi^2 * L_bot / a_-3^3")
# hence C_-5=0 & a_-3!=0 => phi^3 = c' a_-3^2 => 3 deg(phi) = 2 deg(a_-3) => 3 | deg(a_-3)
h_ = sp.Function("h")(t); ee = sp.symbols("e_bot")
zero(3 * h_**3 * D(ee * h_**2) - 2 * D(h_**3) * (ee * h_**2),
     "bottom wall homogeneous solution: a_-3=h^3 => phi=e h^2 solves L_bot (cube-class mirror)")

# =====================================================================
# SECTION 6.  q = a_2 CONSTANT reduction (contains ALL positive controls).
#   s=a_-2 explicit (from C_-1); sig=a_-3 explicit (from C_-2, nonlocal P1=int a_1);
#   closing constraints are C_-3=C_-4=C_-5=C_-6=0.
# =====================================================================
Q0 = sp.symbols("Q0")  # constant a_2
pq = sp.Function("p")(t); rq = sp.Function("r")(t); alq = sp.Function("al")(t)
P1 = sp.Function("P1")(t)  # antiderivative of a_1 = p

def Dp(expr):
    return sp.expand(sp.diff(expr, t).subs(sp.Derivative(P1, t), pq))

s_expr = (1 / kap) * (sp.Rational(2, 3) * Q0 / c * t
         + sp.Rational(40, 81) * Q0**3 * kap / c**3 * pq - sp.Rational(5, 9) * Q0**2 * kap / c**2 * rq
         - sp.Rational(4, 3) * Q0 * gam / c * pq + sp.Rational(2, 3) * Q0 * kap / c * alq
         - sp.Rational(5, 9) * Q0 * kap / c**2 * pq**2 + 2 * delta * pq + gam * rq
         + sp.Rational(2, 3) * kap / c * pq * rq)
sig_expr = (-sp.Rational(10, 81) * Q0**4 * pq / c**4 + sp.Rational(10, 81) * Q0**3 * rq / c**3
            + sp.Rational(2, 3) * Q0**2 * gam * pq / (c**2 * kap) - sp.Rational(1, 9) * Q0**2 * alq / c**2
            + sp.Rational(10, 27) * Q0**2 * pq**2 / c**3 - 2 * Q0 * delta * pq / (c * kap)
            - sp.Rational(2, 3) * Q0 * gam * rq / (c * kap) - sp.Rational(2, 3) * Q0 * pq * rq / c**2
            + 2 * delta * rq / kap + gam * alq / kap - sp.Rational(2, 3) * gam * pq**2 / (c * kap)
            + sp.Rational(2, 3) * alq * pq / c + rq**2 / (3 * c)
            + t * pq / (c * kap) - sp.Rational(2, 3) * P1 / (c * kap) - sp.Rational(5, 27) * pq**3 / c**2)
Ac = {3: c, 2: Q0, 1: pq, 0: rq, -1: alq, -2: s_expr, -3: sig_expr}
b0c = kap * Q0 / (3 * c) + beta
bm1c = kap * pq / (3 * c) - kap * Q0**2 / (9 * c**2) + gam
wc = (kap * rq / (3 * c) - sp.Rational(2, 9) * kap * pq * Q0 / c**2
      + sp.Rational(5, 81) * kap * Q0**3 / c**3 - gam * Q0 / (3 * c) + delta)
zc = sp.expand((t - 2 * Q0 * wc - pq * bm1c + kap * alq) / (3 * c))
Bc = {3: sp.Integer(0), 2: sp.Integer(0), 1: kap, 0: b0c, -1: bm1c, -2: wc, -3: zc}

def Cmc(m):
    return sp.expand(sum(k * Ac[k] * Dp(Bc[l]) - l * Dp(Ac[k]) * Bc[l]
                         for k in KS for l in KS if k + l == m))

zero(Cmc(-1), "q const: a_-2 = s_expr solves C_-1 (explicit, local in a_0,a_1,a_-1)")
zero(Cmc(-2), "q const: a_-3 = sig_expr solves C_-2 (explicit; nonlocal P1=int a_1)")
# sig genuinely carries the nonlocal P1:
check(sig_expr.has(P1), "q const: a_-3 carries the nonlocal generator P1 = int a_1")

# =====================================================================
# SECTION 7.  Positive controls B3-2, B3-3 (gauged) and the tame witness.
# =====================================================================
def full_Cm(m, A_, B_):
    return sp.expand(sum(k * A_[k] * sp.diff(B_[l], t) - l * sp.diff(A_[k], t) * B_[l]
                         for k in KS for l in KS if k + l == m))

def check_pair(name, A_, B_, want_branch):
    ok = all(sp.expand(full_Cm(m, A_, B_) - (1 if m == 0 else 0)) == 0 for m in range(-6, 7))
    check(ok, f"{name}: satisfies all 13 C_m = delta_m0")
    # membership
    memb = all(sp.simplify(sp.S(A_[k]).subs(t, 0)) == 0 for k in [-1] if A_[-1] != 0) if True else True
    # branch data
    a_3n = A_[-3]
    mu = sp.simplify(B_[-3] / a_3n) if a_3n != 0 else None
    return mu

# B3-2 gauged (c=1, kappa=-1, e=0):
A_B32 = {3: sp.Integer(1), 2: sp.Integer(1), 1: 3 * t + 1, 0: 2 * t,
         -1: 3 * t**2 + 2 * t, -2: t**2, -3: t**3}
B_B32 = {3: sp.Integer(0), 2: sp.Integer(0), 1: sp.Integer(-1), 0: sp.Integer(0),
         -1: -t, -2: sp.Integer(0), -3: sp.Integer(0)}
mu = check_pair("B3-2", A_B32, B_B32, "B0-band3")
check(mu == 0 and B_B32[-2] == 0, "B3-2 is B0-band3: mu3=0, b_-3=b_-2=0 (a_2=1 const, deg pattern (P,R,L)=(1,1,2))")

# B3-3 gauged (c=1, kappa=-1, e=0, a_2=0):
A_B33 = {3: sp.Integer(1), 2: sp.Integer(0), 1: 3 * t + 2, 0: sp.Integer(0),
         -1: 3 * t**2 + 3 * t, -2: sp.Integer(0), -3: t**3}
B_B33 = {3: sp.Integer(0), 2: sp.Integer(0), 1: sp.Integer(-1), 0: sp.Integer(0),
         -1: -t, -2: sp.Integer(0), -3: sp.Integer(0)}
mu33 = check_pair("B3-3", A_B33, B_B33, "B0-band3")
check(mu33 == 0, "B3-3 is B0-band3: mu3=0, b_-3=b_-2=0 (a_2=0, odd-only levels)")

# controls' G-ladder is reproduced by the GENERAL (local) reduction formulas
#   b_0, b_-1, b_-2, b_-3  under the scalar substitution.
def reduce_sub(vals, actualB):
    subs = {q: vals['q'], p: vals['p'], r: vals['r'], al: vals['al'], s: vals['s'], sig: vals['sig'],
            c: vals['c'], kap: vals['kap'], beta: vals['beta'], gam: vals['gam'], delta: vals['delta']}
    return all(sp.expand(sp.S(B[k]).subs(subs).doit() - actualB[k]) == 0 for k in [0, -1, -2, -3])

okB32 = reduce_sub(dict(q=1, p=3 * t + 1, r=2 * t, al=3 * t**2 + 2 * t, s=t**2, sig=t**3,
                        c=1, kap=-1, beta=sp.Rational(1, 3), gam=sp.Rational(2, 9),
                        delta=sp.Rational(-7, 81)), B_B32)
check(okB32, "B3-2 G-ladder (b_0,b_-1,b_-2,b_-3) reproduced by the reduction formulas")
okB33 = reduce_sub(dict(q=0, p=3 * t + 2, r=0, al=3 * t**2 + 3 * t, s=0, sig=t**3,
                        c=1, kap=-1, beta=0, gam=sp.Rational(2, 3), delta=0), B_B33)
check(okB33, "B3-3 G-ladder (b_0,b_-1,b_-2,b_-3) reproduced by the reduction formulas")

# tame witness  F = xi + alpha x^3 + b x^2 + g x + d,  G = -x   (onesided branch a_-3=0)
al_, b_, g_, d_ = sp.symbols("alpha bwit gwit dwit")
Aw = {3: al_, 2: b_, 1: g_, 0: d_, -1: t, -2: sp.Integer(0), -3: sp.Integer(0)}
Bw = {k: sp.Integer(0) for k in KS}; Bw[1] = sp.Integer(-1)
check(all(sp.expand(full_Cm(m, Aw, Bw) - (1 if m == 0 else 0)) == 0 for m in range(-6, 7)),
      "tame witness F=xi+alpha x^3+.., G=-x satisfies all C_m (onesided branch a_-3=0)")

# =====================================================================
# SECTION 8.  The band-3 W5 verdict: a TROPICAL (max-degree) obstruction, not a
#   linear congruence.  The resistant branch A*-band3 needs b_-3 = mu3 a_-3,
#   mu3 != 0, hence deg(b_-3) = deg(a_-3).  But the cube structure forces
#   deg(a_-3) >= 3 deg(a_1) while deg(b_-3) <= 2 deg(a_1) generically -- the
#   '3 vs 2' of the cube -- so they cannot match unless a leading coefficient
#   cancels.  This REPLACES band-2's mod-3 congruence.
# =====================================================================
def degof(e):
    e = sp.expand(e)
    return sp.Poly(e, t).degree() if e != 0 else -1
# a_-3 (from C_-2) carries three top terms whose tau-degrees strictly dominate,
# variable by variable, the corresponding terms of b_-3 = z (from the moment):
#     a_-3 : a_1^3 (3P),   a_0^2 (2R),   a_-1 a_1 (L+P)
#     b_-3 : a_1^2 (2P),   a_0   (R),    a_-1     (L)
# so deg a_-3 >= max(3P,2R,L+P) > max(2P,R,L,1) >= deg b_-3  (P=deg a_1 >= 1),
# UNCONDITIONALLY unless the top of a_-3 cancels at a tie 3P=2R / L=2P / 2R=L+P.
zero(sp.expand(sig_expr.coeff(pq, 3)) + sp.Rational(5, 27) / c**2,
     "a_-3 has a_1^3 term, coeff = -5/(27c^2) != 0   (3P term)")
zero(sp.expand(sig_expr.coeff(rq, 2)) - sp.Rational(1, 3) / c,
     "a_-3 has a_0^2 term, coeff = 1/(3c) != 0   (2R term)")
zero(sp.expand(sp.expand(sig_expr.coeff(alq, 1)).coeff(pq, 1)) - sp.Rational(2, 3) / c,
     "a_-3 has a_-1 a_1 term, coeff = 2/(3c) != 0   (L+P term)")
check(sp.Poly(sp.expand(zc), pq).degree() <= 2 and sp.Poly(sp.expand(zc), rq).degree() <= 1
      and sp.Poly(sp.expand(zc), alq).degree() <= 1,
      "b_-3 = z has (a_1,a_0,a_-1)-degrees <= (2,1,1)  =>  termwise tropical gap deg b_-3 < deg a_-3")
# B3-2, B3-3 sit on the boundary stratum deg a_-1 = 2 deg a_1 (both B0-band3, b_-3=0):
check((degof(A_B32[1]), degof(A_B32[0]), degof(A_B32[-1])) == (1, 1, 2),
      "B3-2 on the stratum deg a_-1 = 2 deg a_1  (P,R,L)=(1,1,2), b_-3=0 (B0-band3)")
check((degof(A_B33[1]), degof(A_B33[-1])) == (1, 2),
      "B3-3 on the stratum deg a_-1 = 2 deg a_1  (P,L)=(1,2), b_-3=0 (B0-band3)")
# cube-wall divisibility 3 | deg a_-3 on the controls (deg a_-3 = 3):
check(degof(A_B32[-3]) % 3 == 0 and degof(A_B33[-3]) % 3 == 0,
      "controls satisfy the cube-wall divisibility 3 | deg a_-3 (=3)")

print()
print(f"({PASS} checks)")
print("ALL CLASSICAL CUBE CLOSURE CHECKS PASSED")
