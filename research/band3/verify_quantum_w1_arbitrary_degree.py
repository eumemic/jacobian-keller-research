#!/usr/bin/env python3
"""Exact algebra certificate for the arbitrary-degree W1 obstruction.

Run: python research/band3/verify_quantum_w1_arbitrary_degree.py
Ends: ALL QUANTUM W1 ARBITRARY-DEGREE CHECKS PASSED
"""
import sympy as sp

E = sp.symbols("E")
LEVELS = range(-3, 4)


def sh(f, n):
    return sp.expand(sp.sympify(f).subs(E, E + n))


def falling(n):
    return sp.prod(E-j for j in range(n))


def q_m(A, B, m):
    """Stipulated convention Q_m=sum_(k+l=m)(b_l^[k]a_k-a_k^[l]b_l)."""
    return sp.expand(sum(sh(B[l], k)*A[k] - sh(A[k], l)*B[l]
                         for k in LEVELS for l in LEVELS if k+l == m))


def potential(A, B):
    return sp.expand(sum(sh(A[k], j-k)*sh(B[-k], j)
                         - sh(B[k], j-k)*sh(A[-k], j)
                         for k in range(1, 4) for j in range(k)))


def k3(c):
    return sp.expand(sum(sh(a3, j-3)*sh(c, j) for j in range(3)))


def h2(v):
    return sp.expand(sum(sh(b2, j-2)*sh(v, j) for j in range(2)))


def ev(name, n):
    return sp.symbols(f"{name}_{n}")


def check_zero(value, label):
    residual = sp.factor(sp.expand(value))
    if residual != 0:
        raise AssertionError(f"{label}: residual {residual}")
    print("PASS", label)


def check(condition, label):
    if not condition:
        raise AssertionError(label)
    print("PASS", label)


def solve_fixed(eq, variable, substitutions, pivot, label):
    reduced = sp.expand(eq.subs(substitutions))
    check_zero(reduced.coeff(variable)-pivot, label + " fixed rational pivot")
    rest = sp.expand(reduced-pivot*variable)
    check(not rest.has(variable), label + " is linear in selected evaluation")
    substitutions[variable] = sp.factor(-rest/pivot)


print("--- 0. stipulated formula, W1 top, full identities, and membership ---")
a3 = sp.expand(E*(E-2)*(E-4))
b2 = sp.expand((E-1)*(E-4))
check_zero(sh(b2, 3)*a3-sh(a3, 2)*b2, "Q_5 wall for W1")

# Generic polynomials make these direct polynomial identities rather than sampled values.
def poly(name, degree):
    cs = sp.symbols(f"{name}_0:{degree+1}")
    return sp.expand(sum(cs[j]*E**j for j in range(degree+1)))

A = {k: sp.Integer(0) for k in LEVELS}
B = {k: sp.Integer(0) for k in LEVELS}
A[3], B[2] = a3, b2
A[2], A[1], A[0] = poly("a2", 5), poly("a1", 5), poly("a0", 5)
A[-1] = falling(1)*poly("am1", 5)
A[-2] = falling(2)*poly("am2", 5)
A[-3] = falling(3)*poly("am3", 5)
B[1], B[0] = poly("b1", 6), poly("b0", 6)
B[-1] = falling(1)*poly("bm1", 6)
B[-2] = falling(2)*poly("bm2", 6)
B[-3] = falling(3)*poly("bm3", 5)
check_zero(B[3], "gauge b3=0")
for k in (1, 2, 3):
    check_zero(sp.Poly(A[-k], E).rem(sp.Poly(falling(k), E)).as_expr(),
               f"a_-{k} has Weyl membership factor (E)_{k}")
    check_zero(sp.Poly(B[-k], E).rem(sp.Poly(falling(k), E)).as_expr(),
               f"b_-{k} has Weyl membership factor (E)_{k}")
G = potential(A, B)
R_generic = sp.expand(sh(A[1], -1)*B[-1]-sh(B[1], -1)*A[-1]
                      + sh(A[2], -2)*B[-2]+sh(A[2], -1)*sh(B[-2], 1))
check_zero(G-(R_generic+k3(B[-3])-h2(A[-2])),
           "G=R+K3[b_-3]-H2[a_-2] with the stipulated shifts")
check_zero(G.subs(E, 0), "full Weyl membership gives G(0)=0")
check_zero(q_m(A, B, 0)-(sh(G, 1)-G),
           "Q0=(T-1)G from stipulated Q_m formula")

print("\n--- 1. arbitrary-degree boundary evaluation certificate ---")
# Work in independent evaluation variables. Membership is imposed exactly at the
# required boundary values; a rejected weakened system guards against omissions.
a2 = lambda n: ev("a2", n)
a1 = lambda n: ev("a1", n)
a0 = lambda n: ev("a0", n)
am1 = lambda n: ev("am1", n)
am2 = lambda n: ev("am2", n)
b1 = lambda n: ev("b1", n)
b0 = lambda n: ev("b0", n)
bm1 = lambda n: ev("bm1", n)
bm2 = lambda n: ev("bm2", n)

membership = {am1(0): 0, am2(0): 0, am2(1): 0,
              bm1(0): 0, bm2(0): 0, bm2(1): 0}
check(len(membership) == 6, "all six middle-level boundary membership equations are present")
weakened = dict(membership)
del weakened[bm2(1)]
R_at_0_raw = a1(-1)*bm1(0)-b1(-1)*am1(0)+a2(-2)*bm2(0)+a2(-1)*bm2(1)
check(sp.expand(R_at_0_raw.subs(weakened)) != 0,
      "missing-boundary guard rejects omission of b_-2(1)=0")
check_zero(R_at_0_raw.subs(membership), "membership gives R(0)=0")

# Canonical non-filler residual: G1 plus the solved b_-2 half of G2.
def R(n):
    return sp.expand(a1(n-1)*bm1(n)-b1(n-1)*am1(n)
                     + a2(n-2)*bm2(n)+a2(n-1)*bm2(n+1))

lam_eight = sp.expand(
    a1(2)*bm1(3)-b1(2)*am1(3)+a2(1)*bm2(3)+a2(2)*bm2(4)
    -a1(3)*bm1(4)+b1(3)*am1(4)-a2(2)*bm2(4)-a2(3)*bm2(5)
    +a1(4)*bm1(5)-b1(4)*am1(5)+a2(3)*bm2(5)+a2(4)*bm2(6))
check_zero((R(3)-R(4)+R(5)-R(0)).subs(membership)-lam_eight,
           "lambda0(R) is the displayed eight boundary terms after cancellations")

# Evaluation form of Q_m. Only values touched by the certificate are generated.
def value(item, n):
    if item == 0:
        return sp.Integer(0)
    if item == a3:
        return sp.expand(a3.subs(E, n))
    if item == b2:
        return sp.expand(b2.subs(E, n))
    return ev(item, n)

AE = {3: a3, 2: "a2", 1: "a1", 0: "a0", -1: "am1", -2: "am2", -3: "am3"}
BE = {3: 0, 2: b2, 1: "b1", 0: "b0", -1: "bm1", -2: "bm2", -3: "bm3"}

def q_eval(m, n):
    return sp.expand(sum(value(BE[l], n+k)*value(AE[k], n)
                         - value(AE[k], n+l)*value(BE[l], n)
                         for k in LEVELS for l in LEVELS if k+l == m))

subs = dict(membership)
# Five Q4 evaluations: four displayed b1 values and the compatibility relation.
solve_fixed(q_eval(4, 0), b1(0), subs, -3, "Q4(0)")
solve_fixed(q_eval(4, 2), b1(2), subs, 3, "Q4(2)")
solve_fixed(q_eval(4, 1), b1(4), subs, 3, "Q4(1)")
solve_fixed(q_eval(4, 3), b1(6), subs, -3, "Q4(3)")
check_zero(subs[b1(0)]+sp.Rational(2,3)*a2(0)+sp.Rational(4,3)*a2(2), "displayed b1(0)")
check_zero(subs[b1(2)]+sp.Rational(2,3)*a2(4), "displayed b1(2)")
check_zero(subs[b1(4)]-sp.Rational(2,3)*a2(1), "displayed b1(4)")
check_zero(subs[b1(6)]-sp.Rational(4,3)*a2(3)-sp.Rational(2,3)*a2(5), "displayed b1(6)")
q44 = sp.factor(q_eval(4, 4).subs(subs))
check_zero(q44+10*(a2(1)-a2(4)), "Q4(4) gives a2(4)=a2(1)")

# Preserve the relation as an ideal generator while performing the fixed-pivot
# sequence. The normal form is computed after all substitutions.
relation = a2(1)-a2(4)
for n, variable, pivot in [(0, a1(2), -4), (1, a1(1), -2),
                            (3, a1(5), 2), (4, a1(4), 4)]:
    solve_fixed(q_eval(3, n), variable, subs, pivot, f"Q3({n})")
q32 = sp.factor(q_eval(3, 2).subs(subs))
check_zero(q32+sp.Rational(1, 3)*(a2(0)-a2(5))*relation,
           "exact residual Q3(2)=-1/3(a2(0)-a2(5))(a2(1)-a2(4))")
check_zero(sp.rem(q32, relation, a2(1)), "Q3(2) vanishes modulo a2(1)-a2(4)")

for n, variable, pivot in [(0, a0(0), 4), (1, bm1(4), 3),
                            (2, a0(4), 2), (3, a0(5), 2)]:
    solve_fixed(q_eval(2, n), variable, subs, pivot, f"Q2({n})")
for n, variable, pivot in [(0, am1(2), -4), (1, bm2(4), 3),
                            (2, am1(4), 2), (3, am1(5), 2)]:
    solve_fixed(q_eval(1, n), variable, subs, pivot, f"Q1({n})")

lam_reduced = sp.factor(lam_eight.subs(subs))
check_zero(sp.rem(lam_reduced, relation, a2(1)),
           "final normal form lambda0(R)=0 modulo cascade and membership")

print("\n--- 2. retained b0 kernel and direct filler annihilators ---")
kappa = sp.symbols("kappa")
B_shifted = dict(B)
B_shifted[0] = sp.expand(B[0]+kappa)
check_zero(q_m(A, B_shifted, 3)-q_m(A, B, 3), "constant b0 is an operator kernel for Q3")
check_zero(potential(A, B_shifted)-potential(A, B), "constant b0 does not enter G or R")
for m in range(1, 5):
    check_zero(q_m(A, B_shifted, m)-q_m(A, B, m),
               f"constant-b0 shift cancels in Q{m}")

# Independent values make these checks valid for every admissible filler, rather
# than only for a bounded polynomial ansatz.  The formulas are built directly
# from the displayed W1 a3,b2 and lambda0(f)=f(3)-f(4)+f(5)-f(0).
c = lambda n: ev("c", n)
v = lambda n: ev("v", n)

def k3_eval(n):
    return sp.expand(sum(a3.subs(E, n+j-3)*c(n+j) for j in range(3)))


def h2_eval(n):
    return sp.expand(sum(b2.subs(E, n+j-2)*v(n+j) for j in range(2)))


c_membership = {c(0): 0, c(1): 0, c(2): 0}
v_membership = {v(0): 0, v(1): 0}
lam_k3 = k3_eval(3)-k3_eval(4)+k3_eval(5)-k3_eval(0)
lam_h2 = h2_eval(3)-h2_eval(4)+h2_eval(5)-h2_eval(0)
check_zero(lam_k3.subs(c_membership),
           "lambda0(K3[c])=0 for arbitrary c with c(0)=c(1)=c(2)=0")
check_zero(lam_h2.subs(v_membership),
           "lambda0(H2[v])=0 for arbitrary v with v(0)=v(1)=0")
check_zero((lam_k3-lam_h2).subs({**c_membership, **v_membership}),
           "Im(Phi) is contained in ker(lambda0) directly")

print("\n--- 3. unit bridge and final direct obstruction ---")
check_zero(sh(E, 1)-E-1, "(T-1)E=1")
n = sp.symbols("n", integer=True, positive=True)
check_zero(sp.binomial(n, n-1)-n,
           "arbitrary-degree leading coefficient of (E+1)^n-E^n is n")
constant = sp.symbols("constant")
check(sp.solve(sp.Eq((E+constant).subs(E, 0), 0), constant) == [0],
      "G(0)=0 kills the additive constant in G=E+constant")
check_zero(E.subs(E, 3)-E.subs(E, 4)+E.subs(E, 5)-E.subs(E, 0)-4,
           "lambda0(E)=4")
check_zero(sp.rem(4-lam_reduced, relation, a2(1))-4,
           "lambda0(E-R)=4 after the cascade normal form")
check(4 != 0, "direct annihilator excludes E-R from Im(Phi)")

print("\nScope: normalized W1 datum only; no AP-family, W2, full Band-3, DC1, or JC2 closure; no Weyl pair constructed.")
print("ALL QUANTUM W1 ARBITRARY-DEGREE CHECKS PASSED")
