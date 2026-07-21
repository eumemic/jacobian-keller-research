#!/usr/bin/env python3
"""Exact algebra certificate for the arbitrary-degree, symbolic-r AP lambda obstruction.

Extends the fixed normalized W1 (r=0) theorem of quantum-w1-arbitrary-degree.py
across the ENTIRE step-2 arithmetic-progression exotic family

    a_3 = (E-r)(E-r-2)(E-r-4),   b_2 = (E-r-1)(E-r-4),   b_3 = 0

with the wall parameter r kept symbolic.  The functional

    lambda_r(f) = f(r+3) - f(r+4) + f(r+5) - f(0)

annihilates Im Phi(r) for every r, and lambda_r(R) = 0 on the positive-cascade
variety for every r (a degree-free certificate with r-independent rational
pivots).  Since lambda_r(E) = r+4, this gives lambda_r(E-R) = r+4, nonzero for
every r != -4.  Hence NO polynomial cascade of any coefficient degree can satisfy
Q_4=Q_3=Q_2=Q_1=0, genuine Weyl membership, and Q_0=1 for any AP member with
r != -4.

The single exceptional value r = -4 is exactly the W2 member; there the functional
degenerates (lambda_{-4}(E)=0, Im Phi jumps to codimension three) and the
obstruction reduces to the moment-unit slope question.  That reduction, and the
finite-degree infeasibility already certified for W2 at d=1,2, are recorded here
honestly and NOT claimed as an arbitrary-degree closure.

Companion memo: quantum-ap-lambda.md.  Overnight base commit: e4e704f.

Run:  python research/band3/verify_quantum_ap_lambda.py
Ends: ALL QUANTUM AP LAMBDA CHECKS PASSED
"""
import sympy as sp

E, r = sp.symbols("E r")
LEVELS = range(-3, 4)


def sh(f, n):
    return sp.expand(sp.sympify(f).subs(E, E + n))


def falling(n):
    return sp.prod(E - j for j in range(n))


def poly(name, degree):
    cs = sp.symbols(f"{name}_0:{degree+1}")
    return sp.expand(sum(cs[j] * E**j for j in range(degree + 1)))


def q_m(A, B, m):
    """Stipulated convention Q_m = sum_(k+l=m) (b_l^[k] a_k - a_k^[l] b_l)."""
    return sp.expand(sum(sh(B[l], k) * A[k] - sh(A[k], l) * B[l]
                         for k in LEVELS for l in LEVELS if k + l == m))


def potential(A, B):
    return sp.expand(sum(sh(A[k], j - k) * sh(B[-k], j)
                         - sh(B[k], j - k) * sh(A[-k], j)
                         for k in range(1, 4) for j in range(k)))


def check(condition, label):
    if not condition:
        raise AssertionError(label)
    print("PASS", label)


def check_zero(value, label):
    residual = sp.factor(sp.expand(value))
    if residual != 0:
        raise AssertionError(f"{label}: residual {residual}")
    print("PASS", label)


# =====================================================================
print("--- 0. AP wall, W1/W2 identification, and stipulated identities ---")
# =====================================================================
a3 = sp.expand((E - r) * (E - r - 2) * (E - r - 4))
b2 = sp.expand((E - r - 1) * (E - r - 4))
check_zero(sh(b2, 3) * a3 - sh(a3, 2) * b2,
           "AP top {r,r+2,r+4}: b_2 solves the Q_5 wall for ALL symbolic r")
check_zero(a3.subs(r, 0) - E * (E - 2) * (E - 4), "W1 is the AP member r=0 (a_3)")
check_zero(b2.subs(r, 0) - (E - 1) * (E - 4), "W1 is the AP member r=0 (b_2)")
check_zero(a3.subs(r, -4) - E * (E + 2) * (E + 4), "W2 is the AP member r=-4 (a_3)")
check_zero(b2.subs(r, -4) - E * (E + 3), "W2 is the AP member r=-4 (b_2)")

# a_3, b_2 evaluated at the r-shifted node r+j are r-INDEPENDENT constants:
for j in range(-2, 8):
    check_zero(a3.subs(E, r + j) - j * (j - 2) * (j - 4),
               f"a_3(r+{j}) = {j*(j-2)*(j-4)} is r-independent")
    check_zero(b2.subs(E, r + j) - (j - 1) * (j - 4),
               f"b_2(r+{j}) = {(j-1)*(j-4)} is r-independent")

# Generic polynomials with symbolic r make the decomposition identities exact.
A = {k: sp.Integer(0) for k in LEVELS}
B = {k: sp.Integer(0) for k in LEVELS}
A[3], B[2] = a3, b2
A[2], A[1], A[0] = poly("a2", 5), poly("a1", 5), poly("a0", 5)
A[-1] = falling(1) * poly("am1", 5)
A[-2] = falling(2) * poly("am2", 5)
A[-3] = falling(3) * poly("am3", 5)
B[1], B[0] = poly("b1", 6), poly("b0", 6)
B[-1] = falling(1) * poly("bm1", 6)
B[-2] = falling(2) * poly("bm2", 6)
B[-3] = falling(3) * poly("bm3", 5)
check_zero(B[3], "gauge b_3 = 0")
for k in (1, 2, 3):
    check_zero(sp.Poly(A[-k], E).rem(sp.Poly(falling(k), E)).as_expr(),
               f"a_-{k} carries the Weyl membership factor (E)_{k}")
    check_zero(sp.Poly(B[-k], E).rem(sp.Poly(falling(k), E)).as_expr(),
               f"b_-{k} carries the Weyl membership factor (E)_{k}")


def k3(c):
    return sp.expand(sum(sh(a3, j - 3) * sh(c, j) for j in range(3)))


def h2(v):
    return sp.expand(sum(sh(b2, j - 2) * sh(v, j) for j in range(2)))


G = potential(A, B)
R_generic = sp.expand(sh(A[1], -1) * B[-1] - sh(B[1], -1) * A[-1]
                      + sh(A[2], -2) * B[-2] + sh(A[2], -1) * sh(B[-2], 1))
check_zero(G - (R_generic + k3(B[-3]) - h2(A[-2])),
           "G = R + K_3[b_-3] - H_2[a_-2]  (symbolic r, stipulated shifts)")
check_zero(G.subs(E, 0), "full Weyl membership gives G(0)=0  (symbolic r)")
check_zero(q_m(A, B, 0) - (sh(G, 1) - G),
           "Q_0 = (T-1)G from the stipulated Q_m formula (symbolic r)")
check_zero((sh(E, 1) - E) - 1, "(T-1)E = 1")
check_zero((E.subs(E, r + 3) - E.subs(E, r + 4) + E.subs(E, r + 5) - E.subs(E, 0))
           - (r + 4), "lambda_r(E) = r+4  (symbolic r)")


# =====================================================================
print("\n--- 1. filler annihilation: Im Phi(r) subset ker(lambda_r) for ALL r ---")
# =====================================================================
# Generic degree-5 C,V with membership built in provide a symbolic-r regression.
# The adjacent nodewise formulas are the degree-free proof for arbitrary fillers.
Cgen = poly("C", 5)
Vgen = poly("V", 5)
cfill = sp.expand(falling(3) * Cgen)
vfill = sp.expand(falling(2) * Vgen)


def lam_r(f):
    f = sp.sympify(f)
    return sp.expand(f.subs(E, r + 3) - f.subs(E, r + 4)
                     + f.subs(E, r + 5) - f.subs(E, 0))


check_zero(lam_r(k3(cfill)),
           "lambda_r(K_3[(E)_3 C]) = 0 for generic C, symbolic r")
check_zero(lam_r(h2(vfill)),
           "lambda_r(H_2[(E)_2 V]) = 0 for generic V, symbolic r")
check_zero(lam_r(k3(cfill) - h2(vfill)),
           "generic degree-5 regression: Phi(r) lies in ker(lambda_r)")
# displayed boundary-value formulas (a_3(r+j), b_2(r+j) are r-independent):
kc = k3(cfill)
hv = h2(vfill)
check_zero(kc.subs(E, r + 3) - 3 * cfill.subs(E, r + 4), "K_3[c](r+3) = 3 c(r+4)")
check_zero(kc.subs(E, r + 4) - (3 * cfill.subs(E, r + 4) - 3 * cfill.subs(E, r + 6)),
           "K_3[c](r+4) = 3 c(r+4) - 3 c(r+6)")
check_zero(kc.subs(E, r + 5) + 3 * cfill.subs(E, r + 6), "K_3[c](r+5) = -3 c(r+6)")
check_zero(hv.subs(E, r + 3) + 2 * vfill.subs(E, r + 4), "H_2[v](r+3) = -2 v(r+4)")
check_zero(hv.subs(E, r + 4) - (-2 * vfill.subs(E, r + 4) - 2 * vfill.subs(E, r + 5)),
           "H_2[v](r+4) = -2 v(r+4) - 2 v(r+5)")
check_zero(hv.subs(E, r + 5) + 2 * vfill.subs(E, r + 5), "H_2[v](r+5) = -2 v(r+5)")
# Both Phi blocks vanish at E=0 too, so ev_0 also annihilates Im Phi(r):
check_zero(k3(cfill).subs(E, 0), "ev_0 annihilates K_3[(E)_3 C]")
check_zero(h2(vfill).subs(E, 0), "ev_0 annihilates H_2[(E)_2 V]")
check_zero(lam_r(E) - (r + 4),
           "lambda_r(E)=r+4: independent of ev_0 on Im Phi(r) exactly when r!=-4")


# =====================================================================
print("\n--- 2. R(0)=0 and the r-shifted boundary reduction of lambda_r(R) ---")
# =====================================================================
# Offset evaluation symbols: value of a coefficient at the r-shifted node r+j.
# The wall values a_3(r+j)=j(j-2)(j-4), b_2(r+j)=(j-1)(j-4) are r-INDEPENDENT
# constants, so the whole certificate below has r-independent rational pivots.
def off(name, j):
    return sp.symbols(f"{name}_{j}".replace("-", "m"))


def a3o(j):
    return sp.Integer(j * (j - 2) * (j - 4))


def b2o(j):
    return sp.Integer((j - 1) * (j - 4))


a2o = lambda j: off("a2", j)
a1o = lambda j: off("a1", j)
a0o = lambda j: off("a0", j)
am1o = lambda j: off("am1", j)
am2o = lambda j: off("am2", j)
am3o = lambda j: off("am3", j)
b1o = lambda j: off("b1", j)
b0o = lambda j: off("b0", j)
bm1o = lambda j: off("bm1", j)
bm2o = lambda j: off("bm2", j)
bm3o = lambda j: off("bm3", j)

# R(0): evaluation of R at the ABSOLUTE node 0.  Membership at absolute 0,1 zeros
# a_-1(0)=b_-1(0)=b_-2(0)=b_-2(1)=0, giving R(0)=0 for EVERY r.
am1_abs0, bm1_abs0, bm2_abs0, bm2_abs1 = sp.symbols("am1_abs0 bm1_abs0 bm2_abs0 bm2_abs1")
a1_absm1, b1_absm1, a2_absm2, a2_absm1 = sp.symbols("a1_absm1 b1_absm1 a2_absm2 a2_absm1")
R_at_0_raw = (a1_absm1 * bm1_abs0 - b1_absm1 * am1_abs0
              + a2_absm2 * bm2_abs0 + a2_absm1 * bm2_abs1)
abs_membership = {am1_abs0: 0, bm1_abs0: 0, bm2_abs0: 0, bm2_abs1: 0}
# missing-boundary guard: dropping b_-2(1)=0 destroys the R(0)=0 conclusion.
weak = dict(abs_membership)
del weak[bm2_abs1]
check(sp.expand(R_at_0_raw.subs(weak)) != 0,
      "missing-boundary guard rejects omission of b_-2(1)=0")
check_zero(R_at_0_raw.subs(abs_membership), "membership gives R(0)=0 for every r")


def R_off(n):
    """R evaluated at the r-shifted node r+n, in offset symbols."""
    return sp.expand(a1o(n - 1) * bm1o(n) - b1o(n - 1) * am1o(n)
                     + a2o(n - 2) * bm2o(n) + a2o(n - 1) * bm2o(n + 1))


# lambda_r(R) = R(r+3)-R(r+4)+R(r+5)-R(0); R(0)=0.  Interior level-2 products
# cancel, leaving exactly eight signed boundary products.
lam_eight = sp.expand(
    a1o(2) * bm1o(3) - b1o(2) * am1o(3) + a2o(1) * bm2o(3)
    - a1o(3) * bm1o(4) + b1o(3) * am1o(4)
    + a1o(4) * bm1o(5) - b1o(4) * am1o(5) + a2o(4) * bm2o(6))
check_zero((R_off(3) - R_off(4) + R_off(5)) - lam_eight,
           "lambda_r(R) is the displayed eight boundary terms after cancellations")


# =====================================================================
print("\n--- 3. r-uniform cascade certificate: lambda_r(R) = 0 for ALL r ---")
# =====================================================================
def q_off(m, n):
    """Q_m at the r-shifted node r+n, in offset symbols (r-independent pivots)."""
    AE = {3: a3o, 2: a2o, 1: a1o, 0: a0o, -1: am1o, -2: am2o, -3: am3o}
    BE = {3: (lambda j: sp.Integer(0)), 2: b2o, 1: b1o, 0: b0o,
          -1: bm1o, -2: bm2o, -3: bm3o}
    return sp.expand(sum(BE[l](n + k) * AE[k](n) - AE[k](n + l) * BE[l](n)
                         for k in LEVELS for l in LEVELS if k + l == m))


def solve_fixed(eq, variable, subs, pivot, label):
    reduced = sp.expand(eq.subs(subs))
    check_zero(reduced.coeff(variable) - pivot, label + " fixed r-independent pivot")
    rest = sp.expand(reduced - pivot * variable)
    check(not rest.has(variable), label + " is linear in the selected evaluation")
    subs[variable] = sp.factor(-rest / pivot)


subs = {}
# Q_4 solves b_1 at offsets 0,2,4,6 with r-independent pivots -3,3,3,-3
solve_fixed(q_off(4, 0), b1o(0), subs, -3, "Q4(r+0)")
solve_fixed(q_off(4, 2), b1o(2), subs, 3, "Q4(r+2)")
solve_fixed(q_off(4, 1), b1o(4), subs, 3, "Q4(r+1)")
solve_fixed(q_off(4, 3), b1o(6), subs, -3, "Q4(r+3)")
# Q_4 compatibility relation at offset 4: a_2(r+1) = a_2(r+4)
relation = a2o(1) - a2o(4)
check_zero(sp.factor(q_off(4, 4).subs(subs)) + 10 * relation,
           "Q4(r+4) gives the relation a_2(r+4)=a_2(r+1)  (r-independent)")
# Q_3 solves a_1 at offsets 0,1,3,4 with pivots -4,-2,2,4
for n, variable, pivot in [(0, a1o(2), -4), (1, a1o(1), -2),
                           (3, a1o(5), 2), (4, a1o(4), 4)]:
    solve_fixed(q_off(3, n), variable, subs, pivot, f"Q3(r+{n})")
check_zero(sp.factor(q_off(3, 2).subs(subs))
           + sp.Rational(1, 3) * (a2o(0) - a2o(5)) * relation,
           "Q3(r+2) residual = -1/3 (a_2(r)-a_2(r+5)) (a_2(r+1)-a_2(r+4))")
# Q_2 solves a_0, b_-1 at offsets 0,1,2,3 with pivots 4,3,2,2
for n, variable, pivot in [(0, a0o(0), 4), (1, bm1o(4), 3),
                           (2, a0o(4), 2), (3, a0o(5), 2)]:
    solve_fixed(q_off(2, n), variable, subs, pivot, f"Q2(r+{n})")
# Q_1 solves a_-1, b_-2 at offsets 0,1,2,3 with pivots -4,3,2,2
for n, variable, pivot in [(0, am1o(2), -4), (1, bm2o(4), 3),
                           (2, am1o(4), 2), (3, am1o(5), 2)]:
    solve_fixed(q_off(1, n), variable, subs, pivot, f"Q1(r+{n})")

lam_reduced = sp.factor(lam_eight.subs(subs))
check_zero(sp.rem(lam_reduced, relation, a2o(1)),
           "FINAL: lambda_r(R) = 0 modulo the cascade relation, for EVERY r")
# No membership at r-shifted nodes was used: the eight-term reduction is closed
# inside the r-shifted evaluations, with only the Q_4(r+4) relation.
used = set().union(*[sp.sympify(val).free_symbols for val in subs.values()])
check(am1_abs0 not in used and bm2_abs1 not in used,
      "the eight-term reduction uses NO absolute-node membership (self-contained)")


# =====================================================================
print("\n--- 4. main obstruction: lambda_r(E-R)=r+4 != 0 for all r != -4 ---")
# =====================================================================
lam_E_minus_R = sp.expand((r + 4) - lam_reduced)  # lambda_r(E) - lambda_r(R)
check_zero(sp.rem(lam_E_minus_R, relation, a2o(1)) - (r + 4),
           "lambda_r(E-R) = r+4 after the cascade normal form")
check_zero((r + 4).subs(r, 0) - 4,
           "at r=0 this is 4, reproducing the committed normalized-W1 theorem")
check(sp.solve(sp.Eq(r + 4, 0), r) == [-4],
      "lambda_r(E-R)=r+4 vanishes ONLY at r=-4; unique exceptional locus")
print("   => for every r != -4: E-R not in Im Phi(r), so Q_0=1 is impossible at ANY degree.")


# =====================================================================
print("\n--- 5. the UNIQUE exceptional locus r=-4 (the W2 member): honest ledger ---")
# =====================================================================
# 5a. The r!=-4 obstruction degenerates: lambda_{-4}(E)=lambda_{-4}(E-R)=0.
check_zero((r + 4).subs(r, -4), "lambda_{-4}(E) = 0: the functional is degenerate at r=-4")
check_zero(sp.rem(lam_E_minus_R, relation, a2o(1)).subs(r, -4),
           "lambda_{-4}(E-R) = 0: the r!=-4 obstruction yields NO contradiction at r=-4")

# 5b. Im Phi(-4) = D*F[E], D=E(E-1)(E+1) (codimension THREE); reuse filler image.
D = sp.expand(E * (E - 1) * (E + 1))
a3_w2 = sp.expand(a3.subs(r, -4))
b2_w2 = sp.expand(b2.subs(r, -4))


def k3_w2(cc):
    return sp.expand(sum(sh(a3_w2, j - 3) * sh(cc, j) for j in range(3)))


def h2_w2(vv):
    return sp.expand(sum(sh(b2_w2, j - 2) * sh(vv, j) for j in range(2)))


check(all(sp.Poly(k3_w2(falling(3) * E**j), E).rem(sp.Poly(D, E)).is_zero
          and sp.Poly(h2_w2(falling(2) * E**j), E).rem(sp.Poly(D, E)).is_zero
          for j in range(6)),
      "finite regression j=0..5: W2 basis fillers are divisible by D=E(E-1)(E+1)")
# annihilators of D*F[E] are exactly ev_0, ev_1, ev_-1 (a poly is in D*F[E] iff it
# vanishes at 0,1,-1).  ev_0(E-R)=-R(0)=0 already.

# 5c. The certificate still forces R(1)+R(-1)=0 at r=-4 (rigorous).
#     At r=-4: node r+3=-1, r+4=0, r+5=1.  R(0)=0, so
#     lambda_{-4}(R) = R(-1) - R(0) + R(1) - R(0) = R(1)+R(-1).
check_zero(sp.rem(lam_reduced.subs(r, -4), relation, a2o(1)),
           "at r=-4 the certificate gives R(1)+R(-1) = 0  (rigorous, arbitrary degree)")

# 5d. R(1) = G(1) = slope = Q_0(0) (fillers vanish at +-1 at r=-4): the residual
#     value E-R must miss Im Phi(-4) is EXACTLY R(1)=1, i.e. the moment slope = 1.
Gw = potential(A, B).subs(r, -4)
Rw = R_generic.subs(r, -4)
Q0w = q_m(A, B, 0).subs(r, -4)
check_zero(sp.expand(Rw.subs(E, 1) - Gw.subs(E, 1)),
           "r=-4: R(1) = G(1) (the two filler blocks vanish at E=1)")
check_zero(sp.expand(Gw.subs(E, 1) - Q0w.subs(E, 0)),
           "r=-4: slope G(1) = Q_0(0)")
# so E-R in Im Phi(-4)  <=>  R(1)=1  <=>  Q_0(0)=1 -- CONSISTENT with Q_0=1.
R1s, Rm1s, al, be, ga = sp.symbols("R1s Rm1s alpha beta gamma")
mu_EmR = sp.expand(al * (0) + be * (1 - R1s) + ga * (-1 - Rm1s))
check_zero(mu_EmR.subs(Rm1s, -R1s) - (be - ga) * (1 - R1s),
           "every Im Phi(-4) annihilator sends E-R to (beta-gamma)(1-R(1))")

# 5e. The r-shifted certificate that PINS lambda_r(R) for r!=-4 does NOT pin R(1)
#     at r=-4: its relaxed shifted-evaluation normal form for R(1)=R(r+5) is
#     nonzero.  This proves only that this certificate does not force R(1)=0; it
#     does not prove value-one achievability or a free modulus in the actual
#     polynomial cascade.  The explicit W2-DECISIVE family supplies that proof.
R1_reduced = sp.factor(sp.rem(R_off(5).subs(subs), relation, a2o(1)))
lam_eight_reduced = sp.rem(lam_eight.subs(subs), relation, a2o(1))
check(R1_reduced != 0 and bool(R1_reduced.free_symbols),
      "the relaxed certificate leaves a nonzero symbolic R(1) normal form (no value-one inference)")
check(sp.expand(lam_eight_reduced) == 0 and R1_reduced != 0,
      "contrast: the alternating combination lambda_r(R) is pinned to 0, but this certificate does not pin R(1)")
print("   => at r=-4 central completion requires the moment-unit slope=1,")
print("      which this functional certificate does NOT resolve; see W2-DECISIVE.")

# 5f. Finite-degree closure at r=-4 (bounded degree only).  The committed
#     verify_quantum_exotic_cokernel.py certifies W2 infeasibility at d=1,2:
#     there the positive cascade forces R in D*F[E] (slope=R(1)=0 != 1).  We
#     reproduce the d=1 forcing compactly and self-containedly.
def falling_(n):
    return sp.prod(E - j for j in range(n))


def clean_solve_min(A_, B_, m, lkey, memb, name, raw_deg):
    cs = list(sp.symbols(f"{name}_0:{raw_deg+1}"))
    unknown = sp.expand(falling_(memb) * sum(cs[j] * E**j for j in range(raw_deg + 1)))
    trial = dict(B_); trial[lkey] = unknown
    eqs = sp.Poly(q_m(A_, trial, m), E).all_coeffs()
    M, rhs = sp.linear_eq_to_matrix(eqs, cs)
    conds = [sp.expand(nv.dot(rhs)) for nv in M.T.nullspace()]
    conds = [cc for cc in conds if cc != 0]
    indep = sp.zeros(0, len(cs)); selrhs = []
    for i in range(M.rows):
        cand = indep.col_join(M[i, :])
        if cand.rank() > indep.rank():
            indep = cand; selrhs.append(rhs[i])
    if indep.rows == 0:
        vals = [sp.Integer(0)] * len(cs)
    else:
        solu, _ = indep.gauss_jordan_solve(sp.Matrix(selrhs))
        vals = [x.subs({t: 0 for t in solu.free_symbols}) for x in solu]
    res = sp.expand(unknown.subs(dict(zip(cs, vals))))
    return res, conds


d = 1
a2f = poly("z2", d); a1f = poly("z1", d); a0f = poly("z0", d)
Af = {3: a3_w2, 2: a2f, 1: a1f, 0: a0f,
      -1: falling_(1) * poly("zm1", d), -2: falling_(2) * poly("zm2", d),
      -3: falling_(3) * poly("zm3", d)}
Bf = {k: sp.Integer(0) for k in LEVELS}
Bf[2] = b2_w2
Bf[-3] = sp.expand(sp.symbols("mu3") * Af[-3])
condsf = []
for m, lkey, name, memb, rawd in [(4, 1, "wb1", 0, d + 3), (3, 0, "wb0", 0, 2 * d + 2),
                                   (2, -1, "wbm1", 1, 2 * d + 3), (1, -2, "wbm2", 2, 2 * d + 4)]:
    Bf[lkey], cc = clean_solve_min(Af, Bf, m, lkey, memb, name, rawd)
    condsf += cc
Af0, Bf0 = dict(Af), dict(Bf); Af0[-2] = 0; Bf0[-3] = 0
Rf = potential(Af0, Bf0)
solve_vars = list(sp.symbols("z2_0 z2_1 z1_0 z1_1 z0_0 z0_1 zm1_0 zm1_1"))
branches = sp.solve(condsf, solve_vars, dict=True)
check(bool(branches), "W2 d=1: positive cascade has solved branch(es)")
forced = all(sp.expand(Rf.subs(sol).subs(E, 1)) == 0 for sol in branches)
check(forced, "W2 d=1: positive cascade FORCES slope R(1)=0 != 1  => infeasible at d=1")
print("   d=2 infeasibility is certified in the committed verify_quantum_exotic_cokernel.py.")
print("   Arbitrary-degree closure at r=-4 remains OPEN (a moment-unit / DC1 question).")


print("\nScope: closes the step-2 AP exotic family a_3=(E-r)(E-r-2)(E-r-4), b_2=(E-r-1)(E-r-4)")
print("at ARBITRARY coefficient degree for every r != -4 (unique exceptional locus r=-4 = W2,")
print("where central completion reduces to the moment unit).  W2 is exactly excluded by the")
print("combined slope+tail certificate at raw cap d=3, but remains open at arbitrary degree.")
print("Non-AP deg>=6 exotic tops are out")
print("of scope.  No Weyl pair or counterexample is constructed.")
print("ALL QUANTUM AP LAMBDA CHECKS PASSED")
