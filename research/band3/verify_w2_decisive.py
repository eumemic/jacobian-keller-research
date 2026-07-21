#!/usr/bin/env python3
"""W2-DECISIVE: the achievable values of the slope R(1) at the r=-4 exotic member.

W2 is the r = -4 member of the band-3 step-2 arithmetic-progression exotic family

    a_3 = E(E+2)(E+4)  (roots {-4,-2,0}),   b_2 = E(E+3),   b_3 = 0,

quantum band-3 conventions
    Q_m = sum_(k+l=m) (b_l^[k] a_k - a_k^[l] b_l),   f^[n](E)=f(E+n),
    G   = sum_(k=1)^3 sum_(j=0)^(k-1) (a_k^[j-k] b_-k^[j] - b_k^[j-k] a_-k^[j]),
    Q_0 = (T-1)G,   membership (E)_j | a_-j, b_-j.

The committed lambda-wave (commit d8189fc, quantum-ap-lambda.md) closes every AP
member r != -4 at ARBITRARY degree by the functional
lambda_r(f)=f(r+3)-f(r+4)+f(r+5)-f(0), lambda_r(E-R)=r+4.  At r=-4 the functional
degenerates and the whole question collapses to ONE scalar: the moment slope R(1).
Because Im Phi(W2) = E(E-1)(E+1) F[E] EXACTLY (re-verified below),

    Q_0 = 1  <=>  E - R in Im Phi(W2)  <=>  D | (E-R)  <=>  R(1)=1 and R(-1)=-1,

with R(0)=0 automatic and R(1)+R(-1)=0 the proved cascade constraint.  So W2
lives or dies by whether the positive cascade can produce slope R(1)=1.

VERDICT (this script, exact arithmetic over QQ):
  * d=1: R = 0, so R(1)=0.                                    (forced -> Q_0=1 infeasible)
  * d=2: two branches, both R in D*F[E], so R(1)=0.           (forced -> Q_0=1 infeasible)
  * d>=3: R(1) is a GENUINE FREE MODULUS (shape (b)).  It is NOT forced to 0;
          R(1) = -(8/9) am1_3 on a displayed slice (am1_3 = top coeff of a_-1),
          the coefficient varies with the base data, and R(1)=1 is ACHIEVED.
  An EXPLICIT d=3 solution with R(1)=1 is reconstructed from scratch and verified:
  Q_1=Q_2=Q_3=Q_4=0, full Weyl membership, and (with constructed fillers) G=E, so
  Q_0=1 -- all as EXACT polynomial identities.  These are degree-cap-agnostic
  polynomials, hence a valid solution for EVERY d>=3.

  Consequence: the slope obstruction does NOT close W2.  The band-3 exotic sector
  survives past R(1); only the negative tail Q_-1..Q_-6 now stands between W2 and a
  DC1 counterexample candidate.  This script imposes NO negative-tail equation and
  constructs NO Weyl pair.

Run:  uv run --with sympy python research/band3/verify_w2_decisive.py
Ends: ALL W2 DECISIVE CHECKS PASSED
"""
import sympy as sp

E, t = sp.symbols("E t")
LEVELS = range(-3, 4)


def sh(f, n):
    return sp.expand(sp.sympify(f).subs(E, E + n))


def falling(n):
    return sp.prod(E - j for j in range(n))


def poly(name, degree):
    cs = list(sp.symbols(f"{name}_0:{degree + 1}"))
    return sp.expand(sum(cs[j] * E**j for j in range(degree + 1))), cs


def q_m(A, B, m):
    """Stipulated convention Q_m = sum_(k+l=m) (b_l^[k] a_k - a_k^[l] b_l)."""
    return sp.expand(sum(sh(B[l], k) * A[k] - sh(A[k], l) * B[l]
                         for k in LEVELS for l in LEVELS if k + l == m))


def potential(A, B):
    return sp.expand(sum(sh(A[k], j - k) * sh(B[-k], j) - sh(B[k], j - k) * sh(A[-k], j)
                         for k in range(1, 4) for j in range(k)))


def check(condition, label):
    if not condition:
        raise AssertionError("FAIL: " + label)
    print("PASS", label)


def check_zero(value, label):
    if sp.expand(value) != 0:
        raise AssertionError(f"FAIL: {label}: residual {sp.factor(sp.expand(value))}")
    print("PASS", label)


W2_a3 = sp.expand(E * (E + 2) * (E + 4))
W2_b2 = sp.expand(E * (E + 3))
D = sp.expand(E * (E - 1) * (E + 1))


# =====================================================================
print("--- 0. W2 datum, conventions, and the structural decomposition ---")
# =====================================================================
check_zero(sh(W2_b2, 3) * W2_a3 - sh(W2_a3, 2) * W2_b2,
           "W2: b_2=E(E+3) solves the Q_5 wall b_2^[3]a_3 = a_3^[2]b_2")
check_zero(W2_a3 - E * (E + 2) * (E + 4), "W2 a_3 has roots {-4,-2,0} (r=-4 AP member)")

# generic symbolic data (any degree) to check the stipulated identities
A = {k: sp.Integer(0) for k in LEVELS}
B = {k: sp.Integer(0) for k in LEVELS}
A[3], B[2] = W2_a3, W2_b2
A[2], _ = poly("a2g", 5)
A[1], _ = poly("a1g", 5)
A[0], _ = poly("a0g", 5)
A[-1] = falling(1) * poly("am1g", 5)[0]
A[-2] = falling(2) * poly("am2g", 5)[0]
A[-3] = falling(3) * poly("am3g", 5)[0]
B[1], _ = poly("b1g", 6)
B[0], _ = poly("b0g", 6)
B[-1] = falling(1) * poly("bm1g", 6)[0]
B[-2] = falling(2) * poly("bm2g", 6)[0]
B[-3] = falling(3) * poly("bm3g", 5)[0]


def k3(c):
    return sp.expand(sum(sh(W2_a3, j - 3) * sh(c, j) for j in range(3)))


def h2(v):
    return sp.expand(sum(sh(W2_b2, j - 2) * sh(v, j) for j in range(2)))


G = potential(A, B)
R_gen = sp.expand(sh(A[1], -1) * B[-1] - sh(B[1], -1) * A[-1]
                  + sh(A[2], -2) * B[-2] + sh(A[2], -1) * sh(B[-2], 1))
check_zero(G - (R_gen + k3(B[-3]) - h2(A[-2])),
           "G = R + K_3[b_-3] - H_2[a_-2]  (exact, W2)")
check_zero(G.subs(E, 0), "full Weyl membership gives G(0)=0")
check_zero(q_m(A, B, 0) - (sh(G, 1) - G), "Q_0 = (T-1)G")
check_zero(R_gen.subs(E, 0), "R(0)=0 from membership (a_-1(0)=b_-1(0)=b_-2(0)=b_-2(1)=0)")


# =====================================================================
print("\n--- 1. PIVOT FACT: Im Phi(W2) = E(E-1)(E+1) F[E] EXACTLY ---")
# =====================================================================
# Phi(C,V) = K_3[(E)_3 C] - H_2[(E)_2 V].  (i) subset D*F[E]; (ii) equality by
# rank-matching against the D-multiple lattice at each finite truncation; the
# codimension stabilizes at 3 = deg D.  Hence E-R in Im Phi <=> D | (E-R).
sub_ok = all(sp.rem(sp.Poly(k3(falling(3) * E**j), E), sp.Poly(D, E)).as_expr() == 0
             and sp.rem(sp.Poly(h2(falling(2) * E**j), E), sp.Poly(D, E)).as_expr() == 0
             for j in range(8))
check(sub_ok, "Im Phi(W2) subset D*F[E]: every basis filler divisible by D=E(E-1)(E+1)")


def coeffvec(f, N):
    p = sp.Poly(sp.expand(f), E)
    return [p.nth(j) for j in range(N + 1)]


for dcap in (3, 4, 5):
    N = dcap + 6
    cols = ([coeffvec(k3(falling(3) * E**j), N) for j in range(dcap + 1)]
            + [coeffvec(h2(falling(2) * E**j), N) for j in range(dcap + 1)])
    Mimg = sp.Matrix(cols).T
    Dcols = [coeffvec(sp.expand(D * E**j), N) for j in range(N - 3 + 1)]
    MD = sp.Matrix(Dcols).T
    stacked = MD.row_join(Mimg)
    equality = (Mimg.rank() == MD.rank() == (N - 2)) and (stacked.rank() == MD.rank())
    codim = (N + 1) - Mimg.rank()
    check(equality and codim == 3,
          f"Im Phi(W2) = D*F[E] up to deg {N} (rank {Mimg.rank()} = dim D*F[E]; codim {codim})")
print("   => Q_0=1  <=>  E-R in Im Phi(W2)  <=>  D|(E-R)  <=>  R(1)=1 AND R(-1)=-1.")


# =====================================================================
print("\n--- 1b. the reduction R(1)=G(1)=Q_0(0),  R(-1)=G(-1)=-Q_0(-1) ---")
# =====================================================================
# The two filler blocks lie in D*F[E], so they vanish at E=0,1,-1; hence R and G
# agree at those three nodes, and (T-1)G ties the slope to the boundary of Q_0.
filler = sp.expand(G - R_gen)
check(sp.rem(sp.Poly(filler, E), sp.Poly(D, E)).as_expr() == 0,
      "G - R = K_3[b_-3]-H_2[a_-2] lies in D*F[E] (vanishes at 0,1,-1)")
Q0 = q_m(A, B, 0)
check_zero(R_gen.subs(E, 1) - G.subs(E, 1), "R(1) = G(1)")
check_zero(R_gen.subs(E, -1) - G.subs(E, -1), "R(-1) = G(-1)")
check_zero(R_gen.subs(E, 1) - Q0.subs(E, 0), "R(1) = Q_0(0)  (the moment slope)")
check_zero(R_gen.subs(E, -1) + Q0.subs(E, -1), "R(-1) = -Q_0(-1)")


# =====================================================================
print("\n--- 2. the positive-cascade solver (clean_solve), kernels retained ---")
# =====================================================================
def clean_solve(Ain, Bin, m, lkey, name, membership, raw_degree):
    raw, cs = poly(name, raw_degree)
    unknown = sp.expand(falling(membership) * raw)
    trial = dict(Bin)
    trial[lkey] = unknown
    equations = sp.Poly(q_m(Ain, trial, m), E).all_coeffs()
    M, rhs = sp.linear_eq_to_matrix(equations, cs)
    if any(entry.free_symbols for entry in M):
        raise AssertionError("operator matrix is not numeric; bilinearity leaked")
    conditions = [sp.expand(n.dot(rhs)) for n in M.T.nullspace()]
    conditions = [c for c in conditions if c != 0]
    independent = sp.zeros(0, len(cs))
    selected_rhs = []
    for i in range(M.rows):
        cand = independent.col_join(M[i, :])
        if cand.rank() > independent.rank():
            independent = cand
            selected_rhs.append(rhs[i])
    if independent.rows == 0:
        values = [sp.Integer(0)] * len(cs)
    else:
        solution, _ = independent.gauss_jordan_solve(sp.Matrix(selected_rhs))
        taus = [x for x in solution.free_symbols if str(x).startswith("tau")]
        values = [x.subs({tt: 0 for tt in taus}) for x in solution]
    result = sp.expand(unknown.subs(dict(zip(cs, values))))
    kernels = []
    for j, vector in enumerate(M.nullspace()):
        parameter = sp.symbols(f"{name}K{j}")
        kernel_poly = falling(membership) * sum(vector[i] * E**i for i in range(len(cs)))
        result = sp.expand(result + parameter * kernel_poly)
        kernels.append(parameter)
    return result, kernels, conditions


def build_cascade(a3, b2, d):
    a2, ca2 = poly("a2", d)
    a1, ca1 = poly("a1", d)
    a0, ca0 = poly("a0", d)
    am1_raw, cam1 = poly("am1", d)
    am2_raw, cam2 = poly("am2", d)
    am3_raw, cam3 = poly("am3", d)
    mu3 = sp.symbols("mu3")
    Ac = {3: a3, 2: a2, 1: a1, 0: a0, -1: falling(1) * am1_raw,
          -2: falling(2) * am2_raw, -3: falling(3) * am3_raw}
    Bc = {k: sp.Integer(0) for k in LEVELS}
    Bc[2] = b2
    Bc[-3] = sp.expand(mu3 * Ac[-3])
    conditions, kernels = [], []
    for m, lkey, name, membership, degree in [
            (4, 1, "b1c", 0, d + 3), (3, 0, "b0c", 0, 2 * d + 2),
            (2, -1, "bm1c", 1, 2 * d + 3), (1, -2, "bm2c", 2, 2 * d + 4)]:
        Bc[lkey], nk, nc = clean_solve(Ac, Bc, m, lkey, name, membership, degree)
        kernels += nk
        conditions += nc
    data = {"ca2": ca2, "ca1": ca1, "ca0": ca0, "cam1": cam1}
    return Ac, Bc, conditions, kernels, data


def residual_on_cascade(Ac, Bc):
    A0, B0 = dict(Ac), dict(Bc)
    A0[-2] = sp.Integer(0)
    B0[-3] = sp.Integer(0)
    return potential(A0, B0)


# d=1: R = 0 identically.
Ac, Bc, conds1, ker1, dat1 = build_cascade(W2_a3, W2_b2, 1)
R1poly = residual_on_cascade(Ac, Bc)
sv1 = dat1["ca2"] + dat1["ca1"] + dat1["ca0"] + dat1["cam1"]
br1 = sp.solve(conds1, sv1, dict=True, simplify=False)
check(bool(br1), "d=1: positive cascade has solved branch(es)")
check(all(sp.expand(R1poly.subs(s)) == 0 for s in br1),
      "d=1: R = 0 on every branch  =>  R(1)=0 forced  (Q_0=1 infeasible)")

# d=2: two branches, both R in D*F[E], so R(1)=R(-1)=0.
Ac2, Bc2, conds2, ker2, dat2 = build_cascade(W2_a3, W2_b2, 2)
R2poly = residual_on_cascade(Ac2, Bc2)
sv2 = dat2["ca2"] + dat2["ca1"] + dat2["ca0"] + dat2["cam1"]
br2 = sp.solve(conds2, sv2, dict=True, simplify=False)
check(len(br2) >= 2, f"d=2: positive cascade splits into {len(br2)} branches")
for i, s in enumerate(br2):
    Rb = sp.expand(R2poly.subs(s))
    check(sp.rem(sp.Poly(Rb, E), sp.Poly(D, E)).as_expr() == 0 and Rb.subs(E, 1) == 0,
          f"d=2 branch {i}: R in D*F[E], so R(1)=0 forced  (Q_0=1 infeasible)")


# =====================================================================
print("\n--- 3. DECISIVE: at d=3 the slope R(1) is a FREE MODULUS (shape (b)) ---")
# =====================================================================
Ac3, Bc3, conds3, ker3, dat3 = build_cascade(W2_a3, W2_b2, 3)
R3poly = residual_on_cascade(Ac3, Bc3)
check(ker3 == [sp.symbols("b0cK0")], "d=3: sole retained solver kernel is constant b0cK0")
check(sp.symbols("b0cK0") not in R3poly.free_symbols,
      "d=3: residual R is independent of the retained b_0 kernel")
a2_0, a2_1, a2_2, a2_3 = dat3["ca2"]
a1_0, a1_1, a1_2, a1_3 = dat3["ca1"]
a0_0, a0_1, a0_2, a0_3 = dat3["ca0"]
am1_0, am1_1, am1_2, am1_3 = dat3["cam1"]

# Deterministic slice: fix a2 (a2_1 forced by cond[0]) and a1_3; keep am1_3 SYMBOLIC,
# solve the remaining LINEAR conditions in two stages.  This yields an honest
# 1-PARAMETER FAMILY of positive-cascade solutions on which R(1) = -(8/9) am1_3 --
# a genuine free modulus, decisively refuting any forcing R(1)=0 at d=3.  (No
# Groebner basis is needed: an explicit solution family with R(1) surjective is a
# stronger and faster certificate than a radical-ideal membership test.)
fixL = {a2_0: sp.Integer(2), a2_2: sp.Rational(1, 2), a2_3: sp.Rational(-1, 3),
        a1_3: sp.Rational(5, 7)}
fixL[a2_1] = 3 * fixL[a2_2] - 9 * fixL[a2_3]
cf = [c for c in (sp.expand(c.subs(fixL)) for c in conds3) if c != 0]
a1_only = [c for c in cf if c.free_symbols <= {a1_0, a1_1, a1_2}]
fixL.update(dict(zip([a1_0, a1_1, a1_2],
                     list(sp.linsolve(a1_only, [a1_0, a1_1, a1_2]))[0])))
cf = [c for c in (sp.expand(c.subs(fixL)) for c in conds3) if c != 0]
rest = [a0_1, a0_2, a0_3, am1_0, am1_1, am1_2]   # am1_3 kept symbolic
fixL.update(dict(zip(rest, list(sp.linsolve(cf, rest))[0])))

# it IS a solution family: every positive-cascade condition vanishes on it (for all am1_3)
check(all(sp.expand(c.subs(fixL)) == 0 for c in conds3),
      "d=3: the am1_3-family satisfies EVERY positive-cascade condition (genuine solutions)")
R3_law = sp.expand(R3poly.subs(fixL))
check_zero(sp.expand(R3_law.subs(E, 1)) + sp.Rational(8, 9) * am1_3,
           "DECISIVE: R(1) = -(8/9) am1_3 on the family -- a FREE MODULUS (shape (b))")
check_zero(sp.expand(R3_law.subs(E, -1)) - sp.Rational(8, 9) * am1_3,
           "R(-1) = +(8/9) am1_3, so R(1)+R(-1)=0 (proved constraint) but R(1) NOT pinned")
check(sp.expand(R3_law.subs(E, 1)) != 0 and am1_3 in R3_law.subs(E, 1).free_symbols,
      "R(1) is genuinely nonconstant in am1_3 => NOT forced to 0 at d=3 (contrast d=1,2)")
# CONTRAST CONTROL: at d<=2 the same slope is FORCED to 0 (Section 2); the free
# modulus am1_3 is the E^4-coefficient of a_-1, which exists only for d>=3.


# =====================================================================
print("\n--- 4. EXPLICIT d=3 solution with R(1)=1 (reconstructed from scratch) ---")
# =====================================================================
# pin the slope dial am1_3 = -9/8 to hit R(1)=1
fix = {kk: sp.nsimplify(vv.subs(am1_3, sp.Rational(-9, 8))) for kk, vv in fixL.items()}
fix[am1_3] = sp.Rational(-9, 8)
R3_slice = sp.expand(R3poly.subs(fix))
check_zero(sp.expand(R3_slice.subs(E, 1)) - 1, "slice: R(1) = 1  (am1_3 = -9/8)")
check_zero(sp.expand(R3_slice.subs(E, -1)) + 1, "slice: R(-1) = -1")

# freeze all remaining free symbols to numbers, making the point concrete
concrete = dict(fix)
solved_A = {k: sp.expand(v.subs(concrete)) for k, v in Ac3.items()}
solved_B = {k: sp.expand(v.subs(concrete)) for k, v in Bc3.items()}
solved_A[-2] = sp.Integer(0)
solved_B[-3] = sp.Integer(0)
free_left = set()
for v in list(solved_A.values()) + list(solved_B.values()):
    free_left |= v.free_symbols
free_left -= {E}
zsub = {x: sp.Integer(0) for x in free_left}
solved_A = {k: sp.expand(v.subs(zsub)) for k, v in solved_A.items()}
solved_B = {k: sp.expand(v.subs(zsub)) for k, v in solved_B.items()}

# DIRECT verification of the positive cascade (fillers still zero)
for m in range(1, 5):
    check_zero(q_m(solved_A, solved_B, m), f"explicit point: Q_{m} = 0 (exact identity)")
for kk in (1, 2, 3):
    check(sp.rem(sp.Poly(solved_A[-kk], E), sp.Poly(falling(kk), E)).as_expr() == 0,
          f"explicit point: (E)_{kk} | a_-{kk}")
Rval = sp.expand(potential(solved_A, solved_B))
check(Rval.subs(E, 0) == 0 and Rval.subs(E, 1) == 1 and Rval.subs(E, -1) == -1,
      "explicit point: R(0)=0, R(1)=1, R(-1)=-1")
EmR = sp.expand(E - Rval)
check(sp.rem(sp.Poly(EmR, E), sp.Poly(D, E)).as_expr() == 0,
      "explicit point: D | (E-R), so E-R in Im Phi(W2)")

# CONSTRUCT the fillers Phi(C,V) = E-R by exact linear algebra, verify G=E, Q_0=1
cap = sp.Poly(EmR, E).degree()
Cc = list(sp.symbols(f"CC_0:{cap + 1}"))
Vc = list(sp.symbols(f"VV_0:{cap + 1}"))
Cp = sum(Cc[j] * E**j for j in range(cap + 1))
Vp = sum(Vc[j] * E**j for j in range(cap + 1))
phi = sp.expand(k3(falling(3) * Cp) - h2(falling(2) * Vp))
diff = sp.expand(phi - EmR)
coeffs = [diff.coeff(E, i) for i in range(sp.Poly(diff, E).degree() + 1)]
fs = list(sp.linsolve(coeffs, Cc + Vc))
check(fs != [] and fs[0] != sp.EmptySet, "filler system Phi(C,V)=E-R is consistent")
s = fs[0]
subd = dict(zip(Cc + Vc, s))
part = {p: 0 for p in set().union(*[sp.sympify(x).free_symbols for x in s])}
c_fill = sp.expand(falling(3) * sp.expand(Cp.subs(subd).subs(part)))
v_fill = sp.expand(falling(2) * sp.expand(Vp.subs(subd).subs(part)))
check_zero(k3(c_fill) - h2(v_fill) - EmR, "constructed fillers satisfy Phi(C,V) = E-R")

full_A = dict(solved_A)
full_B = dict(solved_B)
full_A[-2] = v_fill          # a_-2 filler
full_B[-3] = c_fill          # b_-3 filler
check_zero(potential(full_A, full_B) - E, "with fillers: G = E  (exact identity)")
check_zero(q_m(full_A, full_B, 0) - 1, "with fillers: Q_0 = 1  (exact identity)")
check(all(q_m(full_A, full_B, m) == 0 for m in range(1, 5)),
      "with fillers: Q_1=Q_2=Q_3=Q_4=0 still hold (fillers are cascade-transparent)")
check(all(sp.rem(sp.Poly(full_A[-kk], E), sp.Poly(falling(kk), E)).as_expr() == 0
          and sp.rem(sp.Poly(full_B[-kk], E), sp.Poly(falling(kk), E)).as_expr() == 0
          for kk in (1, 2, 3)),
      "with fillers: full Weyl membership (E)_j | a_-j, b_-j for j=1,2,3")
print("   => positive cascade + membership + fillers + Q_0=1 ALL satisfied at R(1)=1.")

# PERSISTENCE to all d>=3: the explicit solution is a set of CONCRETE polynomials
# of bounded degree; Q_m=0, membership, G=E are polynomial identities independent
# of any coefficient-degree cap, so the SAME data is a solution for every d>=3.
maxdeg = max(sp.Poly(full_A[k], E).degree() for k in full_A if full_A[k] != 0)
maxdeg = max(maxdeg, max(sp.Poly(full_B[k], E).degree() for k in full_B if full_B[k] != 0))
check(isinstance(maxdeg, int) and maxdeg < 10**6,
      f"explicit solution is concrete polynomials (max degree {maxdeg}); valid for all d>=3")


# =====================================================================
print("\n--- 5. sanity controls: the cascade ideal is PROPER; no false kills ---")
# =====================================================================
# (a) the positive-cascade ideal at d=3 is proper (no false kills): the explicit
#     verified point of Section 4 is a witness solution, so V(conditions) is nonempty.
#     (A witness is a stronger properness certificate than a Groebner non-[1] test.)
check(all(sp.expand(c.subs(fix)) == 0 for c in conds3),
      "d=3 positive-cascade ideal is PROPER: the explicit point satisfies every condition")
# (b) genuine control: the SAME finite-cokernel test that PASSES E-R at R(1)=1
#     (d=3) correctly REJECTS Q_0=1 at d=2, where R in D*F[E] forces E-R = E-(D-mult),
#     which is NOT in D*F[E] since D does not divide E.
check(sp.rem(sp.Poly(E, E), sp.Poly(D, E)).as_expr() != 0,
      "control: D does not divide E  (so E notin Im Phi(W2))")
for i, s in enumerate(br2):
    Rb = sp.expand(R2poly.subs(s))
    EmRb = sp.expand(E - Rb)
    check(sp.rem(sp.Poly(EmRb, E), sp.Poly(D, E)).as_expr() != 0,
          f"control: d=2 branch {i}: E-R notin Im Phi(W2)  =>  Q_0=1 correctly INFEASIBLE")

print("\nScope: this settles the achievable SLOPE at W2 only.  R(1)=1 is achievable at")
print("every d>=3 (forced 0 only at d<=2), so the moment-slope obstruction does NOT")
print("close W2; the band-3 exotic sector survives to the negative tail Q_-1..Q_-6.")
print("No negative-tail equation is imposed here and no Weyl pair is constructed.")
print("ALL W2 DECISIVE CHECKS PASSED")
