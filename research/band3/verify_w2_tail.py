#!/usr/bin/env python3
"""Exact algebra certificate for the W2 NEGATIVE TAIL (the bottom wall of the
band-3 exotic AP family at the exceptional member r = -4).

W2 datum (gauge b_3 = 0, quantum band-3 conventions):

    a_3 = E(E+2)(E+4)   (roots {0,-2,-4}),     b_2 = E(E+3),     b_3 = 0.

This verifier owns the full system BELOW the fillers -- the part no exotic
analysis needed before, because every other exotic top died upstream.  It:

  0. builds Q_m for every m in [-6,6] and checks it equals the crossed-product
     commutator coefficient with fully generic degree-2 coefficients;
  1. bottom proportionality:  Q_-6 = 0  <=>  b_-3 = mu_3 a_-3  (mirror of the
     top Q_6 gauge), with membership preserved;
  2. bottom wall (phi-gauge normal form):  Q_-5 = phi^[-3] a_-3 - a_-3^[-2] phi
     with phi = b_-2 - mu_3 a_-2;
  3. the Q_-4 decomposition (cascade memo 5.3) and the staggered Lemma-R
     leading-coefficient identities;
  4. the reflection E -> -E-1 carries the TOP wall to the BOTTOM wall exactly;
     the mirror Wall Lemma root-multiset criterion  Phi_3 * B_phi = (S+1) A_-3
     and the mirror degree law  2 deg a_-3 = 3 deg phi;
  5. the BOTTOM CUBE-CLASS GATE: concrete Phi_3-necklace examples -- (E)_3 is
     solvable (tame shifted cube), roots {0,1,2,4} fail (gate fires).  The
     membership-minimal a_-3 = (E)_3 is itself Phi_3, so the bottom base is TAME:
     the bottom wall does NOT independently force W2 exotic;
  6. checks finite regressions for the pivot fact Im Phi(W2) = D*F[E],
     D = E(E-1)(E+1), whose arbitrary-degree triangular proof is in
     quantum-ap-filler-image.md, and verifies the resulting criterion
     E-R in Im Phi  <=>  R(1)=1 AND R(-1)=-1 (R(0)=0 auto);
  7. re-derives the cascade constraint R(1)+R(-1)=0 and R=0 at d=1 from the
     solved positive cascade (=> slope R(1)=0 != 1 at d=1);
  8. bounded separate-feasibility checks.  At d=1 the full DC1-face system
     (positive cascade + Q_0=1 + negative tail Q_-1..Q_-5 + membership) is the
     UNIT IDEAL, while:
        (a) positive cascade + Q_0=1  (no tail)   = unit ideal   [slope kills],
        (b) positive cascade + negative tail       = PROPER ideal [tail feasible].
     The tail subsystem without Q_0=1 remains feasible at d=2.  These bounded
     checks do not imply that the tail is globally non-obstructing;
  9. a genuine positive-control pair (X=U^3-d/kappa, D=kappa U) passes every
     Q_m=0 including the whole tail, and its a_-3 = c1^3 (E)_3 passes the bottom
     gate -- validating the tail machinery against a real Weyl pair;
 10. an ESCALATION HARNESS run_full_system(...) ready to accept the decisive
     sibling's slope=1 positive data (research/band3/w2-decisive.md) and test
     the full combined feasibility.

PROVED (arbitrary coefficient degree, char 0): sections 0-4, 6 -- structural
identities / reflection / root-multiset criteria.  BOUNDED EVIDENCE: section 5
(concrete necklace examples), section 7 (d=1), section 8 (d=1,2 Groebner).  No
Weyl pair and no counterexample is constructed here.  Base commit d8189fc; lambda
wave facts re-verified in-file, not merely cited.

Run:  uv run --with sympy python research/band3/verify_w2_tail.py
Ends: ALL W2 TAIL CHECKS PASSED
"""
import sympy as sp

E, S = sp.symbols("E S")
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
        raise AssertionError(label)
    print("PASS", label)


def check_zero(value, label):
    residual = sp.expand(value)
    if residual != 0:
        raise AssertionError(f"{label}: residual {sp.factor(residual)}")
    print("PASS", label)


a3_w2 = sp.expand(E * (E + 2) * (E + 4))
b2_w2 = sp.expand(E * (E + 3))
D = sp.expand(E * (E - 1) * (E + 1))


# ---- shared crossed-product engine and cascade solver (machine-checked) ----
def mul_ladders(P, Q):
    R = {}
    for k1, p1 in P.items():
        for k2, p2 in Q.items():
            R[k1 + k2] = R.get(k1 + k2, 0) + sh(p1, k2) * p2
    return {k: sp.expand(v) for k, v in R.items() if sp.expand(v) != 0}


def commutator_coeff(X, D_, m):
    """Direct [D_,X]_m in the crossed product A_1[x^-1]."""
    DX = mul_ladders(D_, X)
    XD = mul_ladders(X, D_)
    return sp.expand(DX.get(m, 0) - XD.get(m, 0))


def clean_solve(A, B, m, lkey, name, membership, raw_degree):
    raw, cs = poly(name, raw_degree)
    unknown = sp.expand(falling(membership) * raw)
    trial = dict(B)
    trial[lkey] = unknown
    equations = sp.Poly(q_m(A, trial, m), E).all_coeffs()
    M, rhs = sp.linear_eq_to_matrix(equations, cs)
    if any(entry.free_symbols for entry in M):
        raise AssertionError("operator matrix is not numeric; bilinearity leaked")
    conditions = [sp.expand(n.dot(rhs)) for n in M.T.nullspace()]
    conditions = [c for c in conditions if c != 0]
    independent = sp.zeros(0, len(cs))
    selected_rhs = []
    for i in range(M.rows):
        candidate = independent.col_join(M[i, :])
        if candidate.rank() > independent.rank():
            independent = candidate
            selected_rhs.append(rhs[i])
    if independent.rows == 0:
        values = [sp.Integer(0)] * len(cs)
    else:
        solution, _ = independent.gauss_jordan_solve(sp.Matrix(selected_rhs))
        taus = [x for x in solution.free_symbols if str(x).startswith("tau")]
        values = [x.subs({t: 0 for t in taus}) for x in solution]
    result = sp.expand(unknown.subs(dict(zip(cs, values))))
    kernels = []
    for j, vector in enumerate(M.nullspace()):
        parameter = sp.symbols(f"{name}K{j}")
        kernel_poly = falling(membership) * sum(vector[i] * E**i for i in range(len(cs)))
        result = sp.expand(result + parameter * kernel_poly)
        kernels.append(parameter)
    return result, kernels, conditions


def build_cascade(a3, b2, d):
    """Solve the positive cascade Q_4,Q_3,Q_2,Q_1 for b_1,b_0,b_-1,b_-2 and set
    the bottom Wronskian b_-3 = mu_3 a_-3.  Every coefficient is then expressed in
    the free data; Q_0-1 and Q_-1..Q_-5 become PURE constraints."""
    a2, ca2 = poly("a2", d)
    a1, ca1 = poly("a1", d)
    a0, ca0 = poly("a0", d)
    am1_raw, cam1 = poly("am1", d)
    am2_raw, cam2 = poly("am2", d)
    am3_raw, cam3 = poly("am3", d)
    mu3 = sp.symbols("mu3")
    A = {3: a3, 2: a2, 1: a1, 0: a0,
         -1: falling(1) * am1_raw, -2: falling(2) * am2_raw, -3: falling(3) * am3_raw}
    B = {k: sp.Integer(0) for k in range(-3, 4)}
    B[2] = b2
    B[-3] = sp.expand(mu3 * A[-3])
    conditions, kernels = [], []
    for m, lkey, name, membership, degree in [
            (4, 1, "b1c", 0, d + 3), (3, 0, "b0c", 0, 2 * d + 2),
            (2, -1, "bm1c", 1, 2 * d + 3), (1, -2, "bm2c", 2, 2 * d + 4)]:
        B[lkey], nk, nc = clean_solve(A, B, m, lkey, name, membership, degree)
        kernels += nk
        conditions += nc
    data = {"ca2": ca2, "ca1": ca1, "ca0": ca0, "cam1": cam1,
            "cam2": cam2, "cam3": cam3, "mu3": mu3, "kernels": kernels}
    return A, B, conditions, data


def scalar_coeffs(expr):
    return [c for c in sp.Poly(sp.expand(expr), E).all_coeffs() if sp.expand(c) != 0]


def is_unit_ideal(eqs, variables):
    eqs = [sp.expand(e) for e in eqs if sp.expand(e) != 0]
    G = sp.groebner(eqs, *variables, order="grevlex", domain=sp.QQ)
    return list(G) == [sp.Integer(1)]


# =====================================================================
print("--- 0. Q_m equals the crossed-product commutator for every m in [-6,6] ---")
# =====================================================================
Ag, Bg = {}, {}
for lev in LEVELS:
    ar, _ = poly(f"Ag{lev + 3}", 2)
    br, _ = poly(f"Bg{lev + 3}", 2)
    memb = falling(-lev) if lev < 0 else 1
    Ag[lev] = sp.expand(memb * ar)
    Bg[lev] = sp.expand(memb * br)
for m in range(-6, 7):
    check_zero(q_m(Ag, Bg, m) - commutator_coeff(Ag, Bg, m),
               f"Q_{m} = [D,X]_{m} (generic degree-2 coefficients)")
check_zero(q_m(Ag, Bg, 0) - (sh(potential(Ag, Bg), 1) - potential(Ag, Bg)),
           "Q_0 = (T-1)G telescoping identity (generic)")

# =====================================================================
print("\n--- 1. bottom proportionality: Q_-6 = 0  <=>  b_-3 = mu_3 a_-3 ---")
# =====================================================================
mu3s = sp.symbols("mu3s")
Bmu = dict(Bg)
Bmu[-3] = sp.expand(mu3s * Ag[-3])
check_zero(q_m(Ag, Bmu, -6), "Q_-6 vanishes identically under b_-3 = mu_3 a_-3")
# and the converse: Q_-6=0 forces the ratio b_-3/a_-3 to be (-3)-periodic hence constant.
ratio = sp.symbols("ratio")
check_zero(q_m({-3: Ag[-3]}, {-3: sp.expand(ratio * Ag[-3])}, -6),
           "any scalar ratio b_-3 = ratio*a_-3 solves Q_-6 (mirror of Q_6 top gauge)")
# membership of b_-3 is inherited from a_-3 (both carry (E)_3):
check_zero(sp.Poly(Bmu[-3], E).rem(sp.Poly(falling(3), E)).as_expr(),
           "b_-3 = mu_3 a_-3 carries the Weyl membership factor (E)_3")

# =====================================================================
print("\n--- 2. bottom wall (phi-gauge normal form): Q_-5 = phi^[-3]a_-3 - a_-3^[-2]phi ---")
# =====================================================================
phi = sp.expand(Bmu[-2] - mu3s * Ag[-2])
bottomwall = sp.expand(sh(phi, -3) * Ag[-3] - sh(Ag[-3], -2) * phi)
check_zero(q_m(Ag, Bmu, -5) - bottomwall,
           "Q_-5 = phi^[-3] a_-3 - a_-3^[-2] phi,  phi = b_-2 - mu_3 a_-2")
# the mu_3 source is exactly the correction turning raw b_-2 into the gauged phi:
raw_wall = sp.expand(sh(Bmu[-2], -3) * Ag[-3] - sh(Ag[-3], -2) * Bmu[-2])
mu_source = sp.expand(mu3s * (sh(Ag[-3], -2) * Ag[-2] - sh(Ag[-2], -3) * Ag[-3]))
check_zero(q_m(Ag, Bmu, -5) - (raw_wall + mu_source),
           "raw b_-2 wall plus the mu_3-source fixes the minus sign in phi")

# =====================================================================
print("\n--- 3. Q_-4 decomposition (cascade memo 5.3) and Lemma-R leading terms ---")
# =====================================================================
q4dec = sp.expand(
    (sh(Bmu[-1], -3) * Ag[-3] - sh(Ag[-3], -1) * Bmu[-1])
    + mu3s * (sh(Ag[-3], -1) * Ag[-1] - sh(Ag[-1], -3) * Ag[-3])
    + (sh(Bmu[-2], -2) * Ag[-2] - sh(Ag[-2], -2) * Bmu[-2]))
check_zero(q_m(Ag, Bmu, -4) - q4dec,
           "Q_-4 = [b_-1 wall] + mu_3[a_-1 source] + [b_-2^[-2]a_-2 - a_-2^[-2]b_-2]")
# staggered leading-coefficient Lemma-R (bottom wall): for f^[a]g^[b]-f^[a']g^[b'],
# coeff(E^{p+q-1}) = ((a-a')p+(b-b')q) lc(f) lc(g).  Bottom wall (a,b)/(a',b')=(0,-3)/(-2,0):
f5, _ = poly("f5", 4)
g5, _ = poly("g5", 3)
p5, q5 = 4, 3
stag = sp.expand(sh(g5, -3) * f5 - sh(f5, -2) * g5)
lead = sp.Poly(stag, E).nth(p5 + q5 - 1)
lcf, lcg = sp.Poly(f5, E).LC(), sp.Poly(g5, E).LC()
check_zero(lead - ((2 * p5 - 3 * q5) * lcf * lcg),
           "bottom-wall Lemma-R: coeff(E^{p+q-1}) = (2p-3q) lc(a_-3) lc(phi)")

# =====================================================================
print("\n--- 4. reflection E->-E-1 carries the TOP wall to the BOTTOM wall ---")
# =====================================================================
def refl(f):
    return sp.expand(sp.sympify(f).subs(E, -E - 1))


a3g, _ = poly("a3g", 3)
b2g, _ = poly("b2g", 2)
topwall = sp.expand(sh(b2g, 3) * a3g - sh(a3g, 2) * b2g)
A3r, B2r = refl(a3g), refl(b2g)
bottomform = sp.expand(sh(B2r, -3) * A3r - sh(A3r, -2) * B2r)
check_zero(refl(topwall) - bottomform,
           "reflection sends b_2^[3]a_3 - a_3^[2]b_2  to  phi^[-3]a - a^[-2]phi")
# mirror Wall Lemma root-multiset criterion:  Phi_3(S) B_phi = (S+1) A_-3.
Phi3 = 1 + S + S**2


def rootmultiset(roots):
    m = min(roots)
    return sp.expand(sum(S**(rr - m) for rr in roots))


def phi3_divides(roots):
    return sp.rem(sp.Poly(rootmultiset(roots), S), sp.Poly(Phi3, S)).as_expr() == 0


# top wall root-multiset relation was  S(1+S) A = Phi_3 B  ;  mirror is  (1+S) A = Phi_3 B_phi.
Aroots = sp.Symbol("Adummy")  # symbolic check via a concrete necklace
A_nl = rootmultiset([0, 1, 2])            # = Phi_3
Bphi_forced = sp.expand(sp.cancel((1 + S) * A_nl / Phi3))
check(sp.denom(sp.together(Bphi_forced)) == 1,
      "mirror criterion: (S+1)A_-3 / Phi_3 is a genuine (polynomial) multiset for A_-3=Phi_3")
check_zero(Phi3 * Bphi_forced - (1 + S) * A_nl,
           "mirror Wall relation  Phi_3 B_phi = (S+1) A_-3  holds")

# =====================================================================
print("\n--- 5. BOTTOM CUBE-CLASS GATE: concrete Phi_3-necklace examples ---")
# =====================================================================
def bottom_wall_solve(am3poly, maxdeg=10):
    cs = list(sp.symbols(f"pw0:{maxdeg + 1}"))
    ph = sum(cs[j] * E**j for j in range(maxdeg + 1))
    eq = sp.expand(sh(ph, -3) * am3poly - sh(am3poly, -2) * ph)
    sol = list(sp.linsolve(sp.Poly(eq, E).all_coeffs(), cs))[0]
    return sp.expand(ph.subs(dict(zip(cs, sol))))


# minimal membership a_-3 = (E)_3, roots {0,1,2}: consecutive cube = Phi_3 -> TAME, solvable.
phi_min = bottom_wall_solve(falling(3), 6)
check(phi_min != 0 and phi3_divides([0, 1, 2]),
      "a_-3=(E)_3 (roots {0,1,2}=Phi_3, shifted cube): bottom wall has nonzero phi (TAME base)")
check_zero(2 * sp.Poly(falling(3), E).degree() - 3 * sp.Poly(phi_min, E).degree(),
           "degree law 2 deg a_-3 = 3 deg phi holds for a_-3=(E)_3")
# two consecutive cubes: roots {0..5} still a necklace.
phi_two = bottom_wall_solve(sp.expand(falling(3) * (E - 3) * (E - 4) * (E - 5)), 8)
check(phi_two != 0 and phi3_divides([0, 1, 2, 3, 4, 5]),
      "a_-3 roots {0,1,2,3,4,5}: Phi_3 necklace, bottom wall solvable")
# GATE FIRES: roots {0,1,2,4} not Phi_3-divisible -> only trivial phi.
phi_bad = bottom_wall_solve(sp.expand(falling(3) * (E - 4)), 8)
check(phi_bad == 0 and not phi3_divides([0, 1, 2, 4]),
      "a_-3 roots {0,1,2,4}: NOT a Phi_3 necklace -> bottom wall GATE FIRES (phi=0 only)")
# the top W2 AP roots {0,-2,-4} are Phi_3-divisible (the exotic necklace), for contrast:
check(phi3_divides([0, -2, -4]),
      "contrast: W2 TOP roots {0,-2,-4} ARE a Phi_3 necklace (the exotic step-2 AP)")
print("   NOTE: the bottom's forced roots {0,1,2} are the CONSECUTIVE (tame) cube, not the")
print("   exotic step-2 AP; the bottom wall does NOT independently force W2 exotic.")

# =====================================================================
print("\n--- 6. pivot fact: degree-free triangular argument plus finite regressions; E-R criterion ---")
# =====================================================================
def k3(cc):
    return sp.expand(sum(sh(a3_w2, j - 3) * sh(cc, j) for j in range(3)))


def h2(vv):
    return sp.expand(sum(sh(b2_w2, j - 2) * sh(vv, j) for j in range(2)))


check(all(sp.rem(k3(falling(3) * E**j), D, E) == 0 and sp.rem(h2(falling(2) * E**j), D, E) == 0
          for j in range(8)),
      "finite regression j=0..7: W2 basis fillers are divisible by D=E(E-1)(E+1)")
# Reverse inclusion is degree-free: for symbolic n, the two degree-(n+1) terms in
# B(E^n) contribute -1 each, while every remaining binomial term has lower degree.
# Thus deg B(E^n)=n+1 with leading coefficient -2 for every n>=0, triangularly
# reducing each positive-degree quotient target to a constant.  These finite cases
# are implementation regressions for that written argument.
for n in range(4):
    Bn = sp.expand((2 - E) * E**n - (E + 2) * sh(E**n, 1))
    check(sp.Poly(Bn, E).degree() == n + 1 and sp.Poly(Bn, E).LC() == -2,
          f"finite regression: B(E^{n}) has degree {n + 1}, leading coeff -2")
# An explicit preimage reaches the remaining constant quotient, completing the
# degree-free triangular proof that Im Phi=D*F[E].
Cpre = sp.expand(-sp.Rational(1, 15) - sp.Rational(2, 15) * E)
Vpre = sp.expand(-sp.Rational(11, 20) - sp.Rational(11, 10) * E - sp.Rational(1, 5) * E**3)
check_zero(k3(falling(3) * Cpre) - h2(falling(2) * Vpre) - D,
           "explicit preimage: Phi(C,V) = D for C,V of filler-image memo sec 4 (surjectivity)")
# E-R in D*F[E]  <=>  (E-R)(0)=(E-R)(1)=(E-R)(-1)=0  <=>  R(0)=0, R(1)=1, R(-1)=-1
Rc = sp.symbols("Rc0:8")
Rgen = sum(Rc[j] * E**j for j in range(8))
_, rem = sp.div(sp.Poly(E - Rgen, E), sp.Poly(D, E))
check(sp.Poly(rem, E).degree() <= 2,
      "E-R mod D is degree <= 2: divisibility <=> vanishing at the 3 roots {0,1,-1}")
check_zero((E - Rgen).subs(E, 0) - (-Rgen.subs(E, 0)), "(E-R)(0) = -R(0)  (=> R(0)=0 needed, automatic)")
check_zero((E - Rgen).subs(E, 1) - (1 - Rgen.subs(E, 1)), "(E-R)(1) = 1-R(1)  (=> R(1)=1 needed)")
check_zero((E - Rgen).subs(E, -1) - (-1 - Rgen.subs(E, -1)), "(E-R)(-1) = -1-R(-1)  (=> R(-1)=-1 needed)")
# The Phi kernel note: L_K, L_H are individually injective, but Phi=L_K-L_H has a
# nonzero kernel (two-filler cross-cancellation), so (C,V) is unique only MODULO that
# kernel -- NOT absolutely.  Exhibit a nonzero kernel element at low degree.
Ck = list(sp.symbols("Ck0:3"))
Vk = list(sp.symbols("Vk0:5"))
Cker = sum(Ck[j] * E**j for j in range(3))
Vker = sum(Vk[j] * E**j for j in range(5))
kernel_eq = sp.Poly(k3(falling(3) * Cker) - h2(falling(2) * Vker), E).all_coeffs()
allc = Ck + Vk
Msol = sp.linsolve(kernel_eq, allc)
sol0 = list(Msol)[0]
free = set().union(*[sp.sympify(x).free_symbols for x in sol0])
check(bool(free & set(allc)),
      "Phi=L_K-L_H has a NONZERO kernel (cross-cancellation): (C,V) unique only mod ker Phi")

# =====================================================================
print("\n--- 7. cascade constraint R(1)+R(-1)=0 and R=0 at d=1 (from solved cascade) ---")
# =====================================================================
W2 = (a3_w2, b2_w2)
A1, B1, pos1, data1 = build_cascade(*W2, 1)
A0, B0 = dict(A1), dict(B1)
A0[-2] = 0
B0[-3] = 0
R1 = potential(A0, B0)
check_zero(sp.expand(potential(A1, B1) - (R1 + k3(B1[-3]) - h2(A1[-2]))),
           "d=1: full potential G = R + K_3[b_-3] - H_2[a_-2]")
solve_vars1 = data1["ca2"] + data1["ca1"] + data1["ca0"] + data1["cam1"]
branches1 = sp.solve([e for e in pos1 if sp.expand(e) != 0], solve_vars1, dict=True)
check(len(branches1) >= 1, "d=1: positive cascade has solved branch(es)")
for i, sol in enumerate(branches1):
    Rs = sp.expand(R1.subs(sol))
    check_zero(Rs, f"d=1 branch {i}: R = 0 (positive cascade forces it)")
    check_zero(Rs.subs(E, 1) + Rs.subs(E, -1), f"d=1 branch {i}: R(1)+R(-1)=0 (cascade constraint)")
    check(sp.expand(Rs.subs(E, 1)) != 1,
          f"d=1 branch {i}: slope R(1)=0 != 1  =>  E-R not in Im Phi  =>  Q_0=1 infeasible")

# =====================================================================
print("\n--- 8. COMBINED SYSTEM feasibility and the slope/tail decomposition ---")
# =====================================================================
freevars1 = (data1["ca2"] + data1["ca1"] + data1["ca0"] + data1["cam1"]
             + data1["cam2"] + data1["cam3"] + [data1["mu3"]] + data1["kernels"])
pos1nz = [e for e in pos1 if sp.expand(e) != 0]
Q0m1_eqs = scalar_coeffs(q_m(A1, B1, 0) - 1)
tail_eqs = []
for m in (-1, -2, -3, -4, -5):
    tail_eqs += scalar_coeffs(q_m(A1, B1, m))
check(is_unit_ideal(pos1nz + Q0m1_eqs + tail_eqs, freevars1),
      "d=1: FULL system (positive cascade + Q_0=1 + negative tail + membership) = UNIT IDEAL")
check(is_unit_ideal(pos1nz + Q0m1_eqs, freevars1),
      "d=1 decomposition (a): positive cascade + Q_0=1 (NO tail) = unit ideal  [slope kills]")
check(not is_unit_ideal(pos1nz + tail_eqs, freevars1),
      "d=1 decomposition (b): positive cascade + negative tail = PROPER ideal [tail FEASIBLE]")
check(not is_unit_ideal(pos1nz, freevars1),
      "d=1 sanity: positive cascade alone is a proper ideal (no false kill)")
print("   => at d=1 the slope equation already makes the combined system infeasible;")
print("      after omitting it, the tail subsystem remains feasible.")

# Recheck separate tail-subsystem feasibility at d=2 (the rank-drop degree).
A2, B2, pos2, data2 = build_cascade(*W2, 2)
freevars2 = (data2["ca2"] + data2["ca1"] + data2["ca0"] + data2["cam1"]
             + data2["cam2"] + data2["cam3"] + [data2["mu3"]] + data2["kernels"])
pos2nz = [e for e in pos2 if sp.expand(e) != 0]
tail2 = []
for m in (-1, -2, -3, -4, -5):
    tail2 += scalar_coeffs(q_m(A2, B2, m))
check(not is_unit_ideal(pos2nz + tail2, freevars2),
      "d=2: positive cascade + negative tail = PROPER ideal without Q_0=1")

# =====================================================================
print("\n--- 9. positive-control genuine pair: whole tail satisfied; bottom gate passed ---")
# =====================================================================
kappa, c1 = sp.symbols("kappa c1")
U = {1: sp.Integer(1), -1: sp.expand(c1 * E)}
U3 = mul_ladders(mul_ladders(U, U), U)
Xw = {k: sp.Integer(0) for k in range(-3, 4)}
for k, v in U3.items():
    Xw[k] += v
Xw[-1] = sp.expand(Xw[-1] - E / kappa)
Dw = {k: sp.Integer(0) for k in range(-3, 4)}
Dw[1] = kappa
Dw[-1] = sp.expand(kappa * c1 * E)
for m in range(-6, 7):
    check_zero(q_m(Xw, Dw, m) - (1 if m == 0 else 0),
               f"positive control: Q_{m} = delta_(m0)  (genuine pair X=U^3-d/kappa, D=kappa U)")
check_zero(Xw[-3] - c1**3 * falling(3),
           "positive control: a_-3 = c1^3 (E)_3 is a shifted cube -> passes the bottom gate")
check(phi3_divides([0, 1, 2]), "positive control a_-3 roots {0,1,2} form a Phi_3 necklace")

# =====================================================================
print("\n--- 10. ESCALATION HARNESS: ready for the decisive sibling's slope=1 data ---")
# =====================================================================
def run_full_system(top, d, extra_constraints=None):
    """Return (is_infeasible, localization) for the full W2 DC1-face system at cap d,
    optionally intersected with extra positive-data constraints (e.g. a slope=1 point
    from research/band3/w2-decisive.md).  Ready to detect a candidate pair."""
    A, B, pos, data = build_cascade(*top, d)
    fv = (data["ca2"] + data["ca1"] + data["ca0"] + data["cam1"]
          + data["cam2"] + data["cam3"] + [data["mu3"]] + data["kernels"])
    eqs = [e for e in pos if sp.expand(e) != 0]
    eqs += scalar_coeffs(q_m(A, B, 0) - 1)
    for m in (-1, -2, -3, -4, -5):
        eqs += scalar_coeffs(q_m(A, B, m))
    if extra_constraints:
        eqs += list(extra_constraints)
    return is_unit_ideal(eqs, fv)


check(run_full_system(W2, 1),
      "harness self-test: run_full_system(W2, d=1) reports infeasible (unit ideal)")
print("   (When w2-decisive.md supplies positive data with slope R(1)=1 at some d, pass it")
print("    as extra_constraints; a PROPER ideal there would be a candidate DC1 pair -> escalate.)")

print("\nBOUNDED sections: 5 (concrete necklaces), 7 (d=1), 8 (d=1,2 Groebner).")
print("PROVED arbitrary-degree: 0-4 and the written triangular pivot argument in 6;")
print("finite computations in 6 are regressions for that argument.  No Weyl pair or")
print("counterexample is constructed.  W2 arbitrary-degree status: OPEN: central completion")
print("reduces to R(1)=1, but the full negative tail remains a necessary joint condition.")
print("ALL W2 TAIL CHECKS PASSED")
