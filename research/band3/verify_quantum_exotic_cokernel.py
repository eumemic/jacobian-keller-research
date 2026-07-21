#!/usr/bin/env python3
"""Exact finite cokernel computations for the W1/W2 exotic systems at d=1,2.

This verifier keeps every kernel parameter emitted by the positive-cascade solver.
It checks finite matrices only; it makes no arbitrary-degree or full-pair claim.

Run: python research/band3/verify_quantum_exotic_cokernel.py
Ends: ALL QUANTUM EXOTIC COKERNEL CHECKS PASSED
"""
import itertools
import sympy as sp

E, r = sp.symbols("E r")


def sh(f, n):
    return sp.expand(sp.sympify(f).subs(E, E + n))


def falling(n):
    return sp.prod(E - j for j in range(n))


def poly(name, degree):
    cs = list(sp.symbols(f"{name}_0:{degree + 1}"))
    return sp.expand(sum(cs[j] * E**j for j in range(degree + 1))), cs


def q_m(A, B, m):
    return sp.expand(sum(sh(B[l], k) * A[k] - sh(A[k], l) * B[l]
                         for k in range(-3, 4) for l in range(-3, 4)
                         if k + l == m))


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
    a2, ca2 = poly("a2", d)
    a1, ca1 = poly("a1", d)
    a0, ca0 = poly("a0", d)
    am1_raw, cam1 = poly("am1", d)
    am2_raw, cam2 = poly("am2", d)
    am3_raw, cam3 = poly("am3", d)
    mu3 = sp.symbols("mu3")
    A = {3: a3, 2: a2, 1: a1, 0: a0,
         -1: falling(1) * am1_raw, -2: falling(2) * am2_raw,
         -3: falling(3) * am3_raw}
    B = {k: sp.Integer(0) for k in range(-3, 4)}
    B[2] = b2
    B[-3] = sp.expand(mu3 * A[-3])
    conditions, kernels = [], []
    for m, lkey, name, membership, degree in [
            (4, 1, "b1c", 0, d + 3), (3, 0, "b0c", 0, 2*d + 2),
            (2, -1, "bm1c", 1, 2*d + 3), (1, -2, "bm2c", 2, 2*d + 4)]:
        B[lkey], new_kernels, new_conditions = clean_solve(
            A, B, m, lkey, name, membership, degree)
        kernels += new_kernels
        conditions += new_conditions
    data = {"ca2": ca2, "ca1": ca1, "ca0": ca0, "cam1": cam1,
            "cam2": cam2, "cam3": cam3, "mu3": mu3}
    return A, B, conditions, kernels, data


def potential(A, B):
    return sp.expand(sum(sh(A[k], j-k) * sh(B[-k], j)
                         - sh(B[k], j-k) * sh(A[-k], j)
                         for k in range(1, 4) for j in range(k)))


def k_block(a, c):
    return sp.expand(sum(sh(a, j-3) * sh(c, j) for j in range(3)))


def h_block(u, v):
    return sp.expand(sum(sh(u, j-2) * sh(v, j) for j in range(2)))


def coefficient_vector(f, degree):
    p = sp.Poly(sp.expand(f), E)
    if p.degree() > degree:
        raise AssertionError(
            f"coefficient cap {degree} truncates polynomial of degree {p.degree()}")
    vector = sp.Matrix([p.nth(j) for j in range(degree + 1)])
    if sp.expand(sum(vector[j] * E**j for j in range(degree + 1)) - p.as_expr()) != 0:
        raise AssertionError("coefficient vector does not reconstruct polynomial")
    return vector


def phi_matrix(a3, b2, d):
    output_degree = d + 6
    columns = [coefficient_vector(k_block(a3, falling(3) * E**j), output_degree)
               for j in range(d + 1)]
    columns += [-coefficient_vector(h_block(b2, falling(2) * E**j), output_degree)
                for j in range(d + 1)]
    return sp.Matrix.hstack(*columns)


def primitive(vector):
    denominator = sp.ilcm(*[x.q for x in vector])
    values = [int(x * denominator) for x in vector]
    divisor = sp.igcd(*values)
    values = [x // divisor for x in values]
    first = next(x for x in values if x)
    return sp.Matrix(values if first > 0 else [-x for x in values])


def check(condition, label):
    if not condition:
        raise AssertionError(label)
    print("PASS", label)


def ideal_equality(generators, certificate, variables):
    """Certify equality of two QQ ideals by mutual Groebner reduction."""
    generator_basis = sp.groebner(
        generators, *variables, order="lex", domain=sp.QQ)
    certificate_basis = sp.groebner(
        certificate, *variables, order="lex", domain=sp.QQ)
    return (all(certificate_basis.reduce(f)[1] == 0 for f in generators)
            and all(generator_basis.reduce(f)[1] == 0 for f in certificate))


def w1_d1_certificate(data):
    a20, a21 = data["ca2"]
    a10, a11 = data["ca1"]
    a00, a01 = data["ca0"]
    m0, m1 = data["cam1"]
    substitution = {
        a21: 0, a10: a20**2, a11: -a20**2/3,
        a01: 0, m0: 0, m1: 0,
    }
    return substitution, [variable-value for variable, value in substitution.items()]


def w1_d2_certificate(data):
    a20, a21, a22 = data["ca2"]
    a10, a11, a12 = data["ca1"]
    a00, a01, a02 = data["ca0"]
    m0, m1, m2 = data["cam1"]
    core = a20**2 - 10*a20*a22 + 24*a22**2
    substitution = {
        a21: -5*a22,
        a10: 8*a12 + core,
        a11: -(18*a12 + core)/3,
        a01: -a12*(7*a20 - 36*a22)/12,
        a02: a12*(a20 - 4*a22)/12,
        m0: -a12*(a20 - 6*a22)*(a20 - 4*a22)/36,
        m1: 2*a12**2/3,
        m2: -a12**2/6,
    }
    return substitution, [variable-value for variable, value in substitution.items()]


def w2_d1_certificate(data):
    a20, a21 = data["ca2"]
    a10, a11 = data["ca1"]
    a00, a01 = data["ca0"]
    m0, m1 = data["cam1"]
    substitution = {
        a21: 0, a10: -a20**2/3, a11: -a20**2/3,
        a01: 0, m0: 0, m1: 0,
    }
    return substitution, [variable-value for variable, value in substitution.items()]


CASCADE_CERTIFICATE_BUILDERS = {
    ("W1", 1): w1_d1_certificate,
    ("W1", 2): w1_d2_certificate,
    ("W2", 1): w2_d1_certificate,
}


W1 = (sp.expand(E*(E-2)*(E-4)), sp.expand((E-1)*(E-4)))
W2 = (sp.expand(E*(E+2)*(E+4)), sp.expand(E*(E+3)))
EXPECTED = {
    ("W1", 1): [[0,0,0,0],[-192,-12,8,-10],[414,-234,-30,5],[-348,450,16,-24],
                 [159,-309,-2,15],[-36,135,0,-2],[3,-33,0,0],[0,3,0,0]],
    ("W1", 2): [[0,0,0,0,0,0],[-192,-12,-72,8,-10,-10],[414,-234,-64,-30,5,-23],
                 [-348,450,-168,16,-24,8],[159,-309,521,-2,15,-19],[-36,135,-306,0,-2,14],
                 [3,-33,116,0,0,-2],[0,3,-30,0,0,0],[0,0,3,0,0,0]],
    ("W2", 1): [[0,0,0,0],[0,-12,0,2],[-18,-18,2,1],[0,-18,0,0],
                 [15,15,-2,-1],[0,27,0,-2],[3,3,0,0],[0,3,0,0]],
    ("W2", 2): [[0,0,0,0,0,0],[0,-12,-24,0,2,2],[-18,-18,-64,2,1,5],
                 [0,-18,-36,0,0,0],[15,15,17,-2,-1,-3],[0,27,54,0,-2,-2],
                 [3,3,44,0,0,-2],[0,3,6,0,0,0],[0,0,3,0,0,0]],
}
EXPECTED_RANK = {("W1", 1): 4, ("W1", 2): 6, ("W2", 1): 4, ("W2", 2): 5}

print("--- 0. symbolic structural identities in the finite band-3 setup ---")
generic_A, generic_B = {}, {}
for level in range(-3, 4):
    a_raw, _ = poly(f"generic_a_{level+3}", 2)
    b_raw, _ = poly(f"generic_b_{level+3}", 2)
    membership = falling(-level) if level < 0 else 1
    generic_A[level] = sp.expand(membership*a_raw)
    generic_B[level] = sp.expand(membership*b_raw)
check(sp.expand(q_m(generic_A, generic_B, 0)
                - (sh(potential(generic_A, generic_B), 1)
                   - potential(generic_A, generic_B))) == 0,
      "generic finite band-3 identity Q_0=(T-1)G")

print("\n--- 1. exact finite Phi matrices and cokernels ---")
matrices = {}
for tag, top in [("W1", W1), ("W2", W2)]:
    for d in (1, 2):
        M = phi_matrix(*top, d)
        matrices[tag, d] = M
        check(M == sp.Matrix(EXPECTED[tag, d]), f"{tag} d={d}: exact rational Phi matrix")
        check(M.rank() == EXPECTED_RANK[tag, d],
              f"{tag} d={d}: rank {M.rank()}, cokernel dimension {M.rows-M.rank()}")
        left = [primitive(x) for x in M.T.nullspace()]
        check(all((x.T*M).is_zero_matrix for x in left),
              f"{tag} d={d}: {len(left)} exact primitive left-null functionals")

print("\n--- 2. AP-parameter rank loci of the finite matrices ---")
a3r = sp.expand((E-r)*(E-r-2)*(E-r-4))
b2r = sp.expand((E-r-1)*(E-r-4))
for d, expected_gcd in [(1, sp.Integer(18)), (2, 36*(r+4))]:
    M = phi_matrix(a3r, b2r, d)
    minors = [sp.factor(M[list(rows), :].det())
              for rows in itertools.combinations(range(M.rows), M.cols)]
    minors = [z for z in minors if z != 0]
    gcd = minors[0]
    for z in minors[1:]:
        gcd = sp.gcd(gcd, z)
    check(sp.expand(gcd-expected_gcd) == 0,
          f"AP d={d}: gcd of maximal minors is {sp.factor(expected_gcd)}")
check(phi_matrix(a3r, b2r, 1).subs(r, -4).rank() == 4,
      "AP d=1: no rank drop at W2 (indeed maximal-minor gcd is constant)")
check(phi_matrix(a3r, b2r, 2).subs(r, -4).rank() == 5,
      "AP d=2: unique rank-drop locus r=-4, with rank 5")

print("\n--- 3. solved-cascade residuals, with all solver kernels retained ---")
for tag, top in [("W1", W1), ("W2", W2)]:
    for d in (1, 2):
        A, B, conditions, kernels, data = build_cascade(*top, d)
        check(kernels == [sp.symbols("b0cK0")],
              f"{tag} d={d}: sole operator-kernel parameter is retained constant b0cK0")
        nonfill_A, nonfill_B = dict(A), dict(B)
        nonfill_A[-2] = 0
        nonfill_B[-3] = 0
        R = potential(nonfill_A, nonfill_B)
        check(sp.expand(potential(A, B) - (
                  R + k_block(A[3], B[-3]) - h_block(B[2], A[-2]))) == 0,
              f"{tag} d={d}: full potential is R+K_3[c]-H_2[v]")
        check(sp.symbols("b0cK0") not in R.free_symbols,
              f"{tag} d={d}: canonical non-filler residual is independent of retained b0 kernel")
        solve_vars = data["ca2"] + data["ca1"] + data["ca0"] + data["cam1"]
        if (tag, d) in CASCADE_CERTIFICATE_BUILDERS:
            named_solution, certificate = CASCADE_CERTIFICATE_BUILDERS[tag, d](data)
            check(ideal_equality(conditions, certificate, solve_vars),
                  f"{tag} d={d}: graph ideal equals the positive-cascade ideal over QQ")
            check(all(sp.expand(f.subs(named_solution)) == 0 for f in conditions),
                  f"{tag} d={d}: named graph parametrization satisfies every condition")
        branches = sp.solve(conditions, solve_vars, dict=True, simplify=False)
        check(bool(branches), f"{tag} d={d}: positive cascade has solved branches")
        left = matrices[tag, d].T.nullspace()
        target_degree = d + 6
        for branch_number, solution in enumerate(branches):
            check(all(sp.expand(condition.subs(solution)) == 0 for condition in conditions),
                  f"{tag} d={d} branch {branch_number}: all positive-cascade residuals vanish")
            solved_A = {level: sp.expand(value.subs(solution))
                        for level, value in A.items()}
            solved_B = {level: sp.expand(value.subs(solution))
                        for level, value in B.items()}
            check(all(q_m(solved_A, solved_B, m) == 0 for m in range(1, 5)),
                  f"{tag} d={d} branch {branch_number}: Q_1,...,Q_4 vanish directly")
            solved_R = sp.expand(R.subs(solution))
            target = coefficient_vector(E-solved_R, target_degree)
            values = [sp.factor(f.dot(target)) for f in left]
            check(any(v != 0 for v in values),
                  f"{tag} d={d} branch {branch_number}: [E-R] is nonzero in finite coker Phi")

# Record the displayed residual and functional formulas explicitly.
A, B, conditions, kernels, data = build_cascade(*W1, 1)
w1d1_solution, _ = w1_d1_certificate(data)
B0 = dict(B); A0 = dict(A); A0[-2] = 0; B0[-3] = 0
check(sp.expand(potential(A0, B0).subs(w1d1_solution)) == 0, "W1 d=1: R=0")
w1d1_values = [sp.factor(f.dot(coefficient_vector(E, 7)))
               for f in matrices["W1", 1].T.nullspace()]
check(w1d1_values == [0, sp.Rational(25887, 99251),
                      -sp.Rational(11030, 99251), sp.Rational(1172, 99251)],
      "W1 d=1: displayed finite-cokernel functional evaluations")
A, B, conditions, kernels, data = build_cascade(*W2, 1)
w2d1_solution, _ = w2_d1_certificate(data)
B0 = dict(B); A0 = dict(A); A0[-2] = 0; B0[-3] = 0
check(sp.expand(potential(A0, B0).subs(w2d1_solution)) == 0, "W2 d=1: R=0")
w2d1_values = [sp.factor(f.dot(coefficient_vector(E, 7)))
               for f in matrices["W2", 1].T.nullspace()]
check(w2d1_values == [0, 1, 0, 0],
      "W2 d=1: displayed finite-cokernel functional evaluations")

# W1 d=2: two displayed left-null evaluations cannot vanish simultaneously.
A, B, conditions, kernels, data = build_cascade(*W1, 2)
w1d2_solution, _ = w1_d2_certificate(data)
A0, B0 = dict(A), dict(B); A0[-2] = 0; B0[-3] = 0
R = sp.expand(potential(A0, B0).subs(w1d2_solution))
expected_R = sp.expand(E * data["ca1"][2]**2 / 18 * (
    2*E**3*data["ca2"][0] - 12*E**3*data["ca2"][2]
    - 15*E**2*data["ca2"][0] + 92*E**2*data["ca2"][2]
    + 20*E*data["ca2"][0] - 140*E*data["ca2"][2]
    + 15*data["ca2"][0] - 44*data["ca2"][2]))
check(sp.expand(R-expected_R) == 0, "W1 d=2: displayed solved residual R")
w = data["ca1"][2]**2*(data["ca2"][0]-4*data["ca2"][2])
left = matrices["W1", 2].T.nullspace()
values = [sp.factor(f.dot(coefficient_vector(E-R, 8))) for f in left]
check(sp.expand(values[1]-(sp.Rational(77587,5759959)-sp.Rational(55275,5759959)*w)) == 0,
      "W1 d=2: first target functional is (77587-55275w)/5759959")
check(sp.expand(values[2]-sp.Rational(2,17279877)*(15982*w-22329)) == 0,
      "W1 d=2: second target functional is 2(15982w-22329)/17279877")
check(sp.Rational(77587,55275) != sp.Rational(22329,15982),
      "W1 d=2: the two cokernel equations are incompatible")

# W2 d=2: certify the full two-branch cover without trusting solve ordering.
A, B, conditions, kernels, data = build_cascade(*W2, 2)
a20, a21, a22 = data["ca2"]
a10, a11, a12 = data["ca1"]
a00, a01, a02 = data["ca0"]
m0, m1, m2 = data["cam1"]
certificate = [
    a21 - 3*a22,
    3*a10 + a20**2 - 2*a20*a22,
    3*a11 - 6*a12 + a20**2 - 2*a20*a22,
    3*a01 - 3*a02 - 2*a12*a22,
    2*a02*a20 - 4*a02*a22 - a12**2 + 6*m0 - 6*m2,
    m1,
    a12*a20,
]
groebner_vars = [a21, a10, a11, a01, m0, m1,
                  a20, a12, a22, a02, m2, a00]
check(ideal_equality(conditions, certificate, groebner_vars),
      "W2 d=2: Groebner certificate equals the positive-cascade ideal")

common = {
    a21: 3*a22,
    a10: -a20*(a20-2*a22)/3,
    a11: -a20*(a20-2*a22)/3 + 2*a12,
    a01: a02 + 2*a12*a22/3,
    m0: m2 - a02*(a20-2*a22)/3 + a12**2/6,
    m1: 0,
}
reduced_conditions = [sp.factor(f.subs(common)) for f in conditions]
check(all(f == 0 or sp.rem(f, a12*a20, a12) == 0
          for f in reduced_conditions),
      "W2 d=2: common substitutions leave only the factor a1_2*a2_0")
check(sp.factor(certificate[-1]) == a12*a20,
      "W2 d=2: exhaustive branches are a1_2=0 or a2_0=0")

A0, B0 = dict(A), dict(B)
A0[-2] = 0
B0[-3] = 0
unsolved_R = potential(A0, B0)
branch_a_R = sp.expand(unsolved_R.subs(common).subs(a12, 0))
branch_b_R = sp.expand(unsolved_R.subs(common).subs(a20, 0))
expected_a_R = sp.expand(
    -sp.Rational(2, 3)*E**2*m2*(E-1)*(E+1)*(a20-2*a22))
expected_b_R = sp.expand(
    sp.Rational(2, 3)*E*(E-1)*(E+1)*(2*E*a22*m2+a02*a12))
check(sp.expand(branch_a_R-expected_a_R) == 0,
      "W2 d=2 branch a1_2=0: displayed residual R")
check(sp.expand(branch_b_R-expected_b_R) == 0,
      "W2 d=2 branch a2_0=0: displayed residual R")
check(all(sp.expand(f.subs(common).subs(a12, 0)) == 0 for f in conditions)
      and all(sp.expand(f.subs(common).subs(a20, 0)) == 0 for f in conditions),
      "W2 d=2: both named branches satisfy every positive-cascade condition")

left = matrices["W2", 2].T.nullspace()
def target_values(residual):
    return [sp.factor(f.dot(coefficient_vector(E-residual, 8))) for f in left]

branch_a_values = target_values(branch_a_R)
branch_b_values = target_values(branch_b_R)
expected_b_values = [0, (a02*a12+9)/9, -a02*a12/9, 0]
check(branch_a_values == [0, 1, 0, 0],
      "W2 d=2 branch a1_2=0: complete target-functional vector")
check(all(sp.expand(x-y) == 0 for x, y in zip(branch_b_values, expected_b_values)),
      "W2 d=2 branch a2_0=0: complete target-functional vector")
check(sp.expand(branch_b_values[1]+branch_b_values[2]-1) == 0,
      "W2 d=2 branch a2_0=0: two target functionals sum to 1")

print("\nFinite computation only; unrestricted filler-image intersection remains infinite-dimensional.")
print("ALL QUANTUM EXOTIC COKERNEL CHECKS PASSED")
