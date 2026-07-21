#!/usr/bin/env python3
"""Exact checks for quantum-ap-filler-image.md.

The memo contains the arbitrary-degree proofs. This verifier checks the displayed
finite certificates and runs explicitly labeled regressions of the written
annihilator, filtration, and quotient arguments.

Run: python research/band3/verify_quantum_ap_filler_image.py
Ends: ALL QUANTUM AP FILLER-IMAGE CHECKS PASSED
"""
import sympy as sp

E, r = sp.symbols("E r")


def shift(f, amount):
    """Exact polynomial shift f(E) -> f(E+amount)."""
    return sp.expand(sp.sympify(f).subs(E, E + amount))


def falling(n):
    """The membership factor (E)_n."""
    return sp.prod(E - j for j in range(n))


def divides(divisor, value):
    """Exact membership/divisibility test in the polynomial variable E."""
    value = sp.Poly(sp.expand(value), E)
    divisor = sp.Poly(sp.expand(divisor), E)
    return value.rem(divisor).is_zero


def k3(a, c):
    return sp.expand(sum(shift(a, j - 3) * shift(c, j) for j in range(3)))


def h2(b, v):
    return sp.expand(sum(shift(b, j - 2) * shift(v, j) for j in range(2)))


def ap_top():
    return (sp.expand((E-r)*(E-r-2)*(E-r-4)),
            sp.expand((E-r-1)*(E-r-4)))


def phi(a, b, C, V):
    return sp.expand(k3(a, falling(3)*C) - h2(b, falling(2)*V))


def coefficient_vector(value, cap):
    """Ascending coefficient vector, rejecting any accidental truncation."""
    polynomial = sp.Poly(sp.expand(value), E)
    actual_degree = -1 if polynomial.is_zero else polynomial.degree()
    if actual_degree > cap:
        raise ValueError(f"coefficient cap {cap} truncates degree {actual_degree}")
    result = sp.Matrix([polynomial.nth(j) for j in range(cap + 1)])
    reconstructed = sum(result[j] * E**j for j in range(cap + 1))
    if sp.expand(reconstructed - polynomial.as_expr()) != 0:
        raise AssertionError("coefficient vector reconstruction failed")
    return result


def phi_matrix(d):
    a, b = ap_top()
    cap = d + 6
    columns = [coefficient_vector(phi(a, b, E**j, 0), cap)
               for j in range(d + 1)]
    columns += [coefficient_vector(phi(a, b, 0, E**j), cap)
                for j in range(d + 1)]
    return sp.Matrix.hstack(*columns)


def check(condition, label):
    if not condition:
        raise AssertionError(label)
    print("PASS", label)


def check_zero(value, label):
    residual = sp.expand(value)
    if residual != 0:
        raise AssertionError(f"{label}: residual {residual}")
    print("PASS", label)


print("--- 0. exact helpers and truncation guard ---")
check_zero(shift(E**3 + 2*E, -2) - ((E-2)**3 + 2*(E-2)),
           "exact shift helper")
check(divides(falling(3), falling(3)*(E**4-r*E+1)),
      "exact membership helper accepts a multiple of (E)_3")
check(not divides(falling(3), falling(2)*(E+1)),
      "exact membership helper rejects a nonmember")
try:
    coefficient_vector(E**4 + 1, 3)
except ValueError:
    print("PASS coefficient truncation guard rejects degree overflow")
else:
    raise AssertionError("coefficient truncation guard did not fire")


print("\n--- 1. written-proof annihilator regressions ---")
a3, b2 = ap_top()
columns = [("C", j, phi(a3, b2, E**j, 0)) for j in range(9)]
columns += [("V", j, phi(a3, b2, 0, E**j)) for j in range(9)]
for kind, degree, column in columns:
    check_zero(column.subs(E, 0),
               f"[REGRESSION of written proof] {kind}_{degree}: evaluation at 0 annihilates Phi")
    lam = (column.subs(E, r+3) - column.subs(E, r+4)
           + column.subs(E, r+5) - column.subs(E, 0))
    check_zero(lam,
               f"[REGRESSION of written proof] {kind}_{degree}: lambda_r annihilates Phi")
check_zero((r+3)-(r+4)+(r+5)-(r+4),
           "lambda_r(E)=r+4 (independent from evaluation at 0 when r!=-4)")


print("\n--- 2. exact cap-4 base certificate ---")
M4 = phi_matrix(4)
rows1 = [2, 3, 4, 5, 6, 7, 8, 9, 10]
cols1 = [0, 1, 2, 3, 4, 5, 6, 7, 9]
rows2 = [1, 3, 4, 5, 6, 7, 8, 9, 10]
cols2 = [0, 1, 3, 4, 5, 6, 7, 8, 9]
D1 = sp.factor(M4.extract(rows1, cols1).det())
D2 = sp.factor(M4.extract(rows2, cols2).det())
expected_D1 = -648*(r+4)**2*(r+7)*(3*r**2+18*r+14)
expected_D2 = -648*(r+4)*(r**2+8*r+6)*(r**2+8*r+18)
check_zero(D1-expected_D1, "exact displayed cap-4 minor D1")
check_zero(D2-expected_D2, "exact displayed cap-4 minor D2")
check_zero(sp.factor(sp.gcd(D1, D2))-648*(r+4),
           "gcd(D1,D2)=648(r+4)")
A = sp.cancel(D1/(r+4))
B = sp.cancel(D2/(r+4))
bezout = ((15*r**3+58*r**2-406*r-2176)*A
          -(45*r**3+219*r**2-1284*r-7564)*B)
check_zero(bezout-23379840, "exact displayed Bezout identity")


print("\n--- 3. cap-3 anomaly certificate and filtered regressions ---")
q = r**2 + 8*r + 6
M3 = phi_matrix(3)
generic_minor = M3.subs(r, 0).extract(range(1, 9), range(8)).det()
check_zero(generic_minor + 13302123264,
           "[EXACT FINITE CERTIFICATE] cap-3 rank is generically eight")
kernel = sp.Matrix([
    -(3*r+16)/6, sp.Rational(2, 3), 0, 0,
    -(11*r+58)/2, (7*r+86)/4, -(7*r+34)/4, 1,
])
check(kernel[7] == 1,
      "[EXACT FINITE CERTIFICATE] displayed cap-3 q-kernel is nonzero")
check(all(sp.rem(sp.Poly(entry, r), sp.Poly(q, r)) == 0 for entry in M3*kernel),
      "[EXACT FINITE CERTIFICATE] displayed cap-3 source vector is a kernel modulo q")
rank_minor_rows = [2, 4, 5, 6, 7, 8, 9]
rank_minor_cols = [0, 1, 2, 3, 4, 5, 6]
rank_minor = sp.factor(M3.extract(rank_minor_rows, rank_minor_cols).det())
check_zero(rank_minor-972*(r+4)**2,
           "[EXACT FINITE CERTIFICATE] displayed cap-3 7x7 minor")
check(sp.gcd(sp.Poly(q, r), sp.Poly(r+4, r)).is_one,
      "q=0 does not meet r=-4, so the displayed 7x7 minor is nonzero")
check(sp.factor(sp.discriminant(q, r)) == 40,
      "q-roots require sqrt(10) (and q has no rational root)")

for d in range(4, 10):
    new_C = phi(a3, b2, E**d, 0)
    check(sp.Poly(new_C, E).degree() == d+6
          and sp.Poly(new_C, E).nth(d+6) == 3,
          f"[REGRESSION of filtered proof] C_{d} column has degree {d+6}, LC 3")
    inherited = []
    for j in range(d):
        inherited += [phi(a3, b2, E**j, 0), phi(a3, b2, 0, E**j)]
    check(all(sp.Poly(value, E).degree() <= d+5 for value in inherited),
          f"[REGRESSION of filtered proof] cap {d-1} columns have degree <= {d+5}")


print("\n--- 4. r=-4 quotient identities and arbitrary-degree proof regressions ---")
D = E*(E-1)*(E+1)
a_w2 = sp.expand(a3.subs(r, -4))
b_w2 = sp.expand(b2.subs(r, -4))
c_coeffs = sp.symbols("c0:6")
v_coeffs = sp.symbols("v0:6")
C_generic = sum(c_coeffs[j]*E**j for j in range(6))
V_generic = sum(v_coeffs[j]*E**j for j in range(6))
K = k3(a_w2, falling(3)*C_generic)
minus_H = -h2(b_w2, falling(2)*V_generic)
K_quotient = ((E-3)*(E-2)*(E-1)*C_generic
              + E*(E-2)*(E+2)*shift(C_generic, 1)
              + (E+1)*(E+2)*(E+3)*shift(C_generic, 2))
B_quotient = (2-E)*V_generic - (E+2)*shift(V_generic, 1)
check_zero(K-D*K_quotient,
           "[REGRESSION of written identity] generic-coefficient K/D quotient formula, including (E-1)")
check_zero(minus_H-D*B_quotient,
           "[REGRESSION of written identity] generic-coefficient (-H)/D quotient formula with explicit sign")
for j in range(9):
    check(divides(D, phi(a_w2, b_w2, E**j, 0))
          and divides(D, phi(a_w2, b_w2, 0, E**j)),
          f"[REGRESSION of written inclusion] W2 basis fillers degree {j} are divisible by D")

for n in range(10):
    Bn = sp.expand((2-E)*E**n - (E+2)*shift(E**n, 1))
    check(sp.Poly(Bn, E).degree() == n+1
          and sp.Poly(Bn, E).nth(n+1) == -2,
          f"[REGRESSION of triangular proof] B(E^{n}) has degree {n+1}, LC -2")

C_one = -sp.Rational(1, 15) - sp.Rational(2, 15)*E
V_one = (-sp.Rational(11, 20) - sp.Rational(11, 10)*E
         - sp.Rational(1, 5)*E**3)
quotient_preimage = sp.cancel(phi(a_w2, b_w2, C_one, V_one)/D)
check_zero(quotient_preimage-1,
           "[EXACT CERTIFICATE] displayed quotient preimage maps to 1")
check_zero(phi(a_w2, b_w2, C_one, V_one)-D,
           "[EXACT CERTIFICATE] displayed fillers map exactly to D")

print("\nALL QUANTUM AP FILLER-IMAGE CHECKS PASSED")
