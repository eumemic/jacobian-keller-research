#!/usr/bin/env python3
"""Exact arithmetic checks for main.tex. Requires SymPy."""
import warnings
import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning
warnings.filterwarnings("ignore", category=SymPyDeprecationWarning)

x, y, z, nu = sp.symbols("x y z nu")
a, b, c = sp.symbols("a b c")
t = x*y
u0 = 1 + t
s = x**2*z


def primitive_integral(poly):
    p = sp.Poly(poly, nu, domain=sp.QQ)
    den = sp.ilcm(*[q.q for q in p.all_coeffs()])
    return sp.Poly(sp.expand(den*p.as_expr()), nu, domain=sp.ZZ).primitive()[1]


def family(E, h):
    A = sp.integrate(nu*E, nu)
    C = 2*sp.integrate(E, nu)
    g = sp.expand(h - s)
    nux = sp.expand(u0*g)
    S = sp.expand(A.subs(nu, nux) + u0*g**2)
    T = sp.expand(C.subs(nu, nux) + 2*g)
    F = [sp.cancel(S/(x**2*g**2)), sp.cancel(T/(x*g)), x*g]
    return A, C, g, nux, F


def jacobian(F):
    return sp.expand(sp.det(sp.Matrix([
        [sp.diff(f, q) for q in (x, y, z)] for f in F
    ])))


def mod_poly(expr, modulus):
    num = sp.together(expr).as_numer_denom()[0]
    return sp.rem(sp.Poly(num, nu), sp.Poly(modulus, nu)).as_expr()


def verify_factorization(Q, prime, factors):
    product = sp.prod(factors)
    lc = int(Q.LC()) % prime
    plc = int(sp.Poly(product, nu).LC()) % prime
    scalar = lc * pow(plc, -1, prime) % prime
    assert sp.Poly(Q.as_expr() - scalar*product, nu, modulus=prime).is_zero
    coeff, actual = sp.factor_list(Q.as_expr(), modulus=prime)
    actual_degrees = sorted(sp.degree(f, nu) for f, exponent in actual for _ in range(exponent))
    expected_degrees = sorted(sp.degree(f, nu) for f in factors)
    assert actual_degrees == expected_degrees


# Shifted-Legendre universal family: exact moments and all three constraints.
r = sp.symbols("r")
for n in range(2, 13):
    Ln = sp.legendre(n, 2*nu - 1)
    assert sp.integrate(Ln, (nu, 0, 1)) == 0
    assert sp.integrate(nu*Ln, (nu, 0, 1)) == 0
    E = 2 - 6*nu + Ln
    A = sp.integrate(nu*E, nu)
    C = 2*sp.integrate(E, nu)
    assert sp.simplify(A.subs(nu, 1)) == -1
    assert sp.simplify(C.subs(nu, 1)) == -2
    assert sp.simplify(E.subs(nu, 1)) == -3
    assert sp.simplify(1 + E.subs(nu, 1) - 2*(E.subs(nu, 1) + 2)) == 0

# n=1 member and several Legendre members: polynomiality, determinant, exact degree.
members = [(2 - 6*nu, 1 - sp.Rational(3, 2)*t)]
for n in range(2, 6):
    members.append((2 - 6*nu + sp.legendre(n, 2*nu - 1), 1 - 2*t))
for E, h in members:
    A, C, g, nux, F = family(E, h)
    assert all(f.is_polynomial(x, y, z) for f in F)
    assert jacobian(F) == -2
    assert sp.degree(2*A - nu*C, nu) == sp.degree(E, nu) + 2

# The degree-zero constraints are inconsistent both for a != 0 and a = 0.
e, aa = sp.symbols("e aa", nonzero=True)
assert sp.solve([e*aa**2/2 + aa**2, 2*e*aa + 2*aa], [e, aa], dict=True) == []
assert sp.solve([e/2 + 1, 2*e + 2], [e], dict=True) == []

# Degree-five map, primitive equation, and rational reconstruction.
E5 = 4 - 15*nu + 10*nu**3
A5, C5, g5, nux5, F5 = family(E5, sp.Integer(1))
assert all(f.is_polynomial(x, y, z) for f in F5)
assert jacobian(F5) == -2
P5 = nu**5 - 5*nu**3 + 4*nu**2 - b*c*nu + 2*a*c**2
assert sp.factor(P5.subs({nu: nux5, a: F5[0], b: F5[1], c: F5[2]})) == 0
grec = (b*c - C5)/2
subs_rec = {
    x: c/grec,
    y: ((nu/grec) - 1)/(c/grec),
    z: (1 - grec)/(c/grec)**2,
}
for f, target in zip(F5, (a, b, c)):
    assert sp.factor(mod_poly(f.subs(subs_rec) - target, P5)) == 0

D5_expected = 16*c**2*(
    3125*a**4*c**6 + 37500*a**3*c**4 - 2500*a**2*b**2*c**4
    - 37125*a**2*b*c**3 - 1875*a**2*c**2 + 800*a*b**3*c**3
    + 7000*a*b**2*c**2 - 2700*a*b*c - 2176*a - 16*b**5*c**3
    - 200*b**4*c**2 + 95*b**3*c + 68*b**2)
assert sp.expand(sp.discriminant(P5, nu) - D5_expected) == 0

targets = [
    (sp.Rational(2, 3), sp.Rational(-1, 5), sp.Rational(3, 4)),
    (sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(-2, 7)),
    (sp.Rational(-3, 5), sp.Rational(4, 9), sp.Rational(1, 6)),
]
Q5_expected = [
    20*nu**5 - 100*nu**3 + 80*nu**2 + 3*nu + 15,
    147*nu**5 - 735*nu**3 + 588*nu**2 + 14*nu + 12,
    270*nu**5 - 1350*nu**3 + 1080*nu**2 - 20*nu - 9,
]
D5_values = [797038565664000, -244036829099277031104,
             17125713993062891250000]
Q5 = []
for target, expected, disc in zip(targets, Q5_expected, D5_values):
    Q = primitive_integral(P5.subs(dict(zip((a, b, c), target))))
    assert sp.expand(Q.as_expr() - expected) == 0
    assert sp.discriminant(Q.as_expr(), nu) == disc
    Q5.append(Q)
verify_factorization(Q5[0], 7, [nu**5 + 2*nu**3 - 3*nu**2 - 3*nu - 1])
assert sp.Poly(nu**5 + 2*nu**3 - 3*nu**2 - 3*nu - 1,
               nu, modulus=7).is_irreducible
verify_factorization(Q5[0], 103,
                     [nu+11, nu+18, nu+21, nu**2-50*nu+40])
assert all(d % p for d in D5_values[:1] for p in (7, 103))
assert not sp.integer_nthroot(abs(D5_values[0]*D5_values[1]), 2)[1]
assert not sp.integer_nthroot(abs(D5_values[0]*D5_values[2]), 2)[1]

# Degree-four map and its exact S4 certificates.
E4 = 3 - 12*nu + 6*nu**2
A4, C4, g4, nux4, F4 = family(E4, 1 - 2*t)
assert all(f.is_polynomial(x, y, z) for f in F4)
assert jacobian(F4) == -2
P4 = nu**4 - 4*nu**3 + 3*nu**2 - b*c*nu + 2*a*c**2
assert sp.factor(P4.subs({nu: nux4, a: F4[0], b: F4[1], c: F4[2]})) == 0
R_expected = [
    20*nu**4 - 80*nu**3 + 60*nu**2 + 3*nu + 15,
    147*nu**4 - 588*nu**3 + 441*nu**2 + 14*nu + 12,
    270*nu**4 - 1080*nu**3 + 810*nu**2 - 20*nu - 9,
]
D4_values = [-54873514800, -466029671480352, 4549902124608000]
R = []
for target, expected, disc in zip(targets, R_expected, D4_values):
    Q = primitive_integral(P4.subs(dict(zip((a, b, c), target))))
    assert sp.expand(Q.as_expr() - expected) == 0
    assert sp.discriminant(Q.as_expr(), nu) == disc
    R.append(Q)
verify_factorization(R[0], 11, [nu**4-4*nu**3+3*nu**2+4*nu-2])
assert sp.Poly(nu**4-4*nu**3+3*nu**2+4*nu-2, nu, modulus=11).is_irreducible
verify_factorization(R[0], 7, [nu+3, nu**3+3*nu+2])
verify_factorization(R[0], 19, [nu+7, nu+8, nu**2+4])
assert all(D4_values[0] % p for p in (7, 11, 19))
assert D4_values[0]*D4_values[2] < 0

print("Shifted-Legendre identities and universal constraints: PASS")
print("Family polynomiality, determinant -2, and exact degrees: PASS")
print("Degree-zero inconsistency: PASS")
print("Quintic equation, reconstruction, and S5 certificates: PASS")
print("Quartic equation and S4 certificates: PASS")
print("ALL CHECKS PASSED")
