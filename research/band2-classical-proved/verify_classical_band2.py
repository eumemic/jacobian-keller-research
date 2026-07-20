#!/usr/bin/env python3
"""Exact sign/equation checks for classical band 2.

Convention: {G,F} = G_xi F_x - G_x F_xi, tau=x*xi.
This checks identities and examples only; it is not a completeness proof.
"""
import sympy as sp

x, xi, t = sp.symbols("x xi tau")
ks = range(-2, 3)
A = {k: sp.Function(f"a{k}")(t) for k in ks}
B = {k: sp.Function(f"b{k}")(t) for k in ks}


def ck_component(m):
    return sp.expand(sum(
        k*A[k]*sp.diff(B[l], t) - l*sp.diff(A[k], t)*B[l]
        for k in ks for l in ks if k+l == m
    ))


def poisson(g, f):
    return sp.diff(g, xi)*sp.diff(f, x) - sp.diff(g, x)*sp.diff(f, xi)


def assert_zero(expr, label):
    if sp.simplify(sp.expand(expr)) != 0:
        raise AssertionError(label + ": " + str(sp.expand(expr)))
    print("PASS", label)


# Directly derive every component from the two-variable Poisson bracket.
f = sum(x**k * A[k].subs(t, x*xi) for k in ks)
g = sum(x**k * B[k].subs(t, x*xi) for k in ks)
direct = sp.expand(poisson(g, f))
formula = sp.expand(sum(x**m * ck_component(m).subs(t, x*xi) for m in range(-4, 5)))
assert_zero(direct-formula, "all nine CK components match direct differentiation")

# Freeze the displayed extreme and cross-term formulas.
expected = {
  4: 2*A[2]*sp.diff(B[2],t)-2*sp.diff(A[2],t)*B[2],
  3: 2*A[2]*sp.diff(B[1],t)-sp.diff(A[2],t)*B[1]
     +A[1]*sp.diff(B[2],t)-2*sp.diff(A[1],t)*B[2],
  2: 2*A[2]*sp.diff(B[0],t)+A[1]*sp.diff(B[1],t)
     -sp.diff(A[1],t)*B[1]-2*sp.diff(A[0],t)*B[2],
  1: 2*A[2]*sp.diff(B[-1],t)+sp.diff(A[2],t)*B[-1]
     +A[1]*sp.diff(B[0],t)-sp.diff(A[0],t)*B[1]
     -A[-1]*sp.diff(B[2],t)-2*sp.diff(A[-1],t)*B[2],
  0: 2*A[2]*sp.diff(B[-2],t)+2*sp.diff(A[2],t)*B[-2]
     +A[1]*sp.diff(B[-1],t)+sp.diff(A[1],t)*B[-1]
     -A[-1]*sp.diff(B[1],t)-sp.diff(A[-1],t)*B[1]
     -2*A[-2]*sp.diff(B[2],t)-2*sp.diff(A[-2],t)*B[2],
 -1: A[1]*sp.diff(B[-2],t)+2*sp.diff(A[1],t)*B[-2]
     +sp.diff(A[0],t)*B[-1]-A[-1]*sp.diff(B[0],t)
     -2*A[-2]*sp.diff(B[1],t)-sp.diff(A[-2],t)*B[1],
 -2: 2*sp.diff(A[0],t)*B[-2]+sp.diff(A[-1],t)*B[-1]
     -A[-1]*sp.diff(B[-1],t)-2*A[-2]*sp.diff(B[0],t),
 -3: 2*sp.diff(A[-1],t)*B[-2]-A[-1]*sp.diff(B[-2],t)
     +sp.diff(A[-2],t)*B[-1]-2*A[-2]*sp.diff(B[-1],t),
 -4: 2*sp.diff(A[-2],t)*B[-2]-2*A[-2]*sp.diff(B[-2],t),
}
for m in range(-4, 5):
    assert_zero(ck_component(m)-expected[m], f"m={m} displayed equation")

# The m=3 homogeneous reduction after b2=lambda*a2.
lam = sp.symbols("lambda")
u = sp.Function("u")(t)
a2 = sp.Function("a2")(t)
a1 = sp.Function("a1")(t)
b1 = lam*a1 + u
m3 = 2*a2*sp.diff(b1,t)-sp.diff(a2,t)*b1 \
     +a1*sp.diff(lam*a2,t)-2*sp.diff(a1,t)*lam*a2
assert_zero(m3-(2*a2*sp.diff(u,t)-sp.diff(a2,t)*u),
            "m=3 reduces to 2*a2*u'-a2'*u")

# Square and nonsquare representatives.
h = t**2+t+1
assert_zero(2*(3*h**2)*sp.diff(5*h,t)-sp.diff(3*h**2,t)*(5*h),
            "square-class representative has nonzero homogeneous solution")
assert_zero(2*t*sp.diff(sp.Integer(0),t)-sp.diff(t,t)*0,
            "nonsquare representative has the zero homogeneous solution")
c = sp.symbols("c0:6")
u_ansatz = sum(c[i]*t**i for i in range(6))
h3_ansatz = sp.Poly(sp.expand(2*t*sp.diff(u_ansatz,t)-u_ansatz), t)
sol = sp.solve(h3_ansatz.all_coeffs(), c, dict=True)
if sol != [{ci: 0 for ci in c}]:
    raise AssertionError("a2=t bounded polynomial kernel was not zero: " + str(sol))
print("PASS a2=t has zero polynomial kernel through degree five")

# Reduced nonsquare localized branch: F=x^2*tau, G=(1/2)x^-2.
floc = x**2*(x*xi)
gloc = sp.Rational(1,2)*x**-2
assert_zero(poisson(gloc, floc)-1,
            "localized nonsquare branch has bracket one")
if sp.rem(sp.Poly(sp.Rational(1,2), t), sp.Poly(t**2, t)) == 0:
    raise AssertionError("localized branch incorrectly passed polynomiality")
print("PASS localized branch fails tau^2 divisibility at degree -2")

# Affine symplectic pair represented in band 1, as a baseline.
faff = 2*x + 5 - 3*xi
gaff = x + 7 - xi
assert_zero(poisson(gaff, faff)-1, "affine symplectic baseline")

print("ALL CLASSICAL BAND-2 CHECKS PASSED")
