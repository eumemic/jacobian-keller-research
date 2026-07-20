#!/usr/bin/env python3
"""Exact regression checks for the full classical band-2 assembly.

These symbolic identities support regression testing only; they are not a
proof of classification or case completeness.
"""
import sympy as sp

x, xi, t = sp.symbols("x xi t")


def poisson(g, f):
    """Convention used by the theorem: {g,f}=g_xi*f_x-g_x*f_xi."""
    return sp.expand(sp.diff(g, xi) * sp.diff(f, x)
                     - sp.diff(g, x) * sp.diff(f, xi))


def assert_zero(expr, label):
    value = sp.simplify(sp.expand(expr))
    if value != 0:
        raise AssertionError(f"{label}: {value}")
    print("PASS", label)


# Band 1: derive all five coefficient identities directly.
ks = range(-1, 2)
A = {k: sp.Function(f"a{k}")(t) for k in ks}
B = {k: sp.Function(f"b{k}")(t) for k in ks}


def component(m):
    return sp.expand(sum(
        k * A[k] * sp.diff(B[l], t) - l * sp.diff(A[k], t) * B[l]
        for k in ks for l in ks if k + l == m
    ))


f = sum(x**k * A[k].subs(t, x * xi) for k in ks)
g = sum(x**k * B[k].subs(t, x * xi) for k in ks)
formula = sum(x**m * component(m).subs(t, x * xi) for m in range(-2, 3))
assert_zero(poisson(g, f) - formula,
            "all five band-1 components match direct differentiation")

expected = {
    2: A[1] * sp.diff(B[1], t) - sp.diff(A[1], t) * B[1],
    1: A[1] * sp.diff(B[0], t) - sp.diff(A[0], t) * B[1],
    0: A[1] * sp.diff(B[-1], t) + sp.diff(A[1], t) * B[-1]
       - A[-1] * sp.diff(B[1], t) - sp.diff(A[-1], t) * B[1],
    -1: sp.diff(A[0], t) * B[-1] - A[-1] * sp.diff(B[0], t),
    -2: sp.diff(A[-1], t) * B[-1] - A[-1] * sp.diff(B[-1], t),
}
for m in range(-2, 3):
    assert_zero(component(m) - expected[m], f"band-1 m={m} identity")

# Check the target shear in the generic oriented band-1 reduction.
a1 = sp.Function("a1")(t)
am1 = sp.Function("am1")(t)
bm1 = sp.Function("bm1")(t)
lam = sp.symbols("lambda")
central = (a1 * sp.diff(bm1, t) + sp.diff(a1, t) * bm1
           - am1 * sp.diff(lam * a1, t)
           - sp.diff(am1, t) * lam * a1)
sheared_bm1 = bm1 - lam * am1
assert_zero(central - sp.diff(a1 * sheared_bm1, t),
            "band-1 target shear gives (a1*(b-1-lambda*a-1))'")

# Reflection Rf(x,xi)=f(xi,x): coefficient map and Poisson sign.
coeffs = {
    -2: t**2 * (t + 2),
    -1: t * (2 * t - 3),
    0: t**2 + 1,
    1: 3 * t + 4,
    2: t**2 - t + 5,
}
poly = sp.expand(sum(x**k * coeffs[k].subs(t, x * xi) for k in range(-2, 3)))
reflected_direct = sp.expand(poly.subs({x: xi, xi: x}, simultaneous=True))
reflected_coefficients = {
    k: sp.expand(t**(-k) * coeffs[-k]) for k in range(-2, 3)
}
reflected_formula = sp.expand(sum(
    x**k * reflected_coefficients[k].subs(t, x * xi)
    for k in range(-2, 3)
))
assert_zero(reflected_direct - reflected_formula,
            "representative reflection maps (Rf)_k=t^(-k)*a_-k")

f_test = x**2 + x * xi + 2 * xi**2 + x

g_test = 3 * x**2 - 2 * x * xi + xi**2 + xi
rf = f_test.subs({x: xi, xi: x}, simultaneous=True)
rg = g_test.subs({x: xi, xi: x}, simultaneous=True)
rbracket = poisson(g_test, f_test).subs({x: xi, xi: x}, simultaneous=True)
assert_zero(poisson(rg, rf) + rbracket,
            "representative reflection reverses the Poisson sign")

# M5 normal form: exact bracket and explicit polynomial inverse.
kappa = sp.symbols("kappa", nonzero=True)
c0, c1, Aconst, beta, lam = sp.symbols("c0 c1 A beta lambda")
U = x + c0 + c1 * xi
F = sp.expand(U**2 - xi / kappa - Aconst)
G = sp.expand(lam * F + kappa * U + beta)
assert_zero(poisson(G, F) - 1, "M5 normal form has bracket one")

normal_F_coefficients = {
    2: sp.Integer(1),
    1: 2 * c0,
    0: 2 * c1 * t + c0**2 - Aconst,
    -1: t * (2 * c0 * c1 - 1 / kappa),
    -2: c1**2 * t**2,
}
normal_G_coefficients = {
    2: lam,
    1: 2 * lam * c0 + kappa,
    0: lam * normal_F_coefficients[0] + kappa * c0 + beta,
    -1: lam * normal_F_coefficients[-1] + kappa * c1 * t,
    -2: lam * normal_F_coefficients[-2],
}
normal_F_ladder = sp.expand(sum(
    x**k * normal_F_coefficients[k].subs(t, x * xi)
    for k in range(-2, 3)
))
normal_G_ladder = sp.expand(sum(
    x**k * normal_G_coefficients[k].subs(t, x * xi)
    for k in range(-2, 3)
))
assert_zero(F - normal_F_ladder,
            "normal-form F has band-2 support and polynomial membership")
assert_zero(G - normal_G_ladder,
            "normal-form G has band-2 support and polynomial membership")

P, Q = sp.symbols("P Q")
Uinv = (Q - lam * P - beta) / kappa
xi_inv = sp.expand(kappa * (Uinv**2 - P - Aconst))
x_inv = sp.expand(Uinv - c0 - c1 * xi_inv)
forward_P = sp.expand(F.subs({x: x_inv, xi: xi_inv}, simultaneous=True))
forward_Q = sp.expand(G.subs({x: x_inv, xi: xi_inv}, simultaneous=True))
assert_zero(forward_P - P, "normal-form inverse recovers F coordinate")
assert_zero(forward_Q - Q, "normal-form inverse recovers G coordinate")

print("ALL FULL CLASSICAL BAND-2 REGRESSION CHECKS PASSED")
