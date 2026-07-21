#!/usr/bin/env python3
"""Exact finite corroboration for `two-filler-cross-cancellation.md`.

The symbolic identities and finite truncated matrices checked here corroborate
the note's arbitrary-degree linear-algebra proof; finite witnesses do not prove
that theorem. This script makes no full-Weyl-pair or remaining-ladder claim.

Run: uv run --with sympy python research/dc1-program/verify_two_filler_cross_cancellation.py
Ends: ALL TWO-FILLER CROSS-CANCELLATION CHECKS PASSED
"""
import sympy as sp

E = sp.symbols("E")


def sh(f, r):
    return sp.expand(sp.sympify(f).subs(E, E + r))


def falling(n):
    return sp.prod(E - i for i in range(n))


def shift_sum(f, n):
    return sp.expand(sum(sh(f, j) for j in range(n)))


def K_block(a, c, k):
    return sp.expand(sum(sh(a, -r) * sh(c, k - r)
                         for r in range(1, k + 1)))


def H_block(u, v, k):
    return sp.expand(sum(sh(u, -r) * sh(v, k - 1 - r)
                         for r in range(1, k)))


def az(expr, label):
    residual = sp.expand(sp.sympify(expr))
    if residual != 0:
        raise AssertionError(f"{label}: residual = {residual}")
    print("PASS", label)


def ok(condition, label):
    if not condition:
        raise AssertionError(f"{label}: FALSE")
    print("PASS", label)


def generic_poly(prefix, degree):
    coeffs = sp.symbols(f"{prefix}0:{degree + 1}")
    return sp.expand(sum(coeffs[j] * E**j for j in range(degree + 1))), coeffs


def coeff_vector(poly, max_degree):
    p = sp.Poly(sp.expand(poly), E)
    return sp.Matrix([p.nth(j) for j in range(max_degree + 1)])


def truncated_matrix(block, top, membership, input_degree, output_degree, k):
    """Matrix from free C coefficients to block(top, (E)_membership C)."""
    columns = []
    for j in range(input_degree + 1):
        image = block(top, sp.expand(falling(membership) * E**j), k)
        columns.append(coeff_vector(image, output_degree))
    return sp.Matrix.hstack(*columns)


print("--- 1. exact shift, degree/leading-coefficient, and telescoping identities ---")
for k in range(2, 6):
    # Algebraically independent symbolic coefficient examples; fixed degrees are
    # finite corroboration of the written arbitrary-degree formulas.
    a, ac = generic_poly(f"a{k}_", k + 1)
    c, cc = generic_poly(f"c{k}_", k + 2)
    u, uc = generic_poly(f"u{k}_", k)
    v, vc = generic_poly(f"v{k}_", k + 3)
    K = K_block(a, c, k)
    H = H_block(u, v, k)

    az(K - shift_sum(sh(a, -k) * c, k),
       f"k={k}: K_k[c] = S_k((T^-k a)c)")
    az(H - shift_sum(sh(u, -(k - 1)) * v, k - 1),
       f"k={k}: H_(k-1)[v] = S_(k-1)((T^-(k-1)u)v)")

    pK = sp.Poly(K, E)
    pH = sp.Poly(H, E)
    ok(pK.degree() == 2 * k + 3 and
       sp.expand(pK.LC() - k * ac[-1] * cc[-1]) == 0,
       f"k={k}: deg K=deg a+deg c and lc K=k lc(a)lc(c)")
    ok(pH.degree() == 2 * k + 3 and
       sp.expand(pH.LC() - (k - 1) * uc[-1] * vc[-1]) == 0,
       f"k={k}: deg H=deg u+deg v and lc H=(k-1)lc(u)lc(v)")

    az(sh(K, 1) - K - (a * sh(c, k) - sh(a, -k) * c),
       f"k={k}: (T-1)K=a(E)c(E+k)-a(E-k)c(E)")
    az(sh(H, 1) - H - (u * sh(v, k - 1) - sh(u, -(k - 1)) * v),
       f"k={k}: analogous telescoping identity for H")

print("\n--- 2. exact action of S_n on finite degree filtrations ---")
for n in range(1, 6):
    for degree in range(0, 7):
        matrix = sp.Matrix.hstack(*[
            coeff_vector(shift_sum(E**j, n), degree) for j in range(degree + 1)
        ])
        ok(matrix.det() == n ** (degree + 1),
           f"S_{n} on F[E]_<={degree}: triangular determinant {n}^{degree + 1}")

print("\n--- 3. truncated image-intersection certificates ---")


def certify_intersection(a, u, k, c_degree, v_degree, label):
    output_degree = sp.degree(a, E) + k + c_degree
    assert output_degree == sp.degree(u, E) + (k - 1) + v_degree
    MK = truncated_matrix(K_block, a, k, c_degree, output_degree, k)
    MH = truncated_matrix(H_block, u, k - 1, v_degree, output_degree, k)
    joined = MK.row_join(-MH)
    kernel = joined.nullspace()
    ok(kernel and any(MK * z[:c_degree + 1, :] != sp.zeros(output_degree + 1, 1)
                      for z in kernel),
       f"{label}: exact truncated matrix has a nonzero common image")
    z = next(z for z in kernel
             if MK * z[:c_degree + 1, :] != sp.zeros(output_degree + 1, 1))
    common = MK * z[:c_degree + 1, :]
    ok(common - MH * z[c_degree + 1:, :] == sp.zeros(output_degree + 1, 1),
       f"{label}: nullspace vector certifies M_K C=M_H V exactly")
    return MK, MH, joined


# A finite certificate independent of the displayed witnesses: increasing free
# degrees for one top pair yields exact nonzero intersections at every truncation.
for cap in range(0, 5):
    a0 = E + 2
    u0 = E**2 + E + 1
    k0 = 4
    cdeg = cap + 5
    vdeg = cap + 5
    MK, MH, joined = certify_intersection(
        a0, u0, k0, cdeg, vdeg, f"k=4 cap={cap}")
    ok(MK.rank() == cdeg + 1 and MH.rank() == vdeg + 1,
       f"k=4 cap={cap}: both truncated admissible maps are injective")

print("   Finite truncated certificates corroborate, but do not prove, the written")
print("   arbitrary-degree finite-codimension/intersection argument.")

print("\n--- 4. required k=2 shifted-power witness ---")
k = 2
a = sp.expand(E * (E - 1))
u = E - 1
c = sp.expand(E * (E - 1) / 2)
v = sp.expand(E * (E - 1)**2)
P = sp.expand(E * (E - 1)**2 * (E - 2))
az(sh(u, k) * a - sh(a, k - 1) * u,
   "k=2 shifted-power top: wall identity checked separately")
az(sp.rem(c, falling(k), E), "k=2 witness: c has (E)_2 membership")
az(sp.rem(v, falling(k - 1), E), "k=2 witness: v has (E)_1 membership")
az(K_block(a, c, k) - P, "k=2 witness: K_2[c]=E(E-1)^2(E-2)")
az(H_block(u, v, k) - P, "k=2 witness: H_1[v]=E(E-1)^2(E-2)")
MK, MH, joined = certify_intersection(a, u, k, 0, 2, "k=2 shifted-power witness")
witness_vector = sp.Matrix([sp.Rational(1, 2), 1, -2, 1])
ok(joined * witness_vector == sp.zeros(joined.rows, 1),
   "k=2 witness: displayed C=1/2 and V=(E-1)^2 give an exact nullspace vector")

print("\n--- 5. exotic k=3 two-block witness; wall and membership separate ---")
k = 3
a = sp.expand(E * (E - 2) * (E - 4))
u = sp.expand((E - 1) * (E - 4))
c = sp.expand(sp.Rational(2, 27) * E * (E - 1) * (E - 2) *
              (E - 4) * (9 * E - 26))
v = sp.expand(sp.Rational(1, 9) * E * (E - 1) * (E - 2) * (E - 4) *
              (9 * E**2 - 44 * E + 107))
# Wall is logically independent of the filler membership checks.
az(sh(u, k) * a - sh(a, k - 1) * u,
   "k=3 exotic top: wall u(E+3)a(E)=a(E+2)u(E)")
az(sp.rem(c, falling(k), E), "k=3 exotic witness: c has (E)_3 membership")
az(sp.rem(v, falling(k - 1), E), "k=3 exotic witness: v has (E)_2 membership")
P = K_block(a, c, k)
az(P - H_block(u, v, k),
   "k=3 exotic witness: exact two-block identity K_3[c]=H_2[v]")
ok(P != 0, "k=3 exotic witness: common output is nonzero")
ok(sp.degree(a, E) + sp.degree(c, E) == sp.degree(u, E) + sp.degree(v, E),
   "k=3 exotic witness: degree resonance holds")
az(k * sp.LC(sp.Poly(a, E)) * sp.LC(sp.Poly(c, E)) -
   (k - 1) * sp.LC(sp.Poly(u, E)) * sp.LC(sp.Poly(v, E)),
   "k=3 exotic witness: leading-coefficient resonance holds")
certify_intersection(a, u, k, 2, 4, "k=3 exotic witness truncation")

print("\nFinite witnesses corroborate the arbitrary-degree written proof; they do not prove it.")
print("No full Weyl pair or compatibility with the remaining ladder equations is asserted.")
print("ALL TWO-FILLER CROSS-CANCELLATION CHECKS PASSED")
