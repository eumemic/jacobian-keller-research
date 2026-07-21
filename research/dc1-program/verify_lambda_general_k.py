#!/usr/bin/env python3
"""Exact certificate for the general-k cokernel-dual (annihilating-functional) theory.

Companion to research/dc1-program/lambda-general-k.md.

The structural question: WHY does an annihilating functional of the two-filler
image Phi exist, and does it exist at every band k?  The answer is the adjoint
(cokernel-dual) picture.  With  S_n = 1 + T + ... + T^(n-1)  a linear automorphism
of F[E] (two-filler memo, commit 0889f8a/e4e704f tree) and

    Im L_K = S_k(q_K F[E]),   q_K = (T^-k a)(E)_k,
    Im L_H = S_{k-1}(q_H F[E]),  q_H = (T^-(k-1) u)(E)_{k-1},

a functional lambda annihilates Im L_K iff the box-k moving sum of its
coefficients is supported on the roots of q_K (equivalently S_k^* lambda lies in
the annihilator of the ideal q_K F[E]).  Ann(Im Phi) is the intersection of the
two transported annihilators; its dimension is codim Im Phi.

This file checks, over exact rationals:
  0. S_n automorphism facts and the adjoint S_n^* ev_x = ev_x+...+ev_{x+n-1}.
  1. The shift forms and the moving-sum annihilator criterion (bands 2-5).
  2. ev_0 in Ann(Im Phi) at every band (coker Phi != 0 arbitrary k), degree-free.
  3. Exact W1 recovery of the e4e704f certificate {ev_0, lambda_0}, codim 2,
     lambda_0(E)=4=common root, lambda_0(R)=0, lambda_0(E-R)=4.
  4. The AP band-3 family for symbolic r and arbitrary degree: lambda_r
     annihilates Im Phi, lambda_r(E)=r+4=common root, escape hatch r=-4.
  5. dim Ann(Im Phi) equals the stable finite cokernel dimension (cross-check).
  6. Bands 4 and 5 exotic tops: codim >= 2, an explicit annihilator kills Im Phi,
     and lambda(E) is NOT k+1 (disambiguates the W1 coincidence 4 = k+1).

Run: uv run --with sympy python research/dc1-program/verify_lambda_general_k.py
Ends: ALL LAMBDA GENERAL K CHECKS PASSED
"""
import sympy as sp

E = sp.symbols("E")


# ---------------------------------------------------------------------------
# primitives
# ---------------------------------------------------------------------------
def sh(f, n):
    return sp.expand(sp.sympify(f).subs(E, E + n))


def falling(n):
    return sp.prod(E - j for j in range(n))


def Kk(a, c, k):
    """K_k[c] = sum_{r=1}^k a(E-r) c(E+k-r)."""
    return sp.expand(sum(sh(a, -r) * sh(c, k - r) for r in range(1, k + 1)))


def Hk(u, v, k):
    """H_{k-1}[v] = sum_{r=1}^{k-1} u(E-r) v(E+k-1-r)."""
    return sp.expand(sum(sh(u, -r) * sh(v, k - 1 - r) for r in range(1, k)))


def Sn(g, n):
    """S_n g = g + Tg + ... + T^{n-1} g."""
    return sp.expand(sum(sh(g, i) for i in range(n)))


def moving_sum(coeffs, n):
    """coefficient function of S_n^* lambda, where lambda = sum coeffs[x] ev_x.
    (S_n^* lambda)(f) = sum_x coeffs[x] (f(x)+...+f(x+n-1))."""
    out = {}
    for x, w in coeffs.items():
        for i in range(n):
            out[x + i] = out.get(x + i, 0) + w
    return {y: w for y, w in out.items() if w != 0}


def integer_roots(poly):
    out = {}
    for root, m in sp.roots(sp.Poly(poly, E)).items():
        if getattr(root, "is_integer", False):
            out[int(root)] = m
    return out


def apply_functional(poly, coeffs):
    """evaluate lambda = sum coeffs[x] ev_x on a polynomial."""
    return sp.expand(sum(w * poly.subs(E, x) for x, w in coeffs.items()))


def qK_qH(a, u, k):
    return sp.expand(sh(a, -k) * falling(k)), sp.expand(sh(u, -(k - 1)) * falling(k - 1))


def annihilator(a, u, k, window, orders=3):
    """Basis of the (integer-supported, order<orders) annihilator of Im Phi via
    the moving-sum criterion.  Returns (basis_vectors, var_index) or None if
    q_K,q_H are not fully integer-rooted (a separate coset analysis then applies).
    var_index maps (support_point, derivative_order) -> column."""
    qK, qH = qK_qH(a, u, k)
    ZK, ZH = integer_roots(qK), integer_roots(qH)
    if sp.Poly(qK, E).degree() != sum(ZK.values()):
        return None
    if sp.Poly(qH, E).degree() != sum(ZH.values()):
        return None
    xs = list(range(-window, window + 1))
    var = {}
    n = 0
    for o in range(orders):
        for x in xs:
            var[(x, o)] = n
            n += 1
    rows = []

    def add(width, Z):
        # coefficient of ev_y^{(o)} in S_width^* lambda must vanish unless o<mult(y)
        for o in range(orders):
            for y in range(-window, window + width + 2):
                if o < Z.get(y, 0):
                    continue
                row = [0] * n
                nz = False
                for i in range(width):
                    key = (y - i, o)
                    if key in var:
                        row[var[key]] += 1
                        nz = True
                if nz:
                    rows.append(row)

    add(k, ZK)
    add(k - 1, ZH)
    return sp.Matrix(rows).nullspace(), var


def vec_to_coeffs(vec, var):
    """order-0 part -> {x: weight} (all our certified annihilators are order 0)."""
    d = {}
    for (x, o), j in var.items():
        if o == 0 and vec[j] != 0:
            d[x] = vec[j]
    return d


def lamE(vec, var):
    """value of the functional (encoded by vec) on f = E."""
    t = 0
    for (x, o), j in var.items():
        if o == 0:
            t += vec[j] * x
        elif o == 1:
            t += vec[j]
    return sp.nsimplify(t)


def true_codim(a, u, k, d):
    """codim of the actual Phi image at coefficient cap d (finite truncation)."""
    cols = []
    for i in range(d + 1):
        cols.append(sp.Poly(Kk(a, sp.expand(falling(k) * E**i), k), E))
    for i in range(d + 1):
        cols.append(sp.Poly(-Hk(u, sp.expand(falling(k - 1) * E**i), k), E))
    Dmax = max(p.degree() for p in cols)
    M = sp.zeros(Dmax + 1, len(cols))
    for j, p in enumerate(cols):
        for i, ci in enumerate(p.all_coeffs()[::-1]):
            M[i, j] = ci
    return (Dmax + 1) - M.rank()


PASS = 0


def check(cond, label):
    global PASS
    if not cond:
        raise AssertionError("FAIL: " + label)
    PASS += 1
    print("PASS", label)


def check_zero(val, label):
    check(sp.expand(val) == 0, label)


# ===========================================================================
print("=== 0. S_n automorphism and the adjoint S_n^* ===")
# generic degree-4 polynomial: deg S_n f = deg f, lc(S_n f) = n lc(f)
gc = sp.symbols("g0:5")
gpoly = sum(gc[i] * E**i for i in range(5))
for n in (2, 3, 4, 5):
    Sg = Sn(gpoly, n)
    P = sp.Poly(Sg, E)
    check(P.degree() == 4, f"deg S_{n}(f) = deg f (generic degree 4)")
    check_zero(P.LC() - n * gc[4], f"lc S_{n}(f) = {n}*lc(f)")

# adjoint identity: for any lambda (finite support) and any g,
#   lambda(S_n g) = (S_n^* lambda)(g),  with S_n^* ev_x = ev_x+...+ev_{x+n-1}.
lam_test = {0: sp.Integer(2), 3: sp.Integer(-1), 5: sp.Integer(7)}
for n in (2, 3, 4, 5):
    lhs = apply_functional(Sn(gpoly, n), lam_test)
    rhs = apply_functional(gpoly, moving_sum(lam_test, n))
    check_zero(lhs - rhs, f"adjoint identity lambda(S_{n} g)=(S_{n}^* lambda)(g)")

# ===========================================================================
print("\n=== 1. shift forms and the moving-sum annihilator criterion ===")
tops = {
    3: (sp.expand(E * (E - 2) * (E - 4)), sp.expand((E - 1) * (E - 4))),        # W1
    4: (sp.expand(E * (E - 2) * (E - 3) * (E - 5)), sp.expand((E - 1) * (E - 3) * (E - 5))),
}
for k, (a, u) in tops.items():
    # shift forms K_k[c]=S_k((T^-k a) c), H_{k-1}[v]=S_{k-1}((T^-(k-1) u) v)
    ctest = sp.expand(E**2 - 2)
    vtest = sp.expand(E + 3)
    check_zero(Kk(a, ctest, k) - Sn(sp.expand(sh(a, -k) * ctest), k),
               f"k={k}: K_k[c] = S_k((T^-k a) c)")
    check_zero(Hk(u, vtest, k) - Sn(sp.expand(sh(u, -(k - 1)) * vtest), k - 1),
               f"k={k}: H_(k-1)[v] = S_(k-1)((T^-(k-1) u) v)")

# criterion: if S_k^* lambda is supported on roots(q_K) and S_(k-1)^* lambda on
# roots(q_H) then lambda kills Phi(C,V) for every admissible C,V.
k = 3
a, u = tops[3]
qK, qH = qK_qH(a, u, k)
lam0 = {3: 1, 4: -1, 5: 1, 0: -1}
check(set(moving_sum(lam0, k)).issubset(set(integer_roots(qK))),
      "W1: S_3^* lambda_0 supported on roots(q_K)")
check(set(moving_sum(lam0, k - 1)).issubset(set(integer_roots(qH))),
      "W1: S_2^* lambda_0 supported on roots(q_H)")
Cc = sp.symbols("C0:4")
Vv = sp.symbols("V0:4")
Cp = sum(Cc[i] * E**i for i in range(4))
Vp = sum(Vv[i] * E**i for i in range(4))
c = sp.expand(falling(k) * Cp)
v = sp.expand(falling(k - 1) * Vp)
Phi = sp.expand(Kk(a, c, k) - Hk(u, v, k))
check_zero(apply_functional(Phi, lam0),
           "W1: lambda_0(Phi(C,V))=0 for generic admissible C,V (criterion works)")

# ===========================================================================
print("\n=== 2. ev_0 in Ann(Im Phi) at every band (coker != 0, degree-free) ===")
# K_k[(E)_k C](0)=0 and H_(k-1)[(E)_(k-1) V](0)=0 for arbitrary C,V, any k.
for k in (2, 3, 4, 5):
    a = sp.prod(E - t for t in range(0, 2 * k, 2))
    u = sp.prod(E - t for t in range(1, 2 * (k - 1), 2)) if k > 2 else (E - 1)
    a, u = sp.expand(a), sp.expand(u)
    Cp = sum(sp.symbols(f"c0:{4}")[i] * E**i for i in range(4))
    Vp = sum(sp.symbols(f"d0:{4}")[i] * E**i for i in range(4))
    c = sp.expand(falling(k) * Cp)
    v = sp.expand(falling(k - 1) * Vp)
    val = sp.expand((Kk(a, c, k) - Hk(u, v, k)).subs(E, 0))
    check_zero(val, f"k={k}: ev_0 annihilates Im Phi (membership forces value 0 at E=0)")
# structural reason: the length-k falling factorial evaluated at 0..k-1 vanishes
for k in (2, 3, 4, 5):
    for m in range(k):
        check_zero(falling(k).subs(E, m), f"k={k}: (E)_k vanishes at E={m}")

# ===========================================================================
print("\n=== 3. exact W1 recovery of the e4e704f certificate ===")
k = 3
a3 = sp.expand(E * (E - 2) * (E - 4))
b2 = sp.expand((E - 1) * (E - 4))
qK, qH = qK_qH(a3, b2, k)
check(sorted(integer_roots(qK)) == [0, 1, 2, 3, 5, 7], "W1 roots(q_K)={0,1,2,3,5,7}")
check(sorted(integer_roots(qH)) == [0, 1, 3, 6], "W1 roots(q_H)={0,1,3,6}")

ns, var = annihilator(a3, b2, k, window=14)
check(len(ns) == 2, "W1: 2 independent point annihilators (= stable codim Im Phi 2)")
# recover span{ev_0, lambda_0}
found = {frozenset(vec_to_coeffs(v, var).items()) for v in ns}
# express: ev_0 and (ev_3-ev_4+ev_5) both lie in Ann; check membership by residual
def in_span(target, basis, var):
    cols = sp.Matrix([[v[j] for j in range(len(v))] for v in basis]).T
    tvec = sp.zeros(cols.rows, 1)
    for (x, o), j in var.items():
        tvec[j] = target.get((x, o), 0)
    sol = cols.solve_least_squares(tvec)
    return sp.expand((cols * sol - tvec).norm()) == 0

check(in_span({(0, 0): 1}, ns, var), "ev_0 in Ann(Im Phi) at W1")
check(in_span({(3, 0): 1, (4, 0): -1, (5, 0): 1}, ns, var),
      "ev_3-ev_4+ev_5 in Ann(Im Phi) at W1")

# lambda_0 and its arithmetic
lam0 = {3: 1, 4: -1, 5: 1, 0: -1}
check_zero(apply_functional(E, lam0) - 4, "lambda_0(E)=4")
# 4 is the common root of a3 and b2
check_zero(a3.subs(E, 4), "a3(4)=0  (common root)")
check_zero(b2.subs(E, 4), "b2(4)=0  (common root)")
# Im Phi subset ker lambda_0 (arbitrary admissible c,v)
cc = sp.expand(falling(3) * sum(sp.symbols("p0:4")[i] * E**i for i in range(4)))
vv = sp.expand(falling(2) * sum(sp.symbols("q0:4")[i] * E**i for i in range(4)))
check_zero(apply_functional(Kk(a3, cc, 3), lam0), "lambda_0(K_3[c])=0 (arbitrary c)")
check_zero(apply_functional(Hk(b2, vv, 3), lam0), "lambda_0(H_2[v])=0 (arbitrary v)")

# residual: lambda_0(R)=0 on the displayed d=2 residual, lambda_0(E-R)=4
a1_2, a2_0, a2_2 = sp.symbols("a1_2 a2_0 a2_2")
Rd2 = E * a1_2**2 / 18 * (2 * E**3 * a2_0 - 12 * E**3 * a2_2 - 15 * E**2 * a2_0
                          + 92 * E**2 * a2_2 + 20 * E * a2_0 - 140 * E * a2_2
                          + 15 * a2_0 - 44 * a2_2)
check_zero(apply_functional(Rd2, lam0), "lambda_0(R_{d=2})=0 (residual on solution variety)")
check_zero(apply_functional(sp.expand(E - Rd2), lam0) - 4,
           "lambda_0(E-R)=4 != 0  => E-R not in Im Phi (W1 obstruction)")

# ===========================================================================
print("\n=== 4. AP band-3 family: symbolic r, arbitrary degree ===")
r = sp.symbols("r")
a3r = sp.expand((E - r) * (E - r - 2) * (E - r - 4))
b2r = sp.expand((E - r - 1) * (E - r - 4))
Cp = sum(sp.symbols("s0:4")[i] * E**i for i in range(4))
Vp = sum(sp.symbols("t0:4")[i] * E**i for i in range(4))
cR = sp.expand(falling(3) * Cp)
vR = sp.expand(falling(2) * Vp)
KR = Kk(a3r, cR, 3)
HR = Hk(b2r, vR, 3)


def lam_r(P):
    return sp.expand(P.subs(E, r + 3) - P.subs(E, r + 4) + P.subs(E, r + 5) - P.subs(E, 0))


check_zero(lam_r(KR), "lambda_r(K_3[c])=0 for symbolic r, arbitrary c")
check_zero(lam_r(HR), "lambda_r(H_2[v])=0 for symbolic r, arbitrary v")
check_zero(lam_r(E) - (r + 4), "lambda_r(E)=r+4")
check_zero(a3r.subs(E, r + 4), "a3(r+4)=0  (common root is r+4)")
check_zero(b2r.subs(E, r + 4), "b2(r+4)=0  (common root is r+4)")
check_zero((r + 4).subs(r, -4), "escape hatch: common root r+4 hits 0 exactly at r=-4")
# codimension 2 off r=-4, 3 at r=-4
for rv in (-3, -2, -1, 0, 1, 2, 3):
    a = sp.expand(a3r.subs(r, rv))
    u = sp.expand(b2r.subs(r, rv))
    res = annihilator(a, u, 3, window=max(14, abs(rv) + 12))
    if res is None:
        continue
    ns, var = res
    expected = 3 if rv == -4 else 2
    check(len(ns) == expected, f"AP r={rv}: codim = {len(ns)} (expected {expected})")
res = annihilator(sp.expand(a3r.subs(r, -4)), sp.expand(b2r.subs(r, -4)), 3, window=14)
check(res is not None and len(res[0]) == 3, "AP r=-4: codim jumps to 3 (Im Phi = E(E-1)(E+1)F[E])")

# ===========================================================================
print("\n=== 5. dim Ann(Im Phi) equals the stable finite cokernel dimension ===")
for name, a, u, k, expect in [
    ("W1", sp.expand(E * (E - 2) * (E - 4)), sp.expand((E - 1) * (E - 4)), 3, 2),
    ("AP r=1", sp.expand((E - 1) * (E - 3) * (E - 5)), sp.expand((E - 2) * (E - 5)), 3, 2),
    ("AP r=-4", sp.expand((E + 4) * (E + 2) * E), sp.expand((E + 3) * E), 3, 3),
]:
    ad = len(annihilator(a, u, k, window=16)[0])
    codims = [true_codim(a, u, k, d) for d in (3, 4, 5)]
    check(ad == expect and codims[-1] == expect and codims[-2] == expect,
          f"{name}: dim Ann={ad} = stable finite codim {codims} (converges to {expect})")

# ===========================================================================
print("\n=== 6. bands 4 and 5: codim>=2, explicit kill, lambda(E) != k+1 ===")
band_tops = [
    (4, sp.expand(E * (E - 2) * (E - 3) * (E - 5)),
        sp.expand((E - 1) * (E - 3) * (E - 5)), [3, 5]),        # common roots 3,5
    (5, sp.expand(E * (E - 2) * (E - 4) * (E - 6) * (E - 8)),
        sp.expand((E - 1) * (E - 3) * (E - 6) * (E - 8)), [6, 8]),
]
for k, a, u, common in band_tops:
    # verify wall Q_{2k-1}=0: u(E+k)a(E)=a(E+k-1)u(E)
    check_zero(sh(u, k) * a - sh(a, k - 1) * u, f"k={k}: top wall Q_(2k-1)=0 holds")
    res = annihilator(a, u, k, window=8 * k + 6)
    check(res is not None, f"k={k}: q_K,q_H fully integer-rooted (tame)")
    ns, var = res
    check(len(ns) >= 2, f"k={k}: codim Im Phi = {len(ns)} >= 2 (method has ammunition)")
    # common roots are shared by a and u (wall signature)
    for cr in common:
        check_zero(a.subs(E, cr), f"k={k}: a({cr})=0 (common root)")
        check_zero(u.subs(E, cr), f"k={k}: u({cr})=0 (common root)")
    # explicit annihilator kills Im Phi on generic admissible C,V
    Cp = sum(sp.symbols("m0:3")[i] * E**i for i in range(3))
    Vp = sum(sp.symbols("w0:3")[i] * E**i for i in range(3))
    c = sp.expand(falling(k) * Cp)
    v = sp.expand(falling(k - 1) * Vp)
    Phi = sp.expand(Kk(a, c, k) - Hk(u, v, k))
    for bvec in ns:
        coeffs = vec_to_coeffs(bvec, var)
        check_zero(apply_functional(Phi, coeffs),
                   f"k={k}: annihilator basis functional kills Im Phi (generic C,V)")
    # lambda(E) values over the basis are NOT equal to k+1 -> the W1 '4' is not k+1
    lamvals = sorted({lamE(bvec, var) for bvec in ns})
    check((k + 1) not in lamvals,
          f"k={k}: no basis functional has lambda(E)=k+1={k+1}; values={lamvals} "
          f"(the W1 value 4 = common root, coincidentally k+1 only at k=3)")

# ===========================================================================
print(f"\n{PASS} checks passed.")
print("ALL LAMBDA GENERAL K CHECKS PASSED")
