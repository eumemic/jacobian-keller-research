#!/usr/bin/env python3
"""Independent exact checker for the explicit A_3 Dixmier-counterexample witness.

Certifies the formulas in
  archive-import/provisional/weyl-a3/endomorphism-formulas.md
by recomputing everything from the map F alone, in pure Python with exact
Fraction/integer arithmetic (no sympy, no floating point, no external
dependencies).  Run with:  python3 verify_a3_witness.py

Checks performed (all exact):
  1.  The archived phi(x_i) formulas equal F_i.
  2.  det JF = -2 identically.
  3.  JF^{-1} := -adj(JF)/2 satisfies JF*JF^{-1} = JF^{-1}*JF = I,
      and det JF^{-1} = -1/2.
  4.  All 9 archived phi(d_j) coefficient polynomials equal (JF^{-1})_{kj},
      entry by entry.
  5.  The vector fields Y_j = sum_k (JF^{-1})_{kj} d/dx_k pairwise commute
      (vector-field bracket identities, coefficientwise).
  6.  Self-test of an independent normal-ordered Weyl-algebra implementation
      (defining relations, a Leibniz identity, associativity instances).
  7.  In that Weyl algebra: [phi(d_i), phi(x_j)] = delta_ij,
      [phi(d_i), phi(d_j)] = 0, [phi(x_i), phi(x_j)] = 0, as full operator
      identities (this includes the direct Weyl-algebra confirmation of the
      commutativity checked as vector fields in 5).
  8.  The fiber of F over (-1/4, 0, 0) contains the three distinct points
      (0,0,-1/4), (1,-3/2,13/2), (-1,3/2,13/2); hence F is not injective.
  9.  F is the telescoping-family member (E,h) = (2-3v, 2-3t) of the
      in-repo paper (paper/main.tex, Theorem "Telescoping family"), i.e.
      S = x^2 g^2 F_1, T = x g F_2, F_3 = x g as polynomial identities,
      and the three scalar polynomiality conditions hold at a = h(0) = 2.

A successful run ends with the line:  ALL A3 WITNESS CHECKS PASSED
"""

from fractions import Fraction
from math import comb
import sys

# ----------------------------------------------------------------------
# Exact multivariate polynomials over Q in x, y, z.
# Representation: dict mapping exponent triples (i,j,k) -> nonzero Fraction.
# ----------------------------------------------------------------------

class Poly:
    __slots__ = ("c",)

    def __init__(self):
        self.c = {}

    @staticmethod
    def const(v):
        p = Poly()
        v = Fraction(v)
        if v:
            p.c[(0, 0, 0)] = v
        return p

    @staticmethod
    def var(i):
        p = Poly()
        e = [0, 0, 0]
        e[i] = 1
        p.c[tuple(e)] = Fraction(1)
        return p

    @staticmethod
    def _coerce(o):
        if isinstance(o, Poly):
            return o
        if isinstance(o, (int, Fraction)):
            return Poly.const(o)
        return NotImplemented

    def __add__(self, o):
        o = Poly._coerce(o)
        if o is NotImplemented:
            return NotImplemented
        r = Poly()
        r.c = dict(self.c)
        for k, v in o.c.items():
            nv = r.c.get(k, Fraction(0)) + v
            if nv:
                r.c[k] = nv
            elif k in r.c:
                del r.c[k]
        return r

    __radd__ = __add__

    def __neg__(self):
        r = Poly()
        r.c = {k: -v for k, v in self.c.items()}
        return r

    def __sub__(self, o):
        o = Poly._coerce(o)
        if o is NotImplemented:
            return NotImplemented
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        o = Poly._coerce(o)
        if o is NotImplemented:
            return NotImplemented
        r = Poly()
        for k1, v1 in self.c.items():
            for k2, v2 in o.c.items():
                k = (k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2])
                nv = r.c.get(k, Fraction(0)) + v1 * v2
                if nv:
                    r.c[k] = nv
                elif k in r.c:
                    del r.c[k]
        return r

    __rmul__ = __mul__

    def __truediv__(self, n):
        return self * Fraction(1, n)

    def __pow__(self, n):
        r = Poly.const(1)
        for _ in range(n):
            r = r * self
        return r

    def diff(self, i):
        r = Poly()
        for k, v in self.c.items():
            if k[i]:
                kk = list(k)
                kk[i] -= 1
                r.c[tuple(kk)] = v * k[i]
        return r

    def eval_at(self, pt):
        # Exact evaluation of the polynomial at a rational point (Fraction
        # arithmetic only; unrelated to the Python builtin eval).
        s = Fraction(0)
        for k, v in self.c.items():
            t = v
            for i in range(3):
                t *= Fraction(pt[i]) ** k[i]
            s += t
        return s

    def is_zero(self):
        return not self.c

    def __eq__(self, o):
        o = Poly._coerce(o)
        if o is NotImplemented:
            return NotImplemented
        return self.c == o.c

    def __hash__(self):
        return hash(frozenset(self.c.items()))


x, y, z = Poly.var(0), Poly.var(1), Poly.var(2)

# ----------------------------------------------------------------------
# Normal-ordered Weyl algebra A_3.
# Element: dict mapping d-multi-index (a1,a2,a3) -> Poly coefficient,
# representing  sum_alpha  p_alpha(x,y,z) * d1^a1 d2^a2 d3^a3.
# Multiplication via the Leibniz rule:
#   d^alpha q = sum_{gamma <= alpha} binom(alpha,gamma) (D^gamma q) d^{alpha-gamma}
# ----------------------------------------------------------------------

def w_clean(A):
    return {k: v for k, v in A.items() if not v.is_zero()}


def w_from_poly(p):
    return w_clean({(0, 0, 0): p})


def w_d(j):
    e = [0, 0, 0]
    e[j] = 1
    return {tuple(e): Poly.const(1)}


def w_add(A, B):
    r = dict(A)
    for k, v in B.items():
        r[k] = r.get(k, Poly()) + v
    return w_clean(r)


def w_neg(A):
    return {k: -v for k, v in A.items()}


def w_sub(A, B):
    return w_add(A, w_neg(B))


def w_mul(A, B):
    r = {}
    for al, p in A.items():
        for be, q in B.items():
            for g0 in range(al[0] + 1):
                dq0 = q
                for _ in range(g0):
                    dq0 = dq0.diff(0)
                for g1 in range(al[1] + 1):
                    dq1 = dq0
                    for _ in range(g1):
                        dq1 = dq1.diff(1)
                    for g2 in range(al[2] + 1):
                        dq2 = dq1
                        for _ in range(g2):
                            dq2 = dq2.diff(2)
                        coef = comb(al[0], g0) * comb(al[1], g1) * comb(al[2], g2)
                        idx = (al[0] - g0 + be[0],
                               al[1] - g1 + be[1],
                               al[2] - g2 + be[2])
                        term = p * dq2 * coef
                        r[idx] = r.get(idx, Poly()) + term
    return w_clean(r)


def w_comm(A, B):
    return w_sub(w_mul(A, B), w_mul(B, A))


def w_is_zero(A):
    return not w_clean(A)


def w_eq(A, B):
    return w_is_zero(w_sub(A, B))


# ----------------------------------------------------------------------
# Check harness
# ----------------------------------------------------------------------

FAILED = []


def check(name, ok):
    print(("PASS  " if ok else "FAIL  ") + name)
    if not ok:
        FAILED.append(name)


# ----------------------------------------------------------------------
# The map F (recomputation source: these three lines only).
# ----------------------------------------------------------------------

F1 = (1 + x * y) ** 3 * z + y ** 2 * (1 + x * y) * (4 + 3 * x * y)
F2 = y + 3 * x * (1 + x * y) ** 2 * z + 3 * x * y ** 2 * (4 + 3 * x * y)
F3 = 2 * x - 3 * x ** 2 * y - x ** 3 * z
F = [F1, F2, F3]

# --- 1. archived phi(x_i) formulas ------------------------------------

arch_phix = [
    (x * y + 1) * (x ** 2 * y ** 2 * z + 3 * x * y ** 3 + 2 * x * y * z + 4 * y ** 2 + z),
    3 * x ** 3 * y ** 2 * z + 9 * x ** 2 * y ** 3 + 6 * x ** 2 * y * z + 12 * x * y ** 2 + 3 * x * z + y,
    -x * (x ** 2 * z + 3 * x * y - 2),
]
for i in range(3):
    check(f"archived phi(x{i+1}) equals F{i+1}", F[i] == arch_phix[i])

# --- 2. Jacobian and determinant --------------------------------------

JF = [[F[i].diff(j) for j in range(3)] for i in range(3)]


def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


detJF = det3(JF)
check("det JF = -2 identically", detJF == Poly.const(-2))

# --- 3. inverse Jacobian ----------------------------------------------


def adjugate3(M):
    C = [[None] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            r = [k for k in range(3) if k != i]
            c = [k for k in range(3) if k != j]
            minor = M[r[0]][c[0]] * M[r[1]][c[1]] - M[r[0]][c[1]] * M[r[1]][c[0]]
            C[j][i] = minor * ((-1) ** (i + j))   # adjugate = transpose of cofactors
    return C


A = adjugate3(JF)
JFinv = [[A[i][j] * Fraction(-1, 2) for j in range(3)] for i in range(3)]


def mat_mul(M, N):
    return [[sum((M[i][k] * N[k][j] for k in range(3)), Poly())
             for j in range(3)] for i in range(3)]


def is_identity(M):
    return all(M[i][j] == (1 if i == j else 0) for i in range(3) for j in range(3))


check("JF * JF^{-1} = I", is_identity(mat_mul(JF, JFinv)))
check("JF^{-1} * JF = I", is_identity(mat_mul(JFinv, JF)))
check("det JF^{-1} = -1/2", det3(JFinv) == Poly.const(Fraction(-1, 2)))

# --- 4. archived phi(d_j) coefficients vs (JF^{-1})_{kj} ---------------
# arch_a[k][j] = archived coefficient of d_{k+1} in phi(d_{j+1}).

arch_a = [[None] * 3 for _ in range(3)]
arch_a[0][0] = x ** 3 * (3 * x ** 3 * y * z + 9 * x ** 2 * y ** 2 + 3 * x ** 2 * z + 3 * x * y - 4)
arch_a[1][0] = 3 * x * (x ** 3 * y * z + 3 * x ** 2 * y ** 2 + x ** 2 * z + x * y - 1)
arch_a[2][0] = (-9 * x ** 5 * y * z ** 2 - 45 * x ** 4 * y ** 2 * z - 9 * x ** 4 * z ** 2
                - 54 * x ** 3 * y ** 3 - 30 * x ** 3 * y * z - 27 * x ** 2 * y ** 2
                + 9 * x ** 2 * z + 21 * x * y + 1)
arch_a[0][1] = -x ** 2 * (3 * x ** 4 * y ** 2 * z + 9 * x ** 3 * y ** 3 + 6 * x ** 3 * y * z
                          + 12 * x ** 2 * y ** 2 + 3 * x ** 2 * z - x * y - 3) / 2
arch_a[1][1] = -(3 * x ** 4 * y ** 2 * z + 9 * x ** 3 * y ** 3 + 6 * x ** 3 * y * z
                 + 12 * x ** 2 * y ** 2 + 3 * x ** 2 * z - 2) / 2
arch_a[2][1] = (9 * x ** 5 * y ** 2 * z ** 2 + 45 * x ** 4 * y ** 3 * z + 18 * x ** 4 * y * z ** 2
                + 54 * x ** 3 * y ** 4 + 75 * x ** 3 * y ** 2 * z + 9 * x ** 3 * z ** 2
                + 81 * x ** 2 * y ** 3 + 21 * x ** 2 * y * z + 6 * x * y ** 2
                - 6 * x * z - 16 * y) / 2
arch_a[0][2] = -(x * y + 1) ** 2 * (3 * x ** 4 * y ** 2 * z + 9 * x ** 3 * y ** 3
                                    + 6 * x ** 3 * y * z + 12 * x ** 2 * y ** 2
                                    + 3 * x ** 2 * z - x * y - 1) / 2
arch_a[1][2] = -3 * (x * y + 1) ** 2 * (x ** 2 * y ** 2 * z + 3 * x * y ** 3
                                        + 2 * x * y * z + 4 * y ** 2 + z) / 2
arch_a[2][2] = (9 * x ** 5 * y ** 4 * z ** 2 + 45 * x ** 4 * y ** 5 * z + 36 * x ** 4 * y ** 3 * z ** 2
                + 54 * x ** 3 * y ** 6 + 165 * x ** 3 * y ** 4 * z + 54 * x ** 3 * y ** 2 * z ** 2
                + 189 * x ** 2 * y ** 5 + 216 * x ** 2 * y ** 3 * z + 36 * x ** 2 * y * z ** 2
                + 222 * x * y ** 4 + 117 * x * y ** 2 * z + 9 * x * z ** 2
                + 89 * y ** 3 + 21 * y * z) / 2

arch_match = True
for k in range(3):
    for j in range(3):
        if not (JFinv[k][j] == arch_a[k][j]):
            arch_match = False
            print(f"      mismatch at (JF^-1)[{k+1}][{j+1}]")
check("all 9 archived phi(d_j) coefficients equal (JF^{-1})_{kj}", arch_match)

# --- 5. vector-field brackets [Y_i, Y_j] = 0 ---------------------------
# Y_j = sum_k JFinv[k][j] d/dx_k ;  [Y_i,Y_j] = sum_l (Y_i(a_lj) - Y_j(a_li)) d/dx_l.

vf_ok = True
for i in range(3):
    for j in range(i + 1, 3):
        for l in range(3):
            s = Poly()
            for k in range(3):
                s = s + JFinv[k][i] * JFinv[l][j].diff(k) \
                      - JFinv[k][j] * JFinv[l][i].diff(k)
            if not s.is_zero():
                vf_ok = False
                print(f"      nonzero bracket coefficient: [Y{i+1},Y{j+1}], d/dx{l+1}")
check("vector fields Y_1, Y_2, Y_3 pairwise commute (bracket identities)", vf_ok)

# --- 6. Weyl-algebra implementation self-test --------------------------

Wx = [w_from_poly(Poly.var(i)) for i in range(3)]
Wd = [w_d(i) for i in range(3)]

st = all(w_eq(w_comm(Wd[i], Wx[j]), w_from_poly(Poly.const(1 if i == j else 0)))
         for i in range(3) for j in range(3))
st = st and all(w_is_zero(w_comm(Wd[i], Wd[j])) for i in range(3) for j in range(3))
st = st and all(w_is_zero(w_comm(Wx[i], Wx[j])) for i in range(3) for j in range(3))
check("Weyl self-test: defining relations [d_i,x_j]=delta_ij, others 0", st)

# d1 * x^2 = x^2 d1 + 2x   and   d1^2 * x = x d1^2 + 2 d1
lhs = w_mul(w_d(0), w_from_poly(x ** 2))
rhs = w_add({(1, 0, 0): x ** 2}, w_from_poly(2 * x))
ok_leib = w_eq(lhs, rhs)
lhs = w_mul({(2, 0, 0): Poly.const(1)}, w_from_poly(x))
rhs = w_add({(2, 0, 0): x}, {(1, 0, 0): Poly.const(2)})
ok_leib = ok_leib and w_eq(lhs, rhs)
check("Weyl self-test: Leibniz normal-ordering identities", ok_leib)

# associativity on nontrivial elements involving second-order terms
E1 = w_add({(1, 0, 1): x * y + 1}, w_from_poly(z ** 2 + 3))
E2 = w_add({(0, 2, 0): x ** 2 - y}, {(1, 0, 0): z})
E3 = w_add({(0, 1, 1): y * z}, w_from_poly(x - 5))
assoc = w_eq(w_mul(w_mul(E1, E2), E3), w_mul(E1, w_mul(E2, E3)))
check("Weyl self-test: associativity instance with order-2 operators", assoc)

# --- 7. the endomorphism relations in the Weyl algebra -----------------

phi_x = [w_from_poly(F[i]) for i in range(3)]
phi_d = [w_clean({(1, 0, 0): JFinv[0][j], (0, 1, 0): JFinv[1][j], (0, 0, 1): JFinv[2][j]})
         for j in range(3)]

ok_dx = all(w_eq(w_comm(phi_d[i], phi_x[j]),
                 w_from_poly(Poly.const(1 if i == j else 0)))
            for i in range(3) for j in range(3))
check("[phi(d_i), phi(x_j)] = delta_ij  (9 operator identities)", ok_dx)

ok_dd = all(w_is_zero(w_comm(phi_d[i], phi_d[j]))
            for i in range(3) for j in range(3))
check("[phi(d_i), phi(d_j)] = 0  (direct Weyl-algebra computation)", ok_dd)

ok_xx = all(w_is_zero(w_comm(phi_x[i], phi_x[j]))
            for i in range(3) for j in range(3))
check("[phi(x_i), phi(x_j)] = 0", ok_xx)

# --- 8. the triple fiber over (-1/4, 0, 0) -----------------------------

pts = [(Fraction(0), Fraction(0), Fraction(-1, 4)),
       (Fraction(1), Fraction(-3, 2), Fraction(13, 2)),
       (Fraction(-1), Fraction(3, 2), Fraction(13, 2))]
target = (Fraction(-1, 4), Fraction(0), Fraction(0))
fib_ok = all(tuple(F[i].eval_at(p) for i in range(3)) == target for p in pts)
fib_ok = fib_ok and len(set(pts)) == 3
check("F maps the three distinct points to (-1/4, 0, 0): F not injective", fib_ok)

# first coordinates of the fiber points are pairwise distinct, so x1 is not
# constant on a fiber of F; the same holds for the second coordinates.
sep_ok = (len({p[0] for p in pts}) == 3 and len({p[1] for p in pts}) == 3
          and len({p[2] for p in pts}) == 2)
check("fiber points separated in each coordinate (x1,x2 fully; x3 partially)", sep_ok)

# --- 9. telescoping-family membership (paper/main.tex, Thm 'family') ---
# member (E,h) = (2 - 3v, 2 - 3t):  A(v) = v^2 - v^3, C(v) = 4v - 3v^2,
# g = h(t) - s, u = 1 + t, v = ug, S = A(ug) + u g^2, T = C(ug) + 2g.

t_, s_ = x * y, x ** 2 * z
u_ = 1 + t_
g_ = 2 - 3 * t_ - s_
nu = u_ * g_
Anu = nu ** 2 - nu ** 3
Cnu = 4 * nu - 3 * nu ** 2
S_ = Anu + u_ * g_ ** 2
T_ = Cnu + 2 * g_
fam = (S_ == x ** 2 * g_ ** 2 * F1) and (T_ == x * g_ * F2) and (F3 == x * g_)
check("F equals the telescoping member (E,h) = (2-3v, 2-3t)", fam)

# scalar conditions at a = h(0) = 2 with E(v) = 2 - 3v, h'(0) = -3:
a0 = Fraction(2)
A_a = a0 ** 2 - a0 ** 3          # A(a)
C_a = 4 * a0 - 3 * a0 ** 2       # C(a)
E_a = 2 - 3 * a0                 # E(a)
cond = (A_a == -a0 ** 2) and (C_a == -2 * a0) \
       and (1 + E_a + Fraction(-3) * (E_a + 2) / a0 == 0)
check("family conditions A(a)=-a^2, C(a)=-2a, third identity, at a=2", cond)

# ----------------------------------------------------------------------

print()
if FAILED:
    print(f"{len(FAILED)} CHECK(S) FAILED:")
    for name in FAILED:
        print("  - " + name)
    sys.exit(1)
print("ALL A3 WITNESS CHECKS PASSED")
