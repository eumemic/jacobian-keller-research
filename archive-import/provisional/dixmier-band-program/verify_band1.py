# STATUS: PROVISIONAL RESEARCH ARCHIVE.
# This checker passed in the archive review environment, but passing checks do not
# independently validate the associated mathematical theorem or literature claims.

#!/usr/bin/env python3
"""
verify_prong2.py -- machine verification for Milestone 2 of the DC1 program.

Implements the crossed-product presentation A1[x^-1] = (+)_k x^k C[E], E = x d/dx,
with multiplication (x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E), and checks:

  1. Engine sanity: [d,x]=1, d^3 = x^-3 E(E-1)(E-2)  (falling-factorial membership).
  2. Exact quantization of the tower bottom layer:
       [C(N)+2NW, A(N)+N^2 W] = 2 N^2 W   with A'=nu*E(nu), C'=2E(nu),
     for a fully symbolic E(nu) of degree 3 -- no quantum corrections.
  3. The classical dictionary: {x^l b(tau), x^k a(tau)} = x^{k+l}(k a b' - l a' b),
     tau = x*xi (so the quantum system is a difference deformation of the classical).
  4. Band-1 rigidity steps (Theorem P3):
     (i)  b1(E+1)a1(E) = a1(E+1)b1(E) forces b1 proportional to a1 (nullspace dim 1);
     (ii) with b_k = lambda/mu * a_k and constant middle terms, the degree-0 bracket
          equals (lambda-mu)(V(E)-V(E+1)), V(E)=a1(E-1)a_{-1}(E)  [telescoping],
          so V is linear and the classification follows;
     (iii) instances: an affine symplectic pair has bracket 1 and lies in A1;
           the polar pair (x, d + 4/x) has bracket 1 and fails membership.
Expected output: all PASS, then "ALL PRONG-2 CHECKS PASSED".
"""
import warnings; warnings.filterwarnings("ignore")
import sympy as sp

E, nu, lam, mu, alpha, beta = sp.symbols("E nu lambda mu alpha beta")
RES = []
def check(name, ok):
    RES.append((name, bool(ok)))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

# ---------------- crossed-product engine: element = {k: expr in E} ----------------
def clean(A):
    return {k: sp.expand(v) for k, v in A.items() if sp.expand(v) != 0}

def mul(A, B):
    C = {}
    for a, f in A.items():
        for b, g in B.items():
            C[a + b] = C.get(a + b, 0) + f.subs(E, E + b) * g
    return clean(C)

def add(A, B):
    C = dict(A)
    for k, v in B.items():
        C[k] = C.get(k, 0) + v
    return clean(C)

def neg(A): return {k: -v for k, v in A.items()}
def brk(A, B): return add(mul(A, B), neg(mul(B, A)))
def scal(c, A): return clean({k: c * v for k, v in A.items()})

X1 = {1: sp.Integer(1)}          # x
D1 = {-1: E}                     # d = x^{-1} E

def in_A1(A):
    """x^{-j} c(E) lies in A1 iff E(E-1)...(E-j+1) divides c."""
    for k, v in A.items():
        if k < 0:
            ff = sp.prod([E - i for i in range(-k)])
            if sp.rem(sp.Poly(v, E), sp.Poly(ff, E)).as_expr() != 0:
                return False
    return True

# ---------------- Section 1: engine sanity ----------------
print("Section 1: crossed-product engine")
check("[d, x] = 1", brk(D1, X1) == {0: sp.Integer(1)})
d3 = mul(mul(D1, D1), D1)
check("d^3 = x^-3 E(E-1)(E-2)  (falling factorial)",
      d3 == {-3: sp.expand(E * (E - 1) * (E - 2))})
check("membership: d^2 in A1, but x^-1(E+4) not in A1",
      in_A1(mul(D1, D1)) and not in_A1({-1: E + 4}))

# ---------------- Section 2: exact quantization of the bottom layer ----------------
print("Section 2: the tower bottom quantizes exactly")
e0, e1, e2, e3 = sp.symbols("e0 e1 e2 e3")
Epol = e0 + e1 * nu + e2 * nu**2 + e3 * nu**3
Apol = sp.integrate(nu * Epol, nu)          # A' = nu E
Cpol = 2 * sp.integrate(Epol, nu)           # C' = 2 E
def poly_of_N2(p):
    return clean({k[0]: c for k, c in sp.Poly(p, nu).as_dict().items()})
N2W = mul(mul(X1, X1), D1)                   # N^2 W  (W to the right)
NW  = mul(X1, D1)
S_t = add(poly_of_N2(Apol), N2W)             # S~ = A(N) + N^2 W
T_t = add(poly_of_N2(Cpol), scal(2, NW))     # T~ = C(N) + 2 N W
check("[T~, S~] = 2 N^2 W identically in e0..e3 (no quantum corrections)",
      brk(T_t, S_t) == scal(2, N2W))

# ---------------- Section 3: the classical dictionary ----------------
print("Section 3: classical limit of the graded bracket")
x, xi, tau = sp.symbols("x xi tau")
a0c, a1c, a2c, b0c, b1c, b2c = sp.symbols("a0c a1c a2c b0c b1c b2c")
aC = a0c + a1c * tau + a2c * tau**2
bC = b0c + b1c * tau + b2c * tau**2
ok = True
for (k, l) in [(2, -1), (1, 1), (-1, 3), (-2, -1)]:
    f = x**l * bC.subs(tau, x * xi)
    g = x**k * aC.subs(tau, x * xi)
    pb = sp.expand(sp.diff(f, xi) * sp.diff(g, x) - sp.diff(f, x) * sp.diff(g, xi))
    rhs = sp.expand((x**(k + l) * (k * aC * sp.diff(bC, tau)
                                   - l * sp.diff(aC, tau) * bC)).subs(tau, x * xi))
    ok &= sp.simplify(pb - rhs) == 0
check("{x^l b, x^k a} = x^{k+l}(k a b' - l a' b) for four (k,l) pairs", ok)

# ---------------- Section 4: band-1 rigidity (Theorem P3) ----------------
print("Section 4: band-1 rigidity steps")
# (i) proportionality from the m=+-2 equations, generic a1 of degree 3
a1gen = sp.Rational(3, 7) + 2 * E - sp.Rational(5, 2) * E**2 + E**3
c0, c1, c2, c3, c4 = sp.symbols("c0:5")
b1gen = c0 + c1 * E + c2 * E**2 + c3 * E**3 + c4 * E**4
eqn = sp.expand(b1gen.subs(E, E + 1) * a1gen - a1gen.subs(E, E + 1) * b1gen)
sols = sp.solve([sp.Poly(eqn, E).coeff_monomial(E**i) for i in range(9)],
                [c0, c1, c2, c3, c4], dict=True)
sol_space_ok = False
if sols:
    b1sol = sp.expand(b1gen.subs(sols[0]))
    free = [s for s in (c0, c1, c2, c3, c4) if s not in sols[0]]
    ratio = sp.simplify(sp.cancel(b1sol / a1gen))
    sol_space_ok = (len(free) == 1 and ratio.free_symbols <= set(free))
check("m=2 equation forces b1 = lambda * a1 (solution space is the line C*a1)",
      sol_space_ok)
# (ii) telescoping identity for the m=0 component, fully symbolic a1, a_{-1} of deg 2
p0, p1, p2, q0, q1, q2 = sp.symbols("p0 p1 p2 q0 q1 q2")
a1s = p0 + p1 * E + p2 * E**2
am1s = q0 + q1 * E + q2 * E**2
Xel = clean({1: a1s, 0: alpha, -1: am1s})
Del = clean({1: lam * a1s, 0: beta, -1: mu * am1s})
B = brk(Del, Xel)
V = a1s.subs(E, E - 1) * am1s
tele = sp.expand((lam - mu) * (V - V.subs(E, E + 1)))
check("degree-0 bracket component == (lambda-mu)(V(E)-V(E+1)), V = a1(E-1)a_{-1}(E)",
      sp.expand(B.get(0, 0) - tele) == 0 and
      all(B.get(m, 0) == 0 for m in (1, -1, 2, -2)))
# (iii) instances
Xaff = clean({1: sp.Integer(2), 0: sp.Integer(5), -1: -3 * E})   # 2x + 5 - 3d
Daff = clean({1: sp.Integer(1), 0: sp.Integer(7), -1: -E})       # x + 7 - d  (da-cb=1)
check("affine symplectic pair: bracket 1 and lies in A1",
      brk(Daff, Xaff) == {0: sp.Integer(1)} and in_A1(Xaff) and in_A1(Daff))
Dpol = {-1: E + 4}                                               # d + 4/x
check("polar pair (x, d+4/x): bracket 1 in A1[x^-1] but fails A1-membership",
      brk(Dpol, X1) == {0: sp.Integer(1)} and not in_A1(Dpol))

fails = [n for n, ok in RES if not ok]
print()
print("ALL PRONG-2 CHECKS PASSED" if not fails else f"FAILURES: {fails}")
import sys; sys.exit(0 if not fails else 1)
