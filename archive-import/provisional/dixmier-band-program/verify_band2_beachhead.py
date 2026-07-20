# STATUS: PROVISIONAL RESEARCH ARCHIVE.
# This checker passed in the archive review environment, but passing checks do not
# independently validate the associated mathematical theorem or literature claims.

#!/usr/bin/env python3
"""
verify_band2_beachhead.py -- machine verification for the JC2 attack's first strike.

CLASSICAL SIDE (the plane C^2_{x,xi}, tau = x*xi, Poisson bracket {f,g} = f_x g_xi - f_xi g_x):
  J1 (band-1 rigidity, classical): pairs {D,X} = const with Laurent support in
     x-degrees {-1,0,1} along tau-levels collapse exactly as in the quantum P3:
     (i)   m=+-2 Wronskian equations force b1 || a1  (nullspace = C*a1);
     (ii)  the m=0 component telescopes to a multiple of (a1*a_{-1})' so a1*a_{-1}
           is linear, forcing the affine-symplectic classification;
     (iii) instances: affine symplectic pair has constant bracket and is polynomial;
           the polar pair (x, xi + c/x) has constant bracket but is not polynomial.
  J2 (band-2 bifurcation, classical): in the m=3 equation, b1 = lam2*a1 is always a
     particular solution, and the homogeneous solutions h satisfy h^2 = c*a2:
     they exist over C[tau] iff a2 is a perfect square (up to scalar).

QUANTUM SIDE (crossed product A1[x^-1] = (+) x^k C[E]):
  J2q: the m=3 homogeneous equation h(E+2)a2(E) = a2(E+1)h(E) has polynomial
     solutions iff a2(E) = c*h(E)h(E+1): the falling-factorial deformation of
     "a2 is a square". The m=4 equation still forces b2 || a2.

Expected output: all PASS, then "ALL BEACHHEAD CHECKS PASSED".
"""
import warnings; warnings.filterwarnings("ignore")
import sympy as sp

x, xi, tau, E = sp.symbols("x xi tau E")
lam, mu, lam2, alpha, beta, kappa = sp.symbols("lambda mu lambda2 alpha beta kappa")
RES = []
def check(name, ok):
    RES.append((name, bool(ok)))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

def PB(f, g):
    return sp.expand(sp.diff(f, x) * sp.diff(g, xi) - sp.diff(f, xi) * sp.diff(g, x))

def x_components(expr, kmin, kmax):
    """Decompose a Laurent expression in x with coefficients in tau = x*xi."""
    comp = {}
    e = sp.expand(expr.subs(xi, tau / x))
    for k in range(kmin, kmax + 1):
        c = sp.expand(e.coeff(x, k))
        if c != 0:
            comp[k] = sp.expand(c)
    resid = sp.expand(e - sum(x**k * v for k, v in comp.items()))
    assert resid == 0, "unexpected x-degrees"
    return comp

# ------------------ J1 classical band-1 ------------------
print("Section 1: classical band-1 rigidity (Theorem J1)")
# (i) Wronskian: generic a1 of degree 3, solve a1*b1' - a1'*b1 = 0 for b1 of degree <= 4
a1g = sp.Rational(3, 7) + 2 * tau - sp.Rational(5, 2) * tau**2 + tau**3
c0, c1, c2, c3, c4 = sp.symbols("c0:5")
b1g = c0 + c1 * tau + c2 * tau**2 + c3 * tau**3 + c4 * tau**4
w = sp.expand(a1g * sp.diff(b1g, tau) - sp.diff(a1g, tau) * b1g)
sols = sp.solve([sp.Poly(w, tau).coeff_monomial(tau**i) for i in range(8)],
                [c0, c1, c2, c3, c4], dict=True)
ok = False
if sols:
    b1sol = sp.expand(b1g.subs(sols[0]))
    free = [s for s in (c0, c1, c2, c3, c4) if s not in sols[0]]
    ratio = sp.cancel(b1sol / a1g)
    ok = len(free) == 1 and ratio.free_symbols <= set(free)
check("m=2 Wronskian forces b1 = lambda*a1 (solution space is C*a1)", ok)

# (ii) telescoping of the m=0 component, fully symbolic coefficients of degree 2
p0, p1, p2, q0, q1, q2 = sp.symbols("p0 p1 p2 q0 q1 q2")
a1s = (p0 + p1 * tau + p2 * tau**2)
am1s = (q0 + q1 * tau + q2 * tau**2)
Xc = a1s.subs(tau, x * xi) * x + alpha + am1s.subs(tau, x * xi) / x
Dc = lam * a1s.subs(tau, x * xi) * x + beta + mu * am1s.subs(tau, x * xi) / x
comp = x_components(PB(Dc, Xc), -4, 4)
V = sp.expand(a1s * am1s)
tele = sp.expand((lam - mu) * sp.diff(V, tau))
check("m=0 component == (lambda-mu) d/dtau (a1*a_{-1}); components m=+-1,+-2 vanish",
      sp.expand(comp.get(0, 0) - tele) == 0 and
      all(comp.get(m, 0) == 0 for m in (1, -1, 2, -2)))

# (iii) instances
Xa = 2 * x - 3 * xi + 5
Da = -x + xi + 7                      # bc - ad = (-3)(-1) - (2)(1) = 1 in this convention
pol = xi + 4 / x
check("affine symplectic pair: {D,X} = 1 and polynomial",
      PB(Da, Xa) == 1)
check("polar pair (x, xi+4/x): {D,X} constant but D has a pole (tau does not divide 4)",
      sp.simplify(PB(pol, x)) == -1 and sp.simplify(pol * x - (x * xi + 4)) == 0)

# ------------------ J2 classical band-2 bifurcation ------------------
print("Section 2: classical band-2 m=3 bifurcation (Lemma J2)")
# derive the m=3 component from the actual bracket with X = a2 x^2 + a1 x, D = b2 x^2 + b1 x
r0, r1, s0, s1, t0, t1 = sp.symbols("r0 r1 s0 s1 t0 t1")
a2s = r0 + r1 * tau
b1u = s0 + s1 * tau + t0 * tau**2
Xc2 = a2s.subs(tau, x * xi) * x**2 + a1s.subs(tau, x * xi) * x
Dc2 = lam2 * a2s.subs(tau, x * xi) * x**2 + b1u.subs(tau, x * xi) * x
m3 = x_components(PB(Dc2, Xc2), 0, 6).get(3, 0)
target = -sp.expand(2 * a2s * sp.diff(b1u, tau) - sp.diff(a2s, tau) * b1u
                    - lam2 * (2 * sp.diff(a1s, tau) * a2s - a1s * sp.diff(a2s, tau)))
check("m=3 component == 2 a2 b1' - a2' b1 - lambda2(2 a1' a2 - a1 a2')  [sign-fixed]",
      sp.expand(m3 - target) == 0)
check("b1 = lambda2*a1 is always a particular solution",
      sp.expand(target.subs({s0: lam2 * p0, s1: lam2 * p1, t0: lam2 * p2})) == 0)
w2 = 3 + 2 * tau + tau**2                       # h = w2, a2 = w2^2: square sector
homog = lambda a2f, h: sp.expand(2 * a2f * sp.diff(h, tau) - sp.diff(a2f, tau) * h)
check("square sector: h with h^2 = a2 solves the homogeneous m=3 equation",
      homog(sp.expand(w2**2), w2) == 0)
hcl = c0 + c1 * tau + c2 * tau**2 + c3 * tau**3
sols_c = sp.solve([sp.Poly(homog(tau**2 + 1, hcl), tau).coeff_monomial(tau**i)
                   for i in range(5)], [c0, c1, c2, c3], dict=True)
check("non-square a2 = tau^2+1: only the zero homogeneous solution (deg <= 3)",
      sols_c in ([], [{c0: 0, c1: 0, c2: 0, c3: 0}]))

# ------------------ J2q quantum band-2 bifurcation ------------------
print("Section 3: quantum band-2 m=3 bifurcation and m=4 proportionality")
hq = E**2 + 1
a2q = sp.expand(hq * hq.subs(E, E + 1))          # a2(E) = h(E) h(E+1): shifted square
homq = lambda a2f, h: sp.expand(h.subs(E, E + 2) * a2f - a2f.subs(E, E + 1) * h)
check("shifted-square sector: a2 = h(E)h(E+1) makes h solve the quantum homogeneous eq",
      homq(a2q, hq) == 0)
d0, d1, d2, d3 = sp.symbols("d0:4")
hgen = d0 + d1 * E + d2 * E**2 + d3 * E**3
sols_q = sp.solve([sp.Poly(homq(E**2 + 1, hgen), E).coeff_monomial(E**i)
                   for i in range(6)], [d0, d1, d2, d3], dict=True)
check("non-shifted-square a2 = E^2+1: only the zero homogeneous solution (deg <= 3)",
      sols_q in ([], [{d0: 0, d1: 0, d2: 0, d3: 0}]))
# m=4: b2(E+2) a2(E) - a2(E+2) b2(E) = 0 forces b2 || a2 for generic a2
a2gen = sp.Rational(1, 3) + 2 * E - E**2 + E**3
b2g = c0 + c1 * E + c2 * E**2 + c3 * E**3 + c4 * E**4
eq4 = sp.expand(b2g.subs(E, E + 2) * a2gen - a2gen.subs(E, E + 2) * b2g)
sols4 = sp.solve([sp.Poly(eq4, E).coeff_monomial(E**i) for i in range(8)],
                 [c0, c1, c2, c3, c4], dict=True)
ok4 = False
if sols4:
    b2sol = sp.expand(b2g.subs(sols4[0]))
    free4 = [s for s in (c0, c1, c2, c3, c4) if s not in sols4[0]]
    ratio4 = sp.cancel(b2sol / a2gen)
    ok4 = len(free4) == 1 and ratio4.free_symbols <= set(free4)
check("quantum m=4 forces b2 = lambda2*a2 (2-periodic rational ratio is constant)", ok4)

fails = [n for n, ok in RES if not ok]
print()
print("ALL BEACHHEAD CHECKS PASSED" if not fails else f"FAILURES: {fails}")
import sys; sys.exit(0 if not fails else 1)
