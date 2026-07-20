#!/usr/bin/env python3
"""
Thirteen exact symbolic checks supporting the audited J3 band-2 argument,
classically and quantumly.

The checks cover selected reduced ladder components and product identities for
generic coefficient polynomials of degree at most two, two bounded finite
examples, and 2-periodicity for polynomials of degree at most four. They do not
check the full unreduced coefficient systems, the general homogeneous lemmas,
rational periodicity, arbitrary-degree branches, localization or membership,
or proof completeness. Passing is computational support only, not peer review.

Expected output: all PASS, then "ALL J3 CHECKS PASSED".
"""
import warnings; warnings.filterwarnings("ignore")
import sympy as sp

x, xi, tau, E = sp.symbols("x xi tau E")
lam2, mu2, alpha, gamma = sp.symbols("lambda2 mu2 alpha gamma")
RES = []
def check(name, ok):
    RES.append((name, bool(ok)))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

# ---------------- classical machinery ----------------
def PB(f, g):
    return sp.expand(sp.diff(f, xi) * sp.diff(g, x) - sp.diff(f, x) * sp.diff(g, xi))

def x_components(expr, kmin, kmax):
    comp, e = {}, sp.expand(expr.subs(xi, tau / x))
    for k in range(kmin, kmax + 1):
        c = sp.expand(e.coeff(x, k))
        if c != 0:
            comp[k] = c
    assert sp.expand(e - sum(x**k * v for k, v in comp.items())) == 0
    return comp

p = sp.symbols("p0:3"); q = sp.symbols("q0:3"); r = sp.symbols("r0:3"); s = sp.symbols("s0:3")
a2c  = p[0] + p[1]*tau + p[2]*tau**2
a1c  = q[0] + q[1]*tau + q[2]*tau**2
am1c = r[0] + r[1]*tau + r[2]*tau**2
am2c = s[0] + s[1]*tau + s[2]*tau**2

print("Section 1: classical reduced system and invariants")
lift = lambda f, k: (f.subs(tau, x*xi)) * x**k
Xc = lift(a2c, 2) + lift(a1c, 1) + alpha + lift(am1c, -1) + lift(am2c, -2)
Dc = (lift(lam2*a2c, 2) + lift(lam2*a1c, 1) + (lam2*alpha + gamma)
      + lift(mu2*am1c, -1) + lift(mu2*am2c, -2))
comp = x_components(PB(Dc, Xc), -6, 6)
U = sp.expand(sp.diff(a2c, tau)*am1c + 2*a2c*sp.diff(am1c, tau))
L = sp.expand(sp.diff(am2c, tau)*a1c + 2*am2c*sp.diff(a1c, tau))
moment = sp.expand(sp.diff(2*a2c*am2c + a1c*am1c, tau))
check("components m=+-2,+-3,+-4 vanish for the reduced pair",
      all(comp.get(m, 0) == 0 for m in (2, 3, 4, -2, -3, -4)))
check("m=+1 component == (mu2-lambda2)*U", sp.expand(comp.get(1, 0) - (mu2-lam2)*U) == 0)
check("m=-1 component == (mu2-lambda2)*L", sp.expand(comp.get(-1, 0) - (mu2-lam2)*L) == 0)
check("m=0 component == (mu2-lambda2)*(2 a2 a_{-2} + a1 a_{-1})'",
      sp.expand(comp.get(0, 0) - (mu2-lam2)*moment) == 0)
check("invariants: a_{-1}*U == (a2 a_{-1}^2)'  and  a1*L == (a_{-2} a1^2)'",
      sp.expand(am1c*U - sp.diff(a2c*am1c**2, tau)) == 0 and
      sp.expand(a1c*L - sp.diff(am2c*a1c**2, tau)) == 0)
# no nonzero polynomial a2 solves U=0 for the sample nonconstant a_{-1} = 1+tau
c = sp.symbols("c0:5")
a2u = sum(c[i]*tau**i for i in range(5))
Ueq = sp.expand(sp.diff(a2u, tau)*(1+tau) + 2*a2u)
solsU = sp.solve([sp.Poly(Ueq, tau).coeff_monomial(tau**i) for i in range(6)],
                 list(c), dict=True)
check("U=0 with a_{-1}=1+tau has only the zero polynomial solution a2 (deg<=4)",
      solsU in ([], [{ci: 0 for ci in c}]))

# ---------------- quantum machinery ----------------
print("Section 2: quantum reduced system and invariants")
def clean(A): return {k: sp.expand(v) for k, v in A.items() if sp.expand(v) != 0}
def mul(A, B):
    C = {}
    for a_, f in A.items():
        for b_, g in B.items():
            C[a_+b_] = C.get(a_+b_, 0) + f.subs(E, E+b_)*g
    return clean(C)
def add(A, B):
    C = dict(A)
    for k, v in B.items(): C[k] = C.get(k, 0) + v
    return clean(C)
def brk(A, B): return add(mul(A, B), {k: -v for k, v in mul(B, A).items()})

a2q  = p[0] + p[1]*E + p[2]*E**2
a1q  = q[0] + q[1]*E + q[2]*E**2
am1q = r[0] + r[1]*E + r[2]*E**2
am2q = s[0] + s[1]*E + s[2]*E**2
Xq = clean({2: a2q, 1: a1q, 0: alpha, -1: am1q, -2: am2q})
Dq = clean({2: lam2*a2q, 1: lam2*a1q, 0: lam2*alpha + gamma,
            -1: mu2*am1q, -2: mu2*am2q})
B = brk(Dq, Xq)
Uq = sp.expand(am1q.subs(E, E+2)*a2q - a2q.subs(E, E-1)*am1q)
Lq = sp.expand(am2q.subs(E, E+1)*a1q - a1q.subs(E, E-2)*am2q)
W2 = sp.expand(a2q.subs(E, E-2)*am2q); W1 = sp.expand(a1q.subs(E, E-1)*am1q)
momq = sp.expand(W2.subs(E, E+2) - W2 + W1.subs(E, E+1) - W1)
check("quantum components m=+-2,+-3,+-4 vanish for the reduced pair",
      all(B.get(m, 0) == 0 for m in (2, 3, 4, -2, -3, -4)))
check("quantum m=+1 component == (mu2-lambda2)*U_q",
      sp.expand(B.get(1, 0) - (mu2-lam2)*Uq) == 0)
check("quantum m=-1 component == (mu2-lambda2)*L_q",
      sp.expand(B.get(-1, 0) - (mu2-lam2)*Lq) == 0)
check("quantum m=0 component == (mu2-lambda2)*(W2(E+2)-W2(E)+W1(E+1)-W1(E))",
      sp.expand(B.get(0, 0) - (mu2-lam2)*momq) == 0)
G = sp.expand(a2q*am1q.subs(E, E+1)*am1q.subs(E, E+2))
H = sp.expand(am2q*a1q.subs(E, E-1)*a1q.subs(E, E-2))
check("invariants: G(E)-G(E-1) == a_{-1}(E+1)*U_q(E)  and  H(E+1)-H(E) == a1(E-1)*L_q(E)",
      sp.expand(G - G.subs(E, E-1) - am1q.subs(E, E+1)*Uq) == 0 and
      sp.expand(H.subs(E, E+1) - H - a1q.subs(E, E-1)*Lq) == 0)
a2uq = sum(c[i]*E**i for i in range(5))
Uq0 = sp.expand((E+3)*a2uq - a2uq.subs(E, E-1)*(E+1))   # a_{-1} = E+1
solsUq = sp.solve([sp.Poly(Uq0, E).coeff_monomial(E**i) for i in range(6)],
                  list(c), dict=True)
check("quantum U_q=0 with a_{-1}=E+1 has only the zero polynomial a2 (deg<=4)",
      solsUq in ([], [{ci: 0 for ci in c}]))
pg = sum(c[i]*E**i for i in range(5))
sols2p = sp.solve([sp.Poly(sp.expand(pg.subs(E, E+2) - pg), E).coeff_monomial(E**i)
                   for i in range(5)], list(c), dict=True)
ok2p = False
if sols2p:
    psol = sp.expand(pg.subs(sols2p[0]))
    free = [ci for ci in c if ci not in sols2p[0]]
    ok2p = len(free) == 1 and sp.diff(psol, E) == 0
check("2-periodic polynomials of degree <= 4 are constants", ok2p)

fails = [n for n, ok in RES if not ok]
print()
print("ALL J3 CHECKS PASSED" if not fails else f"FAILURES: {fails}")
import sys; sys.exit(0 if not fails else 1)
