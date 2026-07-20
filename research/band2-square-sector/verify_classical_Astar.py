#!/usr/bin/env python3
"""Exact checks for classical-Astar.md (branch A* closure + constant-h square sector).

Conventions: M4_proof_memo.md ({G,F} = G_xi F_x - G_x F_xi, tau = x*xi).
Identity checks are proofs of the displayed identities; the bounded-degree
sweeps are regression evidence only (the memo's degree arguments are the proofs).
Run:  python3 verify_classical_Astar.py   (requires sympy)
"""
import sympy as sp

t = sp.Symbol('t')
kappa, mutilde, A, e, alpha = sp.symbols('kappa mutilde A e alpha')
OK = True


def check(name, cond):
    global OK
    print(("PASS" if cond else "FAIL"), name, flush=True)
    OK = OK and cond


d = lambda f: sp.diff(f, t)

# ---------- Lemma 2.1: parametrization + two first integrals (s free) ----------
p = sp.Function('p')(t); v = sp.Function('v')(t); w = sp.Function('w')(t); s = sp.Function('s')(t)
P1 = sp.Function('P1')(t); P2 = sp.Function('P2')(t); V1 = sp.Function('V1')(t)
prim = [(d(P1), p), (d(P2), p**2), (d(V1), v)]

aa = {2: sp.Integer(1), 1: p, 0: p**2/4 + 2*v/kappa - A,
      -1: (p*v + 2*w - t - e)/kappa, -2: s}
bb = {2: sp.Integer(0), 1: kappa, 0: kappa*p/2, -1: v, -2: w}


def C(m):
    tot = 0
    for k in range(-2, 3):
        l = m - k
        if -2 <= l <= 2:
            tot += k*aa[k]*d(bb[l]) - l*d(aa[k])*bb[l]
    return sp.expand(tot)


check("C4..C1 == 0 and C0 == 1 under the cascade parametrization",
      all(sp.simplify(C(m)) == 0 for m in (4, 3, 2, 1)) and sp.simplify(C(0)) == 1)

Phi = kappa*p*w + v**2 + sp.Rational(1, 2)*kappa*p*(t+e) - kappa**2*s - sp.Rational(1, 2)*kappa*P1
check("Phi' == kappa*C(-1) identically (s free)",
      sp.simplify(d(Phi).subs(prim) - kappa*C(-1)) == 0)

I2 = (sp.Rational(1, 2)*kappa*p*P1 - sp.Rational(1, 4)*kappa*(t+e)*p**2
      + 2*v*w + (t+e)*v - sp.Rational(1, 4)*kappa*P2 - 2*V1)
check("I2' == kappa*C(-2) - p'*Phi identically",
      sp.simplify(d(I2).subs(prim) - kappa*C(-2) + d(p)*Phi) == 0
      or sp.simplify(d(I2).subs(prim) - kappa*C(-2) - d(p)*Phi + 2*d(p)*Phi) == 0)

# machine says residual of (I2' - kappa*C(-2) - p'Phi) is -2p'Phi, i.e. I2' = kappa*C(-2) - p'Phi:
check("I2' - kappa*C(-2) + p'*Phi == 0 (canonical form)",
      sp.simplify(d(I2).subs(prim) - kappa*C(-2) + d(p)*Phi) == 0)

# ---------- leading coefficient of the I2 p-block ----------
allP = all(sp.simplify(sp.Rational(1, 2)/(Pn+1) - sp.Rational(1, 4) - sp.Rational(1, 4)/(2*Pn+1)
           + sp.Rational(Pn**2, 2*(Pn+1)*(2*Pn+1))) == 0 for Pn in range(1, 12))
check("I2 p-block leading coefficient == -P^2/(2(P+1)(2P+1)) for P=1..11", allP)

# ---------- nonconstant h one-liner ----------
check("h = alpha*(t+e) makes (t+e)/h the constant 1/alpha",
      sp.simplify((t+e)/(alpha*(t+e)) - 1/alpha) == 0)

# ---------- S1: (2w+t)w' = 2w has no poly solution with t^2 | w, deg <= 8 ----------
s1_empty = True
for W in range(2, 9):
    cs = sp.symbols(f'c2:{W+1}')
    wp = sum(cs[i]*t**(i+2) for i in range(W-1))
    eqs = sp.Poly(sp.expand((2*wp + t)*d(wp) - 2*wp), t).all_coeffs()
    for sol in sp.solve(eqs, cs, dict=True):
        if any(sp.simplify(val) != 0 for val in sol.values()):
            s1_empty = False
check("S1 ODE (2w+t)w' = 2w: only w = 0 with t^2|w for deg <= 8", s1_empty)

# ---------- regression sweeps (small profiles; evidence only) ----------
def sweep(Pd, Vd, Wd):
    pc = sp.symbols(f'p0:{Pd+1}')
    vc = sp.symbols(f'v1:{Vd+1}')
    wc = sp.symbols(f'w2:{Wd+1}')
    pe = sum(pc[i]*t**i for i in range(Pd+1))
    ve = sum(vc[i]*t**(i+1) for i in range(Vd)) if Vd >= 1 else sp.Integer(0)
    we = sum(wc[i]*t**(i+2) for i in range(Wd-1)) if Wd >= 2 else sp.Integer(0)
    P1e = sp.integrate(pe, (t, 0, t)); P2e = sp.integrate(pe**2, (t, 0, t)); V1e = sp.integrate(ve, (t, 0, t))
    se = (kappa*pe*we + ve**2 + sp.Rational(1, 2)*kappa*(t*pe - P1e))/kappa**2
    am1e = (pe*ve + 2*we - t)/kappa
    Cm3 = sp.expand(2*d(am1e)*we - am1e*d(we) + d(se)*ve - 2*se*d(ve))
    Cm4 = sp.expand(2*d(se)*we - 2*se*d(we))
    I2e = sp.expand(sp.Rational(1, 2)*kappa*pe*P1e - sp.Rational(1, 4)*kappa*t*pe**2
                    + 2*ve*we + t*ve - sp.Rational(1, 4)*kappa*P2e - 2*V1e)
    eqs = []
    for expr in (I2e, Cm3, Cm4):
        eqs += sp.Poly(sp.expand(kappa**3*expr), t).all_coeffs()
    sols = sp.solve([sp.Eq(x, 0) for x in eqs], list(pc)+list(vc)+list(wc), dict=True)
    bad_astar, bad_a0, bad_b0 = [], [], []
    for sol in sols:
        pp = sp.simplify(pe.subs(sol)); vv = sp.simplify(ve.subs(sol))
        ww = sp.simplify(we.subs(sol)); ss = sp.simplify(se.subs(sol))
        if ww != 0 and ss != 0:
            bad_astar.append(sol)
        elif ww != 0 and ss == 0:
            bad_a0.append(sol)
        elif ww == 0:
            p_const = pp == 0 or sp.degree(pp, t) <= 0
            v_lin = vv == 0 or sp.degree(vv, t) == 1
            if not (p_const and v_lin):
                bad_b0.append(sol)
    return bad_astar, bad_a0, bad_b0


for prof in [(0, 2, 3), (1, 2, 3), (1, 3, 4)]:
    ba, b0_, bb0 = sweep(*prof)
    check(f"sweep P<={prof[0]} V<={prof[1]} W<={prof[2]}: A* empty", len(ba) == 0)
    check(f"sweep P<={prof[0]} V<={prof[1]} W<={prof[2]}: A0 empty", len(b0_) == 0)
    check(f"sweep P<={prof[0]} V<={prof[1]} W<={prof[2]}: B0 tame only", len(bb0) == 0)

print()
print("ALL CLASSICAL ASTAR CHECKS PASSED" if OK else "SOME CLASSICAL ASTAR CHECKS FAILED", flush=True)
raise SystemExit(0 if OK else 1)
