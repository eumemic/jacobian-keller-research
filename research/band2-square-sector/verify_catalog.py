#!/usr/bin/env python3
"""
verify_catalog.py  --  M5 band-2 square-sector exact checks and bounded catalog.

Run:  uv run --with sympy python research/band2-square-sector/verify_catalog.py

This single script is the reproducible backing for:
  * upstream-verification.md   (trust-but-verify of the classical & quantum partial cascades)
  * tame-catalog.md            (bounded tame band-2 catalog + branch residence)

Frozen conventions (M4 memo):
  classical  {G,F} = G_xi F_x - G_x F_xi,  tau = x*xi,  {xi,x}=1  (Keller = {G,F}=1)
  band-2     F = sum_{k=-2}^2 x^k a_k(tau),  G = sum x^k b_k(tau)
  membership tau | a_{-1},b_{-1} ;  tau^2 | a_{-2},b_{-2}
  quantum    A_1[x^-1], (x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E), Q_m ladder, Q_m=delta_{m0}
  gauge      C_4 => b_2 = lambda a_2 ; replace G -> G - lambda F  (b_2 -> 0)
  square     a_2 = h^2 != 0 ; after gauge b_1 = kappa h
  branches   w := gauged b_{-2}, s := a_{-2};  B0: w=0 ;  A0: w!=0,s=0 ;  A: w!=0,s!=0
  A*         h constant, kappa != 0, gauged b_{-2} != 0 AND a_{-2} != 0   (the resistant branch)

Everything below is exact (sympy, rational/symbolic).  No claim of arbitrary-degree
exhaustion of A*; the tame result is proved for the band-2 single-quadratic-shear
normal form (all parameters symbolic) and confirmed by a bounded enumeration.
"""
import sympy as sp

x, xi, t = sp.symbols("x xi tau")
u, v = sp.symbols("u v")
E = sp.symbols("E")
FAILS = []

def ok(cond, label):
    if cond:
        print("PASS", label)
    else:
        print("FAIL", label)
        FAILS.append(label)

# ---------- classical tools ----------
def PB(G, F):
    return sp.diff(G, xi)*sp.diff(F, x) - sp.diff(G, x)*sp.diff(F, xi)

def PBuv(G, F):
    return sp.diff(G, v)*sp.diff(F, u) - sp.diff(G, u)*sp.diff(F, v)

def band_decomp(P, N=16):
    Q = sp.expand(sp.expand(P).subs(xi, t/x))
    Q = sp.expand(Q * x**N)
    out = {}
    for (e,), c in sp.Poly(Q, x).as_dict().items():
        k = e - N
        c = sp.expand(c)
        if c != 0:
            out[k] = c
    return out

def is_band2(P):
    return all(-2 <= k <= 2 for k in band_decomp(P))

def Cm_eval(m, av, bv):
    tot = sp.Integer(0)
    for k in range(-2, 3):
        l = m - k
        if -2 <= l <= 2:
            ak = av.get(k, sp.Integer(0)); bl = bv.get(l, sp.Integer(0))
            tot += k*ak*sp.diff(bl, t) - l*sp.diff(ak, t)*bl
    return sp.expand(tot)

def classify(F, G):
    """Gauge-reduce (kill b2) and return branch profile for a band-2 Keller pair."""
    aF = band_decomp(F); bG = band_decomp(G)
    a2 = sp.expand(aF.get(2, 0))
    p = {'a2': a2, 'suppF': sorted(aF), 'suppG': sorted(bG)}
    if a2 == 0:
        p['square_sector'] = False; p['branch'] = 'not-square (a2=0)'
        return p
    p['square_sector'] = True
    b2 = sp.expand(bG.get(2, 0)); lam = sp.cancel(b2 / a2)
    assert sp.simplify(sp.diff(lam, t)) == 0, "C_4 violated: lambda not constant"
    hh = sp.sqrt(a2)
    p['h'] = sp.simplify(hh); p['h_constant'] = (sp.simplify(sp.diff(hh, t)) == 0)
    p['lambda'] = sp.simplify(lam)
    bG_g = {k: sp.expand(bG.get(k, 0) - lam*aF.get(k, 0)) for k in set(bG) | set(aF)}
    assert sp.expand(bG_g.get(2, 0)) == 0, "gauge failed to kill b2"
    kappa = sp.simplify(sp.cancel(sp.expand(bG_g.get(1, 0)) / hh))
    s = sp.expand(aF.get(-2, 0)); w = sp.expand(bG_g.get(-2, 0))
    p['kappa'] = kappa; p['a_minus2'] = s; p['gauged_b_minus2'] = w
    if w == 0:
        p['branch'] = 'B0'
    elif s == 0:
        p['branch'] = 'A0'
    else:
        p['branch'] = 'A'; p['mutilde'] = sp.simplify(sp.cancel(s/w)); p['mu'] = sp.simplify(sp.cancel(w/s))
    p['in_Astar'] = (p['h_constant'] and kappa != 0 and w != 0 and s != 0)
    return p

# ---------- quantum tools ----------
def sh(f, r):
    return sp.expand(f.subs(E, E + r))

def qmul(X, D):
    out = {}
    for a, f in X.items():
        for b, g in D.items():
            k = a + b
            out[k] = sp.expand(out.get(k, 0) + sh(f, b)*g)
    return {k: val for k, val in out.items() if sp.expand(val) != 0}

def qcomm(D, X):
    DX = qmul(D, X); XD = qmul(X, D)
    return {k: sp.expand(DX.get(k, 0) - XD.get(k, 0))
            for k in set(DX) | set(XD) if sp.expand(DX.get(k, 0) - XD.get(k, 0)) != 0}

def Qm(m, a, b):
    tot = sp.Integer(0)
    for k in range(-2, 3):
        l = m - k
        if -2 <= l <= 2:
            tot += sh(b.get(l, sp.Integer(0)), k)*a.get(k, sp.Integer(0)) \
                 - sh(a.get(k, sp.Integer(0)), l)*b.get(l, sp.Integer(0))
    return sp.expand(tot)

kap, lam_, A0c, ec, beta, gam, mu_ = sp.symbols("kappa lambda A e beta gamma mu")

def poly(name, deg):
    cs = sp.symbols(f"{name}0:{deg+1}")
    return sum(cs[i]*t**i for i in range(deg+1))

print("="*72)
print("SECTION 1  classical nine C_m equations (direct 2-var Poisson vs displayed)")
print("="*72)
ok(sp.simplify(PB(xi, x) - 1) == 0, "base pair orientation {xi,x}=1")
A = {k: sp.Function(f"a{k}")(t) for k in range(-2, 3)}
B = {k: sp.Function(f"b{k}")(t) for k in range(-2, 3)}
Fg = sum(x**k * A[k].subs(t, x*xi) for k in range(-2, 3))
Gg = sum(x**k * B[k].subs(t, x*xi) for k in range(-2, 3))
def Cm_abstract(m):
    return sp.expand(sum(k*A[k]*sp.diff(B[l], t) - l*sp.diff(A[k], t)*B[l]
                         for k in range(-2, 3) for l in range(-2, 3) if k + l == m))
formula = sp.expand(sum(x**m * Cm_abstract(m).subs(t, x*xi) for m in range(-4, 5)))
ok(sp.expand(PB(Gg, Fg) - formula) == 0, "all nine C_m match direct 2-var Poisson expansion")
disp = {
 4: 2*A[2]*B[2].diff(t) - 2*A[2].diff(t)*B[2],
 3: 2*A[2]*B[1].diff(t) - A[2].diff(t)*B[1] + A[1]*B[2].diff(t) - 2*A[1].diff(t)*B[2],
 2: 2*A[2]*B[0].diff(t) + A[1]*B[1].diff(t) - A[1].diff(t)*B[1] - 2*A[0].diff(t)*B[2],
 1: 2*A[2]*B[-1].diff(t) + A[2].diff(t)*B[-1] + A[1]*B[0].diff(t) - A[0].diff(t)*B[1]
    - A[-1]*B[2].diff(t) - 2*A[-1].diff(t)*B[2],
 0: 2*A[2]*B[-2].diff(t) + 2*A[2].diff(t)*B[-2] + A[1]*B[-1].diff(t) + A[1].diff(t)*B[-1]
    - A[-1]*B[1].diff(t) - A[-1].diff(t)*B[1] - 2*A[-2]*B[2].diff(t) - 2*A[-2].diff(t)*B[2],
 -1: A[1]*B[-2].diff(t) + 2*A[1].diff(t)*B[-2] + A[0].diff(t)*B[-1] - A[-1]*B[0].diff(t)
     - 2*A[-2]*B[1].diff(t) - A[-2].diff(t)*B[1],
 -2: 2*A[0].diff(t)*B[-2] + A[-1].diff(t)*B[-1] - A[-1]*B[-1].diff(t) - 2*A[-2]*B[0].diff(t),
 -3: 2*A[-1].diff(t)*B[-2] - A[-1]*B[-2].diff(t) + A[-2].diff(t)*B[-1] - 2*A[-2]*B[-1].diff(t),
 -4: 2*A[-2].diff(t)*B[-2] - 2*A[-2]*B[-2].diff(t),
}
ok(all(sp.expand(Cm_abstract(m) - disp[m]) == 0 for m in range(-4, 5)),
   "nine displayed C_m equations match")

print("="*72)
print("SECTION 2  classical kappa!=0 cascade (3.1)-(3.4), residual (4.1), family (5.1), kappa=0")
print("="*72)
h = poly("h", 3); pp = poly("p", 3); vv = poly("v", 3); ww = poly("w", 3); ss = poly("s", 3)
a2 = h**2
av = {2: a2, 1: h*pp, 0: pp**2/4 + 2*h*vv/kap - A0c,
      -1: pp*vv/kap + 2*h*ww/kap - (t+ec)/(kap*h), -2: ss}
bv = {2: sp.Integer(0), 1: kap*h, 0: kap*pp/2 + beta, -1: vv, -2: ww}
ok(sp.simplify(Cm_eval(3, av, bv)) == 0, "C_3 vanishes (b1=kappa h)")
ok(sp.simplify(Cm_eval(2, av, bv)) == 0, "C_2 vanishes -> (3.1) a1=hp, b0=kappa p/2+beta")
ok(sp.simplify(Cm_eval(1, av, bv)) == 0, "C_1 vanishes -> (3.2) a0")
central = 2*h**2*ww + h*pp*vv - kap*h*av[-1]
ok(sp.simplify(central - (t+ec)) == 0 and sp.simplify(sp.diff(central, t) - 1) == 0,
   "(3.3)-(3.4) central integral = t+e  => h | t+e")
ok(sp.simplify(Cm_eval(0, av, bv) - 1) == 0, "C_0 = 1 under integrated cascade")
# residual (4.1)
a1_, a0_, am1_ = av[1], av[0], av[-1]; b0_ = bv[0]
R = {
 -1: a1_*sp.diff(ww, t) + 2*sp.diff(a1_, t)*ww + sp.diff(a0_, t)*vv - am1_*sp.diff(b0_, t)
     - 2*ss*sp.diff(kap*h, t) - sp.diff(ss, t)*kap*h,
 -2: 2*sp.diff(a0_, t)*ww + sp.diff(am1_, t)*vv - am1_*sp.diff(vv, t) - 2*ss*sp.diff(b0_, t),
 -3: 2*sp.diff(am1_, t)*ww - am1_*sp.diff(ww, t) + sp.diff(ss, t)*vv - 2*ss*sp.diff(vv, t),
 -4: 2*sp.diff(ss, t)*ww - 2*ss*sp.diff(ww, t),
}
ok(all(sp.expand(Cm_eval(m, av, bv) - R[m]) == 0 for m in (-1, -2, -3, -4)),
   "residual system (4.1) R_{-1..-4} == C_{-1..-4}")

# family (5.1)
c0, c1 = sp.symbols("c0 c1")
def family(hpoly, epar):
    ht = hpoly.subs(t, x*xi); qt = ((t+epar)/hpoly).subs(t, x*xi)
    Y = x*ht; Z = x**-1*qt; U = Y + c0 + c1*Z
    return sp.expand(U**2 - Z/kap - A0c), sp.expand(kap*U + beta), Y, Z
Hs = sp.symbols("H0 H1 H2"); hgen = Hs[0] + Hs[1]*t + Hs[2]*t**2
F5, G5, Y5, Z5 = family(hgen, ec)
ok(sp.simplify(PB(Z5, Y5) - 1) == 0, "family (5.1): {Z,Y}=1")
ok(sp.simplify(PB(G5, F5) - 1) == 0, "family (5.1): {G,F}=1 for generic h")
dd_ = sp.symbols("d", nonzero=True)
Fc, Gc, _, _ = family(sp.Integer(1)*dd_, ec)
avc = band_decomp(Fc); bvc = band_decomp(Gc); qd = (t+ec)/dd_
exp_a = {2: dd_**2, 1: 2*c0*dd_, 0: c0**2 + 2*c1*dd_*qd - A0c, -1: 2*c0*c1*qd - qd/kap, -2: c1**2*qd**2}
exp_b = {1: kap*dd_, 0: kap*c0 + beta, -1: kap*c1*qd, -2: sp.Integer(0)}
ok(all(sp.simplify(avc.get(k, 0) - val) == 0 for k, val in exp_a.items())
   and all(sp.simplify(bvc.get(k, 0) - val) == 0 for k, val in exp_b.items()),
   "family (5.1) coefficient table matches; b_{-2}=0 (branch B0)")
alp = sp.symbols("alpha", nonzero=True)
Fnc, _, _, _ = family(alp*(t+ec), ec); anc = band_decomp(Fnc)
ok(sp.simplify(anc.get(-2, 0) - c1**2/alp**2) == 0 and sp.Poly(anc.get(-2, 0), t).degree() == 0,
   "nonconstant-h family: a_{-2}=c1^2/alpha^2 constant (fails tau^2 | a_{-2})")
ok(sp.simplify(anc.get(-1, 0).subs(c1, 0) + 1/(alp*kap)) == 0,
   "nonconstant-h, c1=0: a_{-1}=-1/(alpha kappa) constant (fails tau | a_{-1})")
# kappa=0
hk = poly("hk", 2); vk = poly("vk", 2); wk = poly("wk", 2)
avk = {2: hk**2, -1: vk}; bvk = {2: sp.Integer(0), 1: sp.Integer(0), 0: beta, -1: vk, -2: wk}
avk[1] = sp.Function('a1k')(t); avk[0] = sp.Function('a0k')(t); avk[-2] = sp.Function('am2k')(t)
ok(sp.simplify(Cm_eval(1, avk, bvk) - 2*hk*sp.diff(hk*vk, t)) == 0, "kappa=0: C_1 = 2 h (h v)'")
avk0 = dict(avk); avk0[-1] = sp.Integer(0)
bvk0 = {2: sp.Integer(0), 1: sp.Integer(0), 0: beta, -1: sp.Integer(0), -2: wk}
ok(sp.simplify(Cm_eval(0, avk0, bvk0) - 2*sp.diff(hk**2*wk, t)) == 0,
   "kappa=0,v=0: C_0 = 2 (h^2 w)' = 1  => h^2 | t+2e0 excludes nonconstant h")

print("="*72)
print("SECTION 3  quantum ladder algebra, lemmas, midpoint, tail (4.1)-(4.6)")
print("="*72)
aq = {k: sp.Function(f"qa{k}")(E) for k in range(-2, 3)}
bq = {k: sp.Function(f"qb{k}")(E) for k in range(-2, 3)}
ccm = qcomm(bq, aq)
ok(all(sp.expand(ccm.get(m, 0) - Qm(m, aq, bq)) == 0 for m in range(-4, 5)),
   "quantum Q_m formula matches abstract crossed-product commutator")
hq = sp.Function("h")(E); a2q = sp.expand(hq*sh(hq, 1))
ok(sp.expand(Qm(4, {2: a2q}, {2: lam_*a2q})) == 0, "ladder-4: b2=lambda a2 => Q_4=0")
ok(sp.expand(Qm(3, {2: a2q, 1: sp.Function("qa1")(E)}, {2: sp.Integer(0), 1: kap*hq})) == 0,
   "ladder-3: b1=kappa h (b2=0) => Q_3=0")
pq = sp.Function("p")(E); a1q = sp.expand(hq*pq); Bq = sp.Function("B")(E)
Q2 = Qm(2, {2: a2q, 1: a1q}, {2: sp.Integer(0), 1: kap*hq, 0: Bq})
ok(sp.expand(Q2 - (a2q*(sh(Bq, 2) - Bq) + kap*(sh(hq, 1)*a1q - hq*sh(a1q, 1)))) == 0,
   "ladder-2 (3.1) form; r=a1/h polynomial => h|a1")
Bm = sp.Function("Bm")(E)
lhs32 = sh(Bm, 1) + Bm - (kap*pq + gam)
ok(sp.expand((sh(lhs32, 1) - lhs32) - ((sh(Bm, 2) - Bm) - kap*(sh(pq, 1) - pq))) == 0,
   "quantum midpoint (3.2): B^{[1]}+B=kappa p+gamma  =>  ladder-2 difference eq")
a0q = sp.Function("a0")(E); uq = sp.Function("u")(E); sq = sp.Function("s")(E)
vq = sp.Function("v")(E); wq = sp.Function("w")(E)
aq2 = {2: a2q, 1: a1q, 0: a0q, -1: uq, -2: sq}
bq2 = {2: sp.Integer(0), 1: kap*hq, 0: Bq, -1: vq, -2: wq}
tail = {
 1: a2q*sh(vq, 2) - sh(a2q, -1)*vq + sh(Bq, 1)*a1q - a1q*Bq + kap*hq*a0q - kap*sh(a0q, 1)*hq,
 -4: sh(wq, -2)*sq - sh(sq, -2)*wq,
 -3: sh(wq, -1)*uq + sh(vq, -2)*sq - sh(uq, -2)*wq - sh(sq, -1)*vq,
 -2: wq*a0q + sh(vq, -1)*uq + sh(Bq, -2)*sq - sh(a0q, -2)*wq - sh(uq, -1)*vq - sq*Bq,
 -1: sh(wq, 1)*a1q + vq*a0q + sh(Bq, -1)*uq + kap*sh(hq, -2)*sq
     - sh(a1q, -2)*wq - sh(a0q, -1)*vq - uq*Bq - kap*sh(sq, 1)*hq,
 0: sh(wq, 2)*a2q + sh(vq, 1)*a1q + kap*sh(hq, -1)*uq
    - sh(a2q, -2)*wq - sh(a1q, -1)*vq - kap*sh(uq, 1)*hq,
}
ok(all(sp.expand(Qm(m, aq2, bq2) - tail[m]) == 0 for m in (1, -4, -3, -2, -1, 0)),
   "reduced quantum tail (4.1)-(4.6) all match Q_m")
bq3 = dict(bq2); bq3[-2] = mu_*sq
ok(sp.expand(Qm(-4, aq2, bq3)) == 0, "negative proportionality (4.7): w=mu s => Q_{-4}=0")
# explicit quantum boundary family with [D,X]=1
c0q, d0q, Aq = sp.symbols("c0q d0q Aq"); hlin = c0q*E + d0q
Xex = {2: sp.expand(hlin*sh(hlin, 1)), 0: sp.Integer(1), -1: Aq}
Dex = {1: kap*hlin}
Aval = sp.solve(sp.Eq(Qm(0, Xex, Dex), 1), Aq)[0]; Xex[-1] = Aval
full = qcomm(Dex, Xex)
ok(all(sp.expand(full.get(m, 0) - (1 if m == 0 else 0)) == 0 for m in range(-4, 5))
   and sp.expand(Aval) != 0 and sp.degree(sp.Poly(Aval, E)) == 0,
   "explicit quantum boundary family: [D,X]=1, u=const (fails E|u, polar)")

print("="*72)
print("SECTION 4  generic single-quadratic-shear normal form avoids A*")
print("="*72)
al, be, ga, de, aa, bb, cc, dd, mu, nu = sp.symbols("al be ga de a b c d mu nu")
X1 = al*x + be*xi; Xi1 = ga*x + de*xi
def nf_P():
    Xs, Xis = X1, sp.expand(Xi1 + mu*X1**2)
    return sp.expand(aa*Xs + bb*Xis), sp.expand(cc*Xs + dd*Xis)
def nf_Q():
    Xs, Xis = sp.expand(X1 + nu*Xi1**2), Xi1
    return sp.expand(aa*Xs + bb*Xis), sp.expand(cc*Xs + dd*Xis)
for nm, (F, G) in [("P-shear", nf_P()), ("Q-shear", nf_Q())]:
    aF = band_decomp(F); bG = band_decomp(G)
    a2 = sp.expand(aF.get(2, 0)); b2 = sp.expand(bG.get(2, 0))
    s = sp.expand(aF.get(-2, 0)); wraw = sp.expand(bG.get(-2, 0))
    lam = sp.cancel(b2/a2); wg = sp.cancel(wraw - lam*s)
    ok(sp.simplify(wg) == 0 and sp.simplify(sp.diff(a2, t)) == 0,
       f"normal form {nm}: a2 const (h constant) AND gauged b_(-2) == 0 => B0 (never A*)")
# kappa forced nonzero in square sector (uses ad-bc=1)
aF = band_decomp(nf_P()[0]); bG = band_decomp(nf_P()[1])
lamP = sp.cancel(bG.get(2, 0)/aF.get(2, 0))
b1gP = sp.cancel(bG.get(1, 0) - lamP*aF.get(1, 0))
ok(sp.simplify(b1gP.subs(dd, (1 + bb*cc)/aa) + al/bb) == 0,
   "P-shear: gauged b1 = -alpha/b ; a2=b*mu*alpha^2!=0 => kappa != 0 (kappa=0 unreachable)")

print("="*72)
print("SECTION 5  tame catalog entries (square data + branch residence)")
print("="*72)
catalog = []
def entry(label, F, G, expect):
    F = sp.expand(F); G = sp.expand(G)
    assert sp.simplify(PB(G, F) - 1) == 0, f"{label}: not a Keller pair"
    p = classify(F, G)
    for key, val in expect.items():
        got = p.get(key)
        okk = (str(got) == str(val)) if key == 'branch' or key == 'square_sector' \
              else (sp.simplify((got if got is not None else 0) - val) == 0)
        assert okk, f"{label}: {key} expected {val} got {got}"
    ok(not p.get('in_Astar', False), f"catalog [{label}] branch={p.get('branch')} "
       f"(h={p.get('h')}, kappa={p.get('kappa')}, lambda={p.get('lambda')}) NOT in A*")
    catalog.append((label, p))

# E1 affine baseline (band 1, a2=0)
entry("E1 affine baseline", 2*x - 3*xi + 5, x - xi + 7, {'square_sector': False})
# E2 swap o P(x^2): one-sided square sector (a_{-2}=0)
entry("E2 swap o P(x^2)", xi + x**2, -x,
      {'square_sector': True, 'h': 1, 'lambda': 0, 'kappa': -1, 'a_minus2': 0,
       'gauged_b_minus2': 0, 'branch': 'B0'})
# E3 memo family (5.1), h=1,e=0,c1=1,kappa=2 (two-sided a_{-2}!=0, b_{-2}=0 by construction)
kv = sp.Integer(2); Zc = xi; Uc = x + Zc
entry("E3 memo family (5.1) h=1,e=0,c1=1,kappa=2", Uc**2 - Zc/kv, kv*Uc,
      {'square_sector': True, 'h': 1, 'lambda': 0, 'kappa': 2, 'a_minus2': t**2,
       'gauged_b_minus2': 0, 'branch': 'B0'})
# E4 KEY entry: both raw negative extremes nonzero, gauge kills b_{-2}
X1c = x + xi; Sc = (X1c, xi + X1c**2)
entry("E4 A2[[1,1],[1,2]] o P(X^2) o A1[[1,1],[0,1]]", 1*Sc[0] + 1*Sc[1], 1*Sc[0] + 2*Sc[1],
      {'square_sector': True, 'h': 1, 'lambda': 2, 'kappa': -1, 'a_minus2': t**2,
       'gauged_b_minus2': 0, 'branch': 'B0'})

print("="*72)
print("SECTION 6  bounded enumeration: no sampled band-2 pair in A*")
print("="*72)
SL2 = [(u, v), (v, -u), (2*u, v/2), (u/2, 2*v), (u + v, v), (u, u + v),
       (u - v, v), (u, -u + v), (2*u + v, u + v), (u + v, u + 2*v)]
QUAD = []
for c2 in (1, -1, 2, sp.Rational(1, 2)):
    for c1 in (0, 1):
        QUAD.append((u, v + c2*u**2 + c1*u))
        QUAD.append((u + c2*v**2 + c1*v, v))
for g in SL2 + QUAD:
    assert sp.simplify(PBuv(g[1], g[0]) - 1) == 0, f"non-symplectic generator {g}"
GENS = SL2 + QUAD
def apply_gen(pair, g):
    F, G = pair
    return sp.expand(g[0].subs({u: F, v: G})), sp.expand(g[1].subs({u: F, v: G}))

seen = set(); all_pairs = [(x, xi)]; frontier = [(x, xi)]
for _ in range(3):                    # words up to length 3
    nxt = []
    for pr in frontier:
        for g in GENS:
            npr = apply_gen(pr, g)
            key = (sp.srepr(npr[0]), sp.srepr(npr[1]))
            if key in seen:
                continue
            seen.add(key); nxt.append(npr); all_pairs.append(npr)
    frontier = nxt
tally = {}; astar = 0; both_nonzero = 0; gauge_work = 0; nsq = 0
for (F, G) in all_pairs:
    if not (is_band2(F) and is_band2(G)):
        continue
    assert sp.simplify(PB(G, F) - 1) == 0, "enumerated pair not Keller"
    p = classify(F, G)
    if not p['square_sector']:
        continue
    nsq += 1; tally[p['branch']] = tally.get(p['branch'], 0) + 1
    if p['in_Astar']:
        astar += 1
    aF = band_decomp(F); bG = band_decomp(G)
    if sp.expand(aF.get(-2, 0)) != 0 and sp.expand(bG.get(-2, 0)) != 0:
        both_nonzero += 1
        if p['gauged_b_minus2'] == 0:
            gauge_work += 1
print(f"  enumerated {len(all_pairs)} distinct tame pairs; {nsq} band-2 square-sector")
print(f"  branch tally: {tally}")
print(f"  both raw negative extremes nonzero: {both_nonzero}; gauge collapses b_(-2)->0 in {gauge_work}")
ok(astar == 0 and tally.get('B0', 0) == nsq and nsq > 0,
   "bounded enumeration: every sampled band-2 square-sector pair is B0; ZERO in A*")
ok(both_nonzero > 0 and gauge_work == both_nonzero,
   "gauge does real work: all pairs with both raw extremes nonzero collapse to w=0 (B0)")
# two-quadratic-shear sweep
two = 0; two_astar = 0
for g1 in QUAD:
    for aff in SL2:
        for g2 in QUAD:
            F, G = apply_gen(apply_gen(apply_gen((x, xi), g1), aff), g2)
            if is_band2(F) and is_band2(G):
                p = classify(F, G)
                if p['square_sector']:
                    two += 1
                    if p['in_Astar']:
                        two_astar += 1
ok(two_astar == 0 and two > 0, f"two-quadratic-shear sweep: {two} band-2 square pairs, none in A*")
# Tested affine-cubic-affine words do not reach the square sector.
cub = (u, v + u**3); cub_sq = 0
for a1 in SL2:
    for a2m in SL2:
        F, G = apply_gen(apply_gen(apply_gen((x, xi), a1), cub), a2m)
        if is_band2(F) and is_band2(G) and sp.expand(band_decomp(F).get(2, 0)) != 0:
            cub_sq += 1
ok(cub_sq == 0, "tested affine-cubic-affine sweep yields no band-2 square-sector pair")

print("="*72)
if FAILS:
    print("DEFECTS FOUND:")
    for fl in FAILS:
        print("  -", fl)
    raise SystemExit("verify_catalog.py FAILED: " + "; ".join(FAILS))
print("Astar hits = 0 in the generic single-shear class and stated bounded sweeps")
print("ALL CATALOG CHECKS PASSED")
