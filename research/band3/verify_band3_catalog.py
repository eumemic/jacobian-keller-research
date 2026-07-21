#!/usr/bin/env python3
"""
verify_band3_catalog.py -- exact backing for research/band3/band3-tame-catalog.md

Run:  uv run --with sympy python research/band3/verify_band3_catalog.py

Ground truth for the BAND-3 tame catalog, classical AND quantum. Everything is
exact (sympy). The structure-lemma identities (Sections 2,3,7) are PROOFS for
the displayed normal-form class at all parameters; the bounded enumeration
(Section 6) is corroboration, never an unbounded classification.

Frozen conventions (inherited from band-2, commits f637b1a / 84978b9):
  classical  {G,F} = G_xi F_x - G_x F_xi,  tau = x*xi,  {xi,x} = 1  (Keller = {G,F}=1)
  band-3     F = sum_{k=-3}^3 x^k a_k(tau),  G = sum x^k b_k(tau)
  membership tau | a_{-1},b_{-1} ; tau^2 | a_{-2},b_{-2} ; tau^3 | a_{-3},b_{-3}
  C_m        = sum_{k+l=m}(k a_k b_l' - l a_k' b_l) = delta_{m0},  m in [-6,6]
  gauge      C_6 => b_3 = lambda a_3 ; replace G -> G - lambda F  (b_3 -> 0)
  cube       a_3 = h^3 (scalar-cube) => C_5 gives b_2 = kappa h^2
  quantum    A_1[x^-1], (x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E), E=x*d
  membership x^{-j} c(E) in A_1 iff E(E-1)...(E-j+1) | c(E)
  Q_m        = sum_{i+l=m}[ b_l(E+i) a_i(E) - a_i(E+l) b_l(E) ] = delta_{m0}
  shifted    a_3 = c h(E)h(E+1)h(E+2) => Q_5 gives b_2 = kappa h(E)h(E+1)
"""
import sympy as sp

x, xi, t = sp.symbols("x xi tau")
u, v = sp.symbols("u v")
E = sp.symbols("E")
FAILS = []

def ok(cond, label):
    print("PASS" if cond else "FAIL", label)
    if not cond:
        FAILS.append(label)

# ---------------- classical tools ----------------
def PB(G, F):
    return sp.diff(G, xi)*sp.diff(F, x) - sp.diff(G, x)*sp.diff(F, xi)

def PBuv(G, F):
    return sp.diff(G, v)*sp.diff(F, u) - sp.diff(G, u)*sp.diff(F, v)

def band_decomp(P, N=48):
    """x-level decomposition of a Laurent-in-x, poly-in-tau element."""
    Q = sp.expand(sp.expand(P).subs(xi, t/x))
    Q = sp.expand(Q * x**N)
    out = {}
    for (e,), c in sp.Poly(Q, x).as_dict().items():
        k = e - N
        c = sp.expand(c)
        if c != 0:
            out[k] = c
    return out

def band_of(P):
    d = band_decomp(P)
    return max(abs(k) for k in d) if d else 0

def is_band3(P):
    return all(-3 <= k <= 3 for k in band_decomp(P))

def Cm_eval(m, av, bv):
    tot = sp.Integer(0)
    for k in range(-3, 4):
        l = m - k
        if -3 <= l <= 3:
            ak = av.get(k, sp.Integer(0)); bl = bv.get(l, sp.Integer(0))
            tot += k*ak*sp.diff(bl, t) - l*sp.diff(ak, t)*bl
    return sp.expand(tot)

def apply_gen(pair, g):
    F, G = pair
    return sp.expand(g[0].subs({u: F, v: G})), sp.expand(g[1].subs({u: F, v: G}))

def reflect(P):
    return sp.expand(P.subs({x: xi, xi: x}, simultaneous=True))

def orient(F, G):
    """Bring a genuine band-3 pair to a_3 != 0 via pair-exchange / reflection."""
    for _ in range(4):
        if sp.expand(band_decomp(F).get(3, 0)) != 0:
            return F, G
        if sp.expand(band_decomp(G).get(3, 0)) != 0:
            F, G = G, -F; continue
        if sp.expand(band_decomp(F).get(-3, 0)) != 0:
            F, G = reflect(F), -reflect(G); continue
        if sp.expand(band_decomp(G).get(-3, 0)) != 0:
            F, G = reflect(G), -reflect(F); continue
        return None
    return None

def classify3(F, G):
    """Gauge-reduce a band-3 pair (kill b_3) and return the profile D1/D2 read residence from."""
    aF = band_decomp(F); bG = band_decomp(G)
    a3 = sp.expand(aF.get(3, 0)); am3 = sp.expand(aF.get(-3, 0))
    p = {'suppF': sorted(aF), 'suppG': sorted(bG), 'a3': a3, 'a_m3': am3}
    if a3 == 0:
        p['oriented'] = False; p['branch'] = 'unoriented (a3=0)'
        return p
    p['oriented'] = True
    lam = sp.cancel(sp.expand(bG.get(3, 0)) / a3)
    assert sp.simplify(sp.diff(lam, t)) == 0, "C_6 violated: lambda not constant"
    p['lambda'] = sp.simplify(lam)
    bG_g = {k: sp.expand(bG.get(k, 0) - lam*aF.get(k, 0)) for k in set(bG) | set(aF)}
    assert sp.expand(bG_g.get(3, 0)) == 0
    p['gauged_b3'] = sp.expand(bG_g.get(3, 0))
    p['gauged_b2'] = sp.expand(bG_g.get(2, 0))
    p['gauged_b_m2'] = sp.expand(bG_g.get(-2, 0))
    p['gauged_b_m3'] = sp.expand(bG_g.get(-3, 0))
    p['b2'] = sp.expand(bG.get(2, 0))
    p['a2'] = sp.expand(aF.get(2, 0))
    if am3 == 0:
        p['branch'] = 'onesided-top'
    elif p['gauged_b_m3'] == 0:
        p['branch'] = 'B0-band3 (collapse, mu=lambda)'
    else:
        p['branch'] = 'Astar-band3 (RESISTANT, mu!=lambda)'
    p['resistant'] = (am3 != 0 and p['gauged_b_m3'] != 0)
    return p

# ---------------- quantum tools ----------------
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
    for i in range(-3, 4):
        l = m - i
        if -3 <= l <= 3:
            tot += sh(b.get(l, sp.Integer(0)), i)*a.get(i, sp.Integer(0)) \
                 - sh(a.get(i, sp.Integer(0)), l)*b.get(l, sp.Integer(0))
    return sp.expand(tot)

def qlin(a, b, X, D):
    res = {}
    for k, val in X.items(): res[k] = sp.expand(res.get(k, 0) + a*val)
    for k, val in D.items(): res[k] = sp.expand(res.get(k, 0) + b*val)
    return {k: val for k, val in res.items() if val != 0}

def qpow(P, n):
    R = {0: sp.Integer(1)}
    for _ in range(n):
        R = qmul(R, P)
    return R


print("="*72)
print("SECTION 1  classical thirteen band-3 equations C_m match direct Poisson")
print("="*72)
ok(sp.simplify(PB(xi, x) - 1) == 0, "base orientation {xi,x}=1")
A = {k: sp.Function(f"a{k}")(t) for k in range(-3, 4)}
B = {k: sp.Function(f"b{k}")(t) for k in range(-3, 4)}
Fg = sum(x**k * A[k].subs(t, x*xi) for k in range(-3, 4))
Gg = sum(x**k * B[k].subs(t, x*xi) for k in range(-3, 4))
def Cm_abstract(m):
    return sp.expand(sum(k*A[k]*sp.diff(B[l], t) - l*sp.diff(A[k], t)*B[l]
                         for k in range(-3, 4) for l in range(-3, 4) if k + l == m))
formula = sp.expand(sum(x**m * Cm_abstract(m).subs(t, x*xi) for m in range(-6, 7)))
ok(sp.expand(PB(Gg, Fg) - formula) == 0,
   "all thirteen C_m (m=-6..6) match direct 2-var Poisson expansion")

print("="*72)
print("SECTION 2  top cascade: Wronskian C_6, cube-class C_5, bottom C_{-6}")
print("="*72)
hh = sp.Function('h')(t); kap = sp.symbols('kappa'); lam = sp.symbols('lambda')
# C_6 Wronskian => b_3 = lambda a_3
a3g = sp.Function('A3')(t)
ok(sp.expand(Cm_eval(6, {3: a3g}, {3: lam*a3g})) == 0,
   "C_6 = 3(a3 b3' - a3' b3) = 0  <=  b3 = lambda a3 (top Wronskian, W1)")
# C_5 after gauge (b3=0):  3 a3 b2' - 2 a3' b2 ; cube class a3=h^3 => b2=kappa h^2 solves it
a3cube = hh**3; b2sq = kap*hh**2
C5 = sp.expand(Cm_eval(5, {3: a3cube, 2: sp.Function('a2')(t)},
                          {3: sp.Integer(0), 2: b2sq}))
# the a2/b3 cross-terms vanish since b3=0; isolate the a3,b2 Wronskian-type part
C5_core = sp.expand(3*a3cube*sp.diff(b2sq, t) - 2*sp.diff(a3cube, t)*b2sq)
ok(C5_core == 0, "C_5 core 3 a3 b2' - 2 a3' b2 = 0 for a3=h^3, b2=kappa h^2 (cube class)")
# the invariant: C_5core = 0 & b2!=0  <=>  (b2^3/a3^2)' = 0  <=> a3 = scalar-cube
Aa = sp.Function('A')(t); Bb = sp.Function('B')(t)
core = 3*Aa*sp.diff(Bb, t) - 2*sp.diff(Aa, t)*Bb
dinv = sp.together(sp.diff(Bb**3/Aa**2, t) * Aa**3 / Bb**2)
ok(sp.simplify(dinv - core) == 0,
   "invariant: d/dt(b2^3/a3^2) * a3^3/b2^2 = (3 a3 b2' - 2 a3' b2)  =>  a3 is scalar*cube")
# bottom Wronskian C_{-6} => b_{-3} = mu a_{-3}
am3g = sp.Function('Am3')(t); mu = sp.symbols('mu')
ok(sp.expand(Cm_eval(-6, {-3: am3g}, {-3: mu*am3g})) == 0,
   "C_{-6} = 3(a_{-3}' b_{-3} - a_{-3} b_{-3}') = 0  <=  b_{-3} = mu a_{-3} (bottom Wronskian)")

print("="*72)
print("SECTION 3  CLASSICAL single-shear normal form (symbolic parameters in this class)")
print("  F,G both C-linear combos of X (linear) and Psi (sheared); the whole")
print("  |k|>=2 tower of G is lambda times that of F  =>  ONE gauge kills it.")
print("="*72)
al, be, ga, de, aa, bb, cc, dd, mu3, mu2, mu1 = sp.symbols("al be ga de a b c d mu3 mu2 mu1")
X1 = al*x + be*xi; Xi1 = ga*x + de*xi
# most general SINGLE triangular shear in the X-direction, degree up to 3:
Psi = sp.expand(Xi1 + mu3*X1**3 + mu2*X1**2 + mu1*X1)
F = sp.expand(aa*X1 + bb*Psi); G = sp.expand(cc*X1 + dd*Psi)
aF = band_decomp(F); bG = band_decomp(G)
# All-parameter theorem, with no division: b*b_k = d*a_k at every tower level.
tower_levels = (2, 3, -2, -3)
tower_ok = all(sp.expand(bb*bG.get(k, 0) - dd*aF.get(k, 0)) == 0
               for k in tower_levels)
ok(tower_ok, "single shear: b*b_k = d*a_k at every |k|=2,3 tower level (all parameters)")
# Oriented ratio/gauge check, with b != 0 encoded in the symbol assumption.
b_nz = sp.symbols('b_oriented', nonzero=True)
lam_nf = dd/b_nz
gauged = {k: sp.cancel(bG.get(k, 0).subs(bb, b_nz)
                       - lam_nf*aF.get(k, 0).subs(bb, b_nz))
           for k in tower_levels}
ok(all(gauged[k] == 0 for k in tower_levels),
   "oriented b!=0: lambda=d/b and G->G-(d/b)F kills every |k|=2,3 tower level")
# b=0 edge: use symplectic A2=I.  Then G carries the tower; after
# (F,G)->(G,-F), the first output has nonzero level 3 and ratio 0 kills the other tower.
b0_sub = {aa: 1, bb: 0, cc: 0, dd: 1, al: 1, be: 1, ga: 0, de: 1,
          mu3: 1, mu2: 1, mu1: 0}
F_b0 = sp.expand(F.subs(b0_sub)); G_b0 = sp.expand(G.subs(b0_sub))
F_ex, G_ex = G_b0, -F_b0
a_ex, b_ex = band_decomp(F_ex), band_decomp(G_ex)
ok(a_ex.get(3, 0) != 0 and all(sp.expand(b_ex.get(k, 0)) == 0 for k in tower_levels),
   "b=0 pair exchange: exchanged first output carries level 3; oriented lambda=0 gauge kills the other tower")
# (i) pure cubic (mu2=mu1=0): odd-only support; a2 identically 0
Fpc = F.subs({mu2: 0, mu1: 0}); apc = band_decomp(Fpc)
ok(sp.expand(apc.get(2, 0)) == 0 and sp.expand(apc.get(0, 0)) == 0 and set(apc) <= {-3, -1, 1, 3},
   "pure cubic shear => odd-only support {-3,-1,1,3}, a_2 = a_0 = 0")
# membership: a_{-3} divisible by tau^3, a_{-2} by tau^2 (genuine polynomial pair)
ok(sp.simplify(sp.cancel(sp.expand(aF.get(-3, 0))/t**3)).is_polynomial(t) and
   sp.simplify(sp.cancel(sp.expand(aF.get(-2, 0))/t**2)).is_polynomial(t),
   "membership: tau^3 | a_{-3}, tau^2 | a_{-2} (genuine C[x,xi] pair)")

print("="*72)
print("SECTION 4  blow-up law for displayed normalized opposite-direction shears")
print("="*72)
def shP(d): return (u, v + u**d)
def shQ(d): return (u + v**d, v)
prod_ok = True
rows = []
for a in (1, 2, 3):
    for b in (1, 2, 3):
        Fb, Gb = apply_gen(apply_gen((x, xi), shP(a)), shQ(b))
        band = max(band_of(Fb), band_of(Gb))
        rows.append((a, b, band))
        if band != a*b:
            prod_ok = False
ok(prod_ok, "displayed monic opposite-direction two-shear families have band=a*b  " + str(rows))
# quadratic+affine words realize only bands that are products of 2 -> never 3
cs = [1, -1, 2]
QUAD = [(u, v + c*u**2) for c in cs] + [(u + c*v**2, v) for c in cs]
SL2 = [(u, v), (v, -u), (u+v, v), (u, u+v), (u-v, v), (u, v-u),
       (-v, u), (-u, -v), (2*u, v/2), (u/2, 2*v)]
for g in SL2:
    assert sp.simplify(PBuv(g[1], g[0]) - 1) == 0, f"non-symplectic {g}"
for g in QUAD:
    assert sp.simplify(PBuv(g[1], g[0]) - 1) == 0, f"non-symplectic {g}"
seen = set(); frontier = [(x, xi)]; bands = set()
for _ in range(4):
    nf = []
    for pr in frontier:
        for g in QUAD + SL2:
            npr = apply_gen(pr, g)
            key = (sp.srepr(npr[0]), sp.srepr(npr[1]))
            if key in seen:
                continue
            seen.add(key); nf.append(npr)
            bands.add(max(band_of(npr[0]), band_of(npr[1])))
    frontier = nf
ok(3 not in bands and bands <= {1, 2, 4, 8, 16},
   f"quadratic+affine words (len<=4) realize bands {sorted(bands)} -- NEVER band 3")

print("="*72)
print("SECTION 5  classical catalog entries (profile + cube data + gauge tower)")
print("="*72)
catalog = []
def entry(label, F, G, expect):
    F = sp.expand(F); G = sp.expand(G)
    assert sp.simplify(PB(G, F) - 1) == 0, f"{label}: not a Keller pair"
    assert is_band3(F) and is_band3(G), f"{label}: not band 3"
    Fo, Go = orient(F, G)
    p = classify3(Fo, Go)
    for key, val in expect.items():
        got = p.get(key)
        match = (str(got) == str(val)) if key == 'branch' \
                else (sp.simplify((got if got is not None else 0) - val) == 0)
        assert match, f"{label}: {key} expected {val} got {got}"
    ok(not p.get('resistant', False),
       f"catalog [{label}] branch={p['branch']} lambda={p.get('lambda')} "
       f"a3={p.get('a3')} a2={p.get('a2')} NOT resistant")
    catalog.append((label, p))

# B3-1: minimal one-sided elementary cubic shear (xi -> xi + x^3):  support G = {-1,3}
entry("B3-1 elementary cubic shear (x, xi+x^3)  [minimal, one-sided]",
      x, xi + x**3,
      {'a3': 1, 'a_m3': 0, 'a2': 0, 'lambda': 0, 'branch': 'onesided-top'})
# B3-2: A2=[[1,1],[1,2]] o [xi -> xi + X^3 + X^2] o A1=[[1,1],[0,1]]  (two-sided, EVEN levels)
X1c = x + xi; Psic = sp.expand(xi + X1c**3 + X1c**2)
entry("B3-2 two-sided EVEN-level  A2.[cubic+quad shear].A1  [E4-band3 analogue]",
      1*X1c + 1*Psic, 1*X1c + 2*Psic,
      {'a3': 1, 'a2': 1, 'a_m3': t**3, 'lambda': 2, 'b2': 2,
       'gauged_b2': 0, 'gauged_b_m2': 0, 'gauged_b_m3': 0,
       'branch': 'B0-band3 (collapse, mu=lambda)'})
# B3-3: two-sided ODD-only  A2=[[2,1],[1,1]] o [pure cubic shear] o A1=[[1,1],[0,1]]
X1d = x + xi; Psid = sp.expand(xi + X1d**3)
entry("B3-3 two-sided ODD-only  A2.[pure cubic shear].A1  [odd support -3,-1,1,3]",
      2*X1d + 1*Psid, 1*X1d + 1*Psid,
      {'a3': 1, 'a2': 0, 'a_m3': t**3, 'lambda': 1,
       'gauged_b2': 0, 'gauged_b_m2': 0, 'gauged_b_m3': 0,
       'branch': 'B0-band3 (collapse, mu=lambda)'})

print("="*72)
print("SECTION 6  bounded enumeration: no enumerated word reaches the resistant branch")
print("="*72)
CUBE = [(u, v + c*u**3) for c in cs] + [(u + c*v**3, v) for c in cs]
for g in CUBE:
    assert sp.simplify(PBuv(g[1], g[0]) - 1) == 0, f"non-symplectic {g}"
GENS = SL2 + CUBE + QUAD
seen = set(); allp = [(x, xi)]; frontier = [(x, xi)]
for _ in range(3):          # words up to length 3
    nf = []
    for pr in frontier:
        for g in GENS:
            npr = apply_gen(pr, g)
            key = (sp.srepr(npr[0]), sp.srepr(npr[1]))
            if key in seen:
                continue
            seen.add(key); nf.append(npr); allp.append(npr)
    frontier = nf
genuine = 0; onesided = 0; twosided = 0; resistant = 0; orient_fail = 0
for (F, G) in allp:
    if not (is_band3(F) and is_band3(G)):
        continue
    if max(band_of(F), band_of(G)) != 3:
        continue
    genuine += 1
    assert sp.simplify(PB(G, F) - 1) == 0, "enumerated pair not Keller"
    o = orient(F, G)
    if o is None:
        orient_fail += 1; continue
    p = classify3(*o)
    if p['a_m3'] == 0:
        onesided += 1
    elif p['gauged_b_m3'] == 0:
        twosided += 1
    else:
        resistant += 1
print(f"  {len(allp)} distinct tame pairs; {genuine} genuine band-3; orient_fail={orient_fail}")
print(f"  one-sided top: {onesided}; two-sided collapse (mu=lambda): {twosided}; RESISTANT: {resistant}")
ok(genuine == 1216 and orient_fail == 0 and resistant == 0
   and onesided == 984 and twosided == 232,
   "enumeration: 1216 genuine band-3, ALL oriented, 984 one-sided + 232 collapse, ZERO resistant")
# two quadratic shears never reach genuine band 3
two_b3 = 0
for g1 in QUAD:
    for aff in SL2:
        for g2 in QUAD:
            Fb, Gb = apply_gen(apply_gen(apply_gen((x, xi), g1), aff), g2)
            if max(band_of(Fb), band_of(Gb)) == 3:
                two_b3 += 1
ok(two_b3 == 0, "two-quadratic-shear (quad-affine-quad) words: ZERO reach genuine band 3")

print("="*72)
print("SECTION 7  QUANTUM band 3: Q_m, shifted-cube cascade, single-shear-origin")
print("="*72)
# Q_m match abstract crossed-product commutator
aq = {k: sp.Function(f"qa{k}")(E) for k in range(-3, 4)}
bq = {k: sp.Function(f"qb{k}")(E) for k in range(-3, 4)}
ccm = qcomm(bq, aq)
ok(all(sp.expand(ccm.get(m, 0) - Qm(m, aq, bq)) == 0 for m in range(-6, 7)),
   "quantum Q_m (m=-6..6) match abstract crossed-product commutator")
# Q_6 Wronskian/Casoratian => b_3 = lambda a_3
hq = sp.Function('h')
a3q = hq(E)*hq(E+1)*hq(E+2)               # shifted cube
ok(sp.expand(Qm(6, {3: a3q}, {3: lam*a3q})) == 0,
   "Q_6 = b3(E+3)a3 - a3(E+3)b3 = 0  <=  b3 = lambda a3 (quantum top Wronskian)")
# Q_5 after gauge (b3=0): b2(E+3)a3 - a3(E+2)b2 = 0 ; shifted-cube => b2=kappa h(E)h(E+1)
b2q = kap*hq(E)*hq(E+1)                     # shifted square
Q5core = sp.expand(sh(b2q, 3)*a3q - sh(a3q, 2)*b2q)
ok(Q5core == 0,
   "Q_5 core b2(E+3)a3 - a3(E+2)b2 = 0 for a3=h(E)h(E+1)h(E+2), b2=kappa h(E)h(E+1)")
# concrete h's for the shifted-cube cascade
q5_all = True
for hpoly in [E, E+5, E**2 - 1, sp.Integer(3), 2*E - 7]:
    hf = sp.Lambda(E, hpoly)
    a3v = hf(E)*hf(E+1)*hf(E+2); b2v = kap*hf(E)*hf(E+1)
    if sp.expand(sh(b2v, 3)*a3v - sh(a3v, 2)*b2v) != 0:
        q5_all = False
ok(q5_all, "Q_5 shifted-cube cascade holds for several concrete h(E)")

# quantum single cubic shear normal form (ladder), all parameters symbolic
xX = {1: sp.Integer(1)}; dD = {-1: E}          # x at level 1 ; d = x^{-1} E at level -1
mu = sp.symbols('mu')
X1q = qlin(al, be, xX, dD)
Xi1q = qlin(ga, de, xX, dD)
Psiq = dict(Xi1q)
for k, val in qpow(X1q, 3).items():
    Psiq[k] = sp.expand(Psiq.get(k, 0) + mu*val)
Psiq = {k: val for k, val in Psiq.items() if val != 0}
Xq = qlin(aa, bb, X1q, Psiq)
Dq = qlin(cc, dd, X1q, Psiq)
# All-parameter theorem, division-free at both quantum tower levels.
qtower_levels = (3, -3)
tower_q = all(sp.expand(bb*Dq.get(k, 0) - dd*Xq.get(k, 0)) == 0
              for k in qtower_levels)
ok(tower_q, "quantum single shear: b*b_k = d*a_k at both |k|=3 tower levels (all parameters)")
# Oriented ratio/gauge check, with b != 0 encoded in the symbol assumption.
qb_nz = sp.symbols('qb_oriented', nonzero=True)
lam_q = dd/qb_nz
qgauged = {k: sp.cancel(Dq.get(k, 0).subs(bb, qb_nz)
                        - lam_q*Xq.get(k, 0).subs(bb, qb_nz))
            for k in qtower_levels}
ok(all(qgauged[k] == 0 for k in qtower_levels),
   "quantum oriented b!=0: lambda=d/b and D->D-(d/b)X kills both |k|=3 tower levels")
a3Q = sp.expand(Xq.get(3, 0)); am3Q = sp.expand(Xq.get(-3, 0))
# b=0 edge with symplectic A2=I: D carries the tower.  Pair exchange
# (X,D)->(D,-X) makes its level-3 coefficient first and leaves ratio 0.
qb0_sub = {aa: 1, bb: 0, cc: 0, dd: 1, al: 1, be: 1, ga: 0, de: 1, mu: 1}
X_b0 = {k: sp.expand(val.subs(qb0_sub)) for k, val in Xq.items()}
D_b0 = {k: sp.expand(val.subs(qb0_sub)) for k, val in Dq.items()}
X_ex = D_b0
D_ex = {k: -val for k, val in X_b0.items()}
ok(X_ex.get(3, 0) != 0 and all(sp.expand(D_ex.get(k, 0)) == 0 for k in qtower_levels),
   "quantum b=0 pair exchange: exchanged first output carries level 3; oriented lambda=0 gauge kills the other tower")
# quantum membership: a_{-3} divisible by falling factorial E(E-1)(E-2)
ff = E*(E-1)*(E-2)
ok(sp.simplify(sp.cancel(am3Q/ff)).is_polynomial(E),
   "quantum membership: E(E-1)(E-2) | a_{-3}  (falling-factorial, genuine A_1 pair)")
# concrete quantum pair [D,X]=1 : pure cubic shear X=x, D=d+mu x^3
Xex = {1: sp.Integer(1)}
Dex = {-1: E, 3: mu}
comm = qcomm(Dex, Xex)
ok(all(sp.expand(comm.get(m, 0) - (1 if m == 0 else 0)) == 0 for m in range(-6, 7)),
   "concrete quantum pair X=x, D=d+mu*x^3: [D,X]=1 exactly (support D = {-1,3})")
# quantum Fourier phi is an automorphism (orientation tool, genuinely automorphic, no sign flip)
#   phi: x -> -d, d -> x   on ladder: check [phi(D),phi(X)] = phi([D,X])
phi_x = {-1: -E}       # -d
phi_d = {1: sp.Integer(1)}   # x
ok(qcomm(phi_d, phi_x) == {0: sp.Integer(1)},
   "quantum Fourier phi (x->-d, d->x): [phi d, phi x] = 1 (genuine automorphism)")

print("="*72)
if FAILS:
    print("DEFECTS FOUND:")
    for fl in FAILS:
        print("  -", fl)
    raise SystemExit("verify_band3_catalog.py FAILED: " + "; ".join(FAILS))
print("Structure lemma holds symbolically throughout the stated single-shear normal forms;")
print("resistant branch A*-band3 reached by ZERO tame words in the stated enumeration.")
print("ALL BAND3 CATALOG CHECKS PASSED")
