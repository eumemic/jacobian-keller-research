#!/usr/bin/env python3
"""Exact verifier for astar-band3.md: closure of the mu != lambda resistant
branch A*-band3, both faces.

Builds on Wave A (commit 99fe6ee, all verifiers green):
  - classical-band3-cascade.md / verify_classical_band3.py  (Theorem A, wall, moment)
  - quantum-band3-cascade.md   / verify_quantum_band3.py     (Q_5 wall, inhom. tail)
  - band3-tame-catalog.md       (A*-band3 = mu != lambda; counterexample-or-nothing)
  - band-k-weapons.md           (W1..W6; W3-obstruction at band 3)

Conventions (frozen, identical to the Wave A memos):
  CLASSICAL: C[x,xi], tau = x xi, {G,F} = G_xi F_x - G_x F_xi, {xi,x}=1,
    F = sum_{k=-3}^{3} x^k a_k(tau),  G = sum x^k b_k(tau),
    C_m = sum_{k+l=m}(k a_k b_l' - l a_k' b_l) = delta_{m0},  m in [-6,6];
    membership tau^j | a_{-j}, b_{-j}.
  QUANTUM: A_1[x^{-1}] = (+)_k x^k C[E], (x^a f)(x^b g)=x^{a+b} f(E+b) g(E),
    f^[r](E)=f(E+r), Q_m = sum_{k+l=m}(b_l^[k] a_k - a_k^[l] b_l) = delta_{m0};
    membership E(E-1)...(E-r+1) | a_{-r}, b_{-r}.

Every displayed algebraic identity is machine-verified.  The reflection reduction
+ Theorem A closure of branch A*-I is a PROOF (the degree/membership content of
Theorem A itself is verified in verify_classical_band3.py section 6 and re-checked
here); the bounded-degree sweeps are regression corroboration only.

A successful run ends 'ALL ASTAR BAND3 CHECKS PASSED'.
Run:  uv run --with sympy python research/band3/verify_astar_band3.py
"""
import sympy as sp

# =====================================================================
#  Machinery
# =====================================================================
x, xi, t = sp.symbols("x xi tau")
E = sp.symbols("E")
KS = range(-3, 4)
PASS = 0


def check(cond, label):
    global PASS
    if cond is not True and cond != True:  # noqa: E712
        raise AssertionError("FAIL: " + label)
    PASS += 1
    print("PASS", label)


def zero(expr, label):
    check(sp.simplify(sp.expand(expr)) == 0, label)


# --- classical ---
def Dt(f):
    return sp.diff(f, t)


def Cm(m, A, B):
    return sp.expand(sum(k * A[k] * Dt(B[l]) - l * Dt(A[k]) * B[l]
                         for k in KS for l in KS if k + l == m))


def poisson(g, f):
    return sp.diff(g, xi) * sp.diff(f, x) - sp.diff(g, x) * sp.diff(f, xi)


# --- quantum ---
def sh(f, r):
    return sp.expand(sp.sympify(f).subs(E, E + r))


def Qm(a, b, m, K=3):
    return sp.expand(sum(sh(b[l], k) * a[k] - sh(a[k], l) * b[l]
                         for k in range(-K, K + 1) for l in range(-K, K + 1)
                         if k + l == m))


def fall(r):
    return sp.prod([E - i for i in range(r)]) if r > 0 else sp.Integer(1)


print("=" * 70)
print("SECTION 0.  Sanity: C_m and Q_m machinery agree with direct brackets")
print("=" * 70)
A = {k: sp.Function(f"a{k}")(t) for k in KS}
B = {k: sp.Function(f"b{k}")(t) for k in KS}
Ff = sum(x**k * A[k].subs(t, x * xi) for k in KS)
Gf = sum(x**k * B[k].subs(t, x * xi) for k in KS)
zero(sp.expand(poisson(Gf, Ff)) - sp.expand(sum(x**m * Cm(m, A, B).subs(t, x * xi)
                                                for m in range(-6, 7))),
     "classical: all 13 C_m equal the direct 2-variable Poisson bracket")
Aq = {k: sum(sp.Symbol(f'A{k+3}_{i}') * E**i for i in range(3)) for k in KS}
Bq = {k: sum(sp.Symbol(f'B{k+3}_{i}') * E**i for i in range(3)) for k in KS}


def direct_comm(m, a, b, K=3):
    return sp.expand(sum(sh(b[l], k) * a[k] - sh(a[k], l) * b[l]
                         for k in range(-K, K + 1) for l in range(-K, K + 1)
                         if k + l == m))


for m in [6, 5, 0, -5, -6]:
    zero(Qm(Aq, Bq, m) - direct_comm(m, Aq, Bq),
         f"quantum: Q_{m} equals the direct crossed-product commutator")

print()
print("=" * 70)
print("SECTION 1 (CLASSICAL).  The reflection lemma R:(x,xi)->(xi,x)")
print("=" * 70)
# (RF)_j = tau^{-j} a_{-j};  {RG,RF} = -{G,F} o R  (Jacobian of R is -1).
mu = sp.Symbol("mutilde")
RF = Ff.subs({x: xi, xi: x}, simultaneous=True)
RG = Gf.subs({x: xi, xi: x}, simultaneous=True)
RF_rebuilt = sum(x**j * (t**(-j) * A[-j]).subs(t, x * xi) for j in KS)
zero(RF - RF_rebuilt, "reflection: (RF)_j = tau^{-j} a_{-j}  (generic coefficients)")
zero(sp.simplify(poisson(RG, RF))
     - sp.simplify(-poisson(Gf, Ff).subs({x: xi, xi: x}, simultaneous=True)),
     "reflection: {RG,RF} = -{G,F} o R  (bracket sign flips)")
# concrete genuine pair B3-2 (from band3-tame-catalog.md) reflects to a genuine pair.
pa = {3: 1, 2: 1, 1: 3 * t + 1, 0: 2 * t, -1: 3 * t**2 + 2 * t, -2: t**2, -3: t**3}
pb = {3: 2, 2: 2, 1: 6 * t + 1, 0: 4 * t, -1: 6 * t**2 + 3 * t, -2: 2 * t**2, -3: 2 * t**3}


def build(coef):
    return sum(x**k * (coef[k].subs(t, x * xi) if hasattr(coef[k], "subs") else coef[k])
               for k in KS)


FB, GB = build(pa), build(pb)
zero(poisson(GB, FB) - 1, "B3-2 is a genuine band-3 Keller pair {G,F}=1")
zero(poisson(GB.subs({x: xi, xi: x}, simultaneous=True),
             FB.subs({x: xi, xi: x}, simultaneous=True)) + 1,
     "B3-2 reflected: {RG,RF} = -1  (genuine reflected pair)")

print()
print("=" * 70)
print("SECTION 2 (CLASSICAL).  A*-band3 setup and the bottom wall split")
print("=" * 70)
# Gauge b3 = 0 (single gauge spent on top).  A* = surviving b_-3 = mutilde a_-3,
# mutilde != 0, a_-3 != 0.  C_-6 is solved by the proportionality.
Bm = dict(B)
Bm[-3] = mu * A[-3]
zero(Cm(-6, A, Bm), "A* setup: b_-3 = mutilde a_-3 solves C_-6 identically")
# Bottom wall C_-5 = -(3 a_-3 u2' - 2 a_-3' u2), u2 = b_-2 - mutilde a_-2.
u2 = B[-2] - mu * A[-2]
Ltil = 3 * A[-3] * Dt(u2) - 2 * Dt(A[-3]) * u2
zero(Cm(-5, A, Bm) + Ltil,
     "bottom wall: C_-5 = -Ltil[u2],  u2 = b_-2 - mutilde a_-2,"
     "  Ltil[u]=3 a_-3 u' - 2 a_-3' u")
# cube integrating factor for Ltil (deg a_-3 >= 3 always, since tau^3 | a_-3):
uu = sp.Function("u")(t)
zero(Dt(uu**3 / A[-3]**2) - uu**2 * (3 * A[-3] * Dt(uu) - 2 * Dt(A[-3]) * uu) / A[-3]**3,
     "cube integrating factor: (u^3/a_-3^2)' = u^2 Ltil[u]/a_-3^3")
# => A*-I (a_-3 non-cube): forced u2 = 0, i.e. b_-2 = mutilde a_-2  (Prop 3.2 reflected).
#    A*-II (a_-3 = ctil hhat^3, tau|hhat): u2 = etil hhat^2.

print()
print("=" * 70)
print("SECTION 3 (CLASSICAL).  THEOREM: branch A*-I (a_-3 non-cube) is EMPTY")
print("=" * 70)
# The reflection maps the A* pair to a genuine gauged Keller pair.  Define
#   H = mutilde*RF - RG.   Then H_j = tau^{-j}(mutilde a_-j - b_-j),
#   {H,RF} = -{RG,RF} = 1,   and H_3 = 0 (uses b_-3 = mutilde a_-3).
H = mu * RF - RG
H_rebuilt = sum(x**j * ((t**(-j)) * (mu * A[-j] - B[-j])).subs(t, x * xi) for j in KS)
zero(H - H_rebuilt, "A*-I: H := mutilde RF - RG has coeffs H_j = tau^{-j}(mutilde a_-j - b_-j)")
zero(((t**(-3)) * (mu * A[-3] - B[-3])).subs(B[-3], mu * A[-3]),
     "A*-I: H_3 = 0 (reflected top gauged, from b_-3 = mutilde a_-3)")
# {H,RF} = mutilde{RF,RF} - {RG,RF} = -{RG,RF}; with {RG,RF}=-{G,F}oR and {G,F}=1:
zero(sp.simplify(poisson(mu * RF - RG, RF))
     - sp.simplify(-poisson(RG, RF)),
     "A*-I: {H,RF} = -{RG,RF}  (so = 1 on a genuine pair)")
# (RF)_3 = tau^{-3} a_-3 is a nonzero polynomial (tau^3 | a_-3), and
#   a_-3 is a scalar cube  <=>  (RF)_3 is a scalar cube  (tau^3 is a cube):
gg = sp.Function("g")(t)
zero((t**3 * gg**3) - (t * gg)**3,
     "cube-equivalence: a_-3 = tau^3 g^3 <=> (RF)_3 = g^3 (tau^3 is a cube)")
# Theorem A (verify_classical_band3.py sec 6, PROVED): a genuine gauged band-3
# pair with top coefficient nonzero and NON-cube is empty.  Re-verify its rung
# reductions here so the A*-I closure is self-contained.  (a3 plays the role of (RF)_3.)
a3g = sp.Function("a3")(t)
Bc = {k: sp.Function(f"b{k}")(t) for k in KS}
Bc[3] = sp.Integer(0)
Ag = {3: a3g, 2: sp.Function("a2")(t), 1: sp.Function("a1")(t), 0: sp.Function("a0")(t),
      -1: sp.Function("am1")(t), -2: sp.Function("am2")(t), -3: sp.Function("am3")(t)}
zero(Cm(5, Ag, Bc) - (3 * a3g * Dt(Bc[2]) - 2 * Dt(a3g) * Bc[2]),
     "Theorem A rung: gauged C_5 = 3 a3 b2' - 2 a3' b2 (non-cube => b2=0)")
Bc2 = dict(Bc); Bc2[2] = sp.Integer(0)
zero(Cm(4, Ag, Bc2) - (3 * a3g * Dt(Bc2[1]) - Dt(a3g) * Bc2[1]),
     "Theorem A rung: gauged C_4|b2=0 = 3 a3 b1' - a3' b1 (non-cube => b1=0)")
Bc3 = dict(Bc2); Bc3[1] = sp.Integer(0)
zero(Cm(3, Ag, Bc3) - 3 * a3g * Dt(Bc3[0]), "Theorem A rung: gauged C_3 = 3 a3 b0'")
b0c = sp.Symbol("b0c"); Bc4 = dict(Bc3); Bc4[0] = b0c
zero(Cm(2, Ag, Bc4) - (3 * a3g * Dt(Bc4[-1]) + Dt(a3g) * Bc4[-1]),
     "Theorem A rung: gauged C_2 = 3 a3 b_-1' + a3' b_-1 (deg a3>=1 => b_-1=0)")
Bc5 = dict(Bc4); Bc5[-1] = sp.Integer(0)
zero(Cm(1, Ag, Bc5) - (3 * a3g * Dt(Bc5[-2]) + 2 * Dt(a3g) * Bc5[-2]),
     "Theorem A rung: gauged C_1 = 3 a3 b_-2' + 2 a3' b_-2 (=> b_-2=0)")
Bc6 = dict(Bc5); Bc6[-2] = sp.Integer(0)
zero(Cm(0, Ag, Bc6) - 3 * Dt(a3g * Bc6[-3]),
     "Theorem A endgame: gauged C_0 = 3 (a3 b_-3)' = 1 contradicts tau^3|b_-3")
print("   => (RF, H) is a genuine gauged band-3 pair; if (RF)_3 = tau^-3 a_-3 is")
print("      NOT a cube, Theorem A makes it EMPTY, contradicting its existence.")
print("      THEOREM: A*-I is empty; hence A* forces BOTH extremes to be scalar cubes.")

print()
print("=" * 70)
print("SECTION 4 (CLASSICAL).  A*-I reflected cascade + positive parametrization")
print("=" * 70)
# In A*-I, u2 = 0 => b_-2 = mutilde a_-2.  Then C_-4 is the reflected Theorem-A M-rung:
BI = dict(Bm); BI[-2] = mu * A[-2]
u1 = mu * A[-1] - B[-1]
Mtil = 3 * A[-3] * Dt(u1) - Dt(A[-3]) * u1
zero(Cm(-4, A, BI) - Mtil,
     "A*-I: C_-4 = Mtil[u1], u1 = mutilde a_-1 - b_-1, Mtil[u]=3 a_-3 u' - a_-3' u")
# Positive cascade (constant-h slice a3=c, gauge b3=0, wall e=0): the section-7.3
# parametrization solves C_5,C_4,C_3,C_2 identically.
c, k1, beta, gam = sp.symbols("c kappa1 beta gamma")
a2, a1, a0 = sp.Function("a2")(t), sp.Function("a1")(t), sp.Function("a0")(t)
am1, am2, am3 = sp.Function("am1")(t), sp.Function("am2")(t), sp.Function("am3")(t)
bm2, bm3 = sp.Function("bm2")(t), sp.Function("bm3")(t)
b1 = k1
b0 = k1 * a2 / (3 * c) + beta
bm1 = k1 * a1 / (3 * c) - k1 * a2**2 / (9 * c**2) + gam
Ap = {3: c, 2: a2, 1: a1, 0: a0, -1: am1, -2: am2, -3: am3}
Bp = {3: sp.Integer(0), 2: sp.Integer(0), 1: b1, 0: b0, -1: bm1, -2: bm2, -3: bm3}
for m in [5, 4, 3, 2]:
    zero(Cm(m, Ap, Bp), f"constant-h e=0 positive cascade: C_{m} = 0 identically")
# Moment W4:  C_0 = M',  M = sum_{k>=1} k(a_k b_-k - a_-k b_k);  membership => M = tau.
Mmom = sum(k * (Ap[k] * Bp[-k] - Ap[-k] * Bp[k]) for k in [1, 2, 3])
zero(Cm(0, Ap, Bp) - Dt(Mmom), "moment (W4): C_0 = M'  with M the level-weighted moment")
# M(0) = 0 for ANY genuine pair: every product in M carries a membership-vanishing
# factor (tau^j | a_-j, b_-j).  Check on a fully generic membership representative.
Ag0 = {3: t**2 + 3, 2: t + 1, 1: 5, 0: 7, -1: t * (t + 2), -2: t**2 * (t - 1), -3: t**3 * (t + 4)}
Bg0 = {3: t - 2, 2: 9, 1: t**3, 0: 4, -1: t * (3 * t + 1), -2: t**2 * (t + 5), -3: t**3 * (2 * t - 1)}
Mrep = sum(k * (Ag0[k] * Bg0[-k] - Ag0[-k] * Bg0[k]) for k in [1, 2, 3])
check(sp.simplify(Mrep.subs(t, 0)) == 0,
      "moment M(0)=0 under genuine membership  (=> M(tau)=tau)")

print()
print("=" * 70)
print("SECTION 5 (CLASSICAL).  Doubly-minimal A* (a_-3 = ctil tau^3): closure data")
print("=" * 70)
# In A*-II with hhat = tau (a_-3 = ctil tau^3, the minimal cube) the bottom wall
# u2 = etil tau^2, so b_-2 = mutilde a_-2 + etil tau^2, and C_-5 = 0 is built in.
etil = sp.Symbol("etil")
Adm = {3: sp.Integer(1), 2: a2, 1: a1, 0: a0, -1: am1, -2: am2, -3: t**3}
b1d = k1
b0d = k1 * a2 / 3
bm1d = k1 * a1 / 3 - k1 * a2**2 / 9 + gam
bm2d = mu * am2 + etil * t**2
bm3d = mu * t**3
# a_-1 eliminated by the moment M = tau:
am1d = sp.expand((3 * mu * t**3 + 2 * a2 * bm2d + a1 * bm1d - t) / k1)
Adm[-1] = am1d
Bdm = {3: sp.Integer(0), 2: sp.Integer(0), 1: b1d, 0: b0d, -1: bm1d, -2: bm2d, -3: bm3d}
zero(sum(k * (Adm[k] * Bdm[-k] - Adm[-k] * Bdm[k]) for k in [1, 2, 3]) - t,
     "doubly-minimal: moment M = tau holds by construction (defines a_-1)")
zero(Cm(-5, Adm, Bdm), "doubly-minimal: C_-5 = 0 identically (bottom wall built in)")
zero(Cm(-6, Adm, Bdm), "doubly-minimal: C_-6 = 0 identically (b_-3 = mutilde a_-3)")


# --- exact bounded-degree sweep: no A* pair (regression corroboration) ---
def sweep_classical(mu_v, k1_v, degs, drop=()):
    P2, P1, P0, Q2 = degs
    syms = []

    def poly(name, dmn, dmx):
        cs = [sp.Symbol(f"{name}{j}") for j in range(dmn, dmx + 1)]
        syms.extend(cs)
        return sum(cs[j - dmn] * t**j for j in range(dmn, dmx + 1))
    A2 = poly("A2", 0, P2); A1 = poly("A1", 0, P1); A0 = poly("A0", 0, P0)
    M2 = poly("M2", 2, Q2)
    et_, gm_ = sp.symbols("et_ gm_"); syms += [et_, gm_]
    M3 = t**3
    B1_ = sp.Integer(k1_v); B0_ = B1_ * A2 / 3
    Bm1_ = B1_ * A1 / 3 - B1_ * A2**2 / 9 + gm_
    Bm3_ = mu_v * M3; Bm2_ = mu_v * M2 + et_ * t**2
    M1 = sp.expand((3 * mu_v * M3 + 2 * A2 * Bm2_ + A1 * Bm1_ - t) / B1_)
    AA = {3: sp.Integer(1), 2: A2, 1: A1, 0: A0, -1: M1, -2: M2, -3: M3}
    BB = {3: sp.Integer(0), 2: sp.Integer(0), 1: B1_, 0: B0_, -1: Bm1_, -2: Bm2_, -3: Bm3_}
    eqs = []
    for m in [mm for mm in (1, -1, -2, -3, -4) if mm not in drop]:
        e = sp.expand(sp.numer(sp.together(Cm(m, AA, BB))))
        if e != 0:
            eqs += sp.Poly(e, t).all_coeffs()
    eqs.append(sp.Poly(sp.expand(sp.numer(sp.together(M1))), t).eval(0))  # tau | a_-1
    sol = sp.solve(eqs, syms, dict=True)
    hits = 0
    for s in sol:
        Asub = {k: (AA[k].subs(s) if hasattr(AA[k], "subs") else AA[k]) for k in KS}
        Bsub = {k: (BB[k].subs(s) if hasattr(BB[k], "subs") else BB[k]) for k in KS}
        if all(sp.simplify(Cm(m, Asub, Bsub) - (1 if m == 0 else 0)) == 0 for m in range(-6, 7)):
            hits += 1
    return hits


for degs in [(0, 1, 1, 2), (1, 1, 1, 2), (1, 2, 2, 3)]:
    check(sweep_classical(sp.Integer(2), 1, degs) == 0,
          f"classical A* sweep degs(P2,P1,P0,Q2)={degs}, mu=2: no genuine A* pair")
# redundancy: dropping C_-3 (and C_1) still empty (mirrors band-2 A* using neither)
check(sweep_classical(sp.Integer(2), 1, (1, 2, 2, 3), drop=(-3,)) == 0,
      "classical A* still empty after dropping C_-3 (redundant, as in band 2)")
check(sweep_classical(sp.Integer(2), 1, (1, 2, 2, 3), drop=(1,)) == 0,
      "classical A* still empty after dropping C_1 (redundant)")
# spurious-candidate destruction: the a2=a-2=b-2=0 profile 'solution' fails C_0.
a1v = -1 / (sp.Symbol("M11") * mu)
Asp = {3: 1, 2: 0, 1: a1v, 0: sp.Symbol("A00"), -1: sp.Symbol("M11") * t, -2: 0, -3: t**3}
Bsp = {3: 0, 2: 0, 1: -1 / sp.Symbol("M11"), 0: 0,
       -1: (-1 / sp.Symbol("M11")) * a1v / 3 - 1 / (3 * sp.Symbol("M11")**2 * mu),
       -2: 0, -3: mu * t**3}
zero(Cm(0, Asp, Bsp) - (1 + 9 * mu * t**2),
     "spurious pole-free candidate has C_0 = 1 + 9 mutilde tau^2 (fails unless mu=0):"
     " NO A* counterexample")

print()
print("=" * 70)
print("SECTION 6 (QUANTUM).  A*-band3 on the DC1 face: setup, inhom. tail, sweep")
print("=" * 70)
mu3 = sp.Symbol("mu3")
Bqm = dict(Bq); Bqm[-3] = mu3 * Aq[-3]
zero(Qm(Aq, Bqm, -6), "quantum A* setup: b_-3 = mu3 a_-3 solves Q_-6")
# positive cascade (trivial-cube top a3=1, gauge b3=0, wall b2=0): Q_4 = b1^[3] - b1.
a2q, a1q, a0q = (sum(sp.Symbol(f"{n}_{j}") * E**j for j in range(3)) for n in ("A2q", "A1q", "A0q"))
b1q = sum(sp.Symbol(f"B1q_{j}") * E**j for j in range(3))
aq = {3: sp.Integer(1), 2: a2q, 1: a1q, 0: a0q, -1: sp.Integer(0), -2: sp.Integer(0), -3: sp.Integer(0)}
bq = {3: sp.Integer(0), 2: sp.Integer(0), 1: b1q, 0: sp.Integer(0), -1: sp.Integer(0),
      -2: sp.Integer(0), -3: sp.Integer(0)}
zero(Qm(aq, bq, 4) - (sh(b1q, 3) - b1q),
     "quantum positive cascade: Q_4|(a3=1,b2=0) = b1^[3] - b1  (=> b1 constant)")
# the INHOMOGENEOUS negative tail: Q_-5, Q_-4 carry a mu3-proportional source
# (the genuinely new band-3 feature; cross-coupling lambda3 != mu3).
pw = sh(Bqm[-2], -3) * Aq[-3] - sh(Aq[-3], -2) * Bqm[-2]
src = sh(Aq[-3], -2) * Aq[-2] - sh(Aq[-2], -3) * Aq[-3]
zero(Qm(Aq, Bqm, -5) - (pw + mu3 * src),
     "quantum Q_-5 = [bottom wall in b_-2] + mu3*[a_-3,a_-2 source]  (inhomogeneous)")
d1 = sh(Bqm[-1], -3) * Aq[-3] - sh(Aq[-3], -1) * Bqm[-1]
d2 = sh(Aq[-3], -1) * Aq[-1] - sh(Aq[-1], -3) * Aq[-3]
d3 = sh(Bqm[-2], -2) * Aq[-2] - sh(Aq[-2], -2) * Bqm[-2]
zero(Qm(Aq, Bqm, -4) - (d1 + mu3 * d2 + d3),
     "quantum Q_-4 = [b_-1,a_-3 stagger] + mu3*[source] + [b_-2,a_-2 square]")
# reflection subtlety: the naive E -> -E-1 reflection breaks A_1 membership when
# a3 is constant (the reflected x^-3 coefficient is not falling-factorial divisible).
check(sp.rem(sp.Poly(sp.Integer(1).subs(E, -E - 1), E), sp.Poly(fall(3), E)) != 0,
      "quantum reflection needs a falling-factorial twist: naive E->-E-1 on a3=1 "
      "breaks membership (so the classical reflection+Theorem-A route does NOT "
      "transcribe; quantum Theorem A is itself open, quantum-band3 sec 3.5)")


def sweep_quantum(mu_v, ct, k1, dP, dN):
    syms = []

    def g(name, dm, fac=1):
        cs = [sp.Symbol(f"{name}{j}") for j in range(dm + 1)]
        syms.extend(cs)
        return fac * sum(cs[j] * E**j for j in range(dm + 1))
    a2_ = g("qa2", dP); a1_ = g("qa1", dP); a0_ = g("qa0", dP)
    am1_ = g("qP", dN, fall(1)); am2_ = g("qQ", dN, fall(2))
    b0_ = g("qD", dP + 1); bm1_ = g("qM", dP + 1); bm2_ = g("qF", dN, fall(2))
    am3_ = ct * fall(3)
    aa = {3: sp.Integer(1), 2: a2_, 1: a1_, 0: a0_, -1: am1_, -2: am2_, -3: am3_}
    bb = {3: sp.Integer(0), 2: sp.Integer(0), 1: sp.Integer(k1), 0: b0_, -1: bm1_,
          -2: bm2_, -3: mu_v * am3_}
    eqs = []
    for m in [4, 3, 2, 1, -1, -2, -3, -4, -5]:
        e = sp.expand(Qm(aa, bb, m))
        if e != 0:
            eqs += sp.Poly(e, E).all_coeffs()
    e0 = sp.expand(Qm(aa, bb, 0) - 1)
    if e0 != 0:
        eqs += sp.Poly(e0, E).all_coeffs()
    sol = sp.solve(eqs, syms, dict=True)
    genuine, nonzero_b0 = 0, 0
    for s in sol:
        asub = {k: (aa[k].subs(s) if hasattr(aa[k], "subs") else aa[k]) for k in range(-3, 4)}
        bsub = {k: (bb[k].subs(s) if hasattr(bb[k], "subs") else bb[k]) for k in range(-3, 4)}
        if all(sp.simplify(Qm(asub, bsub, m) - (1 if m == 0 else 0)) == 0 for m in range(-6, 7)):
            genuine += 1
            if any(sp.simplify(asub[j]) != 0 for j in (1, -1)):
                nonzero_b0 += 1
    return genuine, nonzero_b0


for degs in [(0, 1), (1, 1), (1, 2)]:
    g, _ = sweep_quantum(sp.Integer(2), 1, 1, *degs)
    check(g == 0, f"quantum A* sweep dP,dN={degs}, mu3=2: no genuine A* pair")
# validation: with mu3 = 0 the SAME machinery finds genuine B0 pairs (sweep is sound).
_, nz = sweep_quantum(sp.Integer(0), 1, 1, 1, 1)
check(nz >= 1, "quantum sweep validation: mu3=0 (B0) yields genuine nonzero pairs")

print()
print(f"({PASS} checks)")
print("ALL ASTAR BAND3 CHECKS PASSED")
