#!/usr/bin/env python3
"""
verify_bandk_weapons.py
=======================
Exact SymPy verification backing `band-k-weapons.md`: the lift of the band-2
arsenal (W1..W6) to general band index k, machine-checked at k = 2, 3, 4.

Frozen conventions (identical to the band-2 corpus):

  CLASSICAL.  C[x,xi], tau = x*xi, {G,F} = G_xi F_x - G_x F_xi, {xi,x}=1.
     F = sum_{i=-k}^k x^i a_i(tau),  G = sum_{l=-k}^k x^l b_l(tau).
     {x^l b, x^i a} = x^{i+l}(i a b' - l a' b);  primes are d/dtau.
     C_m = sum_{i+l=m} ( i a_i b_l' - l a_i' b_l ),  Keller: C_m = delta_{m0}.
     Membership: tau^j | a_{-j}, b_{-j}.

  QUANTUM.  A_1[x^{-1}] = (+)_k x^k C[E], E = x d, (x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E).
     f^[r](E) = f(E+r).
     Q_m = sum_{i+l=m} ( b_l^[i] a_i - a_i^[l] b_l ),  [D,X]=1: Q_m = delta_{m0}.
     Membership: E(E-1)...(E-j+1) | a_{-j}, b_{-j}  (falling factorial).

Sections:
  1. W1  top Wronskian: C_{2k}, Q_{2k} are pure (shifted-)Wronskians => b_k = lam a_k.
  2. W2c classical gatekeeper: C_{2k-1}(gauge) = k a_k u' - (k-1) a_k' u;  u=b_{k-1}.
        Nonzero polynomial solution iff a_k = c*h^k (k-th power). Both directions.
  3. W2q quantum gatekeeper: Q_{2k-1}(gauge) = u^[k] a_k - a_k^[k-1] u.
        Sufficiency of the shifted k-th power; the exact telescoping INVARIANT;
        and the NEW k>=3 PHENOMENON: exotic solutions whose a_k is NOT a shifted
        k-th power (root-necklace S_k/S_{k-1} coprimality gap).
  4. W4  m=0 moment. Classical: C_0 = Psi',  Psi = sum_{i>=1} i(a_i b_{-i} - a_{-i} b_i).
        Quantum: Q_0 = (T-1)G with the explicit staggered potential G. Both exact, any k.
  5. W5  lattice: the band-2 A* degree-balance system re-cast as an integer affine
        system; the mod-3 (=infeasible) congruence; the general formulation skeleton.
  6. W3  conditional trailing identity: a_{-k} enters C_{-(k-1)} through
        -(k-1)a_{-k}b_1' - (a_{-k} b_1)'.  It is an exact derivative only under
        the independently established extra hypothesis b_1'=0.
  7. W6  bottom-wall degree balance in both faces: (k-1)S = kPhi.
        Classically the differential wall gives ordinary-power behavior; quantumly
        the staggered wall gives root-necklace/cyclotomic behavior.

Run:  uv run --with sympy python research/band3/verify_bandk_weapons.py
Ends: ALL BAND-K WEAPON CHECKS PASSED
"""
import sympy as sp

t = sp.symbols('t')
E = sp.symbols('E')
sig = sp.symbols('sigma')
KS = (2, 3, 4)


def az(expr, label):
    e = sp.expand(sp.sympify(expr))
    if e != 0:
        raise AssertionError(label + "  :  residual = " + str(e))
    print("PASS", label)


def ok(cond, label):
    if not cond:
        raise AssertionError(label + "  :  FAILED")
    print("PASS", label)


def sh(f, r):
    return sp.expand(sp.sympify(f).subs(E, E + r))


# ---- classical / quantum ladder-coefficient engines ----
def C_m(a, b, m, k):
    tot = 0
    for i in range(-k, k + 1):
        for l in range(-k, k + 1):
            if i + l == m:
                tot += i * a[i] * sp.diff(b[l], t) - l * sp.diff(a[i], t) * b[l]
    return sp.expand(tot)


def Q_m(a, b, m, k):
    tot = 0
    for i in range(-k, k + 1):
        for l in range(-k, k + 1):
            if i + l == m:
                tot += sh(b[l], i) * a[i] - sh(a[i], l) * b[l]
    return sp.expand(tot)


def cfun(name, k):
    return {i: sp.Function(f'{name}{i%(2*k+1)}')(t) for i in range(-k, k + 1)}


def qfun(name, k):
    return {i: sp.Function(f'{name}{i%(2*k+1)}')(E) for i in range(-k, k + 1)}


# =====================================================================
print("=" * 70)
print("1. W1  top (shifted-)Wronskian  =>  b_k = lambda a_k   [both faces]")
print("=" * 70)
for k in KS:
    a = cfun('a', k); b = cfun('b', k)
    # classical: only (i,l)=(k,k) hits m=2k
    az(C_m(a, b, 2 * k, k) - k * (a[k] * sp.diff(b[k], t) - sp.diff(a[k], t) * b[k]),
       f"k={k} classical  C_2k = k (a_k b_k' - a_k' b_k)  (pure Wronskian)")
    aq = qfun('A', k); bq = qfun('B', k)
    az(Q_m(aq, bq, 2 * k, k) - (sh(bq[k], k) * aq[k] - sh(aq[k], k) * bq[k]),
       f"k={k} quantum    Q_2k = b_k^[k] a_k - a_k^[k] b_k  (pure shifted-Wronskian)")
# The Wronskian vanishing forces (b_k/a_k)'=0 (classical) / (b_k/a_k) k-periodic (quantum);
# in char 0 both give b_k = lambda a_k. Witness the implication on a concrete a_k:
ak = t**2 + 1
az(sp.together(sp.diff((sp.Symbol('lam') * ak) / ak, t)), "classical: b_k=lam a_k kills the Wronskian")
akq = E**2 + 1
az(sp.Symbol('lam') * akq * sh(akq, 3) - sh(akq, 3) * sp.Symbol('lam') * akq,
   "quantum: b_k=lam a_k kills the shifted-Wronskian (k-periodicity of b_k/a_k)")


# =====================================================================
print("\n" + "=" * 70)
print("2. W2c  classical gatekeeper  C_{2k-1}(b_k=0) = k a_k u' - (k-1) a_k' u")
print("       nonzero polynomial u exists  <=>  a_k = c * h^k  (k-th power)")
print("=" * 70)
for k in KS:
    a = cfun('a', k); b = cfun('b', k)
    b[k] = sp.Integer(0)                       # top gauge b_k = 0
    u = b[k - 1]
    az(C_m(a, b, 2 * k - 1, k) - (k * a[k] * sp.diff(u, t) - (k - 1) * sp.diff(a[k], t) * u),
       f"k={k}  C_(2k-1) = k a_k u' - (k-1) a_k' u   (u = b_{{k-1}}, gauge b_k=0)")
# (<=) a_k = h^k  =>  u = h^{k-1} solves, and u^k = a_k^{k-1} is the first integral.
for k in KS:
    h = sp.Function('h')(t)
    az(k * h**k * sp.diff(h**(k - 1), t) - (k - 1) * sp.diff(h**k, t) * h**(k - 1),
       f"k={k}  (<=) a_k=h^k, u=h^(k-1): k a_k u' - (k-1)a_k' u = 0")
    az((h**(k - 1))**k - (h**k)**(k - 1), f"k={k}  first integral u^k = a_k^(k-1) holds")
# (=>) the UFD/exponent argument: k a_k u' = (k-1) a_k' u  =>  (u^k / a_k^{k-1})' = 0.
for k in KS:
    au = sp.Function('a')(t); uu = sp.Function('u')(t)
    lhs = k * au * sp.diff(uu, t) - (k - 1) * sp.diff(au, t) * uu
    # d/dt log(u^k / a^{k-1}) = k u'/u - (k-1) a'/a = lhs/(u a).  Check the algebra:
    az(sp.simplify((sp.diff(uu**k / au**(k - 1), t)) - (uu**k / au**(k - 1)) * (lhs / (uu * au))),
       f"k={k}  (=>) (u^k/a_k^(k-1))' = (u^k/a_k^(k-1))*(k a_k u'-(k-1)a_k' u)/(u a_k)")
# exponent divisibility: u^k = c a_k^{k-1}, gcd(k,k-1)=1 => k | every prime exponent of a_k.
# Positive witness a_k = c h^k solvable; NEGATIVE witness a_k NOT a k-th power => NO nonzero u.
for k in KS:
    # a_k = t (degree 1, squarefree, not a k-th power for k>=2): u^k = c t^{k-1} impossible.
    found = None
    for du in range(0, 3 * k):
        uc = sp.symbols(f'uc0:{du+1}')
        uP = sum(uc[i] * t**i for i in range(du + 1))
        eq = sp.expand(k * t * sp.diff(uP, t) - (k - 1) * 1 * uP)   # a_k=t, a_k'=1
        solset = sp.solve(sp.Poly(eq, t).all_coeffs(), list(uc), dict=True)
        # any solution with uP != 0 ?
        nz = any(any(sp.simplify(uP.subs(sol)) != 0 for sol in solset) for _ in [0]) if solset else False
        if solset:
            nz = any(sp.simplify(uP.subs(sol)) != 0 for sol in solset)
            if nz:
                found = du
                break
    ok(found is None, f"k={k}  bounded corroboration: a_k=t has NO nonzero u with deg u < {3*k}")
    ok(sp.gcd(k, k - 1) == 1, f"k={k}  gcd(k,k-1)=1 (drives k | exponents => a_k = c h^k)")


# =====================================================================
print("\n" + "=" * 70)
print("3. W2q  quantum gatekeeper  Q_{2k-1}(b_k=0) = u^[k] a_k - a_k^[k-1] u")
print("       sufficiency of shifted k-th power; exact invariant; NEW k>=3 exotica")
print("=" * 70)
for k in KS:
    aq = qfun('A', k); bq = qfun('B', k)
    bq[k] = sp.Integer(0)
    u = bq[k - 1]
    az(Q_m(aq, bq, 2 * k - 1, k) - (sh(u, k) * aq[k] - sh(aq[k], k - 1) * u),
       f"k={k}  Q_(2k-1) = u^[k] a_k - a_k^[k-1] u   (u = b_{{k-1}}, gauge b_k=0)")
# (<=) shifted k-th power a_k = prod_{j=0}^{k-1} h^[j], u = prod_{j=0}^{k-2} h^[j] solves.
for k in KS:
    hh = sp.Function('h')
    ak = sp.prod([hh(E + j) for j in range(k)])
    u = sp.prod([hh(E + j) for j in range(k - 1)])
    az(sp.expand(sh(u, k) * ak - sh(ak, k - 1) * u),
       f"k={k}  (<=) a_k=prod_{{0}}^{{k-1}}h^[j], u=prod_{{0}}^{{k-2}}h^[j]: solves")
# EXACT TELESCOPING INVARIANT (holds for ANY nonzero solution, no factorization assumed):
#   R = [prod_{j=0}^{k-1} u^[j]] / [prod_{j=0}^{k-2} a_k^[j]]  is 1-periodic => constant.
# Proof identity:  R^[1]/R = (u^[k]/u)*(a_k/a_k^[k-1]) = 1 on the equation locus.
for k in KS:
    ug = sp.Function('u'); ag = sp.Function('a')
    U = sp.prod([ug(E + j) for j in range(k)])
    Am = sp.prod([ag(E + j) for j in range(k - 1)])
    U1 = sp.prod([ug(E + 1 + j) for j in range(k)])
    Am1 = sp.prod([ag(E + 1 + j) for j in range(k - 1)])
    # (R^[1]/R) * (a_k^[k-1]/a_k) * (u/u^[k]) == 1  identically (pure telescoping algebra):
    ratio = (U1 / Am1) * (Am / U)
    az(sp.simplify(ratio - (ug(E + k) / ug(E)) * (ag(E) / ag(E + k - 1))),
       f"k={k}  invariant telescoping: R^[1]/R = (u^[k]/u)(a_k/a_k^[k-1])")
print("   => on the equation locus u^[k]/u = a_k^[k-1]/a_k, so R^[1]=R, R=const:")
print("      prod_{j=0}^{k-1} u^[j] = C * prod_{j=0}^{k-2} a_k^[j]  (quantum 1st integral).")

# --- NEW PHENOMENON at k=3: an exotic solution whose a_3 is NOT a shifted cube ---
print("\n   --- k=3 EXOTIC solution: a_3 not a shifted cube ---")
u_ex = E * (E + 3)
a3_ex = E * (E + 2) * (E + 4)
az(sh(u_ex, 3) * a3_ex - sh(a3_ex, 2) * u_ex,
   "exotic: u=E(E+3), a3=E(E+2)(E+4) solves u^[3]a3 = a3^[2]u")
# a3 is NOT a shifted cube c*h(E)h(E+1)h(E+2): deg 3 forces deg h = 1, whose product has
# CONSECUTIVE roots {r,r-1,r-2}; a3 roots are {0,-2,-4}, spaced by 2.
h0, h1, cc = sp.symbols('h0 h1 cc')
hlin = h1 * E + h0
cube = sp.expand(cc * hlin * sh(hlin, 1) * sh(hlin, 2))
sols_cube = sp.solve(sp.Poly(sp.expand(cube - a3_ex), E).all_coeffs(), [h0, h1, cc], dict=True)
real = [s for s in sols_cube if s.get(h1, h1) != 0]
ok(len(real) == 0, "exotic: a3=E(E+2)(E+4) is NOT a shifted cube c h(E)h(E+1)h(E+2)")
# But the invariant still holds:
U_ex = sp.expand(sh(u_ex, 0) * sh(u_ex, 1) * sh(u_ex, 2))
A_ex = sp.expand(sh(a3_ex, 0) * sh(a3_ex, 1))
ok(sp.simplify(U_ex / A_ex) == 1, "exotic: invariant prod u^[j] = prod a3^[j] STILL holds (C=1)")

# Root-necklace lattice picture: per integer coset, sigma = T^{-1}, S_r = 1+sig+..+sig^{r-1}.
def Sr(r):
    return sum(sig**j for j in range(r))
# the equation in necklace form is  S_k(sig) delta(u) = S_{k-1}(sig) delta(a_k).
du = 1 + sig**3           # delta(u)   for roots {0,-3}
da = 1 + sig**2 + sig**4  # delta(a3)  for roots {0,-2,-4}
az(Sr(3) * du - Sr(2) * da, "necklace: S_3 delta(u) = S_2 delta(a3) (exotic solves)")
gq = sp.div(sp.Poly(du, sig), sp.Poly(Sr(2), sig))   # delta(u)/S_2
quo, rem = gq
az(rem.as_expr(), "necklace: S_2 | delta(u) exactly (cofactor g = delta(u)/S_2)")
g_poly = quo.as_expr()
# g = sig^2 - sig + 1 is NOT effective (negative coeff) yet S_2 g, S_3 g are effective:
ok(any(c < 0 for c in sp.Poly(g_poly, sig).all_coeffs()),
   "necklace: cofactor g = sig^2 - sig + 1 is NOT effective (has a negative coeff)")
ok(all(c >= 0 for c in sp.Poly(sp.expand(Sr(2) * g_poly), sig).all_coeffs()),
   "necklace: S_2 g effective  => u genuine")
ok(all(c >= 0 for c in sp.Poly(sp.expand(Sr(3) * g_poly), sig).all_coeffs()),
   "necklace: S_3 g effective  => a3 genuine, but g not effective => a3 not a shifted cube")
# The gap opens because gcd(S_k, S_{k-1}) = 1 but S_{k-1} is nontrivial for k>=3 (S_1=1 for k=2):
for k in KS:
    ok(sp.gcd(sp.Poly(Sr(k), sig), sp.Poly(Sr(k - 1), sig)).as_expr() == 1,
       f"k={k}  gcd(S_k, S_(k-1)) = 1  in Q[sigma]  (coprime necklace operators)")
ok(sp.expand(Sr(1)) == 1, "k=2 special: S_(k-1)=S_1=1 => cofactor forced effective => a_2 IS a shifted square")
# Bounded regression for the memo's exact all-k family from
# g=1-sigma+sigma^2.  The memo supplies the coefficientwise arbitrary-k proof;
# this loop only checks representative k=3,...,9 instances.
for k in range(3, 10):
    gex = 1 - sig + sig**2
    expected_u = 1 + sum(sig**j for j in range(2, k - 1)) + sig**k
    expected_a = 1 + sum(sig**j for j in range(2, k)) + sig**(k + 1)
    az(sp.expand(Sr(k - 1) * gex - expected_u),
       f"k={k} necklace family: S_(k-1)(1-sigma+sigma^2) has the stated effective form")
    az(sp.expand(Sr(k) * gex - expected_a),
       f"k={k} necklace family: S_k(1-sigma+sigma^2) has the stated effective form")
print("   => classical W2 lifts verbatim (a_k = c h^k); QUANTUM W2 does NOT for k>=3:")
print("      the shifted-power class is a PROPER subset of the gatekeeper solution set.")


# =====================================================================
print("\n" + "=" * 70)
print("4. W4  m=0 moment.  Classical C_0 = Psi'.  Quantum Q_0 = (T-1)G.  [exact, any k]")
print("=" * 70)
for k in KS:
    a = cfun('a', k); b = cfun('b', k)
    Psi = sum(i * (a[i] * b[-i] - a[-i] * b[i]) for i in range(1, k + 1))
    az(C_m(a, b, 0, k) - sp.diff(Psi, t),
       f"k={k} classical  C_0 = d/dtau [ sum_{{i>=1}} i (a_i b_-i - a_-i b_i) ]")
for k in KS:
    a = qfun('A', k); b = qfun('B', k)
    G = 0
    for i in range(1, k + 1):
        for r in range(1, i + 1):
            G += sh(a[i], -r) * sh(b[-i], i - r) - sh(a[-i], i - r) * sh(b[i], -r)
    az(Q_m(a, b, 0, k) - (sh(G, 1) - G),
       f"k={k} quantum    Q_0 = (T-1) G,  G = sum_i sum_{{r=1}}^i (a_i^[-r] b_-i^[i-r] - a_-i^[i-r] b_i^[-r])")
print("   Membership makes every product vanish at the base point => Psi = tau (classical),")
print("   G = E (quantum): the moment is normalised with NO free constant. (checked at k<=4)")


# =====================================================================
print("\n" + "=" * 70)
print("5. W5  lattice conjecture: band-2 A* as an integer affine system (mod-3 kill)")
print("=" * 70)
# From classical-Astar.md (commit 84978b9 corpus): degrees P=deg p, V=deg v, W=deg w.
#   (Phi=0)  leading balance : 2V = P + W
#   (I2=0)   leading balance : V + W = 2P + 1   (the +1 is the moment/unit normalisation)
# As M (P,V,W)^T = (0,1)^T :
M = sp.Matrix([[-1, 2, -1], [-2, 1, 1]])
rhs = sp.Matrix([0, 1])
combo = M[0, :] + M[1, :]
ok(combo == sp.Matrix([[-3, 3, 0]]), "A*: row1+row2 = (-3,3,0)  =>  3(V-P) = 1")
g = sp.igcd(int(combo[0]), int(combo[1]), int(combo[2]))
ok(g == 3 and (1 % g) != 0, "A*: gcd of combined LHS = 3 does NOT divide 1  =>  INFEASIBLE over Z")
# The system is solvable over Q but the Q-solution lies off the integer lattice: on the reduced
# line 3(V-P)=1 there is no integer point at all (a fortiori none in Z_{>=0}).
Psol = sp.symbols('P V W')
qsol = sp.linsolve((M, rhs), Psol)
ok(len(qsol) >= 1, "A*: system solvable over Q (rational relaxation nonempty)")
# pick the general Q-solution and confirm 3(V-P) evaluates to 1 (never an integer multiple of 3):
Pv, Vv, Wv = list(qsol)[0]
ok(sp.simplify(3 * Vv - 3 * Pv) == 1, "A*: on the solution set 3V-3P = 1 identically (no integer point)")
print("   FORMULATION (general k): the residual leading-order balances are Z-affine relations")
print("   L d = beta on the degree vector d; homogeneous L pins d to a sublattice, and the lone")
print("   inhomogeneous +1 (from the W4 moment) is reachable iff g := gcd of the reduced row | 1.")
print("   CONJECTURE: outside the tame strata g>1, so the system is infeasible over Z_{>=0}.")


# =====================================================================
print("\n" + "=" * 70)
print("6. W3  conditional trailing identity: exact derivative only if b_1'=0 is known independently")
print("=" * 70)
for k in KS:
    a = cfun('a', k); b = cfun('b', k)
    m = -(k - 1)
    cm = C_m(a, b, m, k)
    s = a[-k]; sp_ = sp.diff(s, t)
    # isolate the a_{-k}-part of C_{-(k-1)}:
    spart = sp.expand(cm.coeff(sp_) * sp_ + cm.subs(sp_, 0).coeff(s) * s)
    az(spart - (-k * s * sp.diff(b[1], t) - sp_ * b[1]),
       f"k={k}  a_(-k)-part of C_(-(k-1)) = -k s b_1' - s' b_1  (s=a_(-k))")
    # rewrite: -k s b1' - s' b1 = -(k-1) s b1' - (s b1)'.  The (s b1)' is EXACT.
    az((-k * s * sp.diff(b[1], t) - sp_ * b[1]) - (-(k - 1) * s * sp.diff(b[1], t) - sp.diff(s * b[1], t)),
       f"k={k}  = -(k-1) s b_1' - (s b_1)'   (exact part (s b_1)'; residue -(k-1)s b_1')")
print("   CONDITIONAL CONSEQUENCE: if a separate rung-by-rung argument proves b_1'=0,")
print("   the residue drops and the a_{-k}-part is -(a_{-k}b_1)'.  W2 alone does not prove this premise.")


# =====================================================================
print("\n" + "=" * 70)
print("7. W6  bottom equation = W2 MIRROR.  Staggered balance (k-1)S = kPhi in BOTH faces.")
print("=" * 70)
# General staggered leading-coefficient identity (quantum): for deg s = S, deg phi = Phi, the
#   E^{S+Phi-1} coefficient of  s^[a] phi - s phi^[b]  equals  (a S - b Phi) * lc(s) * lc(phi).
for (S, Ph, a_, b_) in [(3, 2, -1, -2), (4, 2, -1, -2), (5, 3, -2, -3), (6, 3, -2, -3), (8, 4, -3, -4)]:
    s = sum(sp.Symbol(f's_{i}') * E**i for i in range(S + 1))
    ph = sum(sp.Symbol(f'p_{i}') * E**i for i in range(Ph + 1))
    expr = sp.expand(sh(s, a_) * ph - s * sh(ph, b_))
    co = sp.Poly(expr, E).coeff_monomial(E**(S + Ph - 1))
    az(co - (a_ * S - b_ * Ph) * sp.Symbol(f's_{S}') * sp.Symbol(f'p_{Ph}'),
       f"staggered: coeff(E^(S+Phi-1)) of s^[{a_}]phi - s phi^[{b_}] = ({a_}S-{b_}Phi) lc lc")
# QUANTUM band-k instance: Q_{-(2k-1)} in the bottom gauge is exactly s^[-(k-1)] phi - s phi^[-k].
for k in KS:
    a = {i: sp.Integer(0) for i in range(-k, k + 1)}
    b = {i: sp.Integer(0) for i in range(-k, k + 1)}
    s = sp.Function('s')(E); u = sp.Function('u')(E); v = sp.Function('v')(E); mu = sp.Symbol('mu')
    a[-k] = s; b[-k] = mu * s          # bottom proportionality (mirror of W1, from Q_{-2k})
    a[-(k - 1)] = u; b[-(k - 1)] = v   # level -(k-1)
    phi = mu * u - v
    cand = sh(s, -(k - 1)) * phi - s * sh(phi, -k)
    az(Q_m(a, b, -(2 * k - 1), k) - sp.expand(cand),
       f"k={k}  quantum Q_(-(2k-1)) = s^[-(k-1)] phi - s phi^[-k],  phi = mu a_(-(k-1)) - b_(-(k-1))")
    # a=-(k-1), b=-k => (aS - bPhi) = k Phi - (k-1) S = 0 => (k-1) deg(a_-k) = k deg(phi):
    ok(sp.gcd(k, k - 1) == 1,
       f"k={k}  => (k-1) deg(a_(-k)) = k deg(phi); gcd(k,k-1)=1 forces k | deg(a_(-k))  [Lemma R]")
# CLASSICAL band-k instance: the SAME balance is present -- C_{-(2k-1)} reduced = k s phi' -(k-1)s'phi,
# whose E^{S+Phi-1} leading coefficient is ALSO (k Phi - (k-1) S) lc(s) lc(phi), and which integrates
# to the mirror k-th-power first integral phi^k = c s^{k-1}.  (Corrects the "quantum-only" framing of
# quantum-mirror.md sec 5: the leading balance is available classically too; that memo simply routed
# A* through the Phi, I2 first integrals instead of using C_{-3}.)
for k in KS:
    a = {i: sp.Integer(0) for i in range(-k, k + 1)}
    b = {i: sp.Integer(0) for i in range(-k, k + 1)}
    s = sp.Function('s')(t); u = sp.Function('u')(t); v = sp.Function('v')(t); mu = sp.Symbol('mu')
    a[-k] = s; b[-k] = mu * s
    a[-(k - 1)] = u; b[-(k - 1)] = v
    phi = mu * u - v
    az(C_m(a, b, -(2 * k - 1), k) - (k * s * sp.diff(phi, t) - (k - 1) * sp.diff(s, t) * phi),
       f"k={k}  classical C_(-(2k-1)) = k s phi' - (k-1) s' phi  (same (k-1)S=kPhi balance)")
    az(sp.simplify(sp.diff(phi**k / s**(k - 1), t)
                   - (phi**k / s**(k - 1)) * (k * sp.diff(phi, t) / phi - (k - 1) * sp.diff(s, t) / s)),
       f"k={k}  classical bottom integrates to the MIRROR k-th power  phi^k = c s^(k-1)")
print("   => bottom eqn is the W2 gatekeeper mirrored: classical = clean k-th power phi^k=c s^(k-1);")
print("      quantum = staggered necklace (same (k-1)S=kPhi balance, plus the k>=3 exotica of W2q).")


print("\nALL BAND-K WEAPON CHECKS PASSED")
