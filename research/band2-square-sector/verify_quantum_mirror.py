#!/usr/bin/env python3
"""
verify_quantum_mirror.py
========================
Exact SymPy verification backing the memo `quantum-mirror.md`:
the QUANTUM shifted-square sector of two-sided band 2, resistant branch

        h = 1,   kappa != 0,   mu != 0,   s = a_{-2} != 0,

is EMPTY (no genuine Weyl-membership pair) at ARBITRARY degree.

Convention (quantum memo `quantum-shifted-square-sector-partial.md`):
    A_1[x^{-1}] = (+)_{k in Z} x^k C[E],   (x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E),
    f^[r](E) = f(E+r),
    ladder-m coefficient of [D,X]:  Q_m = sum_{k+l=m} ( b_l^[k] a_k - a_k^[l] b_l ),
    [D,X] = 1  <=>  Q_m = delta_{m0}.

Branch data after the audited positive cascade (memo (4.x)):
    a_2 = 1, a_1 = p, a_0, a_{-1}=u, a_{-2}=s ;
    b_2 = 0, b_1 = kappa, b_0 = B, b_{-1}=v, b_{-2}=w = mu*s ;
    memberships  E | u,v   and   E(E-1) | s,w.

The script verifies, in order:
  0. Q_m reproduce the memo tail equations (4.1),(4.3)-(4.6) for general h.
  1. The three exact arbitrary-degree integrals:
        midpoint      B^[1]+B = kappa p + gamma            (ladder 2)
        ladder-1      kappa a_0 - (v^[1]+v) - (B^2-gamma B)/kappa = A
        central       w^[1]+w + p^[-1]v - kappa u = E       (point conds => const = 0)
  2. Structural lemmas (leading coefficients) driving the emptiness proof:
        Q_{-3} = s^[-1] phi - s phi^[-2],  phi = mu u - v ;  coeff(E^{S+F-1}) = s_S phi_F (2F-S)
        Q_{-2} = s*Psi + (v^[-1]u - u^[-1]v),  coeff(E^{U+V-1}) of bracket = (U-V) u_U v_V
        Q_{-1} = mu(s^[1]p - p^[-2]s) + kappa(s-s^[1]) + v(a0-a0^[-1]) + u(B^[-1]-B),
                 coeff(E^{S+P-1}) of first term = (S+2P) s_S p_P
        coeff(E^{A0-1}) of (a0 - a0^[-2]) = 2 A0 lc(a0)
  3. Case closers:
        phi=0 collapse ;  P=0 case (a) ;  P=0 case (b) ;  Q_{-2} T1 dominance ;
        locus L (S=3P,V=2P) killed by Q_{-1} leading term.
  4. Bounded-degree EXHAUSTIVE emptiness (Groebner + Rabinowitsch saturation of
     kappa*mu != 0 and s != 0): a curated grid, all EMPTY.

Run:  uv run --with sympy python verify_quantum_mirror.py
Ends: ALL QUANTUM MIRROR CHECKS PASSED
"""
import sympy as sp

E = sp.symbols('E')
kappa, gamma, mu, Aconst = sp.symbols('kappa gamma mu A')


def sh(f, r):
    return sp.expand(f.subs(E, E + r))


def poly(name, deg):
    cs = sp.symbols(f'{name}_0:{deg+1}')
    return sum(cs[i] * E**i for i in range(deg + 1)), list(cs)


def Qm(a, b, m):
    """a,b: dicts k-> coeff poly, k in -2..2. Ladder-m coefficient of [D,X]."""
    return sp.expand(sum(
        sh(b[l], k) * a[k] - sh(a[k], l) * b[l]
        for k in range(-2, 3) for l in range(-2, 3) if k + l == m))


def az(expr, label):
    if sp.expand(expr) != 0:
        raise AssertionError(label + " : " + str(sp.expand(expr)))
    print("PASS", label)


# =====================================================================
# 0. Q_m reproduce the memo's displayed tail equations (general h)
# =====================================================================
print("--- Part 0: tail equations match the audited quantum memo ---")
h, _ = poly('h', 3)
p, _ = poly('p', 3)
a0, _ = poly('a0', 3)
u, _ = poly('u', 3)        # a_{-1}
s, _ = poly('s', 3)        # a_{-2}
B, _ = poly('B', 3)        # b_0
v, _ = poly('v', 3)        # b_{-1}
w, _ = poly('w', 3)        # b_{-2}
a = {2: sp.expand(h * sh(h, 1)), 1: sp.expand(h * p), 0: a0, -1: u, -2: s}
b = {2: sp.Integer(0), 1: sp.expand(kappa * h), 0: B, -1: v, -2: w}

memo41 = (sh(h, 1) * sh(v, 2) - sh(h, -1) * v + p * (sh(B, 1) - B) + kappa * (a0 - sh(a0, 1)))
memo43 = (sh(w, -1) * u + sh(v, -2) * s - sh(u, -2) * w - sh(s, -1) * v)
memo44 = (w * a0 + sh(v, -1) * u + sh(B, -2) * s - sh(a0, -2) * w - sh(u, -1) * v - s * B)
memo45 = (sh(w, 1) * a[1] + v * a0 + sh(B, -1) * u + kappa * sh(h, -2) * s
          - sh(a[1], -2) * w - sh(a0, -1) * v - u * B - kappa * sh(s, 1) * h)
memo46 = (sh(w, 2) * a[2] + sh(v, 1) * a[1] + kappa * sh(h, -1) * u
          - sh(a[2], -2) * w - sh(a[1], -1) * v - kappa * sh(u, 1) * h)
memo42 = (sh(w, -2) * s - sh(s, -2) * w)

az(Qm(a, b, 1) - h * memo41, "Q_1 = h * memo(4.1)")                # ladder 1 carries a factor h
az(Qm(a, b, -3) - memo43,    "Q_{-3} = memo(4.3)")
az(Qm(a, b, -2) - memo44,    "Q_{-2} = memo(4.4)")
az(Qm(a, b, -1) - memo45,    "Q_{-1} = memo(4.5)")
az(Qm(a, b, 0) - memo46, "Q_0 = memo(4.6) LHS  (central normalization = 1)")
az(Qm(a, b, -4) - memo42,    "Q_{-4} = memo(4.2)")
az(Qm(a, b, 4), "Q_4 = 0 identically (gauge b2=0)")
az(Qm(a, b, 3), "Q_3 = 0 identically (b1 = kappa h)")


# =====================================================================
# 1. Exact arbitrary-degree integrals  (branch h = 1)
# =====================================================================
print("\n--- Part 1: exact integrated identities (h=1) ---")
a1 = {2: sp.Integer(1), 1: p, 0: a0, -1: u, -2: s}
b1 = {2: sp.Integer(0), 1: kappa, 0: B, -1: v, -2: w}

# midpoint rewrite  p (B^[1]-B) = (1/kappa)(T-1)(B^2 - gamma B)  using B^[1]+B = kappa p + gamma
az(sh(B**2 - gamma*B, 1) - (B**2 - gamma*B) - (sh(B, 1)-B)*(sh(B, 1)+B-gamma),
   "(T-1)(B^2-gamma B) = (B^[1]-B)(B^[1]+B-gamma)")

# ladder-1 integral: under midpoint, memo(4.1)|_{h=1} = -(T-1)[ kappa a0 - (v^[1]+v) - (B^2-gamma B)/kappa ]
p_mid = (sh(B, 1) + B - gamma) / kappa
memo41_h1 = (sh(v, 2) - v + p_mid*(sh(B, 1) - B) + kappa*(a0 - sh(a0, 1)))
Lpot = kappa*a0 - (sh(v, 1) + v) - (B**2 - gamma*B)/kappa
az(memo41_h1 + (sh(Lpot, 1) - Lpot),
   "ladder-1 integral: kappa a0 - (v^[1]+v) - (B^2-gamma B)/kappa = const (A)")

# central telescoping:  Q_0 = (T-1) G,  G = w^[1]+w + p^[-1] v - kappa u
G = sh(w, 1) + w + sh(p, -1)*v - kappa*u
az((sh(G, 1) - G) - Qm(a1, b1, 0),
   "central: Q_0 = (T-1)[ w^[1]+w + p^[-1]v - kappa u ]  => G = E + c0")

# point conditions => c0 = 0.  With s=E(E-1)sig, v=E nu, w=mu s :  G(0)=0.
sig_, _ = poly('sig', 2)
nu_, _ = poly('nu', 2)
s_m = sp.expand(E*(E-1)*sig_)
v_m = sp.expand(E*nu_)
w_m = mu*s_m
u_m = sp.expand((mu*(sh(s_m, 1) + s_m) + sh(p, -1)*v_m - E)/kappa)   # from central, solved for u
G_m = sh(w_m, 1) + w_m + sh(p, -1)*v_m - kappa*u_m
az(G_m.subs(E, 0), "point conditions s(0)=s(1)=v(0)=0 => G(0)=0 => c0=0")
az(G_m - E, "central identity holds exactly: w^[1]+w+p^[-1]v-kappa u = E")
az(u_m.subs(E, 0), "membership E|u is automatic (u(0)=0)")


# =====================================================================
# 2. Structural lemmas (leading coefficients)
# =====================================================================
print("\n--- Part 2: structural leading-coefficient lemmas ---")

# Q_{-3} = s^[-1] phi - s phi^[-2] ;  coeff(E^{S+F-1}) = s_S phi_F (2F-S)   (F = deg phi)
for (S, F) in [(2, 1), (3, 1), (4, 2), (2, 2), (5, 2), (6, 3)]:
    sc = sp.symbols(f'sA_0:{S+1}'); ss = sum(sc[i]*E**i for i in range(S+1))
    fc = sp.symbols(f'fA_0:{F+1}'); ff = sum(fc[i]*E**i for i in range(F+1))
    Q = sp.expand(sh(ss, -1)*ff - ss*sh(ff, -2))
    az(sp.Poly(Q, E).nth(S+F-1) - sc[S]*fc[F]*(2*F - S),
       f"Q_-3 leading (S={S},F={F}): coeff = s_S phi_F (2F-S)  => phi!=0 forces S=2F")

# Q_{-2} regroup and bracket coeff
pg, _ = poly('pg', 3); a0g, _ = poly('a0g', 3); ug, _ = poly('ug', 3); sg, _ = poly('sg', 3)
Bg, _ = poly('Bg', 3); vg, _ = poly('vg', 3); wg = mu*sg
ag = {2: sp.Integer(1), 1: pg, 0: a0g, -1: ug, -2: sg}
bg = {2: sp.Integer(0), 1: kappa, 0: Bg, -1: vg, -2: wg}
Psi = mu*(a0g - sh(a0g, -2)) + (sh(Bg, -2) - Bg)
az(Qm(ag, bg, -2) - (sg*Psi + (sh(vg, -1)*ug - sh(ug, -1)*vg)),
   "Q_-2 = s*Psi + (v^[-1]u - u^[-1]v),  Psi = mu(a0-a0^[-2])+(B^[-2]-B)")
for (Ud, Vd) in [(3, 3), (4, 2), (2, 4), (5, 3)]:
    uc = sp.symbols(f'uB_0:{Ud+1}'); uu = sum(uc[i]*E**i for i in range(Ud+1))
    vc = sp.symbols(f'vB_0:{Vd+1}'); vv = sum(vc[i]*E**i for i in range(Vd+1))
    az(sp.Poly(sp.expand(sh(vv, -1)*uu - sh(uu, -1)*vv), E).nth(Ud+Vd-1) - (Ud-Vd)*uc[Ud]*vc[Vd],
       f"(v^[-1]u-u^[-1]v) coeff(E^{Ud+Vd-1}) = (U-V)u_U v_V   (=0 iff U=V) [U={Ud},V={Vd}]")

# Q_{-1} regroup and leading of s^[1]p - p^[-2]s
az(Qm(ag, bg, -1) - (mu*(sh(sg, 1)*pg - sh(pg, -2)*sg) + kappa*(sg - sh(sg, 1))
                     + vg*(a0g - sh(a0g, -1)) + ug*(sh(Bg, -1) - Bg)),
   "Q_-1 = mu(s^[1]p - p^[-2]s) + kappa(s-s^[1]) + v(a0-a0^[-1]) + u(B^[-1]-B)")
for (S, P) in [(6, 2), (3, 1), (4, 2), (9, 3), (5, 2)]:
    sc = sp.symbols(f'sC_0:{S+1}'); ss = sum(sc[i]*E**i for i in range(S+1))
    pc = sp.symbols(f'pC_0:{P+1}'); pp = sum(pc[i]*E**i for i in range(P+1))
    az(sp.Poly(sp.expand(sh(ss, 1)*pp - sh(pp, -2)*ss), E).nth(S+P-1) - (S+2*P)*sc[S]*pc[P],
       f"(s^[1]p - p^[-2]s) coeff(E^{S+P-1}) = (S+2P)s_S p_P (!=0, P>=1) [S={S},P={P}]")

# (a0 - a0^[-2]) leading = 2 A0 lc(a0)
for A0 in [1, 2, 3, 4]:
    ac = sp.symbols(f'aC_0:{A0+1}'); aa = sum(ac[i]*E**i for i in range(A0+1))
    az(sp.Poly(sp.expand(aa - sh(aa, -2)), E).nth(A0-1) - 2*A0*ac[A0],
       f"(a0-a0^[-2]) coeff(E^{A0-1}) = 2 A0 lc(a0) [A0={A0}]")


# =====================================================================
# 3. Case closers
# =====================================================================
print("\n--- Part 3: case closers ---")

# helper: build reduced data for chosen degrees (p possibly constant)
def build_reduced(dp, dv, ds, pconst=False):
    if pconst:
        pp = sp.symbols('p0'); Bdeg = 0
    else:
        pc = sp.symbols(f'p1_0:{dp+1}'); pp = sum(pc[i]*E**i for i in range(dp+1)); Bdeg = dp
    nu = sp.symbols(f'n1_0:{dv}'); vv = sp.expand(E*sum(nu[i]*E**i for i in range(dv)))
    sig = sp.symbols(f'g1_0:{ds-1}'); ss = sp.expand(E*(E-1)*sum(sig[i]*E**i for i in range(ds-1)))
    bc = sp.symbols(f'B1_0:{Bdeg+1}'); BB = sum(bc[i]*E**i for i in range(Bdeg+1))
    solB = sp.solve(sp.Poly(sp.expand(sh(BB, 1)+BB-(kappa*pp+gamma)), E).all_coeffs(), bc, dict=True)[0]
    BB = sp.expand(BB.subs(solB))
    a0v = sp.expand((sh(vv, 1)+vv)/kappa + (BB**2-gamma*BB)/kappa**2 + Aconst/kappa)
    uu = sp.expand((mu*(sh(ss, 1)+ss)+sh(pp, -1)*vv - E)/kappa)
    aa = {2: sp.Integer(1), 1: pp, 0: a0v, -1: uu, -2: ss}
    bb = {2: sp.Integer(0), 1: kappa, 0: BB, -1: vv, -2: mu*ss}
    return dict(p=pp, v=vv, s=ss, B=BB, a0=a0v, u=uu, sig=list(sig), nu=list(nu),
                Q=lambda m: Qm(aa, bb, m))

# self-consistency: derived B,a0,u make Q_2=0, Q_1=0, Q_0=1 hold IDENTICALLY, so the
# full system [D,X]=1 is faithfully reduced to Q_{-1}=Q_{-2}=Q_{-3}=0 (plus w=mu s).
Msc = build_reduced(3, 3, 4)
az(Msc['Q'](2),     "reduction faithful: Q_2 = 0 identically (midpoint => ladder 2)")
az(Msc['Q'](1),     "reduction faithful: Q_1 = 0 identically (ladder-1 integral)")
az(Msc['Q'](0) - 1, "reduction faithful: Q_0 = 1 identically (central integral)")
az(Msc['Q'](-4),    "reduction faithful: Q_{-4} = 0 identically (w = mu s)")

# --- phi = 0  =>  v^[-1]u - u^[-1]v = 0, so Q_{-2}=s*Psi ---
uT, _ = poly('uT', 3)
az(sh(mu*uT, -1)*uT - sh(uT, -1)*(mu*uT),
   "phi=0 (v=mu u): v^[-1]u - u^[-1]v = 0  => Q_-2 = s*Psi => Psi=0")
#   Psi=0 forces deg a0 <= P (RHS B-B^[-2] has degree P-1); then u polynomial gives deg u = S-P,
#   and the Q_-1 leading term (S+2P)mu s_S p_P dominates -> contradiction (checked generally above).
#   For P=0: Psi=0 => a0 2-periodic => a0 const => v const => v=0 => u=0 => central mu(s^[1]+s)=E,
#   impossible because deg(s^[1]+s) = deg s = S >= 2 while RHS E has degree 1.
for ds in (2, 3, 4):
    sig0 = sp.symbols(f'sig0x{ds}_0:{ds-1}')
    s0 = sp.expand(E*(E-1)*sum(sig0[i]*E**i for i in range(ds-1)))
    d = sp.Poly(mu*(sh(s0, 1) + s0) - E, E).degree()   # equals ds >= 2, so != degree-1 poly E
    if d != ds or ds < 2:
        raise AssertionError(f"mu(s^[1]+s)-E degree {d} != s-degree {ds}")
    print(f"PASS central-with-u=v=0 impossible: deg(mu(s^[1]+s)-E)={d}=S>=2 (S={ds})")

# --- P=0, case (a): V=S. coeff(E^{2S-1}) of Q_-2 = 4 mu S s_S v_S / kappa (!=0) ---
for S in (2, 3):
    M = build_reduced(0, S, S, pconst=True)   # p const, deg v = deg s = S
    az(sp.Poly(M['Q'](-2), E).nth(2*S-1) - 4*mu*S*M['sig'][S-2]*M['nu'][S-1]/kappa,
       f"P=0 case(a) S={S}: coeff(E^{2*S-1}) Q_-2 = 4 mu S s_S v_S/kappa (!=0)")

# --- P=0, case (b): mu*p0 = kappa  =>  phi = (mu/kappa)(mu(s^[1]+s)-E), deg phi = S => Q_-3: S=2S ---
p0 = kappa/mu
sigb = sp.symbols('gb_0:2'); sb = sp.expand(E*(E-1)*sum(sigb[i]*E**i for i in range(2)))
nub = sp.symbols('nb_0:2'); vb = sp.expand(E*sum(nub[i]*E**i for i in range(2)))
ub = sp.expand((mu*(sh(sb, 1)+sb)+p0*vb - E)/kappa)
az(mu*ub - vb - (mu/kappa)*(mu*(sh(sb, 1)+sb) - E),
   "P=0 case(b): phi = (mu/kappa)(mu(s^[1]+s)-E), deg phi = deg s => Q_-3 forces S=2S, contra")

# --- Q_-2 T1 dominance (2P>V): top coeff = 2 mu s_S p_P^2 (!=0) ---
M = build_reduced(2, 1, 3)   # P=2, V=1, S=3 ; 2P=4>V=1 ; P+V=S
az(sp.Poly(M['Q'](-2), E).degree() - 6, "Q_-2 degree = S+2P-1 = 6 (T1 dominant)")
az(sp.Poly(M['Q'](-2), E).LC() - 2*mu*M['sig'][1]*M['p'].as_poly(E).nth(2)**2,
   "Q_-2 top coeff = 2 mu s_S p_P^2 (!=0) when 2P>V")

# --- Locus L (S=3P, V=2P, lc(a0)=0, 2 mu s_S+p_P v_V=0): Q_-1 coeff(E^{4P-1}) != 0 ---
P, S, V = 2, 6, 4
M = build_reduced(P, V, S)
sS = M['sig'][S-2]; pP = M['p'].as_poly(E).nth(P); vV = M['nu'][V-1]
solL = sp.solve([2*mu*sS + pP*vV, pP**2/4 + 2*vV/kappa], [vV, sS], dict=True)[0]
c_top = sp.simplify(sp.Poly(sp.expand(M['Q'](-1).subs(solL)), E).nth(S + P - 1))   # E^{S+P-1}=E^{4P-1}=E^7
# hand value: coeff(E^{4P-1}) of Q_-1 in locus L = (S+2P) mu s_S p_P = 5 P mu s_S p_P,
# and with the two L-relations s_S = kappa p_P^3/(16 mu):  = 5 kappa p_P^4 / 8
az(c_top - 5*kappa*pP**4/8, "locus L (P=2): Q_-1 coeff(E^7) = 5 kappa p_P^4/8 (!=0) => L empty")


# =====================================================================
# 4. Bounded-degree EXHAUSTIVE emptiness (Groebner + Rabinowitsch saturation)
# =====================================================================
print("\n--- Part 4: bounded-degree exhaustive emptiness (kappa*mu!=0, s!=0) ---")

def reduced_eqs(dp, dv, ds):
    pc = sp.symbols(f'P_0:{dp+1}'); pp = sum(pc[i]*E**i for i in range(dp+1))
    nu = sp.symbols(f'N_0:{dv}'); vv = sp.expand(E*sum(nu[i]*E**i for i in range(dv))) if dv >= 1 else sp.Integer(0)
    sig = sp.symbols(f'G_0:{ds-1}'); ss = sp.expand(E*(E-1)*sum(sig[i]*E**i for i in range(ds-1))) if ds >= 2 else sp.Integer(0)
    bc = sp.symbols(f'B_0:{dp+1}'); BB = sum(bc[i]*E**i for i in range(dp+1))
    solB = sp.solve(sp.Poly(sp.expand(sh(BB, 1)+BB-(kappa*pp+gamma)), E).all_coeffs(), bc, dict=True)[0]
    BB = sp.expand(BB.subs(solB))
    a0v = sp.expand((sh(vv, 1)+vv)/kappa + (BB**2-gamma*BB)/kappa**2 + Aconst/kappa)
    uu = sp.expand((mu*(sh(ss, 1)+ss)+sh(pp, -1)*vv - E)/kappa)
    aa = {2: sp.Integer(1), 1: pp, 0: a0v, -1: uu, -2: ss}
    bb = {2: sp.Integer(0), 1: kappa, 0: BB, -1: vv, -2: mu*ss}
    eqs = []
    for m in (-1, -2, -3):
        for c in sp.Poly(Qm(aa, bb, m), E).all_coeffs():
            num = sp.expand(sp.fraction(sp.together(c))[0])
            if num != 0:
                eqs.append(num)
    unknowns = list(pc) + list(nu) + list(sig) + [kappa, gamma, mu, Aconst]
    return eqs, unknowns, list(sig)

def empty_branch(dp, dv, ds):
    """True iff resistant branch (kappa*mu!=0, s!=0) empty at these degrees."""
    eqs, unk, sig = reduced_eqs(dp, dv, ds)
    if not sig:
        return True
    y = sp.symbols('y_slack')
    for c in sig:                       # s!=0  <=>  some sig coeff != 0
        polys = list(eqs) + [sp.expand(kappa*mu*c*y - 1)]
        Gb = sp.groebner(polys, *(unk + [y]), order='grevlex')
        if list(Gb.exprs) != [sp.Integer(1)]:
            return False                # nonempty branch found
    return True

grid = [(dp, dv, ds) for ds in (2, 3) for dp in (0, 1, 2) for dv in (1, 2, 3)]
grid += [(0, 2, 4), (1, 2, 4), (2, 3, 4)]     # a few ds=4 profiles
allok = True
for (dp, dv, ds) in grid:
    e = empty_branch(dp, dv, ds)
    allok = allok and e
    print(f"  (dp={dp},dv={dv},ds={ds}): {'EMPTY' if e else 'NONEMPTY!!'}")
if not allok:
    raise AssertionError("a bounded profile was NONEMPTY")
print("bounded exhaustive grid: all EMPTY")

print("\nALL QUANTUM MIRROR CHECKS PASSED")
