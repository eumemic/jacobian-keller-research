#!/usr/bin/env python3
"""
verify_quantum_completion.py
============================
Exact SymPy verification backing `quantum-completion.md`: the COMPLETE
classification of the quantum band-2 shifted-square sector.

Convention (as in verify_quantum_mirror.py):
    A_1[x^{-1}] = (+)_k x^k C[E],  (x^a f)(x^b g) = x^{a+b} f(E+b) g(E),
    f^[r](E) = f(E+r),   E = x*partial,  partial = x^{-1}E,
    Q_m = sum_{k+l=m}( b_l^[k] a_k - a_k^[l] b_l ),  [D,X]=1 <=> Q_m = delta_{m0}.

Sector: a_2 = h h^[1] (shifted square), gauge b_2=0, b_1 = kappa h; genuine
memberships E|u,v and E(E-1)|s,w  (u=a_-1, v=b_-1, s=a_-2, w=b_-2).

Branches (constant h => normalize h=1):
  d.  h nonconstant                    -> EMPTY  (general-h central integral)
  a.  kappa = 0                        -> EMPTY
  c.  s = 0, w != 0                    -> EMPTY
  *.  s != 0, w != 0  (mu != 0)        -> EMPTY   [proved in quantum-mirror.md]
  b.  s != 0, w = 0   (mu = 0)         -> exactly the tame family (c = 1/kappa^2)
  base. s = 0, w = 0                   -> tame family with c_1 = 0
The surviving pairs are exactly the tame automorphism family
  U = x + c0 + c1 partial,  X = U^2 - partial/kappa - A,  D = lambda X + kappa U + beta.

Run:  uv run --with sympy python verify_quantum_completion.py
Ends: ALL QUANTUM COMPLETION CHECKS PASSED
"""
import sympy as sp

E = sp.symbols('E')
kappa, gamma, A, c, lam, mu, alpha = sp.symbols('kappa gamma A c lambda mu alpha')


def sh(f, r):
    return sp.expand(f.subs(E, E + r))


def poly(name, deg):
    cs = sp.symbols(f'{name}_0:{deg+1}')
    return sum(cs[i]*E**i for i in range(deg+1)), list(cs)


def Qm(a, b, m):
    return sp.expand(sum(sh(b[l], k)*a[k] - sh(a[k], l)*b[l]
                         for k in range(-2, 3) for l in range(-2, 3) if k+l == m))


def az(expr, label):
    if sp.expand(expr) != 0:
        raise AssertionError(label + " : " + str(sp.expand(expr)))
    print("PASS", label)


def midB(p, dp):
    bc = sp.symbols(f'MB_0:{dp+1}'); B = sum(bc[i]*E**i for i in range(dp+1))
    s0 = sp.solve(sp.Poly(sp.expand(sh(B, 1)+B-(kappa*p+gamma)), E).all_coeffs(), bc, dict=True)[0]
    return sp.expand(B.subs(s0))


def clearden(eqs):
    return [sp.expand(sp.fraction(sp.together(e))[0]) for e in eqs if sp.expand(e) != 0]


# ==========================================================================
# BRANCH d: nonconstant h -> EMPTY, via the general-h central integral
# ==========================================================================
print("--- Branch d: h nonconstant is EMPTY (general-h central integral) ---")
h, _ = poly('h', 3); p, _ = poly('p', 3); u, _ = poly('u', 3); v, _ = poly('v', 3); w, _ = poly('w', 3)
a2 = sp.expand(h*sh(h, 1)); a1 = sp.expand(h*p)
a = {2: a2, 1: a1, 0: poly('a0', 2)[0], -1: u, -2: poly('s', 2)[0]}
b = {2: sp.Integer(0), 1: kappa*h, 0: poly('B', 2)[0], -1: v, -2: w}
G = sp.expand(sh(h, -1)*(h*sh(w, 1) + sh(h, -2)*w + sh(p, -1)*v - kappa*u))
az(sh(G, 1) - G - Qm(a, b, 0),
   "general-h central integral: Q_0 = (T-1)[ h^[-1]( h w^[1] + h^[-2] w + p^[-1] v - kappa u ) ]")
# G = E + const  =>  h^[-1] | (E + const)  =>  deg h <= 1 (h affine when nonconstant).
# h affine, h = alpha(E - rho): bracket = (E+const)/h^[-1] = 1/alpha (constant).
rho = sp.symbols('rho'); haff = alpha*(E - rho)
# verify (E+const)/h^[-1] is the constant 1/alpha when h^[-1] | (E+const):
#   pick const so that h^[-1] | (E+const): h^[-1] = alpha(E-1-rho), root E=1+rho => const = -(1+rho)
const = -(1 + rho)
quot = sp.cancel((E + const)/sh(haff, -1))
az(quot - 1/alpha, "h affine => (E+const)/h^[-1] = 1/alpha (constant)")
# then kappa u = h w^[1] + h^[-2] w + p^[-1] v - 1/alpha; evaluate at E=0 with w(0)=w(1)=v(0)=0.
# membership: w = E(E-1) wsig, v = E vsig, giving w(0)=w(1)=v(0)=0.
wsig, _ = poly('ws', 2); vsig, _ = poly('vs', 2)
wm = sp.expand(E*(E-1)*wsig); vm = sp.expand(E*vsig)
ku = sp.expand(haff*sh(wm, 1) + sh(haff, -2)*wm + sh(p, -1)*vm - 1/alpha)  # = kappa*u
az(ku.subs(E, 0) - (-1/alpha), "kappa*u(0) = -1/alpha  (w(0)=w(1)=v(0)=0)")
print("  => u(0) = -1/(kappa*alpha) != 0 violates E|u.  BRANCH d EMPTY.\n")

# ==========================================================================
# BRANCH a: kappa = 0 -> EMPTY  (b1=0, b2=0)
# ==========================================================================
print("--- Branch a: kappa = 0 is EMPTY ---")
hA, _ = poly('hA', 3); vA, _ = poly('vA', 3); wA, _ = poly('wA', 3)
a2A = sp.expand(hA*sh(hA, 1))
aA = {2: a2A, 1: poly('a1A', 3)[0], 0: poly('a0A', 2)[0], -1: poly('uA', 2)[0], -2: poly('sA', 2)[0]}
Bconst = sp.symbols('betaA')
bA = {2: sp.Integer(0), 1: sp.Integer(0), 0: Bconst, -1: vA, -2: wA}
az(Qm(aA, bA, 2), "kappa=0: Q_2 = 0 with B const (B 2-periodic => B const)")
# Q_1 = h(h^[1] v^[2] - h^[-1] v); vanishing => h^[-1] v is 2-periodic => constant
az(Qm(aA, bA, 1) - hA*(sh(hA, 1)*sh(vA, 2) - sh(hA, -1)*vA), "kappa=0: Q_1 = h(h^[1]v^[2]-h^[-1]v)")
#   nonconstant h => h^[-1]v const => v=0; h=1 => v 2-periodic & E|v => v=0.  Take v=0.
bA0 = {2: sp.Integer(0), 1: sp.Integer(0), 0: Bconst, -1: sp.Integer(0), -2: wA}
Q0A = Qm(aA, bA0, 0)
PsiA = sp.expand(sh(a2A, -1)*sh(wA, 1) + sh(a2A, -2)*wA)   # = h^[-1](h w^[1] + h^[-2] w)
az(sh(PsiA, 1) - PsiA - Q0A, "kappa=0,v=0: Q_0 = (T-1)[ h^[-1](h w^[1]+h^[-2]w) ]")
az(PsiA - sp.expand(sh(hA, -1)*(hA*sh(wA, 1) + sh(hA, -2)*wA)), "  potential factors through h^[-1]")
#   h=1: w^[1]+w = E => w = E/2 - 1/4, but w(0) must be 0 (E(E-1)|w): contradiction.
w0, w1 = sp.symbols('w0 w1'); wlin = w0 + w1*E
solw = sp.solve(sp.Poly(sp.expand(sh(wlin, 1)+wlin-E), E).all_coeffs(), [w0, w1], dict=True)[0]
az(wlin.subs(solw) - (E/2 - sp.Rational(1, 4)), "kappa=0,h=1: Q_0 forces w = E/2 - 1/4")
if wlin.subs(solw).subs(E, 0) == 0:
    raise AssertionError("w(0) should be nonzero")
print("  h=1: w(0) = -1/4 != 0 violates E(E-1)|w;  nonconstant h: deg(h w^[1]+h^[-2]w)=1+W>=3 != 0.")
print("  BRANCH a EMPTY (all h).\n")

# ==========================================================================
# BRANCH c: s = 0, w != 0 -> EMPTY
# ==========================================================================
print("--- Branch c: s = 0, w != 0 is EMPTY ---")
uc, _ = poly('uc', 3); vc, _ = poly('vc', 3)
az(sh(lam*vc, -1)*vc - sh(vc, -1)*(lam*vc), "Q_-2|s=0 = v^[-1]u - u^[-1]v = 0 => u = lambda v")
# central (s=0): w^[1]+w + p^[-1]v - kappa u = E, w = c u u^[-1], u = lam v.
# LHS = v * bracket with bracket = c lam^2 (v^[1]+v^[-1]) + p^[-1] - kappa lam, so the
# equation reads v*bracket = E, hence v | E, hence v linear (with E|v).
pc3, _ = poly('pc', 3); vgen, _ = poly('vc2', 3)
u_ = lam*vgen; w_ = sp.expand(c*u_*sh(u_, -1))
central = sp.expand(sh(w_, 1) + w_ + sh(pc3, -1)*vgen - kappa*u_)
bracket = sp.expand(c*lam**2*(sh(vgen, 1) + sh(vgen, -1)) + sh(pc3, -1) - kappa*lam)
az(central - vgen*bracket, "central (s=0,u=lam v) = v * (polynomial) => v | E => v linear")
# finite system (v linear, u linear, w quadratic, p affine): Groebner EMPTY (w!=0)
def branch_c_empty(dp, dv):
    p, pcf = poly('P', dp); nu, ncf = poly('N', dv-1); v = sp.expand(E*nu)
    u = lam*v; wv = sp.expand(c*u*sh(u, -1)); B = midB(p, dp)
    a0 = sp.expand((sh(v, 1)+v)/kappa + (B**2-gamma*B)/kappa**2 + A/kappa)
    aa = {2: sp.Integer(1), 1: p, 0: a0, -1: u, -2: sp.Integer(0)}
    bb = {2: sp.Integer(0), 1: kappa, 0: B, -1: v, -2: wv}
    eqs = clearden(sp.Poly(Qm(aa, bb, 0)-1, E).all_coeffs()
                   + sp.Poly(Qm(aa, bb, -1), E).all_coeffs()
                   + sp.Poly(Qm(aa, bb, -2), E).all_coeffs())
    unk = pcf + ncf + [c, lam, gamma, A]
    for nci in ncf:
        y = sp.Symbol('yy')
        Gg = sp.groebner(eqs + [kappa*c*lam*nci*y - 1], *(unk+[y]), order='grevlex')
        if list(Gg.exprs) != [sp.Integer(1)]:
            return False
    return True
for (dp, dv) in [(1, 1), (1, 2), (2, 2)]:
    if not branch_c_empty(dp, dv):
        raise AssertionError(f"branch c NONEMPTY at ({dp},{dv})")
    print(f"PASS branch c EMPTY (w!=0, v!=0) at (dp={dp},dv={dv})")
# zero-poly sub-cases: u=0 => w=c u u^[-1]=0 contradicts w!=0.  v=0 (u!=0):
#   central w^[1]+w-kappa u = E with w=c u u^[-1] gives u*(c(u^[1]+u^[-1])-kappa)=E => u|E => u=u1 E;
#   then E^2 coeff is 2 c u1^2 = 0, impossible (c!=0,u1!=0).
u1 = sp.symbols('u1'); uLin = u1*E; wLin = sp.expand(c*uLin*sh(uLin, -1))
lhs = sp.expand(sh(wLin, 1) + wLin - kappa*uLin)     # central with v=0
az(sp.Poly(lhs, E).nth(2) - 2*c*u1**2, "branch c v=0 sub-case: central E^2-coeff = 2 c u1^2 (=0 impossible)")
print("  BRANCH c EMPTY (u=0 => w=0 contra; v=0 => 2 c u1^2=0 contra; v!=0 => finite, Groebner empty).\n")

# ==========================================================================
# BRANCH b: w = 0, s != 0 -> exactly the tame family (c = 1/kappa^2)
# ==========================================================================
print("--- Branch b: w = 0, s != 0 is exactly the tame family ---")
# Q_-3|w=0 => s v^[-2] = s^[-1] v => s = c v v^[-1]
vb, _ = poly('vb', 3)
sinv = sp.expand(c*vb*sh(vb, -1))
az(sh(sinv, -1)*vb - sh(vb, -2)*sinv, "Q_-3|w=0: s = c v v^[-1] solves s v^[-2] = s^[-1] v")
# leading coeff of Q_-2(b) = -(P/kappa) lc(v)^2 lc(p) (c kappa^2 - 1)  (P>=1) => c = 1/kappa^2
def build_b(P, V, cval=c):
    pp, pcf = poly('pb', P); nu, ncf = poly('nb', V-1); v = sp.expand(E*nu)
    s = sp.expand(cval*v*sh(v, -1)); B = midB(pp, P); uu = sp.expand((sh(pp, -1)*v - E)/kappa)
    a0 = sp.expand((sh(v, 1)+v)/kappa + (B**2-gamma*B)/kappa**2 + A/kappa)
    aa = {2: sp.Integer(1), 1: pp, 0: a0, -1: uu, -2: s}
    bb = {2: sp.Integer(0), 1: kappa, 0: B, -1: v, -2: sp.Integer(0)}
    return aa, bb, pcf, ncf
for (P, V) in [(1, 1), (1, 2), (2, 2), (2, 3)]:
    aa, bb, pcf, ncf = build_b(P, V)
    lc2 = sp.Poly(Qm(aa, bb, -2), E).LC()
    az(lc2 - (-sp.Rational(1, 1)*P/kappa * ncf[V-1]**2 * pcf[P] * (c*kappa**2 - 1)),
       f"lc(Q_-2(b)) = -(P/kappa) lc(v)^2 lc(p)(c kappa^2-1)  [P={P},V={V}] => c=1/kappa^2")
# at c = 1/kappa^2, lc(Q_-1(b)) = (P/2) lc(p) != 0 for P>=1 => P>=1 empty
for (P, V) in [(1, 1), (1, 2), (2, 1), (2, 2), (3, 2)]:
    aa, bb, pcf, ncf = build_b(P, V, cval=1/kappa**2)
    lc1 = sp.Poly(Qm(aa, bb, -1), E).LC()
    az(lc1 - sp.Rational(1, 2)*P*pcf[P], f"c=1/kappa^2: lc(Q_-1(b)) = (P/2) lc(p) != 0 [P={P},V={V}] => P>=1 empty")
# P=0: Q_-2(b) = [(E-1)v - E v^[-1]]/kappa = 0 => v/E is 1-periodic => v linear
p0 = sp.symbols('p0'); vP0, _ = poly('vP0', 3)
uP0 = sp.expand((p0*vP0 - E)/kappa)  # p const, w=0 => B const => s(B^[-2]-B)=0
Q2_P0 = sp.expand(sh(vP0, -1)*uP0 - sh(uP0, -1)*vP0)   # s-term drops (B const)
az(Q2_P0 - sp.expand(((E-1)*vP0 - E*sh(vP0, -1))/kappa), "P=0: Q_-2(b) = [(E-1)v - E v^[-1]]/kappa => v linear")
# P=0, v = v1 E (linear), s = c v v^[-1]: the finite system Q_-1(b)=0 forces c = 1/kappa^2 (v1!=0).
v1b, p0b = sp.symbols('v1b p0b')
vL = v1b*E; sL = sp.expand(c*vL*sh(vL, -1)); uLb = sp.expand((p0b*vL - E)/kappa)
Bb = midB(p0b, 0)
a0b = sp.expand((sh(vL, 1)+vL)/kappa + (Bb**2-gamma*Bb)/kappa**2 + A/kappa)
ab = {2: sp.Integer(1), 1: p0b, 0: a0b, -1: uLb, -2: sL}
bb = {2: sp.Integer(0), 1: kappa, 0: Bb, -1: vL, -2: sp.Integer(0)}
eqsb = clearden(sp.Poly(Qm(ab, bb, -1), E).all_coeffs() + sp.Poly(Qm(ab, bb, -2), E).all_coeffs())
solb = sp.solve(eqsb, [c, p0b, v1b, gamma, A], dict=True)
# solutions must have c=1/kappa^2 or v1b=0 (v=0, excluded by s!=0)
ok = all((sd.get(c, None) == 1/kappa**2) or (sd.get(v1b, None) == 0) or (v1b in sd and sd[v1b] == 0)
         or (sp.simplify(sd.get(c, c) - 1/kappa**2) == 0) for sd in solb)
if not (solb and ok):
    raise AssertionError(f"branch b P=0 did not force c=1/kappa^2: {solb}")
print("PASS branch b P=0 finite system forces c = 1/kappa^2 (or v=0, excluded by s!=0)")
print("  BRANCH b = tame family: P=0, v linear, c=1/kappa^2 (verified genuine below).\n")

# ==========================================================================
# BRANCH * (resistant, s!=0, w!=0): EMPTY  [proved in quantum-mirror.md];
# regression confirmation here at small degree.
# ==========================================================================
print("--- Branch * (s!=0,w!=0): EMPTY [see quantum-mirror.md]; regression ---")
def resistant_empty(dp, dv, ds):
    pp, pcf = poly('pr', dp); nu, ncf = poly('nr', dv-1); v = sp.expand(E*nu)
    sig, scf = poly('gr', ds-2); s = sp.expand(E*(E-1)*sig); wv = mu*s
    B = midB(pp, dp); uu = sp.expand((mu*(sh(s, 1)+s) + sh(pp, -1)*v - E)/kappa)
    a0 = sp.expand((sh(v, 1)+v)/kappa + (B**2-gamma*B)/kappa**2 + A/kappa)
    aa = {2: sp.Integer(1), 1: pp, 0: a0, -1: uu, -2: s}
    bb = {2: sp.Integer(0), 1: kappa, 0: B, -1: v, -2: wv}
    eqs = clearden(sum([sp.Poly(Qm(aa, bb, m), E).all_coeffs() for m in (-1, -2, -3)], []))
    unk = pcf + ncf + scf + [mu, gamma, A]
    for sci in scf:
        y = sp.Symbol('yr')
        Gg = sp.groebner(eqs + [kappa*mu*sci*y - 1], *(unk+[y]), order='grevlex')
        if list(Gg.exprs) != [sp.Integer(1)]:
            return False
    return True
for (dp, dv, ds) in [(1, 1, 2), (1, 2, 2), (2, 2, 3)]:
    if not resistant_empty(dp, dv, ds):
        raise AssertionError(f"resistant NONEMPTY at ({dp},{dv},{ds})")
    print(f"PASS resistant EMPTY at (dp={dp},dv={dv},ds={ds})")
print()

# ==========================================================================
# TAME FAMILY (normal form): U=x+c0+c1 d, X=U^2 - d/kappa - A, D=lam X + kappa U + beta
# ==========================================================================
print("--- Tame family: [D,X]=1, band 2, w=0, s=(1/kappa^2) v v^[-1], memberships ---")
c0, c1, beta = sp.symbols('c0 c1 beta')
def mul(P, Q):
    R = {}
    for aexp, f in P.items():
        for bexp, g in Q.items():
            R[aexp+bexp] = sp.expand(R.get(aexp+bexp, 0) + sh(f, bexp)*g)
    return {k: v for k, v in R.items() if v != 0}
def add(*Ps):
    R = {}
    for P in Ps:
        for k, v in P.items():
            R[k] = sp.expand(R.get(k, 0) + v)
    return {k: v for k, v in R.items() if v != 0}
def scal(sg, P): return {k: sp.expand(sg*v) for k, v in P.items()}
def commut(P, Q): return add(mul(P, Q), scal(-1, mul(Q, P)))
U = {1: sp.Integer(1), 0: c0, -1: c1*E}; d = {-1: E}
X = add(mul(U, U), scal(-1/kappa, d), {0: -A})
D = add(scal(lam, X), scal(kappa, U), {0: beta})
DX = commut(D, X)
az((DX.get(0, 0) - 1), "tame family: [D,X] = 1 (ladder 0)")
for k in DX:
    if k != 0:
        az(DX[k], f"tame family: [D,X] ladder {k} = 0")
D0 = add(scal(kappa, U), {0: beta})   # reduced gauge b2=0
az(X.get(2, 0) - 1, "tame family: a_2 = 1")
az(D0.get(2, 0), "tame family (reduced): b_2 = 0")
az(D0.get(1, 0) - kappa, "tame family (reduced): b_1 = kappa")
az(D0.get(-2, 0), "tame family (reduced): w = b_-2 = 0")
sT = X.get(-2, 0); vT = D0.get(-1, 0)
az(sT - c1**2*E*(E-1), "tame family: s = a_-2 = c1^2 E(E-1)")
az(sT - (1/kappa**2)*vT*sh(vT, -1), "tame family: s = (1/kappa^2) v v^[-1]  (c = 1/kappa^2)")
az(sT.subs(E, 0), "tame family membership: s(0)=0"); az(sT.subs(E, 1), "tame family membership: s(1)=0")
az(vT.subs(E, 0), "tame family membership: v(0)=0")
az(X.get(-1, 0).subs(E, 0), "tame family membership: a_-1(0)=0")
print("  Tame family is a genuine band-2 pair; U=x+c0+c1 d is an affine symplectic")
print("  generator, so (X,D) generate A_1 (tame automorphism image).\n")

print("ALL QUANTUM COMPLETION CHECKS PASSED")
