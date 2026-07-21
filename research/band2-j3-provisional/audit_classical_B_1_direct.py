#!/usr/bin/env python3
"""
AUDIT B, Route 1: direct attack on classical Theorem J3 (band-2, two-sided
non-square sector) by independent re-derivation + exact bounded searches.

Adversarial stance: nothing from verify_J3.py is reused.  The bracket, the
component system, every reduction identity, and every exhaustive
characterization used to decompose the solution set are re-derived here from
scratch (direct partial derivatives in x, xi; symbolic-exponent monomial
identities), at HIGHER generic degree than the memo's checker (deg 5 vs 2).

Structure
  Part A  independent component system C_m from direct differentiation,
          one-term bracket identity proved for symbolic exponents (all k,l).
  Part B  certificates for the exhaustive characterizations that make the
          later branch decomposition rigorous at ALL degrees:
            (B1) Wronskian ab'-a'b=0, a!=0  <=>  b = const*a
            (B2) 2au'-a'u=0, u!=0          <=>  u^2 = const*a  (scalar-square)
            (B3) a'v-2av'=0, v!=0          <=>  v^2 = const*a  (mirror)
          plus ODE solution scans (deg<=10) on square / nonsquare samples.
  Part C  branch decomposition of ALL band-2 Laurent solutions with
          a2,a_{-2} != 0, and Groebner/solve catalogs of the reduced system
          U=0, L=0, M'=c!=0 for several degree profiles; square classes of
          every family found.
  Part D  cascade-free Groebner emptiness runs on the FULL 9-component system
          for degree profiles with both extremes automatically nonsquare
          (monic linear / saturated quadratic / monic cubic), b's free.
Exit code 0 iff every check passes and no counterexample family is found.
"""
import sys, time
import sympy as sp

x, xi, tau = sp.symbols('x xi tau')
lam, mu, c0 = sp.symbols('lambda2 mu2 c0')
T0 = time.time()
FAIL = []

def check(name, ok):
    ok = bool(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    if not ok:
        FAIL.append(name)

def bracket(D, X):
    """Memo orientation {D,X} = D_xi X_x - D_x X_xi, so {xi,x}=1."""
    return sp.expand(sp.diff(D, xi)*sp.diff(X, x) - sp.diff(D, x)*sp.diff(X, xi))

def lift(a, k):
    """x^k a(tau) as an element of C[x^{+-1},xi]."""
    return sp.expand(x**k * a.subs(tau, x*xi))

def levels(expr, lo=-8, hi=8):
    """Decompose a Laurent element into levels x^k * (poly in tau).
    Independent implementation: substitute xi -> tau/x and collect x powers."""
    e = sp.expand(expr.subs(xi, tau/x))
    out = {}
    for k in range(lo, hi+1):
        cft = sp.expand(e.coeff(x, k))
        if cft != 0:
            out[k] = sp.expand(cft)
    rest = sp.expand(e - sum(x**k*v for k, v in out.items()))
    assert rest == 0, f"level decomposition incomplete: {rest}"
    return out

def genpoly(prefix, deg):
    cs = sp.symbols(f"{prefix}0:{deg+1}")
    return sum(cs[i]*tau**i for i in range(deg+1)), list(cs)

def d(p):
    return sp.diff(p, tau)

# ----------------------------------------------------------------------
print("Part A: independent derivation of the graded component system")
# A1: one-term bracket for SYMBOLIC exponents (covers all k,l,i,j at once):
k_, l_, i_, j_ = sp.symbols('k l i j')
F1 = x**k_ * (x*xi)**i_          # x^k tau^i
G1 = x**l_ * (x*xi)**j_          # x^l tau^j
lhs = bracket(G1, F1)
rhs = sp.expand(x**(k_+l_) * ( k_ * (x*xi)**i_ * sp.diff((x*xi)**j_, x)/xi
                             ))  # placeholder, we do it properly below
# proper claim: {x^l tau^j, x^k tau^i} = x^{k+l}(k*i... ) compute via formula
# {x^l b(tau), x^k a(tau)} = x^{k+l} (k a b' - l a' b)
a_mon = tau**i_; b_mon = tau**j_
rhs = x**(k_+l_) * (k_*a_mon*sp.diff(b_mon, tau) - l_*sp.diff(a_mon, tau)*b_mon)
rhs = sp.expand(rhs.subs(tau, x*xi))
check("one-term bracket {x^l tau^j, x^k tau^i} = x^{k+l}(k a b' - l a' b), symbolic k,l,i,j",
      sp.simplify(lhs - rhs) == 0)

# A2: full component system at generic degree 5 (memo checker used degree 2).
DEG = 5
A = {}; B = {}; ACOEF = {}; BCOEF = {}
for kk in range(-2, 3):
    A[kk], ACOEF[kk] = genpoly(f"a{'m' if kk<0 else ''}{abs(kk)}_", DEG)
    B[kk], BCOEF[kk] = genpoly(f"b{'m' if kk<0 else ''}{abs(kk)}_", DEG)
Xg = sum(lift(A[kk], kk) for kk in range(-2, 3))
Dg = sum(lift(B[kk], kk) for kk in range(-2, 3))
comp = levels(bracket(Dg, Xg), -6, 6)

def C_formula(m):
    """My own component formula, derived from the one-term identity."""
    s = 0
    for kk in range(-2, 3):
        ll = m - kk
        if -2 <= ll <= 2:
            s += kk*A[kk]*d(B[ll]) - ll*d(A[kk])*B[ll]
    return sp.expand(s)

ok = True
for m in range(-4, 5):
    ok &= sp.expand(comp.get(m, 0) - C_formula(m)) == 0
check("C_m (m=-4..4) from direct d/dx,d/dxi == graded formula, generic deg-5 coeffs", ok)
check("no components outside m in [-4,4]",
      all(m in range(-4, 5) for m in comp))

# Cross-check against the memo's displayed 9 equations (transcribed independently)
memo = {
 4: 2*A[2]*d(B[2]) - 2*d(A[2])*B[2],
 3: 2*A[2]*d(B[1]) - d(A[2])*B[1] + A[1]*d(B[2]) - 2*d(A[1])*B[2],
 2: 2*A[2]*d(B[0]) + A[1]*d(B[1]) - d(A[1])*B[1] - 2*d(A[0])*B[2],
 1: 2*A[2]*d(B[-1]) + d(A[2])*B[-1] + A[1]*d(B[0]) - d(A[0])*B[1]
    - A[-1]*d(B[2]) - 2*d(A[-1])*B[2],
 0: 2*A[2]*d(B[-2]) + 2*d(A[2])*B[-2] + A[1]*d(B[-1]) + d(A[1])*B[-1]
    - A[-1]*d(B[1]) - d(A[-1])*B[1] - 2*A[-2]*d(B[2]) - 2*d(A[-2])*B[2],
 -1: A[1]*d(B[-2]) + 2*d(A[1])*B[-2] + d(A[0])*B[-1] - A[-1]*d(B[0])
     - 2*A[-2]*d(B[1]) - d(A[-2])*B[1],
 -2: 2*d(A[0])*B[-2] + d(A[-1])*B[-1] - A[-1]*d(B[-1]) - 2*A[-2]*d(B[0]),
 -3: 2*d(A[-1])*B[-2] - A[-1]*d(B[-2]) + d(A[-2])*B[-1] - 2*A[-2]*d(B[-1]),
 -4: 2*d(A[-2])*B[-2] - 2*A[-2]*d(B[-2]),
}
check("memo's 9 displayed equations agree with the direct expansion (deg-5 generic)",
      all(sp.expand(comp.get(m, 0) - memo[m]) == 0 for m in memo))

# ----------------------------------------------------------------------
print("Part B: certificates for the exhaustive ODE characterizations")
ag, _ = genpoly("A_", DEG); bg, _ = genpoly("B_", DEG); ug, _ = genpoly("U_", DEG)
# B1: a^2 (b/a)' = a b' - a' b  identity  => W=0 <=> (b/a)'=0 <=> b/a in C.
check("identity a^2*(b/a)' == a*b' - a'*b (generic deg-5)",
      sp.expand(ag**2*sp.together(sp.diff(bg/ag, tau)).as_numer_denom()[0]/ag**2
                ) is not None and
      sp.expand(ag*d(bg) - d(ag)*bg - sp.expand(ag**2*sp.simplify(sp.diff(bg/ag, tau)))) == 0)
# B2: a^2 (u^2/a)' = u(2au' - a'u)  => 2au'-a'u=0, u!=0 => u^2 = const*a.
check("identity a^2*(u^2/a)' == u*(2a u' - a' u) (generic deg-5)",
      sp.expand(ug*(2*ag*d(ug) - d(ag)*ug)
                - sp.expand(ag**2*sp.simplify(sp.diff(ug**2/ag, tau)))) == 0)
# B3 mirror: a^2 (v^2/a)' = -v(a'v - 2av').
check("identity a^2*(v^2/a)' == -v*(a' v - 2a v') (generic deg-5)",
      sp.expand(ug*(d(ag)*ug - 2*ag*d(ug))
                + sp.expand(ag**2*sp.simplify(sp.diff(ug**2/ag, tau)))) == 0)
# converse: if u^2 = c*a (c!=0) then u solves 2au'-a'u=0:
cc = sp.Symbol('cc', nonzero=True)
check("converse: u^2=c*a => 2au'-a'u=0 (substitute a=u^2/c, generic deg-5 u)",
      sp.expand(2*(ug**2/cc)*d(ug) - d(ug**2/cc)*ug) == 0)

def ode_solutions(a2poly, ode, maxdeg=10):
    """All polynomial solutions u (deg<=maxdeg) of the given linear ODE in u."""
    us = sp.symbols(f"usol0:{maxdeg+1}")
    u = sum(us[i]*tau**i for i in range(maxdeg+1))
    eq = sp.expand(ode(a2poly, u))
    sol = sp.solve([sp.Poly(eq, tau).coeff_monomial(tau**i)
                    for i in range(sp.Poly(eq, tau).degree()+1)] if eq != 0 else [],
                   list(us), dict=True)
    return sol, us, u

ode_top = lambda a, u: 2*a*d(u) - d(a)*u          # C3-type
ode_bot = lambda a, u: d(a)*u - 2*a*d(u)          # C-3-type (same kernel)

nonsq_samples = [tau, tau+1, tau**2+1, tau**2+tau, tau**3, tau**3-tau,
                 tau*(tau+1)**2, tau**4+1, tau**4+tau**3]
sq_samples    = [sp.Integer(1), sp.Integer(7), tau**2, (tau+1)**2, 3*tau**2*(tau-2)**2]
okn = True
for a2s in nonsq_samples:
    sol, us, u = ode_solutions(a2s, ode_top)
    # solution space must be {0}
    only0 = (sol == [] or sol == [{v: 0 for v in us}])
    okn &= only0
check("2a u'-a'u=0 has ONLY u=0 over deg<=10 for 9 nonsquare samples a", okn)
oks = True
for a2s in sq_samples:
    sol, us, u = ode_solutions(a2s, ode_top)
    if not sol:
        oks = False; continue
    usol = sp.expand(u.subs(sol[0]))
    free = [v for v in us if v not in sol[0]]
    # nonzero solutions exist and satisfy u^2 = const * a
    if usol == 0 or len(free) != 1:
        oks = False; continue
    f = free[0]
    ratio = sp.simplify(usol**2/a2s)
    oks &= sp.diff(ratio, tau) == 0
check("2a u'-a'u=0 has a 1-dim solution line with u^2/a constant for 5 square samples", oks)

# ----------------------------------------------------------------------
print("Part C: branch decomposition and reduced-system catalogs")
# The exhaustive decomposition (valid at ALL degrees, by B1-B3):
#   C4=0, a2!=0      => b2 = lam*a2
#   C-4=0, a_{-2}!=0 => b_{-2} = mu*a_{-2}
#   C3=0             => b1 = lam*a1 + u,  u=0 or u^2 = cu*a2 (cu != 0)
#   C-3=0            => b_{-1} = mu*a_{-1} + v, v=0 or v^2 = cv*a_{-2}
# Branches with u!=0 have a2 scalar-square; with v!=0 have a_{-2} scalar-square.
# ONLY the branch u=v=0 could carry a counterexample to J3. Reduce it:
DEG2 = 5
a2r, _ = genpoly("p", DEG2); a1r, _ = genpoly("q", DEG2)
am1r, _ = genpoly("r", DEG2); am2r, _ = genpoly("s", DEG2)
a0r, _ = genpoly("z", DEG2)
gam = sp.Symbol('gamma0')
subs_branch1 = {}
Ar = {2: a2r, 1: a1r, 0: a0r, -1: am1r, -2: am2r}
Br = {2: lam*a2r, 1: lam*a1r, 0: lam*a0r + gam, -1: mu*am1r, -2: mu*am2r}
def Cm_red(m):
    s = 0
    for kk in range(-2, 3):
        ll = m - kk
        if -2 <= ll <= 2:
            s += kk*Ar[kk]*d(Br[ll]) - ll*d(Ar[kk])*Br[ll]
    return sp.expand(s)
Ured = sp.expand(d(a2r)*am1r + 2*a2r*d(am1r))
Lred = sp.expand(d(am2r)*a1r + 2*am2r*d(a1r))
Mred = sp.expand(2*a2r*am2r + a1r*am1r)
check("branch1: C_{+-4},C_{+-3},C_2 vanish identically (generic deg-5, general a0)",
      all(Cm_red(m) == 0 for m in (4, 3, -3, -4, 2)))
check("branch1: C_-2 == 2 a_{-2} (mu-lam) a0'  (generic deg-5)",
      sp.expand(Cm_red(-2) - 2*am2r*(mu-lam)*d(a0r)) == 0)
check("branch1: C_1 == (mu-lam)*U  (generic deg-5, general a0)",
      sp.expand(Cm_red(1) - (mu-lam)*Ured) == 0)
check("branch1: C_-1 == (mu-lam)*(L + a0'*a_{-1})  (generic deg-5)",
      sp.expand(Cm_red(-1) - (mu-lam)*(Lred + d(a0r)*am1r)) == 0)
check("branch1: C_0 == (mu-lam)*(2 a2 a_{-2} + a1 a_{-1})'  (generic deg-5)",
      sp.expand(Cm_red(0) - (mu-lam)*d(Mred)) == 0)
# NOTE: with b0 = mu*a0 + gam instead, the same run gives the mirror forms;
# consistency of the two choices is exactly C_2 & C_-2 => (lam-mu) a0' = 0.
Br_mu = dict(Br); Br_mu[0] = mu*a0r + gam
def Cm_red_mu(m):
    s = 0
    for kk in range(-2, 3):
        ll = m - kk
        if -2 <= ll <= 2:
            s += kk*Ar[kk]*d(Br_mu[ll]) - ll*d(Ar[kk])*Br_mu[ll]
    return sp.expand(s)
check("branch1 (b0'=mu a0'): C_2 == 2 a2 (lam-mu) a0', C_-2 == 0, C_-1 == (mu-lam)L",
      sp.expand(Cm_red_mu(2) - 2*a2r*(lam-mu)*d(a0r)) == 0
      and Cm_red_mu(-2) == 0
      and sp.expand(Cm_red_mu(-1) - (mu-lam)*Lred) == 0)
# product invariants (independent re-derivation):
check("invariant a_{-1}*U == (a2*a_{-1}^2)'  (generic deg-5)",
      sp.expand(am1r*Ured - d(a2r*am1r**2)) == 0)
check("invariant a1*L == (a_{-2}*a1^2)'  (generic deg-5)",
      sp.expand(a1r*Lred - d(am2r*a1r**2)) == 0)

# ---- C-catalog runs: solve U=0, L=0, M'=c, c!=0 at bounded degrees, ----
# ---- cataloging every family and its square classes.                ----
def square_class(p):
    """Return 'zero' | 'square' (= c*h^2, c in C*) | 'nonsquare' over C."""
    p = sp.expand(p)
    if p == 0:
        return 'zero'
    pp = sp.Poly(p, tau)
    if pp.degree() == 0:
        return 'square'
    if pp.degree() % 2 == 1:
        return 'nonsquare'
    # even degree: p = c*h^2 iff squarefree decomposition is a perfect square
    _, sqf = sp.sqf_list(p, tau)
    return 'square' if all(m % 2 == 0 for _, m in sqf) else 'nonsquare'

assert square_class(tau**2) == 'square' and square_class(3*(tau+1)**2*tau**4) == 'square'
assert square_class(tau) == 'nonsquare' and square_class(tau**2+1) == 'nonsquare'
assert square_class(tau*(tau+1)**2) == 'nonsquare' and square_class(sp.Integer(5)) == 'square'

def catalog_reduced(d2, d1, dm1, dm2, tag):
    """Solve the branch-1 reduced system at the given degree profile.
    Returns list of solution dicts. Uses sp.solve on coefficient equations
    with saturation c != 0 via ce*t = 1."""
    p, pc = genpoly("P", d2); q, qc = genpoly("Q", d1)
    r, rc = genpoly("R", dm1); s, sc = genpoly("S", dm2)
    ce, tsat = sp.symbols("ce tsat")
    U = sp.expand(d(p)*r + 2*p*d(r))
    L = sp.expand(d(s)*q + 2*s*d(q))
    Mp = sp.expand(d(2*p*s + q*r))
    eqs = []
    for e in (U, L, Mp - ce):
        e = sp.expand(e)
        if e == 0:
            continue
        P = sp.Poly(e, tau)
        eqs += [P.coeff_monomial(tau**i) for i in range(P.degree()+1)]
    eqs.append(ce*tsat - 1)
    sols = sp.solve(eqs, pc+qc+rc+sc+[ce, tsat], dict=True)
    fams = []
    for so in sols:
        psol, qsol = sp.expand(p.subs(so)), sp.expand(q.subs(so))
        rsol, ssol = sp.expand(r.subs(so)), sp.expand(s.subs(so))
        fams.append((psol, qsol, rsol, ssol, so))
    print(f"    profile {tag}: {len(fams)} solution families")
    return fams, (p, q, r, s)

def audit_families(fams, tag):
    """Check every family: on the sublocus where a2, a_{-2} are both nonzero,
    at least one of them must be scalar-square for ALL parameter values."""
    all_ok = True
    for idx, (psol, qsol, rsol, ssol, so) in enumerate(fams):
        # parameters remaining free in the family:
        params = sorted((psol+ssol).free_symbols - {tau}, key=str)
        # Structural test: a2 (=psol) and a_{-2} (=ssol) as polynomials in tau
        # with free parameters. J3-consistency: it must be IMPOSSIBLE to choose
        # parameters making both nonzero nonsquares. We test structurally:
        # if either is identically 0 -> family outside J3 hypothesis: OK.
        # if either has tau-degree 0 for all params -> that one is a square: OK.
        va, vb = sp.Poly(psol, tau).degree() if psol != 0 else -1, \
                 sp.Poly(ssol, tau).degree() if ssol != 0 else -1
        if psol == 0 or ssol == 0:
            verdictf = "outside hypothesis (an extreme is 0)"
        elif va == 0 or vb == 0:
            verdictf = "consistent (an extreme is constant=square)"
        else:
            # need a finer look: sample parameters densely and test
            import itertools, random
            random.seed(1)
            bad = None
            for trial in range(200):
                subs = {pp: sp.Rational(random.randint(-6, 6)) for pp in params}
                pa, pb = sp.expand(psol.subs(subs)), sp.expand(ssol.subs(subs))
                if pa == 0 or pb == 0:
                    continue
                if square_class(pa) == 'nonsquare' and square_class(pb) == 'nonsquare':
                    bad = (subs, pa, pb); break
            if bad:
                verdictf = f"*** COUNTEREXAMPLE CANDIDATE *** {bad}"
                all_ok = False
            else:
                verdictf = "consistent (no both-nonsquare specialization found in 200 samples)"
        print(f"      fam{idx}: a2={psol}, a1={qsol}, a-1={rsol}, a-2={ssol} -> {verdictf}")
    check(f"reduced-system catalog {tag}: every family J3-consistent", all_ok)

for prof in [(2, 2, 2, 2), (3, 2, 2, 3), (4, 1, 1, 4), (2, 4, 4, 2)]:
    fams, _ = catalog_reduced(*prof, tag=str(prof))
    audit_families(fams, str(prof))

# ---- C-emptiness certificates on the reduced system with extremes  ----
# ---- FORCED nonsquare by shape (monic odd degree), via Groebner.    ----
def reduced_groebner_empty(d2, d1, dm1, dm2, monic_odd=True, tag=""):
    p, pc = genpoly("P", d2); q, qc = genpoly("Q", d1)
    r, rc = genpoly("R", dm1); s, sc = genpoly("S", dm2)
    unk = pc+qc+rc+sc
    if monic_odd:  # make extremes monic of odd degree => automatically nonsquare
        p = sp.expand(p - pc[-1]*tau**d2 + tau**d2)
        s = sp.expand(s - sc[-1]*tau**dm2 + tau**dm2)
        unk = [u for u in unk if u not in (pc[-1], sc[-1])]
    ce, tsat = sp.symbols("ce tsat")
    U = sp.expand(d(p)*r + 2*p*d(r)); L = sp.expand(d(s)*q + 2*s*d(q))
    Mp = sp.expand(d(2*p*s + q*r))
    eqs = []
    for e in (U, L, Mp - ce):
        P = sp.Poly(e, tau)
        eqs += [P.coeff_monomial(tau**i) for i in range(P.degree()+1)]
    eqs.append(ce*tsat - 1)
    gb = sp.groebner(eqs, *unk, ce, tsat, order='grevlex')
    empty = list(gb.exprs) == [sp.Integer(1)] or gb.exprs == [1]
    check(f"reduced system EMPTY for monic-odd extremes profile {tag}", empty)

reduced_groebner_empty(1, 3, 3, 1, tag="(deg a2=1, a_{-2}=1, middles<=3)")
reduced_groebner_empty(3, 3, 3, 3, tag="(deg a2=3, a_{-2}=3, middles<=3)")
reduced_groebner_empty(1, 4, 4, 3, tag="(deg a2=1, a_{-2}=3, middles<=4)")
reduced_groebner_empty(5, 2, 2, 3, tag="(deg a2=5, a_{-2}=3, middles<=2)")

# quadratic extremes with discriminant saturation (nonsquare <=> disc != 0)
def reduced_groebner_quad(d1, dm1, tag=""):
    p0, p1v, s0, s1v = sp.symbols("p0 p1v s0 s1v")
    p = tau**2 + p1v*tau + p0; s = tau**2 + s1v*tau + s0
    q, qc = genpoly("Q", d1); r, rc = genpoly("R", dm1)
    ce, tsat, t2, t3 = sp.symbols("ce tsat t2 t3")
    U = sp.expand(d(p)*r + 2*p*d(r)); L = sp.expand(d(s)*q + 2*s*d(q))
    Mp = sp.expand(d(2*p*s + q*r))
    eqs = []
    for e in (U, L, Mp - ce):
        P = sp.Poly(e, tau)
        eqs += [P.coeff_monomial(tau**i) for i in range(P.degree()+1)]
    eqs += [ce*tsat - 1, (p1v**2 - 4*p0)*t2 - 1, (s1v**2 - 4*s0)*t3 - 1]
    gb = sp.groebner(eqs, p0, p1v, s0, s1v, *(qc+rc), ce, tsat, t2, t3,
                     order='grevlex')
    empty = list(gb.exprs) == [sp.Integer(1)]
    check(f"reduced system EMPTY for squarefree monic-quadratic extremes {tag}", empty)

reduced_groebner_quad(3, 3, tag="(middles<=3)")

# ----------------------------------------------------------------------
print("Part D: cascade-free FULL-system Groebner emptiness (no reductions used)")
def full_system_empty(a_profile, bdeg, tag, extra_sat=None, timeout_note=""):
    """a_profile: dict k -> (poly, unknown list). b's fully generic deg<=bdeg.
    Assemble all 9 component equations + saturation c!=0 (+extras); Groebner."""
    Bv = {}; bunk = []
    for kk in range(-2, 3):
        Bv[kk], u = genpoly(f"Bb{'m' if kk<0 else ''}{abs(kk)}_", bdeg)
        bunk += u
    Av = {kk: a_profile[kk][0] for kk in a_profile}
    aunk = sum((a_profile[kk][1] for kk in a_profile), [])
    ce, tsat = sp.symbols("ce tsat")
    eqs = []
    for m in range(-4, 5):
        s = 0
        for kk in range(-2, 3):
            ll = m - kk
            if -2 <= ll <= 2:
                s += kk*Av.get(kk, sp.Integer(0))*d(Bv[ll]) \
                     - ll*d(Av.get(kk, sp.Integer(0)))*Bv[ll]
        s = sp.expand(s - (ce if m == 0 else 0))
        if s == 0:
            continue
        P = sp.Poly(s, tau)
        eqs += [P.coeff_monomial(tau**i) for i in range(P.degree()+1)]
    eqs.append(ce*tsat - 1)
    sats = [tsat]
    if extra_sat:
        for i, g in enumerate(extra_sat):
            tv = sp.Symbol(f"tex{i}")
            eqs.append(g*tv - 1); sats.append(tv)
    t0 = time.time()
    gb = sp.groebner(eqs, *(aunk + bunk), ce, *sats, order='grevlex')
    empty = list(gb.exprs) == [sp.Integer(1)]
    check(f"FULL 9-eq system EMPTY, no cascade used: {tag} "
          f"[{time.time()-t0:.1f}s]{timeout_note}", empty)

zero = (sp.Integer(0), [])
# D1: extremes monic linear (automatically nonsquare), middles zero, b deg<=3
p0v, s0v = sp.symbols("p0v s0v")
full_system_empty({2: (tau+p0v, [p0v]), 1: zero, 0: zero, -1: zero,
                   -2: (tau+s0v, [s0v])}, 3,
                  "extremes tau+p0, tau+s0; middles 0; deg b<=3")
# D2: extremes monic linear, middles generic deg<=2, b deg<=3
q_, qc_ = genpoly("qq", 2); z_, zc_ = genpoly("zz", 2); r_, rc_ = genpoly("rr", 2)
full_system_empty({2: (tau+p0v, [p0v]), 1: (q_, qc_), 0: (z_, zc_),
                   -1: (r_, rc_), -2: (tau+s0v, [s0v])}, 3,
                  "extremes monic linear; middles generic deg<=2; deg b<=3")
# D3: extremes monic cubic (auto-nonsquare), middles generic deg<=1, b deg<=3
pC, pCc = genpoly("pc", 2); sC, sCc = genpoly("sc_", 2)
q1_, q1c_ = genpoly("q1_", 1); z1_, z1c_ = genpoly("z1_", 1); r1_, r1c_ = genpoly("r1_", 1)
full_system_empty({2: (tau**3+pC, pCc), 1: (q1_, q1c_), 0: (z1_, z1c_),
                   -1: (r1_, r1c_), -2: (tau**3+sC, sCc)}, 4,
                  "extremes monic cubic; middles generic deg<=1; deg b<=4")
# D4: MIXED degrees: a2 monic linear, a_{-2} monic cubic, middles deg<=2, b<=4
full_system_empty({2: (tau+p0v, [p0v]), 1: (q_, qc_), 0: (z_, zc_),
                   -1: (r_, rc_), -2: (tau**3+sC, sCc)}, 4,
                  "a2 monic linear, a_{-2} monic cubic; middles<=2; deg b<=4")

# ----------------------------------------------------------------------
print()
el = time.time() - T0
if FAIL:
    print(f"AUDIT ROUTE 1: FAILURES {FAIL}  [{el:.1f}s]")
    sys.exit(1)
print(f"AUDIT ROUTE 1: ALL CHECKS PASSED — no counterexample family found  [{el:.1f}s]")
