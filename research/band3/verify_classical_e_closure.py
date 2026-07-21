#!/usr/bin/env python3
"""Exact verifier for the classical band-3  e != 0  MIXED-SECTOR closure memo
(classical-e-nonzero-closure.md).

Sector: the cubic sector a_3 = c h^3 (gauge b_3 = 0) with a NONZERO wall constant
e (gauged b_2 = e h^2 != 0), constant h.  Normalise h = 1, c = 1 (a genuine
diagonal scaling, classical-hard-branches.md/classical-cube-closure.md).  Wave B
(classical-hard-branches.md, commit ebfc64d) left this sector open: it proved
"tame => e = 0" (counterexample-or-nothing), reduced the positive cascade, and
pinned the W3 obstruction to the single residue (4e/3) a_-3 a_2', concluding that
only ONE clean first integral (the moment M = tau) was available and that the
band-3 analogue of band-2's SECOND integral I_2 was the precise missing step.

This script separates exact symbolic verification from bounded exploration:

  EXACT SECTIONS 0-4 check the 13 C_m formula, the e != 0 reduction, Phi' = C_-1,
             and I_2' = C_-2 - (2/3) a_2' Phi, together with the Euler residues.
             The multiplier is unique only up to an additive constant:
             m'=(2/3)a_2', hence m=(2/3)a_2+constant. They also check the two
             conserved quantities.
  EXACT PART OF SECTION 5 checks the trailing 2x2 determinant and the exact
             leading coefficients in one a_2-dominant specialization. The solve
             is conditional on det != 0; this script does not prove denominator
             cancellation, polynomiality, or coefficient membership.
  BOUNDED PART OF SECTION 5 explores generic leading monomials only for
             Q=1..3, P,R=0..3, L=1..3. It sets conserved constants CP=CI=0 and
             does not resolve cancellation/tie loci. It is corroboration, not an
             arbitrary-degree proof and not a sector-emptiness theorem.
  SECTION 6 checks the a_2-constant residue switch-off, one bounded survivor, and
             one exact bounded Groebner emptiness box.

A successful run prints 'ALL CLASSICAL E CLOSURE CHECKS PASSED'. Exact checks and
bounded exploration are labelled separately in the output.
"""
import sympy as sp

t = sp.symbols("tau")
KS = range(-3, 4)
PASS = 0


def D(f):
    return sp.diff(f, t)


def Cm(m, A, B):
    return sp.expand(sum(k*A[k]*D(B[l]) - l*D(A[k])*B[l]
                         for k in KS for l in KS if k + l == m))


def poisson(g, f, x, xi):
    return sp.diff(g, xi)*sp.diff(f, x) - sp.diff(g, x)*sp.diff(f, xi)


def check(cond, label):
    global PASS
    if cond is not True and cond != True:  # noqa: E712
        raise AssertionError("FAIL: " + label)
    PASS += 1
    print("PASS", label)


def zero(expr, label):
    check(sp.simplify(sp.expand(expr)) == 0, label)


# =====================================================================
# SECTION 0.  The 13 C_m match the direct two-variable Poisson bracket.
# =====================================================================
x, xi = sp.symbols("x xi")
Ag = {k: sp.Function(f"a{k}")(t) for k in KS}
Bg = {k: sp.Function(f"b{k}")(t) for k in KS}
f_full = sum(x**k*Ag[k].subs(t, x*xi) for k in KS)
g_full = sum(x**k*Bg[k].subs(t, x*xi) for k in KS)
formula = sp.expand(sum(x**m*Cm(m, Ag, Bg).subs(t, x*xi) for m in range(-6, 7)))
zero(poisson(g_full, f_full, x, xi) - formula,
     "all 13 C_m match the direct 2-variable Poisson bracket (generic coefficients)")

# =====================================================================
# SECTION 1.  The e != 0 reduction (h = 1, c = 1).
#   Free F-data: q=a_2, p=a_1, r=a_0, al=a_-1 and trailing s=a_-2, sig=a_-3.
# =====================================================================
e, k1, beta, gam, delta = sp.symbols("e kappa1 beta gamma delta")
q = sp.Function("q")(t); p = sp.Function("p")(t); r = sp.Function("r")(t)
al = sp.Function("al")(t); s = sp.Function("s")(t); sig = sp.Function("sig")(t)
A = {3: sp.Integer(1), 2: q, 1: p, 0: r, -1: al, -2: s, -3: sig}

# wall + mixed coupling (classical-hard-branches.md 4.1/§5, ebfc64d):
b2 = e
b1 = k1 + sp.Rational(2, 3)*e*q
b0 = (2*e*p + k1*q)/3 - e*q**2/9 + beta
bm1 = sp.Rational(1, 3)*(2*e*r + k1*p - sp.Rational(2, 3)*e*p*q
                         - sp.Rational(1, 3)*k1*q**2 + sp.Rational(4, 27)*e*q**3) + gam
# b_-2 EXPLICIT: the C_1 right-hand side is a total tau-derivative:
bm2 = sp.expand(
    -sp.Rational(1, 3)*(-sp.Rational(4, 9)*e*p*q**2 + sp.Rational(1, 3)*e*p**2
                        + sp.Rational(7, 81)*e*q**4 + sp.Rational(2, 3)*e*q*r - 2*e*al
                        + gam*q + sp.Rational(2, 3)*k1*p*q - sp.Rational(5, 27)*k1*q**3
                        - k1*r) + delta)
# b_-3 from the moment M = tau:
bm3 = sp.expand((t - 2*(q*bm2 - s*e) - (p*bm1 - al*b1))/3)
B = {3: sp.Integer(0), 2: b2, 1: b1, 0: b0, -1: bm1, -2: bm2, -3: bm3}

zero(Cm(3, A, B), "e!=0: b_0 solves C_3")
zero(Cm(2, A, B), "e!=0: b_-1 solves C_2")
zero(Cm(1, A, B), "e!=0: b_-2 is EXPLICIT (C_1 RHS is a total tau-derivative)")
check(not bm2.has(sp.Integral), "e!=0: b_-2 carries no surviving integral (fully local)")
# b_-2 genuinely depends on a_-1 (new vs e=0): the -2e a_-1 term.
check(bm2.has(al), "e!=0: b_-2 depends on a_-1 (the e-coupling; absent when e=0)")
# moment identity C_0 = M':
Mmom = sum(kk*(A[kk]*B[-kk] - A[-kk]*B[kk]) for kk in [1, 2, 3])
zero(Cm(0, A, B) - D(Mmom), "e!=0: C_0 = M', moment M = 3(a3 b-3)+2(a2 b-2 - a-2 b2)+(a1 b-1 - a-1 b1)")
zero(Cm(0, A, B) - 1, "e!=0: C_0 = 1 identically (b_-3 chosen so M = tau)")

print(f"\n--- Section 0-1 (reduction) done ({PASS} checks) ---\n")

# =====================================================================
# SECTION 2.  Phi : Phi' = C_-1, nonlocal Q1 = int a_2.
# =====================================================================
Q1 = sp.Function("Q1")(t)


def Dq(expr):
    return sp.expand(sp.diff(expr, t).subs(sp.Derivative(Q1, t), q))


Phi = (
    # ---- e = 0 part (equals cube-closure Phi_1 at c=1) ----
    2*delta*p - sp.Rational(5, 3)*delta*q**2 - sp.Rational(4, 3)*gam*p*q
    + sp.Rational(14, 27)*gam*q**3 + gam*r + sp.Rational(2, 3)*k1*al*q
    - sp.Rational(5, 9)*k1*p**2*q + sp.Rational(40, 81)*k1*p*q**3 + sp.Rational(2, 3)*k1*p*r
    - sp.Rational(22, 243)*k1*q**5 - sp.Rational(5, 9)*k1*q**2*r - k1*s + t*q - Q1/3
    # ---- e-part ----
    + e*(sp.Rational(2, 3)*al*p - sp.Rational(4, 9)*al*q**2 - sp.Rational(4, 27)*p**3
         + sp.Rational(14, 27)*p**2*q**2 - sp.Rational(70, 243)*p*q**4 - sp.Rational(8, 9)*p*q*r
         + sp.Rational(91, 2187)*q**6 + sp.Rational(28, 81)*q**3*r + sp.Rational(2, 3)*q*s
         + sp.Rational(1, 3)*r**2 - 2*sig))
Phi = sp.expand(Phi)
Cm1 = sp.expand(Cm(-1, A, B))
zero(Dq(Phi) - Cm1, "Phi' = C_-1 exactly (the first integral; nonlocal Q1 = int a_2)")

# Phi|_{e=0} equals cube-closure Phi_1 (commit ebfc64d) exactly:
Phi1_cc = (2*delta*p + gam*r - k1*s - sp.Rational(5, 3)*delta*q**2 - sp.Rational(4, 3)*gam*p*q
           + sp.Rational(2, 3)*k1*al*q + sp.Rational(2, 3)*k1*p*r + t*q - Q1/3
           + sp.Rational(14, 27)*gam*q**3 - sp.Rational(5, 9)*k1*p**2*q - sp.Rational(5, 9)*k1*q**2*r
           + sp.Rational(40, 81)*k1*p*q**3 - sp.Rational(22, 243)*k1*q**5)
zero(Phi.subs(e, 0) - Phi1_cc,
     "Phi|_{e=0} = cube-closure Phi_1 (honest e!=0 generalisation of the W3-replacement)")

# Phi carries a_-2 and a_-3 LINEARLY (trailing coefficients):
check(sp.expand(Phi.coeff(sig, 2)) == 0 and sp.expand(Phi.coeff(s, 2)) == 0,
      "Phi is linear in a_-2 and a_-3 (coeff of a_-3 = -2e, coeff of a_-2 = (2e/3)a_2 - kappa1)")
zero(Phi.coeff(sig, 1) + 2*e, "Phi: coefficient of a_-3 is -2e")
zero(Phi.coeff(s, 1) - (sp.Rational(2, 3)*e*q - k1), "Phi: coefficient of a_-2 is (2e/3)a_2 - kappa1")

print(f"\n--- Section 2 (Phi' = C_-1) done ({PASS} checks) ---\n")

# =====================================================================
# SECTION 3.  THE SECOND INTEGRAL I_2 (the band-3 analogue of band-2's I_2).
#   I_2' = C_-2 - (2/3) a_2' Phi ; nonlocals Q1 = int a_2, Q2 = int a_2^2, P1 = int a_1.
# =====================================================================
Q2 = sp.Function("Q2")(t); P1 = sp.Function("P1")(t)


def Dqq(expr):
    return sp.expand(sp.diff(expr, t).subs({sp.Derivative(Q1, t): q,
                     sp.Derivative(Q2, t): q**2, sp.Derivative(P1, t): p}))


I2 = (
    # ---- e = 0 part ----
    -2*delta*p*q + sp.Rational(10, 27)*delta*q**3 + 2*delta*r + gam*al - sp.Rational(2, 3)*gam*p**2
    + sp.Rational(2, 3)*gam*p*q**2 - sp.Rational(7, 81)*gam*q**4 - sp.Rational(2, 3)*gam*q*r
    + sp.Rational(2, 3)*k1*al*p - sp.Rational(1, 9)*k1*al*q**2 - sp.Rational(5, 27)*k1*p**3
    + sp.Rational(10, 27)*k1*p**2*q**2 - sp.Rational(10, 81)*k1*p*q**4 - sp.Rational(2, 3)*k1*p*q*r
    + sp.Rational(22, 2187)*k1*q**6 + sp.Rational(10, 81)*k1*q**3*r + sp.Rational(1, 3)*k1*r**2
    - k1*sig + t*p - sp.Rational(1, 3)*t*q**2 - sp.Rational(2, 3)*P1 + sp.Rational(2, 9)*Q1*q
    + sp.Rational(1, 9)*Q2
    # ---- e-part ----
    + e*(-sp.Rational(4, 9)*al*p*q + sp.Rational(4, 81)*al*q**3 + sp.Rational(2, 3)*al*r
         + sp.Rational(20, 81)*p**3*q - sp.Rational(56, 243)*p**2*q**3 - sp.Rational(4, 9)*p**2*r
         + sp.Rational(14, 243)*p*q**5 + sp.Rational(4, 9)*p*q**2*r + sp.Rational(2, 3)*p*s
         - sp.Rational(26, 6561)*q**7 - sp.Rational(14, 243)*q**4*r - sp.Rational(2, 9)*q*r**2
         - sp.Rational(2, 3)*q*sig))
I2 = sp.expand(I2)
Cm2 = sp.expand(Cm(-2, A, B))
target = sp.expand(Cm2 - sp.Rational(2, 3)*D(q)*Phi)
zero(Dqq(I2) - target,
     "I_2' = C_-2 - (2/3) a_2' Phi  (THE second integral; nonlocals Q1, Q2=int a_2^2, P1=int a_1)")


# --- the multiplier (2/3)a_2 is FORCED by the W3 residue ---
def euler(L, u):
    return sp.expand(sp.diff(L, u) - D(sp.diff(L, D(u))))


zero(euler(Cm2, sig) + sp.Rational(4, 3)*e*D(q),
     "W3 obstruction: Euler_a-3(C_-2) = -(4e/3) a_2' != 0 (a_-3 not an exact derivative in C_-2)")
zero(euler(target, sig),
     "the multiplier (2/3)a_2 ABSORBS it: Euler_a-3(C_-2 - (2/3)a_2' Phi) = 0")
zero(euler(target, s),
     "and Euler_a-2(C_-2 - (2/3)a_2' Phi) = 0 (the a_-2 residue is absorbed too)")
# the multiplier value m is FORCED: Euler_a-3(C_-2 - (m)' Phi) = -(4e/3)a_2' + 2e*m'
# vanishes iff m'=(2/3)a_2', i.e. m=(2/3)a_2+constant.
# (Phi's a_-3 coefficient is -2e.)
mgen = sp.Function("m")(t)
zero(euler(sp.expand(Cm2 - D(mgen)*Phi), sig) - (-sp.Rational(4, 3)*e*D(q) + 2*e*D(mgen)),
     "for ANY multiplier m: Euler_a-3(C_-2 - m' Phi) = -(4e/3)a_2' + 2e m'; "
     "=0 forces m=(2/3)a_2 + constant")

# I_2 carries a_-2, a_-3 linearly too:
zero(I2.coeff(sig, 1) + (k1 + sp.Rational(2, 3)*e*q),
     "I_2: coefficient of a_-3 is -(kappa1 + (2e/3)a_2) = -b_1")
zero(I2.coeff(s, 1) - sp.Rational(2, 3)*e*p, "I_2: coefficient of a_-2 is (2e/3)a_1")

print(f"\n--- Section 3 (the SECOND integral I_2) done ({PASS} checks) ---\n")

# =====================================================================
# SECTION 4.  On a solution: Phi = const, and J := I_2 + (2/3) Phi_0 a_2 = const.
# =====================================================================
Phi0 = sp.Symbol("Phi0")
# on solution C_-1 = 0 => Phi' = 0 => Phi = Phi0 ; then C_-2 = 0 => I_2' = -(2/3)Phi0 a_2'
zero(sp.expand(Dqq(I2 + sp.Rational(2, 3)*Phi0*q)
               - (target + sp.Rational(2, 3)*Phi0*D(q))),
     "J := I_2 + (2/3)Phi0 a_2 has J' = C_-2 - (2/3)a_2'(Phi - Phi0);"
     " so on Phi=Phi0 & C_-1=C_-2=0, J' = 0 => J = const")

print(f"\n--- Section 4 (two conserved quantities) done ({PASS} checks) ---\n")

# =====================================================================
# SECTION 5. Conditional trailing solve and bounded degree exploration.
#   Exact: Phi and I_2 are linear in (a_-2,a_-3), with the determinant below.
#   The rational solve applies only when det != 0; polynomiality/membership are
#   not checked here. The later scan is bounded generic-leading-monomial evidence.
# =====================================================================
CP, CI = sp.symbols("CP CI")
Ps, PS = Phi.coeff(sig, 1), Phi.coeff(s, 1)
P0 = sp.expand(Phi - Ps*sig - PS*s)
Is, IS = I2.coeff(sig, 1), I2.coeff(s, 1)
I0 = sp.expand(I2 - Is*sig - IS*s)
det = sp.expand(Ps*IS - PS*Is)
zero(det - (-sp.Rational(4, 3)*e**2*p + sp.Rational(4, 9)*e**2*q**2 - k1**2),
     "the [a_-2,a_-3] system determinant = -(4/3)e^2 a_1 + (4/9)e^2 a_2^2 - kappa1^2")
sig_num = sp.expand(IS*(CP - P0) - PS*(CI - I0))   # a_-3 = sig_num/det
bm3_num = sp.expand(det*t - 2*(det*q*bm2 - e*(-Is*(CP - P0) + Ps*(CI - I0)))
                    - det*(p*bm1 - al*b1))          # 3 det b_-3


def leaddata(expr, Q, P, R, L):
    """Return degree data after a bounded generic-leading-monomial substitution.

    CP and CI are set to zero. This substitution does not classify coefficient
    cancellation/tie loci and does not test divisibility by the determinant.
    """
    lcq, lcp, lcr, lcl = sp.symbols("lcq lcp lcr lcl")
    sub = {q: lcq*t**Q, p: lcp*t**P, r: lcr*t**R, al: lcl*t**L,
           Q1: lcq*t**(Q + 1)/(Q + 1), Q2: lcq**2*t**(2*Q + 1)/(2*Q + 1),
           P1: lcp*t**(P + 1)/(P + 1), CP: 0, CI: 0}
    ex = sp.expand(expr.subs(sub))
    if ex == 0:
        return (-1, sp.Integer(0))
    pol = sp.Poly(ex, t)
    return (pol.degree(), pol.LC())


# Exact q-dominant specialization certificates: deg a_-3=6Q and deg b_-3=5Q.
# These are regime-specific and conditional on the nonzero determinant.
dd, lcdet = leaddata(det, 2, 1, 1, 1)
ds, lcsig = leaddata(sig_num, 2, 1, 1, 1)
db, lcbm3 = leaddata(bm3_num, 2, 1, 1, 1)
lcq = sp.symbols("lcq")
rep_scope = ("[REPRESENTATIVE (Q,P,R,L)=(2,1,1,1); CP=CI=0; "
             "generic leading monomials; ties/cancellation omitted] ")
check(ds - dd == 12, rep_scope + "deg a_-3 = 12")
check(db - dd == 10, rep_scope + "deg b_-3 = 10 < deg a_-3")
zero(sp.simplify(lcsig/lcdet) + sp.Rational(13, 2187)*lcq**6,
     rep_scope + "leading coeff of a_-3 is -13 lc(a_2)^6/2187 != 0")
zero(sp.simplify(lcbm3/lcdet) + sp.Rational(25, 243)*e*lcq**5,
     rep_scope + "leading coeff of b_-3 is -25 e lc(a_2)^5/243 != 0")

# BOUNDED EXPLORATION ONLY: Q=1..3, P,R=0..3, L=1..3.
# Generic leading monomials, CP=CI=0; cancellation/tie loci are omitted.
gap_fail = []
a3_not6Q = []
for Q in range(1, 4):
    for P in range(0, 4):
        for R in range(0, 4):
            for L in range(1, 4):
                dda, _ = leaddata(det, Q, P, R, L)
                dsa, _ = leaddata(sig_num, Q, P, R, L)
                dba, _ = leaddata(bm3_num, Q, P, R, L)
                deg_a3 = dsa - dda
                deg_b3 = dba - dda
                if deg_a3 >= 0 and deg_b3 >= deg_a3:
                    gap_fail.append((Q, P, R, L, deg_a3, deg_b3))
                # a_-3 top is max(6Q, 2R, ...): check it is >= 6Q and gap
                if deg_a3 < 6*Q:
                    a3_not6Q.append((Q, P, R, L, deg_a3))
check(gap_fail == [],
      "[BOUNDED GENERIC SCAN Q=1..3,P,R=0..3,L=1..3; CP=CI=0; ties omitted] "
      "no sampled degree has deg b_-3 >= deg a_-3")
check(a3_not6Q == [],
      "[BOUNDED GENERIC SCAN; not an arbitrary-degree proof] no sampled degree has "
      "deg a_-3 < 6 deg a_2")

print(f"\n--- Section 5 (exact conditional certificates + bounded generic scan) done ({PASS} checks) ---\n")

# =====================================================================
# SECTION 6.  a_2 = const (Q = 0) residual stratum + Groebner consistency.
# =====================================================================
# When a_2 is constant, b_1=kappa1+(2e/3)a_2 is constant and this particular W3
# residue vanishes. That does not establish that the full stratum is unobstructed.
# We check only the residue switch-off:
zero(euler(Cm2, sig).subs(D(q), 0), "a_2 constant => the particular W3 residue -(4e/3)a_2' switches off")

# One bounded Q=0 survivor exists outside [CUBE]'s three displayed tie equations.
# This does not classify the Q=0 stratum or prove that its global survivor set is finite.
P_, R_, L_ = 0, 1, 1
dda0, _ = leaddata(det, 0, P_, R_, L_)
dsa0, _ = leaddata(sig_num, 0, P_, R_, L_)
dba0, _ = leaddata(bm3_num, 0, P_, R_, L_)
off_cube = not ((3*P_ == 2*R_) or (L_ == 2*P_) or (2*R_ == L_ + P_))
check((dsa0 - dda0) <= (dba0 - dda0) and off_cube,
      "a_2=const stratum: (P,R,L)=(0,1,1) survives the gap and is OFF {3P=2R,L=2P,2R=L+P}"
      " (one bounded survivor; no global classification claimed)")

# Independent bounded corroboration: re-verify one SMALL exact Groebner emptiness
# box for e != 0. This does not follow from, or establish, an unbounded tropical kill.
def _pc(z):
    z = sp.expand(z)
    return [] if z == 0 else sp.Poly(z, t).all_coeffs()


def _gp(n, d, lo=0):
    cs = [sp.Symbol(f"{n}_{i}") for i in range(d + 1)]
    return sum(cs[i]*t**i for i in range(lo, d + 1)), cs[lo:]


def enonzero_box_empty(degs):
    da2, da1, da0, dam1, dam2, dam3 = degs
    a2f, ca2 = _gp("A2", da2)
    a1f, ca1 = _gp("A1", da1)
    a0f, ca0 = _gp("A0", da0)
    am1f, cam1 = _gp("Am1", dam1, 1)
    am2f, cam2 = _gp("Am2", dam2, 2)
    am3f, cam3 = _gp("Am3", dam3, 3)
    Ad = {3: sp.Integer(1), 2: a2f, 1: a1f, 0: a0f, -1: am1f, -2: am2f, -3: am3f}
    ee, kk = sp.Symbol("e_"), sp.Symbol("k1_")
    b1v = kk + sp.Rational(2, 3)*ee*a2f
    b0f, cb0 = _gp("B0", da2 + 2)
    bm1g, cbm1 = _gp("Bm1", dam1 + 3, 1)
    bm2g, cbm2 = _gp("Bm2", dam2 + 4, 2)
    bm3g, cbm3 = _gp("Bm3", dam3 + 5, 3)
    Bd = {3: sp.Integer(0), 2: ee, 1: b1v, 0: b0f, -1: bm1g, -2: bm2g, -3: bm3g}
    cons = []
    for m in range(-6, 7):
        cons += _pc(sp.expand(Cm(m, Ad, Bd) - (1 if m == 0 else 0)))
    cons = [z for z in (sp.expand(c) for c in cons) if z != 0]
    w = sp.Symbol("w_")
    cons = cons + [ee*w - 1]  # Rabinowitsch: e != 0
    gens = ([w, ee, kk] + list(ca2) + list(ca1) + list(ca0) + list(cam1)
            + list(cam2) + list(cam3) + list(cb0) + list(cbm1) + list(cbm2) + list(cbm3))
    Gb = sp.groebner(cons, *gens, order='grevlex')
    return list(Gb.exprs) == [sp.Integer(1)]


check(enonzero_box_empty((1, 0, 0, 1, 2, 3)) is True,
      "[COMPUTED] e!=0 box deg(a2..a-3)=(1,0,0,1,2,3): Groebner = (1), no pair"
      " (bounded computation only)")

print(f"\n--- Section 6 (a_2=const stratum + Groebner corroboration) done ({PASS} checks) ---\n")
print()
print(f"({PASS} checks)")
print("ALL CLASSICAL E CLOSURE CHECKS PASSED")
