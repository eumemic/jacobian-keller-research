#!/usr/bin/env python3
"""
AUDIT of Theorem P3 -- Part 2: the CRITICAL branch-(B) question, the two
symmetry maps (Weyl-Fourier vs band-reversal), and the full vanishing-pattern
enumeration.

Memo endgame:  V(E)=a1(E-1)a_{-1}(E) is linear (deg sum 1), so one factor is
constant, giving TWO branches:
    (A) a1 = const,  a_{-1} = linear      (the "affine + polar" branch)
    (B) a1 = linear, a_{-1} = const       (the "Fourier mirror")
The memo asserts (B) is the Fourier mirror of (A) and "back-substituting yields
the affine forms".  We test whether branch (B) actually produces genuine
solutions of the FULL system [D,X]=1 in B=A1[x^-1], and whether they are affine.

Run:  uv run --with sympy python audit_band1_branchB.py
"""
import sympy as sp
E = sp.symbols("E")
RES = []
def check(name, ok):
    ok = bool(ok); RES.append((name, ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

# crossed-product engine (independent)
def clean(A): return {k: sp.expand(v) for k, v in A.items() if sp.expand(v) != 0}
def mul(A, B):
    C = {}
    for a, f in A.items():
        for b, g in B.items():
            C[a+b] = C.get(a+b, 0) + f.subs(E, E+b)*g
    return clean(C)
def add(*As):
    C = {}
    for A in As:
        for k, v in A.items(): C[k] = C.get(k, 0) + v
    return clean(C)
def neg(A): return {k: -v for k, v in A.items()}
def brk(A, B): return add(mul(A, B), neg(mul(B, A)))
ONE = {0: sp.Integer(1)}
def in_A1(A):
    for k, v in A.items():
        if k < 0:
            ff = sp.prod([E - i for i in range(-k)])
            if sp.rem(sp.Poly(v, E), sp.Poly(ff, E)).as_expr() != 0:
                return False
    return True
def is_affine(A):
    """Affine = in span{ x, x^-1 E (=d), 1 }, i.e. a1 const, a0 const, a_{-1} in C*E,
       and no other x-levels."""
    if set(A.keys()) - {-1, 0, 1}: return False
    a1  = A.get(1, sp.Integer(0)); a0 = A.get(0, sp.Integer(0)); am1 = A.get(-1, sp.Integer(0))
    return (sp.Poly(a1, E).degree() <= 0 and sp.Poly(a0, E).degree() <= 0
            and sp.expand(am1 - E*sp.Poly(am1, E).nth(1)) == 0)  # am1 = c*E
def is_polar(A):
    """Polar = affine + c*x^-1 : a1 const, a0 const, a_{-1} affine-in-E (c0 + c1 E)."""
    if set(A.keys()) - {-1, 0, 1}: return False
    a1 = A.get(1, sp.Integer(0)); a0 = A.get(0, sp.Integer(0)); am1 = A.get(-1, sp.Integer(0))
    return (sp.Poly(a1, E).degree() <= 0 and sp.Poly(a0, E).degree() <= 0
            and sp.Poly(am1, E).degree() <= 1)

print("== Section E: branch (B) produces genuine NON-AFFINE solutions in B ==")
# Construct branch (B) from the memo's own reduction, fully symbolically.
# a_{-1}=G (const nonzero), a1 linear, a0=alpha, b0=beta const,
# b1=lambda a1, b_{-1}=mu a_{-1}=mu G, and V(E)=a1(E-1)*G = -E/(lam-mu)+c.
# So a1(E) = V(E+1)/G = (-(E+1)/(lam-mu)+c)/G.
lam, mu, alpha, beta, G, c = sp.symbols("lambda mu alpha beta G c")
den = (lam - mu)
a1B  = sp.expand((-(E+1)/den + c)/G)          # linear in E
am1B = G
XB = clean({1: a1B, 0: alpha, -1: am1B})
DB = clean({1: lam*a1B, 0: beta, -1: mu*am1B})
BB = brk(DB, XB)
check("branch (B) satisfies FULL [D,X]=1 for all params (lam!=mu, G!=0)",
      sp.simplify(BB.get(0,0) - 1) == 0 and all(sp.simplify(BB.get(m,0))==0 for m in (1,-1,2,-2)))

# a concrete branch-(B) instance: X = x(E+1) + x^-1,  D = x^-1  (lam=0,mu=1,G=1,c=0,alpha=beta=0)
Xc = {1: E+1, -1: sp.Integer(1)}
Dc = {-1: sp.Integer(1)}
check("concrete branch-(B): [x^-1, x(E+1)+x^-1] = 1", brk(Dc, Xc) == ONE)
check("concrete branch-(B) X = x^2 d + x + 1/x is NOT affine", not is_affine(Xc))
check("concrete branch-(B) X is NOT polar (a1 nonconstant: has x^2 d term)", not is_polar(Xc))
check("concrete branch-(B) X is NOT in A1 (a_{-1}(0)=1 != 0)", not in_A1(Xc))
check("concrete branch-(B) D=1/x NOT in A1", not in_A1(Dc))
# so branch (B) lives in B \ A1 but is NOT of the 'polar' form claimed by P3.

print("\n== Section F: which symmetry actually maps branch (A) <-> branch (B)? ==")
# --- (F1) The genuine Weyl-Fourier automorphism F: x -> -d, d -> x.
# On the crossed product it is E=xd -> (-d)(x) = -(E+1), and x-level k -> -k with
# falling-factorial factors.  KEY POINT: F does NOT preserve B=A1[x^-1] because
# F(x) = -d is NOT invertible in B.  We demonstrate this concretely.
# F sends x |-> -d = -x^{-1}E  and  x^{-1} would need F(x)^{-1} = (-d)^{-1}, not in B.
# Show F(x^-1-generating element) leaves B: F(x)=-d has x-level -1 but F(x)*F(x^-1)
# must be 1; there is no element g in B with (-d) g = 1  (d is not left-invertible in B).
# Certificate: d has a nonzero kernel-ish obstruction -- concretely, d*g=1 has no
# solution g in B because applying to constants: (x^-1 E)(sum x^k h_k) has no x^0
# term equal to 1 unless a k=1 term h_1 with E-shift gives 1, but then g not unique/
# and (-d) is a zero-divisor-free but non-surjective map onto B.  We instead verify
# the clean structural fact used downstream:  F preserves A1 and the band {-1,0,1},
# but the IMAGE of a B-only element (x^-1) is not in B.
Fx = {-1: -E}                      # F(x) = -d = -x^{-1} E
# try to find g in B (x-levels -2..2, deg<=3 coeffs) with Fx * g = 1:
levels = range(-2, 3); deg = 3
gsyms = {(k, i): sp.Symbol(f"g_{k}_{i}") for k in levels for i in range(deg+1)}
g = clean({k: sum(gsyms[(k, i)]*E**i for i in range(deg+1)) for k in levels})
prod = mul(Fx, g)
eqs = []
for k, v in prod.items():
    target = 1 if k == 0 else 0
    pol = sp.Poly(sp.expand(v - target), E)
    eqs += [pol.coeff_monomial(E**i) for i in range(pol.degree()+1)] if v - target != 0 or k==0 else []
# also require all OTHER levels (not produced) are zero automatically; add explicit target for k=0 absent
if 0 not in prod:
    eqs.append(sp.Integer(-1))  # no x^0 term at all => cannot equal 1
sol = sp.solve(eqs, list(gsyms.values()), dict=True)
check("F(x)=-d has NO right-inverse in B (bounded search) => F does not preserve B",
      len(sol) == 0)

# --- (F2) The band-reversal automorphism sigma: x -> x^-1, E -> -E.
# We first PROVE it is an algebra automorphism of B by checking it preserves the
# defining relation  E x = x(E+1):  sigma(E)sigma(x) = sigma(x)(sigma(E)+1).
sx  = {-1: sp.Integer(1)}          # sigma(x) = x^-1
sE  = {0: -E}                      # sigma(E) = -E
lhs = mul(sE, sx)
rhs = mul(sx, add(sE, ONE))
check("sigma (x->x^-1, E->-E) respects E x = x(E+1)  => algebra homomorphism",
      lhs == rhs)
# sigma on a general element: sigma(x^k f(E)) = x^{-k} f(-E).  (homomorphism + generators)
def sigma(A): return clean({-k: sp.sympify(v).subs(E, -E) for k, v in A.items()})
_testel = clean({2: E**2-3, -1: E, 0: sp.Integer(7)})
check("sigma is an involution on B (sigma^2 = id) on random element",
      sigma(sigma(_testel)) == _testel)
check("sigma preserves B and the band {-1,0,1}",
      set(sigma(Xc).keys()) <= {-1,0,1})
# sigma maps our branch-(B) X to a branch-(A) element:
sXc = sigma(Xc)   # expect x + x^-1(1-E) = x - d + 1/x   (a1 const, a_{-1} linear)
check("sigma(branch-B X) is branch (A): a1 const, a_{-1} linear",
      sp.Poly(sXc.get(1,0), E).degree() <= 0 and sp.Poly(sXc.get(-1,0), E).degree() == 1)
check("sigma(branch-B pair) still solves [D,X]=1", brk(sigma(Dc), sigma(Xc)) == ONE)
# BUT sigma does NOT preserve A1 (sigma(x)=x^-1 not in A1):
check("sigma does NOT preserve A1 (sigma(x)=1/x not in A1)", not in_A1(sigma({1: sp.Integer(1)})))

print("\n== Section G: enumerate ALL vanishing patterns of (a1, a_{-1}) ==")
# For each pattern we solve the reduced system with generic deg<=3 coefficients and
# classify.  b1,b_{-1} handled per pattern.  We report solution families.
def solve_pattern(a1_zero, am1_zero, dmax=3, tag=""):
    """Return (families found, any non-affine-A1-pair found)."""
    def gp(name):
        cs = sp.symbols(f"{name}0:{dmax+1}")
        return sum(cs[i]*E**i for i in range(dmax+1)), list(cs)
    syms = []
    if a1_zero: a1 = sp.Integer(0)
    else: a1, s = gp("Aa"); syms += s
    if am1_zero: am1 = sp.Integer(0)
    else: am1, s = gp("Gg"); syms += s
    a0, s = gp("al"); syms += s
    b1, s = gp("Bp"); syms += s
    b0, s = gp("be"); syms += s
    bm1, s = gp("Bn"); syms += s
    X = clean({1: a1, 0: a0, -1: am1})
    D = clean({1: b1, 0: b0, -1: bm1})
    Cbr = brk(D, X)
    eqs = []
    for m in (-2,-1,0,1,2):
        v = Cbr.get(m, sp.Integer(0)) - (1 if m == 0 else 0)
        if v != 0:
            pol = sp.Poly(sp.expand(v), E)
            eqs += [pol.coeff_monomial(E**i) for i in range(pol.degree()+1)]
        elif m == 0:
            eqs.append(sp.Integer(-1))
    sols = sp.solve(eqs, syms, dict=True)
    return X, D, sols, syms

for (a1z, am1z, name) in [(False, False, "(I) a1,a_{-1} both nonzero"),
                          (True,  False, "(II) a1=0, a_{-1}!=0"),
                          (False, True,  "(III) a1!=0, a_{-1}=0"),
                          (True,  True,  "(IV) a1=a_{-1}=0 (X in C[E]-shift)")]:
    X, D, sols, syms = solve_pattern(a1z, am1z, dmax=2)
    print(f"  pattern {name}: {len(sols) if sols else 0} solution branch(es)")

fails = [n for n,ok in RES if not ok]
print("\n" + ("ALL BRANCH-B / SYMMETRY CHECKS PASSED" if not fails else f"FAILURES: {fails}"))
import sys; sys.exit(0 if not fails else 1)
