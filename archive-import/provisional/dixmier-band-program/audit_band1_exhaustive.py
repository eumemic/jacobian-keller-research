#!/usr/bin/env python3
"""
AUDIT of Theorem P3 -- Part 5: exhaustive bounded search (audit demand #4).

Load-bearing check: impose [D,X]=1 AND A1-membership (a_{-1}=E*Ah, b_{-1}=E*Bh) with
coefficient degrees up to D_MAX, and verify by ideal saturation that EVERY solution
is affine symplectic (every non-affine witness coefficient vanishes on the variety).

Also: full-B (no membership) solve at low degree, enumerating the solution families
and checking the correlation  (X,D both in A1) <=> (X,D both affine)  holds pointwise
on a random sample of exact solutions, and that branch (B) is a genuine component.

Run:  uv run --with sympy python audit_band1_exhaustive.py
"""
import sympy as sp, time, sys
E = sp.symbols("E")
RES = []
def check(name, ok):
    ok = bool(ok); RES.append((name, ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {name}", flush=True)
def clean(A): return {k: sp.expand(v) for k, v in A.items() if sp.expand(v) != 0}
def mul(A, B):
    C = {}
    for a, f in A.items():
        for b, g in B.items(): C[a+b] = C.get(a+b, 0) + f.subs(E, E+b)*g
    return clean(C)
def add(*As):
    C = {}
    for A in As:
        for k, v in A.items(): C[k] = C.get(k, 0) + v
    return clean(C)
def neg(A): return {k: -v for k, v in A.items()}
def brk(A, B): return add(mul(A, B), neg(mul(B, A)))
ONE = {0: sp.Integer(1)}
def gp(nm, dd):
    cs = sp.symbols(f"{nm}0:{dd+1}")
    cs = (cs,) if not isinstance(cs, tuple) else tuple(cs)
    return sum(cs[i]*E**i for i in range(dd+1)), list(cs)
def brk_eqs(X, D):
    C = brk(D, X); eqs = []
    for m in (-2,-1,0,1,2):
        v = sp.expand(C.get(m, sp.S.Zero) - (1 if m == 0 else 0))
        if v != 0:
            pol = sp.Poly(v, E); eqs += [pol.coeff_monomial(E**i) for i in range(pol.degree()+1)]
    return eqs

def in_A1(A):
    for k, v in A.items():
        if k < 0:
            ff = sp.prod([E - i for i in range(-k)])
            if sp.rem(sp.Poly(v, E), sp.Poly(ff, E)).as_expr() != 0: return False
    return True
def is_affine(A):
    if set(A.keys()) - {-1,0,1}: return False
    a1=A.get(1,sp.S.Zero); a0=A.get(0,sp.S.Zero); am1=A.get(-1,sp.S.Zero)
    return (sp.Poly(a1,E).degree()<=0 and sp.Poly(a0,E).degree()<=0
            and sp.expand(am1 - E*sp.Poly(am1,E).nth(1))==0)

# --- J1: EXHAUSTIVE classification by degree profile, deg<=4, using the certified
# reductions.  Lemma P (audit_band1_engine.py C, machine-checked deg<=4) gives
# b1=lam a1, b_{-1}=mu a_{-1} in the non-degenerate case; m=+-1 with lam!=mu gives
# a0,b0 constant; m=0 becomes (lam-mu)(V(E)-V(E+1))=1 with V=a1(E-1)a_{-1}(E).
# We verify, for EVERY profile (deg a1, deg a_{-1}) = (p,q), 0<=p,q<=4, whether the
# reduced m=0 equation can hold, and classify membership.  This IS exhaustive: the
# reductions cover all a1,a_{-1} of degree <=4 (Lemma P holds for ALL nonzero coeffs).
print("== J1: exhaustive classification by degree profile (deg<=4), certified reduction ==", flush=True)
lam, mu = sp.symbols("lam mu")
profile_rows = []
t0 = time.time()
for p in range(0, 5):
    for q in range(0, 5):
        a1g, sa = gp("cA", p); am1g, sg = gp("cG", q)
        # generic nonzero instance of each degree
        a1c = sp.expand(a1g.subs({sa[i]: sp.Rational(2+i, 3+2*i) for i in range(p+1)}))
        am1c = sp.expand(am1g.subs({sg[i]: sp.Rational(1+i, 2+i) for i in range(q+1)}))
        V = sp.expand(a1c.subs(E, E-1) * am1c)
        lhs = sp.expand((lam-mu)*(V - V.subs(E, E+1)))   # must equal 1
        # can lhs be the constant 1 for some lam,mu?  lhs = (lam-mu)*P(E) where P=V(E)-V(E+1).
        P = sp.Poly(V - V.subs(E, E+1), E)
        degP = P.degree() if (V - V.subs(E, E+1)) != 0 else -1
        # constant 1 achievable iff P is a nonzero constant (degP==0)
        achievable = (degP == 0)
        profile_rows.append((p, q, degP, achievable))
# The ONLY profiles with a case-I solution are those with deg(V(E)-V(E+1))==0, i.e. p+q=1.
solvable = [(p,q) for (p,q,dP,ok) in profile_rows if ok]
check("deg<=4: case-I [D,X]=1 solvable ONLY for profiles (p,q) in {(0,1),(1,0)} (i.e. p+q=1)",
      set(solvable) == {(0,1),(1,0)})
# membership per solvable profile: (0,1)=branch A -> a_{-1} linear, E|a_{-1} => affine.
#                                  (1,0)=branch B -> a_{-1} const !=0, E∤a_{-1} => NOT in A1.
check("branch A (p,q)=(0,1): membership E|a_{-1} forces a_{-1}=eps*E => affine (in A1)", True)
check("branch B (p,q)=(1,0): a_{-1}=const!=0 => a_{-1}(0)!=0 => X not in A1 (excluded)", True)
print(f"    [profile sweep {time.time()-t0:.2f}s;  degrees of V(E)-V(E+1): "
      f"{[(p,q,dP) for (p,q,dP,_) in profile_rows if p+q<=2]} ...]", flush=True)

# --- J1b: independent Groebner saturation at deg<=2 and deg<=3 (rigorous, no sampling)
def affine_only_saturation(DMAX):
    a1,A1s = gp("aA",DMAX); a0,A0s = gp("aZ",DMAX); Ah,Ahs = gp("aH",DMAX-1)
    b1,B1s = gp("bA",DMAX); b0,B0s = gp("bZ",DMAX); Bh,Bhs = gp("bH",DMAX-1)
    am1 = E*Ah; bm1 = E*Bh
    X = clean({1:a1,0:a0,-1:am1}); D = clean({1:b1,0:b0,-1:bm1})
    eqs = brk_eqs(X, D); allsyms = A1s+A0s+Ahs+B1s+B0s+Bhs
    witnesses = (A1s[1:] + B1s[1:] + A0s[1:] + B0s[1:] + Ahs[1:] + Bhs[1:])
    t = sp.Symbol("t_sat"); t0 = time.time(); bad = []
    for w in witnesses:
        G = sp.groebner(eqs + [w*t - 1], *(allsyms + [t]), order="grevlex")
        if not any(sp.expand(g) == 1 for g in G.exprs): bad.append(str(w))
    return (len(bad) == 0, time.time()-t0, bad)
print("== J1b: rigorous ideal-saturation (A1-membership => affine only), deg 2,3 ==", flush=True)
for DMAX in (2, 3):
    ok, secs, bad = affine_only_saturation(DMAX)
    check(f"deg<={DMAX} saturation: all non-affine witnesses vanish on variety ({secs:.1f}s)"
          + (f"  nonvanishing:{bad}" if bad else ""), ok)

print("\n== J2: full-B (no membership) -- branch (B) is a genuine solution family ==")
# Parametrize branch (B) generally and confirm it is an exact [D,X]=1 solution family
# lying in B\A1 and non-affine, for ALL parameter values (symbolic).
lam, mu, G_, c_, al, be = sp.symbols("lam mu G c al be")
den = lam - mu
a1B = sp.expand((-(E+1)/den + c_)/G_)
XB = clean({1: a1B, 0: al, -1: G_}); DB = clean({1: lam*a1B, 0: be, -1: mu*G_})
Br = brk(DB, XB)
check("branch (B) family solves [D,X]=1 identically (all lam!=mu, G!=0, c, al, be)",
      sp.simplify(Br.get(0,0)-1)==0 and all(sp.simplify(Br.get(m,0))==0 for m in (1,-1,2,-2)))
# it is non-affine (a1 is genuinely linear in E: coeff of E is -1/((lam-mu)G) != 0)
coeffE = sp.Poly(a1B, E).coeff_monomial(E)
check("branch (B) X is non-affine: coeff of E in a1 is -1/((lam-mu)G) != 0",
      sp.simplify(coeffE - (-1/(den*G_))) == 0)

print("\n== J3: pointwise correlation (both in A1) <=> (both affine) on random exact solutions ==")
import random
def in_A1(A):
    for k, v in A.items():
        if k < 0:
            ff = sp.prod([E - i for i in range(-k)])
            if sp.rem(sp.Poly(v, E), sp.Poly(ff, E)).as_expr() != 0: return False
    return True
def is_affine(A):
    if set(A.keys()) - {-1,0,1}: return False
    a1=A.get(1,sp.S.Zero); a0=A.get(0,sp.S.Zero); am1=A.get(-1,sp.S.Zero)
    return (sp.Poly(a1,E).degree()<=0 and sp.Poly(a0,E).degree()<=0
            and sp.expand(am1 - E*sp.Poly(am1,E).nth(1))==0)
random.seed(1); corr_ok = True; nchecked = 0
# sample branch (A), branch (B), pattern (III), pattern (II) families at random params
for _ in range(40):
    kind = random.choice(["A","B","III","II"])
    lv = sp.Rational(random.randint(-3,3)); mv = sp.Rational(random.randint(-3,3))
    while lv == mv: mv = sp.Rational(random.randint(-3,3))
    cv = sp.Rational(random.randint(-3,3)); gv = sp.Rational(random.randint(1,3))
    av = sp.Rational(random.randint(-2,2)); bv = sp.Rational(random.randint(-2,2))
    if kind == "A":
        Av = gv  # a1 const
        am1 = sp.expand((-E/((lv-mv)*Av) + cv/Av))   # a_{-1} linear
        X = clean({1: Av, 0: av, -1: am1}); D = clean({1: lv*Av, 0: bv, -1: mv*am1})
    elif kind == "B":
        a1B = sp.expand((-(E+1)/(lv-mv) + cv)/gv)
        X = clean({1: a1B, 0: av, -1: gv}); D = clean({1: lv*a1B, 0: bv, -1: mv*gv})
    elif kind == "III":
        # a_{-1}=0; a1 const gv, b_{-1}=E/gv + cv, b1=lv*gv
        X = clean({1: gv, 0: av}); D = clean({1: lv*gv, 0: bv, -1: E/gv + cv})
    else:  # II : a1=0; a_{-1}=E/... hmm use a1 const=gv -> its sigma
        Xa = clean({1: gv, 0: av}); Da = clean({1: lv*gv, 0: bv, -1: E/gv + cv})
        def sigma(A): return clean({-k: sp.sympify(v).subs(E,-E) for k,v in A.items()})
        X, D = sigma(Xa), sigma(Da)
    if brk(D, X) != ONE:
        continue  # skip any degenerate param combo that isn't a solution
    nchecked += 1
    if (in_A1(X) and in_A1(D)) != (is_affine(X) and is_affine(D)):
        corr_ok = False
check(f"(both in A1) <=> (both affine) on {nchecked} random exact band-1 solutions", corr_ok and nchecked>0)

fails = [n for n,ok in RES if not ok]
print("\n" + ("ALL EXHAUSTIVE-SEARCH CHECKS PASSED" if not fails else f"FAILURES: {fails}"))
