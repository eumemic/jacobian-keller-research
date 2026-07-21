#!/usr/bin/env python3
"""
verify_band5_comparison.py
==========================
Exact SymPy verification backing `band5-comparison.md`: the BAND-5 comparison
against band-3 (prime k=3) and band-4 (composite k=4), isolating what drives the
quantum "exotic slack" -- primality of the necklace operator S_k = (S^k-1)/(S-1)
versus bare k >= 3.

Frozen conventions (identical to the band-3 corpus; commits 99fe6ee, ebfc64d,
b9f9cf3, 9fa9f74, 050a4c0):
    A_1[x^{-1}] = (+)_k x^k C[E],  (x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E),
    f^[r](E) = f(E+r),  E = x d,  ladder-m coeff of [D,X]:
    Q_m = sum_{k+l=m} ( b_l^[k] a_k - a_k^[l] b_l ),  [D,X]=1 <=> Q_m = delta_{m0}.
    Genuine A_1 membership: E(E-1)...(E-r+1) | a_{-r}, b_{-r}  (falling factorial).
    Gauge b_k = 0 (spent on the top, Q_{2k}).

Root-necklace calculus (S = shift on root positions; sigma = T^{-1} = S^{-1}):
    S_r = 1 + s + ... + s^{r-1}  (necklace operator);  Phi_k the k-th cyclotomic.

Sections:
  0. Engine at k=5: Q_m == direct crossed-product commutator; top proportionality
     Q_10; the wall Q_9 = b_4^[5] a_5 - a_5^[4] b_4 in gauge b_5=0.
  1. The k=5 necklace lemma: WALL-M reduction  S*(1+S+S^2+S^3)*A = Phi_5(S)*B, i.e.
     S_5 delta(u) = S_4 delta(a_5); coprimality gcd(S_5,S_4)=1; shifted-5th-power
     sufficiency; the exotic gap via the universal cofactor g = 1 - sigma + sigma^2.
  2. PRIMALITY DISSECTION: S_5 = Phi_5 IRREDUCIBLE (5 prime) vs S_4 = Phi_2*Phi_4
     REDUCIBLE (4 composite) vs S_3 = Phi_3 irreducible (3 prime).  Independent
     k=4 wall Q_7 = b_3^[4]a_4 - a_4^[3]b_3 and its necklace S_4 delta = S_3 delta.
  3. ENUMERATION of minimal-degree exotic tops: deg a_k = k, single coset.
     k=3 -> 1 exotic ({0,2,4}); k=4 -> 4; k=5 -> 13 (stable).  Cofactor structure;
     reflection pairing; step-2 AP {0,2,..,2(k-1)} admissible IFF k ODD (cofactor
     Phi_{2k}); non-cyclotomic cofactors at k=5; the composite Phi_2-obstruction.
  4. MOMENT-UNIT test at k=5: W4 potential G band-agnostic (Q_0=(T-1)G); G(0)=0
     under membership; slope = G(1) = const coeff Q_0.  Positive cascade Q_8..Q_1
     forward-solves b_3..b_{-4}; {cascade} u {Q_0=1} INFEASIBLE (Groebner=[1]) while
     {Q_0=0} feasible -- the moment carries no unit -- across the exotic families,
     d=1 and d=2, a5(0)=0 and a5(0)!=0.  Positive control (genuine band-5 pair) +
     no-spurious-conditions guard.

Run:  uv run --with sympy python research/dc1-program/verify_band5_comparison.py
Ends: ALL BAND5 COMPARISON CHECKS PASSED
"""
import sympy as sp
import time
from itertools import combinations

E, S, sig = sp.symbols('E S sigma')


# ---------------------------------------------------------------- primitives
def sh(f, r):
    return sp.expand(sp.sympify(f).subs(E, E + r))


def poly(name, deg):
    if deg < 0:
        return sp.Integer(0), []
    cs = list(sp.symbols(f'{name}_0:{deg+1}'))
    return sp.expand(sum(cs[i] * E**i for i in range(deg + 1))), cs


def falling(r):
    return sp.prod([E - i for i in range(r)]) if r > 0 else sp.Integer(1)


def Qm(a, b, m, K):
    return sp.expand(sum(sh(b[l], k) * a[k] - sh(a[k], l) * b[l]
                         for k in range(-K, K + 1) for l in range(-K, K + 1) if k + l == m))


def gauged(K):
    return {k: sp.Integer(0) for k in range(-K, K + 1)}


def Sr(r):
    """necklace operator S_r = 1 + sigma + ... + sigma^{r-1}."""
    return sum(sig**j for j in range(r))


def SrS(r):
    return sum(S**j for j in range(r))


def az(expr, label):
    e = sp.expand(sp.sympify(expr))
    if e != 0:
        raise AssertionError(label + "  :  residual = " + str(e))
    print("PASS", label)


def ok(cond, label):
    if not cond:
        raise AssertionError(label + "  :  FAILED")
    print("PASS", label)


# =====================================================================
print("=" * 72)
print("0. Engine at k=5: Q_m == commutator; top proportionality Q_10; the wall Q_9")
print("=" * 72)
K = 5
Ag = {k: poly(f'A{k+5}', 2)[0] for k in range(-K, K + 1)}
Bg = {k: poly(f'B{k+5}', 2)[0] for k in range(-K, K + 1)}


def direct_commutator(m, K):
    return sp.expand(sum(sh(Bg[l], k) * Ag[k] - sh(Ag[k], l) * Bg[l]
                         for k in range(-K, K + 1) for l in range(-K, K + 1) if k + l == m))


for m in range(-10, 11):
    az(direct_commutator(m, K) - Qm(Ag, Bg, m, K), f"k=5  Q_{m} = direct crossed-product commutator")

# top Q_10: only (5,5) survives -> b_5^[5] a_5 - a_5^[5] b_5 (pure shifted-Wronskian)
az(Qm(Ag, Bg, 10, K) - (sh(Bg[5], 5) * Ag[5] - sh(Ag[5], 5) * Bg[5]),
   "k=5  Q_10 = b_5^[5] a_5 - a_5^[5] b_5  => b_5 = lambda_5 a_5 (5-periodicity), gauge b_5=0")

# the wall Q_9 in gauge b_5 = 0: only (5,4) survives -> b_4^[5] a_5 - a_5^[4] b_4
Bg0 = dict(Bg); Bg0[5] = sp.Integer(0)
az(Qm(Ag, Bg0, 9, K) - (sh(Bg0[4], 5) * Ag[5] - sh(Ag[5], 4) * Bg0[4]),
   "k=5  Q_9 = b_4^[5] a_5 - a_5^[4] b_4   (THE BAND-5 WALL, gauge b_5=0)")


# =====================================================================
print("\n" + "=" * 72)
print("1. The k=5 necklace lemma  S_5 delta(u) = S_4 delta(a_5)  and the exotic gap")
print("=" * 72)
# WALL-M: from b4(E+5)a5(E) = a5(E+4)b4(E), equate root multisets (S = +1 on positions):
#   S^{-5}B + A = S^{-4}A + B  <=>  (1-S^{-4})A = (1-S^{-5})B  <=> (xS^5)  S(S^4-1)A=(S^5-1)B
#   <=> cancel (S-1):  S*(1+S+S^2+S^3)*A = Phi_5(S)*B.
Phi5 = 1 + S + S**2 + S**3 + S**4
S4S = 1 + S + S**2 + S**3
# check the algebraic reduction on the multiset operators (identity of Laurent polys).
# Raw balance from equating root multisets of  b4(E+5)a5(E) = a5(E+4)b4(E):
#   S^{-5}B + A = S^{-4}A + B  <=>  (1 - S^{-4})A = (1 - S^{-5})B.  Multiply by S^5:
#   (S^5 - S)A = (S^5 - 1)B  <=>  S(S-1)(1+S+S^2+S^3)A = (S-1)Phi_5 B  <=> (cancel S-1)
#   S*(1+S+S^2+S^3)*A = Phi_5*B.
Avar, Bvar = sp.symbols('Avar Bvar')
raw = sp.expand((1 - S**(-4)) * Avar - (1 - S**(-5)) * Bvar)   # the raw balance (=0 form)
az(sp.expand(S**5 * raw - (S - 1) * (S * S4S * Avar - Phi5 * Bvar)),
   "k=5  WALL-M: S^5*(raw balance) = (S-1)*(S*(1+S+S^2+S^3)*A - Phi_5*B)  [cancel S-1 => WALL-M]")
# necklace form in sigma = S^{-1}:  S_5(sigma) delta(u) = S_4(sigma) delta(a_5)
az(Sr(5) * Sr(4) - Sr(4) * Sr(5), "k=5  necklace operators commute (S_5 S_4 = S_4 S_5)")
ok(sp.gcd(sp.Poly(Sr(5), sig), sp.Poly(Sr(4), sig)).as_expr() == 1,
   "k=5  gcd(S_5, S_4) = 1 in Q[sigma]  (coprime necklace operators)")

# shifted-5th-power (cube) sufficiency: a5 = prod_{j=0}^4 h^[j], u = prod_{j=0}^3 h^[j] solves the wall
h = sp.Function('h')
a5cube = sp.prod([h(E + j) for j in range(5)]); b4cube = sp.prod([h(E + j) for j in range(4)])
az(sh(b4cube, 5) * a5cube - sh(a5cube, 4) * b4cube,
   "k=5  shifted-5th-power sufficiency: a5=prod_0^4 h^[j], u=prod_0^3 h^[j] solves the wall")

# THE EXOTIC GAP.  Necklace equation is  S_5 delta(u) = S_4 delta(a5).  A solution
# forces delta(u)=S_4 g, delta(a5)=S_5 g for a cofactor g in Z[sigma]; a nonzero u
# exists iff BOTH S_4 g and S_5 g are effective.  The universal weapons-cofactor
# g = 1 - sigma + sigma^2 (works for every k>=3) is NON-effective yet gives both effective:
g = 1 - sig + sig**2
Sg5 = sp.expand(Sr(5) * g); Sg4 = sp.expand(Sr(4) * g)
ok(any(c < 0 for c in sp.Poly(g, sig).all_coeffs()), "k=5  universal cofactor g=1-sigma+sigma^2 is NON-effective")
ok(all(c >= 0 for c in sp.Poly(Sg5, sig).all_coeffs()), "k=5  S_5 g effective => a_5 genuine  (S_5 g = " + str(Sg5) + ")")
ok(all(c >= 0 for c in sp.Poly(Sg4, sig).all_coeffs()), "k=5  S_4 g effective => u=b_4 genuine (S_4 g = " + str(Sg4) + ")")
print("   => the exotic (non-shifted-5th-power) branch is NON-EMPTY at k=5, exactly as at k=3.")

# the WALL-DEG law: evaluate S_5 delta(u) = S_4 delta(a5) at sigma=1: 5|B|... (k-1)|A|=k|B|
ok(Sr(5).subs(sig, 1) == 5 and Sr(4).subs(sig, 1) == 4,
   "k=5  WALL-DEG: (k-1) deg a_5 = k deg b_4, i.e. 4 deg a_5 = 5 deg b_4 (5 | deg a_5)")


# =====================================================================
print("\n" + "=" * 72)
print("2. PRIMALITY DISSECTION: S_k = (S^k-1)/(S-1) irreducible <=> k prime")
print("=" * 72)
ok(sp.factor(Sr(3)) == sig**2 + sig + 1, "k=3 prime : S_3 = Phi_3 = 1+s+s^2 IRREDUCIBLE")
ok(sp.factor(Sr(5)) == sig**4 + sig**3 + sig**2 + sig + 1, "k=5 prime : S_5 = Phi_5 IRREDUCIBLE")
ok(sp.factor(Sr(4)) == (sig + 1) * (sig**2 + 1), "k=4 comp. : S_4 = Phi_2*Phi_4 = (1+s)(1+s^2) REDUCIBLE")
# S_k always has all-1 coefficients regardless of primality; primality controls FACTORIZATION only:
for k in (3, 4, 5):
    ok(all(c == 1 for c in sp.Poly(Sr(k), sig).all_coeffs()),
       f"k={k}  S_k = 1+sigma+...+sigma^(k-1) has all-1 coeffs (primality touches factorization, not shape)")

# Independent k=4 wall (do NOT coordinate with the sibling; derive from scratch):
K4 = 4
A4 = {k: poly(f'C{k+4}', 2)[0] for k in range(-K4, K4 + 1)}
B4 = {k: poly(f'D{k+4}', 2)[0] for k in range(-K4, K4 + 1)}
B4[4] = sp.Integer(0)
az(Qm(A4, B4, 7, K4) - (sh(B4[3], 4) * A4[4] - sh(A4[4], 3) * B4[3]),
   "k=4  Q_7 = b_3^[4] a_4 - a_4^[3] b_3  (independent band-4 wall, gauge b_4=0)")
ok(sp.gcd(sp.Poly(Sr(4), sig), sp.Poly(Sr(3), sig)).as_expr() == 1,
   "k=4  necklace S_4 delta(u)=S_3 delta(a_4); gcd(S_4,S_3)=1")
# composite obstruction: divisibility 'S_4 | A' = (Phi_2|A AND Phi_4|A), a CONJUNCTION:
Phi2, Phi4 = 1 + S, 1 + S**2
partial = 0
for combo in combinations(range(1, 13), 3):
    roots = [0] + list(combo)
    A = sum(S**r for r in roots)
    if sp.rem(A, Phi4, S) == 0 and sp.rem(A, Phi2, S) != 0:
        partial += 1
ok(partial > 0, f"k=4  {partial} deg-4 tops are Phi_4-divisible but NOT Phi_2-divisible => NO wall (composite-only obstruction)")
print("   => at k prime the wall-divisibility is a SINGLE irreducible condition Phi_k|A;")
print("      at k composite it is a CONJUNCTION Phi_d|A (d|k, d>1) -- new obstructions appear.")


# =====================================================================
print("\n" + "=" * 72)
print("3. ENUMERATION: minimal-degree (deg a_k = k) exotic tops, single coset")
print("=" * 72)
def enumerate_exotic(k, span):
    Phik = SrS(k); Skm1 = SrS(k - 1)
    exotic = []; cube = []
    for combo in combinations(range(1, span + 1), k - 1):
        roots = [0] + list(combo)
        A = sum(S**r for r in roots)
        q, rem = sp.div(sp.expand(A), Phik, S)
        if sp.expand(rem) != 0:
            continue
        cof = sp.Poly(q, S); nc = any(cof.nth(i) < 0 for i in range(cof.degree() + 1))
        qb, remb = sp.div(sp.expand(S * Skm1 * A), Phik, S); Bp = sp.Poly(qb, S)
        if sp.expand(remb) != 0 or not all(Bp.nth(i) >= 0 for i in range(Bp.degree() + 1)):
            continue
        (exotic if nc else cube).append((roots, q.as_expr()))
    return exotic, cube


counts = {}
for k, span in [(3, 24), (4, 24), (5, 24)]:
    ex, cu = enumerate_exotic(k, span)
    counts[k] = len(ex)
    print(f"   k={k}: {len(ex)} exotic + {len(cu)} cube (span<={span}, stable)")
ok(counts[3] == 1, "k=3: exactly 1 minimal exotic top up to translation ({0,2,4}) -- the unique step-2 AP")
ok(counts[4] == 4, "k=4: exactly 4 minimal exotic tops up to translation")
ok(counts[5] == 13, "k=5: exactly 13 minimal exotic tops up to translation (13x richer than k=3, both PRIME)")
# stability: k=5 count unchanged from span 20 to 30
ok(len(enumerate_exotic(5, 20)[0]) == 13 and len(enumerate_exotic(5, 30)[0]) == 13,
   "k=5: exotic count STABLE at 13 (span 20 == span 30) -- the both-effective condition bounds the Newton span")
# the counts 1, 4, 13 match (3^(k-2)-1)/2 -- a smooth bare-k law BLIND to primality
# (k=4 composite and k=5 prime both fit).  k=6 -> 40 verified off-verifier (~10 s, span 22).
ok(all(counts[k] == (3 ** (k - 2) - 1) // 2 for k in (3, 4, 5)),
   "k=3,4,5: exotic count = (3^(k-2)-1)/2 = 1,4,13 (CONJECTURAL bare-k law; k=6->40 off-verifier; primality-blind)")

# the unique k=3 exotic has cofactor Phi_6; among the 13 at k=5 the cofactors are NOT all cyclotomic:
ex5, _ = enumerate_exotic(5, 24)


def _is_cyclo(poly_expr):
    """True iff poly_expr equals some cyclotomic Phi_n (monic, degree <= 8 suffices here)."""
    d = sp.Poly(sp.expand(poly_expr), S).degree()
    for n in range(1, 60):
        cn = sp.cyclotomic_poly(n, S)
        if sp.Poly(cn, S).degree() == d and sp.expand(poly_expr - cn) == 0:
            return True
    return False


def has_noncyclotomic_factor(c):
    """True iff some monic irreducible factor of c is not any cyclotomic Phi_n."""
    for fac, _mult in sp.factor_list(c)[1]:
        if not _is_cyclo(fac):
            return True
    return False


ok(any(has_noncyclotomic_factor(c) for _, c in ex5),
   "k=5: some exotic cofactors are NON-cyclotomic (e.g. S^3-S^2+1) -- the class is tropical, not cyclotomic")

# step-2 AP {0,2,...,2(k-1)}: wall-admissible IFF k ODD; when admissible cofactor = Phi_{2k}
for k in (3, 4, 5, 6, 7):
    roots = list(range(0, 2 * k, 2))
    A = sum(S**r for r in roots)
    q, rem = sp.div(sp.expand(A), SrS(k), S)
    admissible = sp.expand(rem) == 0
    ok(admissible == (k % 2 == 1),
       f"k={k} ({'odd' if k%2 else 'even'}): step-2 AP wall-admissible = {admissible}  (== k odd)")
    if admissible:
        ok(sp.expand(q.as_expr() - sp.cyclotomic_poly(2 * k, S)) == 0,
           f"k={k}: admissible step-2 AP has cofactor EXACTLY Phi_{2*k}")
print("   => the k=3 'step-2 AP' exotic lifts to every ODD k (cofactor Phi_{2k}); EVEN k obstructs it")
print("      via Phi_2-divisibility failure.  Count: k=3->1, k=4->4, k=5->13 (grows with k, not primality).")


# =====================================================================
print("\n" + "=" * 72)
print("4. MOMENT-UNIT TEST at k=5: the exotic branch is killed at Q_0 (no unit)")
print("=" * 72)


def potential_G(A, B, K):
    return sp.expand(sum(sh(A[k], j - k) * sh(B[-k], j) - sh(B[k], j - k) * sh(A[-k], j)
                         for k in range(1, K + 1) for j in range(0, k)))


# W4 potential band-agnostic at k=5:
az(Qm(Ag, Bg, 0, K) - (sh(potential_G(Ag, Bg, K), 1) - potential_G(Ag, Bg, K)),
   "k=5  Q_0 = (T-1)G  (W4 telescoping lemma, band-agnostic potential, generic band-5)")
# G(0)=0 under membership:
Am = {k: poly(f'Em{k+5}', 2)[0] for k in range(0, K + 1)}
Bm = {k: poly(f'Fm{k+5}', 2)[0] for k in range(0, K + 1)}
for k in range(1, K + 1):
    Am[-k] = sp.expand(falling(k) * poly(f'amt{k}', 2)[0])
    Bm[-k] = sp.expand(falling(k) * poly(f'bmt{k}', 2)[0])
Gm = potential_G(Am, Bm, K)
az(Gm.subs(E, 0), "k=5  G(0) = 0 identically under membership  =>  Q_0=1 <=> G=E")
az(sp.Poly(sp.expand(sh(Gm, 1) - Gm), E).nth(0) - Gm.subs(E, 1),
   "k=5  slope := const coeff of Q_0 = G(1)  (the 'moment unit' slot)")
# membership-protection of the level-5 (bottom) block: its slope contribution is the
# SINGLE surviving product mu5 * a5(0) * a_-5(5)  (all a_-5(0..4) killed by E^underline5):
a5g, _ = poly('gA', 5)
am5t, _ = poly('gB', 3)
am5 = sp.expand(falling(5) * am5t)
mu5g = sp.symbols('mu5g')
bm5 = sp.expand(mu5g * am5)
G5blk = sp.expand(sum(sh(a5g, j - 5) * sh(bm5, j) for j in range(5)))   # gauge b_5=0
az(G5blk.subs(E, 0), "k=5  level-5 block G5(0) = 0 (a_-5(0..4) killed by membership)")
az(G5blk.subs(E, 1) - mu5g * a5g.subs(E, 0) * am5.subs(E, 5),
   "k=5  level-5 block G5(1) = mu5 * a5(0) * a_-5(5)  (the ONLY surviving bottom slope term)")


# ---- positive cascade solver (generalized to band 5) ----
def clean_solve(A, B, m, lkey, name, memb, bdeg_raw, K):
    braw, bc = poly(f'{name}c_', bdeg_raw)
    bnew = sp.expand(falling(memb) * braw)
    Bt = dict(B); Bt[lkey] = bnew
    eq = sp.expand(Qm(A, Bt, m, K))
    coeffs = sp.Poly(eq, E).all_coeffs() if eq != 0 else [sp.Integer(0)]
    M, rhs = sp.linear_eq_to_matrix(coeffs, bc)
    if not all(len(e.free_symbols) == 0 for e in M):
        raise AssertionError("operator matrix not numeric (bilinearity leaked m=%d)" % m)
    conds = [sp.expand(n.dot(rhs)) for n in M.T.nullspace()]
    conds = [c for c in conds if c != 0]
    Mred = sp.zeros(0, len(bc)); rhs_red = []
    for i in range(M.shape[0]):
        cand = Mred.col_join(M[i, :])
        if cand.rank() > Mred.rank():
            Mred = cand; rhs_red.append(rhs[i])
    if Mred.shape[0] == 0:
        bc_vals = [sp.Integer(0)] * len(bc)
    else:
        sol, _ = Mred.gauss_jordan_solve(sp.Matrix(rhs_red))
        tau = [s for s in sol.free_symbols if str(s).startswith('tau')]
        bc_vals = [x.subs({t: 0 for t in tau}) for x in sol]
    bsol = sp.expand(bnew.subs(dict(zip(bc, bc_vals))))
    ker = []
    for idx, kv in enumerate(M.nullspace()):
        gk = sp.symbols(f'{name}K{idx}')
        kp = sp.expand(falling(memb) * sum(kv[i] * E**i for i in range(len(bc))))
        bsol = sp.expand(bsol + gk * kp); ker.append(gk)
    return bsol, ker, conds


def build_cascade5(a5, b4, d):
    a4, ca4 = poly('a4', d); a3, ca3 = poly('a3', d); a2, ca2 = poly('a2', d)
    a1, ca1 = poly('a1', d); a0, ca0 = poly('a0', d)
    negs = {}; cneg = {}
    for r in range(1, 6):
        raw, csr = poly(f'am{r}', d)
        negs[-r] = sp.expand(falling(r) * raw); cneg[-r] = csr
    mu5 = sp.symbols('mu5')
    A = {5: a5, 4: a4, 3: a3, 2: a2, 1: a1, 0: a0}
    for r in range(1, 6):
        A[-r] = negs[-r]
    B = gauged(5); B[4] = b4; B[-5] = sp.expand(mu5 * negs[-5])
    plan = [(8, 3, 'b3', 0), (7, 2, 'b2', 0), (6, 1, 'b1', 0), (5, 0, 'b0', 0),
            (4, -1, 'bm1', 1), (3, -2, 'bm2', 2), (2, -3, 'bm3', 3), (1, -4, 'bm4', 4)]
    pos = []; ks = []
    for (m, l, nm, mem) in plan:
        bb, kk, cc = clean_solve(A, B, m, l, nm, mem, 4 * d + 5, 5)
        B[l] = bb; pos += cc; ks += kk
    allvars = sum([list(ca4), list(ca3), list(ca2), list(ca1), list(ca0)], [])
    for r in range(1, 6):
        allvars += list(cneg[-r])
    allvars += [mu5] + ks
    return A, B, pos, allvars


def q0_conditions(A, B, target):
    ex = sp.expand(Qm(A, B, 0, 5) - target)
    return [sp.expand(c) for c in sp.Poly(ex, E).all_coeffs() if sp.expand(c) != 0] if ex != 0 else []


def b4_from_top(roots):
    A = sum(S**r for r in roots)
    qb, remb = sp.div(sp.expand(S * S4S * A), Phi5, S)
    assert sp.expand(remb) == 0
    Bp = sp.Poly(qb, S)
    br = [i for i in range(Bp.degree() + 1) if Bp.nth(i) != 0]
    return sp.expand(sp.prod([E - i for i in br])), br


# the exotic families (representatives of distinct cofactor types + a translate with a5(0)!=0):
FAMILIES = [
    ("{0,2,3,4,6} cofactor Phi_6  (universal g; a5(0)=0)", [0, 2, 3, 4, 6]),
    ("{0,2,4,6,8} step-2 AP cofactor Phi_10   (a5(0)=0)", [0, 2, 4, 6, 8]),
    ("{0,1,4,7,8} self-reflective cofactor Phi_12", [0, 1, 4, 7, 8]),
    ("{0,1,3,4,7} NON-cyclotomic cofactor S^3-S^2+1", [0, 1, 3, 4, 7]),
    ("{1,3,4,5,7} translate, a5(0)!=0 (bottom enters unit eqn)", [1, 3, 4, 5, 7]),
]
for tag, roots in FAMILIES:
    a5 = sp.expand(sp.prod([E - r for r in roots]))
    b4, _ = b4_from_top(roots)
    az(sh(b4, 5) * a5 - sh(a5, 4) * b4, f"wall OK: {tag}")
    A, B, pos, allvars = build_cascade5(a5, b4, 1)
    GU = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(1)), *allvars, order='grevlex'))
    GH = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(0)), *allvars, order='grevlex'))
    ok(GU == [sp.Integer(1)] and GH != [sp.Integer(1)],
       f"MOMENT KILL d=1: {tag}: {{cascade}}u{{Q_0=1}} INFEASIBLE, {{Q_0=0}} feasible => unit is the killer")

# one d=2 instance (free-degree beyond the base):
a5 = sp.expand(sp.prod([E - r for r in [0, 2, 3, 4, 6]])); b4, _ = b4_from_top([0, 2, 3, 4, 6])
A, B, pos, allvars = build_cascade5(a5, b4, 2)
GU = list(sp.groebner(pos + q0_conditions(A, B, sp.Integer(1)), *allvars, order='grevlex'))
ok(GU == [sp.Integer(1)], "MOMENT KILL d=2: {0,2,3,4,6}: {cascade}u{Q_0=1} INFEASIBLE (free-degree 2)")


# ---- positive control: a GENUINE band-5 pair, plus no-spurious-conditions guard ----
def mul(P, Q):
    R = {}
    for k1, p1 in P.items():
        for k2, p2 in Q.items():
            R[k1 + k2] = R.get(k1 + k2, 0) + sh(p1, k2) * p2
    return {k: sp.expand(v) for k, v in R.items() if sp.expand(v) != 0}


kap = sp.symbols('kappa')
U = {1: sp.Integer(1), -1: sp.expand(E)}
U5 = mul(mul(mul(mul(U, U), U), U), U)
Xc = gauged(5)
for k, v in U5.items():
    Xc[k] = Xc.get(k, 0) + v
Xc[-1] = sp.expand(Xc[-1] - E / kap)
Dc = gauged(5); Dc[1] = kap; Dc[-1] = sp.expand(kap * E)
ok(all(sp.expand(Qm(Xc, Dc, m, 5) - (1 if m == 0 else 0)) == 0 for m in range(-10, 11)),
   "positive control U=x+d, X=U^5-d/kappa, D=kappa U: [D,X]=1 (all Q_m=delta_{m0}) -- engine validated at band 5")
# guard: clean_solve reconstructs the genuine b_1 = kappa with NO spurious conditions
Av = dict(Xc)
for r in range(1, 6):
    Av.setdefault(-r, sp.Integer(0)); Av.setdefault(r, sp.Integer(0))
Bv = gauged(5); Bv[-1] = Dc[-1]
rb1, kb1, cv1 = clean_solve(Av, Bv, 6, 1, 'ctl', 0, 6, 5)
ok(len(cv1) == 0 and len(kb1) >= 1 and sp.expand(rb1.subs({kb1[0]: kap}) - kap) == 0,
   "positive control: clean_solve emits NO spurious conditions and reconstructs b_1=kappa (kill is real, not artifact)")


print("\nALL BAND5 COMPARISON CHECKS PASSED")
