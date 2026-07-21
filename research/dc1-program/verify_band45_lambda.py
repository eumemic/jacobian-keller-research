#!/usr/bin/env python3
"""
verify_band45_lambda.py
=======================
Exact SymPy verification backing `band45-lambda.md`: the ANNIHILATOR TEMPLATE
applied to the band-4 and band-5 minimal exotic tops.

The overnight W1 breakthrough (`b7e85e8..e4e704f`) killed the normalized band-3
W1 datum by an annihilating functional `lambda_0` of the two-filler image `Phi`
with `lambda_0(E-R) != 0`.  This file establishes, DEGREE-FREE and UNIFORMLY, the
closed-form annihilator for every band-4/5 exotic top, and reports the intrinsic
identity of `lambda`:  it is built from the top's necklace COFACTOR
`C = A / S_k`  (A = root necklace of a_k, S_k = 1+S+..+S^{k-1}), namely

        hat_lambda(t) = t^k C(t) - C(1),      lambda = sum_p [t^p] hat_lambda * ev_p.

The two block-annihilation conditions are the two generating-function support
identities
        hat_lambda * S_k     = t^k A       - C(1) S_k       (supp {0..k-1} u {k+rho})
        hat_lambda * S_{k-1} = t^{k-1} B   - C(1) S_{k-1}   (supp {0..k-2} u {k-1+sig})
where B = t S_{k-1} C is the wall root-necklace of u=b_{k-1}.  Both are proved
here as EXACT polynomial identities in t, and the resulting
    lambda(K_k[(E)_k C_poly]) = 0,   lambda(H_{k-1}[(E)_{k-1} V_poly]) = 0
are checked as SYMBOLIC identities in the free coefficients of C_poly,V_poly
(arbitrary degree).  Consequently  Im Phi subset ker lambda  degree-free.

Closed forms proved:  C(1)=1;  lambda(1)=0;  lambda(E) = (k+1)/2 + (sum rho)/k > 0
for every normalized top (all roots >= 0), so  E not in Im Phi  arbitrary degree.

Reflection E -> -E-1: the reflected top's cofactor is the REVERSAL of C; self-dual
tops have palindromic C; the pair {0,1,3,6}~{0,3,5,6} maps by reversal.

The residual step  lambda(R) = 0  (which upgrades  E not in Im Phi  to the full
moment-unit kill  E - R not in Im Phi  <=>  Q_0 = 1 impossible) is CASCADE-
DEPENDENT and is the single open arbitrary-degree step (exactly as flagged for
general k in lambda-general-k.md).  It is verified here as BOUNDED evidence: for
every top, after forward-solving the positive cascade at free degree d, lambda(G)
= lambda(R) reduces to the constant 0 modulo the positive-cascade ideal, while
lambda(E) != 0 -- so lambda(E-R) = lambda(E) != 0 and Q_0=1 is infeasible at that
degree.  These are per-degree Groebner facts, not an arbitrary-degree theorem.

Conventions (frozen, identical to verify_band4_experiment.py /
verify_band5_comparison.py):
    A_1[x^{-1}] = (+)_k x^k C[E],  (x^a f)(x^b g) = x^{a+b} f(E+b) g(E),
    f^[r](E)=f(E+r),  Q_m = sum_{k+l=m} (b_l^[k] a_k - a_k^[l] b_l),
    membership E(E-1)..(E-r+1) | a_{-r},b_{-r};  gauge b_k=0.
    closure-form potential
      G = sum_{k=1}^K sum_{j=0}^{k-1} [ a_k^[j-k] b_{-k}^[j] - b_k^[j-k] a_{-k}^[j] ],
      Q_0 = (T-1)G,  G(0)=0 under membership  =>  Q_0=1 <=> G=E.

Run:  uv run --with sympy python research/dc1-program/verify_band45_lambda.py
Ends: ALL BAND45 LAMBDA CHECKS PASSED
"""
import sympy as sp

E, S, t = sp.symbols('E S t')


# ----------------------------------------------------------------- primitives
def sh(f, s):
    return sp.expand(sp.sympify(f).subs(E, E + s))


def poly(name, deg):
    if deg < 0:
        return sp.Integer(0), []
    cs = list(sp.symbols(f'{name}_0:{deg+1}'))
    return sp.expand(sum(cs[i] * E**i for i in range(deg + 1))), cs


def falling(k):
    return sp.prod([E - i for i in range(k)]) if k > 0 else sp.Integer(1)


def a_of_roots(roots):
    return sp.expand(sp.prod([E - r for r in roots]))


def az(expr, label):
    e = sp.expand(sp.sympify(expr))
    if e != 0:
        raise AssertionError(label + "  :  residual = " + str(e))
    print("PASS", label)


def ok(cond, label):
    if not cond:
        raise AssertionError(label + "  :  FALSE")
    print("PASS", label)


def Qm(A, B, m, K):
    return sp.expand(sum(sh(B[l], k) * A[k] - sh(A[k], l) * B[l]
                         for k in range(-K, K + 1) for l in range(-K, K + 1) if k + l == m))


def potential_G(A, B, K):
    return sp.expand(sum(sh(A[k], j - k) * sh(B[-k], j) - sh(B[k], j - k) * sh(A[-k], j)
                         for k in range(1, K + 1) for j in range(0, k)))


# --------------------------------------------------- necklace / cofactor / lambda
def cofactor(roots, K):
    """C(t) = A(t)/S_K(t);  A = sum t^rho,  S_K = 1+t+..+t^{K-1}."""
    A = sum(t**r for r in roots)
    SK = sum(t**j for j in range(K))
    C, rem = sp.div(sp.expand(A), sp.expand(SK), t)
    assert sp.expand(rem) == 0, f"top {roots}: not S_{K}-divisible"
    return sp.expand(C)


def bnecklace(roots, K):
    """B(t) = t S_{K-1} C ; roots of u=b_{K-1} = support of B (effective tops)."""
    C = cofactor(roots, K)
    Skm1 = sum(t**j for j in range(K - 1))
    B = sp.expand(t * Skm1 * C)
    Bp = sp.Poly(B, t)
    assert all(Bp.nth(i) >= 0 for i in range(Bp.degree() + 1)), "u necklace not effective"
    br = [i for i in range(Bp.degree() + 1) for _ in range(int(Bp.nth(i)))]
    return B, br


def lam_weights(roots, K):
    """hat_lambda(t) = t^K C(t) - C(1) ;  weights w_p = [t^p] hat_lambda."""
    C = cofactor(roots, K)
    hat = sp.expand(t**K * C - C.subs(t, 1))
    hp = sp.Poly(hat, t)
    return {i: int(hp.nth(i)) for i in range(hp.degree() + 1) if hp.nth(i) != 0}, hat


def lam(f, w):
    return sp.expand(sum(c * sp.sympify(f).subs(E, p) for p, c in w.items()))


def Kblock(aK, c, K):
    """K_K[c] = sum_{j=0}^{K-1} a_K(E+j-K) c(E+j)."""
    return sp.expand(sum(sh(aK, j - K) * sh(c, j) for j in range(K)))


def Hblock(bKm1, v, K):
    """H_{K-1}[v] = sum_{j=0}^{K-2} b_{K-1}(E+j-(K-1)) v(E+j)."""
    return sp.expand(sum(sh(bKm1, j - (K - 1)) * sh(v, j) for j in range(K - 1)))


def reversal(C):
    p = sp.Poly(C, t)
    d = p.degree()
    return sp.expand(sum(p.nth(i) * t**(d - i) for i in range(d + 1)))


def is_cyclotomic(pexpr):
    d = sp.Poly(sp.expand(pexpr), t).degree()
    for n in range(1, 60):
        cn = sp.cyclotomic_poly(n, t)
        if sp.Poly(cn, t).degree() == d and sp.expand(pexpr - cn) == 0:
            return True
    return False


# ---------------------------------------------------- exotic tops (from corpus)
TOPS = {
    4: [('{0,2,3,5}', [0, 2, 3, 5]),   # cofactor Phi_6, self-dual
        ('{0,1,3,6}', [0, 1, 3, 6]),   # reflection pair with {0,3,5,6}
        ('{0,3,5,6}', [0, 3, 5, 6]),   # reflection pair with {0,1,3,6}
        ('{0,3,6,9}', [0, 3, 6, 9])],  # step-3 AP, self-dual
    5: [('{0,2,3,4,6}', [0, 2, 3, 4, 6]),   # universal cofactor Phi_6 (g=1-s+s^2)
        ('{0,2,4,6,8}', [0, 2, 4, 6, 8]),   # step-2 AP, cofactor Phi_10
        ('{0,1,3,4,7}', [0, 1, 3, 4, 7])],  # NON-cyclotomic cofactor S^3-S^2+1
}


# ============================================================================
print("=" * 74)
print("0. CONVENTION ECHO:  Q_0=(T-1)G ;  G(0)=0 under membership  (bands 4,5)")
print("=" * 74)
for K in (4, 5):
    Ag = {k: poly(f'A{k+K}_', 2)[0] for k in range(-K, K + 1)}
    Bg = {k: poly(f'B{k+K}_', 2)[0] for k in range(-K, K + 1)}
    az(Qm(Ag, Bg, 0, K) - (sh(potential_G(Ag, Bg, K), 1) - potential_G(Ag, Bg, K)),
       f"k={K}: Q_0 = (T-1) G   (closure-form telescoping potential)")
    Am = {k: poly(f'Cm{k+K}_', 2)[0] for k in range(0, K + 1)}
    Bm = {k: poly(f'Dm{k+K}_', 2)[0] for k in range(0, K + 1)}
    for r in range(1, K + 1):
        Am[-r] = sp.expand(falling(r) * poly(f'amt{K}{r}', 2)[0])
        Bm[-r] = sp.expand(falling(r) * poly(f'bmt{K}{r}', 2)[0])
    az(potential_G(Am, Bm, K).subs(E, 0), f"k={K}: G(0)=0 under membership  =>  Q_0=1 <=> G=E")


# ============================================================================
print("\n" + "=" * 74)
print("1. THE DECOMPOSITION  G = R + K_K[b_-K] - H_{K-1}[a_-(K-1)]  (generic identity)")
print("=" * 74)
# R = sum_{i=1}^{K-1} K_i[b_-i] - sum_{i=1}^{K-2} H_i[a_-i], with each block using
#    the level-i coefficient.  Verified as a polynomial identity in generic
#    symbolic coefficients (gauge b_K=0), i.e. degree-free structural fact.
for K in (4, 5):
    A = {k: poly(f'ax{k+K}_', 2)[0] for k in range(-K, K + 1)}
    B = {k: poly(f'bx{k+K}_', 2)[0] for k in range(-K, K + 1)}
    B[K] = sp.Integer(0)                      # gauge

    def Ki(i, c):
        return sp.expand(sum(sh(A[i], j - i) * sh(c, j) for j in range(i)))

    def Hi(i, v):
        return sp.expand(sum(sh(B[i], j - i) * sh(v, j) for j in range(i)))

    KK = Ki(K, B[-K])                          # K_K[b_-K]
    HH = Hi(K - 1, A[-(K - 1)])                # H_{K-1}[a_-(K-1)]
    R = sp.expand(sum(Ki(i, B[-i]) for i in range(1, K))
                  - sum(Hi(i, A[-i]) for i in range(1, K - 1)))
    az(potential_G(A, B, K) - (R + KK - HH),
       f"k={K}: G == R + K_K[b_-K] - H_{K-1}[a_-(K-1)]  (R = sum K_i[b_-i] - sum H_i[a_-i])")


# ============================================================================
print("\n" + "=" * 74)
print("2. THE COFACTOR ANNIHILATOR  (DEGREE-FREE, uniform over every exotic top)")
print("=" * 74)
lamE_table = {}
for K, tops in TOPS.items():
    Skm1_t = sum(t**j for j in range(K - 1))
    SK_t = sum(t**j for j in range(K))
    for tag, roots in tops:
        aK = a_of_roots(roots)
        C = cofactor(roots, K)
        Aroot = sum(t**r for r in roots)
        Bnk, br = bnecklace(roots, K)
        u = a_of_roots(br)

        # 2a. cofactor facts
        az(C.subs(t, 1) - 1, f"k={K} {tag}: C(1)=1  (=> lambda(1)=0)")
        # wall solved by u=b_{K-1}
        wall = (sh(u, K) * aK - sh(aK, K - 1) * u)
        az(wall, f"k={K} {tag}: u=b_{K-1} (roots {br}) solves the wall b^[K]a - a^[K-1]b = 0")
        az(Bnk - sp.expand(t * Skm1_t * C), f"k={K} {tag}: necklace B = t*S_{K-1}*C  (wall root-necklace)")

        # 2b. the annihilator + the TWO generating-function support identities (THE proof)
        w, hat = lam_weights(roots, K)
        # hat * S_K = t^K A - C(1) S_K  (support {0..K-1} u {K+rho})
        az(sp.expand(hat * SK_t) - sp.expand(t**K * Aroot - C.subs(t, 1) * SK_t),
           f"k={K} {tag}: hat_lambda * S_K = t^K A - C(1) S_K  (K-block support identity)")
        # hat * S_{K-1} = t^{K-1} B - C(1) S_{K-1}  (support {0..K-2} u {K-1+sigma})
        az(sp.expand(hat * Skm1_t) - sp.expand(t**(K - 1) * Bnk - C.subs(t, 1) * Skm1_t),
           f"k={K} {tag}: hat_lambda * S_{K-1} = t^(K-1) B - C(1) S_{K-1}  (H-block support identity)")

        # 2c. resulting block annihilation as SYMBOLIC identities (arbitrary degree):
        Cpoly, cc = poly(f'Cc{K}', 4)
        Vpoly, vv = poly(f'Vv{K}', 4)
        cadm = sp.expand(falling(K) * Cpoly)
        vadm = sp.expand(falling(K - 1) * Vpoly)
        az(lam(Kblock(aK, cadm, K), w),
           f"k={K} {tag}: lambda(K_K[(E)_K C])=0 for GENERIC deg-4 C  (degree-free annihilation)")
        az(lam(Hblock(u, vadm, K), w),
           f"k={K} {tag}: lambda(H_{K-1}[(E)_{K-1} V])=0 for GENERIC deg-4 V  (degree-free annihilation)")

        # 2d. lambda(1)=0 ; lambda(E) closed form and positivity
        az(sum(c for c in w.values()), f"k={K} {tag}: lambda(1) = sum w_p = 0")
        lamE_direct = sum(c * p for p, c in w.items())
        lamE_form = sp.Rational(K + 1, 2) + sp.Rational(sum(roots), K)
        az(lamE_direct - lamE_form,
           f"k={K} {tag}: lambda(E)=sum w_p*p = (K+1)/2 + (sum rho)/K = {lamE_form}")
        ok(lamE_form > 0, f"k={K} {tag}: lambda(E) = {lamE_form} > 0  (=> E not in Im Phi, arbitrary degree)")
        lamE_table[(K, tag)] = lamE_form

        # 2e. ev_0 is the always-present (non-obstructing) annihilator
        az(Kblock(aK, cadm, K).subs(E, 0), f"k={K} {tag}: ev_0(K_K[c])=0 (membership)")
        az(Hblock(u, vadm, K).subs(E, 0), f"k={K} {tag}: ev_0(H_{K-1}[v])=0 (membership); ev_0(E)=0 non-obstructing")


# ============================================================================
print("\n" + "=" * 74)
print("3. REFLECTION  E -> -E-1 :  C*(refl) = reversal(C) ; self-dual => palindromic")
print("=" * 74)


def reflect(roots):
    M = max(roots)
    return sorted(M - x for x in roots)


for K, tops in TOPS.items():
    for tag, roots in tops:
        rr = reflect(roots)
        C = cofactor(roots, K)
        Cr = cofactor(rr, K)
        az(Cr - reversal(C), f"k={K} {tag} -> {rr}: cofactor of reflected top = reversal(C)")
        selfdual = (rr == sorted(roots))
        if selfdual:
            az(C - reversal(C), f"k={K} {tag}: SELF-DUAL => C palindromic => lambda reflection-symmetric")
        # lambda(E) transforms as (K+1)/2 + M - (sum rho)/K
        lamE_r = sp.Rational(K + 1, 2) + sp.Rational(sum(rr), K)
        az(lamE_r - (sp.Rational(K + 1, 2) + max(roots) - sp.Rational(sum(roots), K)),
           f"k={K} {tag}: lambda_refl(E) = (K+1)/2 + max(rho) - (sum rho)/K")
ok(sorted([tuple(reflect([0, 1, 3, 6])), tuple([0, 3, 5, 6])])[0] == (0, 3, 5, 6),
   "band-4 reflection PAIR: {0,1,3,6} <-> {0,3,5,6} (halves the work)")


# ============================================================================
print("\n" + "=" * 74)
print("4. THE NON-CYCLOTOMIC TEST:  mechanism does NOT depend on cyclotomic structure")
print("=" * 74)
Cnc = cofactor([0, 1, 3, 4, 7], 5)
ok(sp.expand(Cnc - (t**3 - t**2 + 1)) == 0, "k=5 {0,1,3,4,7}: cofactor C = S^3 - S^2 + 1")
ok(not is_cyclotomic(Cnc), "k=5 {0,1,3,4,7}: C = S^3-S^2+1 is NON-cyclotomic (not any Phi_n)")
# but its reflection cofactor is also non-cyclotomic and the annihilator works identically:
Cncr = cofactor(reflect([0, 1, 3, 4, 7]), 5)
ok(not is_cyclotomic(Cncr), "k=5 {0,1,3,4,7} reflects to another NON-cyclotomic top (S^3-S+1)")
print("   => the sharpest test passes: the annihilator template is cyclotomic-independent.")


# ============================================================================
print("\n" + "=" * 74)
print("5. BOUNDED RESIDUAL KILL:  lambda(G)=lambda(R) --> CONST 0 mod positive cascade")
print("   (per-degree Groebner; NOT an arbitrary-degree theorem -- see band45-lambda.md)")
print("=" * 74)


def clean_solve(A, B, m, lkey, name, memb, bdeg_raw, K):
    braw, bc = poly(f'{name}c_', bdeg_raw)
    bnew = sp.expand(falling(memb) * braw)
    Bt = dict(B)
    Bt[lkey] = bnew
    eq = sp.expand(Qm(A, Bt, m, K))
    coeffs = sp.Poly(eq, E).all_coeffs() if eq != 0 else [sp.Integer(0)]
    M, rhs = sp.linear_eq_to_matrix(coeffs, bc)
    conds = [sp.expand(n.dot(rhs)) for n in M.T.nullspace()]
    conds = [c for c in conds if c != 0]
    Mred = sp.zeros(0, len(bc))
    rhs_red = []
    for i in range(M.shape[0]):
        cand = Mred.col_join(M[i, :])
        if cand.rank() > Mred.rank():
            Mred = cand
            rhs_red.append(rhs[i])
    if Mred.shape[0] == 0:
        bc_vals = [sp.Integer(0)] * len(bc)
    else:
        sol, _ = Mred.gauss_jordan_solve(sp.Matrix(rhs_red))
        tau = [s for s in sol.free_symbols if str(s).startswith('tau')]
        bc_vals = [x.subs({tt: 0 for tt in tau}) for x in sol]
    bsol = sp.expand(bnew.subs(dict(zip(bc, bc_vals))))
    ker = []
    for idx, kv in enumerate(M.nullspace()):
        g = sp.symbols(f'{name}K{idx}')
        kp = sp.expand(falling(memb) * sum(kv[i] * E**i for i in range(len(bc))))
        bsol = sp.expand(bsol + g * kp)
        ker.append(g)
    return bsol, ker, conds


def build_cascade(aK, uK, d, K):
    """Forward-solve the positive cascade for band K; returns (A,B,pos,allvars)."""
    A = {K: aK}
    allc = []
    for i in range(K - 1, -1, -1):
        p, c = poly(f'a{i}_', d)
        A[i] = p
        allc += c
    cneg = {}
    for r in range(1, K + 1):
        raw, cs = poly(f'am{r}_', d)
        A[-r] = sp.expand(falling(r) * raw)
        cneg[-r] = cs
        allc += cs
    mu = sp.symbols(f'mu{K}')
    B = {k: sp.Integer(0) for k in range(-K, K + 1)}
    B[K - 1] = uK
    B[-K] = sp.expand(mu * A[-K])
    # solve Q_{2K-2}, ..., Q_1 for b_{K-2},...,b_{-(K-1)}
    plan = []
    for l in range(K - 2, -(K), -1):
        m = l + K            # pair (K, l): Q_{K+l}
        mem = max(0, -l)
        plan.append((m, l, mem))
    pos = []
    ks = []
    cap = (2 * d + K + 2)
    for (m, l, mem) in plan:
        bb, kk, cc = clean_solve(A, B, m, l, f'b{l}_', mem, cap, K)
        B[l] = bb
        pos += cc
        ks += kk
    allvars = allc + [mu] + ks
    return A, B, pos, allvars


DEG = {4: [1, 2], 5: [1, 2]}     # d=2 only for a representative subset (below)
D2_SUBSET = {'{0,2,3,5}', '{0,3,6,9}', '{0,2,3,4,6}'}
for K, tops in TOPS.items():
    for tag, roots in tops:
        aK = a_of_roots(roots)
        _, br = bnecklace(roots, K)
        uK = a_of_roots(br)
        w, _ = lam_weights(roots, K)
        lamE = sum(c * p for p, c in w.items())
        degs = [1] + ([2] if tag in D2_SUBSET else [])
        for d in degs:
            A, B, pos, allvars = build_cascade(aK, uK, d, K)
            G = potential_G(A, B, K)
            lamG = lam(G, w)
            gb = sp.groebner(pos, *allvars, order='grevlex')
            red = sp.expand(gb.reduce(sp.expand(lamG))[1])
            ok(red == 0,
               f"k={K} {tag} d={d}: lambda(G)=lambda(R) reduces to CONST 0 mod positive cascade "
               f"(lambda(E)={lamE} != 0  =>  lambda(E-R)={lamE}  =>  E-R not in Im Phi: Q_0=1 killed)")
            # cascade must itself be feasible (not the unit ideal) -- no false kill
            ok(list(gb) != [sp.Integer(1)],
               f"k={K} {tag} d={d}: positive cascade is FEASIBLE (proper ideal) -- kill is real, not vacuous")


print("\n" + "=" * 74)
print("SUMMARY  lambda(E) = (k+1)/2 + (sum rho)/k  per top:")
for (K, tag), v in lamE_table.items():
    print(f"   k={K:>1}  {tag:<11}  lambda(E) = {v}")
print("=" * 74)
print("\nALL BAND45 LAMBDA CHECKS PASSED")
