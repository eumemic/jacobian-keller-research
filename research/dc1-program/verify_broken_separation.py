#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# verify_broken_separation.py
#
# Exact-SymPy certificate for research/dc1-program/broken-separation.md.
#
# THE BROKEN SEPARATION CLASSES.  The band-3 shifted-cube wall (gauge b_3=0,
# a_3 = h h^[1] h^[2], b_2 = kappa h h^[1]) is closed at ARBITRARY degree for
# 2-separated h (gcd(h,h^[1])=gcd(h,h^[2])=1), including the diff-3 class -- see
# shifted-power-residuals.md.  The last shifted-power walls are the two classes
# where the Q_4-rung coprimality genuinely BREAKS:
#
#     diff-1:  h = (E-r)(E-r-1),   gcd(h,h^[1]) = (E-r) != 1,   r free,
#     diff-2:  h = (E-r)(E-r-2),   gcd(h,h^[2]) = (E-r) != 1,   r free,
#
# plus the multiple-root h = (E-r)^2 (checked below: it is CUBE-SEPARATED, hence
# already closed by the arbitrary-degree theorem -- NOT a separate broken class).
#
# WHAT THIS FILE ESTABLISHES (all re-derived in-file from the ladder engine):
#
#  TASK 1 -- the exact DEGRADED Q_4 forcing, symbolic in r:
#    * diff-1:  Q_4=0 forces  (E-r) | a_2  AND the coupled correction
#               a_2(r-1) + a_2(r+1) = 0.  The clean  h h^[1] | a_2  FAILS.
#    * diff-2:  Q_4=0 forces the clean PROPER-FACTOR divisibility
#               (E-r)(E-r-1) | a_2.  The clean  h h^[1] | a_2  FAILS.
#    * mult-root (E-r)^2: cube-separated => clean Q_4 forces h h^[1] | a_2 and
#      the affine kill runs => CLOSED arbitrary degree (covered, not broken).
#
#  TASK 2 -- the cascade Q_3, Q_2 with degraded shapes:
#    * the general Q_3 identity (NO shape assumption) is re-derived;
#    * a_3's root multiset acquires DOUBLE nodes (diff-1: at r-1 and r; diff-2:
#      at r) -- the source of the derivative-node criterion of Task 3;
#    * a GENERIC Q_4 solution does NOT extend to Q_3 (linear inconsistency over
#      the function field -- Q_3 constrains the a_2 freedom FURTHER);
#    * tested Gröbner normal forms show ideal membership at the gcd node and
#      ideal nonmembership at a clean-extra node. This does NOT decide radical
#      membership, geometric restoration, or alternative proper-factor mechanisms.
#
#  TASK 3 -- the multiplicity-extended moment/adjoint criterion:
#    * the jet divisibility criterion a_3 | P  <=>  P^(j)(nu)=0 for j < m_nu;
#    * the exact derivative-node (ev_rho') equations at the double nodes;
#    * the Lemma-P moment slope  G(1) = sum_i (a_i(0)b_-i(i) - a_-i(i)b_i(0))
#      re-derived for the DEGENERATE tops (G(0)=0 by membership).
#
#  TASK 4 -- bounded emptiness certificates per class (positive cascade
#            Q_4=Q_3=Q_2=Q_1=0 + Q_0=1 + membership):
#    * committed (exact SymPy; displayed coefficient domain ZZ[r,kappa], with r,kappa
#      treated as coefficients): cap d=2 GENERIC-FIBER UNIT IDEAL over QQ(r,kappa),
#      BOTH classes -- not a uniform certificate on every parameter specialization.
#    * HEAVY (msolve, exact QQ, SKIP if absent): a grid of sampled (r,kappa)
#      covering the special r-loci at cap d=2 AND d=3 (and d=4), BOTH classes.
#
# Conventions frozen from the DC1 corpus:
#   A_1[x^-1] = (+)_k x^k C[E],  E = x d,  (x^a f)(x^b g) = x^{a+b} f(E+b) g(E),
#   f^[n](E) = f(E+n),  T f = f^[1],  S_n = 1+T+...+T^{n-1},
#   Q_m = sum_{k+l=m} ( b_l^[k] a_k - a_k^[l] b_l ),  [D,X]=1 <=> Q_m = delta_{m0},
#   membership (E)_j = E(E-1)...(E-j+1) | a_-j, b_-j,  gauge b_k = 0,
#   G = sum_{k>=1} sum_{j=0}^{k-1} ( a_k^[j-k] b_-k^[j] - b_k^[j-k] a_-k^[j] ),
#   Q_0 = (T-1)G.
#
# EVIDENCE TIERS:
#   * committed = reproduced by the DEFAULT run (exact SymPy over QQ; symbolic
#     parameter calculation displayed over ZZ[r,kappa] and interpreted on the
#     generic fiber QQ(r,kappa)): the degraded Q_4 forcing, the Q_3 cascade
#     structure, the multiplicity criterion + moment slope, the mult-root closure,
#     and cap-d=2 generic-fiber emptiness for both classes.
#   * HEAVY (set HEAVY=1) = msolve corroboration on a sampled (r,kappa) grid at
#     cap d=2,3,4.  If msolve is absent or HEAVY is unset these print SKIP and the
#     final status line does NOT claim them.
#
# Run:   uv run --with sympy python research/dc1-program/verify_broken_separation.py
#        HEAVY=1 uv run --with sympy python research/dc1-program/verify_broken_separation.py
# The final summary distinguishes no-skip success from executed-check success with skips.
# ---------------------------------------------------------------------------
import sympy as sp
import os
import time
import shutil
import subprocess
import tempfile

E = sp.symbols('E')
kappa, r = sp.symbols('kappa r')

HEAVY = os.environ.get('HEAVY', '') == '1'
HAVE_MSOLVE = shutil.which('msolve') is not None

CHECKS = []
SKIPS = []


def check(name, cond):
    ok = bool(cond)
    CHECKS.append((name, ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}", flush=True)
    return ok


def skip(name, why):
    SKIPS.append((name, why))
    print(f"[SKIP] {name}  ({why})", flush=True)


def sh(f, n):                                   # f^[n](E) = f(E+n)
    return sp.expand(sp.sympify(f).subs(E, E + n))


def gen(name, d):                               # generic degree-d poly + coeff list
    cs = sp.symbols(f'{name}_0:{d+1}')
    if not isinstance(cs, (list, tuple)):
        cs = (cs,)
    return sum(cs[i] * E**i for i in range(d + 1)), list(cs)


def divides(fac, c):                            # fac | c over C[E]
    c = sp.expand(sp.sympify(c))
    if c == 0:
        return True
    if sp.expand(fac) == 0:
        return False
    return sp.rem(sp.Poly(c, E), sp.Poly(sp.expand(fac), E)) == 0


def prodsh(h, m):                               # prod_{j=0}^{m-1} h^[j]
    return sp.expand(sp.prod([sh(h, j) for j in range(m)])) if m > 0 else sp.Integer(1)


def gcd_deg(a, b):
    return sp.gcd(sp.Poly(sp.expand(a), E), sp.Poly(sp.expand(b), E)).degree()


# ---- crossed-product ladder engine (dicts index -> C[E]) ----
def cp_mul(P, Q):
    R = {}
    for k, fk in P.items():
        for l, gl in Q.items():
            R[k + l] = sp.expand(R.get(k + l, 0) + sh(fk, l) * gl)
    return {m: v for m, v in R.items() if sp.expand(v) != 0}


def cp_sub(P, Q):
    R = dict(P)
    for k, v in Q.items():
        R[k] = sp.expand(R.get(k, 0) - v)
    return {m: v for m, v in R.items() if sp.expand(v) != 0}


def commutator(D, X):                           # [D,X] = DX - XD
    return cp_sub(cp_mul(D, X), cp_mul(X, D))


def Qm(X, D, m):                                # ladder-m coeff of [D,X]
    tot = 0
    for k, ak in X.items():
        l = m - k
        if l in D:
            tot += sh(D[l], k) * ak - sh(ak, l) * D[l]
    return sp.expand(tot)


def Gpot(X, D, K=3):                            # closed-form moment potential
    G = 0
    for k in range(1, K + 1):
        for j in range(0, k):
            G += sh(X.get(k, 0), j - k) * sh(D.get(-k, 0), j) \
               - sh(D.get(k, 0), j - k) * sh(X.get(-k, 0), j)
    return sp.expand(G)


def msolve_unit(eqs, allc, tag, timeout=300):
    """True iff ideal(eqs) is the unit ideal over QQ (reduced GB = [1])."""
    varstr = ", ".join(str(v) for v in allc)
    # msolve's characteristic-zero parser is most reliable on integral
    # polynomials.  Clear numeric denominators (notably the r=1/2 grid point)
    # equation-by-equation; multiplying by a nonzero rational preserves the ideal.
    cleared = [sp.fraction(sp.together(sp.expand(e)))[0] for e in eqs]
    body = ",\n".join(str(e).replace('**', '^').replace(' ', '') for e in cleared)
    ms = f"{varstr}\n0\n{body}\n"
    if '**' in ms:
        raise ValueError("msolve input must use '^' not '**'")
    t = time.time()
    with tempfile.TemporaryDirectory(prefix="broken-separation-") as tmp:
        path = os.path.join(tmp, f"{tag}.ms")
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(ms)
        try:
            rr = subprocess.run(['msolve', '-g', '2', '-f', path],
                                capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return None
    if rr.returncode != 0:
        raise RuntimeError(f"msolve[{tag}] failed with status {rr.returncode}: {rr.stderr.strip()}")
    out = rr.stdout.strip()
    if not out:
        raise RuntimeError(f"msolve[{tag}] produced empty output; stderr: {rr.stderr.strip()}")
    dt = time.time() - t
    tail = out.splitlines()[-1]
    if not tail.rstrip().endswith(']:'):
        raise RuntimeError(f"msolve[{tag}] malformed output tail: {tail[:200]!r}")
    unit = tail.rstrip().endswith('[1]:')
    print(f"        msolve[{tag}]: nvars={len(allc)} neqs={len(eqs)} "
          f"time={dt:.1f}s GB_tail={tail[:24]} UNIT={unit}", flush=True)
    return unit


# sector shapes shared by several sections -----------------------------------
def sector_top(hc, ksym=kappa):
    """a_3 = h h^[1] h^[2],  b_2 = kappa h h^[1] (the Q_5 wall)."""
    hc = sp.expand(hc)
    h1, h2 = sh(hc, 1), sh(hc, 2)
    return sp.expand(hc * h1 * h2), sp.expand(ksym * hc * h1)


# ===========================================================================
print("\n=== §0  crossed-product engine (re-derived, self-contained) ===")
# ===========================================================================
Xg = {k: gen(f'a{k}', 2)[0] for k in range(-3, 4)}
Dg = {k: gen(f'b{k}', 2)[0] for k in range(-3, 4)}
Cg = commutator(Dg, Xg)
check("Q_m == [D,X]_m for all m in [-6,6] (generic band-3 data)",
      all(sp.expand(Cg.get(m, 0) - Qm(Xg, Dg, m)) == 0 for m in range(-6, 7)))
check("Q_0 = (T-1)G with the closed-form staggered potential (generic band-3)",
      sp.expand(sh(Gpot(Xg, Dg), 1) - Gpot(Xg, Dg) - Qm(Xg, Dg, 0)) == 0)
check("[d, x] = 1 (reference pair (x,d))",
      commutator({-1: E}, {1: sp.Integer(1)}) == {0: sp.Integer(1)})


# ===========================================================================
print("\n=== §1  the shifted-cube wall + the general Q_4, Q_3 identities ===")
# ===========================================================================
# The Q_5 wall b_2^[3] a_3 = a_3^[2] b_2 is solved by b_2 = kappa h h^[1] for any h
# (re-derived here at symbolic h; the 1-dimensionality is in the sibling memos).
for hc, tag in [(E - r, "deg1"), ((E - r) * (E - r - 1), "diff-1"),
                ((E - r) * (E - r - 2), "diff-2"), ((E - r)**2, "mult-root")]:
    a3, b2 = sector_top(hc)
    check(f"Q_5 wall: b_2 = kappa h h^[1] solves b_2^[3]a_3 = a_3^[2]b_2  [{tag}]",
          sp.expand(sh(b2, 3) * a3 - sh(a3, 2) * b2) == 0)

# The GENERAL Q_4 identity (gauge b_3=0, a_3=h h^[1]h^[2], b_2=kappa h h^[1]),
# NO shape assumption on a_2,b_1:
#   Q_4 = h^[1]h^[2]( b_1^[3] h - h^[3] b_1 ) + kappa( h^[2]h^[3] a_2 - h h^[1] a_2^[2] ).
hh = sp.expand((E - r) * (E - r - 1))
h1, h2, h3 = sh(hh, 1), sh(hh, 2), sh(hh, 3)
a3, b2 = sector_top(hh)
a2G, _ = gen('a2', 3); b1G, _ = gen('b1', 3)
X4 = {3: a3, 2: a2G, 1: sp.Integer(0)}
D4 = {3: sp.Integer(0), 2: b2, 1: b1G}
Q4_id = sp.expand(h1 * h2 * (sh(b1G, 3) * hh - h3 * b1G)
                  + kappa * (h2 * h3 * a2G - hh * h1 * sh(a2G, 2)))
check("general Q_4 = h^[1]h^[2](b_1^[3]h - h^[3]b_1) + kappa(h^[2]h^[3]a_2 - h h^[1]a_2^[2])",
      sp.expand(Qm(X4, D4, 4) - Q4_id) == 0)

# The GENERAL Q_3 identity (NO shape assumption on a_2,b_1,a_1):
#   Q_3 = a_3(b_0^[3]-b_0) + (b_1^[2]a_2 - a_2^[1]b_1) + kappa h^[1](h^[2]a_1 - h a_1^[2]).
a1G, _ = gen('a1', 3); b0G, _ = gen('b0', 3); a0G, _ = gen('a0', 3)
X3 = {3: a3, 2: a2G, 1: a1G, 0: a0G}
D3 = {3: sp.Integer(0), 2: b2, 1: b1G, 0: b0G}
Q3_id = sp.expand(a3 * (sh(b0G, 3) - b0G) + (sh(b1G, 2) * a2G - sh(a2G, 1) * b1G)
                  + kappa * h1 * (h2 * a1G - hh * sh(a1G, 2)))
check("general Q_3 = a_3(b_0^[3]-b_0) + (b_1^[2]a_2 - a_2^[1]b_1) + kappa h^[1](h^[2]a_1 - h a_1^[2])",
      sp.expand(Qm(X3, D3, 3) - Q3_id) == 0)


# ===========================================================================
print("\n=== §2  TASK 1: the exact DEGRADED Q_4 forcing per class (symbolic r) ===")
# ===========================================================================
# Q_4=0 is LINEAR in the joint coefficient vector (a_2, b_1) (the top b_2 a_2 term
# is linear in a_2, the b_1 terms linear in b_1).  So the "generic Q_4 solution"
# is the FULL solution space of a linear system -- reliable, no component dropping.
def q4_forcing(hc, dcap=4):
    hc = sp.expand(hc)
    a3, b2 = sector_top(hc)
    a2, a2c = gen('a2', dcap); b1, b1c = gen('b1', dcap)
    X = {3: a3, 2: a2, 1: sp.Integer(0)}; D = {3: sp.Integer(0), 2: b2, 1: b1}
    sols = sp.solve(sp.Poly(Qm(X, D, 4), E).all_coeffs(), b1c + a2c, dict=True)
    if len(sols) != 1:
        raise RuntimeError(f"Q_4 not a single linear component: {len(sols)}")
    a2s = sp.expand(a2.subs(sols[0]))
    return a2s

# diff-1: (E-r) | a_2  AND  a_2(r-1) + a_2(r+1) = 0;  clean h h^[1] | a_2 FAILS.
hd1 = sp.expand((E - r) * (E - r - 1))
a2d1 = q4_forcing(hd1)
check("DIFF-1 Q_4 forcing: (E-r) | a_2 (the gcd factor divides)",
      divides(E - r, a2d1))
check("DIFF-1 Q_4 forcing: coupled correction a_2(r-1) + a_2(r+1) = 0 (EXACT, symbolic r)",
      sp.expand(a2d1.subs(E, r - 1) + a2d1.subs(E, r + 1)) == 0)
check("DIFF-1 Q_4: clean h | a_2 FAILS and h h^[1] | a_2 FAILS (genuinely degraded)",
      (not divides(hd1, a2d1)) and (not divides(sp.expand(hd1 * sh(hd1, 1)), a2d1)))
# and a_2'(r) is NOT forced to 0 (the double-root condition of the clean case is lost)
check("DIFF-1 Q_4: a_2'(r) is FREE (clean double-root condition a_2'(r)=0 is lost)",
      sp.expand(sp.diff(a2d1, E).subs(E, r)) != 0)

# diff-2: clean PROPER-FACTOR divisibility (E-r)(E-r-1) | a_2;  h h^[1] | a_2 FAILS.
hd2 = sp.expand((E - r) * (E - r - 2))
a2d2 = q4_forcing(hd2)
check("DIFF-2 Q_4 forcing: (E-r)(E-r-1) | a_2 (clean proper factor: gcd node + midpoint)",
      divides(sp.expand((E - r) * (E - r - 1)), a2d2))
check("DIFF-2 Q_4: the divisor (E-r)(E-r-1) is NOT h and NOT a factor of h "
      "(h roots r,r+2; forced roots r,r+1 -- no coupling term)",
      (not divides(hd2, a2d2)) and gcd_deg(hd2, sp.expand((E - r) * (E - r - 1))) == 1)
check("DIFF-2 Q_4: clean h h^[1] | a_2 FAILS (loses vanishing at r-1 and r+2)",
      not divides(sp.expand(hd2 * sh(hd2, 1)), a2d2))

# mult-root h=(E-r)^2: it is CUBE-SEPARATED, so the clean theorem applies:
hmr = sp.expand((E - r)**2)
check("MULT-ROOT (E-r)^2 is CUBE-SEPARATED: gcd(h,h^[j])=1 for j=1,2,3 "
      "(zero root-gap is not in {1,2,3}) -- NOT a broken class",
      all(gcd_deg(hmr, sh(hmr, j)) == 0 for j in (1, 2, 3)))
a2mr = q4_forcing(hmr)
check("MULT-ROOT Q_4 forces the CLEAN h h^[1] | a_2 (cube-separated, no slack)",
      divides(sp.expand(hmr * sh(hmr, 1)), a2mr))
check("MULT-ROOT affine kill: h^[-1] (deg 2) does NOT divide E => no nonconstant top "
      "(closed at arbitrary degree by the shifted-cube theorem)",
      not divides(sh(hmr, -1), E))


# ===========================================================================
print("\n=== §3  TASK 2: the cascade Q_3, Q_2 with the degraded shapes ===")
# ===========================================================================
# (3a) a_3's root MULTISET acquires DOUBLE nodes in the broken classes -- the
#      source of the derivative-node criterion (Task 3).  For h with a repeated
#      root-difference in {1,2}, two of the three shifted copies share a root.
def a3_multiset(hc):
    a3 = prodsh(sp.expand(hc), 3)
    return sp.roots(sp.Poly(a3, E))     # {root: multiplicity}
mult_d1 = a3_multiset((E - r) * (E - r - 1))
mult_d2 = a3_multiset((E - r) * (E - r - 2))
mult_mr = a3_multiset((E - r)**2)
check("DIFF-1 a_3 root multiset has DOUBLE nodes at r-1 and r "
      "(a_3 = (E-r+2..)... with mult-2 at r-1,r)",
      mult_d1.get(r - 1) == 2 and mult_d1.get(r) == 2 and max(mult_d1.values()) == 2)
check("DIFF-2 a_3 root multiset has a DOUBLE node at r (mult-2 at r; rest simple)",
      mult_d2.get(r) == 2 and sorted(mult_d2.values()) == [1, 1, 1, 1, 2])
check("MULT-ROOT a_3 = (h h^[1]h^[2]) is a perfect square (all nodes double) "
      "but cube-separated => clean theorem still applies",
      set(mult_mr.values()) == {2})

# (3b) a GENERIC Q_4 solution does NOT extend to Q_3.  Two-stage over the field:
#      Q_4 linear -> (a2g,b1g) with free theta; substitute; Q_3=0 is then linear
#      in (a1,b0) over QQ(theta,r,kappa) and is INCONSISTENT (no polynomial b_0,a_1).
def generic_no_q3_extension(hc):
    hc = sp.expand(hc); a3, b2 = sector_top(hc)
    a2, a2c = gen('a2', 5); b1, b1c = gen('b1', 5)
    s1 = sp.solve(sp.Poly(Qm({3: a3, 2: a2, 1: sp.Integer(0)},
                             {3: sp.Integer(0), 2: b2, 1: b1}, 4), E).all_coeffs(),
                  b1c + a2c, dict=True)[0]
    a2g = sp.expand(a2.subs(s1)); b1g = sp.expand(b1.subs(s1))
    a1, a1c = gen('a1', 6); b0, b0c = gen('b0', 6)
    q3 = Qm({3: a3, 2: a2g, 1: a1, 0: sp.Integer(0)},
            {3: sp.Integer(0), 2: b2, 1: b1g, 0: b0}, 3)
    s2 = sp.solve(sp.Poly(q3, E).all_coeffs(), a1c + b0c, dict=True)
    return len(s2) == 0
check("DIFF-1: a GENERIC Q_4 solution does NOT extend to Q_3 (Q_3 constrains a_2 further; symbolic r,kappa)",
      generic_no_q3_extension((E - r) * (E - r - 1)))
check("DIFF-2: a GENERIC Q_4 solution does NOT extend to Q_3 (Q_3 constrains a_2 further; symbolic r,kappa)",
      generic_no_q3_extension((E - r) * (E - r - 2)))

# (3c) At tested specializations, Gröbner reduction proves ideal membership for
#      the gcd-node evaluation and ideal nonmembership for the clean-extra
#      evaluation. Ideal nonmembership does not imply radical nonmembership or
#      geometric non-forcing, so restoration remains unresolved.
def restoration_status(hc, r0, k0):
    hc = sp.expand(hc); a3, b2 = sector_top(hc, sp.Integer(k0))
    a3 = sp.expand(a3.subs(r, r0)); b2 = sp.expand(b2.subs(r, r0))
    a2, a2c = gen('a2', 4); b1, b1c = gen('b1', 4); a1, a1c = gen('a1', 5); b0, b0c = gen('b0', 5)
    X = {3: a3, 2: a2, 1: a1, 0: sp.Integer(0)}; D = {3: sp.Integer(0), 2: b2, 1: b1, 0: b0}
    eqs = sp.Poly(Qm(X, D, 4), E).all_coeffs() + sp.Poly(Qm(X, D, 3), E).all_coeffs()
    G = sp.groebner([sp.expand(e) for e in eqs if sp.expand(e) != 0], *(a2c + b1c + a1c + b0c),
                    order='grevlex')
    gcd_forced = G.reduce(sp.expand(a2.subs(E, r0)))[1] == 0            # a_2(r)=0
    extra_forced = G.reduce(sp.expand(a2.subs(E, r0 - 1)))[1] == 0      # a_2(r-1)=0 (clean-extra)
    return gcd_forced, extra_forced
for hc, tag in [((E - r) * (E - r - 1), "diff-1"), ((E - r) * (E - r - 2), "diff-2")]:
    for r0, k0 in [(2, 1), (0, 2)]:
        gf, ef = restoration_status(hc, r0, k0)
        check(f"{tag} IDEAL TEST [r={r0},kappa={k0}]: gcd-node evaluation is in the ideal={gf}; "
              f"clean-extra evaluation is not in the ideal={not ef}", gf and (not ef))
print("   => ideal nonmembership is inconclusive about radical membership and geometric restoration;")
print("      restoration and alternative proper-factor mechanisms remain open.")


# ===========================================================================
print("\n=== §4  TASK 3: the multiplicity-extended moment/adjoint criterion ===")
# ===========================================================================
# (4a) the jet divisibility criterion, extending ev_rho to the jet (ev_rho, ev_rho', ...):
#      a_3 | P  <=>  P^(j)(nu) = 0 for all j < m_nu, at every node nu of a_3.
def jet_divides(a3, P):
    a3p = sp.Poly(sp.expand(a3), E); P = sp.expand(P)
    for nu, m in sp.roots(a3p).items():
        for j in range(m):
            if sp.expand(sp.diff(P, E, j).subs(E, nu)) != 0:
                return False
    return True
# agree with polynomial remainder on random data, INCLUDING at a double node:
a3t = sp.expand(E * (E - 1)**2 * (E - 3))
for Pt in [sp.expand(E * (E - 1)**2 * (E - 3) * (2 * E + 5)),   # divisible
           sp.expand(E * (E - 1) * (E - 3) * (2 * E + 5))]:      # NOT (only order 1 at 1)
    check(f"jet criterion == polynomial-remainder divisibility (double node at 1) [{Pt!=a3t*0}]",
          jet_divides(a3t, Pt) == divides(a3t, Pt))

# (4b) the EXACT derivative-node (ev_rho') equations Q_3=0 imposes at the double
#      nodes.  Q_3 = a_3(b_0^[3]-b_0) + Src_3, Src_3 = (b_1^[2]a_2 - a_2^[1]b_1)
#      + kappa h^[1](h^[2]a_1 - h a_1^[2]); a_3 | Src_3 requires, at a double node
#      nu, BOTH Src_3(nu)=0 and Src_3'(nu)=0 (the derivative-node equation).
def src3(hc, a2, b1, a1):
    hc = sp.expand(hc); h1, h2 = sh(hc, 1), sh(hc, 2)
    return sp.expand((sh(b1, 2) * a2 - sh(a2, 1) * b1) + kappa * h1 * (h2 * a1 - hc * sh(a1, 2)))
# symbolic check: with GENERIC a_2,b_1,a_1 the Src_3 is a genuine polynomial and the
# a_3-divisibility is exactly the jet conditions at a_3's nodes (double at r-1,r for diff-1).
hcs = sp.expand((E - r) * (E - r - 1)); a3s = prodsh(hcs, 3)
a2s, _ = gen('za2', 3); b1s, _ = gen('zb1', 3); a1s, _ = gen('za1', 3)
S = src3(hcs, a2s, b1s, a1s)
# the derivative-node equation at nu=r is ev_r'(Src_3)=0 -- a genuine, nonzero linear functional
ev_r = sp.expand(S.subs(E, r))
ev_r_prime = sp.expand(sp.diff(S, E).subs(E, r))
check("DIFF-1 Q_3 requires a derivative-node equation ev_r'(Src_3)=0 at the DOUBLE node r "
      "(a genuine extra linear functional, not present at simple nodes)",
      ev_r_prime != 0 and ev_r != 0 and ev_r_prime != ev_r)
# and the ordinary node ev_{r-2}(Src_3)=0 at the simple node r-2 (a3 has simple root r-2):
ev_rm2 = sp.expand(S.subs(E, r - 2))
check("DIFF-1 Q_3 also requires the ordinary node ev_{r-2}(Src_3)=0 at the SIMPLE node r-2",
      ev_rm2 != 0)

# (4c) the Lemma-P moment slope, re-derived for the DEGENERATE tops.  Membership
#      gives G(0)=0, and Q_0=1 forces G=E, so the slope G(1) must equal 1:
#         G(1) = sum_{i=1}^3 ( a_i(0) b_-i(i) - a_-i(i) b_i(0) ).
def slope_check(hc):
    hc = sp.expand(hc); a3, b2 = sector_top(hc)
    a2, _ = gen('sa2', 2); a1, _ = gen('sa1', 2); a0, _ = gen('sa0', 2)
    b1, _ = gen('sb1', 2); b0, _ = gen('sb0', 2)
    am1 = sp.expand(E * gen('sam1', 2)[0]); am2 = sp.expand(E * (E - 1) * gen('sam2', 2)[0])
    am3 = sp.expand(E * (E - 1) * (E - 2) * gen('sam3', 2)[0])
    bm1 = sp.expand(E * gen('sbm1', 2)[0]); bm2 = sp.expand(E * (E - 1) * gen('sbm2', 2)[0])
    bm3 = sp.expand(E * (E - 1) * (E - 2) * gen('sbm3', 2)[0])
    A = {3: a3, 2: a2, 1: a1, 0: a0, -1: am1, -2: am2, -3: am3}
    B = {3: sp.Integer(0), 2: b2, 1: b1, 0: b0, -1: bm1, -2: bm2, -3: bm3}
    G = Gpot(A, B)
    slope = sum(A[i].subs(E, 0) * B[-i].subs(E, i) - A[-i].subs(E, i) * B[i].subs(E, 0)
                for i in (1, 2, 3))
    return sp.expand(G.subs(E, 0)) == 0, sp.expand(G.subs(E, 1) - slope) == 0
for hc, tag in [((E - r) * (E - r - 1), "diff-1"), ((E - r) * (E - r - 2), "diff-2")]:
    g0, gs = slope_check(hc)
    check(f"{tag} moment slope: G(0)=0 (membership) and G(1)=sum_i(a_i(0)b_-i(i)-a_-i(i)b_i(0)) "
          f"(Lemma P, symbolic r) -- Q_0=1 needs this slope = 1", g0 and gs)
print("   => the finite-cap COKERNEL of this slope/adjoint system is exactly the")
print("      bounded emptiness of §5 (E - R not in Im Phi at cap d).")


# ===========================================================================
print("\n=== §5  TASK 4: bounded emptiness certificates per class ===")
# ===========================================================================
# The FULL positive cascade Q_4=Q_3=Q_2=Q_1=0 with Q_0=1 and genuine membership
# (E)_j | a_-j, b_-j, gauge b_3=0, top a_3 = h h^[1]h^[2], wall b_2 = kappa h h^[1].
# Emptiness = the ideal of coefficient equations is the UNIT ideal (1 in ideal).
def build_sector(hc, Dcap, rval=None, kval=None):
    hc = sp.expand(hc)
    if rval is not None:
        hc = sp.expand(hc.subs(r, rval))
    ks = kappa if kval is None else sp.Integer(kval)
    h1, h2 = sh(hc, 1), sh(hc, 2)
    a3 = sp.expand(hc * h1 * h2); b2 = sp.expand(ks * hc * h1)
    a2, a2c = gen('a2', Dcap + 2); a1, a1c = gen('a1', Dcap + 2); a0, a0c = gen('a0', Dcap + 2)
    b1, b1c = gen('b1', Dcap + 3); b0, b0c = gen('b0', Dcap + 3)
    am1 = sp.expand(E * gen('am1', Dcap)[0]); am1c = gen('am1', Dcap)[1]
    am2 = sp.expand(E * (E - 1) * gen('am2', Dcap)[0]); am2c = gen('am2', Dcap)[1]
    am3 = sp.expand(E * (E - 1) * (E - 2) * gen('am3', Dcap)[0]); am3c = gen('am3', Dcap)[1]
    bm1 = sp.expand(E * gen('bm1', Dcap)[0]); bm1c = gen('bm1', Dcap)[1]
    bm2 = sp.expand(E * (E - 1) * gen('bm2', Dcap)[0]); bm2c = gen('bm2', Dcap)[1]
    bm3 = sp.expand(E * (E - 1) * (E - 2) * gen('bm3', Dcap)[0]); bm3c = gen('bm3', Dcap)[1]
    A = {3: a3, 2: a2, 1: a1, 0: a0, -1: am1, -2: am2, -3: am3}
    B = {3: sp.Integer(0), 2: b2, 1: b1, 0: b0, -1: bm1, -2: bm2, -3: bm3}
    allc = a2c + a1c + a0c + b1c + b0c + am1c + am2c + am3c + bm1c + bm2c + bm3c
    eqs = []
    for m in [4, 3, 2, 1]:
        eqs += sp.Poly(Qm(A, B, m), E).all_coeffs()
    eqs += sp.Poly(Qm(A, B, 0) - 1, E).all_coeffs()
    return [sp.expand(e) for e in eqs if sp.expand(e) != 0], allc


def sector_empty_sympy(hc, Dcap, rval=None, kval=None):
    eqs, allc = build_sector(hc, Dcap, rval, kval)
    Gb = sp.groebner(eqs, *allc, order='grevlex')
    return list(Gb) == [sp.Integer(1)]


# (5a) COMMITTED, GENERIC FIBER (r,kappa): cap d=2 unit ideal after treating
#      r,kappa as coefficient parameters, both classes. SymPy displays the
#      coefficient domain as ZZ[r,kappa], but this calculation certifies the
#      generic fiber over QQ(r,kappa), not uniform emptiness at every specialization.
t0 = time.time()
check("DIFF-1 bounded (cap d=2, exact SymPy; displayed domain ZZ[r,kappa], generic fiber QQ(r,kappa)): sector EMPTY",
      sector_empty_sympy((E - r) * (E - r - 1), 2))
check("DIFF-2 bounded (cap d=2, exact SymPy; displayed domain ZZ[r,kappa], generic fiber QQ(r,kappa)): sector EMPTY",
      sector_empty_sympy((E - r) * (E - r - 2), 2))
print(f"   (§5a generic-fiber (r,kappa) SymPy Groebner time: {time.time() - t0:.1f}s)", flush=True)

# (5b) COMMITTED, SPECIFIC (r,kappa): cap d=2 UNIT IDEAL at integer/special r
#      (covers r-values off the generic denominator locus), both classes.
t0 = time.time()
GRID = [(0, 1), (1, 1), (-1, 2), (2, 1)]
for hc, tag in [((E - r) * (E - r - 1), "diff-1"), ((E - r) * (E - r - 2), "diff-2")]:
    all_empty = all(sector_empty_sympy(hc, 2, rval, kval) for rval, kval in GRID)
    check(f"{tag} bounded (cap d=2, exact SymPy over QQ): sector EMPTY at all "
          f"(r,kappa) in {GRID}", all_empty)
print(f"   (§5b specific-(r,kappa) SymPy grid time: {time.time() - t0:.1f}s)", flush=True)

# (5c) the multiple-root class needs NO emptiness search: it is cube-separated
#      (§2), hence closed at arbitrary degree by the shifted-cube theorem.
check("MULT-ROOT (E-r)^2 closed ARBITRARY DEGREE (cube-separated; §2) -- no bounded "
      "search needed (h forced constant by the clean theorem)",
      all(gcd_deg((E - r)**2, sh((E - r)**2, j)) == 0 for j in (1, 2, 3)))

# (5d) HEAVY: msolve, exact QQ, sampled (r,kappa) grid at cap d=2, d=3, d=4, both
#      classes.  Each is instant (0-dimensional/unit).  SKIP cleanly if absent.
HEAVY_GRID = [(0, 1), (0, 2), (1, 1), (-1, 1), (2, 1), (sp.Rational(1, 2), 1)]
if HEAVY and HAVE_MSOLVE:
    for Dcap in (2, 3, 4):
        for hc, tag in [((E - r) * (E - r - 1), "diff1"), ((E - r) * (E - r - 2), "diff2")]:
            allok = True
            timed_out = []
            for rval, kval in HEAVY_GRID:
                eqs, allc = build_sector(hc, Dcap, rval, kval)
                u = msolve_unit(eqs, allc, f"{tag}_D{Dcap}_r{sp.nsimplify(rval)}_k{kval}".replace('/', '_'))
                if u is None:
                    timed_out.append((rval, kval))
                else:
                    allok = allok and u
            if timed_out:
                skip(f"RESIDUAL HEAVY cap d={Dcap} grid [{tag}]", f"msolve timeout at {timed_out}")
            if len(timed_out) < len(HEAVY_GRID):
                check(f"RESIDUAL HEAVY (cap d={Dcap}, msolve exact QQ): sector EMPTY at every "
                      f"completed (r,kappa) grid point  [{tag}]", allok)
else:
    why = "set HEAVY=1" if not HEAVY else "msolve not found"
    skip("RESIDUAL HEAVY: msolve emptiness grid at cap d=2,3,4 (diff-1, diff-2)", why)


# ===========================================================================
print("\n=== §6  bookkeeping (see broken-separation.md §6) ===")
# ===========================================================================
print("   Band-3 shifted-cube broken-separation status after this memo:")
print("     TASK 1 (degraded Q_4): EXACT per class, symbolic r --")
print("        diff-1: (E-r)|a_2 + coupled correction a_2(r-1)+a_2(r+1)=0;")
print("        diff-2: clean proper-factor (E-r)(E-r-1)|a_2 (no coupling);")
print("        mult-root (E-r)^2: cube-separated => CLOSED arbitrary degree.")
print("     TASK 2 (cascade): general Q_3 identity; a_3 double nodes; generic Q_4")
print("        does NOT extend to Q_3; tested normal forms establish only ideal")
print("        membership/nonmembership, not radical or geometric restoration status.")
print("     TASK 3 (adjoint): jet/derivative-node criterion at double nodes;")
print("        Lemma-P moment slope for the degenerate tops.")
print("     TASK 4 (bounded): cap d=2 EMPTY on the generic fiber and tested")
print("        specific (r,kappa) fibers, both classes; no uniform parameter claim")
print("        (committed); msolve cap d=2,3,4 grid (HEAVY).")
print("   OPEN (not claimed): diff-1, diff-2 at ARBITRARY degree (bounded only).")

# ===========================================================================
print("\n" + "=" * 72)
n_pass = sum(1 for _, ok in CHECKS if ok)
n_tot = len(CHECKS)
n_skip = len(SKIPS)
print(f"{n_pass}/{n_tot} checks passed; {n_skip} skipped", flush=True)
if SKIPS:
    print("SKIPPED (not claimed in the final status):")
    for name, why in SKIPS:
        print(f"   SKIP: {name}  ({why})")
if n_pass == n_tot:
    if n_skip:
        print("ALL EXECUTED BROKEN SEPARATION CHECKS PASSED; OPTIONAL CHECKS SKIPPED")
    else:
        print("ALL BROKEN SEPARATION CHECKS PASSED; NO SKIPS")
else:
    print("SOME CHECKS FAILED:")
    for name, ok in CHECKS:
        if not ok:
            print("   FAIL:", name)
    raise SystemExit(1)
