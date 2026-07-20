#!/usr/bin/env python3
"""
audit_classical_A_identities.py  --  Referee audit (auditor A) of the CLASSICAL
half of provisional Theorem J3 (band-2 two-sided non-square sector).

This script certifies, by exact symbolic computation with FULLY SYMBOLIC
polynomial coefficients (several degree profiles, not spot checks):

  Sec 0. Orientation audit: the memo bracket {G,F} = G_xi F_x - G_x F_xi vs the
         bracket used in verify_J3.py (which is its negative).
  Sec 1. From-scratch derivation of all nine C_m from the two-variable Poisson
         bracket, compared to the M4 memo's displayed equations, per profile.
  Sec 2. Cascade reduction identities (steps m=+-4, +-3, +-2, +-1, 0), each as
         a polynomial identity in all symbolic coefficients, with the exact
         residual terms displayed (b0 kept generic).
  Sec 3. Product-invariant identities of the endgame.
  Sec 4. Lemma certificates:
         (W)  Wronskian lemma  a b' - a' b = 0, a != 0  =>  b = const * a:
              - the rational identity W = a^2 (b/a)',
              - top-coefficient identity (e-d) a_d b_e (kills deg b != deg a),
              - bounded-exhaustive kernel-dimension certificate via minors +
                Rabinowitsch saturation (all a with lead != 0, d = 1..4).
         (J2) 2 a u' - a' u = 0: complete polynomial solution set:
              - identity u*(2au'-a'u) = a*(u^2)' - a'*u^2  (so u^2 = c*a),
              - converse: a = d h^2 admits u = e h,
              - top-coefficient identity (2e-d) a_d u_e (kills deg u != d/2,
                hence odd-degree a has only u = 0 -- complete for ALL u),
              - bounded-exhaustive certificates: d=2 kernel nonzero iff
                discriminant = 0 (exact ideal computation, saturated at lead),
                d=4 square-parametrization direction + generic rank.

No tau-divisibility is used anywhere (matching the theorem's claim that the
statement lives in the localized algebra C[x^{+-1}, xi]).

Exit code 0 iff every check passes.
"""
import sys
from itertools import combinations

import sympy as sp

x, xi, tau, z = sp.symbols("x xi tau z")
lam, mu = sp.symbols("lambda2 mu2")
Dt = lambda f: sp.diff(f, tau)

FAILED = []


def check(name, ok):
    print(("  [PASS] " if ok else "  [FAIL] ") + name)
    if not ok:
        FAILED.append(name)


def gpoly(prefix, deg):
    """Generic polynomial with fresh symbolic coefficients. deg None -> zero poly."""
    if deg is None:
        return sp.Integer(0), ()
    cs = sp.symbols(f"{prefix}0:{deg + 1}")
    return sum(cs[i] * tau**i for i in range(deg + 1)), cs


def PB_memo(G, F):
    """Frozen M4-memo orientation: {G,F} = G_xi F_x - G_x F_xi (so {xi,x}=1)."""
    return sp.expand(sp.diff(G, xi) * sp.diff(F, x) - sp.diff(G, x) * sp.diff(F, xi))


def PB_j3script(G, F):
    """The orientation coded in verify_J3.py: PB(f,g)=f_x g_xi - f_xi g_x."""
    return sp.expand(sp.diff(G, x) * sp.diff(F, xi) - sp.diff(G, xi) * sp.diff(F, x))


def x_components(expr, lo=-8, hi=8):
    e = sp.expand(expr.subs(xi, tau / x))
    comp = {}
    for k in range(lo, hi + 1):
        c = sp.expand(e.coeff(x, k))
        if c != 0:
            comp[k] = c
    assert sp.expand(e - sum(x**k * v for k, v in comp.items())) == 0, \
        "component extraction incomplete"
    return comp


def Cm(m, A, B):
    """Memo formula: C_m = sum_{k+l=m} (k a_k b_l' - l a_k' b_l)."""
    s = sp.Integer(0)
    for k in range(-2, 3):
        l = m - k
        if -2 <= l <= 2:
            s += k * A[k] * Dt(B[l]) - l * Dt(A[k]) * B[l]
    return sp.expand(s)


# ----------------------------------------------------------------------------
# Degree profiles: (a2,a1,a0,am1,am2 ; b2,b1,b0,bm1,bm2). None = zero poly.
# Symbolic-coefficient identities specialize downward, so the largest profile
# already subsumes the others; several are run anyway as redundancy, including
# mismatched b-degrees and vanishing middle coefficients.
# ----------------------------------------------------------------------------
PROFILES = [
    ("P1 all deg 3",            (3, 3, 3, 3, 3), (3, 3, 3, 3, 3)),
    ("P2 a deg 2 / b deg 4",    (2, 2, 2, 2, 2), (4, 4, 4, 4, 4)),
    ("P3 mixed (4,1,2,0,3)",    (4, 1, 2, 0, 3), (3, 2, 4, 1, 3)),
    ("P4 a1=a-1=0",             (3, None, 2, None, 3), (3, 3, 3, 3, 3)),
]


def run_profile(tag, adeg, bdeg):
    print(f"Profile {tag}: a-degs {adeg}, b-degs {bdeg}")
    a2, _ = gpoly("p", adeg[0]);  a1, _ = gpoly("q", adeg[1])
    a0, a0c = gpoly("g", adeg[2]); am1, _ = gpoly("r", adeg[3])
    am2, _ = gpoly("s", adeg[4])
    b2, _ = gpoly("B2_", bdeg[0]); b1, _ = gpoly("B1_", bdeg[1])
    b0, b0c = gpoly("B0_", bdeg[2]); bm1, _ = gpoly("Bm1_", bdeg[3])
    bm2, _ = gpoly("Bm2_", bdeg[4])
    A = {2: a2, 1: a1, 0: a0, -1: am1, -2: am2}
    B = {2: b2, 1: b1, 0: b0, -1: bm1, -2: bm2}

    # --- Sec 0/1: from-scratch derivation & orientation ---
    lift = lambda f, k: f.subs(tau, x * xi) * x**k
    F = sum(lift(A[k], k) for k in range(-2, 3))
    G = sum(lift(B[k], k) for k in range(-2, 3))
    comp = x_components(PB_memo(G, F))
    check("direct {G,F} (memo orientation) x-components == memo C_m formulas, "
          "components confined to |m|<=4",
          set(comp) <= set(range(-4, 5)) and
          all(sp.expand(comp.get(m, 0) - Cm(m, A, B)) == 0 for m in range(-4, 5)))
    comp_flip = x_components(PB_j3script(G, F))
    check("verify_J3.py bracket == MINUS the memo bracket (orientation flag)",
          all(sp.expand(comp_flip.get(m, 0) + Cm(m, A, B)) == 0 for m in range(-4, 5)))

    # --- Sec 2: cascade reduction identities (memo orientation throughout) ---
    W = lambda a, b: sp.expand(a * Dt(b) - Dt(a) * b)

    # m=+-4
    check("C_4 == 2*(a2 b2' - a2' b2)  [Wronskian form]",
          sp.expand(Cm(4, A, B) - 2 * W(a2, b2)) == 0)
    check("C_-4 == -2*(a-2 b-2' - a-2' b-2)  [mirror Wronskian form]",
          sp.expand(Cm(-4, A, B) + 2 * W(am2, bm2)) == 0)

    # m=+3 with b2 -> lam*a2, b1 generic; u := b1 - lam*a1
    B3 = dict(B); B3[2] = lam * a2
    u = b1 - lam * a1
    check("C_3|_{b2=lam a2} == 2 a2 u' - a2' u  with u = b1 - lam a1 (b1 generic)",
          sp.expand(Cm(3, A, B3) - (2 * a2 * Dt(u) - Dt(a2) * u)) == 0)

    # m=-3 with b-2 -> mu*a-2, b-1 generic; v := b-1 - mu*a-1  (mirror ORIENTATION)
    Bm3 = dict(B); Bm3[-2] = mu * am2
    v = bm1 - mu * am1
    check("C_-3|_{b-2=mu a-2} == -(2 a-2 v' - a-2' v)  with v = b-1 - mu a-1 "
          "(same homogeneous form as (H3); sign is a global factor)",
          sp.expand(Cm(-3, A, Bm3) + (2 * am2 * Dt(v) - Dt(am2) * v)) == 0)

    # m=+-2 with proportional b's at levels 1,2 (resp. -1,-2); b0 generic
    B2s = dict(B); B2s[2] = lam * a2; B2s[1] = lam * a1
    check("C_2|_{b2=lam a2, b1=lam a1} == 2 a2 (b0' - lam a0')  "
          "[(1,1) terms cancel identically]",
          sp.expand(Cm(2, A, B2s) - 2 * a2 * (Dt(b0) - lam * Dt(a0))) == 0)
    check("(1,1) cancellation alone: a1 b1' - a1' b1 == 0 under b1 = lam a1",
          sp.expand(a1 * Dt(lam * a1) - Dt(a1) * lam * a1) == 0)
    Bm2s = dict(B); Bm2s[-2] = mu * am2; Bm2s[-1] = mu * am1
    check("C_-2|_{b-2=mu a-2, b-1=mu a-1} == -2 a-2 (b0' - mu a0')",
          sp.expand(Cm(-2, A, Bm2s) + 2 * am2 * (Dt(b0) - mu * Dt(a0))) == 0)

    # m=+-1: all four outer proportionalities, b0 generic; exact residuals
    Bp = dict(B); Bp[2] = lam * a2; Bp[1] = lam * a1
    Bp[-1] = mu * am1; Bp[-2] = mu * am2
    U = sp.expand(Dt(a2) * am1 + 2 * a2 * Dt(am1))
    L = sp.expand(Dt(am2) * a1 + 2 * am2 * Dt(a1))
    check("C_1|_{4 props} == (mu-lam)*U + a1*(b0' - lam a0')   (b0 generic)",
          sp.expand(Cm(1, A, Bp) - ((mu - lam) * U + a1 * (Dt(b0) - lam * Dt(a0)))) == 0)
    check("C_-1|_{4 props} == (mu-lam)*L - a-1*(b0' - mu a0')  (b0 generic)",
          sp.expand(Cm(-1, A, Bp) - ((mu - lam) * L - am1 * (Dt(b0) - mu * Dt(a0)))) == 0)

    # m=0: full telescoping; b0, a0 play no role at all
    M = 2 * a2 * am2 + a1 * am1
    check("C_0|_{4 props} == (mu-lam) * d/dtau(2 a2 a-2 + a1 a-1)   "
          "(memo orientation; prose's (lam-mu) matches the flipped bracket)",
          sp.expand(Cm(0, A, Bp) - (mu - lam) * Dt(M)) == 0)
    c0 = Cm(0, A, B)
    check("C_0 contains neither b0 nor a0 structurally",
          all(sp.diff(c0, s_) == 0 for s_ in (list(b0c) + list(a0c))))

    # --- Sec 3: endgame product invariants ---
    check("a-1 * U == (a2 a-1^2)'   and   a1 * L == (a-2 a1^2)'",
          sp.expand(am1 * U - Dt(a2 * am1**2)) == 0 and
          sp.expand(a1 * L - Dt(am2 * a1**2)) == 0)
    print()


# ----------------------------------------------------------------------------
# Sec 4: lemma certificates
# ----------------------------------------------------------------------------
def lemma_certificates():
    print("Lemma certificates (bounded-exhaustive, symbolic coefficients)")

    # (J2) key identity: u*(2au' - a'u) == a*(u^2)' - a'*u^2  == a^2 (u^2/a)' * (1/a) ...
    a_, ac = gpoly("A", 6)
    u_, uc = gpoly("U", 6)
    check("J2 identity: u*(2 a u' - a' u) == a*(u^2)' - a'*u^2  (deg 6 generic)",
          sp.expand(u_ * (2 * a_ * Dt(u_) - Dt(a_) * u_)
                    - (a_ * Dt(u_**2) - Dt(a_) * u_**2)) == 0)
    # (W) rational identity a^2 (b/a)' == a b' - a' b
    b_, bc = gpoly("Bb", 6)
    check("W identity: a^2 * (b/a)' == a b' - a' b  (deg 6 generic)",
          sp.simplify(a_**2 * sp.diff(b_ / a_, tau) - (a_ * Dt(b_) - Dt(a_) * b_)) == 0)
    # J2 converse: a = d h^2 admits u = e h
    d_, e_ = sp.symbols("d e")
    h_, hc = gpoly("H", 4)
    check("J2 converse: a = d h^2, u = e h solves 2 a u' - a' u = 0 (h deg 4 generic)",
          sp.expand(2 * (d_ * h_**2) * Dt(e_ * h_) - Dt(d_ * h_**2) * (e_ * h_)) == 0)

    # Top-coefficient identities (complete degree control, all e)
    ok_all = True
    for d in range(1, 5):
        A_, Ac = gpoly(f"Ad{d}_", d)
        for e in range(0, 9):
            Ue, Uec = gpoly(f"Ue{d}_{e}_", e)
            expr = sp.expand(2 * A_ * Dt(Ue) - Dt(A_) * Ue)
            top = expr.coeff(tau, d + e - 1) if d + e - 1 >= 0 else 0
            ok_all &= sp.expand(top - (2 * e - d) * Ac[d] * Uec[e]) == 0
    check("J2 top coeff of (2 a u' - a' u) at tau^{d+e-1} == (2e-d)*a_d*u_e "
          "for d=1..4, e=0..8  => nonzero u forces deg u = d/2; odd d => u=0",
          ok_all)
    ok_all = True
    for d in range(1, 5):
        A_, Ac = gpoly(f"Wd{d}_", d)
        for e in range(0, 9):
            Be, Bec = gpoly(f"We{d}_{e}_", e)
            expr = sp.expand(A_ * Dt(Be) - Dt(A_) * Be)
            top = expr.coeff(tau, d + e - 1) if d + e - 1 >= 0 else 0
            ok_all &= sp.expand(top - (e - d) * Ac[d] * Bec[e]) == 0
    check("W top coeff of (a b' - a' b) at tau^{d+e-1} == (e-d)*a_d*b_e "
          "for d=1..4, e=0..8  => nonzero kernel b forces deg b = deg a",
          ok_all)

    # (W) bounded-exhaustive kernel dimension: for every a of degree exactly d
    # (lead invertible), kernel of b |-> a b' - a' b restricted to deg b <= d is
    # exactly the line C*a: rank == d for ALL specializations with a_d != 0.
    for d in range(1, 5):
        A_, Ac = gpoly(f"WK{d}_", d)
        bsyms = sp.symbols(f"wb{d}_0:{d + 1}")
        Bgen = sum(bsyms[i] * tau**i for i in range(d + 1))
        expr = sp.expand(A_ * Dt(Bgen) - Dt(A_) * Bgen)
        rows = sp.Poly(expr, tau).degree() if expr != 0 else 0
        Mrows = []
        for i in range(2 * d):
            Mrows.append([sp.expand(sp.diff(expr.coeff(tau, i), bs)) for bs in bsyms])
        Mmat = sp.Matrix(Mrows)
        # a itself is in the kernel:
        kerA = sp.Matrix([[Ac[i]] for i in range(d + 1)])
        in_ker = all(sp.expand(v) == 0 for v in (Mmat * kerA))
        # all d x d minors; saturate at a_d: Groebner of (minors, 1 - z*a_d) == (1)
        minors = []
        rows_idx = range(Mmat.rows); cols_idx = range(Mmat.cols)
        for rs in combinations(rows_idx, d):
            for cs_ in combinations(cols_idx, d):
                minors.append(sp.expand(Mmat[list(rs), list(cs_)].det()))
        gb = sp.groebner(minors + [1 - z * Ac[d]], *(list(Ac) + [z]), order="grevlex")
        check(f"W kernel certificate d={d}: a in kernel; rank == d for ALL a with "
              f"a_d != 0 (minors ideal + Rabinowitsch == (1)) => kernel = C*a",
              in_ker and list(gb.exprs) == [sp.Integer(1)])

    # (J2) even-degree certificates.
    # d=2, e=1: kernel nonzero  <=>  disc(a) = 0   (on a_2 != 0)
    A_, Ac = gpoly("JA2_", 2)
    usyms = sp.symbols("ju0:2")
    Ugen = usyms[0] + usyms[1] * tau
    expr = sp.expand(2 * A_ * Dt(Ugen) - Dt(A_) * Ugen)
    Mmat = sp.Matrix([[sp.diff(expr.coeff(tau, i), us) for us in usyms] for i in range(3)])
    minors = [sp.expand(Mmat[list(rs), :].det()) for rs in combinations(range(3), 2)]
    disc = sp.expand(Ac[1]**2 - 4 * Ac[2] * Ac[0])
    # direction 1: disc=0 (with a2 != 0) => all minors vanish: substitute the
    # double-root parametrization a = a2*(tau+t0)^2
    t0 = sp.Symbol("t0")
    subs_sq = {Ac[0]: Ac[2] * t0**2, Ac[1]: 2 * Ac[2] * t0}
    dir1 = all(sp.expand(m.subs(subs_sq)) == 0 for m in minors)
    # ... and the kernel vector is u = tau + t0:
    dir1b = sp.expand(expr.subs(subs_sq).subs({usyms[0]: t0, usyms[1]: 1})) == 0
    # direction 2: minors=0 & a2 != 0 => disc = 0: disc in radical of
    # (minors) saturated at a2.
    gb = sp.groebner(minors + [1 - z * Ac[2]], *(list(Ac) + [z]), order="grevlex")
    red = gb.reduce(disc)[1]
    dir2 = (red == 0) or (gb.reduce(sp.expand(disc**2))[1] == 0)
    check("J2 d=2 exhaustive: kernel nonzero <=> disc==0 (squares), for ALL "
          "a2-lead != 0; kernel vector is the square root", dir1 and dir1b and dir2)

    # d=4, e=2: squares have kernel (parametrized direction) + generic rank full
    A4c = sp.symbols("A4_0:5")
    A4 = sum(A4c[i] * tau**i for i in range(5))
    u2 = sp.symbols("ju2_0:3")
    U2 = u2[0] + u2[1] * tau + u2[2] * tau**2
    expr = sp.expand(2 * A4 * Dt(U2) - Dt(A4) * U2)
    M4m = sp.Matrix([[sp.diff(expr.coeff(tau, i), us) for us in u2] for i in range(6)])
    pp, qq, cc = sp.symbols("pp qq cc")
    sq_subs = {A4c[4]: cc, A4c[3]: 2 * pp * cc, A4c[2]: cc * (pp**2 + 2 * qq),
               A4c[1]: 2 * pp * qq * cc, A4c[0]: cc * qq**2}   # a = cc*(t^2+pp*t+qq)^2
    Msq = M4m.subs(sq_subs)
    kvec = sp.Matrix([[qq], [pp], [1]])
    dir1 = all(sp.expand(v) == 0 for v in (Msq * kvec))
    dir2 = M4m.rank() == 3   # generic coefficients: kernel zero
    # a fixed nonsquare quartic has zero kernel:
    spec = {A4c[0]: 1, A4c[1]: 1, A4c[2]: 0, A4c[3]: 0, A4c[4]: 1}
    dir3 = M4m.subs(spec).rank() == 3
    check("J2 d=4: scalar-square quartics admit kernel u = sqrt; generic and "
          "sample nonsquare quartics have zero kernel (converse is the exact "
          "identity u^2 = c*a)", dir1 and dir2 and dir3)
    print()


def main():
    print("=" * 78)
    print("AUDIT A (classical J3): forward cascade identities and lemma certificates")
    print("=" * 78)
    for tag, adeg, bdeg in PROFILES:
        run_profile(tag, adeg, bdeg)
    lemma_certificates()
    if FAILED:
        print("FAILURES:", *FAILED, sep="\n  ")
        sys.exit(1)
    print("ALL AUDIT-A IDENTITY/LEMMA CHECKS PASSED")


if __name__ == "__main__":
    main()
