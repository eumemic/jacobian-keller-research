#!/usr/bin/env python3
"""Affine-quotient obstruction for the explicit telescoping Keller maps F3, F4, F5.

Milestone: research/jc2-attack-memo.md, Section 4 ("affine-quotient obstruction").

Question.  For F in {F3, F4, F5} (the paper's degree-3, S4, and S5 members of the
telescoping family, all with det Jac = -2), do there exist
  * an affine-linear submersion pi = (l1, l2), l_i = r_i x + s_i y + t_i z + u_i,
    with 2x3 linear-part matrix of rank 2, and
  * polynomials G1, G2 in two variables with deg G_i <= deg F,
such that pi o F = G o pi identically and det Jac(G) is a nonzero constant?

Method (all exact, over QQ; conclusions valid over CC).

  Lemma 1 (fiber-derivative criterion).  Fix pi with rank-2 linear part M and let
  v != 0 span ker M.  For a polynomial L in C[x,y,z]:
      L in C[l1, l2]  <=>  d_v L == 0  (directional derivative along v),
  because (l1, l2, lam) is an affine coordinate system for any affine lam with
  lam(v) = 1, and d_v = d/d(lam) in those coordinates.  Consequently, for fixed pi,
  a polynomial G with pi o F = G o pi exists iff
      d_v(l_i o F) == 0   for i = 1, 2,
  and then G_i is unique with deg G_i = deg(l_i o F) <= deg F automatically; the
  translation parts u_i are irrelevant (d_v kills constants).  This removes G
  (the linear level of the two-level structure) entirely.

  Lemma 2 (gauge reduction).  Target-affine changes (pi, G) -> (A pi + b, G')
  with A in GL2, G' = (A,b) o G o (A,b)^{-1} preserve existence, degrees, and
  det Jac(G).  Writing M = A R with R the reduced row-echelon form, every rank-2
  pi is gauge-equivalent to exactly one of three charts:
      C12: l1 = x + p z,  l2 = y + q z     (kernel v = (-p, -q, 1))
      C13: l1 = x + p y,  l2 = z           (kernel v = (-p, 1, 0))
      C23: l1 = y,        l2 = z           (kernel v = (1, 0, 0))
  So the full search is three charts with at most two scalar unknowns each.

  Per chart we expand E_i = d_v(l_i o F) in x, y, z and collect the coefficient
  ideal I in QQ[p, q].  A generator that is a nonzero rational constant (a "unit
  certificate") proves I = (1), hence no solutions in any field extension of QQ;
  we also compute a Groebner basis of I over QQ and require it to be [1].
  Rank < 2 linear parts are excluded by the definition of submersion.

Cross-checks.
  * Sanity gate: each member is polynomial with det Jac = -2 (as in paper/verify.py),
    with the expected component degrees and monomial leading forms.
  * Positive controls: two synthetic maps with known affine quotients confirm the
    pipeline finds (pi, G), verifies pi o F = G o pi exactly, and distinguishes
    constant from non-constant det Jac(G).
  * Lemma-free sampled check: for sampled chart parameters, direct subalgebra
    membership l_i o F in QQ[l1, l2] is refuted by exact linear algebra mod a
    prime (infeasibility mod p implies infeasibility over QQ), independently of
    Lemma 1.

Outcome A (a quotient found for a Keller member) would be flagged loudly as a
JC2 counterexample candidate and the script exits with code 2 demanding manual
triple-verification (collision descent).  Outcome B prints the final line
    OBSTRUCTION CERTIFIED: no affine quotient for F3, F4, F5

Run:  uv run --with sympy python run_obstruction.py
Writes certificates.json next to this file.  Deterministic (fixed RNG seed).
"""

import json
import os
import random
import sys

import sympy as sp

x, y, z, nu = sp.symbols("x y z nu")
p, q = sp.symbols("p q")
w1, w2 = sp.symbols("w1 w2")
t = x * y

PRIME = 1000003
RNG_SEED = 20260720

CERT = {"prime": PRIME, "rng_seed": RNG_SEED, "members": {}, "controls": {}, "modp_crosscheck": {}}


# ----------------------------------------------------------------------------
# The telescoping construction (paper/main.tex Section 2; cf. paper/verify.py).
# ----------------------------------------------------------------------------

def family(E, h):
    """Return [F1, F2, F3] for the telescoping member with data (E, h)."""
    A = sp.integrate(nu * E, nu)          # A(nu) = int_0^nu xi E(xi) dxi
    C = 2 * sp.integrate(E, nu)           # C(nu) = 2 int_0^nu E(xi) dxi
    u0 = 1 + t
    s = x**2 * z
    g = sp.expand(h - s)
    nux = sp.expand(u0 * g)
    S = sp.expand(A.subs(nu, nux) + u0 * g**2)
    T = sp.expand(C.subs(nu, nux) + 2 * g)
    return [sp.expand(sp.cancel(S / (x**2 * g**2))),
            sp.expand(sp.cancel(T / (x * g))),
            sp.expand(x * g)]


def jacobian_det(F):
    return sp.expand(sp.det(sp.Matrix([[sp.diff(f, v) for v in (x, y, z)] for f in F])))


MEMBERS = {
    # name: (E, h, expected component degrees)
    "F3": (2 - 6 * nu, 1 - sp.Rational(3, 2) * t, (7, 6, 4)),    # announced degree-3 member
    "F4": (3 - 12 * nu + 6 * nu**2, 1 - 2 * t, (12, 11, 4)),      # S4 member
    "F5": (4 - 15 * nu + 10 * nu**3, sp.Integer(1), (17, 16, 4)),  # S5 member
}


def leading_data(f):
    P = sp.Poly(f, x, y, z)
    d = P.total_degree()
    top = [(m, sp.nsimplify(c)) for m, c in P.terms() if sum(m) == d]
    return d, top


def sanity_gate():
    built = {}
    for name, (E, h, expected_degs) in MEMBERS.items():
        F = family(E, h)
        assert all(f.is_polynomial(x, y, z) for f in F), name + ": not polynomial"
        dj = jacobian_det(F)
        assert dj == -2, name + ": det Jac != -2"
        info = {"detJac": str(dj), "components": []}
        degs = []
        for i, f in enumerate(F, 1):
            d, top = leading_data(f)
            degs.append(d)
            entry = {"degree": d,
                     "leading_form": [{"monomial": list(m), "coeff": str(c)} for m, c in top]}
            info["components"].append(entry)
            if i <= 2:
                # leading form is a single monomial divisible by x*y*z
                assert len(top) == 1, name + (": F%d leading form not a monomial" % i)
                assert all(e >= 1 for e in top[0][0]), \
                    name + (": F%d leading monomial not divisible by xyz" % i)
        assert tuple(degs) == expected_degs, name + ": unexpected component degrees"
        assert degs[0] > degs[1] > degs[2], name + ": degrees not strictly decreasing"
        built[name] = F
        CERT["members"][name] = {"E": str(E), "h": str(h), "sanity": info}
        print("[gate] %s: polynomial, det Jac = -2, degrees %s, monomial leading forms OK"
              % (name, degs))
    return built


# ----------------------------------------------------------------------------
# Chart framework (Lemmas 1 and 2).
# ----------------------------------------------------------------------------

CHARTS = {
    "C12": dict(rows=[[1, 0, p], [0, 1, q]], v=(-p, -q, 1), unknowns=(p, q),
                section={x: w1, y: w2, z: 0}),
    "C13": dict(rows=[[1, p, 0], [0, 0, 1]], v=(-p, 1, 0), unknowns=(p,),
                section={x: w1, y: 0, z: w2}),
    "C23": dict(rows=[[0, 1, 0], [0, 0, 1]], v=(1, 0, 0), unknowns=(),
                section={x: 0, y: w1, z: w2}),
}


def chart_equations(F, chart):
    """E_i = d_v(l_i o F), i = 1, 2, as expanded polynomials in x, y, z over QQ[p,q]."""
    ch = CHARTS[chart]
    v = ch["v"]

    def dv(f):
        return sp.expand(v[0] * sp.diff(f, x) + v[1] * sp.diff(f, y) + v[2] * sp.diff(f, z))

    return [dv(sp.expand(r[0] * F[0] + r[1] * F[1] + r[2] * F[2])) for r in ch["rows"]]


def analyze_chart(F, chart):
    """Return a record deciding feasibility of the chart's coefficient ideal."""
    ch = CHARTS[chart]
    eqs = chart_equations(F, chart)
    gens = set()
    units = []
    for i, eq in enumerate(eqs, 1):
        if eq == 0:
            continue
        for m, c in sp.Poly(eq, x, y, z).terms():
            c = sp.expand(sp.nsimplify(c))
            if c != 0:
                gens.add(c)
                if c.is_Rational:
                    units.append({"equation": i, "monomial": list(m), "value": str(c)})
    gens = sorted(gens, key=sp.default_sort_key)
    rec = {"chart": chart, "n_generators": len(gens), "unit_certificates": units}
    if not ch["unknowns"]:
        # No parameters: feasible iff both E_i are identically zero.
        infeasible = len(gens) > 0
        rec["groebner_basis"] = ["1"] if infeasible else []
        rec["infeasible"] = infeasible
        rec["solutions"] = []
        return rec, gens, eqs
    if not gens:
        rec["groebner_basis"] = []
        rec["infeasible"] = False
        rec["solutions"] = ["all parameter values"]
        return rec, gens, eqs
    gb = sp.groebner(gens, *ch["unknowns"], order="grevlex")
    gbex = list(gb.exprs)
    rec["groebner_basis"] = [str(e) for e in gbex]
    rec["infeasible"] = gbex == [sp.Integer(1)]
    if units:
        assert rec["infeasible"], chart + ": unit certificate present but Groebner basis != [1]"
    rec["solutions"] = []
    if not rec["infeasible"]:
        sols = sp.solve(gens, list(ch["unknowns"]), dict=True)
        rec["solutions"] = sols
    return rec, gens, eqs


def reconstruct_quotient(F, chart, sol):
    """Given a feasible chart parameter value, rebuild (pi, G) and verify exactly."""
    ch = CHARTS[chart]
    rows = [[sp.sympify(e).subs(sol) for e in r] for r in ch["rows"]]
    l1 = rows[0][0] * x + rows[0][1] * y + rows[0][2] * z
    l2 = rows[1][0] * x + rows[1][1] * y + rows[1][2] * z
    L1 = sp.expand(rows[0][0] * F[0] + rows[0][1] * F[1] + rows[0][2] * F[2])
    L2 = sp.expand(rows[1][0] * F[0] + rows[1][1] * F[1] + rows[1][2] * F[2])
    G1 = sp.expand(L1.subs(ch["section"], simultaneous=True))
    G2 = sp.expand(L2.subs(ch["section"], simultaneous=True))
    ok = (sp.expand(G1.subs({w1: l1, w2: l2}, simultaneous=True) - L1) == 0 and
          sp.expand(G2.subs({w1: l1, w2: l2}, simultaneous=True) - L2) == 0)
    jac = sp.expand(sp.diff(G1, w1) * sp.diff(G2, w2) - sp.diff(G1, w2) * sp.diff(G2, w1))
    jac_const_nonzero = jac.is_constant() and jac != 0
    return (l1, l2), (G1, G2), ok, jac, jac_const_nonzero


# ----------------------------------------------------------------------------
# Positive controls: the machinery must FIND quotients when they exist.
# ----------------------------------------------------------------------------

def run_controls():
    ctrl1 = [x + (y + z)**2, y + x**3, z - x**3]
    # pi = (x, y+z), G = (w1 + w2^2, w2), J(G) = 1  -> full Outcome-A shape
    ctrl2 = [x + (y + z)**2, (y + z)**3 - z + x**5, z - x**5]
    # pi = (x, y+z), G = (w1 + w2^2, w2^3), J(G) = 3 w2^2 -> semiconjugation only
    for name, F, expect_jconst in [("CTRL1", ctrl1, True), ("CTRL2", ctrl2, False)]:
        found = []
        crec = {}
        for chart in CHARTS:
            rec, gens, _ = analyze_chart(F, chart)
            crec[chart] = {"groebner_basis": rec["groebner_basis"],
                           "infeasible": rec["infeasible"]}
            if not rec["infeasible"]:
                for sol in rec["solutions"]:
                    pi, G, ok, jac, jcn = reconstruct_quotient(F, chart, sol)
                    assert ok, name + ": reconstructed identity failed"
                    found.append((chart, sol, pi, G, jac, jcn))
                    crec[chart].setdefault("quotients", []).append(
                        {"solution": {str(k): str(v) for k, v in sol.items()},
                         "pi": [str(pi[0]), str(pi[1])],
                         "G": [str(G[0]), str(G[1])],
                         "JacG": str(jac),
                         "JacG_nonzero_constant": bool(jcn)})
        assert len(found) == 1, name + ": expected exactly one quotient, got %d" % len(found)
        chart, sol, pi, G, jac, jcn = found[0]
        assert chart == "C12" and sol == {p: 0, q: 1}, name + ": wrong chart/solution"
        assert jcn == expect_jconst, name + ": Jacobian constancy misclassified"
        CERT["controls"][name] = crec
        print("[control] %s: quotient found in %s at %s, G = %s, Jac(G) = %s (%s) -- as designed"
              % (name, chart, sol, G, jac,
                 "nonzero constant" if jcn else "NOT constant"))


# ----------------------------------------------------------------------------
# Lemma-free cross-check: sampled subalgebra membership refuted mod PRIME.
# ----------------------------------------------------------------------------

def rat_mod(r, pr):
    r = sp.Rational(r)
    return (int(r.p) % pr) * pow(int(r.q) % pr, -1, pr) % pr


def rank_modp(rows, pr):
    rows = [r[:] for r in rows]
    m = len(rows)
    n = len(rows[0]) if m else 0
    rank = 0
    for col in range(n):
        piv = next((i for i in range(rank, m) if rows[i][col] % pr), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        inv = pow(rows[rank][col], -1, pr)
        rows[rank] = [(v * inv) % pr for v in rows[rank]]
        for i in range(m):
            if i != rank and rows[i][col]:
                f = rows[i][col]
                rows[i] = [(a - f * b) % pr for a, b in zip(rows[i], rows[rank])]
        rank += 1
        if rank == m:
            break
    return rank


def membership_infeasible_modp(L, l1c, l2c, pr, rng):
    """Certify that no polynomial identity L = sum c_ab l1^a l2^b (a+b <= deg L) holds.

    Sampled linear system in the c_ab at random integer points, reduced mod pr.
    Any exact rational solution of the full identity would satisfy the sampled
    system mod pr, so inconsistency mod pr is a valid infeasibility certificate.
    (Consistency of the sample would be inconclusive; we then fall back to the
    caller's other test.)"""
    P = sp.Poly(L, x, y, z)
    d = P.total_degree()
    basis = [(a, s_ - a) for s_ in range(d + 1) for a in range(s_ + 1)]
    nrows = len(basis) + 20
    A, rhs = [], []
    for _ in range(nrows):
        X, Y, Z = (rng.randint(-40, 40) for _ in range(3))
        l1v = (l1c[0] * X + l1c[1] * Y + l1c[2] * Z) % pr
        l2v = (l2c[0] * X + l2c[1] * Y + l2c[2] * Z) % pr
        A.append([pow(l1v, a, pr) * pow(l2v, b, pr) % pr for a, b in basis])
        rhs.append(rat_mod(P(X, Y, Z), pr))
    rA = rank_modp(A, pr)
    rAb = rank_modp([row + [v] for row, v in zip(A, rhs)], pr)
    return dict(deg=d, nbasis=len(basis), nrows=nrows, rankA=rA, rankAug=rAb,
                infeasible=bool(rAb > rA))


def run_modp_crosscheck(built):
    rng = random.Random(RNG_SEED)
    plan = [("C12", [(0, 0), (1, 2), (-3, 5)]),
            ("C13", [(0,), (2,), (-7,)]),
            ("C23", [()])]
    for name, F in built.items():
        recs = []
        for chart, samples in plan:
            for smp in samples:
                if chart == "C12":
                    pv, qv = smp
                    l1c, l2c = (1, 0, pv), (0, 1, qv)
                    tests = [("l2oF", F[1] + qv * F[2]), ("l1oF", F[0] + pv * F[2])]
                elif chart == "C13":
                    (pv,) = smp
                    l1c, l2c = (1, pv, 0), (0, 0, 1)
                    tests = [("l2oF", F[2]), ("l1oF", F[0] + pv * F[1])]
                else:
                    l1c, l2c = (0, 1, 0), (0, 0, 1)
                    tests = [("l2oF", F[2]), ("l1oF", F[1])]
                joint_infeasible = False
                detail = []
                for lbl, L in tests:
                    r = membership_infeasible_modp(L, l1c, l2c, PRIME, rng)
                    detail.append({"test": lbl, **r})
                    if r["infeasible"]:
                        joint_infeasible = True
                        break
                assert joint_infeasible, \
                    "%s %s %s: sampled membership unexpectedly consistent" % (name, chart, smp)
                recs.append({"chart": chart, "sample": list(smp), "tests": detail,
                             "joint_infeasible": True})
        CERT["modp_crosscheck"][name] = recs
        print("[mod-p] %s: all %d sampled chart parameters refuted independently of Lemma 1"
              % (name, len(recs)))


# ----------------------------------------------------------------------------
# Main.
# ----------------------------------------------------------------------------

def main():
    print("Affine-quotient obstruction computation (memo Section 4 milestone)")
    print("sympy", sp.__version__)
    built = sanity_gate()
    run_controls()

    outcome_A = []
    for name, F in built.items():
        for chart in CHARTS:
            rec, gens, _ = analyze_chart(F, chart)
            CERT["members"][name].setdefault("charts", {})[chart] = {
                k: rec[k] for k in ("n_generators", "unit_certificates",
                                    "groebner_basis", "infeasible")}
            if rec["infeasible"]:
                nunits = len(rec["unit_certificates"])
                print("[main] %s %s: INFEASIBLE (%d generators, %d unit certificates, GB = [1])"
                      % (name, chart, rec["n_generators"], nunits))
                continue
            # ---- Outcome A branch (never reached if unit certificates exist) ----
            print("!" * 78)
            print("!! %s %s: chart FEASIBLE -- investigating candidate quotient" % (name, chart))
            for sol in rec["solutions"]:
                pi, G, ok, jac, jcn = reconstruct_quotient(F, chart, sol)
                print("!!   solution %s  pi = %s" % (sol, pi))
                print("!!   G = %s" % (G,))
                print("!!   identity pi o F == G o pi: %s ; Jac(G) = %s" % (ok, jac))
                CERT["members"][name]["charts"][chart].setdefault("quotients", []).append(
                    {"solution": {str(k): str(v) for k, v in sol.items()},
                     "pi": [str(pi[0]), str(pi[1])], "G": [str(G[0]), str(G[1])],
                     "identity_verified": bool(ok), "JacG": str(jac),
                     "JacG_nonzero_constant": bool(jcn)})
                if ok and jcn:
                    outcome_A.append((name, chart, sol, pi, G, jac))
            print("!" * 78)

    run_modp_crosscheck(built)
    write_certificates()

    if outcome_A:
        print("=" * 78)
        print("JC2 COUNTEREXAMPLE CANDIDATE: an affine quotient with constant nonzero")
        print("Jacobian was found for a Keller member. THIS REQUIRES INDEPENDENT")
        print("TRIPLE-VERIFICATION, including descent of exact collisions of F to")
        print("collisions of G (F noninjective => G noninjective if a colliding pair")
        print("separates under pi). Do not announce before manual verification.")
        for item in outcome_A:
            print("  ", item)
        print("=" * 78)
        sys.exit(2)

    print("All charts infeasible for all three members; the coefficient ideals are")
    print("unit ideals over QQ, hence empty over CC. No polynomial G exists for any")
    print("rank-2 affine pi (any deg G, no Jacobian condition needed); a fortiori none")
    print("with deg G_i <= deg F and det Jac(G) a nonzero constant.")
    print("OBSTRUCTION CERTIFIED: no affine quotient for F3, F4, F5")


def write_certificates():
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "certificates.json")
    with open(path, "w") as fh:
        json.dump(CERT, fh, indent=1, default=str)
    print("[out] wrote", path)


if __name__ == "__main__":
    main()
