#!/usr/bin/env python3
"""
sieve_dc1_candidate.py
======================
The CHAR-P SIEVE for DC1 counterexample candidates.

DC1 asks: is every endomorphism of the first Weyl algebra A_1 = C<x,d | [d,x]=1>
an automorphism?  An endomorphism is a pair (X, D) with [D, X] = 1 (it is
automatically injective since A_1 is simple; only surjectivity is in question).

This module implements a prototype central-reduction pipeline and a conditional
tame-obstruction candidate screen.  The geometric theorem needed to turn the
screen into an exclusion theorem is not established here:

  Given a pair (X, D) with rational (or integer) coefficients and [D, X] = 1,
  for each prime p not dividing a denominator:

    1. reduce mod p to a pair in A_1 tensor F_p;
    2. form the two central elements  X^p, D^p  in  Z = F_p[x^p, d^p];
       (that they are central is the BKK bridge; it is CHECKED exactly here and
        it FAILS for pairs that are not genuine Weyl pairs mod p);
    3. read them as a polynomial self-map  psi_p = (P, Q) : A^2_{F_p} -> A^2_{F_p}
       of the center  Spec Z = A^2,  with  u = x^p,  v = d^p;
    4. verify  Jac(psi_p) = 1  exactly  (the Keller / Jacobian-1 property);
    5. compute the geometric degree  d_p = [F_p(u,v) : F_p(P,Q)]  exactly, by
       elimination resultants;
    6. emit the tame-obstruction verdict.

If a future finite tame-extension theorem establishes all required geometric and
reduction hypotheses, a verified degree-2 center map at p > 2 would be an
OBSTRUCTION_CANDIDATE.  This module does not establish that theorem, so no output
is a theorem-level exclusion.

CONVENTIONS
-----------
A_1 tensor F_p is modeled by NORMAL-ORDERED polynomials: an element is a dict
    { (i, j) : coeff }   meaning   sum coeff * x^i * d^j   (x on the left).
The single relation d*x = x*d + 1 gives the normal-ordering product
    (x^i1 d^j1)(x^i2 d^j2) = sum_{k>=0} C(j1,k) C(i2,k) k!  x^(i1+i2-k) d^(j1+j2-k).
The coefficients C(j1,k) C(i2,k) k! are integers, so reduction mod p is exact.

The center of A_1 tensor F_p is  Z = F_p[x^p, d^p]  (proved in charp-sieve.md):
a normal-ordered monomial x^i d^j is central iff p | i and p | j.  Hence a
central element maps bijectively to a polynomial in u = x^p, v = d^p.

NOTHING here proves DC1.  The prototype only reports non-theorem obstruction
candidates or inconclusive outcomes; it never excludes a hypothetical
counterexample.  Bounded prime sweeps are corroboration only.

Run the self-demo:   uv run --with sympy python sieve_dc1_candidate.py
"""

from math import comb, factorial
import sympy as sp

u, v, a, b = sp.symbols('u v a b')


# =====================================================================
# 1. Exact Weyl-algebra arithmetic over Q (normal ordered), and mod p.
# =====================================================================
def w_add(A, B):
    out = dict(A)
    for k, c in B.items():
        out[k] = out.get(k, 0) + c
    return {k: c for k, c in out.items() if c != 0}


def w_scal(A, s):
    return {k: c * s for k, c in A.items() if c * s != 0}


def w_mul(A, B):
    """Normal-ordered product over an exact (Q) coefficient ring."""
    out = {}
    for (i1, j1), c1 in A.items():
        for (i2, j2), c2 in B.items():
            base = c1 * c2
            if base == 0:
                continue
            for k in range(0, min(j1, i2) + 1):
                coef = comb(j1, k) * comb(i2, k) * factorial(k)
                key = (i1 + i2 - k, j1 + j2 - k)
                out[key] = out.get(key, 0) + base * coef
    return {k: c for k, c in out.items() if c != 0}


def w_comm(A, B):
    return w_add(w_mul(A, B), w_scal(w_mul(B, A), -1))


def w_pow(A, n, reducer=None):
    """A^n; optional per-step reducer (e.g. mod p) to keep coefficients small."""
    R = {(0, 0): sp.Integer(1)}
    for _ in range(n):
        R = w_mul(R, A)
        if reducer is not None:
            R = reducer(R)
    return R


def elt(*terms):
    """Build an element from (coeff, i, j) triples. coeff may be int/Rational."""
    out = {}
    for c, i, j in terms:
        c = sp.nsimplify(c) if not isinstance(c, sp.Basic) else c
        out[(i, j)] = out.get((i, j), 0) + c
    return {k: c for k, c in out.items() if c != 0}


X_STD = elt((1, 1, 0))   # x
D_STD = elt((1, 0, 1))   # d


# ---- mod-p handling ----
def bad_primes(*elts):
    """Primes dividing any coefficient denominator: the bridge is invalid there."""
    dens = set()
    for A in elts:
        for c in A.values():
            r = sp.Rational(c)
            if r.q != 1:
                dens.add(int(r.q))
    ps = set()
    for d0 in dens:
        for f in sp.factorint(d0):
            ps.add(int(f))
    return ps


def reduce_mod_p(A, p):
    out = {}
    for k, c in A.items():
        r = sp.Rational(c)
        if r.q % p == 0:
            raise ValueError(f"coefficient {c} has denominator divisible by p={p}")
        val = (int(r.p) * pow(int(r.q), -1, p)) % p
        if val:
            out[k] = val
    return out


def is_central_mod_p(A, p):
    return all(i % p == 0 and j % p == 0 for (i, j) in A)


def central_to_uv(A, p):
    """Central mod-p element -> sympy polynomial in u=x^p, v=d^p (balanced ints)."""
    assert is_central_mod_p(A, p)
    e = sp.Integer(0)
    for (i, j), c in A.items():
        cc = c % p
        cc = cc if cc <= p // 2 else cc - p     # balanced representative
        e += cc * u**(i // p) * v**(j // p)
    return sp.expand(e)


# =====================================================================
# 2. The central-reduction bridge:  (X, D) mod p  |->  psi_p = (P, Q).
# =====================================================================
class BridgeFailure(Exception):
    pass


def induced_center_map(X, D, p):
    """
    Return (P, Q) with P = X^p, Q = D^p expressed in u=x^p, v=d^p over F_p.
    Raise BridgeFailure if X^p or D^p is not central (i.e. the bridge is
    undefined: the pair is not a genuine Weyl pair mod p).
    """
    Xp = w_pow(reduce_mod_p(X, p), p, reducer=lambda A: {k: c % p for k, c in A.items() if c % p})
    Dp = w_pow(reduce_mod_p(D, p), p, reducer=lambda A: {k: c % p for k, c in A.items() if c % p})
    if not is_central_mod_p(Xp, p):
        raise BridgeFailure(f"X^p not central mod p={p}: exponents {sorted(Xp)}")
    if not is_central_mod_p(Dp, p):
        raise BridgeFailure(f"D^p not central mod p={p}: exponents {sorted(Dp)}")
    return central_to_uv(Xp, p), central_to_uv(Dp, p)


# =====================================================================
# 3. Jacobian and geometric degree of a plane self-map over F_p.
# =====================================================================
def jacobian_mod_p(P, Q, p):
    """Return Jac(P,Q) = P_u Q_v - P_v Q_u  reduced mod p, as a balanced-int poly."""
    J = sp.expand(sp.diff(P, u) * sp.diff(Q, v) - sp.diff(P, v) * sp.diff(Q, u))
    return _reduce_expr_mod_p(J, p)


def _reduce_expr_mod_p(expr, p):
    expr = sp.expand(expr)
    if expr == 0:
        return sp.Integer(0)
    P = sp.Poly(expr, u, v)
    out = sp.Integer(0)
    for (i, j), c in P.terms():
        cc = int(c) % p
        cc = cc if cc <= p // 2 else cc - p
        if cc:
            out += cc * u**i * v**j
    return sp.expand(out)


def _u_degree_mod_p(res, keepvar, p):
    """res in ZZ[a,b,keepvar]; reduce coeffs mod p; degree in keepvar (or -oo/None)."""
    res = sp.expand(res)
    if res == 0:
        return None
    P = sp.Poly(res, keepvar)
    deg = -1
    for (e,), coeff in P.terms():
        cpoly = sp.Poly(coeff, a, b)
        reduced = sum((int(c) % p) * a**m[0] * b**m[1] for m, c in cpoly.terms())
        if sp.expand(reduced) != 0:
            deg = max(deg, e)
    return deg if deg >= 0 else None


def reconcile_elimination_degrees(du, dv):
    """Return a positive common degree, or None if missing, zero, or different."""
    return du if du is not None and dv is not None and du == dv and du > 0 else None


def geometric_degree(P, Q, p):
    """
    Geometric degree d = [F_p(u,v):F_p(P,Q)] via elimination resultants.

    Returns (d, detail) where detail records the two eliminations. A geometric
    degree is returned only when both elimination degrees are available and agree.
    Missing or inconsistent eliminations return d=None and require hand inspection;
    this routine never substitutes the maximum of inconsistent values.
    """
    Pa, Qb = sp.Poly(sp.expand(P) - a, v), sp.Poly(sp.expand(Q) - b, v)
    rv = sp.resultant(Pa, Qb)                        # eliminate v -> poly in u
    Pa2, Qb2 = sp.Poly(sp.expand(P) - a, u), sp.Poly(sp.expand(Q) - b, u)
    ru = sp.resultant(Pa2, Qb2)                      # eliminate u -> poly in v
    dv = _u_degree_mod_p(rv, u, p)
    du = _u_degree_mod_p(ru, v, p)
    available = (dv is not None and du is not None and dv > 0 and du > 0)
    d = reconcile_elimination_degrees(du, dv)
    consistent = d is not None
    return d, {"elim_v_deg_u": dv, "elim_u_deg_v": du,
               "available": available, "consistent": consistent}


# =====================================================================
# 4. The tame cyclic obstruction verdict  (the P2 screen).
# =====================================================================
def tame_obstruction_verdict(
        d, p, *, weyl_relation_verified=False,
        denominator_good_prime_verified=False,
        center_powers_central_verified=False, jac_is_one=None,
        degree_consistent=None, global_good_reduction_established=False,
        external_theorem_established=False):
    """Return a non-theorem verdict for the conditional tame-obstruction screen.

    Every local gate defaults to unverified.  A degree-2 result can be labeled an
    OBSTRUCTION_CANDIDATE only when the caller explicitly supplies all local
    evidence.  Global good reduction and the external theorem are recorded but are
    not established or consumed as theorem-level exclusion evidence by this code.
    """
    local_evidence = {
        "weyl_relation_verified": weyl_relation_verified is True,
        "denominator_good_prime_verified": denominator_good_prime_verified is True,
        "center_powers_central_verified": center_powers_central_verified is True,
        "jacobian_one_verified": jac_is_one is True,
        "consistent_elimination_degree_verified": degree_consistent is True,
    }
    flags = {
        "d": d,
        "p": p,
        "excluded": False,
        "obstruction_candidate": False,
        "local_evidence": local_evidence,
        "global_good_reduction_established": global_good_reduction_established,
        "external_theorem_established": external_theorem_established,
    }
    if jac_is_one is False:
        flags["regime"] = "non-keller"
        flags["verdict"] = "INCONCLUSIVE: induced center map does not have Jacobian exactly one."
        return flags
    missing = [name for name, verified in local_evidence.items() if not verified]
    if missing:
        flags["regime"] = "unverified-local-evidence"
        flags["missing_local_evidence"] = missing
        flags["verdict"] = ("INCONCLUSIVE: obstruction screening requires explicit verified "
                            "local evidence for " + ", ".join(missing) + ".")
        return flags
    if d is None:
        flags["regime"] = "undefined"
        flags["verdict"] = ("INCONCLUSIVE: both elimination degrees must be available "
                            "and equal before assigning a geometric degree.")
        return flags
    if d == 1:
        flags["regime"] = "automorphism"
        flags["verdict"] = ("INCONCLUSIVE FOR DC1: the sampled center map has degree 1; "
                            "this is consistent with automorphy but does not prove it.")
        return flags
    if p <= d:
        flags["regime"] = "possibly-wild"
        flags["verdict"] = (f"INCONCLUSIVE: degree {d} is not prime to p={p}; an "
                            "Artin-Schreier-type wild quotient is not ruled out.")
        return flags
    flags["regime"] = "tame-candidate"
    if d == 2:
        flags["obstruction_candidate"] = True
        flags["verdict"] = (f"OBSTRUCTION_CANDIDATE (non-theorem): degree 2 at p={p}>2 with "
                            "all local prerequisites verified. Global good reduction and "
                            "the external finite tame-extension theorem remain unestablished.")
    else:
        flags["verdict"] = (f"INCONCLUSIVE: degree {d} at p={p}>d is a conditional tame "
                            "regime; full monodromy and the required global geometry are "
                            "unavailable.")
    return flags


# =====================================================================
# 5. Whole-candidate sieve.
# =====================================================================
def aggregate_obstruction_verdicts(verdicts):
    """Aggregate local verdicts without upgrading mixed evidence.

    A candidate label requires at least one local candidate and no sampled local
    failure or missing prerequisite. Any mixed set remains inconclusive.
    """
    verdicts = list(verdicts)
    required_local_evidence = {
        "weyl_relation_verified",
        "denominator_good_prime_verified",
        "center_powers_central_verified",
        "jacobian_one_verified",
        "consistent_elimination_degree_verified",
    }
    candidate_flags = [
        verdict.get("obstruction_candidate", False) is True
        and verdict.get("excluded", False) is False
        for verdict in verdicts
    ]
    all_local_evidence_verified = [
        set(verdict.get("local_evidence", {})) == required_local_evidence
        and all(verdict["local_evidence"][key] is True
                for key in required_local_evidence)
        for verdict in verdicts
    ]
    no_local_failure = all(
        verified and verdict.get("regime") not in {
            "non-keller", "unverified-local-evidence", "undefined"
        }
        for verdict, verified in zip(verdicts, all_local_evidence_verified)
    )
    return {
        "excluded": False,
        "obstruction_candidate": bool(verdicts) and any(candidate_flags) and no_local_failure,
        "mixed_or_incomplete_local_evidence": bool(verdicts) and not no_local_failure,
    }


def aggregate_prime_entries(prime_entries):
    """Aggregate requested-prime entries, failing closed on attempted reductions.

    Denominator primes are outside the good-reduction domain and are skipped. Any
    other entry without a complete verified verdict counts as incomplete evidence.
    """
    verdicts = []
    for entry in prime_entries:
        status = entry.get("status")
        if status == "skipped":
            continue
        verdicts.append(entry.get("verdict", {}) if status == "ok" else {})
    return aggregate_obstruction_verdicts(verdicts)


def sieve_candidate(X, D, primes, name="candidate", verbose=True):
    """
    Screen a pair (X, D) [dicts over Q with [D,X]=1] against P2 over `primes`.
    Returns a report dict. Raises ValueError if [D,X] != 1 over Q.
    """
    comm = w_comm(D, X)
    comm_ok = (comm == {(0, 0): sp.Integer(1)})
    report = {"name": name, "commutator": comm, "commutator_ok": comm_ok, "primes": {}}
    if not comm_ok:
        report["verdict"] = "REJECTED PRECONDITION: [D,X] != 1 over Q (not an endomorphism)."
        if verbose:
            _print_header(report)
            print(f"  [D,X] = {_fmt(comm)}   !=  1     ->  {report['verdict']}")
        return report

    bad = bad_primes(X, D)
    for p in primes:
        entry = {}
        if p in bad:
            entry["status"] = "skipped"
            entry["reason"] = "p divides a coefficient denominator"
            report["primes"][p] = entry
            continue
        try:
            P, Q = induced_center_map(X, D, p)
        except BridgeFailure as e:
            entry["status"] = "bridge-failure"
            entry["reason"] = str(e)
            report["primes"][p] = entry
            continue
        J = jacobian_mod_p(P, Q, p)
        jac_is_one = (sp.expand(J - 1) == 0)
        d, ddetail = geometric_degree(P, Q, p)
        verdict = tame_obstruction_verdict(
            d, p,
            weyl_relation_verified=comm_ok,
            denominator_good_prime_verified=(p not in bad),
            center_powers_central_verified=True,
            jac_is_one=jac_is_one,
            degree_consistent=ddetail["consistent"],
            global_good_reduction_established=False,
            external_theorem_established=False,
        )
        entry.update({
            "status": "ok",
            "P": P, "Q": Q, "jacobian": J, "jac_is_one": jac_is_one,
            "geom_degree": d, "degree_detail": ddetail,
            "local_checks": verdict["local_evidence"],
            "global_good_reduction_established": False,
            "external_theorem_established": False,
            "verdict": verdict,
        })
        report["primes"][p] = entry

    report["screen_degrees"] = sorted({e.get("geom_degree") for e in report["primes"].values()
                                       if e.get("status") == "ok"} - {None})
    aggregate = aggregate_prime_entries(report["primes"].values())
    report.update(aggregate)
    degs = report["screen_degrees"]
    if report["obstruction_candidate"]:
        report["verdict"] = ("OBSTRUCTION_CANDIDATE (non-theorem) at a sampled prime. "
                             "The external finite tame-extension theorem and its geometric "
                             "and reduction hypotheses are not established here.")
    elif degs == [1]:
        report["verdict"] = ("INCONCLUSIVE FOR DC1: every valid sampled center map has "
                             "geometric degree 1. This is consistent with automorphy but "
                             "does not prove it.")
    else:
        report["verdict"] = (f"INCONCLUSIVE: sampled geometric degrees {degs}; no "
                             "theorem-level exclusion is implemented.")
    if verbose:
        _print_report(report)
    return report


# ---- pretty printing ----
def _fmt(A):
    if not A:
        return "0"
    parts = []
    for (i, j), c in sorted(A.items()):
        mon = "1" if (i == 0 and j == 0) else f"x^{i} d^{j}"
        parts.append(f"({c})*{mon}")
    return " + ".join(parts)


def _print_header(report):
    print("=" * 72)
    print(f"CANDIDATE: {report['name']}")
    print("=" * 72)


def _print_report(report):
    _print_header(report)
    print(f"  [D,X] = {_fmt(report['commutator'])}   (=1 ? {report['commutator_ok']})")
    for p, e in sorted(report["primes"].items()):
        if e["status"] == "skipped":
            print(f"  p={p:>3}: SKIPPED ({e['reason']})")
        elif e["status"] == "bridge-failure":
            print(f"  p={p:>3}: BRIDGE FAILURE -- {e['reason']}")
        else:
            v = e["verdict"]
            print(f"  p={p:>3}: psi_p = ( u -> {e['P']},   v -> {e['Q']} )")
            print(f"        Jac = {e['jacobian']}  (=1 ? {e['jac_is_one']});  "
                  f"geom.deg = {e['geom_degree']}  [{v['regime']}]")
            print(f"        local checks = {e['local_checks']}")
            print("        global good reduction established = "
                  f"{e['global_good_reduction_established']}")
            print(f"        -> {v['verdict']}")
    print("-" * 72)
    print(f"  SCREEN VERDICT: {report['verdict']}")
    print()


# =====================================================================
# 6. A small library of named candidates for the self-demo.
# =====================================================================
CANDIDATES = {
    # tame automorphism images (must all pass: geom degree 1, Jac 1)
    "identity":        (X_STD, D_STD),
    "shear d->d+x^2":  (X_STD, elt((1, 0, 1), (1, 2, 0))),
    "rotation":        (elt((1, 0, 1)), elt((-1, 1, 0))),
    "band2 X=x+d^2":   (elt((1, 1, 0), (1, 0, 2)), elt((1, 0, 1))),
    "composite":       (elt((1, 0, 1)), w_add(w_pow(elt((1, 0, 1)), 2), elt((-1, 1, 0)))),
    # a rational-coefficient tame pair (exercises bad-prime skipping at p=2)
    "half-shear":      (X_STD, elt((1, 0, 1), (sp.Rational(1, 2), 2, 0))),
    # deliberately broken pairs (must fail loudly)
    "BROKEN [D,X]=2":  (X_STD, elt((2, 0, 1))),
    "BROKEN [D,X]=1+x": (X_STD, elt((1, 0, 1), (1, 1, 1))),
}


def demo():
    primes = (3, 5, 7, 11)
    for name, (X, D) in CANDIDATES.items():
        sieve_candidate(X, D, primes, name=name)


if __name__ == "__main__":
    demo()
