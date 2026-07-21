#!/usr/bin/env python3
"""
verify_charp_sieve.py
=====================
Exact SymPy verification backing `charp-sieve.md` and `sieve_dc1_candidate.py`:
the Tsuchimoto / BKK central-reduction bridge and the P2 tame cyclic obstruction,
mechanized as a screen for DC1 counterexample candidates.

Everything is exact arithmetic (integers mod p, symbolic resultants over Z then
reduced mod p). Bounded prime ranges are corroboration only; the load-bearing
facts (center theorem, Jacobian = 1 on each tested pair, geometric degree of each
tested map, the verdict logic) are checked as exact identities.

Sections:
  1. Weyl arithmetic soundness: d*x = x*d + 1, associativity, [D_std, X_std] = 1.
  2. Center theorem: a normal-ordered monomial is central mod p  <=>  p | i, p | j;
     and X^p, D^p ARE central for genuine [D,X]=1 pairs (the bridge is defined).
  3. The bridge on tame pairs: psi_p computed, Jac(psi_p) = 1 exactly, geom.deg 1.
     Exact regression on the small-prime (p=3) Artin-Schreier correction term.
  4. Broken pairs fail loudly: [D,X] != 1 rejected; a non-Weyl pair makes X^p / D^p
     NON-central, so induced_center_map raises BridgeFailure.
  5. Geometric degree engine: identity->1, (u^2,v)->2, (u^3,v)->3, wild
     Artin-Schreier (u, v - v^p)->p; the wild map has Jac = 1 (a Keller map that
     is NOT an automorphism -- the exact loophole P2 isolates).
  6. Tame-obstruction verdict logic: the full truth table
     (automorphism / possibly-wild / tame; the degree-2 EXCLUSION).
  7. Degree bound independent of p: deg_{u,v} psi_p <= bideg(X,D) across primes.
  8. End-to-end: tame candidates NOT excluded; a synthetic degree-2 tame-prime
     map IS excluded (exercising the exclusion path the sieve exists to fire).

Run:  uv run --with sympy python research/dc1-program/verify_charp_sieve.py
Ends: ALL CHARP SIEVE CHECKS PASSED
"""
import sympy as sp

import sieve_dc1_candidate as S
from sieve_dc1_candidate import (
    elt, w_mul, w_add, w_comm, w_pow, X_STD, D_STD,
    reduce_mod_p, is_central_mod_p, central_to_uv,
    induced_center_map, jacobian_mod_p, geometric_degree,
    tame_obstruction_verdict, sieve_candidate, BridgeFailure, CANDIDATES,
)

u, v = sp.symbols('u v')
FAILS = []


def ok(cond, label):
    if cond:
        print("PASS", label)
    else:
        print("FAIL", label)
        FAILS.append(label)


def raises(fn, exc, label):
    try:
        fn()
    except exc:
        print("PASS", label)
        return
    except Exception as e:  # noqa
        print("FAIL", label, "-- wrong exception:", repr(e))
        FAILS.append(label)
        return
    print("FAIL", label, "-- no exception raised")
    FAILS.append(label)


# =====================================================================
print("=" * 70)
print("1. Weyl arithmetic soundness")
print("=" * 70)
# d*x = x*d + 1  i.e.  (d)(x) - (x)(d) = 1
lhs = w_add(w_mul(D_STD, X_STD), S.w_scal(w_mul(X_STD, D_STD), -1))
ok(lhs == {(0, 0): 1}, "d*x - x*d = 1")
# [D_std, X_std] = 1
ok(w_comm(D_STD, X_STD) == {(0, 0): 1}, "[D_std, X_std] = 1")
# associativity spot check on x, d, and E=x*d
E_elt = elt((1, 1, 1))
A, B, C = X_STD, D_STD, E_elt
ok(w_mul(w_mul(A, B), C) == w_mul(A, w_mul(B, C)), "associativity (x*d)*E = x*(d*E)")
# d^2 * x^2 normal order:  d^2 x^2 = x^2 d^2 + 4 x d + 2
d2x2 = w_mul(elt((1, 0, 2)), elt((1, 2, 0)))
ok(d2x2 == {(2, 2): 1, (1, 1): 4, (0, 0): 2}, "d^2 x^2 = x^2 d^2 + 4 x d + 2")

# =====================================================================
print("\n" + "=" * 70)
print("2. Center theorem and bridge-well-definedness")
print("=" * 70)
for p in (3, 5, 7):
    # monomial central  <=>  p | i and p | j
    good = elt((1, p, 0), (1, 0, p), (1, p, p))
    ok(is_central_mod_p(reduce_mod_p(good, p), p), f"p={p}: x^p, d^p, x^p d^p are central")
    bad = elt((1, 1, 0))
    ok(not is_central_mod_p(reduce_mod_p(bad, p), p), f"p={p}: x is NOT central")
    bad2 = elt((1, p, 1))
    ok(not is_central_mod_p(reduce_mod_p(bad2, p), p), f"p={p}: x^p d is NOT central (p nmid j)")
# X^p, D^p central for genuine [D,X]=1 pairs (the BKK bridge is defined)
tame_pairs = [(k, CANDIDATES[k]) for k in
              ("identity", "shear d->d+x^2", "rotation", "band2 X=x+d^2", "composite")]
for name, (X, D) in tame_pairs:
    for p in (3, 5, 7, 11):
        Xp = w_pow(reduce_mod_p(X, p), p, reducer=lambda Adict: {k: c % p for k, c in Adict.items() if c % p})
        Dp = w_pow(reduce_mod_p(D, p), p, reducer=lambda Adict: {k: c % p for k, c in Adict.items() if c % p})
        ok(is_central_mod_p(Xp, p) and is_central_mod_p(Dp, p),
           f"[{name}] p={p}: X^p and D^p are central (bridge defined)")
# central_to_uv round-trips balanced representatives
for p in (5, 7):
    cen = reduce_mod_p(elt((1, p, 0), (-1, 0, p), (2, p, p)), p)
    poly = central_to_uv(cen, p)
    ok(sp.expand(poly - (u - v + 2 * u * v)) == 0, f"p={p}: central_to_uv gives u - v + 2 u v (balanced)")

# =====================================================================
print("\n" + "=" * 70)
print("3. The bridge on tame pairs: psi_p, Jac = 1, geom.deg = 1")
print("=" * 70)
for name, (X, D) in tame_pairs:
    for p in (3, 5, 7, 11):
        P, Q = induced_center_map(X, D, p)
        J = jacobian_mod_p(P, Q, p)
        ok(sp.expand(J - 1) == 0, f"[{name}] p={p}: Jac(psi_p) = 1 exactly")
        d, det = geometric_degree(P, Q, p)
        ok(d == 1 and det["consistent"], f"[{name}] p={p}: geometric degree = 1 (automorphism)")
# exact regression: the p=3 small-prime Artin-Schreier correction
P3, Q3 = induced_center_map(*CANDIDATES["shear d->d+x^2"], 3)
ok(sp.expand(P3 - u) == 0 and sp.expand(Q3 - (u**2 + v - 1)) == 0,
   "regression: shear d->d+x^2 gives (u, u^2+v-1) at p=3 (AS correction -1)")
P5, Q5 = induced_center_map(*CANDIDATES["shear d->d+x^2"], 5)
ok(sp.expand(Q5 - (u**2 + v)) == 0, "regression: shear d->d+x^2 gives (u, u^2+v) at p=5 (no correction)")

# =====================================================================
print("\n" + "=" * 70)
print("4. Broken pairs fail loudly")
print("=" * 70)
# [D,X] = 2 : rejected at the precondition gate
rep = sieve_candidate(*CANDIDATES["BROKEN [D,X]=2"], (5, 7), name="BROKEN2", verbose=False)
ok(rep["commutator_ok"] is False and "REJECTED" in rep["verdict"], "[D,X]=2 rejected (precondition)")
# [D,X] = 1 + x : rejected at the precondition gate
rep = sieve_candidate(*CANDIDATES["BROKEN [D,X]=1+x"], (5, 7), name="BROKEN1x", verbose=False)
ok(rep["commutator_ok"] is False and "REJECTED" in rep["verdict"], "[D,X]=1+x rejected (precondition)")
# a non-Weyl pair fed straight to the bridge: X^p or D^p is NON-central -> BridgeFailure
Xbad, Dbad = elt((1, 1, 0)), elt((1, 0, 1), (1, 1, 1))  # [D,X] = 1 + x
raises(lambda: induced_center_map(Xbad, Dbad, 5), BridgeFailure,
       "non-Weyl pair (X=x, D=d+xd): induced_center_map raises BridgeFailure at p=5")
raises(lambda: induced_center_map(Xbad, Dbad, 7), BridgeFailure,
       "non-Weyl pair (X=x, D=d+xd): induced_center_map raises BridgeFailure at p=7")
# confirm the failure is genuinely non-centrality of D^p
Dp = w_pow(reduce_mod_p(Dbad, 5), 5, reducer=lambda Adict: {k: c % 5 for k, c in Adict.items() if c % 5})
ok(not is_central_mod_p(Dp, 5), "the non-Weyl D^p is genuinely NON-central (mod 5)")

# =====================================================================
print("\n" + "=" * 70)
print("5. Geometric-degree engine + the wild Artin-Schreier loophole")
print("=" * 70)
ok(geometric_degree(u, v, 5)[0] == 1, "geom.deg(identity) = 1")
ok(geometric_degree(u**2, v, 5)[0] == 2, "geom.deg((u^2, v)) = 2")
ok(geometric_degree(u**3, v, 7)[0] == 3, "geom.deg((u^3, v)) = 3")
ok(geometric_degree(u + v**2, v, 5)[0] == 1, "geom.deg((u+v^2, v)) = 1 (total deg 2 but automorphism)")
for p in (3, 5, 7):
    d, _ = geometric_degree(u, v - v**p, p)
    ok(d == p, f"geom.deg(Artin-Schreier (u, v - v^{p})) = {p}")
    J = jacobian_mod_p(u, v - v**p, p)
    ok(sp.expand(J - 1) == 0, f"Artin-Schreier (u, v - v^{p}) has Jac = 1 (Keller, NOT an automorphism)")

# =====================================================================
print("\n" + "=" * 70)
print("6. Tame-obstruction verdict truth table")
print("=" * 70)
# automorphism: never excluded
for p in (3, 5, 7):
    ver = tame_obstruction_verdict(1, p)
    ok(ver["excluded"] is False and ver["regime"] == "automorphism", f"verdict(d=1, p={p}) = automorphism, not excluded")
# tame degree 2 with p > 2: EXCLUDED
for p in (3, 5, 7, 11):
    ver = tame_obstruction_verdict(2, p)
    ok(ver["excluded"] is True and ver["regime"] == "tame", f"verdict(d=2, p={p}) = EXCLUDED (tame Kummer)")
# degree 2 at p = 2: NOT excluded (p <= d, possibly wild)
ver = tame_obstruction_verdict(2, 2)
ok(ver["excluded"] is False and ver["regime"] == "possibly-wild", "verdict(d=2, p=2) = possibly-wild, NOT excluded")
# degree p at p (the Artin-Schreier degree): possibly-wild, not excluded
for p in (3, 5, 7):
    ver = tame_obstruction_verdict(p, p)
    ok(ver["excluded"] is False and ver["regime"] == "possibly-wild", f"verdict(d={p}, p={p}) = possibly-wild")
# tame degree 3 with p > 3: NOT excluded by the degree screen (S_3 monodromy allowed)
ver = tame_obstruction_verdict(3, 7)
ok(ver["excluded"] is False and ver["regime"] == "tame", "verdict(d=3, p=7) = tame, not excluded by degree screen")
# undefined degree
ver = tame_obstruction_verdict(None, 5)
ok(ver["excluded"] is False and ver["regime"] == "undefined", "verdict(d=None) = inconclusive")

# =====================================================================
print("\n" + "=" * 70)
print("7. Geometric degree of psi_p is bounded independently of p")
print("=" * 70)
def bideg(A):
    return max((i + j) for (i, j) in A) if A else 0
for name, (X, D) in tame_pairs:
    bound = max(bideg(X), bideg(D))
    for p in (3, 5, 7, 11, 13):
        P, Q = induced_center_map(X, D, p)
        du = sp.Poly(P, u, v).total_degree() if P != 0 else 0
        dv = sp.Poly(Q, u, v).total_degree() if Q != 0 else 0
        ok(max(du, dv) <= bound,
           f"[{name}] p={p}: deg_uv(psi_p) = {max(du, dv)} <= bideg(X,D) = {bound}")

# =====================================================================
print("\n" + "=" * 70)
print("8. End-to-end sieve verdicts")
print("=" * 70)
for name in ("identity", "shear d->d+x^2", "rotation", "band2 X=x+d^2", "composite", "half-shear"):
    rep = sieve_candidate(*CANDIDATES[name], (3, 5, 7, 11), name=name, verbose=False)
    ok(rep["excluded"] is False and rep["screen_degrees"] == [1],
       f"[{name}]: NOT excluded, reduces to automorphism at all sampled primes")
# half-shear must SKIP p=2 (denominator) and pass elsewhere
rep = sieve_candidate(*CANDIDATES["half-shear"], (2, 3, 5, 7), name="half-shear", verbose=False)
ok(rep["primes"][2]["status"] == "skipped", "half-shear: p=2 skipped (denominator)")

# The exclusion path fires on a synthetic tame-prime degree-2 map (the object the
# sieve exists to reject). No genuine [D,X]=1 pair yields this -- that IS P2 --
# so we exercise the verdict engine directly on a hand-supplied center map.
d_syn, _ = geometric_degree(u**2, v, 5)      # a genuine geometric-degree-2 map
ver = tame_obstruction_verdict(d_syn, 5)
ok(ver["excluded"] is True, "synthetic degree-2 center map at p=5: EXCLUDED (exclusion path fires)")

# a synthetic 'candidate' whose center map is degree 2 would be excluded end-to-end;
# demonstrate the aggregate flag by monkey-checking the aggregate rule directly:
fake_report_excluded = any(tame_obstruction_verdict(2, p)["excluded"] for p in (5, 7, 11))
ok(fake_report_excluded is True, "aggregate rule: a degree-2 map is excluded at tame primes 5,7,11")

# =====================================================================
print("\n" + "=" * 70)
if FAILS:
    print(f"{len(FAILS)} CHECK(S) FAILED:")
    for f in FAILS:
        print("  -", f)
    raise SystemExit(1)
print("ALL CHARP SIEVE CHECKS PASSED")
