#!/usr/bin/env python3
"""Exact-algebra certificate for the W2 COMBINED-SYSTEM VERDICT.

Question (the last gate before a DC1 counterexample candidate at band 3):
does ANY slope-1 datum simultaneously satisfy, at a common raw-degree cap d,

    positive cascade  Q_4=Q_3=Q_2=Q_1=0
    slope / moment unit  Q_0 = 1        (<=> R(1)=1, R(-1)=-1)
    full negative tail  Q_-1=Q_-2=Q_-3=Q_-4=Q_-5 = 0
    membership  (E)_j = E(E-1)...(E-j+1) | a_-j, b_-j   (j=1,2,3)

for the exceptional exotic top  a_3=E(E+2)(E+4), b_2=E(E+3), b_3=0  (r=-4, "W2")?

Both tail branches are covered (they are exhaustive: Q_-6=0 forces the ratio
b_-3/a_-3 to be (-3)-periodic hence constant when a_-3 != 0):
    branch A:  a_-3 = (E)_3 * am3_raw (generic),  b_-3 = mu3 * a_-3
    branch B:  a_-3 = 0,                          b_-3 = (E)_3 * bm3_raw (free filler)

SELF-CONTAINED VERDICT (this file, exact over QQ): INFEASIBLE at d=3, BOTH
branches -- exact identities checked directly against the freshly reconstructed
FULL systems prove that both ideals are UNIT. Optional msolve runs can additionally
reproduce the reported d=4 exact-QQ and mod-65003 results. If msolve is absent
those checks are explicit SKIP, not successes.

SCOPE (honest): this proves no encoded solution at raw cap d=3 in the fixed W2
orientation/gauge.  A successful optional msolve run also checks d=4.  It does
NOT prove W2 dead at arbitrary degree, and it does not execute the documentary
five-prime d=3/d=4 sweeps or d=5 mod-65003 report.  See w2-verdict.md.

Method.  sympy alone cannot Groebner the raw 62-eq/26-var d=3 system in a sane
budget.  The committed companion witness gives sparse multipliers h_i for the
freshly reconstructed FULL equations f_i.  This verifier parses all coefficients
in QQ and checks the polynomial identity sum_i h_i*f_i = 1 by exact coefficient
collection; that independently proves the ideal is the whole QQ polynomial ring.
The identity was generated with Singular, but Singular is neither trusted nor
needed here.  A feasibility-preserving linear reducer remains only for its d=2
control.  Separate msolve 0.10.1 computations were reported for d=3/d=4 over QQ
and five primes; no outputs are committed, so those results remain documentary
unless an optional command is run here.

Controls (no false kill):
  - The explicit w2-decisive Section-4 slope-1 datum satisfies the cascade, Q_0=1
    and all six memberships but FAILS the tail (all five Q_-m != 0); its raw free
    degrees fit cap d for every d>=3, witnessing cascade+slope FEASIBLE in branch B.
  - d<=2: cascade+slope UNIT and FULL UNIT in both branches; cascade+tail is checked
    PROPER at d=1 in both branches.  The sibling tail verifier separately checks
    cascade+tail at d=2.

Run:  uv run --with sympy python research/band3/verify_w2_verdict.py
The final summary lists every executed PASS and absent optional-tool SKIP.
"""
import sympy as sp
import time
from pathlib import Path

E = sp.symbols("E")
LEVELS = range(-3, 4)
_T0 = time.time()
_STATUSES = []


def sh(f, n):
    return sp.expand(sp.sympify(f).subs(E, E + n))


def falling(n):
    return sp.prod(E - j for j in range(n))


def poly(name, degree):
    cs = list(sp.symbols(f"{name}_0:{degree + 1}"))
    return sp.expand(sum(cs[j] * E**j for j in range(degree + 1))), cs


def q_m(A, B, m):
    """Stipulated convention Q_m = sum_(k+l=m) (b_l^[k] a_k - a_k^[l] b_l)."""
    return sp.expand(sum(sh(B[l], k) * A[k] - sh(A[k], l) * B[l]
                         for k in LEVELS for l in LEVELS if k + l == m))


def potential(A, B):
    return sp.expand(sum(sh(A[k], j - k) * sh(B[-k], j) - sh(B[k], j - k) * sh(A[-k], j)
                         for k in range(1, 4) for j in range(k)))


def record(status, label):
    _STATUSES.append((status, label))
    print(f"{status} [{time.time() - _T0:6.1f}s] {label}")


def check(condition, label):
    if not condition:
        record("FAIL", label)
        raise AssertionError("FAIL " + label)
    record("PASS", label)


def check_zero(value, label):
    if sp.expand(value) != 0:
        record("FAIL", label)
        raise AssertionError(f"FAIL {label}: residual {sp.factor(sp.expand(value))}")
    record("PASS", label)


def skip(label):
    record("SKIP", label)


# ---- crossed-product engine (identical to the sibling tail verifier) ----
def mul_ladders(P, Q):
    R = {}
    for k1, p1 in P.items():
        for k2, p2 in Q.items():
            R[k1 + k2] = R.get(k1 + k2, 0) + sh(p1, k2) * p2
    return {k: sp.expand(v) for k, v in R.items() if sp.expand(v) != 0}


def commutator_coeff(X, Dop, m):
    return sp.expand(mul_ladders(Dop, X).get(m, 0) - mul_ladders(X, Dop).get(m, 0))


def clean_solve(A, B, m, lkey, name, membership, raw_degree):
    """Solve Q_m = 0 linearly for B[lkey] = (E)_membership * (free poly).  Returns
    the solved coefficient (plus kernel-parameter freedom) and the compatibility
    conditions on the remaining free data (the positive-cascade constraints)."""
    raw, cs = poly(name, raw_degree)
    unknown = sp.expand(falling(membership) * raw)
    trial = dict(B)
    trial[lkey] = unknown
    equations = sp.Poly(q_m(A, trial, m), E).all_coeffs()
    M, rhs = sp.linear_eq_to_matrix(equations, cs)
    if any(entry.free_symbols for entry in M):
        raise AssertionError("operator matrix is not numeric; bilinearity leaked")
    conditions = [c for c in (sp.expand(n.dot(rhs)) for n in M.T.nullspace()) if c != 0]
    independent = sp.zeros(0, len(cs))
    selected_rhs = []
    for i in range(M.rows):
        candidate = independent.col_join(M[i, :])
        if candidate.rank() > independent.rank():
            independent = candidate
            selected_rhs.append(rhs[i])
    if independent.rows == 0:
        values = [sp.Integer(0)] * len(cs)
    else:
        solution, _ = independent.gauss_jordan_solve(sp.Matrix(selected_rhs))
        taus = [x for x in solution.free_symbols if str(x).startswith("tau")]
        values = [x.subs({t: 0 for t in taus}) for x in solution]
    result = sp.expand(unknown.subs(dict(zip(cs, values))))
    kernels = []
    for j, vector in enumerate(M.nullspace()):
        parameter = sp.symbols(f"{name}K{j}")
        kernel_poly = falling(membership) * sum(vector[i] * E**i for i in range(len(cs)))
        result = sp.expand(result + parameter * kernel_poly)
        kernels.append(parameter)
    return result, kernels, conditions


def positive_cascade(a3, b2, d):
    """Solve Q_4..Q_1 for b_1,b_0,b_-1,b_-2 (independent of a_-3,b_-3)."""
    a2, ca2 = poly("a2", d)
    a1, ca1 = poly("a1", d)
    a0, ca0 = poly("a0", d)
    am1_raw, cam1 = poly("am1", d)
    am2_raw, cam2 = poly("am2", d)
    A = {3: a3, 2: a2, 1: a1, 0: a0, -1: falling(1) * am1_raw, -2: falling(2) * am2_raw}
    B = {k: sp.Integer(0) for k in range(-3, 4)}
    B[2] = b2
    A[-3] = sp.Integer(0)
    B[-3] = sp.Integer(0)
    conditions, kernels = [], []
    for m, lkey, name, membership, degree in [
            (4, 1, "b1c", 0, d + 3), (3, 0, "b0c", 0, 2 * d + 2),
            (2, -1, "bm1c", 1, 2 * d + 3), (1, -2, "bm2c", 2, 2 * d + 4)]:
        B[lkey], nk, nc = clean_solve(A, B, m, lkey, name, membership, degree)
        kernels += nk
        conditions += nc
    return A, B, conditions, ca2 + ca1 + ca0 + cam1 + cam2, kernels


def build_branch(a3, b2, d, branch):
    A, B, pos, base_vars, kernels = positive_cascade(a3, b2, d)
    if branch == "A":
        am3_raw, cam3 = poly("am3", d)
        mu3 = sp.symbols("mu3")
        A[-3] = sp.expand(falling(3) * am3_raw)
        B[-3] = sp.expand(mu3 * A[-3])
        extra = cam3 + [mu3]
    elif branch == "B":
        bm3_raw, cbm3 = poly("bm3", d)
        A[-3] = sp.Integer(0)
        B[-3] = sp.expand(falling(3) * bm3_raw)
        extra = cbm3
    else:
        raise ValueError(branch)
    return A, B, pos, base_vars + extra + kernels


def scalar_coeffs(expr):
    return [c for c in sp.Poly(sp.expand(expr), E).all_coeffs() if sp.expand(c) != 0]


def systems(a3, b2, d, branch):
    A, B, pos, fv = build_branch(a3, b2, d, branch)
    posnz = [e for e in pos if sp.expand(e) != 0]
    slope = scalar_coeffs(q_m(A, B, 0) - 1)
    tail = []
    for m in (-1, -2, -3, -4, -5):
        tail += scalar_coeffs(q_m(A, B, m))
    return {"cascade": posnz, "cascade+tail": posnz + tail,
            "cascade+slope": posnz + slope, "FULL": posnz + slope + tail}, fv


def clear_denoms(e):
    num, _ = sp.fraction(sp.together(sp.expand(e)))
    return sp.expand(num)


def is_unit(eqs, vs, domain=sp.QQ):
    eqs = [clear_denoms(e) for e in eqs if sp.expand(e) != 0]
    if not eqs:
        return False
    G = sp.groebner(eqs, *vs, order="grevlex", domain=domain)
    return list(G.exprs) == [sp.Integer(1)]


def linear_reduce(eqs, fv):
    """Feasibility-preserving elimination: while some equation is c*v+rest with c a
    nonzero rational constant and rest free of v, set v=-rest/c and substitute.
    Preserves V(.)=empty exactly.  Returns (residual eqs, residual vars)."""
    eqs = [sp.expand(e) for e in eqs if sp.expand(e) != 0]
    vs = list(fv)
    changed = True
    while changed:
        changed = False
        best = None
        for e in eqs:
            fsyms = e.free_symbols & set(vs)
            if not fsyms:
                continue
            # fv order is intentional: set iteration made this reducer nondeterministic.
            for v in (v_ for v_ in vs if v_ in fsyms):
                if sp.Poly(e, v).degree() != 1:
                    continue
                coeff = e.coeff(v, 1)
                rest = sp.expand(e - coeff * v)
                if coeff == 0 or coeff.free_symbols or v in rest.free_symbols:
                    continue
                score = len(sp.Add.make_args(rest))
                candidate = (score, vs.index(v), e, v, sp.expand(-rest / coeff))
                if best is None or candidate[:2] < best[:2]:
                    best = candidate
        if best is not None:
            _, _, e0, v0, sol = best
            eqs = [sp.expand(ee.subs(v0, sol)) for ee in eqs if ee is not e0]
            eqs = [ee for ee in eqs if sp.expand(ee) != 0]
            vs.remove(v0)
            changed = True
    return eqs, vs


def load_qq_certificate(branch, variables):
    """Load and validate the sparse multipliers for sum_i h_i f_i = 1."""
    if branch not in {"A", "B"}:
        raise AssertionError(f"unknown certificate branch {branch!r}")
    path = Path(__file__).with_name("w2_d3_qq_certificates.txt")
    sections = {"A": {}, "B": {}}
    section = None
    seen_headers = set()
    for line_number, raw in enumerate(path.read_text(encoding="ascii").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") or line.endswith("]"):
            if not (line.startswith("[") and line.endswith("]")):
                raise AssertionError(f"malformed certificate section at line {line_number}")
            section = line[1:-1]
            if section not in sections:
                raise AssertionError(f"unknown certificate section {section!r}")
            if section in seen_headers:
                raise AssertionError(f"duplicate certificate section {section!r}")
            seen_headers.add(section)
            continue
        if section is None or ":" not in line:
            raise AssertionError(f"malformed certificate entry at line {line_number}")
        index_text, multiplier = line.split(":", 1)
        try:
            index = int(index_text)
        except ValueError as exc:
            raise AssertionError(
                f"non-integer certificate index at line {line_number}"
            ) from exc
        if index < 0:
            raise AssertionError(f"negative certificate index at line {line_number}")
        if index in sections[section]:
            raise AssertionError(
                f"duplicate equation index {index} in certificate section {section}"
            )
        if not multiplier.strip():
            raise AssertionError(f"empty certificate multiplier at line {line_number}")
        sections[section][index] = sp.Poly(
            multiplier,
            *variables,
            domain=sp.QQ,
        ).as_expr()
    if seen_headers != set(sections):
        missing = sorted(set(sections) - seen_headers)
        raise AssertionError(f"missing certificate section(s): {missing}")
    if not sections[branch]:
        raise AssertionError(f"empty d=3 QQ certificate for branch {branch}")
    return sections[branch]


def verify_qq_unit_certificate(eqs, fv, branch):
    """Check a Nullstellensatz identity exactly by sparse QQ coefficient collection."""
    eqs = [sp.expand(e) for e in eqs if sp.expand(e) != 0]
    expected_vars = {"A": 26, "B": 25}
    if len(eqs) != 62 or len(fv) != expected_vars[branch]:
        raise AssertionError(
            f"unexpected d=3 branch-{branch} shape: {len(eqs)} equations, "
            f"{len(fv)} variables"
        )
    ys = sp.symbols(f"y0:{len(fv)}")
    substitution = dict(zip(fv, ys))
    multipliers = load_qq_certificate(branch, ys)
    if max(multipliers) >= len(eqs):
        raise AssertionError("certificate equation index is out of range")
    identity = sum(h * eqs[i].subs(substitution) for i, h in multipliers.items()) - 1
    return sp.Poly(identity, *ys, domain=sp.QQ).is_zero, len(eqs), len(ys), len(multipliers)


a3_w2 = sp.expand(E * (E + 2) * (E + 4))
b2_w2 = sp.expand(E * (E + 3))
W2 = (a3_w2, b2_w2)

# =====================================================================
print("--- 0. machinery: Q_m equals the crossed-product commutator, Q_0=(T-1)G ---")
# =====================================================================
Ag, Bg = {}, {}
for lev in LEVELS:
    ar, _ = poly(f"Ag{lev + 3}", 2)
    br, _ = poly(f"Bg{lev + 3}", 2)
    memb = falling(-lev) if lev < 0 else 1
    Ag[lev] = sp.expand(memb * ar)
    Bg[lev] = sp.expand(memb * br)
for m in range(-6, 7):
    check_zero(q_m(Ag, Bg, m) - commutator_coeff(Ag, Bg, m),
               f"Q_{m} = [D,X]_{m} (generic degree-2 coefficients)")
check_zero(q_m(Ag, Bg, 0) - (sh(potential(Ag, Bg), 1) - potential(Ag, Bg)),
           "Q_0 = (T-1)G telescoping identity")

# =====================================================================
print("\n--- 1. explicit witness: the w2-decisive Section-4 slope-1 datum ---")
# =====================================================================
R = sp.Rational
Ax = {
    3: a3_w2,
    2: sp.expand(-E**3 / 3 + E**2 / 2 + R(9, 2) * E + 2),
    1: sp.expand(R(5, 7) * E**3 + R(107, 63) * E**2 + R(118, 63) * E + R(10, 9)),
    0: sp.expand(-R(775, 5103) * E**3 - R(92545, 3402) * E**2 - R(277597, 10206) * E),
    -1: sp.expand(-R(9, 8) * E**4 - R(219830, 11907) * E),
    -2: sp.expand(R(1, 5) * E**5 - R(1, 5) * E**4 + R(967027, 2755620) * E**3
                  - R(53597707, 8573040) * E**2 + R(455302607, 77157360) * E),
    -3: sp.Integer(0),
}
Bx = {
    3: sp.Integer(0),
    2: b2_w2,
    1: sp.expand(-R(2, 9) * E**2 + R(8, 9) * E + R(4, 3)),
    0: sp.expand(R(263, 567) * E**2 + R(179, 567) * E),
    -1: sp.expand(-R(256, 5103) * E**2 - R(31130, 1701) * E),
    -2: sp.expand(-R(3, 4) * E**3 + R(2747993, 2571912) * E**2 - R(819059, 2571912) * E),
    -3: sp.expand(R(2, 15) * E**4 - R(5, 12) * E**3 + R(19, 60) * E**2 - R(1, 30) * E),
}
for m in (4, 3, 2, 1):
    check_zero(q_m(Ax, Bx, m), f"witness: positive cascade Q_{m} = 0")
check_zero(q_m(Ax, Bx, 0) - 1, "witness: Q_0 = 1 (slope R(1)=1 -- slope IS achievable)")
for j in (1, 2, 3):
    check_zero(sp.rem(Ax[-j], falling(j), E), f"witness: membership (E)_{j} | a_-{j}")
    check_zero(sp.rem(Bx[-j], falling(j), E), f"witness: membership (E)_{j} | b_-{j}")
tail_witness = [q_m(Ax, Bx, m) for m in (-1, -2, -3, -4, -5)]
check(all(t != 0 for t in tail_witness),
      "witness: the tail FAILS at this slope-1 point (all Q_-1..Q_-5 != 0)")
# The raw free quotient polynomials fit cap d=3, so this same point embeds in
# every larger raw cap; cascade+slope is therefore feasible for all d>=3.
degfit = (sp.Poly(Ax[2], E).degree() <= 3 and sp.Poly(Ax[-1], E).degree() <= 4
          and sp.Poly(Ax[-2], E).degree() <= 5)
check(degfit, "witness: coeff degrees fit cap d=3 (hence all d>=3): cascade+slope feasible d>=3")

# =====================================================================
print("\n--- 2. d<=2 controls: slope kills alone; tail is separately feasible ---")
# =====================================================================
for d in (1, 2):
    for br in ("A", "B"):
        sysd, fv = systems(*W2, d, br)
        check(is_unit(sysd["cascade+slope"], fv),
              f"d={d} br {br}: cascade+Q_0=1 = UNIT (slope alone kills at d<=2)")
        check(is_unit(sysd["FULL"], fv), f"d={d} br {br}: FULL = UNIT (infeasible)")
# tail feasible WITHOUT the slope (no false kill); checked at d=1 (positive-dim
# Groebner is instant there; the sibling verify_w2_tail.py checks it at d=1,2).
for br in ("A", "B"):
    sysd, fv = systems(*W2, 1, br)
    check(not is_unit(sysd["cascade+tail"], fv),
          f"d=1 br {br}: cascade+tail = PROPER (tail feasible without the slope)")

# =====================================================================
print("\n--- 3. reducer validation at d=2: reduced-unit == direct-unit ---")
# =====================================================================
for br in ("A", "B"):
    sysd, fv = systems(*W2, 2, br)
    red, vs = linear_reduce(sysd["FULL"], fv)
    direct = is_unit(sysd["FULL"], fv)
    reduced = is_unit(red, vs)
    check(direct == reduced and reduced,
          f"d=2 br {br}: linear_reduce preserves the UNIT verdict (direct={direct}, reduced={reduced})")

# =====================================================================
print("\n--- 4. THE d=3 VERDICT: FULL system is the UNIT IDEAL, both branches (exact QQ) ---")
# =====================================================================
for br in ("A", "B"):
    sysd, fv = systems(*W2, 3, br)
    valid, n_eqs, n_vars, n_nonzero = verify_qq_unit_certificate(sysd["FULL"], fv, br)
    print(f"    d=3 br {br}: exact identity uses {n_nonzero}/{n_eqs} multipliers"
          f" in QQ[y0,...,y{n_vars - 1}]")
    check(valid,
          f"d=3 br {br}: committed QQ identity proves FULL = UNIT  =>  INFEASIBLE")

# =====================================================================
print("\n--- 5. optional d=4 confirmation of the reported full-system verdict ---")
# =====================================================================
# The raw d=4 system (73 eqs / ~32 vars) is beyond a quick SymPy Groebner.
# Reproduce the reported exact-QQ result and one modular check opportunistically
# when msolve is on PATH; otherwise record explicit SKIPs. The reported five-prime
# sweep is documentary and is not executed here; see w2-verdict.md.
import shutil
import subprocess
import tempfile
import os


def msolve_verdict(eqs, fv, char):
    """Return UNIT only for msolve's complete [-1] output; errors propagate."""
    xs = sp.symbols(f"y0:{len(fv)}")
    sub = dict(zip(fv, xs))
    with tempfile.NamedTemporaryFile("w", suffix=".ms", delete=False) as f:
        inp = f.name
        f.write(",".join(str(x) for x in xs) + "\n" + str(char) + "\n")
        f.write(",\n".join(sp.sstr(sp.expand(e.subs(sub))).replace("**", "^").replace(" ", "")
                           for e in eqs if sp.expand(e) != 0) + "\n")
    out = inp + ".out"
    try:
        proc = subprocess.run(["msolve", "-f", inp, "-o", out], check=True,
                              capture_output=True, text=True, timeout=1200)
        if not os.path.exists(out):
            raise RuntimeError("msolve completed without creating its output file")
        with open(out, encoding="utf-8") as f:
            result = f.read().strip()
        if result == "[-1]":
            return "UNIT"
        raise RuntimeError(
            "msolve returned an unrecognized/non-unit result; "
            f"stdout={proc.stdout[-500:]!r}, stderr={proc.stderr[-500:]!r}, "
            f"output={result[:500]!r}"
        )
    finally:
        for p_ in (inp, out):
            if os.path.exists(p_):
                os.remove(p_)


if shutil.which("msolve"):
    for br in ("A", "B"):
        sysd, fv = systems(*W2, 4, br)
        check(msolve_verdict(sysd["FULL"], fv, 0) == "UNIT",
              f"d=4 br {br}: FULL = UNIT over QQ via msolve  =>  INFEASIBLE (exact)")
        check(msolve_verdict(sysd["FULL"], fv, 65003) == "UNIT",
              f"d=4 br {br}: FULL = UNIT mod 65003 via msolve")
else:
    skip("d=4 branches A/B: exact QQ via optional msolve (msolve not on PATH)")
    skip("d=4 branches A/B: mod 65003 via optional msolve (msolve not on PATH)")

# These computations are reported in the memo but are not encoded by this verifier.
skip("documentary five-prime d=3 sweep (no committed inputs/outputs)")
skip("documentary five-prime d=4 sweep (no committed inputs/outputs)")
skip("documentary d=5 mod 65003 computation (no committed inputs/outputs)")

passed = [label for status, label in _STATUSES if status == "PASS"]
skipped = [label for status, label in _STATUSES if status == "SKIP"]
print("\n" + "=" * 70)
print("EXECUTED RESULT: exact SymPy checks passed, including FULL d=3 over QQ")
print("on both tail branches. No encoded solution exists at raw cap d=3 in")
print("the fixed normalized W2 orientation/gauge. Arbitrary degree remains OPEN.")
print(f"STATUS: {len(passed)} PASS, {len(skipped)} SKIP; any failure exits immediately")
for label in skipped:
    print(f"  SKIP: {label}")
print("=" * 70)
print(f"\n(total {time.time() - _T0:.1f}s)")
print("W2 VERDICT EXECUTED CHECKS PASSED; OPTIONAL/DOCUMENTARY ITEMS SKIPPED ABOVE")
