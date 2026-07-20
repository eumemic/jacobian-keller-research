# Affine-quotient obstruction for the explicit telescoping Keller maps

> **PROVISIONAL RESEARCH MEMO — NOT PEER REVIEWED.** This document records the
> outcome of the first milestone selected in `research/jc2-attack-memo.md`,
> Section 4 (affine-quotient obstruction). All computations are exact (symbolic
> over Q); the run is deterministic and rerunnable (`run_obstruction.py`). The
> result is a **negative (obstruction) result — Outcome B** of the milestone's
> falsifiable success criterion. Nothing here bears on the truth of JC2 itself.

## 1. Setting and question

Work over C. Let F3, F4, F5 be the three explicit members of the telescoping
family of `paper/main.tex` (Theorem 2.1 there; construction as in
`paper/verify.py`), namely `F_{E,h}` for

| name | E(nu)                    | h(t)        | geometric degree | component degrees |
|------|--------------------------|-------------|------------------|-------------------|
| F3   | `2 - 6*nu`               | `1 - (3/2)*t` | 3              | (7, 6, 4)         |
| F4   | `3 - 12*nu + 6*nu**2`    | `1 - 2*t`   | 4                | (12, 11, 4)       |
| F5   | `4 - 15*nu + 10*nu**3`   | `1`         | 5                | (17, 16, 4)       |

with `t = x*y`, `s = x**2*z`, `g = h(t) - s`, `nu = (1+t)*g`,
`A(nu) = int_0^nu xi*E(xi) dxi`, `C(nu) = 2*int_0^nu E(xi) dxi`, and

```
F = ( (A(nu) + (1+t)*g**2) / (x**2*g**2),  (C(nu) + 2*g) / (x*g),  x*g ).
```

**Sanity gate (checked).** Each of F3, F4, F5 is a polynomial map with
`det Jac = -2`, with the component degrees listed above (this is the sanity
check of `paper/verify.py`, re-run here, plus exact degree/leading-form data).
Geometric degrees and monodromy are taken from the paper and are not re-derived
here.

**Question (memo Section 4, exact form treated).** For F in {F3, F4, F5}, do
there exist

* an affine-linear submersion `pi = (l1, l2)`, `l_i = r_i*x + s_i*y + t_i*z + u_i`
  with 2x3 linear-part matrix of rank 2 (rank < 2 is excluded by the definition
  of submersion), and
* polynomials `G1, G2` in two variables with `deg G_i <= deg F`,

such that `pi o F = G o pi` identically and `det Jac(G)` is a nonzero constant?
(The attack memo writes `J(G) = 1`; any nonzero constant is accepted here, since
a target-linear rescaling normalizes the constant.) Coefficients are allowed in
C throughout; `deg F` denotes the maximum total degree of the components
(17, 12, 7 respectively — the result below is independent of this reading).

## 2. Result

**Theorem (obstruction; Outcome B — computer-checked, exact).**
For each F in {F3, F4, F5} there is **no** pair `(pi, G)` with `pi` affine-linear
of linear-part rank 2 and `G = (G1, G2)` polynomial — **of any degrees, with no
condition on Jac(G) at all** — such that `pi o F = G o pi`. In particular there
is none with `deg G_i <= deg F` and `det Jac(G)` a nonzero constant, so the
affine-quotient descent route of `jc2-attack-memo.md` Section 4 is closed for
these three maps. The infeasibility holds over C (the relevant coefficient
ideals contain 1 over Q).

Since no semiconjugation exists, the milestone's further steps (saturation by
`det Jac(G) != 0`, normalization `J(G) = 1`, and descent of exact collisions of
F to collisions of G) are moot; there is nothing to descend.

**What is *not* claimed.** The theorem is stated for these three explicit maps
in their given coordinate presentation. It says nothing about: nonlinear
polynomial submersions `pi` (equivalently, affine `pi` after conjugating F by a
nonlinear automorphism of A^3); the invariant-graph test or flat-specialization
test of the attack memo (tests 1 and 3); other members of the telescoping
family (although Section 5 notes the leading-form mechanism looks uniform in
the family, this has been verified only for the three members above); or the
truth of JC2.

## 3. Method

The computation is organized around the two-level structure noted in the attack
memo — the semiconjugacy equation is linear in G for fixed pi — but the linear
level collapses completely:

**Lemma 1 (fiber-derivative criterion).** Fix `pi = (l1, l2)` with rank-2 linear
part M, and let `v != 0` span ker M. For a polynomial L in C[x,y,z]:

```
L in C[l1, l2]   <=>   d_v L == 0,
```

where `d_v = v1*d/dx + v2*d/dy + v3*d/dz`.
*Proof.* Choose an affine form `lam` with `lam(v) = 1`; then `(l1, l2, lam)` is
an invertible affine coordinate system (the linear parts of l1, l2 span the
2-plane of forms annihilating v, and lam is outside it), and in these
coordinates `d_v` is the partial derivative in the third coordinate. So
`d_v L = 0` iff L is a polynomial in the first two coordinates, i.e.
`L = G(l1, l2)` for a (unique) polynomial G. QED

Consequently, for fixed pi, a polynomial G with `pi o F = G o pi` exists iff

```
E_i := d_v(l_i o F) == 0   for i = 1, 2,
```

and then `deg G_i = deg(l_i o F) <= deg F` **automatically** (invertible affine
substitutions preserve total degree), so the milestone's degree bound is not an
active constraint. The translation parts `u_i` are irrelevant (`d_v` kills
constants). G is eliminated from the search entirely.

**Lemma 2 (gauge reduction to three charts).** Target-affine changes
`(pi, G) -> (A*pi + b, (A,b) o G o (A,b)^{-1})`, A in GL2(C), preserve
existence, polynomiality, degrees, and `det Jac(G)`. Writing M = A*R with R the
reduced row-echelon form of the rank-2 linear part, every admissible pi is
gauge-equivalent to exactly one of

```
C12:  l1 = x + p*z,  l2 = y + q*z    (ker spanned by v = (-p, -q, 1)),  p, q in C
C13:  l1 = x + p*y,  l2 = z          (ker spanned by v = (-p, 1, 0)),   p in C
C23:  l1 = y,        l2 = z          (ker spanned by v = (1, 0, 0)).
```

So the entire problem is: for each map and each chart, decide whether the
polynomial identities `E_1 = E_2 = 0` (coefficientwise in x, y, z) have a
solution in the at most two scalar parameters. The coefficients generate an
ideal I in Q[p, q]; the chart is infeasible over C iff 1 in I. A generator that
is itself a nonzero rational constant (**unit certificate**) proves this by
inspection; a Groebner basis equal to [1] proves it in general.

**Why it fails — leading forms (by-hand proof, machine-confirmed).**
The exact leading forms of the components are:

| map | LF(F_1)            | LF(F_2)           | LF(F_3) |
|-----|--------------------|-------------------|---------|
| F3  | `2*x^3*y^3*z`      | `6*x^3*y^2*z`     | `-x^3*z` |
| F4  | `(3/2)*x^6*y^4*z^2`| `4*x^6*y^3*z^2`   | `-x^3*z` |
| F5  | `-2*x^9*y^5*z^3`   | `-5*x^9*y^4*z^3`  | `-x^3*z` |

In each map the component degrees are strictly decreasing, and the leading
forms of the first two components are **single monomials divisible by x·y·z**.
Suppose `l o F = G(l1, l2)` with l having x-coefficient r != 0. Since
`deg F_1 > deg F_2 > deg F_3`, the top form of `l o F` is `r*LF(F_1)`, a
monomial whose linear factors are x, y, z. But the top form of `G(l1, l2)` is
`G_D(l1_lin, l2_lin)` (D = deg G >= 1, and any binary form factors into linear
forms), a product of linear forms lying in the 2-dimensional span of `l1_lin,
l2_lin`. By unique factorization in C[x,y,z], x, y, z would all lie in that
2-plane — impossible. Hence `r_1 = r_2 = 0`; the same argument with `LF(F_2)`
gives `s_1 = s_2 = 0`; then the linear-part matrix has only its z-column and
rank <= 1, contradicting rank 2. This is the conceptual reason; the certified
result below does not depend on this argument (the machine computation is
primary and confirms it).

## 4. Computation and exact results

`run_obstruction.py` computes `E_1, E_2` exactly for each (map, chart), extracts
all coefficient generators in Q[p, q], records every unit certificate, and
computes a Groebner basis (grevlex, over Q). Results (full data in
`certificates.json`):

| map | chart | #generators | #unit certificates | sample unit certificate | Groebner basis |
|-----|-------|------------:|-------------------:|-------------------------|----------------|
| F3  | C12   | 23 | 5  | coeff of `x^3*y^3` in E_1 is `2`        | `[1]` |
| F3  | C13   | 18 | 7  | coeff of `x^3*y^2*z` in E_1 is `6`      | `[1]` |
| F3  | C23   | 6  | 8  | constant term of E_2 is `1`             | `[1]` |
| F4  | C12   | 61 | 17 | coeff of `x^6*y^4*z` in E_1 is `3`      | `[1]` |
| F4  | C13   | 48 | 15 | coeff of `x^6*y^3*z^2` in E_1 is `6`    | `[1]` |
| F4  | C23   | 16 | 16 | constant term of E_2 is `1`             | `[1]` |
| F5  | C12   | 71 | 30 | coeff of `x^9*y^5*z^2` in E_1 is `-6`   | `[1]` |
| F5  | C13   | 59 | 19 | coeff of `x^9*y^4*z^3` in E_1 is `-10`  | `[1]` |
| F5  | C23   | 18 | 19 | constant term of E_2 is `1`             | `[1]` |

Every cell contains at least one generator that is a nonzero rational constant
independent of (p, q), so every coefficient ideal is the unit ideal; the
Groebner bases confirm this redundantly. Hence no chart admits any parameter
value in C, which by Lemmas 1 and 2 proves the Theorem.

The sample certificates are exactly the ones predicted by the leading-form
analysis: for C12 the listed coefficient is `d/dz` applied to `LF(F_1)`
(e.g. for F5, `d/dz(-2*x^9*y^5*z^3)` contributes `-6*x^9*y^5*z^2`, and no other
term of `E_1 = d_v(F_1 + p*F_3)` can produce that monomial for any p, q); for
C13 it is `d/dy` applied to `LF(F_1)`; for C23, `E_2 = d/dx(F_3)` has constant
term `h(0) = 1`.

A near-miss worth recording: in chart C13 with `p = 0` (`pi = (x, z)`), the
**second** equation alone is satisfiable for F5 — indeed `l2 o F5 = x*g =
x - x^3*z` lies in C[x, z] — and the obstruction there comes only from the
first coordinate (`F5_1` genuinely involves y). So the obstruction is not
vacuously uniform across coordinates; the sampled cross-check below detects
this case (membership test for `l2 o F` passes, for `l1 o F` fails).

## 5. Cross-checks

1. **Positive controls (the pipeline can find quotients).** Two synthetic maps
   with known affine quotients are run through the identical pipeline:
   * `CTRL1 = (x + (y+z)^2, y + x^3, z - x^3)`: the pipeline finds chart C12,
     `(p, q) = (0, 1)`, i.e. `pi = (x, y+z)`, reconstructs `G = (w1 + w2^2, w2)`
     via the chart section, verifies `pi o F = G o pi` exactly, and computes
     `Jac(G) = 1` (nonzero constant — the full Outcome-A shape).
   * `CTRL2 = (x + (y+z)^2, (y+z)^3 - z + x^5, z - x^5)`: quotient found with
     `G = (w1 + w2^2, w2^3)`, and `Jac(G) = 3*w2^2` correctly classified as
     **not** constant.
   Both controls are asserted in the script; a machinery regression aborts the
   run before any obstruction claim is printed.

2. **Lemma-free sampled membership test.** Independently of Lemma 1, for
   sampled chart parameters ((p,q) in {(0,0), (1,2), (-3,5)} for C12, p in
   {0, 2, -7} for C13, and C23), the direct question "`is l_i o F a polynomial
   in (l1, l2)`?" is posed as a linear system in the unknown coefficients of
   `G_i` (basis `l1^a*l2^b`, `a + b <= deg(l_i o F)`), evaluated at random
   integer points, and shown **inconsistent modulo the prime 1000003**
   (rank of the augmented matrix exceeds rank of the coefficient matrix).
   Inconsistency mod p of a specialization certifies inconsistency of the exact
   rational system, hence non-membership. All 21 sampled (map, chart, parameter)
   cases are refuted (7 per map; fixed RNG seed 20260720). This samples the
   parameter space only — the covering statement rests on Section 4 — but it
   exercises none of the derivative/gauge reasoning.

## 6. Status, and what remains open

* **Established (exact, rerunnable):** the Theorem of Section 2 — no affine
  two-dimensional semiconjugate quotient exists for F3, F4, F5, in the strong
  degree-free, Jacobian-condition-free form. This is the milestone's Outcome B
  ("proved obstruction for these members"), strengthened from "no Keller
  quotient of bounded degree" to "no polynomial semiconjugation at all".
* **Not established / open routes** (in the attack memo's terms): the
  invariant-graph test (nonlinear graphs `z = h(x,y)`); nonlinear polynomial
  submersions `pi`; flat specializations; the same question for other family
  members or after nonlinear coordinate changes of A^3. The leading-form
  mechanism (top forms of the first two components are monomials divisible by
  x·y·z of strictly separated degrees) appears to be uniform in the telescoping
  family — for `deg E = n`, the expected top forms are
  `x^(3n)*y^(n+2)*z^n` and `x^(3n)*y^(n+1)*z^n` — but this has been verified
  here only for the three explicit members and is **conjectural** beyond them.
* The obstruction is consistent with the attack memo's Section 2 heuristic that
  the third variable is an essential state variable for the telescoping
  mechanism; it converts that heuristic into a theorem for the affine-linear
  descent route on these members.

## 7. Reproduction

```
cd research/jc2-affine-quotient-obstruction
uv run --with sympy python run_obstruction.py
```

Runtime ~2.5 s (sympy 1.14.0, Python 3.12). The run is deterministic, writes
`certificates.json` (sanity data, all unit certificates, Groebner bases,
control records, mod-p cross-check records), and ends with

```
OBSTRUCTION CERTIFIED: no affine quotient for F3, F4, F5
```

Had a quotient with nonzero-constant `Jac(G)` been found for a Keller member,
the script would instead print a loud JC2-candidate warning (demanding
independent triple-verification and collision descent) and exit with code 2.

## Files

* `run_obstruction.py` — the complete exact computation (this memo's evidence).
* `certificates.json` — machine-written certificate data from the last run.
* `../jc2-attack-memo.md` Section 4 — the milestone this executes.
* `../../paper/main.tex`, `../../paper/verify.py` — construction and baseline
  checks for F3, F4, F5.

---

## Independent verification (appended by a separate verifying agent, 2026-07-20)

**Verdict: ENDORSED.** The Theorem of Section 2 (Outcome B, in the strong
degree-free, Jacobian-condition-free form) was re-certified by an independently
written script, `verify_obstruction_independent.py` (same directory; run with
`uv run --with sympy python verify_obstruction_independent.py`; ~2 s; all
arithmetic exact). The verification deliberately used **different routes** than
`run_obstruction.py`; the fiber-derivative machinery of Lemma 1 was *not* used
in the primary re-certification. All checks passed. Findings:

1. **Reconstruction of F3, F4, F5 from `paper/main.tex` (independent code).**
   The three maps were rebuilt from the paper's construction with separately
   written antiderivatives and exact multivariate division (remainder asserted
   0, rather than `cancel`). The paper's polynomiality criterion (c1)–(c3) was
   verified for each (E, h); F4's E was confirmed to equal the paper's Legendre
   recipe `2 - 6*nu + L_2(2*nu - 1)`. `det Jac = -2` was recomputed by an own
   cofactor expansion; component degrees (7,6,4), (12,11,4), (17,16,4) and the
   leading-form table of Section 3 were confirmed exactly. The expanded
   polynomials **agree term-by-term** with the output of the builder's
   `family()` for all three members, and the builder's (E, h) data match the
   paper's.

2. **Structural facts re-derived.** (a) *Chart completeness*: the RREF of a
   rank-2 2x3 matrix has pivot columns {1,2}, {1,3} or {2,3}, whose general
   shapes are exactly C12, C13, C23; confirmed by an own Fraction-arithmetic
   RREF on 603 rank-2 samples (random plus adversarial), each with exact gauge
   factor M = A*R, det A != 0, and the kernel formulas checked symbolically.
   No chart is missing. (b) *Lemma 1* was spot-checked on 5 random affine
   rank-2 pi in both directions, with the converse decided by explicit affine
   coordinate inversion, not by the lemma. (c) *Linearity in G*: for fixed pi
   the semiconjugacy is a linear system in G's coefficients; on a randomly
   gauged copy of CTRL1 an own exact-Fraction Gaussian elimination solved it
   and (pi, G) was verified with Jac(G) = 1; for F3 with
   pi = (x+2y+3z+1, 4x+5y-6z-1) the exact system over Q is inconsistent.

3. **Re-certification of all 9 (map, chart) cells by a derivative-free route.**
   Substitution route: complete (l1, l2) to an invertible affine coordinate
   system, substitute the inverse, and use the from-scratch criterion
   "membership in C[l1, l2] iff no monomial contains the third coordinate".
   Certification per chart used deliberately different engines than the
   builder's: C23 by direct inspection (all obstructing coefficients are
   explicit nonzero rationals); C13 by an **own univariate gcd over Fraction**
   (gcd of all generators = 1 in Q[p], no Groebner engine involved); C12 by
   unit-constant generators (nonzero rationals independent of p, q — valid for
   all complex parameters), with the unit sitting at the leading-form position
   (coefficient of X^a Y^b Z^c for (a,b,c) = the LF(F_1) exponents, value = the
   LF coefficient), confirming the Section 3 leading-form collapse from
   scratch. Additionally, a third flavour for C12: point-evaluation
   (fiber-constancy) differences `l_i o F(P+v) - l_i o F(P)` at explicit
   integer points give ideals in Q[p,q] with Groebner basis [1] for all three
   maps; and for C23 explicit integer pairs (0,0,-1), (1,0,-1) in one pi-fiber
   split under F_3 (values 0 vs 2, exact). The Section 4 near-miss (F5, C13,
   p = 0, second equation satisfiable) was confirmed, including that it does
   *not* occur for F3, F4.

4. **Positive controls and artifact audit.** CTRL1/CTRL2 were pushed through
   the independent substitution pipeline: unique feasible chart C12 at
   (p,q) = (0,1), G reconstructed and verified exactly, Jac(G) classification
   (1 vs 3*w2^2) correct. The memo's Section 4 table, the sample unit
   certificates, and `certificates.json` were all re-derived by an own
   implementation of the builder's fiber-derivative route and agree; Groebner
   bases were recomputed in **both grevlex and lex** and equal [1] in every
   parametric cell. A rerun of `run_obstruction.py` reproduced
   `certificates.json` byte-identically (deterministic as claimed).

**Nits found (no bearing on the result):**

* Section 1, "Question": the parenthetical "(17, 12, 7 respectively)" lists
  the maximum component degrees in reverse order relative to the enumeration
  {F3, F4, F5}; it should read (7, 12, 17). Immaterial — the theorem proved is
  degree-free.
* Section 4 table semantics: "#generators" counts *distinct* coefficient
  values (deduplicated across both equations) while "#unit certificates"
  counts term *occurrences*; hence rows where units exceed generators (e.g.
  F3/C23: 6 vs 8) are consistent, not an error.

**Scope agreed.** The memo's claim boundaries are accurate: the result covers
only affine-linear pi for these three explicit members in their given
coordinates; monodromy/geometric-degree data are quoted from the paper, not
re-derived; the family-uniformity remark is correctly labelled conjectural;
nothing is claimed about JC2 itself. The mod-p cross-check is correctly
described as sampling parameters only. No overclaims were found.

* `verify_obstruction_independent.py` — the independent verification script
  (this section's evidence).
