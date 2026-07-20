# AUDIT (classical half) — auditor A — Theorem J3, band-2 two-sided non-square sector

```
STATUS: AUDIT COMPLETE (classical half only; quantum half NOT audited here)
VERDICT: CONFIRMED WITH MINOR REPAIRS
         The classical Theorem J3 is TRUE as stated. Every step of the
         written argument is certified below; three presentation-level
         defects were found and repaired (none mathematical): a bracket-
         orientation/sign inconsistency with the frozen M4 conventions, a
         missing (trivial) const = 0 branch in the endgame, and an implicit
         band bound on D. An exhaustive bounded search (all coefficient
         degree profiles deg a_k <= 3) found no violating solution and is
         consistent only with the theorem.
```

Audited artifact: `theorem_J3_band2_nonsquare.md` (classical Statement + Proof
(classical) sections). Conventions frozen by
`../band2-classical-proved/M4_proof_memo.md` (Sections 2–3). Referee stance:
attempt to refute.

Audit scripts (this directory; run with a Python that has SymPy ≥ 1.12):

- `audit_classical_A_identities.py` — forward cascade identities from scratch
  (memo orientation), 4 degree profiles with fully symbolic coefficients;
  Wronskian and J2 lemma certificates (bounded-exhaustive). **All checks pass.**
- `audit_classical_A_search.py` — exhaustive bounded search deg a_k ≤ 3
  (Gröbner emptiness certificates per stratum + independent full-system
  Gröbner + random exact sweeps + generic-rank obstruction + positive
  controls). **All checks pass; no violating solution exists.**

Both scripts exit 0. The claims below cite the specific check names.

---

## 0. Statement as audited

X = Σ_{k=−2}^{2} x^k a_k(τ), D = Σ_{l=−2}^{2} x^l b_l(τ), a_k, b_l ∈ C[τ],
{D,X} ∈ C* in the memo orientation {G,F} = G_ξF_x − G_xF_ξ, with a₂ ≠ 0 and
a₋₂ ≠ 0. Claim: at least one of a₂, a₋₂ is a scalar-square c·h², c ∈ C*,
h ∈ C[τ] (constants count as squares; "non-square" is therefore in particular
nonconstant).

Two statement-level observations (repairs R1, R2 below):

- The displayed D omits the range of l. The proof requires supp(D) ⊆ [−2,2];
  the "band width ≤ 2 on both sides" phrasing intends this, but the display
  should say Σ_{l=−2}^{2}.
- "No membership/divisibility hypothesis" is legitimate and verified: the
  band-≤2 elements of C[x^{±1}, ξ] are exactly the Σ x^k a_k(τ) with
  a_k ∈ C[τ] free (x^kτ^j = x^{k+j}ξ^j, so no divisibility is imposed once x
  is inverted and ξ is not). Every audited step uses only that a_k, b_l are
  polynomials in τ; τ-divisibility is used nowhere (loophole hunt, §3.4).

## 1. Per-step certificate table

Orientation note for the whole table: all component formulas below are in the
**memo orientation**. The theorem prose (and `verify_J3.py`) use the opposite
bracket; see repair R2. Zero-equalities are orientation-independent.

| Step | Claim in the theorem | Status | Certificate |
|---|---|---|---|
| C_m derivation | nine memo equations, coefficient of x^m of {D,X}, support ⊆ [−4,4] | CERTIFIED | `direct {G,F} (memo orientation) x-components == memo C_m formulas` — from-scratch two-variable differentiation, 4 symbolic profiles P1–P4 (all-deg-3; deg-2 a vs deg-4 b; mixed (4,1,2,0,3)/(3,2,4,1,3); a₁=a₋₁=0) |
| m=±4 | Wronskians vanish ⇒ b₂ = λ₂a₂, b₋₂ = μ₂a₋₂, λ₂, μ₂ ∈ C | CERTIFIED | identity `C_4 == 2(a₂b₂′−a₂′b₂)`, `C_-4 == −2(a₋₂b₋₂′−a₋₂′b₋₂)` (4 profiles); Lemma W hand proof (§2.1) + machine: rational identity `a²(b/a)′ = ab′−a′b`, top-coefficient identity `(e−d)a_d b_e` for d=1..4, e=0..8 (kills deg b ≠ deg a for *all* degrees), and bounded-exhaustive kernel certificate: for every a with lead ≠ 0, d = 1..4, the kernel of b ↦ ab′−a′b on deg ≤ d is exactly C·a (all d×d minors + Rabinowitsch at lead = unit ideal). Constant-ratio-in-C(τ) inference and polynomial conclusion verified (§2.1); b₂ = 0 (λ₂ = 0) included |
| m=±3, Lemma J2 | b₂-substitution leaves exactly 2a₂u′ − a₂′u = 0, u = b₁−λ₂a₁; nonzero u ⟺ a₂ scalar-square; nonsquare a₂ ⇒ b₁ = λ₂a₁ with the same λ₂ | CERTIFIED | identity `C_3|_{b₂=λa₂} == 2a₂u′ − a₂′u` (b₁ generic, 4 profiles). **Mirror orientation derived from scratch**: `C_-3|_{b₋₂=μa₋₂} == −(2a₋₂v′ − a₋₂′v)`, v = b₋₁−μ₂a₋₁ — the same homogeneous form (global −1 factor irrelevant for = 0), so Lemma J2 applies verbatim on the mirror side with a₋₂. Lemma J2 both directions: identity `u(2au′−a′u) = a(u²)′ − a′u²` ⇒ Wronskian(a,u²)=0 ⇒ u² = c·a (Lemma W); converse a = dh² ⇒ u = eh solves (symbolic h, deg 4); d=2 exhaustive: kernel nonzero ⟺ disc = 0 on lead ≠ 0 (ideal computation, both directions, kernel vector = the square root); d=4: square-parametrization direction + generic/specimen rank; top-coefficient identity `(2e−d)a_d u_e` (d=1..4, e=0..8) ⇒ deg u = d/2 forced, odd-degree a₂ has u = 0 with no degree bound needed. "Same λ₂" is by construction (u is defined with the λ₂ from m=4) |
| m=±2 | (1,1) terms cancel identically; remainder 2a₂(b₀′−λ₂a₀′) = 0 and (memo sign) −2a₋₂(b₀′−μ₂a₀′) = 0; ⇒ b₀′ = λ₂a₀′ = μ₂a₀′ | CERTIFIED | identities `C_2|_{b₂=λa₂,b₁=λa₁} == 2a₂(b₀′−λa₀′)` and `C_-2|… == −2a₋₂(b₀′−μa₀′)` with b₀ generic; separate check of the bare (1,1) cancellation. Last inference is the domain step a·g = 0, a ≠ 0 ⇒ g = 0 in C[τ] (hand, §2.2). Corollary (λ₂−μ₂)a₀′ = 0 |
| m=±1 | collapse to (λ₂−μ₂)U = 0 and (λ₂−μ₂)L = 0 (theorem's sign; memo sign (μ₂−λ₂)) | CERTIFIED | exact-residual identities with b₀ **generic**: `C_1|_{4 props} == (μ−λ)U + a₁(b₀′−λa₀′)` and `C_-1|_{4 props} == (μ−λ)L − a₋₁(b₀′−μa₀′)`; the prose's "using b₀′ = λ₂a₀′, resp. μ₂a₀′" is load-bearing and correct — auxiliary check: using only b₀′=λ₂a₀′ leaves C₋₁ = (μ−λ)(L + a₋₁a₀′), which the second m=−2 relation kills. Both relations hold simultaneously, so the collapse is exact |
| m=0 | surviving terms telescope to (λ₂−μ₂)·(2a₂a₋₂+a₁a₋₁)′ = const ≠ 0 (memo sign (μ₂−λ₂)); hence λ₂ ≠ μ₂, a₀, b₀ constant, M linear nonconstant | CERTIFIED | identity `C_0|_{4 props} == (μ−λ)·d/dτ(2a₂a₋₂+a₁a₋₁)` — full C_0 with all cross terms, 4 profiles; structural check that C_0 contains **no** b₀ or a₀ terms at all (the (0,0) pair contributes k·(…) − l·(…) = 0). λ₂ = μ₂ branch closed: if λ₂ = μ₂ then D = λ₂X + γ and ALL nine components vanish identically (machine check), so {D,X} = 0 ∉ C*. Then a₀′ = 0 from (λ₂−μ₂)a₀′ = 0, b₀′ = λ₂a₀′ = 0, and M′ = const/(μ₂−λ₂) ∈ C* ⇒ M linear nonconstant |
| Endgame invariants | a₋₁U = (a₂a₋₁²)′, a₁L = (a₋₂a₁²)′ | CERTIFIED | polynomial identities, 4 profiles |
| Endgame dichotomy | U = 0 ⇒ a₂a₋₁² = const; a₋₁ = 0 forced; mirror a₁ = 0; then M = 2a₂a₋₂ linear ⇒ deg a₂ + deg a₋₂ = 1 ⇒ one factor constant = square, contradiction | CERTIFIED (with repair R3) | hand degree argument (§2.3) incl. the const = 0 branch missing from the prose; machine: top-coefficient identity coeff_{τ^{d₂+d₋₂−1}}((2a₂a₋₂)′) = 2(d₂+d₋₂)·lead·lead ≠ 0 for d₂,d₋₂ ≥ 1 (all nine degree pairs ≤ 3), so M′ constant is impossible with both extremes nonconstant; plus the Part B Gröbner emptiness certificates below, which re-prove the whole endgame ideal-theoretically on every deg ≤ 3 stratum |
| Edge case (task item 6) | a₂a₋₁² = nonzero const with a₋₁ constant nonzero ⇒ a₂ constant | NOT A LOOPHOLE | under the frozen definition constants ARE scalar-squares (c·1²), so a₂ constant directly contradicts "a₂ non-square". The definition ("constants included") makes this branch airtight |

## 2. Hand arguments the computation cannot reach (with proofs)

**2.1 Lemma W (Wronskian ⇒ constant ratio).** a, b ∈ C[τ], a ≠ 0,
ab′ − a′b = 0 ⇒ b = λa with λ ∈ C. *Proof.* In C(τ), (b/a)′ = (ab′−a′b)/a² = 0.
The kernel of d/dτ on C(τ) is C: write b/a = P/Q, gcd(P,Q) = 1, Q monic; then
P′Q = PQ′, so Q | Q′P, so Q | Q′; since char 0, deg Q′ = deg Q − 1 unless
Q′ = 0, forcing Q′ = 0, Q = 1; then P′ = 0 and (char 0) P ∈ C. So b = λa as
polynomials (λ = 0 allowed, covering b = 0). ∎ — The machine complements this
with the bounded-exhaustive kernel certificates (all a of degree 1..4 with
lead ≠ 0) and the top-coefficient identity that excludes kernel elements of
any other degree; the hand proof covers all degrees.

**2.2 Domain steps.** C[τ] is an integral domain: a·g = 0, a ≠ 0 ⇒ g = 0.
Used at m=±2 (divide by 2a₂, −2a₋₂), in the endgame (a₂a₋₁² = 0, a₋₁ ≠ 0 ⇒
a₂ = 0), and to bound deg b₀ (b₀′ = λ₂a₀′ ⇒ deg b₀ ≤ max(deg a₀, 0)).

**2.3 Endgame degree dichotomy (repaired form).** From λ₂ ≠ μ₂ and
(μ₂−λ₂)U = 0: U = 0, hence (a₂a₋₁²)′ = 0, hence a₂a₋₁² = κ ∈ C. If a₋₁ ≠ 0:
either κ = 0, and then a₂ = 0 (domain), contradicting a₂ ≠ 0 — **this branch
is missing from the prose (repair R3)** — or κ ≠ 0, and then
deg a₂ + 2 deg a₋₁ = 0, so a₂ is a nonzero constant, i.e. a scalar-square,
contradicting non-squareness. Hence a₋₁ = 0; mirror gives a₁ = 0. Then
M = 2a₂a₋₂ with M′ a nonzero constant, so deg a₂ + deg a₋₂ = 1; both are
nonzero, and the factor of degree 0 is a nonzero constant = scalar-square —
contradiction. (Equivalently: both non-squares are nonconstant, so
deg M ≥ 2 > 1.) ∎

## 3. Loophole hunt (task item 7)

**3.1 a₁ = 0 or a₋₁ = 0 from the start.** No escape. Profile P4
(a₁ = a₋₁ = 0) passes every identity; the m=±3 step then reads
2a₂b₁′ − a₂′b₁ = 0 directly, i.e. u = b₁, and nonsquare a₂ still forces
b₁ = 0 = λ₂a₁. The endgame only simplifies (U = L = 0 automatic, M = 2a₂a₋₂),
and the same degree contradiction closes it. Also covered by the Part B/C
searches (middle profiles include the zero polynomial).

**3.2 λ₂ = μ₂ consistency when some intermediate vanishes.** Impossible under
the hypotheses: once both extremes are nonsquare, the four proportionalities
are forced with no homogeneous freedom, and λ₂ = μ₂ makes D = λ₂X + (b₀ with
b₀′ = λ₂a₀′) = λ₂X + γ, for which **all nine** components vanish identically
(machine check `lam=mu branch`), i.e. {D,X} = 0 ∉ C*. No intermediate
vanishing (a₁ = 0, a₋₁ = 0, a₀′ = 0, U = 0, L = 0 …) changes this: C_0 =
(μ₂−λ₂)M′ is an identity, so C_0 ∈ C* forces λ₂ ≠ μ₂ unconditionally.

**3.3 Bracket constant vs normalization G → G − λF.** The written proof never
uses the normalization; and {G−λF, F} = {G,F}, so it could be used freely
without touching the constant. The constant's exact nonzero value is used
nowhere except as "≠ 0" at m=0; rescaling D rescales it harmlessly. No
interaction.

**3.4 Secret τ-divisibility.** None. Audited step by step: Lemma W, Lemma J2,
the domain steps and the degree dichotomies use only that coefficients are
polynomials in τ. The identity checks run with fully generic a_k, b_l ∈ C[τ]
(no divisibility imposed), and the search parts A–D impose none. Contrast
with the M4 polynomial theorem, whose final contradiction *does* invoke
τ² | b₋₂; J3's classical proof is genuinely divisibility-free, as claimed.
(Boundary remark: the statement is for C[x^{±1}, ξ] — ξ **not** inverted —
which is exactly what makes a_k ∈ C[τ] polynomials; nothing here is claimed
for C[x^{±1}, ξ^{±1}], where coefficients would be Laurent in τ.)

**3.5 Wider D.** If D were allowed band 3, the cascade would need the m=±5
equations first; the statement's "band ≤ 2 on both sides" excludes this, but
the display should fix the range (repair R1).

## 4. Exhaustive bounded search (deg a_k ≤ 3)

All in `audit_classical_A_search.py`; every check passed, no violating
solution found. The search is *graded*: the bilinear system is linear in the
b-side, and the four outer equations have machine-certified complete solution
sets, which linearizes the hunt without loss. Two normalizations are used,
both legitimate because they preserve the hypotheses and the scalar-square
class of every coefficient: **N1** — X → tX composed with the diagonal
symplectic scaling (x,ξ) → (ρx, ρ⁻¹ξ) sends a_k → tρ^k a_k; solving
tρ² = 1/lead(a₂), tρ⁻² = 1/lead(a₋₂) over C makes both extremes monic;
**N2** — D → D − λ₂X (the memo's allowed normalization, bracket-preserving)
sets λ₂ = 0 in the reduced variant.

- **Part A (fully independent exhaustion — no cascade input whatsoever).**
  The complete bilinear system {τ-coeffs of C_m (m≠0), τ-coeffs ≥ 1 of C_0,
  [C_0]₀ = c ≠ 0} with **all ten coefficient polynomials undetermined** —
  extremes monic (N1) on each of the 9 nonsquare strata, middles a₁, a₀,
  a₋₁ fully symbolic of degree ≤ 3, all five b_l fully symbolic of degree
  ≤ 3, disc ≠ 0 and c ≠ 0 by Rabinowitsch inverses — yields the **unit
  ideal** (Gröbner) on every stratum. This is machine quantifier
  elimination over the whole stratum: emptiness for every coefficient
  choice, every b of degree ≤ 3, every c ≠ 0, with no lemma grading at all.
  Together with the lemma-certified degree bounds on b (which exclude
  b-degrees above 3 given a-degrees ≤ 3), Part A alone re-proves the
  deg ≤ 3 exhaustion independently of Part B.
- **Part B (the deg ≤ 3 exhaustion, lemma-graded).** The strata (deg a₂, deg a₋₂) ∈
  {1, 2-nonsquare (disc ≠ 0), 3}² exhaust all nonzero nonsquare extremes of
  degree ≤ 3 (degree 0 is a square; odd degrees are never scalar-squares;
  degree-2 scalar-squares are exactly disc = 0). Middles a₁, a₀, a₋₁ fully
  symbolic of degree ≤ 3 (subsuming every lower/zero profile). After the
  certified outer reductions b±₂ = (λ₂,μ₂)a±₂, b±₁ = (λ₂,μ₂)a±₁ — complete
  for b of **arbitrary** degree by Lemma W / Lemma J2 — both a "full" variant
  (λ₂, μ₂ free, b₀ undetermined to degree 4, complete since C_2 = 0 bounds
  deg b₀ ≤ deg a₀) and a "reduced" variant (N2, collapse identities + domain
  step) run Gröbner + Rabinowitsch(c, disc) per stratum: **all 9 strata give
  the unit ideal in both variants** (18 emptiness certificates). Every
  solution of the C_m system with both extremes nonzero nonsquare of degree
  ≤ 3 is thereby excluded — for all coefficient values, all b degrees, all
  c ≠ 0.
- **Part C (random exact sweep + positive controls).** 1125 systems — one
  per degree profile (deg a₂, deg a₋₂) ∈ {1,2,3}² × (deg a₁, a₀, a₋₁) ∈
  {absent,0,1,2,3}³ — with random integer coefficients, both extremes
  verified nonsquare via square-free multiplicity parity (field-independent
  in char 0), b free to degree 4, plus 146 extra fully random draws with b
  free to degree 5: the constant functional [C_0]₀ vanished identically on
  the exact QQ-solution space of the remaining constraints in every single
  case. Positive controls: the machinery *does* report c ≠ 0 solutions for
  the memo's localized near-miss (a₂ = τ, a₋₂ = 0 — outside J3's
  hypotheses) and for the band-1 pair X = x, D = ξ — so the sweep's
  emptiness is informative, not vacuous.
- **Part D (generic-rank obstruction).** With a₂, a₋₂ fully symbolic (random
  middles), appending the c-row to the constraint matrix does not raise the
  rank over the rational-function field: no solution with c ≠ 0 at the
  generic point of the stratum. Completed on the low strata within the
  per-stratum time budget (90 s); the heavier deg-3 symbolic ranks time out
  and are skipped — pure redundancy, since Parts A/B cover those strata
  ideal-theoretically at every point, not just generically.

Conclusion of the search: **no violating solution exists at deg a_k ≤ 3**;
every outcome matches the theorem's prediction (solutions with c ≠ 0 exist
only when an extreme vanishes or a square appears).

## 5. Repairs (all presentation-level; proofs above)

- **R1 (statement).** Display D = Σ_{l=−2}^{2} x^l b_l(τ); the band bound on
  D is used by the proof (it is intended by "band width ≤ 2 on both sides").
- **R2 (orientation/sign bookkeeping).** The theorem file cites the memo's
  frozen conventions but writes the collapses with the sign of the *opposite*
  bracket: machine-certified, `verify_J3.py`'s `PB(D,X) = D_xX_ξ − D_ξX_x`
  is the **negative** of the memo's {D,X} = D_ξX_x − D_xX_ξ, and in memo
  orientation the exact identities are C_1 = (μ₂−λ₂)U, C_{−1} = (μ₂−λ₂)L,
  C_0 = (μ₂−λ₂)·(2a₂a₋₂+a₁a₋₁)′ — the prose's (λ₂−μ₂) prefactors match the
  flipped bracket. No logical damage anywhere (only = 0 and ∈ C* are used),
  but under the M4 memo's "fixed bracket orientation" policy the file should
  either flip the prefactors or state its orientation. Likewise the prose's
  "2a₋₂(b₀′ − μ₂a₀′) = 0" is, in memo orientation, −2a₋₂(b₀′ − μ₂a₀′) = 0
  (equivalent).
- **R3 (endgame branch).** "If a₋₁ ≠ 0, then a₂ = const/a₋₁² …" silently
  assumes const ≠ 0. Add: if the constant is 0, then a₂a₋₁² = 0 with
  a₋₁ ≠ 0 gives a₂ = 0 in the domain C[τ], contradicting a₂ ≠ 0. (§2.3.)
- **R4 (dependency note).** The README flags Lemma J2 as a missing
  dependency. For the classical half this blocker is discharged: Lemma J2 is
  fully reconstructed and proved here (identity u(2au′−a′u) = a(u²)′ − a′u²
  plus Lemma W), and machine-certified exhaustively at bounded degree. J1/P3
  is not used by the classical Theorem J3 proof (only by the "Consequences"
  section, which is outside this audit's scope).

## 6. What this audit does NOT cover

The quantum half (Theorem J3q), the "Consequences and the shape of the
induction" section, and the M5 program. The `verify_J3.py` checker itself was
re-run (13/13 PASS) but is computational support only; note that its classical
section checks the *reduced* pair (weak direction) with a₀ constant, whereas
this audit certifies the *forward* cascade with a₀ generic — the direction
that carries the proof.
