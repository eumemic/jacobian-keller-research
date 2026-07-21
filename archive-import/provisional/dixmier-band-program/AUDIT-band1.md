# Hostile audit of Theorem P3 (quantum band-1 rigidity)

**Auditor role:** hostile referee. **Target:** `band1-rigidity-milestone.md` §3
(Theorem P3). **Date:** 2026-07-20. **Environment:** sympy 1.14 via
`uv run --with sympy`; exact arithmetic; no git commits.

## VERDICT: CONFIRMED_WITH_REPAIRS

The **load-bearing** content of P3 is **correct** and is now proved cleanly at
arbitrary Bernstein degree:

- **(a) A₁-membership ⇒ affine.** Every band-1 pair `X,D ∈ A₁` with `[D,X]=1` is
  an affine symplectic pair. — **CONFIRMED** (repaired proof below; the memo's
  own proof of this is valid *for the A₁ conclusion* but reaches it through a
  flawed sub-argument).
- **(c) Generation.** Such a pair generates `A₁`; hence a pair generating a
  *proper* subalgebra of `A₁` must have a generator with `ad(E)`-support outside
  `{−1,0,1}`. — **CONFIRMED**.

Two defects were found; neither touches (a)/(c), but both must be repaired
because they are stated as theorem content:

- **(b) The B∖A₁ classification is WRONG/incomplete.** P3 claims the non-A₁
  solutions "form the one-parameter polar families such as `(x, ∂+c/x)`." This is
  **false**: there is a genuine non-degenerate, **non-affine, non-polar** family
  in `A₁[x⁻¹]∖A₁`, e.g.

      X = x²∂ + x + x⁻¹,   D = x⁻¹      with   [D,X] = 1.

  This is the "branch (B)" the memo itself names but mishandles.
- **The "Fourier automorphism" is mis-identified.** The map used to fold the
  branches/degenerate cases, `x ↦ −∂, ∂ ↦ x`, is the genuine Weyl-Fourier
  automorphism of `A₁` — but it does **not** preserve `B = A₁[x⁻¹]` (it sends
  `x↦−∂`, and `∂` is not invertible in `B`). The correct band-symmetry is the
  band-reversal `σ: x ↦ x⁻¹, E ↦ −E`, which **is** an automorphism of `B` but
  **not** of `A₁`. Using the wrong map is why branch (B) was wrongly reported to
  "back-substitute to affine forms."

**Consequence for band-2.** Band-2's `a₂=0` sector must cite **(a)/(c)** — "band-1
A₁-pairs are affine, hence generate A₁, hence are tame / never a DC1
counterexample" — and must **not** cite the B∖A₁ classification. The clean citable
statement is written to
[`../../../research/band2-square-sector/quantum-band1.md`](../../../research/band2-square-sector/quantum-band1.md).

**Novelty caveat (from the literature check).** P3's true content is essentially a
known result. Bavula–Levandovskyy (2020) prove Dixmier's Problem 1 for elements
that are sums of ≤2 homogeneous `ad(E)`-components, in *exactly* this crossed-
product localization; Han–Tan (2024) extend it. Band-1 (≤3 components) is a mild
instance. The memo's own hedge ("presumably accessible to Dixmier's ad-semisimple
analysis — literature check pending") is confirmed. See §5.

---

## 1. Per-step certificate table

Every row is machine-checked by the named script (all end with a PASS banner).

| # | Audit demand / claim | Method | Result |
|---|---|---|---|
| 1 | Independent crossed-product engine; re-derive all 5 component equations of `[D,X]=1` | `audit_band1_engine.py` §A,§B — engine cross-checked vs. naive operator model on 200 random monomials; each `m∈{−2..2}` component matched a hand-derived closed form | **PASS** |
| 1 | Certify degree `±2` periodicity ⇒ proportionality (`b₁=λa₁`) | `..._engine.py` §C — kernel of `b₁↦b₁(E+1)a₁−a₁(E+1)b₁` is exactly `C·a₁` for generic `a₁` of each degree 0–4 and 8 special `a₁`; plus the analytic 1-periodic-rational argument | **PASS** |
| 1 | Certify degree-0 telescoping `[D,X]₀=(λ−μ)(V−V^{[1]})`, `V=a₁(E−1)a₋₁(E)` | `..._engine.py` §D | **PASS** |
| 2 | **CRITICAL** branch (B) (`a₁` linear, `a₋₁` const): trace back-substitution | `audit_band1_branchB.py` §E — branch (B) satisfies the **full** system for all params; concrete `X=x²∂+x+1/x, D=1/x` has `[D,X]=1`, is non-affine, non-polar, and `∉A₁` | **GAP FOUND — repaired** |
| 3 | Enumerate all vanishing patterns of `(a₁,a₋₁)` and close each | `audit_band1_search.py` §H1–H4 (cases I–IV) + `audit_band1_classification.py` §I1 (degree-drop lemma closes II/III at arbitrary degree) | **PASS** |
| 3 | Verify the Fourier map negates `ad(E)`-degrees and maps the band presentation | `..._branchB.py` §F — genuine Fourier `F` has **no** right-inverse for `F(x)=−d` in `B` (⇒ `F` doesn't preserve `B`); band-reversal `σ` **does**, is an involution, swaps branch (A)↔(B), and does **not** preserve `A₁` | **DEFECT FOUND — repaired** |
| 4 | Exhaustive bounded search, coefficient degrees ≤ 4; classify every family | `audit_band1_exhaustive.py` §J1 (degree-profile sweep: case-I solvable only for `(deg a₁,deg a₋₁)∈{(0,1),(1,0)}`), §J1b (rigorous ideal-saturation: A₁-membership ⇒ affine only, deg ≤ 2 in 0.7 s and deg ≤ 3 in 52 s), §J3 (correlation `both∈A₁ ⇔ both affine` on 40 exact solutions) | **PASS** |
| 5 | Literature (Dixmier 1968, Joseph 1975, …) — novelty framing | §5 below | **Prior art found** |
| — | Claim (a): A₁ solutions = affine symplectic pairs | repaired proof §3; `audit_band1_classification.py` §I3, `..._exhaustive.py` §J1b | **CONFIRMED** |
| — | Claim (b): B∖A₁ solutions = polar families | §2 | **REFUTED (incomplete)** |
| — | Claim (c): proper-subalgebra generator leaves the band | §3 Cor.; generation identity machine-checked | **CONFIRMED** |

Scripts (all alongside this file):
`audit_band1_engine.py`, `audit_band1_branchB.py`, `audit_band1_search.py`,
`audit_band1_classification.py`, `audit_band1_exhaustive.py`.

---

## 2. Issue 1 — the branch-(B) gap (audit demand #2)

The memo's endgame: `V(E)=a₁(E−1)a₋₁(E)` is linear, so (deg sum 1) one factor is
constant, giving two branches — **(A)** `a₁` const, `a₋₁` linear, and **(B)** `a₁`
linear, `a₋₁` const — and it asserts (B) is "its Fourier mirror" and
"back-substituting yields the affine forms."

**This is where P3 breaks.** Back-substituting branch (B) does **not** yield affine
forms. Take `λ=0, μ=1, a₋₁≡1, a₀=b₀=0`, and `V(E)=E` (so `a₁(E)=V(E+1)=E+1`):

    X = x·(E+1) + x⁻¹ = x²∂ + x + 1/x,     D = x⁻¹ = 1/x.

Direct Weyl-algebra computation (not the crossed product): `DX = x∂+1+x⁻²`,
`XD = x∂+x⁻²`, so `[D,X]=1`. Support `= {−1,1} ⊆ {−1,0,1}` (band-1). But
`X` contains `x²∂` — **not affine**, and not of the claimed form
`X = δx + α − ε∂ + c·x⁻¹/δ′`. Branch (B) is **non-degenerate** by the memo's own
definition (`a₁≠0` and `a₋₁≠0`), so it is not covered by the "degenerate subcases."
`audit_band1_branchB.py` §E verifies the whole branch-(B) family solves the full
system for all parameters.

So the claim "the B∖A₁ solutions form the polar families" is **false**: branch (B)
(and the analogous non-affine families in patterns (II)/(III), e.g. `X=x²∂`,
`D=λx²∂+x⁻¹`) are genuine band-1 `[D,X]=1` pairs in `B∖A₁` that are neither affine
nor polar.

**Why (a) still survives.** Branch (B) has `a₋₁ = const ≠ 0`, so `a₋₁(0)≠0` and
`X∉A₁`. Membership *excludes* branch (B) from the A₁ classification. The memo's
error is conflating the (correct) A₁ statement with a (false) full-`B` statement;
the "back-substituting yields affine" is true only after imposing membership, i.e.
only in branch (A).

## 3. Repaired proof of the load-bearing statement

The clean, arbitrary-degree proof is in
[`quantum-band1.md`](../../../research/band2-square-sector/quantum-band1.md)
(Theorem 1 + Corollary 2). Its two engines:

- **Lemma P** (periodicity ⇒ proportionality) — the memo's `m=±2` step, made
  rigorous: `r(E+1)=r(E)` for a rational `r` ⇒ `r` const, so `b₁=λa₁` for **any**
  nonzero `a₁`. (Machine-checked, `..._engine.py` §C.)
- **Lemma D** (degree drop of the twisted Wronskian) — the missing tool. For
  nonzero `f,g` of degrees `p,q`,
  `W(f,g)=g(E+1)f(E)−f(E−1)g(E)` has degree **exactly** `p+q−1` with leading
  coefficient `(p+q)·lc(f)·lc(g)≠0`. Hence `W=1` forces `p+q=1`. (Machine-checked
  for all `p,q≤4`, `..._classification.py` §I1.)

Lemma D closes the degenerate patterns (II)/(III) — which the memo waved through
"via the Fourier automorphism" — **without** any automorphism: the `m=0` equation
there is literally `W=1`, so `deg a₁+deg b₋₁=1`, and membership (`b₋₁(0)=0`) then
forces the nonconstant factor to be `a₁`, i.e. affine. Case (I) reduces (Lemma P +
telescoping) to `V` linear ⇒ branches (A)/(B), and membership kills (B). Case (IV)
is impossible (`[D,X]₀=0`). Full case tree: `..._search.py` §H, `..._classification.py`.

**Corollary 2 (generation).** For `X=ax+b∂+e`, `D=cx+d∂+f`, `ad−bc=1`:
`dX−bD=x+const` and `−cX+aD=∂+const`; with `1=[D,X]` these give `x,∂∈C⟨X,D⟩`, so
`⟨X,D⟩=A₁`. (Machine-checked.) Thus no band-1 A₁-pair generates a proper
subalgebra — claim (c). ∎

## 4. Issue 2 — the Fourier vs. band-reversal mix-up (audit demand #3)

`E=x∂`. Genuine Weyl-Fourier `F: x↦−∂, ∂↦x` gives `F(E)=(−∂)(x)=−(E+1)`, and it
does negate `ad(E)`-degree and preserve the symmetric band `{−1,0,1}` — **inside
`A₁`**. But `F(x)=−∂` is **not invertible in `B=A₁[x⁻¹]`**, so `F(B)=A₁[∂⁻¹]≠B`:
`F` cannot be used to fold degenerate cases that live in `B`. `..._branchB.py` §F
verifies (bounded search) that `F(x)=−∂` has no right-inverse in `B`.

The correct symmetry is the **band-reversal** `σ: x↦x⁻¹, E↦−E`, i.e.
`σ(x^k f(E))=x^{−k}f(−E)`. It respects `Ex=x(E+1)` (checked), is an involution of
`B`, sends band `{−1,0,1}` to itself, and swaps branch (A)↔(B)
(`σ(x²∂+x+1/x)=x−∂+1/x`, branch A). Crucially `σ(x)=x⁻¹∉A₁`, so **σ does not
preserve `A₁`** — which is exactly why it may not be used inside the A₁
classification, and why the two branches are genuinely distinct there (one in `A₁`
after membership, one not). Bookkeeping repaired: `E↦−E` (not `E↦−E−1`) is the
band-reversal used; `E↦−E−1` belongs to the genuine (non-`B`-preserving) Fourier.

## 5. Literature check (audit demand #5) — affects novelty only

- **Dixmier, "Sur les algèbres de Weyl", Bull. SMF 96 (1968) 209–242.** Poses the
  six problems; **Problem 1** is exactly: `[P,Q]=1 ⇒ ⟨P,Q⟩=A₁` (every endomorphism
  is an automorphism) — i.e. the generation claim (c). Introduces the
  `ad`-semisimple/nilpotent element analysis.
- **Joseph, "The Weyl Algebra — Semisimple and Nilpotent Elements", Amer. J. Math.
  97 (1975) 597–615.** Solves Dixmier's Problems 3 and 6; develops the
  `ad`-semisimple machinery underlying the `ad(E)`-grading used here.
- **Bavula–Levandovskyy, "A remark on the Dixmier conjecture", Canad. Math. Bull.
  63 (2020) 6–12** (arXiv:1812.00042). Uses **precisely** this localization
  `B=S⁻¹A₁=K(H)[X,X⁻¹;σ]`, `σ(H)=H−1`, `Z`-graded by `ad(H)` (`H=pq=E`), and proves
  Dixmier's Problem 1 for `P,Q` each a **sum of ≤2 homogeneous components**, no
  degree bound. Their `p=q` step ("`Q=λP_p`; then `(P,Q−λP)` is obtained from
  `(Y,X)` by automorphisms") is the memo's `b₁=λa₁` step.
- **Han–Tan, "Some progress in the Dixmier conjecture for A₁", Comm. Algebra 52
  (2024).** `spec(ad(pq))=Z`; A₁ is the `ad(pq)`-graded algebra (= our band
  grading); proves generation when an element has **no negative-support
  component**, and extends the ≤2-component result.
- **Guccione–Guccione–Valqui (2025), Comm. Algebra 53(3) 1307**, "Number of
  homogeneous components of counterexamples to the Dixmier conjecture" — directly
  about the band/component count.

**Assessment.** Band-1 = "supports in `ad(E)`-levels `{−1,0,1}`" = sum of ≤3
homogeneous components. Bavula–Levandovskyy already prove generation for ≤2
components; Han–Tan push further. Theorem 1/Cor 2 here are a ≤3-component instance
and **should not be framed as new mathematics** — only as a self-contained,
audited induction base recorded in the program's own crossed-product notation. The
"difference deformation of JC2" framing (QK vs CK) is the memo's genuine
contribution; the band-1 rigidity fact itself is prior art.

## 6. Reproduce

```sh
cd archive-import/provisional/dixmier-band-program
for s in engine branchB search classification exhaustive; do
  uv run --with sympy python audit_band1_$s.py
done
```

Expected: each ends `... CHECKS PASSED`. (The deg-3 ideal-saturation in
`audit_band1_exhaustive.py` §J1b takes ~1 min.)
