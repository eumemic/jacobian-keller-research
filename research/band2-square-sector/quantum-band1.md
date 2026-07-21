# Quantum band-1 rigidity (audited, citable)

**Status: hostile-audited and repaired.** This file states the band-1 facts that
the band-2 assembly may cite, in the form that survives the audit in
[`../../archive-import/provisional/dixmier-band-program/AUDIT-band1.md`](../../archive-import/provisional/dixmier-band-program/AUDIT-band1.md).
It supersedes, for citation purposes, the classification prose of Theorem P3 in
`band1-rigidity-milestone.md` §3 (whose B∖A₁ classification and "Fourier"
bookkeeping were found to be incomplete/mis-stated — see the audit).

> **Do not cite** "the only band-1 solutions in `A₁[x⁻¹]` are affine + polar."
> That is false (there is a non-affine family — see the Remark). **Cite**
> Theorem 1 (A₁-membership ⇒ affine) and Corollary 2 (generation), which are
> proved below at arbitrary Bernstein degree.

## Conventions

`A₁ = C⟨x, ∂⟩`, `[∂,x]=1`. Localize `B := A₁[x⁻¹] = ⊕_{k∈Z} x^k C[E]`,
`E = x∂`, with the crossed-product law

    (x^a f(E))·(x^b g(E)) = x^{a+b} f(E+b) g(E),      ∂ = x⁻¹E.

`ad(E)` acts by `x^k f(E) ↦ k·x^k f(E)`; the "x-level" / "band" grading is the
`ad(E)`-eigenspace grading (identical to Dixmier's `ad(pq)`-grading, `pq=E`).
Membership: `x^{−j} c(E) ∈ A₁ ⇔ E(E−1)⋯(E−j+1) | c(E)`; in particular for
`j=1`, `x⁻¹c(E) ∈ A₁ ⇔ c(0)=0`.

"Band-1" means support (set of nonzero x-levels) `⊆ {−1,0,1}`. Write

    X = x a₁(E) + a₀(E) + x⁻¹ a₋₁(E),   D = x b₁(E) + b₀(E) + x⁻¹ b₋₁(E).

An **affine symplectic pair** is `X = a x + b ∂ + e`, `D = c x + d ∂ + f` with
`ad − bc = 1` (equivalently a₁,a₀ const, a₋₁ = b·E, etc.).

## The five component equations of [D,X]=1

Expanding `[D,X]=DX−XD` by x-level, `[D,X]=1` is equivalent to (`Δf := f(E+1)−f(E)`):

    m=+2:  b₁(E+1)a₁(E) − a₁(E+1)b₁(E) = 0
    m=+1:  a₁·Δb₀ − b₁·Δa₀ = 0
    m= 0:  b₋₁(E+1)a₁(E) − a₁(E−1)b₋₁(E) + b₁(E−1)a₋₁(E) − a₋₁(E+1)b₁(E) = 1
    m=−1:  a₋₁·(b₀(E−1)−b₀(E)) + b₋₁·(a₀(E)−a₀(E−1)) = 0
    m=−2:  b₋₁(E−1)a₋₁(E) − a₋₁(E−1)b₋₁(E) = 0

(Re-derived independently and machine-checked, `audit_band1_engine.py` §B.)

## Two lemmas

**Lemma P (periodicity ⇒ proportionality).** If `a₁ ≠ 0` and
`b₁(E+1)a₁(E)=a₁(E+1)b₁(E)`, then `b₁ = λ a₁` for a unique `λ ∈ C`. *Proof.* The
ratio `r=b₁/a₁ ∈ C(E)` satisfies `r(E+1)=r(E)`; a 1-periodic rational function has
a shift-invariant (hence empty) finite pole set, so `r` is a 1-periodic
polynomial, i.e. a constant. ∎ Likewise `m=−2` gives `b₋₁=μ a₋₁` when `a₋₁≠0`.

**Lemma D (degree drop of the twisted Wronskian).** For nonzero `f (deg p)`,
`g (deg q)`,

    W(f,g)(E) := g(E+1)f(E) − f(E−1)g(E)

has **degree exactly `p+q−1`** with leading coefficient `(p+q)·lc(f)·lc(g) ≠ 0`
(and `W=0` iff `p=q=0`, i.e. both constant). Same for the mirror
`W'(f,g):=f(E−1)g(E)−g(E+1)f(E)`. *Proof.* The top terms `E^{p+q}` cancel; the
`E^{p+q−1}` coefficient computes to `(p+q)lc(f)lc(g)`. (Machine-verified for all
`p,q ≤ 4`, `audit_band1_classification.py` §I1.) ∎

**Corollary.** If `W(f,g)=1` (a nonzero constant) then `p+q=1`: one of `f,g` is a
nonzero constant and the other is exactly linear.

## Theorem 1 (band-1 A₁-rigidity)

*Let `X, D ∈ A₁` (genuine Weyl-algebra elements) with `[D,X]=1`, both band-1.
Then `(X,D)` is an affine symplectic pair.*

**Proof.** Membership at level `−1` reads `a₋₁(0)=0` and `b₋₁(0)=0`. Split on the
top/bottom coefficients.

*Case I: `a₁≠0` and `a₋₁≠0`.* By Lemma P, `b₁=λa₁`, `b₋₁=μa₋₁`. If `λ=μ` then
`D−λX ∈ C[E]` and `[D,X]=[D−λX, X]` has zero level-0 component (a level-0 element
commutes with the level-0 part of `X`), contradicting `[D,X]=1`; so `λ≠μ`. Then
`m=±1` give `Δb₀=λΔa₀` and `Δb₀=μΔa₀`, whence `(λ−μ)Δa₀=0`, so `a₀=α`, `b₀=β`
are constant. Now `m=0` telescopes:

    [D,X]₀ = (λ−μ)(V(E) − V(E+1)) = 1,     V(E) := a₁(E−1)·a₋₁(E),

so `ΔV` is the nonzero constant `−1/(λ−μ)` and `V` is **exactly linear**. Since
`deg V = deg a₁ + deg a₋₁ = 1`, either (A) `a₁` const, `a₋₁` linear, or (B) `a₁`
linear, `a₋₁` const. In case (B), `a₋₁ = const ≠ 0`, so `a₋₁(0) ≠ 0` and `X ∉ A₁`
— excluded. In case (A), `a₁ = A` const and membership `a₋₁(0)=0` forces
`a₋₁ = ε·E`, giving `X = A·x + α + ε∂` — affine; and `D = λA·x + β + με∂` —
affine.

*Case II/III: exactly one of `a₁,a₋₁` is 0.* Say `a₋₁=0` (Case III; `X=xa₁+a₀`).
Then `X∈A₁` automatically and `m=0` is `W(a₁,b₋₁)=1` (with `b₁=λa₁`). By the
Corollary, `deg a₁ + deg b₋₁ = 1`. If `a₁` is linear then `b₋₁` is constant, but
membership `b₋₁(0)=0` forces `b₋₁=0`, contradicting `W=1`; so `a₁ = A` is constant
and `b₋₁ = E/A` (from `AΔb₋₁=1` and `b₋₁(0)=0`) — affine. Case II (`a₁=0`) is
identical with `W'(b₁,a₋₁)=1`.

*Case IV: `a₁=a₋₁=0`.* Then `X=a₀(E)∈C[E]` and `[D,X]₀=0≠1` — impossible.

All cases give an affine symplectic pair. ∎

## Corollary 2 (generation; Dixmier's Problem 1 in band 1)

*Every band-1 pair `X,D ∈ A₁` with `[D,X]=1` generates `A₁`. Hence no such pair
generates a proper subalgebra; equivalently, a pair generating a proper
subalgebra of `A₁` must have a generator with `ad(E)`-support outside `{−1,0,1}`.*

**Proof.** By Theorem 1, `X=ax+b∂+e`, `D=cx+d∂+f`, `ad−bc=1`. Then
`dX−bD = x + (be−df)` and `−cX+aD = ∂ + (ce−af)`, and `1=[D,X]`, so `x,∂` lie in
`C⟨X,D⟩`. ∎ (Machine-checked in `audit_band1_classification.py`/`_exhaustive.py`.)

## Remark (what is NOT true — the correct B∖A₁ picture)

`B=A₁[x⁻¹]` **does** contain non-affine band-1 pairs; they are exactly the
solutions killed by membership in the cases above. Representatives:

    branch (B):     X = x²∂ + x + x⁻¹,        D = x⁻¹            ([D,X]=1)
    pattern (III):  X = x²∂,                  D = λx²∂ + x⁻¹     ([D,X]=1, any λ)

These have `[D,X]=1` and band-1 support but are **not** affine and **not** of the
"polar" form `affine + c·x⁻¹`; each lies in `B∖A₁` (e.g. `x⁻¹ ∉ A₁`). They are
irrelevant to DC1/Dixmier (which concerns genuine `A₁` pairs), which is why
Theorem 1 — not any full-`B` classification — is the citable fact. The
band-reversing automorphism `σ: x↦x⁻¹, E↦−E` (an automorphism of `B`, **not** of
`A₁`; the naive Weyl "Fourier" `x↦−∂` does *not* preserve `B`) swaps branch (A)
and branch (B).

## Provenance / prior art

The setup and the `m=±2 ⇒ b₁=λa₁` step are exactly the localization
`B = S⁻¹A₁ = K(H)[X,X⁻¹;σ]` and the degree case-analysis of
**Bavula–Levandovskyy, "A remark on the Dixmier conjecture", Canad. Math. Bull.
63 (2020) 6–12** (arXiv:1812.00042), who prove Dixmier's Problem 1 for elements
that are sums of ≤2 homogeneous `ad(E)`-components. **Han–Tan, "Some progress in
the Dixmier conjecture for A₁", Comm. Algebra 52 (2024)** extend this and give a
positive-support criterion. Theorem 1 above is the ≤3-component (band-1) instance
and is not claimed as new mathematics; it is recorded here in a self-contained,
audited form for the band-2 induction base.

## Verification

    uv run --with sympy python ../../archive-import/provisional/dixmier-band-program/audit_band1_engine.py
    uv run --with sympy python ../../archive-import/provisional/dixmier-band-program/audit_band1_branchB.py
    uv run --with sympy python ../../archive-import/provisional/dixmier-band-program/audit_band1_search.py
    uv run --with sympy python ../../archive-import/provisional/dixmier-band-program/audit_band1_classification.py
    uv run --with sympy python ../../archive-import/provisional/dixmier-band-program/audit_band1_exhaustive.py
