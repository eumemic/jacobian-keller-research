# The full quantum band-2 theorem

**ASSEMBLY OF INDEPENDENTLY PROVED AND MACHINE-VERIFIED COMPONENTS — NOT PEER REVIEWED**

This memo assembles the quantum band-2 program into a single theorem, mirroring
the classical assembly `../band2-classical-full/full_classical_band2_theorem.md`
(commit `f637b1a`) one floor up the quantization ladder. Every component is
proved in its own memo with its own exact verifier; this document contributes
only the composition, the orientation bookkeeping, and the tameness of the
undoing operations. The driver `verify_quantum_band2_all.py` re-runs the entire
verification chain and checks the assembly-specific facts; a successful run
ends `ALL QUANTUM BAND-2 ASSEMBLY CHECKS PASSED`.

## Theorem and exact scope

Work in the first Weyl algebra `A₁ = C⟨x, ∂⟩`, `[∂, x] = 1`, with `E = x∂` and
the ladder presentation

```
X = Σ_{k=-2}^{2} x^k a_k(E),      D = Σ_{l=-2}^{2} x^l b_l(E),
```

where genuine `A₁`-membership means the falling-factorial divisibilities
`E | a₋₁, b₋₁` and `E(E−1) | a₋₂, b₋₂`. Suppose `[D, X] = 1` and both ladder
supports are contained in `[−2, 2]` (an extreme coefficient may vanish).

> **Theorem.** `(X, D)` is the image of `(x, ∂)` under a **tame automorphism**
> of `A₁` (a composition of affine symplectic automorphisms and triangular
> automorphisms `exp(ad p(x))`, `exp(ad q(∂))`). In particular
> `C⟨X, D⟩ = A₁`: no such pair generates a proper subalgebra.

This is the band-2 case of the Dixmier conjecture for `A₁` in the fixed
`E`-ladder band filtration. It does **not** claim DC1 itself: a general Weyl
pair need not admit this band bound, exactly as the classical counterpart
(`f637b1a`) does not claim JC2.

## Components consumed (each pinned and separately verified)

| # | Statement | Source | Verifier |
|---|---|---|---|
| C1 | Band-1: every `A₁` band-1 pair is affine symplectic and generates `A₁` (audited repair of provisional P3; mild instance of Bavula–Levandovskyy 2020, Han–Tan 2024) | `quantum-band1.md` (this dir); audit `../../archive-import/provisional/dixmier-band-program/AUDIT-band1.md` | `audit_band1_*.py` |
| C2 | `a₂ = 0` reduction: the sector partition Z / O1 / O2 / O3 is exhaustive; Z collapses to band-1; O1–O3 route via pair exchange and/or the Fourier automorphism `φ` to an equivalent genuine `A₁` pair with top coefficient nonzero | `quantum-a2-zero.md` | `verify_quantum_a2zero.py` |
| C3 | M4q: a genuine band-2 pair with `a₂ ≠ 0` has `a₂ = c·h(E)h(E+1)` a shifted square (no hypothesis on `a₋₂`) | `quantum-M4.md` | `verify_quantum_M4.py` |
| C4 | Shifted-square sector: `h` is forced constant, `κ ≠ 0`, `b₋₂ = 0`, and the pair is a member of the tame family `U = x + c₀ + c₁∂`, `X = U² − ∂/κ − A`, `D = λX + κU + β` (resistant branch: `quantum-mirror.md`, commit `ad43ab5`, repaired `9e70871`, audit-endorsed; all remaining branches: `quantum-completion.md`) | `quantum-mirror.md` + `quantum-completion.md` | `verify_quantum_mirror.py`, `verify_quantum_completion.py` |

## Assembly

**Step 0 (band 1).** If all four extremes `a₂, a₋₂, b₂, b₋₂` vanish, both
members are band-1 and C1 gives an affine symplectic pair — tame, generating.

**Step 1 (orientation).** Otherwise, by C2's orientation lemma, applying
pair exchange `(X, D) ↦ (D, −X)` and/or the Fourier automorphism
`φ: x ↦ −∂, ∂ ↦ x` (a genuine automorphism of `A₁`, so `[φD, φX] = φ([D,X])`
exactly) produces an equivalent genuine band-2 pair with top coefficient
nonzero. Both operations preserve `[·,·] = 1`, membership, and the band.

**Step 2 (square class).** C3 applies to the oriented pair: its top is
`c·h(E)h(E+1)`. Over `C`, absorb `c` by `h ↦ √c·h`, so `a₂ = h h⁽¹⁾`.

**Step 3 (classification).** C4 classifies the oriented pair: after the
reversible normalizations (gauge `D ↦ D − λX`, diagonal scaling, additive
constants) it is the displayed tame family, which is the image of `(x, ∂)`
under an explicit composition of affine symplectic and triangular
automorphisms; C4's verifiers exhibit the generation certificate
(`x` and `∂` recovered from `X, D` by explicit polynomial expressions).

**Step 4 (undoing).** The operations of Steps 1–3 must themselves be tame:

- gauge, scaling, additive constants: affine symplectic / triangular — tame
  by definition;
- pair exchange: `(D, −X)` is the image of `(X, D)` under composition with
  `φ⁻¹` acting on the pair `(x, ∂) ↦ (−∂, x)`… precisely, if `(D, −X)` is the
  image of `(x, ∂)` under a tame `ψ`, then `(X, D) = (ψ∘φ)(x, ∂ )` up to the
  sign conventions checked in the driver;
- the Fourier automorphism `φ` is itself tame, via the standard triangular
  factorization
  `φ = exp(ad(−x²/2)) ∘ exp(ad(−∂²/2)) ∘ exp(ad(−x²/2))`
  (the SL₂ identity shear(1)·shear(−1)·shear(1) = rotation; checked exactly
  on generators in the driver).

Composing, the original `(X, D)` is the image of `(x, ∂)` under a tame
automorphism of `A₁`, and generation follows. ∎

## What is new here versus the components

Only Step 4's tameness bookkeeping (in particular the exact triangular
factorization of `φ` and the pair-exchange sign conventions) and the
end-to-end composition. Everything else is a citation to a separately
verified memo.

## Relation to the classical theorem

The statement, proof shape, and even the normal form deform exactly:
classical `(0.1)` in `f637b1a` is the `ħ → 0` limit of the C4 family, the
classical orientation lemma's determinant-(−1) reflection is replaced by the
genuinely automorphic `φ` (no sign flip), and the classical membership
divisor `τ^j` deforms to the falling factorial `E(E−1)⋯(E−j+1)` — the
mechanism that closes M4q where polynomiality closed M4. The band-2 floor of
the width induction now stands on both faces:

- **JC2 face:** every band-2 Keller pair is a tame polynomial automorphism
  pair (`f637b1a`; independent square-sector proof in `classical-Astar.md`).
- **DC1 face:** every band-2 Weyl pair is a tame automorphism image
  (this memo).

## Verification

```sh
python3 research/band2-square-sector/verify_quantum_band2_all.py
```

runs, in order: the assembly-specific checks (φ is an automorphism; the
triangular factorization of φ on generators; pair-exchange bookkeeping; a
random tame round-trip), then the component verifiers
(`verify_quantum_a2zero.py`, `verify_quantum_M4.py`,
`verify_quantum_mirror.py`, `verify_quantum_completion.py`, and the band-1
audit scripts), failing loudly if any sub-verifier fails.

## Status and residuals

- Not peer reviewed; the components carry their own audit statuses
  (`quantum-mirror.md` audit-endorsed and repaired upstream at `9e70871`;
  band-1 audited CONFIRMED WITH REPAIRS; `quantum-completion.md`,
  `quantum-M4.md`, `quantum-a2-zero.md`, and this assembly await their
  hostile passes).
- The `B∖A₁` (localized/polar) landscape is deliberately out of scope —
  the audit of P3 showed the localized band-1 classification is richer than
  the retracted "affine + polar" claim, and no statement here depends on it.
- Next floor of the induction: band 3, where the cross-coupling of
  proportionality constants (flagged in the archived milestone) is the first
  genuinely new phenomenon.
