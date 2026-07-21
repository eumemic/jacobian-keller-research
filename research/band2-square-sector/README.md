# Band-2 square sector: classical closure and the quantum band-2 assembly

**INDEPENDENTLY DERIVED — REPAIRED AFTER HOSTILE AUDIT — EXACT CHECKS — NOT PEER REVIEWED**

This directory contains M5 results downstream of the audited partial cascades
as frozen at commit `91a053a`. The classical memo closes the historical
classical branch A\* and classifies the classical square sector (in parallel,
commit `44c66d5` independently rewrote the classical partial into a complete
direct-substitution proof with the same conclusion). The quantum side has
since been carried to a full fixed-band-2 assembly and repaired after a hostile
audit. Exact scripts support identities and bounded regressions; arbitrary-
degree completeness comes from the written proofs. This is not peer review and
does not claim DC1.

## Results

| Result | File | Status |
|---|---|---|
| Classical branch A\* is empty (mod-3 degree obstruction on two new first integrals); constant-h branches B0 = exactly tame, A0 empty; nonconstant h excluded by a one-line membership argument; **complete classical square-sector classification** | `classical-Astar.md` | Proof written; identities machine-verified; degree arguments are prose proofs; awaiting hostile audit |
| Gauged quantum branch (`h = 1, κμs ≠ 0`) is empty at arbitrary degree via three exact difference integrals + Lemma R; **no DC1 counterexample in this branch** | `quantum-mirror.md` | Narrow structural proof passed hostile audit after cancellation-aware exposition and scope repairs (`9e70871`); not peer reviewed |
| Remaining shifted-square branches, including the exact base branch; explicit recovery identities and tame word ⇒ **shifted-square sector = tame family** | `quantum-completion.md` | Written arbitrary-degree proof repaired after hostile audit; exact identities plus clearly labeled bounded Gröbner regressions |
| **M4q**: a genuine band-2 Weyl pair with `a₂ ≠ 0` has `a₂` a shifted square (no hypothesis on `a₋₂`; falling-factorial membership endgame) | `quantum-M4.md` | Written proof with exact symbolic checks; repaired package, not peer reviewed |
| `a₂ = 0` reduction: exhaustive Z/O1/O2/O3 partition; Z closed via band 1; O1–O3 routed by pair exchange / genuine Fourier automorphism φ | `quantum-a2-zero.md` | Retrospective status repaired; exact orientation checks, written partition supplies completeness |
| Audited band-1 statement (Theorem 1 + Cor. 2): `A₁` band-1 pairs are affine symplectic and generate `A₁`; repaired from provisional P3, whose localized "affine + polar" classification is retracted (counterexample `X = x²∂ + x + x⁻¹, D = x⁻¹`); essentially prior art (Bavula–Levandovskyy 2020, Han–Tan 2024) | `quantum-band1.md` | Hostile-audited: CONFIRMED WITH REPAIRS (`../../archive-import/provisional/dixmier-band-program/AUDIT-band1.md`) |
| **Assembly — the full fixed quantum band-2 theorem**: orientation → M4q → shifted-square classification → explicit tame undoing | `quantum-band2-theorem.md` | Repaired after hostile audit; written proof supplies completeness, driver reruns exact and bounded checks |
| Classical tame catalog: the generic single-quadratic-shear normal form lands in B0 after gauge; bounded 7,874-pair enumeration, 1,120 square-sector examples, all B0 | `tame-catalog.md` | Exact symbolic normal-form identity + bounded enumeration; no quantum conclusion and no unbounded tame-word classification claimed |
| Independent re-derivation of the load-bearing identities in the two partial memos as frozen at `91a053a` | `upstream-verification.md` | Historical-snapshot identity audit; machine-verified |

## Combined statement (classical, complete)

Every polynomial Keller pair `F, G ∈ C[x,ξ]`, `{G,F} = 1`, with both supports
in x-levels `[-2,2]` and `a₂ ≠ 0`, is a tame automorphism pair. (Inputs: the
proved M4 square theorem; the `../band2-m5-partial/` cascade as audited at
commit `91a053a`; the closures in this directory. The full classical band-2
theorem, including `a₂ = 0`, is `../band2-classical-full/` at `f637b1a`.)

## Combined statement (quantum, repaired fixed-band theorem)

Every genuine `A₁` pair `[D, X] = 1` with both ladder supports in `[−2, 2]`
is the image of `(x, ∂)` under a tame automorphism of `A₁`; in particular it
generates `A₁`. This is the band-2 case of the Dixmier conjecture in the band
filtration (see `quantum-band2-theorem.md` for exact scope). The band bound is
essential: this neither proves nor claims DC1. The package remains internal and
not peer reviewed despite the hostile-audit repairs.

## Still open

- Independent external review of the repaired fixed-band theorem.
- The localized (`B∖A₁`) polar landscape — richer than previously claimed
  (see the band-1 audit); deliberately out of scope for the theorems above.
- Band 3: cross-coupled proportionality constants — the next floor of the
  width induction on both faces.

## Verification

```sh
python3 research/band2-square-sector/verify_quantum_band2_all.py   # full quantum chain, end to end
python3 research/band2-square-sector/verify_classical_Astar.py     # ALL CLASSICAL ASTAR CHECKS PASSED
python3 research/band2-square-sector/verify_catalog.py             # ALL CATALOG CHECKS PASSED
```

Scripts check identities and run bounded regression sweeps; the arbitrary-
degree claims rest on the written proofs in the memos.
