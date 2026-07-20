# Band-2 square sector: resistant-branch closures and tame catalog

**INDEPENDENTLY DERIVED — MACHINE-CHECKED — NOT PEER REVIEWED — NOT YET CROSS-AUDITED**

This directory contains the session's M5 results downstream of the audited
partial cascades in `../band2-m5-partial/` (as of commit `91a053a`). Together
they close **both resistant branches** that package left open, classify the
classical square sector completely, and prove no tame pair ever inhabited the
resistant loci. Note: in parallel, upstream commit `44c66d5` rewrote the
classical partial into its own complete proof of the classical square sector
by direct substitution — independent of `classical-Astar.md` here, with
agreeing conclusions; the quantum closure and the tame-catalog results in this
directory are not covered by that rewrite.

## Results

| Result | File | Status |
|---|---|---|
| Classical branch A\* is empty (mod-3 degree obstruction on two new first integrals); constant-h branches B0 = exactly tame, A0 empty; nonconstant h excluded by a one-line membership argument; **complete classical square-sector classification** ("every band-2 polynomial Keller pair with a₂ ≠ 0 is a tame automorphism pair") | `classical-Astar.md` | Proof written; identities machine-verified (`verify_classical_Astar.py`); degree arguments are prose proofs; awaiting hostile audit |
| Quantum resistant branch (h = 1, κμs ≠ 0) is empty at arbitrary degree via three exact difference integrals + Lemma R (a staggered-shift leading-coefficient rigidity with no differential analogue); **no DC1 counterexample in this branch** | `quantum-mirror.md` | Complete structural proof, machine-verified (`verify_quantum_mirror.py`, 50 checks); awaiting hostile audit |
| No tame pair lands in A\* or its quantum mirror: both extreme weights of any band-2 tame pair originate in a single quadratic shear, so the b₂-killing gauge kills b₋₂ simultaneously; 7,874-pair enumeration: 1,120 square-sector tames, all branch B0. Hence any A\* inhabitant would have been a JC2 counterexample outright — the emptiness proofs above rule that out | `tame-catalog.md` | Symbolic mechanism proof (all parameters) + bounded enumeration (`verify_catalog.py`) |
| Independent re-derivation of every load-bearing identity in both upstream partial memos: no defects found | `upstream-verification.md` | Machine-verified |

## Combined statement (classical, complete)

Every polynomial Keller pair `F, G ∈ C[x,ξ]`, `{G,F} = 1`, with both supports
in x-levels `[-2,2]` and `a₂ ≠ 0`, is a tame automorphism pair. (Inputs: the
proved M4 square theorem; the audited `../band2-m5-partial/` cascade; the
closures in this directory.)

## Combined statement (quantum, partial)

In the shifted-square sector with constant `h` (`a₂ = 1` after scaling), the
resistant branch `κ ≠ 0, μ ≠ 0, s ≠ 0` is empty; the surviving constant-h
families are the tame swapped-shear and Z-generator families.

## Still open

- Quantum: the `κ = 0` branch, the remaining constant-h sub-branches
  (`s ≠ 0, μ = 0`; `s = 0, w ≠ 0`), and the nonconstant-h ("affine") lemma —
  the general-h central telescoping identity recorded in `quantum-mirror.md`
  §7 is the suggested lever.
- Classical and quantum sectors with `a₂ = 0` (top level ≤ 1), which route
  through band-1 analysis — required for the full band-2 JC2/DC1 fragments.
- Hostile cross-audit of `classical-Astar.md` and `quantum-mirror.md`
  (the tame-catalog agent ran both verifiers read-only and found them
  consistent, but no dedicated referee pass has happened yet).

## Verification

```sh
python3 research/band2-square-sector/verify_classical_Astar.py   # ALL CLASSICAL ASTAR CHECKS PASSED
python3 research/band2-square-sector/verify_quantum_mirror.py    # ALL QUANTUM MIRROR CHECKS PASSED
python3 research/band2-square-sector/verify_catalog.py           # ALL CATALOG CHECKS PASSED
```

Scripts check identities and run bounded regression sweeps; the arbitrary-
degree claims rest on the written proofs in the memos.
