# Band-2 square sector: classical closure and a quantum branch exclusion

**INDEPENDENTLY DERIVED — HOSTILE-AUDITED IN PART — EXACT CHECKS — NOT PEER REVIEWED**

This directory contains M5 results downstream of the audited partial cascades
as frozen at commit `91a053a`. The classical memo closes the historical
classical branch A\* and classifies the classical square sector. The quantum
memo proves emptiness only for the distinct gauged branch
`h = 1, κμs ≠ 0`; it does not classify the full constant-`h` quantum sector.
The tame catalog concerns the classical A\* branch only and proves no quantum
tame exclusion. In parallel, commit `44c66d5` independently rewrote the
classical partial into a complete direct-substitution proof with the same
conclusion.

## Results

| Result | File | Status |
|---|---|---|
| Classical branch A\* is empty (mod-3 degree obstruction on two new first integrals); constant-h branches B0 = exactly tame, A0 empty; nonconstant h excluded by a one-line membership argument; **complete classical square-sector classification** ("every band-2 polynomial Keller pair with a₂ ≠ 0 is a tame automorphism pair") | `classical-Astar.md` | Proof written; identities machine-verified (`verify_classical_Astar.py`); degree arguments are prose proofs; awaiting hostile audit |
| Gauged quantum branch (`h = 1, κμs ≠ 0`) is empty at arbitrary degree via three exact difference integrals + Lemma R; **no DC1 counterexample in this branch** | `quantum-mirror.md` | Narrow structural proof passed hostile audit after cancellation-aware exposition and scope repairs; verifier emits 49 exact `PASS` assertions plus one aggregate bounded-grid check; not peer reviewed |
| Classical tame catalog: the generic single-quadratic-shear normal form lands in B0 after gauge; a bounded 7,874-pair enumeration found 1,120 square-sector examples, all B0 | `tame-catalog.md` | Exact symbolic normal-form identity + bounded enumeration (`verify_catalog.py`); no quantum conclusion and no unbounded tame-word classification claimed |
| Independent re-derivation of the load-bearing identities in the two partial memos as frozen at commit `91a053a` | `upstream-verification.md` | Historical-snapshot identity audit; machine-verified |

## Combined statement (classical, complete)

Every polynomial Keller pair `F, G ∈ C[x,ξ]`, `{G,F} = 1`, with both supports
in x-levels `[-2,2]` and `a₂ ≠ 0`, is a tame automorphism pair. (Inputs: the
proved M4 square theorem; the `../band2-m5-partial/` cascade as audited at
commit `91a053a`; the closures in this directory.)

## Combined statement (quantum, partial)

After the stated gauge and positive cascade, the shifted-square branch
`h = 1, κ ≠ 0, μ ≠ 0, s ≠ 0` is empty. No conclusion is drawn here about the
other constant-`h` subbranches, the nonconstant-`h` sector, quantum tame
classification, or DC1.

## Still open

- Quantum: the `κ = 0` branch, the remaining constant-h sub-branches
  (`s ≠ 0, μ = 0`; `s = 0, w ≠ 0`), and the nonconstant-h ("affine") lemma —
  the general-h central telescoping identity recorded in `quantum-mirror.md`
  §7 is the suggested lever.
- Classical and quantum sectors with `a₂ = 0` (top level ≤ 1), which route
  through band-1 analysis — required for the full band-2 JC2/DC1 fragments.
- Independent hostile audit of `classical-Astar.md` remains open. The narrow
  theorem in `quantum-mirror.md` has passed a hostile audit after the
  cancellation-aware and scope repairs recorded here; it remains non-peer-reviewed.

## Verification

```sh
python3 research/band2-square-sector/verify_classical_Astar.py   # ALL CLASSICAL ASTAR CHECKS PASSED
python3 research/band2-square-sector/verify_quantum_mirror.py    # ALL QUANTUM MIRROR CHECKS PASSED
python3 research/band2-square-sector/verify_catalog.py           # ALL CATALOG CHECKS PASSED
```

Scripts check identities and run bounded regression sweeps; the arbitrary-
degree claims rest on the written proofs in the memos.
