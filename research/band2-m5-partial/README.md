# Band-2 M5 classical theorem and quantum partial cascade package

**INDEPENDENTLY DERIVED AND AUDITED — NOT PEER REVIEWED — SECTOR-SCOPED RESULTS**

This package records independently derived results for a specifically oriented classical square sector and a quantum shifted-square sector related to, but not exhaustively produced by, the earlier J3 obstruction. It also preserves an adversarial audit of the supplied M5 checker.

## Status and scope

- The classical memo now gives a complete arbitrary-degree classification of the oriented `a2`-square sector `a2 = h^2 != 0` when both entries' ladder supports are contained in `[-2,2]` and the genuine polynomial membership conditions hold. After its stated normalizations, every member has the displayed parameterized normal form (0.1).
- This classical result is sector-scoped by itself and does not establish orientation or case coverage beyond the oriented `a2`-square sector. The later [`../band2-classical-full/`](../band2-classical-full/) package supplies the independent band-1 and orientation bridge and assembles this result with M4 into the full support-contained classical band-2 theorem over `C`.
- The quantum shifted-square-sector memo remains partial. Its arbitrary-degree constant-`h` branch with nonzero positive and negative homogeneous parameters remains unresolved.
- The checker audit preserves the supplied checker's sign defects and fatal completeness gaps. The old checker did not prove the classical classification; the new classical proof closes its resistant branch independently, without bounded-degree computational exhaustion.
- No classification checker has been added, and completeness is not inferred from symbolic checks or computation.
- This sector memo proves neither JC2 nor DC1. Its use in the later full classical band-2 assembly still does not prove JC2 because that theorem retains the fixed support bound. No counterexample to JC2/DC1 or to the unresolved quantum classification was found.

## Contents

- [`classical-square-sector-partial.md`](classical-square-sector-partial.md): complete proof and classification of the oriented classical `a2`-square sector under the memo's exact hypotheses. The path is retained for stable links.
- [`quantum-shifted-square-sector-partial.md`](quantum-shifted-square-sector-partial.md): independent quantum finite-difference cascade and unresolved nonlinear tail.
- [`checker-audit.md`](checker-audit.md): historical independent audit of the supplied checker, with a current note distinguishing that checker's failure from the later independent classical proof.

These notes preserve distinct bracket/commutator conventions stated within each memo. Claims should be read under the convention and hypotheses of the relevant file.
