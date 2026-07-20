# Band-2 M5 partial cascade package

**INDEPENDENTLY DERIVED AND AUDITED — NOT PEER REVIEWED — PARTIAL RESULTS ONLY**

This package records independently derived classical and quantum results for the remaining square/shifted-square sectors after the earlier J3 obstruction. It also records an adversarial audit of the supplied M5 checker. The derivations establish useful arbitrary-degree cascade reductions, but they do not close the resistant branches and therefore do not classify band 2.

## Status and scope

- The classical square-sector memo proves top and central cascade identities, including an arbitrary-degree divisibility consequence, and gives an exact residual negative system. Its resistant arbitrary-degree branch with both negative extremes surviving remains unresolved.
- The quantum shifted-square-sector memo proves the top homogeneous parameter, a ladder-2 divisibility result for the nonzero-parameter branch, the quantum midpoint equation, and an exact reduced finite-difference tail. The arbitrary-degree constant-`h` branch with nonzero positive and negative homogeneous parameters remains unresolved.
- The checker audit finds useful exact regression identities and examples, but also sign defects and fatal completeness gaps. The supplied checker is not included here and is not endorsed as a classification certificate.
- Exact symbolic checks, bounded-degree searches, and displayed examples support only the identities or examples tested. They do not prove exhaustion or completeness.
- No counterexample was found.

In particular, this package proves neither a full classical nor a full quantum band-2 classification. It does not prove JC2 or DC1. The classical and quantum resistant arbitrary-degree branches described in the memos remain open.

## Contents

- [`classical-square-sector-partial.md`](classical-square-sector-partial.md): independent classical derivation, exact residual branch split, proved family-specific membership conclusions, and unresolved arbitrary-degree branch.
- [`quantum-shifted-square-sector-partial.md`](quantum-shifted-square-sector-partial.md): independent quantum finite-difference cascade and unresolved nonlinear tail.
- [`checker-audit.md`](checker-audit.md): independent audit explaining why the supplied checker cannot certify classification or completeness.

These notes preserve distinct bracket/commutator conventions stated within each memo. Claims should be read under the convention and hypotheses of the relevant file.
