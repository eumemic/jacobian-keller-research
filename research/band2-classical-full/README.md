# Full classical band-2 theorem package

**INDEPENDENTLY ASSEMBLED AND INTERNALLY AUDITED — NOT PEER REVIEWED — BAND-SCOPED**

This package assembles the proved classical M4 non-square theorem and M5
square-sector classification with an independent band-1 proof and explicit
pair-exchange/reflection orientation lemmas.

## Result

Over `C`, if `F,G in C[x,xi]` satisfy `{G,F}=1` and both have ladder support
contained in `[-2,2]` for `t=x*xi`, including the genuine polynomial membership
conditions at negative ladder degrees, then `(F,G)` is a polynomial
automorphism.

This does not prove JC2: the theorem assumes the fixed band-2 support bound.

## Files

- [`full_classical_band2_theorem.md`](full_classical_band2_theorem.md): exact
  statement, band-1 proof, orientation lemmas, M4+M5 assembly, and inverse.
- [`verify_full_classical_band2.py`](verify_full_classical_band2.py): exact
  SymPy regression checks for the band-1 identities and shear, representative
  reflection mapping/sign examples, and the normal-form bracket/inverse. The
  script is support, not proof completeness.

## Run the verifier

From the repository root:

```sh
python3 research/band2-classical-full/verify_full_classical_band2.py
```

A successful run ends with
`ALL FULL CLASSICAL BAND-2 REGRESSION CHECKS PASSED`.
