# Band-2 J3/M4 audited package

**INDEPENDENTLY AUDITED — NOT PEER REVIEWED**

This directory contains a conservative statement and self-contained proof of a narrow classical Laurent and quantum localized obstruction. The independent audit found the core argument correct after repairs now incorporated in the theorem. The result requires both entries to have ladder support contained in `[-2,2]`, requires both extreme coefficients of `X` to be nonzero, and proves only that those two extremes cannot both be non-square classically or non-shifted-square quantumly.

It does not cover one-sided sectors, show that each extreme is square, classify band-2 pairs, establish invertibility, or prove JC2/DC1. Corollaries for genuine polynomial Keller and Weyl pairs require membership in the unlocalized algebra and the same common band presentation.

## Files

- `theorem_J3_band2_nonsquare.md`: audited statement, conventions, self-contained essential lemmas, proofs, and scope limits.
- `verify_J3.py`: 13 exact SymPy checks of selected reduced identities and finite examples.

## Run the checker

From the repository root:

```sh
python3 research/band2-j3-provisional/verify_J3.py
```

A successful run ends with `ALL J3 CHECKS PASSED`. The checks provide computational support for selected reductions only. They do not test the full unreduced systems, prove the arbitrary-degree lemmas, establish branch completeness, check localization or Weyl membership, or certify proof completeness.
