# Band-2 J3/M4 provisional package

**PROVISIONAL — INDEPENDENT AUDIT PENDING — NOT PEER REVIEWED**

This directory stages a supplied argument and its symbolic checker for inspection. Inclusion does not assert that the argument is correct, complete, novel, or publishable, and it must not be represented as an established theorem or a proof of JC2 or DC1.

The argument refers to earlier items labeled `Lemma J2`, `Lemma J2q`, and `J1/P3`, but those dependencies were not present in the supplied archive. The package is therefore not self-contained; this is an explicit audit blocker rather than something repaired in staging.

## Files

- `theorem_J3_band2_nonsquare.md`: provisional mathematical argument, lightly edited only to make its audit status explicit and to avoid presenting a claimed consequence as established.
- `verify_J3.py`: exact SymPy checks for selected identities used by the argument.

## Run the checker

From the repository root:

```sh
python3 research/band2-j3-provisional/verify_J3.py
```

A successful run ends with `ALL J3 CHECKS PASSED`. Passing the checker is computational support only. It does not establish the general lemmas invoked by the prose, certify every inference, replace proof review, or constitute peer review or independent validation.
