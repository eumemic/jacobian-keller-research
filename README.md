# Jacobian, Keller, and JC2 research materials

This staging repository separates the current conservative manuscript from provisional research notes. Inclusion here is for review and reproducibility; it is not a claim of peer review, publication, priority, or mathematical validation.

## Status

| Material | Location | Status |
|---|---|---|
| Telescoping Keller maps manuscript | `paper/` | Current conservative manuscript; author and release metadata remain unresolved. The included PDF was built from the included source before staging. |
| Direct JC2 attack memo | `research/jc2-attack-memo.md` | Provisional, non-peer-reviewed research directions; no claimed result. |
| JC2 frontier source audit | `research/jc2-frontier-audit.md` | Time-bounded literature audit with qualifications and primary-source links. |
| Classical band-2 non-square theorem | `research/band2-classical-proved/` | Proved for polynomial Keller pairs over an algebraically closed characteristic-zero field; independent derivation and exact identity checker included. The square sector remains open. |
| Stronger classical/quantum J3 claim | `research/band2-j3-provisional/` | **PROVISIONAL — INDEPENDENT AUDIT PENDING.** Not an established theorem or claimed result. |
| Imported session artifacts | `archive-import/` | Deduplicated historical/provisional materials with a disposition manifest; not current results. |

## Verification

Requirements: Python 3 and SymPy.

Run the manuscript checks from the repository root:

```sh
python3 paper/verify.py
```

A successful run ends with `ALL CHECKS PASSED`.

Run the independently derived classical band-2 identity checks:

```sh
python3 research/band2-classical-proved/verify_classical_band2.py
```

A successful run ends with `ALL CLASSICAL BAND-2 CHECKS PASSED`.

Run the stronger provisional J3 symbolic checks separately:

```sh
python3 research/band2-j3-provisional/verify_J3.py
```

A successful run ends with `ALL J3 CHECKS PASSED`.

These scripts use exact symbolic arithmetic for the identities they test. Passing them provides computational support only; it is not a proof, proof review, peer review, or independent validation of the surrounding mathematical arguments.

To rebuild the manuscript PDF (with a TeX installation providing the packages used by `main.tex`):

```sh
cd paper
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

## Attribution and release metadata

The manuscript itself contains cautious attribution and does not assert exact priority for overlapping work. Its author name, affiliation, email, acknowledgments, and date are still placeholders and must be resolved before publication. The provisional materials should not be cited as established results without independent review.

No license file is included because no explicit license was present in the supplied source materials. A repository owner must choose and approve a license before public release; absent a license, normal copyright restrictions apply.
