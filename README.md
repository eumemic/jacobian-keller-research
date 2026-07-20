# Jacobian, Keller, and JC2 research materials

This living repository separates the current conservative manuscript from provisional research notes. Inclusion here is for review and reproducibility; it is not a claim of peer review, priority, or mathematical validation.

## Status

| Material | Location | Status |
|---|---|---|
| Telescoping Keller maps manuscript | `paper/` | Current conservative manuscript; author and release metadata remain unresolved. The included PDF was built from the included source before staging. |
| Direct JC2 attack memo | `research/jc2-attack-memo.md` | Provisional, non-peer-reviewed research directions; no claimed result. |
| JC2 frontier source audit | `research/jc2-frontier-audit.md` | Time-bounded literature audit with qualifications and primary-source links. |
| Classical band-2 non-square theorem | `research/band2-classical-proved/` | Proved for polynomial Keller pairs over an algebraically closed characteristic-zero field; independent derivation and exact identity checker included. An oriented `a2`-square sector is classified separately under the M5 memo's exact support-containment and membership hypotheses; full orientation/case coverage remains open. |
| Audited classical/quantum band-2 obstruction (J3) | [`research/band2-j3-provisional/`](research/band2-j3-provisional/) | Core Laurent/localized theorem independently audited and repaired; not peer reviewed. With both entries' ladder supports contained in `[-2,2]` and `a2*a-2 != 0` for the designated member, it excludes only the case in which both extremes lie outside the relevant (shifted-)square classes; it does not classify band 2 or prove JC2/DC1. |
| Classical band-2 M5 sector theorem / quantum partial cascade | [`research/band2-m5-partial/`](research/band2-m5-partial/) | Independently derived and audited, not peer reviewed. The oriented classical `a2`-square sector is completely classified at arbitrary degree when both entries' ladder supports are contained in `[-2,2]` and genuine membership holds; this is not full orientation/case coverage. Quantum remains partial. No full band-2 theorem or proof of JC2/DC1; the flawed supplied checker is not the source of the classical proof. |
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

Run the audited J3 package's selected symbolic checks separately:

```sh
python3 research/band2-j3-provisional/verify_J3.py
```

A successful run ends with `ALL J3 CHECKS PASSED`.

These scripts use exact symbolic arithmetic for the identities they test. Passing them provides computational support only; it is not by itself a proof, proof review, peer review, or completeness verification. The J3 written argument has undergone an independent theorem audit, while remaining non-peer-reviewed. The M5 package deliberately includes no classification checker: its checker audit explains why the supplied script could not certify arbitrary-degree exhaustion, while the later oriented classical `a2`-square theorem closes its stated sector by an independent written proof. The quantum cascade remains partial.

To rebuild the manuscript PDF (with a TeX installation providing the packages used by `main.tex`):

```sh
cd paper
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

## Attribution and release metadata

The manuscript itself contains cautious attribution and does not assert exact priority for overlapping work. Its author name, affiliation, email, acknowledgments, and date are still placeholders and must be resolved before publication. The provisional materials should not be cited as established results without independent review.

No license file is included because no explicit license was present in the supplied source materials. A repository owner must choose and approve a license before public release; absent a license, normal copyright restrictions apply.
