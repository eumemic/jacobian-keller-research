# Jacobian, Keller, and JC2 research materials

This living repository separates the current conservative manuscript from provisional research notes. Inclusion here is for review and reproducibility; it is not a claim of peer review, priority, or mathematical validation.

## Status

| Material | Location | Status |
|---|---|---|
| Telescoping Keller maps manuscript | `paper/` | Current conservative manuscript; author and release metadata remain unresolved. The included PDF was built from the included source before staging. |
| Direct JC2 attack memo | `research/jc2-attack-memo.md` | Provisional, non-peer-reviewed research directions; no claimed result. |
| JC2 frontier source audit | `research/jc2-frontier-audit.md` | Time-bounded literature audit with qualifications and primary-source links. |
| Full classical band-2 theorem | [`research/band2-classical-full/`](research/band2-classical-full/) | Independently assembled and internally audited, not peer reviewed. Over `C`, every genuine polynomial Keller pair whose two ladder supports are contained in `[-2,2]` is a polynomial automorphism. The assembly cites the proved M4 and M5 sector arguments, adds independent band-1 and orientation proofs, and does not claim JC2. |
| Classical band-2 M4 non-square theorem | [`research/band2-classical-proved/`](research/band2-classical-proved/) | Proved for polynomial Keller pairs over an algebraically closed characteristic-zero field; independent derivation and exact identity checker included. It supplies the non-square exclusion used by the full classical assembly. |
| Audited classical/quantum band-2 obstruction (J3) | [`research/band2-j3-provisional/`](research/band2-j3-provisional/) | Core Laurent/localized theorem independently audited and repaired; not peer reviewed. With both entries' ladder supports contained in `[-2,2]` and `a2*a-2 != 0` for the designated member, it excludes only the case in which both extremes lie outside the relevant (shifted-)square classes; it does not classify band 2 or prove JC2/DC1. |
| Classical band-2 M5 sector theorem / quantum partial cascade | [`research/band2-m5-partial/`](research/band2-m5-partial/) | Independently derived and audited, not peer reviewed. The oriented classical `a2`-square sector is completely classified at arbitrary degree when both entries' ladder supports are contained in `[-2,2]` and genuine membership holds; it supplies the square-sector result used by the full classical assembly. Quantum remains partial; no DC1 result is claimed. |
| Imported session artifacts | `archive-import/` | Deduplicated historical/provisional materials with a disposition manifest; not current results. |

## Verification

Requirements: Python 3 and SymPy.

Run the manuscript checks from the repository root:

```sh
python3 paper/verify.py
```

A successful run ends with `ALL CHECKS PASSED`.

Run the full classical band-2 assembly's exact regression checks:

```sh
python3 research/band2-classical-full/verify_full_classical_band2.py
```

A successful run ends with `ALL FULL CLASSICAL BAND-2 REGRESSION CHECKS PASSED`.

Run the underlying M4 classical identity checks separately:

```sh
python3 research/band2-classical-proved/verify_classical_band2.py
```

A successful run ends with `ALL CLASSICAL BAND-2 CHECKS PASSED`.

Run the audited J3 package's selected symbolic checks separately:

```sh
python3 research/band2-j3-provisional/verify_J3.py
```

A successful run ends with `ALL J3 CHECKS PASSED`.

These scripts use exact symbolic arithmetic for the identities they test. Passing them provides computational support only; it is not by itself a proof, proof review, peer review, or completeness verification. The full classical theorem's completeness comes from its written band-1/orientation split and the cited arbitrary-degree M4 and M5 proofs, not from its regression script. The J3 written argument has undergone an independent theorem audit, while remaining non-peer-reviewed. The M5 package's classical result comes from its written proof rather than bounded-degree computation. The quantum cascade remains partial.

To rebuild the manuscript PDF (with a TeX installation providing the packages used by `main.tex`):

```sh
cd paper
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

## Attribution and release metadata

The manuscript itself contains cautious attribution and does not assert exact priority for overlapping work. Its author name, affiliation, email, acknowledgments, and date are still placeholders and must be resolved before publication. The provisional materials should not be cited as established results without independent review.

No license file is included because no explicit license was present in the supplied source materials. A repository owner must choose and approve a license before public release; absent a license, normal copyright restrictions apply.
