# Archive import manifest

## Review scope

This import was reviewed from a 14-file research-session ZIP. The ZIP contained 14 regular files, no directories, no hidden entries, no symbolic links, and no absolute or parent-traversal paths. Total uncompressed size: 445,149 bytes. Archive SHA-256: `9a58b3254c1382585460f0dfc10ecd8959ba49d4f1a351106e15803cdddf086a`.

Only sanitized, nonduplicate material with continuing research value is staged here. A `provisional` label means that a checker may reproduce symbolic identities while the accompanying mathematical interpretation, completeness, and literature positioning still require independent review.

## File disposition

| Original relative path | Bytes | SHA-256 | Classification | Staged path or omission reason |
|---|---:|---|---|---|
| `README.md` | 3,376 | `5f7ab0aed35bab68f13110ee4c8096c1c134dd27f48a53bbed4b9badc6be2d39` | private/internal | Omitted: session packaging notes, model/vendor attribution, posting instructions, contact suggestions, author placeholders, and unresolved checklist items. |
| `big_game_jc2_attack.md` | 7,758 | `bc5b60651b2b515c080f0457184b81092c028d153d3b3064f2e1fb7a80417031` | private/internal | Omitted: time-sensitive competitive strategy, monitoring instructions, unsupported novelty language, and source claims explicitly marked for re-verification. The conservative JC2 memo already in the repository is the public-safe replacement, not an exact duplicate. |
| `degree5_S5_etale_map.md` | 5,847 | `eda43cc78f81bcfa45e5300dad44d452c5949da6da928e3c4048d04ad373d48c` | superseded-but-historically-useful | Omitted: session/conversation references and claims superseded by the more conservative current paper, which supplies cleaner exact certificates and qualified attribution. |
| `etale_note.pdf` | 336,679 | `4decf0619428e8048b89045f9a9b6d432f19bf467cb42b09edc0e247fbb33e7a` | generated/build artifact | Omitted: compiled binary of the superseded draft. |
| `etale_note.tex` | 31,082 | `9f961ce7152e163b09c16232cbe3c8492c7f953722c1747553be2967ce1807d9` | superseded-but-historically-useful | Omitted: superseded manuscript with author/contact placeholders, unresolved bibliography TODOs, model attribution, unsupported novelty claims, and a known incorrect degree-zero explanation. The current paper corrects and narrows these points. |
| `obstruction_theorems.md` | 6,377 | `25f02bb7ad1d985d1e30d6dfd3c459bdd360a56d6fd0bcb7b6ef6b27467bc05b` | unsafe/false | Omitted: presents the cyclic obstruction and degree-two exclusion as new and proved, while the current paper replaces the degree-two argument with the classical Galois-case result and does not retain the broader novelty claim. It also describes representative symbolic checks as verification of a general planar theorem. Requires expert resolution before public inclusion. |
| `prong2_DC1_opening.md` | 4,277 | `c4eeb0029e828f02e06338d8b0d77af054a79106d9363ec9d31c1cf4f0a10376` | provisional | `provisional/dixmier-band-program/dc1-opening.md`; added a status banner and replaced one session-log phrase. |
| `prong2_milestone2_band_rigidity.md` | 7,966 | `0e1145878ae1fd19ce63070bdfcdc0d0a55fc5ae38fc33d9b2d975082642b961` | provisional | `provisional/dixmier-band-program/band1-rigidity-milestone.md`; added a status banner. Claims of equivalence, completeness, and novelty remain unreviewed. |
| `theorem_J3_band2_nonsquare.md` | 7,476 | `287c249e55fe571294945cd3169528f8b2c4c62f68fe1c25dbe5c8a6f53055e0` | duplicate | Omitted: byte-identical to the already staged provisional copy and to the two-file comparison ZIP. |
| `verify_J3.py` | 6,704 | `5ce540801ab6d5cbd500ae1b986c03cbbb90bcd5d1ad0391b3d8c1589bf2df62` | duplicate | Omitted: byte-identical to the already staged provisional copy and to the two-file comparison ZIP. |
| `verify_all.py` | 10,758 | `95407e83b64b33d4b4a6994658c098a88326af4ea8fab23c7d20712d5101d6df` | superseded-but-historically-useful | Omitted: all checks pass, but the script accompanies the superseded manuscript and includes numerical root work. The current paper has a narrower exact checker without numerical root counts. |
| `verify_band2_beachhead.py` | 7,329 | `890f1eca44f5b03b7d294b15bd6b1c1c67454a7c851846d451928531de3b1ab1` | provisional | `provisional/dixmier-band-program/verify_band2_beachhead.py`; added a status header. All 11 archived checks passed during import review. |
| `verify_prong2.py` | 6,679 | `cd732233da26b5bc1225f901220e17c3625a4ac81bb6ad894a292f208297faf7` | provisional | `provisional/dixmier-band-program/verify_band1.py`; added a status header. All 9 archived checks passed during import review. |
| `weyl_A3_endomorphism.md` | 2,841 | `af8afcbedbac7b28dcacc32e437431279f2cf2fe8e5a5abfca253c62922e50db` | provisional | `provisional/weyl-a3/endomorphism-formulas.md`; added a status banner. Preserved as explicit formulas; the non-surjectivity interpretation needs independent review. |

## Duplicate and version relationships

- `theorem_J3_band2_nonsquare.md` and `verify_J3.py` are exact duplicates of the existing provisional repository files and of the comparison ZIP entries.
- `etale_note.tex`, `etale_note.pdf`, and `verify_all.py` form an older manuscript/checker bundle. The current paper is not byte-identical: it narrows claims, corrects the degree-zero discussion, replaces the claimed-new degree-two obstruction with the classical Galois-case argument, adds an all-degrees Legendre construction, and uses cleaner exact monodromy certificates.
- `degree5_S5_etale_map.md` is an earlier working-note version of material incorporated more conservatively into the current paper.
- `big_game_jc2_attack.md` is not a version of the existing conservative JC2 memo; it is a more speculative and internally strategic branch and is intentionally omitted.

## Security and privacy scan

No credential, token, password, login, private-key, authorization-header, or obvious secret pattern was found. No hidden archive entries or private filesystem paths were found. No personal contact data was present beyond unresolved author/address/email placeholders. Internal/session references and model/vendor attribution were confined to omitted source files. The staged files contain no raw archives, caches, TeX auxiliaries, compiled binaries, login details, conversation logs, or local filesystem paths.

## Verification status

During import review, the archived scripts reported:

- `verify_prong2.py`: 9/9 checks passed.
- `verify_band2_beachhead.py`: 11/11 checks passed.
- `verify_J3.py`: 13/13 checks passed; omitted as a duplicate.
- `verify_all.py`: 31/31 checks passed; omitted with its superseded manuscript bundle.

These executions reproduce the scripts' tested identities only; they do not establish literature novelty, proof completeness, or the truth of untested generalizations.
