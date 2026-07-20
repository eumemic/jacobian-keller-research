# Trust-but-verify: the M5 band-2 square-sector partial cascades

**INDEPENDENT RE-VERIFICATION — EXACT SYMPY — NOT A CLASSIFICATION CERTIFICATE**

This memo records an independent, from-scratch symbolic re-derivation of every
load-bearing identity in the two upstream partial cascades

- `research/band2-m5-partial/classical-square-sector-partial.md`
- `research/band2-m5-partial/quantum-shifted-square-sector-partial.md`

against the frozen conventions of `research/band2-classical-proved/M4_proof_memo.md`.
All checks are reproduced by `verify_catalog.py` (Sections 1–3), run with exact
arithmetic. No arbitrary-degree exhaustion is claimed; this verifies the identities
the memos actually assert, not the resistant branches they explicitly leave open.

## Verdict (report defects loudly)

**No defects were found in either partial-cascade memo.** Every equation, integration,
reduction, family construction, and membership statement I re-derived independently
matches the upstream text exactly. Both memos are also honest about their gaps
(the classical A\* / quantum `h=1, κμs≠0` resistant branches), and those honestly-labelled
gaps are genuine — I did not close them and neither did they.

The only sign defects on record in this campaign are those the checker-audit
(`research/band2-m5-partial/checker-audit.md`) found in the **supplied M5 checker**,
which is a *separate, rejected* artifact and is **not** one of the two memos verified
here. Those defects do not touch the two partial cascades. I re-derived everything
in the single self-consistent orientation `{G,F}=G_ξ F_x − G_x F_ξ`, `{ξ,x}=1`
(classical) and `Q_m` from the crossed-product commutator `[D,X]` (quantum); all
signs close.

## Conventions frozen and independently re-checked

- Classical Poisson `{G,F}=G_ξ F_x − G_x F_ξ`, `τ=xξ`, base pair `{ξ,x}=1`.
- Band-2 decomposition `F=Σ_{k=−2}^2 x^k a_k(τ)`, obtained by the substitution `ξ→τ/x`
  and collecting powers of `x`. Membership `τ|a_{−1},b_{−1}`, `τ²|a_{−2},b_{−2}`
  is automatic here because the objects checked are genuine polynomials.
- Quantum `A_1[x^{−1}]`, `(x^a f(E))(x^b g(E))=x^{a+b} f(E+b) g(E)`,
  ladder coefficient `Q_m=Σ_{k+l=m}(b_l^{[k]}a_k − a_k^{[l]}b_l)`, `Q_m=δ_{m0}`.

## Classical (Sections 1–2 of the verifier)

| upstream claim | independent check | result |
|---|---|---|
| nine `C_m` are the `x^m`-components of `{G,F}` | direct 2-variable differentiation of `F,G` vs the coefficient sum, all `m∈[−4,4]` | **match** |
| the nine displayed `C_m` (§1) incl. all cross-terms `a_2'b_{−1}` etc. | compare each displayed formula to `C_m` | **match** |
| top: `C_4⇒b_2=λa_2`, gauge `b_2→0`, `C_3⇒b_1=κh` | substitute, confirm `C_3≡0` with `b_1=κh`, `a_2=h²`, `b_2=0` | **match** |
| (3.1) `a_1=hp`, `b_0=κp/2+β` | confirm `C_2≡0` under (3.1) | **match** |
| (3.2) `a_0=p²/4+2hv/κ−A` | confirm `C_1≡0` under (3.2) | **match** |
| (3.3)-(3.4) central integral `2h²w+hpv−κh·a_{−1}=t+e`, hence `h∣t+e` | confirm derivative `=1`, integral `=t+e`, and `C_0≡1` | **match** |
| exact residual system (4.1) `R_{−1..−4}` | confirm `R_m` equals `C_m` for `m=−1,−2,−3,−4` under the integrated forms | **match** |
| canonical family (5.1) `F=U²−Z/κ−A`, `G=λF+κU+β` | `{Z,Y}=1` and `{G,F}=1` for **generic** `h(τ)` (rational identity) | **match** |
| (5.1) coefficient table (reduced gauge) | full `a_k,b_l` table for constant `h=d`; in particular `b_{−2}=0` (branch B0) | **match** |
| nonconstant `h=α(τ+e)` fails membership | `a_{−2}=c_1²/α²` (constant, fails `τ²∣·`); if `c_1=0` then `a_{−1}=−1/(ακ)` (constant, fails `τ∣·`) | **match** |
| `κ=0` exclusion (§7) | `C_1=2h(hv)'`; with `v=0`, `C_0=2(h²w)'=1`, so `h²∣t+2e_0` kills nonconstant `h` | **match** |

The central first integral was verified two ways: as the exact antiderivative
`(2h²w+hpv−κh a_{−1})'=1`, and by confirming `C_0=1` on the fully substituted
coefficient vector. The divisibility conclusion `h∣τ+e` therefore follows from
coefficient polynomiality alone, exactly as the memo states, with no degree bound.

## Quantum (Section 3 of the verifier)

| upstream claim | independent check | result |
|---|---|---|
| `Q_m` ladder formula | build `[D,X]=DX−XD` in the abstract crossed product for generic `a_k(E),b_l(E)`; compare each `x^m`-coefficient to `Q_m` | **match** |
| ladder-4 `b_2=λa_2` (`a_2=h(E)h(E+1)`) | `Q_4=0` under `b_2=λa_2` | **match** |
| ladder-3 `b_1=κh` (gauge `b_2=0`) | `Q_3=0` under `b_1=κh` | **match** |
| ladder-2 (3.1) `a_2(B^{[2]}−B)+κ(h^{[1]}a_1−h a_1^{[1]})=0` ⇒ `h∣a_1` | `Q_2` equals the displayed form; `r=a_1/h` polynomial ⇒ `a_1=hp` | **match** |
| quantum midpoint (3.2) `B^{[1]}+B=κp+γ` (NOT the naive `κp/2+β`) | applying `(T−1)` recovers the ladder-2 difference equation `B^{[2]}−B=κ(p^{[1]}−p)` | **match** |
| exact reduced tail (4.1)-(4.6) | each displayed equation matches `Q_1,Q_{−4},Q_{−3},Q_{−2},Q_{−1},Q_0` (with the staggered shifts) | **match** |
| negative proportionality (4.7) `w=μs` | `Q_{−4}=0` under `w=μs` | **match** |
| §5 quantum boundary family, `[D,X]=1`, `u` constant | constructed explicitly (linear `h`, `p=0`, `v=s=w=0`); solved `A=−1/(κc)`; confirmed `[D,X]=1` and `u` constant nonzero (fails `E∣u`, polar) | **match** |

The quantum midpoint check is the sharp one: the memo is correct that
`b_0=κp/2+β` is *generally false* for nonconstant `p`; the right statement is the
finite-difference `(T+1)B=κp+γ`, and my check confirms this is the exact ladder-2
consequence. The staggered-shift tail equations (4.1)-(4.6) are reproduced verbatim
from the commutator, including the `κh^{[−2]}s` / `κs^{[1]}h` type terms that are
easy to mis-shift by hand.

## What is *not* verified (unchanged open status)

- The classical resistant branch **A\*** (`h` const, `κ≠0`, gauged `b_{−2}≠0`,
  `a_{−2}≠0`) and the quantum counterpart `h=1, κμs≠0` remain open. My verification
  confirms the *identities and reductions*; it does not solve or exclude those
  branches at arbitrary degree. (But see `tame-catalog.md`: no *tame* pair reaches
  A\*, which reframes what closing A\* means.)
- No JC2 / DC1 claim. No completeness / exhaustion claim.

## Reproduce

```
uv run --with sympy python research/band2-square-sector/verify_catalog.py
```

Sections 1–3 are the classical+quantum re-verification above; the script prints a
`PASS` line per row and aborts with a `DEFECTS FOUND` block if any identity fails.
It currently passes all rows and ends `ALL CATALOG CHECKS PASSED`.
