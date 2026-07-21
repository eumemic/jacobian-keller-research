# The classical band-3 e ≠ 0 mixed sector: the second first integral and conditional tropical certificates

**INDEPENDENTLY DERIVED — MACHINE-CHECKED IDENTITIES — NOT PEER REVIEWED — BAND-SCOPED**

This memo attacks the one branch of the classical cubic sector `a_3 = c h³`
(gauge `b_3 = 0`) that Wave B named **counterexample-or-nothing** and left open:
the **`e ≠ 0` mixed sector** — gauged wall constant `b_2 = e h² ≠ 0`, constant `h`.
It **produces the missing second first integral** that
`research/band3/classical-hard-branches.md` (commit `ebfc64d`, cited [HARD]) pinned
as *"the precise missing step"* — the band-3 analogue of band-2's `I₂`
(`research/band2-square-sector/classical-Astar.md`, commit `84978b9`, cited [B2A]).
It then derives a **conditional trailing solve off an explicit determinant locus** and
records an `a_2`-dominant tropical degree gap. A bounded generic-leading-monomial
scan corroborates that gap in a small range but does not close the arbitrary-degree
sector. The determinant, cancellation/tie, polynomiality, membership, and constant-
`a_2` loci remain explicit.

Conventions are frozen and identical to Wave A/B: over `C`, `τ = xξ`,
`{G,F} = G_ξ F_x − G_x F_ξ`, `F = Σ_{k=-3}^{3} x^k a_k(τ)`, `G = Σ x^k b_k(τ)`,
membership `τ^j | a_{-j}, b_{-j}`, ladder `C_m = Σ_{k+l=m}(k a_k b_l' − l a_k' b_l) = δ_{m0}`,
primes `d/dτ`. Diagonal scaling normalizes `c = 1` (genuine, [HARD] §3;
`classical-cube-closure.md` §0, cited [CUBE]) and `h = 1`. Every displayed algebraic
identity is machine-checked by `verify_classical_e_closure.py` (same directory,
31 checks); a successful run ends `ALL CLASSICAL E CLOSURE CHECKS PASSED`. The
stated leading data are checked only at the representative
`(Q,P,R,L)=(2,1,1,1)`, with conserved constants zero and no tie or cancellation;
the wider bounded leading-monomial and Gröbner sweeps are **corroboration only** and do not prove
polynomiality, membership, cancellation/tie cases, or arbitrary-degree emptiness. Framework verifiers (cascade,
weapons) are green at commit `b9f9cf3`.

---

## 0. Headline

> **The band-3 analogue of band-2's `I₂` exists in the `e ≠ 0` sector.** With the
> two nonlocal-augmented first integrals
> ```
>   Φ' = C_{-1}            (nonlocal Q₁ = ∫a_2),
>   I₂' = C_{-2} − (2/3) a_2' Φ   (nonlocals Q₁, Q₂ = ∫a_2², P₁ = ∫a_1),
> ```
> the multiplier **`m=(2/3) a_2 + constant`** is forced up to its additive constant,
> and its derivative exactly absorbs the W3 obstruction residue
> `(4e/3) a_{-3} a_2'` that [HARD] §5 proved blocks the naive
> trailing integral: `Euler_{a_{-3}}(C_{-2}) = −(4e/3) a_2' ≠ 0`, but
> `Euler_{a_{-3}}(C_{-2} − (2/3)a_2' Φ) = 0` (and the `a_{-2}` residue vanishes too).
> `Φ|_{e=0}` is **exactly** [CUBE]'s `Φ₁`, so `Φ` is the honest `e ≠ 0` generalization.
> *(PROVED, machine-checked, §2–§3.)*
>
> **Conditional consequence and bounded corroboration.** `Φ` and `I₂` are **linear**
> in the trailing pair `(a_{-2}, a_{-3})`. They determine that pair only where
> `det = −(4/3)e²a_1+(4/9)e²a_2²−κ₁²` is nonzero, and the resulting rational
> expressions still require denominator cancellation, polynomiality, and relevant
> membership checks. At the representative `(Q,P,R,L)=(2,1,1,1)`, with conserved
> constants zero and under the stated no-tie/no-cancellation assumptions, the checked
> leading data are `deg a_{-3}=12`, `deg b_{-3}=10`, with respective leading
> coefficients `−13 lc(a_2)⁶/2187` and `−25e lc(a_2)⁵/243`. A generic-leading-
> monomial scan over `Q=1..3`, `P,R=0..3`, `L=1..3`, with conserved constants set
> to zero and cancellation ties omitted, finds the same strict gap throughout that
> bounded box. This is evidence, not an arbitrary-degree emptiness theorem.
>
> **Shape verdict: tropical evidence, not a derived modulus.** The `a_2`-dominant
> certificate has a `6 : 5` degree signature, versus the `e=0` `3 : 2` signature
> ([CUBE] §7). No Band-2-shaped congruence lattice is derived here; wider regimes,
> cancellation/tie loci, and the determinant locus remain open.

---

## 1. The `e ≠ 0` reduction (h = 1, c = 1)

By [HARD] §5 the positive cascade determines, with free F-data
`q = a_2, p = a_1, r = a_0, al = a_{-1}` and trailing `s = a_{-2}, sig = a_{-3}`
(verifier §1):
```
b_2 = e,   b_1 = κ₁ + (2e/3) a_2,
b_0 = (2e a_1 + κ₁ a_2)/3 − e a_2²/9 + β,
b_{-1} = (1/3)[ 2e a_0 + κ₁ a_1 − (2e/3) a_1 a_2 − (κ₁/3) a_2² + (4e/27) a_2³ ] + γ.
```
`C_1 = 0` integrates **exactly** (its RHS is a total `τ`-derivative) to an explicit
`b_{-2}`, which — unlike the `e = 0` case — **depends on `a_{-1}`** (the `−2e a_{-1}`
term); this `e`-coupling is the structural novelty. The moment `C_0 = M' = 1`
(`M = 3(a_3 b_{-3}) + 2(a_2 b_{-2} − a_{-2} b_2) + (a_1 b_{-1} − a_{-1} b_1)`) gives
`b_{-3}`, so **all of `b_0..b_{-3}` are explicit** and the residual is
`C_{-1} = ⋯ = C_{-6} = 0` plus memberships. *(Verifier §1, all identities.)*

## 2. Φ : the first integral of `C_{-1}` (the honest `e ≠ 0` W3-replacement)

The Euler operators of `C_{-1}` vanish except `Euler_{a_2}(C_{-1}) = −1/3` (a
*constant* density — exactly [CUBE]'s obstruction), so a single nonlocal generator
`Q₁ = ∫a_2` restores integrability:

> **Lemma 2.1.** `Φ' = C_{-1}` identically (verifier §2), where `Φ = Φ_{e=0} + e·Φ_e`,
> `Φ_{e=0} = 2δ a_1 + γ a_0 − κ₁ a_{-2} − (5δ/3)a_2² − (4γ/3)a_1 a_2 + (2κ₁/3)a_{-1}a_2
>  + (2κ₁/3)a_1 a_0 + τ a_2 − Q₁/3 + (14γ/27)a_2³ − (5κ₁/9)a_1²a_2 − (5κ₁/9)a_2²a_0
>  + (40κ₁/81)a_1 a_2³ − (22κ₁/243)a_2⁵`
> is **exactly [CUBE]'s `Φ₁`** (verifier §2), and
> `Φ_e = (2/3)a_{-1}a_1 − (4/9)a_{-1}a_2² − (4/27)a_1³ + (14/27)a_1²a_2² − (70/243)a_1 a_2⁴
>  − (8/9)a_1 a_2 a_0 + (91/2187)a_2⁶ + (28/81)a_2³a_0 + (2/3)a_2 a_{-2} + (1/3)a_0² − 2 a_{-3}`.

`Φ` is **linear** in the trailing coefficients: coefficient of `a_{-3}` is `−2e`,
coefficient of `a_{-2}` is `(2e/3)a_2 − κ₁` (verifier §2). This is the exact `e ≠ 0`
generalization of [CUBE] Lemma 3.1 / [B2A] Lemma 2.1's `Φ`.

## 3. `I₂` : the second first integral (the missing step)

[HARD] §5 pinned the obstruction: `a_{-3}` enters `C_{-2}` with the residue
`(4e/3) a_{-3} a_2'` — non-exact for `e≠0` — obstructing the raw, gauge-free
Band-2-shaped total derivative of `C_{-2}`; *"with a single balance there is no
congruence to close on;
producing the second integral is the precise missing step."* The resolution mirrors
[B2A]'s `I₂' = κ C_{-2} − p'Φ`:

> **Lemma 3.1 (the second first integral).** With the multiplier **`(2/3) a_2`**,
> ```
>   I₂' = C_{-2} − (2/3) a_2' Φ          (verifier §3),
> ```
> where `I₂` is the explicit polynomial in `{a_2,a_1,a_0,a_{-1},a_{-2},a_{-3}}` and
> the nonlocal generators `Q₁ = ∫a_2`, `Q₂ = ∫a_2²`, `P₁ = ∫a_1` displayed in the
> verifier. `I₂` is linear in the trailing pair: coefficient of `a_{-3}` is
> `−b_1 = −(κ₁ + (2e/3)a_2)`, coefficient of `a_{-2}` is `(2e/3)a_1`.

**Why the multiplier is `(2/3) a_2` up to a constant, and why it works.** For *any* multiplier `m`,
```
   Euler_{a_{-3}}(C_{-2} − m' Φ) = −(4e/3) a_2'  +  2e·m'
```
(verifier §3; the `2e` is `Φ`'s `a_{-3}`-coefficient `−2e` with sign). This vanishes
**iff `m' = (2/3) a_2'`, i.e. `m = (2/3) a_2 + constant`** — the multiplier is
*forced only up to an additive constant*, and its derivative
absorbs the W3 residue **exactly**. The same combination kills the `a_{-2}` residue:
`Euler_{a_{-2}}(C_{-2} − (2/3)a_2' Φ) = 0` (verifier §3). Both deep trailing
coefficients then sit inside exact derivatives, so `I₂` exists (with the three
nonlocal generators — one level up from [B2A]'s `∫a_1, ∫a_1², ∫b_{-1}`, matching the
cube's extra free level, exactly as [CUBE]'s `Q₁ = ∫a_2` was one level up from
[B2A]'s `∫a_1`).

**The two conserved quantities.** On a genuine pair (`C_{-1} = C_{-2} = 0`):
`Φ = Φ_0` (const), and then `I₂' = −(2/3)Φ_0 a_2'`, so
```
   J := I₂ + (2/3) Φ_0 a_2 = const                 (verifier §4).
```
Unlike [B2A], where memberships forced `Φ_0=0`, the identities here retain
`Φ_0=Φ(0)` as an unrestricted conserved scalar (it may vanish). Any later degree
analysis may use only that `Φ` and `J` are constant.

## 4. `Φ, I₂` conditionally determine the trailing pair off `det=0`

`Φ = Φ_0` and `I₂ = J − (2/3)Φ_0 a_2` are **linear in `(a_{-2}, a_{-3})`**; the `2×2`
coefficient matrix has determinant (verifier §5)
```
   det = −(4/3) e² a_1 + (4/9) e² a_2² − κ₁².
```
Where this polynomial determinant is **nonzero**, Cramer's rule determines
`a_{-2}` and `a_{-3}` as rational expressions in the upper data
`(a_2,a_1,a_0,a_{-1})` and the scalars. This is the determination the W3 obstruction
denied `C_{-2}` alone, but it is conditional: the `det=0` locus requires separate
analysis, and away from it one must still prove denominator cancellation,
polynomiality, and the relevant membership analogues before treating the rational
solve as a genuine trailing polynomial pair. (For `e=0` the determinant becomes
`−κ₁²` and reduces to [CUBE] §2's `κ≠0` lever.)

## 5. Conditional tropical certificates and bounded exploration

Write `Q=deg a_2`, `P=deg a_1`, `R=deg a_0`, `L=deg a_{-1}`. The following two
claims must be kept separate.

> **Exact leading-coefficient certificates at a representative `a_2`-dominant
> specialization (machine-checked).** At `(Q,P,R,L)=(2,1,1,1)`, provided the
> determinant has the expected nonzero leading term and no competing term ties or
> cancels it, the Cramer solve gives `deg a_{-3}=6Q` with leading coefficient
> `−13 lc(a_2)⁶/2187≠0`; the moment then gives `deg b_{-3}=5Q` with leading
> coefficient `−25e lc(a_2)⁵/243≠0`.

Thus, in this regime and subject to the conditional solve and polynomiality, the
bottom Wronskian proportionality `b_{-3}=μ_3a_{-3}` is incompatible with the strict
`6Q:5Q` degree gap. This is a regime-specific obstruction, not a theorem covering
all degree orderings.

> **Bounded generic-leading-monomial exploration (corroboration only).** The verifier
> substitutes one generic leading monomial for each upper coefficient and scans
> exactly `Q=1..3`, `P,R=0..3`, `L=1..3`. It sets the conserved constants `CP,CI`
> to zero. Within those substitutions it reports no case with
> `deg b_{-3}≥deg a_{-3}` and no case with `deg a_{-3}<6Q`.

This scan is not an arbitrary-degree proof. It does not enumerate cancellation/tie
loci, does not prove that division by `det` yields polynomials, and does not impose
all coefficient memberships after the rational solve. In particular, it cannot
establish emptiness of the full `e≠0`, `deg a_2≥1` sector.

**Bounded consistency only.** [HARD] found the degree box
`(deg a_2,a_1,a_0,a_{-1},a_{-2},a_{-3})=(1,1,1,2,2,3)` empty by Gröbner
computation. The verifier independently finds the smaller box `(1,0,0,1,2,3)`
empty. These are exact bounded computations, not consequences of an unbounded
emptiness theorem.

## 6. Open loci and limitations

- **`det=0` locus.** For `e≠0`,
  `a_1=a_2²/3−3κ₁²/(4e²)` makes the trailing matrix singular; the two-integral
  Cramer lever does not apply.
- **Denominator cancellation and polynomiality.** When `det` is nonconstant and
  nonzero, divisibility of the Cramer numerators by `det` is not proved. Nor are
  the classical membership conditions (and their Weyl-side analogues, if used)
  proved for the resulting expressions.
- **Cancellation and degree-tie loci.** Competing leading monomials can tie or
  cancel in `det`, the Cramer numerators, or `b_{-3}`. The bounded generic scan
  intentionally does not resolve these coefficient loci.
- **`a_2=const` (`Q=0`).** Then `b_1=κ₁+(2e/3)a_2` is constant and the W3 residue
  `−(4e/3)a_2'` vanishes. The verifier only exhibits the bounded survivor
  `(P,R,L)=(0,1,1)`, which lies off the three `e=0` tie equations; it does not
  classify the constant-`a_2` stratum or prove a finite global survivor list.
- **Leading-coefficient vanishing.** Any locus on which the displayed
  `a_2`-dominant coefficients cease to control the top degree requires a separate
  descent and remains open.
- **onesided-top (`a_{-3}=0`).** `C_{-6}` is vacuous and routes toward band `≤2`
  (classical-band3-cascade §9).
- **Nonconstant-`h`, `e≠0` cross-branch.** Out of scope ([HARD] §6); the available
  divisibility observations do not close it.

## 7. Claim disposition

**Proved (machine-checked identities, arbitrary degree):**
- the `e ≠ 0` reduction: `b_0, b_{-1}` solve `C_3, C_2`; `b_{-2}` explicit (`C_1` RHS a
  total derivative), now `a_{-1}`-dependent; `b_{-3}` from the moment; `C_0 = 1` (§1).
- **`Φ' = C_{-1}`** (nonlocal `Q₁ = ∫a_2`); `Φ|_{e=0} = ` [CUBE]'s `Φ₁` exactly (§2).
- **`I₂' = C_{-2} − (2/3) a_2' Φ`** — the missing second first integral — with
  `m=(2/3)a_2+constant` uniquely forced up to its additive constant and `m'`
  absorbing the W3 residue
  (`Euler_{a_{-3}}(C_{-2}) = −(4e/3)a_2'`, `Euler_{a_{-3}}(C_{-2} − (2/3)a_2'Φ) = 0`,
  `Euler_{a_{-2}}(⋯) = 0`); nonlocals `Q₁, Q₂ = ∫a_2², P₁ = ∫a_1` (§3).
- the two conserved quantities `Φ = const`, `J = I₂ + (2/3)Φ_0 a_2 = const` (§4); the
  `2×2` determinant and the trailing-pair determination off `det = 0` (§4–§5).

**Proved conditionally / in a specified regime:**
- off the zero locus of `det`, the two linear identities give a rational Cramer
  determination of `(a_{-2},a_{-3})`; polynomiality and membership are additional
  unresolved conditions (§4);
- at the representative `a_2`-dominant specialization `(Q,P,R,L)=(2,1,1,1)`,
  under the no-tie/no-cancellation assumptions, the exact leading data are
  `deg a_{-3}=6Q` (lc `−13lc(a_2)⁶/2187`) and `deg b_{-3}=5Q`
  (lc `−25e lc(a_2)⁵/243`), yielding a specialization-specific proportionality
  obstruction.

**Computed (bounded, corroboration only):**
- the generic-leading-monomial scan covers exactly `Q=1..3`, `P,R=0..3`,
  `L=1..3`, sets `CP=CI=0`, and omits cancellation/tie analysis (§5);
- the `e≠0` Gröbner box `(1,0,0,1,2,3)` is `(1)`; [HARD]'s separate
  `(1,1,1,2,2,3)` box is also bounded evidence.

**Open (precisely delimited):**
- `det=0`; denominator cancellation and polynomiality; classical membership and
  any relevant membership analogues used in a later transfer; cancellation and degree-tie loci;
- `a_2=const` (`Q=0`) and all wider degree regimes not covered by the
  `a_2`-dominant certificate;
- onesided-top (`a_{-3}=0`) and nonconstant-`h`, `e≠0`.

**Not claimed:** emptiness of the full `e≠0`, `deg a_2≥1` sector; a full band-3
or A\*-band3 theorem; a Band-2-shaped congruence lattice; the quantum mirror; JC2;
DC1; or any unbounded completeness from the bounded scan.

## 8. Verification

```
uv run --with sympy python research/band3/verify_classical_e_closure.py
```
Exact SymPy; a successful run ends `ALL CLASSICAL E CLOSURE CHECKS PASSED`. Its
exact symbolic checks certify the 13 `C_m` formula, reduction and moment, `Φ'=C_{-1}`,
`Φ|_{e=0}=Φ₁`, **`I₂'=C_{-2}−(2/3)a_2'Φ`**, the multiplier derivative condition
(`m=(2/3)a_2+constant`), the conserved quantities, determinant, and the leading
coefficients only at representative `(Q,P,R,L)=(2,1,1,1)`, with `CP=CI=0` and
ties/cancellation omitted. Separately, it performs a bounded generic-
leading-monomial exploration over `Q=1..3`, `P,R=0..3`, `L=1..3`, with conserved
constants zero and cancellation ties omitted, plus one small Gröbner box. Those
computations are corroboration only and do not certify arbitrary-degree emptiness,
denominator cancellation, polynomiality, or membership.
