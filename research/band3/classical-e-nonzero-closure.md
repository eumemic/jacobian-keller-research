# The classical band-3 e ≠ 0 mixed sector: the second first integral, and closure of A\*

**INDEPENDENTLY DERIVED — MACHINE-CHECKED IDENTITIES — NOT PEER REVIEWED — BAND-SCOPED**

This memo attacks the one branch of the classical cubic sector `a_3 = c h³`
(gauge `b_3 = 0`) that Wave B named **counterexample-or-nothing** and left open:
the **`e ≠ 0` mixed sector** — gauged wall constant `b_2 = e h² ≠ 0`, constant `h`.
It **produces the missing second first integral** that
`research/band3/classical-hard-branches.md` (commit `ebfc64d`, cited [HARD]) pinned
as *"the precise missing step"* — the band-3 analogue of band-2's `I₂`
(`research/band2-square-sector/classical-Astar.md`, commit `84978b9`, cited [B2A]).
It then uses it to **close the resistant branch A\*-band3 for `deg a_2 ≥ 1`** by a
tropical degree gap, and delimits sharply what remains.

Conventions are frozen and identical to Wave A/B: over `C`, `τ = xξ`,
`{G,F} = G_ξ F_x − G_x F_ξ`, `F = Σ_{k=-3}^{3} x^k a_k(τ)`, `G = Σ x^k b_k(τ)`,
membership `τ^j | a_{-j}, b_{-j}`, ladder `C_m = Σ_{k+l=m}(k a_k b_l' − l a_k' b_l) = δ_{m0}`,
primes `d/dτ`. Diagonal scaling normalizes `c = 1` (genuine, [HARD] §3;
`classical-cube-closure.md` §0, cited [CUBE]) and `h = 1`. Every displayed algebraic
identity is machine-checked by `verify_classical_e_closure.py` (same directory,
31 checks); a successful run ends `ALL CLASSICAL E CLOSURE CHECKS PASSED`. Degree /
membership statements are **written arguments**; the bounded box and Gröbner sweeps
are **regression corroboration only**, labelled as such. Framework verifiers (cascade,
weapons) are green at commit `b9f9cf3`.

---

## 0. Headline

> **The band-3 analogue of band-2's `I₂` exists in the `e ≠ 0` sector.** With the
> two nonlocal-augmented first integrals
> ```
>   Φ' = C_{-1}            (nonlocal Q₁ = ∫a_2),
>   I₂' = C_{-2} − (2/3) a_2' Φ   (nonlocals Q₁, Q₂ = ∫a_2², P₁ = ∫a_1),
> ```
> the multiplier **`(2/3) a_2`** is *forced by, and exactly absorbs*, the W3
> obstruction residue `(4e/3) a_{-3} a_2'` that [HARD] §5 proved blocks the naive
> trailing integral: `Euler_{a_{-3}}(C_{-2}) = −(4e/3) a_2' ≠ 0`, but
> `Euler_{a_{-3}}(C_{-2} − (2/3)a_2' Φ) = 0` (and the `a_{-2}` residue vanishes too).
> `Φ|_{e=0}` is **exactly** [CUBE]'s `Φ₁`, so `Φ` is the honest `e ≠ 0` generalization.
> *(PROVED, machine-checked, §2–§3.)*
>
> **Consequence — A\*-band3 (`e ≠ 0`) is empty for `deg a_2 ≥ 1`.** `Φ` and `I₂` are
> **linear** in the trailing pair `(a_{-2}, a_{-3})`; the `2×2` solve (invertible off
> the locus `det = 0`) **determines `a_{-3}`** — the determination the W3 obstruction
> denied `C_{-2}` alone. Its leading coefficient is a **nonzero** multiple of the data
> (`−13 lc(a_2)⁶/2187` when `a_2` dominates), giving `deg a_{-3} = 6·deg a_2`, whereas
> the moment-determined `b_{-3}` has `deg b_{-3} = 5·deg a_2` (leading coefficient
> `−25 e lc(a_2)⁵/243 ≠ 0`). So **`deg b_{-3} < deg a_{-3}`**: the bottom Wronskian
> `C_{-6}` (which forces `b_{-3} = μ_3 a_{-3}`) is unsatisfiable for `μ_3 ≠ 0`. The
> resistant branch **A\*-band3 is empty at arbitrary degree whenever `deg a_2 ≥ 1`**.
> *(PROVED modulo the termwise-domination degree argument, §5; box-corroborated.)*
>
> **Modulus verdict: `other` (tropical), and steeper than `e = 0`.** As in [CUBE]
> the band-3 obstruction is a **max-degree gap**, not a linear congruence — but the
> `e ≠ 0` gap is **`6 : 5`** (from the `e`-wall's `a_2`-power growth) versus the
> `e = 0` gap `3 : 2` ([CUBE] §7). No residue-class modulus exists.

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
`(4e/3) a_{-3} a_2'` — non-exact for `e ≠ 0` — so no gauge-free trailing integral of
`C_{-2}` exists, and *"with a single balance there is no congruence to close on;
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

**Why the multiplier is `(2/3) a_2`, and why it works.** For *any* multiplier `m`,
```
   Euler_{a_{-3}}(C_{-2} − m' Φ) = −(4e/3) a_2'  +  2e·m'
```
(verifier §3; the `2e` is `Φ`'s `a_{-3}`-coefficient `−2e` with sign). This vanishes
**iff `m' = (2/3) a_2'`, i.e. `m = (2/3) a_2`** — the multiplier is *forced*, and it
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
Unlike [B2A], where memberships forced `Φ_0 = 0`, here `Φ_0 = Φ(0)` is a genuine
nonzero constant in `a_2(0), a_1(0), a_0(0)`; the degree argument uses only that
`Φ, J` are **constant**, i.e. their positive-degree parts vanish.

## 4. `Φ, I₂` determine the trailing pair — the W3 lever restored

`Φ = Φ_0` and `I₂ = J − (2/3)Φ_0 a_2` are **linear in `(a_{-2}, a_{-3})`**; the `2×2`
coefficient matrix has determinant (verifier §5)
```
   det = −(4/3) e² a_1 + (4/9) e² a_2² − κ₁².
```
Off the locus `det ≡ 0`, the two integrals **determine `a_{-2}` and `a_{-3}`** from
the upper data `(a_2, a_1, a_0, a_{-1})` and the scalars. This is exactly the
determination the W3 obstruction denied `C_{-2}` alone: **`I₂` is what makes the
trailing coefficient `a_{-3}` solvable in the `e ≠ 0` sector.** (For `e = 0` the
determinant degenerates to `−κ₁²` and this reduces to [CUBE] §2's `κ ≠ 0` lever.)

## 5. The tropical gap: A\*-band3 (`e ≠ 0`, `deg a_2 ≥ 1`) is empty

Write `Q = deg a_2, P = deg a_1, R = deg a_0, L = deg a_{-1}`. Solving the `2×2` for
the leading behaviour (verifier §5):

> **Leading-coefficient certificates (`a_2` dominant, PROVED — machine-checked).**
> `deg a_{-3} = 6Q` with leading coefficient `−13 lc(a_2)⁶/2187 ≠ 0`; and the
> moment-determined `deg b_{-3} = 5Q` with leading coefficient `−25 e lc(a_2)⁵/243 ≠ 0`.

So `a_{-3}` grows **one `a_2`-power faster than `b_{-3}`** — the cube's `6 : 5`
signature (from `Φ`'s `a_2⁶` and `I₂`'s `a_2⁷` versus the moment's `a_2⁵`). Across all
regimes the general behaviour is `deg a_{-3} = max(6Q, 2R, …)` with a nonzero leading
coefficient, and **`deg b_{-3} < deg a_{-3}` for every `deg a_2 ≥ 1`** (verifier §5:
zero failures over the box `Q ≤ 3, P,R,L ≤ 3`; the arbitrary-degree statement is the
termwise-domination written argument, exactly as [CUBE] §7).

> **Theorem 5.1 (A\*-band3, `e ≠ 0`, `deg a_2 ≥ 1`: EMPTY).** The bottom Wronskian
> `C_{-6} = 0` with `a_{-3} ≠ 0` forces `b_{-3} = μ_3 a_{-3}` for a **constant** `μ_3`
> ([HARD] §1). With `μ_3 ≠ 0` this requires `deg b_{-3} = deg a_{-3}` and proportional
> polynomials — impossible since `deg b_{-3} < deg a_{-3}`. Hence **no A\*-band3 pair
> has `deg a_2 ≥ 1`**.

*Proof.* `a_{-3}` is determined (§4, `det ≠ 0`) with nonzero leading coefficient, and
`deg b_{-3} < deg a_{-3}` for every `deg a_2 ≥ 1` (verifier §5 box; termwise argument).
Two polynomials of different degree are not proportional, so `b_{-3} = μ_3 a_{-3}`
fails for `μ_3 ≠ 0`. ∎

> **Corollary 5.2 (`deg a_2 ≥ 1`, `a_{-3} ≠ 0`: EMPTY off `lc(b_{-3}) = 0`).** The gap
> also forces the `μ_3 = 0` (B0) alternative `b_{-3} ≡ 0`; but the moment-determined
> `b_{-3}` has **nonzero** leading coefficient `−25 e lc(a_2)⁵/243` (a nonzero multiple
> of `e`, since `e ≠ 0`), so `b_{-3} ≢ 0`. Contradiction. Hence the **entire** `e ≠ 0`
> sector with `deg a_2 ≥ 1` and `a_{-3} ≠ 0` is empty **off the finite union of
> leading-coefficient loci where `lc(b_{-3})` vanishes** (where `b_{-3}` drops degree
> and the argument descends).

Theorem 5.1 is unconditional on those loci (the gap only widens when `lc(b_{-3})`
vanishes, since `a_{-3}`'s leading coefficient is robustly nonzero); Corollary 5.2's
descent is the residual coefficient problem, exactly the shape [CUBE] §7 met at `e = 0`.

**Consistency.** The [HARD] §5 Gröbner box `deg(a_2,a_1,a_0,a_{-1},a_{-2},a_{-3}) =
(1,1,1,2,2,3)` has `deg a_2 = 1`, so it is covered by Theorem 5.1; [HARD] found it
`Gröbner = (1)` (empty), and the verifier re-derives a small `e ≠ 0` box
`(1,0,0,1,2,3)` as `(1)` independently (§6). The theorem **explains** the bounded
emptiness as an arbitrary-degree phenomenon for `deg a_2 ≥ 1`.

## 6. What remains (sharp delimitation)

The tropical gap closes A\*-band3 for `deg a_2 ≥ 1`. Precisely delimited residue:

- **`a_2 = const` (`Q = 0`) stratum.** With `a_2` constant `b_1 = κ₁ + (2e/3)a_2` is
  **constant**, so the W3 residue `−(4e/3)a_2'` **switches off** (verifier §6) and the
  sector is unobstructed there — the `e ≠ 0` analogue of [CUBE]'s `e = 0` slice, with
  shifted constant `κ₁' = κ₁ + (2e/3)a_2`. A **finite** set of degree loci survives the
  gap (the `Q = 0` survivors — related to [CUBE]'s `{3P=2R, L=2P, 2R=L+P}` but enlarged
  by the `b_2 = e` coupling; verifier §6 exhibits survivors off those three, e.g.
  `(P,R,L) = (0,1,1)`). Bounded-empty ([HARD] §5 Gröbner). **Not closed at arbitrary
  degree here.**
- **Leading-coefficient cancellation loci** (`deg a_2 ≥ 1`): the finite union where
  `lc(b_{-3}) = 0` (Corollary 5.2) — a coefficient problem, not a congruence.
- **`det ≡ 0` locus** `a_1 = (a_2² )/3 − 3κ₁²/(4e²)` (a_1 tied to a_2): the `2×2` lever
  degenerates; separate handling required.
- **onesided-top (`a_{-3} = 0`)**: `C_{-6}` is vacuous; routes toward band `≤ 2`
  (classical-band3-cascade §9).
- **Nonconstant-`h`, `e ≠ 0` cross-branch**: out of scope by mandate ([HARD] §6 open
  item); the divisibility engine bites there (`h | a_2`, `rad(h) | a_2/h`) but the
  `τ`-order kill is not carried through.

## 7. Claim disposition

**Proved (machine-checked identities, arbitrary degree):**
- the `e ≠ 0` reduction: `b_0, b_{-1}` solve `C_3, C_2`; `b_{-2}` explicit (`C_1` RHS a
  total derivative), now `a_{-1}`-dependent; `b_{-3}` from the moment; `C_0 = 1` (§1).
- **`Φ' = C_{-1}`** (nonlocal `Q₁ = ∫a_2`); `Φ|_{e=0} = ` [CUBE]'s `Φ₁` exactly (§2).
- **`I₂' = C_{-2} − (2/3) a_2' Φ`** — the missing second first integral — with the
  multiplier `(2/3)a_2` forced by and absorbing the W3 residue
  (`Euler_{a_{-3}}(C_{-2}) = −(4e/3)a_2'`, `Euler_{a_{-3}}(C_{-2} − (2/3)a_2'Φ) = 0`,
  `Euler_{a_{-2}}(⋯) = 0`); nonlocals `Q₁, Q₂ = ∫a_2², P₁ = ∫a_1` (§3).
- the two conserved quantities `Φ = const`, `J = I₂ + (2/3)Φ_0 a_2 = const` (§4); the
  `2×2` determinant and the trailing-pair determination off `det = 0` (§4–§5).

**Proved modulo the termwise-domination degree argument (leading coefficients
machine-checked, gap box-verified):**
- `deg a_{-3} = 6·deg a_2` (lc `−13 lc(a_2)⁶/2187`), `deg b_{-3} = 5·deg a_2`
  (lc `−25 e lc(a_2)⁵/243`), hence `deg b_{-3} < deg a_{-3}` for `deg a_2 ≥ 1` (§5).
- **Theorem 5.1: A\*-band3 (`e ≠ 0`, `deg a_2 ≥ 1`) is EMPTY**, and **Corollary 5.2**:
  the full sector with `deg a_2 ≥ 1`, `a_{-3} ≠ 0` is empty off the finite
  `lc(b_{-3}) = 0` loci (§5).
- **Modulus verdict = `other` (tropical `6 : 5`)**, steeper than the `e = 0` gap
  `3 : 2` — no residue-class modulus in the `e ≠ 0` sector (§0, §5).

**Computed (bounded, corroboration only):**
- `e ≠ 0` Gröbner box `(1,0,0,1,2,3)` is `(1)` (§6); consistency with [HARD]'s
  `(1,1,1,2,2,3)` box (both are `deg a_2 = 1`, covered by Theorem 5.1).

**Open (precisely delimited):**
- the `a_2 = const` (`Q = 0`) stratum's finite surviving degree loci (bounded-empty);
- the `deg a_2 ≥ 1` leading-coefficient loci `lc(b_{-3}) = 0` and the `det ≡ 0` locus;
- onesided-top (`a_{-3} = 0`, routes to band `≤ 2`);
- the nonconstant-`h`, `e ≠ 0` cross-branch ([HARD] §6, out of scope).

**Not claimed:** a full band-3 theorem, the whole `e ≠ 0` sector at `deg a_2 = 0`,
the `A*`-band3 quantum mirror, JC2, DC1, or any statement beyond §1–§6.

## 8. Verification

```
uv run --with sympy python research/band3/verify_classical_e_closure.py
```
Exact SymPy, 31 checks; ends `ALL CLASSICAL E CLOSURE CHECKS PASSED`. It certifies:
the 13 `C_m` vs the direct Poisson bracket; the `e ≠ 0` reduction and moment; `Φ' =
C_{-1}` and `Φ|_{e=0} = Φ₁`; **`I₂' = C_{-2} − (2/3)a_2' Φ`** with the forced-multiplier
Euler-residue identities; the two conserved quantities; the `2×2` determinant and the
leading-coefficient certificates `deg a_{-3} = 6·deg a_2`, `deg b_{-3} = 5·deg a_2`
with their exact leading coefficients; the gap over the box; the `a_2`-constant
residue switch-off; and a small `e ≠ 0` Gröbner emptiness box. The degree/membership
arguments of §5–§6 are the proofs; the box and Gröbner sweeps are corroboration.
