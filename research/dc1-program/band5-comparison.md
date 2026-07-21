# Band-5 comparison: what drives the exotic slack — primality of the cyclotomic, or bare `k ≥ 3`?

**INDEPENDENTLY DERIVED — MACHINE-CHECKED IDENTITIES (k = 5, with k = 3, 4 controls) — NOT PEER REVIEWED — BAND-SCOPED**

This memo is the **band-5 comparison arm** of the DC1 campaign. `k = 5` is the
second prime band index (after `k = 3`): the necklace operator
`S_5 = (S^5−1)/(S−1) = Φ_5` is irreducible, exactly as `S_3 = Φ_3` is. Its
purpose is to isolate the driver of the quantum **exotic slack** — the
`k ≥ 3` phenomenon by which the top-wall admits non-shifted-`k`-th-power tops
(`band-k-weapons.md` W2q, `quantum-band3-cascade.md` §3, commit `9fa9f74`) — by
holding *primality* fixed (`k = 3, 5` prime) and comparing against the sibling's
*composite* band `k = 4` (`S_4 = Φ_2Φ_4` reducible; derived here **independently**,
not coordinated with the sibling).

> **One-line verdict.** Primality of `S_k` does **not** drive the exotic slack.
> *Existence* of exotica is a **bare-`k ≥ 3`** fact (the universal cofactor
> `g = 1 − σ + σ²` works for every `k`, prime or composite); *richness* grows
> smoothly with `k` — the number of minimal-degree exotic tops is
> `1, 4, 13, 40` at `k = 3, 4, 5, 6`, matching `(3^{k−2}−1)/2`, a law **blind to
> primality** (`k = 4` composite and `k = 5` prime both fit). What primality *does*
> control is the **factorization** of the top-divisibility condition (a single
> irreducible `Φ_k | A` when `k` is prime; a conjunction `Φ_d | A`, `d ∣ k`, `d>1`,
> when `k` is composite) and, through **parity**, which specific low-complexity
> families are admissible (the step-2 AP exotic exists iff `k` is *odd*).
> **And the moment-unit kill is indifferent to all of this:** at `k = 5`, every
> one of the (13× richer) exotic families is killed at `Q_0` by the same
> band-agnostic, membership-protected mechanism that killed the unique `k = 3`
> exotic. That indifference is the load-bearing evidence that the **moment-unit
> principle is a candidate law**, not a band-3 coincidence.

Conventions are frozen exactly as in the band-3 corpus (`quantum-band3-cascade.md`
`99fe6ee`, `quantum-exotic-branch.md` `ebfc64d`, `quantum-exotic-closure.md`
`9fa9f74`, hardening `050a4c0`):

- `A₁[x⁻¹] = ⊕_k x^k ℂ[E]`, `E = x∂`, `(x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E)`,
  `f^[r](E) = f(E+r)`.
- ladder-`m` coefficient `Q_m = Σ_{k+l=m}(b_l^[k] a_k − a_k^[l] b_l)`,
  `[D,X] = 1 ⇔ Q_m = δ_{m0}` (band 5: `m ∈ [−10, 10]`).
- genuine `A₁` membership: `E(E−1)⋯(E−r+1) | a_{−r}, b_{−r}` (falling factorial).
- root-necklace calculus: `S` = shift `+1` on root positions, `σ = T⁻¹ = S⁻¹`,
  `S_r = 1 + σ + ⋯ + σ^{r−1}` (necklace operator), `Φ_n` the `n`-th cyclotomic.

Everything displayed is checked exactly by
[`verify_band5_comparison.py`](verify_band5_comparison.py) (72 `PASS` lines, ~30 s;
ends `ALL BAND5 COMPARISON CHECKS PASSED`).

---

## 1. The `k = 5` wall `Q_9` and its necklace lemma

**Top proportionality and gauge.** With `a₅ ≠ 0`, only `(k,l) = (5,5)` reaches
`m = 10`:
```
   Q_10 = b₅^[5] a₅ − a₅^[5] b₅ ,
```
a pure shifted-Wronskian (`W1`), so `b₅/a₅` is 5-periodic hence constant:
`b₅ = λ₅ a₅`. Subtracting `λ₅ X` from `D` realizes the gauge `b₅ = 0`, preserving
`[D,X]=1`, band-5 support, and membership. **[PROVED, `verify §0`.]**

**The wall.** In the gauge `b₅ = 0`, the only surviving pair at `m = 9` is
`(5,4)`:
```
   Q_9 = b₄^[5] a₅ − a₅^[4] b₄ .                    (WALL, band-5 J2q)
```
So `Q_9 = 0` is the staggered homogeneous equation
`b₄(E+5)·a₅(E) = a₅(E+4)·b₄(E)` — a `+5` shift on `b₄` against a `+4` shift on
`a₅`. This is the `k = 5` instance of `W2q` (`u := b_{k−1} = b₄`,
`Q_{2k−1} = u^[k]a_k − a_k^[k−1]u`). **[PROVED identity, `verify §0`.]**

### 1.1 The necklace reduction

Equating leading coefficients is automatic; equating **root multisets** (`A, B`
the position data of `a₅, b₄`) gives, following the band-3 derivation verbatim
with `k = 5`:
```
   S⁻⁵B + A = S⁻⁴A + B  ⇔  (1−S⁻⁴)A = (1−S⁻⁵)B
   ⇔ (×S⁵)  S(S⁴−1)A = (S⁵−1)B  ⇔ (cancel S−1)   S·(1+S+S²+S³)·A = Φ₅(S)·B.   (WALL-M)
```
Since `Φ₅` is a non-zero-divisor in `ℂ[S,S⁻¹]`, `B = S(1+S+S²+S³)A/Φ₅` is
**uniquely determined** by `A`, coset by coset. Equivalently in `σ = S⁻¹`:
```
   S_5(σ)·δ(b₄) = S_4(σ)·δ(a₅) ,     S_5 = Φ_5 = 1+σ+σ²+σ³+σ⁴,  S_4 = 1+σ+σ²+σ³.
```
The two necklace operators are **coprime**, `gcd(S_5, S_4) = 1` in `ℚ[σ]` (their
roots are primitive 5th vs `{2,4}`-th roots of unity — disjoint). Summing
(WALL-M) at `σ = 1` gives the degree law `(k−1) deg a₅ = k deg b₄`, i.e.
**`4 deg a₅ = 5 deg b₄`** (`5 | deg a₅`). **[PROVED structure, `verify §1`.]**

> **Wall Lemma (band-5).** Fix `a₅ ≠ 0`, gauge `b₅ = 0`. A nonzero `b₄` solving
> (WALL) exists **iff**, in every `mod-ℤ` coset, (i) `Φ₅(S) | A` and (ii) the
> forced quotient `B = S(1+S+S²+S³)A/Φ₅` is effective (a genuine multiset). When it
> exists, `b₄` is unique up to a scalar (1-dimensional freedom).

### 1.2 Shifted-5th-power sufficiency, and the exotic gap

If `a₅ = c·h h^[1] h^[2] h^[3] h^[4]` (a **shifted 5th power**), then
`b₄ = κ·h h^[1] h^[2] h^[3]` solves (WALL) — the cofactor `C := A/Φ₅` is effective.
The wall, however, only needs the *image* `B = S·S_4·C` effective, and the map
`C ↦ S·S_4·C` can send a **non-effective** `C` to an effective `B`. That is the
exotic gap. As at every `k ≥ 3`, the **universal cofactor** `g = 1 − σ + σ²`
(`band-k-weapons.md` §W2q, `9fa9f74`) witnesses it at `k = 5`:
```
   g = 1 − σ + σ²   (non-effective),
   S_5 g = 1 + σ² + σ³ + σ⁴ + σ⁶     (effective ⇒ a₅ genuine),
   S_4 g = 1 + σ² + σ³ + σ⁵          (effective ⇒ b₄ genuine).
```
So the exotic (non-shifted-5th-power) branch is **non-empty** at `k = 5`, exactly
as at `k = 3`. `S_5 g` names the top `a₅` with roots `{0,2,3,4,6}` — the
smallest-span `k = 5` exotic (below). **[PROVED, `verify §1`.]**

### 1.3 Structural comparison of the walls (`k = 3` prime, `k = 4` composite, `k = 5` prime)

| `k` | wall `Q_{2k−1}` | necklace | `S_k = (S^k−1)/(S−1)` | factorization |
|---|---|---|---|---|
| 3 (prime) | `b₂^[3]a₃ − a₃^[2]b₂` | `S_3 δ(u) = S_2 δ(a₃)` | `Φ_3 = 1+σ+σ²` | **irreducible** |
| 4 (comp.) | `b₃^[4]a₄ − a₄^[3]b₃` | `S_4 δ(u) = S_3 δ(a₄)` | `Φ_2Φ_4 = (1+σ)(1+σ²)` | **reducible** |
| 5 (prime) | `b₄^[5]a₅ − a₅^[4]b₄` | `S_5 δ(u) = S_4 δ(a₅)` | `Φ_5 = 1+σ+σ²+σ³+σ⁴` | **irreducible** |

`S_k` **always** has all-1 coefficients (it is `1+σ+⋯+σ^{k−1}` by definition);
primality touches only its **factorization**, hence only the *shape of the
divisibility condition* on the top, never the raw operator. **[PROVED / computed,
`verify §2`.]** Independently derived `k = 4` wall `Q_7 = b₃^[4]a₄ − a₄^[3]b₃` and
its necklace `S_4 δ(u) = S_3 δ(a₄)` are machine-checked here from scratch.

---

## 2. Primality dissection: what factorization actually controls

**Prime `k` (3, 5):** wall-divisibility is the **single irreducible** condition
`Φ_k | A`. The minimal cube is the consecutive block `{0,1,…,k−1}` (`A = Φ_k`).

**Composite `k` (4):** `S_4 = Φ_2·Φ_4`, so wall-divisibility is a **conjunction**
`Φ_2 | A` **and** `Φ_4 | A`. This creates a phenomenon **absent at prime `k`**:
tops divisible by one factor but not the other. Concretely, **9** degree-4 tops of
span `≤ 12` are `Φ_4`-divisible but **not** `Φ_2`-divisible — they satisfy "half"
the necklace and admit **no** wall solution. (The step-2 AP `{0,2,4,6}` is the
smallest: `Φ_4 | A` but `A(−1) = 4 ≠ 0`.) **[Computed, `verify §2`.]**

> **This is the sharpest thing primality controls:** at prime `k` the top passes a
> single irreducible cyclotomic gate; at composite `k` it must pass **every**
> `Φ_d`, `d ∣ k`, `d > 1`, simultaneously — and the partial-divisibility tops that
> fail are a composite-only obstruction with no prime analogue.

But — decisively — this factorization does **not** control the *size* or *nature*
of the exotic class, as §3 shows.

---

## 3. Enumeration of minimal-degree exotic tops (Task 2)

Minimal genuine degree is `deg a_k = k` (from `4 deg a₅ = 5 deg b₄`, `5 | deg a₅`).
For single-coset tops (translate `min → 0`), the wall-admissible **exotic**
(non-shifted-power) tops are the cofactors `C = A/S_k` with `C(1)=1`, `C`
non-effective, and both `A = S_k C`, `B = S·S_{k−1}·C` effective. Exhaustive scan
(bounded span, count **stable** — the both-effective condition bounds the Newton
polytope of `C`):

| `k` | minimal exotic tops (up to translation) | `(3^{k−2}−1)/2` |
|---|---|---|
| 3 (prime) | **1** — only `{0,2,4}` (cofactor `Φ_6`) | 1 |
| 4 (comp.) | **4** | 4 |
| 5 (prime) | **13** | 13 |
| 6 (comp.) | **40** (off-verifier, span 22, ~10 s) | 40 |

**[Computed / exhaustive within bounded span; `verify §3` covers `k = 3,4,5`.]**

Two facts jump out.

**(a) Richness is a bare-`k` law, blind to primality.** `k = 3` and `k = 5` are
both prime, yet have `1` vs `13` minimal exotic tops. The counts `1,4,13,40`
interpolate `(3^{k−2}−1)/2` across *both* primes and composites — a smooth law that
does **not** see primality. The exotic slack grows because larger `k` gives the
"both-effective" convolution more Newton-polytope room, not because of any
arithmetic of `S_k`. *(The closed form `(3^{k−2}−1)/2` is **CONJECTURAL** — four
matching data points `k = 3..6`, no proof; presented as a numerical law, not a
theorem.)*

**(b) The exotic class is *tropical*, not cyclotomic.** At `k = 3` the unique
cofactor `Φ_6` is cyclotomic, which could tempt a "cofactors are cyclotomic"
reading. **`k = 5` refutes it:** among the 13 cofactors, several are **not
cyclotomic** — e.g. `{0,1,3,4,7}` has cofactor `S³−S²+1`, irreducible but with no
root of unity. The defining condition is purely Newton-polytope effectivity of
`S_k C` and `S·S_{k−1}·C`, indifferent to whether `C` factors into cyclotomics.
**[Computed, `verify §3`.]**

### 3.1 The AP/tiling motif: step-2 AP admissible **iff `k` odd**

The `k = 3` closure showed the degree-3 exotic branch **is exactly** the step-2 AP
class `{r,r+2,r+4}` (`quantum-exotic-closure.md` Fact 1, `9fa9f74`), via a
domino/triomino tiling of `Φ_3·B = S(1+S)·A`. At `k = 5` the analogue tiles
`Φ_5·B = S·S_4·A` — 4 pentominoes (from `Φ_5·b₄`) against 5 tetrominoes (from
`S·S_4·a₅`) covering the same 20 cells — but the AP is **no longer the whole
branch**: `{0,2,3,4,6}` (differences `2,1,1,2`) is exotic and is not an AP. So the
`k = 3` "AP = branch" collapse is *special to `k = 3`*.

What *does* lift cleanly is the step-2 AP itself, and its admissibility is a clean
**parity** theorem:

> **Step-2 AP theorem.** The step-2 AP `{0,2,…,2(k−1)}` is wall-admissible **iff
> `k` is odd**. Indeed its cofactor is `A/S_k = (S^k+1)/(S+1)`, a polynomial iff
> `S = −1` is a root of `S^k+1`, i.e. iff `(−1)^k = −1`, i.e. `k` odd. For `k` an
> **odd prime** the cofactor is exactly `Φ_{2k}`.

Verified `k = 3, 5, 7` admissible (cofactors `Φ_6, Φ_10, Φ_14`) and `k = 4, 6`
obstructed — the even-`k` failure is precisely the composite `Φ_2`-divisibility
gap of §2 (even spacing kills `A(−1)`). **[PROVED (parity argument) + `verify §3`.]**
So the `k = 3` exotic motif propagates to `k = 5, 7, …` as an *odd*-`k` family, not
a *prime*-`k` family — parity, coarser than primality, is the relevant invariant.

---

## 4. The moment-unit test at `k = 5` (Task 3)

The campaign's candidate uniform mechanism is the **moment-unit principle**: the
central potential `G` (the `W4` moment, `Q_0 = (T−1)G`) is band-agnostic, and
membership pins `G(0) = 0`, so `Q_0 = 1 ⇔ G = E ⇔` the slope `G(1) = 1`; on every
quantum-resistant branch the structure forces slope `0`, so the unit in `[D,X]=1`
is unrealizable. We instantiate it at `k = 5`.

**Band-agnostic core (degree-free).** At `k = 5`:
```
   Q_0 = (T−1)G ,   G = Σ_{k=1}^{5} Σ_{j=0}^{k−1} ( a_k^[j−k] b_{−k}^[j] − b_k^[j−k] a_{−k}^[j] ) ,
```
an exact identity for generic band-5 coefficients; under genuine membership
`G(0) = 0` **identically** (every term of `G(0)` carries a value
`a_{−k}(j)` or `b_{−k}(j)`, `0 ≤ j ≤ k−1`, killed by `E^{underline k}`); and the
**slope** `:=` constant coefficient of `Q_0` `= G(1) − G(0) = G(1)`.
**[PROVED, `verify §4`, band-agnostic — this is the `W4` weapon at `k = 5`.]**

**The kill.** In gauge `b₅ = 0`, the wall gives `b₄`, and the positive cascade
`Q_8, Q_7, Q_6, Q_5, Q_4, Q_3, Q_2, Q_1` forward-solves
`b₃, b₂, b₁, b₀, b_{−1}, b_{−2}, b_{−3}, b_{−4}` from the free data
`(a₄,a₃,a₂,a₁,a₀)` and the membership-protected negatives (with `b_{−5} = μ₅ a_{−5}`
from `Q_{−10}`). We then test the moment. For **every** exotic family tried — a
representative of each distinct cofactor type — at free-degree `d = 1`:
```
   {positive cascade} ∪ {Q_0 = 1}   is INFEASIBLE   (Gröbner basis = [1]),
   {positive cascade} ∪ {Q_0 = 0}   is  FEASIBLE    (Gröbner ≠ [1]).
```
So the **unit is exactly the killer** — the moment carries no unit — and the
exotic branch is empty at that instance. Verified families (`verify §4`):

| top (roots) | cofactor `C` | `a₅(0)` | `d = 1` kill |
|---|---|---|---|
| `{0,2,3,4,6}` | `Φ_6` (universal `g`) | 0 | INFEASIBLE / feasible |
| `{0,2,4,6,8}` | `Φ_10` (step-2 AP) | 0 | INFEASIBLE / feasible |
| `{0,1,4,7,8}` | `Φ_12` (self-reflective) | 0 | INFEASIBLE / feasible |
| `{0,1,3,4,7}` | `S³−S²+1` (**non-cyclotomic**) | 0 | INFEASIBLE / feasible |
| `{1,3,4,5,7}` | `Φ_6`, translated | `−420 ≠ 0` | INFEASIBLE / feasible |

The last row has `a₅(0) ≠ 0`, so the bottom `(μ₅, a_{−5})` **enters** the unit
equation (via `G₅(1) = μ₅ a₅(0) a_{−5}(5)`, the membership-protected level-5 term)
— and still cannot rescue the unit. `{0,2,3,4,6}` is additionally verified at
`d = 2` (INFEASIBLE). **[Computed per instance — each Gröbner `= [1]` is an exact
emptiness *proof* for that top and free-degree; not an arbitrary-degree theorem.]**

**No false kill.** A genuine band-5 pair — `U = x + ∂`, `X = U⁵ − ∂/κ`, `D = κU`,
`[D,X] = 1` (all `Q_m = δ_{m0}`, `a₅ = 1`, `a_{−5} = c·E(E−1)(E−2)(E−3)(E−4)`) —
validates the engine at band 5, and the solver reconstructs its `b₁ = κ` with
**no spurious conditions**. The infeasibility above is a real kill, not a pipeline
artifact; and the `{Q_0 = 0}` feasibility already certifies the solver finds real
solutions when they exist. **[`verify §4`.]**

**The point.** The kill mechanism is **identical** to `k = 3`
(`quantum-exotic-branch.md` §3, `quantum-exotic-closure.md`, `9fa9f74`): the wall
`Q_9` is too weak to force a shifted 5th power, but the moment `Q_0` forbids
everything the wall let through. Crucially it does **not care** that `k = 5` has
`13×` more exotic tops, non-cyclotomic cofactors, or a richer AP structure —
because the kill lives at the moment, driven by the band-agnostic potential and
membership-protection of the extreme, **not** by the wall/top structure. That is
why the moment-unit principle survives intact across the band-3 → band-5 jump
despite the exotic class exploding.

---

## 5. Synthesis across `k` (Task 4): the sharpest supportable conjectures

Two prime bands (`k = 3, 5`, both killed) plus the composite control (`k = 4`) now
constrain the general-`k` picture. The raw material:

- **Existence** of exotica: universal, `k ≥ 3` (cofactor `g = 1 − σ + σ²`). *Not*
  a primality effect.
- **Richness**: `1, 4, 13, 40` at `k = 3..6`, `≈ (3^{k−2}−1)/2`, primality-blind.
- **Nature**: tropical (Newton-polytope), not cyclotomic (non-cyclotomic cofactors
  at `k = 5`).
- **Which families**: parity- and factorization-gated (step-2 AP iff `k` odd;
  composite `k` adds partial-divisibility obstructions).
- **Fate**: killed at `Q_0` in every reached instance, uniformly, by the
  band-agnostic membership-protected moment.

### 5.1 The moment-unit conjecture (sharpest supportable form)

> **Moment-unit principle (conjecture).** For every band `k ≥ 2`, in gauge
> `b_k = 0` with a nonzero level-`k` extreme, the `W4` potential `G` (band-agnostic
> closed form) satisfies `G(0) = 0` under genuine `A₁` membership, and on every
> quantum-resistant branch — i.e. every top **outside** the shifted-`k`-th-power
> (tame) sector — the positive cascade `Q_{2k−2}, …, Q_1` forces the moment slope
> `G(1) = 0`. Hence `Q_0 = 1` (the unit of `[D,X] = 1`) is **unrealizable** off the
> tame sector: *the unit lives only where the top is a shifted `k`-th power.*

Support ledger: the `G`-identity and `G(0) = 0` are **PROVED band-agnostic**
(`W4`, all `k`). The slope-forced-to-`0` is **PROVED degree-free at `k = 2`**
(`quantum-M4.md`, `84978b9`), **closed at `k = 3`** for `d = 1` generic `r` and
per-instance to `d ≤ 4` / `deg a₃ ≤ 9` (`quantum-exotic-closure.md`, `9fa9f74`),
and now **verified per-instance at `k = 5`** across all cofactor-types at `d = 1`
(+ `d = 2`). The uniform-in-`(k, free-degree)` slope-`= 0` statement is the **open
general-`k` theorem** — the `W5` lattice/rank infeasibility applied to the exotic
residual (`band-k-weapons.md` §W5, band-`k` modulus open). The band-5 evidence
**sharpens** the conjecture by removing the last plausible escape hatch: one might
have feared the kill was an artifact of the `k = 3` exotic class being a *single*
AP; band 5 shows it holds against a `13×`-richer, non-cyclotomic, mixed-parity
exotic class — so the kill is genuinely a property of the *moment*, not of the
*top*.

### 5.2 The exotic-class structural conjecture (sharpest supportable form)

> **Exotic-class structure (conjecture).** At minimal degree `deg a_k = k`, the
> wall-admissible exotic tops (single coset, up to translation) are exactly the
> cofactors `C = A/S_k` (`C(1) = 1`, `C` non-effective) for which **both** `S_k C`
> and `S·S_{k−1}·C` are effective — a purely **tropical** condition on the Newton
> polytope of `C`, independent of its cyclotomic factorization. Their number is
> finite (the both-effective condition bounds the polytope) and grows as
> `(3^{k−2}−1)/2` `= 1,4,13,40,…` (**numerically conjectural**, `k = 3..6`),
> a bare-`k` law **independent of the primality of `k`**. The **shape** of the
> class is gated by the factorization of `S_k`: the universal member `g = 1−σ+σ²`
> exists for all `k ≥ 3`; the step-2 AP member (cofactor `Φ_{2k}` for `k` odd
> prime) exists iff `k` is odd; and composite `k` admits partial-`Φ_d`-divisibility
> obstructions with no prime analogue.

### 5.3 Answer to the driving question

> **Does primality of `(S^k−1)/(S−1)` control the shape of the exotic class? — No,
> not the slack itself.** Primality (⇔ irreducibility of `S_k`) controls the
> *factorization* of the top-divisibility gate and, through parity, *which*
> low-complexity families appear; it leaves the *existence*, the *count*, and the
> *tropical nature* of the exotic class untouched. Those are **bare-`k ≥ 3`**
> phenomena. And the campaign's real target — whether the exotic branch survives —
> is decided *downstream of all of this*, at the moment `Q_0`, by a mechanism that
> is blind to primality, parity, count, and cofactor-type alike. The exotic slack
> is a red herring for DC1 in exactly the sense band 3 first suggested: **the wall
> is loose, the moment is tight, and the moment does not care how loose the wall
> got.**

---

## 6. Rigor ledger

**PROVED (band-agnostic / structural, machine-checked at `k = 5` with `k = 3,4` controls):**
- the wall `Q_9 = b₄^[5]a₅ − a₅^[4]b₄` and top `Q_10` proportionality (`§1`, `verify §0`);
- the necklace reduction (WALL-M) `S·S_4·A = Φ_5·B`, `gcd(S_5,S_4)=1`, the degree
  law `4 deg a₅ = 5 deg b₄`, shifted-5th-power sufficiency, and the exotic gap via
  the universal cofactor `g = 1−σ+σ²` (`§1`, `verify §1`);
- `S_k` irreducible ⇔ `k` prime; the composite-`k` conjunction-of-divisibilities
  and its partial-divisibility obstruction (`§2`, `verify §2`);
- the step-2 AP admissibility theorem (admissible **iff `k` odd**; cofactor `Φ_{2k}`
  for odd prime `k`) (`§3.1`, `verify §3`);
- the `W4` moment core at `k = 5`: `Q_0 = (T−1)G`, `G(0) = 0` under membership,
  `slope = G(1)` (`§4`, `verify §4`) — degree-free.

**COMPUTED (exact, bounded — corroboration or per-instance proof, not
arbitrary-degree):**
- the minimal-exotic-top **counts** `1, 4, 13` (`k = 3,4,5`, stable span-scan;
  `k = 6 → 40` off-verifier) and the presence of **non-cyclotomic** cofactors at
  `k = 5` (`§3`, `verify §3`);
- the **moment kill**: `{cascade} ∪ {Q_0 = 1}` infeasible (Gröbner `= [1]`) while
  `{Q_0 = 0}` feasible, for 5 exotic families at `d = 1` and one at `d = 2` — each
  Gröbner `= [1]` an **exact emptiness proof** for that top and free-degree (`§4`,
  `verify §4`);
- the positive control (genuine band-5 pair; no spurious conditions) (`§4`, `verify §4`).

**CONJECTURED:**
- the closed form `(3^{k−2}−1)/2` for the minimal-exotic-top count (4 data points,
  no proof — a numerical law, possibly spurious beyond `k = 6`);
- the **moment-unit principle** at general `k` (§5.1): the uniform-in-`(k,
  free-degree)` slope-`= 0` statement — the `W5` band-`k` modulus, open;
- the **exotic-class structural** conjecture (§5.2).

**NOT claimed:** any band-5 theorem; any DC1/JC2 statement (no counterexample is
produced — the exotic branch is empty everywhere the argument reaches, consistent
with DC1); closure of the `k = 5` exotic branch uniform in free degree or across
all 13 families for degree `> 1`; anything about `deg a₅ > 5` exotic tops.

---

## 7. Verification

```sh
uv run --with sympy python research/dc1-program/verify_band5_comparison.py
```
Exact SymPy; 72 `PASS` lines, ~30 s; ends `ALL BAND5 COMPARISON CHECKS PASSED`.
`§0` engine (`Q_m` = crossed-product commutator, all `m ∈ [−10,10]`; `Q_10`, `Q_9`),
`§1` the necklace lemma + exotic gap, `§2` the primality dissection + independent
`k = 4` wall, `§3` the exotic-top enumeration (counts, cofactor nature, step-2 AP
parity), `§4` the moment-unit kill + positive control. The arbitrary-`k`
completeness of the THEOREM items is the written structural argument (necklace
divisibility, tiling, telescoping, membership), not the finite-`k` computation;
the counts and the moment kill are exact but bounded (span / free-degree).
```
Citations pinned: quantum-band3-cascade.md `99fe6ee`, quantum-exotic-branch.md
`ebfc64d`, quantum-exotic-closure.md + hardening `9fa9f74`/`050a4c0`,
band-k-weapons.md `9fa9f74`, quantum-M4.md `84978b9`.
```
