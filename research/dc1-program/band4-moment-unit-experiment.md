# The band-4 moment-unit stress test: the composite cyclotomic wall opens a wider exotic class, and the moment kills all of it

**INDEPENDENTLY DERIVED AND MACHINE-VERIFIED — NOT PEER REVIEWED — BAND-SCOPED**

DC1 campaign, Wave α, **BAND-4 EXPERIMENT**. This memo runs the decisive
moment-unit stress test one band above the band-3 exotic closure
(`../band3/quantum-exotic-branch.md`, `ebfc64d`;
`../band3/quantum-exotic-closure.md`, `9fa9f74`). The target is DC1 (every
endomorphism of `A₁` is an automorphism), attacked via the band filtration; the
candidate uniform mechanism under test is the **moment-unit principle** — the
central potential `G` with `Q₀ = (T−1)G` is band-agnostic (`W4`,
`../band3/band-k-weapons.md`), membership pins `G(0) = 0`, so `Q₀ = 1 ⇔ G = E` is
a **slope** statement, and in every quantum-resistant branch so far the branch
structure forces slope `0`, making the unit in `[D,X] = 1` unrealizable.

Band 4 is the sharpest test yet because the necklace factor
```
   (S⁴−1)/(S−1) = S₄ = 1+S+S²+S³ = (1+S)(1+S²) = Φ₂ · Φ₄
```
is **composite** — unlike band 3's irreducible `Φ₃`. The pre-registered
hypothesis was that composite cyclotomics might open *wider* exotic classes (root
data divisible by "either factor pattern"). This memo settles that hypothesis and
runs the experiment.

> **Verdict.** The moment-unit principle **SURVIVES at band 4.** Every minimal
> exotic top (the `Q₇`-wall-admissible, non-shifted-4th-power `a₄`) is **KILLED at
> `Q₀`**: `{positive cascade} ∪ {Q₀ = 1}` is infeasible (exact Gröbner `= [1]`)
> while `{positive cascade} ∪ {Q₀ = 0}` is feasible — the moment **unit** is the
> exact killer, precisely as at band 3. **No exotic top admits `Q₀ = 1`; no DC1
> counterexample candidate arises.** The composite cyclotomic changes the *breadth*
> of the exotic class (band 4's minimal exotic tops are **4** distinct necklaces in
> **3** reflection classes — NOT a single AP family as at band 3), but it changes
> **nothing** about the kill: the same slope-forcing structure closes every branch.

The "either factor" form of the hypothesis is **false**: because
`gcd(Φ₃, Φ₂Φ₄) = 1`, the wall forces the **full** `S₄ | A`, never just one factor
(machine-verified). The genuinely wider phenomenon is a larger *effectivity gap*
in the wall's forced quotient, not a relaxed divisibility.

Everything below is checked exactly by
[`verify_band4_experiment.py`](verify_band4_experiment.py) (ends
`ALL BAND4 EXPERIMENT CHECKS PASSED`).

## 0. Setup and the descent (frozen conventions)

`A₁[x⁻¹] = ⊕_k x^k ℂ[E]`, `(x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E)`,
`f^[r](E) = f(E+r)`, `E = x∂`. Band-4 supports:
`X = Σ_{k=−4}^{4} x^k a_k(E)`, `D = Σ_{l=−4}^{4} x^l b_l(E)`. The ladder-`m`
coefficient of `[D,X]` is
```
   Q_m = Σ_{k+l=m} ( b_l^[k] a_k − a_k^[l] b_l ),    [D,X]=1 ⇔ Q_m = δ_{m0},  m ∈ [−8,8].
```
Genuine `A₁` membership: `E(E−1)⋯(E−r+1) | a_{−r}, b_{−r}` for `r = 1..4`.
`verify §0` checks `Q_m` equals the direct crossed-product commutator for every
`m ∈ [−8,8]` with generic degree-2 coefficients.

**Top and gauge.** `Q₈ = b₄^[4] a₄ − a₄^[4] b₄`; assuming `a₄ ≠ 0`, `b₄/a₄` is
4-periodic, hence constant, so `b₄ = λ₄ a₄`. Gauge `D ↦ D − λ₄ X` to set
`b₄ = 0`. In this gauge the "new" `b`-coefficient at level `Q_m` is `b_{m−4}`,
isolated by the **descent operator** `L_m[b] := b^[4] a₄ − a₄^[m−4] b` (the
`(4, m−4)` pair), with `Q_m = L_m[b_{m−4}] + (lower pairs)` for
`m = 7,6,…,0` (`verify §0`).

**The moment.** For band `K = 4`, `Q₀ = (T−1)G` with the band-agnostic potential
(`W4`; `verify §0` checks it in both the `W4` and the closure forms, and their
equality)
```
   G = Σ_{k=1}^{4} Σ_{j=0}^{k−1} ( a_k^[j−k] b_{−k}^[j] − b_k^[j−k] a_{−k}^[j] ).
```
Under membership `G(0) = 0` **identically**, so `Q₀ = 1 ⇔ G = E`, and the single
number that decides it is the **slope** `= const. coeff. of Q₀ = G(1)`
(`verify §0`). The experiment is: *can the constrained bilinear `G` have
`G(1) = 1`?*

## 1. The `Q₇` wall — the k = 4 necklace lemma (composite `Φ₂Φ₄`)

In gauge `b₄ = 0` the only pair at `m = 7` is `(4,3)`:
```
   Q₇ = b₃^[4] a₄ − a₄^[3] b₃        (verify §1).
```
So `Q₇ = 0` is the staggered homogeneous **wall** `b₃(E+4)·a₄(E) = a₄(E+3)·b₃(E)`.

### 1.1 Necklace reduction

Per `mod-ℤ` coset, with `A, B` the root-multiset data of `a₄, b₃` (`S` = shift
`+1`): equating root multisets and cancelling the common `(S−1)` gives
```
   S₄(S)·B = S·S₃(S)·A ,     S₄ = 1+S+S²+S³ = Φ₂Φ₄,   S₃ = 1+S+S² = Φ₃.   (WALL-M)
```
The forced quotient is `B = S·S₃·A / S₄`. Since **`gcd(S₃, S₄) = 1`** in
`ℚ[S]` (`verify §1`), `S₄ | S₃A ⇔ S₄ | A`: a nonzero `b₃` exists **iff**, in
every coset, (i) `S₄ = Φ₂Φ₄` divides `A`, and (ii) the forced quotient
`B = S·S₃·A/S₄` is effective (nonnegative). Summing (WALL-M) at `S = 1`:
`4·deg b₃ = 3·deg a₄` (`verify §1`), so the minimal nonzero top is
`deg a₄ = 4`, `deg b₃ = 3`.

> **The "either factor" hypothesis is FALSE.** One might hope the composite
> `S₄ = Φ₂Φ₄` lets `A` be divisible by only *one* factor (`Φ₂` **or** `Φ₄`). It
> cannot: (WALL-M) plus `gcd(Φ₃, Φ₂Φ₄) = 1` forces the **full** `S₄ | A`. The
> complete deg-4 scan (`verify §2`) confirms every wall-admissible top has both
> `Φ₂ | A` and `Φ₄ | A`. The composite structure does **not** relax the
> divisibility condition (i).

### 1.2 Shifted-4th-power sufficiency; the effectivity gap

If `a₄ = c·h h^[1] h^[2] h^[3]` (shifted 4th power) then
`b₃ = κ·h h^[1] h^[2]` solves the wall (`verify §1`, symbolic `h`), with cofactor
`C := A/S₄` effective. The **exotic** tops are exactly those where `C` is *not*
effective (not a shifted 4th power) yet `B = S·S₃·C` *is* effective — the map
`C ↦ S·S₃·C` can carry a non-effective `C` to an effective `B`. This is where
band 4 genuinely differs from band 3: the longer blocks `S₄` (span 4) and
`S₃` (span 3) admit a **larger** effectivity gap than band 3's `Φ₃`/`S₂`.

## 2. Enumeration of the minimal exotic tops — no AP collapse

A complete scan of degree-4 single-coset tops (`verify §2`; span bounded by the
12-cell tiling below) gives **exactly 5** wall-admissible tops up to translation:
one shifted 4th power `{0,1,2,3}` and **four EXOTIC** tops:

| exotic `a₄` (roots) | `b₃` (roots) | cofactor `A/S₄` | reflection class |
|---|---|---|---|
| `{0,2,3,5}` | `{1,3,5}` | `Φ₆ = 1−S+S²` | palindromic |
| `{0,1,3,6}` | `{1,2,6}` | `1−S²+S³` | pair with `{0,3,5,6}` |
| `{0,3,5,6}` | `{1,5,6}` | `1−S+S³` | pair with `{0,1,3,6}` |
| `{0,3,6,9}` | `{1,5,9}` | `Φ₆·Φ₁₂` | palindromic (step-3 AP) |

So the four exotic tops form **3 reflection classes** (`verify §2`; reflection is
the quantum Fourier involution `E ↦ −E−1`). This is the central structural
finding:

> **No AP collapse at band 4.** At band 3 the degree-3 exotic branch collapsed to
> the *single* step-2 AP family `{r, r+2, r+4}`
> (`../band3/quantum-exotic-closure.md` Fact 1). At band 4 the minimal exotic
> branch does **not** collapse to one AP family: **3 of the 4 exotic tops are
> non-arithmetic-progressions** (`verify §2`). Only `{0,3,6,9}` is an AP (step 3).
> The composite `Φ₂Φ₄` opens a **wider** exotic class — this is the true content
> of the "composite ⇒ wider" hypothesis (via the effectivity gap of §1.2, **not**
> via relaxed divisibility).

Note `{0,2,3,5}` re-uses the band-3 exotic cofactor `Φ₆ = 1−S+S²`, and the
step-3 AP `{0,3,6,9}` has the cyclotomic cofactor `Φ₆·Φ₁₂` (`verify §2`).

### 2.1 The tiling picture (k = 4: triomino/tetromino)

(WALL-M) reads, cell by cell, as a multiset tiling: `S₄·B` places a **tetromino**
(span-4 block) at each root of `b₃` — **3 tetrominoes** — and `S·S₃·A` places a
**triomino** (span-3 block) at each root of `a₄` shifted by `1` — **4 triominoes**
— covering the same 12-cell multiset (`verify §2`). This is the band-4 analogue
of band 3's "2 triominoes = 3 dominoes" (`../band3/quantum-exotic-closure.md`
§2). The band-3 tiling had a **unique** nontrivial solution (the AP); the band-4
tiling **3 tetrominoes = 4 triominoes** has **four** (three reflection classes) —
the combinatorial source of the wider class. Because `lcm(3,4) = 12` admits no
proper sub-block splitting, the tiled region is a single block, bounding the span
(max exotic span `= 9`, the step-3 AP) and making the finite scan complete.

## 3. THE EXPERIMENT — the moment kills every exotic top

For each exotic top we fix `a₄` and its wall solution `b₃` (normalized
`κ₃ = 1`, i.e. `b₃ ≠ 0` — the genuinely exotic branch), take **free lower data**
`a₃, a₂, a₁, a₀` and membership-carrying negatives `a_{−1}, a_{−2}, a_{−3},
a_{−4}` (degree `d`), with `b_{−4} = μ₄ a_{−4}` from `Q_{−8}`, and forward-solve
the **positive cascade**
```
   Q₆ → b₂,  Q₅ → b₁,  Q₄ → b₀,  Q₃ → b_{−1},  Q₂ → b_{−2},  Q₁ → b_{−3}
```
subject to solvability conditions on the free data. This cascade is **solvable**
(`verify §3`: explicit feasible point for `{0,2,3,5}`; `{cascade} ∪ {Q₀ = 0}`
feasible for every top) — so, exactly as at band 3, the positive side alone does
**not** kill the exotic top. The kill is downstream, at `Q₀`.

> **Result (per-instance emptiness proofs).** For **every** minimal exotic top,
> `{positive cascade} ∪ {Q₀ = 1}` is **infeasible** (exact Gröbner `= [1]`), while
> `{positive cascade} ∪ {Q₀ = 0}` is **feasible**. Verified at free-degree
> `d = 1` and `d = 2` for all four exotic tops (`{0,2,3,5}, {0,1,3,6}, {0,3,5,6},
> {0,3,6,9}`), and for integer-shift instances `r ∈ {0,1,−1,2}` of the step-3 AP
> family `{r,r+3,r+6,r+9}` at `d = 1` (`verify §3`). The **moment unit** is
> exactly the killer.

Each `Gröbner = [1]` is an exact per-instance emptiness *proof* (not mere
evidence): the ideal generated by the cascade conditions and the `Q₀ = 1`
coefficient system is the unit ideal, so the variety is empty over any field
containing `ℚ`. The contrast with `{Q₀ = 0}` (feasible) localizes the obstruction
to the `+1`.

### 3.1 The `b₃ = 0` sub-branch

Setting `κ₃ = 0` (`b₃ = 0`) for an exotic `a₄` is a separate degenerate branch.
The clean structural fact is the **`L₀` degree obstruction** (`verify §4`,
symbolic in `deg b`): for `deg a₄ = 4`,
```
   L₀[b] = b^[4] a₄ − a₄^[−4] b   has degree deg b + 3 with leading coefficient
   4(4 + deg b)·lc(a₄)·lc(b) ≠ 0,
```
so `L₀[b]` is either `0` (`b = 0`) or of degree `≥ 3` — it can **never** equal the
nonzero constant `1`. Hence a *fully collapsed* `b₃ = 0` branch dies at
`Q₀ = L₀[b_{−4}] = 1`, arbitrary degree — the direct analogue of band 3's
Lemma 2 (`../band3/quantum-exotic-branch.md` §2), whose obstruction was `3(3+q)`.

**Honest caveat (a genuine band-4 novelty).** Unlike band 3, the `b₃ = 0`
collapse is **incomplete** at band 4: the operator `L₆` (the `(4,2)` pair,
shift `+2`) has a **nontrivial kernel** for the exotic top, because e.g.
`{0,2,3,5}` has root data `A = (1+S²)(1+S³)`, and `b₂ = (E−2)(E−5)` solves
`b₂^[4] a₄ = a₄^[2] b₂` (`verify §4`). So `b₂` is not forced to `0` and the clean
"collapse then `L₀`" route does not finish the `b₃ = 0` branch by itself. We
therefore close it, too, by direct Gröbner: for every exotic top,
`{b₃ = 0 cascade} ∪ {Q₀ = 1}` is infeasible and `∪ {Q₀ = 0}` feasible
(`verify §4`).

### 3.2 No counterexample; anti-false-kill control

No exotic top admits `Q₀ = 1`, so **no DC1/JC2 counterexample candidate arises**
and none is claimed; the result is consistent with DC1. To rule out a pipeline
artifact, `verify §5` runs a **genuine band-4 pair** through the same cascade:
```
   U = x + ∂,   X = U⁴ − ∂,   D = U,     [D,X] = −[U,∂] = 1     (a₄ = 1, trivial 4th power),
```
with membership `E(E−1)(E−2)(E−3) | a_{−4}` holding. The solver reconstructs `D`
exactly (`b₁ = 1`, `b_{−1} = E`, all else `0`), every emitted condition involves
only the `b`-kernel freedom (no false `a`-obstruction), and `Q₀ = 1` **holds** on
the reconstructed pair. The pipeline detects genuine feasibility — the exotic
infeasibility is real.

## 4. The kill mechanism vs band 3 — the composite changes nothing

The band-4 gatekeeper now reads, extending the band-3 table
(`../band3/quantum-exotic-closure.md` §7):

| face / sector | mechanism at `Q₀` |
|---|---|
| band-3 quantum exotic `b₂ = 0` | `L₀[b_{−3}] = 1` impossible: `deg ≥ 2` or `0` (obstruction `3(3+q)`) |
| band-3 quantum exotic `b₂ ≠ 0` | moment slope `G(1)` forced to `0`: the moment carries no unit |
| **band-4 quantum exotic `b₃ = 0`** | `L₀[b_{−4}] = 1` impossible: `deg ≥ 3` or `0` (obstruction `4(4+q)`); collapse incomplete, Gröbner-closed |
| **band-4 quantum exotic `b₃ ≠ 0`** | moment slope `G(1)` forced to `0`: `{cascade} ∪ {Q₀=1}` Gröbner `= [1]` |

The killing equation is the **same** — the `W4` moment `Q₀`, whose unit is the
slope `G(1)` — and the reason is the **same** membership-protected mechanism: a
membership-protected extreme lets the moment carry only its natural low-order
content (`E`), never a residual unit. The composite cyclotomic `Φ₂Φ₄`:

- **widens the wall's exotic input** (4 tops / 3 reflection classes vs band 3's
  single AP family) — via the effectivity gap, not divisibility;
- **does not widen the divisibility** (full `S₄ | A` forced; "either factor" is
  false);
- **does not weaken the kill** (every wider branch dies at the identical moment
  unit).

So the moment-unit principle is **band-agnostic** in exactly the sense the
campaign conjectured: the wall's arithmetic gets richer with `k`, but the `W4`
central integral `Q₀ = (T−1)G` with `G(0) = 0` forbids the unit regardless. This
is the band-4 rung on the DC1 face, and it is the first datum showing the
principle survives a **composite** cyclotomic wall.

## 5. Status of claims (proved / verified / residual)

**PROVED (arbitrary degree, machine-checked identities):**
- the `Q_m` engine and gauge/descent isolation `Q_m = L_m[b_{m−4}] + lower`
  (§0); `Q₀ = (T−1)G` with the band-agnostic `G`, `G(0) = 0` under membership,
  and `slope = G(1)` (§0);
- the `Q₇` wall necklace reduction `S₄B = S·S₃·A`, the composite factorization
  `S₄ = Φ₂Φ₄`, `gcd(Φ₃, Φ₂Φ₄) = 1` forcing full `S₄ | A` (**"either factor"
  false**), shifted-4th-power sufficiency, and the `4 deg b₃ = 3 deg a₄` degree
  law (§1);
- the `b₃ = 0` `L₀` degree obstruction `4(4+q)` and the `L₆`-kernel obstruction
  to a clean collapse (§3.1).

**VERIFIED (exact, bounded — per-instance emptiness PROOFS, not arbitrary-degree):**
- the complete degree-4 single-coset enumeration: exactly 5 realizable tops, 4
  exotic in 3 reflection classes, **no AP collapse**; the `3-tetromino =
  4-triomino` tiling (§2);
- the moment-unit kill: `{cascade} ∪ {Q₀ = 1}` infeasible (Gröbner `= [1]`) vs
  `{Q₀ = 0}` feasible, for all 4 exotic tops at `d = 1, 2`, the step-3 AP
  instances `r ∈ {0,±1,2}` at `d = 1`, and the `b₃ = 0` sub-branch of every
  exotic top (§3, §3.1);
- the positive control: a genuine band-4 pair reproduced with `Q₀ = 1` holding —
  no false kill (§3.2).

**RESIDUAL (honestly open, precisely delimited):**
- the **uniform-in-`r` and uniform-in-`d`** closure. Each instance is a Gröbner
  proof, but no single closed-form slope certificate uniform in the shift `r` or
  the free-degree `d` is extracted here (the symbolic-`r` elimination that closed
  band 3 at `d = 1` is computationally much heavier at band 4 and is not carried
  out). This is the same status frontier as the band-3 closure's uniform-degree
  residual (`../band3/quantum-exotic-closure.md` §8).
- the **non-single-coset and `deg a₄ ≥ 8`** exotic tops (multi-coset products and
  higher-degree necklaces), untreated here — the band-4 analogue of band 3's
  `deg a₃ ≥ 6` non-AP residual.

**NOT claimed:** any DC1/JC2 statement (no counterexample is produced, consistent
with DC1); a full band-4 theorem; closure of the exotic branch uniform in `r`/`d`;
any statement about non-exotic (shifted-4th-power) tops beyond the wall/moment
reduction.

## 6. Verification

```sh
uv run --with sympy python research/dc1-program/verify_band4_experiment.py
```
Exact SymPy; §0 engine (`Q_m` = commutator, `m ∈ [−8,8]`; descent operators;
`Q₀ = (T−1)G`; `G(0) = 0`; `slope = G(1)`), §1 the `Q₇` wall (composite `Φ₂Φ₄`,
full `S₄ | A`, sufficiency, degree law), §2 the complete deg-4 enumeration (5
tops, 4 exotic / 3 reflection classes, no AP collapse, the tetromino/triomino
tiling, `Φ₆`/`Φ₆Φ₁₂` cofactors), §3 the experiment (`{cascade} ∪ {Q₀ = 1}`
infeasible vs `{Q₀ = 0}` feasible, `d = 1, 2` + AP instances), §4 the `b₃ = 0`
sub-branch (`L₀` obstruction, `L₆` kernel, Gröbner kill), §5 the genuine
positive control. A successful run ends `ALL BAND4 EXPERIMENT CHECKS PASSED`.
