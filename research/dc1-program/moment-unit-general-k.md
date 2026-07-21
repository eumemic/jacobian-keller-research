# The moment-unit principle at arbitrary band `k` (the DC1 induction step)

**INDEPENDENTLY DERIVED — MACHINE-CHECKED IDENTITIES (k = 2,3,4,5,6,7) — NOT PEER
REVIEWED — BAND-SCOPED.**

This memo attempts the **moment-unit theorem at symbolic band index `k`** — the
statement that would furnish the quantum resistant-branch rung of the width
(band-filtration) induction toward **DC1** (every endomorphism of `A₁` is an
automorphism). It **proves the four band-agnostic pillars of the mechanism** (the
reformulation, the moment-slope formula, the wall/necklace classification, the
tropical filler skeleton) as arbitrary-`k` theorems, isolates the **single
residual lemma** on which a full band-`k` theorem hangs, and **machine-verifies
that residual through band 5** — crucially, past the composite-cyclotomic
threshold `k = 4` (`S₄ = Φ₂Φ₄`) where the band-3 proof's central structural
crutch (the AP collapse) provably fails. The honest deliverable is a **precise
reduction**: *moment-unit at band `k` holds provided Lemma S (the low-end
slope-forcing certificate)*, with Lemma S proved at `k = 3`, proved for a whole
`k = 4` translation family at once, and verified at `k = 4, 5`.

Conventions are frozen exactly as in the quantum corpus: `A₁[x⁻¹] = ⊕_k x^k ℂ[E]`,
`E = x∂`, `(x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E)`, `f^[r](E) = f(E+r)`,
```
   X = Σ_{i=-k}^{k} x^i a_i(E),   D = Σ_{l=-k}^{k} x^l b_l(E),
   Q_m = Σ_{i+l=m} ( b_l^[i] a_i − a_i^[l] b_i ) = δ_{m0}   ⇔   [D,X] = 1,
```
membership `E(E−1)⋯(E−r+1) | a_{−r}, b_{−r}`, gauge `b_k = 0` (top proportionality
`W1`), bottom proportionality `b_{−k} = μ_k a_{−k}` (from `Q_{−2k}`). Citations are
commit-pinned: band-2 assembly `84978b9`; band-3 reconnaissance `99fe6ee`;
band-`k` weapons + Wave-B exotic branch `ebfc64d`; induction-weapon hardening
`050a4c0`; Wave-C exotic closure `9fa9f74` (current HEAD). `[WEAP]` =
`research/band3/band-k-weapons.md`; `[CASC]` = `research/band3/quantum-band3-cascade.md`;
`[EXO]` = `research/band3/quantum-exotic-branch.md`; `[CLO]` =
`research/band3/quantum-exotic-closure.md`.

Everything displayed is checked exactly by
[`verify_moment_unit_k.py`](verify_moment_unit_k.py) (81 `PASS` lines, ~30 s;
ends `ALL MOMENT UNIT K CHECKS PASSED`).

---

## 0. The sharpened target

> **Moment-Unit Theorem at band `k` (target).** Let `[D,X] = 1` in `A₁[x⁻¹]` have
> supports in `[−k,k]` and genuine `A₁` membership. Gauge `b_k = 0`. Suppose
> `a_k ≠ 0` solves the top wall `Q_{2k−1} = 0` but is **not** a shifted `k`-th
> power `c·h h^[1]⋯h^[k−1]` — i.e. `a_k` lies in the **exotic (root-necklace)**
> class, its cofactor `g` non-effective (`[WEAP]` §W2q). Then the central equation
> `Q₀ = 1` is **unsatisfiable**: the potential's slope `G(1)` is forced to `0`.

Combined with the shifted-power sectors (closed at `k ≤ 3`; the band-2/3 pattern
`[CASC]` §6, `[EXO]` §2) and the `a_k = 0` orientation reductions to lower bands,
the Moment-Unit Theorem at all `k` is the quantum band-`k` induction step; the
band-filtration then gives DC1. **This memo does not prove the target at all `k`;
it proves the reduction and verifies the residual through `k = 5`.** No
DC1/JC2 statement is claimed and no counterexample is produced (consistent with
DC1).

The tool `[WEAP]` §W4q already proves, band-agnostically, that `Q₀ = (T−1)G` with
the explicit staggered potential
```
   G = Σ_{i=1}^{k} Σ_{r=1}^{i} ( a_i^[−r] b_{−i}^[i−r] − a_{−i}^[i−r] b_i^[−r] ),   (W4q)
```
and that membership forces `G(0) = 0`, whence `Q₀ = 1 ⇔ G = E`. So the entire
target is the **single number** `G(1)` — the *slope* of the central potential,
equal to the constant coefficient of `Q₀` (`verify §0–§1`). The theorem asserts
this slope cannot be `1` on the exotic variety; equivalently it is forced to `0`.

---

## 1. Route 1 — Lemma P: the general moment-slope formula (PROVED, any `k`)

The band-2/3 experience was that "membership protects the deepest levels." We make
this exact and band-agnostic. It is the load-bearing structural theorem of the
memo.

> **Lemma P (moment-slope formula).** Under genuine membership, the slope of the
> `W4q` potential collapses from a double sum to **one product per side per level**:
> ```
>    G(1) = Σ_{i=1}^{k} [ a_i(0) · b_{−i}(i)  −  a_{−i}(i) · b_i(0) ].
> ```
> In the gauge `b_k = 0` with `b_{−k} = μ_k a_{−k}`, the level-`k` term is
> `μ_k · a_k(0) · a_{−k}(k)`, so
> ```
>    G(1) = Σ_{i=1}^{k−1} [ a_i(0) b_{−i}(i) − a_{−i}(i) b_i(0) ]  +  μ_k a_k(0) a_{−k}(k).
> ```

*Proof.* Evaluate `(W4q)` at `E = 1`. The level-`i` block's first sum is
`Σ_{r=1}^{i} a_i(1−r) · b_{−i}(1+i−r)`. Membership gives `b_{−i}(n) = 0` for
`n ∈ {0,1,…,i−1}`; the argument `n = 1+i−r` lands in that killed range for every
`r = 2,…,i` and escapes it (`n = i`) only at `r = 1`. So the first sum contributes
the **single** term `a_i(0) b_{−i}(i)`. The second sum
`Σ_{r=1}^{i} a_{−i}(1+i−r) b_i(1−r)` is annihilated identically by membership on
`a_{−i}` except at `r = 1`, contributing `a_{−i}(i) b_i(0)`. Summing `i = 1..k`
with the block sign gives the formula. The gauged form uses `b_k(0) = 0` and
`b_{−k}(k) = μ_k a_{−k}(k)`. ∎

This is a genuine arbitrary-`k` theorem (the membership argument is uniform in `k`);
`verify §1` confirms it as an identity at `k = 2,3,4,5`, both raw and gauged. It
**recovers and generalizes** the band-3 fact `[CLO]` §3 / `[EXO]` §3.2
(`G₃(1) = μ₃ a₃(0) a_{−3}(3)`) and exhibits it as the `i = k` instance of a clean
general law. Two consequences the formula makes transparent:

- **Membership protection is exact and total at the deep levels.** Every level-`i`
  bottom coefficient `a_{−i}, b_{−i}` enters the slope through its **first
  membership-liberated value** `(i)` only — the entire double ladder `r = 2..i`
  is killed. This is precisely the `a_{−k}(1) = … = a_{−k}(k−1) = 0` protection
  the band-2/3 corpus observed, now proved to isolate one value per coefficient.
- **The slope is bilinear top×bottom, level-diagonal.** `G(1)` pairs the value of
  each top coefficient at `0` with the value of the mirror bottom coefficient at
  `i`. No cross-level products survive. This is what makes the tropical skeleton
  (Route 3) triangular.

---

## 2. Route 2 — the wall/necklace classification and the AP-collapse's failure

The exotic branch is governed by the top wall `Q_{2k−1} = u^[k] a_k − a_k^[k−1] u`
(`u := b_{k−1}`). Passing to root multisets (`δ(f^[r]) = S^{−r}δ(f)`, `S` = shift
on root positions) reduces the wall to the **necklace identity**
```
   S_k(S) · δ(u) = S · S_{k−1}(S) · δ(a_k),     S_r(S) := 1 + S + ⋯ + S^{r−1}   (WALL-N)
```
(`verify §2`; this is `[WEAP]` §W2q / `[CASC]` §3.1, re-derived here for general
`k`). Because `gcd(S_k, S·S_{k−1}) = 1` in `ℚ[S]`, a nonzero `u` exists **iff**
`δ(a_k) = S_k · g`, `δ(u) = S_{k−1}·g` (up to overall shift) for a cofactor
`g ∈ ℤ[S]` with **both** `S_k g` and `S_{k−1} g` effective. *Shifted `k`-th power*
= `g` itself effective; *exotic* = `g` non-effective but the two products
effective. This is the exact necklace slack.

**Three arbitrary-`k` facts (`verify §2`):**

1. **The `S_k` prime/composite dichotomy.** `S_k = ∏_{d|k, d>1} Φ_d`. For **prime**
   `k`, `S_k = Φ_k` is a single cyclotomic (`S_2=Φ_2, S_3=Φ_3, S_5=Φ_5, S_7=Φ_7`);
   for **composite** `k` it factors (`S_4 = Φ₂Φ₄`, `S_6 = Φ₂Φ₃Φ₆`). The composite
   case is the new slack structure flagged for the campaign — and the moment-unit
   kill (Route 4) is verified to survive it at `k = 4`.
2. **A universal exotic cofactor.** `g = 1 − S + S²` is non-effective yet makes
   `S_k g` and `S_{k−1} g` effective for **every** `k ≥ 3` (`[WEAP]` §W2q). It
   yields the canonical exotic top `δ(a_k) = {0, 2, 3, …, k−1, k+1}` with wall
   solution `δ(u) = {1, 3, 4, …, k−1, k+1}` (`verify §2`, wall solved directly in
   `E` for `k = 3..7`). This is one explicit exotic family in every band.
3. **The AP-collapse is a `k = 3` accident.** The band-3 tiling theorem `[CLO]`
   §2 (domino/triomino) collapsed the *entire* degree-3 exotic branch to a
   **single** one-parameter family `{r, r+2, r+4}`. This does **not** generalize:
   at minimal degree `deg a_k = k`, the realizable exotic tops number
   ```
        k = 2 : 0    (none — band 2 is clean, no exotic sector at all)
        k = 3 : 1    ({0,2,4}, the AP)
        k = 4 : 4    ({0,1,3,6}, {0,2,3,5}, {0,3,5,6}, {0,3,6,9})
        k = 5 : 13   (incl. the step-4 AP {0,4,8,12,16})
        k = 6 : 41
   ```
   (single-coset enumeration, min-root normalized to `0`, `verify §2` for `k ≤ 5`;
   `k = 6` off-verifier). So `k = 2`'s
   emptiness explains why band 2 needed no exotic argument; `k = 3`'s uniqueness is
   what made the one-parameter AP proof possible; and **for `k ≥ 4` the exotic base
   is a finite multi-family set, not one AP**. The proof technique must therefore
   run family-by-family (which it does — Route 4) rather than on a single AP.

*(Rigor note: the count is an exact single-coset scan whose value is
**bound-sensitive** — the max-root cutoff must exceed `~4k` before it stabilizes
(the naive `3k` cutoff misses `{0,4,8,12,16}` and undercounts `k = 5` as `12`);
`verify §2` scans to `5k`, checked stable off-verifier at `6k, 8k`. The counts
corroborate finiteness and the multi-family phenomenon, not a closed classification
of every exotic top at every degree; the multi-coset and higher-degree strata are
not enumerated here. The `k = 2` emptiness and `k = 3` uniqueness match the proved
band-2/3 results; the `k = 4, 5` counts `4, 13` independently agree with the
sibling band-4/band-5 experiments of this wave (`verify_band4_experiment.py`,
`verify_band5_comparison.py` / `band5-comparison.md`).)*

---

## 3. Route 3 — the tropical filler skeleton (Lemma-R blocks, PROVED any `k`)

At `Q₀ = (T−1)G` the free **negative** data enters linearly through exactly **two
fillers** — `a_{−(k−1)}` (in the level-`(k−1)` block) and `μ_k a_{−k}` (in the
level-`k` block) — because every other negative `a_{−1},…,a_{−(k−2)}` is already
constrained by the positive cascade `Q_{2k−2},…,Q_1` (it pairs with a nonzero
`b_{m+j}` at some `m ≥ 1`), while `a_{−(k−1)}` and `a_{−k}` first surface at `Q₀`.
This is the general-`k` version of `[CLO]` §4's "two fillers `a_{−2}`, `μ₃ a_{−3}`."

> **Lemma R (block leading coefficients, any `k`).** The two filler blocks cannot
> self-cancel their tops:
> ```
>   level-k   K_k[c] := Σ_{r=1}^{k}   a_k^[−r] c^[k−r] :   lead = k·lc(a_k)·lc(c),      deg = deg a_k + deg c;
>   level-(k−1) filler := −Σ_{r=1}^{k−1} b_{k−1}^[−r] a_{−(k−1)}^[k−1−r] :   lead = −(k−1)·lc(b_{k−1})·lc(a_{−(k−1)}).
> ```
> Both leading coefficients are nonzero in characteristic 0.

*Proof.* Each of the `k` (resp. `k−1`) shifted products has the same top-degree
term `lc(a_k)lc(c)E^{deg a_k + deg c}` (resp. with `b_{k−1}, a_{−(k−1)}`), all of one
sign, so they add rather than cancel; the count is the block length `k` (resp.
`k−1`). `verify §3` checks the degree and coefficient symbolically for
`k = 2,…,6`. ∎

Because neither filler block self-cancels, the highest coefficients of `Q₀ = 1`
(which no lower block reaches) force the **leading coefficient of each filler to
vanish**, then the next, cascading downward — the degree-free *skeleton* of the
kill (`[CLO]` §4, here general). What the skeleton alone does not finish is the
**low end**, where the positive cascade couples in; that is Lemma S.

---

## 4. Route 4 — the slope-forcing residual, and the honest reduction

Assemble Routes 1–3. Gauge `b_k = 0`; the wall fixes `b_{k−1} = u`; the positive
cascade `Q_{2k−2},…,Q_1` forward-solves `b_{k−2},…,b_{−(k−1)}` (subject to
solvability conditions on the positive middle `a_{k−1},…,a_0`); `b_{−k} = μ_k a_{−k}`.
Then `Q₀ − c` (slope `c` symbolic) is a polynomial in `E` whose coefficient system,
together with the two fillers, decides the unit. Routes 1–3 reduce the decision to:

> **Lemma S (low-end slope-forcing certificate).** On the positive-cascade variety,
> after the tropical skeleton annihilates the filler leading coefficients, the
> residual coupling of the slope `c` to the positive data forces **`c = 0`**.

**Lemma S is the one open lemma.** Granting it, the moment-unit theorem at band `k`
follows: `G = E` (slope `1`) is impossible, so `Q₀ = 1` is unsatisfiable and the
exotic branch is empty. Everything else in the chain (L1 reformulation `[WEAP]`
§W4q; L2 Lemma P §1; L3 wall classification §2; L4 Lemma R §3) is a proved
arbitrary-`k` theorem.

**Status of Lemma S — proved / verified (this memo, `verify §4–§5`):**

- **`k = 3`: PROVED.** The AP class `{r,r+2,r+4}`, `r` **symbolic**, `d = 1`, solves
  through `r`-independent pivots to the sole residual `c = 0` (`[CLO]` §5, Theorem 1);
  per-instance through `d ≤ 4`, `deg a₃ ≤ 9`. Reproduced here (`verify §4`).
- **`k = 4`: PROVED for a whole translation family; verified for all minimal
  families.** For the universal family `{r, r+2, r+3, r+5}` with `r` **symbolic**
  and `d = 1`, the positive cascade solves through **`r`-independent pivots** and
  the **sole residual is `c`** ⇒ slope forced to `0` for the entire family at once
  (`verify §5`, 10 pivots, all `r`-free). All **four** minimal exotic families are
  killed at `d = 1` (`{cascade} ∪ {Q₀=1}` Gröbner `= [1]`, `{Q₀=0}` feasible — the
  **unit** is the killer), together with integer translates `t ∈ {−1,0,1,2}` of each
  (16 instances), and the universal family at `d = 2` (`verify §4`).
- **`k = 5`: VERIFIED (prime cyclotomic `Φ₅`).** The universal exotic top
  `{0,2,3,4,6}` at `d = 1` gives `{cascade} ∪ {Q₀=1}` Gröbner `= [1]`, `{Q₀=0}`
  feasible (`verify §4`).
- **Positive control.** The genuine band-4 tame pair `U = x+∂`, `X = U⁴−∂`, `D = U`
  (`[D,X]=1`, `a₄ = 1` trivial cube) is **admitted** by the same pipeline with no
  spurious infeasibility (`verify §6`) — the kills above are real, not artifacts.

**Why this is the right evidence.** The band-3 closure `[CLO]` left exactly this
"uniform-in-free-degree" certificate as its residual (§8 there). The worry was that
its proof leaned on the AP collapse — a `k = 3` accident (§2 above). This memo
shows the concern is unfounded at the next two bands: the identical
moment-carries-no-unit kill holds at `k = 4` **past the composite-cyclotomic
threshold** (`S₄ = Φ₂Φ₄`, four exotic families, all killed) and at `k = 5` (prime
`Φ₅`), and the `r`-independent-pivot mechanism that powered the `k = 3` symbolic
proof **transfers verbatim** to the `k = 4` universal family. The mechanism is
uniform in `k`; only a **closed-form** certificate uniform in the free degree `d`
(and across the finitely many families per band) is missing.

**Why the pivots are `r`-independent (the mechanism, `verify §5`).** Translation
`E ↦ E − r` is an equivariance of the *positive* cascade (positive data is free of
membership), so the pivots — leading coefficients of the forward-solve steps — are
translation-invariant, hence `r`-free. Membership is pinned at `0` (not `r`), so it
breaks the equivariance only in the **bottom**; yet Lemma P routes the entire
bottom into the single protected products `a_{−i}(i)`, which the tropical skeleton
(Route 3) annihilates. The surviving slope residual is therefore a rational
identity in `r` that reduces to `c = 0` — the family is killed uniformly in `r`.
This is the general reason the `k = 3` AP proof and the `k = 4` universal-family
proof both terminate at "sole residual `= c`."

---

## 5. The reduction as a DC1 roadmap (the named finite list of lemmas)

> **Moment-unit at band `k` holds provided the following lemmas.**
>
> | # | Lemma | Status |
> |---|---|---|
> | L1 | `Q₀ = (T−1)G`, `G(0)=0`, so `Q₀=1 ⇔ G=E`, `slope = G(1)` | **PROVED any `k`** (`[WEAP]` §W4q; `verify §0–§1`) |
> | L2 | Moment-slope formula (Lemma P): `G(1) = Σ_i[a_i(0)b_{−i}(i) − a_{−i}(i)b_i(0)]`, level-`k` term `μ_k a_k(0)a_{−k}(k)` | **PROVED any `k`** (§1; `verify §1`) |
> | L3 | Wall necklace `(WALL-N)`: exotic = non-effective cofactor `g`; minimal degree `k`; finite multi-family for `k ≥ 4` | **PROVED any `k`** (structure); enumeration **computed** (§2; `verify §2`) |
> | L4 | Tropical filler skeleton (Lemma R): filler block leads `k·lc·lc`, `−(k−1)·lc·lc` ≠ 0 | **PROVED any `k`** (§3; `verify §3`) |
> | L5 | Low-end slope-forcing certificate (Lemma S): residual forces `c = 0` | **OPEN in general**; PROVED `k=3` + `k=4` universal family; VERIFIED `k=4` (all families, `d≤2`), `k=5` (§4; `verify §4–§5`) |
>
> With L1–L5 for a given `k`, the exotic sector at band `k` is empty. Combined with
> the shifted-`k`-th-power sector (a genuine cube gatekeeper, closed `k ≤ 3`) and
> the `a_k = 0` orientation reductions to band `< k`, this is the quantum band-`k`
> induction step; the band filtration then yields DC1.

**The single obstruction to DC1 (quantum face), precisely.** It is **L5, uniform in
the free degree `d`**. All the band-agnostic scaffolding (L1–L4) is proved; the
per-band exotic base is finite (L3); the kill is verified through band 5. What
resists is a **closed-form slope certificate** — a `d`-uniform reason the low-end
residual is a nonzero multiple of the slope `c`. The tropical skeleton (L4)
reduces the fillers but does not, in closed form, carry the positive-cascade
coupling at the bottom degree. This is the exact `d`-uniform frontier `[CLO]` §8
named at band 3, now shown to be the **only** general-`k` gap and shown to persist
verbatim across the prime/composite cyclotomic divide.

---

## 6. Rigor ledger

**PROVED (arbitrary `k`, structural, machine-checked identities `k = 2..5/7`):**
- L1 the reformulation `Q₀=(T−1)G`, `G(0)=0`, `slope = G(1) =` const. coeff. `Q₀`
  (`[WEAP]` §W4q re-derived; `verify §0–§1`);
- **L2 Lemma P**, the moment-slope formula and its membership collapse to one
  product per side per level, gauged form `μ_k a_k(0)a_{−k}(k)` (§1; `verify §1`) —
  the central new theorem;
- L3 the wall necklace reduction `(WALL-N)`, the `S_k = ∏_{d|k}Φ_d`
  prime/composite dichotomy, the universal cofactor `g = 1−S+S²` for all `k ≥ 3`,
  and minimal exotic degree `= k` (§2; `verify §2`);
- **L4 Lemma R**, the two filler-block leading coefficients `k·lc·lc`,
  `−(k−1)·lc·lc` (§3; `verify §3`).

**PROVED per family / instance (exact Gröbner / slope-residual decision):**
- **L5 (Lemma S)** at `k = 3` (AP, `r` symbolic, `d = 1`; reproduced) and at
  `k = 4` for the universal translation family (`r` symbolic, `d = 1`, sole residual
  `= c`, `r`-independent pivots) (§4; `verify §4–§5`).

**VERIFIED (exact, bounded — per-instance emptiness proofs, not a `d`-uniform proof):**
- `{cascade} ∪ {Q₀=1}` empty (Gröbner `= [1]`) while `{Q₀=0}` feasible, at:
  `k = 3` `{0,2,4}` `d=1`; `k = 4` all four minimal families + 12 integer translates
  `d=1`, universal `d=2`; `k = 5` universal `d=1` (§4; `verify §4`);
- the minimal-exotic multi-family counts `0,1,4,13,41` for `k = 2..6` (§2; `verify §2` for `k ≤ 5`);
- the band-4 tame positive control admitted with no false kill (§4; `verify §6`).

**CONJECTURE (the residual, precisely delimited):**
- **L5 uniform in the free degree `d`** (and across the finitely many exotic
  families per band). Verified `k ≤ 5` at low `d`; no closed-form `d`-uniform
  certificate is extracted. This is the sole general-`k` gap; it is the band-`k`
  incarnation of `[CLO]` §8 and `[WEAP]` §W5 (the "moment-unit-unrealizable"
  weapon / lattice modulus), here shown to be the only obstruction and to persist
  across the prime/composite cyclotomic divide.

**AUDIT-SENSITIVE POINTS (called out honestly):**
- The exotic-family enumeration (§2) is a **single-coset, bounded-root scan**;
  finiteness and the multi-family phenomenon are its content, not a full
  classification of exotic tops at all degrees/cosets.
- The Gröbner kills are **per-instance** (fixed integer top, fixed free-degree `d`).
  Each is an exact emptiness proof for that instance; none is an arbitrary-`d`
  proof. The `r`-symbolic results (`k = 3`, `k = 4` universal) are the only ones
  covering an infinite family, and only at `d = 1`.
- The positive control's forward solver emits one condition; it forces a
  **periodicity-kernel parameter to 0** (consistent with the tame `b_2 = 0`), and
  setting the `b_1`-kernel parameter to `1` recovers `D` — i.e. genuine feasibility
  is detected, so the infeasibilities are real (`verify §6`).

**NOT claimed:** the moment-unit theorem at all `k`; any band-`k` classification; the
shifted-power sector beyond `k ≤ 3`; a `d`-uniform Lemma S; JC2; DC1; any
counterexample (none is produced — consistent with DC1).

---

## 7. What genuinely differs at general `k` (vs the band-3 proof)

1. **The AP collapse is gone.** Band 3's single-family exotic base was a `k = 3`
   accident; `k ≥ 4` has a finite **multi-family** exotic base (4 at `k=4`, 13 at
   `k=5`). The proof runs family-by-family; the `r`-symbolic pivot method still
   closes each family (verified for the `k=4` universal family). *(New, sharp.)*
2. **Composite cyclotomic slack appears and is survived.** `S₄ = Φ₂Φ₄` opens a
   qualitatively new necklace slack absent at prime bands; the moment still carries
   no unit there (all four `k=4` families killed). *(New.)*
3. **Lemma P is the uniform protection theorem** the band-2/3 corpus was reaching
   for: the deep-level membership protection is exactly "one liberated value
   `a_{−i}(i)` per bottom coefficient," proved for all `k`. *(New, load-bearing.)*
4. **L1–L4 are band-agnostic theorems**; only L5's `d`-uniformity resists — the
   same frontier at every band, now isolated as the unique general-`k` gap.

---

## 8. Verification

```sh
uv run --with sympy python research/dc1-program/verify_moment_unit_k.py
```
Exact SymPy; 81 `PASS` lines (~30 s). §0 engine (`Q_m` = commutator; `Q₀=(T−1)G`
for `k=2..5`); §1 **Lemma P** (`G(0)=0`, `slope=G(1)`, the moment-slope formula,
raw and gauged, `k=2..5`); §2 the wall necklace (`S_k` prime/composite
factorization, universal cofactor `k=3..7`, minimal-exotic counts `k=2..5`); §3
**Lemma R** filler-block leading coefficients (`k=2..6`); §4 the **moment-unit
kill** (`k=3` reproduced; `k=4` all four families + translates + `d=2`; `k=5`
prime-cyclotomic) with `{Q₀=0}` feasible; §5 the **`r`-symbolic slope certificate**
for the `k=4` universal family; §6 the band-4 positive control. A successful run
ends `ALL MOMENT UNIT K CHECKS PASSED`.
