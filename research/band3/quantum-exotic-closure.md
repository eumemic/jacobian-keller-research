# Closing the quantum exotic `b₂ ≠ 0` sub-case: the moment carries no unit

**INDEPENDENTLY DERIVED — MACHINE-CHECKED IDENTITIES — NOT PEER REVIEWED — BAND-SCOPED**

This memo attacks the one residual gap left by the quantum exotic branch memo
(`quantum-exotic-branch.md`, commit `ebfc64d`, §7 "CONJECTURED": *the
arbitrary-degree closure of the `b₂ ≠ 0` half*). That gap is the load-bearing
input for the quantum band-3 floor: `astar-band3.md` (`ebfc64d`, §6) routes
quantum **A\*-I** into exactly this question — *"is a wall-admitting
non-shifted-cube top killed downstream?"* — so closing it firms the quantum
band-3 induction rung on the Dixmier (DC1) face.

Conventions are frozen, identical to every sibling quantum memo: `A₁[x⁻¹] =
⊕_k x^k ℂ[E]`, `(x^a f)(x^b g) = x^{a+b} f^[b] g`, `f^[r](E) = f(E+r)`,
`E = x∂`, ladder-`m` coefficient `Q_m = Σ_{k+l=m}(b_l^[k]a_k − a_k^[l]b_l)`,
`[D,X]=1 ⇔ Q_m = δ_{m0}` (`m ∈ [−6,6]`), membership
`E(E−1)⋯(E−r+1) | a_{−r}, b_{−r}`, gauge `b₃ = 0`. Citations: the cascade
`quantum-band3-cascade.md` = [CASC] (current HEAD `b9f9cf3`; the Telescoping
Lemma, the `Q₅` wall, and the bottom proportionality are unchanged in form from
`99fe6ee`); the Wave-B exotic memo `quantum-exotic-branch.md` = [EXO]
(`ebfc64d`, unchanged); `band-k-weapons.md` = [WEAP]; `astar-band3.md` = [AST].
Everything displayed is checked exactly by
[`verify_quantum_exotic_closure.py`](verify_quantum_exotic_closure.py) (66
`PASS` lines, ~15 s; ends `ALL QUANTUM EXOTIC CLOSURE CHECKS PASSED`).

## 0. Verdict

> **The exotic `b₂ ≠ 0` kill is the statement that *the moment carries no
> unit*.** After the reduction `Q₀ = 1 ⇔ G = E` ([CASC] §4), the unit is exactly
> the **slope** of the central potential `G`, i.e. `G(1) =` the constant
> coefficient of `Q₀`. On the positive-cascade variety this slope is **forced to
> `0`** — proved for the whole degree-3 (= AP) branch at `d = 1` and verified at
> higher degree — so `G = E` (slope `1`) is impossible. This is the quantum
> incarnation of the classical Theorem-A step-6 principle ([EXO] §3.2: "the moment
> can carry only `τ`/`E`, never a residual unit against a membership-protected
> extreme") and of the [WEAP] §W5 "moment-unit-unrealizable" weapon.
>
> Two exact structural facts sharpen the branch and make the kill degree-free at
> the base degree:
>
> - **(Structural.)** At `deg a₃ = 3` the realizable exotic tops are **exactly**
>   the step-2 arithmetic-progression class `{r, r+2, r+4}` — every other
>   `Φ₃`-divisible non-cube 3-root multiset fails `b₂`-effectivity, i.e. is *not*
>   in the `b₂ ≠ 0` branch at all. **So the whole degree-3 exotic branch IS the
>   AP class.** (`verify §1`.)
> - **(Main, PROVED.)** For the AP top `{r, r+2, r+4}` with `r` **symbolic** and
>   `d = 1` free data, the moment slope is forced to `0` — the branch is **empty**
>   for generic `r` at once (`verify §4`), with the small integer values
>   `r ∈ {0,−4,1,−1,2,3}` confirmed directly (`verify §5`). This closes the `d = 1`
>   slice of the entire degree-3 exotic branch.

No DC1/JC2 counterexample is produced; the exotic branch is empty everywhere the
argument reaches, consistent with DC1. The **residual** is honestly stated in §8:
the *uniform-in-free-degree* certificate (`d ≥ 2` is proved per instance but not
by one closed form), and the non-AP realizable exotic tops that appear only at
`deg a₃ ≥ 6`.

## 1. The moment-unit reformulation (degree-free)

By the Telescoping Lemma ([CASC] §4, re-derived in `verify §0` as an identity for
generic band-3 coefficients),
```
   Q₀ = (T−1)G ,   G = Σ_{k=1}^{3} Σ_{j=0}^{k−1} ( a_k^[j−k] b_{−k}^[j] − b_k^[j−k] a_{−k}^[j] ) ,
```
and under genuine membership `G(0) = 0` **identically** (`verify §0`: every term
of `G(0)` carries a value `a_{−k}(j)` or `b_{−k}(j)` with `0 ≤ j ≤ k−1`, killed
by `E^{underline k} | a_{−k}, b_{−k}`). A polynomial with `(T−1)G` constant is
affine, so
```
   Q₀ = 1  ⇔  G = E .                                          (UNIT)
```
The single number that decides (UNIT) is the **slope** of `G`. Because `G(0)=0`,
```
   slope := constant coefficient of Q₀  =  G(1) − G(0)  =  G(1)   (verify §0).
```
So the exotic kill is: *can the constrained bilinear form `G` have `G(1) = 1`?*
Writing the kill this way — **"the moment carries no unit" ⇔ `G(1)` is forced to
`0`** — is the exact quantum mirror of the classical statement that the band-3
moment `M = τ` carries only order-`1` content and never a residual unit ([EXO] §4
table, row "classical non-cube"; [CASC]/`classical-band3-cascade.md` §6).

In gauge `b₃ = 0` the potential splits by level, `G = G₁ + G₂ + G₃` with
```
   G₁ = a₁^[−1] b_{−1} − b₁^[−1] a_{−1},
   G₂ = a₂^[−2] b_{−2} + a₂^[−1] b_{−2}^[1] − b₂^[−2] a_{−2} − b₂^[−1] a_{−2}^[1],
   G₃ = a₃^[−3] b_{−3} + a₃^[−2] b_{−3}^[1] + a₃^[−1] b_{−3}^[2]     (b₃ = 0),
```
and `b_{−3} = μ₃ a_{−3}` from `Q_{−6}` ([CASC] §5.1). The positive cascade
`Q₄,Q₃,Q₂,Q₁ = 0` forward-solves `b₁, b₀, b_{−1}, b_{−2}` from `(a₃,b₂,a₂,a₁,a₀)`
subject to solvability conditions on `(a₂,a₁,a₀)` (and `a_{−1}` enters at `Q₁`
through the `(−1,2)` pair `b₂^[−1] a_{−1}`, [EXO] §3.1). The free negative data
`a_{−1}, a_{−2}, a_{−3}, μ₃` then enters `Q₀` **only** — the exotic kill uses the
positive cascade `Q₄..Q₁` and the moment `Q₀`, nothing lower (this is the strong
form: fewer equations, [EXO] §3).

## 2. Structural: the degree-3 exotic branch is the step-2 AP class

The `Q₅` Wall Lemma ([CASC] §3.1) admits a nonzero `b₂` iff, coset by coset, the
root multiset `A` of `a₃` is `Φ₃(S)`-divisible **and** the forced quotient
`B = S(1+S)A/Φ₃` is effective (a genuine multiset). Non-shifted-cube = the
cofactor `A/Φ₃` is *not* effective.

> **Fact 1.** The realizable tops at `deg a₃ = 3` are **exactly** `{0,1,2}` (the
> consecutive cube) and `{0,2,4}` (the exotic counterexample). Hence **the
> degree-3 exotic `b₂ ≠ 0` branch is precisely the step-2 AP class
> `{r, r+2, r+4}`** (`a₃ = (E−r)(E−r−2)(E−r−4)`, `b₂ = (E−r−1)(E−r−4)`, wall-solved
> for all `r`, `verify §4/§5`).

*Proof.* A `Φ₃`-divisible degree-3 top has **distinct** roots (a repeated root
gives `A = 2S^i + S^j`, and `A(ω) = 2ω^i + ω^j ≠ 0` for a primitive cube root
`ω` since `|2ω^i| = 2 ≠ 1 = |ω^j|`), so `A = 1 + S^a + S^b`, `0 < a < b`, after
translation. By `WALL-DEG` ([CASC] §3.3) `deg b₂ = 2`, so `B` has exactly two
roots and `B·Φ₃ = S(1+S)A`, i.e.
```
   (S^p + S^q)(1+S+S²) = (S+S²)(1+S^a+S^b)
                       = {1,2} ⊎ {a+1,a+2} ⊎ {b+1,b+2}   (exponent multiset).
```
The right side is **three consecutive-pair "dominoes"**; the left is **two
consecutive-triple "triominoes"**. A single triomino (span 3) can meet two
dominoes only if their left ends differ by `≤ 2`; since the dominoes sit at
`1 < a+1 < b+1`, tiling forces `a+1 ≤ 3` **and** `b+1 ≤ (a+1)+2`, i.e.
`a ≤ 2, b ≤ a+2`. The finitely many candidates `(a,b) ∈
{(1,2),(1,3),(2,3),(2,4)}` reduce by `Φ₃`-divisibility (distinct residues
`mod 3`): `(1,3)` gives residues `{0,1,0}`, `(2,3)` gives `{0,2,0}` — both
rejected; `(1,2)` is the cube, `(2,4)` the exotic. ∎ (`verify §1` corroborates by
exhaustive scan `b ≤ 12`; extended to `b ≤ 40` off-verifier — no further
solution, confirming the bound.)

Every other `Φ₃`-divisible non-cube 3-root multiset (e.g. `{0,1,5}`, `{0,4,8}`,
`{1,5,6}`) has a **non-effective** `B` — no nonzero `b₂` — so it is not in the
`b₂ ≠ 0` branch at all.

This is genuinely new: it collapses "the exotic branch" at the base degree to a
one-parameter family, the exact class the wall counterexample `E(E−2)(E−4)`
([CASC] §3.2) generates. **Caveat (honest, `verify §1`):** the collapse is
degree-3-specific. At `deg a₃ = 6` there are **22** realizable non-AP exotic tops
(e.g. `{0,1,2,3,5,7}`, `{0,1,2,4,6,8}`), so for `deg a₃ ≥ 6` the AP class is a
proper subclass — see §8.

## 3. Membership-protection of the bottom (degree-free)

The level-3 block `G₃ = Σ_{j=0}^{2} a₃^[j−3] b_{−3}^[j]` with `b_{−3} = μ₃ a_{−3}`
is membership-protected in its contribution to the slope (`verify §2`):
```
   G₃(0) = 0        (a_{−3}(0)=a_{−3}(1)=a_{−3}(2)=0),
   G₃(1) = μ₃ · a₃(0) · a_{−3}(3)      (the only surviving bottom term).
```
So the **entire negative bottom** `(μ₃, a_{−3})` enters the slope `G(1)` through
the single product `μ₃ a₃(0) a_{−3}(3)`. For any AP top with `0` a root
(`r ∈ {0,−2,−4}`, in particular the Wave-A witnesses `W1 = {0,2,4}`,
`W2 = {0,−2,−4}`) this is `0`: **the bottom is absent from the unit equation
entirely**, and the kill lives in the positive + middle data. For `a₃(0) ≠ 0`
(e.g. `r = 1,−1,2,3`) the bottom enters — and still cannot rescue the unit
(§6). This is the exact structural fact [EXO] §3.2 flagged, here isolated as a
one-line membership identity.

## 4. The tropical filler-annihilation skeleton (degree-free)

The free negatives enter `Q₀ = (T−1)G` linearly as two **fillers**: `a_{−2}` (in
`G₂`) and `μ₃ a_{−3}` (in `G₃`; `a_{−3}` occurs only through the product
`μ₃ a_{−3}`, so the relaxation `μ₃ a_{−3} ↦` free is exact, [EXO] §3.2). Their
leading coefficients are pinned by the Lemma-R staggered identities ([CASC] §5.3),
verified symbolically in degree (`verify §3`):
```
   G₃-block K₃[c] = a₃^[−3]c + a₃^[−2]c^[1] + a₃^[−1]c^[2] :
        deg = deg a₃ + deg c,   lead = 3·lc(a₃)·lc(c) ≠ 0 ;
   G₂-filler  −b₂^[−2]a_{−2} − b₂^[−1]a_{−2}^[1] :
        deg = deg b₂ + deg a_{−2},  lead = −2·lc(b₂)·lc(a_{−2}) ≠ 0 .
```
Because `K₃` cannot self-cancel its top (all three shifts share the sign of
`lc(a₃)lc(c)`), the highest `Q₀`-coefficients — which no lower block reaches —
force the **top coefficient of `a_{−3}` to vanish**, then the next, cascading
downward; the `G₂`-filler does the same for `a_{−2}`. This is the degree-free
*skeleton* of the kill: the top of `Q₀ = 1` triangularly annihilates the fillers
from above. (It is the same "termwise degree domination" the classical sibling
used, `classical-cube-closure.md` §7 at `99fe6ee`, here driving fillers to `0`
rather than clashing two extremes.) What the skeleton alone does **not** finish
is the low end, where the positive data couples in; that is completed next.

## 5. Main theorem: the AP class carries no unit (r symbolic, d = 1)

Run the positive cascade for the AP top with `r` **symbolic** and `d = 1` free
data (all of `a₂,a₁,a₀` degree ≤ 1, negatives degree ≤ 1 modulo membership),
relax `μ₃ a_{−3} ↦` free, and eliminate the fillers `(a_{−2}, μ₃ a_{−3})` and
the cascade solvability from the coefficient system of `Q₀ = c` (slope `c`
symbolic). The elimination leaves (`verify §4`):

> **Theorem 1 (moment carries no unit, AP class, `d = 1`).** With `r` symbolic,
> the positive cascade solves through **`r`-independent** pivots (`verify §4`
> asserts every pivot is free of `r`; the `W1` walk-through §4b exhibits them), and
> after the filler elimination the **sole residual condition is `c = 0`** — a
> rational identity in `r`. Hence the moment slope
> `G(1)` is forced to `0`, `G = E` is impossible, and the **degree-3 exotic
> `b₂ ≠ 0` branch has no `d = 1` member for generic `r`.** The finitely many
> integer values where the elimination could degenerate are covered directly: `r ∈
> {0,−4,1,−1,2,3}` each give `{cascade} ∪ {Q₀=1}` Gröbner-`= [1]` (`verify §5`).
> [PROVED, symbolic `r` + explicit small values.]

The mechanism is transparent for `W1` (`r = 0`), spelled out in `verify §4b`.
The positive cascade solves to `a₂ = a₂(0)` const, `a₀` const, `a_{−1} = 0`,
`a₁ = a₂(0)² − (a₂(0)²/3)E`; then `Q₀`, after these, depends only on the four
fillers `u₀,u₁ = ` coefficients of `a_{−2}` and `p₀,p₁ = μ₃·(`coefficients of
`a_{−3})`. The top coefficients form a triangular chain
```
   [E⁶] = 21 p₁              ⇒ p₁ = 0,
   [E⁵] = 9(2p₀ − 15p₁)      ⇒ p₀ = 0,
   [E⁴] = −5(2u₁ + 27p₀ − 57p₁) ⇒ u₁ = 0,
   [E³] = −8u₀ + 40u₁ + 336p₀ − 441p₁ ⇒ u₀ = 0,
```
annihilating every filler, after which the constant coefficient (the slope) is
```
   [E⁰] = −8u₀ − 16u₁ = 0 ,
```
**stranded at `0`** — the unit `1` has nothing to stand on. This is precisely the
skeleton of §4 running all the way to the bottom, with the membership-freed slope
left as `0`.

## 6. Verification across the class and beyond `d = 1`

The kill is not `W1`-specific and extends past the `d = 1` proof (each an exact
Gröbner decision `= [1]`, hence a per-instance emptiness **proof**, not mere
evidence):

- **Whole AP class, `d = 1` (`verify §5`).** For `r ∈ {0, −4, 1, −1, 2, 3}`
  (both `a₃(0) = 0` and `a₃(0) ≠ 0` cases): `{cascade} ∪ {Q₀ = 1}` is
  **infeasible** (`Gröbner = [1]`) while `{cascade} ∪ {Q₀ = 0}` is **feasible** —
  the **unit** is exactly the killer.
- **`W1`, `W2` at `d = 2` (`verify §5`, Gröbner `= [1]`); `W1` at `d = 3, 4`
  (`verify §5`, slope residual a nonzero multiple of `c` ⇒ slope forced to `0`).**
  Notably the controlling combination `w` survives **only at `d = 2`**: at
  `d = 1, 3, 4` the residual collapses directly to "slope `= 0`", so `d = 2` is the
  lone free-degree needing the two-rung `w`-certificate below.
- **The exact `d = 2` certificate (`verify §5`).** Eliminating the fillers leaves
  two residuals on the positive data: a homogeneous one that is a nonzero multiple
  of the **controlling combination**
  ```
     w := a₁₂² (a₂₀ − 4 a₂₂)          (a₁₂ = lc(a₁), a₂₀ = a₂(0), a₂₂ = lc(a₂)),
  ```
  forcing `w = 0`, and one carrying the slope `c` that couples `c` to `w`; on
  `{w = 0}` it forces `c = 0`. (This is [EXO] §3.2's `{8w = 0, 7w = 9}` in
  slope-normalized form: the homogeneous rung forces `w = 0`, so the slope cannot
  be `1`.)
- **Arbitrary top-degree AP: `deg a₃ = 3, 6, 9` (`verify §5b`).** The step-2 AP
  top `{0,2,…,2(3g−1)}` (`3g` roots, wall-admissible iff `3 | (# roots)`, [CASC]
  §3.1) at `d = 1` gives `{cascade} ∪ {Q₀ = 1}` infeasible for `g = 1,2,3`. The
  mechanism scales in the top degree.
- **Positive control (`verify §6`).** On the genuine `b₂ = 0` tame pair
  `U = x+∂, X = U³−∂, D = U` (`[D,X]=1`), the same solver reproduces `D`
  (`b₁ = 1, b_{−1} = E`) with **no spurious conditions** — the infeasibility above
  is a real kill, not a pipeline artifact.

## 7. Why the moment, not the wall — and the classical mirror

The wall `Q₅` is too weak to force a shifted cube (that is the whole necklace gap,
[CASC] §3.2); the operative equation is the moment `Q₀`. The three incarnations of
the band-3 gatekeeper now read (extending [EXO] §4):

| face / sector | mechanism at `Q₀` / `C₀` |
|---|---|
| classical non-cube ([CASC-cl] §6) | `M = τ` vs `τ³ \| b_{−3}`: unit forces order 1, membership needs ≥ 3 |
| quantum exotic `b₂ = 0` ([EXO] §2) | `L₀[b_{−3}] = 1` impossible: `deg ≥ 2` or `0`, never a nonzero const |
| **quantum exotic `b₂ ≠ 0` (this memo)** | **moment slope `G(1)` forced to `0`: the moment carries no unit** |

All three are the same principle — *a membership-protected extreme lets the
moment carry only its natural low-order content (`τ`, `E`), never a residual
unit*. The `b₂ ≠ 0` case is the one with no positive collapse, so the unit's
absence must be read directly off the constant coefficient of `(T−1)G`, which is
what Theorem 1 does. The [WEAP] §W5 "moment-unit-unrealizable" weapon is thereby
instantiated for the band-3 exotic residual: the `+1` there is this slope.

## 8. Claim disposition (proved / verified / residual)

**PROVED (arbitrary parameter, machine-checked):**
- the moment-unit reformulation `Q₀ = 1 ⇔ G = E`, `slope = G(1) = ` const. coeff
  of `Q₀`, and `G(0) = 0` under membership (§1, `verify §0`) — degree-free
  identities;
- **Fact 1:** the degree-3 realizable exotic branch is exactly the step-2 AP class
  `{r,r+2,r+4}` (§2, `verify §1`) — the exotic base collapses to one parameter;
- the bottom membership-protection `G₃(1) = μ₃ a₃(0) a_{−3}(3)` (§3, `verify §2`);
- the Lemma-R block leading-coefficient identities driving the filler-annihilation
  cascade (§4, `verify §3`), symbolic in degree;
- **Theorem 1:** for the AP class with `r` **symbolic** and `d = 1` free data
  (`r`-independent cascade pivots; sole residual `c = 0`), the moment slope is
  forced to `0` — the whole degree-3 exotic branch's `d = 1` slice is **empty** for
  generic `r`, with `r ∈ {0,−4,1,−1,2,3}` confirmed directly (§5, `verify §4`+§5).

**PROVED per instance (exact Gröbner / slope-residual decision):**
- `{cascade} ∪ {Q₀ = 1}` empty for `W1, W2` at `d = 2` (Gröbner `= [1]`), `W1` at
  `d = 3, 4` (slope forced to `0`), the AP class `r ∈ {0,±1,±2,3}` at `d = 1`, and
  the AP tops of `deg a₃ = 3, 6, 9` at `d = 1`; `{cascade} ∪ {Q₀ = 0}` feasible —
  the unit is the killer (§6);
- the exact `d = 2` slope certificate via `w = a₁₂²(a₂₀ − 4a₂₂)`, and the fact
  that `w` survives **only** at `d = 2` (`d = 1,3,4` collapse directly) (§6).

**RESIDUAL (honestly open, precisely delimited):**
- the **uniform-in-free-degree** closure. `d = 1` is proved (Theorem 1: generic
  `r` + the small integer values); `d = 2, 3, 4` are proved per instance; but the
  *controlling combination* `w` is degree-dependent (`d = 1, 3, 4` collapse
  directly — sole residual a multiple of the slope `c` — while `d = 2` carries a
  genuine `w`), and no single closed-form certificate uniform in `d` is extracted
  here. This is the same status frontier as the classical sibling's
  leading-coefficient loci (`classical-cube-closure.md` §7, `99fe6ee`) — the
  tropical skeleton (§4) reduces the problem but does not carry the uniform-`d`
  proof in closed form.
- the **non-AP realizable exotic tops at `deg a₃ ≥ 6`** (22 exist at degree 6,
  `verify §1`). Fact 1 makes the branch = AP only at degree 3; the higher-degree
  non-AP exotic tops are not treated here.

**NOT claimed:** a full band-3 theorem; any DC1/JC2 statement (no counterexample
is produced, consistent with DC1, [EXO] §5); closure of the `b₂ ≠ 0` half uniform
in free degree; the non-AP higher-degree exotic sector.

## 9. Relation to the induction and to astar-band3

[AST] §6 reduced **quantum A\*-I** to "is the reflected non-shifted-cube-necklace
top killed downstream?" — the open question this memo attacks. Theorem 1 answers
it in the affirmative at the base degree for the whole AP family (`d = 1`, generic
`r` + the small integer values), and the per-instance proofs extend it through
`d = 4` and `deg a₃ = 9`. The
band-3 floor on the DC1 face therefore stands, at the exotic gatekeeper, on: the
`Q₅` wall structure (necklace), the central integral `Q₀`, and now the
*moment-carries-no-unit* closure of the exotic residual — with the honest residual
being the uniform free-degree statement and the higher-degree non-AP tops.

## 10. Verification

```sh
uv run --with sympy python research/band3/verify_quantum_exotic_closure.py
```
Exact SymPy; §0 engine (`Q_m` = commutator; `Q₀ = (T−1)G`; `G(0)=0`;
`slope = G(1)`), §1 the degree-3 exotic = AP structural scan (+ the degree-6
non-AP residual), §2 the bottom membership-protection, §3 the Lemma-R block
leading coefficients, §4 the main theorem (`r` symbolic, `d = 1`) + the explicit
`W1` triangular annihilation, §5 the class/degree verification and the exact
`d = 2` certificate, §6 the positive control. A successful run prints 66 `PASS`
lines and ends `ALL QUANTUM EXOTIC CLOSURE CHECKS PASSED`.
