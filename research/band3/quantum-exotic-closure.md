# Moment-unit obstructions and certificates in proved quantum-exotic slices

**INDEPENDENTLY DERIVED — MACHINE-CHECKED IDENTITIES — NOT PEER REVIEWED — BAND-SCOPED**

This memo advances the residual gap left by the quantum exotic branch memo
(`quantum-exotic-branch.md`, commit `ebfc64d`, §7 "CONJECTURED": *the
arbitrary-degree closure of the `b₂ ≠ 0` half*). It proves a generic-`r` statement
for the degree-3 AP family at free-data degree `d = 1` and selected exact instances,
but does not close exceptional `r` loci, the gap uniformly in free degree, or
higher-degree non-AP tops.
`astar-band3.md` (`ebfc64d`, §6) routes quantum **A\*-I** through this question;
the advances here are therefore bounded/sliced evidence for that induction rung,
not a completed quantum band-3 floor or a DC1 result.

Conventions are frozen, identical to every sibling quantum memo: `A₁[x⁻¹] =
⊕_k x^k ℂ[E]`, `(x^a f)(x^b g) = x^{a+b} f^[b] g`, `f^[r](E) = f(E+r)`,
`E = x∂`, ladder-`m` coefficient `Q_m = Σ_{k+l=m}(b_l^[k]a_k − a_k^[l]b_l)`,
`[D,X]=1 ⇔ Q_m = δ_{m0}` (`m ∈ [−6,6]`), membership
`E(E−1)⋯(E−r+1) | a_{−r}, b_{−r}`, gauge `b₃ = 0`. Citations: the cascade
`quantum-band3-cascade.md` = [CASC] (current HEAD `b9f9cf3`; the Telescoping
Lemma, the `Q₅` wall, and the bottom proportionality are unchanged in form from
`99fe6ee`); the Wave-B exotic memo `quantum-exotic-branch.md` = [EXO]
(`ebfc64d`, unchanged); `band-k-weapons.md` = [WEAP]; `astar-band3.md` = [AST].
The verifier exactly checks identities and finite computations specifically
attributed to it. The displayed arbitrary-degree leading formulas and Fact 1
additionally rely on the written proofs below. The verifier is
[`verify_quantum_exotic_closure.py`](verify_quantum_exotic_closure.py) (53 `PASS`
lines; ends `ALL QUANTUM EXOTIC CLOSURE CHECKS PASSED`).

## 0. Verdict

> **In the proved slices, the exotic `b₂ ≠ 0` obstruction is that *the moment
> cannot carry the required unit*.** After the reduction `Q₀ = 1 ⇔ G = E`
> ([CASC] §4), the unit occupies the **slope** slot of the central potential `G`,
> i.e. `G(1) =` the constant coefficient of `Q₀`. The verifier symbolically forces
> this slope to `0` for the degree-3 AP family at `d = 1`, and supplies slope
> certificates for W1 at `d = 2,3,4`; other listed exact instances establish only
> `Q₀ = 1` infeasibility unless stated otherwise. Thus `G = E` (slope `1`) is
> impossible exactly in the documented slices. This is the quantum
> incarnation of the classical Theorem-A step-6 principle ([EXO] §3.2: "the moment
> can carry only `τ`/`E`, never a residual unit against a membership-protected
> extreme") and of the [WEAP] §W5 "moment-unit-unrealizable" weapon.
>
> Two structural facts sharpen the degree-3, `d = 1` slice with symbolic AP
> translation parameter `r`, subject, for the moment-obstruction statement, to
> genericity conditions:
>
> - **(Structural.)** At `deg a₃ = 3` the realizable exotic tops are **exactly**
>   the step-2 arithmetic-progression class `{r, r+2, r+4}` — every other
>   `Φ₃`-divisible non-cube 3-root multiset fails `b₂`-effectivity, i.e. is *not*
>   in the `b₂ ≠ 0` branch at all. **So the whole degree-3 exotic branch IS the
>   AP class.** (`verify §1`.)
> - **(Main, PROVED generically.)** For the AP top `{r, r+2, r+4}` with symbolic
>   `r` and `d = 1` free data, rational elimination forces the moment slope to `0`
>   for generic `r` (`verify §4`). The six values
>   `r ∈ {0,−4,1,−1,2,3}` are separately checked (`verify §5`), but they do not
>   exhaust all exceptional rank or denominator loci. Unchecked exceptional `r`
>   values remain open; no entire-slice or entire-branch closure is claimed.

No DC1/JC2 counterexample is produced; the exotic branch is empty everywhere the
argument reaches, consistent with DC1. The **residual** is honestly stated in §8:
the *uniform-in-free-degree* certificate (selected instances at `d = 2,3,4` are
proved), and all higher-degree non-AP realizable exotic tops.

## 1. The moment-unit reformulation (degree-free)

By the Telescoping Lemma ([CASC] §4, checked in `verify §0` from the stipulated
`Q_m` convention as an identity for generic band-3 coefficients),
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
degree-3-specific. In the verifier's bounded normalized-root window — integer
roots, `0` included, the remaining roots in `1..14` (equivalently max root
`≤ 14`) — `deg a₃ = 6` has **exactly 22** realizable non-AP exotic tops (e.g.
`{0,1,2,3,5,7}`, `{0,1,2,4,6,8}`). Thus the AP class is already a proper
subclass at degree 6; the residual is all higher-degree non-AP realizable tops,
not merely these 22 bounded-search examples — see §8.

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

## 4. Isolated filler blocks and the cross-cancellation gap (degree-free)

The free negatives enter `Q₀ = (T−1)G` linearly as two **fillers**: `a_{−2}` (in
`G₂`) and `μ₃ a_{−3}` (in `G₃`; `a_{−3}` occurs only through the product
`μ₃ a_{−3}`, so the relaxation `μ₃ a_{−3} ↦` free is exact, [EXO] §3.2). Their
leading coefficients follow from the displayed leading-term expansion and the
Lemma-R staggered identities ([CASC] §5.3); `verify §3` corroborates them only at
the finite degree pairs listed there:
```
   G₃-block K₃[c] = a₃^[−3]c + a₃^[−2]c^[1] + a₃^[−1]c^[2] :
        deg = deg a₃ + deg c,   lead = 3·lc(a₃)·lc(c) ≠ 0 ;
   G₂-filler  −b₂^[−2]a_{−2} − b₂^[−1]a_{−2}^[1] :
        deg = deg b₂ + deg a_{−2},  lead = −2·lc(b₂)·lc(a_{−2}) ≠ 0 .
```
These formulas prove that neither isolated filler block self-cancels. They do
not triangularly annihilate either filler in the combined equation: the two
blocks can have equal output degree and cross-cancel, and solved cascade blocks
may also contribute at relevant degrees. In fact, with unrestricted filler
degree, the two admissible filler images have an infinite-dimensional
intersection; see
[`two-filler-cross-cancellation.md`](../dc1-program/two-filler-cross-cancellation.md).
Filler elimination is therefore valid only in the explicitly computed
fixed-degree systems below, not as a degree-free theorem.

## 5. Main theorem: the AP class carries no unit (r symbolic, d = 1)

Run the positive cascade for the AP top with `r` **symbolic** and `d = 1` free
data (all of `a₂,a₁,a₀` degree ≤ 1, negatives degree ≤ 1 modulo membership),
relax `μ₃ a_{−3} ↦` free, and eliminate the fillers `(a_{−2}, μ₃ a_{−3})` and
the cascade solvability from the coefficient system of `Q₀ = c` (slope `c`
symbolic). The elimination leaves (`verify §4`):

> **Theorem 1 (generic moment obstruction, AP class, `d = 1`).** With `r` symbolic,
> rational elimination through the displayed pivots leaves the residual condition
> `c = 0` wherever those pivots and denominators are valid. Hence the moment slope
> `G(1)` is forced to `0`, and `G = E` is impossible for generic `r`. The values
> `r ∈ {0,−4,1,−1,2,3}` are checked individually by
> `{cascade} ∪ {Q₀=1}` Gröbner-`= [1]` (`verify §5`). Those six checks do not
> enumerate all exceptional rank or denominator values; all other exceptional `r`
> loci remain open. [PROVED for generic `r`; six exact specializations checked.]

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
**stranded at `0`** — the unit `1` has nothing to stand on. This is an exact
instance-specific annihilation chain in the fixed `W1`, `d=1` system; §4 explains
why it does not extend to an unrestricted degree-free filler theorem.

## 6. Verification across the class and beyond `d = 1`

The obstruction is not `W1`-specific and extends past the `d = 1` proof. Each
listed result is an exact per-instance Gröbner or slope-residual decision, hence a
proof only within the stated finite ansatz, not mere evidence:

- **Six checked AP specializations at `d = 1` (`verify §5`).** For
  `r ∈ {0, −4, 1, −1, 2, 3}` (both `a₃(0) = 0` and `a₃(0) ≠ 0` cases):
  `{cascade} ∪ {Q₀ = 1}` is **infeasible** (`Gröbner = [1]`) while
  `{cascade} ∪ {Q₀ = 0}` is **feasible** — the **unit** is exactly the obstruction.
- **Exact degree-3 `d = 2` infeasibility.** This verifier proves `Q₀ = 1`
  infeasible for W1 and W2 (`verify §5`, Gröbner `= [1]`); the backing Wave-B
  verifier `verify_quantum_exotic.py` §4 additionally proves the AP cases
  `r = 1, −1, 3` at `d = 2`. Only W1 has the slope certificate below. For W1 at
  `d = 3, 4`, this verifier finds a slope residual that is a nonzero multiple of
  `c`, hence forces slope `0`.
- **The exact W1 `d = 2` certificate (`verify §5`).** Eliminating the fillers leaves
  two residuals on the positive data: a homogeneous one that is a nonzero multiple
  of the **controlling combination**
  ```
     w := a₁₂² (a₂₀ − 4 a₂₂)          (a₁₂ = lc(a₁), a₂₀ = a₂(0), a₂₂ = lc(a₂)),
  ```
  forcing `w = 0`, and one carrying the slope `c` that couples `c` to `w`; on
  `{w = 0}` it forces `c = 0`. (This is [EXO] §3.2's `{8w = 0, 7w = 9}` in
  slope-normalized form: the homogeneous rung forces `w = 0`, so the slope cannot
  be `1`.)
- **Selected AP top degrees: `deg a₃ = 3, 6, 9` (`verify §5b`).** The step-2 AP
  top `{0,2,…,2(3g−1)}` (`3g` roots, wall-admissible iff `3 | (# roots)`, [CASC]
  §3.1) at `d = 1` gives `{cascade} ∪ {Q₀ = 1}` infeasible for `g = 1,2,3`.
  These checks establish neither symbolic slope `0` nor `Q₀ = 0` feasibility for
  the degree-6 and degree-9 cases.
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
| **quantum exotic `b₂ ≠ 0` (proved slices here)** | **slope forced to `0` in the certified W1/AP `d=1` slices; `Q₀=1` infeasible in the other listed exact slices** |

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
- the Lemma-R degree and leading-coefficient identities for each isolated filler
  block (§4), derived from the displayed arbitrary-degree leading-term expansion
  and corroborated at selected degree pairs by `verify §3`; these identities do
  not imply unrestricted filler annihilation;
- **Theorem 1:** for the AP class with `d = 1` free data, rational elimination
  forces the moment slope to `0` for generic `r`; the six values
  `r ∈ {0,−4,1,−1,2,3}` are confirmed directly, while all other exceptional rank
  or denominator loci remain open (§5, `verify §4`+§5).

**PROVED per instance (exact Gröbner / slope-residual decision):**
- at `d = 1`, for degree-3 AP `r ∈ {0,−4,1,−1,2,3}`, `{cascade} ∪ {Q₀ = 1}` is
  empty and `{cascade} ∪ {Q₀ = 0}` is feasible; this is direct corroboration of
  the symbolic-`r` slope-zero theorem;
- at `d = 2`, `{cascade} ∪ {Q₀ = 1}` is empty for W1 and W2 in this verifier, and
  for AP `r = 1,−1,3` in the backing Wave-B verifier; W1 additionally has the exact
  slope certificate via `w = a₁₂²(a₂₀ − 4a₂₂)`, while no slope-zero or `Q₀ = 0`
  claim is made here for the other `d = 2` AP cases;
- W1 at `d = 3,4` has slope forced to `0` by a sole nonzero multiple of `c`;
- the AP tops of degrees `3,6,9` at `d = 1` have `{cascade} ∪ {Q₀ = 1}` empty;
  for the degree-6 and degree-9 checks, neither slope zero nor `Q₀ = 0` feasibility
  is established.

**RESIDUAL (honestly open, precisely delimited):**
- exceptional `r` loci in the degree-3 AP `d = 1` family, beyond the six values
  individually checked; generic rational elimination does not classify every
  possible rank or denominator specialization;
- the **uniform-in-free-degree** closure. `d = 1` has a generic-`r` theorem plus
  six direct specializations; selected
  instances at `d = 2, 3, 4` are proved. For W1, the slope residual is direct at
  `d = 1,3,4` while `d = 2` carries the controlling combination `w`; the other
  `d = 2` cases are only infeasibility decisions. No single closed-form certificate
  uniform in `d` is extracted here. This is the same status frontier as the classical sibling's
  leading-coefficient loci (`classical-cube-closure.md` §7, `99fe6ee`) — the
  tropical skeleton (§4) reduces the problem but does not carry the uniform-`d`
  proof in closed form.
- **all non-AP realizable exotic tops at `deg a₃ ≥ 6`**. In the verifier's
  bounded normalized-root window (integer roots, `0` included, remaining roots
  in `1..14`, equivalently max root `≤ 14`), exactly 22 occur at degree 6
  (`verify §1`). Those 22 establish that the AP class is not exhaustive; they do
  not exhaust the residual. Higher-degree non-AP exotic tops are not treated here.

**NOT claimed:** a full band-3 theorem; any DC1/JC2 statement (no counterexample
is produced, consistent with DC1, [EXO] §5); closure of the entire degree-3
`d = 1` parameter family, including unchecked exceptional `r` loci; closure of
the `b₂ ≠ 0` half uniform in free degree; the non-AP higher-degree exotic sector.

## 9. Relation to the induction and to astar-band3

[AST] §6 reduced **quantum A\*-I** to "is the reflected non-shifted-cube-necklace
top killed downstream?" — the open question this memo advances. Theorem 1 answers
it generically in `r` for the degree-3 AP family at `d = 1`, with six direct
specializations but unchecked exceptional loci left open; selected per-instance
proofs cover W1/W2 plus AP `r=1,−1,3` at
`d = 2`, W1 at `d = 3,4`, and AP top degrees `3,6,9` at `d = 1`. These are exact advances toward the exotic
gatekeeper, not a completed band-3 induction rung: uniform free degree and all
higher-degree non-AP tops remain open. No DC1/JC2 or full-band-3 conclusion is
claimed.

## 10. Verification

```sh
uv run --with sympy python research/band3/verify_quantum_exotic_closure.py
```
Exact SymPy; §0 checks the telescoping identity from the stipulated `Q_m`
convention, membership `G(0)=0`, and the slope formula; it does not independently
re-derive or validate the crossed-product engine. §1 checks the degree-3 exotic =
AP structural scan and a bounded degree-6 non-AP search (normalized integer roots,
`0` included, remaining roots in `1..14`, max root `≤14`); §2 checks bottom
membership-protection; §3 corroborates the Lemma-R block leading coefficients at
selected degree pairs; §4 checks the generic-`r`, `d = 1` rational elimination and
explicit W1 triangular annihilation, without exhausting exceptional `r` loci; §5 the
selected class/degree instances and exact `d = 2` certificate; §6 the separate
positive control. A successful run prints 53 `PASS` lines and ends
`ALL QUANTUM EXOTIC CLOSURE CHECKS PASSED`.
