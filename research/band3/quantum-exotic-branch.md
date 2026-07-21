# The quantum band-3 exotic branch: moment obstructions in proved slices

**INDEPENDENTLY DERIVED AND MACHINE-VERIFIED — NOT PEER REVIEWED — BAND-SCOPED**

QUANTUM EXOTIC BRANCH assault (Wave B). This memo decides the fate of the
**necklace gap** exposed by the band-3 reconnaissance: the quantum `Q₅` wall is
strictly weaker than its classical sibling, admitting *non-shifted-cube* tops
`a₃` (the exact refutation of the shifted-cube conjecture, `quantum-band3-cascade.md`
§3, commit `99fe6ee`; `band-k-weapons.md` W2q, `99fe6ee`). Band 3 is the first band where the quantum
theory is *looser* than the classical one at the gatekeeper, and the open
question it left was: **does a wall-admitting non-shifted-cube top survive the
full band-3 system + genuine `A₁` membership, or is it killed downstream?**

> **Verdict (scoped).** For `b₂ = 0`, the verifier fixes the degree-3 W1 top and
> checks the prerequisite collapse kernels for polynomial ansätze through degree
> `≤ 7`; once those collapse hypotheses hold, the final `L₀` leading-degree
> obstruction is valid for arbitrary `deg b₋₃`. This is not an arbitrary-top or
> arbitrary-free-degree proof of emptiness for the whole `b₂ = 0` exotic
> sub-branch. For `b₂ ≠ 0`, the `Q₀` moment-unit obstruction proves the generic-`r`
> part of the degree-3 AP family at free-data degree `d = 1`; six exact
> specializations are also checked, while all other exceptional rank/denominator
> loci remain open. Further bounded/sliced cases cover W1/W2 and AP `r=1,−1,3`
> at `d = 2`, W1 at `d = 3,4`, and AP top degrees `3,6,9` at `d = 1`. A separate
> arbitrary-degree theorem excludes `Q₀=1` for the fixed normalized W1 datum under
> the positive cascade and genuine membership; see
> [`quantum-w1-arbitrary-degree.md`](quantum-w1-arbitrary-degree.md). Uniform closure
> outside normalized W1 remains open, including every other AP datum and all
> higher-degree non-AP tops. No DC1/JC2 counterexample, full band-3 theorem, or unconditional
> induction rung is produced.

The displayed finite computations are checked exactly by
[`verify_quantum_exotic.py`](verify_quantum_exotic.py) (ends
`ALL QUANTUM EXOTIC CHECKS PASSED`, 49 exact checks, ~2 s); the conditional
arbitrary-degree `L₀` conclusion uses the written leading-term proof in §2.

## 0. Setup, the branch, and the two sub-cases

Conventions frozen exactly as in the quantum band-3 cascade (`quantum-band3-cascade.md`,
`99fe6ee`): `A₁[x⁻¹] = ⊕_k x^k ℂ[E]`, `(x^a f)(x^b g) = x^{a+b} f^[b] g`,
`f^[r](E) = f(E+r)`, ladder-`m` coefficient
`Q_m = Σ_{k+l=m}(b_l^[k] a_k − a_k^[l] b_l)`, and `[D,X]=1 ⇔ Q_m = δ_{m0}`
(`m ∈ [−6,6]`). Genuine `A₁` membership: `E(E−1)⋯(E−r+1) | a_{−r}, b_{−r}`.
`verify §0` checks the stipulated `Q_m` convention and isolates, in the gauge
`b₃ = 0` (from `Q₆`, `quantum-band3-cascade.md` §2), the
**descent operator**
```
   L_m[b] := b^[3] a₃ − a₃^[m−3] b        (the (3, m−3) pair of Q_m),
```
so that `Q_m = L_m[b_{m−3}] + (lower pairs)` for `m = 5,4,3,2,1,0`.

**The exotic branch.** Fix an exotic top: `a₃ ≠ 0` solving the wall
`b₂^[3] a₃ = a₃^[2] b₂` for some `b₂`, with `a₃` **not** a shifted cube
`c·h h^[1] h^[2]`. The two Wave-A witnesses (`verify §1`, both certified
non-shifted-cube by an exact ∅-solve):
```
   W1:  a₃ = E(E−2)(E−4)   (roots {0,2,4}),   b₂ = (E−1)(E−4);
   W2:  a₃ = E(E+2)(E+4)   (roots {0,−2,−4}), b₂ = E(E+3).
```
By the Wall Lemma (`quantum-band3-cascade.md` §3, `99fe6ee`) the wall solution
`b₂` is unique up to a scalar `κ₂` (freedom exactly 1-dimensional, `κ₂ = 0`
included). This splits the branch into two exhaustive sub-cases, treated separately:

- **`b₂ = 0`** (`κ₂ = 0`): §2 — bounded collapse checks for W1 through ansatz
  degree `≤ 7`, followed conditionally by an arbitrary-`deg b₋₃` `L₀` obstruction.
- **`b₂ ≠ 0`** (`κ₂ ≠ 0`, the genuinely new sub-case): §3 — empty in the
  proved bounded/sliced cases at `Q₀`; uniform closure remains open.

By the diagonal-plus-rescaling symmetry `x ↦ ρx`, `D ↦ ρ³D`, `X ↦ ρ⁻³X` with
`ρ⁵ = 1/κ₂` (which fixes `a₃` and sends `κ₂ ↦ 1`), we normalize `b₂ = (E−1)(E−4)`
etc. in the `b₂ ≠ 0` sub-case; the `b₂ ≠ 0` exotic subcase is exactly `κ₂ ≠ 0`.

## 1. The descent operators and their kernels (structure)

For the fixed W1 top, the bounded kernel checks govern the freedom in the tested
ansätze (the `Q₅` row is the Wave-A wall, `quantum-band3-cascade.md` §3):

| `Q_m` | operator `L_m[b]` | scoped kernel statement |
|---|---|---|
| `Q₅` | `b^[3]a₃ − a₃^[2] b` | 1-dimensional wall freedom `κ₂` |
| `Q₄` | `b^[3]a₃ − a₃^[1] b` | W1: trivial in ansätze through degree 7 |
| `Q₃` | `a₃(b^[3] − b)` | W1: constants in ansätze through degree 5 |
| `Q₂` | `b^[3]a₃ − a₃^[−1] b` | W1: trivial in ansätze through degree 7 |
| `Q₁` | `b^[3]a₃ − a₃^[−2] b` | W1: trivial in ansätze through degree 7 |
| `Q₀` | `b^[3]a₃ − a₃^[−3] b` | conditional written W1 leading-degree result (Lemma 2) |

The necklace comparison suggests why W1 has only the wall freedom at `m=5` and
constants at `m=3`; it motivates and corroborates the bounded kernel pattern.
It is not an unbounded exotic-top kernel theorem, and no family-wide claim is
made for every other positive rung.

## 2. Sub-branch `b₂ = 0`: bounded collapse checks and a conditional `L₀` obstruction

*(For the fixed W1 degree-3 top, the checked collapse is followed by a pure
leading-degree endgame. The scopes of those two steps are different.)*

With `b₂ = 0`, the positive cascade collapses for the W1 top in the polynomial
ansätze checked by `verify §2` (kernel searches through degree `≤ 7`; the `L₃`
constant-kernel check uses degrees `1..5`):
```
   Q₄ = L₄[b₁] = 0  ⇒ b₁ = 0;    Q₃ = a₃(b₀^[3]−b₀) = 0 ⇒ b₀ = const;
   Q₂ = L₂[b₋₁] = 0 ⇒ b₋₁ = 0;   Q₁ = L₁[b₋₂] = 0 ⇒ b₋₂ = 0.
```
(For W1, the stated kernels are machine-checked only in those bounded polynomial
ansätze; the necklace discussion in §1 is motivation, not an unbounded kernel
proof.) Assuming this collapse, with `b₂ = b₁ = b₋₁ = b₋₂ = 0` and `b₀`
constant, every off-diagonal pair of `Q₀` vanishes and
```
   Q₀ = L₀[b₋₃] = b₋₃^[3] a₃ − a₃^[−3] b₋₃  =  1 .
```

> **Lemma 2 (conditional degree obstruction).** For the fixed monic degree-3 W1
> top used by the verifier, and symbolically for arbitrary `q = deg b`, `L₀[b]`
> has degree `q + 2` with leading coefficient `3(3 + q)·lc(b) ≠ 0`.
> Hence, once the preceding collapse is available, `L₀[b]` is either `0` (`b = 0`)
> or of degree `≥ 2`; it can **never** equal the nonzero constant `1`.

*Proof.* The `E^{q+3}` terms of `b^[3]a₃` and `a₃^[−3]b` (`q := deg b`) are both
`lc(b)lc(a₃)E^{q+3}` and cancel. The `E^{q+2}` coefficient is, by the staggered
leading-coefficient identity (`quantum-band3-cascade.md` §5.3, `99fe6ee`, shifts
`(0,3)/(−3,0)`), `((0−(−3))·3 + (3−0)·q)·lc = 3(3+q)·lc(a₃)lc(b) ≠ 0` in
characteristic 0. `verify §2` checks the degree and this coefficient symbolically
for `q = 0..7`. ∎

Therefore, **conditional on the displayed collapse hypotheses**, `Q₀ = 1` is
impossible for arbitrary `deg b₋₃`. **[PROVED conditional `L₀` obstruction;
bounded-verified collapse for W1.]** This does not prove the entire `b₂ = 0`
exotic sub-branch empty for arbitrary top or free-data degree. No membership is
needed for the conditional `L₀` step itself.

## 3. Sub-branch `b₂ ≠ 0`: bounded/sliced kills at `Q₀` (the moment unit)

This is the genuinely new sub-case, with **no positive collapse**. We proceed
constructively — *attempting to build a pair, verifying to destruction* — and
find exact obstructions in the slices stated below. Those finite calculations
are not themselves an arbitrary-degree theorem; the later separate result
[`quantum-w1-arbitrary-degree.md`](quantum-w1-arbitrary-degree.md) supplies one only
for the fixed normalized W1 datum.

### 3.1 The positive cascade is solvable

With `b₂ = (E−1)(E−4)` (normalized) and free lower `X`-coefficients
`a₂, a₁, a₀` (and negatives with membership), `Q₄, Q₃, Q₂, Q₁` forward-solve for
`b₁, b₀, b₋₁, b₋₂` subject to a set of **solvability conditions** on
`(a₂, a₁, a₀)` (each forced operator raises degree by 2, so its image has
codimension ≤ 2 and the source must lie in it). These conditions
are **satisfiable**: `verify §3` exhibits an explicit point on the positive
solution variety of W1,
```
  a₂ = E(E−5),  a₁ = E²−14E+32,  a₀ = 3E − E²/3,
  a₋₁ = E(−2/3 + 2E/3 − E²/6),
```
for which `Q₄ = Q₃ = Q₂ = Q₁ = 0` hold **exactly** with the reconstructed
`b₁, b₀, b₋₁, b₋₂`. So the positive side alone does **not** kill the exotic top:
the naive "wall ⇒ collapse ⇒ empty" route of band 2 / classical band 3 fails
here, precisely as the reconnaissance warned (`band-k-weapons.md` W2q, `99fe6ee`:
"since `u = b_{k−1}` sits at a positive ladder level, this branch is fully
admissible at the gatekeeper stage"). The kill is **downstream, at `Q₀`.**

Note `Q₁` couples the trailing `a₋₁` through the `(−1,2)` pair `b₂^[−1] a₋₁`
(nonzero because `b₂ ≠ 0`), so `a₋₁` is already constrained by the positive
cascade — the first place the negative tail enters. `Q₀` then brings in
`a₋₂`, `a₋₃`, and `μ₃` (via `b₋₃ = μ₃ a₋₃`, `Q₋₆`, `quantum-band3-cascade.md` §5.1).

### 3.2 `Q₀ = 1` is infeasible — the moment-unit obstruction

`Q₀ = 1` is the central integral `G = E` (`quantum-band3-cascade.md` §4,
`99fe6ee`: `Q₀ = (T−1)G`, `G(0)=0` by membership, so `Q₀ = 1 ⇔ G = E`). On the
positive-solution variety, `Q₀ − 1` is a polynomial in `E` whose coefficient
system, together with the free negative data `(a₋₂, a₋₃, μ₃)`, is
**infeasible**:

> **Result 3 (bounded-verified + exact certificate).** For both witnesses W1, W2,
> the system {positive cascade solvability} ∪ {`Q₀ = 1`} has **no solution**
> (`verify §3`, Gröbner basis `= [1]` at free-degree `d = 2`; `d = 1` is checked
> in `verify_quantum_exotic_closure.py`, which also supplies W1 slope certificates
> at `d = 3,4`). Removing the unit — replacing `Q₀ = 1` by `Q₀ = 0` — makes the system
> **feasible** (`verify §3`, Gröbner `≠ [1]`). Hence the obstruction is exactly
> the **moment unit**.

**The exact certificate (W1, `verify §3`).** Eliminating *all* free negative data
(`a₋₂`, `a₋₃`, `μ₃`; the relaxation `μ₃·a₋₃` ↦ free, which is exact because
`μ₃` and `a₋₃` are independent) from the `Q₀ = 1` coefficient system leaves,
on the positive data, **two residual conditions**
```
   8 w = 0        and        7 w = 9,        w := a₁₂² (a₂₀ − 4 a₂₂),
```
where `a₁₂ = lc(a₁)`, `a₂₀, a₂₂` are the constant/leading coefficients of `a₂`.
The higher-degree coefficients of `Q₀` force `w = 0`; the constant coefficient
carries the moment unit and demands `7 w = 9`, i.e. `w = 9/7 ≠ 0`. These are
**contradictory** (`verify §3`, Gröbner of the residual `= [1]`). For W2 the
elimination collapses even more directly to `0 = 1`. This is a *rank/consistency*
infeasibility, not a congruence one — the moment unit (the `1` of `[D,X]=1`, here
in denominator-cleared form `9`) simply cannot be produced once compatibility with
the higher rungs forces the controlling combination `w` to vanish.

**Why the bottom cannot rescue the unit (structure).** The level-3 part of `G`,
`P₃ = Σ_{j=0}^{2} a₃^[j−3] b₋₃^[j]` with `b₋₃ = μ₃ a₋₃`, is *membership-protected*:
its contribution to the constant coefficient of `Q₀ = (T−1)G` is
`P₃(1) − P₃(0) = μ₃·a₃(0)·a₋₃(3)` (all other `a₋₃(1), a₋₃(2), a₋₃(0)` vanish by
`E(E−1)(E−2) | a₋₃`). For the witnesses (`a₃(0) = 0`) this is `0`, so the bottom
data is *absent from the unit equation entirely*; for the additional AP slices with `a₃(0) ≠ 0` the bottom enters the unit equation
but the corresponding exact eliminations still cannot absorb it. This is the
quantum incarnation, within the proved slices, of the classical Theorem A step-6
mechanism (`classical-band3-cascade.md` §6, `99fe6ee`) — the moment can carry only
`τ`/`E`, never a residual unit against a membership-protected extreme — and of the
`W5` "moment-unit-unrealizable" principle (`band-k-weapons.md` §W5, `99fe6ee`; the
`+1` there is this same unit).

### 3.3 Robustness across the exotic class

`verify §4` confirms the `Q₀` kill is not witness-specific:

- **Step-2 arithmetic-progression tops** `{r, r+2, r+4}` (`a₃ = (E−r)(E−r−2)(E−r−4)`,
  `b₂ = (E−r−1)(E−r−4)`) for `r = 0, 1, −1, 3`: all solve the wall, all give
  `{positive} ∪ {Q₀ = 1}` infeasible (`d = 2`). Note `r = 1, −1, 3` have
  `a₃(0) ≠ 0`, so the bottom **does** enter the unit equation there — and still
  cannot save it.
- **A degree-6 exotic top** `{0,2,4,6,8,10}`: root multiset `Φ₃`-divisible
  (wall-admissible), cofactor `A/Φ₃` has a negative coefficient (**not** a shifted
  cube), `b₂` at roots `{1,4,7,10}` solves the wall — and `{positive} ∪ {Q₀ = 1}`
  is infeasible (`d = 1`). This is one exact case past `deg a₃ = 3`; no uniform
  higher-top-degree mechanism is proved.

### 3.4 Validation of the constructive machinery

`verify §5` guards against a false kill: the forward solver is run on the genuine
band-3 positive control `U = x + ∂`, `X = U³ − ∂`, `D = U` (`[D,X] = 1`,
`quantum-band3-cascade.md` §6, `99fe6ee`; this is a `b₂ = 0` *tame* pair). The
solver reproduces its `D` exactly (`b₁ = 1`, `b₋₁ = E`) and emits **no spurious
conditions** — so the pipeline detects real feasibility, and the exotic
infeasibility above is genuine, not an artifact. (The full pipeline's conditions
were also checked to vanish identically at this real pair.)

## 4. Scoped `Q₀` mechanisms in the checked sectors

The kill lives at the **same equation** in every incarnation of the band-3
gatekeeper, though the mechanism differs:

| face / sector | where the top is pinned | how `Q₀`/`C₀` kills |
|---|---|---|
| classical, non-cube (`classical-band3-cascade.md` §6) | wall forces `b₂ = 0`, collapse | `M = τ` vs `τ³ \| b₋₃`: order 1 vs ≥ 3 |
| quantum exotic, `b₂ = 0` (§2) | W1 collapse checked in bounded ansätze | conditional on collapse, `L₀[b₋₃] = 1` impossible for arbitrary `deg b₋₃` |
| quantum exotic, `b₂ ≠ 0` (§3) | **no collapse** | moment **unit** unrealizable (`w = 0` vs `7w = 9`) |

The band-2 wall could kill at the wall itself (membership dispatched the collapsed
tail one rung up). Band 3 is the first band where the quantum wall is too weak to
collapse the branch, and the gatekeeper's real work is done by the **moment** `Q₀`
— the `W4` central integral. In the cases actually computed, this corrects the
naive wall-only gatekeeper: the wall (`Q₅`) does not force a shifted cube, while
the moment (`Q₀`) kills the proved bounded/sliced exotic cases. Whether it
forbids every `b₂ ≠ 0` top uniformly remains open.

## 5. Relation to DC1 / JC2 (no counterexample)

`[D,X] = 1` gives an algebra endomorphism `φ: A₁ → A₁`, `x ↦ X`, `∂ ↦ D`,
automatically injective (`A₁` simple). DC1 asserts `φ` is surjective (an
automorphism) and remains open. The established bridge gives **JC2 ⇒ DC1**;
therefore a genuine DC1 counterexample would refute JC2, but equivalence is not
claimed here. If an exotic pair existed and generated `A₁`, its endomorphism would be an
automorphism and hence tame (Dixmier's theorem). The bounded/scoped blow-up law
recorded in `band3-tame-catalog.md` supplies useful checks on the catalogued tame
examples, but it is not used here as an unbounded classification of every tame
band-3 word. Accordingly §5 draws no global exclusion from that catalog. The
computed exotic slices are empty, so they produce no counterexample; a hypothetical
non-generating exotic pair would be a non-surjective endomorphism and hence a DC1
counterexample. **No counterexample to DC1/JC2 is produced, and none is claimed.**
The exotic wall witness remains a counterexample only to the *shifted-cube
conjecture for the wall equation*, as already recorded in Wave A.

## 6. Scoped gatekeeper advances (induction status)

> **Proved and computed scope.** Let `[D,X] = 1` be a band-3 pair with `a₃ ≠ 0`,
> gauge `b₃ = 0`, and suppose `a₃` is exotic/non-shifted-cube at the `Q₅` wall.
>
> - **`b₂ = 0` half:** for W1, collapse kernels are checked in bounded polynomial
>   ansätze through degree `≤ 7` (with the `L₃` check through degree `5`); conditional
>   on collapse, Lemma 2 excludes `Q₀ = 1` for arbitrary `deg b₋₃`.
> - **`b₂ ≠ 0` half:** the `Q₀` moment-unit obstruction proves the generic-`r`
>   part of the degree-3 AP family at `d = 1`; six exact specializations are also
>   checked; all other exceptional rank/denominator loci remain open. Selected
>   exact instances additionally cover W1/W2 at `d = 2`, W1 at `d = 3,4`, and AP
>   top degrees `3,6,9` at `d = 1` (§3 and `quantum-exotic-closure.md`). The
>   separate theorem [`quantum-w1-arbitrary-degree.md`](quantum-w1-arbitrary-degree.md)
>   excludes `Q₀=1` for the fixed normalized W1 datum, under the positive cascade
>   and genuine membership, at arbitrary coefficient degree.

These results supply exact advances toward the quantum band-3 gatekeeper. They
do **not** restore the induction step unconditionally: uniform closure outside
normalized W1 remains open, including other AP data and all higher-degree non-AP
realizable tops. The W1 theorem does not classify scalar multiples or opposite orientations,
W2, other AP parameters, higher-degree/non-AP tops, or all Band 3. No full band-3,
DC1, or JC2 theorem is claimed.

## 7. Status of claims (proved / computed / conjectured)

**PROVED (with scope explicit):**
- the descent-operator isolation `Q_m = L_m[b_{m−3}] + lower` (§0–§1); the listed
  bounded kernel checks should not be read as an unbounded tame-word or exotic-top
  classification;
- for the fixed W1 top, the collapse kernels in the verifier's stated bounded
  polynomial ansätze; and, **conditional on that collapse**, Lemma 2's `Q₀`
  leading-degree obstruction for arbitrary `deg b₋₃` (§2). Neither statement is an
  arbitrary-top/free-degree emptiness proof for the whole `b₂ = 0` sub-branch.

**VERIFIED (exact, bounded — corroboration with an exact certificate, not an
arbitrary-degree proof):**
- the positive cascade `Q₄..Q₁` is solvable for the exotic top (explicit witness
  point, §3.1);
- `{positive} ∪ {Q₀ = 1}` is **infeasible** in the stated exact slices. The
  homogeneous system `{positive} ∪ {Q₀ = 0}` is also verified feasible for
  W1/W2 at `d = 2` here and for the degree-3 `d = 1` instances listed in
  `quantum-exotic-closure.md`; no homogeneous-feasibility claim is made for the
  other listed slices. In the former cases this isolates the **moment unit** as
  the checked obstruction (§3.2–3.3);
- the exact residual certificate `{8w = 0, 7w = 9}` (W1) and direct W2 Gröbner
  infeasibility (§3.2);
- the pipeline reproduces the genuine positive control with no spurious conditions
  (§3.4).

**OPEN (residual gap):**
- unconditional closure of the `b₂ = 0` half beyond the fixed-top bounded collapse
  checks, including arbitrary exotic tops and arbitrary free-data degree;
- closure of the `b₂ ≠ 0` half uniformly in free degree outside the fixed
  normalized W1 datum, and for all higher-degree non-AP realizable tops. The exact
  finite certificates localize the obstruction to the `Q₀` moment unit in their
  proved slices, while the separate W1 boundary certificate is degree-free only
  for that normalized datum. No complete higher-degree non-AP classification is
  supplied; the uniform branch is not closed.

**NOT claimed:** settlement of DC1 or JC2 (§5); a full band-3 theorem; unconditional
arbitrary-top/free-degree closure of the `b₂ = 0` half; closure of the `b₂ ≠ 0`
half at arbitrary degree outside the fixed normalized W1 theorem; any
classification of scalar multiples or opposite orientations of W1, W2, other AP
parameters, or higher-degree/non-AP tops; any statement about non-exotic
(shifted-cube) tops beyond the gatekeeper reduction.

## 8. Verification

```sh
uv run --with sympy python research/band3/verify_quantum_exotic.py
```
runs §0 (operator isolation from the stipulated `Q_m` convention in gauge `b₃=0`), §1
(both wall witnesses; non-shifted-cube certificates), §2 (the fixed-W1 `b₂=0`
slice: bounded collapse-kernel checks plus direct corroboration of the `L₀` degree
formula for `deg b₋₃=0..7`; the memo's displayed leading-term expansion supplies
the arbitrary-degree formula, conditional on the boundedly checked collapse), §3 (the `b₂≠0` sub-branch: explicit positive
solution; `Q₀=1` infeasible
vs `Q₀=0` feasible via Gröbner; the exact `{8w=0, 7w=9}` certificate; the
positive-control validation), §4 (the exotic AP class and a degree-6 exotic top).
A successful run prints 49 `PASS` lines and ends
`ALL QUANTUM EXOTIC CHECKS PASSED` (~2 s).
