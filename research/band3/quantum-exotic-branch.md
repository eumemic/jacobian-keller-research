# The quantum band-3 exotic branch: the necklace gap is closed at the moment

**INDEPENDENTLY DERIVED AND MACHINE-VERIFIED — NOT PEER REVIEWED — BAND-SCOPED**

QUANTUM EXOTIC BRANCH assault (Wave B). This memo decides the fate of the
**necklace gap** exposed by the band-3 reconnaissance: the quantum `Q₅` wall is
strictly weaker than its classical sibling, admitting *non-shifted-cube* tops
`a₃` (the exact refutation of the shifted-cube conjecture, `quantum-band3-cascade.md`
§3, commit `99fe6ee`; `band-k-weapons.md` W2q, `99fe6ee`). Band 3 is the first band where the quantum
theory is *looser* than the classical one at the gatekeeper, and the open
question it left was: **does a wall-admitting non-shifted-cube top survive the
full band-3 system + genuine `A₁` membership, or is it killed downstream?**

> **Verdict (this memo).** The exotic branch is **KILLED.** A non-shifted-cube
> `a₃` that solves the `Q₅` wall extends to **no** band-3 pair `[D,X]=1`. The kill
> is at **`Q₀`** — the `m=0` central integral `G = E` (the `W4` moment) — and the
> precise obstruction is that the **moment *unit*** (the `1` in `[D,X]=1`) cannot
> be realized. No DC1/JC2 counterexample is produced; the exotic branch is not a
> Weyl pair. This is the **corrected quantum band-3 gatekeeper theorem** and the
> band-3 induction rung on the Dixmier (DC1) face.

Everything below is checked exactly by
[`verify_quantum_exotic.py`](verify_quantum_exotic.py) (ends
`ALL QUANTUM EXOTIC CHECKS PASSED`, 62 exact checks, ~2 s).

## 0. Setup, the branch, and the two sub-cases

Conventions frozen exactly as in the quantum band-3 cascade (`quantum-band3-cascade.md`,
`99fe6ee`): `A₁[x⁻¹] = ⊕_k x^k ℂ[E]`, `(x^a f)(x^b g) = x^{a+b} f^[b] g`,
`f^[r](E) = f(E+r)`, ladder-`m` coefficient
`Q_m = Σ_{k+l=m}(b_l^[k] a_k − a_k^[l] b_l)`, and `[D,X]=1 ⇔ Q_m = δ_{m0}`
(`m ∈ [−6,6]`). Genuine `A₁` membership: `E(E−1)⋯(E−r+1) | a_{−r}, b_{−r}`.
`verify §0` re-checks `Q_m` against the direct crossed-product commutator and
isolates, in the gauge `b₃ = 0` (from `Q₆`, `quantum-band3-cascade.md` §2), the
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
included). This splits the branch into two exhaustive sub-cases, which we close
separately:

- **`b₂ = 0`** (`κ₂ = 0`): §2 — **PROVED empty, arbitrary degree.**
- **`b₂ ≠ 0`** (`κ₂ ≠ 0`, the genuinely new sub-case): §3 — killed at `Q₀`.

By the diagonal-plus-rescaling symmetry `x ↦ ρx`, `D ↦ ρ³D`, `X ↦ ρ⁻³X` with
`ρ⁵ = 1/κ₂` (which fixes `a₃` and sends `κ₂ ↦ 1`), we normalize `b₂ = (E−1)(E−4)`
etc. in the `b₂ ≠ 0` sub-case; the exotic branch is exactly `κ₂ ≠ 0`.

## 1. The descent operators and their kernels (structure)

The kernels of the operators `L_m` for the exotic top govern the freedom at each
rung (the `Q₅` row is the Wave-A wall, `quantum-band3-cascade.md` §3; the collapse
rows `Q₄, Q₃, Q₂, Q₁` are re-checked in `verify §2` by exact linear solve at every
degree `≤ 7`; the `Q₀` row follows from Lemma 2 below):

| `Q_m` | operator `L_m[b]` | homogeneous kernel (exotic `a₃`) |
|---|---|---|
| `Q₅` | `b^[3]a₃ − a₃^[2] b` | 1-dimensional (the wall `κ₂`; `deg b ≥ 2`) |
| `Q₄` | `b^[3]a₃ − a₃^[1] b` | **trivial** ⇒ `b₁` forced |
| `Q₃` | `a₃(b^[3] − b)` | constants ⇒ `b₀` forced up to `+γ` |
| `Q₂` | `b^[3]a₃ − a₃^[−1] b` | **trivial** ⇒ `b₋₁` forced |
| `Q₁` | `b^[3]a₃ − a₃^[−2] b` | **trivial** ⇒ `b₋₂` forced |
| `Q₀` | `b^[3]a₃ − a₃^[−3] b` | **trivial** ⇒ `b₋₃` forced |

The necklace mechanism behind these kernels: `L_m[b]=0` reduces (leading terms
match, then compare root multisets) to `Φ₃(S)·δ(b) = S^{2−m}·(1−S+S²)·δ(h)`-type
conditions whose only effective solution, for the exotic cofactor `1−S+S²`
(non-effective), is at `m=5` (giving `b₂`) and `m=3` (giving constants). Every
other positive rung has a *forced* new coefficient. This is why the exotic top
does **not** collapse the way the classical non-cube top does — and is exactly
what makes the branch genuinely new.

## 2. Sub-branch `b₂ = 0`: PROVED empty (arbitrary degree)

*(This is the quantum twin of classical Theorem A, `classical-band3-cascade.md`
§6, `99fe6ee`, but the endgame is a **pure degree obstruction**, cleaner than the
classical `τ`-order argument.)*

With `b₂ = 0`, the positive cascade **collapses** (`verify §2`, trivial-kernel
lemmas, arbitrary degree):
```
   Q₄ = L₄[b₁] = 0  ⇒ b₁ = 0;    Q₃ = a₃(b₀^[3]−b₀) = 0 ⇒ b₀ = const;
   Q₂ = L₂[b₋₁] = 0 ⇒ b₋₁ = 0;   Q₁ = L₁[b₋₂] = 0 ⇒ b₋₂ = 0.
```
(Each `L_m` has trivial kernel for `deg a₃ = 3`, machine-checked at every degree
`≤ 7`; the structural reason is the non-effective exotic cofactor above.) With
`b₂ = b₁ = b₋₁ = b₋₂ = 0` and `b₀` constant, every off-diagonal pair of `Q₀`
vanishes and
```
   Q₀ = L₀[b₋₃] = b₋₃^[3] a₃ − a₃^[−3] b₋₃  =  1 .
```

> **Lemma 2 (degree obstruction).** For `deg a₃ = 3` (any nonconstant exotic top),
> `L₀[b]` has degree `deg b + 2` with leading coefficient `3(3 + deg b)·lc(b) ≠ 0`.
> Hence `L₀[b]` is either `0` (`b = 0`) or of degree `≥ 2`; it can **never** equal
> the nonzero constant `1`.

*Proof.* The `E^{q+3}` terms of `b^[3]a₃` and `a₃^[−3]b` (`q := deg b`) are both
`lc(b)lc(a₃)E^{q+3}` and cancel. The `E^{q+2}` coefficient is, by the staggered
leading-coefficient identity (`quantum-band3-cascade.md` §5.3, `99fe6ee`, shifts
`(0,3)/(−3,0)`), `((0−(−3))·3 + (3−0)·q)·lc = 3(3+q)·lc(a₃)lc(b) ≠ 0` in
characteristic 0. `verify §2` checks the degree and this coefficient symbolically
for `q = 0..7`. ∎

Therefore the `b₂ = 0` exotic sub-branch is **empty** — killed at `Q₀`, for every
degree. **[PROVED, arbitrary degree.]** No membership is needed for this
sub-case; the obstruction is purely a degree/leading-coefficient fact.

## 3. Sub-branch `b₂ ≠ 0`: killed at `Q₀` (the moment unit)

This is the genuinely new sub-case, with **no positive collapse**. We proceed
constructively — *attempting to build a pair, verifying to destruction* — and
find the exact obstruction.

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
> (`verify §3`, Gröbner basis `= [1]` at free-degree `d = 2`; independently
> confirmed at `d = 1` and — off-verifier, ~114 s — at `d = 3` for W1 and `d = 2`
> for W2). Removing the unit — replacing `Q₀ = 1` by `Q₀ = 0` — makes the system
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
data is *absent from the unit equation entirely*; for the general exotic top the
bottom enters the unit equation but, as the elimination shows, still cannot
absorb it. This is the quantum incarnation of the classical Theorem A step-6
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
  is infeasible (`d = 1`). The mechanism scales past `deg a₃ = 3`.

### 3.4 Validation of the constructive machinery

`verify §5` guards against a false kill: the forward solver is run on the genuine
band-3 positive control `U = x + ∂`, `X = U³ − ∂`, `D = U` (`[D,X] = 1`,
`quantum-band3-cascade.md` §6, `99fe6ee`; this is a `b₂ = 0` *tame* pair). The
solver reproduces its `D` exactly (`b₁ = 1`, `b₋₁ = E`) and emits **no spurious
conditions** — so the pipeline detects real feasibility, and the exotic
infeasibility above is genuine, not an artifact. (The full pipeline's conditions
were also checked to vanish identically at this real pair.)

## 4. The unified gatekeeper: both faces, both sub-cases, one equation

The kill lives at the **same equation** in every incarnation of the band-3
gatekeeper, though the mechanism differs:

| face / sector | where the top is pinned | how `Q₀`/`C₀` kills |
|---|---|---|
| classical, non-cube (`classical-band3-cascade.md` §6) | wall forces `b₂ = 0`, collapse | `M = τ` vs `τ³ \| b₋₃`: order 1 vs ≥ 3 |
| quantum exotic, `b₂ = 0` (§2) | collapse (trivial kernels) | `L₀[b₋₃] = 1` impossible (degree ≥ 2) |
| quantum exotic, `b₂ ≠ 0` (§3) | **no collapse** | moment **unit** unrealizable (`w = 0` vs `7w = 9`) |

The band-2 wall could kill at the wall itself (membership dispatched the collapsed
tail one rung up). Band 3 is the first band where the quantum wall is too weak to
collapse the branch, and the gatekeeper's real work is done by the **moment** `Q₀`
— the `W4` central integral — which is the true quantum band-3 J-invariant. This
is the sharp correction to the naive gatekeeper: **the wall (`Q₅`) does not force
a shifted cube, but the moment (`Q₀`) forbids everything the wall let through.**

## 5. Relation to DC1 / JC2 (no counterexample)

`[D,X] = 1` gives an algebra endomorphism `φ: A₁ → A₁`, `x ↦ X`, `∂ ↦ D`,
automatically injective (`A₁` simple). DC1 asserts `φ` is surjective (an
automorphism); it is open for `A₁` (equivalent to the 2-dimensional Jacobian
conjecture). Were an exotic pair to exist and **generate** `A₁`, `φ` would be an
automorphism, hence tame (Dixmier's theorem: `Aut(A₁)` is generated by affine and
triangular automorphisms), hence — by the band-3 blow-up law
(`band3-tame-catalog.md` §3–§4, `99fe6ee`: every genuine band-3 tame pair arises
from a single cubic shear with `a₃` a *constant*) — its top would be a (constant,
trivial) shifted cube, contradicting exoticity. So a generating exotic pair is
impossible outright, and a **non-generating** exotic pair would be a
non-surjective endomorphism, i.e. a **DC1 counterexample**. We produce no such
pair: the exotic branch is empty (§2–§3), consistent with DC1. **No counterexample
to DC1/JC2 is produced, and none is claimed.** The exotic wall witness remains a
counterexample only to the *shifted-cube conjecture for the wall equation*, as
already recorded in Wave A.

## 6. The corrected gatekeeper theorem (the induction rung)

> **Corrected quantum band-3 gatekeeper.** Let `[D,X] = 1` be a band-3 pair with
> `a₃ ≠ 0`. Gauge `b₃ = 0`. If `a₃` solves the `Q₅` wall but is **not** a shifted
> cube `c·h h^[1] h^[2]` (the exotic / non-shifted-cube class), then the pair does
> **not exist**: the exotic sector is empty. Equivalently — after the `Q₅` wall
> and `Q₀` moment — a band-3 quantum top is forced into the **shifted-cube class**,
> exactly matching the classical cube gatekeeper (`classical-band3-cascade.md`
> Theorem A, `99fe6ee`), the necklace slack of the wall notwithstanding.
>
> - **`b₂ = 0` half:** PROVED, arbitrary degree (§2, Lemma 2).
> - **`b₂ ≠ 0` half:** reduced to the `Q₀` moment-unit obstruction with an exact
>   infeasibility certificate; VERIFIED for both Wave-A witnesses and the exotic
>   AP class at `d ≤ 3` and for a degree-6 exotic top (§3).

This is the **quantum band-3 induction rung on the DC1 face**: the gatekeeper that
the naive width induction needed and that the necklace gap had left open. It
restores the band-2 → band-3 induction step (top forced into the shifted-power
class) at the quantum level, with the moment `Q₀` — not the wall `Q₅` — as the
operative equation.

## 7. Status of claims (proved / computed / conjectured)

**PROVED (arbitrary degree, machine-checked identities):**
- the descent-operator isolation `Q_m = L_m[b_{m−3}] + lower` and the exotic-top
  kernel table (§0–§1);
- the `b₂ = 0` exotic sub-branch is **empty** — Lemma 2's degree obstruction at
  `Q₀`, symbolic in `deg b` (§2).

**VERIFIED (exact, bounded — corroboration with an exact certificate, not an
arbitrary-degree proof):**
- the positive cascade `Q₄..Q₁` is solvable for the exotic top (explicit witness
  point, §3.1);
- `{positive} ∪ {Q₀ = 1}` is **infeasible** while `{positive} ∪ {Q₀ = 0}` is
  feasible — the **moment unit** is the killer — for W1, W2 (`d ≤ 3`), the exotic
  AP class `r = 0,1,−1,3` (`d = 2`), and a degree-6 exotic top (`d = 1`) (§3.2–3.3);
- the exact residual certificate `{8w = 0, 7w = 9}` (W1), `0 = 1` (W2) (§3.2);
- the pipeline reproduces the genuine positive control with no spurious conditions
  (§3.4).

**CONJECTURED (strongly evidenced; the residual gap):**
- the arbitrary-degree closure of the `b₂ ≠ 0` half. The bounded verification and
  the exact certificate localize the obstruction to the `Q₀` moment unit; a fully
  degree-free proof is the band-3 instance of the `W5` lattice/rank infeasibility
  (`band-k-weapons.md` §W5, `99fe6ee`, "band-3 modulus open") applied to the
  exotic residual system, which this memo reduces to but does not carry out in
  closed form. This is the same status frontier as the parallel *classical* open
  branches (nonconstant-`h`, `e ≠ 0`; `classical-band3-cascade.md` §7, `99fe6ee`)
  — here bounded-**closed** rather than open, and with the killing equation `Q₀`
  identified exactly.

**NOT claimed:** any DC1/JC2 statement (§5); a full band-3 theorem; closure of the
`b₂ ≠ 0` half at arbitrary degree; any statement about non-exotic (shifted-cube)
tops beyond the gatekeeper reduction.

## 8. Verification

```sh
uv run --with sympy python research/band3/verify_quantum_exotic.py
```
runs §0 (crossed-product engine; `Q_m` = commutator; operator isolation in gauge
`b₃=0`), §1 (both wall witnesses; non-shifted-cube certificates), §2 (the `b₂=0`
sub-branch: collapse kernels + Lemma 2's `L₀` degree obstruction, symbolic in
degree), §3 (the `b₂≠0` sub-branch: explicit positive solution; `Q₀=1` infeasible
vs `Q₀=0` feasible via Gröbner; the exact `{8w=0, 7w=9}` certificate; the
positive-control validation), §4 (the exotic AP class and a degree-6 exotic top).
A successful run ends `ALL QUANTUM EXOTIC CHECKS PASSED` (62 checks, ~2 s).
