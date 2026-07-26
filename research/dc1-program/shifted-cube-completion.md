# The shifted-cube wall at band 3: the `κ = 0 ∧ b₁ = 0` hole is CLOSED, and closed far wider than the hole

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED — ONE NEW ARBITRARY-DEGREE THEOREM (REPAIRING FIVE CORPUS SECTORS) +
TWO REFUTATIONS + NAMED GAPS**

This memo closes the residual named in [`band3-sectors.md`](band3-sectors.md) §4.1
and §9-open-2 — *the sub-branch `κ = 0` **and** `b₁ = 0`, open in **every**
shifted-cube class including the ones the corpus already reports closed* — and
reports what happened when the mission's route into diff-2 was pushed.

The one-line summary:

> **TARGET A: the hole is not a wall. It is a hole in a floor that was already
> solid, and the repair is much bigger than the hole.**
>
> With gauge `b₃ = 0`, the `Q₅` wall `b₂ = κ h h^{[1]}` at `κ = 0`, and `b₁ = 0`,
> the operator `D` has **band ≤ 0**. That single structural fact — *not* anything
> about `h` — collapses the whole ladder:
>
> ```
> Q₃ = a₃ (T³−1) b₀                                    ⇒  b₀ constant
> Q₂ = b₋₁^{[3]}a₃ − a₃^{[-1]}b₋₁      (b₀ const)      ⇒  b₋₁ = 0
> Q₁ = b₋₂^{[3]}a₃ − a₃^{[-2]}b₋₂      (b₋₁ = 0)       ⇒  b₋₂ = 0
> Q₀ = b₋₃^{[3]}a₃ − a₃^{[-3]}b₋₃  = 1 (b₋₂ = 0)       ⇒  CONTRADICTION
> ```
>
> **Every coefficient of `X` except the top drops out identically at every rung.**
> The kill is one line of leading-coefficient arithmetic: for nonzero `φ, a` of
> degrees `n, A`,
>
> ```
> [E^{N−1}]( φ^{[k]} a − a^{[−j]} φ )  =  lc(φ) lc(a) · ( k n + j A ),    N = n + A,
> ```
>
> and `k n + j A` is **strictly positive** whenever `n ≥ 1`. Membership
> `(E)_j | b₋ⱼ` forces `n ≥ j`, so each rung kills its `b₋ⱼ` outright, and the
> moment unit `Q₀ = 1` has nowhere left to come from.
>
> > **NONPOSITIVE-D EXCLUSION THEOREM (arbitrary degree, arbitrary band `k ≥ 2`, ARBITRARY nonzero top
> > `a_k`).** There is **no** pair with `band X = k ≥ 2`, `band D ≤ 0`,
> > membership-valid negative tail, and `[D,X] = 1`.
>
> The sector `κ = 0 ∧ b₁ = 0` is therefore **EMPTY at arbitrary degree, in every
> shifted-cube class at once** — cube-separated, 2-separated, diff-1, diff-2,
> diff-3, `(E−r)²`, and the **A\*-band3 constant-top corner `κ₂ = 0 ∧ b₁ = 0`**
> (not all `κ₂ = 0`; see §5). The corpus's cube-separated / 2-separated / diff-1
> closures are **REPAIRED**, not broken. No new wall.
>
> **TARGET B: the mission's route into diff-2 is REFUTED with a nonzero-branch
> witness (`κ ≠ 0`, `λ ≠ 0`).** On the branch `b₁(r+1) = 0`, the two analyzed
> clean `Q₂` nodes (at `r` and `r+1`) do force a genuine new degree-free
> condition `(E−r−1) | a₁` (with `κ ≠ 0` and `a₂(r−1) ≠ 0`, i.e. on the intended
> surviving branch), so the [`band3-sectors.md`](band3-sectors.md) §5.1 witness
> with `a₁ = −3` does not extend past `Q₃` for that specific `a₁`. The published
> Q₄, Q₃ witness was correct for its stated Q₄, Q₃ scope; the new content is
> only that it does not extend one more rung with `a₁ = −3`. Repairing the
> witness restores it completely: at symbolic `(r, κ)` and any scale `λ`,
> ```
> a₂ = 3λ(E−r)(E−r−1),  b₁ = 2κλ(E−r−1)²,  a₁ = −3λ²(E−r−1)²,  b₀ = 0,  a₀ = c
> ```
> solves `Q₆ = Q₅ = Q₄ = Q₃ = Q₂ = Q₁ = 0` **exactly**. Interpreted as a point of
> the intended surviving branch (`b₁(r+1) = 0` with `a₂(r−1) ≠ 0`) it requires
> **both** `κ ≠ 0` and `λ ≠ 0`; at `κ = 0` the entire `b₁` vanishes and at `λ = 0`
> the whole family collapses to zero. With both nonzero (e.g. `κ = λ = 1`)
> `a₂(r−1) = a₂(r+2) = 6λ ≠ 0` and the refined family shows the entire positive
> cascade fails to close the diff-2 branch. The kill must come from `Q₀ = 1` +
> membership + the negative tail. This memo localises exactly what is missing:
> **two scalar covectors** `C₁, C₃` on the tail.

Exact certificate:
[`verify_shifted_cube_completion.py`](verify_shifted_cube_completion.py). Every
load-bearing fact in the Nonpositive-D Exclusion Theorem (the crossed-product
ladder engine, `Q₀ = (T−1)G`, `G(0) = 0` under membership, and the
`(T^n−1)` kernel) is **re-derived in file**, not cited. The shifted-cube gauge
`b₃=0` and wall shape `b₂=κhh^{[1]}` are frozen corpus inputs, not independently
proved here; the verifier's `Q₅` computation is only the bounded concrete
cross-check stated in §6. Historical run counts (`100/100`, `112/112`) predate
the 2026-07-26 audit edits and are **superseded** unless re-executed on this
checkout; no new counts are asserted here.

Conventions frozen from the corpus: `A₁[x^{-1}] = ⊕_k x^k C[E]`, `E = x∂`,
`(x^a f)(x^b g) = x^{a+b} f(E+b) g(E)`, `f^{[n]}(E) = f(E+n)`, `T f = f^{[1]}`,
`Q_m = ∑_{k+l=m}[b_l^{[k]}a_k − a_k^{[l]}b_l]`, `[D,X]=1 ⇔ Q_m=δ_{m0}`,
membership `(E)_j | a_{-j},b_{-j}`, gauge `b₃=0`, `Q₀=(T−1)G`, `G(0)=0`.
Sector: `a₃ = h h^{[1]}h^{[2]}`, `b₃ = 0`, `b₂ = κ h h^{[1]}`;
`diff-1: h = (E−r)(E−r−1)`, `diff-2: h = (E−r)(E−r−2)`.

---

## 0. What the hole was

[`band3-sectors.md`](band3-sectors.md) §4.1 named it precisely, and the naming was
correct:

> *In the sub-branch `κ = 0` and `b₁ = 0`, `Q₃ = a₃(T³−1)b₀` only forces `b₀`
> constant and puts **no** condition on `a₂`. This gap is `h`-independent: the
> corpus's cube-separated `κ = 0` chain
> ([`shifted-power-residuals.md`](shifted-power-residuals.md) §1.2,
> "`Q₄ ⇒ b₁ = c h`, `Q₃ ⇒ h h^{[1]} | a₂`") silently assumes `c ≠ 0` in the second
> step.*

The verifier confirms the diagnosis before repairing it: with `b₃ = b₂ = b₁ = 0`
and coefficients entered as **abstract `sympy.Function`s**, `Q₃` literally does
not contain `a₂` (verifier `§5`, degree-free). So the chain that closes the
cube-separated, 2-separated and diff-1 classes really does have a gap at its
first rung, in **all** of them.

The mission asked whether this is a repair or a **new wall**. It is a repair —
and the reason is that the corpus was looking at the wrong invariant. The
sub-branch is not "a shifted-cube sector with a degenerate constant"; it is the
sector **`band D ≤ 0`**, and that is killable without ever mentioning `h`.

---

## 1. The ladder collapse (degree-free, and `X`-tail-free)

Gauge `b₃ = 0`, wall `b₂ = κ h h^{[1]}` with `κ = 0`, and `b₁ = 0` give
`b_l = 0` for every `l ≥ 1`. So `D = b₀ + x^{-1}b₋₁ + x^{-2}b₋₂ + x^{-3}b₋₃`.

For `k + l = m` with `l ≤ 0` and `|k| ≤ 3`, the surviving pairs are few, and
`b₀` being **constant** annihilates every `s = 0` term (`b₀^{[k−j]}a_{k−j} −
a_{k−j}b₀ = 0`). Verifier `§1` derives, with **all** coefficients of `X` and `D`
entered as undetermined functions (no polynomial ansatz, no degree cap, `k = 3, 4, 5`):

```
Q_k       = a_k (T^k − 1) b₀
Q_{k−j}   = b_{−j}^{[k]} a_k − a_k^{[−j]} b_{−j}      (b₀ const, b_{−1}…b_{−(j−1)} = 0),  j = 1..k
```

> **Note what is absent.** `a_{k−1}, …, a_0, a_{−1}, …, a_{−k}` appear in **none**
> of these rungs. The `κ = 0 ∧ b₁ = 0` sector is decided entirely by
> `(a_k, b₀, b_{−1}, …, b_{−k})`. In particular the corpus's worry — that `Q₃`
> constrains no `a₂` — is correct and **irrelevant**: nothing needs to.

A negative control confirms the collapse is not an artefact: restoring `b₁ ≠ 0`
makes `Q₂` acquire extra terms immediately.

---

## 2. The kill: one leading-coefficient identity, at symbolic degree

For a nonzero `P` of degree `p`, write `τ(P) = [E^{p−1}]P / [E^p]P` (minus the sum
of the roots). Two rules suffice, and both are **derived by sympy with the degree
left as a free symbol** — verifier `§2` expands in `u = 1/E`:

```
τ(P^{[s]}) = τ(P) + s·deg P            (from (1+su)^p = 1 + p s u + O(u²)),
τ(P·Q)     = τ(P) + τ(Q).
```

Both are then **instance-validated against real polynomial algebra** (degrees 1–6,
shifts `−3..3`) with negative controls: the degree-blind rule `τ ↦ τ + s` fails.

Apply them to a rung. `A = φ^{[k]}a` and `B = a^{[−j]}φ` have the *same* degree
`N = n + A_deg` and the *same* leading coefficient, so the `E^N` terms cancel
identically, and

> **KEY FORMULA (symbolic `n, A, k, j`; verifier `§2`).**
> ```
> [E^{N−1}]( φ^{[k]} a − a^{[−j]} φ )  =  lc(φ)·lc(a)·( k·n + j·A ).
> ```

Instance-validated on real polynomials for `(k,j) ∈ {(3,1),(3,2),(3,3),(4,2),(5,4),(2,1)}`
× `deg φ ∈ 1..4` × `deg a ∈ 0..4`, including the `E^N` cancellation; the
sign-flipped variant `k n − j A` fails on the same instances.

**Now the chain (verifier `§3`).** Let `k ≥ 2`.

1. **`Q_{k−j} = 0` for `j = 1..k−1`.** If `b₋ⱼ ≠ 0` then membership `(E)_j | b₋ⱼ`
   gives `n = deg b₋ⱼ ≥ j ≥ 1`, so `N ≥ 1` and the `E^{N−1}` coefficient of the
   rung must vanish: `lc·lc·(k n + j A) = 0`. But `k n + j A > 0` strictly
   (`n ≥ 1`, `A ≥ 0`, `k, j ≥ 1`). **Contradiction ⇒ `b₋ⱼ = 0`.**
2. **`Q₀ = 1`.** With `b₋₁ = ⋯ = b₋₍ₖ₋₁₎ = 0`, `Q₀ = b₋ₖ^{[k]}a_k − a_k^{[−k]}b₋ₖ`.
   If `b₋ₖ = 0` this reads `0 = 1`. Otherwise `n ≥ k ≥ 2`, so `N = n + A ≥ 2` and
   `N − 1 ≥ 1`: the `E^{N−1}` coefficient must **vanish** (the right-hand side `1`
   lives in degree `0` only), i.e. `k n + k A = 0`. **Contradiction.**

> ### NONPOSITIVE-D EXCLUSION THEOREM (arbitrary degree; arbitrary `k ≥ 2`; **arbitrary** nonzero top `a_k`)
> There is no pair `(X, D)` in `A₁[x^{-1}]` with `band X = k ≥ 2`, `band D ≤ 0`,
> membership-valid negative tail `(E)_j | b₋ⱼ`, and `[D,X] = 1`.
>
> **Machine scope: degree-free and `k`-free in the arithmetic; the rung collapse
> is machine-derived with abstract functions at `k = 3, 4, 5`, and the `k`-uniform
> pattern is visible in the derivation (the `s = 0` term is killed by `b₀` constant
> for every `k`).**

### 2.1 Why `k = 1` survives, and why membership is the load-bearing hypothesis

Both boundary cases are **machine-exhibited**, not argued (verifier `§4`):

- **`k = 1` must survive, and does.** `D = ∂ = x^{-1}E`, `X = x` has `band D = −1 ≤ 0`,
  `band X = 1`, genuine membership `E | b₋₁`, and the verifier **computes**
  `[D,X] = 1` in the crossed product. The escape is exactly `N = 1`: with
  `deg b₋₁ = 1`, `deg a₁ = 0`, the `E^{N−1}` coefficient **is** the constant term,
  and `k n + j A = 1` is the moment unit rather than a contradiction. A theorem
  that killed band 1 would be false; this one does not.
- **Membership is what does the killing.** `X = x³`, `D = c + x^{-3}(E/3)`
  satisfies `[D,X] = 1` **exactly** — verifier computes it — with `band D ≤ 0`,
  `band X = 3`. It is not a counterexample because `(E)_3 = E(E−1)(E−2)` does not
  divide `E/3`, so `D ∉ A₁`. This is precisely the `deg b₋₃ = 1 < 3` escape the
  key formula predicts, and it certifies that Nonpositive-D Exclusion Theorem's hypothesis set is
  **tight**: drop membership and the sector is nonempty *with the moment unit*.

### 2.2 Bounded cross-checks and non-vacuity

Verifier `§4`, exact SymPy over `ℚ`, cap `d ∈ {1,2}`, the **full** cascade
`Q_m = δ_{m0}` for `m ∈ [−6,6]` plus membership, at six different tops:

| top `h` | `d = 1` | `d = 2` |
|---|---|---|
| diff-1 `E(E−1)` | unit ideal ⇒ EMPTY | unit ideal ⇒ EMPTY |
| diff-2 `E(E−2)` | unit ideal ⇒ EMPTY | unit ideal ⇒ EMPTY |
| cube-separated `E(2E−1)` | unit ideal ⇒ EMPTY | unit ideal ⇒ EMPTY |
| diff-3 `E(E−3)` | unit ideal ⇒ EMPTY | unit ideal ⇒ EMPTY |
| double root `E²` | unit ideal ⇒ EMPTY | unit ideal ⇒ EMPTY |
| constant top `h = 1` | unit ideal ⇒ EMPTY | unit ideal ⇒ EMPTY |

These are an **independent cross-check** of Nonpositive-D Exclusion Theorem, not its evidence.
Non-vacuity is by **explicit point**, never by a normal form (house rule): drop
`Q₀ = 1` and `(a₃ = h h^{[1]}h^{[2]}, b₀ = 5, all else 0)` is a point of the sector
with every `Q_m = 0` (commutator computed `= 0`); drop membership and §2.1's
`(x³, c + x^{-3}E/3)` is a point *with* `Q₀ = 1`.

---

## 3. What Nonpositive-D Exclusion Theorem repairs

The theorem never mentions `h`. Every class whose `κ = 0 ∧ b₁ = 0` sub-branch was
open inherits the closure verbatim (verifier `§5`):

| corpus claim | what was actually missing | status now |
|---|---|---|
| cube-separated `h` closed at arbitrary degree ([`shifted-power-residuals.md`](shifted-power-residuals.md) §1.2) | the `κ = 0` chain's step `Q₃ ⇒ h h^{[1]} \| a₂` assumed `c ≠ 0` in `b₁ = c h` | **REPAIRED** |
| 2-separated `h`, incl. the whole diff-3 class ([`shifted-power-residuals.md`](shifted-power-residuals.md) §2.3) | same silent `c ≠ 0` | **REPAIRED** |
| diff-1 closed at arbitrary degree ([`band3-sectors.md`](band3-sectors.md) §4) | stated for `κ ≠ 0`, and for `κ = 0` with `b₁ ≠ 0` | **REPAIRED — diff-1 is now closed for ALL `κ` and ALL `b₁`** |
| diff-2 `κ = 0` sub-branch | same | **REPAIRED** (diff-2's `κ ≠ 0` branch is still open — §4) |
| A\*-band3 constant top with gauged `κ₂ = 0` and `b₁ = 0` | not previously isolated | **CLOSED** as a special case (`a_k` arbitrary) |

> *(AUDIT SCOPE, 2026-07-25: only the third branch — `κ = 0 ∧ b₁ = 0` — is proved
> in THIS verifier, by Nonpositive-D Exclusion Theorem. The other two branches (`κ ≠ 0`, and `κ = 0 ∧
> b₁ ≠ 0`) are **cited** from [`band3-sectors.md`](band3-sectors.md) §4/§4.1 and are
> not re-derived here, so the "no side condition" statement is the combination of
> this arm with that memo, not a self-contained result of this file.)*
>
> **The diff-1 class is now closed at arbitrary degree with no side condition at
> all.** [`band3-sectors.md`](band3-sectors.md) §4's theorem carried the caveat
> "*for `κ ≠ 0`, and for `κ = 0` with `b₁ ≠ 0`*". That caveat can be dropped.

**This is not a new wall.** The mission asked to say so loudly if it were; it is
not, and the finding is the opposite — a corpus-wide repair from a single
`h`-independent theorem.

---

## 4. TARGET B — diff-2's surviving branch: the positive cascade does not close it

The mission's first Target-B bullet was: *"on the branch `b₁(r+1) = 0`, what do
the deeper rungs force?"* The supported answer here is **one new condition from
the two analyzed clean `Q₂` nodes, plus an explicit refined family showing that
the cascade `Q₂, Q₁` does not close the branch**.

### 4.1 The new rung (degree-free) — the two analyzed clean `Q₂` nodes

Verifier `§6` re-derives the diff-2 node table with `a₂, a₁, a₀, b₁, b₀` as
undetermined functions, then finds that the tail-free node set of `Q₂` is
**exactly `{r−1, r, r+1, r+2}`** (exclusivity machine-tested over offsets
`−5..5`). The two nodes actually used for the displayed forcing are `r` and
`r+1`:

```
Q₂(r)   = a₁(r)b₁(r+1) − a₁(r+1)b₁(r) − a₂(r)(b₀(r) − b₀(r+2)),
Q₂(r+1) = a₁(r+1)b₁(r+2) − a₁(r+2)b₁(r+1) − a₂(r+1)(b₀(r+1) − b₀(r+3)).
```

The other two tail-free nodes `r−1` and `r+2` are **unclassified forcing data**
in this memo — they are not exploited for the displayed forcing and no claim is
made about what they add.

On the `Q₄` locus (`a₂(r) = a₂(r+1) = 0`, `b₁(r) = b₁(r+2) = κa₂(r−1)/3`) **and**
the surviving branch `b₁(r+1) = 0`, the two used nodes collapse to
`∓(κ/3)·a₂(r−1)·a₁(r+1)`. With `κ ≠ 0` and `a₂(r−1) ≠ 0` — the defining property
of the intended surviving branch —

> **`a₁(r+1) = 0`, i.e. `(E−r−1) | a₁`.**  *(from the two ANALYZED clean nodes,
> arbitrary degree)*

A negative control records that the same two nodes say nothing about `a₁(r)` or
`a₁(r+2)`.

### 4.2 The published witness, the specific `a₁ = −3` non-extension, and the refined family

[`band3-sectors.md`](band3-sectors.md) §5.1's witness has `a₁ = −3`, so
`a₁(r+1) = −3 ≠ 0`. The verifier re-derives that witness independently
(it does solve `Q₅ = Q₄ = Q₃ = 0`) and then **shows `Q₂ ≠ 0` on it**.

> **HISTORY (scoped, audit 2026-07-26).**
> The [`band3-sectors.md`](band3-sectors.md) §5.1 witness was **correct for its
> stated `Q₄, Q₃` scope**. It is not smaller than reported at that scope. What
> is new here is only that it does not extend one more rung with the specific
> choice `a₁ = −3`: the two analyzed clean `Q₂` nodes force `(E−r−1) | a₁`, and
> `a₁ = −3` violates that. The member of the scaled family that does extend has
> `a₁ = −3λ²(E−r−1)²`.

The extending witness is exact:

> **REFINED WITNESS (verifier `§6`; symbolic `r`, `κ`, scale `λ`, constant `c`).**
> ```
> h  = (E−r)(E−r−2),   a₃ = h h^{[1]}h^{[2]},   b₃ = 0,   b₂ = κ h h^{[1]},
> a₂ = 3λ(E−r)(E−r−1), b₁ = 2κλ(E−r−1)²,  a₁ = −3λ²(E−r−1)²,  a₀ = c,  b₀ = 0.
> ```
> `Q₆ = Q₅ = Q₄ = Q₃ = Q₂ = Q₁ = 0` **exactly**, for every `r`, `κ`, `λ`, `c`.
> Interpreted as a point of the intended surviving branch (`b₁(r+1) = 0` with
> `a₂(r−1) ≠ 0`, i.e. `h h^{[1]} ∤ a₂`), this requires **both `κ ≠ 0` and
> `λ ≠ 0`**: at `κ = 0` the whole `b₁` vanishes (so "the surviving branch
> `b₁(r+1) = 0`" is vacuous), and at `λ = 0` the family collapses to zero (so
> `a₂(r−1) = 0`, not the branch's defining `≠ 0`). At the concrete
> specialization `κ = λ = 1` (both nonzero, machine-checked) `a₂(r−1) = 6 ≠ 0`
> and the clean divisibility is genuinely absent.
>
> The symbolic equality `a₂(r−1) = a₂(r+2) = 6λ` is **not by itself a
> nonvanishing test** — it collapses at `λ = 0`; the branch interpretation needs
> `λ ≠ 0` (and `κ ≠ 0`, as above) to hold.
>
> The structural reason `Q₂` vanishes on the refined family: `a₁` and `b₁`
> become **proportional** (both multiples of `(E−r−1)²`), so the `Q₂` middle
> term `b₁^{[1]}a₁ − a₁^{[1]}b₁` is identically zero.
>
> **SCOPE GUARD (machine-checked).** This is **not** a candidate Weyl pair. With
> the zero negative tail shown, `Q₀ = 0 ≠ 1`, so `[D,X] ≠ 1`. It solves the
> *positive cascade only*. Nothing in this memo constructs a genuine pair, and
> nothing here was (or needed to be) run through
> [`sieve_dc1_candidate.py`](sieve_dc1_candidate.py).

> ### TARGET B VERDICT (machine-checked; supported scope)
> **The "push to `Q₂` and `Q₁`" route is REFUTED by the explicit refined family
> with `κ ≠ 0` and `λ ≠ 0`, which solves the entire positive cascade
> `Q₅…Q₁ = 0` with `a₂(r−1) ≠ 0`.** Supported statements:
> (i) the two analyzed clean `Q₂` nodes force `(E−r−1) | a₁`;
> (ii) the explicit refined family shows that `Q₂` and `Q₁` do **not** close
> the branch.
> **NOT claimed** as universal statements: "`Q₂` adds exactly one new
> condition" (the two extra tail-free nodes `r−1, r+2` are unclassified) or
> "`Q₁` adds none" (demonstrated on the refined family only, not universally).
> Any kill must use `Q₀ = 1`, membership, and the negative tail.

### 4.3 The moment-unit route, localised to two covectors

The mission's second Target-B bullet asked for the correct `G`-divisor. There is
none available beyond `h^{[-1]}`, and the verifier says precisely **why and by how
much** (`§7`, degree-free). Splitting the potential,

```
G = ( five terms carrying a₃^{[-1]}, a₃^{[-2]}, a₃^{[-3]}, b₂^{[-1]}, b₂^{[-2]} )  +  R,
R = a₁^{[-1]}b₋₁ − b₁^{[-1]}a₋₁ + a₂^{[-2]}b₋₂ + a₂^{[-1]}b₋₂^{[1]},
```

the five listed terms are **automatically** divisible by `h^{[-1]}`, and `R` is
not — verified, with the negative control that none of `a₂^{[-1]}, a₂^{[-2]},
a₁^{[-1]}, b₁^{[-1]}` is `h^{[-1]}`-divisible even on the refined branch shapes.
Since `h^{[-1]} = (E−r−1)(E−r−3)` has two **distinct simple** roots for diff-2,

> **`h^{[-1]} | G` ⟺ `C₁ = C₃ = 0`, where (on the `Q₄` locus)**
> ```
> C₁ = a₁(r)  b₋₁(r+1) − b₁(r)  a₋₁(r+1) + a₂(r−1)·b₋₂(r+1),
> C₃ = a₁(r+2)b₋₁(r+3) − b₁(r+2)a₋₁(r+3) + a₂(r+2)·b₋₂(r+4).
> ```

and then the diff-1 endgame runs verbatim, and in its shortest form: `Q₀ = 1` gives
`(T−1)G = 1 = (T−1)E`, so `G − E` is in `ker(T−1) =` constants (now degree-free,
§2), and `G(0) = 0` under membership gives `G = E`. With `h^{[-1]} | G` that reads
`h^{[-1]} | E` — impossible, since `deg h^{[-1]} = 2 > 1` and `E ≠ 0`. *(No claim
about `M(0)` is needed: the degree count alone closes it.)*

> **CONDITIONAL CLOSURE (arbitrary degree).** *diff-2 closes at arbitrary degree
> **modulo** the two scalar conditions `C₁ = C₃ = 0`.* The residual is real, not
> vacuous: an explicit tail assignment makes `C₁ ≠ 0` (verified). **This replaces
> a vague open item with two named scalars — the sharpest available formulation of
> what diff-2 still needs.**

### 4.4 Bounded backstop, both engines, and a much higher cap

Verifier `§7`, `r = 0`, `κ = 1`, integer systems, full cascade `Q_m = δ_{m0}`
(`m ∈ [−6,6]`) + membership. Before any load-bearing call, the verifier
unconditionally tests denominator-safe serialization: coefficient-only rational
denominators are cleared through the expanded numerator, denominators involving
solver unknowns are rejected, and the serialized body contains neither `/` nor
`**`. It also parses the complete comment-stripped msolve record exactly: only
`[-1]` is empty, and only `[1,N,-1,[]]` for a positive integer `N` is recognized
as nonempty; malformed, prefixed, suffixed, concatenated, garbage, and unknown
records yield no verdict. Live controls then check a known unit ideal, the
**complex-only** ideal `(x²+1)`, and a real-rooted feasible ideal.

| cap `d` | vars | diff-1 | diff-2 | engine / tier |
|---|---|---|---|---|
| `d = 1` | 22 | EMPTY | EMPTY | sympy, rerun on this checkout |
| `d = 2` | 33 | EMPTY | EMPTY | sympy, rerun on this checkout |
| `d = 3` | 44 | EMPTY | EMPTY | sympy, rerun on this checkout |
| `d = 4` | 55 | EMPTY | EMPTY | **both** engines, historical `HEAVY=1`; not rerun |
| `d = 6` | 77 | EMPTY | EMPTY | historical msolve char-0; not rerun |
| `d = 8` | 99 | EMPTY | EMPTY | historical msolve char-0; not rerun |
| `d = 10` | 121 | EMPTY | EMPTY | historical msolve char-0; not rerun |
| `d = 12` | 143 | EMPTY | EMPTY | historical msolve char-0; not rerun |

The `d = 4..12` rows are preserved as historical evidence, not current-run results;
live msolve and the `HEAVY` tier were unavailable or skipped in this audit. The
omitted historical timings were environment-specific and are not evidence; the
record is the exact engine, cap, algebraic verdict, and rerun status.

**Historically this triples the corpus's bounded cap**
([`band3-sectors.md`](band3-sectors.md) §5.2 reached `d = 4` under `HEAVY`).
`[-1]` is msolve's *empty over `C̄`* verdict, and the complex-only parser control
is what licenses that interpretation when the live msolve leg is executed.

**Historical engine finding, recorded for the corpus.** On these rigid sectors
msolve was **two orders of magnitude faster** than sympy at `d = 4` (0.5 s vs
~50 s) — the same engine-selection lesson as the A\*-band3 situation, where
sympy stalls and msolve flies ([`band3-sectors.md`](band3-sectors.md) §6).
House rule 6 still applies: benchmark the available exact engines rather than
assuming one universally dominates.

Non-vacuity at the same encoding: dropping `Q₀ = 1` leaves an **explicit point**
(`a₃, b₂` as given, all else `0`, commutator computed `= 0`). *(Recorded negative
engine result: the `Q₀`-dropped ideal is positive-dimensional and msolve does
**not** finish it at `d = 6` within 240 s — which is exactly why the non-vacuity
tier here is an explicit point and not a Gröbner verdict.)*

---

## 5. What this changes in the band-3 ledger

| class | status before | **status now** |
|---|---|---|
| cube-separated / 2-separated `h` (incl. diff-3, `(E−r)²`) | "closed", with an unnoticed `c ≠ 0` assumption in the `κ = 0` chain | **closed, arbitrary degree, gap REPAIRED** |
| **diff-1** `(E−r)(E−r−1)` | closed for `κ ≠ 0`, and `κ = 0 ∧ b₁ ≠ 0` | **CLOSED at arbitrary degree** (this arm supplies the third branch `κ = 0 ∧ b₁ = 0` by Nonpositive-D Exclusion Theorem; the other two are cited from `band3-sectors.md`) |
| **diff-2** `(E−r)(E−r−2)` | `κ = 0` sub-branch open; `κ ≠ 0` branch open, restoration refuted | `κ = 0 ∧ b₁ = 0` **CLOSED by Nonpositive-D Exclusion Theorem**; `κ = 0 ∧ b₁ ≠ 0` **closed** by direct cascade; `κ ≠ 0` branch **still open**, now reduced to the two covectors `C₁, C₃`; positive cascade **exhausted** on the refined family (with `κ ≠ 0`, `λ ≠ 0`); bounded-empty `d ≤ 4` both engines |
| `κ = 0 ∧ b₁ = 0` sub-branch (all classes) | named open gap | **CLOSED at arbitrary degree by NONPOSITIVE-D EXCLUSION THEOREM** |
| A\*-band3 constant-top corner `κ₂ = 0 ∧ b₁ = 0` (NOT all `κ₂ = 0`) | not isolated | **CLOSED** (special case of NONPOSITIVE-D EXCLUSION THEOREM; the existing `κ₂ = 0` tame witness in `astar-band3.md` §6 has positive band-one `D` with `b₁ ≠ 0`, so no contradiction with Nonpositive-D Exclusion Theorem) |
| A\*-band3 general `κ₂ ≠ 0` at arbitrary degree | open | **UNCHANGED — still open** |
| general band `k ≥ 2`, `band D ≤ 0` | not previously stated | **CLOSED at arbitrary degree by NONPOSITIVE-D EXCLUSION THEOREM** |

**No Weyl pair and no counterexample is constructed; DC1/JC2 untouched.**

### 5.1 A\*-band3 `κ₂ = 0`: only the corner `b₁ = 0` is closed by Nonpositive-D Exclusion Theorem

Nonpositive-D Exclusion Theorem closes only the sub-branch `κ = 0 ∧ b₁ = 0` (equivalently, at the
constant-top A\*-band3, `κ₂ = 0 ∧ b₁ = 0`). It does **not** close all
`κ₂ = 0`. The existing `κ₂ = 0` tame witness recorded in
[`astar-band3.md`](astar-band3.md) §6 (`U = x + 2∂`, `X = U³ − ∂`, `D = U`) has
`κ₂ = b₂ = 0` but **positive band-one `D` with `b₁ ≠ 0`** (its `D` is `U`, of
band 1), so it falls outside Nonpositive-D Exclusion Theorem's `band D ≤ 0 ∧ b₁ = 0` hypothesis and
there is **no contradiction** between the tame witness and Nonpositive-D Exclusion Theorem.

---

## 6. Honest ledger

> **Note (audit 2026-07-26):** the `Q₂` tail-free node set is **exactly
> `{r−1, r, r+1, r+2}`**, of which only `{r, r+1}` are used for the displayed
> forcing (exclusivity now machine-tested over offsets `−5..5`). The two extra
> tail-free equations at `r−1` and `r+2` are **unclassified** in this memo — they
> are not exploited and no claim is made about what they add.

**Proved (arbitrary degree; degree-free proof objects in the verifier — abstract
`sympy.Function` coefficients with no cap, plus a symbolic-degree series
derivation):**
1. Engine `Q_m = [D,X]_m`, `Q₀ = (T−1)G`, `G(0) = 0` under membership, `(T^n−1)`
   kernel = constants (`§0`), each with a negative control.
2. The `band D ≤ 0` ladder collapse: `Q_k = a_k(T^k−1)b₀` and
   `Q_{k−j} = b₋ⱼ^{[k]}a_k − a_k^{[−j]}b₋ⱼ`, with **every other coefficient of `X`
   dropping out identically** (`§1`; abstract Functions, `k = 3, 4, 5`).
3. The trace rules `τ(P^{[s]}) = τ(P) + s·deg P`, `τ(PQ) = τ(P)+τ(Q)` derived at
   **symbolic degree**, and the KEY FORMULA
   `[E^{N−1}](φ^{[k]}a − a^{[−j]}φ) = lc·lc·(kn + jA)` at symbolic `n, A, k, j`
   (`§2`), each instance-validated with negative controls. As a by-product the
   same formula **upgrades `(T^n−1)P = 0 ⇒ P constant` from the bounded `deg ≤ 3`
   probe of `§0` to a degree-free fact** (take `φ = P`, `a = 1`, `k = n`, `j = 0`:
   `[E^{p−1}](P^{[n]} − P) = lc(P)·n·p ≠ 0` for `p ≥ 1`), so no step of Nonpositive-D Exclusion Theorem
   rests on a capped kernel computation.
4. **NONPOSITIVE-D EXCLUSION THEOREM**: `band X = k ≥ 2`, `band D ≤ 0`, membership, `[D,X] = 1` is
   impossible (`§3`). Hence the `κ = 0 ∧ b₁ = 0` sector is empty in every
   shifted-cube class and, as a special case, in the A\*-band3 constant-top
   corner `κ₂ = 0 ∧ b₁ = 0` (§5). The existing `κ₂ = 0` tame witness has
   positive band-one `D` (`b₁ ≠ 0`), so Nonpositive-D Exclusion Theorem does not touch it.
5. **diff-1 closure at arbitrary degree with no side condition** (§3 above + the
   corpus's `κ ≠ 0` and `κ = 0, b₁ ≠ 0` results).
6. **diff-2, new rung**: on the surviving branch `b₁(r+1) = 0` with `κ ≠ 0` and
   `a₂(r−1) ≠ 0`, the two analyzed clean `Q₂` nodes (at `r` and `r+1`) force
   `(E−r−1) | a₁` (`§6`, degree-free). No universal "exactly one condition"
   claim: the two extra tail-free nodes are unclassified.
7. **diff-2 `G`-split**: `G − R` is `h^{[-1]}`-divisible term-by-term, and
   `h^{[-1]} | G ⟺ C₁ = C₃ = 0`; given that, diff-2 closes at arbitrary degree
   (`§7`, degree-free).

**Refuted (machine-checked, explicit exact witnesses; nonzero-branch scope):**
- **The mission's Target-B route (with `κ ≠ 0` and `λ ≠ 0`).** The deeper rungs
  `Q₂, Q₁` do **not** close the diff-2 surviving branch: an explicit family at
  symbolic `(r, κ, λ, c)`, interpreted as a point of the intended surviving
  branch under both `κ ≠ 0` and `λ ≠ 0` (concrete specialization `κ = λ = 1`
  machine-checked), solves the **entire positive cascade** `Q₆ … Q₁ = 0` with
  `a₂(r−1) ≠ 0`. *(Audit phrasing note: on this witness `Q₆` and `Q₁` are
  **vacuous** rather than solved — `Q₆`'s only pair dies on the gauge `b₃ = 0`,
  and `Q₁`'s four pairs die term-by-term on `b₀ = 0`, constant `a₀` and the
  zero tail. The real content is `Q₅…Q₂`.)*
- **The specific `a₁ = −3` extension of [`band3-sectors.md`](band3-sectors.md)
  §5.1's witness.** The published §5.1 family is correct at its stated `Q₄, Q₃`
  scope. What is new here is that it does not extend past `Q₃` with the
  specific `a₁ = −3` (which violates the new `Q₂` forcing `a₁(r+1) = 0`); the
  member that does extend is `a₁ = −3λ²(E−r−1)²`. This is a non-extension of a
  specific choice, not a shrinkage of the §5.1 family at its published scope.

**Bounded / finite evidence (exact scope stated):**
- **Instance sweep of Nonpositive-D Exclusion Theorem's conclusion** at `k = 3`: over 105 combinations
  `(j ∈ {1,2,3}) × (deg of the b₋ⱼ cofactor ∈ 0..4) × (deg a₃ ∈ 0..6)` with
  membership imposed, the rung's `E^{N−1}` coefficient matches the key formula and
  is nonzero, so the rung is never `0` and never `1`. Paired with the negative
  control that **without** membership the rung *can* equal `1`.
- `κ = 0 ∧ b₁ = 0` full-cascade + `Q₀ = 1` + membership emptiness at cap `d ≤ 2`
  for six tops (exact SymPy over `ℚ`, unit ideal), as a cross-check of Nonpositive-D Exclusion Theorem.
- diff-1 / diff-2 full sector at `r = 0, κ = 1`: `d ≤ 3` rerun on this checkout
  (SymPy); the `d = 4` **both-engine** result (msolve char-0 `[-1]`, SymPy unit
  ideal) is historical `HEAVY=1` evidence and was not rerun in this audit.
- **The `§1` `Q₅` computation with `h = (E−1)(E−4)` and `b₂` degree cap 5 is a
  BOUNDED CONCRETE CROSS-CHECK only** — it does not establish an
  arbitrary-`h`/degree necessity proof or classify the whole solution space of
  `Q₅`; the verifier's check label is relabelled accordingly.
- Non-vacuity everywhere by **explicit points**, not normal forms.

**Open / NOT claimed:**
1. **diff-2 at arbitrary degree with `κ ≠ 0`, branch `b₁(r+1) = 0`** — reduced to
   the two tail covectors `C₁ = C₃ = 0` (§4.3). Deriving those from `Q₋₁, Q₋₂`
   (or refuting them) is now the whole of the residual.
2. A\*-band3 general `κ₂ ≠ 0` at arbitrary degree; composite tame-word escape.
3. Imbalanced coprime walls; general-`k` negative tail with `band D > 0`; **W2**;
   radical forcing at coupling widths `k = 4, 5`.
4. Nonpositive-D Exclusion Theorem's rung collapse is machine-derived at `k = 3, 4, 5`; the `k`-uniform
   statement is the written argument (the `s = 0` term dies because `b₀` is
   constant, for every `k`), not a symbolic-`k` machine identity.

No Weyl pair, no counterexample; DC1/JC2 untouched.

---

## 7. Verification

```sh
uv run --with sympy python research/dc1-program/verify_shifted_cube_completion.py
HEAVY=1 uv run --with sympy python research/dc1-program/verify_shifted_cube_completion.py
```

Historical run counts (default `100/100`, `HEAVY=1` `112/112`) predate the
2026-07-26 audit edits and are **superseded** unless reproduced on this
checkout; no new run-count claim is asserted here.

`§0` engine; `§1` the `band D ≤ 0` collapse (degree-free, `k = 3,4,5`); `§2` the
trace functional at symbolic degree + key formula; `§3` NONPOSITIVE-D EXCLUSION THEOREM; `§4` controls
(band 1 survives; membership load-bearing; bounded cross-checks; explicit-point
non-vacuity); `§5` what it repairs; `§6` diff-2's new rung, the cut witness and the
refined witness; `§7` the two covectors + bounded emptiness (both engines).
Runtime is environment-dependent. The msolve evidence uses denominator-guarded
polynomial serialization and exact whole-record verdict parsing, so unknown solver
output cannot become PASS. The final banner distinguishes *all checks passed* from
*all executed checks passed; optional checks skipped*, and lists every skip with
its reason.
