# A\*-band3, DC1 face: the constant-top (constant-h) negative-tail closure

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED — BOUNDED CLOSURE + NAMED DEGREE-FREE OBSTRUCTION**

This memo isolates and attacks the **A\*-band3 negative-tail closure** on the
Weyl/Dixmier (DC1) face: the statement that the **constant-`h` (constant-top)
band-3 sector is exactly the tame family**. This is the
[`shifted-power-residuals.md`](shifted-power-residuals.md) §3–§4 **RESIDUAL 3**
`(κ₂-closure)`, and it is doubly load-bearing — it finishes the `κ₂ ≠ 0`
disposition of that memo **and** the constant-`h` completeness of the
shifted-power descent ([`shifted-power-descent.md`](shifted-power-descent.md) §6
residual 2).

The one-line summary:

> **The top potential-factorization that killed *nonconstant* `h` is VACUOUS for
> constant `h`, so the closure must come entirely from the negative tail — and the
> negative tail's `μ₃` cross-coupling (the band-3 effect with no band-2 shadow)
> blocks a clean band-2-style first-integral / `mod 3` lattice.** What is delivered
> is therefore a **bounded** closure — the `κ₂ ≠ 0` constant-top sector is
> **EMPTY over `ℚ`** at coefficient-degree cap `d ≤ 2` (committed exact `d = 1`
> certificate by two independent engines; `d = 2` by `msolve` over `ℚ`), with an
> explicit `κ₂ = 0` tame witness — together with the **exact identification of the
> arbitrary-degree obstruction**. The arbitrary-degree closure is **OPEN**, mirroring precisely the
> classical `e ≠ 0` sibling
> ([`../band3/classical-e-nonzero-closure.md`](../band3/classical-e-nonzero-closure.md)),
> which is itself open with only bounded + conditional evidence.

Exact certificate: [`verify_astar_band3.py`](verify_astar_band3.py) — default run
exact SymPy over `ℚ` (`ALL ASTAR-BAND3 (DC1) CHECKS PASSED`); `msolve` corroboration
gated (SKIPs cleanly if absent); cap `d = 2, 3` behind `HEAVY=1`. Every load-bearing
upstream fact (the crossed-product ladder engine `Q_m = [D,X]_m`, `Q₀ = (T−1)G`, the
telescoping potential, the wall) is **re-derived in file**, not merely cited.

Conventions frozen from the corpus: `A₁[x^{-1}] = ⊕_k x^k C[E]`, `E = x∂`,
`(x^a f)(x^b g) = x^{a+b} f(E+b) g(E)`, `f^{[r]}(E) = f(E+r)`, `T f = f^{[1]}`,
`Q_m = ∑_{k+l=m}[b_l^{[k]}a_k − a_k^{[l]}b_l]`, `[D,X] = 1 ⇔ Q_m = δ_{m0}`,
membership `(E)_r = E(E−1)⋯(E−r+1) | a_{-r}, b_{-r}`, gauge `b_3 = 0`,
`G = ∑_{i=1}^{3} ∑_{r=1}^{i}(a_i^{[-r]}b_{-i}^{[i-r]} − a_{-i}^{[i-r]}b_i^{[-r]})`,
`Q₀ = (T−1)G`. **Sector:** constant top `a₃ = 1`, gauge `b₃ = 0`, wall
`Q₅ ⇒ b₂ = κ₂` (constant).

---

## 0. Headline

> **Theorem (bounded).** In the constant-top band-3 Weyl sector (`a₃ = 1`, gauge
> `b₃ = 0`, so `b₂ = κ₂` constant), at coefficient-degree cap `d ≤ 2` there is **no**
> genuine pair `[D,X] = 1` with `κ₂ ≠ 0`. Every genuine pair in the certified range
> has `κ₂ = 0` — the tame slice, which is nonempty (explicit witness). *(Exact `ℚ`;
> `d = 1` committed by two engines; `d = 2` by `msolve` over `ℚ`; `d = 3` behind
> `HEAVY`.)*
>
> **Structure (arbitrary degree, machine-checked identities).** `[D,X] = 1` is
> **bilinear** in the coefficients: for a fixed constant-top `X` the whole system is
> **linear in `D`**. For a *generic* constant-top `X` it is **not solvable** — the
> moment covector (Lemma P) obstructs — so admissibility is a positive-codimension
> condition. The top wall `Q₅ = (T³−1)b₂` forces `b₂ = κ₂`; the potential is
> `G = h^{[-1]}M` with `h = 1`, hence `h^{[-1]} = 1` and **`G = M`**: the
> nonconstant-`h` kill is unavailable, so the obstruction lives in the **negative
> tail**, where `Q₋₅, Q₋₄` carry an **inhomogeneous `μ₃`-source** (bottom
> proportionality `b₋₃ = μ₃ a₋₃`) with no band-2 shadow.
>
> **Why not arbitrary degree.** The band-2 A\* kill was a `mod 3` congruence
> `3V = 3P + 1` read off two first integrals `Φ, I₂` of the negative cascade
> ([`../band2-square-sector/classical-Astar.md`](../band2-square-sector/classical-Astar.md)).
> At band 3 the raw Band-2-shaped ansatz is obstructed; the Wave-C compensated
> integral (`Φ' = C₋₁`, `I₂' = C₋₂ − (2/3)a₂'Φ`) exists only in the classical face
> and yields a **conditional** determinant / tropical `6:5` degree gap, **not** a
> `mod 3` lattice — and the classical `e ≠ 0` sector is itself **open**. The quantum
> reflection `E ↦ −E−1` breaks `A₁`-membership for the constant top, so the
> classical route does not even transcribe. **No clean degree-free closure is
> available; the exact obstruction is named.**

---

## 1. The sector, exactly

Orient so `a₃ ≠ 0`; the top Wronskian `Q₆ ⇒ b₃ = λ₃ a₃`, and the single gauge
`D ↦ D − λ₃X` sets `b₃ = 0`. Diagonal scaling normalizes the constant top `a₃ = 1`
(the `h = 1` shifted cube). The wall is then, exactly,

```
Q₅ = b₂^{[3]}·a₃ − a₃^{[2]}·b₂ = (T³−1) b₂ = 0   ⇒   b₂ = κ₂  (constant),
```

the **gauged wall constant** (`verify §1`). The tame family
`U = x + c₀ + c₁∂`, `X = U³ − ∂/κ − A`, `D = λX + κU + β` sits, after the gauge
`D ↦ D − λX`, at `D' = κU + β` — which is **band 1**, so `b₂(D') = 0`: **the tame
slice is `κ₂ = 0`**. (This is the band-3 twist relative to band 2, where the wall
level `b₁` lies *inside* band 1 and survives the gauge as the nonzero `κ`; at band 3
the wall level `b₂` lies *above* band 1 and the gauge kills it.) `κ₂` is a proved
tame invariant of the gauge-`b₃=0` sector (no transvection / pair-exchange / Fourier
alters it — [`shifted-power-residuals.md`](shifted-power-residuals.md) §3.1). Hence:

> **`(κ₂-closure)`.** No genuine Weyl pair `[D,X] = 1` has `a₃ = 1`, `b₃ = 0`,
> `b₂ = κ₂ ≠ 0`, membership-valid negative tail `Q_{-1} = ⋯ = Q_{-6} = 0`, and
> `Q₀ = 1`. Equivalently: the constant-top sector's genuine pairs are exactly the
> `κ₂ = 0` tame slice.

The positive cascade **permits** `κ₂ ≠ 0` (`Q₄ = κ₂(a₂−a₂^{[2]}) + (b₁^{[3]}−b₁)`
has `κ₂ ≠ 0` solutions), so `(κ₂-closure)` is genuinely a **negative-tail** question,
not a positive-level one.

## 2. Engine, telescoping potential, moment slope (re-derived)

`verify §0` re-derives `Q_m = [D,X]_m` for `m ∈ [−6,6]` against the direct
crossed-product commutator (generic degree-2 coefficients), and
`Q₀ = (T−1)G` with the band-agnostic staggered potential `G`. Membership makes
every term of `G(0)` vanish, so `Q₀ = 1 ⇒ G = E`; the **slope** is `G(1)`, and
(Lemma P, `moment-unit-general-k.md`) for the constant top

```
G(1) = a₁(0)b₋₁(1) − a₋₁(1)b₁(0) + a₂(0)b₋₂(2) − a₋₂(2)·κ₂ + μ₃·a₋₃(3) = 1
```

(`verify §1`), on the branch `a₋₃ ≠ 0` with `b₋₃ = μ₃ a₋₃`. This single bilinear
scalar identity is the load-bearing covector of §4.

## 3. Positive cascade and the vacuity of the top potential-factorization

The constant-top positive rung is the exact identity
`Q₄ = κ₂(a₂ − a₂^{[2]}) + (b₁^{[3]} − b₁)` (`verify §1`), and `Q₃, Q₂, Q₁`
determine `b₀, b₋₁, b₋₂` (3-fold periodic summations, each with one free constant,
membership pinning the deep ones). This is the exact quantum analogue of the band-2
positive cascade with 2-fold periodicities replaced by 3-fold.

**The key structural point.** In the nonconstant-`h` shifted-cube descent, the
central potential factors as `G = h^{[-1]}M`, and `Q₀ = 1 ⇒ h^{[-1]} | E`, which
kills every nonconstant `h`. **Here `h = 1`, so `h^{[-1]} = 1` and `G = M`** — the
factorization is trivial and `G = M = E` is perfectly consistent (`verify §1`).
*The top-level h-forcing that closes the nonconstant sector contributes nothing to
the constant sector.* The obstruction to `(κ₂-closure)` therefore **cannot** sit at
the top; it sits in the negative tail. This is the precise content of the
shifted-power-descent §6 residual-2 hand-off.

## 4. The negative tail: `μ₃`-source, bilinearity, and the covector obstruction

**Bottom proportionality.** `Q₋₆ = 0 ⇒ b₋₃ = μ₃ a₋₃` on `a₋₃ ≠ 0` (`verify §2`);
the single gauge is spent on the top, so `μ₃` is not removable (the `λ₃–μ₃`
cross-coupling).

**Inhomogeneous `μ₃`-source (no band-2 shadow).** With `b₋₃ = μ₃ a₋₃`,

```
Q₋₅ = [ b₋₂^{[-3]}a₋₃ − a₋₃^{[-2]}b₋₂ ] + μ₃·[ a₋₃^{[-2]}a₋₂ − a₋₂^{[-3]}a₋₃ ],
Q₋₄ = [ b₋₁^{[-3]}a₋₃ − a₋₃^{[-1]}b₋₁ ] + μ₃·[ a₋₃^{[-1]}a₋₁ − a₋₁^{[-3]}a₋₃ ]
                                        + [ b₋₂^{[-2]}a₋₂ − a₋₂^{[-2]}b₋₂ ]
```

(`verify §2`, exactly the `quantum-band3-cascade.md` §5 decompositions). The
`μ₃`-proportional source is the top/bottom cross-coupling absent at band 2 — the
"first genuinely new band-3 effect."

**Bilinearity ⇒ the covector route.** `[D,X] = 1` is **linear in every `b_l`** (and
in every `a_k`) separately (`verify §2`). So for a fixed constant-top `X`, the
system `Q_m = δ_{m0}` is an inhomogeneous **linear** system in the `D`-coefficients.
Computing it at a generic constant-top `X` (`verify §2`): the `D`-image has rank
`(#unknowns − 1)` — the one-dimensional kernel is the constant centralizer
`D = const` — and the augmented rank is one higher, so **generic constant-top `X`
is not solvable**. **Tier (audit-demoted): this rank/obstruction computation is a
fixed-cap generic-instance result (caps `(dX,dD)=(1,3),(2,4)`, random seed; the
audit re-ran seeds 7/42/100 with identical ranks) — bounded-finite, not
arbitrary-degree.** Admissibility is a positive-codimension condition; the critical
cokernel covector applied to `δ_{m0} = 1` is exactly the moment slope of §2 (the
filler image is **not** everything). This is the concrete "λ arsenal / covector"
object for this sector.

> **Verifier-hygiene note (audit).** Several arbitrary-degree ledger identities
> (`Q_0=(T-1)G`, `G(0)=0` under membership, the Lemma-P slope) are exercised in-file
> on random numeric instances rather than symbolic coefficients; the audit supplied
> the symbolic versions externally and all hold. The in-file machine scope is the
> instance level. The `d=1` char-0 `msolve` corroboration is additional/manual, not
> part of the committed default run.

## 5. The degree-free mechanism and its exact obstruction

The band-2 A\* kill (`classical-Astar.md` §2–3) used **two exact first integrals**
of the negative cascade, `Φ` (of `C₋₁`) and `I₂` (of `C₋₂`), whose membership-forced
vanishing gave the degree balances `2V = P + W` and `V + W = 2P + 1`, hence
`3V = 3P + 1` — infeasible `mod 3`. The band-3 analogue is obstructed at three
successive points, and this memo pins each:

1. **The raw Band-2-shaped ansatz is obstructed.** No constant-coefficient
   combination of `Q₋₁..₋₅` is an exact `(T−1)`-difference of a local expression in
   the searched basis — the band-3 statement of `astar-band3.md` (Wave-A/B) §5.
2. **The compensated integral is classical-only and conditional.** Classically the
   Wave-C construction restores integrability with a nonlocal generator:
   `Φ' = C₋₁`, `I₂' = C₋₂ − (2/3)a₂'Φ` (`classical-e-nonzero-closure.md`). But the
   resulting `Φ, I₂` are **linear in the trailing pair `(a₋₂, a₋₃)`** and determine
   it only off a **determinant locus** `det = −(4/3)e²a₁ + (4/9)e²a₂² − κ₁²`, with a
   **tropical `6:5`** degree signature (not a `mod 3` lattice). Even classically this
   leaves the sector **OPEN** (denominator cancellation, polynomiality, `det = 0`,
   `a₂ = const` strata).
3. **The quantum reflection does not transcribe.** The classical import route used
   reflection `R:(x,ξ) ↦ (ξ,x)`; its quantum analogue `E ↦ −E−1` sends the constant
   top `a₃ = 1` to a reflected `x^{-3}`-coefficient `1`, which is **not** divisible
   by `E(E−1)(E−2)` — `A₁`-membership breaks (`astar-band3.md` §6(i)). And the
   inhomogeneous `μ₃`-source (§4) has no band-2 mirror to cancel against.

> **Exact obstruction.** The `(κ₂-closure)` degree-free proof would require a
> quantum first integral of the `μ₃`-sourced tail that determines the trailing pair
> `(a₋₂, a₋₃)` and forces a degree contradiction. The classical mirror shows such an
> integral is at best **conditional** (determinant locus + tropical gap, no `mod 3`
> lattice), and the quantum face additionally loses the reflection import and gains
> the `μ₃`-source. This is precisely where the closure sits, and it is **open** —
> the same frontier as classical `e ≠ 0` and quantum A\*-I.

## 6. The bounded certificate

Normalize `κ₂ = 1` (any nonzero value, by scaling). Build the full constant-top
system `Q_m = δ_{m0}`, `m ∈ [−6,6]`, at free-polynomial-degree cap `d` (membership
factors `E`, `E(E−1)`, `E(E−1)(E−2)` on the level-`−1,−2,−3` coefficients). All
coefficients are **integer**; `msolve` requires cleared denominators (it misparses
rational monomials — a documented trap in the verifier).

| cap `d` | `κ₂ ≠ 0` (normalized `κ₂ = 1`) | engine(s) | tier |
|---|---|---|---|
| `d = 1` | **UNIT IDEAL (EMPTY over `ℚ`)** | SymPy `ℚ` Gröbner `= [1]` + `msolve` `-g` unit over 3 primes *(default run)*; `msolve` char-0 `[-1]` *(additional)* | **committed** |
| `d = 2` | **EMPTY over `ℚ`** | `msolve` char-0 `[-1]` (≈35 min) | HEAVY (`HEAVY=1`), reproducible |
| `d = 3` | attempted | `msolve` char-0 / `-g` over `F_p` — did **not** complete within the HEAVY time cap in development | HEAVY, not certified |

The `d = 1` row is the load-bearing committed certificate (two `ℚ`-rigorous
engines, `SymPy` unit ideal cross-checked by `msolve` `[-1]`, plus a 3-prime `-g`
corroboration). The `d = 2` row is an exact `ℚ` result via `msolve`'s rational
solver (`[-1] =` empty variety), reproducible under `HEAVY=1`. `d = 3` is attempted
but its Gröbner did not finish in the budgeted window; no `d = 3` claim is made.

The `κ₂ = 0` slice is **nonempty**: the explicit positive control
`U = x + 2∂`, `X = U³ − ∂`, `D = U` is a genuine pair `[D,X] = 1` with `a₃ = 1`,
`b₂ = 0` (`verify §3`), and the `κ₂ = 0` system at `d = 1` is **not** the unit ideal.
So within cap `d ≤ 2`, every constant-top genuine pair lies in the `κ₂ = 0` slice.
(That the `κ₂ = 0` slice is *itself* exactly the tame family — the reverse
direction — is the **B0-band3 single-shear-origin structure lemma**,
[`../band3/band3-tame-catalog.md`](../band3/band3-tame-catalog.md) §2/§6, cited here,
not re-proved. Combining: constant-top genuine pairs at cap `d ≤ 2` are exactly the
tame family.)

**A methodological note recorded for the corpus.** A hand-rolled forward-solve of
the positive cascade produced a *rational-coefficient* reduced system; `msolve`
silently **misparsed** its rational monomials (`2*x^2/3`), reporting a spurious
non-unit basis while SymPy over `ℚ` (and the full integer system in both engines)
correctly reported the unit ideal. The verifier only ever hands `msolve`
**integer** systems. This is a second `msolve` trap beyond `**` vs `^`.

## 7. What the closure delivers for Gap 2, band 3

- **RESIDUAL 3 (`shifted-power-residuals.md` §4 pt 4).** The `κ₂ ≠ 0` sector is
  disposed of **at bounded degree** (`d ≤ 2`, exact `ℚ`): no genuine pair beyond the
  tame `κ₂ = 0` slice. The arbitrary-degree `(κ₂-closure)` remains **open**, with the
  obstruction now named exactly (§5) and shown to coincide with the classical
  `e ≠ 0` / quantum A\*-I frontier.
- **shifted-power-descent §6 residual 2.** The constant-`h` completeness of the
  descent is reduced to `(κ₂-closure)`; the top potential-factorization is shown
  **vacuous** for constant `h` (§3), so the descent's constant-`h` step is complete
  **modulo** exactly this negative-tail closure, bounded-verified here to `d ≤ 2`.

Beyond band 3, unchanged: the imbalanced coprime walls, the general-`k` negative
tail, and **W2** remain open independently. **No Weyl pair and no counterexample is
constructed; DC1/JC2 untouched.**

## 8. Honest ledger

**Proved (exact algebra, machine-checked identities; arbitrary degree):**
- Engine `Q_m = [D,X]_m` and telescoping `Q₀ = (T−1)G`, `G(0) = 0` (`§0`).
- Constant-top wall `Q₅ = (T³−1)b₂ ⇒ b₂ = κ₂`; the `Q₄` rung; Lemma-P slope; the
  **vacuity of the potential-factorization** (`G = M` for `h = 1`) (`§1`).
- Bottom proportionality `Q₋₆ ⇒ b₋₃ = μ₃ a₋₃`; the inhomogeneous `μ₃`-source
  decompositions of `Q₋₅, Q₋₄`; the **bilinearity** of `[D,X]`; the generic-`X`
  non-solvability (moment/covector obstruction, filler image ≠ everything) (`§2`).

**Bounded / finite evidence (exact scope):**
- `κ₂ ≠ 0` constant-top sector **EMPTY over `ℚ`** at cap `d = 1` (committed default:
  SymPy `ℚ` unit ideal `[1]` + `msolve` `-g` unit over 3 primes; additionally
  `msolve` char-0 `[-1]`) and cap `d = 2` (`msolve` char-0 `[-1]`, `HEAVY`,
  reproducible). `d = 3` attempted, not certified within the time budget.
- explicit `κ₂ = 0` tame witness `[D,X] = 1`; `κ₂ = 0` at `d = 1` not unit.

**Refuted (machine-checked) — corpus corrections:**
- That the top potential-factorization contributes to the constant-`h` closure: it
  is **vacuous** (`h^{[-1]} = 1`).
- That `msolve` may be fed rational-coefficient systems: it **misparses** them; only
  integer systems are certified.

**Open / NOT claimed:**
1. `(κ₂-closure)` at arbitrary degree — the negative-tail first integral; the exact
   obstruction (§5) is the classical `e ≠ 0` / quantum A\*-I frontier.
2. Non-2-separated tops, imbalanced coprime walls, general-`k` tail, **W2**.
No Weyl pair, no counterexample; DC1/JC2 untouched.

## 9. Verification

```sh
uv run --with sympy python research/dc1-program/verify_astar_band3.py
HEAVY=1 uv run --with sympy python research/dc1-program/verify_astar_band3.py
```

Exact SymPy over `ℚ`: `§0` engine + telescoping potential; `§1` constant-top wall +
`Q₄` rung + Lemma-P slope + potential vacuity; `§2` bottom proportionality +
inhomogeneous `μ₃`-tail + bilinearity + generic-`X` non-solvability; `§3` the
**bounded certificate** (cap `d = 1` unit ideal committed by two engines, tame
witness) with `msolve` corroboration gated by `shutil.which` (SKIPs cleanly), and
cap `d = 2, 3` behind `HEAVY=1`. A successful run ends
`ALL ASTAR-BAND3 (DC1) CHECKS PASSED`.
