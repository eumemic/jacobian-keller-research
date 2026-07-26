# A\*-band3, DC1 face: the constant-top (constant-h) negative-tail closure

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED — BOUNDED CLOSURE + NAMED DEGREE-FREE OBSTRUCTION**

This memo isolates and attacks the **A\*-band3 negative-tail closure** on the
Weyl/Dixmier (DC1) face: whether the **constant-`h` (constant-top)** band-3 sector
reduces to the known tame family. This is the
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
> is therefore a **bounded** closure — the normalized `κ₂ = 1` fiber is
> **EMPTY after base change to an algebraic closure** at coefficient-degree cap `d ≤ 2`
> (an exact `ℚ` certificate at `d = 1`, finite-field corroboration, and a recorded exact
> characteristic-zero `msolve` result at `d = 2`). Since the scaling below sends
> `κ₂ ↦ s⁵κ₂`, every nonzero `κ₂` over an algebraically closed characteristic-zero field
> can be normalized by choosing `s⁵ = κ₂⁻¹`; this yields bounded emptiness of that sector
> over the algebraic closure, not a rational-coordinate fiber-by-fiber computation. There is an
> explicit `κ₂ = 0` tame witness — together with the **exact identification of the
> arbitrary-degree obstruction**. The arbitrary-degree closure is **OPEN**, mirroring precisely the
> classical `e ≠ 0` sibling
> ([`../band3/classical-e-nonzero-closure.md`](../band3/classical-e-nonzero-closure.md)),
> which is itself open with only bounded + conditional evidence.

Exact certificate: [`verify_astar_band3.py`](verify_astar_band3.py) — default exact
SymPy over `ℚ`; `msolve` corroboration is optional and cap `d = 2, 3` is behind
`HEAVY=1`. Runtime is environment-dependent. The final summary reports either all checks
passed with no skips or all executed checks passed with explicit skips. Every load-bearing
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
> `b₃ = 0`, so `b₂ = κ₂` constant), the normalized fiber `κ₂ = 1` is empty after base
> change to an algebraic closure at coefficient-degree cap `d ≤ 2`. Under the diagonal
> variable/operator scaling of §6, `κ₂ ↦ s⁵κ₂`; hence over an algebraically closed
> characteristic-zero field any `κ₂ ≠ 0` can be normalized by choosing
> `s⁵ = κ₂⁻¹`, and normalized-fiber emptiness transports to every nonzero sector over
> that algebraic closure. This is not generally a rational-coordinate equivalence because
> the fifth root need not lie in `ℚ`. The `κ₂ = 0` slice is nonempty and contains the
> explicit tame witness. *(At normalized `κ₂=1`, `d = 1` has an exact SymPy certificate
> over `ℚ` plus finite-field corroboration; the recorded `d = 2` result is an exact
> characteristic-zero `msolve` computation; `d = 3` was attempted but is not certified.)*
>
> **Structure.** `[D,X] = 1` is **bilinear** in the coefficients: for a fixed
> constant-top `X` the whole system is **linear in `D`** (an arbitrary-degree
> identity). At the tested fixed degree caps, generic sampled constant-top `X` is
> **not solvable**, so the sampled admissible locus has positive codimension. The rank
> calculation proves existence of a cokernel functional at those samples; it does not
> identify that functional with the Lemma-P moment slope or prove generic nonsolvability
> at arbitrary degree. The top wall `Q₅ = (T³−1)b₂` forces `b₂ = κ₂`; the potential is
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
> `mod 3` lattice — and the classical `e ≠ 0` sector is itself **open**. Genuine
> Weyl Fourier preserves `A₁` membership and reverses support, but exits the
> displayed top-wall/gauge chart, so its chart parameter is not directly
> comparable; formal coefficient reflection is not Fourier. **No clean degree-free closure is
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
`D ↦ D − λX`, at `D' = κU + β` — which is **band 1**, so `b₂(D') = 0`: the
**displayed tame family lies in `κ₂ = 0`**. (This is the band-3 twist relative to
band 2, where the wall level `b₁` lies *inside* band 1 and survives the gauge as the
nonzero `κ`; at band 3 the wall level `b₂` lies *above* band 1 and the displayed gauge
kills it.) The existing corpus verifies `κ₂` only under the displayed direct tame
generators; invariance under arbitrary composite tame words remains open
([`shifted-power-residuals.md`](shifted-power-residuals.md) §3.1,
[`shifted-power-descent.md`](shifted-power-descent.md) §6).

> **`(κ₂-closure)` (open at arbitrary degree).** The target exclusion is that no
> genuine Weyl pair `[D,X] = 1` has `a₃ = 1`, `b₃ = 0`, `b₂ = κ₂ ≠ 0`, a
> membership-valid negative tail, and `Q₀ = 1`. This memo proves that exclusion only
> at the stated coefficient caps. Even a full exclusion would show `κ₂=0`; it would
> not by itself classify every `κ₂=0` pair as tame.

The positive cascade **permits** `κ₂ ≠ 0` (`Q₄ = κ₂(a₂−a₂^{[2]}) + (b₁^{[3]}−b₁)`
has `κ₂ ≠ 0` solutions), so `(κ₂-closure)` is genuinely a **negative-tail** question,
not a positive-level one.

> **Cross-corpus (audit 2026-07-26): what the Nonpositive-D Exclusion Theorem
> of [`shifted-cube-completion.md`](shifted-cube-completion.md) closes here, and
> what it does NOT.**
> The Nonpositive-D Exclusion Theorem closes the constant-top corner
> **`κ₂ = 0 ∧ b₁ = 0`** (as a special case of "band `D ≤ 0` with band `X = k ≥ 2`
> + membership + `[D,X] = 1` is impossible"). It does **not** close all `κ₂ = 0`:
> the explicit `κ₂ = 0` tame witness of §6 (`U = x + 2∂`, `X = U³ − ∂`, `D = U`)
> has `κ₂ = 0` **but positive band-one `D`** — its `D = U` is band 1, so
> `b₁(D) ≠ 0`, and the witness falls outside the theorem's
> `band D ≤ 0 ∧ b₁ = 0` hypothesis. There is no contradiction between the tame
> witness and the Nonpositive-D Exclusion Theorem. The general
> `κ₂ ≠ 0` question at arbitrary degree therefore **remains open**.

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
is not solvable**. **Tier (audit-demoted): the in-file rank/obstruction computation is a
fixed-cap generic-instance result (`(dX,dD)=(1,3)`, seed 11; an external audit
re-ran seeds 7/42/100 at that cap with identical ranks) — bounded-finite, not
arbitrary-degree.** At those samples, admissibility is a positive-codimension condition and the filler
image is **not** everything. The rank computation supplies some nonzero cokernel
functional, but the verifier does not construct it or prove proportionality to the
moment-slope functional of §2.

> **Verifier-hygiene note (audit).** Several arbitrary-degree ledger identities
> (`Q_0=(T-1)G`, `G(0)=0` under membership, the Lemma-P slope) are exercised in-file
> on random numeric instances rather than symbolic coefficients; the audit supplied
> the symbolic versions externally and all hold. The in-file machine scope is the
> instance level. In this older verifier, the `d=1` char-0 `msolve` corroboration is
> additional/manual, not part of the committed default run. The newer
> [`verify_band3_sectors.py`](verify_band3_sectors.py) conditionally executes its msolve legs,
> reports solver identity, and offers `--require-msolve`; when msolve is unavailable those
> newer legs are explicitly `NOT RERUN/SKIPPED`. This does not retroactively change the
> execution status recorded here.

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
3. **Genuine Fourier reverses the chart rather than breaking membership.** Weyl
   Fourier sends `a_i` to `a'_{−i}=(-1)^i(E)_i a_i(-E−1)`; in particular, the
   constant top `a₃=1` maps to `−(E)_3`, so the reflected negative-band coefficient
   satisfies `A₁` membership. The move reverses support and exits the displayed
   top-wall/gauge chart, however, so `κ₂` is not directly comparable there. A bare
   coefficient reflection `E↦−E−1` is not Weyl Fourier. The inhomogeneous
   `μ₃`-source (§4) still has no band-2 mirror to cancel against.

> **Exact obstruction.** The `(κ₂-closure)` degree-free proof would require a
> quantum first integral of the `μ₃`-sourced tail that determines the trailing pair
> `(a₋₂, a₋₃)` and forces a degree contradiction. The classical mirror shows such an
> integral is at best **conditional** (determinant locus + tropical gap, no `mod 3`
> lattice), while genuine Fourier reverses/exits the displayed chart rather than
> importing that closure into it. Composite tame escape and arbitrary-degree tail
> closure remain open, alongside the `μ₃`-source frontier — the same frontier as
> classical `e ≠ 0` and quantum A\*-I.

## 6. The bounded certificate

Work on the normalized fiber `κ₂ = 1`. Under

```text
x ↦ sx,   ∂ ↦ s⁻¹∂,   X ↦ s⁻³X,   D ↦ s³D,
```

one has `κ₂ ↦ s⁵κ₂`. Thus normalization of a nonzero parameter chooses
`s⁵ = κ₂⁻¹`, which is available over `C̄` (or any algebraically closed
characteristic-zero field) but not generally over `ℚ`. Build the full constant-top system
`Q_m = δ_{m0}`, `m ∈ [−6,6]`, at free-polynomial-degree cap `d` (membership factors
`E`, `E(E−1)`, `E(E−1)(E−2)` on the level-`−1,−2,−3` coefficients). All coefficients
in this normalized fiber are **integer**; `msolve` requires cleared denominators (it
misparses rational monomials — a documented trap in the verifier).

| cap `d` | `κ₂ ≠ 0` (normalized `κ₂ = 1`) | engine(s) | tier |
|---|---|---|---|
| `d = 1` | **UNIT IDEAL at `κ₂=1` (hence EMPTY after base change)** | SymPy `ℚ` Gröbner `= [1]` + `msolve` `-g` unit over 3 primes *(default run)*; `msolve` char-0 `[-1]` *(additional)* | **committed** |
| `d = 2` | **historically recorded EMPTY at `κ₂=1` over `C̄`** | historical `msolve` char-0 `[-1]` transcript (≈35 min) | HEAVY leg implemented, but not completed in the current checkout audit |
| `d = 3` | attempted | `msolve` char-0 / `-g` over `F_p` — did **not** complete within the HEAVY time cap in development | HEAVY, not certified |

The `d = 1` row has a load-bearing SymPy-over-`ℚ` unit-ideal certificate for the
normalized `κ₂=1` ideal. A unit ideal over `ℚ` remains a unit ideal after base change, so
this proves emptiness of that normalized fiber over `C̄`; default finite-field `msolve -g`
is corroboration. A characteristic-zero `msolve [-1]` result is recorded as an
additional/manual computation rather than reproduced by the default verifier path. The
`d = 2` row records a historical exact characteristic-zero transcript for the same normalized
fiber via `msolve`'s rational solver (`[-1] =` empty variety). The HEAVY leg is implemented,
but it did not complete in the current checkout audit, so this is not a presently rerun result.
Scaling with `s⁵=κ₂⁻¹` then transports normalized-fiber emptiness to every nonzero `κ₂`
over the algebraic closure. It does not assert a rational-coordinate equivalence or an
independent computation for each nonzero rational fiber. `d = 3` is attempted but its
Gröbner did not finish in the budgeted window; no `d = 3` claim is made.

The `κ₂ = 0` slice is **nonempty**: the explicit positive control
`U = x + 2∂`, `X = U³ − ∂`, `D = U` is a genuine pair `[D,X] = 1` with `a₃ = 1`,
`b₂ = 0` (`verify §3`), and the `κ₂ = 0` system at `d = 1` is **not** the unit ideal.
Thus over an algebraically closed characteristic-zero field, within cap `d ≤ 2`, every
constant-top genuine pair lies in the `κ₂ = 0` slice.
The cited tame catalog proves the displayed single-shear-origin family lies on this
slice; it explicitly does **not** classify arbitrary tame words or all genuine pairs.
Therefore no reverse implication `κ₂=0 ⇒ tame` is claimed here, even within the cap.

**A methodological note recorded for the corpus.** A hand-rolled forward-solve of
the positive cascade produced a *rational-coefficient* reduced system; `msolve`
silently **misparsed** its rational monomials (`2*x^2/3`), reporting a spurious
non-unit basis while SymPy over `ℚ` (and the full integer system in both engines)
correctly reported the unit ideal. The verifier only ever hands `msolve`
**integer** systems. This is a second `msolve` trap beyond `**` vs `^`.

## 7. What the closure delivers for Gap 2, band 3

- **RESIDUAL 3 (`shifted-power-residuals.md` §4 pt 4).** The normalized `κ₂=1` fiber
  is disposed of **at bounded degree** (`d ≤ 2`) after base change to an algebraic closure:
  the exact `ℚ` certificate at that fiber remains valid after base change, and
  `κ₂ ↦ s⁵κ₂` with `s⁵=κ₂⁻¹` transports the result to every nonzero `κ₂` over an
  algebraically closed characteristic-zero field. Thus every pair in the certified cap over
  that field must lie on `κ₂ = 0`, a slice containing the displayed tame family but not
  classified here. The arbitrary-degree `(κ₂-closure)` remains **open**, with the
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
  decompositions of `Q₋₅, Q₋₄`; and the **bilinearity** of `[D,X]` (`§2`).

**Bounded / finite evidence (exact scope):**
- Generic sampled constant-top `X` is nonsolvable at the implemented cap
  `(dX,dD)=(1,3)`; rank/augmented-rank proves a cokernel obstruction at each tested
  seed but does not
  identify it with the moment slope or extend to arbitrary degree.
- The normalized `κ₂=1` constant-top fiber is empty after base change at cap `d = 1`
  (SymPy `ℚ` unit ideal `[1]`; `msolve -g` over three primes is corroboration;
  additionally `msolve` char-0 `[-1]`) and has a historical cap-`d = 2` transcript
  (`msolve` char-0 `[-1]`); the implemented HEAVY rerun did not complete in the current
  checkout audit. Because the variable/operator scaling sends `κ₂` to
  `s⁵κ₂`, choosing `s⁵=κ₂⁻¹` transports this to every nonzero `κ₂` over an
  algebraically closed characteristic-zero field, not generally by rational coordinates.
  `d = 3` attempted, not certified within the time budget.
- explicit `κ₂ = 0` tame witness `[D,X] = 1`; `κ₂ = 0` at `d = 1` not unit.

**Refuted (machine-checked) — corpus corrections:**
- That the top potential-factorization contributes to the constant-`h` closure: it
  is **vacuous** (`h^{[-1]} = 1`).
- That `msolve` may be fed rational-coefficient systems: it **misparses** them; only
  integer systems are certified.

**Open / NOT claimed:**
1. `(κ₂-closure)` at arbitrary degree for general `κ₂ ≠ 0` — the negative-tail
   first integral; the exact obstruction (§5) is the classical `e ≠ 0` /
   quantum A\*-I frontier. The Nonpositive-D Exclusion Theorem of
   [`shifted-cube-completion.md`](shifted-cube-completion.md) closes only the
   corner `κ₂ = 0 ∧ b₁ = 0`; the existing `κ₂ = 0` tame witness of §6 has
   positive band-one `D` (`b₁ ≠ 0`), so there is no contradiction with that theorem.
2. Whether `κ₂` is invariant under arbitrary composite tame words, and whether every
   `κ₂=0` constant-top pair is tame; the cited single-shear catalog does not classify these.
3. Non-2-separated tops, imbalanced coprime walls, general-`k` negative tail
   with `band D > 0`, **W2**, radical forcing at coupling widths `k = 4, 5`.
No Weyl pair, no counterexample; DC1/JC2 untouched.

## 9. Verification

```sh
uv run --with sympy python research/dc1-program/verify_astar_band3.py
HEAVY=1 uv run --with sympy python research/dc1-program/verify_astar_band3.py
```

Exact SymPy over `ℚ`: `§0` engine + telescoping potential; `§1` constant-top wall +
`Q₄` rung + Lemma-P slope + potential vacuity; `§2` bottom proportionality +
inhomogeneous `μ₃`-tail + bilinearity + capped sampled generic-`X` nonsolvability;
`§3` the **bounded certificate** (committed exact SymPy unit ideal at cap `d = 1`,
tame witness) with `msolve` corroboration gated by `shutil.which` (SKIPs cleanly), and
cap `d = 2, 3` behind `HEAVY=1`. Runtime varies with the environment. The final banner is
PASS/SKIP-aware: it says all checks passed only when none were skipped, and otherwise says
all executed checks passed while listing the skips.
