# Band-3 tame catalog: the single-shear-origin structure lemma, classical and quantum

**EXACT SYMPY — NORMAL-FORM STRUCTURE LEMMA (PROVED, ALL PARAMETERS) + BOUNDED ENUMERATION — NOT A FULL CLASSIFICATION**

Backing script: `verify_band3_catalog.py` (same directory), all arithmetic exact.
A successful run ends `ALL BAND3 CATALOG CHECKS PASSED`.

Conventions frozen from band 2 (classical assembly `f637b1a`, quantum assembly
`84978b9`, catalog `verify_catalog.py`): classical `{G,F}=G_ξ F_x−G_x F_ξ`,
`τ=xξ`, `{ξ,x}=1`, `F=Σ_{k=-3}^{3}x^k a_k(τ)`, `G=Σ x^k b_k(τ)`, membership
`τ^j | a_{-j}, b_{-j}`, coefficient equations `C_m=Σ_{k+l=m}(k a_k b_l' − l a_k' b_l)=δ_{m0}`
for `m∈[−6,6]`. Quantum `A_1[x^{-1}]=⊕_k x^k C[E]`, `E=x∂`,
`(x^a f(E))(x^b g(E))=x^{a+b}f(E+b)g(E)`, membership
`x^{-j}c(E)∈A_1 ⟺ E(E−1)⋯(E−j+1) | c(E)`,
`Q_m=Σ_{i+l=m}[b_l(E+i)a_i(E) − a_i(E+l)b_l(E)]=δ_{m0}`.

---

## 0. Strategic headline

Band 2's load-bearing discovery was structural: every band-2 tame pair's two
*extremes* (`a_{±2}`) originate in a **single quadratic shear** `μX²`, so one
gauge `G→G−λF` kills `b_2` and `b_{−2}` at once and the resistant branch A\* is
never reached (`tame-catalog.md`, `f637b1a`). Band 3 sharpens this into a
**tower** statement and adds a genuinely new phenomenon.

Three exact structural facts, all machine-verified for **all parameters**:

1. **Blow-up law for the displayed two-shear families.** Composing the
   normalized monic, nonzero nonlinear `x`- and `ξ`-shears used in `verify` §4
   produces band **exactly `a·b`**: the displayed extreme coefficient is nonzero.
   Thus two quadratics in those opposite-direction families have band 4, while
   the bounded quadratic/affine enumeration reaches no band-3 word. Zero or
   affine shear terms, more general interleavings, and possible cancellations
   require separate treatment; no unbounded tame-word theorem is inferred.

2. **Single-shear-origin, upgraded to a tower.** In the normal form
   `Φ = A₂ ∘ S ∘ A₁` with `A₁,A₂` affine symplectic and `S` a **single**
   triangular shear (one coordinate stays linear), `F` and `G` are both
   `C`-linear combinations of two functions `X` (linear, band 1) and
   `Ψ` (sheared, band 3). The **entire weight-`|k|≥2` tower** of `F` and `G`
   lives in `Ψ`, so the all-parameter statement is the division-free identity
   `b·b_k=d·a_k` for every `|k|≥2`. If `b≠0`, this becomes the common ratio
   `b_k=(d/b)a_k`, and `G→G−(d/b)F` annihilates **`b_3, b_2, b_{−2}, b_{−3}`
   simultaneously**. If `b=0`, symplecticity `ad−bc=1` forces `d≠0`, so the
   opposite output `G` carries the tower; exchanging `(F,G)` for `(G,−F)` gives
   an orientation with nonzero first tower coefficient before applying the
   ratio/gauge statement. `verify` §3, §7 (classical and quantum, all parameters).

3. **Cross-coupling of proportionality constants → the resistant locus.** On
   the two-sided branch `a_3a_{−3}≠0`, the top Wronskian `C_6` gives
   `b_3=λ a_3` and the bottom Wronskian `C_{−6}` gives
   `b_{−3}=μ a_{−3}`. The **A\*-band3 resistant branch** is exactly
   `μ≠λ` (gauged `b_{−3}=(μ−λ)a_{−3}≠0`). The single-shear normal form forces
   `μ=λ` because both extremes come from one `Ψ`, and no word in the stated
   bounded enumeration reaches `μ≠λ`. The heuristic that distinct positive and
   negative cubic sources should blow the band up to at least 9 is proved only
   for the displayed nondegenerate two-shear families; it is not an unbounded
   classification. Thus resistant-locus emptiness remains conjectural and must
   be settled from the coefficient equations.

Everything here is band-scoped and claims neither JC2 nor DC1.

---

## 1. The thirteen coefficient equations and the top/bottom cascade

`verify` §1 checks that the thirteen `C_m` (`m=−6..6`) match the direct
two-variable Poisson expansion of `{G,F}` for the fully generic band-3 pair.
The two extreme rungs and the rung one down are the load-bearing ones.

**Top Wronskian `C_6` (weapon W1).** On the branch `a_3≠0`,
`C_6 = 3(a_3 b_3' − a_3' b_3) = 0 ⟹ (b_3/a_3)' = 0 ⟹ b_3 = λ a_3`, `λ∈C`.
Gauge `G→G−λF` sets `b_3=0`. If `a_3=0`, `C_6=0` is vacuous and belongs
to a separate vanishing branch. `verify` §2.

**Cube-class `C_5` (band-3 M4 analogue).** On the same `a_3≠0` branch, after
the gauge, `C_5` reduces to `3 a_3 b_2' − 2 a_3' b_2 = 0`. Its integrating
invariant is `b_2^3/a_3^2`, with the denominator-safe cleared identity
`a_3^3·d/dτ(b_2^3/a_3^2) = b_2^2(3 a_3 b_2' − 2 a_3' b_2)` exactly
(`verify` §2). Hence, for `b_2≠0`,

```
C_5 = 0  ⟺  b_2^3 = c·a_3^2  ⟺  a_3 is a scalar·cube  (a_3 = c h^3),  b_2 = κ h^2.
```

This is the **cube class**: where band 2 had `a_2=h²` (square) with `b_1=κh`,
band 3 has `a_3=h³` (cube) with `b_2=κh²`, one weight-level lower. Proving `a_3`
is forced to be a cube even when `b_2=0` (the full band-3 M4 theorem, needing
the whole cascade) is D1/D2's job; here it is recorded as the cascade's shape.

**Bottom Wronskian `C_{−6}`.** If `a_{−3}≠0`,
`C_{−6} = 3(a_{−3}' b_{−3} − a_{−3} b_{−3}') = 0 ⟹ b_{−3} = μ a_{−3}`, `μ∈C`.
If `a_{−3}=0`, this extreme equation is vacuous and belongs to a separate
vanishing branch. Thus `(λ,μ)` is the **cross-coupling** datum on the two-sided
branch `a_3a_{−3}≠0`. `verify` §2.

---

## 2. The structure lemma (proved for the normal-form class)

> **Structure Lemma (band-3 single-shear-origin).** Let
> `Φ = A₂ ∘ S ∘ A₁` where `A₁,A₂` are affine symplectic (SL₂) and `S` is a
> *single* triangular shear `Ξ ↦ Ξ + p(X)` in the `A₁`-coordinates,
> `deg p ≤ 3`. Write `X = A₁`-image of the unsheared coordinate (linear,
> band 1) and `Ψ = Ξ + p(X)` (band `deg p`). With `A₂=[[a,b],[c,d]]`,
> `F = aX + bΨ`, `G = cX + dΨ`. Then:
>
> (i) **Tower proportionality.** For every `|k|≥2`, `(X)_k=0` (X is band 1), so
> `a_k = b·Ψ_k` and `b_k = d·Ψ_k`; hence the division-free identity
> `b·b_k=d·a_k` holds for all parameters and every `|k|≥2`.
>
> (ii) **One oriented gauge kills the tower.** If `b≠0`, the common ratio is
> `λ=d/b`, and `G→G−λF` annihilates `b_3, b_2, b_{−2}, b_{−3}` simultaneously.
> If `b=0`, symplecticity forces `d≠0`; the tower is then carried by `G`, and
> pair exchange `(F,G)→(G,−F)` supplies an orientation with nonzero first tower
> coefficient to which the ratio and gauge apply.
>
> (iii) **λ = μ after orientation.** The oriented top and bottom proportionality
> constants coincide, so gauged `b_{−3}=0`: the pair is **never** in the
> resistant branch.
>
> (iv) **Support.** If `p` is a pure cubic `μX³`, the support is odd-only
> `{−3,−1,1,3}` (`a_2=a_0=a_{−2}=0`); a lower-degree tail in `p` adds the even
> levels. Genuine polynomial membership `τ^j | a_{−j}` holds automatically.

*Proof of (i)–(iii).* `X` is a nonzero linear form in `x,ξ`, hence weight `±1`
only: its `x^k`-coefficient vanishes for `|k|≥2`. `Ψ = Ξ + p(X)` carries all the
weight-`≥2` content. Since `F,G` are the *same* two functions `X,Ψ` combined
with constant rows `(a,b)` and `(c,d)`, every weight-`|k|≥2` coefficient obeys
`b·b_k=d·a_k`, without division. When `b≠0` this gives the ratio `d/b`,
independent of `k`. When `b=0`, `ad−bc=1` gives `d≠0`; thus `F` has no tower,
`G` carries it, and pair exchange gives the nonzero-first-output orientation used
for the ratio and tower-killing gauge. `∎`

`verify` §3 proves (i)–(iv) symbolically for the fully generic
`α,β,γ,δ,a,b,c,d,μ₃,μ₂,μ₁` (the most general single shear of degree ≤ 3),
including the membership divisibilities. **This is a proof for the normal-form
class at all parameters, not a bounded check.** What is *not* claimed here is
that every unbounded tame band-3 word reduces to this normal form — exactly the
scope caveat of the band-2 catalog.

**Why this is richer than band 2.** Band 2's lemma coupled the two *extremes*
`a_{±2}`. Band 3 couples an entire two-level tower `{a_{±2}, a_{±3}}` of `G` to
that of `F` through the single `Ψ`, and the cascade `C_5` ties the middle level
`b_2` to the top `a_3` via the cube class. The "cross-coupling of proportionality
constants" flagged as band 3's first new phenomenon is exactly the forced
coincidence `λ_top = μ_bottom` (and the `C_5`-mediated middle relation) inside
the single-shear form.

---

## 3. The interleaving question: two quadratic shears never reach band 3

The mandate asks which compositions of two quadratic shears stay inside band 3.
The blow-up law answers it completely.

- **Opposite directions in the displayed nondegenerate families** (`x`-shear
  then `ξ`-shear, or vice versa), with nonzero leading coefficients and degrees
  `a,b≥2`: band `=a·b`, because the computed extreme coefficient is nonzero.
  Two quadratics in this family give band 4. If a leading coefficient vanishes,
  the effective degree drops; longer words may permit cancellations not covered
  by this two-shear calculation. `verify` §4.
- **Same direction**: two `ξ`-shears `Ξ↦Ξ+p₁(X)` then `Ξ↦Ξ+p₂(X)` compose
  additively (`Ξ↦Ξ+p₁+p₂`), so the band is `max(deg p₁, deg p₂)`, never new
  weight `±3` from two quadratics.
- **Enumeration**: over affine SL₂ + quadratic shears, words of length ≤ 4
  realize bands `{1,2,4,8,16}` and **never 3** (`verify` §4); the
  quad-affine-quad sweep produces **zero** genuine band-3 pairs (`verify` §6).

**Bounded conclusion.** No quadratic–quadratic interleaving tested in §4 or §6
lands in band 3. The exact single-shear lemma explains every displayed catalog
entry and forces `λ=μ` there. It does **not** prove that every unbounded tame
band-3 word contains, or reduces to, one cubic shear; zero/extreme degeneracies,
longer interleavings, and cancellations remain outside this catalog's theorem.

---

## 4. The catalog

Gauge-reduced data per entry (`λ` = top constant `b_3/a_3`; `b_2` before gauge;
`gauged b_{−3}` = `(μ−λ)a_{−3}`; branch = residence in the case tree). Every
entry is machine-verified `{G,F}=1`, band 3, oriented to `a_3≠0`, and classified
by `verify` §5. **`B0-band3`** = collapse (`μ=λ`, gauged `b_{−3}=0`);
**`A*-band3`** = resistant (`μ≠λ`).

| # | tame construction | supp `F` / supp `G` | `a_3` | `a_2` | `a_{−3}` | `λ` | `b_2` | gauged `b_{−3}` | branch |
|---|---|---|---|---|---|---|---|---|---|
| B3-1 | elementary cubic shear `(x, ξ+x³)` | `{1}` / `{−1,3}` | `1`* | `0` | `0` | `0` | `0` | — | one-sided top |
| B3-2 | `A₂[[1,1],[1,2]] ∘ (Ξ↦Ξ+X³+X²) ∘ A₁[[1,1],[0,1]]` | `{−3..3}` / `{−3..3}` | `1` | `1` | `τ³` | `2` | `2` | `0` | **B0-band3** |
| B3-3 | `A₂[[2,1],[1,1]] ∘ (Ξ↦Ξ+X³) ∘ A₁[[1,1],[0,1]]` | `{−3,−1,1,3}` / same | `1` | `0` | `τ³` | `1` | `0` | `0` | **B0-band3** |

*B3-1 is oriented by pair exchange so that the cube sits in `a_3`; before
orientation it is `(x, ξ+x³)` with the cube in `b_3`.

**Full coefficient profiles** (exact, `verify` §5 / `profiles`):

- **B3-1** `(x, ξ+x³)`: `a_1=1`; `b_3=1, b_{−1}=τ`. The minimal band-3 pair;
  the cubic shear `ξ↦ξ+x³` creates only the `+3` extreme, so it is **one-sided**
  (`a_{−3}=0`). No resistant question arises.

- **B3-2** (two-sided, **even** levels): the E4-band3 analogue.
  ```
  a_3=1,  a_2=1,  a_1=3τ+1,  a_0=2τ,  a_{−1}=3τ²+2τ,  a_{−2}=τ²,  a_{−3}=τ³
  b_3=2,  b_2=2,  b_1=6τ+1,  b_0=4τ,  b_{−1}=6τ²+3τ,  b_{−2}=2τ²,  b_{−3}=2τ³
  ```
  Both negative extremes present *before* gauging (`a_{−3}=τ³`, `b_{−3}=2τ³`),
  the raw profile that A\*-band3 would need. But `b_{−3}=λ·a_{−3}=2τ³`, and the
  gauge `G→G−2F` sends the **entire tower** `b_3,b_2,b_{−2},b_{−3}→0`. Cube data:
  `a_3=1=1³` (`h=1`), `b_2=2=2·1²` (`κ=2`). Branch **B0-band3**.

- **B3-3** (two-sided, **odd**-only `{−3,−1,1,3}`):
  ```
  a_3=1,  a_1=3τ+2,  a_{−1}=3τ²+3τ,  a_{−3}=τ³
  b_3=1,  b_1=3τ+1,  b_{−1}=3τ²+2τ,  b_{−3}=τ³   (a_2=a_0=a_{−2}=b_2=b_0=b_{−2}=0)
  ```
  Pure cubic shear → no even levels. `λ=1`, gauge kills `b_3` and `b_{−3}`
  together. Branch **B0-band3**.

**No entry in A\*-band3** — none exists in the stated enumeration (§5).

### Gauge-reduced data for D1/D2's case tree

D1/D2's memos may not exist yet. The residence datum they will read off is the
triple, after gauge `G→G−λF`:

```
(a_3,  a_2 = gauged b_2 free-vs-cube,  gauged b_{−2},  gauged b_{−3}).
```

For every catalog entry and every one of the 1216 enumerated genuine band-3
pairs (§5): gauged `b_{−3}=0`. The branch split is:
`onesided-top` (`a_{−3}=0`) · `B0-band3` (`a_{−3}≠0`, gauged `b_{−3}=0`) ·
`A*-band3` (`a_{−3}≠0`, gauged `b_{−3}≠0`, **conjecturally empty**). This mirrors
band 2's `B0 / A0 / A*` split one weight-level up, with the extra middle-tier
coupling `b_2 = κ h²` from `C_5`.

---

## 5. Bounded enumeration

Words of length ≤ 3 over 10 affine SL₂ generators + 6 cubic shears + 6 quadratic
shears produced **4620** distinct tame pairs. Band-max histogram
`{1,2,3,4,6,8,9,12,18,27}` — every value a product of shear degrees, confirming
the blow-up law (`27=3³`, `9=3²`, `8=2³`, …). Among them, **1216** are genuine
band-3 (band exactly 3). After orientation to `a_3≠0` (pair exchange /
reflection, band-2 §3 orientation lemma, transferring verbatim), all 1216
classify:

```
one-sided top (a_{-3}=0):            984
two-sided collapse  (μ=λ, B0-band3): 232
RESISTANT           (μ≠λ, A*-band3):   0
```

- **Rigorous for all parameters:** the structure-lemma identities (§2), the
  cascade identities (§1), the blow-up law `band=a·b` for the displayed shear
  compositions (§3).
- **Enumerated (corroboration only):** `0` of `1216` genuine band-3 tame words
  reach the resistant branch; `0` two-quadratic words reach band 3.
- **Not claimed:** exhaustion over unbounded word length; exclusion of *all*
  tame automorphisms from A\*-band3; the full band-3 M4 (`a_3` a cube when
  `b_2=0`); any JC2/DC1 statement.

The bounded enumeration suggests A\*-band3 is empty; this memo does **not**
promote that to a theorem. The written proof (from the coefficient equations,
the band-2 `classical-Astar.md` route lifted one level) is D1/D2's deliverable —
this catalog supplies the positive controls (B0-band3 nonempty) and the
conjectural counterexample-or-nothing locus (A\*-band3).

---

## 6. Quantum band 3

Everything deforms one floor up the quantization ladder (`verify` §7).

- **`Q_m` match** the abstract crossed-product commutator for `m=−6..6`.
- **Top Wronskian/Casoratian `Q_6`**: `b_3(E+3)a_3 − a_3(E+3)b_3 = 0 ⟹ b_3=λ a_3`.
- **Shifted-cube `Q_5`**: after gauge, `b_2(E+3)a_3 − a_3(E+2)b_2 = 0`. The
  **shifted-cube class** solves it:
  ```
  a_3 = c·h(E)h(E+1)h(E+2)   ⟹   b_2 = κ·h(E)h(E+1),
  ```
  verified for several concrete `h(E)` (`E`, `E+5`, `E²−1`, const, `2E−7`).
  This is the exact deformation of the classical cube `a_3=h³ ⟹ b_2=κh²`: the
  three literal factors of `h³` become the three shifted factors
  `h(E)h(E+1)h(E+2)`, and the two of `h²` become `h(E)h(E+1)`.
- **Single-shear-origin (quantum), all parameters**: for the ladder normal form
  `A₂ ∘ (∂↦∂+μX³) ∘ A₁`, the division-free identity `b·b_k=d·a_k` holds for
  `|k|=3`. In an orientation with `b≠0`, `λ_top=μ_bottom=d/b` and the gauge kills
  `b_{±3}`. If `b=0`, symplecticity forces `d≠0`; pair exchange makes the
  tower-carrying output first and supplies such an orientation. Membership is the
  **falling factorial**: `a_{−3}` is divisible by `E(E−1)(E−2)` automatically
  (`verify` §7) — the deformation of `τ³ | a_{−3}`, the mechanism that closes
  the quantum extreme where polynomiality closed the classical one.
- **Orientation** uses the Fourier automorphism `φ: x↦−∂, ∂↦x`, genuinely
  automorphic (`[φ∂,φx]=1`, no sign flip), verified in `verify` §7.

**Quantum catalog entries** (`[D,X]=1` exact):

| # | construction | supp `X` / supp `D` | `a_3` | `a_{−3}` (falling fact.) | `λ` | branch |
|---|---|---|---|---|---|---|
| Q3-1 | `X=x`, `D=∂+μx³` | `{1}` / `{−1,3}` | — | — | — | one-sided (`b_3=μ`, `b_{−1}=E`) |
| Q3-2 | `A₂[[1,1],[1,2]] ∘ (∂↦∂+X³) ∘ A₁[[1,1],[0,1]]` | `{−3,−1,1,3}` / same | `1` | `E(E−1)(E−2)` | `2` | **B0-band3** (gauge kills `b_{±3}`) |

Q3-2 profile: `a_1=3E+4`, `a_{−1}=E(3E+2)`, `a_{−3}=E(E−1)(E−2)`;
`b_3=2`, `b_1=6E+7`, `b_{−1}=3E(2E+1)`, `b_{−3}=2E(E−1)(E−2)`; gauge `D→D−2X`
kills `b_3` and `b_{−3}` simultaneously. The bottom coefficient `a_{−3}` is
*exactly* the falling factorial `E(E−1)(E−2)` — genuine `A_1`-membership with no
slack. Same structural verdict as classical: single cubic shear ⟹ `λ=μ` ⟹ no
resistant branch.

---

## 7. What is proved vs computed vs conjectured

**Proved (exact, all parameters / arbitrary degree in the stated class):**
- the thirteen `C_m` and their quantum `Q_m` (§1, §6);
- top/bottom Wronskian proportionalities and the cube / shifted-cube cascade
  identities (§1, §6);
- the single-shear-origin structure lemma for `A₂∘S∘A₁`, classical and quantum:
  `b·b_k=d·a_k` at tower levels for all parameters, plus `λ=μ` and the
  tower-killing gauge after `b≠0` orientation (including the `b=0` pair exchange)
  (§2, §6);
- the blow-up law `band=a·b` for the displayed nondegenerate
  opposite-direction two-shear families (§3).

**Computed (bounded, corroboration only):**
- 1216 genuine band-3 tame pairs, `0` resistant, `0` two-quadratic pairs in
  band 3 (§5).

**Conjectured (counterexample-or-nothing, to be settled by D1/D2 from the
coefficient equations, lifting `classical-Astar.md` / `quantum-mirror.md`):**
- **A\*-band3 is empty**: no genuine band-3 Keller/Weyl pair has gauged
  `b_{−3}≠0` (i.e. `μ≠λ`). Positive control: B0-band3 is nonempty (B3-2, B3-3,
  Q3-2). A motivating heuristic is that independent positive and negative cubic
  sources produce band 9 in the displayed two-shear models; no general
  source-decomposition or unbounded-word obstruction is proved here.
- the **full band-3 M4**: `a_3` is a cube (resp. shifted cube) even when `b_2=0`.

No JC2, DC1, or unbounded-word claim is made anywhere in this memo.

---

## Reproduce

```
uv run --with sympy python research/band3/verify_band3_catalog.py
```

Sections: (1) thirteen `C_m`; (2) top/bottom cascade + cube invariant; (3)
classical structure lemma at generic parameters; (4) blow-up law + powers-of-2
enumeration; (5) three classical catalog entries; (6) 4620-pair bounded
enumeration (1216 genuine, 0 resistant) + two-quadratic sweep; (7) quantum
`Q_m`, shifted-cube cascade, quantum structure lemma, falling-factorial
membership, Fourier automorphism. The script ends `ALL BAND3 CATALOG CHECKS PASSED`.
