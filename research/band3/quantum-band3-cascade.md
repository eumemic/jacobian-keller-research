# The quantum band-3 cascade: top proportionality, the Q₅ wall, telescoping, and the negative tail

**INDEPENDENTLY DERIVED AND MACHINE-VERIFIED — NOT PEER REVIEWED — BAND-SCOPED**

QUANTUM MIRROR assault, band 3. This memo opens the next floor of the width
induction on the Dixmier (DC1) face, one rung above the assembled quantum
band-2 theorem (`../band2-square-sector/quantum-band2-theorem.md`, commit
`84978b9`). It derives the complete `Q_m` system for `E`-ladder supports in
`[-3,3]` (`m ∈ [-6,6]`), machine-verified against the crossed-product
commutator with generic coefficients, and pushes the descent as far as the
top proportionality, the `Q₅` **wall**, the `m=0` telescoping identity, and the
first negative-tail rigidities. It is the independent quantum sibling of the
classical band-3 work (developed in parallel; the point is an independent
derivation, not a shared one).

Everything below is checked exactly by
[`verify_quantum_band3.py`](verify_quantum_band3.py) (ends
`ALL QUANTUM BAND3 CHECKS PASSED`).

Every ratio statement below is restricted to the branch where its denominator
extreme is nonzero. A vanishing extreme is treated as a separate branch rather
than divided through.

## 0. Headline

> The naive mirror **breaks at band 3.** The band-2 wall `Q₃` forced the top to
> be a *shifted square* `a₂ = c·h h⁽¹⁾` — cleanly, with no arithmetic residue.
> Its band-3 analogue `Q₅` does **not** force a *shifted cube*
> `a₃ = c·h h⁽¹⁾ h⁽²⁾`. The exact condition is strictly weaker: in each
> `mod-ℤ` coset the root multiset of `a₃` must be divisible by the cyclotomic
> `Φ₃(S)=1+S+S²`, together with a realizability (nonnegativity) side condition —
> and the shifted cube is only the special sub-case where the *cofactor* is
> itself effective. The polynomial `a₃ = E(E−2)(E−4)` (roots `{0,2,4}`, an
> arithmetic progression of step 2, **not** a shifted cube) solves the wall with
> `b₂ = (E−1)(E−4)`. This is the first genuinely new band-3 phenomenon, and it
> refutes the conjecture that the wall alone forces a shifted cube.

What survives the mirror intact: the top proportionality (`Q₆`, now a **3**-fold
periodicity), the `m=0` telescoping to a local potential (in closed form for
*every* band), and the falling-factorial membership machinery. What is new:
the `Φ₃`/sesqui structure of the wall, the forced `2/3` degree ratio
`2 deg a₃ = 3 deg b₂`, and the **cross-coupling of the two proportionality
constants** `λ₃` (top) and `μ₃` (bottom) — the effect the archived milestone
flagged as the band-3 obstruction, here made exact.

## 1. Setup and conventions (frozen)

Work in `A₁[x⁻¹] = ⊕_{k∈ℤ} x^k ℂ[E]`, `E = x∂`, `∂ = x⁻¹E`, with
```
(x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E),     f^[r](E) := f(E+r),
```
and the band-3 presentation
```
X = Σ_{k=-3}^{3} x^k a_k(E),      D = Σ_{l=-3}^{3} x^l b_l(E).
```
The ladder-`m` coefficient of `[D,X] = DX − XD` is
```
Q_m = Σ_{k+l=m} ( b_l^[k] a_k − a_k^[l] b_l ),        [D,X] = 1  ⇔  Q_m = δ_{m0},
```
for `m ∈ [-6,6]`. **Genuine `A₁` membership** (quantum polynomiality):
`x⁻ʳ c(E) ∈ A₁` iff the falling factorial `E^{underline r} = E(E−1)⋯(E−r+1)`
divides `c(E)`. In band 3:
```
E | a_-1, b_-1 ;    E(E−1) | a_-2, b_-2 ;    E(E−1)(E−2) | a_-3, b_-3.
```
`verify` §0 checks `Q_m` equals the direct crossed-product commutator for every
`m ∈ [-6,6]` with fully generic degree-2 coefficients, and records the closed
forms of `Q₆, Q₅, Q₀, Q₋₅, Q₋₆`.

`T` denotes `f ↦ f^[1]`. `S` denotes the shift acting on **root positions**
(`S·t^γ = t^{γ+1}`), used only in the root-multiset calculus of §3. `Φ₃(S) =
1+S+S²`.

## 2. `Q₆` — top proportionality (3-periodicity) and the gauge

Assume `a₃≠0`. Only `(k,l)=(3,3)` contributes to `Q₆`:
```
Q₆ = b₃^[3] a₃ − a₃^[3] b₃ .
```
`Q₆ = 0` reads `b₃^[3]/a₃^[3] = b₃/a₃`, so the rational function `r = b₃/a₃`
satisfies `r^[3] = r`: it is **3-periodic**. A rational function with a nonzero
additive period is constant (a pole would generate an infinite orbit spaced by
3, but a rational function has finitely many poles — the rational periodicity
lemma used verbatim at band 2, `../band2-square-sector/quantum-M4.md`
`Q₄` step). Hence
```
b₃ = λ₃ a₃ ,     λ₃ ∈ ℂ.
```
`verify` §1 certifies the periodicity engine exactly as an identity of rational
functions, `r^[3] − r = (b₃^[3]a₃ − a₃^[3]b₃)/(a₃ a₃^[3])`, and records the
period-3 leading-coefficient rigidity for the *polynomial* case (`3-step` finite
difference of a degree-`d` polynomial has degree `d−1` and leading coefficient
`3d·lc`, nonzero for `d ≥ 1`).

**Gauge.** Subtracting `λ₃X` from `D` realizes `b₃ = 0`. The gauge `D ↦ D − λ₃X`
preserves `[D,X]=1` (since `[X,X]=0`), band-3 support, the generated subalgebra,
and — because `b_{-r}, a_{-r}` both carry `E^{underline r}` — genuine membership
of every gauged coefficient (`verify` §1, all `m` and all `r`). We henceforth
work in the gauge `b₃ = 0`.

## 3. `Q₅` — THE QUANTUM WALL (the band-3 J2q)

In the gauge `b₃ = 0`, the only surviving pair at `m=5` is `(k,l)=(3,2)`:
```
Q₅ = b₂^[3] a₃ − a₃^[2] b₂        (verify §0).
```
So `Q₅ = 0` is the **staggered homogeneous equation**
```
        b₂(E+3) · a₃(E) = a₃(E+2) · b₂(E) .                    (WALL)
```
This is the band-3 analogue of the band-2 J2q equation
`b₁(E+2)a₂(E) = a₂(E+1)b₁(E)` (which forced `a₂ = c·h h⁽¹⁾`). **The staggering
is asymmetric — a `+3` shift on `b₂` against a `+2` shift on `a₃` — and that
single fact is what makes band 3 genuinely new.** In band 2 the shifts were `+2`
against `+1`, and the clean invariant `r = a₂/(b₁ b₁^[1])` with `r^[1]=r` closed
it. Here `deg a₃` and `deg b₂` are not in a `1:1` ratio (see §3.3), so no such
single-product invariant exists.

### 3.1 The exact reduction to a cyclotomic divisibility

Assume `a₃ ≠ 0` and take a **nonzero** `b₂` solving (WALL). Both sides are
nonzero polynomials; equating leading coefficients is automatic
(`lc(b₂)lc(a₃)` both sides), and equating **root multisets** gives, with
`A, B` the root-multiset Laurent data of `a₃, b₂` (positions grouped by
`mod-ℤ` coset; `S` = shift by `+1`):
```
   (roots b₂ − 3) ⊎ (roots a₃) = (roots a₃ − 2) ⊎ (roots b₂)
   ⇔  S⁻³B + A = S⁻²A + B
   ⇔  S(1+S)·A = (1+S+S²)·B .                                  (WALL-M)
```
(Derivation: `(1−S⁻²)A = (1−S⁻³)B`; multiply by `S³` and cancel the common
non-zero-divisor `(S−1)`.) Cosets do not mix under `S`, so (WALL-M) **decouples
coset by coset**. Within one coset, `Φ₃(S)=1+S+S²` is a non-zero-divisor in
`ℂ[S,S⁻¹]`, so
```
   B = S(1+S)·A / Φ₃(S)                                        (B-FORMULA)
```
is **uniquely determined** by `A`. Consequently:

> **Wall Lemma (band-3 J2q).** Fix `a₃ ≠ 0` and the gauge `b₃=0`. A nonzero
> `b₂` solving (WALL) exists **iff**, in every `mod-ℤ` coset `c`,
>  (i) the root multiset `A_c` of `a₃` is divisible by `Φ₃(S)=1+S+S²`, **and**
>  (ii) the forced quotient `B_c = S(1+S)A_c/Φ₃(S)` has nonnegative
>  coefficients (is a genuine multiset).
> When it exists, `b₂` is unique up to a scalar `κ₂` (freedom **exactly
> 1-dimensional**, `κ₂ = 0` included). Roots of `b₂` lie only in cosets where
> `a₃` has roots (elsewhere `A_c=0` forces `B_c=0`).

Because `gcd(Φ₃, S(1+S)) = 1`, condition (i) is equivalently `Φ₃ | A_c`
directly. `verify` §2 machine-checks the 1-dimensionality by direct linear solve
(cube top of degree 6 and the counterexample top of degree 3 each give a
single-scalar solution space), and reconstructs `b₂` from `A` via (B-FORMULA)
for representative `a₃`, confirming it solves (WALL).

### 3.2 Shifted cube is **sufficient but not necessary** — the refutation

**Sufficiency (mirror survives one way).** If `a₃ = c·h h⁽¹⁾ h⁽²⁾` (shifted
cube), then `b₂ = κ·h h⁽¹⁾` solves (WALL): in root data
`A = (1+S⁻¹+S⁻²)η = S⁻²Φ₃(S)η` (with `η =` roots of `h`), and (B-FORMULA) gives
`B = S(1+S)S⁻²η = (1+S⁻¹)η`, i.e. `b₂ = κ·h h⁽¹⁾`. `verify` §2 checks this at
symbolic `deg h = 1,2,3`.

**Non-necessity (mirror breaks the other way).** Take
```
   a₃ = E(E−2)(E−4)      (roots {0,2,4}),     b₂ = (E−1)(E−4).
```
Then `b₂^[3] a₃ = (E+2)(E−1)·E(E−2)(E−4) = a₃^[2] b₂` **exactly** (`verify` §2),
so `a₃` solves the wall. But `a₃` is **not** a shifted cube: a degree-1 `h`
(the only possibility, since `deg a₃ = 3`) with `c·h h⁽¹⁾ h⁽²⁾ = a₃` has no
solution (`verify` §2 solves the coefficient system and returns ∅). In root
data, `A = 1+S²+S⁴ = Φ₃(S)·(1−S+S²)`: divisible by `Φ₃`, but the cofactor
`C = 1−S+S²` has a **negative** coefficient, so `A` is not `Φ₃·(effective)`,
i.e. not a shifted cube. Yet `B = S(1+S)C = S+S⁴ ≥ 0`, a genuine multiset —
hence a valid `b₂`. **The shifted-cube conjecture for the band-3 wall is false.**

The precise gap: shifted cube ⇔ the *cofactor* `C := A/Φ₃` is effective
(`C ≥ 0`); the wall only needs the *image* `B = S(1+S)C` to be effective
(`B ≥ 0`). The map `C ↦ S(1+S)C` can send a non-effective `C` to an effective
`B` (as `1−S+S² ↦ S+S⁴`). This is the exact source of the new solutions.

### 3.3 The forced `2/3` degree ratio (a clean Lemma-R at the top)

Summing (WALL-M) over positions gives `2|A| = 3|B|`, i.e.
```
   2 deg a₃ = 3 deg b₂ .                                       (WALL-DEG)
```
This is also visible as a **staggered leading-coefficient** identity, the
band-3 top instance of the Lemma-R family (§6): writing
`Q₅ = a₃^[0] b₂^[3] − a₃^[2] b₂^[0]` with `p = deg a₃`, `q = deg b₂`, the
`E^{p+q}` terms cancel identically and
```
   coeff(E^{p+q−1}) of Q₅  =  (3q − 2p)·lc(a₃)·lc(b₂)          (verify §2),
```
so `Q₅ = 0` forces `3q = 2p`. In particular `3 | deg a₃`, and the degree ratio
`deg b₂ : deg a₃ = 2 : 3` is the quantum incarnation of the classical square
sector's `2/3`-power bookkeeping. (No cancellation caveat here: this is one
coefficient of the equation and `lc(a₃)lc(b₂) ≠ 0`.)

### 3.4 Realizability is a genuine second condition

`Φ₃`-divisibility of `A` alone is **not** sufficient. The top
`a₃ = (E−1)(E−5)(E−6)` (roots `{1,5,6}`) has `A = S+S⁵+S⁶ = Φ₃(S)·(S−S²+S⁴)`,
divisible by `Φ₃`, yet (B-FORMULA) gives `B = S²−S⁴+S⁵+S⁶` with a **negative**
coefficient — no effective multiset, hence **no nonzero `b₂`**. `verify` §2
confirms by direct solve (only `b₂ = 0` for `deg b₂ ≤ 8`). So the Wall Lemma's
two conditions are both necessary, and the bounded corroboration (`verify` §7)
finds that among single-coset 3-root tops `{0,a,b}` exactly `{0,1,2}` (the
consecutive cube `h=E`) and `{0,2,4}` (the counterexample) admit a `b₂`.

### 3.5 What the wall does and does not settle

The Wall Lemma classifies the **homogeneous freedom of one equation** `Q₅`. It
does **not** assert that a wall-admitting non-cube top extends to a genuine
band-3 pair `[D,X]=1`. Whether it does is the band-3 analogue of the band-2
question "is `h` forced constant?" — which band 2 answered *yes* (`h` constant,
`../band2-square-sector/quantum-completion.md` §4, and nonconstant-`h`
killed by membership). At band 3 the corresponding statement — *"is `a₃` forced
to a trivial cube (constant), so that the counterexample tops are killed
downstream?"* — is **open** and is the sharpest single question this cascade
exposes. The wall counterexample shows it cannot be settled by `Q₅` alone.

## 4. `Q₀` — the `m=0` telescoping (W4/W6 generalization, closed form, any band)

For any band `K`, the `m=0` commutator coefficient is a total finite difference
with an explicit **local** potential:

> **Telescoping Lemma.** `Q₀ = (T−1)G` with
> ```
>   G = Σ_{k=1}^{K} Σ_{j=0}^{k-1} ( a_k^[j−k] b_{-k}^[j]  −  b_k^[j−k] a_{-k}^[j] ) .
> ```

*Proof (exact).* Pair the `m=0` summands `(k,−k)` for `k>0` with `(−k,k)`. The
identity `g_k^[1] − g_k = a_k b_{-k}^[k] − a_k^[-k] b_{-k}` holds for
`g_k = Σ_{j=0}^{k-1} a_k^[j−k] b_{-k}^[j]` (telescoping in `j`), and the
mirror identity closes the `(−k,k)` block with a sign. Summing over `k=1..K`
gives the displayed `G`. `verify` §3 checks `Q₀ = G^[1] − G` **identically**
for generic band-3 coefficients (`K=3`). □

For `K=3` this is a finite ℂ-combination of shifted products of the coefficient
polynomials — no rational functions, no `(T+1)⁻¹`. It **specializes to band 2**:
with `K=2`, gauge `b₂=0`, `a₂ = h h⁽¹⁾`, `a₁ = h p`, `b₁ = κh`,
```
   G = h⁽⁻¹⁾( h w⁽¹⁾ + h⁽⁻²⁾ w + p⁽⁻¹⁾ v − κ u ),
```
exactly the potential of `../band2-square-sector/quantum-completion.md` §2
(`verify` §3). This is the general-`h` central telescoping that memo recorded
as having "no one-line closed form" — the closed form is the band-agnostic `G`
above.

**Point condition.** Under genuine membership, `G(0) = 0` **identically**
(`verify` §3): every term of `G(0)` carries a factor `b_{-k}(j)` or `a_{-k}(j)`
with `0 ≤ j ≤ k−1`, and `E^{underline k} | b_{-k}, a_{-k}` makes each such
value vanish. Therefore `Q₀ = 1` (i.e. `(T−1)G = 1`, forcing `G = E + c`) fixes
`c = 0`:
```
   G = E .                                                     (CENTRAL)
```
This is the band-3 form of the band-2 central integral
`w⁽¹⁾+w+p⁽⁻¹⁾v−κu = E` — the same mechanism (`W4`/`W6`), one band up.

## 5. The negative tail: bottom mirror and the cross-coupling of `λ₃, μ₃`

### 5.1 Bottom proportionality and bottom wall

Assume `a₋₃≠0`. By the same 3-periodicity argument reflected to the bottom,
`Q₋₆ = b₋₃^[−3]a₋₃ − a₋₃^[−3]b₋₃ = 0` gives
```
   b₋₃ = μ₃ a₋₃ ,     μ₃ ∈ ℂ                                   (verify §4).
```
If `a₋₃=0`, `Q₋₆=0` is vacuous and this `μ₃` gauge is unavailable; that is a
separate vanishing branch. On the nonzero branch, direct substitution in `Q₋₅`
shows that the correct gauged bottom variable is
```
   ϕ := b₋₂ − μ₃a₋₂,
   Q₋₅ = ϕ^[−3]a₋₃ − a₋₃^[−2]ϕ.                         (BOTTOM-WALL)
```
Indeed, expanding the right side gives the displayed `b₋₂` wall plus
`μ₃[a₋₃^[−2]a₋₂−a₋₂^[−3]a₋₃]`, fixing the minus sign in `ϕ`. The quantum Fourier
reflection `E↦−E−1` (`../band2-square-sector/quantum-a2-zero.md` §2, a genuine
automorphism, no sign flip) flips shift signs — `(P^[−3])|_φ=(P|_φ)^[+3]` and
`(P^[−2])|_φ=(P|_φ)^[+2]` (`verify` §4) — carrying (BOTTOM-WALL) to the top-wall
form. Hence, when `ϕ≠0`, **`a₋₃` and `ϕ` obey the same `Φ₃` necklace criterion**
as `a₃,b₂`, with mirror degree law `2 deg a₋₃=3 deg ϕ`. No such law for raw
`b₋₂` follows when `μ₃a₋₂` is present.

### 5.2 The first genuinely new phenomenon: `λ₃ ≠ μ₃`

At band 2 the single gauge `D ↦ D − λX` handled the top, and the bottom kept one
free proportionality constant (`w = μs`, `../band2-square-sector/quantum-mirror.md`).
On the two-sided branch `a₃a₋₃≠0`, at band 3 the same single gauge is spent on
the **top** (`b₃ = 0` via `λ₃`), so the **bottom** proportionality
`b₋₃ = μ₃ a₋₃` survives with `μ₃` **independent of `λ₃`**. The two constants
couple through the middle equations. In raw
`b₋₂`, `Q₋₅` is inhomogeneous with a `μ₃`-source, while the substitution
`ϕ=b₋₂−μ₃a₋₂` converts it exactly to the homogeneous (BOTTOM-WALL) equation.
This
is precisely the "cross-coupling of proportionality constants" the archived
milestone flagged as the first new band-3 effect (quoted in the band-2 theorem
residuals, `84978b9`): the top and bottom cannot both be gauged flat, and the
residue `μ₃` contaminates the bottom wall. There is no band-2 analogue.

### 5.3 Lemma-R bottom-wall balance: same degrees, different structure

The classical and quantum bottom walls have the same degree balance. For the
quantum staggered difference `f^[a]g^[b] − f^[a']g^[b']`, with `p = deg f` and
`q = deg g`, cancellation of the top terms leaves
```
   coeff(E^{p+q−1}) = ( (a−a')p + (b−b')q )·lc(f)·lc(g)        (verify §5).
```
The classical differential wall has the corresponding leading coefficient and
forces the same relation; for band `k` both give `(k−1)S=kΦ`. The distinction is
not a stronger quantum degree law: the classical differential equation yields
ordinary-power behavior, whereas the quantum shifted equation retains the
root-necklace/cyclotomic structure described in §3. Instances proved exactly
(`verify` §5), each firing **when the named term dominates**:

| equation / term | staggered shifts `(a,b)/(a',b')` | obstruction `(a−a')p+(b−b')q` | forced degree relation |
|---|---|---|---|
| top wall `Q₅`, `a₃^[0]b₂^[3] − a₃^[2]b₂^[0]` | `(0,3)/(2,0)` | `−2p+3q` | `2 deg a₃ = 3 deg b₂` |
| gauged bottom wall `Q₋₅`, `a₋₃^[0]ϕ^[−3] − a₋₃^[−2]ϕ^[0]` | `(0,−3)/(−2,0)` | `2p−3q` | `2 deg a₋₃ = 3 deg ϕ` (when `ϕ≠0`) |
| `Q₋₄`, `a₋₃^[0]b₋₁^[−3] − a₋₃^[−1]b₋₁^[0]` | `(0,−3)/(−1,0)` | `p−3q` | `deg a₋₃ = 3 deg b₋₁` |
| `Q₋₄`, `μ₃(a₋₃^[−1]a₋₁^[0] − a₋₃^[0]a₋₁^[−3])` | `(−1,0)/(0,−3)` | `−p+3q` | `deg a₋₃ = 3 deg a₋₁` |

and the middle piece of `Q₋₄` is the **J2/J2q "square classification one rung
down"** (`W2`): `b₋₂^[−2]a₋₂ − a₋₂^[−2]b₋₂`, whose isolated vanishing makes
`b₋₂/a₋₂` 2-periodic, hence `b₋₂ = ν a₋₂` (`verify` §5, the rational
2-periodicity engine). In the full system this piece is coupled to the
staggered `b₋₁`/`μ₃` terms, so the clean `ν` proportionality only emerges after
a case split — the residual work of the band-3 negative tail.

`verify` §5 records the exact `Q₋₄` decomposition
```
   Q₋₄ = [ b₋₁^[−3]a₋₃ − a₋₃^[−1]b₋₁ ]
        + μ₃·[ a₋₃^[−1]a₋₁ − a₋₁^[−3]a₋₃ ]
        + [ b₋₂^[−2]a₋₂ − a₋₂^[−2]b₋₂ ] .
```

## 6. Positive control

`verify` §6 exhibits a **genuine** band-3 pair, verifying the machinery accepts
real pairs:
```
   U = x + c₁∂ ,     X = U³ − ∂/κ ,     D = κU ,
   [D,X] = κ[U,X] = κ·(1/κ) = 1     (all Q_m = δ_{m0}, machine-checked),
```
with every membership holding, `a₃ = 1` (trivial shifted cube `h=1`), and
`a₋₃ = c₁³ E(E−1)(E−2)`. Its `b₂ = 0` (the wall solution with `κ₂ = 0`). A
genuine pair with **nonconstant** top would require the wall satisfied
nontrivially and the whole tail consistent — exactly the open question of §3.5.

## 7. Case tree (with band-2 reductions cited)

Writing the descent in the gauge `b₃ = 0`:

1. **`Q₆`** ⇒ `b₃ = λ₃ a₃`; gauge `b₃ = 0`. *(Proved, §2.)*
2. **`Q₅` wall** ⇒ either `b₂ = 0`, or `a₃` is `Φ₃`-compatible per coset and
   `b₂ = κ₂·g(a₃)` (1-dimensional). *(Proved, §3; shifted-cube **not** forced.)*
3. **`Q₋₆`** ⇒ if `a₋₃≠0`, then `b₋₃ = μ₃ a₋₃` (independent scalar);
   if `a₋₃=0`, the extreme equation is vacuous and routes to a separate
   vanishing branch. *(Proved, §5.1.)*
4. **`Q₄, Q₃, Q₂`** (positive middle) — inhomogeneous cascades coupling
   `b₁, b₀` to `(a₃, a₂, b₂)`. When `a₃` is a nonzero **constant** (trivial
   cube, normalized `a₃=1` by diagonal scaling), the wall `Q₅ = a₃(b₂^[3]−b₂)`
   forces `b₂` constant `=κ₂`, and `Q₄ = κ₂(a₂ − a₂^[2]) + (b₁^[3] − b₁)`
   (`verify` §7): the band-3 positive cascade, **structurally parallel** to band 2
   (`../band2-square-sector/quantum-completion.md` §1, `84978b9`) but with
   **3-fold** periodicities `b₁^[3]=b₁` replacing the 2-fold ones — the natural
   band-3 tame sector (the positive control of §6 lives here). For nonconstant
   `a₃` these are the genuinely new residual systems. *(Routed / open — this
   memo does not solve the positive middle.)*
5. **`Q₀`** ⇒ `G = E` (central integral). *(Proved, §4.)*
6. **`Q₋₁ … Q₋₄`** (negative middle) — the staggered rigidities of §5.3, plus
   the `μ₃`-cross-coupling. *(Partial; degree-forcing lemmas proved, closure
   open.)*
7. **`Q₋₅` bottom wall** ⇒ the pair `(a₋₃,ϕ)` with
   `ϕ=b₋₂−μ₃a₋₂` obeys the reflected `Φ₃` necklace law.
   *(Proved structure, §5; closure open.)*

**Band-2 reductions consumed.** A pair drops to band at most 2 only when its
level-3 extremes vanish (after any stated normalization), not merely when their
coefficient polynomials are constant. Such genuinely lower-band sectors are
closed by the assembled quantum band-2 theorem
(`../band2-square-sector/quantum-band2-theorem.md`, `84978b9`): the
shifted-square sector (`quantum-mirror.md` `ad43ab5` + `quantum-completion.md`),
the non-shifted-square kill (`quantum-M4.md`), the `a₂=0` reduction and Fourier
orientation (`quantum-a2-zero.md`), and band-1 rigidity (`quantum-band1.md`).
A genuine band-3 sector has, after orientation, a **nonzero level-3 extreme**;
a nonzero constant extreme remains genuine band 3 and includes the tame sector
above. Nonconstant extremes instead identify the wall-novel/open residual
behavior, where the `Φ₃`/sesqui wall and the `λ₃–μ₃` cross-coupling have no
band-2 shadow.

## 8. Status of claims

- **Proved, arbitrary degree (structural, machine-checked identities):**
  the `Q_m` system (§0/§1); `Q₆ ⇒ b₃ = λ₃a₃` (§2); the Wall Lemma and its exact
  cyclotomic reduction (WALL-M)/(B-FORMULA), the 1-dimensional freedom, the
  `2 deg a₃ = 3 deg b₂` degree law, shifted-cube **sufficiency** (§3.1–3.3); the
  **refutation** of shifted-cube necessity via `E(E−2)(E−4)` and the
  realizability gap `{1,5,6}` (§3.2, §3.4); the Telescoping Lemma `Q₀=(T−1)G`
  with the band-agnostic closed-form `G`, its band-2 specialization, and
  `G(0)=0 ⇒ G=E` (§4); on the branch `a₋₃≠0`, the bottom proportionality
  `b₋₃ = μ₃a₋₃`, the reflection shift-flip, the `Q₋₅`/`Q₋₄` decompositions, and
  the staggered leading-coefficient Lemma-R identities (§5).
- **Computed / corroboration only:** the bounded wall-admissibility sweep and
  the genuine positive-control witness (`verify` §6–§7).
- **Open (the frontier this cascade exposes):**
  1. whether the full band-3 system + membership forces `a₃` to a *trivial*
     (constant) cube — i.e. whether the wall's non-cube tops are killed
     downstream (the band-3 analogue of "constant `h`"); **not settled by `Q₅`**;
  2. the closure of the `μ₃`-cross-coupled negative tail (§5), hence any band-3
     tame classification;
  3. any global DC1/JC2 statement — **not claimed**; band 3 is not closed here.
- **No counterexample to DC1 is produced.** The `E(E−2)(E−4)` construction and
  other wall witnesses solve `Q₅` only. They are not asserted to satisfy the
  remaining `Q_m`, are not Weyl pairs, and refute only shifted-cube necessity for
  the isolated wall equation. The §6 positive control is separately verified as
  a genuine Weyl pair.

## 9. Verification

```sh
uv run --with sympy python research/band3/verify_quantum_band3.py
```
runs §0 (crossed-product engine; `Q_m` = commutator, all `m ∈ [-6,6]`; closed
forms), §1 (`Q₆` proportionality, period-3 rigidity, gauge legitimacy +
membership), §2 (the wall: shifted-cube sufficiency at symbolic degree; the
multiset reconstruction; the **counterexample** and shifted-cube refutation;
1-dimensionality; the degree law; the realizability failure), §3 (`Q₀=(T−1)G`,
the band-2 specialization, `G(0)=0`), §4 (bottom mirror, reflection,
cross-coupling), §5 (the Lemma-R staggered identities and the `Q₋₄`
decomposition), §6 (the genuine positive control), §7 (bounded corroboration).
A successful run ends `ALL QUANTUM BAND3 CHECKS PASSED`.

## 10. Relation to the classical band-3 sibling and the induction

The classical band-3 sibling works the `{G,F}=1`, `τ=xξ` face in parallel. Its
top `C₆` gives `b₃ = λ₃a₃` by `(b₃/a₃)' = 0`, and its `C₅` wall proves the
ordinary-power alternative: a nonzero wall variable forces `a₃=c·h³`. The
quantum `Q₅` wall has the **same degree balance** `2 deg a₃=3 deg b₂`, but its
shifted equation retains root-necklace/cyclotomic cofactors absent from the
classical differential equation. Thus `E(E−2)(E−4)` is a purely quantum
**isolated-wall** phenomenon: it solves `Q₅`, but it is not asserted to satisfy
the full cascade, genuine membership, or generation, and is not a Weyl pair or
a DC1 counterexample. Quantum is not stronger at the level of this degree law;
the solution structures differ (ordinary powers classically, necklaces
quantumly), exactly as in the corrected W6 comparison.

The band-3 floor of the width induction therefore stands, on the DC1 face, on:
top proportionality ✓, the wall structure ✓ (with its new `Φ₃`/sesqui form),
the central integral ✓, and an **open** negative tail whose new content is the
`λ₃–μ₃` cross-coupling. The lattice/Diophantine flavor of the classical mod-3
kill (`3V=3P+1`) reappears here as the `2 deg a₃ = 3 deg b₂` ratio and the
`Φ₃`-divisibility — the arithmetic the induction will have to carry.
