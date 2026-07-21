# The quantum band-3 shifted-cube sector: positive cascade, the S₃ midpoint, and h-constant forcing

**INDEPENDENTLY DERIVED AND MACHINE-VERIFIED — NOT PEER REVIEWED — BAND-SCOPED**

QUANTUM MIRROR assault, band 3, the **non-exotic** half of the wall dichotomy.
This memo works the shifted-cube sector

```text
gauge b_3 = 0,   a_3 = h(E) h(E+1) h(E+2) != 0   (a genuine shifted cube),
```

the *sufficient* branch of the `Q₅` Wall Lemma (`quantum-band3-cascade.md` §3.2,
commit `050a4c0`). It is independent of the exotic (non-cube) `Φ₃`-necklace tops
`E(E−2)(E−4)`, W1, W2 handled by the sibling exotic memos
(`quantum-exotic-*.md`, `quantum-w1-arbitrary-degree.md`, commit `e4e704f`). It
pushes Wave A/B of the cascade memo: the positive-middle cascades `Q₄,Q₃,Q₂` and
the h-constant question.

Everything below is checked exactly by
[`verify_quantum_shifted_cube.py`](verify_quantum_shifted_cube.py) (ends
`ALL QUANTUM SHIFTED CUBE CHECKS PASSED`).

Conventions are frozen from the cascade memo:
`A₁[x⁻¹]=⊕_k x^k ℂ[E]`, `E=x∂`, `(x^a f)(x^b g)=x^{a+b} f(E+b) g(E)`,
`f^[r](E)=f(E+r)`, `T f = f^[1]`,
`Q_m = Σ_{k+l=m}(b_l^[k] a_k − a_k^[l] b_l)`, `[D,X]=1 ⇔ Q_m=δ_{m0}`,
genuine membership `E^{underline r}=E(E−1)⋯(E−r+1) | a_{−r}, b_{−r}`,
`Q₀=(T−1)G` with the closed-form potential
`G = Σ_{k=1}^3 Σ_{j=0}^{k−1}(a_k^[j−k] b_{−k}^[j] − b_k^[j−k] a_{−k}^[j])`.

Write `S_n = 1 + T + ⋯ + T^{n−1}`, so `T^n − 1 = (T−1)S_n`.

## 0. Headline

> Band 2 forced `h` constant in the shifted-**square** sector by a general-`h`
> central telescoping `G = h^[-1](h w^[1] + h^[-2] w + p^[-1] v − κu) = E` plus a
> one-line membership evaluation (`quantum-completion.md` §5, commit `b9f9cf3`).
> **The band-3 shifted-cube sector reproduces exactly this mechanism, one weight
> level up, with `S₂ = 1+T` replaced by the quantum midpoint `S₃ = 1+T+T²`.** The
> positive cascade forces the divisibilities `h h^[1] | a₂`, `h | a₁`, `h | b₁`
> (the band-3 analogues of band-2's `h | a₁`); the central telescoping then
> factors as `G = h^[-1] M` with a *closed form* `M`; and `Q₀ = 1` forces
> `h^[-1] | E`, hence `deg h ≤ 1`, with `deg h = 1` killed by `M(0) = 0` versus
> `M = 1/α`. **For cube-separated `h` and every `κ`, `h` is forced constant.**
> The constant-`h` sector contains the tame family
> `X = U³ − ∂/κ − A`, `D = λX + κU + β`.

The one genuinely new obstacle relative to band 2 is that the divisibility step
uses the coprimalities `gcd(h, h^[j]) = 1` for `j=1,2,3` (**cube-separation**);
tops whose roots differ by `1,2,3` are a documented residual, killed at bounded
degree but not by the arbitrary-degree argument here.

## 1. The sector and the wall

Assume `a₃ = h h^[1] h^[2] ≠ 0` and the top gauge `b₃ = 0`. The `Q₅` wall
`b₂^[3] a₃ = a₃^[2] b₂` is solved, for a shifted cube, by

```text
   b₂ = κ · h h^[1],           κ ∈ ℂ,
```

and the wall solution space is exactly `1`-dimensional (the Wall Lemma;
`verify` §1 checks `b₂ = κ h h^[1]` at symbolic `deg h = 1,2,3` and the
`1`-dimensionality at concrete cubes). This is the defining shifted-cube reduction
`κ·h³ ⟹ κh²` deformed into `κ·h h^[1] h^[2] ⟹ κ·h h^[1]`
(`band3-tame-catalog.md` §6, commit `050a4c0`).

## 2. The positive cascade and the quantum S₃ midpoint

In the gauge `b₃=0` with `b₂ = κ h h^[1]`, the three positive-middle equations
each carry a *staggered* operator on the next unknown `b₁, b₀, b₋₁`, all staggered
against `a₃ = h h^[1] h^[2]`:

```text
   Q₄ :  b₁^[3] a₃ − a₃^[1] b₁   (stagger +3 vs +1)  =  h^[1]h^[2]( b₁^[3] h − h^[3] b₁ ),
   Q₃ :  b₀^[3] a₃ − a₃    b₀   (stagger +3 vs  0)  =  a₃ ( b₀^[3] − b₀ ),
   Q₂ :  b₋₁^[3] a₃ − a₃^[-1] b₋₁ (stagger +3 vs −1) =  h h^[1]( b₋₁^[3] h^[2] − h^[-1] b₋₁ ).
```

(`verify` §2 checks the `Q₂` operator factorization exactly.)

### 2.1 `Q₄ → b₁`: the divisibilities `h h^[1] | a₂`, `h | b₁`, and the `S₃` midpoint

`Q₄ = 0` reads, in gauge,

```text
   h^[1]h^[2]( b₁^[3] h − h^[3] b₁ ) = κ ( h h^[1] a₂^[2] − h^[2]h^[3] a₂ ).      (Q4)
```

The left side is divisible by `h^[1]h^[2]`; therefore so is the right. Using
cube-separation `gcd(h,h^[1])=gcd(h,h^[2])=gcd(h,h^[3])=1`:

- `h^[2] |` RHS and `h^[2] | κ h^[2]h^[3] a₂` give `h^[2] | κ h h^[1] a₂^[2]`,
  hence `h^[2] | a₂^[2]`, i.e. **`h | a₂`**;
- `h^[1] |` RHS and `h^[1] | κ h h^[1] a₂^[2]` give `h^[1] | κ h^[2]h^[3] a₂`,
  hence **`h^[1] | a₂`**.

So `h h^[1] | a₂`; write `a₂ = h h^[1] g`. Substituting into (Q4) and cancelling
`h^[1]h^[2]` gives `b₁^[3] h − h^[3] b₁ = κ h h^[3](g − g^[2])`, so
`h^[3] b₁ = h( b₁^[3] − κ h^[3](g−g^[2]) )` and, by `gcd(h,h^[3])=1`,
**`h | b₁`**; write `b₁ = h β`. Then (`verify` §2)

```text
   Q₄ = a₃ h^[3] [ (T³−1) β − κ (T²−1) g ].
```

Because `T³−1 = (T−1)S₃` and `T²−1 = (T−1)S₂`, `Q₄ = 0` is the **quantum band-3
midpoint**

```text
   S₃ β = κ S₂ g + γ₁,     i.e.   β^[2] + β^[1] + β = κ(g^[1] + g) + γ₁,       (MID)
```

with a free constant `γ₁`. This is band 2's midpoint `B^[1] + B = κp + γ`
(`quantum-completion.md` §2) with the `2`-fold sum `S₂` promoted to the `3`-fold
sum `S₃` — the exact deformation the cascade memo predicted (§7, "3-fold
periodicities replacing the 2-fold ones").

### 2.2 `Q₃ → b₀`: the divisibility `h | a₁`

With `a₂ = h h^[1] g`, `b₁ = h β`, `b₂ = κ h h^[1]`, the equation `Q₃ = 0` reads
`a₃(b₀^[3]−b₀) + h h^[1]h^[2](β^[2]g − g^[1]β) + κ h^[1](h^[2]a₁ − h a₁^[2]) = 0`.
The first two terms carry the full factor `a₃ = h h^[1] h^[2]`, so
`a₃ | κ h^[1](h^[2]a₁ − h a₁^[2])`, i.e. `h h^[2] | κ(h^[2]a₁ − h a₁^[2])`;
`h | κ h^[2] a₁` with `gcd(h,h^[2])=1` gives **`h | a₁`**. Writing `a₁ = h p`,

```text
   Q₃ = a₃ [ (T³−1) b₀ + (β^[2] g − g^[1] β) + κ (p − p^[2]) ].
```

Here `T³−1 = (T−1)S₃` is surjective on `ℂ[E]` with kernel the constants (indeed
`S₃` is a linear automorphism, `../dc1-program/two-filler-cross-cancellation.md`,
commit `0889f8a`), so `b₀` **always exists** and is determined up to one additive
constant — precisely the "constant `b₀` freedom" of the reconstruction
(`quantum-exotic-cokernel.md` §1, commit `e4e704f`). No divisibility beyond
`h | a₁` is imposed by `Q₃`.

### 2.3 `Q₂ → b₋₁`

`Q₂ = 0` has leading operator `h h^[1]( b₋₁^[3] h^[2] − h^[-1] b₋₁ )` on the first
negative unknown `b₋₁` (membership `E | b₋₁`), inhomogeneously coupled to `a₀`,
the solved `b₀,b₁`, and the constant-`b₀` freedom. It is not needed for the
h-constant forcing below (which uses only `Q₄, Q₃, Q₀`); for the `κ = 0` route it
supplies the last divisibility (§6).

### 2.4 The three band-3 divisibilities

The clean output of the positive cascade is the triple

```text
   h h^[1] | a₂  (a₂ = h h^[1] g),    h | a₁  (a₁ = h p),    h | b₁  (b₁ = h β),
```

the band-3 analogue of band-2's single `h | a₁`. `verify` §2b confirms the
*necessity* of `h h^[1] | a₂` and `h | b₁` at the concrete cube-separated top
`h = E(2E−1)` (roots `{0, 1/2}`): a fully generic `a₂` admits a polynomial `b₁`
solving `Q₄ = 0` **iff** `h h^[1] | a₂`.

## 3. The general-h central identity `G = h^[-1] M`

Substituting the three divisibilities into the closed-form potential collapses
every one of the nine terms of `G` to a common left factor `h^[-1]`:

```text
   G = h^[-1] · M,
   M =  p^[-1] b₋₁ − β^[-1] a₋₁
      + h^[-2] g^[-2] b₋₂ + h g^[-1] b₋₂^[1]
      + h^[-3]h^[-2] b₋₃ + h^[-2]h b₋₃^[1] + h h^[1] b₋₃^[2]
      − κ h^[-2] a₋₂ − κ h a₋₂^[1].
```

This is a **term-by-term** factorization, hence degree-independent: e.g.
`a₁^[-1] b₋₁ = (h p)^[-1] b₋₁ = h^[-1] p^[-1] b₋₁`, and likewise for all nine
terms. `verify` §3 checks `G = h^[-1] M` **identically** at `deg h = 1, 2, 3`
(cube-separated instances of the same manifest factoring). The `a₃`- and
`b₂`-blocks factor automatically from the cube/wall structure; the residual block
needs exactly the divisibilities of §2.4. This is the band-3 form of the band-2
identity `G = h^[-1](h w^[1] + h^[-2] w + p^[-1] v − κ u)`
(`quantum-completion.md` §2).

**Membership evaluation.** Every term of `M` carries a negative coefficient
`a₋_k` or `b₋_k` evaluated at one of `0,1,…,k−1`, all of which vanish by
`E^{underline k} | a₋_k, b₋_k`. Hence

```text
   M(0) = 0        (identically under membership).                            (M0)
```

`verify` §4 checks `M(0)=0` at symbolic `h`, and the companion `G(0)=0` for fully
generic band-3 data.

## 4. h-constant forcing (the affine kill)

**Theorem (shifted-cube h-forcing).** Let `char F = 0`, `a₃ = h h^[1] h^[2] ≠ 0`
with `h` **cube-separated** (`gcd(h,h^[j])=1`, `j=1,2,3`), gauge `b₃=0`, genuine
membership. If polynomial data satisfy the positive cascade and `Q₀ = 1`, then
`h` is constant.

*Proof.* By §2 the divisibilities hold (for `κ ≠ 0` from `Q₄,Q₃`; for `κ = 0`
from `Q₄,Q₃,Q₂`, §6), so `G = h^[-1] M` with `M` a polynomial. Now
`Q₀ = (T−1)G = 1 = (T−1)E`, so `G − E` is `1`-periodic; over characteristic zero
the leading term of `(T−1)q` is `deg(q)·lc(q)·E^{deg q −1}`, nonzero for
`deg q > 0`, so a `1`-periodic polynomial is constant and `G = E + c`. Membership
gives `G(0) = h^[-1](0)M(0) = 0` (indeed `G(0)=0` for all band-3 data), so `c = 0`
and

```text
   G = h^[-1] M = E.
```

Thus `h^[-1] | E`, a degree-`1` polynomial. If `deg h ≥ 2` then `deg h^[-1] ≥ 2`
and `h^[-1] ∤ E` — contradiction. If `deg h = 1`, write `h = α(E−ρ)`,
`h^[-1] = α(E−1−ρ)`; `h^[-1] | E` forces the root `1+ρ = 0`, so `ρ = −1`,
`h^[-1] = α E`, and `α E · M = E` gives `M = 1/α` (constant). But `M(0) = 0` by
(M0), so `1/α = 0` — impossible. Hence `deg h = 0`: `h` is constant. ∎

`verify` §5 machine-checks each step: the arbitrary-`n` leading coefficient `n`
of `(E+1)^n − E^n`; that `h^[-1] = αE` divides `E` while `α(E−1−ρ)` with `ρ≠−1`
does not; that `αE·(1/α) = E`; and that a degree-`2` `h^[-1]` does not divide `E`.
Together with (M0) and `G = h^[-1] M` these certify the theorem.

**Remark (why it is cleaner than band 2).** Band 2 evaluated
`κu = h w^[1] + h^[-2] w + p^[-1] v − 1/α` at `E=0`, where the single term `−κu`
survived to give `κu(0) = −1/α` contradicting `u(0)=0`. In band 3 **all nine**
terms of `M` already vanish at `E=0`, so the contradiction is the bare
`M(0) = 0 = 1/α`.

## 5. Corroboration of the forcing

Every *cube-separated* nonconstant top (any degree, any `κ`) is excluded by the
theorem of §4 itself — no computation needed: `deg h = 1` tops such as `E+1`,
`2E−1` are killed by the `M(0)=0` versus `1/α` contradiction, and
`deg h ≥ 2` separated tops such as `E(2E−1)` by the degree bound `h^[-1] ∤ E`.
The bounded Gröbner evidence in `verify` §8 is therefore reserved for the
**non-cube-separated** residual, which the theorem does *not* reach (§7).

## 6. The `κ = 0` branch folds in

If `κ = 0` then `b₂ = 0` and `Q₄` becomes `b₁^[3] a₃ = a₃^[1] b₁`, i.e.
`(b₁/h)^[3] = b₁/h`: the ratio is `3`-periodic, hence constant, so **`b₁ = c h`**
(`verify` §6). Then `Q₃ = 0` forces `h h^[1] | a₂` and `Q₂ = 0` forces `h | a₁`
(both checked in `verify` §6). These are exactly the divisibilities of §2.4 with
`β = c`, `κ = 0`, so `G = h^[-1] M` still holds and the affine kill of §4 applies
verbatim. **The `κ=0` branch also forces `h` constant** for cube-separated `h`;
there is no separate surviving `κ=0` sector (contrast band 2, where `κ=0` was a
distinct — and empty — branch).

## 7. Exceptional locus: non-cube-separated tops

The divisibility derivation of §2 uses `gcd(h,h^[j]) = 1` for `j=1,2,3`. If two
roots of `h` differ by `1,2,3` (necessarily `deg h ≥ 2`), it genuinely fails: for
`h = E(E−1)` (roots `0,1`), `verify` §7b exhibits a generic `a₂` **not** divisible
by `h h^[1]` for which `Q₄ = 0` still has a polynomial `b₁` (the divisibility is
false there). So the closed-form factorization `G = h^[-1] M` is not available and
the arbitrary-degree proof does **not** cover these tops.

They remain a **residual**, but the bounded evidence is uniformly negative: exact
Gröbner elimination at cap `D=2` (`verify` §8) finds the positive cascade `+ Q₀=1`
`+` membership **inconsistent** (unit ideal) for

```text
   h = E(E−1)  (roots differ by 1),   h = E(E−2)  (by 2),
   h = E(E−3)  (by 3),                h = (E−1)(E−3)  (by 2).
```

This is corroboration only; it is **not** an arbitrary-degree exclusion for the
integer-spaced-root tops. It is the single explicit gap between this memo and a
full shifted-cube closure.

## 8. The constant-h sector and the tame family

When `h` is constant (normalize `h = 1`, `a₃ = 1` by diagonal scaling), the wall
gives `b₂ = κ₂` constant and the positive cascade is the `3`-fold-periodic tame
system, e.g.
`Q₄ = κ₂(a₂ − a₂^[2]) + (b₁^[3] − b₁)` (the `S₃` midpoint; `verify` §7, matching
the cascade memo `verify` §7). The sector contains the explicit **tame family**

```text
   U = x + c₀ + c₁∂,     X = U³ − (1/κ)∂ − A,     D = λX + κU + β,
   c₀,c₁,A,κ,λ,β ∈ ℂ,   κ ≠ 0,
```

with `[U,X] = 1/κ`, hence `[D,X] = κ[U,X] = 1` (`verify` §7, all parameters). Its
reduced data are `a₃ = 1`, `a₋₃ = c₁³ E(E−1)(E−2)` (an exact falling factorial,
genuine membership with no slack), and after the gauge `D ↦ D − λX` the entire top
tower `b₃, b₂, b₋₃` vanishes at once — the **B0-band3 collapse**
(`band3-tame-catalog.md` §6, Q3-2; commit `050a4c0`). The cascade memo's §6
positive control `X = U³ − ∂/κ`, `D = κU` is the `c₀=A=λ=β=0` sub-case (`verify`
§7). Note the gauged wall constant is `κ₂ = 0`: the tame family lives in the
gauged-`b₂ = 0` branch, the band-3 image of band 2's gauged-`b₂ = 0` shifted
square.

**Scope of the constant-h classification.** This memo proves the tame family is a
genuine sub-sector and pins the positive-cascade structure. It does **not** prove
the constant-`h` sector is *exactly* the tame family: that requires the closure of
the `μ₃`-cross-coupled negative tail and the exclusion of the resistant locus
`A*-band3` (`quantum-band3-cascade.md` §5, §7; `band3-tame-catalog.md` §7,
conjecturally empty), which are open one level down as well. The bounded probe
for constant-`h` pairs with `κ₂ ≠ 0` was inconclusive (Gröbner too costly at the
attempted caps); no claim is made there.

## 9. Honest ledger — where this sits in the band-3 rung

The full quantum band-3 rung on the DC1 face decomposes as:

| piece | who supplies it | status |
|---|---|---|
| `Q₆ ⇒ b₃ = λ₃a₃`, gauge `b₃=0` | cascade memo `050a4c0` | proved |
| `Q₅` Wall Lemma: `b₂=0` or `a₃` is `Φ₃`-compatible | cascade memo `050a4c0` | proved |
| **shifted-cube tops `a₃ = c·h h^[1] h^[2]`** | **this memo** | **h constant forced (cube-separated, any κ); tame family** |
| exotic (non-cube) tops `E(E−2)(E−4)`, W1, W2 | sibling exotic memos `e4e704f` | W1 normalized datum killed; others scoped |
| `a₃ = 0` orientation ⇒ band ≤ 2 | Fourier/`a₂=0` routing → band-2 theorem `84978b9` | proved (band 2) |
| negative-tail / `μ₃` cross-coupling closure | open | open |

**This memo supplies, and only supplies:**
- the shifted-cube positive cascade `Q₄→b₁`, `Q₃→b₀`, `Q₂→b₋₁`, with the quantum
  `S₃` midpoint and the three divisibilities `h h^[1] | a₂`, `h | a₁`, `h | b₁`
  (arbitrary degree, machine-checked identities; necessity for cube-separated `h`);
- the general-`h` central identity `G = h^[-1] M` and, via `Q₀=1`, the forcing of
  `h` **constant** for cube-separated `h` and every `κ` (arbitrary degree);
- the constant-`h` tame family as a genuine sub-sector, with the `B0-band3`
  gauge collapse.

**It does not supply:** the exotic-branch closure; the `a₃ = 0` reduction; an
arbitrary-degree exclusion for non-cube-separated tops (only bounded evidence);
the complete constant-`h` classification (negative tail); or any DC1/JC2 claim.
No Weyl pair beyond the tame family, and no counterexample, is constructed.

## 10. Status of claims

- **Proved, arbitrary degree (machine-checked identities + coprimality/degree
  prose):** the wall solution `b₂ = κ h h^[1]` and its `1`-dimensionality; the
  `Q₄,Q₃,Q₂` operator forms and the `S₃` midpoint; the divisibilities
  `h h^[1]|a₂, h|a₁, h|b₁` for **cube-separated** `h` (any `κ`, including the
  `κ=0` fold-in); the central identity `G = h^[-1] M`; `M(0)=0` and `G(0)=0`;
  and the affine kill forcing `h` constant.
- **Exact symbolic checks:** all displayed identities; the tame family
  `[D,X]=1` with its coefficients and memberships; the `B0-band3` gauge collapse;
  the section-6 positive control.
- **Bounded computation (corroboration only):** unit-ideal emptiness at cap `D=2`
  for the non-cube-separated tops `E(E−1), E(E−2), E(E−3), (E−1)(E−3)` (the
  residual the arbitrary-degree theorem does not reach).
- **Open / not claimed:** arbitrary-degree exclusion for non-cube-separated tops;
  the full constant-`h` classification; the exotic branch; `a₃=0` reduction;
  DC1/JC2.

## 11. Verification

```sh
uv run --with sympy python research/band3/verify_quantum_shifted_cube.py
```

runs §0 (`Q_m` = direct commutator, sector data), §1 (wall + `1`-dim), §2
(cascade operators, `S₃` midpoint, `Q₂` operator, divisibility necessity), §3
(`G = h^[-1] M` at `deg h = 1,2`), §4 (`M(0)=0`, `G(0)=0`), §5 (the affine kill
pieces), §6 (`κ=0` fold-in), §7 (tame family, memberships, gauge collapse,
positive control, `S₃` midpoint), §7b (divisibility failure for the non-separated
`E(E−1)`), §8 (bounded emptiness). A successful run ends
`ALL QUANTUM SHIFTED CUBE CHECKS PASSED`.
