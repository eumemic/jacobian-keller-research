# Closure of the classical resistant branch A\*, and the complete constant-h square sector

**INDEPENDENTLY DERIVED — MACHINE-CHECKED IDENTITIES — NOT PEER REVIEWED**

This memo closes the residual branch left open by
`research/band2-m5-partial/classical-square-sector-partial.md` (cited below as
[UP]; all [UP] section and equation numbers refer to the revision at commit
`91a053a` — the file was subsequently rewritten at commit `44c66d5` into a
complete direct-substitution proof of the same sector, derived in parallel and
independently of this memo; the two arguments share the cascade but close the
resistant branch by different mechanisms, and their conclusions agree), namely
branch A\* (both negative extremes surviving), and with it
completes the classification of the **constant-h** classical square sector.
Combined with [UP]'s arbitrary-degree cascade and a one-line strengthening of
its nonconstant-h analysis (Section 6 below), the result is:

> **Theorem (classical square sector, complete).** Let `F, G ∈ C[x,ξ]` with
> `{G,F} = 1`, both supported in x-levels `[-2,2]` in the `τ = xξ`
> presentation, and `a₂ ≠ 0` (equivalently, by the proved M4 theorem, `a₂ = h²`
> a nonzero square). Then, after the normalizations of [UP] (subtract the top
> multiple of `F` from `G`; diagonal scaling; additive constants), the pair is
> a member of the tame family [UP,(5.1)] with constant `h` and `e = 0`.
> In particular every such pair is an automorphism pair: the plane map
> `(F,G)` is a tame polynomial automorphism of `C²`.

All bracket conventions, the nine coefficient equations `C_m`, and the
membership divisibilities (`τ | a₋₁, b₋₁`, `τ² | a₋₂, b₋₂`) are those of
`research/band2-classical-proved/M4_proof_memo.md`. Everything below is over
`C`. The verifier `verify_classical_Astar.py` (same directory) machine-checks
every displayed identity and runs regression sweeps; a successful run ends
`ALL CLASSICAL ASTAR CHECKS PASSED`.

## 1. Setup and normalizations

By [UP] §7, every genuine polynomial square-sector pair has `κ ≠ 0`. By [UP]
§3, in the gauge `b₂ = 0` with `b₁ = κh`:

```
a₁ = h p,   b₀ = (κ/2) p + β,   a₀ = p²/4 + 2hv/κ − A,
a₋₁ = pv/κ + 2hw/κ − (t+e)/(κh),   h | t+e,
```

with `v = b₋₁`, `w = b₋₂`, `s = a₋₂`. Section 6 below disposes of nonconstant
`h`; from here on `h` is constant, normalized `h = 1` by the diagonal
symplectic scaling `(x,ξ) → (ρx, ρ⁻¹ξ)` (which multiplies `a₂ = h²` by `ρ²`),
and `β = 0` by adding a constant to `G` (only `b₀'` occurs in the equations).

Membership `τ | a₋₁` evaluated at `t = 0`, using `v(0) = w(0) = 0`
(memberships on `b₋₁, b₋₂`), gives `a₋₁(0) = −e/κ`, hence

```
e = 0.
```

With these substitutions the equations `C₄, C₃, C₂, C₁` vanish identically and
`C₀ = 1` identically (machine-checked); the surviving system is
`C₋₁ = C₋₂ = C₋₃ = C₋₄ = 0` in the unknowns `p, v, w, s ∈ C[t]` with scalars
`κ ≠ 0, A`, and memberships `t | v`, `t² | w`, `t² | s`.

Write `P₁ = ∫₀ᵗ p`, `P₂ = ∫₀ᵗ p²`, `V₁ = ∫₀ᵗ v`, and `P = deg p`,
`V = deg v`, `W = deg w` (degrees of nonzero polynomials).

## 2. Two exact first integrals (all branches at once)

**Lemma 2.1.** With `s` free (no proportionality assumption), the following
are polynomial identities in the parametrized unknowns (machine-checked):

```
Φ  := κpw + v² + (κ/2)tp − κ²s − (κ/2)P₁        satisfies  Φ' = κ C₋₁ ;
I₂ := (κ/2)pP₁ − (κ/4)tp² + 2vw + tv − (κ/4)P₂ − 2V₁
                                                 satisfies  I₂' = κ C₋₂ − p'Φ .
```

Consequently, on any solution, `Φ` is the constant `Φ(0) = −κ²s(0) = 0`
(membership `t²|s`), and then `I₂` is the constant `I₂(0) = 0`. Thus

```
(Φ = 0):   κ²s = κpw + v² + (κ/2)(tp − P₁),
(I₂ = 0):  (κ/2)pP₁ − (κ/4)tp² + 2vw + tv = (κ/4)P₂ + 2V₁,
```

and the system reduces to these two identities together with `C₋₃ = C₋₄ = 0`.
Note `(Φ=0)` **determines `s`** from `(p,v,w)`; the conditions `s(0)=s'(0)=0`
are automatic from the memberships on `v, w`.

Both integrals degenerate helpfully at the top: for `P ≥ 1` the polynomial
`(κ/2)(tp − P₁)` has degree exactly `P+1` with leading coefficient
`(κ/2)·P/(P+1)·lc(p) ≠ 0`, and the p-block
`(κ/2)pP₁ − (κ/4)tp² − (κ/4)P₂` of `I₂` has degree exactly `2P+1` with leading
coefficient `−κ·P²/(2(P+1)(2P+1))·lc(p)² ≠ 0` (both machine-checked; the
second follows from `2(2P+1) − (P+1)(2P+1) − (P+1) = −2P²`).

Branches, following [UP]: **B0** (`w = 0`), **A0** (`w ≠ 0, s = 0`),
**A\*** (`w ≠ 0, s = μ̃w, μ̃ ≠ 0`; the proportionality is `C₋₄`).

## 3. Branch A\* is empty

Assume `w ≠ 0`, `s = μ̃ w`, `μ̃ ≠ 0`. Memberships give `W ≥ 2`, and `V ≥ 1`
when `v ≠ 0`.

**Case `p` constant, `p = κμ̃`.** Here `(Φ=0)` reads `v² = 0` (the `w` terms
cancel and `tp − P₁ = 0`), so `v = 0`; the system collapses to the single ODE
`(2w + t)w' = 2w`. For `W ≥ 2`, `deg(2ww') = 2W−1 > W = deg(2w − tw')` and
`lc(2ww') = 2W·lc(w)² ≠ 0`: no cancellation is possible, so `w = 0`,
a contradiction. (Regression: exhaustive solve for `W ≤ 8` finds only `0`.)

**Case `p` constant, `p = c ≠ κμ̃`.** `(Φ=0)` gives
`w = −v²/(κ(c−κμ̃))`, so `v = 0` would force `w = 0`; assume `v ≠ 0`. The
p-block of `(I₂=0)` vanishes identically for constant `p`, leaving
`2vw + tv = 2V₁`. Differentiating: `tv' − v = 6v²v'/(κ(c−κμ̃))`. If `V ≥ 2`
the right side has degree `3V−1 > V` and leading coefficient
`6V·lc(v)³/(κ(c−κμ̃)) ≠ 0`, impossible; if `V = 1`, the left side vanishes
identically (`v = v₁t`), forcing `v₁ = 0`. Contradiction.

**Case `P ≥ 1, v = 0`.** `(Φ=0)` reads `κ(p−κμ̃)w = (κ/2)(P₁ − tp)`; the left
side has degree `P + W ≥ P + 2`, the right side degree exactly `P+1`.
Contradiction.

**Case `P ≥ 1, v ≠ 0`.** Two forced degree balances:

*From `(Φ=0)` with `s = μ̃w`:* the term degrees are `P+W` (coefficient
`κ·lc(p)lc(w) ≠ 0`), `2V`, `W`, `P+1`. Since `P+W` strictly exceeds `W` and
`P+1`, the only possible partner for it is `v²`:

```
2V = P + W,      lc(v)² = −κ·lc(p)·lc(w).
```

*From `(I₂=0)`:* the v-block `2vw + tv − 2V₁` has top degree exactly `V+W`
(coefficient `2·lc(v)lc(w) ≠ 0`, since `V+W > V+1`), and the p-block has top
degree exactly `2P+1` with nonzero coefficient (Lemma 2.1). These two blocks
must cancel each other:

```
V + W = 2P + 1.
```

Substituting `W = 2P+1−V` into `2V = P+W` gives `3V = 3P + 1`, which is
impossible in integers. ∎

**Theorem 3.1.** Branch A\* contains no polynomial pair. (Note the proof never
uses `C₋₃`: the branch dies of `Φ`, `I₂`, and memberships alone.)

## 4. Branch A0 is empty

Assume `w ≠ 0`, `s = 0`. Equation `C₋₄` is vacuous; `C₋₃` reads
`2a₋₁'w − a₋₁w' = 0`, i.e. `(a₋₁²/w)' = 0`, so `a₋₁² = c·w`, `c ∈ C`.

**Subcase `a₋₁ = 0`** (`c = 0` collapses here unless `w=0`): then
`pv + 2w − t = 0`, so `w = (t−pv)/2`, and membership `t²|w` forces
`p(0)v'(0) = 1`. `(Φ=0)` becomes `v² − (κ/2)p²v + (κ/2)(2tp − P₁) = 0`.
For `p = 0` this gives `v = 0`, then `w = t/2`, violating `t²|w`. For `p`
a nonzero constant, `t ∈ (v)` forces `v = v₁t`, and the `t²` coefficient
gives `v₁ = 0`, then `(κ/2)pt = 0`, contradiction. For `P ≥ 1`: comparing
top degrees in the displayed quadratic (`2V` vs `2P+V` vs `P+1`) forces
`V = 2P` with `lc(v) = (κ/2)lc(p)²`; then in `(I₂=0)` the v-block has top
degree `V + W = V + (P+V) = 5P` (here `W = deg(t−pv) = P+V = 3P`), while the
p-block has top degree `2P+1 < 5P`; the v-block leading coefficient
`2·lc(v)lc(w) ≠ 0` is unmatched. Contradiction.

**Subcase `a₋₁ ≠ 0`, `c ≠ 0`, `w = a₋₁²/c`.** Let `L = deg a₋₁ ≥ 1`
(membership `t|a₋₁`), so `W = 2L`. In `κa₋₁ = pv + 2w − t`, the right side's
top candidates are `P+V` and `2L`; since the left side has degree `L < 2L`,
the tops must cancel: `P + V = 2L` (in particular `v ≠ 0` or `p v = 0` makes
this subcase collapse to the previous ones; for `v = 0` see below). Then in
`(Φ=0)` (with `s = 0`) the term degrees are `P+2L = 2P+V`, `2V`, `P+1`.
If `V < 2P`, the top `2P+V` has coefficient `κ·lc(p)lc(w) ≠ 0` and no partner
(`2P+V = P+1` forces `P+V=1`, impossible for `P,V ≥ 1`). If `V > 2P`, the
top `2V` needs `2V = P+1`, i.e. `V = (P+1)/2 > 2P`, forcing `P < 1/3`,
impossible. So `V = 2P`, and the two candidates share degree `4P`; then in
`(I₂=0)` the v-block top degree is `V + W = V + P + V = 5P` against the
p-block's `2P+1 < 5P` for `P ≥ 1`: unmatched, contradiction. For `P = 0`
(`p` constant): `(Φ=0)` gives `pw = −v²/κ`; `p = 0` forces `v = 0`, and
`p ≠ 0` gives `w = −v²/(κp)`, whence `a₋₁ = γv` and
`(κγ−p)v − 2v²/(κp) − t = 0`, whose degree-`2V` top is unmatched for `V ≥ 1`.

**Subcase `v = 0`.** `(I₂=0)` reduces to its p-block, forcing `p` constant;
`(Φ=0)` then gives `pw = 0`, so `p = 0`; then `κa₋₁ = 2w − t` and
`a₋₁² = cw` give `(2w−t)² = κ²cw`, impossible by degrees for `W ≥ 2`
(and `c = 0` gives `w = t/2`, violating `t²|w`). ∎

**Theorem 4.1.** Branch A0 contains no polynomial pair.

## 5. Branch B0 is exactly the tame family

Assume `w = 0`. `C₋₄` is vacuous; `C₋₃` reads `s'v − 2sv' = 0`.

If `v = 0`: `(I₂=0)` reduces to its p-block, forcing `p` constant, and then
`(Φ=0)` gives `κ²s = 0`. The resulting pair is
`F = (x + p/2)² − A − ξ/κ`, `G = κ(x + p/2)`, visibly a tame automorphism
pair (solve `G` for `x + p/2`, then `F` for `ξ`).

If `v ≠ 0`: `(s/v²)' = 0` gives `s = c v²`, `c ∈ C`, and `(Φ=0)` reads
`(κ²c − 1)v² = (κ/2)(tp − P₁)`.

*If `c = 1/κ²`:* `tp = P₁`, so `tp' = 0`, `p` constant; the p-block of
`(I₂=0)` vanishes and `tv = 2V₁` forces `v = v₁t`; then `s = c v₁² t²`
satisfies `t²|s` automatically. These are precisely the constant-h members of
the family [UP,(5.1)] with `e = 0` (parameters `c₀ = p/2`, `c₁ = v₁/κ`):
tame, as proved in [UP] §5.

*If `c ≠ 1/κ²`:* then `v² = (κ/2)(tp−P₁)/(κ²c−1)`. For `p` constant the right
side vanishes, forcing `v = 0`, done above. For `P ≥ 1`: degrees give
`2V = P+1`; in `(I₂=0)` (with `w = 0`) the v-block is `tv − 2V₁`, of top
degree `V+1` with coefficient `(V−1)/(V+1)·lc(v)`, which vanishes only for
`V = 1`. For `V ≥ 2` the balance against the p-block forces `V + 1 = 2P + 1`,
i.e. `V = 2P`, contradicting `2V = P+1` (`P = 1/3`). For `V = 1` the v-block
vanishes identically, forcing the p-block to vanish, i.e. `P = 0`,
a contradiction with `P ≥ 1`. ∎

**Theorem 5.1.** Branch B0 consists exactly of the constant-h, `e = 0`
members of the tame family [UP,(5.1)] and their `v = 0` degenerations. All
are tame automorphism pairs.

## 6. Nonconstant h: a one-line exclusion

[UP,(3.4)] proves `h | t+e` at arbitrary degree, so a nonconstant `h` is
affine: `h = α(t+e)`, `α ≠ 0`. Then in [UP,(3.3)] the term `(t+e)/(κh)` is
the **constant** `1/(κα)`, so, using `v(0) = w(0) = 0`,

```
a₋₁(0) = −1/(κα) ≠ 0,
```

violating the membership `τ | a₋₁`. Hence **no genuine polynomial pair has
nonconstant `h`** — with no analysis of the residual system required. (This
strengthens [UP] §5, which proved membership failure only for the members of
the canonical family (5.1).)

## 7. Assembly and scope

- `κ = 0`: no genuine pairs ([UP] §7, proved).
- `κ ≠ 0`, `h` nonconstant: no genuine pairs (§6).
- `κ ≠ 0`, `h` constant: `e = 0` forced; branches A\* (§3) and A0 (§4) empty;
  B0 exactly tame (§5).

This proves the Theorem stated at the head of the memo. **Not covered here:**
the sectors with `a₂ = 0` (top level ≤ 1 on the `F` side), which the sector
tree routes through lower-band analysis, and the quantum (Weyl) mirror. The
full band-2 classical claim "every band-2 Keller pair is an automorphism
pair" therefore remains open exactly in the `a₂ = 0` sectors.

## 8. Verification

`verify_classical_Astar.py` checks, exactly in sympy: the parametrized
vanishing of `C₄..C₁` and `C₀ = 1`; both integral identities of Lemma 2.1;
the two leading-coefficient formulas; the `S1` ODE exclusion for `W ≤ 8`;
and per-branch bounded-degree regression sweeps (A\*, A0 empty; B0 tame only).
The degree-balance arguments of §§3–5 are proofs, not computations; the
sweeps are regression evidence only.
