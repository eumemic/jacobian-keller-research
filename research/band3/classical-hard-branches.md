# The classical band-3 hard branches: the e ≠ 0 mixed sector and nonconstant h

**INDEPENDENTLY DERIVED — MACHINE-CHECKED IDENTITIES — NOT PEER REVIEWED — BAND-SCOPED**

This memo attacks the two branches of the classical cubic sector `a_3 = c h^3`
(gauge `b_3 = 0`) that Wave A (`research/band3/classical-band3-cascade.md`, commit
`99fe6ee`, all verifiers green) left **open** and where it proved the band-2
weapons fail:

- **(I) the `e ≠ 0` mixed sector** — gauged `b_2 = e h^2 ≠ 0`;
- **(II) nonconstant `h`** — `h` a nonconstant cube root of `a_3`.

Every displayed algebraic identity is machine-verified by `verify_classical_hard.py`
(same directory); a successful run ends `ALL CLASSICAL HARD CHECKS PASSED`. The
degree / divisibility / membership contradictions are **written proofs**; the
bounded Gröbner sweeps are **regression corroboration only**, labelled as such.

Conventions are frozen and identical to Wave A: over `C`, `τ = xξ`,
`{G,F} = G_ξ F_x − G_x F_ξ`, `F = Σ_{k=-3}^{3} x^k a_k(τ)`, `G = Σ x^k b_k(τ)`,
membership `τ^j | a_{-j}, b_{-j}`, and `C_m = Σ_{k+ℓ=m}(k a_k b_ℓ' − ℓ a_k' b_ℓ) = δ_{m0}`.
Primes are `d/dτ`. We work in the oriented, gauged cubic frame of Wave A §7:
`a_3 = c h^3` (`c ≠ 0`), `b_3 = 0`, `b_2 = e h^2`, and (Wave A 4.1)
`b_1 = κ_1 h + (2e/3c) a_2/h` with `e ≠ 0 ⇒ h | a_2`. Diagonal scaling normalizes
`c = 1` throughout the proofs (it multiplies `a_3` by `ρ^3`, `h` by `ρ`).

**Headline results (new, this memo).**

1. **The bottom is a cube (reflection).** The anti-symplectic reflection
   `R: f(x,ξ) ↦ f(ξ,x)` turns Theorem A into a statement about the trailing
   coefficient: `a_{-3} ≠ 0` forces `a_{-3} = ĉ h̄^3` with `τ | h̄`, and the entire
   **bottom** cascade is the mirror of the top — a bottom wall `b_{-2} = μ_3 a_{-2} + ê h̄^2`
   and a bottom mixed coupling with `ê ≠ 0 ⇒ h̄ | a_{-2}`. The cubic sector is
   **two-sided reflection-symmetric**, with data `(e, κ_1)` on top and `(ê, κ̂_1, μ_3)`
   on the bottom (`μ_3` the A\*-band3 datum). *(PROVED, §1.)*

2. **`e ≠ 0` is a counterexample-or-nothing locus.** Every single-shear tame
   band-3 pair has gauged `u_2 = b_2 − λ a_2 = 0`, i.e. **`e = 0`** — the tower
   proportionality of the structure lemma forces it. So the positive controls
   realize exactly `e = 0` (with `κ_1` sweeping `C^*`); `e ≠ 0` is the band-3 analogue of the
   resistant branch A\*, now at the **level-2/level-3 proportionality mismatch**
   rather than the `±3` one. *(PROVED, §2.)* Bounded Gröbner: **empty** in the
   B3-2 degree box. *(COMPUTED, §5.)*

3. **Nonconstant `h`, `e = 0`: EMPTY (complete proof).** The band-2 one-liner
   provably fails at band 3 (Wave A §7.1); the honest replacement is a two-step
   negative-cascade argument. Polynomiality of the *determined* coefficients forces
   `h^2 | a_2` then `h | a_1`; reducing the moment `M = τ` modulo `h` gives the
   divisibility **`h | τ`** Wave A said was not produced — hence `h = στ`. At that
   fixed point every moment term has `τ`-order `≥ 2` while `M = τ` has order `1`, a
   contradiction; the complementary `κ_1 = 0` case dies by a Theorem-A endgame. So
   **branch (II) with `e = 0` is empty at arbitrary degree**. *(PROVED, §3–§4;
   independent bounded Gröbner agrees, §5.)*

---

## 1. The bottom is a cube: reflection symmetry of the cubic sector

`R: f(x,ξ) ↦ f(ξ,x)` is **anti-symplectic**: `{RG, RF} = −R{G,F}` (verifier §1),
so `(RF, RG)` is again a Keller pair (bracket `−1`, rescalable to `+1`), and in
band form
```
(RF)_m = τ^{-m} a_{-m},        (RG)_m = τ^{-m} b_{-m}.
```

**Proposition 1.1 (bottom cube).** If `a_{-3} ≠ 0` then `a_{-3} = ĉ h̄^3` for some
`h̄ ∈ C[τ]`, `ĉ ∈ C^*`, with `τ | h̄`.

*Proof.* The top coefficient of `RF` is `(RF)_3 = τ^{-3} a_{-3}`, a nonzero
polynomial (membership `τ^3 | a_{-3}`). Apply **Theorem A** (Wave A §6) to the
Keller pair `(RF, RG)`: if `(RF)_3` were not a scalar cube, that pair would be
empty — but it exists (it is `R` of our pair). Hence `τ^{-3} a_{-3} = ĉ ĥ^3`, so
`a_{-3} = ĉ (τ ĥ)^3 = ĉ h̄^3` with `h̄ = τ ĥ`, and `τ | h̄`. ∎

The reflection carries the entire top cascade (Wave A §2–§4) to the bottom
(verifier §1):

- **`C_{-6}` (bottom Wronskian):** `b_{-3} = μ_3 a_{-3}`, `μ_3 ∈ C`.
- **`C_{-5}` (bottom wall):** with `b_{-2} = μ_3 a_{-2} + û_2`, `C_{-5}` reduces
  **exactly** to the reflected 2/3-power wall `−(3 a_{-3} û_2' − 2 a_{-3}' û_2)`;
  since `a_{-3} = ĉ h̄^3` is a cube, `û_2 = ê h̄^2`, so
  `b_{-2} = μ_3 a_{-2} + ê h̄^2` — the mirror of `b_2 = e h^2`.
- **`C_{-4}` (bottom mixed coupling):**
  `b_{-1} = μ_3 a_{-1} + κ̂_1 h̄ + (2ê/3ĉ) a_{-2}/h̄`, so **`ê ≠ 0 ⇒ h̄ | a_{-2}`** —
  the mirror of `e ≠ 0 ⇒ h | a_2`.

So the cubic sector is completely reflection-symmetric. The datum is
`(e, κ_1 | ê, κ̂_1 | μ_3)`; the top gauge sets the top `λ_3 = 0` but leaves `μ_3`
free, and `μ_3 ≠ 0` (gauged `b_{-3} ≠ 0`) is precisely the resistant branch
**A\*-band3**. The `moment` `M = τ` (Wave A §5) is the one two-sided first
integral coupling `(e, κ_1)` to `(ê, κ̂_1, μ_3)`.

> **Consequence (structural).** A band-3 cubic pair carries *two independent* wall
> constants `e` (top) and `ê` (bottom), and both negative extremes are cubes
> divisible by the membership power `τ^3`. Every known pair has `e = ê = 0` and
> `μ_3 = λ_3 = 0` (single-shear); the "hard" branches are exactly the deviations
> `e ≠ 0`, `μ_3 ≠ 0`, and nonconstant `h`/`h̄`.

## 2. The `e ≠ 0` mixed sector is counterexample-or-nothing

**Proposition 2.1 (tame ⇒ e = 0).** Every single-shear-origin tame band-3 pair
has gauged `e = 0`.

*Proof.* The structure lemma (`band3-tame-catalog.md`, `99fe6ee`, §2) proves the
**tower proportionality** `b_k = λ a_k` for all `|k| ≥ 2` with a *single* constant
`λ = d/b`. At `k = 2` this is `b_2 = λ a_2`, so the gauged wall datum
`u_2 = b_2 − λ_3 a_2 = b_2 − λ a_2 = 0` (the gauge is `λ_3 = λ`). Hence `e = 0`. ∎

The verifier exhibits this on the two positive controls: catalog **B3-2** (even
levels, `λ = 2`) and **B3-3** (odd-only, `λ = 1`) both satisfy `b_2 − λ a_2 = 0`,
and the single gauge `G → G − λF` annihilates the whole tower
`b_3, b_2, b_{-2}, b_{-3}` simultaneously (verifier §2). The tame family
`F = ξ + αx^3 + βx^2 + γx + δ`, `G = −x` realizes `e = 0` with `b_1 = −1 = κ_1 h`,
`h = α^{1/3}`, so `κ_1 = −α^{-1/3}` sweeps **every nonzero value** as `α` varies:
the tame positive controls realize exactly `(e, κ_1) ∈ {0} × C^*`, never `e ≠ 0`.

Therefore `e ≠ 0` is the band-3 analogue of the resistant A\* branch: it is the
proportionality **mismatch between level 2 and level 3** (`b_2 ≠ λ_3 a_2`,
equivalently the top and level-2 constants differ), exactly as A\*-band3 is the
mismatch between level `+3` and level `−3` (`μ_3 ≠ λ_3`). By the blow-up law
(`band3-tame-catalog.md` §3) a genuine such mismatch would need independent shear
sources and hence band `≥ 9`; no band-3 tame word reaches it. The sector is
**empty-or-rigid**; §5 gives its residual system and bounded-empty evidence.

## 3. The divisibility engine (both branches)

Work at `c = 1`, `a_3 = h^3`, gauge `b_3 = 0`. The forward cascade
`C_5, C_4, C_3, C_2, C_1, C_0` integrates each `b`-level in closed form; the
**level-`(−j)` operator** is `3 a_3 b_{-j}' + j a_3' b_{-j} = 3h^3 b_{-j}' + 3j h^2 h' b_{-j}`,
whose homogeneous solution is `b_{-j} ∝ h^{-j}`. Two facts drop out immediately:

**(A) Integration-constant killing.** For **nonconstant `h`** and `j ≥ 1`, the
homogeneous solution `h^{-j}` is not a polynomial, so the level-`(−j)` integration
constant must vanish. (This removes the additive freedoms `γ, δ, ε` at levels
`−1, −2, −3`; the level `0, 1, 2` constants `β, κ_1, e` survive since `h^0, h, h^2`
are polynomial.)

**(B) Closed forms and forced divisibility.** Machine-checked (verifier §3):
```
b_0    = (2e/3) a_1/h + (κ_1/3) a_2/h^2 − (e/9) a_2^2/h^4 + β                (gen. e)
b_0    = (κ_1/3) a_2/h^2 + β                                                (e = 0)
b_{-1} = (κ_1/3) a_1/h^2 − (κ_1/9) a_2^2/h^5 + γ/h                          (e = 0)
```
Polynomiality forces, for `e = 0`, `κ_1 ≠ 0`:
```
b_0  polynomial  ⇒  h^2 | a_2 ;
h^2 · b_{-1} = (κ_1/3) a_1 + h(γ − (κ_1/9)(a_2/h^2)^2)   ⇒   h | a_1 .
```
The first is immediate; the second because `h^2 | [(κ_1/3)a_1 + h(·)]` gives
`h | (κ_1/3)a_1`. For `e ≠ 0` the cascade also bites: with `a_2 = h ã_2`
(Wave A), `9h^2(b_0 − β) = 6e h a_1 + 3κ_1 h ã_2 − e ã_2^2`, so
`h | e ã_2^2 ⇒ rad(h) | ã_2` (verifier §6). The engine is the same on the bottom
via reflection (§1). **This is why nonconstant `h`/`h̄` is the obstructed
direction**: each downward integration demands another power of `h` divide the
free data, which finitely many coefficients cannot sustain — made precise next.

## 4. Nonconstant `h` in the `e = 0` sector is EMPTY (complete proof)

**Theorem 4.1 (`h | τ`).** In the `e = 0` cubic sector with `κ_1 ≠ 0`, a genuine
polynomial pair with nonconstant `h` forces `h | τ`; hence `h = στ`, `σ ≠ 0`.

*Proof.* By §3(B), `h^2 | a_2` and `h | a_1`. The moment (Wave A §5, verifier §3)
is `M = 3(a_3 b_{-3} − a_{-3}b_3) + 2(a_2 b_{-2} − a_{-2}b_2) + (a_1 b_{-1} − a_{-1}b_1) = τ`.
Here `b_3 = 0` and `b_2 = 0` (`e = 0`), so
```
M = 3 h^3 b_{-3} + 2 a_2 b_{-2} + a_1 b_{-1} − a_{-1}(κ_1 h).
```
Every term on the right is divisible by `h`: `h^3 | 3h^3 b_{-3}`,
`h^2 | a_2 ⇒ h | 2 a_2 b_{-2}`, `h | a_1 ⇒ h | a_1 b_{-1}`, and `h | κ_1 h a_{-1}`.
Thus `h | M = τ`. Since `deg h ≥ 1`, `h | τ` forces `h = στ` (degree exactly 1). ∎

This is the honest band-3 replacement for the band-2 one-liner (`h | τ + e`,
`classical-Astar.md` §6). Wave A §7.1 noted the band-2 argument fails because the
central identity is polynomial without producing a divisibility; Theorem 4.1
recovers the divisibility **from the negative cascade** (`b_{-1}` polynomiality),
exactly the "genuine negative-cascade argument" Wave A flagged as required.

**Proposition 4.2 (`b_1 = 0` sub-case).** In the `e = 0` sector with `κ_1 = 0`
(so `b_1 = b_2 = b_3 = 0`, `G` supported in `x`-levels `≤ 0`), a nonconstant-`h`
pair is empty.

*Proof (verifier §4, Theorem-A endgame).* `C_3 = 3h^3 b_0' ⇒ b_0` const.
`C_2 = 3h^2(h b_{-1})' ⇒ h b_{-1} = γ ⇒ b_{-1} = γ/h`; nonconstant `h` ⇒ `γ = 0`,
`b_{-1} = 0`. `C_1 = 3h(h^2 b_{-2})' ⇒ b_{-2} = δ/h^2 ⇒ δ = 0`, `b_{-2} = 0`.
`C_0 = 3(h^3 b_{-3})' = 1 ⇒ h^3 b_{-3} = τ/3`; membership `τ^3 | b_{-3}` makes the
left side have `τ`-order `≥ 3`, contradicting the `τ`-order-1 right side (and
`b_{-3} = 0` gives `0 = τ/3`). ∎

**Theorem 4.3 (τ-order kill of `h = στ`).** In the `e = 0` sector with `κ_1 ≠ 0`,
the case `h = στ` (equivalently `a_3 = c τ^3`) is empty.

*Proof (verifier §4).* With `h = στ`, the forced divisibilities become
`τ^2 | a_2` (from `h^2 | a_2`) and `τ | a_1` (from `h | a_1`), while membership
gives `τ | a_{-1}`, `τ | b_{-1}`, `τ^2 | b_{-2}`, `τ^3 | b_{-3}`, and
`a_3 = c τ^3`, `b_1 = κ_1 στ`. Read off the `τ`-order of each moment term:
```
3 a_3 b_{-3} : ord ≥ 3 + 3 = 6      a_1 b_{-1} : ord ≥ 1 + 1 = 2
2 a_2 b_{-2} : ord ≥ 2 + 2 = 4      a_{-1} b_1 : ord ≥ 1 + 1 = 2
```
Every term has `τ`-order `≥ 2`, so `ord_τ(M) ≥ 2`. But `M = τ` has `τ`-order `1`.
Contradiction (verifier checks `M(0) = M'(0) = 0`). ∎

**Corollary 4.4 (branch (II), `e = 0`: EMPTY).** The nonconstant-`h`, `e = 0`
cubic sector contains **no genuine polynomial pair**, at arbitrary degree.

*Proof.* `κ_1 = 0` is Prop. 4.2; `κ_1 ≠ 0` gives `h = στ` (Thm 4.1), excluded by
Thm 4.3. ∎

This **completely closes branch (II) in the `e = 0` sector** — a full theorem, not
a bounded check. The mechanism is the honest replacement Wave A called for: the
band-2 one-liner is superseded by a *two-step* negative-cascade argument
(`h | τ` from `b_{-1}` polynomiality, then a `τ`-order count on the moment), and
the fixed point `a_3 = c τ^3` — whose striking profile is `τ^{|k|} | a_k` for
**every** `k ≠ 0` (reflection-symmetric "double membership") — dies precisely
because that profile pushes every moment term to `τ`-order `≥ 2`. The independent
bounded Gröbner sweep (§5) agrees, as a redundant cross-check.

## 5. Residual systems and bounded corroboration

**The `e ≠ 0` residual system (constant `h`, `h = 1`, `c = 1`).** The forward
cascade determines (verifier §5)
```
b_2 = e,   b_1 = κ_1 + (2e/3) a_2,
b_0 = (2e a_1 + κ_1 a_2)/3 − e a_2^2/9 + β,
b_{-1} = (1/3)[2e a_0 + κ_1 a_1 − (2e/3) a_1 a_2 − (κ_1/3) a_2^2 + (4e/27) a_2^3] + γ,
```
`b_{-2}` from `C_1`, `b_{-3}` from the moment `C_0 = 1`. The residual is
`C_{-1} = C_{-2} = C_{-3} = 0` together with the reflected bottom data of §1
(`C_{-4}, C_{-5}, C_{-6}`) and memberships `τ|b_{-1}, τ^2|b_{-2}, τ^3|b_{-3}`,
`τ|a_{-1}, τ^2|a_{-2}, τ^3|a_{-3}`. Two free lower `F`-functions (`a_2, a_1`)
survive into the negative cascade — the extra work band 3 demands over band 2.

**W3 is obstructed (verifier §5).** The trailing coefficient `a_{-3}` enters `C_{-2}`
as
```
C_{-2} ⊃ d/dτ[(μ_3 a_1 − b_1) a_{-3}] + 2 a_{-3}(μ_3 a_1' − b_1'),
```
with residue `2 a_{-3}(μ_3 a_1' − b_1')`; at `μ_3 = 0` this is `(4e/3) a_{-3} a_2' ≠ 0`
for `e ≠ 0`. So no gauge-free W3 first integral exists — confirming Wave A §8's
obstruction. The available first integral is the moment `M = τ` alone, which
carries `a_{-3}, b_{-3}` linearly; the missing second integral is why unbounded
closure stays open.

**Lattice framework (W5) status.** The band-2 A\* kill (`band-k-weapons.md` §W5,
`99fe6ee`) needs **two** independent leading-degree balances `L d = β` to expose
the `gcd ∤ m` congruence (there: `2V = P+W` from `Φ` and `V+W = 2P+1` from `I_2`,
giving `3V = 3P+1`). Band 3's `e ≠ 0` sector supplies **only one** clean balance
here — the moment `M = τ` — because the second band-2 integral `I_2` had a
*trailing-coefficient-free* structure that the W3 obstruction destroys at `e ≠ 0`
(the residue `(4e/3)a_{-3}a_2'` is exactly the term that spoils exactness). With a
single balance there is no congruence to close on, so the lattice framework does
**not** fire in this sector; producing the second integral is the precise missing
step. This is a sharper localization than Wave A's "W3 obstructed": the
obstruction is now pinned to the one residue that a second lattice balance would
need to vanish.

**Bounded Gröbner corroboration (COMPUTED, regression only).** With the closed
forms substituted, the residual becomes a polynomial system in the `F`-data and
scalars; a Gröbner basis over `Q` decides emptiness of each bounded box.

- `e ≠ 0`, constant `h`: adjoining `e·w = 1` (Rabinowitsch) in the B3-2 degree box
  `deg(a_2,a_1,a_0,a_{-1},a_{-2},a_{-3}) = (1,1,1,2,2,3)` yields Gröbner basis `= (1)`:
  **no `e ≠ 0` pair in the box.** (Verifier gates this behind `HARD_FULL=1`; it is
  the one genuinely-new emptiness datum, since `e ≠ 0` has no unbounded proof yet.)
- `a_3 = c τ^3` (`h = στ`), `e = 0`: small boxes yield Gröbner basis `= (1)` — an
  **independent cross-check of Theorem 4.3**, which already closes this case at
  arbitrary degree. Runs by default (fast).

These are bounded and claim nothing beyond their boxes.

## 6. Claim disposition

**Proved (machine-checked identities + written degree/divisibility/membership
arguments), at arbitrary degree:**
- Reflection anti-symplecticity and `(RF)_m = τ^{-m} a_{-m}`; the **bottom cube**
  `a_{-3} = ĉ h̄^3` (`τ|h̄`) and the reflection-symmetric bottom cascade
  (bottom wall `b_{-2} = μ_3 a_{-2} + ê h̄^2`, bottom coupling with `ê ≠ 0 ⇒ h̄|a_{-2}`). §1.
- **tame ⇒ e = 0**: `e ≠ 0` contains no single-shear tame pair; it is the
  level-2/level-3 proportionality-mismatch (A\*-type) locus. §2.
- The divisibility engine: closed forms `b_0, b_{-1}`; `e=0,κ_1≠0 ⇒ h^2|a_2, h|a_1`;
  integration-constant killing; the `e≠0` divisibility `rad(h) | a_2/h`. §3.
- **Nonconstant-`h`, `e = 0`: EMPTY at arbitrary degree** (Cor. 4.4). Via
  **Theorem 4.1** (`h | τ`, hence `h = στ`), **Proposition 4.2** (`κ_1 = 0`
  Theorem-A endgame), and **Theorem 4.3** (the `τ`-order kill of `h = στ`: every
  moment term has `τ`-order `≥ 2` ≠ 1). §4.
- The `e≠0` residual system and the **W3 obstruction** (no gauge-free trailing
  integral for `e ≠ 0`). §5.

**Computed (bounded, regression corroboration only — redundant given the proofs):**
- `e ≠ 0`, constant `h`: no pair in the B3-2 degree box (Gröbner `= (1)`). §5.
- `a_3 = c τ^3`, `e = 0`: no pair in small boxes (Gröbner `= (1)`) — an independent
  cross-check of Theorem 4.3, which already proves this at arbitrary degree. §5.

**Open (precisely delimited):**
- **Unbounded** closure of `e ≠ 0` (constant `h`): empty-or-rigid. The reflection
  structure and moment are in hand; the missing ingredient is a *second* negative
  first integral (W3 is obstructed), the band-3 analogue of band-2's `I_2`.
  Bounded-empty in the B3-2 box (§5).
- The `e ≠ 0` **nonconstant-`h`** cross-branch (both deviations at once): §3 shows
  the divisibility engine bites here too (`h | a_2` and `rad(h) | a_2/h`), but the
  `τ`-order kill has not been carried through — the `b_1 = κ_1 h + (2e/3)a_2/h`
  term reduces to `(2e/3)ã_2 ≢ 0 (mod h)`, so `h | τ` is not immediate. Not closed.

**Not claimed:** a full band-3 theorem, the A\*-band3 (`μ_3 ≠ 0`) closure, JC2,
DC1, the quantum mirror, or any unbounded-degree statement beyond §1–§4.

## 7. Verification

```
uv run --with sympy python research/band3/verify_classical_hard.py
```
Exact SymPy. The script checks: reflection anti-symplecticity and the bottom
cube/wall/coupling reductions (§1); the tame `⇒ e=0` collapse on B3-2, B3-3 and
the tame witness (§2); the general-`h` closed forms and forced divisibilities,
including the `h|τ` divisibility identity (§3–§4); the `b_1=0` Theorem-A endgame
(§4); the `e≠0` residual closed forms, the moment first integral, and the W3
obstruction residue (§5); and the smallest Gröbner emptiness box of each branch
(§5, corroboration). A successful run prints `ALL CLASSICAL HARD CHECKS PASSED`.
These checks certify the algebra; the degree/divisibility/membership arguments of
§1–§5 are the proofs.
