# The shifted-power descent: band-reduction Gap 2, resolved at band 3 (cube-separated)

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED — PARTIAL RESULT WITH NAMED GAPS**

Architecture step 3, [`band-reduction.md`](band-reduction.md) §9 **Gap 2**. The
Band Normalization Theorem forces a tame-minimal non-generating pair onto a
`W(k,q)` wall; on the **shifted-power** branch (effective cofactor) the top is a
consecutive shifted `k`-th power `a_k = c·h(E)h(E+1)⋯h(E+k-1)`. `W(k,q)`
constrains only the top **two** ladder rungs. Gap 2 asks whether the
shifted-power structure **propagates** to sub-leading bands and licenses a
genuine band-lowering tame move — the difference-operator analogue of the
Abhyankar–Moh / Newton-polygon descent.

This memo settles Gap 2 **at band 3 for cube-separated `h`**, and states exactly
what is uniform in `k` and what remains open. The one-line summary:

> **The descent is not a leading-form move.** The nonconstant shifted cube sits
> at the Dixmier-stuck leading form `(x²ξ)³` `(a,b)=(3,2)`, where **no single
> transvection lowers `(n+m,k)`**. The reduction is entirely **sub-leading**: the
> positive cascade + the moment unit `Q₀=1` **force `h` constant** (arbitrary
> degree, cube-separated, **conditional** on the full sub-leading ansatz — §3),
> after which an **explicit two-move tame
> word** collapses the pair to **band 1** (instance-verified). What resolves the gap is the
> *propagation* (h-forcing), not a symbol-level cancellation.

Exact certificate: [`verify_shifted_power_descent.py`](verify_shifted_power_descent.py)
(69 checks, ends `ALL SHIFTED POWER DESCENT CHECKS PASSED`, ~1 s). Every
load-bearing upstream fact (crossed-product engine, `Q₀=(T-1)G`, the `W(k,q)`
wall, the h-constant forcing `G=h^{[-1]}M`) is **re-derived in file**, not merely
cited.

Conventions frozen from the corpus: `A₁[x^{-1}]=⊕_k x^k C[E]`, `E=x∂`,
`(x^a f)(x^b g)=x^{a+b}f(E+b)g(E)`, `f^{[n]}(E)=f(E+n)`, `T f=f^{[1]}`,
`S_n=1+T+⋯+T^{n-1}`, `Q_m=∑_{k+l=m}[b_l^{[k]}a_k-a_k^{[l]}b_l]`,
`[D,X]=1 ⇔ Q_m=δ_{m0}`, membership `(E)_j=E(E-1)⋯(E-j+1)|a_{-j},b_{-j}`, gauge
`b_k=0`, `G=∑_{k≥1}∑_{j=0}^{k-1}(a_k^{[j-k]}b_{-k}^{[j]}-b_k^{[j-k]}a_{-k}^{[j]})`,
`Q₀=(T-1)G`. Invariant: `(n+m, k)`, `n=deg X`, `m=deg D` (Bernstein degree
`max_k(k+2\,deg\,a_k)`), `k=max(band X, band D)`, ordered lexicographically.

---

## 0. What Gap 2 is, and the shape of the answer

`W(k,q)` (`band-reduction.md` §5) gives, for `band X=k`, `band D=q`, the
unconditional top equation `Q_{k+q}=b_q^{[k]}a_k-a_k^{[q]}b_q=0`. On the balanced
wall `q=k-1`, gauge `b_k=0` (`W1`), the shifted-power branch is

```
a_k = c·∏_{j=0}^{k-1} h^{[j]}   (consecutive k-th power),   b_{k-1} = κ·∏_{j=0}^{k-2} h^{[j]}.
```

The classical intuition (`band-reduction.md` §9 Gap 2) is the **positive
control** `X = Uᵏ - ∂/κ`, `D = κU`, which "should" telescope down. The gap is
that `W(k,q)` pins only two rungs; propagation to lower bands is unproved.

The resolution has two independent parts, and it is important to keep them
separate:

- **(P) Propagation / h-forcing** — an *arbitrary-degree theorem* (cube-separated
  `h`): the next rung `Q_{2k-2}` forces the sub-leading coefficients to inherit
  shifted-power divisibility, the central potential collapses to `G=h^{[-1]}M`,
  and `Q₀=1` forces **`h` constant** (band 3, §3 below; re-verified from
  [`../band3/quantum-shifted-cube.md`](../band3/quantum-shifted-cube.md)).
- **(R) Reduction / the tame move** — an *explicit tame word*: once `h` is
  constant the pair is in the classical tame family, and a two-move word
  (gauge + one mirror transvection) collapses it to **band 1** (§4). This part is
  verified on genuine `[D,X]=1` data; its degree-free completeness for *every*
  constant-`h` pair rests on the (open) constant-`h` negative-tail classification.

The gap's own phrasing — "shifted-power is *expected* to reduce by the classical
telescoping tower" — is **corrected**: the reduction is real, but it is not a
telescoping cancellation of leading forms (§4.3 refutes that); it is the h-forcing
(P) that does the work.

---

## 1. The engine and the wall (re-derived)

Verifier `§0` re-derives `Q_m=[D,X]_m` for `m∈[-6,6]` and `Q₀=(T-1)G` with the
closed-form staggered potential on generic band-3 data; and confirms the
tame-move generators (transvections `exp(ad p(X))`, `exp(ad q(D))`,
pair-exchange) preserve `[D,X]=1`. Verifier `§1` re-derives, for `k=3,4` and
`deg h=1,2`:

- `W1`: `Q_{2k}=b_k^{[k]}a_k-a_k^{[k]}b_k=0` under `b_k=λa_k` (top Wronskian);
- the wall `W(k,k-1)`: `Q_{2k-1}=b_{k-1}^{[k]}a_k-a_k^{[k-1]}b_{k-1}=0` is solved
  by the shifted power `a_k=∏h^{[j]}`, `b_{k-1}=κ∏h^{[j]}`;
- the necklace form of the wall for the consecutive power
  (`(σ^k-1)S_{k-1}=(σ^{k-1}-1)S_k`, both `=(σ-1)S_kS_{k-1}`), the `g=1` instance
  of the cofactor parametrization `δ(a_k)=S_k g`, `δ(b_{k-1})=S_{k-1}g`.

---

## 2. Task 1 — the next rung `Q_{k+q-1}=Q_{2k-2}` and h-compatibility

Balanced `q=k-1`, so `k+q-1=2k-2`. In gauge `b_k=0`, only three `(i,l)` with
`i+l=2k-2` and `i,l∈[-k,k]` contribute, and `b_k=0` kills `(k-2,k)`:

```
Q_{2k-2} = [ b_{k-2}^{[k]} a_k − a_k^{[k-2]} b_{k-2} ]   +   [ b_{k-1}^{[k-1]} a_{k-1} − a_{k-1}^{[k-1]} b_{k-1} ].
```

Verifier `§2` checks this explicit four-term form against the engine at `k=3,4`.
Substituting the shifted powers `a_k=∏_{j<k}h^{[j]}`, `b_{k-1}=κ∏_{j<k-1}h^{[j]}`
and using **shift-separation** `gcd(h,h^{[j]})=1`, the rung forces the
**h-compatibility divisibilities**

```
∏_{j=0}^{k-2} h^{[j]}  |  a_{k-1}          (the (k-1)-fold shifted product),
∏_{j=0}^{k-3} h^{[j]}  |  b_{k-2}          (the (k-2)-fold shifted product),
```

i.e. **each sub-leading coefficient is a shifted-power correction times
lower-order slack**: `a_{k-1}=∏_{j<k-1}h^{[j]}·g`, `b_{k-2}=∏_{j<k-2}h^{[j]}·β`.
This is the exact analogue, one band up, of the band-3 shifted-cube divisibilities
`h h^{[1]}|a₂`, `h|b₁`.

**k=3 (arbitrary degree).** `Q₄` collapses to
`a₃h^{[3]}[(T³-1)β - κ(T²-1)g]`, and via `T^n-1=(T-1)S_n` this is the **quantum
`S₃` midpoint** `S₃β = κS₂g + γ₁` (`band-reduction`/`quantum-shifted-cube` §2).
Verifier `§2` confirms the collapse and the midpoint symbolically at generic `g,β`.

**Necessity (machine-checked).** At the concrete shift-separated top `h=E(2E-1)`
(roots `{0,1/2}`), solving `Q₄=0` (`k=3`) / `Q₆=0` (`k=4`) for a *fully generic*
sub-leading pair `(a_{k-1},b_{k-2})` returns a locus on which the divisibilities
above hold and the sub-leading data is nonzero — the divisibility is **forced**,
not assumed (verifier `§2`, `k=3` and `k=4`).

> **Status.** The explicit `Q_{2k-2}` form and the `k=3` `S₃`-midpoint collapse
> are arbitrary-degree identities. The **h-compatibility divisibility** is the
> shift-separated coprimality argument (degree-free in prose, as in
> `quantum-shifted-cube` §2), machine-verified as *necessity* at a concrete
> separated `h` for `k=3` **and** `k=4`. General-`k` uniformity is the written
> coprimality argument, not a per-`k` computation.

---

## 3. Tasks 2/3 (Propagation) — `h` is forced constant (band 3, cube-separated)

This is the arbitrary-degree half of the descent, re-verified in file from
[`../band3/quantum-shifted-cube.md`](../band3/quantum-shifted-cube.md). Feeding
the §2 divisibilities into the closed-form potential factors it term-by-term:

```
G = h^{[-1]} · M,     M = (nine explicit terms, each carrying a negative-band coefficient).
```

Membership `(E)_j|a_{-j},b_{-j}` makes every term of `M` vanish at `E=0`, so
`M(0)=0`. Then `Q₀=(T-1)G=1=(T-1)E` forces `G=E+c`; `G(0)=0` gives `G=h^{[-1]}M=E`,
so `h^{[-1]}|E`. A degree-`≥2` `h^{[-1]}` cannot divide the linear `E`; a degree-1
`h=α(E-ρ)` needs `ρ=-1`, forcing `M=1/α`, which contradicts `M(0)=0`. **Hence `h`
is constant.** Verifier `§3` checks `G=h^{[-1]}M` and `M(0)=0` for `h` **linear**
(`E-ρ`) **and** `h` **quadratic** (cube-separated `(E-1/3)(E-9/4)`), plus every
step of the affine kill.

> **Conditionality (audit-flagged).** The factorization `G=h^{[-1]}M` *requires
> the full sub-leading shifted-power ansatz*, in particular `h | a₁` (the
> `a₁=h·p` substitution imported from
> [`../band3/quantum-shifted-cube.md`](../band3/quantum-shifted-cube.md)) — with
> generic `a₁` the divisibility of `G` by `h^{[-1]}` fails (independently
> confirmed in audit). What §2 machine-forces from the top rung is only the top
> layer (`h h^{[1]} | a₂`, `h | b₁`). So "`Q₀=1 ⇒ h` constant" is proved
> arbitrary-degree **given** the full ansatz; deriving `h | a₁` from the ladder
> in-file (rather than importing it) is the remaining step of the propagation.

The same forcing seen from the `U`-side is a two-liner: the shifted power is the
top of `U^k` with `U=x·h(E)`, and

```
[∂, x·h(E)] = (E+1)h(E) − E·h(E−1) = 1   ⟺   h = 1.
```

Verifier `§3` checks `[∂, x·h]` for `h∈{1, E-ρ, E²-ρ}`: it equals `1` **iff**
`h` is the constant `1`. So the shifted-power block `U^k` is a genuine Weyl
building block only when `U` is affine — exactly the ladder h-forcing, restated.

> **This is the resolution of Gap 2's "propagation."** The shifted-power top does
> propagate to the sub-leading bands (§2), and the propagation, combined with the
> moment unit, **eliminates every nonconstant cube-separated top**. Arbitrary
> degree, cube-separated `h`, any `κ` (including the `κ=0` fold-in) — **modulo
> the conditionality note above** (`h | a₁` imported, not derived in-file).

---

## 4. Task 2 (Reduction) — the explicit band-3 tame move

### 4.1 The move

With `h` constant the pair lands in the classical tame family
`U=x+c₀+c₁∂`, `X=U³-∂/κ-A`, `D=λX+κU+β`. The **band-lowering tame word** is

```
(i)   D  ↦  D' = D − λX              gauge  exp(ad ½λX²): kills b₃, b₂, b₋₃ at once
                                      (the B0-band3 collapse; gauged wall const κ₂=0)  ⇒  band D' = 1;
(ii)  X  ↦  X'' = X − ((D'−β)/κ)³     mirror exp(ad q(D')), q'(D')=−U³ (q=−(D'−β)⁴/4κ³)  ⇒  band X'' = 1.
```

Move (i) drops `band D` from 3 to 1 (`D'=κU+β` is band 1); move (ii) cancels the
`U³` top (`(D'-β)/κ = U`) leaving `X''=-∂/κ-A`, band 1. The result
`(X'',D')=(-∂/κ-A, κU+β)` is **band 1 — affine symplectic — and generates `A₁`**
(band-1 rigidity `P3`). Both moves are transvections, hence preserve `[D,X]=1`.

### 4.2 Verified on genuine data (verifier `§4`)

| datum | start `(n+m,k)` | after gauge | after mirror | `[D,X]` |
|---|---|---|---|---|
| positive control `x³-∂/κ, κx` | `(4,3)` | — | `(2,1)` | `=1` throughout |
| full tame family `c₁≠0` | `(6,3)` | `(4,3)` | `(2,1)` | `=1` throughout |

Each step is checked to **strictly lower `(n+m,k)`** lexicographically **and**
preserve `[D,X]=1`; the mirror move's `((D'-β)/κ)` is confirmed to recover `U`
exactly, and `X''=X-U³` to be band 1. The general uniform move (any `k`) is §5.

### 4.3 Refutation — the descent is *not* a single leading-form move

The nonconstant shifted cube has Dixmier exponents `(a,b)=(3,2)` (leading form
`(x²ξ)³`), **mutually non-dividing**. Verifier `§4` confirms that for `D` of
band 2, `D²` already has a nonzero **band-4** top, so any `q'(D)` reaching band 3
(needs `deg≥2` in `D`) simultaneously injects band 4 and **raises** `band X`.
Hence **no single transvection lowers `(n+m,k)`** from `(3,2)`: the classical
"cancel the leading forms" step **stalls**. The reduction of §4.1 is available
only *after* the sub-leading h-forcing (§3) has collapsed the top to a constant —
which is precisely why Gap 2 is a genuine (difference-operator) Abhyankar–Moh
step and not a one-line symbol computation.

---

## 5. Task 4 — general `(k,q)`: the uniform move and the imbalanced walls

**Uniform reducing move (any `k`).** For the shifted `k`-th power positive control
`X=Uᵏ-∂/κ`, `D=κU` (`U=x`, `h≡1`), the single mirror transvection

```
X ↦ X − (D/κ)^k = −∂/κ           (band k ⤳ band 1)
```

lowers `(n+m,k)` from `(k+1,k)` to `(2,1)` and preserves `[D,X]=1`. Verified
`k=3,4,5` (verifier `§5`). The `Q_{2k-2}` h-forcing rung (§2, verified `k=3,4`)
plus the `[∂,U]=1` forcing (§3) are the uniform structural inputs; the balanced
wall's descent is uniform in `k` *modulo* the two open items below.

**Imbalanced walls `q<k-1`, step-`d` `(d=gcd(k,q))` shifted powers.** These are
the same twisted-Wronskian necklace with a different `(k,q)` (`band-reduction`
§5). For the **coprime** case `d=1`:

- `(k,1)`: `band D=1`, and `1|k`, so the positive control reduces by the single
  mirror transvection above (verified `k=3,4,5`);
- `(k,q≥2)` with `q∤k`: a single transvection cannot band-match (`qj≠k`), so it
  **stalls** — this is the **composite-move escape** (`band-reduction.md` §9
  **Gap 1**, the classical DC1 core), *not* closed here. Verified `(5,2),(5,3),
  (7,2),(7,3)` fall in this class.

So the shifted-power descent closes the **balanced** wall (where the constant-`h`
gauge collapses `band D` all the way to 1); the imbalanced coprime walls with
`q∤k` reduce only **modulo Gap 1**.

---

## 6. Task 5 — consequence bookkeeping (no overclaim)

The band-3 top wall splits into three genuinely different objects that **share
the leading form `(x²ξ)³`** (symbol blindness, `band-reduction` §4;
verifier `§6`):

| branch | band-3 top `a₃` | status |
|---|---|---|
| **shifted-cube** (this memo) | `c·h h^{[1]} h^{[2]}` (consecutive roots) | **cube-separated `h`: `h` forced constant (§3), tame family descends to band 1 (§4)**; residuals below |
| **exotic AP `r≠-4`** | `E(E+r)(E+2r)`, `r≠-4` | closed by the cofactor functional `λ_r` ([`../band3/quantum-ap-lambda.md`](../band3/quantum-ap-lambda.md)) |
| **W2 singular hatch** | `E(E+2)(E+4)` (`r=2` anchored) | **OPEN** — Joint Filler Covector Lemma ([`w2-joint-theorem.md`](w2-joint-theorem.md)) |

**If the shifted-cube descent closed completely**, band 3 would read
`{exotic AP r≠-4 (closed) + W2 (open) + shifted-power (closed)}`. It does **not**
close completely. What remains for full band-3 shifted-power closure:

1. **Non-cube-separated tops** (two roots of `h` differing by `1,2,3`,
   `deg h≥2`): the §2 divisibility derivation genuinely fails and `G=h^{[-1]}M` is
   unavailable. Only **bounded evidence** exists (unit ideal at cap `D=2` for
   `E(E-1),E(E-2),E(E-3),(E-1)(E-3)`, `quantum-shifted-cube` §8). This is a **new
   sub-wall**, open at arbitrary degree.
2. **Constant-`h` with `κ₂≠0`** (gauged `b₂` at band 2, leading form `(3,2)`,
   `p=x`): the §4 reduction is proved for the **tame family** (gauged `κ₂=0`,
   `band D→1`), not for *every* constant-`h` pair. Deciding whether the
   constant-`h` sector is exactly the tame family needs the negative-tail
   `A*-band3` closure ([`../band3/quantum-band3-cascade.md`](../band3/quantum-band3-cascade.md)
   §5,§7), **open**.

And **W2 remains open** independently: the shifted-power descent does not touch
the singular hatch (its top is not a shifted cube — non-consecutive roots;
verifier `§6`). Full band-3 closure still requires the W2 joint filler lemma.

---

## 7. Honest ledger — proved / bounded / refuted / open

**Proved (exact algebra, machine-checked identities; arbitrary degree where
stated):**
- Engine `Q_m=[D,X]_m`, `Q₀=(T-1)G`; tame-move invariance of `[D,X]=1` (`§0`).
- `W(k,q)` wall and its shifted-power solution; `W1`; the consecutive-power
  necklace identity `(σ^k-1)S_{k-1}=(σ^{k-1}-1)S_k` (`§1`, `k=3,4`, `deg h=1,2`).
- **Task 1:** the explicit `Q_{2k-2}` four-term rung (`k=3,4`); the `k=3` `S₃`
  midpoint collapse `S₃β=κS₂g+γ₁` (arbitrary degree); the h-compatibility
  divisibilities `∏_{j<k-1}h^{[j]}|a_{k-1}`, `∏_{j<k-2}h^{[j]}|b_{k-2}` as
  *necessity* at a concrete shift-separated `h` (`k=3` **and** `k=4`) — the
  degree-free content being the coprimality argument (`§2`).
- **Propagation (Tasks 2/3):** the identity `G=h^{[-1]}M` with `M(0)=0` for
  cube-separated `h` (re-verified for `h` linear and quadratic) and the affine
  kill. **Tier (audit-corrected): the identity is arbitrary-degree, but the
  conclusion "`Q₀=1 ⇒ h` constant" is CONDITIONAL on the full sub-leading
  ansatz including `h | a₁`, of which only the top layer is machine-forced in
  §2 (see the §3 conditionality note).** The `U`-side forcing
  `[∂,x·h]=1 ⇔ h` constant is arbitrary-degree by a trivial degree argument;
  machine coverage is the three instances `h ∈ {1, E-ρ, E²-ρ}` (`§3`).

**Bounded / finite evidence (exact scope; includes AUDIT-DEMOTED items):**
- **Reduction move (Task 2) — audit-demoted from "proved":** the explicit tame
  word (gauge + mirror transvection) strictly lowers `(n+m,k)` to `(2,1)` and
  preserves `[D,X]=1`, machine-verified at a **single rational parameter point**
  (`c₀=2,c₁=3,A=5,λ=7,β=11`, `κ` symbolic) on the positive control and the tame
  family `c₁≠0` (`§4`). Instance-level, not arbitrary-degree.
- **General `k` (Task 4) — audit-demoted from "proved":** the uniform move
  `X↦X-(D/κ)^k` for the positive control, machine scope `k ∈ {3,4,5}` (`k` never
  symbolic); `k`-uniformity is by inspection (`(D/κ)^k = U^k` = top of `X`), not
  a machine identity (`§5`).
- h-compatibility necessity is machine-checked at the single concrete
  shift-separated top `h=E(2E-1)` for `k=3,4` (the general statement is the
  written coprimality argument, not a per-degree computation); likewise the
  `Q_{2k-2}` four-term rung and the necklace identity are machine-instanced at
  `k=3,4` (arbitrary-`k`-true by telescoping, by inspection).
- The tame-move reduction covers the tame family only; "every genuine
  constant-`h` pair is in the tame family" is **not** established here.
- Non-cube-separated emptiness is the sibling's cap-`D=2` Gröbner (imported,
  bounded).

**Refuted (machine-checked):**
- The gap's expectation that the shifted power reduces by a *leading-form*
  (single-transvection) telescoping cancellation. The `(3,2)` leading form is
  Dixmier-stuck: no single transvection lowers `(n+m,k)` (band-3 `D²` reaches
  band 4). The descent is sub-leading (the h-forcing) (`§4.3`).

**Open / NOT claimed:**
1. **Non-cube-separated shifted-cube tops** (roots differing by `1,2,3`): no
   arbitrary-degree forcing — a new sub-wall (only bounded evidence).
2. **Constant-`h` with `κ₂≠0`** (leading form `(3,2)`, `p=x`): the full
   constant-`h` classification (negative-tail `A*-band3`) is open; §4 covers the
   tame family only.
3. **Imbalanced coprime walls** `q<k-1`, `q∤k`, `q≥2`: reduce only modulo the
   composite-move escape (Gap 1 = DC1 core), untouched here.
4. **General-`k` shifted-power descent** beyond the balanced positive control:
   the rung-by-rung h-forcing is verified `k=3,4`; general `k` is the written
   argument. The full negative tail at any `k` is untouched.
5. **W2** (the singular hatch) remains open independently; band-3 closure still
   needs its joint filler lemma. No Weyl pair and no counterexample is
   constructed; DC1/JC2 untouched.

---

## 8. Verification

```sh
uv run --with sympy python research/dc1-program/verify_shifted_power_descent.py
```

Exact SymPy, 69 checks (~1 s): `§0` engine + tame moves + invariant; `§1`
`W(k,q)` wall + necklace; `§2` the `Q_{2k-2}` rung + h-compatibility (`k=3,4`);
`§3` the h-constant forcing `G=h^{[-1]}M` (`h` linear/quadratic) + affine kill +
`[∂,x·h]=1⇔h` const; `§4` the explicit tame move on genuine data + the `(3,2)`
reducedness refutation; `§5` the uniform move (`k=3,4,5`) + imbalanced walls;
`§6` band-3 branch distinctness. A successful run ends
`ALL SHIFTED POWER DESCENT CHECKS PASSED`.
