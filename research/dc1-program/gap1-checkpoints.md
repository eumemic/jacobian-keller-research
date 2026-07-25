# Walls as checkpoints: band-reduction Gap 1 and a conditional primitive-degree floor for the W2 leading-form stratum

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / BOUNDED MACHINE REGRESSIONS — NOT PEER
REVIEWED — EXPLORATORY, WITH A CONDITIONAL PAPER THEOREM AND NAMED GAPS**

[`band-reduction.md`](band-reduction.md) §9 **Gap 1** — the composite-move escape —
is the campaign's classical core: reducedness is proved only against *single*
transvections, so a tame *word* could temporarily raise the invariant `(n+m,k)`
and come back down outside the classified stratum. This memo formulates the
"walls as checkpoints" attack, checks bounded instances of the generator and
wall-data arithmetic, explains why the proposed local-minimum implication is not
established, and gives a conditional paper proof of a primitive-degree floor.
For any genuine pair attaining that floor, the cited inputs exclude tame descent
below it; descent from strata above the floor remains open.

Exact bounded regression support:
[`verify_gap1_checkpoints.py`](verify_gap1_checkpoints.py). It reports check and
search counts dynamically; runtime is machine-dependent. `HEAVY=1` adds the
optional depth-4 search leg. Base commit `201c2f6`. Every load-bearing
upstream fact (the ladder engine, `Q_0=(T-1)G`, the `W(k,q)` wall and its cofactor
parametrisation, the Fourier ladder rule, the leading-monomial rule) is
**re-derived in file**; the Fourier rule is additionally re-derived against a
normal-ordered `A_1` implemented from scratch.

Conventions frozen from the corpus: `A_1[x^{-1}]=⊕_k x^k C[E]`, `E=x∂`,
`(x^a f)(x^b g)=x^{a+b}f(E+b)g(E)`, `f^{[n]}(E)=f(E+n)`,
`Q_m=∑_{k+l=m}[b_l^{[k]}a_k-a_k^{[l]}b_l]`, `[D,X]=1 ⟺ Q_m=δ_{m0}`, membership
`(E)_j|a_{-j},b_{-j}`, gauge `b_k=0`. Invariant `(n+m,k)`: `n=deg X`, `m=deg D`
Bernstein, `k=max(width X, width D)`, lexicographic.

---

## 0. Headline

> **Conditional paper statement.** Over `C` (more generally an algebraically closed
> characteristic-zero field), assuming Dixmier's common-primitive leading-symbol
> lemma, every genuine Weyl pair with primitive degree `e≥2` has `a,b≥2`,
> `(a,b)≠(2,2)`, and therefore `n+m≥5e`, with equality exactly when
> `{a,b}={2,3}`. Conditional on the cited tame and fixed-band inputs, a genuine
> floor-attaining pair is tame-orbit minimal in `(n+m,k)`. This does not say that
> every orbit-minimal or minimal counterexample attains equality: descent from
> `a+b>5` to the floor is open. The verifier is bounded regression support, not a
> proof object.**
>
> **TIER (audit, 2026-07-25) — read this before citing anything below.** The
> headline is a **PAPER PROOF, CONDITIONAL** on four inputs not re-derived in
> file: Dixmier's 1968 leading-symbol lemma, Makar-Limanov `Aut(A_1)`=tame, and
> the campaign's band-1 rigidity `P3` + band-2 theorem (`84978b9`). Every
> "(PROVED, degree-free)" tag in items 1–4 below means *paper-proved with
> fixed-degree or single-instance machine support* — no item carries a
> symbolic/degree-free proof object in the verifier, and the former `§5.2` floor
> check was a tautology (now removed). See the §8 ledger for the exact machine
> scope of each line.
>
> 1. **Generator arithmetic (PROVED, degree-free).** `deg(X^s)=s·deg X`,
>    `band(X^s)=s·band(X)` at both extremes (the extreme coefficient of `X^s` is a
>    nonzero product `∏_j a_k^{[jk]}` in the domain `C[E]`). Hence the exact
>    trichotomy for `D↦D-f(X)`, `deg f=s`: `sn>m` ⟹ strict raise to `m'=sn`;
>    `sn<m` ⟹ `m'=m`; `sn=m` ⟹ strict drop **iff** `σ_D=cσ_X^s`. Fourier and
>    pair-exchange preserve `(n+m,k)` exactly; affine `X↦X+bD` sends `n↦max(n,m)`.
>    **Corollary:** a single generator lowers `n+m` only if `a|b` or `b|a`.
> 2. **The proposed local-minimum proof/implication is not established.** A local
>    minimum of the trajectory certifies only that the two *adjacent* letters do
>    not lower the invariant; that observation does not prove that *no* generator
>    lowers it. No committed genuine-pair witness refuting the statement is supplied.
>    What survives is the
>    **last-ascent checkpoint lemma** (§3): the escape must leave the
>    mutually-non-dividing stratum and re-enter below `μ`, and its final move is a
>    *deep* cancellation whose every unit of depth is one more `p`-power condition.
> 3. **Wall-data transformation table (PROVED).** A raising transvection sends
>    `W(k,q)↦W(sk,k)` with cofactor `g↦cyc(k,d)·g = δ(a_k)`, necklace law
>    `δ(b'_{sk})=(1+σ^k+⋯+σ^{(s-1)k})δ(a_k)`; the image wall has `q'|k'`, so **any
>    raising transvection destroys the singular-hatch branch** and lands on the
>    forced-effective (shifted-power) branch. Pair-exchange preserves `g` and swaps
>    `(k,q)`; **Fourier is not a symmetry of the top wall at all** — it maps the
>    top wall to the *bottom* wall (machine-checked refutation: two elements with
>    identical top data have different Fourier tops).
> 4. **Monovariant hunt.** `gcd(k,q)`, cofactor effectivity, and the non-dividing
>    property of `(a,b)` fail in the tested transformations. The over-approximate
>    degree-transition relation admits a short escape from `(9,6)`, so that tested
>    relation is insufficient by itself. This does not refute every conceivable
>    monovariant depending on actual `(n,m,k)` trajectories. The useful invariant is
>    the conditional **primitive degree**
>    `e=deg p`, invariant along every trajectory (§5).
> 5. **FLOOR THEOREM (conditional paper proof, arbitrary degree).** Over the stated
>    characteristic-zero field and assuming Dixmier's lemma, a genuine pair with
>    `e≥2` has `a≥2`, `b≥2`, `(a,b)≠(2,2)`, hence `n+m=(a+b)e≥5e`, with equality
>    iff `{a,b}={2,3}`. The W2 and shifted-cube displayed data are only formal
>    wall/top representatives. The conclusion applies conditionally to any genuine
>    pair having those leading forms; the representatives themselves are not called
>    tame-minimal.
> 6. **Bounded search (consistency).** A bounded, history-aware depth-3 search
>    (`HEAVY`: depth 4) over the disclosed alphabet finds no escape at the formal
>    W2 and shifted-cube data, while finding known reductions on tame controls.
>    Counts are reported dynamically. This is not a proof or realizability test.

---

## 1. The trajectory picture, made precise

The tame group `𝒯` acts on Weyl pairs (`[D,X]=1`) by: affine symplectic, Fourier
`φ:x↦-∂,∂↦x`, pair-exchange `S:(X,D)↦(D,-X)`, transvections `exp(ad p(X)):D↦D-p'(X)`
and mirrors `exp(ad q(D)):X↦X+q'(D)`. All preserve `[D,X]=1` and `⟨X,D⟩`
(re-verified, verifier `§2`).

For a word `w=g_N⋯g_1` write `V_j=g_j⋯g_1V_0` and `I_j=I(V_j)=(n_j+m_j,k_j)`. The
**trajectory** is `I_0,…,I_N`. `V` is **single-move-reduced** if no single
generator lowers `I(V)`; it is **`𝒯`-minimal** if no word does. Gap 1 is exactly

```
   Gap 1  ⟺  every single-move-reduced pair is 𝒯-minimal.
```

Two floors are proved in file before anything else:

- **Degree floor (verifier `§1`).** A Bernstein-degree-1 member is affine
  (`αx+β∂+γ`); then `Y=∂/α` satisfies `[Y,X]=1`, `D-Y` centralises `X`, the
  centraliser of an affine `X` is `C[X]` (machine-checked in a degree-≤4 window and
  transported), so `D=Y+c(X)` and `⟨X,D⟩∋D-c(X)=Y`: **the pair generates**. Since
  generation is a `𝒯`-orbit invariant, a non-generating trajectory has
  `min(n_j,m_j) ≥ 2` at *every* step.
- **Band floor (imported).** Non-generating ⟹ band `≥3` (band-1 rigidity `P3`,
  band-2 theorem `84978b9`).

---

## 2. Exact generator arithmetic (paper proof with bounded regression support)

Verifier `§2` supplies the bounded scope recorded in §8. The arbitrary-degree paper
argument uses the engine identity

```
   extreme-band coeff of X^s  =  ∏_{j=0}^{s-1} a_k^{[jk]}  ≠ 0   in the domain C[E],
```

so `deg(X^s)=s·deg X`, `σ(X^s)=σ(X)^s`, `bandtop(X^s)=s·bandtop(X)`,
`bandbot(X^s)=s·bandbot(X)`, `width(X^s)=s·width(X)` — **no cancellation is ever
possible in a power**. Consequently:

| generator | `(n,m)` | `k` |
|---|---|---|
| scaling `(λX, D/λ)`, translation | unchanged | unchanged |
| pair-exchange `S` | `(m,n)` | unchanged |
| Fourier `φ` | unchanged | unchanged (band support **reversed**) |
| affine `X↦X+bD` | `n↦max(n,m)` (drop only if `n=m` and forms cancel) | `≤max` |
| transvection `D↦D-f(X)`, `deg f=s`, `sn>m` | `m↦sn` (**strict raise**) | `≤max(k_D,s·k_X)` |
| transvection, `sn<m` | unchanged | may still move |
| transvection, `sn=m` | `m'<m` **iff** `σ_D=cσ_X^s` | — |
| mirror `X↦X+f(D)` | same law with `(n,m)`,`(a,b)` swapped | — |

**Corollary (single-move lowering criterion).** A single tame generator strictly
lowers `n+m` only if `n|m` or `m|n`, i.e. `a|b` or `b|a`. This upgrades
`band-reduction.md` §3 (stated there for transvections) to the **full generator
set** — Fourier, exchange, scaling and translation cannot lower `n+m` at all.

**Fourier, exactly.** Derived in file and cross-checked against a from-scratch
normal-ordered `A_1`:

```
   a_i  ↦  a'_{-i} = (-1)^i (E)_i a_i(-E-1)     (i ≥ 0),
   a_{-j}=(E)_j c  ↦  a'_{j} = c(-E-1)          (j > 0),
```

which preserves the Bernstein degree band-by-band and reverses the band support.
A bare coefficient reflection `E↦-E-1` that drops the falling factorial is not
the `A_1` Fourier map. In particular, `a_3=1` maps to `-(E)_3`, which is
membership-valid. As recorded in
[`shifted-power-residuals.md`](shifted-power-residuals.md) §3.1, genuine Fourier
instead reverses/exits the displayed top-wall/gauge chart, so `κ₂` is not directly
comparable there; composite tame escape remains open.

---

## 3. The checkpoint lemma — what is true, and what the proposal got wrong

> **NOT ESTABLISHED by the proposed argument.** *"At every local minimum of the
> trajectory the pair is single-transvection-reduced, hence lies in the classified
> stratum."* A local minimum only certifies that `g_j` did not lower it and
> `g_{j+1}` does not lower it — two letters, not the whole generating set. Thus the
> proposed inference to single-move-reduced does not follow. No explicit committed
> genuine-pair counterexample is given here. (The
> *trajectory-global* minimum of a length-minimal word achieving the orbit minimum
> *is* single-move-reduced, but that is the orbit minimum by definition and carries
> no new information.)

What survives is the following prose lemma; verifier `§3` prints the result and
checks only the same-side collapse arithmetic listed in §8:

> **Checkpoint lemma (last ascent).** Let `V_0` be single-move-reduced with
> `I(V_0)=μ`, and let `w=g_N⋯g_1` satisfy `I(wV_0)<μ`. Put
> `j=max{i≤N : I_i ≥ μ}`. Then
> 1. `j<N` and `g_{j+1}` strictly lowers `I(V_j)` below `μ`;
> 2. hence (§2 corollary) `V_j` has **dividing** Dixmier exponents — the escape
>    must *leave* the mutually-non-dividing stratum and *re-enter* below `μ`;
> 3. the member of `V_j` untouched by `g_{j+1}` has Bernstein degree `<μ₁`, while
>    `n_j+m_j ≥ μ₁`, so the cancellation depth is `≥ (s+1)n_j-μ₁+1`;
> 4. **cancellation-depth lemma:** `D-f(X)` can only remove leading forms at
>    Bernstein degrees divisible by `n=deg X` whose form is the corresponding power
>    of `σ_X=αp^a`. Every further unit of depth is one more `p`-power condition on
>    the sub-leading symbols of `D`;
> 5. reduced words **alternate sides** (consecutive same-side transvections
>    collapse: `(D-f(X))-g(X)=D-(f+g)(X)`), and `V_{j+1}=(X_j,D_j-f(X_j))` is again
>    a genuine Weyl pair (`f(X)` commutes with `X`).

Also checked in the bounded degree model: **non-cancelling words never lower
`n+m`** over the disclosed `7` starting states, length `≤5`, and `s≤4`; the
transition count is reported dynamically. Every escape in that bounded model
contains a cancelling step.

---

## 4. Wall-data transformation table (machine-checked)

Wall `W(k,q)`: `δ(a_k)=cyc(k,d)g`, `δ(b_q)=cyc(q,d)g`, `d=gcd(k,q)`,
`cyc(k,d)=(σ^k-1)/(σ^d-1)`; Dixmier `σ_X=αp^a`, `σ_D=βp^b`. Verifier `§4`, checked
on the W2 hatch (`g=1-σ+σ²`), the tame cube (`g=1`), and the census hatches
`k=4,5`.

| generator | `(k,q)` | `(a,b)` | cofactor `g` | `p` |
|---|---|---|---|---|
| scaling, translation | `(k,q)` | `(a,b)` | `g` | `p` |
| pair-exchange `S` | `(q,k)` | `(b,a)` | `g` | `p` |
| Fourier `φ` | **top wall ↦ bottom wall** — top data is *not* `φ`-covariant | `(a,b)` | — | `p∘φ` |
| `D↦D-cX^s` (`sk>q`) | `(sk,k)`, `d'=k` | `(a,sa)` | `cyc(k,d)g = δ(a_k)` | `p` |
| `X↦X+cD^s` (`sq>k`) | `(q,sq)`, `d'=q` | `(sb,b)` | `cyc(q,d)g = δ(b_q)` | `p` |
| cancelling move (`sn=m`) | `(k,q'')`, `q''` **unconstrained** by the old data | new | new | `p` |

with the exact necklace law for the raising moves

```
   δ(b'_{sk}) = (1 + σ^k + σ^{2k} + ⋯ + σ^{(s-1)k}) · δ(a_k) ,
```

verified identically in a symbolic `a_k` for `(k,s)∈{(3,2),(3,3),(4,2)}`, and
verified to satisfy `W(k,sk)` identically.

**Two structural readings.** (i) Raising moves act on the cofactor by
*multiplication* by `cyc(k,d)`, and always land on a wall with `q'|k'` — i.e. on
the **forced-effective / shifted-power** branch (`band-reduction.md` §6: exotic
requires `q∤k`). So **the singular-hatch branch is destroyed by any raising
transvection and can only be re-created by a cancelling move.** (ii) Fourier
exchanges top and bottom walls, so **no function of the top wall alone can be a
`𝒯`-monovariant** — any candidate must be symmetrised over top and bottom, or the
word class must exclude Fourier.

---

## 5. The monovariant hunt

### 5.1 Four refutations, including a decisive one

| candidate | verdict |
|---|---|
| `d=gcd(k,q)` | **non-decreasing under raising** (`d'=k≥d`), but **refuted globally**: explicit word with `d: 1→3→1` |
| cofactor effectivity | destroyed by raising, restored by the inverse — **refuted** |
| `(a,b)` mutually non-dividing | one raising move gives `(3,2)↦(3,6)` — **refuted** |
| the tested over-approximate `(n,m,k)` transition relation alone | **insufficient by itself**: it permits a short escape that the true algebra blocks. This does not rule out every monovariant depending on actual `(n,m,k)` trajectories |

### 5.2 The winner: the primitive degree `e=deg p`

Dixmier's 1968 leading-symbol lemma (classical input, cited, instance-verified in
file) gives, for `n+m>2`, `σ_X=αp^a`, `σ_D=βp^b` with `p` **primitive**
homogeneous; set `e=deg p`, so `n=ae`, `m=be` and `e | gcd(n,m)`.

> **Invariance (conditional paper proof; bounded verifier support in `§5.1`).** Over
> `C` or an algebraically closed characteristic-zero field, assuming Dixmier's
> lemma, `e` is constant along tame trajectories whose states satisfy `n+m>2`.
>
> *Reason.* A transvection changes one member only, so the other member retains its
> primitive, and Dixmier applied to the new pair controls the changed member.
> Ambient linear symplectic changes of `x,ξ` are a separate action: they carry `p`
> by an invertible linear substitution and preserve `deg p`. For an invertible
> linear recombination of the pair `(X,D)`, if the degrees differ, the higher
> leading form is retained in at least one component; cancellation in the other is
> controlled by Dixmier using that retained member. If the degrees agree, the two
> common leading forms are proportional powers with equal exponent; an invertible
> matrix cannot cancel both, and a retained component fixes the same primitive.
> Pair exchange is the special recombination swapping `(a,b)`. The verifier checks
> only bounded instances of transvections, exchange, and Fourier.

### 5.3 FLOOR THEOREM

> **Theorem (conditional paper proof, arbitrary degree).** Over `C` (or an
> algebraically closed characteristic-zero field), assume Dixmier's
> common-primitive leading-symbol lemma. Let `[D,X]=1` have primitive degree
> `e≥2`. Then `a≥2`, `b≥2`, and `(a,b)≠(2,2)`; hence
> ```
>        n + m = (a+b)·e  ≥  5e ,     with equality iff {a,b}={2,3}.
> ```
>
> *Proof.* If `a=1` then `σ_D=βp^b` is proportional to `σ_X^b=α^bp^b`, so the
> transvection `D↦D₁=D-(β/α^b)X^b` cancels the top: `deg D₁<be`. If `deg D₁=0`
> then `D₁∈C` and `[D₁,X]=0≠1` — contradiction (`[D₁,X]=[D,X]=1`). Otherwise
> `n+deg D₁ ≥ e+1 ≥ 3 > 2`, so the lemma applies to `(X,D₁)`; since `σ_X=αp` with
> `p` primitive, the common primitive is `p` and `deg D₁=b₁e` with `1≤b₁<b`.
> Iterate: `b` strictly decreases, so we reach `b=1`, where the same move gives
> `deg D'<e` and `e|deg D'` forces `deg D'=0` — contradiction. Hence `a≥2`, and by
> the mirror argument `b≥2`. If `a=b=2`, the affine move `D↦D-(β/α)X` gives
> `deg D₁<2e` with `e|deg D₁`, so `deg D₁∈{0,e}`; `0` is excluded as above and `e`
> means `b₁=1`, excluded by the previous step. ∎

> **Corollary (T-minimality at the floor).** A Weyl pair with `e≥2` and
> `{a,b}={2,3}` is **`𝒯`-minimal**. Indeed `e` is orbit-invariant (§5.2) — and the
> orbit cannot leak out through `n+m=2`, since a transvection leaves one member's
> degree in `eℤ_{≥2}` while `n+m=2` forces `n=m=1`, and the linear letters cannot
> lower `n+m` at all — so every state of the orbit has `n+m ≥ 5e = n+m`. Further,
> `e≥2` pairs are never tame images of `(x,∂)`: a degree-raising move out of an
> affine state produces `σ = ` a power of the *linear* symbol of the affine member,
> i.e. `e=1` (machine-checked at `(x,∂)` for `s=2,3,4` on both sides; degree-free in
> general because a Bernstein-degree-1 member has a linear symbol), and `e` is
> invariant thereafter. So by the band-2 theorem the band is `≥3`, which is the
> minimum. Both components of the lexicographic invariant are therefore at their
> orbit floor.

**Conditional consequence for the W2 leading-form stratum.** The W2 hatch
(`σ_X=(x²ξ)³`, `σ_D=(x²ξ)²`) and shifted-cube top data share `e=3`,
`(a,b)=(3,2)`, `n+m=15=5·3`, `k=3`. The displayed representatives are formal,
not Weyl pairs, so no orbit-minimality claim is made about them. Rather, the
corollary says that **if a genuine Weyl pair has this leading-form stratum, no
composite tame word lowers its invariant** (conditional on the stated inputs).

### 5.4 Two consequences the campaign should record

- **The band-`k` census hatches, `k≥4`, are NOT at the floor**: `e=3`,
  `(a,b)=(k,k-1)`, `n+m=(2k-1)·3 > 15`. A composite descent to a band-3
  configuration is *not* excluded for them (they die at the moment unit anyway,
  [`hatch-census.md`](hatch-census.md)).
- **The lowest arithmetically admissible `e≥2` floor case is `e=2`**:
  `n+m=10<15`, with `p~xξ` up to `SL₂`, `σ_X=(xξ)²`, `σ_D=(xξ)³` after choosing
  suitable coordinates. In those `p=xξ` coordinates, both `p²` and `p³` are band
  `0`. A generic ambient `SL₂` change can spread either top across several bands
  in the original fixed ladder coordinates. Band is not `SL₂`-invariant, so
  existing fixed single-band classifications do not automatically transfer; it is
  also wrong to call the stratum intrinsically multi-band. Realizability remains
  open, and candidate `e=1` sectors remain open independently.

---

## 6. Bounded short-word search (verifier `§6`)

Alphabet per node: `oD`,`oX` — the Bernstein-degree-**optimal** transvection over
**all** `f` of degree `≤S` (exact: the `X^j` have pairwise distinct degrees `jn`,
so greedy leading-form elimination attains the minimum over the whole linear span);
`rD1..rDS`,`rX1..rXS` — unit-coefficient raising transvections; pair-exchange;
Fourier (when membership holds). Consecutive same-side transvections are excluded
(they collapse). Cap: invariant sum `≤45`.

| seed | `[D,X]=1`? | start `I` | best reached | words |
|---|---|---|---|---|
| `(x³-∂, x)` (positive control) | yes | `(4,3)` | **`(2,1)` via `oX`** | reported dynamically |
| tame family `c₁≠0` (positive control) | yes | `(6,3)` | **`(2,1)` via `oD·oX`** | reported dynamically |
| **W2 formal wall datum** | **no** — positive cascade `Q_{m>0}=0` and `Q₀=1` hold, the negative tail does not | `(15,3)` | no drop found within the bound | reported dynamically |
| **shifted-cube formal top datum** `a₃=E(E+1)(E+2)`, `b₂=E(E+1)` | not a Weyl pair | `(15,3)` | no drop found within the bound | reported dynamically |

**Not covered** (exact scope): non-unit raising coefficients; transvections of
degree `>S` (`S=2` for the band-3 data, `S=3` for the controls); affine-symplectic
letters with both off-diagonal entries nonzero; words longer than the depth. The
search is a **consistency probe** of §5.3, not a proof — and note that §5.3 *does*
prove what the search only samples, for the first component and for the band.

The W2 datum used throughout is the `w2-decisive.md` §4 data: machine-checked here
to satisfy the positive cascade and `Q₀=1` but **not** `Q_{-1..-5}=0`. It is a
**wall representative**, not a Weyl pair; the §5 theorem is what carries the
weight, since it is a statement about leading forms only.

---

## 7. The precise conditional statement

> **Conditional floor statement.** Over the stated characteristic-zero field,
> assuming Dixmier's lemma and the cited band inputs, a genuine pair `V` with
> `e(V)≥2` and `{a,b}={2,3}` is tame-orbit minimal in the claimed invariant.
>
> **Gap 1 for the remaining strata** reduces to: *for a single-move-reduced pair
> `V` with `e(V)≥2` and `a+b>5`, no tame word reaches a state with the same `e` and
> a strictly smaller `(a+b, k)`.* By §3 such a word must contain a cancelling step
> at a state with dividing exponents, and by §3.4 that step's cancellation depth is
> a stack of `p`-power conditions on the sub-leading symbols. A **proposed
> per-parameter condition** is to test those cancellation equations for each fixed
> `(a,b,s,e)`. This is not a completed or global finite reduction: polynomiality
> and uniform degree control remain open.
>
> **For `e=1` the floor argument is vacuous** (`5e=5` is below everything of
> interest and the exponents are the degrees themselves): the constant-`h`/`κ₂`
> sector of [`shifted-power-residuals.md`](shifted-power-residuals.md) §3, whose
> leading form is `(3,2)` with `p=x` linear, is **untouched by this memo**.

**Relation to the classical descent.** For *automorphism* pairs the greedy descent
never stalls: divisibility `a|b` or `b|a` holds at every stage (Jung–van der Kulk;
Makar-Limanov's theorem `Aut(A_1)` = tame), so no composite escape question arises.
The present memo isolates *why* the endomorphism case differs and *where* it does
not: the classical proof uses invertibility to force divisibility, and what
replaces it here is the primitive-degree invariant `e`. In these terms:

```
   generates  ⟹  automorphism  ⟹ (Makar-Limanov) tame  ⟹  e = 1.
```

Thus `e≥2` describes a candidate-counterexample sector at symbol level, while
`e=1` candidate sectors remain open. The floor theorem says only that equality,
when attained by a genuine pair, is exactly `{a,b}={2,3}`. It does **not** imply
that a minimal counterexample or every orbit-minimal pair attains equality; descent
from `a+b>5` to the floor is the open Gap 1 problem.

---

## 8. Honest ledger

> **AUDIT RE-TIERING (2026-07-25).** An adversarial audit found that **no** item
> below carried a symbolic/degree-free *proof object in the verifier*: each is a
> paper derivation supported by fixed-degree or single-instance machine checks,
> and the `§5.2` floor check was a **tautology** (its `floor_chain` helper never
> read its `e` argument, so the check passed identically at `e=1`, where the
> file's own tame control refutes the encoded conclusion; the check has been
> removed — see `§9`). The mathematics survives; the *tier* does not. Everything
> formerly listed here as "Proved (arbitrary degree)" is re-tiered below as
> **paper proof + bounded machine scope**, which is what it is.

**Paper proofs (arbitrary degree in the argument), with their in-file machine scope:**
- The power law `deg(X^s)=s·deg X`, `σ(X^s)=σ(X)^s`, extreme-band multiplicativity;
  the generator arithmetic table and the single-move lowering criterion (§2).
  *Machine scope: `s∈{2,3}`, three band shapes, coefficient degree `≤2`.*
- The **transvection trichotomy** and its corollary. *Machine scope: one pair
  `(x+3, ∂+x²)`, `s=1..5`; the `s·n<m` row only at `s=1`.*
- The Fourier ladder rule, cross-checked against normal-ordered `A_1`.
  *Machine scope: `i≤3`, `deg f≤3`, `j≤3`. General by linearity (prose).*
- The degree-1 generation floor ⟹ `min(n_j,m_j)≥2` along a non-generating
  trajectory. *Machine scope: centraliser of `x` in a degree-`≤4` window; the
  affine transport and the conclusion are prose.*
- The checkpoint (last-ascent) lemma (§3), same-side collapse, the
  cancellation-depth lemma. *These are **prose lemmas**: statements (i)–(iv) are
  printed, not machine-checked. Real machine content: the same-side collapse
  arithmetic only.*
- The wall-data transformation table (§4), the raising necklace law, and the
  destruction of the singular-hatch branch by raising moves. *Machine scope:
  `(k,q)∈{(3,2),(4,3),(5,4)}`, `s∈{2,3}`; the symbolic-`a_k` identity only at
  `(k,s)∈{(3,2),(3,3),(4,2)}`, `deg a_k≤3`. The universal "any raising
  transvection" is **not** supported by a degree-free object.*
- **`e`-invariance and the FLOOR THEOREM `n+m ≥ 5e` for `e≥2`**, and the resulting
  `𝒯`-minimality of the `{a,b}={2,3}` stratum (§5) — **a paper proof, and
  CONDITIONAL** on four inputs **not** re-derived here: Dixmier's 1968
  leading-symbol lemma (instance-verified only as `{σ_X,σ_D}=0` on one datum),
  Makar-Limanov `Aut(A_1)`=tame, and the campaign's band-1 rigidity `P3` and
  band-2 theorem (`84978b9`) for the band clause. *In-file machine scope: the
  cancellation step at `p=x²ξ`, `b≤3`. `e`-invariance is checked on two data sets
  for transvections `s=1,2,3` (both sides), pair-exchange and Fourier — the
  **affine symplectic letter is not checked at all**, and it is the case with real
  content (both members move).*

**Bounded / finite evidence (exact scope):**
- The non-cancelling-word monotonicity: exhaustive in the degree model over 7
  starts, length `≤5`, `s≤4`.
- The history-aware short-word escape search: depth 3 by default, depth 4 under
  `HEAVY=1`, `S=2` at the band-3 seeds, unit raising coefficients only. Explored
  state/word counts are printed dynamically; runtime is machine-dependent.
- The centraliser-of-`x` computation is a degree-≤4 window.
- Census/wall checks at `k=3,4,5`; `(k,s)∈{(3,2),(3,3),(4,2)}` for the symbolic
  raising-wall identity.

**Disposition of proposed implications and monovariants:**
- The proposed proof of the **local-minimum implication** ("at every local minimum
  the pair is single-move-reduced") is not established: adjacent-letter information
  does not supply the universal single-move claim. No explicit committed genuine-pair
  witness is included, so this memo does not label the implication false/refuted.
- `gcd(k,q)`, cofactor effectivity, and non-dividing `(a,b)` as monovariants
  (explicit words, machine-checked).
- The tested over-approximate `(n,m,k)` transition relation is insufficient by
  itself: it permits an escape that the true algebra blocks. This does not refute
  every conceivable monovariant depending on actual `(n,m,k)` trajectories.
- Top-wall data as a Fourier covariant.

**Open / NOT claimed:**
1. Gap 1 **above** the floor (`e≥2`, `a+b>5`) — the proposed per-parameter
   conditions of §7 are stated, not discharged; polynomiality and uniform degree
   control remain open.
2. Gap 1 for **`e=1`** strata (the `κ₂≠0` constant-`h` sector) — untouched.
3. The **`e=2` case** (`n+m=10`, `p~xξ`) — the lowest arithmetically admissible
   `e≥2` floor case, not a realized stratum. In `p=xξ` coordinates `p²,p³` are band
   0; generic `SL₂` changes may spread them in fixed coordinates, and fixed-band
   results do not automatically transfer.
4. Whether any Weyl pair with `e≥2` **exists** at all: nothing here constructs one
   or rules one out. No Weyl pair, no DC1 counterexample, no DC1/JC2 progress
   beyond the stated leading-form statements.
5. The corrected Fourier calculation preserves membership but reverses/exits the
   displayed `κ₂` chart. Composite tame escape and arbitrary-degree negative-tail
   closure remain open.

---

## 9. Verification

```sh
uv run --with sympy python research/dc1-program/verify_gap1_checkpoints.py
HEAVY=1 uv run --with sympy python research/dc1-program/verify_gap1_checkpoints.py  # optional depth-4 search
```

`§0` engine + normal-ordered `A_1` bridge; `§1` invariant + degree floor; `§2`
generator arithmetic + Fourier rule; `§3` checkpoint lemma + cancellation depth;
`§4` wall-data table; `§5` monovariant hunt, refutations, floor theorem, `e=2`
ladder; `§6` bounded word search with two positive controls; `§7` assembled
statement. Default run ends `ALL EXECUTED CHECKS PASSED; OPTIONAL CHECKS SKIPPED`;
under `HEAVY=1`, `ALL GAP1 CHECKPOINT CHECKS PASSED`.
