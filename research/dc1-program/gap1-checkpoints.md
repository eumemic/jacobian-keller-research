# Walls as checkpoints: band-reduction Gap 1, and a primitive-degree floor that closes it at the W2 stratum

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED — EXPLORATORY, WITH ONE NEW ARBITRARY-DEGREE THEOREM AND NAMED GAPS**

[`band-reduction.md`](band-reduction.md) §9 **Gap 1** — the composite-move escape —
is the campaign's classical core: reducedness is proved only against *single*
transvections, so a tame *word* could temporarily raise the invariant `(n+m,k)`
and come back down outside the classified stratum. This memo formulates the
"walls as checkpoints" attack precisely, machine-checks the exact invariant
arithmetic of every tame generator and the wall-data transformation table,
**refutes the proposed local-minimum lemma as stated**, and then — from the
monovariant hunt — proves a new arbitrary-degree floor that **excludes the escape
at the W2 leading-form stratum**.

Exact certificate: [`verify_gap1_checkpoints.py`](verify_gap1_checkpoints.py)
(75 checks, ~45 s, ends `ALL EXECUTED CHECKS PASSED; OPTIONAL CHECKS SKIPPED`;
`HEAVY=1` adds the depth-4 search leg, ~115 s, then ends
`ALL GAP1 CHECKPOINT CHECKS PASSED`). Base commit `201c2f6`. Every load-bearing
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

> **The primitive degree `e=deg p` is a tame-orbit invariant, and every Weyl pair
> with `e≥2` has `n+m ≥ 5e`. The W2 configuration sits exactly on that floor, so
> no composite tame word escapes from it — Gap 1 is closed at the W2
> leading-form stratum, and open above it.**
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
> 2. **The proposed local-minimum lemma is FALSE as stated (REFUTED).** A local
>    minimum of the trajectory certifies only that the two *adjacent* letters do
>    not lower the invariant — not that *no* generator does. What survives is the
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
> 4. **Monovariant hunt (one winner, four refutations).** `gcd(k,q)`, cofactor
>    effectivity, and the non-dividing property of `(a,b)` are all **refuted** as
>    monovariants by explicit two-letter words. Worse: the **degree-only model
>    `(n,m,k)` admits a length-2 escape** from `(9,6)` to sum `9`, so **no
>    monovariant that depends only on the invariant can ever prove Gap 1**;
>    leading-form data is unavoidable. The winner is the **primitive degree**
>    `e=deg p`, invariant along every trajectory (§5).
> 5. **FLOOR THEOREM (new, arbitrary degree).** For `[D,X]=1` with `e≥2`:
>    `a≥2`, `b≥2`, `(a,b)≠(2,2)`, hence `n+m=(a+b)e ≥ 5e`; and `e≥2` pairs never
>    reach `n+m=2` or band `≤2`. The **W2 hatch and the shifted cube both have
>    `e=3`, `(a,b)=(3,2)`, `n+m=15=5e`** — the floor — and band `3`, the minimum
>    for a non-generating pair. **They are T-minimal: Gap 1's escape is
>    excluded there.**
> 6. **Bounded search (consistency).** An exhaustive depth-3 (HEAVY: depth-4)
>    word search over `{optimal transvection, unit raising transvections `s≤2`,
>    exchange, Fourier}` at the W2 datum and the shifted cube finds **no escape**
>    (538 words each at depth 4), while the same search *does* find the known
>    reductions on both tame positive controls. **No escape witness exists in the
>    searched space.**

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

## 2. Exact generator arithmetic (machine-checked, degree-free)

Verifier `§2`. The degree-free engine identity is

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
(Note: the naive "`E↦-E-1`, drop the falling factorial" reading — which
[`shifted-power-residuals.md`](shifted-power-residuals.md) §3.1 appeals to when it
says a constant top `a_3=1` goes to an `x^{-3}`-coefficient `1` — is **not** the
`A_1` Fourier map; with the falling factorial the image is `-(E)_3`, which *is*
membership-valid. The `κ₂` gauge obstruction there should be re-derived with the
rule above before being relied on.)

---

## 3. The checkpoint lemma — what is true, and what the proposal got wrong

> **REFUTED (as proposed).** *"At every local minimum of the trajectory the pair is
> single-transvection-reduced, hence lies in the classified stratum."* A local
> minimum only certifies that `g_j` did not lower it and `g_{j+1}` does not lower
> it — two letters, not the whole generating set. The implication
> "local minimum ⟹ single-move-reduced" is a quantifier error. (The
> *trajectory-global* minimum of a length-minimal word achieving the orbit minimum
> *is* single-move-reduced, but that is the orbit minimum by definition and carries
> no new information.)

What survives, and is proved (verifier `§3`):

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

Also proved in file: **non-cancelling words never lower `n+m`** (exhaustive over
`7` starting states, length `≤5`, `s≤4`, `13448` transitions) — every escape
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
| **anything depending only on `(n,m,k)`** | **refuted decisively**: the degree-only model, with cancelling moves allowed whenever divisibility holds, escapes from `(9,6)` to sum `9` in two letters. The true algebra blocks that very word (`X''=D`, `D'=D-X`, sum back to `15`) — **so leading-form data is unavoidable in any proof of Gap 1** |

### 5.2 The winner: the primitive degree `e=deg p`

Dixmier's 1968 leading-symbol lemma (classical input, cited, instance-verified in
file) gives, for `n+m>2`, `σ_X=αp^a`, `σ_D=βp^b` with `p` **primitive**
homogeneous; set `e=deg p`, so `n=ae`, `m=be` and `e | gcd(n,m)`.

> **Invariance (PROVED, verifier `§5.1`).** `e` is constant along every tame
> trajectory all of whose states satisfy `n+m>2`; `p` itself is constant up to the
> linear symplectic action of the affine/Fourier letters.
>
> *Reason.* A transvection changes one member only, so the **other** member's
> symbol — hence its primitive `p` — is untouched, and Dixmier's lemma applied to
> the new pair forces the new symbol to be a power of the *same* `p`. Exchange
> swaps `(a,b)`. Affine/Fourier act linearly on `(x,ξ)`, preserving `deg p` (if an
> affine move cancels the top, the lemma re-applies with the untouched member's
> `p`). Machine-checked on the W2 datum and the tame control for transvections
> `s=1,2,3` on both sides, exchange, and Fourier (`x²ξ ↦ xξ²`).

### 5.3 FLOOR THEOREM

> **Theorem (arbitrary degree; classical input: Dixmier's lemma).** Let `[D,X]=1`
> with primitive degree `e≥2`. Then `a≥2`, `b≥2`, and `(a,b)≠(2,2)`; hence
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

**This closes Gap 1 at the W2 stratum.** The W2 hatch (`σ_X=(x²ξ)³`,
`σ_D=(x²ξ)²`) and the shifted cube share `e=3`, `(a,b)=(3,2)`, `n+m=15=5·3`,
`k=3`: **no composite tame word lowers the invariant from either.** (The corollary
is a statement about the *leading-form stratum*, so it applies to *every* genuine
pair with that stratum, not merely to the specific data.)

### 5.4 Two consequences the campaign should record

- **The band-`k` census hatches, `k≥4`, are NOT at the floor**: `e=3`,
  `(a,b)=(k,k-1)`, `n+m=(2k-1)·3 > 15`. A composite descent to a band-3
  configuration is *not* excluded for them (they die at the moment unit anyway,
  [`hatch-census.md`](hatch-census.md)).
- **The absolute minimal non-generating leading-form stratum is `e=2`, not
  `e=3`**: `n+m=10 < 15`, with `p ~ xξ` up to `SL₂` (a primitive binary quadratic
  has two distinct roots, so a double root is excluded — it would be imprimitive),
  `σ_X=(xξ)²`, `σ_D=(xξ)³`. In coordinates where `p=xξ` the top is carried by the
  **band-0** coefficient (`i+d=2, d=2 ⟹ i=0`); in general the degree-4 top spreads
  over the five bands with `i+2d=4` (the band grading is *not* `SL₂`-invariant, so
  the normalisation cannot be assumed). Either way the `e=2` stratum is a
  **multi-band-top** configuration — exactly what `band-reduction.md` §9 **Gap 3**
  explicitly leaves unclassified — while every hatch the campaign has attacked has
  `e=3`, `p=x²ξ`, single-band-dominated. **The floor ladder ranks lowest precisely
  the stratum that Gap 3 leaves open.**

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
| `(x³-∂, x)` (positive control) | yes | `(4,3)` | **`(2,1)` via `oX`** | 72 (depth 2, `S=3`) |
| tame family `c₁≠0` (positive control) | yes | `(6,3)` | **`(2,1)` via `oD·oX`** | 78 (depth 2, `S=3`) |
| **W2 hatch datum** | **no** — positive cascade `Q_{m>0}=0` and `Q₀=1` hold, the negative tail does not | `(15,3)` | `(15,3)` — **no escape** | 172 (depth 3), 538 (depth 4, HEAVY) |
| **shifted-cube top** `a₃=E(E+1)(E+2)`, `b₂=E(E+1)` | top-only representative | `(15,3)` | `(15,3)` — **no escape** | 172 / 538 |

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

> **Gap 1 holds for a pair `V` if `e(V)≥2` and `{a,b}={2,3}`** — unconditionally,
> by §5.3 (given Dixmier's 1968 lemma and the campaign's band floors).
>
> **Gap 1 for the remaining strata** reduces to: *for a single-move-reduced pair
> `V` with `e(V)≥2` and `a+b>5`, no tame word reaches a state with the same `e` and
> a strictly smaller `(a+b, k)`.* By §3 such a word must contain a cancelling step
> at a state with dividing exponents, and by §3.4 that step's cancellation depth is
> a stack of `p`-power conditions on the sub-leading symbols. **A sufficient
> finite condition:** for each `(a,b)` with `a+b>5` and each admissible `s`,
> show that no Weyl pair with `σ_X=αp^a` admits `D` with `σ_D=βp^{sa}` whose
> sub-leading symbols are `p`-powers to depth `(s+1)a-(a+b)+1`. This is a *finite*
> system per `(a,b,s)` once `e` is fixed — but it is **not** carried out here.
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

So `e≥2` is exactly the *candidate-counterexample* condition at the symbol level,
and DC1 decomposes as **(i)** no Weyl pair has `e≥2`, **and** **(ii)** every `e=1`
Weyl pair generates. §5.3 says a minimal counterexample of type (i) has
`{a,b}={2,3}` and `n+m=5e` — which is precisely the configuration the campaign has
been attacking at `e=3`, and (§5.4) *not* the one it has been attacking at `e=2`.

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
- The short-word escape search: depth 3 committed (172 words/seed), depth 4 under
  `HEAVY=1` (538 words/seed), `S=2` at the band-3 seeds, unit raising coefficients
  only.
- The centraliser-of-`x` computation is a degree-≤4 window.
- Census/wall checks at `k=3,4,5`; `(k,s)∈{(3,2),(3,3),(4,2)}` for the symbolic
  raising-wall identity.

**Refuted:**
- The proposed **local-minimum lemma** ("at every local minimum the pair is
  single-move-reduced") — a quantifier error; only the last-ascent form survives.
  *In-file this is a printed observation, not a witness; a quantifier gap shows
  the proposed proof invalid, not the statement false. The audit constructed the
  missing explicit witness, so the refutation stands on that.*
- `gcd(k,q)`, cofactor effectivity, and non-dividing `(a,b)` as monovariants
  (explicit words, machine-checked).
- **Any monovariant depending only on the `(n,m,k)` transition relation of the
  over-approximate degree model** — that model escapes from `(9,6)`. *Corrected
  statement (audit): this refutes arguments using only that relation; it does not
  refute every conceivable `(n,m,k)`-valued monovariant, since the true algebra
  blocks the escaping move.*
- Top-wall data as a Fourier covariant.

**Open / NOT claimed:**
1. Gap 1 **above** the floor (`e≥2`, `a+b>5`) — the finite condition of §7 is
   stated, not discharged.
2. Gap 1 for **`e=1`** strata (the `κ₂≠0` constant-`h` sector) — untouched.
3. The **`e=2` stratum** (`n+m=10`, `p~xξ`, multi-band top) — newly identified as
   the absolute floor of the non-generating ladder, and completely unexamined; it
   is `band-reduction.md` Gap 3 territory.
4. Whether any Weyl pair with `e≥2` **exists** at all: nothing here constructs one
   or rules one out. No Weyl pair, no DC1 counterexample, no DC1/JC2 progress
   beyond the stated leading-form statements.
5. The `κ₂` gauge obstruction of `shifted-power-residuals.md` §3.1 should be
   re-derived with the corrected Fourier rule of §2.

---

## 9. Verification

```sh
uv run --with sympy python research/dc1-program/verify_gap1_checkpoints.py          # 75 checks, ~45 s
HEAVY=1 uv run --with sympy python research/dc1-program/verify_gap1_checkpoints.py  # + depth-4 search, ~115 s
```

`§0` engine + normal-ordered `A_1` bridge; `§1` invariant + degree floor; `§2`
generator arithmetic + Fourier rule; `§3` checkpoint lemma + cancellation depth;
`§4` wall-data table; `§5` monovariant hunt, refutations, floor theorem, `e=2`
ladder; `§6` bounded word search with two positive controls; `§7` assembled
statement. Default run ends `ALL EXECUTED CHECKS PASSED; OPTIONAL CHECKS SKIPPED`;
under `HEAVY=1`, `ALL GAP1 CHECKPOINT CHECKS PASSED`.
