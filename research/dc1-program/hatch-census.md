# The hatch census: one silent hatch per band, and the band-4/5 kill at the moment unit

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER REVIEWED — BAND-SCOPED — MIXED PROVED/BOUNDED**

This memo audits architecture step 2's load-bearing prerequisite — *"one silent
hatch per band"* — and tests the combined feasibility of the band-4 and band-5
hatches (the analogues of W2). It **re-derives** every load-bearing upstream fact
in file rather than citing it, and in doing so it **corrects two committed
claims** (one in `../band3/w2-theory.md`, one method bug in the two sibling
`msolve` verifiers).

Conventions frozen from the band-`k` corpus (`b7e85e8..30d8c59`): `A_1[x^{-1}] =
(+)_k x^k C[E]`, `(x^a f)(x^b g)=x^{a+b} f(E+b)g(E)`, `f^{[n]}(E)=f(E+n)`,
`Q_m=sum_{i+l=m}(b_l^{[i]}a_i-a_i^{[l]}b_i)`, gauge `b_k=0`, `u:=b_{k-1}`,
membership `(E)_j=E(E-1)...(E-j+1) | a_{-j},b_{-j}`, `Q_0=(T-1)G`, `G(0)=0`.

The canonical band-`k` **hatch** (translate the unique common root onto the
anchor `0`):

```text
a_k = prod_{i=0}^{k-1} (E + i(k-1))   (step-(k-1) AP, roots {0,-(k-1),...,-(k-1)^2}),
u   = prod_{j=0}^{k-2} (E + j k)      (step-k AP, roots {0,-k,...,-(k-2)k}).
```

Band 3 is **W2** (`a_3=E(E+2)(E+4)`, `u=E(E+3)`). Band 4: `{0,-3,-6,-9}/{0,-4,-8}`.
Band 5: `{0,-4,-8,-12,-16}/{0,-5,-10,-15}`.

Exact certificate: [`verify_hatch_census.py`](verify_hatch_census.py) (ends
`ALL HATCH CENSUS CHECKS PASSED`). Base commit `30d8c59`.

## 0. Headline

> **The band-`k` hatches (`k>=4`) are killed by the MOMENT UNIT ALONE — a
> stronger and earlier kill than W2 — because the census's "reduces to the same
> slope gate `R(1)=1`" is FALSE for `k>=4`.**
>
> 1. **Census, corrected.** Every hatch (bands 3, 4, 5) has `codim Im Phi = 3`
>    (honest left-nullspace, cross-checked three ways). Band 3 (W2) is
>    **point-complete** — `Ann = {ev_-1, ev_0, ev_1}`, `Im Phi = (E(E-1)(E+1))` —
>    so `Q_0=1 <=> R(1)=1` and W2 survives the moment unit. But bands 4, 5 have
>    **only two point annihilators `{ev_0, ev_1}` plus one NON-POINT annihilator**;
>    `Im Phi != (E(E-1))` and `w2-theory.md`'s "band-4 hatch: `Im Phi=(E(E-1))`,
>    codim 2, principal, reduces to the same slope gate" is **wrong** (true codim
>    `3`). The `gcd`-shortcut `im_phi_gen` it used is unreliable (it assumes the
>    image is principal).
> 2. **The kill.** For the band-4 and band-5 hatches the moment unit `Q_0=1`
>    (positive cascade `∧ Q_0=1 ∧` membership) is the **UNIT IDEAL** on both tail
>    branches, at `d=2` (`msolve` + independent `sympy` cross-check) and `d=3`
>    (band 4: two engines; band 5: `msolve`). Since `cascade+Q_0=1` UNIT `=>` the FULL system
>    (`+` negative tail) is UNIT a fortiori, **the negative tail is not even
>    needed** — the hatch dies at the moment unit.
> 3. **Why (band 4).** The pure slope `R(1)=1` **is** achievable at `d=3`, yet
>    `Q_0=1` is still UNIT — so the killer is the **non-point annihilator**, an
>    obstruction that does not exist at W2. W2 (band 3) is the *weakest* case of
>    the tower: it alone needs the joint slope+tail argument
>    ([`w2-joint-theorem.md`](w2-joint-theorem.md)); `k>=4` is more decisively dead.
> 4. **Existence PROVED, uniqueness BOUNDED.** The step-`(k-1)`/step-`k` hatch is
>    wall-admissible (base-representation collapse) and has a unique common root at
>    the anchor (`gcd(k-1,k)=1`) for **all** `k` — a theorem. Its uniqueness among
>    minimal single-coset tops is **verified over a large window** (bands 3–5), not
>    proved in general.
> 5. **No hatch survived** the moment unit → **no counterexample lead** → the
>    escalation protocol was **not** triggered; no Weyl pair and no DC1
>    counterexample is constructed.

## 1. Task 1 — the census audit

### 1.1 Existence of the hatch (THEOREM, arbitrary `k`)

The top wall `Q_{2k-1}=0` is, on the root necklace (`sigma=T^{-1}`, root at `E=-p
<-> sigma^p`), `S_k(sigma) δ(u) = S_{k-1}(sigma) δ(a_k)`, `S_r=1+sigma+...+sigma^{r-1}`
([`../band3/band-k-weapons.md`](../band3/band-k-weapons.md), re-verified).

> **Theorem (hatch, arbitrary `k`).** For the step-`(k-1)`/step-`k` AP
> `δ(a_k)={i(k-1)}_{i<k}`, `δ(u)={jk}_{j<k-1}`:
> - *wall-admissible:* `S_k δ(u) = S_{k-1} δ(a_k) = S_{k(k-1)}` (the multiset
>   `{0,1,...,k(k-1)-1}`), by the base-`k` / base-`(k-1)` two-digit collapse — a
>   proof for every `k` (checked `k=3..7`).
> - *unique common root at the anchor:* `i(k-1)=jk` with `gcd(k-1,k)=1` forces
>   `k|i`, so `i=0`; the ONLY common root is `E=0`.
> - *exotic* (`a_k` non-consecutive) for `k>=3`.

This is the census's positive claim, and it holds as stated (confirms
`w2-theory.md` §3).

### 1.2 Uniqueness (OBSERVATION, strengthened, not a theorem)

A **hatch** is a wall-admissible exotic top with a *unique* common root, its
root translated onto the anchor. `w2-theory.md` verified uniqueness in a small
window (support `<=6`, coeffs `{-1,0,1}`). We strengthen this to an exhaustive
scan of **all minimal squarefree single-coset wall-admissible tops** (top with
exactly `k` distinct integer roots, `S_k | δ(a_k)`, `S_{k-1}·(δ(a_k)/S_k)`
effective) in a large window:

| band | window `[0,W]` | # admissible tops | # unique-common-root hatches |
|---|---|---|---|
| 3 | `[0,22]` | 2 | **1** (`{0,2,4}/{0,3}` = W2) |
| 4 | `[0,20]` | 5 | **1** (`{0,3,6,9}/{0,4,8}`) |
| 5 | `[0,22]` | 14 | **1** (`{0,4,8,12,16}/{0,5,10,15}`) |

Every non-hatch admissible top has `>=2` common roots (the tame consecutive
"shifted-power" top has `k-1` of them). So *within minimal squarefree
single-coset tops*, the hatch is unique **over this window** — bounded evidence.
A general uniqueness theorem is **open**: it is a nontrivial statement about
`S_k`-divisible `0/1` necklaces with a unique `S_{k-1}(δ(a)/S_k)`-overlap, and we
do not have it.

### 1.3 CORRECTION — the E-pairing does NOT "reduce to just `R(1)=1`" for `k>=4`

`w2-theory.md` §3 states the band-4 hatch has `Im Phi=(E(E-1))` (principal,
squarefree, codim 2) and that `E-R in Im Phi <=> R(1)=1`, "the same slope gate",
making the general-`k` hatch "a bona-fide W2-analogue". **This is wrong for
`k>=4`.** Computing `codim Im Phi` honestly (rank of the image-coefficient matrix
= left-nullspace dimension; `verify §2`), cross-checked by the point-annihilator
count and by explicit ideal membership:

| band | honest `codim Im Phi` | point annihilators | non-point annih? | `Im Phi` |
|---|---|---|---|---|
| 3 (W2) | **3** | `{ev_-1, ev_0, ev_1}` (3, complete) | no | `(E(E-1)(E+1))` |
| 4 | **3** | `{ev_0, ev_1}` (2 only) | **yes (1)** | `⊊ (E(E-1))` |
| 5 | **3** | `{ev_0, ev_1}` (2 only) | **yes (1)** | `⊊ (E(E-1))` |

For `k>=4` the point annihilators **undercount** (`dim 2 < codim 3`), exactly the
"single-block" caveat of
[`lambda-general-k.md`](lambda-general-k.md) §2 — but now for the *sum* `Im Phi`.
`E(E-1) not in Im Phi` (checked), so `Im Phi` is a **codimension-1 subspace of**
`(E(E-1))`, cut out by a **non-point (infinite-support) annihilator**. Hence

```text
Q_0 = 1   (k>=4)   <=>   R(0)=0 (auto)  AND  R(1)=1  AND  lambda_np(E-R)=0,
```

a **third** condition beyond the slope, absent at W2. The `gcd`-shortcut
`im_phi_gen` used by `w2-theory.md` is **unreliable**: it takes `gcd` of the
images and tests principal surjectivity, so it silently assumes `Im Phi` is
principal. On the un-translated band-4 top `{0,3,6,9}/{1,5,9}` it returns
`gcd = E` (claiming codim 1) while the true codim is `2` (`verify §2`). Its
band-4 *hatch* answer `(E(E-1))` is off by the same mechanism — codim 3, not 2.

### 1.4 Reconciliation with `band45-lambda.md`, and the residual-identity dependence

`band45-lambda.md`'s normalized top `{0,3,6,9}` (cofactor `C=Φ_6Φ_12`,
`lambda(E)=7`) is the **same** top as the band-4 hatch, un-translated: its wall
sub is `{1,5,9}` and its unique common root sits at the **max** position `9`, not
at the anchor. Translating that common root onto the anchor is precisely what
produces the hatch. So the two memos analyse the same top at two different
translations relative to the fixed membership anchor — different DC1 problems.

The census criterion "unique common root **at the anchor** ⟺ E-pairing
degenerates" therefore rests on the interplay between the common root and the
anchor: at the anchor the obstructing functional's `lambda(E)` drops toward the
slope, off the anchor it is `lambda(E)=rho* != 0`. But whether an off-anchor
common root **keeps** the obstruction alive is the statement
`lambda(E-R)=rho*-lambda(R) != 0`, i.e. the **residual identity `lambda(R)=0`** —
which is proved only at W1 (`e4e704f`) and is the *open* cascade-dependent step
for every `k>=4` (`band45-lambda.md` §5, `lambda-general-k.md` §6). **So the
"unique common root" criterion of the census is conditional on that open
identity, not a theorem.** What *is* unconditional is the honest codimension and
the moment-unit verdict below.

## 2. Task 2 — the band-4 and band-5 hatches, explicitly

```text
band 4:  a_4 = E(E+3)(E+6)(E+9),   u = b_3 = E(E+4)(E+8),   b_4 = 0.
band 5:  a_5 = E(E+4)(E+8)(E+12)(E+16),  u = b_4 = E(E+5)(E+10)(E+15),  b_5 = 0.
```

Both satisfy the top wall exactly; both have a unique common root at the anchor;
both are exotic; both have `codim Im Phi = 3` (`verify §1,§2`).

**Annihilator silence, corrected.** At each hatch the two *point* annihilators are
silent on `E-R` in the sense that `ev_0(E-R)=0` (automatic) and `ev_1(E-R)=1-R(1)`
vanishes at the achievable slope `R(1)=1`. But the **non-point** annihilator is
**not silent** — it carries the extra obstruction. So the "annihilator silence"
picture of `w2-theory.md` holds only for the point part; the full dual is not
silent for `k>=4`.

## 3. Task 3 — combined feasibility: the hatch dies at the moment unit

We test, exactly as `../band3/w2-verdict.md` did for W2, the systems (both
exhaustive `Q_{-2k}`-branches A: `a_{-k}=(E)_k·am_k`, `b_{-k}=mu_k a_{-k}`; B:
`a_{-k}=0`, `b_{-k}=(E)_k·bm_k`):

```text
cascade          Q_{2k-2}=...=Q_1=0            (control)
cascade + Q_0=1  moment unit                    (the slope + non-point conditions)
FULL             cascade + Q_0=1 + Q_-1..Q_-(2k-1) + membership
```

### 3.1 The escape modulus (when does the slope free up?)

`R(1)=G(1)` is the moment slope. On the positive cascade:

| band | `R(1)` forced `0` at | `R(1)` a free modulus at |
|---|---|---|
| 3 (W2) | `d<=2` (proved, `w2-decisive.md`) | `d=3` (proved) |
| 4 | `d<=2` | `d=3` (pure slope `R(1)=1` achievable — `msolve`; `G(1)` a live modulus, `verify §4`) |
| 5 | `d<=3` | `> 3` (still pinned at `d=3`: pure slope `R(1)=1` is UNIT — `msolve`) |

Band 4's pure slope frees up at the **same** degree as W2 (`d=3`); band 5's is
slower (still pinned at `d=3`). This is the escape-modulus analysis the task asks
for: the slope is a genuine free modulus (band 4 at `d=3`), not a permanent zero.
The certificate checks the degree-free part — that `G(1)` gains free moduli from
`d=2` to `d=3`; the exact `R(1)=1` achievability (band 4) and pinning (band 5) at
`d=3` are the `msolve` results recorded here.

### 3.2 The verdict — UNIT at the moment unit (both bands, both branches)

> **The band-4 and band-5 hatch moment units `cascade+Q_0=1` are the UNIT ideal**
> (both tail branches):
> - `d=2`: exact `QQ`, `msolve` (`^`) with an independent `sympy` cross-check
>   (feasibility-preserving linear reduction) on band 4, both bands.
> - `d=3`, band 4: exact `QQ`, `sympy` **and** `msolve`; band 5: `msolve` (`^`).
>
> **Control (no false kill).** `cascade` alone is FEASIBLE (band 4, `d=2`, `sympy`).
> And for band 4 at `d=3` the **pure slope** `R(1)=1` (the single scalar
> `Q_0(0)=1`) is achievable (`msolve`) — so the moment-unit kill is **not** the
> slope being pinned; the extra killer is the non-point annihilator.

Because `cascade+Q_0=1` is UNIT, the FULL system (which *adds* the negative tail)
is UNIT a fortiori — **the hatch dies before the tail is reached.** This is
strictly stronger and earlier than the W2 kill (which is UNIT only for the
*combined* slope+tail system; `w2-verdict.md`).

### 3.3 Consequence

No band-4 or band-5 hatch admits a slope-`1` datum satisfying the moment unit at
`d<=3` (band 4 to `d=3` exactly; band 5 to `d=3`; `d<=2` both). No candidate DC1
pair materializes, so `[D,X]=1`, the char-`p` sieve
([`sieve_dc1_candidate.py`](sieve_dc1_candidate.py)) and the `A_1`-generation
test were **not** run — there is nothing to feed them. **No counterexample lead;
escalation not triggered.** This is a **bounded-degree** kill (`d<=3`); arbitrary
degree needs the degree-free non-point covector (§4), which is open.

## 4. Task 4 — the uniform pattern for step 2

Reading with the sibling [`w2-joint-theorem.md`](w2-joint-theorem.md), which
localizes W2's kill:

- **Band 3 (W2)** — the annihilator dual is the **symmetric point triple**
  `{ev_-1, ev_0, ev_1}` (roots of `D=E(E-1)(E+1)`). The proved cascade relation
  `R(1)+R(-1)=0` collapses the three conditions to the single **achievable**
  `R(1)=1`, so `Q_0=1` is **feasible**. W2 survives the moment unit; its kill is a
  **joint slope+tail filler obstruction** (a Fredholm gap `=1` between the fillers
  `Q_0=1` demands and the fillers the tail allows — `w2-joint-theorem.md`).
- **Band `k>=4`** — the anchor neighborhood supports only `{ev_0, ev_1}` as point
  annihilators; the third annihilator is **non-point** and its `E-R` pairing is
  **not** governed by any cascade symmetry. `Q_0=1` is therefore **infeasible**,
  and the hatch is killed at the **moment unit alone** — the negative tail plays
  no role.

> **Uniform pattern (corrected).** The singular-hatch obstruction is
> `Ann(Im Phi)` paired against `E-R`. It **strengthens with band index**: at `k=3`
> it is a symmetric point triple that collapses to an achievable slope (needing
> the joint tail to kill), and at `k>=4` it acquires a non-point covector that
> kills `Q_0=1` outright. W2 is the **weakest** rung; the tower gets *more*
> decisively dead, not uniformly-the-same. This **corrects** both
> `w2-theory.md` §5 ("every hatch reduces to the same slope gate `R(1)=1`") and
> `w2-joint-theorem.md` §6.2 (the `{Q_0=1, Q_-1}` generic kill "the natural
> band-`k` conjecture"): for `k>=4` the kill is already at `{Q_0=1}`.

What step 2 needs to make this a *theorem* at arbitrary degree is the **degree-free
non-point covector** `lambda_np` with `Im Phi ⊂ ker lambda_np` and
`lambda_np(E-R) != 0` on the cascade — the `k>=4` analogue of the `lambda_r`
cofactor functional, now genuinely off the point lattice. Its existence
degree-free is open (it inherits the same residual-identity gap as `lambda(R)=0`).

## 5. Method finding — the `msolve` `**` parser bug (re-verify upstream)

`msolve 0.10.1` **misparses Python `**` exponent notation**: it silently corrupts
the ideal and can return a **wrong "feasible" verdict**. With `^` it is correct.

> **Minimal reproducer (`verify §5`).** `{3y0^3-2y0^2 y1+9 y0 y4,
> -9y0^5+12y0^4 y1-4y0^3 y1^2-54y0^3 y4+36y0^2 y1 y4-81y0 y4^2-2916}` is provably
> the **unit ideal** (`y0(3y0^2-2y0 y1+9y4)=0`; either `y0=0` or the factor
> vanishes, and in both cases the quintic reduces to `-2916 != 0`), confirmed by
> `sympy` over `QQ` and `GF(p)`. `msolve` with `**` returns **feasible** (a
> spurious 2-element Gröbner basis); with `^` it returns `[1]` (unit).

Both [`../band3/verify_w2_verdict.py`](../band3/verify_w2_verdict.py) (line 361)
and [`verify_w2_joint.py`](verify_w2_joint.py) (line 570) build the `msolve`
input via `str(...).replace(" ","")` — emitting `**`. The bug is **insidious**:
robustly-unit systems survive it (the W2 `d=3` FULL system returns UNIT under
`**` too, and is independently `sympy`-confirmed), but **delicate or feasibility
verdicts can flip**. So:

- The sibling **`d=3` FULL kills** stand (robustly unit; `sympy`-cross-confirmed
  in-file).
- Their **`msolve`-only** results — the `d=4`/`d=5` verdicts (`w2-verdict.md`),
  and the *feasibility* sub-claims of `w2-joint-theorem.md` §S8 (e.g.
  "`cascade+Q_0=1+Q_-1` is NOT unit", "`R(1) not in sqrt(cascade+tail)`") — should
  be **re-run with the `^` fix** before being treated as load-bearing.

`verify_hatch_census.py` uses `^` throughout and cross-checks every kill with an
independent `sympy` Gröbner (feasibility-preserving linear reduction).

## 6. Claim disposition

**Proved (exact algebra, arbitrary `k` / arbitrary degree as stated):**
- The step-`(k-1)`/step-`k` hatch is wall-admissible (base-rep, all `k`), has a
  unique common root at the anchor (`gcd(k-1,k)=1`, all `k`), exotic `k>=3`.
- `Q_m=[D,X]_m` (`m in [-2k,2k]`) and `Q_0=(T-1)G`, `G(0)=0`, at `k=3,4,5`.

**Proved (exact finite, bands 3,4,5):**
- `codim Im Phi = 3` for every hatch (left-nullspace, stable in truncation);
  band 3 point-complete `{ev_-1,ev_0,ev_1}`, `Im Phi=(E(E-1)(E+1))`; bands 4,5
  carry a non-point annihilator, `E(E-1) not in Im Phi`, so `Im Phi != (E(E-1))`
  (correcting `w2-theory.md`).
- The `gcd`-shortcut `im_phi_gen` undercounts codim (demoed on `{0,3,6,9}/{1,5,9}`).
- **Band 4, 5 hatch moment unit `cascade+Q_0=1` is the UNIT ideal**: `d=2` both
  bands both branches (`msolve` `^`; `sympy` cross-check on band 4 branch A),
  `d=3` band 4 both branches (`sympy` branch A + `msolve`) and band 5 both
  branches (`msolve` `^`). Pure slope `R(1)=1` achievable at `d=3` (band 4,
  `msolve`) — the killer is the non-point annihilator. Cascade control feasible
  (`sympy`; no false kill).
- The `msolve` `**` bug (minimal reproducer, both directions).

**Bounded / observational evidence:**
- Uniqueness of the hatch among minimal squarefree single-coset tops (windows
  `W<=22`, bands 3–5): exactly one per band.
- The full-tail FULL system UNIT follows from `cascade+Q_0=1` UNIT; not separately
  Gröbnered at `d=3` (unnecessary — a superset of an infeasible system).
- Band-5 pure-slope freeness threshold is `>3` (still pinned at `d=3`); the band-5
  `d>3` moment unit is not pushed (it is UNIT already at `d<=3`).

**Open / not claimed:**
- General uniqueness of the hatch (arbitrary `k`, arbitrary window / non-single-
  coset / higher-degree tops).
- The residual identity `lambda(R)=0` for the non-point covector at `k>=4` — the
  cascade-dependent step that would make the moment-unit kill **degree-free** (it
  is only bounded here, `d<=3`).
- Arbitrary-degree status of any hatch; DC1; JC2. No Weyl pair, no counterexample.

**Exceptional loci (explicit):** the hatch exists at every `k>=3` (needs only
`gcd(k-1,k)=1`); branch A requires `a_{-k} != 0` (the `mu_k` gauge), branch B is
`a_{-k}=0`; both are covered and exhaustive. The `msolve` bug bites only
delicate/feasibility verdicts, not robustly-unit ones.

## 7. Verification

```sh
uv run --with sympy python research/dc1-program/verify_hatch_census.py
```

Exact `sympy` (+ `msolve` with the `^` fix where present). §0 crossed-product
engine; §1 hatch existence (base-rep + `gcd`, `k=3..7`); §2 honest `codim Im Phi`
and the `w2-theory` correction; §3 census uniqueness (windowed enumeration); §4
the moment-unit kills (bands 4,5, both branches, `sympy`+`msolve`) with cascade
and pure-slope controls; §5 the `msolve` `**` bug. Ends
`ALL HATCH CENSUS CHECKS PASSED`.
