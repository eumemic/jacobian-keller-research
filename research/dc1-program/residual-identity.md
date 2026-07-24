# The residual identity at W2: the depth-3 tail forces the moment-slope residual W=0

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED — BAND-SCOPED — MIXED PROVED/BOUNDED**

This memo attacks the step-1 lemma of the DC1 program — degree-free slope-forcing at W2:

```text
   positive cascade  AND  Q_-1=Q_-2=Q_-3=0  AND  membership   =>   R(1)=G(1)=0.
```

With the slope gate (`Q_0=1 ⇒ R(1)=1`) this makes the depth-3 tail and the moment unit
mutually exclusive, closing W2 at that degree. The factorization `R(1)=a_2(0)·W`
([`slope-forcing-degree-free.md`](slope-forcing-degree-free.md), re-derived in file)
reduces the lemma to: *the tail forces `a_2(0)·W = 0`*. This memo does four things:

1. **Structural key (PROVED, arbitrary degree).** At GENERAL degree `d`, the depth-3
   ladders `Q_-1,Q_-2,Q_-3` are **linear in the two fillers** — a degree-free
   level-incidence fact. Hence filler elimination is **linear algebra over the cascade
   function field at every degree** — the road around the `d≥4` Gröbner wall.
2. **Product-vs-factor at `d=3` (RESOLVED, and CORRECTING the prior "union" reading).**
   The depth-3 tail forces the **single factor `W=0`** (at `d=3`, `W=-(4/9)am1_3`, so
   `am1_3=0`) — **not** a genuine union `{a_2(0)=0} ∪ {W=0}`, and **not** `a_2(0)=0`.
   `a_2(0)` is unconstrained on `cascade+tail`.
3. **The explicit `d=3` consistency-covector certificate** (with the `8×8` determinant
   saturation the naive elimination drops).
4. **The `d=4` linear-elimination route** — the computation the raw `d=4` Gröbner basis
   OOMs on — reaching a strong sampling verdict that `W` is forced at `d=4` too.

W2 datum, gauge `b_3=0`, quantum band-3 conventions
(`Q_m=sum_(k+l=m)[b_l^[k]a_k - a_k^[l]b_l]`, `f^[n](E)=f(E+n)`,
membership `(E)_j=E(E-1)...(E-j+1) | a_-j, b_-j`):

```text
a_3 = E(E+2)(E+4)   (roots {0,-2,-4}),     b_2 = E(E+3),     b_3 = 0.
branch A:  a_-3=(E)_3 am3,  b_-3=mu_3 a_-3;    branch B:  a_-3=0,  b_-3=(E)_3 C.
```

Exact certificate: [`verify_residual_identity.py`](verify_residual_identity.py) (default run
≈ 5.4 min, `102` checks; `HEAVY=1` ≈ 9.5 min, `111` checks). The all-caps PASS is emitted
only when the exact-`QQ` depth-3 `msolve` payload completes; otherwise the verifier ends
with an explicit SKIP while retaining the exact full-tail SymPy checks. Every load-bearing upstream fact — the
crossed-product engine, `Q_0=(T-1)G`, the slope gate, the both-ends Lemma P, the positive
cascade, the factorization `R(1)=a_2(0)·W` — is **re-derived in file**.

## 0. Headline

> **THE RESIDUAL IS FORCED — AS A SINGLE FACTOR (bounded, `d=3`, exact).** On the W2
> positive cascade the depth-3 tail forces the moment-slope residual `W` to vanish:
>
> ```text
>    Q_-1=Q_-2=Q_-3=0   =>   W = 0        (at d=3: W=-(4/9)am1_3, so am1_3=0),
> ```
>
> machine-checked as `am1_3 ∈ sqrt(cascade + Q_-1..Q_-3)` (msolve `^`, `QQ` + `GF(65003)`),
> and `∈ sqrt(cascade + Q_-1..Q_-5)` by exact **sympy** over `QQ` + two primes (second
> engine). The top boundary value `a_2(0)` is **NOT** forced (explicit witness with
> `a_2(0)≠0`, `R(1)=0`; the slice `{cascade+tail, a_2(0)=1}` is feasible). So
>
> ```text
>    cascade + Q_-1..Q_-3   ⊆   V(W),        a_2(0)  unconstrained,
> ```
>
> and `R(1)=a_2(0)·W=0` via `W=0` — the **second factor**, not `a_2(0)`, and not a union.

This **corrects** [`slope-forcing-degree-free.md`](slope-forcing-degree-free.md) §4/§7,
which recorded "the tail forces the PRODUCT `a_2(0)·am1_3`, a genuine union
`{a_2(0)=0} ∪ {W=0}`" and left "product or factor `a_2(0)=0` — undetermined". Neither
holds: the tail forces the factor **`W`** (`am1_3` at `d=3`) alone, cleanly.

## 1. The structural key — filler-linearity at arbitrary degree (PROVED)

In `Q_m=sum_(k+l=m)[b_l^[k]a_k - a_k^[l]b_l]` a term is **bilinear in the two fillers**
only if it multiplies two filler quantities. Branch-B fillers are `a_-2` (a-slot, level
`-2`) and `b_-3` (b-slot, level `-3`); a term is filler-bilinear iff `k=-2 AND l=-3`, i.e.
`m=-5` **only**. Branch A (fixed gauge `mu_3`, fillers `a_-2` and `a_-3`, `b_-3=mu_3 a_-3`):
bilinear iff `k∈{-2,-3} AND l=-3`, i.e. `m∈{-5,-6}` only. Because this depends only on the
level grid, it is **degree-free**:

> **Lemma (filler-linearity, arbitrary `d`).** For every degree `d`, `Q_-1,Q_-2,Q_-3`
> (indeed `Q_-1..Q_-4` on B, `Q_-1..Q_-4` on A) are **affine-linear** in the fillers; the
> positive cascade `Q_4,Q_3,Q_2,Q_1` is **filler-free**, so the cascade-solved
> `b_1,b_0,b_-1,b_-2` are filler-free. Hence, after the cascade forward-solve, the depth-3
> tail is a linear system `M(free)·x + N(free) = 0` in the fillers `x`, with `M,N`
> rational in the free cascade coordinates.

Machine-confirmed at `d=2,3,4,5` (symbolic generic cascade data, symbolic fillers): every
`Q_-1,Q_-2,Q_-3` coefficient has total degree `≤1` in the filler coefficients, `Q_-5` has a
genuine bilinear (degree-2) term, and `Q_4..Q_1` contain no filler symbol
(`verify §S2`). **Consequence: filler elimination is linear algebra over the cascade
function field at every degree** — the OOMing `d≥4` Gröbner wall is replaced by
determinant/minor linear algebra.

## 2. Product-vs-factor at `d=3`: the tail forces the factor `W` (bounded, exact)

On the branch-B parametrized cascade the free data is 9 coordinates and
`R(1)=-(4/9)a_2(0)·am1_3`; the two fillers are `a_-2=(E)_2 V`, `b_-3=(E)_3 C` (8
coefficients). The decomposition, all machine-checked (`verify §S4`):

- **`W=0` forced.** `am1_3 ∈ sqrt(cascade+tail)`: **sympy** exact over `QQ`, `GF(65003)`,
  `GF(32003)` on the full tail `Q_-1..Q_-5`; **msolve `^`** over `QQ` and `GF(65003)` on
  the depth-3 tail `Q_-1..Q_-3`. Two independent engines, both giving the unit
  Rabinowitsch ideal (a machine-checked radical certificate). Equivalently `W ∈
  sqrt(cascade+Q_-1..Q_-3)`.
- **`a_2(0)` NOT forced.** The slice `{cascade+tail, a_2(0)=1}` is **feasible** (msolve `^`
  non-unit mod `65003` ⇒ feasible over `Q̄`), and an **explicit witness**
  `(a_2(0),am1_3,R(1)) = (2,0,0)` sits on `cascade+tail`. Dually, `{cascade+tail,
  am1_3=1}` is **infeasible** — corroborating `am1_3=0`.

So the tail variety lies in `V(W)=V(am1_3)` with `a_2(0)` free: the forcing is the single
factor `W`, and `R(1)=a_2(0)·W=0` because `W=0`. The "genuine union" reading is refuted.

## 3. The consistency-covector certificate at `d=3` (det-saturated)

The depth-3 tail is `24` scalar equations, all linear in the `8` fillers: `M·x + N = 0`,
with `M` a `24×8` matrix over the 9 free coordinates. `M` has **full column rank `8`**
(cokernel dimension `16`). Pick `8` independent (low-degree) rows `I`; `L_I` is the `8×8`
filler block, `det_I` its determinant, `adj(L_I)` its polynomial adjugate. For each of the
other `16` rows `j`, the **explicit covector**

```text
   mu_j = det_I·e_j - M_j·adj(L_I)·(rows I),      mu_j·M = 0   (identically),
```

annihilates every filler column (machine-checked as a polynomial identity), and its
pairing with the non-filler residual is the consistency condition

```text
   G_j = mu_j·N = det_I·N_j - M_j·adj(L_I)·N_I     (the rank[M|N]=rank M minor).
```

These `G_j` (12 nonzero, total degree `15–16`) **carry the `8×8` determinant saturation**
that the naive filler-elimination drops (the `slope-forcing-degree-free.md` warning). The
residual is located inside the saturated consistency ideal:

```text
   a_2(0)·W  ∈  sqrt( <G_1,…,G_12> : det_I^∞ )         (msolve '^' mod 65003, HEAVY),
```

i.e. `a_2(0)·W` vanishes on `{G_j=0} ∖ V(det_I)` — the det≠0 chart of the elimination
variety. The `det_I=0` locus is covered chart-free by §2 (`am1_3 ∈ sqrt(cascade+tail)`).
This is the finite certificate whose **structure** the degree-free argument must
generalize. The fixed `(a_3,b_2)` root necklace belongs to the joint `Q_0`/`Phi`
operator (degree-free, [`lambda-general-k.md`](lambda-general-k.md) Thm A′), not to
the pure `Q_-1,Q_-2,Q_-3` tail. The pure-tail necklace is algebraic after solving
the cascade; its integer nodes are only the membership-window nodes (`{0,1}` for
`a_-2=(E)_2V` and `{0,1,2}` for `b_-3=(E)_3C`). Its datum-dependent
`(a_2,b_1)` block is the same open block as the joint-covector lemma
([`joint-covector.md`](joint-covector.md) §7).

## 4. The `d=4` linear-elimination route (bounded; the road around the GB wall)

The **raw** `d=4` tail-forcing Rabinowitsch GB is intractable: `msolve` (even with the
`-e` filler-elimination order) is memory-bound and does not terminate; sympy is slower
still. The **linear** route removes the fillers first:

- The depth-3 tail is `30` rows **linear** in the `10` fillers, `M` of **full column rank
  `10`** (§1 at `d=4`) — elimination is linear algebra.
- **`W` forced (strong sampling evidence).** Among many random cascade points with `W≠0`,
  the depth-3 tail is **never solvable** (`rank[M|N] > rank[M]`, tested over `GF(2^31-1)`
  with `QQ` re-confirmation) — so tail-solvability forces `W=0` with probability one. This
  is the `d=4` analogue of §2 and the `d=4` **verdict**.
- The consistency minors are **computable fraction-free mod `p`** (`DomainMatrix` over
  `GF(p)[free]`), degree `17`. The det-saturated Nullstellensatz GB is degree `~18` and is
  **not confirmed in budget** — reported as such, not claimed.

So the linearity delivers a decisive `d=4` W-forcing verdict where the raw GB could not,
while the exact `d=4` Nullstellensatz certificate remains out of budget.

## 5. The mirror and branch A

- **Mirror.** `R(-1)=-R(1)` on the cascade (`Q_0(0)=Q_0(-1)`; the reflected 3-term Lemma-P
  boundary pairing at `E=-1` is filler-independent, and the cascade relation
  `R(1)+R(-1)=0` is verified on the parametrized cascade at `d=1,2,3`, `+4` HEAVY). So the
  slope is a single scalar and the tail forces **both ends** to zero.
- **Branch A per gauge `mu_3 ∈ {1,-2,1/3}`.** The depth-3 tail is linear in the fillers
  `(a_-2,a_-3)`; `am1_3 ∈ sqrt(cascade+Q_-1..Q_-5)` (W forced); `a_2(0)` not forced (slice
  feasible + explicit witnesses `a_2(0)∈{2,4,2}`, `R(1)=0`). Same decomposition as branch B.

## 6. Evidence ledger — proved / bounded / refuted / open

**Proved (arbitrary coefficient degree, char 0, machine-checked identities):**
- `Q_m=[D,X]_m` for `m∈[-6,6]`; `Q_0=(T-1)G` (`§S0`). **In-file machine scope
  (audit): generic degree 2** (the identities are degree-free code-equivalence /
  telescoping; audit independently corroborated at `d=6`).
- Slope gate `D|(E-R) ⇔ R(1)=1,R(-1)=-1`; the factorization `R(1)=a_2(0)·W` (`§S1`,
  symbolic).
- **Filler-linearity at arbitrary `d`** (the structural key): `Q_-1,Q_-2,Q_-3` linear in the
  fillers, the positive cascade filler-free — degree-free level-incidence proof (both
  branches, branch A per fixed gauge `mu_3`); the machine loop covers **branch B at
  `d=2,3,4,5`** and **branch A at `d=3`** (`§S2`,`§S6`; audit independently confirmed both
  branches at `d=6`). Filler elimination is linear algebra over the
  cascade function field at every degree.

**Bounded-in-file, degree-free by structure (AUDIT-RETIERED from "proved"):**
- Both-ends Lemma P `R(1)=Q_0(0)`, `R(-1)=-Q_0(-1)`, both **filler-independent**: the
  in-file machine check is at **fixed `d=4`** (`§S1`); the fact is degree-free (membership
  vanishing) and audit confirmed it independently at `d=5,6`, but per house rule the
  in-file evidence tier is bounded.

**Bounded-finite (exact scope stated):**
- `R(1)=a_2(0)·W` on the parametrized cascade at `d=1,2,3` (`+4` HEAVY); mirror
  `R(-1)=-R(1)` there; `W=-(4/9)am1_3` at `d=3`; `R(1)` a free modulus on the cascade alone
  (`§S3`).
- **`d=3` product-vs-factor:** the depth-3 tail forces the **factor `W=0`**
  (`am1_3 ∈ sqrt(cascade+tail)`); `a_2(0)` NOT forced (msolve slice feasible + explicit
  witness) (`§S4`); both branches (`§S6`). **Engine scope (audit-clarified): the
  FULL-tail (`Q_-1..Q_-5`) certificate is sympy, exact `QQ` + 2 primes; the DEPTH-3
  (`Q_-1..Q_-3`) certificate is msolve-only, exact `QQ` + one prime — both depths carry
  an exact-`QQ` certificate, but the depth-3 statement is single-tool.**
- **`d=3` consistency-covector certificate:** explicit covectors `mu_j` (`mu_j·M=0`
  identically), the det-saturated conditions `G_j` (rank minors, degree `15–16`), and
  `a_2(0)·W ∈ sqrt(<G_j>:det_I^∞)` (msolve `^`, HEAVY) (`§S5`).
- **`d=4` (HEAVY):** depth-3 tail linear, `M` full column rank `10`; `W` forced by sampling
  (no `W≠0` tail-solvable witness); consistency minors computable fraction-free mod `p`
  (`§S7`).

**Refuted / corrected:**
- [`slope-forcing-degree-free.md`](slope-forcing-degree-free.md) §4/§7 —
  *"the tail forces the PRODUCT `a_2(0)·am1_3`, a genuine union `{a_2(0)=0} ∪ {W=0}`"* and
  *"product or factor `a_2(0)=0` — undetermined"*. **Corrected:** the tail forces the single
  factor `W` (`am1_3` at `d=3`); `cascade+tail ⊆ V(W)` with `a_2(0)` unconstrained. Not a
  union, and not `a_2(0)=0`.

**Open / not claimed:**
- The residual identity at **arbitrary positive-data degree** — `W ∈ sqrt(cascade+Q_-1..Q_-3)`
  for all `d`. Exact at `d=3` (both branches); strong **sampling** evidence at `d=4`; no
  degree-free certificate. `W` grows with `d`, so this needs the degree-free covector for
  the algebraic `(a_2,b_1)` necklace block — the same gap as the joint-covector lemma.
- The `d≥4` det-saturated **Nullstellensatz GB** (degree `~18`) — not confirmed in budget.
- Whether the `d=4` W-forcing sampling verdict is an exact identity (it is not a
  Nullstellensatz certificate). No Weyl pair, no DC1 counterexample; all of Band 3, DC1,
  JC2 remain open. The infinite-dimensional `Im L_K ∩ Im L_H`
  ([`two-filler-cross-cancellation.md`](two-filler-cross-cancellation.md)) is untouched.

## 7. Verification

```sh
uv run --with sympy python research/dc1-program/verify_residual_identity.py
# HEAVY (d=4 linear route, covector saturation msolve, branch-A depth-3 msolve):
HEAVY=1 uv run --with sympy python research/dc1-program/verify_residual_identity.py
```

`§S0` engine; `§S1` slope gate + both-ends Lemma P + mirror + factorization; `§S2` the
degree-free filler-linearity structural key (`d=2,3,4,5`, both branches); `§S3`
factorization + mirror + `W`-collapse on the parametrized cascade; `§S4` product-vs-factor
at `d=3` (`W` forced two engines when `msolve` is available, `a_2(0)` free via slice +
witness); `§S5` the det-saturated consistency-covector certificate; `§S6` branch A per
gauge; `§S7` (HEAVY) the `d=4` linear-elimination route. The final status is PASS only
when the exact-`QQ` depth-3 payload runs; otherwise it is SKIP with supporting/full-tail
checks passed.
