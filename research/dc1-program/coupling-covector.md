# The coupling covector at W2: the varying-tops coupling collapses onto a FIXED INTEGER-NODE grid — and the fixed-shape minors still miss `W`

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED — BAND-SCOPED — MIXED PROVED/BOUNDED/REFUTED**

The prize this memo attacks ([`algebraic-covector.md`](algebraic-covector.md) §5,
[`residual-identity.md`](residual-identity.md) §6): a **fixed finite recipe** producing, at
every positive-data degree `d`, a covector annihilating the depth-3 filler columns whose
residual pairing is a **unit multiple of `W`**. That would close W2 arbitrary-degree
(tail ⇒ `R(1)=a_2(0)·W=0` contradicts the slope gate `Q_0=1 ⇒ R(1)=1`).

W2 datum, gauge `b_3=0`, quantum band-3 conventions
(`Q_m=sum_(k+l=m)[b_l^[k]a_k − a_k^[l]b_l]`, `f^[n](E)=f(E+n)`,
membership `(E)_j=E(E−1)...(E−j+1) | a_-j,b_-j`):

```text
a_3 = E(E+2)(E+4)  (roots {0,-2,-4}),   b_2 = E(E+3),   b_3 = 0,
branch B: a_-3 = 0,  fillers  P := a_-2 = (E)_2 V,   S := b_-3 = (E)_3 C.
```

Exact certificate: [`verify_coupling_covector.py`](verify_coupling_covector.py) —
default run **~46 s / 166 checks** (one optional leg skipped), `HEAVY=1` **~115 s / 175 checks** (no skips); runtimes are environment-dependent. Every load-bearing
upstream fact (engine, `Q_0=(T−1)G`, slope gate, both-ends Lemma P, `R(1)=a_2(0)·W`) is
**re-derived in file**. The final status line distinguishes an unskipped run from one with
skipped optional legs.

## 0. Headline

> **THE COUPLING COLLAPSES; THE FIXED-SHAPE MINORS DO NOT CARRY `W`.**
> The "varying tops" obstruction of [`algebraic-covector.md`](algebraic-covector.md) §5 is
> **dissolved**: the depth-3 tail, read at **integer nodes**, is a fixed-shape linear system
> in **finitely many filler VALUES** — `(3k−4) × (2k−1)` for every `k ≥ 3` at **every**
> degree `d`, with entries the cascade datum at **fixed integer nodes**. No algebraic
> necklace, no trace forms, no `gcd(a_2(E),a_2(E−3))` condition: the two blocks are a
> **staircase**, six `(rung,node)` pairs are identically silent, and the cascade puts the
> `b`-tops on the **same** integer grid by a degree-free 3-step recursion anchored at the
> `a_3`-roots `{0,−2,−4}`.
> But the recipe does **not** close: at `d=3` **no** maximal minor of the fixed-shape window
> (`k=3,4,5`) is divisible by `W`, and an **explicit exact rational witness** shows the
> `k=3` window is satisfiable with `W = −8/9 ≠ 0`. So the depth-3 W-forcing is **not** a
> finite-integer-node phenomenon at the tested widths.
>
> **CORRECTED SHARPENING (audit, 2026-07-25; §5 has the full statement).** The earlier
> reading — "the forcing consumes the *polynomiality* of the fillers" — is **wrong**: the
> evaluation map onto the window values is **surjective for `k ≤ d+2`**, so at every tested
> width (`k=3,4,5`, `d=3`) the relaxation is **vacuous** and polynomiality is never probed
> (the audit realized the `k=3` witness values by honest degree-`≤3` polynomials). The real
> content: the `3k−4` node equations are **too few** — the true tail is inconsistent at that
> same point. The first width where polynomiality could matter is **`k = d+3`**, untested,
> and a width growing with the degree is itself an obstacle to a degree-free recipe.

## 1. The staircase collapse of the two blocks (PROVED, degree-free)

The two-term block operators of [`algebraic-covector.md`](algebraic-covector.md) §2, written
in the filler values `P=a_-2`, `S=b_-3` (`verify §S2`, generic data at `d=2,3,4,5`):

```text
   a-block(Q_m) = b_{m+2}(E-2) P(E) - b_{m+2}(E) P(E+m+2),
   b-block(Q_m) = a_{m+3}(E) S(E+m+3) - a_{m+3}(E-3) S(E).
```

At `m=−2` the a-block **degenerates to a pure multiplier** `[b_0(E−2)−b_0(E)]·P(E)`
(no shift), and at `m=−3` the b-block degenerates to `[a_0(E)−a_0(E−3)]·S(E)`. The depth-3
filler dependence is therefore a **staircase**, not a generic 2×2 coupling:

| rung | `P` used at | `S` used at |
|---|---|---|
| `Q_-1` | `E, E+1` | `E, E+2` |
| `Q_-2` | `E` (multiplier) | `E, E+1` |
| `Q_-3` | `E, E-1` | `E` (multiplier) |

This is what makes the node bookkeeping below finite and `d`-independent.

## 2. The membership silence table (PROVED, degree-free)

Two independent membership sources kill terms at integer nodes: the **filler** membership
(`P` vanishes on `{0,1}`, `S` on `{0,1,2}`) and the **cascade** membership `b_-1(0)=0`
(level `−1` carries `(E)_1`). Enumerating the staircase against them (`verify §S3`,
prediction + machine confirmation at `d=2,3,4,5`), exactly **six** `(rung,node)` pairs are
identically `0 = 0` — both filler part **and** residual `N_m`:

```text
   Q_-1(0),  Q_-2(0),  Q_-2(1),  Q_-3(0),  Q_-3(1),  Q_-3(2).
```

These are the tail's **silent point functionals**: they carry no information at all
(sharpening [`algebraic-covector.md`](algebraic-covector.md) §2, which recorded that
membership-window covectors annihilate *a block* — in fact at these six nodes the whole rung,
residual included, is vacuous).

## 3. The integer-node window: the fixed finite shape (PROVED, degree-free)

> **Theorem (integer-node window).** For every `k ≥ 3` and **every** degree `d`, the rows
>
> ```text
>    Q_-1(E=e),  e = 1..k-1 ;   Q_-2(E=e),  e = 2..k ;   Q_-3(E=e),  e = 3..k
> ```
>
> involve **only** the `2k−2` filler VALUES `P(2),…,P(k)` and `S(3),…,S(k+1)`; the augmented
> matrix `[M | N]` has the **`d`-independent** shape `(3k−4) × (2k−1)`, and its entries are
> the cascade datum evaluated at **fixed integer nodes**.

(`verify §S4`: the shape claim is a pure level-incidence enumeration verified for
`k=3..7`; the row identities are checked **exactly** against `Q_m(E=e)` on the real
parametrized cascade at `d=2,3` and, under `HEAVY`, `d=4`, for `k=3,4,5`.)

Consequences.

- **Excess.** `(3k−4) − (2k−2) = k−2` consistency conditions: the maximal minors of `[M|N]`
  (size `2k−1`) all vanish on the tail locus. These minors are a **fixed finite recipe**:
  fixed number of rows/columns, entries a fixed finite list of node values, at every degree.
- **The coupling obstruction is gone.** [`algebraic-covector.md`](algebraic-covector.md) §5
  located the gap in "a single algebraic node cannot annihilate a block because
  `gcd(a_2(E),a_2(E−3))=1`, so the covector must couple terms and rungs across the varying
  tops". The window **is** that coupling, in closed form: Cramer across the staircase at
  integer nodes. Nothing algebraic (no root of the datum) is ever named; the trace-form
  apparatus is not needed for the pure tail.
- The relaxation is **legitimate and strictly conservative**: the true fillers are
  polynomials, so every genuine tail solution gives a solution of the value-relaxed window.
  Hence *window-forcing ⇒ true forcing*. (The converse is exactly what fails; §5.)

## 4. On-shell: the cascade puts the b-tops on the same grid (PROVED, degree-free)

Attack angle 2 asked whether the `b`-tops `(b_1,b_0,b_-1)` become expressible through the
`a`-data on-shell. The answer is **partially yes, on the integer grid** (`verify §S5`,
generic data `d=2,3,4,6`). In `Q_4 = b_2(E+2)a_2 − a_2(E+2)b_2 + b_1(E+3)a_3 − a_3(E+1)b_1`
the top is the **fixed** `a_3 = E(E+2)(E+4)`, whose integer roots `{0,−2,−4}` meet every
residue class mod 3. At the five **anchor** nodes the `b_1`-dependence collapses to a
**single node value**:

| anchor `Q_4(x)=0` | `b_1`-dependence | consequence |
|---|---|---|
| `x=0` | `−15 b_1(0)` | `b_1(0) = (2/3) a_2(0)` |
| `x=−3` | `+3 b_1(0)` | `b_1(0) = (2/3) a_2(−3)` |
| `x=−2` | `+3 b_1(−2)` | `b_1(−2) = −(2/3) a_2(0)` |
| `x=−4` | `−3 b_1(−4)` | `b_1(−4)` fixed by `a_2(−4), a_2(−2)` |
| `x=−1` | `−3 b_1(2)` | `b_1(2) = (2/3)(2 a_2(−1) + a_2(1))` |

and the two `x=0,−3` anchors together give the **degree-free cascade node relation**

```text
   a_2(0) - a_2(-3) = Q_4(0)/10 + Q_4(-3)/2      =>   a_2(-3) = a_2(0)  on the cascade.
```

Off the anchors, `Q_4(E)=0` is a **3-step forward recursion**, e.g.
`15 b_1(4) = 48 b_1(1) − 18 a_2(1) + 4 a_2(3)`. The anchors constrain only the node values
`b_1 ∈ {0,−2,−4,2}`, so `b_1(1)` and `b_1(3)` are free seeds **of the value-level
recursion**.

> **Two audit corrections (2026-07-25).**
> 1. **Count.** There are **six** nodes at which `Q_4` sees exactly one `b_1` value —
>    `x ∈ {0,−1,−2,−3,−4,−5}` (`x=−5` sees `b_1(−2)`), not five. The sixth is redundant for
>    the constrained node-value set `{0,−2,−4,2}`, so nothing downstream changes, but the
>    enumeration (and the verifier's hardcoded anchor list) was incomplete.
> 2. **"Free seeds" is a VALUE-LEVEL statement only.** On the actual *polynomial* cascade
>    `Q_4=0` is a linear system for the `b_1` coefficients of **full rank** (rank `7` of `7`
>    at `d=3`, kernel `0`; now machine-checked in `§S5`), so `b_1` is **uniquely determined**
>    by `a_2`. The freedom lives in the value-level 3-step recursion, not in the pair.

**Verdict on angle 2: the collapse is real but PARTIAL.** On-shell the `b`-tops live on the
*same* fixed integer grid as the `a`-tops — which is why the window of §3 is finite and
degree-free — but the two blocks still do **not** become shifts of one polynomial, so the
AP telescoping of [`lambda-general-k.md`](lambda-general-k.md) Thm C does **not** reappear.
*(The earlier reason given here — "two free seeds per `b`-level survive" — is corrected
above: at the polynomial level there are none.)*

## 5. The verdict: the fixed-shape minors are not `W`-multiples (REFUTED, exact)

At `d=3` the cascade has 9 free coordinates, `R(1) = −(4/9) a_2(0)·am1_3`, so
`W = −(4/9) am1_3` and **`W | X ⟺ am1_3 | X`**. The window matrix `M` is `am1_3`-free
(`W` lives only in the residual `N`), so each maximal minor is a polynomial in `am1_3` whose
`W`-divisibility is decided by its residue at `am1_3=0`.

- **No window minor is a unit multiple of `W`.** At exact rational specializations of the
  other eight free coordinates, every nonzero maximal minor has a **nonzero** residue at
  `am1_3=0`, for `k=3` (1 minor), `k=4` (8 minors) and — under `HEAVY`, at two seeds —
  `k=5` (52 nonzero of 55) (`verify §S7`). So the strict prize form
  `μ·N = u·W` is **not** realized by any fixed-shape window minor at these widths.
- **The `k=3` window does not force `W` at all — explicit witness.** On the `d=3` cascade at

  ```text
   a2_0=-4, a2_2=-73426/34959, a2_3=0, a1_3=0, a0_0=2, a0_2=1, am1_2=3, am1_3=2, b0cK0=-4
  ```

  the `k=3` window system is **solvable**: `det[M|N]=0`, `rank M = rank[M|N] = 4` (full), and
  the verifier exhibits the explicit filler-value solution satisfying every window row —
  while `W = −8/9 ≠ 0` (`verify §S7`). This is a **witness point**, not a normal form, so it
  is a radical-correct refutation: `W ∉ sqrt(I_3)`.
- **The symbolic-node (`Q(E)`) combined-rung adjoint fails the same way.** Treating
  `[Q_-1;Q_-2;Q_-3]` at shifted arguments `E+j`, `j=0..J`, as **one** operator on the filler
  values gives a `3(J+1) × 2(J+3)` matrix, again `d`-independent in shape. At `J=3` it is
  square `12×12` and, at the tested exact specialization, has **rank 11** — a unique
  covector — but its residual pairing has a nonzero `am1_3`-free part, so it is **not** a
  `W`-multiple either (`verify §S8`; the rank statement is bounded to that specialization,
  not claimed as an identity).

**What this means — CORRECTED (audit, 2026-07-25).** The depth-3 tail *does* force `W`
(reproduced below), but that forcing is **not** captured by the ladder identities at any
**tested** finite integer window with the fillers relaxed to free values.

> **The earlier reading of this section — "the forcing consumes the *polynomiality* of the
> fillers" — is WRONG, and is contradicted by this file's own data.** The evaluation map
> (filler coefficients) `↦` (the `2k-2` window values `P(2..k), S(3..k+1)`) is **surjective
> exactly when `k ≤ d+2`**. Every width tested here (`k=3,4,5`) satisfies `k ≤ d+2 = 5` at
> `d=3`, so **the relaxation is vacuous at those widths**: the audit exhibited honest
> polynomials `P=(E)_2 V`, `S=(E)_3 C` with `deg V, deg C ≤ 3` realizing the `k=3` witness
> values exactly. No test in this memo probes polynomiality at all.
>
> **The correct diagnosis:** the chosen `3k-4` node equations are simply **too few**. At the
> very same point the *true* tail is inconsistent (audit: `Q_-1..Q_-3` has rank `M=8` vs
> `rank[M|b]=9`; `Q_-1..Q_-5` has Gröbner basis `[1]`, an exact UNSAT certificate) — the
> window omits the rest of the ladder content.
>
> **The corrected sharpening:** the value-relaxed window is too weak at every width
> `k ≤ d+2`, where the relaxation is vacuous; the **first width at which polynomiality could
> matter is `k = d+3`** (i.e. `k=6` at `d=3`), and **that regime is untested**. Note the
> width would then have to grow with the degree, which is itself the obstacle to a
> degree-free recipe.

## 6. Angle 1: what the consistency minors actually factor into (BOUNDED, `d=2`, exact)

The resultant conjecture (angle 1) predicted `Res_E(top(E), top(E−3))`-type cofactors. At
`d=2` (where `W ≡ 0` on the cascade, so forcing is vacuous but the window content is not)
the factorization is **exactly computable** (`verify §S6`):

- every maximal minor of the `k=3` window (1 minor) and of the `k=4` window (8 minors) has
  the same **single non-linear common factor** of total degree 4 (the two cores agree up to a
  scalar): the window's elimination content is **essentially principal**, stable in `k`;
- the linear common factors are **datum node values** at the fixed integer nodes
  (`a_2(−2)`, `a_2(−3) = a_2(0)`), i.e. **Sylvester-like in the node grid**, not generic
  resultants in `E`;
- and the `a_-1` **top coefficient** `am1_2` — the `d=2` look-alike of the `W`-slot — *is* a
  common factor of every `d=2` window minor. That look-alike is precisely what **fails** to
  persist: at `d=3` the true `W`-slot `am1_3` divides **no** window minor (§5). The naive
  extrapolation from `d=2` is therefore refuted in file.

## 7. Controls (per house rules)

- **W-kill reproduced by RABINOWITSCH, not by a normal form:**
  `am1_3 ∈ sqrt(cascade + Q_-1..Q_-5)` — sympy exact unit Rabinowitsch ideal over `QQ` and
  mod `65003` (`verify §S9`).
- **`a_2(0)` not forced:** an explicit cascade+depth-3-tail point with `a_2(0)=4 ≠ 0`,
  `R(1)=0` (the kill is the **factor** `W`).
- **Non-vacuity:** that same point carries an **explicit filler solution** satisfying every
  depth-3 tail equation, so `cascade + tail` is nonempty.

## 8. Evidence ledger — proved / bounded / refuted / open

**Proved (degree-free, symbolic proof object in file):**
- The factorization `R(1)=a_2(0)·W` (proved on abstract symbols, no degree) (`§S1`).

> **AUDIT RE-TIERING (2026-07-25).** The remaining items below were listed as "Proved
> (arbitrary degree)" but rest on **fixed-degree instance checks**, which the house rule
> forbids as a basis for that tier. Each is stated with its in-file machine scope and the
> audit independently re-proved the staircase and the window incidence degree-free with
> abstract `sympy.Function` symbols — so the mathematics stands and the fix is small; the
> *certification artifact* is what is bounded. Additional scope notes: the "exactly six"
> silence table is **one-sided in file** (the six predicted pairs are confirmed to vanish;
> that no *other* pair vanishes is not tested in file — the audit verified the converse
> externally at generic `d=3`), and the window theorem's `every k, every d` quantifier is
> enumerated only to `k≤7` with the row identities tied to the algebra at `d≤4`, `k≤5`.

**Paper/level-incidence proofs with bounded machine scope (audit-retiered):**
- Engine `Q_m=[D,X]_m` (`m∈[−6,6]`), `Q_0=(T−1)G` (in-file machine scope: generic degree 2);
  slope gate; both-ends Lemma P at `E=1`; filler-independence of `R(1)`; the factorization
  `R(1)=a_2(0)·W` (symbolic, degree-free) (`§S0,§S1`).
- **Staircase collapse:** the two-term block formulas plus the degeneration of the `Q_-2`
  a-block and the `Q_-3` b-block to pure multipliers — in-file machine scope `d=2,3,4,5`,
  proof degree-free by level incidence (`§S2`).
- **Membership silence table:** exactly the six pairs `Q_-1(0), Q_-2(0), Q_-2(1), Q_-3(0),
  Q_-3(1), Q_-3(2)` are identically `0=0`; level-incidence prediction + machine confirmation
  at `d=2,3,4,5` (`§S3`).
- **Integer-node window theorem:** for `k=3..7` the level-incidence shape argument
  (`(3k−4)` rows, `2k−2` live filler values, `[M|N]` of shape `(3k−4)×(2k−1)`, entries at
  fixed integer nodes), verified as **exact row identities** on the cascade at `d=2,3`
  (`+4` HEAVY) for `k=3,4,5` (`§S4`).
- **On-shell node recursion:** the `Q_4` anchors each see exactly one `b_1` node value
  (**six** such nodes, audit-corrected from five);
  `b_1(0)=(2/3)a_2(0)`, `b_1(0)=(2/3)a_2(−3)`, the node relation `a_2(0)−a_2(−3) =
  Q_4(0)/10 + Q_4(−3)/2` (hence `a_2(−3)=a_2(0)` on the cascade),
  `b_1(2)=(2/3)(2a_2(−1)+a_2(1))`, `15 b_1(4)=48 b_1(1)−18a_2(1)+4a_2(3)`, and the free
  seeds `b_1(1), b_1(3)` — machine scope `d=2,3,4,6`, argument degree-free (`§S5`).

**Bounded-finite (exact scope stated):**
- `d=2`: the window minors' common-factor structure (one degree-4 core for `k=3,4`, stable up
  to scalar; linear common factors are datum node values; `am1_2` a common factor) (`§S6`).
- `d=3`: `R(1)=−(4/9)a_2(0)·am1_3`; `M` is `am1_3`-free; the minor sweeps at one exact
  specialization (`k=3,4`) and, under `HEAVY`, two specializations plus `k=5` (`§S7`).
- `d=3`: the combined-rung `Q(E)` adjoint at `J=3` is `12×12` of rank 11 at the tested exact
  specialization, unique covector, pairing not a `W`-multiple (`§S8`).
- Controls: the exact-`QQ` + mod-`p` Rabinowitsch W-kill on `Q_-1..Q_-5`; the `a_2(0)`-free
  witness; the non-vacuity filler solution (`§S9`).

**Refuted (exact, in file):**
- *"A fixed-shape window minor is a unit multiple of `W`."* **False** at `d=3` for
  `k=3,4` (default) and `k=5` (HEAVY): every nonzero maximal minor has a nonzero residue at
  `am1_3=0`. **Scope (audit): at `k=5`, 3 of the 55 maximal minors vanish at the tested
  specialization and are therefore NOT tested — the statement covers the 52 that are
  nonzero there (a minor vanishing at a specialization may still be symbolically nonzero).**
- *"The fixed integer-node window forces `W`."* **False at `k=3`**, by an **explicit exact
  rational witness** on the `d=3` cascade with the window solvable and `W=−8/9 ≠ 0`
  (radical-correct: a point, not a normal form).
- *The `d=2` extrapolation* — that the `a_-1` top coefficient (the `W`-slot look-alike),
  which divides every `d=2` window minor, keeps dividing at higher `d`. **False at `d=3`.**
- The *coupling* framing of [`algebraic-covector.md`](algebraic-covector.md) §5 — "the
  annihilator must couple terms and rungs across the varying tops, and no fixed telescoping
  is available" — is **sharpened, not sustained as an obstruction**: the coupling is
  explicitly available as the fixed-shape integer-node Cramer minors (§3). What blocks the
  prize is not the coupling but the **content** of those minors.

**Open / not claimed:**
- Whether **any** `k` has `W ∈ sqrt(I_k)` for the value-relaxed integer-node window. Refuted
  at `k=3`; at `k=4,5` only the strict "single minor `= u·W`" form is refuted, and the ideal
  question is untested (the msolve/Gröbner Rabinowitsch on the window ideal did not terminate
  in budget at `d=3`, `k=3,4`; reported as such, not claimed).
- **The sharpened residual gap:** since the value-relaxed window is provably too weak at
  `k=3`, any degree-free recipe must use the **finite-difference relations** among the filler
  node values (i.e. `deg P, deg S ≤ d`) — but those relations are *degree-dependent*. A
  degree-free recipe must therefore either (i) find a `k` at which the value-relaxed window
  already forces `W`, or (ii) use a degree-uniform surrogate for polynomiality (e.g. the
  `Δ^{n}`-annihilation of a *fixed* number of leading node differences). Neither is obtained.
- Everything the parent memos leave open: the arbitrary-degree slope forcing
  ([`slope-forcing-degree-free.md`](slope-forcing-degree-free.md) §6,
  [`residual-identity.md`](residual-identity.md) §6,
  [`algebraic-covector.md`](algebraic-covector.md) §6); no Weyl pair; all of Band 3, DC1, JC2.

## 9. Verification

```sh
uv run --with sympy python research/dc1-program/verify_coupling_covector.py
# HEAVY (adds the d=4 window row identities, the k=5 minor sweep and a second seed):
HEAVY=1 uv run --with sympy python research/dc1-program/verify_coupling_covector.py
```

`S0` engine; `S1` slope gate + Lemma P + factorization; `S2` staircase collapse; `S3`
membership silence table; `S4` the integer-node window theorem (shape + exact row identities);
`S5` on-shell node recursion; `S6` `d=2` minor factorization; `S7` the `d=3` verdict + the
explicit `k=3` witness; `S8` combined-rung `Q(E)` adjoint; `S9` controls. Default run
**~46 s, 166 checks, one optional leg skipped**; `HEAVY` **~115 s, 175 checks, no skips**
(runtimes environment-dependent).
The verifier never folds a skipped optional leg into an unqualified all-passed banner.
Sampling and specialization legs do not certify the broader conjectures.
