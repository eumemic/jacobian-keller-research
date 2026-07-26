# The broken separation classes: diff-1 and diff-2 non-cube-separated shifted cubes

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED — PARTIAL RESULT WITH NAMED GAPS**

This memo works the **last shifted-power walls** of band-3 Gap 2: the two
non-cube-separated shifted-cube top classes that the 2-separation theorem
([`shifted-power-residuals.md`](shifted-power-residuals.md) §2) does **not**
reach. The closed shifted-cube class is now `2`-separated `h`
(`gcd(h,h^{[1]}) = gcd(h,h^{[2]}) = 1`), so the genuinely broken classes are
exactly the two where a root-gap of `1` or `2` breaks a coprimality:

```
diff-1:  h = (E-r)(E-r-1),   gcd(h,h^{[1]}) = (E-r) ≠ 1,   r a free parameter,
diff-2:  h = (E-r)(E-r-2),   gcd(h,h^{[2]}) = (E-r) ≠ 1,   r a free parameter,
```

together with the multiple-root top `h = (E-r)²`. The one-line summary:

> **The multiple-root top is not a broken class.** `h = (E-r)²` has root
> multiset `{r,r}`, gap `0 ∉ {1,2,3}`, so `gcd(h,h^{[j]}) = 1` for `j=1,2,3`: it
> is **cube-separated**, and the arbitrary-degree shifted-cube theorem already
> forces `h` constant. Only **diff-1** and **diff-2** genuinely break.
>
> **The degraded `Q₄` forcing is sharp and asymmetric (arbitrary degree,
> symbolic `r`).** For **diff-1**, `Q₄ = 0` forces `(E-r) | a₂` **plus a coupled
> correction** `a₂(r-1) + a₂(r+1) = 0`; the clean `h h^{[1]} | a₂` fails and even
> the double-root condition `a₂'(r)=0` is lost. For **diff-2**, `Q₄ = 0` forces a
> **clean proper-factor divisibility** `(E-r)(E-r-1) | a₂` (gcd node `r` plus the
> midpoint `r+1`, **no coupling**); the clean `h h^{[1]} | a₂` fails (the
> vanishings at `r-1, r+2` are lost).
>
> **Audit update: the classes have opposite outcomes.** The general `Q₃` identity has now
> been used radical-correctly. For diff-1, `Q₄∧Q₃` forces the missing value with exponent
> two and restores `h h^{[1]}|a₂`; those nonzero branches combine with Nonpositive-D Exclusion Theorem of
> [`shifted-cube-completion.md`](shifted-cube-completion.md), which closes
> `κ=0,b₁=0`, to close diff-1 at arbitrary degree with no side condition. For diff-2,
> Nonpositive-D Exclusion Theorem closes only `κ=0,b₁=0`; an exact symbolic positive-cascade family shows
> the surviving `κ≠0` branch remains open at arbitrary degree.
> The older fixed-cap normal-form discussion below is superseded by these radical-correct
> results; see [`band3-sectors.md`](band3-sectors.md).
>
> **Multiplicity-extended adjoint criterion.** The `Q₃` divisibility `a₃ | Src₃`
> becomes the **jet criterion** `a₃ | P ⇔ P^{(j)}(ν)=0, j<m_ν`; at the double
> nodes this adds a genuine **derivative-node equation** `ev'_ν(Src₃)=0` beyond
> the ordinary `ev_ν`. The Lemma-P moment slope `G(1) = Σ_i(a_i(0)b_{-i}(i) −
> a_{-i}(i)b_i(0))` is re-derived for the degenerate tops; its finite-cap
> cokernel is the bounded emptiness below.
>
> **Bounded emptiness, both classes.** Committed exact SymPy: with `r,κ` treated
> as coefficient parameters (displayed domain `ZZ[r,κ]`), the full positive cascade
> `+ Q₀=1 +` membership is the **unit ideal on the cap-`d=2` generic fiber over
> `QQ(r,κ)`**. Separate exact checks cover a grid of specific integer parameter fibers
> — an upgrade over the prior fixed-`r=0` evidence. This is not a uniform certificate
> at every specialization of `r,κ`.
> HEAVY (`msolve`, exact `QQ`): historically, the unit ideal on a sampled `(r,κ)` grid at
> cap `d=2,3,4`. This bounded evidence is no longer the status statement: diff-1 is closed
> in the branches above, while diff-2 remains open at arbitrary degree after restoration
> was refuted.

Exact certificate:
[`verify_broken_separation.py`](verify_broken_separation.py) — default exact SymPy;
the optional `msolve` grid is behind `HEAVY=1`. Runtime and check counts are
environment-dependent. The final summary distinguishes all checks passed from all
executed checks passed with skips. Every load-bearing upstream fact (the crossed-product
ladder engine `Q_m = [D,X]_m`, `Q₀ = (T−1)G`, the `Q₅` wall, the general `Q₄`/`Q₃`
identities) is **re-derived in file**, not cited.

Conventions frozen from the corpus: `A₁[x^{-1}] = ⊕_k x^k C[E]`, `E = x∂`,
`(x^a f)(x^b g) = x^{a+b} f(E+b) g(E)`, `f^{[n]}(E) = f(E+n)`, `T f = f^{[1]}`,
`S_n = 1+T+⋯+T^{n-1}`, `Q_m = ∑_{k+l=m}[b_l^{[k]}a_k − a_k^{[l]}b_l]`,
`[D,X]=1 ⇔ Q_m=δ_{m0}`, membership `(E)_j | a_{-j},b_{-j}`, gauge `b_k=0`,
`Q₀=(T−1)G`. Sector: gauge `b₃=0`, `a₃ = h(E)h(E+1)h(E+2)`, `b₂ = κ·h(E)h(E+1)`
(the `Q₅` wall).

---

## 1. The engine, the wall, and the general rung identities (re-derived)

Verifier `§0` re-derives `Q_m = [D,X]_m` for `m ∈ [−6,6]` and `Q₀ = (T−1)G` with
the closed-form staggered potential on generic band-3 data. Verifier `§1`
re-derives the `Q₅` wall solution `b₂ = κ h h^{[1]}` (any `h`) and the two rung
identities used throughout, **with no shape assumption** on the sub-leading data:

```
Q₄ = h^{[1]}h^{[2]}( b₁^{[3]} h − h^{[3]} b₁ ) + κ( h^{[2]}h^{[3]} a₂ − h h^{[1]} a₂^{[2]} ),
Q₃ = a₃( b₀^{[3]} − b₀ ) + ( b₁^{[2]} a₂ − a₂^{[1]} b₁ ) + κ h^{[1]}( h^{[2]} a₁ − h a₁^{[2]} ).
```

In the clean (2-separated) case `a₂ = h h^{[1]} g`, `b₁ = h β` make the `Q₃`
middle term `b₁^{[2]}a₂ − a₂^{[1]}b₁ = a₃(β^{[2]}g − g^{[1]}β)` divisible by `a₃`;
it is exactly this cancellation that the broken classes **lose**.

---

## 2. Task 1 — the exact degraded `Q₄` forcing per class (symbolic `r`)

`Q₄ = 0` is **linear** in the joint coefficient vector `(a₂, b₁)` (the wall term
`b₂ a₂` is linear in `a₂`), so its full solution space is reliable — no
component dropping. Solving it and reading the divisibility of `a₂` gives, exactly
and symbolically in `r`:

### 2.1 diff-1 — a proper factor **plus** a coupled correction

Writing `u = E − r`, `h h^{[1]}` carries `(u+1)²` (the shared factor `h,h^{[1]}`),
so after cancelling one `(u+1)` the modular conditions `u(u+1)(u+2) |` (the
`a₂`-source) are:

```
u=0, u=-2 (nodes r, r-2):   a₂(r) = 0,                       ⟹  (E-r) | a₂,
u=-1     (node r-1):        a₂(r-1) + a₂(r+1) = 0.           ⟹  a COUPLED correction
```

The middle node's coefficients do **not** vanish, coupling `ev_{r-1}` and
`ev_{r+1}`. So:

> **diff-1 degraded `Q₄` forcing:** `(E-r) | a₂` **and** `a₂(r-1) + a₂(r+1) = 0`.
> The clean `h | a₂`, `h h^{[1]} | a₂`, and the double-root condition `a₂'(r)=0`
> **all fail** — `a₂'(r)` is free.

The clean case would demand `a₂` vanish at `{r-1, r, r, r+1}` (four conditions);
the degraded forcing keeps only `a₂(r)=0` and replaces the two separate
vanishings `a₂(r-1)=a₂(r+1)=0` by their **sum** `=0`. Verifier `§2`.

### 2.2 diff-2 — a clean proper-factor divisibility (no coupling)

Here `h,h^{[1]}` are coprime but `h,h^{[2]}` share `(E-r)`. The same reduction
gives `(u-1)(u+2) |` (the `a₂`-source), and at **each** node one of the two
coefficients vanishes, so the two evaluations **decouple**:

```
u=-2 (node r):   a₂(r) = 0,        u=1 (node r+1):   a₂(r+1) = 0.
```

> **diff-2 degraded `Q₄` forcing:** the **clean** `(E-r)(E-r-1) | a₂` — the gcd
> node `r` plus the *midpoint* `r+1` of the two `h`-roots `r, r+2`. **No coupling
> term.** The clean `h h^{[1]} | a₂` fails: the vanishings at `r-1, r+2` are lost.

Note the forced divisor `(E-r)(E-r-1)` is **not** `h` and shares only `(E-r)` with
`h` — the extra root `r+1` is off `h`. So diff-2 is *cleaner* than diff-1 at
`Q₄` (no coupling), but still genuinely weaker than the 2-separated case.
Verifier `§2`.

### 2.3 the multiple-root top is cube-separated (closed, not broken)

`h = (E-r)²` has `gcd(h,h^{[j]}) = gcd((E-r)²,(E-r-j)²) = 1` for `j=1,2,3` (a
common root needs `j=0`). So it is **cube-separated**: `Q₄` forces the clean
`h h^{[1]} | a₂`, the affine kill runs (`h^{[-1]}` degree `2` cannot divide `E`),
and `h` is forced constant at **arbitrary degree** by the shifted-cube theorem.
It is **covered**, not a separate broken class. Verifier `§2`.

---

## 3. Task 2 — historical fixed-cap probes, superseded by radical-correct restoration tests

### 3.1 double nodes and the general `Q₃`

For the broken classes `a₃ = h h^{[1]}h^{[2]}` has a **repeated** root multiset:

```
diff-1:  {r-2, r-1, r-1, r, r, r+1}   — DOUBLE nodes at r-1 and r,
diff-2:  {r-2, r-1, r, r, r+1, r+2}   — DOUBLE node at r.
```

(The mult-root `(E-r)²` gives a perfect square `a₃`, all nodes double, but stays
cube-separated.) These double nodes are the source of the derivative-node
criterion of §4. Verifier `§3`.

### 3.2 a generic `Q₄` solution does **not** extend to `Q₃`

Substituting the generic `Q₄` solution `(a₂, b₁)` (free parameters `θ`) into
`Q₃`, the verifier's **fixed coefficient ansatz** gives a linear system in
`(a₁,b₀)` that is inconsistent for generic symbolic parameters. Thus `Q₃` cuts the
finite-dimensional `Q₄` solution space at the implemented caps. This is bounded
symbolic-parameter evidence, not a theorem over fillers of arbitrary degree.
Verifier `§3` (both classes; fixed caps, symbolic `r,κ`).

### 3.3 ideal-membership probe for clean divisibility (historical and superseded)

The original fixed-cap probe asked whether `Q₃` cuts the `a₂` freedom back onto the
clean `h h^{[1]} | a₂`. At `(r,κ)∈{(2,1),(0,2)}`, Gröbner normal forms in the
coefficient ideal gave:

```
gcd node    a₂(r)    lies in the tested (Q₄,Q₃) ideal,
clean-extra a₂(r-1) has nonzero normal form modulo that ideal.
```

The second statement alone was inconclusive because geometric forcing is controlled
by `sqrt(I)`, not `I`: for example, `x∉(x²)` although `x` vanishes on every point of
`V(x²)`. The later radical-correct calculation resolves the two classes differently.
For diff-1, the missing evaluation is forced with exponent two, restoring
`h h^{[1]}|a₂`; this closes the nonzero branches, while Nonpositive-D Exclusion Theorem of
[`shifted-cube-completion.md`](shifted-cube-completion.md) closes `κ=0,b₁=0` and
therefore completes diff-1 at arbitrary degree without a side condition. For diff-2,
restoration is refuted on the surviving `κ≠0` branch; Nonpositive-D Exclusion Theorem closes only
`κ=0,b₁=0`, not the full class. See [`band3-sectors.md`](band3-sectors.md).

> **Conclusion (Task 2, updated).** The fixed-cap normal forms were useful diagnostics
> but are no longer the status result. Radical-correct forcing closes the stated diff-1
> branches; an exact `(Q₄,Q₃)`-locus witness refutes diff-2 restoration without
> constructing a Weyl pair or resolving the full diff-2 tail.

---

## 4. Task 3 — the multiplicity-extended moment/adjoint criterion

### 4.1 the jet divisibility criterion

For `a₃ = ∏_ν (E-ν)^{m_ν}`, `a₃ | P ⇔ P^{(j)}(ν) = 0` for all `j < m_ν` at every
node `ν`. This extends the moving-sum/adjoint node functional `ev_ν` to the
**jet** `(ev_ν, ev'_ν, …, ev^{(m_ν-1)}_ν)`. Verifier `§4` checks the jet
criterion agrees with polynomial-remainder divisibility, including at a double
node.

### 4.2 the derivative-node equations `ev'_ν`

Applied to `Q₃ = a₃(b₀^{[3]}−b₀) + Src₃`, `Src₃ = (b₁^{[2]}a₂ − a₂^{[1]}b₁) +
κ h^{[1]}(h^{[2]}a₁ − h a₁^{[2]})`, the divisibility `a₃ | Src₃` becomes: at each
simple node the ordinary `ev_ν(Src₃) = 0`, and at each **double** node the extra

```
ev'_ν(Src₃) = 0        (a genuine linear functional, absent at simple nodes),
```

evaluated with `d/dE`. For diff-1 the node `r` is double and `ev'_r(Src₃) = 0` is
a genuine extra equation distinct from `ev_r(Src₃) = 0`; the node `r-2` is simple
and gives only `ev_{r-2}(Src₃) = 0`. This is the exact multiplicity-extension of
the adjoint criterion the moving-sum program needs at degenerate tops. Verifier
`§4`.

### 4.3 the moment slope for the degenerate tops

Membership `(E)_j | a_{-j}, b_{-j}` gives `G(0) = 0`; `Q₀ = 1 ⇔ G = E` then
requires the slope `G(1) = 1`, and the Lemma-P formula
([`moment-unit-general-k.md`](moment-unit-general-k.md) §1) holds **verbatim** for
the degenerate tops:

```
G(1) = ∑_{i=1}^{3} ( a_i(0) b_{-i}(i) − a_{-i}(i) b_i(0) ),   G(0) = 0.
```

Verifier `§4` re-derives this at symbolic `r` for both classes. The finite-cap
cokernel of the combined slope/adjoint system is exactly the bounded emptiness of
§5 (the `[E − R] ∈ Im Φ` question of
[`../band3/quantum-exotic-cokernel.md`](../band3/quantum-exotic-cokernel.md), one
top-level up). No arbitrary-degree kill follows from the slope alone here — the
degenerate cokernel does not collapse symbolically at the caps attempted.

---

## 5. Task 4 — bounded emptiness certificates per class

The full positive cascade `Q₄=Q₃=Q₂=Q₁=0` with `Q₀=1`, gauge `b₃=0`, wall
`b₂ = κ h h^{[1]}`, and genuine membership. Emptiness = the coefficient ideal is
the **unit ideal**.

**Committed (default run, exact SymPy):**
- cap `d=2` **unit ideal on the generic fiber over `QQ(r,κ)`**, both diff-1 and
  diff-2. SymPy displays coefficient domain `ZZ[r,κ]` because `r,κ` are not Gröbner
  variables; the result is generic-fiber evidence, not uniform emptiness under every
  specialization. This upgrades the prior fixed-`r=0`, all-`κ` evidence of
  [`shifted-power-residuals.md`](shifted-power-residuals.md) §2.4 (Verifier `§5a`);
- cap `d=2` unit ideal at specific `(r,κ) ∈ {(0,1),(1,1),(-1,2),(2,1)}`, both
  classes (Verifier `§5b`).

**HEAVY (`msolve`, exact `QQ`, SKIP if absent):**
- cap `d=2, 3, 4` unit ideal on the sampled grid
  `(r,κ) ∈ {(0,1),(0,2),(1,1),(-1,1),(2,1),(½,1)}`, both classes; each `msolve`
  call is instant (0-dimensional/unit, `≤ 0.1 s`) (Verifier `§5d`).

This bounded computation is **corroboration**, not the source of the later
arbitrary-degree conclusions. The radical-correct argument closes the nonzero
diff-1 branches, and Nonpositive-D Exclusion Theorem of
[`shifted-cube-completion.md`](shifted-cube-completion.md) closes `κ=0,b₁=0`, so
diff-1 is closed without side conditions. For diff-2, Nonpositive-D Exclusion Theorem closes only that
zero branch; the surviving `κ≠0` branch remains open. The cap-`d=4` sampled
emptiness and cap-`d=2` **generic-fiber** emptiness over `QQ(r,κ)` remain valid
bounded evidence.

---

## 6. Bookkeeping — what Gap 2's band-3 instance still needs

**Newly established (this memo):**
1. The multiple-root top `h = (E-r)²` is **cube-separated ⇒ closed arbitrary
   degree** — it is *not* one of the broken classes (a corpus clarification: the
   broken set is exactly diff-1 and diff-2, each a one-parameter family in `r`).
2. The **exact degraded `Q₄` forcing** per class, symbolic in `r`:
   diff-1 `(E-r)|a₂` + coupled `a₂(r-1)+a₂(r+1)=0`; diff-2 clean `(E-r)(E-r-1)|a₂`.
3. The original cascade probes add fixed-cap constraints but were inconclusive about
   radicals. The subsequent radical-correct update resolves restoration: diff-1 restores
   `h h^{[1]}|a₂` and closes in the stated nonzero branches; diff-2 has an exact
   `(Q₄,Q₃)`-locus witness refuting restoration. The latter is not a full-tail point and
   generally has `Q₂≠0`.
4. The **multiplicity-extended adjoint criterion** (jet / `ev'_ν` derivative
   nodes) and the Lemma-P moment slope for the degenerate tops.
5. Bounded emptiness upgraded: cap `d=2` **generic `(r,κ)`** (committed) and cap
   `d≤4` sampled grid (HEAVY), both classes. **Non-vacuity is checked exactly in-file:
   `verify_band3_sectors.py` §7 drops `Q₀=1` at cap `d=2` and obtains a proper
   (non-unit) ideal for both diff-1 and diff-2. Thus each sector is nonempty without
   the moment unit, so the emptiness certificates are meaningful, not vacuous.**

**Still open (precisely delimited):**
- **diff-1:** nothing; it is closed at arbitrary degree with no `κ`/`b₁` side
  condition by combining the restored-divisibility branches with Nonpositive-D Exclusion Theorem of
  [`shifted-cube-completion.md`](shifted-cube-completion.md).
- **diff-2 at arbitrary degree:** only the surviving `κ≠0` branch remains open;
  the `κ=0,b₁=0` arm is closed by Nonpositive-D Exclusion Theorem. A different tail obstruction is
  still needed, so no full diff-2 closure is claimed.
- Everything the 2-separation memo left open is unchanged: **RESIDUAL 3**
  (`κ₂ ≠ 0` constant-`h`, the A\*-band3 negative tail); imbalanced coprime walls;
  general-`k` negative tail; **W2**.

After the radical-correct update, the band-3 shifted-cube Gap-2 ledger reads:

| class | status |
|---|---|
| cube-separated / 2-separated `h` (incl. diff-3, mult-root `(E-r)²`) | **closed, arbitrary degree**; old `c=0` hole repaired by [`shifted-cube-completion.md`](shifted-cube-completion.md), Nonpositive-D Exclusion Theorem |
| **diff-1** `(E-r)(E-r-1)` | **closed, arbitrary degree, no side condition**; restored divisibility handles the nonzero branches and Nonpositive-D Exclusion Theorem handles `κ=0,b₁=0` |
| **diff-2** `(E-r)(E-r-2)` | `κ=0,b₁=0` **closed** by Nonpositive-D Exclusion Theorem; surviving `κ≠0` branch **open** at arbitrary degree |

**No Weyl pair and no counterexample is constructed; DC1/JC2 untouched.**

---

## 7. Honest ledger — proved / bounded / refuted / open

**Proved (exact algebra, machine-checked identities; arbitrary degree where
stated):**
- Engine `Q_m=[D,X]_m`, `Q₀=(T−1)G`; the `Q₅` wall; the general `Q₄`/`Q₃`
  identities with no shape assumption (`§0,§1`).
- **Task 1 (arbitrary degree, symbolic `r`):** the exact degraded `Q₄` forcing —
  diff-1 `(E-r)|a₂` + coupled correction `a₂(r-1)+a₂(r+1)=0` (clean `h h^{[1]}|a₂`
  and `a₂'(r)=0` fail); diff-2 clean `(E-r)(E-r-1)|a₂` (no coupling; clean
  `h h^{[1]}|a₂` fails); mult-root `(E-r)²` cube-separated ⇒ closed (`§2`).
- **Task 2:** the double-node multiset structure and the historical fixed-ansatz
  generic-`Q₄`-does-not-extend-to-`Q₃` inconsistency (`§3`). The old normal forms are
  retained only as bounded diagnostics. The radical-correct update restores the clean
  divisor for diff-1 in the stated nonzero branches and gives an exact diff-2
  `(Q₄,Q₃)`-locus witness refuting restoration.
- **Task 3:** the jet/derivative-node divisibility criterion; the `ev'_ν` double-
  node equation; the Lemma-P moment slope for the degenerate tops (`§4`).

**Bounded / finite evidence (exact scope):**
- diff-1, diff-2 sector emptiness: cap `d=2` committed — the **generic fiber** over
  `QQ(r,κ)` and a specific integer-parameter grid (exact SymPy), with no uniform
  claim over every specialization; cap `d=2,3,4` HEAVY —
  a sampled `(r,κ)` grid (`msolve`, exact `QQ`). This is corroboration, **not** an
  arbitrary-degree exclusion.

**Refuted (machine-checked) — a corpus clarification:**
- The multiple-root top `h = (E-r)²`, which a naive "roots differing by
  `1,2,3`" reading might file as broken, is in fact **cube-separated** (root gap
  `0`), hence **closed at arbitrary degree**. It is not a broken class.

**Open / NOT claimed:**
1. The diff-2 `κ ≠ 0` surviving branch at arbitrary degree. The `κ = 0, b₁ = 0`
   arm is closed by Nonpositive-D Exclusion Theorem of
   [`shifted-cube-completion.md`](shifted-cube-completion.md); the direct
   residual/congruence route remains silent at `κ = 0, b₁ = 0`, but Nonpositive-D Exclusion Theorem
   closes it independently. No full diff-2 closure is claimed.
2. **A\*-band3 general `κ₂ ≠ 0`** at arbitrary degree remains open. Nonpositive-D Exclusion Theorem
   closes only the corner `κ₂ = 0, b₁ = 0`; the existing `κ₂ = 0` tame witness
   has positive band-one `D` (`b₁ ≠ 0`), so no contradiction with Nonpositive-D Exclusion Theorem.
3. Everything inherited open from `shifted-power-residuals.md` §5:
   composite tame-word escape, imbalanced coprime walls, general-`k` tail with
   `band D > 0`, **W2**, and radical forcing at coupling widths `k = 4, 5`.
No Weyl pair, no counterexample; DC1/JC2 untouched.

---

## 8. Verification

```sh
uv run --with sympy python research/dc1-program/verify_broken_separation.py
HEAVY=1 uv run --with sympy python research/dc1-program/verify_broken_separation.py
```

Exact SymPy checks (with domains stated by each section): `§0` engine; `§1` wall +
general `Q₄`/`Q₃` identities; `§2` Task 1 (degraded `Q₄` forcing per class,
symbolic `r`; mult-root closure); `§3` Task 2 (double nodes and the historical
fixed-cap probes, now superseded as status evidence by the radical-correct checks in
`band3-sectors.md`); `§4` Task 3 (jet criterion, `ev'_ν` double-node equations,
moment slope); `§5` Task 4 (cap-`d=2`
generic-fiber + tested specific `(r,κ)` emptiness; HEAVY `msolve` grid at `d=2,3,4`); `§6`
bookkeeping. Runtime is environment-dependent. The final summary distinguishes a run
in which all checks passed from one in which all executed checks passed but optional
checks were skipped.
