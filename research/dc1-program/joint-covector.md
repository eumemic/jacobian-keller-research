# The Joint Filler Covector: explicit covector, necklace support, and the d=3–4 uniform kill

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED — BAND-SCOPED — MIXED PROVED/BOUNDED**

This memo attacks architecture step 1, the **Joint Filler Covector Lemma** (the W2
lemma of [`w2-joint-theorem.md`](w2-joint-theorem.md) §4). It does four things:

1. **Extracts the generic covector explicitly.** For the stacked linear filler
   operator `L = [Phi ; M_-1]` on the W2 slope-`1` cascade locus, the Fredholm gap
   is `1`, so the annihilating covector `lam` (with `lam·L=0`, `lam·rhs != 0`) is
   generically unique up to scale. It is produced two ways in file — the left
   kernel of `L`, and an explicit **pivot/residual formula**
   `lam.rhs = L_j L_I^{-1} rhs_I - rhs_j`, a rational function of the datum.
2. **Confirms the union-of-necklaces closed form** (bounded, at the datum): `lam`'s
   block-0 part `mu0` **must leave** the point triple `Ann(Im Phi)={ev_-1,ev_0,ev_1}`
   (which is silent at slope one) and lives on the **fixed `(a_3,b_2)` root
   necklace** `roots(a_3) ∪ roots(b_2) ∪ roots(q_K) ∪ roots(q_H) = {-4,..,3}`; its
   block-1 part `mu1` is on the datum `(a_2,b_1)` necklace.
3. **Establishes the payload as a UNIFORM d=3–4 kill.** `cascade + Q_0=1 + Q_-1`
   (i.e. `[Phi ; M_-1]`) is the **unit ideal at raw cap `d=3` and `d=4`, on both tail
   branches**, over `QQ` (`msolve` with `^`) with an independent `sympy`
   point-cross-check. So the joint obstruction is present on the **whole** slope-`1`
   cascade locus at those caps — not merely generically.
4. **Corrects a committed claim.** The sibling assertion
   ([`w2-joint-theorem.md`](w2-joint-theorem.md) §3, §S8) that *"a slope-`1`
   sub-locus survives `Q_-1`; the uniform kill needs the full tail `Q_-1..Q_-5`"* is
   a `msolve` `**`-parser artifact (flagged by
   [`hatch-census.md`](hatch-census.md) §5). With `^` the system is the **unit
   ideal**; the surviving sub-locus is **empty at `d=3,4`**, and `Q_-2..Q_-5` are
   not needed there.

W2 datum, gauge `b_3=0`, quantum band-3 conventions
(`Q_m = sum_(k+l=m)[b_l^[k] a_k - a_k^[l] b_l]`, `f^[n](E)=f(E+n)`,
membership `(E)_j=E(E-1)...(E-j+1) | a_-j, b_-j`):

```text
a_3 = E(E+2)(E+4)   (roots {0,-2,-4}),   b_2 = E(E+3)   (roots {0,-3}),   b_3 = 0.
```

Branch B fillers `a_-2=(E)_2 V`, `b_-3=(E)_3 C` (`a_-3=0`);
`Phi(C,V)=K_3[(E)_3 C]-H_2[(E)_2 V]`, `Im Phi = D·F[E]`, `D=E(E-1)(E+1)`.
`M_-1` is the filler-linear part of `Q_-1` — a two-block map with `(a_2,b_1)` in
place of `(a_3,b_2)`.

Exact certificate: [`verify_joint_covector.py`](verify_joint_covector.py) (ends
`ALL JOINT COVECTOR CHECKS PASSED`; default run ~8 s, `HEAVY=1` adds the `d=4`
`QQ` legs, ~20 s). Every load-bearing upstream fact — the crossed-product engine,
`Im Phi = D·F[E]`, the slope gate — is **re-derived in file**, not cited.

## 0. The obstruction, as a covector

After the positive cascade `Q_4=Q_3=Q_2=Q_1=0` solves `b_1,b_0,b_-1,b_-2`, the
moment potential splits `G = R + Phi(C,V)` with `R` the non-filler residual, and
`Q_0=(T-1)G` with membership `G(0)=0` gives

```text
Q_0 = 1   <=>   G = E   <=>   Phi(C,V) = E - R      (block 0),
```

while `Q_-1 = 0` reads `M_-1(C,V) = -N_-1` (block 1), where `N_-1` is the
non-filler part. Stacking, `L : (C,V) -> (Phi(C,V), M_-1(C,V))` with
`rhs = (E-R, -N_-1)`. A **joint filler covector** is `lam=(mu0,mu1)` with

```text
mu0(Phi(C,V)) + mu1(M_-1(C,V)) = 0   for all fillers      (lam·L = 0),
mu0(E-R)      - mu1(N_-1)       != 0                       (lam·rhs != 0),
```

which certifies that no filler completes `Q_0=1` and `Q_-1=0` simultaneously.
The block-0 equation uses `G-E` (whose filler-linear part is exactly `Phi`, so
`D | Phi`), **not** `Q_0-1 = (T-1)(G-E)` — this is what makes the necklace clean
(`verify §S1,§S2`).

## 1. The generic covector is explicit (Fredholm gap 1)

At the explicit branch-B slope-`1` datum (from `w2-decisive.md` §4), `L` has
**full column rank** and Fredholm gap `1` at every tested filler degree, so `lam`
is generically unique up to scale (`verify §S3`). It is extracted two ways:

- **Left kernel.** `lam` is the (up-to-scale-unique) left null vector of `L` with
  nonzero `rhs` pairing; `lam·rhs = 16/45 != 0` at the reference datum under the
  `sympy` null-space normalisation (only the **nonvanishing** is invariant; the
  scale/sign is not).
- **Pivot/residual formula** (the "explicit function of the datum"). Pick `nf`
  independent rows `I` of `L`; `x* = L_I^{-1} rhs_I`; any augmenting row `j` gives

```text
lam.rhs = L_j x* - rhs_j = L_j L_I^{-1} rhs_I - rhs_j   ( = -16/45 here, != 0),
```

  realised by `lam = e_j - L_j L_I^{-1}·(rows I)`. Its **denominator is `det L_I`**;
  where that vanishes one re-pivots (§5).

**Minimal genuine coupling.** Filler degree `dV=3` is the *first* degree where
block 0 (`Q_0=1`) **and** block 1 (`Q_-1`) are each individually feasible yet the
join fails (gap `1`). Lower `dV` gives a spurious gap that is only a filler-
truncation artifact of block 0 (its rank is below `nf`), not the joint obstruction
(`verify §S3`).

## 2. Closed form: the union of the two root necklaces (confirmed, bounded)

> **The point triple is silent; the covector lives on the block necklace.**
> A covector with `mu0` supported on `Ann(Im Phi) = {ev_-1,ev_0,ev_1}` (the roots
> of `D`) and **any** `mu1` has obstruction **`0`** (`verify §S4`). The reason is
> structural — degree-free in the covector argument, with two ingredients in-file
> bounded (see §6 ledger): (i) `mu0(Phi)=0` because `D | Phi`; (ii) on the slope locus
> `mu0(E-R)=0` because `(E-R)(±1)=0`; (iii) `mu1` alone cannot obstruct because
> `Q_-1=0` is feasible on its own. So `lam` **must leave** the point triple.

This is exactly why the complete point-annihilator dual of the W2 filler quotient
([`../band3/w2-theory.md`](../band3/w2-theory.md) §4) is *blind* to the joint kill:
it is `{ev_-1,ev_0,ev_1}`, and every member pairs `E-R` to `(beta-gamma)(1-R(1))=0`
on slope. The **joint** covector is genuinely new — it is coupled (`mu0(Phi) != 0`;
`verify §S3`).

> **Where it lives (confirms the union-of-necklaces conjecture).** A covector with
> `mu0` on the **fixed `(a_3,b_2)` root necklace**
>
> ```text
> N_fix = roots(a_3) ∪ roots(b_2) ∪ roots(q_K) ∪ roots(q_H)
>       = {0,-2,-4} ∪ {0,-3} ∪ {-1,0,1,2,3} ∪ {-1,0,1,2}  =  {-4,-3,-2,-1,0,1,2,3},
> ```
>
> (with `q_K = a_3(E-3)(E)_3`, `q_H = b_2(E-2)(E)_2` the shift-form necklace
> polynomials of [`lambda-general-k.md`](lambda-general-k.md) Thm A′) **couples**:
> a nonzero-obstruction covector exists on `N_fix` and an integer window for `mu1`
> (`verify §S4`). The smaller set `{-4,..,1} = roots(a_3) ∪ roots(b_2) ∪ roots(D)`
> is **insufficient** — the coupled covector needs the shifted membership window
> `{2,3}` of `q_K`.

So the fixed part of `lam` is the moving-sum/adjoint covector on the `(a_3,b_2)`
necklace, as conjectured. The datum-dependent `mu1` sits on the `(a_2,b_1)`
necklace; its exact node support is algebraic (the roots of `a_2,b_1` are not
integers for a generic datum), so `mu1` admits **no fixed integer-node closed form**
— it is represented by the pivot/residual formula. Comparing to the **both-ends
Lemma P**: `lam.rhs = mu0(E-R) - mu1(N_-1)`, and the `mu0(E-R)` term is the pairing
of the necklace covector against the Lemma-P boundary values (which carry `R(1)`,
`R(-1)`).

## 3. Filler-degree-freeness (bounded, `dV=4..10`)

At the reference datum, `L` has **full column rank, gap `1`, and nonzero
obstruction** at every filler degree `dV = 4,5,6,7,8,9,10` (`verify §S5`) — the
datum admits **no tail-clearing filler of any tested degree**. Structurally the
`mu0`/`Phi` block is degree-free by the adjoint criterion (a necklace covector
annihilates `Im Phi` at arbitrary degree, [`lambda-general-k.md`](lambda-general-k.md)
Thm A′); the `mu1` node support grows with `dV` (the `(a_2,b_1)` necklace is
algebraic), so the covector itself is not a single fixed functional across all
degrees, but the gap stays `1`.

## 4. The payload — a UNIFORM d=3–4 kill (the sharp bounded result)

> **`[Phi ; M_-1]` is the unit ideal at `d=3` and `d=4`, both branches.** The
> system `cascade ∧ Q_0=1 ∧ Q_-1 ∧ membership` (fillers included, raw cap `d`) is
> the **unit ideal over `QQ`** at raw cap `d=3` **and** `d=4`, on branch A and
> branch B (`msolve` with `^`, cross-checked mod `65003`; `verify §S7`,
> `d=4` behind `HEAVY=1`). So **no** slope-`1` cascade datum (raw cap `d≤4`) admits
> any filler (raw cap `d≤4`) satisfying `Q_0=1` and `Q_-1=0` — the Joint Filler
> Covector Lemma holds **uniformly**, not just generically, at these caps, and the
> full negative tail `Q_-2..Q_-5` is **not needed**.

**No false kill (controls).** `cascade` alone is feasible, and `{Q_0=1}` alone is
feasible in the fillers at the explicit datum (`sympy`, `verify §S7`); dropping
`Q_-1` breaks the kill. The independent `sympy` cross-check confirms `{Q_0=1,Q_-1}`
is the unit ideal **in the fillers** at the explicit datum on both branches — an
exact `QQ` Gröbner witness of the same obstruction the symbolic `msolve` run
certifies uniformly.

**Branch A (per gauge).** The `d=3,4` unit certificate already covers branch A
(`a_-3=(E)_3 am3`, `b_-3=mu_3 a_-3`). At the datum, for each fixed gauge
`mu_3 ∈ {1,-2,1/3}` the stacked operator has full column rank, gap `1`, and nonzero
obstruction (`verify §S8`).

## 5. The residual sub-locus is empty at d=3–4 (correction)

[`w2-joint-theorem.md`](w2-joint-theorem.md) §3 posited a *"special slope-`1`
sub-locus [that] survives `Q_-1`,"* with the *uniform* kill requiring the full tail
`Q_-1..Q_-5`. That rested on the §S8 `msolve` run **built with Python `**`
exponents**, which `msolve 0.10.1` silently misparses (a wrong "feasible" verdict;
[`hatch-census.md`](hatch-census.md) §5). Re-run with `^` (`verify §S7,§S9`), the
system is the **unit ideal** at `d=3,4` both branches. Hence:

- The **residual sub-locus is empty at `d=3,4`**: `[Phi ; M_-1]` is the uniform
  kill; adjoining `M_-2,...` (or the bilinear `Q_-5`) is unnecessary at these caps.
- The generic covector's denominator `det L_I` **can** still vanish on a proper
  subvariety (where a *particular* pivot degenerates); there one re-pivots to a
  different `lam`, and the obstruction persists (the unit ideal guarantees it). So
  "the covector loses this pivot" and "the obstruction survives" are consistent —
  the sub-locus of *pivot* degeneracy is not a sub-locus of *feasibility*.

Whether a genuine feasibility sub-locus (requiring `Q_-2..`) appears at raw cap
`d>4` is **not settled**.

The `**`→`^` parser bug is reproduced directly in file on the
[`hatch-census.md`](hatch-census.md) §5 minimal reproducer: a provably unit ideal
(`sympy`) that `msolve` returns as unit under `^` and as a spurious non-unit basis
under `**` (`verify §S9`).

> **A second, FLAGGED, `**`→`^` flip (not load-bearing here; needs independent
> re-derivation).** The other `msolve` feasibility sub-claim of
> [`w2-joint-theorem.md`](w2-joint-theorem.md) §S8 — *"`R(1)` is NOT in
> `sqrt(cascade+tail)`"*, the machine-checked pillar of that memo's §2 headline
> that the joint kill is *"not slope-forcing"* and `R(1)` is a *"nonvanishing free
> modulus"* on `cascade+tail` — **also flips under `^`.** Re-running the exact
> Rabinowitsch test with `^` (`cascade + Q_-1..Q_-5 + {1 - t·R(1)}`, branch A,
> `d=3`) returns the **unit ideal** over `QQ` (msolve, 258 s) **and** mod `65003`
> (124 s), i.e. **`R(1) ∈ sqrt(cascade+tail)`** — `R(1)=0` on the `cascade+tail`
> variety. That variety appears **positive-dimensional/feasible** (msolve did not
> terminate parametrising `cascade+tail` alone within 300 s), so the radical
> membership is **non-vacuous**, and it would make *"`R(1)=0` on `cascade+tail`"*
> **TRUE** — the opposite of `w2-joint-theorem.md` §2. **This is a single-method,
> heavy `msolve`-only verdict** on precisely the house's most-audited failure mode;
> it is **not** used anywhere in this memo's results, and it is recorded as a
> **flagged discrepancy to be independently re-derived** (a `sympy`/second-engine
> radical certificate, or an explicit `cascade+tail` witness with `R(1)≠0`) before
> the "not slope-forcing" framing is trusted or overturned. Reproduce:
> `cascade+tail + {1 - t·R(1)}` via `msolve` with `^` (see `verify_w2_joint.py`
> §S8, which now carries the `^` fix at line 570).

## 6. Evidence ledger — proved vs bounded vs refuted

**Proved (arbitrary coefficient degree, char 0, machine-checked identities):**
- `Q_m=[D,X]_m` for `m in [-6,6]`; `Q_0=(T-1)G`, `G=` potential — the
  crossed-product engine as symbolic identities in generic (fully symbolic
  degree-2) coefficients across all seven levels (`§S0`).
- **Slope gate:** `D=E(E-1)(E+1)` divides `E-R` iff `R(1)=1` and `R(-1)=-1`
  (`R(0)=0` automatic) — the elementary root-evaluation fact (`§S1`).

**Bounded / structural (AUDIT-DEMOTED from "proved arbitrary degree"):**
- **Point-triple decoupling:** any `lam` with `mu0 ∈ Ann(Im Phi)={ev_-1,ev_0,ev_1}`
  has obstruction `0` on the slope locus (`§S4`). The covector *argument* is
  degree-free, but two of its three ingredients are only in-file **bounded**:
  (a) `Im Phi(W2) ⊆ D·F[E]` is re-verified at filler degree `jdeg ≤ 5` (`§S1`),
  its arbitrary-degree proof being the cited triangular argument of
  [`../band3/quantum-ap-filler-image.md`](../band3/quantum-ap-filler-image.md);
  (b) `Q_-1`-feasibility-alone is instance-checked at the reference datum (`§S3`);
  and the machine check itself is a single datum at `dV=6`. Structural + bounded,
  not an arbitrary-degree machine identity.

**Bounded / generic evidence (exact, stated scope):**
- The **explicit generic covector** (left kernel and pivot/residual formula),
  `lam.rhs != 0`, at the reference datum (`§S3`).
- **Union-of-necklaces closed form:** `mu0` on `N_fix={-4,..,3}` couples; `{-4,..,1}`
  is insufficient; the point triple decouples — at the datum, `dV=6` (`§S4`).
- **Filler-degree-freeness:** full column rank, gap `1`, `lam.rhs != 0` at
  `dV=4..10` (`§S5`).
- **Generic in positive data:** 6 distinct slope-`1` data, gap `1`, `lam.rhs != 0`
  (`§S6`).
- **Branch A per gauge** `mu_3 ∈ {1,-2,1/3}`: gap `1`, `lam.rhs != 0` (`§S8`).

**Bounded exact certificate (the payload):**
- **`cascade + Q_0=1 + Q_-1` is the UNIT ideal** at raw cap `d=3` (`QQ` + mod
  `65003`) **and** `d=4` (`QQ`, `HEAVY`), on **both** tail branches (`msolve` `^`),
  with a `sympy` `QQ` cross-check in the fillers at the explicit datum (`§S7`). The
  Joint Filler Covector Lemma holds **uniformly** at raw cap `d≤4`.

**Refuted / corrected:**
- The sibling claim *"a slope-`1` sub-locus survives `Q_-1`; the uniform kill needs
  the full tail `Q_-1..Q_-5`"* (`w2-joint-theorem.md` §3, §S8) — a `msolve` `**`
  parser artifact. Corrected: `[Phi ; M_-1]` is the **uniform** kill at `d=3,4`,
  the sub-locus is **empty** there, and `Q_-2..Q_-5` are not needed (`§S7,§S9`).
- A degree-free **constant** value of `lam.rhs`: the obstruction is
  **datum-dependent** (`16/45, -16/15, -32/15, 16/15, 32/15, ...` across the six
  data under a fixed normalisation), so there is no universal constant analogous to
  the W1 `lambda_0(E-R)=4`.

**Flagged (recorded, NOT claimed — needs independent re-derivation):**
- The `w2-joint-theorem.md` §2 pillar *"`R(1)` not in `sqrt(cascade+tail)`"* also
  flips under `^` (`R(1) ∈ sqrt(cascade+tail)` per `msolve` `QQ` + mod `p`;
  `cascade+tail` appears feasible), which would make *"`R(1)=0` on `cascade+tail`"*
  **true** and reopen whether W2's kill is slope-forcing (§5). A single heavy
  `msolve`-only verdict; **flagged**, not used here (see §5).

**Open / not claimed:**
- The Joint Filler Covector Lemma at **arbitrary positive-data degree** (`d>4`):
  the uniform kill is only a bounded `d≤4` certificate; a degree-free proof needs
  the adjoint machinery for the `(a_2,b_1)` necklace, whose nodes are algebraic.
- A **degree-free `lam.rhs` identity** (the analogue of the W1 `lambda_0(R)=0`
  elimination): not obtained; `lam.rhs` is a datum-dependent rational function,
  nonzero on the whole locus at `d≤4` (by the unit certificate) but not reduced to
  a closed nonvanishing invariant.
- Whether a genuine feasibility sub-locus needing `Q_-2..` appears at `d>4`.
- **Filler degree freeness** beyond `dV=10` is evidence, not proof; and the full
  arbitrary-degree W2 joint theorem, all of Band 3, DC1, JC2 remain open. No Weyl
  pair and no counterexample is constructed. The infinite-dimensional
  `Im L_K ∩ Im L_H`
  ([`two-filler-cross-cancellation.md`](two-filler-cross-cancellation.md)) is
  untouched.

## 7. What it hands architecture step 2

The reduction now hands step 2 a **uniform, computable** object at `d≤4`: the joint
cokernel of `[Phi ; M_-1]`, with its covector's fixed part pinned to the `(a_3,b_2)`
root necklace. The single remaining degree-free gap is the `(a_2,b_1)` necklace
adjoint — the `mu1` block — because its nodes are the (algebraic) roots of the
solved `a_2,b_1`. A degree-free description of `mu1` (e.g. as an `S^*`-transported
functional on the datum necklace, à la
[`lambda-general-k.md`](lambda-general-k.md) Thm A′, now with algebraic support)
would lift the bounded `d≤4` kill to arbitrary positive-data degree — the sharp
next computation.

## 8. Verification

```sh
uv run --with sympy python research/dc1-program/verify_joint_covector.py
# optional heavier legs (d=4 joint unit over QQ, both branches):
HEAVY=1 uv run --with sympy python research/dc1-program/verify_joint_covector.py
```

Runs `§S0` (engine, `Q_0=(T-1)G`), `§S1` (`Im Phi ⊆ D·F[E]`, slope gate), `§S2`
(explicit datum, stacked `L`), `§S3` (explicit covector: left kernel + pivot
formula; minimal coupling `dV=3`), `§S4` (point-triple decoupling; the necklace
`N_fix` couples, `{-4,..,1}` insufficient), `§S5` (filler-degree-freeness
`dV=4..10`), `§S6` (six slope-`1` data), `§S7` (the `d=3` — and `HEAVY` `d=4` —
unit certificate over `QQ` + mod `p`, both branches, with `sympy` controls and
cross-check), `§S8` (branch A per gauge), `§S9` (the `**`→`^` `msolve` parser bug
reproducer). A successful run ends `ALL JOINT COVECTOR CHECKS PASSED`.
