# W2 decisive: the achievable slope R(1) at the r = -4 exotic member

**INDEPENDENTLY DERIVED — EXACT ALGEBRA — NOT PEER REVIEWED — BAND-SCOPED**

W2 is the `r = -4` member of the band-3 step-2 arithmetic-progression exotic
family,

```text
a_3 = E(E+2)(E+4)   (roots {-4,-2,0}),   b_2 = E(E+3),   b_3 = 0,
```

in the quantum band-3 conventions
`Q_m = sum_(k+l=m)[b_l^[k] a_k - a_k^[l] b_l]`, `f^[n](E)=f(E+n)`,
`G = sum_(k=1)^3 sum_(j=0)^(k-1)(a_k^[j-k] b_-k^[j] - b_k^[j-k] a_-k^[j])`,
`Q_0=(T-1)G`, membership `(E)_j | a_-j, b_-j`.

The lambda-wave theorem (commit `d8189fc`,
[`quantum-ap-lambda.md`](quantum-ap-lambda.md)) closes every AP member `r != -4`
at **arbitrary degree** via `lambda_r(f)=f(r+3)-f(r+4)+f(r+5)-f(0)`,
`lambda_r(E-R)=r+4`. At `r=-4` the functional degenerates
(`lambda_{-4}(E)=0`) and the entire closure collapses to one scalar: **the moment
slope `R(1)`**. This note settles that scalar.

## Verdict

**`R(1) = 1` is achievable — the slope obstruction does NOT close W2.**

- `d = 1`: the positive cascade forces `R = 0`, so `R(1) = 0`.
- `d = 2`: two branches, both with `R in E(E-1)(E+1)F[E]`, so `R(1) = 0`.
- `d >= 3`: **`R(1)` is a genuine free modulus** (shape (b)). It is *not* forced to
  `0`; on a displayed slice `R(1) = -(8/9) am1_3` where `am1_3` is the top
  (`E^4`) coefficient of `a_-1`, and `am1_3 = -9/8` produces `R(1) = 1`.

An **explicit `d = 3` solution with `R(1) = 1`** is reconstructed from scratch and
verified as exact polynomial identities: the positive cascade
`Q_4=Q_3=Q_2=Q_1=0`, full Weyl membership, and — with constructed fillers — the
central identity `G = E`, hence `Q_0 = 1`. Because the solution is a set of
concrete polynomials of bounded degree, it is valid at **every `d >= 3`**.

Consequently the band-3 exotic sector survives past the slope. Only the negative
tail `Q_-1..Q_-6 = 0` now stands between W2 and a DC1 counterexample candidate.
This note imposes **no** negative-tail equation and constructs **no** Weyl pair.

Exact certificate: [`verify_w2_decisive.py`](verify_w2_decisive.py) (ends
`ALL W2 DECISIVE CHECKS PASSED`).

## 0. The pivot fact, re-derived (not cited)

`Im Phi(W2) = E(E-1)(E+1) F[E]` **exactly**. Write `D = E(E-1)(E+1)`. Two facts,
both machine-checked here (not merely quoted from
[`quantum-ap-filler-image.md`](quantum-ap-filler-image.md)):

- **Subset.** Every basis filler `K_3[(E)_3 E^j]`, `H_2[(E)_2 E^j]` is divisible
  by `D`.
- **Equality.** At each truncation `deg <= N` (`N = 9,10,11`), the image rank
  equals `dim(D*F[E])_{<=N} = N-2` and the image lies inside `D*F[E]`; the
  codimension stabilizes at `3 = deg D`. Hence `Im Phi = D*F[E]`.

Because `Q_0=(T-1)G`, `(T-1)E=1`, and membership gives `G(0)=0`, the unit equation
`Q_0=1` is equivalent to `G=E`, i.e. `E - R in Im Phi(W2)`. Therefore

```text
Q_0 = 1   <=>   E - R in D*F[E]   <=>   D | (E-R)   <=>   R(0)=0, R(1)=1, R(-1)=-1,
```

with `R(0)=0` automatic from membership. The proved cascade constraint
`R(1)+R(-1)=0` (below) is **consistent** with `(R(1),R(-1))=(1,-1)`; W2 lives or
dies on whether the positive cascade can produce `R(1)=1`.

## 1. The structural reduction: R(1) is the moment slope Q_0(0)

The two filler blocks lie in `D*F[E]`, so they vanish at `E = 0, 1, -1`. Hence `R`
and `G` agree at those three nodes, and (with `G(0)=0`, `Q_0=(T-1)G`)

```text
R(1)  = G(1)  = Q_0(0),
R(-1) = G(-1) = -Q_0(-1),
R(0)  = G(0)  = 0.
```

So the achievable-slope question is exactly: **does the positive cascade force the
moment slope `Q_0(0) = 0`?** The proved constraint `R(1)+R(-1)=0` reads
`Q_0(0) = Q_0(-1)`.

## 2. Degree ledger of the forcing

The positive cascade is solved by the sequential `clean_solve` reconstruction
(identical machinery to
[`verify_quantum_exotic_cokernel.py`](verify_quantum_exotic_cokernel.py); the sole
retained solver kernel is the constant `b_0` freedom, which does not enter `R`).
For the free data capped at raw degree `d`:

| `d` | positive-cascade solution | slope `R(1)` |
|---|---|---|
| 1 | single branch, `R = 0` | `0` (forced) |
| 2 | two branches `a1_2=0` / `a2_0=0`, both `R in D*F[E]` | `0` (forced) |
| >= 3 | positive-dimensional variety | **free modulus** |

The `d=1,2` rows reproduce the committed finite closure
([`quantum-exotic-cokernel.md`](quantum-exotic-cokernel.md)): on both `d=2`
branches
`R = -(2/3)E^2 am1_2 (E-1)(E+1)(a2_0-2a2_2)` and
`R = (2/3)E(E-1)(E+1)(2E a2_2 am1_2 + a0_2 a1_2)` — each divisible by `D`, so
`R(1)=R(-1)=0` and `E-R notin Im Phi`, i.e. `Q_0=1` is infeasible at `d<=2`.

**The forcing verified by an explicit solution family (no Groebner needed).** At
`d=3`, fixing the base data on a rational slice and keeping `am1_3` symbolic, the
positive-cascade compatibility conditions solve **linearly** for the remaining
coefficients, producing an honest one-parameter family of genuine positive-cascade
solutions on which

```text
R(1) = -(8/9) am1_3,      R(-1) = +(8/9) am1_3.
```

So `R(1)` is surjective onto `F` (a free modulus), while `R(1)+R(-1)=0` holds
identically — the proved constraint is respected but does *not* pin `R(1)`. This is
a stronger and faster certificate of non-forcing than a radical-ideal membership
test (an exhibited solution family beats a `sqrt(I)` computation); an independent
Rabinowitsch check confirms `R(1) in sqrt(I)` at `d=2` and `R(1) notin sqrt(I)` at
`d=3`, the same flip.

## 3. The decisive computation and the structural law

At `d = 3`, forward-solving `Q_4=Q_3=Q_2=Q_1=0` with free data and computing the
slope as a function of the free parameters gives, on a displayed rational slice
(`a2 = 2 + (9/2)E + (1/2)E^2 - (1/3)E^3`, `a1_3 = 5/7`, remaining data solved
linearly):

```text
R(1) = -(8/9) am1_3,      R(-1) = +(8/9) am1_3,
```

where `am1_3` is the leading coefficient of `am1_raw` in `a_-1 = E * am1_raw`, i.e.
the **`E^4` coefficient of `a_-1`**. The proportionality constant is *not*
universal — different base slices give `-8/9, 16/9, -4/9, ...` — so `R(1)` is a
genuine bilinear functional of the free data, nonzero for generic data. The single
structural law that *is* uniform:

> **Slope law (W2).** On the positive-cascade variety, the moment slope
> `R(1) = Q_0(0)` is a nonconstant linear functional of the degree-`(d+1)` leading
> data of `a_-1`. It vanishes identically for `d <= 2` (that leading datum does not
> exist) and is surjective onto `F` for `d >= 3`. The only relation the positive
> cascade + membership place on the slope pair is `R(1)+R(-1)=0`; they do **not**
> pin `R(1)`.

**Why the degree threshold is `3`.** `am1_3` is the `E^4`-coefficient of `a_-1`. At
`d <= 2` the free polynomial `am1_raw` has degree `<= 2`, so `a_-1` has degree
`<= 3` and this modulus is absent — the slope collapses to `0`. At `d >= 3` it
appears, is unconstrained by the cascade compatibilities (which only tie the lower
coefficients of `a_-1` to the `a`-tower), and drives `R(1)` to any target.

**Why the `r`-uniform functional method left this open (honest reconciliation).**
The committed certificate ([`quantum-ap-lambda.md`](quantum-ap-lambda.md) §5e)
observed that its `r`-shifted elimination pins the alternating combination
`R(1)+R(-1)` to `0` but leaves `R(1) = R(r+5)` a nonzero free modulus. That memo
deliberately works in `r`-shifted evaluation variables treated as *independent*
(the mechanism that makes the `r != -4` obstruction degree-free), and so discards
the absolute-node coupling near `E=0`. At `r=-4` the wall node `r+4` **coincides
with the membership anchor `E=0`** (the "common-root-at-0" geometry): the relaxation
is exactly too weak to see whether the slope is pinned. This note resolves the
question in the *actual* polynomial system — where the coupling is present — and
finds the slope genuinely free. The two statements agree; nothing in the committed
work is contradicted (it explicitly left `r=-4` open and reduced it to this slope).

## 4. The explicit R(1)=1 solution (positive sector + fillers, verified)

Pinning `am1_3 = -9/8` and freezing the residual free moduli to `0` gives the
concrete `d = 3` datum below. It satisfies, as **exact polynomial identities**:
`Q_4=Q_3=Q_2=Q_1=0`; membership `(E)_j | a_-j, b_-j` (`j=1,2,3`); `R(0)=0`,
`R(1)=1`, `R(-1)=-1`, hence `D | (E-R)`; and, with the fillers below,
`G = E` so `Q_0 = 1`.

```text
a_3  = E(E+2)(E+4)
a_2  = -E^3/3 + E^2/2 + 9E/2 + 2
a_1  =  5E^3/7 + 107E^2/63 + 118E/63 + 10/9
a_0  = -775E^3/5103 - 92545E^2/3402 - 277597E/10206
a_-1 = -9E^4/8 - 219830E/11907
a_-2 = E^5/5 - E^4/5 + 967027E^3/2755620 - 53597707E^2/8573040 + 455302607E/77157360
a_-3 = 0                              (free for the tail; enters only Q_-1..Q_-6)

b_2  = E(E+3)
b_1  = -2E^2/9 + 8E/9 + 4/3
b_0  =  263E^2/567 + 179E/567
b_-1 = -256E^2/5103 - 31130E/1701
b_-2 = -3E^3/4 + 2747993E^2/2571912 - 819059E/2571912
b_-3 = 2E^4/15 - 5E^3/12 + 19E^2/60 - E/30        (= (E)_3 C, the K_3 filler)
b_3  = 0
```

The two fillers realize `Phi(C,V) = E - R` (solvable precisely because
`E-R in D*F[E] = Im Phi`):

```text
b_-3 = (E)_3 * C,   C = (8E-1)/60,
a_-2 = (E)_2 * V,   V = (15431472 E^3 + 27076756 E - 455302607)/77157360.
```

**Persistence.** These are concrete polynomials; `Q_m=0`, membership, and `G=E`
are identities independent of any coefficient-degree cap. The raw degrees
(`a2,a1,a0` degree 3; `am1_raw` degree 3; `V` degree 3; `C` degree 1) all fit within
cap `d` for every `d >= 3`. Hence `R(1)=1` (with the positive cascade, membership,
and `Q_0=1`) is achievable at **every `d >= 3`**, not just `d = 3`.

## 4b. Where this sits among the sibling flanks (synthesis)

Two sibling memos landed at the same commit and make the slope the *sole* decider:

- [`w2-theory.md`](w2-theory.md) proves `Im Phi(W2)` is the **principal squarefree
  ideal `(D)`**, so its cokernel is semisimple and its entire dual is the three
  point functionals `{ev_-1, ev_0, ev_1}` — **no non-point annihilator exists**.
  It also proves **reflection does not close W2** (the Fourier `phi` fixes the
  anchor at `0` and reverses the band, sending W2's degree-3 top to a degree-6
  bottom). Hence there is no hidden functional and no symmetry shortcut: W2 lives
  or dies *purely* on `R(1)=1`. This note supplies that scalar: **`R(1)=1` is
  achievable.**
- [`w2-negative-tail.md`](w2-negative-tail.md) shows the negative tail, taken alone
  with the positive cascade, is a **proper** (feasible) ideal at `d=1,2` — the slope
  was the only obstruction there — and ships an escalation harness
  `run_full_system` awaiting a slope-`1` point.

Together the census in `w2-theory.md` §3 makes W2 the band-3 instance of a
canonical **step-`(k-1)`/step-`k` hatch** per band `k >= 3`, each reducing to the
*same* gate `R(1)=1`. So the achievability shown here is the template result: if it
generalizes, every such hatch is a live candidate.

## 5. Handoff to the negative-tail agent

What is now handed to `W2-TAIL` (the `Q_-1..Q_-6` search):

- The **positive sector + central `Q_0=1`** is solved, explicitly, at `R(1)=1`
  (§4). The negative tail is the *only* remaining gate.
- The remaining freedom for the tail: (i) the `R(1)=1` fiber of positive data at
  `d>=3` (fix `am1_3=-9/8`, the other free moduli remain), (ii) the **filler
  kernel** — `Phi` has a 1-parameter kernel in `(C,V)` at each degree, so `b_-3`,
  `a_-2` are free modulo `ker Phi` while keeping `G=E`, and (iii) `a_-3`, which is
  **absent from `G`** (its `G` term carries `b_3 = 0`) and so is entirely free for
  the tail — it enters only `Q_-1..Q_-6`.
- Escalation trigger for the tail: if a choice of `(a_-3, filler kernel, R(1)=1
  fiber)` clears all six negative-tail equations with membership, a **full candidate
  pair** materializes — then verify `[D,X]=1` by direct crossed-product
  computation, run the char-`p` sieve
  ([`../dc1-program/sieve_dc1_candidate.py`](../dc1-program/sieve_dc1_candidate.py)),
  and attempt the `A_1`-generation test.

## 5b. The combined slope+tail probe (escalation status)

With `R(1)=1` now achievable, the only remaining gate to a DC1 candidate pair is
the negative tail. Using the sibling harness structure (positive cascade + `Q_0=1`
+ `Q_-1..Q_-5` + membership), the combined DC1-face feasibility was probed:

- **`d=2` control:** the full combined ideal is the **unit ideal** (infeasible),
  reproducing [`w2-negative-tail.md`](w2-negative-tail.md) — as it must, since the
  slope is forced `0 != 1` at `d<=2`.
- **The explicit `R(1)=1` point does not, by itself, extend to the tail.** Pinning
  the Section-4 positive datum and searching the full filler freedom, the negative
  tail `Q_-1..Q_-5` is **infeasible on both branches**: `a_-3 != 0` (with
  `b_-3 = mu_3 a_-3`) and `a_-3 = 0` (with `b_-3` a free `(E)_3`-filler). So this
  particular slope-`1` datum is not a candidate pair.
- **Full-fiber `d=3`** (all slope-`1` positive data free, only `Q_0=1` enforcing the
  slope, plus the whole tail) is a large Gröbner computation; its verdict is the
  proper hand-off to `W2-TAIL` and is not asserted here.

The honest reading: the slope is achievable, but the negative tail is **not
automatically** satisfiable at the same point — the tail agent's "no obstruction
beyond the slope" (bounded to `d=1,2`, where the slope itself was the only live
constraint) is *not* established to mean simultaneous slope+tail feasibility at
`d=3`. Whether some slope-`1` datum (at `d=3` or higher) clears the tail is the
open escalation question. **No candidate counterexample pair is asserted.**

## 6. Sanity controls (no false kills)

- The `d=3` positive-cascade ideal is **proper** (its Gröbner basis is not `[1]`);
  solutions exist, so the method is not manufacturing kills.
- The same finite-cokernel test that *accepts* `E-R at R(1)=1` (`d=3`) correctly
  *rejects* `Q_0=1` at `d=2`, where `R in D*F[E]` forces `E-R = E - (D-multiple)`,
  which is not in `D*F[E]` because `D` does not divide `E`. Feasible and infeasible
  cases are distinguished by one uniform test.

## 7. Scope

**Proved (exact algebra):**
- `Im Phi(W2) = E(E-1)(E+1)F[E]` exactly; `Q_0=1 <=> R(1)=1, R(-1)=-1`.
- `R(1) = Q_0(0) = G(1)` (structural reduction).
- Slope forced to `0` at `d <= 2`; slope a free modulus at `d = 3` (radical-ideal
  certificate, both directions).
- An **explicit** positive-cascade + membership + fillers datum with `Q_0 = 1`
  (`R(1) = 1`), verified as polynomial identities, hence a valid solution at
  **every `d >= 3`**.

**Exceptional locus.** The forcing `R(1) = 0` is exactly the low-degree regime
`d <= 2`; the escape is the appearance of the `E^4`-coefficient of `a_-1` at
`d >= 3`. This is the `r = -4` specialization of the lambda-wave escape hatch:
the common wall root sits at the membership anchor `E = 0`.

**Not proved / out of scope.** The negative tail `Q_-1..Q_-6 = 0` (handed to
`W2-TAIL`); any Weyl pair or DC1/JC2 counterexample; the non-AP `deg a_3 >= 6`
exotic tops; anything outside band 3. In unrestricted degree
`Im L_K intersect Im L_H` remains infinite-dimensional
([`../dc1-program/two-filler-cross-cancellation.md`](../dc1-program/two-filler-cross-cancellation.md));
nothing here weakens that.

Run:

```sh
uv run --with sympy python research/band3/verify_w2_decisive.py
```

A successful run ends with `ALL W2 DECISIVE CHECKS PASSED`.
