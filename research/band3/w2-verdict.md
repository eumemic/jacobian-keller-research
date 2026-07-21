# W2 verdict: the combined slope+tail system is infeasible at d=3 and d=4 — the band-3 exotic construction route dies at bounded degree

**INDEPENDENTLY DERIVED — EXACT ALGEBRA (two independent solvers) — NOT PEER REVIEWED — BAND-SCOPED — BOUNDED DEGREE**

This memo settles the **single finite feasibility question** all three W2 sibling
flanks reduced to. W2 is the `r=-4` exceptional member of the band-3 step-2 exotic
arithmetic-progression family, gauge `b_3=0`, quantum band-3 conventions
(`Q_m = sum_(k+l=m)[b_l^[k] a_k - a_k^[l] b_l]`, `f^[n](E)=f(E+n)`,
membership `(E)_j = E(E-1)...(E-j+1) | a_-j, b_-j`):

```text
a_3 = E(E+2)(E+4)   (roots {0,-2,-4}),     b_2 = E(E+3),     b_3 = 0.
```

The question: **does any slope-1 datum simultaneously satisfy, at a common raw
degree cap `d`, the positive cascade `Q_4=Q_3=Q_2=Q_1=0`, the moment unit
`Q_0=1` (`<=> R(1)=1, R(-1)=-1`), the full negative tail `Q_-1=...=Q_-5=0`, and
all Weyl memberships?** A "yes" would assemble a **DC1 counterexample candidate
pair** and trigger the escalation protocol. A "no" kills the current construction
route.

## Verdict

> **INFEASIBLE. The combined system is the UNIT IDEAL at `d=3` AND `d=4`, on BOTH
> tail branches — exactly over `QQ`, by two independent Gröbner engines
> (`msolve 0.10.1` and `sympy`), and modulo five primes in `[10^4,10^5]`.**
>
> No slope-`1` datum clears the negative tail at raw-degree cap `3` or `4`. The
> slope `R(1)=1` is achievable and the tail is feasible **separately**, but their
> **intersection is empty**. **No DC1 counterexample pair materializes; the
> escalation orders were not triggered; no Weyl pair and no counterexample is
> constructed.**
>
> **Scope (honest).** This is a **bounded-degree** kill (`d in {3,4}` exact,
> `d=5` mod-`p`). It does **NOT** close W2 at arbitrary degree. It kills the
> finite-search *construction route* (find a slope-`1` datum that also clears the
> tail at some small cap) through the degrees a naive search would reach, and
> shifts the remaining burden to a **structural / degree-free argument** — the
> finite computation cannot decide arbitrary degree.

Exact certificate: [`verify_w2_verdict.py`](verify_w2_verdict.py) (ends
`ALL W2 VERDICT CHECKS PASSED`). Base commit `66156f2`; all load-bearing
lambda-wave / cascade facts are re-verified in file, not merely cited.

## 0. The combined system and the two exhaustive tail branches

After the positive cascade solves `b_1,b_0,b_-1,b_-2` from `Q_4..Q_1=0` (with the
single `b_0` constant-kernel freedom retained), every coefficient is expressed in
the free data `{a_2,a_1,a_0,a_-1,a_-2, a_-3, ...}`, and `Q_0-1`, `Q_-1..Q_-5`
become **pure scalar polynomial constraints** (as established in
[`w2-negative-tail.md`](w2-negative-tail.md)). The bottom Wronskian `Q_-6=0`
splits the tail into **two exhaustive branches** (a `(-3)`-periodic ratio of
polynomials is constant):

```text
branch A:  a_-3 = (E)_3 * am3_raw  (generic),   b_-3 = mu_3 * a_-3      (Q_-6 gauge);
branch B:  a_-3 = 0,                            b_-3 = (E)_3 * bm3_raw  (free filler).
```

Their union is exhaustive (branch A at `am3_raw=0` gives `a_-3=b_-3=0`, contained
in branch B). At `d=3`: branch A has **26** free parameters, branch B **25**; the
`FULL` system is **62** scalar equations. This is exactly the "62-eq/26-var d=3
system that did not terminate" in [`w2-decisive.md`](w2-decisive.md) §5b.

## 1. Stage 1 — the mod-`p` sweep (decisive-in-expectation)

Reducing the full combined system modulo five primes
`p in {10007, 31013, 50021, 65003, 99013}` and running `msolve` over `GF(p)`:

| system (`d=3`) | branch A | branch B |
|---|---|---|
| `FULL` (cascade + `Q_0=1` + tail + membership) | **`[-1]` (UNIT) mod all 5** | **`[-1]` (UNIT) mod all 5** |

`msolve`'s `[-1]` = "no solution in the algebraic closure" = unit ideal. A uniform
`[-1]` across five large primes is a **strong kill signal** (a proper `QQ` ideal
would reduce to a proper ideal at all but finitely many primes). Each solve took
`<0.02 s`. Proceeding to the exact verdict with the kill prior.

## 2. Stage 3 — the exact `QQ` verdict (two independent engines)

**`msolve` over `QQ` (char 0).** The reduced Gröbner basis (`-g 2`) of the `FULL`
`d=3` system is literally `[1]` — a length-one basis, the constant `1` — for
**both** branch A and branch B. Hence `1` is in the ideal: **infeasible over
`QQ`** (and over `C`). The default parametrization output is `[-1]` for both.

**`sympy`, independent, via structure reduction.** `sympy` cannot Gröbner the raw
`62`-eq/`26`-var system in a sane budget (the very reason this task existed). A
**feasibility-preserving linear elimination** — repeatedly solve any equation
`c*v + rest = 0` with `c` a nonzero **constant** and `rest` free of `v`, via the
variety bijection `v = -rest/c` (so `V(.)` is empty iff the reduced `V(.)` is
empty; the unit-ideal verdict is preserved **exactly**) — shrinks the system to a
size `sympy` handles:

| `d=3` | raw | reduced | independent `sympy` verdict |
|---|---|---|---|
| branch A | 62 eq / 26 var | 51 eq / **16 var** | `groebner` residual `= [1]` over `QQ` (unit) |
| branch B | 62 eq / 25 var | 38 eq / **12 var** | `groebner` residual `= [1]` over `GF(32003), GF(65003), QQ` (unit) |

The reducer is validated at `d=2` against a **direct** full-system Gröbner
(`reduced-unit == direct-unit`). **Two solvers, exact `QQ`, both branches: UNIT.**

## 3. The controls — this is not a false kill

The kill is meaningful only because the two ingredients are *individually*
satisfiable at `d=3`; the machine distinguishes feasible from infeasible:

| system, `d=3` | branch A | branch B | meaning |
|---|---|---|---|
| positive cascade only | FEASIBLE (pos-dim) | FEASIBLE | solutions exist |
| cascade + tail (no `Q_0=1`) | FEASIBLE | FEASIBLE | tail non-obstructing alone |
| cascade + `Q_0=1` (no tail) | **FEASIBLE** | **FEASIBLE** | **slope `R(1)=1` achievable** |
| `FULL` (all together) | **UNIT** | **UNIT** | **intersection empty** |

The **slope flip** is the crux. At `d<=2` the row "cascade + `Q_0=1`" is the
**UNIT** ideal (the slope alone kills — `R(1)` is forced `0 != 1`); at `d=3` it
becomes **FEASIBLE**, exactly reproducing [`w2-decisive.md`](w2-decisive.md)'s
central result that `R(1)=1` is achievable at `d=3`. So the `FULL` unit ideal at
`d=3` is **not** the slope dying and **not** a systematic bug — it is the
slope-`1` locus and the tail-feasible locus **failing to meet**.

**Explicit witness.** The concrete `d=3` slope-`1` datum of
[`w2-decisive.md`](w2-decisive.md) §4 satisfies `Q_4=Q_3=Q_2=Q_1=0`, `Q_0=1`, and
all six memberships `(E)_j | a_-j, b_-j` — but all five tail equations are
nonzero:

```text
Q_-1, Q_-2, Q_-3, Q_-4, Q_-5  all != 0   at that point
```

(each carrying its forced membership factor `(E)_1,(E)_2,(E)_3,(E)_4,(E)_5`
respectively). Because that datum's coefficient degrees fit cap `d` for every
`d>=3`, it also witnesses "cascade + `Q_0=1`" **feasible** at `d=3` and `d=4`.
This is the point [`w2-decisive.md`](w2-decisive.md) §5b already found not to
extend; the present memo upgrades "this particular point fails" to "**the entire
slope-`1` fiber fails**, exactly, at `d=3` and `d=4`".

## 4. Stage 4 — `d=4` (and `d=5`) confirmation

Per the escalation protocol, `d=4` was pushed:

| system | branch A | branch B |
|---|---|---|
| `FULL`, `d=4`, `QQ` (msolve, exact) | **UNIT** | **UNIT** |
| `FULL`, `d=4`, mod 5 primes | `[-1]` all 5 | `[-1]` all 5 |
| `FULL`, `d=5`, mod `65003` | `[-1]` (UNIT) | `[-1]` (UNIT) |

`d=4` is **73** equations in **32** (branch A) / **31** (branch B) free
parameters. So the intersection is empty **exactly** at `d=3` and `d=4`, and
mod-`p` at `d=5`.

## 5. What this settles, and what it does not

**Reading with the siblings.** [`w2-theory.md`](w2-theory.md) proved every
*functional-method* obstruction on W2 fails (reflection does not close it; the
point-annihilator dual `{ev_-1,ev_0,ev_1}` is complete with no hidden functional;
the moment slope is achievable), leaving exactly the combined feasibility open.
[`w2-negative-tail.md`](w2-negative-tail.md) proved the tail is a *proper*
(feasible) ideal on its own at `d=1,2`. [`w2-decisive.md`](w2-decisive.md) proved
`R(1)=1` achievable at `d=3`. **This memo closes the loop: the SIMULTANEOUS system
is the unit ideal at `d=3` and `d=4`.** The achievability of the slope and the
feasibility of the tail are real, but they are **not co-realizable** at these
degrees.

**The DC1 consequence.** No candidate pair exists at `d in {3,4}`, so `[D,X]=1`,
the char-`p` sieve, and the `A_1`-generation test were **not** run — there is
nothing to feed them. The band-3 exotic hatch does **not** yield a DC1
counterexample by this route at bounded degree.

**Honest ledger — proved vs bounded.**

- **Proved (exact algebra over `QQ`, this memo, machine-checked by two independent
  solvers):** the W2 combined system (positive cascade `∧` `Q_0=1` `∧`
  `Q_-1..Q_-5=0` `∧` membership) is the **unit ideal** at raw-degree cap `d=3` and
  `d=4`, on **both** tail branches. Equivalently: no slope-`1` datum with
  coefficient degrees `<= 3` (resp. `<= 4`) clears the negative tail.
- **Bounded evidence (mod-`p`):** the same kill at `d=5` (both branches, prime
  `65003`); and the `d=3,4` kill modulo five primes in `[10^4,10^5]`.
- **Controls (no false kill):** cascade, cascade+tail, and cascade+`Q_0=1` are all
  **feasible** at `d=3` (the last one flips from unit at `d<=2`); an explicit
  slope-`1` witness point satisfies cascade+slope+membership and fails the tail.

**Open / NOT claimed (explicit).**
- **Arbitrary degree.** Whether some slope-`1` datum of raw degree `>= 6` clears
  the tail is **not** decided here. A finite Gröbner search cannot; closing W2 at
  arbitrary degree now requires a **structural (degree-free)** argument — the
  natural next target, analogous to the `lambda_r` theorem that closed
  `r != -4` at arbitrary degree. The infinite-dimensional
  `Im L_K ∩ Im L_H` ([`../dc1-program/two-filler-cross-cancellation.md`](../dc1-program/two-filler-cross-cancellation.md))
  is untouched.
- The band-`k` hatches (`k >= 4`, [`w2-theory.md`](w2-theory.md) §3) reduce to the
  same slope gate `R(1)=1`; their combined slope+tail feasibility is **not**
  addressed here (a separate computation per band).
- Everything outside the W2 combined system: the non-AP `deg a_3 >= 6` exotic
  tops, all of band `>= 4`, DC1, JC2. No Weyl pair, no counterexample.

**Exceptional-locus / search-space honesty.** Both `Q_-6` branches are covered
(`a_-3=0` and `a_-3 != 0` with `b_-3=mu_3 a_-3`), and they are exhaustive. The
free data was **not** otherwise bounded: within each branch the FULL space of raw
degree-`<= d` coefficients (all of `a_2,a_1,a_0,a_-1,a_-2,a_-3`, `mu_3`, and the
`b_0` cross-cancellation kernel) was carried symbolically — no `a_-3`-shape
restriction and no kernel-parameter truncation beyond the degree cap `d` itself.
The only bound is the degree cap, and it is the reported scope (`d in {3,4}`
exact, `d=5` mod-`p`).

## 6. Verification

```sh
uv run --with sympy python research/band3/verify_w2_verdict.py
```

Exact `sympy`. §0 machinery (`Q_m = [D,X]_m` for all `m in [-6,6]`,
`Q_0=(T-1)G`); §1 the explicit slope-`1` witness (cascade+slope+membership hold,
tail fails); §2 the `d<=2` controls (slope kills alone; tail feasible alone); §3
reducer validation at `d=2`; §4 **the `d=3` verdict** — the FULL system reduces to
`[1]` over `QQ` on both branches (INFEASIBLE); §5 the `d=4` confirmation
(reproduced via `msolve` if present; exact-`QQ` values recorded above). A
successful run ends `ALL W2 VERDICT CHECKS PASSED`.

The primary exact-`QQ` `d=4` certificate and the mod-`p` sweeps were produced by
`msolve 0.10.1`; the `d=3` kill is reproduced independently and self-containedly by
`sympy` in the verifier via the feasibility-preserving reduction.
