# W2 bounded verdict: no encoded slope+tail solution at raw cap d=3

**INDEPENDENTLY DERIVED — EXACT ALGEBRA — NOT PEER REVIEWED — FIXED ORIENTATION/GAUGE — BOUNDED RAW CAP**

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
pair** and trigger the escalation protocol. A "no" establishes only that no
encoded solution exists at the tested raw cap(s).

## Verdict

> **BOUNDED VERDICT.** In the fixed normalized W2 orientation/gauge, the encoded
> combined system is the unit ideal at raw cap `d=3`, on both tail branches, by
> exact `QQ` polynomial identities checked directly against the reconstructed FULL
> systems. An exact `d=4` result was reported
> externally and can be checked by the verifier only when optional `msolve 0.10.1`
> is installed; without committed run artifacts or a completed optional run, it is
> documentary rather than part of the self-contained verdict.
>
> Thus no encoded slope-`1` datum clears the negative tail at raw cap `d=3`.
> Separate feasibility plus the empty conjunction is established only to the
> extent stated in §3; this is not an arbitrary-degree joint theorem. No DC1
> candidate pair materializes, and no Weyl pair or counterexample is constructed.
>
> **Scope.** The self-contained exact rational result is bounded to raw cap `d=3`.
> The next unverified exact rational raw cap is `d=4`; `d=5` modulo `65003` is a
> reported finite-field computation, not a proof over `QQ`. Nothing here closes
> W2 at arbitrary degree or any other Band-3 top.

Committed certificate: [`verify_w2_verdict.py`](verify_w2_verdict.py), with the
exact rational multiplier witness in
[`w2_d3_qq_certificates.txt`](w2_d3_qq_certificates.txt). Its final status states
exactly what ran, says only that the **executed** checks passed, and lists every
optional or documentary item as `SKIP`. Base commit `66156f2`.

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

## 1. Reported external mod-`p` sweep (documentary, not reproduced here)

The original computation report says the full combined system was reduced modulo
five primes `p in {10007, 31013, 50021, 65003, 99013}` and run through `msolve`
over `GF(p)`:

| system (`d=3`) | branch A | branch B |
|---|---|---|
| `FULL` (cascade + `Q_0=1` + tail + membership) | **`[-1]` (UNIT) mod all 5** | **`[-1]` (UNIT) mod all 5** |

`msolve`'s `[-1]` denotes no solution in the algebraic closure, hence the unit
ideal over that finite field. These timings and outputs are an **external
computation report**: no per-prime input/output artifacts are committed, and the
verifier does not execute this five-prime sweep. Reproduction requires generating
the relevant `FULL` system for each branch/cap in msolve format (variable line,
characteristic line, comma-separated equations), then running
`msolve -f INPUT.ms -o OUTPUT.out` and retaining inputs, tool/version metadata,
and outputs beginning `[-1]`. Do not treat this report as the exact-`QQ` proof.

## 2. Exact `QQ` verdict: committed d=3 check and external cross-checks

**Reported external `msolve` over `QQ` (char 0).** The original report records a
reduced Gröbner basis (`-g 2`) of `[1]` for the `FULL` `d=3` system on both
branches, with default output `[-1]`. No raw msolve input/output is committed, so
this is documentary corroboration rather than a second committed certificate.

**Committed exact rational identity.** `sympy` cannot Gröbner the raw
`62`-eq/`26`-var system in a sane budget (the very reason this task existed).
Instead, [`w2_d3_qq_certificates.txt`](w2_d3_qq_certificates.txt) commits sparse
Nullstellensatz multipliers `h_i` for the freshly reconstructed `FULL` equations
`f_i`. The principal verifier renames the free parameters deterministically to
`y0,...`, parses every multiplier coefficient in `QQ`, and checks

```text
sum_i h_i f_i = 1
```

by exact multivariate coefficient collection. Branch A uses 17 nonzero
multipliers among the 62 equations in `QQ[y0,...,y25]`; branch B uses 14 among 62
in `QQ[y0,...,y24]`. This identity directly implies `1` belongs to each ideal.
The witnesses were generated with Singular, but the generator is not required or
trusted by the verifier: only the independently checked identities matter.

The feasibility-preserving linear reducer remains as a `d=2` control and is
validated there against a **direct** full-system Gröbner
(`reduced-unit == direct-unit`). The committed, self-contained exact result is
therefore: **exact polynomial identities over `QQ`, d=3, both branches: UNIT**.
The reported msolve result is a solver cross-check without committed artifacts.

## 3. The controls — this is not a false kill

The committed controls distinguish several feasible and infeasible systems, but
they do **not** establish every entry of the former d=3-by-branch table. The exact
committed evidence is:

| control | committed evidence | scope |
|---|---|---|
| cascade + tail, no `Q_0=1` | proper ideal | `d=1`, both branches in this verifier; `d=2` in `verify_w2_tail.py` |
| cascade + `Q_0=1`, no tail | explicit exact witness | `d=3`, branch B (`a_-3=0`, free `b_-3`) |
| `FULL` conjunction | unit ideal | `d=3`, both branches |

Thus separate feasibility is established, and the conjunction is empty, but a
claim of separate d=3 feasibility on **both** branches would exceed the committed
controls. Branch A cascade+slope feasibility is not needed for the bounded unit
verdict and is not asserted here.

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
respectively). Because that datum's **raw free-polynomial degrees** fit cap `d` for every
`d>=3`, it witnesses "cascade + `Q_0=1`" feasible in branch B at those caps.
The cap does not bound every resulting coefficient polynomial: membership factors
and the solved cascade formulas raise degrees (for example `a_-1` and solved
`b` coefficients can exceed `d`). This is the point [`w2-decisive.md`](w2-decisive.md) §5b already found not to
extend; the present memo upgrades "this particular point fails" to "**the entire
slope-`1` fiber fails**, exactly, at `d=3`". The analogous `d=4` conclusion remains
an external report unless the optional computation is actually run.

## 4. Optional d=4 reproduction and reported d=5 evidence

The original external report records:

| system | branch A | branch B | publication status |
|---|---|---|---|
| `FULL`, `d=4`, `QQ` (msolve, exact) | **UNIT** | **UNIT** | optionally reproducible by the verifier when `msolve` is on `PATH` |
| `FULL`, `d=4`, mod 5 primes | `[-1]` all 5 | `[-1]` all 5 | documentary; no artifacts committed |
| `FULL`, `d=5`, mod `65003` | `[-1]` | `[-1]` | documentary finite-field evidence only |

`d=4` has 73 equations in 32 (branch A) / 31 (branch B) free parameters. With
`msolve` installed, the committed verifier generates the `d=4` systems and runs
exact char 0 and char 65003 via `msolve -f INPUT -o OUTPUT`; absent `msolve`, both
are explicit `SKIP`. It does not run the five-prime sweeps or `d=5`.

The self-contained exact rational ledger stops at `d=3`; `d=4` is the next
unverified exact rational raw cap unless the optional msolve run completes. The
reported char-65003 `d=5` result is not a `QQ` proof.

## 5. What this settles, and what it does not

**Reading with the siblings.** [`w2-theory.md`](w2-theory.md) showed that reflection does not close W2 and that
the complete dual of the W2 filler quotient is the point-annihilator span
`{ev_-1,ev_0,ev_1}`; no additional functional **in that dual** remains. It did not
rule out every conceivable linear invariant of larger systems. Together with the
achievable moment slope, this left the encoded combined feasibility open.
[`w2-negative-tail.md`](w2-negative-tail.md) proved the tail is a *proper*
(feasible) ideal on its own at `d=1,2`. [`w2-decisive.md`](w2-decisive.md) proved
`R(1)=1` achievable at `d=3`. **This memo closes the self-contained loop at `d=3`:
the simultaneous system is the unit ideal on both branches.** The analogous
`d=4` result is externally reported and only optional here.

**The DC1 consequence.** No encoded candidate pair exists at raw cap `d=3` in
this fixed normalized W2 orientation/gauge, so `[D,X]=1`, the char-`p` sieve, and
the `A_1`-generation test were not run. This says only that no encoded solution
exists at that cap; it does not kill a construction route at arbitrary degree or
address other Band-3 tops.

**Honest ledger — proved vs bounded.**

- **Committed exact `QQ` certificate:** the W2 combined system (positive cascade
  `∧` `Q_0=1` `∧` `Q_-1..Q_-5=0` `∧` membership) is the unit ideal at raw cap
  `d=3`, on both tail branches. The sparse identities are checked directly against
  all 62 freshly reconstructed FULL equations; no d=3 reduction is load-bearing.
- **Exact `QQ` result with optional reproduction:** the same unit result at raw cap
  `d=4`, both branches, is reported from msolve and reruns when optional msolve is
  available; otherwise the verifier marks it `SKIP`.
- **Reported external finite-field evidence:** `d=5` at prime `65003`, and `d=3,4`
  at five listed primes. No machine-checkable artifacts for these sweeps are
  committed; they remain documentary and are not `QQ` proofs.
- **Controls:** cascade+tail is committed feasible at `d=1` in both branches here
  and at `d=2` in the sibling tail verifier. An explicit branch-B slope-`1`
  witness at raw cap `d=3` satisfies cascade+slope+membership and fails the tail.

**Open / NOT claimed (explicit).**
- **Arbitrary degree.** Whether some slope-`1` datum at a higher raw cap clears
  the tail is **not** decided here. The next unverified exact rational raw cap is
  `d=4`; the reported `d=5` modular result is not a `QQ` proof. A finite Gröbner
  search cannot settle arbitrary degree; closing W2 at
  arbitrary degree now requires a **structural (degree-free)** argument — the
  natural next target, analogous to the `lambda_r` theorem that closed
  `r != -4` at arbitrary degree. The infinite-dimensional
  `Im L_K ∩ Im L_H` ([`../dc1-program/two-filler-cross-cancellation.md`](../dc1-program/two-filler-cross-cancellation.md))
  is untouched.
- A canonical common-root hatch family exists for every `k>=3`. Band 4 has a
  checked analogous slope gate, but uniqueness of one hatch per band, reduction
  of all hatches to one gate, and any uniform higher-band induction role are
  conjectural. No higher-band combined slope+tail system is addressed here.
- Everything outside the W2 combined system: the non-AP `deg a_3 >= 6` exotic
  tops, all of band `>= 4`, DC1, JC2. No Weyl pair, no counterexample.

**Exceptional-locus / search-space honesty.** Both `Q_-6` branches are covered
(`a_-3=0` and `a_-3 != 0` with `b_-3=mu_3 a_-3`), and they are exhaustive. The
free data was not otherwise specialized: within each branch, the raw free
polynomials for `a_2,a_1,a_0` and the membership quotients of `a_-1,a_-2,a_-3`
(or `b_-3` in branch B), together with `mu_3` and retained solver kernels, were
carried symbolically up to raw cap `d`. The cap applies to those raw free
polynomials **after factoring out required membership factors**. Consequently
`a_-j` and cascade-solved `b` coefficients can have degrees greater than `d`;
phrases such as "all coefficient degrees <= d" are incorrect. No additional
`a_-3` shape restriction or kernel specialization was imposed.

## 6. Verification

```sh
uv run --with sympy python research/band3/verify_w2_verdict.py
```

Exact SymPy executes §0 machinery, §1 the explicit branch-B slope witness, §2
low-cap controls, §3 reducer validation, and §4 the FULL `d=3` unit test over `QQ`
on both branches. Section 5 runs `d=4` exact `QQ` and mod `65003` only when
`msolve` is on `PATH`; otherwise both are explicit `SKIP`. The five-prime sweeps
and `d=5` report are always listed as documentary `SKIP` because no committed
artifact path is available. The final line is
`W2 VERDICT EXECUTED CHECKS PASSED; OPTIONAL/DOCUMENTARY ITEMS SKIPPED ABOVE`,
preceded by exact PASS/SKIP counts.

The committed self-contained exact certificate consists of direct `QQ` identities
against the reconstructed FULL `d=3` systems.
The exact-`QQ` `d=4` result is optionally reproducible with msolve. The five-prime
and `d=5` statements remain external reports until inputs, outputs, versions, and
commands are committed.
