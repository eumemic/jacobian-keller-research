# Degree-free slope forcing at W2: the a_2(0) factorization of the moment slope

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED — BAND-SCOPED — MIXED PROVED/BOUNDED**

This memo attacks THE new step-1 lemma: degree-free slope-forcing at W2,

```text
   positive cascade  AND  Q_-1=Q_-2=Q_-3=0  AND  membership   =>   R(1)=G(1)=0.
```

If the displayed implication were proved at arbitrary positive-data degree, then
`Q_0=1` forcing `R(1)=1` would close W2 at that scope. It is **not** proved here:
the tail-forcing residual remains open in §6. The bounded `d=3` adjudication
([`slope-forcing-verdict.md`](slope-forcing-verdict.md), commit `968ec09`) established
`R(1) ∈ sqrt(cascade+tail)` at `d=3`, both branches, two engines; on the cascade the
slope collapses to a single monomial. This memo proves only the **degree-free
factorization backbone** and pushes the still-open tail elimination as far as it
currently reaches.

W2 datum, gauge `b_3=0`, quantum band-3 conventions
(`Q_m=sum_(k+l=m)[b_l^[k]a_k - a_k^[l]b_l]`, `f^[n](E)=f(E+n)`,
membership `(E)_j=E(E-1)...(E-j+1) | a_-j,b_-j`):

```text
a_3 = E(E+2)(E+4)   (roots {0,-2,-4}),     b_2 = E(E+3),     b_3 = 0.
branch A:  a_-3=(E)_3 am3,  b_-3=mu_3 a_-3;    branch B:  a_-3=0,  b_-3=(E)_3 C.
```

Exact certificate: [`verify_slope_forcing_degree_free.py`](verify_slope_forcing_degree_free.py)
(ends `ALL SLOPE FORCING DEGREE FREE CHECKS PASSED`; default sympy legs, `HEAVY=1` adds
the `d=4` `msolve` `^` extension and depth trajectory). Every load-bearing upstream fact
— the crossed-product engine, `Q_0=(T-1)G`, the slope gate, the both-ends Lemma P, the
positive cascade — is **re-derived in file**.

## 0. Headline

> **THE FACTORIZATION (PROVED, arbitrary positive-data degree).** On the W2 positive
> cascade the moment slope factors through the **top boundary value `a_2(0)`**:
>
> ```text
>    R(1) = G(1) = a_2(0) · W,                                          (FACTORIZATION)
>    W = b_-2(2) - (2/3) a_-1(1) + (1/4)[(2/3)a_2(1) - b_1(2)] · b_-1(1),
> ```
>
> where `b_1(2), b_-1(1), b_-2(2)` are the cascade-solved boundary `b`-values. In
> particular `a_2(0) | R(1)`, so **`R(1)=0` whenever `a_2(0)=0`**. Consequently the
> degree-free slope-forcing lemma **reduces to a single statement**:
>
> ```text
>    the tail Q_-1,Q_-2,Q_-3 forces  a_2(0) · W = 0.               (RESIDUAL IDENTITY)
> ```

The factorization is two single-point cascade evaluations, both degree-free (they use
only the roots `a_3(0)=b_2(0)=0` and the nonzero pivots `a_3(1)=15`, `b_2(1)=4`,
`b_2(2)=10`), substituted into the both-ends Lemma-P slope formula. It is the exact W2
analogue of the W1 `lambda_0(R)=0` fixed-pivot cascade elimination
([`../band3/quantum-w1-arbitrary-degree.md`](../band3/quantum-w1-arbitrary-degree.md),
commit `e4e704f`) — but here the elimination lands on a **factor of the slope itself**,
not on an annihilator of the filler image.

## 1. The two boundary identities (degree-free)

Evaluate the cascade walls at the top anchor `E=0`. Because `a_3(0)=a_3(-4)=0` and
`b_2(0)=b_2(-3)=0`, almost every term drops and two clean identities survive, valid for
**fully generic coefficients at every degree** (`verify §S2`, checked at `d=2,3,5`):

```text
   Q_4(0) = 10 a_2(0) - 15 b_1(0),
   Q_3(0) =  4 a_1(0) + a_2(0) b_1(2) - a_2(1) b_1(0).
```

`Q_4=0` is `b_1^[3]a_3 - a_3^[1]b_1 + b_2^[2]a_2 - a_2^[2]b_2`; at `E=0` only
`-a_3(1)b_1(0)+b_2(2)a_2(0) = -15 b_1(0)+10 a_2(0)` survives. `Q_3(0)` similarly keeps
the pivot `b_2(1)=4` on `a_1(0)`. On the cascade (`Q_4=Q_3=0`, hence
`Q_4(0)=Q_3(0)=0`) these solve the two top boundary values:

```text
   b_1(0) = (2/3) a_2(0),                                         (a_2(0)-PROPORTIONAL)
   a_1(0) = (a_2(0)/4) [ (2/3) a_2(1) - b_1(2) ].                 (a_2(0)-DIVISIBLE)
```

Both are manifestly divisible by `a_2(0)`. This is the whole of the degree-free content:
`a_2(0)` divides the third Lemma-P term via `b_1(0)`, and divides the first via `a_1(0)`.

## 2. From the boundary identities to the factorization

The both-ends Lemma P ([`w2-joint-theorem.md`](w2-joint-theorem.md) §1;
[`moment-unit-general-k.md`](moment-unit-general-k.md) Lemma P), re-proved in file as an
arbitrary-degree identity with the W2 level-3 drop, gives the filler-independent slope

```text
   R(1) = G(1) = a_1(0) b_-1(1) + a_2(0) b_-2(2) - a_-1(1) b_1(0).
```

Substituting `b_1(0)=(2/3)a_2(0)` and `a_1(0)=(a_2(0)/4)[(2/3)a_2(1)-b_1(2)]`:

```text
   R(1) = a_2(0) · { (1/4)[(2/3)a_2(1) - b_1(2)] b_-1(1) + b_-2(2) - (2/3) a_-1(1) }
        = a_2(0) · W.
```

This identity is machine-checked two ways (`verify §S2, §S3`): as a symbolic identity in
free boundary symbols (`a_1(0), b_1(0)` eliminated via the two walls), and **exactly on
the parametrized cascade** (no `clear_denoms`; `R(1)-a_2(0)W = 0`) at `d=1,2,3` (default),
`d=4` (`HEAVY`), and `d=5` in exploration.

**The mirror.** `R(-1)+R(1)=0` on the cascade (`verify §S3/§S5`, the reflection
`E→-E-1`). So the slope is a single scalar `R(1)=a_2(0)W`, and the tail must annihilate
it.

## 3. What the factorization buys, degree by degree

`W` on the cascade (`verify §S3`):

| `d` | `R(1)` on the cascade | remark |
|---|---|---|
| 1 | `0` | cascade alone forces `R=0` (slope-forcing vacuous) |
| 2 | `0` | cascade alone forces `R(1)=0` (slope-forcing vacuous) |
| 3 | `-(4/9) a_2(0) am1_3` | **first genuine free modulus**; `W=-(4/9)am1_3` |
| 4 | `a_2(0)·W_4`, `W_4` a degree-4 form (13 free coords) | verifier `HEAVY`; `am1_3` (not `am1_4`) + `a_2,a_1,a_0` terms |
| 5 | `a_2(0)·W_5`, `W_5` a degree-5 form (17 free coords) | exploration-only; `am1_3, am1_5` appear, **not** `am1_4` |

(The default verifier checks `R(1)=a_2(0)W` exactly at `d=1,2,3`; `d=4` under `HEAVY`;
the `d=5` row is a recorded exploration observation, not shipped in the verifier — its
`d=5` parametrization build costs ~305 s.)

So the slope-forcing lemma is **vacuous below `d=3`** (the cascade already kills the
slope; there is nothing for the tail to force). The genuine free-modulus slope first
appears at `d=3`, where `W` collapses to `-(4/9)am1_3` — the top coefficient of the
membership quotient of `a_-1`. The clear-denominator normalization of
[`slope-forcing-verdict.md`](slope-forcing-verdict.md) reports this as
`R(1)=-108 a2_0 am1_3`; the exact rational slope is `-(4/9)a_2(0)am1_3` (they differ by
the `clear_denoms` factor `243`; the forcing conclusion is identical). At `d≥4`, `W` is
**not** a single top coefficient — it is a genuine degree-growing polynomial in the free
cascade coordinates. The conjecture "R(1) mod cascade is always a product of an
`a_2`-boundary value and a top coefficient of `a_-1`" is therefore **false at `d≥4`**:
the correct degree-free invariant is the factor `a_2(0)`, and the second factor `W`
grows.

## 4. The tail forcing (bounded; exact scope)

On branch B the depth-3 tail is **filler-linear** (`verify §S5e`): at `d=3` it is
**24 scalar equations, all linear in the 8 fillers** (`a_-2=(E)_2 V`, `b_-3=(E)_3 C`;
4+4 coefficients). The filler map has **full column rank 8**, so its cokernel is
**16-dimensional** — 16 explicit covectors `lam` with `lam·[filler columns]=0`, whose
pairings `lam·N` (N = the non-filler residual) are the elimination conditions on the
positive data. On the variety they cut, `a_2(0)·W = 0`. This is the joint-covector
Fredholm mechanism of [`joint-covector.md`](joint-covector.md), now applied to the
depth-3 tail **without** the `Q_0=1` block.

**Slope forcing, machine-checked (two engines):**

- **`d=3` full tail, sympy EXACT over `QQ` + `GF(65003)` + `GF(32003)`, both branches:**
  `a_2(0)·am1_3 ∈ sqrt(cascade + Q_-1..Q_-5)` — the parametrized-tail Rabinowitsch ideal
  is the unit ideal, a machine-checked Nullstellensatz radical certificate (`verify §S4a`).
- **`d=3` depth-3, `msolve` `^`:** `a_2(0)·am1_3 ∈ sqrt(cascade + Q_-1..Q_-3)` — the
  **graded** refinement (`kmin≤3`; `verify §S4b`). Branch B is fast and reliable (`QQ` +
  `GF(65003)`, ~2 s); branch A's depth-3 `msolve` is memory-flaky (~43 s or OOM), so it
  is attempted over two primes and **downgraded to a note on OOM** — branch A's forcing is
  already certified by the sympy full-tail leg above, and its depth-3 value `kmin=3` is
  the `slope-forcing-verdict.md` `msolve` result. Depth 3 is prohibitively slow in sympy
  (fewer generators → Buchberger takes minutes to reach the unit ideal), so this leg is
  `msolve`; the two engines are complementary.
- **`d=4` FACTORIZATION — confirmed.** `R(1)=a_2(0)·W` holds exactly on the parametrized
  `d=4` cascade (`verify §S3`, `HEAVY`, ~40 s, exact sympy). So the degree-free backbone
  is corroborated at `d=4`.
- **`d=4` TAIL FORCING — not tractable in this session's budget.** The `d=4` parametrized
  Rabinowitsch system (13 free + ~10 filler variables) is beyond both engines here:
  `msolve` `^` is memory-bound (climbs past 6 GB without terminating), and the sympy GB
  does not clear even setup within 2 min (`verify §S6` attempts it and downgrades to a
  note on OOM/timeout, robustly). So `d=4` slope forcing (`a_2(0)·W ∈ sqrt(cascade+tail)`)
  is **not machine-confirmed here** — `d=3` (both branches, `QQ`+two primes) is the
  confirmed tail-forcing cap. The `d`-uniform statement is the §0–§2 factorization, not a
  bounded `msolve` sweep.

**Product vs. `a_2(0)` — the finer mechanism is OPEN.** Whether the tail forces the whole
**product** `a_2(0)·W = 0` (a genuine union `{a_2(0)=0} ∪ {W=0}`) or already forces the
**factor** `a_2(0)=0` is **not resolved here**. The discriminating test — is `a_2(0) ∈
sqrt(cascade+Q_-1..Q_-3)`? — is a non-unit (positive-dimensional) Gröbner computation
that `msolve` does not terminate in budget (it OOMs; `verify §S5d` downgrades it to a
note). The naive filler-elimination (solve 8 fillers, project) is **not** a valid
substitute: it drops the saturation by the `8×8` filler determinant, and indeed reports
even the verified product `a_2(0)·am1_3` as *not* in its (unsaturated) radical — so it
under-approximates the true elimination ideal and cannot settle the union either. What
**is** machine-checked (`verify §S5e`) is the depth-3 tail's Fredholm data — **24
filler-linear equations, filler map full column rank 8, cokernel dimension 16**.

## 5. Depth trajectory

The minimal forcing depth `kmin` (least `k` with `R(1)∈sqrt(cascade+Q_-1..Q_-k)`):

- **`kmin>1` (PROVEN):** an explicit witness on `cascade+Q_-1` with `R(1)≠0` (solve the
  single filler-linear row `Q_-1` at a cascade point with `a_2(0)am1_3≠0`); `verify §S5`
  reproduces `R(1)=-(4/9)·a_2(0)·am1_3 ≠ 0`. So `Q_-1` alone is the **filler Fredholm
  gap** of [`joint-covector.md`](joint-covector.md) §4, not slope-forcing.
- **`kmin>2`:** attempted by the same witness search on `cascade+Q_-1,Q_-2`; when a
  witness is found it proves `kmin>2`, else it remains evidence-only (as in
  `slope-forcing-verdict.md`). Rigorously `1<kmin<=3`; the `d=3` verdict pins `kmin=3`.
- **Depth `3` suffices at `d=3`** (`verify §S4b`): the depth-3 tail is the unit
  Rabinowitsch ideal (`R(1)∈sqrt(cascade+Q_-1..Q_-3)`, `msolve` `^`). Whether the depth
  stays `3` at `d=4` and beyond is **open** here — the `d=4` forcing GB is not tractable
  in budget (§4), so the depth trajectory past `d=3` is not machine-decided.

## 6. The remaining combinatorial identity (open, sharply posed)

The degree-free lemma is now **exactly one statement**, thanks to the factorization:

> **Residual slope-forcing identity (open).** For every positive-data degree `d`, on the
> W2 positive cascade,
>
> ```text
>    Q_-1 = Q_-2 = Q_-3 = 0   =>   a_2(0) · W = 0,
> ```
>
> where `W = b_-2(2) - (2/3)a_-1(1) + (1/4)[(2/3)a_2(1)-b_1(2)]·b_-1(1)` is the explicit
> residual of §0. Equivalently `a_2(0)·W ∈ sqrt(cascade + Q_-1..Q_-3)` at all `d`. (Whether
> the tail forces the whole product or already the factor `a_2(0)` — §4 — is itself open;
> either way `R(1)=a_2(0)W=0`.)

This is a Fredholm/cokernel statement: the 16-covector elimination of §4 must pair the
non-filler residual `N` to a multiple of `a_2(0)·W` at arbitrary degree — the
degree-free analogue of the bounded joint-covector kill (committed at `d=3`; optional/heavy at `d=4`)
([`joint-covector.md`](joint-covector.md) §4), now on the depth-3 tail cokernel and with
`W` in place of the moment-unit `rhs`. What blocks it from being degree-free today is the
same gap as the joint-covector lemma: the covectors' fixed part lives on the `(a_3,b_2)`
root necklace (degree-free, adjoint criterion,
[`lambda-general-k.md`](lambda-general-k.md) Thm A′), but the datum-dependent part sits
on the `(a_2,b_1)` necklace, whose nodes are the **algebraic** roots of the solved
`a_2,b_1`; a degree-free description of that block would close the lemma. `W` growing with
`d` (§3) is exactly this datum-dependence.

## 7. Evidence ledger — proved / bounded / refuted / open

**Proved (arbitrary coefficient degree, char 0, machine-checked identities):**
- `Q_m=[D,X]_m` for `m∈[-6,6]`; `Q_0=(T-1)G` (`§S0`).
- Slope gate `D|(E-R) ⇔ R(1)=1,R(-1)=-1`; both-ends Lemma P
  `R(1)=a_1(0)b_-1(1)+a_2(0)b_-2(2)-a_-1(1)b_1(0)`, filler-independent (`§S1`).
- **The two boundary identities** `Q_4(0)=10a_2(0)-15b_1(0)`,
  `Q_3(0)=4a_1(0)+a_2(0)b_1(2)-a_2(1)b_1(0)`, for fully generic coefficients (`§S2`).
- **THE FACTORIZATION** `R(1)=a_2(0)·W` on the cascade, degree-free; hence `a_2(0)|R(1)`
  and `R(1)=0` whenever `a_2(0)=0` (`§S2`). This is the new degree-free structural
  theorem — the slope is divisible by the top boundary value.

**Bounded-finite evidence (exact scope stated):**
- `R(1)=0` on the cascade alone for `d≤2`; `W=-(4/9)am1_3` at `d=3` (so
  `R(1)=-(4/9)a_2(0)am1_3`); `R(1)=a_2(0)·W` verified exactly on the parametrized cascade
  at `d=1,2,3` (default) and `d=4` (`HEAVY`), plus `d=5` in exploration; `R(1)` a genuine
  free modulus on the cascade (attains `1`) (`§S3`).
- Slope forcing `a_2(0)·W ∈ sqrt(cascade+tail)` at **`d=3`, both branches**: full tail
  over `QQ`+`GF(65003)`+`GF(32003)` (sympy exact Rabinowitsch radical certificate,
  `§S4a`); the depth-3 graded refinement over `QQ`+`GF(65003)` (branch B) / two primes
  (branch A) via `msolve` `^` (`§S4b`). **`d=4` tail forcing is NOT machine-confirmed
  here** — the larger `d=4` Rabinowitsch GB is beyond both engines in budget (`msolve`
  memory-bound >6 GB; sympy GB slow); `§S6` attempts it and downgrades robustly. (The
  `d=4` *factorization* `R(1)=a_2(0)W` **is** confirmed, above.)
- Depth-3 tail is **24 filler-linear equations, filler map rank 8, 16-covector cokernel**
  (machine-checked, `§S5e`); the certificate artifact
  [`w2_slope_forcing_df_cert.txt`](w2_slope_forcing_df_cert.txt) records the collapse
  certificate and the elimination path (12 conditions of total degree 10–12 on the 9 free
  coordinates, after clearing the `8×8` filler determinant). Depth: rigorous `1<kmin<=3`
  (`§S5a–c`; `kmin>1` proven by an explicit witness `R(1)=-16/9≠0`), pinned to `kmin=3` by
  `slope-forcing-verdict.md`.

**Refuted / corrected:**
- The conjecture "R(1) mod cascade is always a product of an `a_2`-boundary value and a
  top coefficient of `a_-1`" — **true at `d=3` only**; false at `d≥4`, where `W` is a
  degree-growing polynomial (not a single `a_-1` coefficient) (`§S3`). The invariant
  factor is `a_2(0)`; the second factor grows.
- Clarification (not a refutation): the exact rational slope at `d=3` is
  `-(4/9)a_2(0)am1_3`; the `-108 a2_0 am1_3` of `slope-forcing-verdict.md` is the
  `clear_denoms`-normalized numerator (factor `243`).

**Open / not claimed:**
- The **residual slope-forcing identity** `a_2(0)·W ∈ sqrt(cascade+Q_-1..Q_-3)` at
  **arbitrary degree** (§6) — machine-confirmed at `d=3` only (the `d=4` forcing GB is not
  tractable in budget). The degree-free covector for the `(a_2,b_1)`-necklace block is not
  obtained; `W`'s datum-dependence is the obstacle.
- Whether the tail forces the **product** `a_2(0)·W=0` (a genuine union) or already the
  **factor** `a_2(0)=0`: the discriminating test `a_2(0) ∈ sqrt(cascade+Q_-1..Q_-3)?` is a
  non-unit positive-dimensional Gröbner basis that `msolve` OOMs on; the naive
  filler-elimination is invalid (drops the `8×8`-determinant saturation). Undetermined.
- Whether `kmin=3` for `d>4`; the explicit minimal Nullstellensatz power `m` (only
  `m≥1`, radical membership, is certified — the elimination Gröbner bases do not
  terminate in a sane sympy budget).
- No Weyl pair, no DC1 counterexample. All of Band 3, DC1, JC2 remain open. The
  infinite-dimensional `Im L_K ∩ Im L_H`
  ([`two-filler-cross-cancellation.md`](two-filler-cross-cancellation.md)) is untouched.

## 8. Verification

```sh
uv run --with sympy python research/dc1-program/verify_slope_forcing_degree_free.py
# d=4 msolve '^' extension + depth trajectory (needs msolve 0.10.1 on PATH):
HEAVY=1 uv run --with sympy python research/dc1-program/verify_slope_forcing_degree_free.py
```

`§S0` engine; `§S1` slope gate + both-ends Lemma P + filler-independence; `§S2` the two
boundary identities + THE FACTORIZATION `R(1)=a_2(0)W`; `§S3` factorization on the
parametrized cascade (`d=1,2,3`; `d=4` `HEAVY`), `W` collapse, slope-is-free control;
`§S4a` slope forcing at `d=3` on the full tail both branches (`QQ`+two primes, sympy
exact); `§S4b` the depth-3 graded refinement (`msolve` `^`, `QQ`+one prime, both
branches); `§S5` depth trajectory (witnesses; rigorous `1<kmin<=3`), the `a_2(0)`-not-forced
union probe (`msolve`, budget-permitting), and `§S5e` the depth-3 Fredholm structure
(24 filler-linear eqs, rank 8, cokernel 16); `§S6` (`HEAVY`) the `d=4` `msolve` `^`
attempt + depth trajectory (memory-bound). Ends
`ALL SLOPE FORCING DEGREE FREE CHECKS PASSED`.
