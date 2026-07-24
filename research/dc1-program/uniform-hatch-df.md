# Uniform Hatch Slope-Forcing: the degree-free backbone (both cokernel pairings)

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED — BAND-SCOPED — MIXED PROVED/BOUNDED**

This memo runs architecture's **unification** target: at the band-`k` hatch
(`k>=4`), does the negative tail (no `Q_0`) force **both** cokernel pairings of the
residual `R` to vanish — `R(1)=0` **and** `lambda_3(R)=0` — at arbitrary degree? A
single such mechanism would collapse steps **1 (slope-forcing)** and **2 (covector)**
into one for the whole tower. The bounded band-4 `d=3` confirmation already exists
([`hatch-slope-forcing.md`](hatch-slope-forcing.md)). The question here is the
**degree-free structure**: what is the band-4 analogue of the W2 `a_2(0)`
factorization ([`slope-forcing-degree-free.md`](slope-forcing-degree-free.md)), and
does the depth-2 kill become a covector statement at *every* degree?

Conventions frozen from the band-`k` corpus (gauge `b_k=0`): `(x^a f)(x^b g)=x^{a+b}
f(E+b)g(E)`, `f^[n](E)=f(E+n)`, `Q_m=sum_{i+l=m}(b_l^[i]a_i - a_i^[l]b_l)`,
`u:=b_{k-1}`, membership `(E)_j=E(E-1)...(E-j+1) | a_-j,b_-j`, `Q_0=(T-1)G`, `G(0)=0`,
`R(1)=Q_0(0)=G(1)`, `lambda_3(R)=lambda_3(G)`. Hatch: `a_k=prod_{i<k}(E+i(k-1))`,
`u=prod_{j<k-1}(E+jk)`. Band 4: `a_4=E(E+3)(E+6)(E+9)`, `u=E(E+4)(E+8)`,
`lambda_3=ev_-2-ev_-4+ev_-5`, `lambda_3(E)=-3`. **Phi fillers** `= a_{-(k-1)}, b_{-k}`.

Exact certificate: [`verify_uniform_hatch_df.py`](verify_uniform_hatch_df.py). The
default run (`~165 s`) prints a PASS/SKIP evidence ledger and ends
`ALL UNIFORM HATCH DF CHECKS PASSED` only when the band-4 `d=3` covector-kill payload
actually ran on `msolve` and returned UNIT; without `msolve` that payload is `SKIP`
and the final line says so. Heavy corroboration (the determinant-saturated covector
certificate, branch A, second prime, `QQ` char-0, band-5 `d=2`) is behind `HEAVY=1`
(note: the full `HEAVY` run takes ~12 minutes — it does not fit a single 600 s budget).
Every load-bearing upstream fact — the crossed-product engine, `Q_0=(T-1)G`,
`lambda_3` and `Im Phi subset ker lambda_3`, the positive cascade — is **re-derived
in file**.

## 0. Headline

> **The W2 single-factor slope factorization does NOT survive to band 4; the
> degree-free backbone is instead (i) a uniform both-ends Lemma-P and (ii) a uniform
> top-wall proportionality, and (iii) the depth-2 tail is filler-linear at every
> degree — so the depth-2 kill of BOTH targets is a consistency-covector statement at
> every degree, machine-certified at band 4 `d=3` (`kmin=2`).**
>
> 1. **BACKBONE-1 — both-ends Lemma-P (PROVED, arbitrary degree, `k=3,4,5`).** On the
>    band-`k` datum,
>    ```text
>        R(1) = Q_0(0) = G(1) = sum_{i=1}^{k-1} [ a_i(0) b_-i(i) - a_-i(i) b_i(0) ],
>    ```
>    filler-independent. (Band 4: `a_1(0)b_-1(1) + a_2(0)b_-2(2) + a_3(0)b_-3(3) -
>    a_-1(1)b_1(0) - a_-2(2)b_2(0)` — the exact `k=4` image of the W2 three-term form.)
>    An analogous point-functional boundary form holds for `lambda_3(R)=lambda_3(G)`.
> 2. **BACKBONE-2 — top-wall proportionality (PROVED, arbitrary degree, `k=3,4,5`).**
>    At the anchor, the top wall keeps only two boundary values:
>    `Q_{2k-2}(0) = c_a·a_{k-1}(0) - c_b·b_{k-2}(0)`, so on the cascade
>    ```text
>        b_{k-2}(0) = rho_k · a_{k-1}(0),   rho_3=2/3,  rho_4=21/80,  rho_5=8/55.
>    ```
>    Hence two Lemma-P terms of **both** targets — the explicit `a_{k-1}(0)·b_-(k-1)(k-1)`
>    and `-a_-(k-2)(k-2)·b_{k-2}(0)=-rho_k a_-(k-2)(k-2)·a_{k-1}(0)` — carry a **shared**
>    factor `a_{k-1}(0)`. This is the *partial* backbone.
> 3. **REFUTED — no single divisor (band 4).** Unlike W2, where `R(1)=a_2(0)·W` factors
>    through a **single** boundary value, at band 4 **no** `a_j(0)` divides `R(1)` or
>    `lambda_3(R)` modulo the cascade (`d=3`, exact `GF(p)` radical membership: each is
>    UNIT-**false**), and no tested boundary value is a common divisor. The
>    remaining Lemma-P terms carry `a_2(0)` (the *next* boundary value), so the slope is
>    a genuine **two-boundary** object; W2's clean factorization is special to `k=3`.
> 4. **BACKBONE-3 — depth-2 filler-linearity (PROVED, arbitrary degree, band 4 & 5).**
>    The depth-2 tail `Q_-1, Q_-2` is **linear** in the Phi fillers `a_{-(k-1)}, b_{-k}`
>    (max total degree `1`; band 4 `d=2,3,4`, band 5 `d=3`), and both targets are
>    Phi-filler-independent. **Therefore the depth-2 kill is, at every degree, the
>    consistency (cokernel) condition of a linear filler system** — a covector
>    statement, not a factorization.
> 5. **THE COVECTOR KILL (BOUNDED, band 4 `d=3`, `kmin=2`).** Both targets lie in
>    `sqrt(cascade + Q_-1 + Q_-2)` (`msolve '^'` UNIT), while `Q_-1` alone forces
>    **neither** (explicit feasible witnesses). The depth-2 filler map is `18×8`, full
>    column rank `8`, **cokernel dimension 10**. Under `HEAVY`, both targets are
>    certified in the **determinant-saturated** consistency-covector ideal
>    `sqrt((cascade + 8 explicit covectors) : det^inf)` — the valid form; naive
>    filler-elimination (dropping the saturation) is invalid.
> 6. **Band 5 is DEGENERATE through `d=3`.** Via the reduce-first route (tractable where
>    the raw `msolve` probe was not), both targets reduce to `0` on the band-5 cascade
>    **alone** at `d=3` — so the tail is not tested there. The genuine non-degenerate
>    band-5 probe is `d>=4`, which is **out of budget** here.

This is a **structural + bounded** result. It proves the degree-free *scaffold* of the
uniform statement (Lemma-P, top wall, filler-linearity) and machine-certifies the
depth-2 both-target kill at band 4 `d=3`; it does **not** prove the arbitrary-degree
kill, and it builds no Weyl pair and settles neither DC1 nor JC2.

## 1. The backbone, derived

### 1.1 Both-ends Lemma-P (BACKBONE-1)

`Q_0=(T-1)G` with `G(0)=0` gives `R(1)=Q_0(0)=G(1)`. Expanding `G(1)` on the band-`k`
datum, every interior/`Phi`-filler contribution cancels and the surviving terms are
the `k-1` **boundary pairings** displayed in §0.1. The verifier checks the identity
against the raw `Q_0(0)` for fully generic coefficients at `k=3,4,5` and `d=2,3`
(`verify §S2`). This is the `k`-uniform generalization of the W2 formula
`R(1)=a_1(0)b_-1(1)+a_2(0)b_-2(2)-a_-1(1)b_1(0)`
([`slope-forcing-degree-free.md`](slope-forcing-degree-free.md) §2).

### 1.2 Top-wall proportionality (BACKBONE-2)

At the anchor `E=0`, `a_k(0)=u(0)=0` kills all but two terms of the top wall:
`Q_{2k-2}(0)=c_a a_{k-1}(0)-c_b b_{k-2}(0)` (an exact identity, `verify §S2`). On the
cascade (`Q_{2k-2}=0`) this pins `b_{k-2}(0)=rho_k a_{k-1}(0)` with
`rho_3=2/3, rho_4=21/80, rho_5=8/55`. This is the direct analogue of W2's
`b_1(0)=(2/3)a_2(0)`.

### 1.3 Why the single factor dies at band 4 (REFUTED), and what survives (PARTIAL)

At W2 there are only two rungs above the slope, so both boundary values `b_1(0),
a_1(0)` are forced `a_2(0)`-proportional and `R(1)=a_2(0)·W` factors through **one**
value. At band 4 the chain is longer: the wall pins `b_2(0)=rho_4 a_3(0)`, but the
next rung `Q_5(0)=0` pins `b_1(0)` to a combination `[120 a_2(0)+ (a_3(0)-divisible)]/
280` that **carries `a_2(0)`**, and the Lemma-P term `a_2(0)b_-2(2)` is `a_2(0)` too.
So exactly the two terms `a_3(0)b_-3(3)` and `-a_-2(2)b_2(0)=-rho_4 a_-2(2)a_3(0)` are
`a_3(0)`-divisible (shared by `lambda_3(R)`, `verify §S2`), but the whole slope is
**not**: the exact `GF(p)` radical test finds `a_3(0), a_2(0), a_1(0)` each a
**non-divisor** of both targets modulo the cascade (`d=3`; `verify §S2`, `a_3(0)`
default, `a_2(0),a_1(0)` under `HEAVY`), with no claim about an arbitrary common factor of the two targets.
**The degree-free invariant is therefore not a factor but the covector structure of
§1.4.**

### 1.4 Depth-2 filler-linearity ⇒ covector statement at every degree (BACKBONE-3)

The tail is bilinear in `(a,b)`; the positive cascade solves `b_{k-2}..b_{-(k-1)}`
linearly in the `a`-data, and the two **Phi fillers** `a_{-(k-1)}, b_{-k}` each enter
every `Q_m` term **once**. Hence `Q_-1, Q_-2` are **linear** in the Phi fillers (max
total degree `1`; machine-checked band 4 `d=2,3,4`, band 5 `d=3`, `verify §S3`), with
coefficients polynomial in the positive data, and the targets are Phi-independent.
Writing the depth-2 rows as `M(y)·f + N(y)=0` (`f`= Phi fillers, `y`= positive data),
**filler-solvability is exactly the vanishing of the cokernel pairings `lam·N`
(`lam` a left-null covector of `M`)** — a set of polynomial conditions on `y`, at
every degree. This is the sense in which the depth-2 kill "is a covector statement at
every degree": the linear structure is degree-uniform; only the *membership of the
targets in the cokernel ideal* is checked per degree (§2).

## 2. The band-4 `d=3` covector kill (bounded)

- **Filler map.** `18` scalar depth-2 rows, `8` Phi fillers, **full column rank 8**
  at a generic point ⇒ **cokernel dimension 10** (`verify §S4`). Both targets are
  free (nonconstant) and **distinct** moduli on the cascade (GB normal forms), and the
  origin lies on `cascade+tail` — the forcing is non-vacuous.
- **The kill (`kmin=2`).** `R(1), lambda_3(R) ∈ sqrt(cascade + Q_-1 + Q_-2)`
  (`msolve '^'` UNIT mod `65003`; branch B default, branch A + mod `32003` + `QQ`
  char-0 under `HEAVY`). `Q_-1` **alone** forces **neither** (explicit feasible
  witnesses with target `=1`, `msolve` GB non-unit on a `0`-dim slice). So `kmin=2`
  for both targets — the depth-2 tail is the minimal forcing depth, and it forces
  **both** cokernel pairings at once.
- **Explicit determinant-saturated covectors (`HEAVY`).** Choosing `8` pivot rows
  (`det = D`), the non-pivot rows give explicit consistency polynomials
  `C_i = D·N_i - M_i·adj(M_I)·N_I`. **The executed certificate constructs and
  certifies with `8` explicit covectors (audit-corrected — the verifier prints "8
  pivot rows, 8 explicit consistency covectors"; the cokernel has dimension 10, so
  8 suffice for the UNIT verdict but do not enumerate the full cokernel).** Then
  ```text
      R(1), lambda_3(R) ∈ sqrt( (cascade + {C_i}) : D^inf )       (msolve '^' UNIT),
  ```
  the **valid** determinant-saturated statement. The naive alternative — solve the
  fillers and substitute — drops the `D`-saturation and under-reports (the W2
  cautionary tale, [`slope-forcing-degree-free.md`](slope-forcing-degree-free.md) §4),
  so keeping the fillers as Rabinowitsch variables (or the explicit `:D^inf`
  saturation) is required.

## 3. Band 5 and general `k`

**Band 5 (reduce-first route).** Forward-substituting the cascade linearly and then
GB-reducing over `GF(p)` makes the band-5 `d=3` verdict tractable (the raw `msolve`
probe was out of budget). The verdict is **degenerate**: both targets reduce to `0`
on the band-5 cascade **alone** at `d=3` (`verify §S6`), extending the known `d<=2`
degeneracy to `d<=3` (consistent with the census's "band-5 `R(1)` forced `0` at
`d<=3`", [`hatch-census.md`](hatch-census.md) §3.1). So the band-5 tail plays no role
at `d=3`; the first non-degenerate band-5 degree is `d>=4`, whose cascade build alone
is `~50 s` and whose reduce+GB exceeds the budget — **open / out of budget**.

**Uniform structural facts (machine witnesses `k=4,5`; `k`-scope per tier notes).**
- Both-ends Lemma-P and top-wall proportionality: `k=3,4,5` (`verify §S2`).
- Depth-2 filler-linearity in the Phi fillers: band 4 (`d=2,3,4`) and band 5 (`d=3`)
  (`verify §S3`).
- `lambda_3` a point functional with `lambda_3(E)=-(k-1)(k-2)/2`, `Im Phi ⊂ ker
  lambda_3` degree-free: re-derived `k=3,4,5` (`verify §S1`); arbitrary-`k` **cited**
  from [`nonpoint-covector.md`](nonpoint-covector.md), not re-derived symbolically in
  `k` here.

## 4. Evidence ledger — proved / bounded / refuted / open

**Proved (arbitrary coefficient degree, machine identities; `k in {3,4,5}` scope):**
- Engine `Q_m=[D,X]_m`, `Q_0=(T-1)G`, `G(0)=0` (`k=3,4,5`).
- `lambda_3` point functional (moving-sum), `lambda_3(E)=-(k-1)(k-2)/2`,
  `lambda_3(1)=1`, `gcd(S_k,S_{k-1})=1`, `Im Phi ⊂ ker lambda_3` in generic degree-`k`
  `C,V`; both targets Phi-filler-independent. **Tier note:** the in-file witnesses are
  fixed low degree; the degree-generality rests on the structural arguments
  (bilinearity; `a_k(0)=u(0)=0`) as in the corpus. Arbitrary-`k` `lambda_3` is CITED,
  not re-derived in `k`.
- **BACKBONE-1** both-ends Lemma-P `R(1)=sum_{i=1}^{k-1}[a_i(0)b_-i(i)-a_-i(i)b_i(0)]`
  (generic coeffs, `k=3,4,5`, `d=2,3`).
- **BACKBONE-2** top-wall identity `Q_{2k-2}(0)=c_a a_{k-1}(0)-c_b b_{k-2}(0)` and the
  proportionality `b_{k-2}(0)=rho_k a_{k-1}(0)` (`k=3,4,5`); the two shared Lemma-P
  terms of both targets are `a_{k-1}(0)`-divisible on the wall.
- **BACKBONE-3** depth-2 filler-linearity: `Q_-1, Q_-2` linear in the Phi fillers,
  targets Phi-independent (band 4 `d=2,3,4`, band 5 `d=3`) — so the depth-2 kill is a
  consistency-covector statement at every degree.

**Bounded-finite evidence (exact scope stated):**
- **Band 4, `d=3` covector kill (`kmin=2`):** `R(1), lambda_3(R) ∈ sqrt(cascade +
  Q_-1 + Q_-2)` — `msolve '^'` UNIT mod `65003` (branch B default; branch A + mod
  `32003` + `QQ` under `HEAVY`); `Q_-1` alone forces neither (explicit witnesses).
- **Band 4, `d=3` cokernel:** `18×8` depth-2 filler map, full column rank `8`,
  cokernel dim `10`; targets free + distinct on the cascade; origin feasible.
- **Band 4, `d=3` determinant-saturated covector certificate (`HEAVY`):** both targets
  `∈ sqrt((cascade + 8 explicit covectors) : det^inf)` — `msolve '^'` UNIT.
- **Band 5, `d<=3` (reduce-first):** both targets reduce to `0` on the cascade alone
  (degenerate); band-5 `d=2` `cascade+tail` also `0` (`HEAVY`).

**Refuted / corrected:**
- The naive generalization "the band-4 slope factors through a single boundary value
  `a_{k-1}(0)`" is **false**: no single `a_j(0)` divides `R(1)` or `lambda_3(R)` mod
  the band-4 cascade (`d=3`, exact `GF(p)`). No arbitrary gcd/common-factor test is claimed. W2's
  `R(1)=a_2(0)·W` is special to `k=3`; for `k>=4` the correct invariant is the
  covector structure, not a factor.

**Open / not claimed:**
- The **arbitrary-degree** kill (both targets in `sqrt(cascade+tail)` at all `d`): the
  covector membership is machine-certified only at band 4 `d=3`. What remains is a
  degree-free description of the `M(y)·f+N(y)` cokernel ideal (the same datum-dependent
  necklace obstruction as [`slope-forcing-degree-free.md`](slope-forcing-degree-free.md)
  §6), now for **two** targets.
- The non-degenerate band-5 probe (`d>=4`) — out of budget.
- Whether the tail forces at depth `<2k-1` uniformly; the `kmin=2` value beyond `d=3`.
- DC1, JC2; no Weyl pair, no counterexample.

**Exceptional loci (explicit):** the hatch and `lambda_3` exist at every `k>=3`
(`gcd(k-1,k)=1`). At `k=3` (W2) the single-factor `a_2(0)` factorization holds
(`lambda_3=ev_-1`, the two targets coincide via `R(1)+R(-1)=0`); for `k>=4` the two
targets are independent and the single factor is gone — W2 is the degenerate base of
the backbone, not the pattern.

## 5. Verification

```sh
uv run --with sympy python research/dc1-program/verify_uniform_hatch_df.py
# heavy: det-saturated covector cert + branch A + mod 32003 + QQ char-0 + band-5 d=2
HEAVY=1 uv run --with sympy python research/dc1-program/verify_uniform_hatch_df.py
```

`S0` engine; `S1` `lambda_3` point functional + `Im Phi ⊂ ker lambda_3` + Phi-filler
independence; `S2` BACKBONE (both-ends Lemma-P, top-wall proportionality, shared
partial divisibility, single-factor refutation); `S3` depth-2 filler-linearity
(general `d`); `S4` band-4 `d=3` controls (cokernel dim `10`, free+distinct targets,
origin); `S5` THE COVECTOR KILL (`msolve` parser validation, both targets in
`sqrt(cascade+Q_-1+Q_-2)`, `kmin=2` witnesses); `S6` band-5 reduce-first degeneracy;
`S7` (`HEAVY`) the determinant-saturated covector certificate + branch A + second
prime + `QQ` + band-5 `d=2`. The ledger reports executed checks as PASS, unavailable
optional legs as SKIP, and ends `ALL UNIFORM HATCH DF CHECKS PASSED` only when the
band-4 `d=3` `msolve` covector-kill payload actually ran and returned UNIT; otherwise
`SKIP -- ... COVECTOR PAYLOAD NOT RUN; SUPPORTING CHECKS PASSED`.
