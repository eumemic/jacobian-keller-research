# The unification probe: the k>=4 hatch tail is slope-forcing in EVERY cokernel direction (band 4, d=3)

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED — BAND-SCOPED — MIXED PROVED/BOUNDED**

This memo runs architecture's **unification probe**: are the `k>=4` hatch tails
*also* slope-forcing? At W2 (band 3) the negative tail `Q_-1..Q_-5` (with **no**
`Q_0` condition) forces the moment slope `R(1)=0` on the `cascade+tail` variety
([`slope-forcing-verdict.md`](slope-forcing-verdict.md)). For the `k>=4` hatches
the cokernel `Ann(Im Phi) = span{ev_0, ev_1, lambda_3}` carries a **second**
nontrivial covector — the point functional `lambda_3`
([`nonpoint-covector.md`](nonpoint-covector.md)) — beyond the slope `ev_1`. The
question is whether the tail forces **both** cokernel pairings of the residual `R`
to vanish:

```text
   on {cascade, tail (no Q_0), membership}:     R(1) = 0   AND   lambda_3(R) = 0 ?
```

If yes, architecture steps **1 (slope-forcing)** and **2 (covector)** collapse into
**one** mechanism — the tail annihilates `E-R` against *every* cokernel direction
beyond `ev_0`, so `Q_0=1` is contradicted in all directions at once.

Conventions frozen from the band-`k` corpus (gauge `b_k=0`): `(x^a f)(x^b g)=
x^{a+b} f(E+b)g(E)`, `f^[n](E)=f(E+n)`, `Q_m=sum_{i+l=m}(b_l^[i]a_i-a_i^[l]b_l)`,
`u:=b_{k-1}`, membership `(E)_j | a_-j,b_-j`, `Q_0=(T-1)G`, `G(0)=0`,
`R(1)=Q_0(0)=G(1)`, `lambda_3(R)=lambda_3(G)`. Hatch:
`a_k=prod_{i<k}(E+i(k-1))`, `u=prod_{j<k-1}(E+jk)`. Band 4:
`a_4=E(E+3)(E+6)(E+9)`, `u=E(E+4)(E+8)`, `lambda_3=ev_-2-ev_-4+ev_-5`,
`lambda_3(E)=-3`. Tail `= Q_-1..Q_-(2k-1)`; `Q_-2k` splits branch A
(`a_-k=(E)_k am_k`, `b_-k=mu_k a_-k`) from branch B (`a_-k=0`, `b_-k=(E)_k C`).

Verifier: [`verify_hatch_slope_forcing.py`](verify_hatch_slope_forcing.py). It
prints a PASS/SKIP/FAIL evidence ledger; without `msolve`, the load-bearing `d=3`
probe is explicitly `SKIP` and the final status does not claim the headline passed.
Heavy corroboration remains behind `HEAVY=1`. Every
load-bearing upstream fact — the crossed-product engine, `Q_0=(T-1)G`, `lambda_3`
and its point-functional closed form, `Im Phi subset ker lambda_3`, the
`Phi`-filler independence of the two targets, the positive cascade — is
**re-derived in file**.

## 0. Headline

> **YES at the band-4 hatch, `d=3`, on BOTH branches (bounded; `msolve '^'` mod
> `65003` + mod `32003` + `QQ` char-0, parser-validated):**
>
> ```text
>    R(1) IN sqrt(cascade+tail)      AND      lambda_3(R) IN sqrt(cascade+tail).
> ```
>
> The negative tail forces **both** cokernel pairings of the residual to vanish.
> Since `lambda_3(E)=-3 != 0` and `ev_1(E)=1 != 0`, on `cascade+tail` one has
> `ev_1(E-R)=1-R(1)=1 != 0` and `lambda_3(E-R)=lambda_3(E)-lambda_3(R)=-3 != 0`, so
> `E-R notin Im Phi` in **every** nontrivial cokernel direction at once, i.e.
> `Q_0=1` is contradicted simultaneously in all of `coker Phi` beyond `ev_0`.
> **Architecture steps 1 and 2 collapse into ONE mechanism** — the tail — at the
> band-4 hatch (`d=3`).
>
> - **Non-trivial:** on the cascade **alone** at `d=3`, `R(1)` and `lambda_3(R)`
>   are **free, non-constant, and independent** moduli (GB normal forms are
>   distinct non-constants). The tail forces both to `0`.
> - **Graded depth:** both targets are forced at the **same** minimal depth
>   `kmin=2` (`cascade+Q_-1..Q_-2` is already unit), and `Q_-1` alone forces
>   **neither** (explicit feasible witnesses with `R(1)=1`, `lambda_3(R)=1`). W2
>   needed depth `3` for its single condition; band 4 forces **two** conditions at
>   depth `2`.
> - **W2 is the base case:** at `k=3`, `lambda_3=ev_-1`, so `lambda_3(R)=R(-1)`,
>   and the cascade relation `R(1)+R(-1)=0` makes the two conditions **coincide** —
>   the W2 slope-forcing result IS the `k=3` instance of the unified principle. For
>   `k>=4` the two conditions are genuinely independent, and both still vanish.
> - **Degenerate low degree:** at `d<=2` (bands 4 and 5) *both* targets are already
>   `0` on the cascade **alone**, so the tail is not tested there; `d=3` (band 4) is
>   the first non-degenerate probe.
> - **No collapse to a monomial.** Unlike W2 (`R(1)=-108 a2_0 am1_3`, a single
>   monomial on the rationally-parametrized cascade), the band-4 `d=3` cascade is
>   **not** linearly parametrizable (4 residual degree-`4/5` conditions on 7 gens),
>   and the forced values `R(1)`, `lambda_3(R)` reduce to genuine **multinomials**
>   mod the cascade (20 terms each). Only their *vanishing on the tail*
>   mirrors W2 — not their cascade-collapsed form.

This is a **bounded** result (`d=3`, band 4). It builds no Weyl pair and settles
neither DC1 nor JC2. Its value is the **mechanism**: the tail, not `Q_0`, does the
killing, and it kills in every cokernel direction.

## 1. What is proved vs. what is bounded

### 1.1 Proved (arbitrary degree, machine identities; `k`-scope per the §5 audit tier notes)

- **Engine** (`k=3,4,5`): `Q_m=[D,X]_m` for `m in [-2k,2k]`, `Q_0=(T-1)G`,
  `G(0)=0` under membership (so `R(1)=Q_0(0)=G(1)`).
- **`lambda_3` is a point functional** with closed form re-derived by the
  moving-sum symbol solve (`k=3,4,5`): the annihilator space is 3-dimensional, every
  symbol `L=P_K/S_k` is a Laurent polynomial (finite support), `gcd(S_k,S_{k-1})=1`
  (`k=3..8`); `lambda_3(E)=-(k-1)(k-2)/2`, `lambda_3(1)=1`.
- **`Im Phi subset ker lambda_3`** degree-free: `lambda_3(K_k[(E)_k C])=0` and
  `lambda_3(H_{k-1}[(E)_{k-1}V])=0` as symbolic identities in generic degree-`k`
  `C,V` (`k=3,4,5`).
- **Both targets are `Phi`-filler-independent** (`k=4,5`): `R(1)=Q_0(0)=G(1)` and
  `lambda_3(R)=lambda_3(G)` do not depend on the level-`-(k-1)` and level-`-k`
  fillers `a_-(k-1), b_-k` — the content of Lemma P together with
  `Im Phi subset ker lambda_3`. (For `k>=4` both targets DO depend on the *middle*
  filler `a_-2`, which is part of `R`; this is a genuine richness beyond W2.)

### 1.2 Bounded-finite evidence (exact scope stated)

- **THE PROBE (band 4, `d=3`, both branches A and B):**
  `R(1) IN sqrt(cascade+tail)` **and** `lambda_3(R) IN sqrt(cascade+tail)`.
  Certified by `msolve 0.10.1 '^'` Rabinowitsch (`cascade + tail + {1 - t*target}`
  is the unit ideal) on the full un-reduced `d=3` system: mod `65003` (default), and
  mod `32003` **and over `QQ`** (char-0, obtained in ~9s per call) under `HEAVY` —
  all four (branch x target) combinations UNIT under every characteristic. The
  `msolve` UNIT parsing is **validated in file** against a known unit ideal `{1}`
  (-> UNIT) and a known feasible ideal `{x(x-1)}` (-> not unit).
- **Controls (band 4, `d=3`):** on the cascade **alone**, `R(1)` and `lambda_3(R)`
  have **non-constant** normal forms modulo the sympy `GF(p)` Gröbner basis of the
  cascade, and are **distinct** (`R(1)-lambda_3(R) != 0` mod cascade). So both are
  genuine free moduli and the tail-forcing is non-vacuous. Feasibility: the
  **origin** lies on `cascade+tail` (both branches) with `R(1)=lambda_3(R)=0`.
- **Graded depth (band 4, `d=3`, branch B):** `kmin(R(1))=kmin(lambda_3(R))=2` —
  `cascade+Q_-1..Q_-2` forces each target (msolve UNIT), while `Q_-1` alone forces
  neither (explicit feasible witness with target `=1`, msolve GB non-unit on a
  0-dim slice of the reduced cascade coordinates).
- **Forced-value structure (band 4, `d=3`):** `R(1)` mod cascade is a 20-term
  multinomial; `lambda_3(R)` mod cascade is a multinomial — **no** collapse to a
  single monomial (contrast W2).
- **Degenerate low degree:** at `d<=2` (bands 4, 5, both branches) both targets are
  already `0` on the cascade **alone** (sympy `GF(p)` linear elimination), and under
  `HEAVY` `msolve '^'` confirms UNIT on `cascade+tail` at band 5 `d=2`.

### 1.3 The second engine, honestly

The load-bearing `d=3` UNIT verdicts rest on `msolve '^'` (parsing validated; two
primes + `QQ` char-0). The **independent sympy Gröbner engine** certifies, over
`GF(p)` exactly and without `msolve`: the whole engine and `lambda_3` layer, the
`d=2` degenerate probe (targets reduce to `0`), the `d=3` free-modulus and
distinctness controls (GB of the cascade), and the feasibility origin. sympy's
Buchberger on the **full** `d=3` band-4 system (15 variables after linear
reduction) does **not** terminate in a sane budget and is reported as such — the
`d=3` radical-membership verdict itself is `msolve`-only, mitigated by the parser
validation, the two primes, and the `QQ` char-0 certificate.

## 2. The unified statement (conjecture) — stated exactly

The `d=3` band-4 evidence, together with the W2 (`k=3`) instance, supports:

> **Uniform Hatch Slope-Forcing Conjecture.** At every band-`k` hatch (`k>=3`), on
> the variety `{cascade (Q_{2k-2}..Q_1=0), tail (Q_-1..Q_-(2k-1)=0, NO Q_0),
> membership}`, both branches, at arbitrary positive-data degree,
>
> ```text
>    R(1) = 0        AND        lambda_3(R) = 0,
> ```
>
> i.e. `lambda(E-R) = lambda(E)` for **every** `lambda in Ann(Im Phi)` beyond
> `ev_0` (the tail annihilates `E-R` against the entire cokernel beyond the point
> `ev_0`). Equivalently, the negative tail is **slope-forcing in every cokernel
> direction**, and `Q_0=1` is contradicted on `cascade+tail` in all directions
> simultaneously.

Status of the conjecture:

- `k=3` (W2): the `R(1)=0` half is the bounded `d<=3` result of
  [`slope-forcing-verdict.md`](slope-forcing-verdict.md); the `lambda_3(R)=0` half
  is **equivalent** to it (`lambda_3=ev_-1`, `lambda_3(R)=R(-1)=-R(1)` on the
  cascade), so the two coincide. **Base case, bounded `d<=3`.**
- `k=4`: **both** halves verified here at `d=3`, both branches (bounded) — `msolve
  '^'` UNIT mod `65003` + mod `32003` + `QQ` char-0, parser-validated, with the sympy
  engine certifying the surrounding controls. At `d<=2` both are degenerate (forced
  on the cascade alone).
- `k=5`: `d<=2` degenerate (both forced on the cascade alone); the non-degenerate
  `d>=3` band-5 tail probe is **open / not run** (budget).
- arbitrary degree, arbitrary `k`: **open**.

Because the conjecture asserts *both* independent cokernel pairings vanish on the
tail, and because at `k=4` they are provably independent free moduli on the cascade,
this is a strictly stronger statement than the W2 slope-forcing lemma — and its
`d=3` band-4 confirmation is the concrete evidence that steps 1 and 2 are one
mechanism.

## 3. Why the low-degree cases are degenerate (and why `d=3` is the real test)

On the positive cascade, the number of free moduli grows with the raw coefficient
cap `d`. At `d<=2` (band 4) and `d<=2` (band 5) the cascade **pins both targets to
`0`** by linear elimination alone — so on `cascade+tail` they are trivially `0`,
and the tail plays no role. The first degree at which `R(1)` and `lambda_3(R)`
become genuine free moduli on the band-4 cascade is `d=3` (matching the census's
"pure slope frees up at `d=3`" for band 4). That is exactly where the probe is
informative, and exactly where both are forced back to `0` by the tail. For band 5
the slope frees up only at `d>3`, so the non-degenerate band-5 probe is out of the
`msolve` budget here and is left open.

## 4. Relation to the moment-unit kill (they are different systems)

The census ([`hatch-census.md`](hatch-census.md)) kills the `k>=4` hatch at the
**moment unit** `cascade+Q_0=1` (no tail needed). That is a different system from
this one (`cascade+tail`, **no** `Q_0`). The two are complementary faces of the
same cokernel geometry:

- **Moment unit** (`cascade + Q_0=1`): demands `E-R in Im Phi`, i.e.
  `R(1)=1 AND lambda_3(R)=lambda_3(E)`; this is infeasible on the cascade (UNIT),
  killing the hatch.
- **Slope-forcing tail** (`cascade + tail`, no `Q_0`): forces `R(1)=0 AND
  lambda_3(R)=0`, i.e. `E-R` is *maximally outside* `Im Phi` (every cokernel
  pairing beyond `ev_0` is non-zero). This is the W2-style mechanism, now shown to
  operate in **both** cokernel directions for `k>=4`.

The unification is that the **tail** — the same object that forces the slope at W2
— is what forces the entire cokernel data for `k>=4`; the extra covector `lambda_3`
that the census needed at the moment unit is *also* slope-forced by the tail.

## 5. Evidence ledger

**Proved (arbitrary degree, exact machine identities; `k`-scope audit-corrected):**
- Engine `Q_m=[D,X]_m`, `Q_0=(T-1)G`, `G(0)=0` (`k=3,4,5`).
- `lambda_3` point functional; `Im Phi subset ker lambda_3` degree-free; both
  targets `Phi`-filler-independent. **Tier note (audit):** the arbitrary-DEGREE
  property of each is genuine (independently reconfirmed in audit at cofactor
  degree 12 / raw `d=3,4`), but the in-file witnesses are fixed low degree
  (engine at generic deg 2; blocks at cofactor deg `k`; filler-independence at
  `d=2`) — the degree-generality rests on the structural arguments (bilinearity;
  `a_k(0)=u(0)=0`), stated in prose.
- `lambda_3` closed form `lambda_3(E)=-(k-1)(k-2)/2`, `lambda_3(1)=1`, every
  annihilator symbol Laurent. **`k`-scope (audit-corrected): machine-exhibited
  only `k in {3,4,5}`** (gcd `k=3..8`); the arbitrary-`k` extension is CITED
  from [`lambda-general-k.md`](lambda-general-k.md) /
  [`nonpoint-covector.md`](nonpoint-covector.md), not re-derived symbolically in
  `k` here. Not load-bearing for the band-4 headline (which uses only the fully
  re-derived `k=4` `lambda_3`).
- The `k=3` reduction *"the two conditions coincide at W2
  (`lambda_3=ev_-1`, `lambda_3(R)=-R(1)` via `R(1)+R(-1)=0` on the cascade)"*:
  the identification `lambda_3=ev_-1` is machine-checked; the coincidence
  relation `R(1)+R(-1)=0` on the `k=3` cascade is **prose/cited**
  (`slope-forcing-verdict.md` / `w2-joint-theorem.md`), not re-derived in-file
  (audit note; not load-bearing for band 4).

**Bounded-finite evidence (exact scope stated):**
- **Band 4, `d=3`, both branches:** `R(1) IN sqrt(cascade+tail)` AND
  `lambda_3(R) IN sqrt(cascade+tail)` — `msolve '^'` UNIT mod `65003` (default) and
  mod `32003` + over `QQ` (char-0) (`HEAVY`), all four branch x target combinations;
  parser validated in file.
- Band 4, `d=3`: free-modulus + distinctness controls (sympy GB of the cascade);
  feasibility origin; graded `kmin(R(1))=kmin(lambda_3(R))=2` with explicit `Q_-1`
  witnesses; forced values are non-monomial (20 terms each mod cascade).
- Band 4/5, `d<=2`: both targets degenerate (`0` on the cascade alone); band 5
  `d=2` `cascade+tail` UNIT (`HEAVY`, `msolve '^'`).

**Refuted / corrected:** nothing prior is refuted. This memo *extends* the W2
slope-forcing verdict (`slope-forcing-verdict.md`) from the single slope covector to
the full cokernel at `k=4`, and gives the sharp reason the W2 result is the `k=3`
instance (`lambda_3=ev_-1`, `lambda_3(R)=-R(1)` on the cascade).

**Open / not claimed:**
- Arbitrary degree at any band; the arbitrary-degree Uniform Hatch Slope-Forcing
  Conjecture.
- The non-degenerate (`d>=3`) band-5 tail probe (budget).
- Independent (non-`msolve`) certification of the `d=3` radical membership
  (sympy Buchberger intractable at 15 variables).
- DC1, JC2; no Weyl pair, no counterexample.

## 6. Verification

```sh
uv run --with sympy python research/dc1-program/verify_hatch_slope_forcing.py
# heavy legs (second prime + QQ char-0 + band 5 d=2 + full graded sweep):
HEAVY=1 uv run --with sympy python research/dc1-program/verify_hatch_slope_forcing.py
```

`S0` engine; `S1` `lambda_3` point functional + `Im Phi subset ker lambda_3` +
`Phi`-filler independence + the W2 base case; `S2` controls (d=2 degenerate, d=3
free moduli, distinctness); `S3` feasibility origin; `S4` THE PROBE (msolve parser
validation, d=2 degenerate, d=3 both branches both targets); `S5` graded depth
(`kmin=2`, `Q_-1` witnesses); `S6` forced-value structure (no monomial collapse);
`S7` (`HEAVY`) band 5 `d=2`. The final evidence ledger reports the executed
supporting checks as PASS, unavailable optional legs as SKIP, and the headline as
`PASS -- ... PAYLOAD PASSED` only when all four default `msolve` radical-membership
probes (both targets, both branches, mod `65003`) actually ran and returned UNIT.
Without `msolve` it ends `SKIP -- ... PAYLOAD NOT RUN; SUPPORTING CHECKS PASSED`.
