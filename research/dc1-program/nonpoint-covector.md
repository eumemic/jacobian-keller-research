# The band-k hatch covector is a POINT functional: correcting the "non-point" census claim

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER REVIEWED — BAND-SCOPED — MIXED PROVED/BOUNDED**

This memo is architecture step 2's "degree-free non-point covector `lambda_np`"
task for the band-`k` hatches, `k>=4`. It **corrects a committed structural
claim**: the census
([`hatch-census.md`](hatch-census.md) §2, §4) and its sibling
([`band45-lambda.md`](band45-lambda.md) §5,
[`lambda-general-k.md`](lambda-general-k.md) §6) assert that each `k>=4` hatch
carries **one NON-POINT (infinite-support) annihilator** of `Im Phi`. **It does
not.** The third annihilator is an explicit **POINT functional** `lambda_3`
(finite support), and this is provable degree-free for every `k`. The census's
"non-point" reading is a **computational artifact** of a too-narrow
point-annihilator search window.

Conventions frozen from the band-`k` corpus (`b7e85e8..e223d21`): `F` char zero,
`A_1[x^{-1}]=(+)_k x^k C[E]`, `(x^a f)(x^b g)=x^{a+b} f(E+b)g(E)`,
`f^{[n]}(E)=f(E+n)`, `T f(E)=f(E+1)`, `S_n=1+T+...+T^{n-1}`,
`Q_m=sum_{i+l=m}(b_l^{[i]}a_i - a_i^{[l]}b_l)`, gauge `b_k=0`, `u:=b_{k-1}`,
membership `(E)_j=E(E-1)...(E-j+1) | a_{-j},b_{-j}`, `Q_0=(T-1)G`, `G(0)=0`.

The hatch (unique common root translated onto the anchor `0`):

```text
band 4:  a_4 = E(E+3)(E+6)(E+9),          u = b_3 = E(E+4)(E+8),          b_4=0.
band 5:  a_5 = E(E+4)(E+8)(E+12)(E+16),   u = b_4 = E(E+5)(E+10)(E+15),   b_5=0.
general: a_k = prod_{i<k}(E+i(k-1)),      u = prod_{j<k-1}(E+jk).
```

Exact certificate: [`verify_nonpoint_covector.py`](verify_nonpoint_covector.py)
(ends `ALL NONPOINT COVECTOR CHECKS PASSED`; heavy `d=3` legs behind `HEAVY=1`).
Base commit `e223d21`.

## 0. Headline

> **The band-`k` hatch obstruction to `Q_0=1` is a POINT covector `lambda_3`, not
> a non-point one. `Ann(Im Phi) = span{ev_0, ev_1, lambda_3}` — three point
> functionals — at every band.**
>
> 1. **`lambda_3` is a point functional (PROVED, degree-free, arbitrary `k`).**
>    The moving-sum/adjoint criterion writes any annihilator's symbol
>    `L(t)=sum_x c_x t^x` with `S_k(t)L(t)=P_K(t)` (finite, supported on
>    `roots(q_K)`) and `S_{k-1}(t)L(t)=P_H(t)` (finite, supported on
>    `roots(q_H)`). Since `gcd(S_k,S_{k-1})=1` in `F[t,t^{-1}]`, `S_k | P_K`, so
>    `L=P_K/S_k` is a **Laurent polynomial** — finite support. Every annihilator
>    is therefore a point functional; there is no room in `coker Phi` (codim 3)
>    for an infinite-support one.
> 2. **The census artifact (identified).** `hatch-census.md`'s `point_annih` scans
>    support `range(-3, 3k+5)`, which starts at `-3`. The genuine `lambda_3` has
>    support at large-negative positions (down to `-(k^2-3k+1)`), **outside** that
>    window, so the scan finds only `{ev_0, ev_1}` and — since `codim=3` — wrongly
>    infers a hidden non-point annihilator. Widening the window recovers the third
>    **point** annihilator (verifier §1).
> 3. **Explicit closed form** (shortest-support representative), with
>    `lambda_3(E) = -(k-1)(k-2)/2` and `lambda_3(1)=1`:
>
>    | `k` | `lambda_3` (support : weights) | `lambda_3(E)` |
>    |---|---|---:|
>    | 3 (W2) | `{-1}` | `-1` |
>    | 4 | `ev_{-2}-ev_{-4}+ev_{-5}` | `-3` |
>    | 5 | `ev_{-3}-ev_{-5}+ev_{-7}-ev_{-10}+ev_{-11}` | `-6` |
>    | 6 | `ev_{-4}-ev_{-6}+ev_{-9}-ev_{-12}+ev_{-14}-ev_{-18}+ev_{-19}` | `-10` |
>    | 7 | 9-term, support `{-29..-5}` | `-15` |
>
> 4. **`Im Phi subset ker(lambda_3)` (PROVED, degree-free, `k=3..7`)** by the two
>    adjoint support identities and symbolic block annihilation in generic `C,V`.
> 5. **Residual (BOUNDED).** On the positive cascade, `lambda_3(R)=0` is forced at
>    filler degree `d<=2` (bands 4,5, both tail branches), so
>    `lambda_3(E-R)=lambda_3(E)=-(k-1)(k-2)/2 != 0` — a clean **single-covector**
>    kill for `d<=2`. **But the identity `lambda_3(R)=0` FAILS at `d=3`**: there
>    `lambda_3(R)` becomes a free modulus (the same freed cascade data that unpins
>    the slope), so the single-covector obstruction is **bounded to `d<=2`** and
>    the `d=3` kill is the **joint** `{R(1)=1, lambda_3(R)=lambda_3(E)}`
>    infeasibility (`Q_0=1` UNIT, reproducing the census verdict).
> 6. **No counterexample.** The hatch still dies (bounded); nothing here builds a
>    Weyl pair or settles DC1/JC2. The moment-unit kill's arbitrary-degree status
>    remains open, now with a **sharp** obstruction: a joint two-covector
>    condition, not one non-point covector.

## 1. Task 1 — `lambda_3` explicit, and it is a point functional

### 1.1 The moving-sum construction (the engine)

By the adjoint criterion ([`lambda-general-k.md`](lambda-general-k.md) Thm A′), a
point functional `lambda=sum_x c_x ev_x` annihilates the two filler blocks iff its
forward moving-sums land on the shifted root necklaces:

```text
B_k c    supported on  roots(q_K),   q_K = a_k(E-k)(E)_k,
B_{k-1}c supported on  roots(q_H),   q_H = u(E-(k-1))(E)_{k-1},
```

`B_n c(y)=sum_{i=0}^{n-1}c_{y-i}`. In generating-function form with
`L(t)=sum_x c_x t^x`, `S_n(t)=1+t+...+t^{n-1}`:

```text
S_k(t) L(t)     = P_K(t)   supported on roots(q_K),
S_{k-1}(t) L(t) = P_H(t)   supported on roots(q_H).
```

For the band-4 hatch, `roots(q_K)={-5,-2,0,1,2,3,4}` (double at `1`),
`roots(q_H)={-5,-1,0,1,2,3}`. Solving `P_K S_{k-1}=P_H S_k` gives a
**3-dimensional** solution space (`= codim Im Phi`), and each basis direction's
`L=P_K/S_k` is a **Laurent polynomial**. The shortest-support representative not
equal to `ev_0` (`L=t^0`) or `ev_1` (`L=t`) is

```text
band 4:  L = (t^3 - t + 1) / t^5   =>   lambda_3 = ev_{-2} - ev_{-4} + ev_{-5},
band 5:  L = (t^8 - t^6 + t^4 - t + 1) / t^{11}
             =>   lambda_3 = ev_{-3} - ev_{-5} + ev_{-7} - ev_{-10} + ev_{-11}.
```

### 1.2 Theorem (all annihilators are point functionals; degree-free, arbitrary `k`)

> **Theorem P.** For every band-`k` hatch, every functional annihilating
> `Im Phi` is a **point** functional (finite support). Equivalently, `coker Phi`
> carries **no** non-point annihilator.

**Proof.** An annihilator has `S_k L=P_K`, `S_{k-1}L=P_H` with `P_K,P_H` finite
Laurent polynomials (supported on the finite sets `roots(q_K)`, `roots(q_H)`).
Then `P_K S_{k-1} = P_H S_k`. In the Laurent ring `F[t,t^{-1}]` (a PID),
`gcd(S_k,S_{k-1})=1` (their non-`1` roots are `k`-th vs `(k-1)`-th roots of unity,
disjoint since `gcd(k,k-1)=1`). Hence `S_k | P_K`, so `L=P_K/S_k in F[t,t^{-1}]`
is a **Laurent polynomial** — finite support. ∎

The verifier checks `gcd(S_k,S_{k-1})=1` and the Laurent-polynomial property of
every basis symbol for `k=3..7`, and the closed forms
`lambda_3(E)=-(k-1)(k-2)/2`, `lambda_3(1)=1`. (Machine scope: `k<=7` for the
symbol solve, `k<=8` for the gcd; the arbitrary-`k` step is this paper argument
together with the *cited* adjoint criterion — see the §6 tier note.)

### 1.3 The census correction, pinned exactly

`hatch-census.md` computes the honest `codim Im Phi = 3` (left-nullspace, stable),
which is **correct**. Its error is the *inference* of a non-point annihilator from

```text
point_annih(a,u,k, list(range(-3, 3k+5)), ...)   ->   dim 2   (only ev_0, ev_1).
```

The window starts at `-3`; `lambda_3` needs support down to `-(k^2-3k+1)`
(`-5` at `k=4`, `-11` at `k=5`). Widening to `range(-(k^2-3k+1)-2, 3k+5)` yields
**dim 3**, all point. The verifier runs both windows side by side (§1): narrow → 2,
wide → 3. Since `{ev_0,ev_1,lambda_3}` are independent on the cokernel basis
`{1,E,E^2}` (determinant `16` at `k=4`, `60` at `k=5`) and `codim=3`, they
**exhaust** `Ann(Im Phi)`. No functional is left over to be non-point.

(For band 5 a naive finite-window scan also throws up two **spurious** large-
coefficient "annihilators" that are truncation artifacts of the finite image set;
the moving-sum symbol solve — which is exact — gives the clean 3-dimensional
answer and is what the verifier uses.)

## 2. Task 2 — `Im Phi subset ker(lambda_3)`, degree-free

Because `lambda_3` is a finite point functional, degree-freeness is the direct
block computation (as in W1 / `band45-lambda.md`), and its structural "why" is the
two support identities. With `q_K=a_k(E-k)(E)_k` one has the shift identity

```text
K_k[(E)_k C](E) = sum_{j=0}^{k-1} (q_K C)(E+j),
```

so `lambda_3(K_k[(E)_k C]) = ((q_K C)(θ) P_K)|_{t=1} = sum_{y in supp(P_K)} P_K[y]
q_K(y) C(y) = 0` because `supp(P_K) subset roots(q_K)`. Identically for the
`H`-block with `q_H, P_H`. The verifier checks, for `k=3..7`:

- `supp(S_k L) subset roots(q_K)` and `supp(S_{k-1}L) subset roots(q_H)` (the two
  support identities), the shift identity above, and
- `lambda_3(K_k[(E)_k C])=0`, `lambda_3(H_{k-1}[(E)_{k-1}V])=0` as **symbolic
  identities** in the free coefficients of generic degree-`k` `C,V` (arbitrary
  degree).

Hence `Im Phi subset ker(lambda_3)` degree-free, and `lambda_3(E)=-(k-1)(k-2)/2 !=
0` gives `E not in Im Phi` at arbitrary filler degree.

## 3. Task 3 — the residual `lambda_3(R)` and the exact boundary

### 3.1 The reduction (structural, exact)

`Q_0=(T-1)G` and `(T-1)E=1`, with `G(0)=0` from membership, give `Q_0=1 <=> G=E`.
Writing `G=R+K_k[b_{-k}]-H_{k-1}[a_{-(k-1)}]` with the fillers admissible, the two
filler blocks lie in `Im Phi subset ker(lambda_3)`, so `lambda_3(G)=lambda_3(R)`
and

```text
Q_0=1  <=>  E-R in Im Phi  <=>  ev_0(E-R)=ev_1(E-R)=lambda_3(E-R)=0
       <=>  R(0)=0 (auto),  R(1)=1,  lambda_3(R)=lambda_3(E).
```

The equivalence is exact because `{ev_0,ev_1,lambda_3}` is a **basis** of the
3-dimensional `coker Phi` (§1.3). The verifier reproduces the census moment-unit
verdict `cascade+Q_0=1 = UNIT` at `d=2` (both bands) and, in the heavy leg, at
`d=3` band 4 **both** via `coeffs(Q_0-1)` and via the 3-covector form — the two
agree, validating the reduction.

### 3.2 Bounded kill (`d<=2`) — the single-covector regime

On the forward-solved positive cascade, at filler degree `d<=2` (bands 4,5, tail
branches A and B), fixed-pivot linear elimination drives

```text
R(0) = 0    and    lambda_3(R) = lambda_3(G) = 0     (mod the positive cascade),
```

so `lambda_3(E-R)=lambda_3(E)=-(k-1)(k-2)/2 != 0`. This is a clean
**single-covector** obstruction: `Q_0=1` is impossible at `d<=2` regardless of the
slope, and it needs no Gröbner engine. (At `d<=2` the slope `R(1)` is itself pinned
to `0`, so `R(1)=1` already fails; the two obstructions are simultaneous.)

### 3.3 The `d=3` boundary — the identity FAILS, the kill goes joint

At `d=3` (band 4) the cascade frees moduli. The verifier's heavy leg shows:

```text
cascade + Q_0=1                                  = UNIT   (hatch dies; = census)
lambda_3(R) mod cascade                          = a FREE, non-constant modulus
```

So `lambda_3(R)=0` is **false** at `d=3`: `lambda_3(R)` acquires the very cascade
moduli that also unpin the slope `R(1)`, and can take the value `lambda_3(E)`. The
single covector no longer obstructs by itself. The kill survives only as the
**joint** infeasibility `{R(1)=1  AND  lambda_3(R)=lambda_3(E)}` on the cascade —
`Q_0=1` is `UNIT` because those two scalar conditions cannot hold together with the
cascade, even though each is individually reachable.

**This is the precise statement of "what blocks the degree-free residual
identity."** At W1 (`e4e704f`) `lambda_0(R)=0` telescoped degree-free by
cascade-block elimination with fixed rational pivots. Here the analogous telescope
closes `lambda_3(R)=0` only for `d<=2`; at `d=3` the same elimination leaves a
non-zero residual that is a genuine function of the freed positive-cascade moduli.
The obstruction is thus **not** a single degree-free covector identity for `k>=4`;
it is a joint slope + `lambda_3` condition — structurally the moment-unit analogue
of W2's joint slope + tail obstruction ([`w2-joint-theorem.md`](w2-joint-theorem.md)),
but with the covector `lambda_3` in place of the negative tail.

## 4. Task 4 — general `k`

**Closed form (proved `k=3..7`, uniform construction all `k`).** `lambda_3` is the
point functional with symbol `L=P_K/S_k` from §1.1; its shortest-support
representative has `2k-5` alternating `+1/-1` weights on support
`subset [-(k^2-3k+1), -(k-2)]`, with

```text
lambda_3(E) = -(k-1)(k-2)/2,     lambda_3(1) = 1.
```

**Degree-free annihilation** `Im Phi subset ker(lambda_3)` is verified for `k=6`
and `k=7` (§2 loop) — the support identities and symbolic block annihilation go
through identically, confirming Theorem P and the construction are not special to
bands 4,5.

**Residual conjecture (general `k`).** On the positive cascade, `lambda_3(R)=0` is
forced for filler degree up to some threshold `d_0(k)` (verified `d_0>=2` at bands
4,5), giving `lambda_3(E-R)=lambda_3(E)=-(k-1)(k-2)/2 != 0` there; beyond `d_0`,
`lambda_3(R)` acquires the freed cascade moduli and the arbitrary-degree kill is the
**joint** `{R(1)=1, lambda_3(R)=lambda_3(E)}` infeasibility. The `k=4` evidence
(`d<=2` single-covector, `d=3` joint-UNIT) is attached in §3.

## 5. Sanity control

`lambda_3(E) = -(k-1)(k-2)/2 != 0` for every `k>=3` (verifier §5). So the covector
is a live obstruction: even where the pure slope `R(1)=1` is achievable, `Q_0=1`
additionally requires `lambda_3(R)=lambda_3(E) != 0`, which is not free of charge.

## 6. Claim disposition

**Proved (exact algebra, arbitrary `k` / arbitrary degree as stated):**
- **Theorem P**: every annihilator of `Im Phi` is a POINT functional. **Tier
  (audit-corrected):** machine-certified at `k=3..7` (moving-sum dimension;
  `gcd(S_k,S_{k-1})=1` spot-checked `k=3..8`); arbitrary `k` rests on the
  (correct, elementary) paper PID/gcd argument of §1.2 **plus** the adjoint
  criterion Thm A′, which is *cited* from
  [`lambda-general-k.md`](lambda-general-k.md), not re-derived symbolically in
  `k` here. The census's "non-point covector" does not exist; its appearance is
  a search-window artifact.
- `lambda_3(E)=-(k-1)(k-2)/2`, `lambda_3(1)=1`, and the moving-sum symbol
  `L=P_K/S_k` closed form (`k=3..7`).
- `Im Phi subset ker(lambda_3)` degree-free via the two support identities and
  symbolic block annihilation in generic degree-`k` `C,V` (`k=3..7`).
- `Q_m=[D,X]_m`, `Q_0=(T-1)G`, `G(0)=0` at `k=4,5`.

**Proved (exact finite):**
- honest `codim Im Phi = 3` (stable left-nullspace) `= dim span{ev_0,ev_1,lambda_3}`;
  the three are independent, hence exhaust `Ann(Im Phi)`. **Scope
  (audit-corrected): the stable-left-nullspace + det-exhaustion route is machine
  run at bands 4,5**; band 3 is covered by the moving-sum dimension (`=3`,
  `k=3..7`) and the prior W2 point-completeness, not by this route.
- the reduction `Q_0=1 <=> {R(0)=0, R(1)=1, lambda_3(R)=lambda_3(E)}`, and its
  agreement with the census moment-unit `UNIT` verdict (`d=2` both bands, sympy;
  `d=3` band 4, `msolve ^`, heavy leg).

**Bounded / per-degree evidence:**
- `lambda_3(R)=0` on the positive cascade at `d<=2` (bands 4,5, both tail
  branches), giving the single-covector kill `lambda_3(E-R)=lambda_3(E) != 0`.
- `d=3` band 4: `lambda_3(R)` is a free modulus (residual identity FAILS) and
  `cascade+Q_0=1` is `UNIT` (joint kill) — heavy leg, `msolve ^` + mod-`p`
  Gröbner.

**Refuted (corrections to committed memos):**
- `hatch-census.md` §0.1/§2/§4 "ONE NON-POINT annihilator `lambda_np`", and the
  matching escape-hatch framing of `band45-lambda.md` §5 and `lambda-general-k.md`
  §6, insofar as they posit a *non-point* obstruction at the hatch: the obstruction
  is the point functional `lambda_3`.
- the census's implied "single-functional degree-free kill is the only missing
  step": for `k>=4` the arbitrary-degree kill is a **joint** two-covector
  condition, since `lambda_3(R)=0` provably fails at `d=3`.

**Open / not claimed:**
- the arbitrary-degree threshold `d_0(k)` and the joint-obstruction
  infeasibility at arbitrary degree (the moment-unit kill remains bounded, `d<=3`).
- general-`k` `codim Im Phi = 3` beyond `k=3..7` (verified windowed).
- DC1, JC2; no Weyl pair, no counterexample.

**Exceptional loci (explicit):** the hatch and `lambda_3` exist at every `k>=3`
(`gcd(k-1,k)=1`); at `k=3` (W2) `lambda_3=ev_{-1}` is the third *point* member of
the already point-complete triple `{ev_{-1},ev_0,ev_1}`, so W2 is the base case of
Theorem P, not an exception.

## 7. Verification

```sh
uv run --with sympy python research/dc1-program/verify_nonpoint_covector.py
# heavy d=3 boundary (msolve + mod-p Groebner):
HEAVY=1 uv run --with sympy python research/dc1-program/verify_nonpoint_covector.py
```

§0 engine (`Q_m=[D,X]`, `Q_0=(T-1)G`, `G(0)=0`); §1 `lambda_3` point functional,
Theorem P, closed forms, honest codim, census-window artifact; §2 degree-free
`Im Phi subset ker(lambda_3)`; §3 the reduction, `d<=2` single-covector kill,
`d=3` joint boundary; §4 general `k` (`k=6,7`); §5 sanity. Ends
`ALL NONPOINT COVECTOR CHECKS PASSED`.
