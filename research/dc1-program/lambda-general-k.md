# The cokernel dual: why the annihilating functional exists, at every band

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER REVIEWED — BAND-SCOPED**

This note answers the structural question behind the W1 template
([`quantum-w1-arbitrary-degree.md`](../band3/quantum-w1-arbitrary-degree.md),
commit `e4e704f`): **why** does an annihilating functional `lambda` of the
two-filler image `Phi` exist, and does one exist at every band `k`? The answer is
an adjoint (cokernel-dual) computation that requires no coefficient-degree bound.
It reproduces the W1 certificate `{ev_0, lambda_0}`, `lambda_0(E-R)=4`, exactly.

What is **proved here degree-free and arbitrary-`k`** is the cokernel-dual
characterization (the adjoint criterion), the primal codimension formula
`codim Im L_K = deg a+k`, the always-present functional `ev_0` (hence
`coker Phi != 0` at every band), and the moving-sum construction of annihilators. What is
**proved arbitrary-degree for the whole band-3 arithmetic-progression (AP)
family** is the explicit second functional `lambda_r`, `Im Phi subset ker
lambda_r`, and `lambda_r(E)=r+4`. What is **only bounded evidence** (bands 3, 4,
5; integer-rooted tame tops) is `codim Im Phi >= 2` in general and the
`lambda(E) = common root` reading. What **remains open** is the residual identity
`lambda(R)=0` beyond W1 and, with it, the full band-`k` exotic closure; the
precise escape hatch is stated loudly in §6. Nothing here constructs a Weyl pair
or settles DC1/JC2.

Exact checks: [`verify_lambda_general_k.py`](verify_lambda_general_k.py) (ends
`ALL LAMBDA GENERAL K CHECKS PASSED`, 88 checks).

Conventions are frozen from the band-`k` corpus (`b7e85e8..e4e704f`): `F` has
characteristic zero, `T f(E)=f(E+1)`, `S_n=1+T+...+T^(n-1)`,
`(E)_j=E(E-1)...(E-j+1)`, membership `(E)_j | a_{-j},b_{-j}`, quantum ladder
`Q_m=sum_(i+l=m)[b_l^{[i]}a_i - a_i^{[l]}b_l]`, gauge `b_k=0`, `u:=b_{k-1}`.

---

## 1. The map, and the two filler images

After the top Wronskian gauge `b_k=0` (W1, [`band-k-weapons.md`](../band3/band-k-weapons.md)),
solving the positive cascade splits the central potential (Lemma P frame,
[`moment-unit-general-k.md`](moment-unit-general-k.md)) as

```text
G = R + K_k[b_-k] - H_{k-1}[a_-(k-1)],
K_k[c]     = sum_(r=1)^k     a(E-r) c(E+k-r),      a := a_k,
H_{k-1}[v] = sum_(r=1)^(k-1) u(E-r) v(E+k-1-r),    u := b_{k-1},
```

with `R` the canonical non-filler residual and, under membership
`c=(E)_k C`, `v=(E)_(k-1) V`,

```text
Phi(C,V) = K_k[(E)_k C] - H_{k-1}[(E)_(k-1) V].
```

Because `Q_0=(T-1)G` and `(T-1)E=1`, the unit equation `Q_0=1` forces `G=E`
(membership gives `G(0)=0`, killing the additive constant), i.e.

```text
Q_0 = 1   <=>   E - R  in  Im Phi.
```

The kill mechanism is a functional `lambda` with `Im Phi subset ker(lambda)` and
`lambda(E-R) != 0`.

**Shift forms (two-filler memo, `0889f8a`).** Re-indexing the sums,

```text
K_k[c]     = S_k((T^-k a) c),          H_{k-1}[v] = S_{k-1}((T^-(k-1) u) v).
```

With membership factors absorbed,

```text
Im L_K = S_k(q_K F[E]),      q_K = (T^-k a)(E)_k     = a(E-k)(E)_k,
Im L_H = S_{k-1}(q_H F[E]),  q_H = (T^-(k-1) u)(E)_(k-1) = u(E-k+1)(E)_(k-1).
```

Since `S_n` is a **linear automorphism** of `F[E]` (degree preserved, leading
coefficient `n`; two-filler memo), `Im L_K` and `Im L_H` are automorphic images
of the ideals `q_K F[E]`, `q_H F[E]`, of codimensions `deg a + k` and
`deg u + (k-1)`.

---

## 2. Theorem A — the cokernel dual (degree-free, arbitrary `k`)

Two claims must be kept apart: the **primal codimensions** (exact, degree-free)
and the **explicit point functionals** that realise annihilators. They are not the
same count — see the caveat below.

> **Theorem A (primal codimension).** For nonzero `a,u` and every `k>=2`,
> `S_k` and `S_{k-1}` are linear automorphisms of `F[E]`, so
>
> ```text
> codim Im L_K = deg q_K = deg a + k,      codim Im L_H = deg q_H = deg u + k-1.
> ```
>
> Hence `Im Phi = Im L_K + Im L_H` has finite codimension, bounded by
> `min(deg a+k, deg u+k-1)`. In particular `coker Phi` is finite-dimensional and
> `codim Im Phi` is a well-defined invariant, computed exactly by finite
> truncation.

**Proof.** `q_K F[E]` has codimension `deg q_K` in `F[E]`; the automorphism `S_k`
preserves codimension, and `Im L_K = S_k(q_K F[E])`. Identically for `H`. A sum of
subspaces has codimension at most the smaller codimension. ∎

Now the adjoint construction — **the "why" an annihilator exists**. Write a point
functional as `lambda = sum_x c_x ev_x`. The shift adjoint translates support,
`T^* ev_x = ev_{x+1}`, hence `S_n^* ev_x = ev_x + ev_{x+1} + ... + ev_{x+n-1}`.

> **Theorem A' (adjoint criterion, degree-free).** For every point functional
> `lambda` and every `C in F[E]` of any degree,
>
> ```text
> lambda(K_k[(E)_k C]) = (S_k^* lambda)(q_K C).
> ```
>
> Consequently `lambda in Ann(Im L_K)` iff `S_k^* lambda` is supported on
> `V(q_K)` (order below multiplicity at each root); likewise `lambda in
> Ann(Im L_H)` iff `S_{k-1}^* lambda` is supported on `V(q_H)`. Any `lambda`
> satisfying **both** support conditions annihilates `Im Phi`.

**Proof.** `lambda(K_k[(E)_k C]) = lambda(S_k(q_K C)) = (S_k^* lambda)(q_K C)` by
the definition of the adjoint; this is an identity in `C`, no degree bound. It
vanishes for all `C` iff `S_k^* lambda` kills the ideal `q_K F[E]`, i.e. is
supported on `V(q_K)`. Same for `H`; `Ann(Im Phi)=Ann(Im L_K) intersect
Ann(Im L_H)`. ∎

**This is the "why".** An annihilating functional is a combination of evaluations
whose forward moving-sum `S_k^*` lands exactly on the shifted root necklace
`V(q_K)` and whose `S_{k-1}^*` lands on `V(q_H)` — the dual-basis functionals at
the two root necklaces, pulled back through the two coprime shift windows. On
integer cosets this is the moving-sum criterion of §2.1. The verifier checks the
automorphism facts (`deg`, `lc=n`), the adjoint identity
`lambda(S_n g)=(S_n^* lambda)(g)`, and the shift forms at `k=2..5`.

**Caveat (do not conflate the two counts).** The point functionals satisfying the
support conditions form a *subspace* of the full annihilator. For a **single
block** this subspace can be strictly smaller than the primal codimension: at W1,
`codim Im L_K = 6` but only `4` annihilators of `Im L_K` are finitely-supported
point functionals (the moving-sum `S_k^*` is not surjective onto functions
supported on `V(q_K)` by finitely-supported coefficients; the missing two carry
infinite support). What is verified case by case (§5, W1 and the AP family) is
that for the **intersection** `Ann(Im Phi)` the point functionals *do* exhaust the
cokernel dual — their dimension equals the truncation-computed `codim Im Phi`. So
the explicit functionals below are honest annihilators (Theorem A'), and give a
lower bound `codim Im Phi >= #(independent point functionals)`; their
*completeness* is a checked fact, not a theorem of §2.

### 2.1 The moving-sum criterion (integer cosets)

When the roots lie on `Z` (the tame/AP case), a point functional `lambda=sum_x c_x
ev_x` has `S_n^* lambda = sum_y (B_n c)(y) ev_y` where `B_n` is convolution with
the box `{0,...,n-1}`, `(B_n c)(y)=sum_{i=0}^{n-1} c_{y-i}`. Theorem A' becomes:

```text
lambda in Ann(Im Phi)   <=>   B_k     c  supported on roots(q_K)
                         and   B_{k-1} c  supported on roots(q_H).
```

Writing the cumulative sum `C(y)=sum_{t<=y} c_t`, this reads `C(y)=C(y-k)` off
`roots(q_K)` and `C(y)=C(y-(k-1))` off `roots(q_H)`: the cumulative profile is
constant along each residue class mod `k` (resp. mod `k-1`) except at the allowed
root jumps. Because `gcd(k,k-1)=1` the two window conditions are rigidly coupled.
The verifier solves this finite linear system for finitely-supported `c`; its
solution space is the space of point annihilators, and in every computed case
(§5) its dimension equals the truncation-computed `codim Im Phi`, so the point
functionals are a complete set of obstructions there.

---

## 3. Theorem B — `ev_0`, and coker `Phi != 0` at every band (degree-free)

> **Theorem B.** For every `k>=2` and all admissible `c=(E)_k C`,
> `v=(E)_(k-1) V`, `K_k[c](0)=0` and `H_{k-1}[v](0)=0`. Hence `ev_0 in
> Ann(Im Phi)`, so `coker Phi != 0` at every band.

**Proof.** `K_k[c](0)=sum_{r=1}^k a(-r) c(k-r)`, and `c(k-r)=(k-r)_k C(k-r)`. The
falling factorial `(k-r)_k=(k-r)(k-r-1)...(k-r-k+1)` is a product of `k` consecutive
integers descending from `k-r`. For `r=1,...,k` we have `k-r in {0,...,k-1}`, so the
factor `(k-r)-(k-r)=0` occurs in the product; thus `(k-r)_k=0` and every term
vanishes. Identically `H_{k-1}[v](0)=0` since `(k-1-r)_{k-1}=0` for `r=1,...,k-1`.
∎

This is the always-present functional the finite-cokernel memo noted
([`quantum-exotic-cokernel.md`](../band3/quantum-exotic-cokernel.md)), now with the
arbitrary-`k` reason. **But `ev_0` alone never obstructs:** membership gives
`R(0)=0`, so `ev_0(E-R)=0-R(0)=0`. The obstruction must come from a *second*,
non-`ev_0` functional. That is the content of §4–§6.

---

## 4. The second functional and the residual pairing — band 3 (arbitrary degree)

For band 3 the wall `Q_5=0` is `u(E+3)a(E)=a(E+2)u(E)`. Its genuine
(non-shifted-power) integer-rooted solutions are exactly the step-2 arithmetic
progressions

```text
a_3 = (E-r)(E-r-2)(E-r-4),      u = b_2 = (E-r-1)(E-r-4)     (r in F),
```

verified by the band-3 wall scan (every genuine top is one of these). Set

```text
lambda_r(f) = f(r+3) - f(r+4) + f(r+5) - f(0).
```

> **Theorem C (band-3 AP, arbitrary degree, all `r`).** For every `r` and every
> admissible `c=(E)_3 C`, `v=(E)_2 V`,
>
> ```text
> lambda_r(K_3[c]) = 0,      lambda_r(H_2[v]) = 0,      so  Im Phi subset ker lambda_r,
> ```
>
> and `lambda_r(E) = r+4`.

**Proof.** By Theorem A' it suffices that `S_3^* lambda_r` be supported on
`roots(q_K)={0,1,2}∪{r+3,r+5,r+7}` and `S_2^* lambda_r` on
`roots(q_H)={0,1}∪{r+3,r+6}`; equivalently, direct block evaluation
([`quantum-ap-filler-image.md`](../band3/quantum-ap-filler-image.md)) gives
`K_3(r+3)=a(r+1)c(r+4)`, `K_3(r+4)=a(r+1)c(r+4)+a(r+3)c(r+6)`,
`K_3(r+5)=a(r+3)c(r+6)`, so the alternating sum at `r+3,r+4,r+5` telescopes to
`0`, and `c(0)=0` by membership; identically for `H_2`. The verifier confirms
`lambda_r(K_3[(E)_3 C])=0` and `lambda_r(H_2[(E)_2 V])=0` as **symbolic identities
in `r` and in the free coefficients of `C,V`** — arbitrary degree, all `r`. And
`lambda_r(E)=(r+3)-(r+4)+(r+5)-0=r+4`. ∎

**The residual pairing.** `lambda_r(E-R)=lambda_r(E)-lambda_r(R)=(r+4)-lambda_r(R)`.
At `r=0` (W1), the positive cascade `Q_4=Q_3=Q_2=Q_1=0` forces `lambda_0(R)=0`
(commit `e4e704f`, re-checked here on the displayed `d=2` residual), whence

```text
lambda_0(E-R) = lambda_0(E) = 4.
```

So the W1 value **`4` is `lambda_0(E) = r+4` at `r=0`, and `r+4` is the common
root of the top `a_3` and the sub-diagonal filler `u=b_2`**: both `a_3(r+4)=0` and
`b_2(r+4)=0` (checked symbolically). It is **not** `k+1`; that `4 = k+1` at `k=3`
is a numerical coincidence, disproved at bands 4–5 in §6.

---

## 5. Exact W1 recovery (sanity, matches `e4e704f`)

At `r=0`: `roots(q_K)={0,1,2,3,5,7}`, `roots(q_H)={0,1,3,6}`. The moving-sum
criterion (Theorem A') gives a 2-dimensional space of point annihilators,

```text
Ann(Im Phi) = span{ ev_0,  ev_3 - ev_4 + ev_5 } = span{ ev_0, lambda_0 },
```

with `S_3^* lambda_0 = ev_3+ev_5+ev_7-ev_0-ev_1-ev_2` (on `roots(q_K)`) and
`S_2^* lambda_0 = ev_3+ev_6-ev_0-ev_1` (on `roots(q_H)`). The verifier confirms
this point-annihilator dimension `2` equals the *stable* finite cokernel
dimension: the truncated `Phi` matrix has codim `4,3,2,2,2,...` at caps
`d=1,2,3,4,...`, converging to `2` (the finite-cokernel memo's `d=1,2` rows `4,3`
are truncation
artifacts of this limit). All of `lambda_0(E)=4`, `Im Phi subset ker lambda_0`,
`lambda_0(R)=0`, `lambda_0(E-R)=4` reproduce the `e4e704f` certificate.

---

## 6. The `k+1` question and the escape hatch (loud)

**Bands 4, 5 (bounded evidence, integer-rooted tame tops).** For the exotic wall
tops

```text
k=4:  a=E(E-2)(E-3)(E-5),           u=(E-1)(E-3)(E-5),        common roots {3,5},
k=5:  a=E(E-2)(E-4)(E-6)(E-8),      u=(E-1)(E-3)(E-6)(E-8),   common roots {6,8},
```

Theorem A' produces `3 >= 2` independent point annihilators (so `codim Im Phi >=
2`); explicit basis functionals annihilate
`Im Phi` on generic admissible `C,V`; and **no basis functional has
`lambda(E)=k+1`** (values `{0,2,3}` at `k=4`, `{0,3,4}` at `k=5`). This kills the
`lambda(E)=k+1` reading and confirms `lambda(E)` tracks *common-root* data, not
the band index. Across every integer-rooted wall-admissible band-3 top (18 of
them) `codim in {2,3}`, **never `<= 1`** — the annihilator method always has
ammunition beyond `ev_0`. The wall itself is the reason: every admissible
`(a_k,u)` shares at least one root (checked), and each shared root feeds the
transported dual construction.

**What would let `E-R in Im Phi` (the counterexample habitat).** The obstruction
is `lambda(E-R) != 0` for some `lambda in Ann(Im Phi)`. Since `ev_0(E-R)=0`
always, it must come from a non-`ev_0` functional. It fails — `E-R` escapes into
`Im Phi` — exactly when **every** functional in `Ann(Im Phi)/<ev_0>` pairs to
zero with `E-R`. Two mechanisms make that happen:

1. **Codimension collapse to 1** (only `ev_0` survives). Not observed in bands
   3–5, but not excluded arbitrary-`k`; a wall top with `codim Im Phi = 1` would
   be immune to this method. **This is the sharpest open target: hunt for a
   wall-admissible top with `codim Im Phi = 1`.**

2. **Common root at the membership anchor**, `rho* = 0`. Then the second
   functional has `lambda(E)=rho*=0`, and (with `lambda(R)=0`) `lambda(E-R)=0`.
   In the band-3 AP family this is exactly `r=-4`: there the common root `r+4`
   hits `0`, `codim` jumps to `3`, and `Im Phi = E(E-1)(E+1)F[E]` — a genuinely
   different image. **`rho* = 0` is the precise escape locus of the `lambda(E)`
   pairing; any exotic top whose top/sub common root sits at the membership point
   is where the `E`-pairing degenerates and the obstruction must be re-derived
   from the full codimension.**

**What remains for full band-`k` exotic closure.** Two gaps, both honestly open:

- `codim Im Phi >= 2` for all wall-admissible tops, arbitrary `k`, arbitrary
  degree (evidence: bands 3–5; the wall forces a shared root, but no degree-free
  proof that this yields a non-`ev_0` annihilator in general).
- `lambda(R) = 0` beyond W1. This is the one cascade-dependent step; it is proved
  only at `r=0`, `k=3`. It needs the general positive cascade
  `Q_{2k-2}=...=Q_1=0`, not pure linear algebra. Until it is proved, the pairing
  `lambda(E-R)=rho*-lambda(R)` is not known to be nonzero for `k>=4` or for
  `r != 0`.

Given both, the theorem would read: *for every band-`k` wall-admissible
non-shifted-power top with common root `rho* != 0`, `coker Phi` carries a
functional `lambda` with `Im Phi subset ker lambda` and `lambda(E-R)=rho* != 0`,
so `Q_0=1` is impossible.* Theorems A–C establish the cokernel-dual scaffold and
close the band-3 AP case for `r != -4` modulo the single residual identity
`lambda_r(R)=0` (proved at `r=0`).

---

## 7. Scope

**Proved degree-free, arbitrary `k`:** Theorem A (primal codimensions
`codim Im L_K = deg a+k`, `codim Im L_H = deg u+k-1`; `coker Phi` finite),
Theorem A' (the adjoint criterion / annihilator construction), and Theorem B
(`ev_0`, `coker Phi != 0`). **Proved arbitrary degree, band 3, all `r`:**
Theorem C (`Im Phi subset ker lambda_r`, `lambda_r(E)=r+4`). **Bounded evidence
(bands 3–5, tame integer-rooted tops):** `codim Im Phi >= 2` (via >=2 independent
point functionals), `lambda(E)=` common root, `lambda(E) != k+1`, and the point
annihilators being complete (their dimension `= ` stable finite cokernel).
**Proved only at W1 (`e4e704f`):** `lambda_0(R)=0`, hence `lambda_0(E-R)=4`.
**Open:** `codim Im Phi >= 2` and `lambda(R)=0` in general; whether a `codim=1`
wall top exists; completeness of point annihilators in general; the full band-`k`
exotic closure; DC1; JC2. No Weyl pair and no counterexample is constructed.

Run:

```sh
uv run --with sympy python research/dc1-program/verify_lambda_general_k.py
```

Ends `ALL LAMBDA GENERAL K CHECKS PASSED`.
