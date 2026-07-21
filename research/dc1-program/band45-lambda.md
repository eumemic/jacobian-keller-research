# The annihilator template at bands 4 and 5: the cofactor functional

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER REVIEWED — BAND-SCOPED**

This note applies the overnight W1 annihilator template
([`../band3/quantum-w1-arbitrary-degree.md`](../band3/quantum-w1-arbitrary-degree.md),
commit `e4e704f`) to the **band-4** and **band-5** minimal exotic tops. It answers
the structural question — *what is `lambda_0` intrinsically?* — with a closed form
that is uniform over every band and every exotic top, and it upgrades the bounded
band-4/5 moment-unit experiments
([`band4-moment-unit-experiment.md`](band4-moment-unit-experiment.md),
[`band5-comparison.md`](band5-comparison.md), corrected at `449cca2`) from
per-ansatz Gröbner kills to a degree-free obstruction on the filler image, modulo
one honestly-open cascade step.

The companion adjoint analysis is the general-`k` sibling
[`lambda-general-k.md`](lambda-general-k.md) (its Theorem A' is the moving-sum
criterion that the closed form below solves in closed form). This note does not
construct a Weyl pair, settle DC1/JC2, or close band 4 or band 5.

Exact checks: [`verify_band45_lambda.py`](verify_band45_lambda.py) — 132 checks,
ends `ALL BAND45 LAMBDA CHECKS PASSED`.

Conventions are frozen from the corpus (`b7e85e8..e4e704f`): `F` characteristic
zero, `A_1[x^{-1}]=(+)_k x^k C[E]`, `(x^a f)(x^b g)=x^{a+b} f(E+b)g(E)`,
`T f(E)=f(E+1)`, `S_n=1+T+...+T^{n-1}`, `(E)_j=E(E-1)...(E-j+1)`, membership
`(E)_j | a_{-j},b_{-j}`, `Q_m=sum_{i+l=m}(b_l^{[i]}a_i-a_i^{[l]}b_l)`, gauge
`b_k=0`, `u:=b_{k-1}`. Necklaces are in `S` (the shift on root positions):
`A(S)=sum_rho S^rho`, `S_k=1+S+...+S^{k-1}`.

---

## 1. The map and the decomposition (degree-free)

After the top Wronskian gauge `b_k=0`, solving the positive cascade splits the
closure-form central potential
`G = sum_{k=1}^K sum_{j=0}^{k-1}(a_k^{[j-k]}b_{-k}^{[j]}-b_k^{[j-k]}a_{-k}^{[j]})`
as

```text
G = R + K_k[b_-k] - H_{k-1}[a_-(k-1)],
K_k[c]     = sum_{j=0}^{k-1} a_k(E+j-k) c(E+j),           a_k := top,
H_{k-1}[v] = sum_{j=0}^{k-2} b_{k-1}(E+j-(k-1)) v(E+j),   b_{k-1} := u,
R          = sum_{i=1}^{k-1} K_i[b_-i] - sum_{i=1}^{k-2} H_i[a_-i]   (interior blocks).
```

Under membership `c=b_-k=(E)_k C`, `v=a_-(k-1)=(E)_{k-1} V`, the two filler blocks
are `Phi(C,V)=K_k[(E)_k C]-H_{k-1}[(E)_{k-1} V]`. Because `Q_0=(T-1)G` and
`(T-1)E=1`, and membership gives `G(0)=0`,

```text
Q_0 = 1   <=>   G = E   <=>   E - R in Im Phi.
```

The verifier checks `G = R + K_k[b_-k] - H_{k-1}[a_-(k-1)]` as an exact polynomial
identity in generic symbolic coefficients at `k=4` and `k=5` (`R` derived here, not
imported). This mirrors the band-3 split `G = R + K_3[b_-3] - H_2[a_-2]`
([`../band3/quantum-exotic-cokernel.md`](../band3/quantum-exotic-cokernel.md)).

---

## 2. What `lambda_0` is: the cofactor functional (the emerging pattern)

Every minimal exotic top is `S_k`-divisible; its **necklace cofactor** is

```text
C(S) = A(S) / S_k        (deg A = max root; the NON-effective quotient that
                          certifies exoticness -- e.g. Phi_6, Phi_10, S^3-S^2+1).
```

> **The pattern.** The annihilating functional is the *cofactor boundary
> functional*
>
> ```text
> hat_lambda(t) = t^k C(t) - C(1),        lambda = sum_p ([t^p] hat_lambda) ev_p.
> ```
>
> It is determined **entirely by the cofactor `C`** — no cascade data, no degree
> bound.

This is not a residue at a pole and not a raw coefficient sum: it is the boundary
pairing dual to the two shifted root necklaces. Concretely, `lambda` is pinned by
the **two generating-function support identities** (both proved exactly for all
seven tops):

```text
hat_lambda * S_k     = t^k A       - C(1) S_k        (K-block:  supp {0..k-1} u {k+rho})
hat_lambda * S_{k-1} = t^{k-1} B   - C(1) S_{k-1}    (H-block:  supp {0..k-2} u {k-1+sigma})
```

where `B(S)=t S_{k-1} C` is the wall root-necklace of `u=b_{k-1}` (the necklace form
of the wall `Q_{2k-1}=0`, [`moment-unit-general-k.md`](moment-unit-general-k.md) §2).
These say exactly that the forward moving-sum `S_k^* lambda` lands on the shifted
top necklace `{k+rho}` and `S_{k-1}^* lambda` on the shifted wall necklace
`{k-1+sigma}` — the moving-sum / adjoint criterion of
[`lambda-general-k.md`](lambda-general-k.md) Theorem A', now solved in closed form.

**Why it annihilates (degree-free).** With `c=(E)_k C_poly`,

```text
lambda(K_k[c]) = sum_n c(n) a_k(n-k) [t^n](hat_lambda * S_k)
              = sum_n c(n) a_k(n-k) ( [n = k+rho] - C(1)[0<=n<=k-1] ).
```

The `[0<=n<=k-1]` terms die because `c(n)=(n)_k C_poly(n)=0` there (membership); the
`[n=k+rho]` terms die because `a_k(n-k)=a_k(rho)=0` (the top's roots). Identically
`lambda(H_{k-1}[v])=0` uses `v(0..k-2)=0` and `b_{k-1}(sigma)=0`. Hence

```text
Im Phi subset ker lambda      (proved, arbitrary degree, every band-4/5 top).
```

The verifier confirms `lambda(K_k[(E)_k C])=0` and `lambda(H_{k-1}[(E)_{k-1} V])=0`
as **symbolic identities in the free coefficients of generic degree-4 `C,V`** — the
degree-free content, not sampled.

### 2.1 The two scalar invariants

```text
C(1) = A(1)/S_k(1) = k/k = 1     =>   lambda(1) = hat_lambda(1) = C(1) - 1 = 0,

lambda(E) = hat_lambda'(1) = k C(1) + C'(1) = (k+1)/2 + (sum rho)/k.
```

`C'(1) = (sum rho)/k - (k-1)/2` because `A'(1)=sum rho` and `S_k'(1)=k(k-1)/2`.
For every **normalized** top (roots `>= 0`, so `sum rho >= 0`),
`lambda(E) >= (k+1)/2 > 0`. Therefore

```text
E not in Im Phi        (proved, arbitrary degree, every normalized exotic top).
```

`lambda(E)` is thus a *weighted mean of the root positions*, `(k+1)/2 + mean(rho)`.
The seven computed values:

| band | top | cofactor `C` | `lambda(E)` |
|---|---|---|---:|
| 4 | `{0,2,3,5}` | `Phi_6 = S^2-S+1` | `5` |
| 4 | `{0,1,3,6}` | `S^3-S^2+1` | `5` |
| 4 | `{0,3,5,6}` | `S^3-S+1` | `6` |
| 4 | `{0,3,6,9}` | `Phi_6 Phi_12` | `7` |
| 5 | `{0,2,3,4,6}` | `Phi_6` (universal `1-s+s^2`) | `6` |
| 5 | `{0,2,4,6,8}` | `Phi_10` | `7` |
| 5 | `{0,1,3,4,7}` | `S^3-S^2+1` (**non-cyclotomic**) | `6` |

The always-present functional `ev_0` (the other annihilator, from `c(0)=v(0)=0`)
has `ev_0(E)=0`, so it never obstructs `E-R`; the obstruction is carried by the
cofactor functional. Two remarks reconcile this with the sibling note's reading:

- The annihilator space `Ann(Im Phi)` is `>= 2`-dimensional; `lambda(E)` is not a
  single number but a functional on that space. The sibling's canonical *second*
  functional has `lambda(E) = `(a common root of `a_k` and `u`); the **cofactor**
  functional here is a different, canonical element with the closed form
  `(k+1)/2 + mean(rho)`. For the two self-dual band-4 tops they coincide with a
  common root (`5 in {3,5}`, `7 in {...}`); for `{0,1,3,6}` they differ
  (`5 not in {1,6}`). Both are honest annihilators; the cofactor one is the uniform
  closed form that generalizes.
- `lambda(E) = k+1` is a band-3 coincidence (`4 = k+1` at `k=3`), already disproved
  in [`lambda-general-k.md`](lambda-general-k.md) §6; the correct reading is the
  root-data mean above.

---

## 3. Reflection `E -> -E-1` (halves the work)

Reflection sends a normalized top with roots `rho`, max `M`, to the top with roots
`{M-rho}`. Its cofactor is the **reversal** of `C`:

```text
C_refl(t) = t^{deg C} C(1/t),        lambda_refl(E) = (k+1)/2 + M - mean(rho).
```

Checked exactly for all seven tops. Consequences used to halve the work:

- **Self-dual tops** (`{0,2,3,5}`, `{0,3,6,9}`, `{0,2,3,4,6}`, `{0,2,4,6,8}`) have
  **palindromic** `C`, so `lambda` is reflection-symmetric and `lambda_refl(E)=lambda(E)`.
- The band-4 **reflection pair** `{0,1,3,6} <-> {0,3,5,6}` maps by reversal
  (`S^3-S^2+1 <-> S^3-S+1`), and `lambda(E)` swaps `5 <-> 6`. Analysing one member
  gives the other for free.

---

## 4. The non-cyclotomic test (the sharpest case)

`{0,1,3,4,7}` at `k=5` has cofactor `C = S^3 - S^2 + 1`, which is **not any
cyclotomic `Phi_n`** (verified against `Phi_1..Phi_59`). The cofactor functional,
its two support identities, `Im Phi subset ker lambda`, and `lambda(E)=6>0` all go
through **identically**. Its reflection `{0,3,4,6,7}` has cofactor `S^3-S+1`, also
non-cyclotomic. Conclusion: **the annihilator mechanism does not depend on
cyclotomic structure** — it is a formal consequence of `A = S_k C` and the wall
necklace `B = t S_{k-1} C`, whatever `C` is.

---

## 5. The residual `lambda(R)` — bounded kill, one open step

The full moment-unit obstruction needs `E - R not in Im Phi`, i.e.
`lambda(E-R)=lambda(E)-lambda(R) != 0`. Since `lambda(E) != 0` is settled degree-free
(§2.1), everything reduces to the single scalar `lambda(R)`.

Because `lambda` annihilates the two filler blocks and the cascade's `b_-k,a_-(k-1)`
satisfy membership, `lambda(G) = lambda(R)`. The verifier forward-solves the
positive cascade `Q_{2k-2}=...=Q_1=0` at free degree `d` and reduces `lambda(G)`
modulo the positive-cascade ideal:

```text
for every band-4/5 top:  lambda(G) = lambda(R)  reduces to the CONSTANT 0
   (d=1 all seven tops; d=2 for {0,2,3,5},{0,3,6,9},{0,2,3,4,6}),
while the cascade ideal is proper (feasible), so the kill is real, not vacuous.
```

Hence at those degrees `lambda(E-R) = lambda(E) in {5,6,7} != 0`, so `Q_0=1` is
infeasible — the moment-unit kill, now pinned to a single functional and sharper
than the full `Q_0=1` Gröbner test of the prior experiments.

**This is bounded (per-degree) evidence, not an arbitrary-degree theorem.** The
residual identity `lambda(R)=0` is the one cascade-dependent step; it is proved
degree-free only at band-3 W1 (`e4e704f`) and, for the whole band-3 AP family, in
the sibling `quantum-ap-lambda.md`. At `k>=4` it is not pure `b`-linear algebra
(the interior `a_i` are themselves cascade-constrained; a bilinear rank test of
`lambda(R) in span{Q_m(p)}` over `Q(a-evals)` is *inconsistent* at finite windows,
so the reduction genuinely needs the joint nonlinear cascade), exactly the escape
hatch flagged in [`lambda-general-k.md`](lambda-general-k.md) §6. Nothing here
proves it beyond the stated `d`.

---

## 6. Scope

**Proved degree-free, arbitrary degree, every band-4/5 exotic top (and, by the
same argument, every `S_k`-divisible top at every band):**

1. `G = R + K_k[b_-k] - H_{k-1}[a_-(k-1)]` (exact identity).
2. the cofactor functional `hat_lambda = t^k C - C(1)` and its two support
   identities `hat_lambda*S_k = t^k A - C(1)S_k`, `hat_lambda*S_{k-1} = t^{k-1}B -
   C(1)S_{k-1}`;
3. `Im Phi subset ker lambda` (symbolic block annihilation, any degree);
4. `C(1)=1`, `lambda(1)=0`, `lambda(E)=(k+1)/2+mean(rho)>0`, hence `E not in Im Phi`;
5. reflection covariance `C_refl = reversal(C)`;
6. cyclotomic-independence (the non-cyclotomic `{0,1,3,4,7}`).

**Bounded evidence (per-degree Gröbner):** `lambda(R)=0` — the full kill
`E-R not in Im Phi`, hence `Q_0=1` impossible — at `d=1` (all seven tops) and `d=2`
(three representatives).

**Open:** `lambda(R)=0` at arbitrary degree for `k>=4` (the single cascade-dependent
step; would upgrade §5 to a band-4/5 arbitrary-degree obstruction); `codim Im Phi
>= 2` in general; whether any wall top has `codim Im Phi = 1` or common root at the
membership anchor `0` (the escape loci); full band-4/5 exotic closure; DC1; JC2. No
Weyl pair and no counterexample is constructed.

Run:

```sh
uv run --with sympy python research/dc1-program/verify_band45_lambda.py
```

Ends `ALL BAND45 LAMBDA CHECKS PASSED`.
