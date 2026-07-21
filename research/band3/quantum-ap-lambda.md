# Arbitrary-degree unit obstruction across the band-3 degree-3 AP exotic family

**INDEPENDENTLY DERIVED — EXACT ALGEBRA — NOT PEER REVIEWED — BAND-SCOPED**

This note extends the fixed normalized-W1 theorem
([`quantum-w1-arbitrary-degree.md`](quantum-w1-arbitrary-degree.md), overnight
commit `e4e704f`) from the single member `r=0` to the **entire** step-2
arithmetic-progression exotic family, with the translation parameter `r` kept
symbolic. Over a characteristic-zero field `F`, for the wall datum

```text
a_3 = (E-r)(E-r-2)(E-r-4),   b_2 = (E-r-1)(E-r-4),   b_3 = 0
```

it proves: **for every `r != -4`, no polynomial cascade of any coefficient
degrees can satisfy `Q_4=Q_3=Q_2=Q_1=0`, genuine Weyl membership, and `Q_0=1`
simultaneously.** The single value `r=-4` is exactly the W2 member; there the
obstructing functional degenerates and the question reduces to the moment-unit
slope. The exact certificate is
[`verify_quantum_ap_lambda.py`](verify_quantum_ap_lambda.py) (ends
`ALL QUANTUM AP LAMBDA CHECKS PASSED`).

It does **not** close `r=-4` at arbitrary degree, higher-degree/non-AP Band-3
tops, DC1, or JC2, and it constructs no Weyl pair or counterexample. A later
bounded verdict, [`w2-verdict.md`](w2-verdict.md), excludes the encoded W2
combined system at raw cap `d=3`; an exact `d=4` exclusion is externally reported
and optionally reproducible. Neither enlarges this theorem's arbitrary-degree
scope.

## 0. Convention and the family

Retain the stipulated convention (identical to the W1 memo)

```text
Q_m = sum_(k+l=m) (b_l^[k] a_k - a_k^[l] b_l),   f^[n](E)=f(E+n),
G   = sum_(k=1)^3 sum_(j=0)^(k-1) (a_k^[j-k] b_-k^[j] - b_k^[j-k] a_-k^[j]),
Q_0 = (T-1)G,
```

with gauge `b_3=0` and Weyl membership `(E)_j | a_-j, b_-j`,
`(E)_j = E(E-1)...(E-j+1)`. The AP top `{r,r+2,r+4}` solves the `Q_5` wall
`b_2^[3] a_3 = a_3^[2] b_2` for **all** symbolic `r`
(`quantum-exotic-closure.md`; re-checked here). Its two named members are

```text
W1 = member r=0:   a_3=E(E-2)(E-4),   b_2=(E-1)(E-4),
W2 = member r=-4:  a_3=E(E+2)(E+4),   b_2=E(E+3).
```

After solving the positive cascade, split the central potential exactly as in
the W1 memo,

```text
G = R + K_3[b_-3] - H_2[a_-2],
R = a_1^[-1] b_-1 - b_1^[-1] a_-1 + a_2^[-2] b_-2 + a_2^[-1] b_-2^[1],
K_3[c] = sum_(j=0)^2 a_3^[j-3] c^[j],   H_2[v] = sum_(j=0)^1 b_2^[j-2] v^[j],
Phi(C,V) = K_3[(E)_3 C] - H_2[(E)_2 V].
```

`R` is the canonical non-filler residual; `Im Phi` is the admissible filler
image. Membership gives `G(0)=0`, and `Q_0=1` is equivalent to `G=E` (the
`1`-periodicity argument of the W1 memo, unchanged), hence `Q_0=1` forces

```text
E - R = Phi(C,V) in Im Phi.
```

The obstruction is an annihilating functional of `Im Phi` that is nonzero on
`E-R`.

## 1. The symbolic-r annihilating functional

Define, for `r in F`,

```text
lambda_r(f) = f(r+3) - f(r+4) + f(r+5) - f(0).
```

**Filler annihilation (holds for every `r`).** Because
`a_3(r+j) = j(j-2)(j-4)` and `b_2(r+j) = (j-1)(j-4)` are `r`-independent
constants, every summand of `Phi` vanishes at `E=0`, and for `c=(E)_3 C`,
`v=(E)_2 V`,

```text
K_3[c](r+3)=3 c(r+4),   K_3[c](r+4)=3 c(r+4)-3 c(r+6),      K_3[c](r+5)=-3 c(r+6),
H_2[v](r+3)=-2 v(r+4),  H_2[v](r+4)=-2 v(r+4)-2 v(r+5),     H_2[v](r+5)=-2 v(r+5),
```

so the alternating values cancel block by block for `C,V` of arbitrary degree.
The verifier checks these identities symbolically in `r` on generic degree-5
polynomials; the displayed nodewise calculation is the degree-free proof. This is
exactly the inclusion of the unrestricted image theorem
[`quantum-ap-filler-image.md`](quantum-ap-filler-image.md):

```text
Im Phi(r) subset ker(ev_0) intersect ker(lambda_r),      lambda_r(E) = r+4.
```

This inclusion — all the obstruction needs — is `r`-uniform and degree-free.

## 2. Arbitrary-degree written proof that `lambda_r(R) = 0`

Membership gives `a_-1(0)=b_-1(0)=0` and `b_-2(0)=b_-2(1)=0`, hence `R(0)=0` for
every `r`. After the two interior level-2 products cancel,

```text
lambda_r(R) = R(r+3)-R(r+4)+R(r+5)
```

is the eight signed boundary products

```text
 a_1(r+2)b_-1(r+3) - b_1(r+2)a_-1(r+3) + a_2(r+1)b_-2(r+3)
-a_1(r+3)b_-1(r+4) + b_1(r+3)a_-1(r+4)
+a_1(r+4)b_-1(r+5) - b_1(r+4)a_-1(r+5) + a_2(r+4)b_-2(r+6).
```

**The decisive structural fact:** the wall values at the `r`-shifted nodes
`r+j` are `r`-independent integers. Evaluating the positive cascade at those
nodes therefore gives elimination steps with **fixed rational pivots that do not
depend on `r`**:

```text
Q_4(r+0),(r+2),(r+1),(r+3):  solve b_1(r+0),(r+2),(r+4),(r+6)  pivots -3, 3, 3, -3;
Q_4(r+4):                    the relation a_2(r+4) = a_2(r+1)   (pivot -10);
Q_3(r+0),(r+1),(r+3),(r+4):  solve a_1(r+2),(r+1),(r+5),(r+4)  pivots -4, -2, 2, 4;
Q_3(r+2):                    residual = -1/3 (a_2(r)-a_2(r+5))(a_2(r+1)-a_2(r+4));
Q_2(r+0),(r+1),(r+2),(r+3):  solve a_0(r+0), b_-1(r+4), a_0(r+4), a_0(r+5)  pivots 4,3,2,2;
Q_1(r+0),(r+1),(r+2),(r+3):  solve a_-1(r+2), b_-2(r+4), a_-1(r+4), a_-1(r+5) pivots -4,3,2,2.
```

Every pivot is a nonzero integer, **independent of `r`**; the whole elimination
is the `r`-translate of the W1 elimination with the same pivots. Substituting the
solved values into the eight boundary products and reducing modulo the single
relation `a_2(r+1)=a_2(r+4)` gives

```text
lambda_r(R) = 0     for EVERY r.
```

No absolute-node membership is used inside this reduction: the eight-term
expression and every relation used live entirely among the `r`-shifted
evaluations, which are treated as independent (this is what makes the argument
degree-free — exactly the mechanism of the W1 certificate, now carried out with
`r` symbolic). The retained constant-`b_0` operator kernel does not enter `G` or
`R`, as in the W1 memo.

## 3. The obstruction for every `r != -4`

Since `lambda_r(E) = (r+3)-(r+4)+(r+5)-0 = r+4` and `lambda_r(R)=0`,

```text
lambda_r(E - R) = r + 4.
```

If a cascade with `Q_0=1` existed, then `E-R in Im Phi(r) subset ker(lambda_r)`
would force `lambda_r(E-R)=0`, i.e. `r+4=0`. Therefore, **for every `r != -4`,
no polynomial data of any coefficient degrees can satisfy the positive cascade,
Weyl membership, and `Q_0=1` for the AP wall `{r,r+2,r+4}`.** At `r=0` the value
is `4`, reproducing the committed W1 theorem.

### Exceptional locus

The elimination pivots are nonzero constants for every `r`, and the filler
inclusion holds for every `r`; the only place the argument can fail is
`lambda_r(E)=r+4=0`. Hence the exceptional set is **exactly** `{ r = -4 }`. The
AP wall is a genuine (non-shifted-cube) exotic top for every `r` — its three
roots `r, r+2, r+4` are always distinct — so there is no degenerate-wall or
shifted-cube coincidence to separate out; `r=-4` is the sole exception.

## 4. The exceptional member `r=-4` (= W2): honest ledger

At `r=-4` the functional degenerates: `lambda_{-4}(E)=0`, so
`lambda_{-4}(E-R)=0` yields no contradiction. Correspondingly the filler image
jumps to codimension three,

```text
Im Phi(-4) = D * F[E],   D = E(E-1)(E+1)
```

(`quantum-ap-filler-image.md`, by its degree-free triangular proof; the lambda
verifier performs a finite basis regression only). The annihilators of `D*F[E]`
are exactly `ev_0, ev_1, ev_-1`, and
`ev_0(E-R)=-R(0)=0` automatically. The `r=-4` specialization of the Section-2
certificate still rigorously gives

```text
R(1) + R(-1) = 0.
```

Using `R(-1) = -R(1)`, every annihilator `alpha ev_0 + beta ev_1 + gamma ev_-1`
sends `E-R` to `(beta - gamma)(1 - R(1))`. But at `r=-4` the two filler blocks
vanish at `E=1`, so

```text
R(1) = G(1) = slope = Q_0(0).
```

Thus `E-R in Im Phi(-4)` is equivalent to `R(1)=1`, i.e. the moment slope equals
`1` — which is **consistent** with `Q_0=1` (then `Q_0(0)=1`). The functional
method that annihilates the `r != -4` obstruction therefore produces no
contradiction at `r=-4`: its relaxed shifted-evaluation normal form pins the
alternating combination `R(1)+R(-1)` to `0` but does not force `R(1)` itself to
zero. That non-forcing statement alone does not prove value-one achievability or
a free modulus in the actual polynomial cascade; the explicit family in
[`w2-decisive.md`](w2-decisive.md) supplies that proof. The arbitrary-degree
closure at `r=-4` remains **open**.

**Finite-degree closure at `r=-4`.** The positive cascade forces the slope
`R(1)=0 != 1` at `d=1` (reproduced compactly in the verifier), and W2 is
infeasible at `d=2` as well (the committed
[`verify_quantum_exotic_cokernel.py`](verify_quantum_exotic_cokernel.py) Gröbner
branch certificate: on both branches `R in D*F[E]`, so `[E-R] != 0` in the finite
cokernel). These are bounded-degree statements only; they do not imply an
arbitrary-degree obstruction at `r=-4`.

## 5. Combined statement and scope

**Proved here (arbitrary coefficient degree, characteristic zero):** the
**normalized degree-3 step-2 AP exotic family** with `b_2 != 0` and top
`{r,r+2,r+4}` is empty — no cascade with the positive wall conditions, genuine
Weyl membership, and `Q_0=1` — for **every `r != -4`**. This is not an
arbitrary-degree exclusion of every Band-3 exotic top. The obstruction is the
single explicit functional `lambda_r(f)=f(r+3)-f(r+4)+f(r+5)-f(0)`, with
`Im Phi(r) subset ker(lambda_r)` and `lambda_r(E-R)=r+4 != 0`.

**Explicit exceptional locus within this normalized AP family:** `r=-4` (W2).
There the functional degenerates and does not resolve the moment-unit slope. W2
is closed at raw caps `d<=2` by earlier finite computation; the later bounded
verdict excludes the encoded combined system at raw cap `d=3`, while exact `d=4`
is externally reported and optionally reproducible. Its arbitrary-degree status
remains open.

**Out of scope / not proved:** the non-AP `deg a_3 >= 6` exotic tops; scalar or
orientation reductions of any datum; the arbitrary-degree `r=-4` case; all of
Band 3, DC1, JC2. In unrestricted degree `Im L_K intersect Im L_H` remains
infinite-dimensional (`../dc1-program/two-filler-cross-cancellation.md`); nothing
here weakens that.

Run:

```sh
python research/band3/verify_quantum_ap_lambda.py
```

A successful run ends with `ALL QUANTUM AP LAMBDA CHECKS PASSED`.
