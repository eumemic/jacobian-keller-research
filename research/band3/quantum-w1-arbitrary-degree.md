# Arbitrary-degree unit obstruction for the normalized W1 wall datum

**INDEPENDENTLY DERIVED — EXACT ALGEBRA — NOT PEER REVIEWED — BAND-SCOPED**

This note proves one narrowly normalized statement over a characteristic-zero field `F`: for the fixed W1 wall datum below, no polynomial cascade satisfying `Q_4=Q_3=Q_2=Q_1=0` and genuine Weyl membership can also satisfy `Q_0=1`, regardless of the coefficient degrees. It does not classify scalar multiples or other orientations of W1, close the full arithmetic-progression family, W2, Band 3, DC1, or JC2, and it constructs no Weyl pair or counterexample. The exact certificate is [`verify_quantum_w1_arbitrary_degree.py`](verify_quantum_w1_arbitrary_degree.py).

## 1. Arbitrary-degree written proof

Work in `F[E]`, with `char(F)=0`. Use the stipulated convention

```text
Q_m = sum_(k+l=m) (b_l^[k] a_k - a_k^[l] b_l),
f^[r](E)=f(E+r),
```

with gauge `b_3=0` and Weyl membership `(E)_j | a_-j,b_-j`. Fix

```text
a_3=E(E-2)(E-4),    b_2=(E-1)(E-4).
```

These satisfy the `Q_5` wall. Suppose the full positive cascade
`Q_4=Q_3=Q_2=Q_1=0` holds as polynomial identities. Split the central potential as

```text
G=R+K_3[b_-3]-H_2[a_-2],
R=a_1^[-1]b_-1-b_1^[-1]a_-1
  +a_2^[-2]b_-2+a_2^[-1]b_-2^[1].
```

Here `R` is the canonical non-filler residual. Membership gives
`a_-1(0)=b_-1(0)=0` and `b_-2(0)=b_-2(1)=0`, hence `R(0)=0`.
Define

```text
lambda_0(f)=f(3)-f(4)+f(5)-f(0).
```

After membership and cancellation of the repeated middle terms,
`lambda_0(R)` is the following eight boundary evaluations:

```text
 a_1(2)b_-1(3)-b_1(2)a_-1(3)+a_2(1)b_-2(3)
-a_1(3)b_-1(4)+b_1(3)a_-1(4)
+a_1(4)b_-1(5)-b_1(4)a_-1(5)+a_2(4)b_-2(6).
```

These are exactly eight signed boundary products; the two repeated interior level-2 products cancel.

Five evaluations of `Q_4=0` give

```text
b_1(0)=-2/3 a_2(0)-4/3 a_2(2),
b_1(2)=-2/3 a_2(4),
b_1(4)= 2/3 a_2(1),
b_1(6)= 4/3 a_2(3)+2/3 a_2(5),
a_2(4)=a_2(1).
```

Next use, in order, `Q_3(0),Q_3(1),Q_3(3),Q_3(4)`, then
`Q_2(0),...,Q_2(3)`, and finally `Q_1(0),...,Q_1(3)`. Every solved evaluation has a fixed nonzero rational pivot; no degree bound or leading-coefficient assumption enters. After those substitutions, the unused equation is

```text
Q_3(2)=-1/3(a_2(0)-a_2(5))(a_2(1)-a_2(4)),
```

so it vanishes by the `Q_4` compatibility relation. Exact nonlinear reduction of the boundary expression itself gives

```text
lambda_0(R) congruent to 0
```

modulo those equations, membership, and `a_2(1)-a_2(4)`. Therefore

```text
lambda_0(R)=0.
```

The sequential reconstruction retains its arbitrary constant `b_0` operator kernel. A constant shift of `b_0` has cancelling shifted contributions in `Q_4,Q_3,Q_2,Q_1`; moreover `b_0` does not occur in `G`, and therefore does not enter `R`. The argument does not silently set this kernel freedom to zero.

It remains to connect this boundary identity to `Q_0`. Full negative-level Weyl membership writes

```text
b_-3=(E)_3 C,    a_-2=(E)_2 V
```

and gives `G(0)=0`. The displayed decomposition is therefore

```text
G=R+Phi(C,V),
Phi(C,V)=K_3[(E)_3 C]-H_2[(E)_2 V].
```

The argument allows arbitrary admissible `C,V`; it does not use bottom proportionality or any lower negative equation. Since `Q_0=(T-1)G` and `(T-1)E=1`, the equation `Q_0=1` makes `G-E` 1-periodic. Over characteristic zero a nonconstant polynomial of degree `n>0` cannot be 1-periodic: the leading term of `p(E+1)-p(E)` is `n lc(p) E^{n-1}`, which is nonzero. Hence `G=E+c`; membership gives `G(0)=0`, so `c=0` and `G=E`. Thus `Q_0=1` would require

```text
E-R in Im Phi.
```

The needed filler statement is proved directly here. For every admissible
`c=(E)_3C` and `v=(E)_2V`, direct substitution of the displayed W1 polynomials into the two filler blocks gives

```text
lambda_0(K_3[c])=0,    lambda_0(H_2[v])=0,
```

because `c(0)=c(1)=c(2)=0`, `v(0)=v(1)=0`, and the remaining values cancel pairwise at `3,4,5`. Consequently

```text
Im Phi subset ker(lambda_0).
```

This inclusion is all the obstruction needs. The separately proved unrestricted filler-image theorem in [`quantum-ap-filler-image.md`](quantum-ap-filler-image.md) gives the stronger equality `Im Phi = ker(ev_0) intersect ker(lambda_0)` at W1, but that result is not load-bearing here.

Now `R(0)=0`, `lambda_0(R)=0`, and

```text
lambda_0(E)=3-4+5-0=4.
```

Therefore `lambda_0(E-R)=4`, so `E-R` is not in `Im Phi`. Hence no polynomial cascade for this normalized W1 wall datum can satisfy `Q_4=...=Q_1=0`, Weyl membership, and `Q_0=1`, regardless of coefficient degrees.

## 2. Exact algebra certificate

The verifier constructs every `Q_m` from the stipulated summation formula, checks the W1 wall, `b_3=0`, polynomial Weyl membership, and `Q_0=(T-1)G`. It then works in exact independent evaluation variables, imposes every boundary membership equation, and includes a negative guard showing that omitting `b_-2(1)=0` destroys the `R(0)` conclusion.

It reproduces all four displayed `b_1` evaluations, `a_2(4)=a_2(1)`, the prescribed `Q_3,Q_2,Q_1` evaluation sequence with fixed rational pivots, the exact residual form of `Q_3(2)`, and the final normal form

```text
lambda_0(R)=0.
```

These are exact symbolic identities modulo the generated cascade and membership relations, not sampled degree calculations. The verifier also checks directly the decomposition `G=R+K_3[b_-3]-H_2[a_-2]`, membership's implication `G(0)=0`, both unrestricted evaluation-variable annihilators `lambda_0(K_3[c])=lambda_0(H_2[v])=0`, `(T-1)E=1`, the symbolic arbitrary-`n` leading coefficient `n` in `(E+1)^n-E^n`, and the elimination of the additive constant by `G(0)=0`. Thus the decisive `Im Phi subset ker(lambda_0)` bridge is checked inside the W1 certificate; the stronger separate image theorem is only corroboration.

## 3. Scope

**Proved here:** for the normalized W1 wall datum displayed above over a characteristic-zero field, no polynomial data of any coefficient degrees can satisfy the positive cascade, genuine Weyl membership, and `Q_0=1` simultaneously.

**Not proved here:** scalar-normalization or orientation reductions for every datum called W1, the unrestricted AP family, W2, the whole exotic branch, Band 3, DC1, or JC2. No Weyl pair and no counterexample is constructed.

Run:

```sh
python research/band3/verify_quantum_w1_arbitrary_degree.py
```

A successful run ends with `ALL QUANTUM W1 ARBITRARY-DEGREE CHECKS PASSED`.
