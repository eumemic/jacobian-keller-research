# Unrestricted filler image for the quantum degree-three AP top

**INDEPENDENTLY DERIVED — EXACT ALGEBRA — NOT PEER REVIEWED — BAND-SCOPED**

This note classifies only the unrestricted polynomial image of one linear filler
map. It does not prove a full cascade obstruction. In an application to the
Band-3 cascade, the cascade-produced residual target `E-R` still has to be tested
for membership in this image. Thus the theorem below does not close Band 3,
DC1, or JC2, and it constructs no Weyl pair or counterexample.

Exact finite certificates and proof regressions are in
[`verify_quantum_ap_filler_image.py`](verify_quantum_ap_filler_image.py).

## 1. Convention and theorem

Let `F` be any characteristic-zero field and let `r in F`. Write

```text
(E)_m = E(E-1)...(E-m+1),
K_3(a,c) = sum_(j=0)^2 a(E+j-3)c(E+j),
H_2(b,v) = sum_(j=0)^1 b(E+j-2)v(E+j),
Phi(C,V) = K_3(a_3,(E)_3 C)-H_2(b_2,(E)_2 V),
a_3=(E-r)(E-r-2)(E-r-4),
b_2=(E-r-1)(E-r-4).
```

**Theorem (unrestricted image of `Phi`).**

1. If `r != -4`, then

   ```text
   Im Phi = {f in F[E] : f(0)=0 and
             f(r+3)-f(r+4)+f(r+5)-f(0)=0}.
   ```

   In particular, `Im Phi` has codimension `2` in `F[E]`.

2. If `r=-4`, then

   ```text
   Im Phi = E(E-1)(E+1) F[E],
   ```

   of codimension `3` in `F[E]`.

These alternatives classify `Phi` only. They are not an unrestricted exception
classification for the nonlinear cascade.

## 2. Proof when `r != -4`

Define

```text
lambda_r(f)=f(r+3)-f(r+4)+f(r+5)-f(0).
```

Every summand in `Phi` vanishes at `E=0`: in the `K_3` block this follows from
the roots at `0,1,2` of the shifted copies of `(E)_3`, and in the `H_2` block
from the roots at `0,1` of the shifted copies of `(E)_2`. For
`c=(E)_3C` and `v=(E)_2V`, direct evaluation at the other three points gives

```text
K_3(a_3,c)(r+3)=a_3(r+1)c(r+4),
K_3(a_3,c)(r+4)=a_3(r+1)c(r+4)+a_3(r+3)c(r+6),
K_3(a_3,c)(r+5)=a_3(r+3)c(r+6),

H_2(b_2,v)(r+3)=b_2(r+2)v(r+4),
H_2(b_2,v)(r+4)=b_2(r+2)v(r+4)+b_2(r+3)v(r+5),
H_2(b_2,v)(r+5)=b_2(r+3)v(r+5).
```

Thus the alternating values cancel separately for both blocks, so evaluation at
zero and `lambda_r` annihilate `Im Phi`. They are independent when `r != -4`,
because

```text
lambda_r(E)=r+4,
```

whereas evaluation at zero sends `E` to zero.

It remains to prove that there are no further conditions. Use coefficient rows
and source columns numbered from zero, with source order

```text
(C_0,...,C_d,V_0,...,V_d)
```

and target order `(1,E,...,E^(d+6))`. At the exact base cap `d=4`, two `9 x 9`
minors are

```text
D1, rows {2,3,4,5,6,7,8,9,10}, cols {0,1,2,3,4,5,6,7,9}:
  -648(r+4)^2(r+7)(3r^2+18r+14),

D2, rows {1,3,4,5,6,7,8,9,10}, cols {0,1,3,4,5,6,7,8,9}:
  -648(r+4)(r^2+8r+6)(r^2+8r+18).
```

Their gcd is `648(r+4)`. More explicitly, with `A=D1/(r+4)` and
`B=D2/(r+4)`, the exact Bézout certificate is

```text
(15r^3+58r^2-406r-2176) A
 -(45r^3+219r^2-1284r-7564) B = 23379840.
```

Therefore at every specialization `r != -4`, at least one minor is nonzero and
the cap-4 image has rank at least `9`. The two annihilators bound its rank by
`9`, which is exactly the dimension of the degree-at-most-10 target satisfying
the two independent conditions. Thus equality holds at the base cap.

For the arbitrary-degree step, pass from cap `d-1` to cap `d`, for `d>=5`.
The new `C_d` column has degree exactly `d+6` and leading coefficient `3`, while
all inherited columns have degree at most `d+5`. The constrained target space
also grows by exactly one dimension. Consequently the new column extends the
previous equality by one dimension, and filtered induction proves equality at
every cap `d>=4`. Every constrained polynomial has bounded degree, so choosing
a sufficiently large cap places it in one of these finite filtered targets.
Taking their union proves the unrestricted statement.

The direct annihilator calculations and leading-coefficient statements are the
written proof. The displayed minors and Bézout identity are exact finite
certificates. Additional basis-column checks in the verifier are regressions of
those arguments, not substitutes for them.

## 3. Conservative cap-3 anomaly

The following is a finite-filtration observation, not part of the unrestricted
exception classification. It records the genuine quadratic anomaly
`q=r^2+8r+6`; no claim about a complete cap-3 maximal-minor gcd is needed here.
At `q=0`, a nonzero source kernel, in the order
`(C_0,...,C_3,V_0,...,V_3)`, is

```text
(-(3r+16)/6, 2/3, 0, 0,
 -(11r+58)/2, (7r+86)/4, -(7r+34)/4, 1).
```

At such a root this kernel gives rank at most `7`, while the `7 x 7` minor
`972(r+4)^2` is nonzero because `q` and `r+4` are coprime. Hence the rank is
exactly `7`. (For comparison, one fixed `8 x 8` minor is nonzero at `r=0`, so
cap `3` has generic rank `8`.) Cap `4` recovers rank `9` at the quadratic roots
by the base Bézout certificate above. The roots of `q` exist individually only over a field
containing `sqrt(10)`; over `Q`, `q` has no root.

## 4. Proof when `r=-4`

Set

```text
D=E(E-1)(E+1).
```

Direct expansion, with the sign on the second filler block retained, gives the
exact quotient identities

```text
K_3(a_3,(E)_3 C)/D
 = (E-3)(E-2)(E-1)C(E)
   + E(E-2)(E+2)C(E+1)
   + (E+1)(E+2)(E+3)C(E+2),

-H_2(b_2,(E)_2 V)/D
 = (2-E)V(E)-(E+2)V(E+1).
```

In particular, the first formula includes the factor `(E-1)`. These formulas
prove `Im Phi subset D F[E]`.

For the reverse inclusion, consider the quotient operator

```text
B(V)=(2-E)V(E)-(E+2)V(E+1).
```

For every `n>=0`, `B(E^n)` has degree `n+1` and leading coefficient `-2`.
Therefore `B` triangularly reduces every positive-degree target to a constant.
The remaining constant is also attained by the full quotient map: an explicit
preimage of `1` is

```text
C=-1/15-(2/15)E,
V=-11/20-(11/10)E-(1/5)E^3.
```

Hence the quotient map is surjective onto `F[E]`, proving
`Im Phi=D F[E]`. This triangular argument is arbitrary-degree; finite scans in
the verifier are regression checks only.
