# Finite cokernel frontier for the W1/W2 exotic systems

**INDEPENDENTLY DERIVED — EXACT FINITE COMPUTATION — NOT PEER REVIEWED — BAND-SCOPED**

This note makes the two-filler cokernel question concrete for the W1/W2 systems
encoded in `verify_quantum_exotic_closure.py`, at the existing free-data caps
`d=1,2`. The exact calculations are over `Q` (and hence apply after extension to
`C`) in the fixed orientation `a_3!=0`, with gauge `b_3=0`, the displayed
normalized nonzero `b_2`, bottom normalization `b_-3=mu_3 a_-3`, and Weyl
membership imposed at every negative level. Here `d` bounds the raw polynomial
factors of the free coefficients after their
membership factors are removed. It does not prove an arbitrary-degree
obstruction, construct a full Weyl pair, or settle DC1/JC2. The unrestricted
admissible filler images still intersect infinitely dimensionally, exactly as proved in
`../dc1-program/two-filler-cross-cancellation.md`; nothing here weakens that
theorem or the isolated-block formulas of Lemma R.

Exact checks are in
[`verify_quantum_exotic_cokernel.py`](verify_quantum_exotic_cokernel.py).
The unrestricted image of `Phi` for the full degree-three AP family is classified
separately in [`quantum-ap-filler-image.md`](quantum-ap-filler-image.md); that
classification still requires testing the cascade-produced target `E-R`.

## 1. Convention and the actual target class

Retain the quantum convention

```text
Q_m = sum_(k+l=m) (b_l^[k] a_k - a_k^[l] b_l),
G = sum_(k=1)^3 sum_(j=0)^(k-1)
      (a_k^[j-k] b_-k^[j] - b_k^[j-k] a_-k^[j]),
Q_0=(T-1)G.
```

After solving `Q_4,Q_3,Q_2,Q_1=0`, split the central potential as

```text
G = R + K_3[c] - H_2[v],
c=b_-3=(E)_3 C,                 v=a_-2=(E)_2 V,
K_3[c]=sum_(j=0)^2 a_3^[j-3] c^[j],
H_2[v]=sum_(j=0)^1 b_2^[j-2] v^[j].
```

Thus the corrected map is

```text
Phi(C,V)=L_K(C)-L_H(V)=K_3[(E)_3 C]-H_2[(E)_2 V].
```

Here `K_3` uses the top `a_3`, while `H_2` uses `u=b_2`. The non-filler
residual is not chosen by setting reconstructed-operator kernel freedoms to
zero. On the positive-cascade solution variety it is

```text
R = G_1 + a_2^[-2]b_-2 + a_2^[-1]b_-2^[1],
```

with the solved `b_1,b_0,b_-1,b_-2`, the constrained `a_-1`, and every kernel
parameter emitted by the sequential `clean_solve` reconstruction retained. In
these four systems that reconstruction has one such parameter, the constant
`b_0` freedom. It contributes neither to `G` nor to `R`, so `R` is independent
of that reconstructed-operator kernel. It still varies with the free positive
data and, for W2 at `d=2`, is displayed branchwise. This limited independence
is checked rather than assumed.

Because `Q_0=1` is equivalent to `G=E`, the obstruction question is

```text
E-R in Im Phi,
```

not whether `R` itself lies in `Im Phi`. Equivalently, its class `[E-R]` must
vanish in the finite coefficient cokernel at the stated cap. It is not the class
of `R` alone that encodes the unit equation.

## 2. Finite coefficient matrices

The two fixed wall pairs are

```text
W1: a_3=E(E-2)(E-4),       b_2=(E-1)(E-4),
W2: a_3=E(E+2)(E+4),       b_2=E(E+3).
```

For cap `d`, use coefficient-ascending bases

```text
(C_0,...,C_d,V_0,...,V_d)
```

in the source and `(1,E,...,E^(d+6))` in the target. The exact matrices are:

```text
M(W1,d=1) =
[ 0    0    0    0]
[-192 -12    8  -10]
[414 -234  -30    5]
[-348 450   16  -24]
[159 -309   -2   15]
[-36  135    0   -2]
[ 3   -33    0    0]
[ 0     3    0    0]

M(W1,d=2) =
[ 0    0    0    0    0    0]
[-192 -12  -72    8  -10  -10]
[414 -234  -64  -30    5  -23]
[-348 450 -168   16  -24    8]
[159 -309  521   -2   15  -19]
[-36  135 -306    0   -2   14]
[ 3   -33  116    0    0   -2]
[ 0     3  -30    0    0    0]
[ 0     0    3    0    0    0]

M(W2,d=1) =
[ 0   0   0   0]
[ 0 -12   0   2]
[-18 -18  2   1]
[ 0 -18   0   0]
[15  15  -2  -1]
[ 0  27   0  -2]
[ 3   3   0   0]
[ 0   3   0   0]

M(W2,d=2) =
[ 0   0    0   0   0   0]
[ 0 -12  -24   0   2   2]
[-18 -18 -64   2   1   5]
[ 0 -18  -36   0   0   0]
[15  15   17  -2  -1  -3]
[ 0  27   54   0  -2  -2]
[ 3   3   44   0   0  -2]
[ 0   3    6   0   0   0]
[ 0   0    3   0   0   0].
```

Their ranks and finite cokernel dimensions are

| top | cap | matrix size | rank | coker dimension |
|---|---:|---:|---:|---:|
| W1 | 1 | `8 x 4` | 4 | 4 |
| W1 | 2 | `9 x 6` | 6 | 3 |
| W2 | 1 | `8 x 4` | 4 | 4 |
| W2 | 2 | `9 x 6` | 5 | 4 |

The coefficient-constant functional is always in the left nullspace because
both membership factors force every filler block to vanish at `E=0`. Complete
left-nullspace bases are generated and checked exactly by the verifier. It also
checks primitive rescalings of those bases; for readability, the residual
evaluations below use SymPy's rational basis normalization.

## 3. Parameter-generic rank and the rank-drop locus

For the AP family

```text
a_3=(E-r)(E-r-2)(E-r-4),
b_2=(E-r-1)(E-r-4),
```

the gcds of all maximal minors are, up to the displayed positive normalization,

```text
d=1: 18,
d=2: 36(r+4).
```

Therefore the finite `d=1` matrix has rank `4` for every specialization of `r`.
The finite `d=2` matrix has generic rank `6` and one rank-drop locus,
`r=-4`, where its rank is `5`. This explains why W2, and not W1, has the extra
finite cokernel functional at `d=2`. These statements concern `Phi` alone; the
positive-cascade solution variety also specializes and branches.

## 4. Solved residual and its class

At `d=1`, the solved positive cascade gives `R=0` for both W1 and W2, with the
constant `b_0` kernel retained. Hence the target class is `[E]`. Exact left-null
functionals evaluate nontrivially:

```text
W1 d=1: (0, 25887/99251, -11030/99251, 1172/99251),
W2 d=1: (0, 1, 0, 0),
```

in the verifier's rational nullspace bases. Thus `[E-R] != 0` in both finite
cokernels.

For W1 at `d=2`, solving the positive conditions without discarding free
parameters gives

```text
R = E a1_2^2 / 18 *
    (2E^3 a2_0 - 12E^3 a2_2 - 15E^2 a2_0 + 92E^2 a2_2
     +20E a2_0 - 140E a2_2 + 15a2_0 - 44a2_2).
```

Set `w=a1_2^2(a2_0-4a2_2)`. Two left-null evaluations of `E-R` are

```text
(77587-55275w)/5759959,
2(15982w-22329)/17279877.
```

They cannot both vanish, since `77587/55275 != 22329/15982`. Therefore
`[E-R] != 0`. This is an exact finite cokernel certificate. Its normalization
differs from the earlier `{8w=0,7w=9}` certificate because the present
functionals act directly on `G`, while the earlier elimination used coefficient
equations for `Q_0=(T-1)G`; both certify the same finite inconsistency.

At W2, `d=2`, the verifier's Gröbner certificate proves that the capped
positive-cascade ideal is exhausted by the two branches `a1_2=0` and `a2_0=0`.
The residual varies over that solution variety, and its displayed formula is
branchwise:

```text
a1_2=0: R=-(2/3) E^2 am1_2 (E-1)(E+1)(a2_0-2a2_2),
a2_0=0: R=(2/3) E(E-1)(E+1)(2E a2_2 am1_2+a0_2 a1_2).
```

On the `a1_2=0` branch one left-null functional evaluates `E-R` to `1`. On the
`a2_0=0` branch, two evaluations are `(a0_2 a1_2+9)/9` and
`-a0_2 a1_2/9`, whose sum is `1`. Hence neither branch has
`E-R in Im Phi`. The rank drop at `r=-4` is handled at the specialization; it
is not inferred from generic rank.

## 5. Scope

These are exact finite computations for the already encoded `d=1,2` systems.
They show concretely that the residual targets miss the finite filler images for
W1 and W2 at those caps, while retaining the sole kernel parameter emitted by
the capped sequential reconstruction (the constant `b_0` freedom). They do not
imply an arbitrary-degree obstruction. In unrestricted degree,
`Im L_K intersect Im L_H` remains infinite-dimensional, and the isolated
Lemma-R formulas remain valid only as isolated-block degree and leading-term
statements.
