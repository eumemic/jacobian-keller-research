# The moment-unit program at arbitrary band `k`

**INDEPENDENTLY DERIVED — EXACT FINITE COMPUTATIONS — NOT PEER REVIEWED.**

This note separates one arbitrary-`k` theorem, Lemma P, from bounded computations
and unresolved reductions in the DC1 band program. It proves neither a general
moment-unit theorem nor a completed band-induction step, DC1, or JC2. The exact
finite checks are in [`verify_moment_unit_k.py`](verify_moment_unit_k.py).

Use `A_1[x^{-1}]=sum_i x^i C[E]`, `E=x∂`, with

```text
(x^a f(E))(x^b g(E)) = x^(a+b) f(E+b)g(E),
Q_m = sum_(i+l=m) (b_l(E+i)a_i(E)-a_i(E+l)b_l(E)).
```

Membership requires `E(E-1)...(E-r+1)` to divide every coefficient at level
`-r`. The orientation and gauge used below have `a_k != 0`, `b_k=0`, and bottom
normalization `b_{-k}=mu_k a_{-k}`. Zero-extreme and opposite-orientation cases
are separate reductions, not consequences of the calculations here.

## 1. Lemma P: the arbitrary-`k` theorem

The band-agnostic potential is

```text
G = sum_(i=1)^k sum_(r=1)^i
      (a_i(E-r)b_{-i}(E+i-r)-a_{-i}(E+i-r)b_i(E-r)),
Q_0=(T-1)G.
```

Membership gives `G(0)=0`. Hence `Q_0=1` would require `G=E`, and its slope is
`G(1)`.

> **Lemma P (moment-slope formula).** Under the displayed membership,
> orientation, gauge, nonzero-extreme, and normalization hypotheses,
>
> ```text
> G(1) = sum_(i=1)^k
>        (a_i(0)b_{-i}(i)-a_{-i}(i)b_i(0)).
> ```
>
> With `b_k=0` and `b_{-k}=mu_k a_{-k}`, this becomes
>
> ```text
> G(1) = sum_(i=1)^(k-1)
>        (a_i(0)b_{-i}(i)-a_{-i}(i)b_i(0))
>        + mu_k a_k(0)a_{-k}(k).
> ```

**Proof.** Evaluate the potential at `E=1`. In the first level-`i` sum,
`b_{-i}(1+i-r)` vanishes by membership for `r=2,...,i`; only `r=1` remains,
giving `a_i(0)b_{-i}(i)`. In the second sum, membership similarly kills every
term except `r=1`, giving `a_{-i}(i)b_i(0)`. Sum over `i`. The gauged formula
uses `b_k=0` and `b_{-k}=mu_k a_{-k}`. ∎

This is a written arbitrary-`k` proof. The verifier's runs at `k=2,3,4,5` are
corroboration, not the source of the theorem.

## 2. Wall and bounded root computations

Writing `u=b_{k-1}`, the top wall becomes

```text
u(E+k)a_k(E)-a_k(E+k-1)u(E)=0.
```

On one integer-root coset its necklace equation is

```text
S_k(S) delta(u) = S S_{k-1}(S) delta(a_k),
S_j=1+S+...+S^(j-1).
```

The factorization `S_k=product_(d|k,d>1) Phi_d` explains a prime/composite
difference in the divisibility gate. The non-effective universal cofactor
`1-S+S^2` supplies an exotic wall solution for every `k>=3`.

The verifier also performs normalized, distinct, single-coset integer-root scans
at degree `k` with `min(root)=0` and `max(root)<=5k`. Within those finite windows
the counts for `k=2,3,4,5` are `0,1,4,13`. These are bounded observations, not a
span bound, a finiteness theorem, or a classification of higher-degree,
repeated-root, noninteger-root, or multiple-coset tops. No `k=6` count is used.

## 3. The Lemma R gap

For the two nominal filler blocks, isolated leading-term calculations give

```text
K_k[c] = sum_(r=1)^k a_k(E-r)c(E+k-r),
lead(K_k[c]) = k lc(a_k)lc(c),
```

and

```text
-sum_(r=1)^(k-1) b_{k-1}(E-r)a_{-(k-1)}(E+k-1-r),
lead = -(k-1)lc(b_{k-1})lc(a_{-(k-1)}).
```

Thus neither filler self-cancels its own leading term in characteristic zero.
This does **not** prove the previously claimed uniform triangular annihilation:

- the two filler blocks can cross-cancel at a common degree; and
- either filler can collide with terms from already solved blocks.

Consequently the proposed L4 reduction is not uniformly proved. A proof must
control the full coefficient system, not merely each filler in isolation.

## 4. Scope of the finite moment systems

For a fixed top, fixed wall solution, and free degree `d`, the verifier
forward-solves the positive cascade and provisions every solved `b` with raw
degree cap

```text
2d + deg(a_k) + 2.
```

It then appends the coefficient equations for `Q_0`. The exact systems cover
selected instances at `k=3,4,5`, including a symbolic translation slice for one
`k=4` family at `d=1`. In those encoded systems, positive cascade plus `Q_0=1`
has unit ideal while the corresponding `Q_0=0` ideal is proper.

These are fixed-top, fixed-degree, forward-solvable polynomial ansätze. A proper
`Q_0=0` ideal only says that this partial system has a common zero; it does not
impose every negative ladder equation and does not construct a complete
homogeneous Weyl pair. The finite checks do not establish arbitrary free degree,
uniform translation outside the symbolic slice, or exceptional pivot and
denominator loci.

The Band-3 input is likewise limited: generic translation parameter `r` at
`d=1`, together with selected exact slices. It is not a completed Band-3 rung or
an arbitrary-degree theorem.

## 5. Residual problems

Lemma S—the low-end slope-forcing certificate—is one important residual only in
the fixed exotic, nonzero-wall, forward-solvable polynomial ansatz. It is not
the only missing ingredient and is not a single obstruction equivalent to DC1.
Independent residuals include:

1. a valid two-filler reduction controlling cross-cancellation and collision
   with solved blocks;
2. exceptional pivot, rank, and denominator loci in the forward solver;
3. arbitrary free degree;
4. higher-degree and multiple-coset top geometry;
5. shifted-power sectors beyond `k<=3`;
6. zero-extreme and opposite-orientation reductions; and
7. realization and induction reductions assembling all positive and negative
   ladder equations into the claimed band step.

Until these are resolved, there is no completed general-`k` rung, no completed
Band-3/4/5 rung supplied by this note, and no unrestricted counterexample
exclusion. If the characteristic-zero bridge is invoked globally, the safe
implication used here is only `JC2 => DC1`.

## 6. Reproduce

```sh
uv run --with sympy python research/dc1-program/verify_moment_unit_k.py
```
