# Two-filler cross-cancellation

**INDEPENDENTLY DERIVED — EXACT ALGEBRA — NOT PEER REVIEWED.**

This note corrects the unrestricted no-cross-cancellation claim for the two
nominal filler blocks. It does not construct a full Weyl pair, settle DC1 or
JC2, or show that these cancellations survive the remaining ladder equations.
Exact finite checks are in
[`verify_two_filler_cross_cancellation.py`](verify_two_filler_cross_cancellation.py).

Let `F` have characteristic zero, `T f(E)=f(E+1)`, and
`S_n=1+T+...+T^(n-1)`. For `k>=2`, set

```text
K_k[c]     = sum_(r=1)^k     a(E-r)c(E+k-r),
H_(k-1)[v] = sum_(r=1)^(k-1) u(E-r)v(E+k-1-r).
```

## Exact shift forms

Changing the summation index gives

```text
K_k[c]     = S_k((T^-k a)c),
H_(k-1)[v] = S_(k-1)((T^-(k-1)u)v).
```

If `f!=0`, then

```text
deg S_n(f)=deg f,              lc(S_n(f))=n lc(f).
```

Thus `S_n:F[E]->F[E]` is injective. On each finite degree filtration it is an
endomorphism with nonzero diagonal `n`, hence surjective; equivalently it is a
linear automorphism of `F[E]`.

The same sums telescope:

```text
(T-1)K_k[c]     = a(E)c(E+k)-a(E-k)c(E),
(T-1)H_(k-1)[v] = u(E)v(E+k-1)-u(E-k+1)v(E).
```

## Admissible images and their intersection

Weyl membership is

```text
c=(E)_k C,                 v=(E)_(k-1)V,
(E)_j=E(E-1)...(E-j+1).
```

Consequently, with

```text
q_K=(T^-k a)(E)_k,         q_H=(T^-(k-1)u)(E)_(k-1),
```

the maps on free polynomial factors have images

```text
Im L_K=S_k(q_K F[E]),      Im L_H=S_(k-1)(q_H F[E]).
```

For a nonzero polynomial `q`, the ideal `qF[E]` has codimension `deg q`.
Automorphisms preserve codimension, so for nonzero `a,u`,

```text
codim Im L_K=deg(a)+k,
codim Im L_H=deg(u)+k-1.
```

The intersection of two finite-codimensional subspaces has codimension at most
the sum of their codimensions. It is therefore finite-codimensional, and hence
infinite-dimensional, in `F[E]`. There are infinitely many nonzero
`P in Im L_K intersect Im L_H`. Injectivity of multiplication by `q_K,q_H` and
of the shift sums makes the corresponding admissible `c,v` unique.

No wall equation was used. The conclusion therefore applies to every wall pair
`(a,u)` for which both polynomials are nonzero. Membership-node equations reveal
some coefficients
triangularly, but only partially; they cannot rule out these common outputs.

## Degree resonance is necessary, not prohibitive

For nonzero `a,c,u,v`, the exact degree and leading-coefficient formulas are

```text
deg K_k[c]=deg(a)+deg(c),       lc K_k[c]=k lc(a)lc(c),
deg H_(k-1)[v]=deg(u)+deg(v),   lc H_(k-1)[v]=(k-1)lc(u)lc(v).
```

Hence every nonzero kernel pair for `K_k[c]-H_(k-1)[v]` obeys

```text
deg(a)+deg(c)=deg(u)+deg(v),
k lc(a)lc(c)=(k-1)lc(u)lc(v).
```

These resonance conditions are compatible with the infinite-dimensional image
intersection. They do not restore triangular annihilation.

## Concrete checks and scope

At `k=2`,

```text
a=E(E-1),                  u=E-1,
c=E(E-1)/2,                v=E(E-1)^2,
K_2[c]=H_1[v]=E(E-1)^2(E-2).
```

The verifier also checks an exact `k=3` two-block identity for the exotic wall
solution `a=E(E-2)(E-4)`, `u=(E-1)(E-4)`, with both memberships imposed. These
finite witnesses corroborate the arbitrary-degree linear-algebra proof above;
they do not prove it.

The result disproves any unrestricted claim that wall plus membership forbids
the two fillers from cross-cancelling. It does not produce a full Weyl pair,
settle DC1/JC2, or establish compatibility with the unsolved ladder equations.
The next object is

```text
coker Phi,       Phi(C,V)=L_K(C)-L_H(V),
```

together with the specific residual produced by the solved cascade blocks. This
cokernel is finite-dimensional; it may vanish depending on `a,u`, and no
nonvanishing claim is made. Nor is any generic-rank formula or assertion that
componentwise multiplication generates all syzygies used here.
