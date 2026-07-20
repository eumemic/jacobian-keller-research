# Tame band-2 catalog and the strategic A\* question

**EXACT SYMPY — PROVED WITHIN STATED BOUNDS — NOT A FULL CLASSIFICATION**

Backing script: `verify_catalog.py` (Sections 4–6). All arithmetic exact.

## Strategic headline: `tameInAStar = FALSE`

**No tame band-2 pair lands in the classical resistant branch A\*.**

A\* is the branch that both upstream memos leave open: after the gauge that kills
`b_2` (`G→G−λF`, `λ=b_2/a_2`), a pair sits in A\* iff

```
a_2 = h² ≠ 0,   h constant,   κ ≠ 0,   gauged b_{−2} ≠ 0,   a_{−2} ≠ 0.
```

Every band-2 tame pair with `a_2≠0` instead has **gauged `b_{−2}=0`**, i.e. it lives
in branch **B0**. It is never in A0, never in A, hence never in A\*.

### Why (the mechanism — this is the load-bearing fact)

A band-2 tame automorphism with nonzero nonlinear part has the normal form

```
Φ = A₂ ∘ S ∘ A₁,   A₁,A₂ affine symplectic (SL₂),   S a single quadratic shear.
```

Write `A₁·(x,ξ)=(X,Ξ)` with `X=αx+βξ`, `Ξ=γx+δξ` (linear, weight ±1). Take
`S=P(μX²)` (the `Q(νΞ²)` case is identical by symmetry). Then the post-shear pair is
`(X, Ξ+μX²)`, and `A₂=[[a,b],[c,d]]` gives

```
F = aX + b(Ξ+μX²),        G = cX + d(Ξ+μX²).
```

The **only** source of weight `±2` (the extreme coefficients `a_{±2}, b_{±2}`) is the
single square term `μX²`, whose weight `+2` piece is `μα²x²` and weight `−2` piece is
`μβ²ξ²=μβ²τ²x^{−2}`. Therefore

```
a_2 = bμα²,   b_2 = dμα²   ⇒  λ = b_2/a_2 = d/b,
a_{−2} = bμβ²τ²,   b_{−2} = dμβ²τ² = (d/b)·a_{−2} = λ·a_{−2}.
```

The **same** ratio `λ=d/b` controls both extremes, because both come from the same
`μX²`. So the gauge `G→G−λF` that annihilates `b_2` **simultaneously** annihilates
`b_{−2}`:

```
gauged b_{−2} = b_{−2} − λ·a_{−2} = 0    (identically, for all parameters).
```

A\* asks for `b_{−2}` to **survive** the gauge that kills `b_2`. For a single-shear
tame pair that is structurally impossible: the two extremes share one origin.
`verify_catalog.py` §4 proves `gauged b_{−2} ≡ 0` symbolically for fully generic
`α,β,γ,δ,a,b,c,d,μ` (and `ν`).

Two corollaries fall out of the same normal form (`verify_catalog.py` §4, §6):

- **`h` is always constant.** `a_2=bμα²` has no `τ`-dependence, so a tame band-2
  square-sector pair is *automatically* in the constant-`h` sector. (Nonconstant `h`
  needs a `x³ξ`-type weight-2 term of total degree 4, which a quadratic shear cannot
  produce, and higher-degree shears leave band 2 — §6 confirms cubic shears never
  reach the square sector.)
- **`κ≠0` is forced.** Using `ad−bc=1`, the gauged `b_1` equals `−α/b` (P-shear) or
  `γ/a` (Q-shear); since `a_2≠0` forces `α≠0` (resp. `γ≠0`), gauged `b_1=κh≠0`. So the
  `κ=0` branch contains **no** tame square-sector pair — consistent with the memo's
  separate proof that `κ=0` admits no genuine polynomial pair with `a_2≠0`.

### Completeness of the normal form (why this covers all tame pairs)

Grade `ℂ[x,ξ]` by weight `w(x^mξ^n)=m−n`; multiplication adds weights and band-2 =
weights in `[−2,2]`. Affine maps keep linear forms in band 1. A degree-`d` shear
applied to a weight-`±1` argument produces weight up to `d`; to stay in band 2 the
shear must be degree ≤ 2, and its argument must be a genuine linear form. After one
quadratic shear the two coordinates are `X` (still weight `±1`) and `Ξ+μX²` (band 2);
any further degree-≥2 shear needs a weight-`±1` argument, but the only such direction
preserved is `X`, and shearing along it merely augments `μX²` — no new structure.
Hence every nonaffine band-2 tame pair reduces to `A₂∘S∘A₁`. This structural argument
is backed empirically in §6 by a bounded enumeration (below): length-≤3 words and all
two-quadratic-shear words produce only B0.

### Bounds under which this is proved

- **Rigorous for all parameters:** the normal-form identity `gauged b_{−2}≡0` (§4).
- **Enumerated:** words of length ≤ 3 over 10 affine SL₂ generators + 16 quadratic
  shears → 7874 distinct tame pairs, **1120** band-2 square-sector, **all B0, 0 in A\***.
  Of these, **186** have *both* raw `a_{−2}≠0` and raw `b_{−2}≠0` (they "look like" A\*
  before gauging) and in **all 186** the gauge collapses `b_{−2}` to 0. A separate
  two-quadratic-shear sweep gives 56 more band-2 square pairs, none in A\*. Cubic
  shears never reach the square sector.
- **Not claimed:** exhaustion over unbounded composition length as a *formal* proof;
  that rests on the weight-grading completeness argument plus the enumeration.

### Strategic consequence for closing A\*

In dimension 2 every polynomial automorphism is tame (Jung–van der Kulk). A band-2
Keller pair is a polynomial map with unit Jacobian; if it is an automorphism it is
tame, and we have just shown tame ⇒ B0. Therefore:

> **A\* (and A0) contain no automorphisms.** Any pair actually sitting in A\* would be a
> polynomial map with unit Jacobian that is *not* an automorphism — i.e. a Jacobian
> Conjecture counterexample, in band 2, constant-`h` square sector.

This tells the campaign how to spend effort on the resistant branch:

- **Do not** search for a tame witness in A\* — there is none. Closing A\* is *not* a
  "classify the tame members" task (the way B0 is); it is an **emptiness** task.
- Proving A\* empty (in the constant-`h` square sector) is **exactly equivalent** to
  ruling out a JC2 counterexample there. Conversely a nonempty A\* *is* a JC2
  counterexample. That is the real stake of the upstream "resistant branch."

## The catalog

Square data per entry (reduced gauge; `w` = gauged `b_{−2}`, `s = a_{−2}`,
`μ̃ = s/w` on branch A). All verified `{G,F}=1` and classified by `verify_catalog.py` §5.

| # | tame construction | `F` | `G` | supp `F`/`G` | `h` | `λ` | `κ` | `s=a_{−2}` | raw `b_{−2}` | `w` (gauged) | branch | in A\*? |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E1 | affine (band 1) | `2x−3ξ+5` | `x−ξ+7` | `{−1,0,1}` | — | — | — | — | — | — | not square (`a_2=0`) | no |
| E2 | `swap ∘ P(x²)` | `x²+ξ` | `−x` | `{−1,2}` / `{1}` | `1` | `0` | `−1` | `0` | `0` | `0` | **B0** (one-sided) | no |
| E3 | memo family (5.1), `h=1,e=0,c₁=1,κ=2` | `(x+ξ)²−ξ/2` | `2x+2ξ` | `{−2,−1,0,2}` / `{−1,1}` | `1` | `0` | `2` | `τ²` | `0` | `0` | **B0** | no |
| E4 | `A₂[[1,1],[1,2]] ∘ P(X²) ∘ A₁[[1,1],[0,1]]` | `x²+2xξ+x+ξ²+2ξ` | `2x²+4xξ+x+2ξ²+3ξ` | `{−2..2}` / `{−2..2}` | `1` | `2` | `−1` | `τ²` | `2τ²` | `0` | **B0** | no |

**E4 is the decisive entry.** It is a genuine two-sided band-2 Keller pair with *both*
negative extremes present before gauging (`a_{−2}=τ²`, `b_{−2}=2τ²`), the exact profile
that would be A\* — except `b_{−2}=λ·a_{−2}=2·τ²`, so the forced gauge `G→G−2F` sends
`b_{−2}→0`. It is B0. The 186 enumerated "both-extremes" pairs all behave this way.

E3 reproduces the upstream canonical family (5.1) (constant `h`, `e=0`): tame, B0,
`b_{−2}=0` by construction — the tame family the memos already knew lives in B0. E2 is
the minimal one-sided square-sector pair. E1 is the affine baseline (out of the square
sector, `a_2=0`).

Branch tally over the whole enumeration: **B0 : 1120,  A0 : 0,  A : 0,  A\* : 0.**

## Reproduce

```
uv run --with sympy python research/band2-square-sector/verify_catalog.py
```

Section 4 = symbolic A\*-exclusion proof; Section 5 = the four catalog entries above;
Section 6 = the bounded enumeration and the "gauge does real work" diagnostic
(186/186 collapse). The script ends `ALL CATALOG CHECKS PASSED`.
