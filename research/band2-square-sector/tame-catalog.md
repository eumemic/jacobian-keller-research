# Tame band-2 catalog and the classical A\* question

**EXACT SYMPY — GENERIC NORMAL-FORM IDENTITY + BOUNDED ENUMERATION — NOT A FULL CLASSIFICATION**

Backing script: `verify_catalog.py` (Sections 4–6). All arithmetic exact.

## Strategic headline for the checked class

**The generic single-quadratic-shear normal form, and every word in the stated
bounded enumeration, avoids the classical resistant branch A\*.**

A\* is the classical branch left open by the historical partial memo at commit
`91a053a`: after the gauge that kills `b_2` (`G→G−λF`, `λ=b_2/a_2`), a pair sits
in A\* iff

```
a_2 = h² ≠ 0,   h constant,   κ ≠ 0,   gauged b_{−2} ≠ 0,   a_{−2} ≠ 0.
```

Every pair in the displayed single-quadratic-shear normal form instead has
**gauged `b_{−2}=0`**, so it lives in branch **B0** rather than A\*. The bounded
word enumeration below found the same behavior in every sampled square-sector
case. No quantum statement is intended.

### Why (the mechanism — this is the load-bearing fact)

Consider the generic single-quadratic-shear normal form

```
Φ = A₂ ∘ S ∘ A₁,   A₁,A₂ affine symplectic (SL₂),   S a single quadratic shear.
```

The identity below is rigorous for all parameters in this class. Reduction of
an arbitrary unbounded tame word to this form is not claimed here.

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

A\* asks for `b_{−2}` to **survive** the gauge that kills `b_2`. In the displayed
affine–single-quadratic-shear–affine class that is structurally impossible: the
two extremes share one origin.
`verify_catalog.py` §4 proves `gauged b_{−2} ≡ 0` symbolically for fully generic
`α,β,γ,δ,a,b,c,d,μ` (and `ν`).

Two corollaries fall out of the same normal form (`verify_catalog.py` §4, §6):

- **`h` is constant in this normal form.** `a_2=bμα²` has no `τ`-dependence.
  The bounded sweep in §6 also found no cubic-shear word in the square sector;
  that observation is not an unbounded classification.
- **`κ≠0` in this normal form.** Using `ad−bc=1`, the gauged `b_1` equals
  `−α/b` (P-shear) or `γ/a` (Q-shear); since `a_2≠0` forces `α≠0` (resp.
  `γ≠0`), gauged `b_1=κh≠0`.

### Scope of the normal-form argument

Grade `ℂ[x,ξ]` by weight `w(x^mξ^n)=m−n`; multiplication adds weights and band 2
means weights in `[−2,2]`. This grading motivates the single-quadratic-shear
model, but it is not by itself a proof that every arbitrarily long tame word
reduces to that model. The exact claim here is therefore limited to the generic
normal form above and the bounded sweeps below.

### Bounds under which this is proved

- **Rigorous for all parameters in the displayed class:** the normal-form identity
  `gauged b_{−2}≡0` (§4).
- **Enumerated:** words of length ≤ 3 over 10 affine SL₂ generators + 16 quadratic
  shears produced 7,874 distinct tame pairs, including **1,120** band-2
  square-sector pairs, all B0 and none in A\*. Of these, **186** have both raw
  negative extremes nonzero, and in all 186 the gauge collapses `b_{−2}` to zero.
  A separate two-quadratic-shear sweep gives 56 more band-2 square pairs, none in
  A\*. The tested affine–cubic-shear–affine words never reach the square sector.
- **Not claimed:** exhaustion over unbounded composition length, exclusion of all
  classical tame automorphisms from A\*, or any quantum tame result.

### Strategic consequence

The exact mechanism shows why a witness search restricted to a single quadratic
shear cannot reach A\*. The bounded enumeration suggests the same obstruction
persists more broadly, but this memo does not promote that evidence to an
unbounded tame-word theorem. The independent written proof in
`classical-Astar.md` establishes that A\* itself is empty; this catalog is not
the completeness argument for that theorem.

## The catalog

Square data per entry (reduced gauge; `w` = gauged `b_{−2}`, `s = a_{−2}`,
`μ̃ = s/w` on branch A). All verified `{G,F}=1` and classified by `verify_catalog.py` §5.

| # | tame construction | `F` | `G` | supp `F`/`G` | `h` | `λ` | `κ` | `s=a_{−2}` | raw `b_{−2}` | `w` (gauged) | branch | in A\*? |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E1 | affine (band 1) | `2x−3ξ+5` | `x−ξ+7` | `{−1,0,1}` | — | — | — | — | — | — | not square (`a_2=0`) | no |
| E2 | `swap ∘ P(x²)` | `x²+ξ` | `−x` | `{−1,2}` / `{1}` | `1` | `0` | `−1` | `0` | `0` | `0` | **B0** (one-sided) | no |
| E3 | memo family (5.1), `h=1,e=0,c₁=1,κ=2` | `(x+ξ)²−ξ/2` | `2x+2ξ` | `{−2,−1,0,2}` / `{−1,1}` | `1` | `0` | `2` | `τ²` | `0` | `0` | **B0** | no |
| E4 | `A₂[[1,1],[1,2]] ∘ P(X²) ∘ A₁[[1,1],[0,1]]` | `x²+2xξ+x+ξ²+2ξ` | `2x²+4xξ+x+2ξ²+3ξ` | `{−2..2}` / `{−2..2}` | `1` | `2` | `−1` | `τ²` | `2τ²` | `0` | **B0** | no |

**E4 illustrates the gauge mechanism.** It is a genuine two-sided band-2 Keller
pair with both negative extremes present before gauging (`a_{−2}=τ²`,
`b_{−2}=2τ²`), the raw profile associated with A\*. But
`b_{−2}=λ·a_{−2}=2·τ²`, so the gauge `G→G−2F` sends `b_{−2}→0`; the pair is B0.
The 186 enumerated both-extremes examples behave the same way.

E3 reproduces the canonical family (5.1) from the partial memo snapshot at commit
`91a053a` (constant `h`, `e=0`): tame, B0, and `b_{−2}=0` by construction. E2 is
the minimal one-sided square-sector pair. E1 is the affine baseline (out of the
square sector, `a_2=0`).

Branch tally within the stated enumeration:
**B0 : 1120, A0 : 0, A : 0, A\* : 0.**

## Reproduce

```
uv run --with sympy python research/band2-square-sector/verify_catalog.py
```

Section 4 proves the gauge identity for the generic single-shear normal form;
Section 5 checks the four catalog entries above; Section 6 runs the bounded
enumeration and the 186/186 gauge-collapse diagnostic. The script ends
`ALL CATALOG CHECKS PASSED`.
