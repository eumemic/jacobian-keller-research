# Band-3 residual sectors, closed radical-correctly: diff-1 falls, diff-2 does not

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED — ONE NEW ARBITRARY-DEGREE CLOSURE + ONE REFUTATION + NAMED GAPS**

This memo settles the **restoration question** left explicitly unresolved by
[`broken-separation.md`](broken-separation.md) §3.3 — *does the cascade
`Q₄ ∧ Q₃` force the clean divisibility `h h^{[1]} | a₂` on the two
non-cube-separated shifted-cube classes?* — and does so **radical-correctly**,
with a certificate on one side and an explicit exact witness on the other. It
then reports what completed of the A\*-band3 `κ₂ ≠ 0` sector at higher caps.

The one-line summary:

> **The two broken classes have OPPOSITE answers, and the reason is visible in a
> single scalar.** Evaluating the rungs at the gcd node `r` gives, for *arbitrary*
> polynomial data:
>
> ```
> diff-1:   Q₃(r) = (κ/2) · a₂(r-1)²          ← a PERFECT SQUARE
> diff-2:   Q₃(r-1) = a₂(r-1) · b₁(r+1)       ← a PRODUCT of two free quantities
> ```
>
> A square puts `a₂(r-1)` in the **radical** of the cascade ideal (exponent 2) —
> so `a₂(r-1)` genuinely vanishes on the variety, while its Gröbner **normal form
> is nonzero**. That nonzero normal form is exactly what the previous pass
> mistook for non-forcing. A product does not: the variety splits, and the branch
> `b₁(r+1) = 0` carries an explicit exact solution with `a₂(r-1) ≠ 0`.
>
> **diff-1 — RESTORATION IS FORCED, and the class CLOSES AT ARBITRARY DEGREE.**
> `Q₄ ∧ Q₃ ⇒ a₂(r-1) = 0`, and then `Q₄` alone delivers the rest of the jet
> (`a₂(r) = a₂(r+1) = a₂'(r) = 0`, `b₁(r) = b₁(r+1) = 0`), i.e. the **full clean
> shapes** `h h^{[1]} | a₂` **and** `h | b₁`. From there `gcd(h,h^{[2]}) = 1`
> gives `h | a₁`, the potential factorises as `G = h^{[-1]}M`, and
> `Q₀ = 1 ⇒ G = E ⇒ h^{[-1]} | E` is impossible for `deg h^{[-1]} = 2`.
> **The diff-1 sector is EMPTY at arbitrary degree with no `κ`/`b₁` side
> condition.** This memo proves the two nonzero branches; Nonpositive-D Exclusion Theorem of
> [`shifted-cube-completion.md`](shifted-cube-completion.md) supplies the remaining
> `κ=0,b₁=0` branch.
>
> **diff-2 — RESTORATION IS REFUTED.** Explicit exact solution of `Q₄ = Q₃ = 0`
> at **symbolic `(r, κ)`**:
> `a₂ = 3(E-r)(E-r-1)`, `b₁ = 2κ(E-r-1)²`, `a₁ = -3`, `b₀ = 0`, with
> `a₂(r-1) = a₂(r+2) = 6 ≠ 0`. The cascade forces **exactly** the degraded
> `(E-r)(E-r-1) | a₂` and nothing more; `Q₃` contributes only the branch
> disjunction `a₂(r-1) · b₁(r+1) = 0`. **The clean-divisibility route to diff-2 is
> dead — a different mechanism is required.** diff-2 remains open at arbitrary
> degree, with bounded emptiness at cap `d ≤ 3` (committed) and `d = 4` (HEAVY).

Exact certificate: [`verify_band3_sectors.py`](verify_band3_sectors.py) — the current
no-`msolve` run recorded **145 passed, 0 failed, 3 skipped in 87.0 s**; an independent
`python -O` rerun recorded the same status in **90.1 s**. The skipped groups were the `HEAVY` cap-`d=4` leg, all
load-bearing `msolve`/live-parser legs, and the Task C `msolve` probe. Historical
solver-backed results below were not rerun in this environment. Every load-bearing upstream
fact (the crossed-product ladder engine `Q_m = [D,X]_m`, `Q₀ = (T−1)G`, `G(0)=0` under
membership, the `Q₅` wall, the
general `Q₄`/`Q₃` identities, `(T^n−1)` surjectivity) is **re-derived in file**,
not cited.

Conventions frozen from the corpus: `A₁[x^{-1}] = ⊕_k x^k C[E]`, `E = x∂`,
`(x^a f)(x^b g) = x^{a+b} f(E+b) g(E)`, `f^{[n]}(E) = f(E+n)`, `T f = f^{[1]}`,
`Q_m = ∑_{k+l=m}[b_l^{[k]}a_k − a_k^{[l]}b_l]`, `[D,X]=1 ⇔ Q_m=δ_{m0}`,
membership `(E)_j | a_{-j},b_{-j}`, gauge `b₃=0`, `Q₀=(T−1)G`, `G(0)=0`.
Sector: `a₃ = h h^{[1]}h^{[2]}`, `b₃ = 0`, `b₂ = κ h h^{[1]}` (the `Q₅` wall),
`diff-1: h = (E-r)(E-r-1)`, `diff-2: h = (E-r)(E-r-2)`.

---

## 0. Method: what "forced" and "not forced" are allowed to mean

House rule 3 of this program says a nonzero Gröbner **normal form** proves
non-membership in the **ideal**, never in the **radical**, and therefore never
establishes geometric non-forcing. `broken-separation.md` §3.3 recorded exactly
that caveat and left the question open.

**This memo's finding is that the caveat was not academic — it was the whole
story for diff-1.** The cascade ideal contains `a₂(r-1)²` but not `a₂(r-1)`. A
normal-form probe therefore *must* return nonzero, and the correct conclusion is
the opposite of the naive reading.

Accordingly, every verdict below is one of:

| verdict | admissible evidence | used for |
|---|---|---|
| **forced** | exponent-`k` membership identity (explicit cofactors) **and** an independent Rabinowitsch unit-ideal certificate | diff-1 |
| **not forced** | an **explicit exact solution** with the value nonzero | diff-2 |
| **empty** | unit ideal (weak Nullstellensatz), with the cap stated | bounded legs |
| **nonempty** | proper ideal (weak Nullstellensatz) or an explicit point | every control |

**Degree-freeness device.** §§3–5 carry no degree cap of any kind: `a₂, b₁, a₁,
b₀` are entered into the verifier as **undetermined functions** (`sympy.Function`),
and every node identity is derived by evaluating explicit `h`-products at
`E = r + j`. Since only shifts and multiplication by explicit polynomials occur,
the identities hold for polynomial data of *arbitrary* degree — indeed for
arbitrary functions. The in-file machine scope of §§3–5 is therefore
**degree-free**, not "verified up to cap `d`".

---

## 1. The engine and the rungs (re-derived, no shape assumption)

Verifier `§0` re-derives `Q_m = [D,X]_m` for `m ∈ [−6,6]` against the direct
crossed-product commutator on generic symbolic band-3 data, `Q₀ = (T−1)G` with
the staggered potential, `G(0) = 0` under membership, and the surjectivity of
`(T^n − 1) : C[E] → C[E]` with kernel the constants (used at every telescoping
step). Each of these carries a **negative control** that fails when the
statement is corrupted (a wrong shift, dropped membership, a wrong wall).

Verifier `§1` re-derives the `Q₅` wall `b₂ = κ h h^{[1]}` for *any* `h`, and the
two rung identities used throughout, with **no shape assumption**:

```
Q₄ = h^{[1]}h^{[2]}( b₁^{[3]} h − h^{[3]} b₁ ) + κ( h^{[2]}h^{[3]} a₂ − h h^{[1]} a₂^{[2]} ),
Q₃ = a₃( b₀^{[3]} − b₀ ) + ( b₁^{[2]} a₂ − a₂^{[1]} b₁ ) + κ h^{[1]}( h^{[2]} a₁ − h a₁^{[2]} ).
```

`§1` also re-derives the cancellation the broken classes are accused of losing:
with `a₂ = h h^{[1]}g`, `b₁ = hβ`, the middle term is `a₃(β^{[2]}g − g^{[1]}β)`.

**Normalisation (verifier `§2`, degree-free).** `E ↦ E + ρ` is a ring
automorphism commuting with `T`, and `Q₄(sκ; a₂, s b₁) = s·Q₄(κ; a₂, b₁)`. So the
`(Q₄,Q₃)` analysis at `(r, κ)` and at `(0, 1)` are in exact bijection. The memo
nonetheless states every result at **symbolic `r` and `κ`**, because §5 needs
membership, which is *not* shift-invariant.

---

## 2. The degree-free node table

Verifier `§3` derives, with `a₂, b₁, a₁, b₀` undetermined:

**diff-1** (`h(r+j)`, `j = −3..2`, equals `[12, 6, 2, 0, 0, 2]`):

```
Q₄(r)     = 12κ a₂(r)
Q₄'(r-1)  = 2κ ( a₂(r-1) + a₂(r+1) )
Q₄(r-3)   = 144 b₁(r) − 72κ a₂(r-1)
Q₄(r+1)   = 72κ a₂(r+1) − 144 b₁(r+1)
Q₄'(r)    = 28κ a₂(r) + 12κ a₂'(r) − 12 b₁(r)
Q₃(r)     = a₂(r) b₁(r+2) − a₂(r+1) b₁(r)
Q₃(r-1)   = a₂(r-1) b₁(r+1) − a₂(r) b₁(r-1)
```

**diff-2** (`h(r+j)`, `j = −3..2`, equals `[15, 8, 3, 0, −1, 0]`):

```
Q₄(r-2)   = −24κ a₂(r)          Q₄(r+1)  = 24κ a₂(r+1)
Q₄(r-3)   = −120κ a₂(r-1) + 360 b₁(r)
Q₄(r+2)   =  120κ a₂(r+2) − 360 b₁(r+2)
Q₄'(r-1)  = −2κ a₂(r-1) + 6κ a₂(r+1) + 6 b₁(r+2)
Q₃(r)     = a₂(r) b₁(r+2) − a₂(r+1) b₁(r)
Q₃(r±1)   = ∓ a₂(r-1) b₁(r+1)   (on the Q₄ locus)
```

Note the `Q₃` rows: `a₃` and `h^{[1]}` both vanish at these nodes, so `Q₃`
there is **purely the bilinear middle term** — no `a₁`, no `b₀`. That is why the
restoration question is decided by two scalars and not by a Gröbner basis.

---

## 3. diff-1 — restoration is FORCED (radical certificate, exponent 2)

On `Q₄ = 0`: `a₂(r) = 0`, `a₂(r+1) = −a₂(r-1)`, `b₁(r) = κ a₂(r-1)/2`. Substituting
into `Q₃(r)`:

> **`Q₃(r) = (κ/2) · a₂(r-1)²`.**

The verifier commits this as an **exact polynomial identity with explicit
integer cofactors** (`§4`, degree-free):

```
144 κ² a₂(r-1)²  =  288κ·Q₃(r) − 24 b₁(r+2)·Q₄(r) + Q₄'(r-1)·Q₄(r-3)
                    + 72κ a₂(r-1)·Q₄'(r-1) − 2κ a₂(r-1)·Q₄(r-3),
```

together with an **independent second certificate** from the node `r-1`:

```
144 κ² a₂(r-1)²  =  −288κ·Q₃(r-1) − 24 b₁(r-1)·Q₄(r)
                    + 72κ a₂(r-1)·Q₄'(r-1) − 2κ a₂(r-1)·Q₄(r+1).
```

Both are checked with perturbation controls that must (and do) fail. Since `κ ≠ 0`
and `char = 0`, `144κ²` is a unit, so `a₂(r-1)² ∈ I` and hence
`a₂(r-1) ∈ √I`. The verifier independently confirms this with **Rabinowitsch**
(`I + (1 − t·a₂(r-1))` is the unit ideal) and pairs it with the control that
`b₁(r-1)` and `b₁(r+2)`, which are genuinely free, are **not** forced.

**And it records the trap explicitly:** in the same node ideal, `a₂(r-1)` has
**nonzero** Gröbner normal form while `a₂(r-1)²` reduces to zero. A normal-form
probe alone would have reported "not forced" — the exact error corrected here.

### 3.1 The full clean shape is restored, not just one node

`Q₄` alone gives `a₂'(r) = a₂(r-1)/2` (so `a₂'(r)` is indeed free at the `Q₄`
level — `broken-separation.md` §2.1 was right about that) and
`b₁(r+1) = κ a₂(r+1)/2`. Once `Q₃` forces `a₂(r-1) = 0`, the whole jet collapses:

```
a₂(r-1) = a₂(r) = a₂'(r) = a₂(r+1) = 0    ⇒   h h^{[1]} | a₂   (clean),
b₁(r)   = b₁(r+1) = 0                     ⇒   h | b₁          (clean).
```

The jet-to-divisibility step is itself degree-free: `P mod h h^{[1]}` has degree
`≤ 3` and the four functionals `(P(r-1), P(r), P'(r), P(r+1))` form an invertible
confluent Vandermonde on that remainder space (verifier checks `det ≠ 0`); same
for `h` with `(P(r), P(r+1))`. Every one of the six functionals is separately
certified by Rabinowitsch, with the non-vacuity control that `b₁(r-1)`, `b₁(r+2)`
are not.

> **Theorem (diff-1 restoration, arbitrary degree, symbolic `r`, `κ ≠ 0`).**
> `Q₄ = Q₃ = 0` ⇒ `h h^{[1]} | a₂` **and** `h | b₁`. The clean shapes are fully
> restored; the `Q₄`-level degradation is repaired by `Q₃` alone.

---

## 4. diff-1 — the class CLOSES at arbitrary degree

With the clean shapes in hand the 2-separation machinery runs, and diff-1 has
the coprimality it needs: `gcd(h, h^{[2]}) = gcd(h, h^{[3]}) = 1` (only
`gcd(h,h^{[1]})` is broken, and no step below uses it). Verifier `§5`:

1. **`h | a₁`.** With the restored shapes the middle term of `Q₃` is `a₃`-divisible,
   so `Q₃ = 0` forces `a₃ | R` with `R = κ h^{[1]}(h^{[2]}a₁ − h a₁^{[2]})`. At the
   double node `r`, `R(r) = 0` automatically and `R'(r) = 2κ a₁(r)`; at the simple
   node `r+1`, `R(r+1) = 12κ a₁(r+1)`. Hence `a₁(r) = a₁(r+1) = 0`, i.e. `h | a₁`.
   *(Degree-free; the `ev'` derivative-node equation of `broken-separation.md` §4.2
   is exactly what does the work at the double node.)*
2. **`G = h^{[-1]} M`.** Pure divisibility bookkeeping, and degree-free: every term
   of `G` carries a factor `F_k^{[j-k]}` with `F₁ = h`, `F₂ = h h^{[1]}`,
   `F₃ = a₃`, `1 ≤ k ≤ 3`, `0 ≤ j < k`, and `h^{[-1]}` divides each of the six.
   The verifier checks all six, plus the control that `h^{[-1]} ∤ h` (so the
   bookkeeping is not vacuous), plus an instance cross-check on generic quotients.
3. **The affine kill.** `Q₀ = 1 ⇒ (T−1)G = 1`, and `G(0) = 0` under membership, so
   `G = E`. Then `h^{[-1]} | E`. But `deg h^{[-1]} = 2 > 1`. **Contradiction.**

> **Theorem (diff-1 closure, ARBITRARY DEGREE).** There is no genuine pair
> `[D,X] = 1` with `a₃ = h h^{[1]}h^{[2]}`, `h = (E-r)(E-r-1)`, `b₃ = 0`,
> `b₂ = κ h h^{[1]}` with `κ ≠ 0`, membership-valid tail and `Q₀ = 1`.
> **The diff-1 broken class is closed at arbitrary degree with no side condition.**
> Sections 3–4 close `κ ≠ 0`, §4.1 closes `κ=0,b₁≠0`, and Nonpositive-D Exclusion Theorem of
> [`shifted-cube-completion.md`](shifted-cube-completion.md) closes the remaining
> `κ=0,b₁=0` branch independently of `h`.

### 4.1 The `κ = 0` branch, and one inherited gap

`κ = 0` (`b₂ = 0`): `Q₄ = h^{[1]}h^{[2]}·W` with `W = h b₁^{[3]} − h^{[3]}b₁`, and
`W(r) = −6b₁(r)`, `W(r+1) = −12b₁(r+1)`, so `h | b₁`; writing `b₁ = hψ` gives
`ψ^{[3]} = ψ`, hence `ψ = c` constant. For `c ≠ 0`, `Q₃ = 0` forces `a₃ | F` with
`F = h^{[2]}a₂ − h a₂^{[1]}`, and `F(r) = 2a₂(r)`, `F(r+1) = 6a₂(r+1)`,
`F(r-2) = −6a₂(r-1)`, `F'(r) = 3a₂(r) + a₂(r+1) + 2a₂'(r)` — the same clean jet.
`Q₂ = 0` then forces `h h^{[1]} | V = h^{[1]}a₁ − h a₁^{[1]}`, and
`V(r+1) = 2a₁(r+1)`, `V(r-1) = −2a₁(r)`, so `h | a₁`. The closure of §4 runs.
**All degree-free (verifier `§5`).**

> **Historical residual, now closed.** This memo correctly identified that when
> `κ=0` and `b₁=0`, `Q₃=a₃(T³−1)b₀` only forces `b₀` constant and does not constrain
> `a₂`; it also correctly located the same silent `c≠0` assumption in the old
> cube-separated chain. Nonpositive-D Exclusion Theorem of
> [`shifted-cube-completion.md`](shifted-cube-completion.md) bypasses that rung:
> `κ=0,b₁=0` gives `band D≤0`, which is impossible for a membership-valid Weyl pair
> with `band X=3`. Thus this inherited, `h`-independent hole is closed at arbitrary
> degree in every shifted-cube class.

---

### 4.2 The symbolic `Q₂` congruence on `κ=0`, `b₁=ch`, `c≠0`

After `b₁=ch` and the restored shape `a₂=h h^{[1]}g`, direct construction of the
crossed-product rung gives

```text
V = h^{[1]}a₁ - h a₁^{[1]},
Q₂ - cV ∈ (h h^{[1]}).
```

The verifier constructs `Q₂` symbolically and divides the difference exactly; this is not a
node sample. It records the congruence needed by the `κ=0,b₁≠0` closure. The complementary
`c=0`, equivalently `b₁=0`, branch is closed independently by Nonpositive-D Exclusion Theorem of
[`shifted-cube-completion.md`](shifted-cube-completion.md).

---

## 5. diff-2 — restoration is REFUTED (explicit exact `(Q₄,Q₃)`-locus witness)

On the `Q₄` locus, `a₂(r) = a₂(r+1) = 0`, `b₁(r) = b₁(r+2) = κ a₂(r-1)/3`, and
`a₂(r+2) = a₂(r-1)` (the tie). Then `Q₃(r) ≡ 0` identically, and the *only*
`Q₃` node conditions free of `a₁, b₀` are

```
Q₃(r-1) = a₂(r-1)·b₁(r+1),      Q₃(r+1) = −a₂(r-1)·b₁(r+1).
```

A **product**, not a power. Rabinowitsch confirms `a₂(r-1)` is **not** in the
radical, while the *same* machinery applied to diff-1 returns forced — so the
two classes are separated by the mathematics, not by the setup.

### 5.1 The witness

> **Explicit exact solution at symbolic `(r, κ)`** (verifier `§6`):
> ```
> h  = (E-r)(E-r-2),    a₃ = h h^{[1]}h^{[2]},   b₃ = 0,   b₂ = κ h h^{[1]},
> a₂ = 3(E-r)(E-r-1),   b₁ = 2κ(E-r-1)²,        a₁ = -3,   b₀ = 0.
> ```
> `Q₄ = 0` and `Q₃ = 0` **exactly**, for every `r` and every `κ`. All coefficients
> are integers. And
> ```
> a₂(r-1) = a₂(r+2) = 6 ≠ 0,      (E-r)(E-r-1) | a₂  ✓,      h h^{[1]} ∤ a₂  ✗.
> ```
> The scaled family `(λ a₂, λ b₁, λ² a₁, 0)` solves `Q₄ = Q₃ = 0` for every `λ`, so
> this is a one-parameter family, not an isolated accident.

> **Corollary (refutation, with exact scope).** For diff-2, the displayed
> `(Q₄,Q₃)`-locus family proves that the cascade does **not** restore the clean
> divisibility `h h^{[1]} | a₂`: every point in the family has `Q₄=Q₃=0`, while
> `h h^{[1]} ∤ a₂`. The witness does not prove that the degraded
> `(E-r)(E-r-1) | a₂` is the strongest universally forced divisor, nor does it
> exclude every possible intermediate factor. It is **not** a genuine polynomial-tail
> or Weyl-pair point; in general its `Q₂` is nonzero. **Only the clean-divisibility
> restoration route is refuted; a different mechanism is required.**

### 5.2 The surviving family, and why it is a genuine wall

The `(Q₄,Q₃)` variety splits into the branch `a₂(r-1) = 0` (which recovers
`(E-r)(E-r-1)(E-r+1)(E-r-2) | a₂`, i.e. the clean divisor, by the tie) and the
branch `b₁(r+1) = 0`, on which `a₂(r-1) = a₂(r+2)` is **free**. On the second
branch the potential factorisation `G = h^{[-1]}M` is **unavailable**: the
verifier checks that for the witness `h^{[-1]} ∤ a₂^{[-1]}`. That is precisely
the step that closes diff-1 in §4, and it is precisely the step diff-2 loses.

Does the survivor extend to a genuine pair? **Not at bounded degree.** Verifier
`§7`, exact SymPy over `ℚ` at `r = 0`, `κ = 1`, full cascade `Q₄=Q₃=Q₂=Q₁=0`,
`Q₀ = 1`, gauge, wall and genuine membership:

| cap `d` | diff-1 | diff-2 | tier |
|---|---|---|---|
| `d = 1` | unit ideal ⇒ **EMPTY** | unit ideal ⇒ **EMPTY** | committed (0.1 s) |
| `d = 2` | unit ideal ⇒ **EMPTY** | unit ideal ⇒ **EMPTY** | committed (0.5 s) |
| `d = 3` | unit ideal ⇒ **EMPTY** | unit ideal ⇒ **EMPTY** | committed (4 s) |
| `d = 4` | unit ideal ⇒ **EMPTY** | unit ideal ⇒ **EMPTY** | `HEAVY=1` (~90 s) |

**Non-vacuity, now in-file** (the previous memo had to import this from an
external audit): dropping `Q₀ = 1` leaves a **proper** ideal at `d = 2` for both
classes, so by the weak Nullstellensatz the sector is **nonempty** without the
moment unit. The emptiness certificates are therefore meaningful.

For diff-1 the bounded table is now an *independent cross-check* of the
arbitrary-degree theorem of §4; for diff-2 it is the only evidence there is.

---

## 6. TASK B — A\*-band3, the constant-top `κ₂ ≠ 0` sector

Verifier `§8` re-derives the sector structure exactly and fast: the constant-top
wall `Q₅ = (T³−1)b₂ ⇒ b₂ = κ₂`; the rung
`Q₄ = (b₁^{[3]} − b₁) + κ₂(a₂ − a₂^{[2]})`; the **vacuity** of the potential
factorisation (`h = 1 ⇒ h^{[-1]} = 1`, `G = M`); the bottom proportionality
`Q₋₆ ⇒ b₋₃ = μ₃a₋₃`; the inhomogeneous `μ₃`-source in `Q₋₅`; and the Lemma-P
moment slope `G(1) = ∑_i(a_i(0)b_{-i}(i) − a_{-i}(i)b_i(0))` with `G(0) = 0`.

There is also an exact scaling normalization over `C̄`. Under

```text
x ↦ sx,   ∂ ↦ s⁻¹∂,   X ↦ s⁻³X,   D ↦ s³D,
```

one has `a_k ↦ s^{k-3}a_k`, `b_k ↦ s^{k+3}b_k`, each rung `Q_m ↦ s^mQ_m`, and
`κ₂ ↦ s⁵κ₂`. Thus `[D,X]=1`, `a₃=1`, the E-degree caps, and negative-band
falling-factorial membership are preserved. For `κ₂≠0`, choose `s⁵=κ₂⁻¹` to normalize
`κ₂` to `1` over `C̄`. Such a fifth root need not lie in `Q`, so this is not generally a
rational normalization.

**Controls, committed:**
- **`κ₂ = 0` slice nonempty, explicit point.** `U = x + 2∂`, `X = U³ − ∂`,
  `D = U`. The verifier *computes* `[D,X]` in the crossed product and checks it
  equals `1` with `a₃ = 1`, `b₃ = 0`, `b₂ = 0` and genuine membership.
- **Strict complete-record msolve parsing.** A result is accepted only when one complete
  balanced bracket record consumes the output (apart from the documented suffix); no prefix
  acceptance is allowed. Deterministic regressions reject `[-1]garbage`,
  `[-1]\nSECOND_RECORD`, `[totally malformed`, and `[]garbage`. When msolve is available,
  live known-unit and known-feasible ideals validate the parser before load-bearing calls.
  Integer-body guards still reject rational literals and Python `**` syntax.

**Engine finding (recorded for the corpus, and it corrects the previous
diagnosis).** SymPy's Gröbner engine does **not** finish the A\*-band3 cap-`d = 1`
system within 20 minutes, although the *rigid* diff-1/diff-2 sectors at `d = 3`
(44 variables) finish in ~4 s; the constant top `a₃ = 1` removes all rigidity.
**msolve does the same cap-`d = 1` system in 0.1 s.** So the binding constraint on
the previous attempt was **engine selection**, not (as previously supposed)
rational-coefficient misparsing — trap #2 is real and is guarded here, but it was
not what made `d = 3` unreachable. Every A\*-band3 emptiness leg in this verifier
is msolve-only.

**Historical msolve run results (not reproduced by the current no-msolve run):**

| leg (cap `d`, `κ₂ = 1`, integer system) | engine / tier | result |
|---|---|---|
| `d = 1` emptiness over `ℚ` | historical characteristic-zero msolve run (0.1 s), not reproduced currently | `[-1]` = **EMPTY over `C̄`** |
| `d = 1` `κ₂ = 0` slice | historical characteristic-zero msolve run (0.3 s), not reproduced currently | `[1,22,-1,[]]` = **NONEMPTY** (positive-dimensional) |
| `d = 1` unit ideal mod `p = 2³⁰+3` | historical msolve `-g` run (0.1 s), not reproduced currently | `[1]` — **corroboration only** |
| `d = 1` positive cascade alone | explicit point, instant | `a₂ = 5`, `b₁ = 7` const ⇒ **NONEMPTY** with `κ₂ = 1` |
| `d = 1` drop `Q₀ = 1` | msolve char-0, `HEAVY=1`, 1200 s cap | slow (positive-dimensional); redundant with the two controls above |
| `d = 2` (branch `a₋₃ = 0` and full) | msolve char-0, `HEAVY=1`, 3000 s cap | **not completed in this session** (>150 s probe; the corpus records ≈35 min) |
| **`d = 3`** (branch `a₋₃ = 0` and full) | msolve char-0, `HEAVY=1`, 3000 s cap | **NOT ACHIEVED** — the leg is implemented and gated; timeout is an optional SKIP normally and a failure with `--require-msolve` |

> **Honest statement of Task B.** The `d = 3` target is **not** met. Historical
> characteristic-zero msolve runs established the cap-`d=1` statements above; finite-field
> output is corroboration only and is not promoted to a characteristic-zero conclusion.
> The newer verifier conditionally reruns every solver leg when msolve is available, reports
> the resolved executable path and version/identity where available, and supports
> `--require-msolve`. Without that flag, unavailable or unsuccessful invoked solver legs are
> explicit optional SKIPs. With it, every solver leg actually invoked fails closed on timeout,
> nonzero exit, missing output, or malformed/unrecognized output; HEAVY-only legs are invoked,
> and therefore required, only when `HEAVY=1`. In the current environment msolve was absent,
> so parser live-validation and all msolve-dependent mathematical legs were explicitly
> `NOT RERUN/SKIPPED`; only the local strict-parser regressions and non-msolve controls ran.

The `F_p` row is **corroboration only**. Unit mod `p` does *not* imply unit over
`ℚ` — `(p x − 1)` is a counterexample — and the verifier says so at the check.

## 7. TASK C — formal constant-top A*-band3, `κ₂=1`, cap `d=1`

In this precisely bounded formal sector, ask whether the Lemma-P moment slope `G(1)` is
forced to `0` by the tail alone (all `Q_m=0` for `m≠0`, membership, and no `Q₀` condition).
A historical characteristic-zero msolve Rabinowitsch run returned the nonempty record
`[1,23,-1,[]]`, so within this cap `G(1)` is not forced to zero. This is a bounded negative,
not a degree-free conclusion and not a claim about the entire A*-band3 sector.

The current no-msolve validation did **not** rerun Task C and reported it as
`NOT RERUN/SKIPPED`. Use `--require-msolve` in an environment where the solver leg is required.

## 8. What this changes in the band-3 ledger

| class | status before | **status now** |
|---|---|---|
| cube-separated / 2-separated `h` (incl. diff-3, `(E-r)²`) | closed, arbitrary degree, with an unnoticed `c=0` hole | **closed, arbitrary degree; hole repaired by [`shifted-cube-completion.md`](shifted-cube-completion.md), Nonpositive-D Exclusion Theorem** |
| **diff-1** `(E-r)(E-r-1)` | degraded `Q₄` forcing; restoration **unresolved**; bounded emptiness only | **CLOSED at arbitrary degree, all three branches**: `κ ≠ 0` by the radical-correct certificate of §4; `κ = 0, b₁ ≠ 0` by direct cascade (§4.1); `κ = 0, b₁ = 0` **by Nonpositive-D Exclusion Theorem** ([`shifted-cube-completion.md`](shifted-cube-completion.md)). The direct residual/congruence route remains silent at `κ = 0, b₁ = 0`, but Nonpositive-D Exclusion Theorem closes it independently. |
| **diff-2** `(E-r)(E-r-2)` | degraded `Q₄` forcing; restoration **unresolved** | `κ = 0, b₁ = 0` **CLOSED by Nonpositive-D Exclusion Theorem**; `κ = 0, b₁ ≠ 0` closed by direct cascade; surviving **`κ ≠ 0` branch open** at arbitrary degree (clean restoration refuted; bounded-empty `d ≤ 4`) |
| `κ = 0 ∧ b₁ = 0` sub-branch | (not previously isolated) | **CLOSED at arbitrary degree by Nonpositive-D Exclusion Theorem**, inherited and `h`-independent |

**No Weyl pair and no counterexample is constructed; DC1/JC2 untouched.**

---

## 9. Honest ledger

**Proved (arbitrary degree; degree-free proof objects in the verifier — no cap
anywhere in §§3–5, data entered as undetermined functions):**
1. Engine `Q_m = [D,X]_m`, `Q₀ = (T−1)G`, `G(0) = 0` under membership,
   `(T^n−1)` surjective with kernel the constants; the `Q₅` wall; the general
   `Q₄`/`Q₃` identities with no shape assumption (`§0, §1`).
2. The `(r, κ)` normalisation is an exact bijection on the `(Q₄,Q₃)` problem (`§2`).
3. The full degree-free node table for both classes (`§3`).
4. **diff-1 restoration:** the exponent-2 identity `144κ²a₂(r-1)² = …` with
   explicit integer cofactors, plus an independent second certificate and an
   independent Rabinowitsch certificate; hence `h h^{[1]} | a₂` and `h | b₁` (`§4`).
5. **diff-1 closure at arbitrary degree with no side condition**: this memo's
   `κ≠0` and `κ=0,b₁≠0` branches combine with Nonpositive-D Exclusion Theorem of
   [`shifted-cube-completion.md`](shifted-cube-completion.md) for `κ=0,b₁=0`.
6. Constant-top structure: wall, `Q₄` rung, vacuity of the potential
   factorisation, `Q₋₆ ⇒ b₋₃ = μ₃a₋₃`, the `μ₃`-source, the Lemma-P slope (`§8`).

**Refuted (machine-checked, with an explicit exact witness):**
- **diff-2 restoration.** `Q₄ ∧ Q₃` does **not** force `h h^{[1]} | a₂`, nor any
  proper factor beyond `(E-r)(E-r-1)`. Witness at symbolic `(r, κ)`, integer
  coefficients, whole scaling family (`§6`).
- **The corpus reading that a nonzero normal form for `a₂(r-1)` indicated
  non-restoration.** For diff-1 the normal form is nonzero *and* the value is
  forced; the ideal contains the square, not the element (`§4`).

**Bounded / finite evidence (exact scope stated):**
- diff-1, diff-2 full-cascade + `Q₀ = 1` + membership emptiness at `r = 0`,
  `κ = 1`, cap `d ≤ 3` committed (exact SymPy over `ℚ`, unit ideal), `d = 4`
  under `HEAVY=1`. Non-vacuity control in-file at `d = 2` (proper ideal without
  `Q₀ = 1`).
- A\*-band3 `κ₂ ≠ 0`: a **historical characteristic-zero msolve run** found cap `d=1`
  empty over `C̄` (`[-1]`), with controls (`κ₂=0` slice nonempty; positive cascade permits
  `κ₂≠0`; explicit tame pair). The finite-field unit is corroboration only. These msolve legs
  were **not rerun** in the current no-msolve environment; `d=2,3` remain uncompleted.
- Task C is narrowed to the formal constant-top A\*-band3 sector, `κ₂=1`, cap `d=1`.
  Historical characteristic-zero msolve found the slope Rabinowitsch system nonempty; the
  current run skipped this solver-dependent leg. No degree-free conclusion is claimed.

**Open / NOT claimed:**
1. **diff-2 at arbitrary degree with `κ ≠ 0`.** The surviving branch remains a
   tail problem; Nonpositive-D Exclusion Theorem of [`shifted-cube-completion.md`](shifted-cube-completion.md)
   closes only `κ = 0, b₁ = 0`, so no full diff-2 closure is claimed.
2. **A\*-band3 general `κ₂ ≠ 0` at arbitrary degree** remains open. Nonpositive-D Exclusion Theorem
   closes only the corner `κ₂ = 0, b₁ = 0`; the existing `κ₂ = 0` tame witness
   has positive band-one `D` (`b₁ ≠ 0`), so no contradiction with Nonpositive-D Exclusion Theorem.
   Composite tame-word escape; whether every `κ₂ = 0` pair is tame; all open.
3. Imbalanced coprime walls; general-`k` negative tail with `band D > 0`; **W2**;
   radical forcing at coupling widths `k = 4, 5`.

No Weyl pair, no counterexample; DC1/JC2 untouched.

---

## 10. Verification

```sh
uv run --with sympy python research/dc1-program/verify_band3_sectors.py
uv run --with sympy python research/dc1-program/verify_band3_sectors.py --require-msolve
HEAVY=1 uv run --with sympy python research/dc1-program/verify_band3_sectors.py --require-msolve
```

The ordinary mode treats unavailable or unsuccessful invoked msolve legs as optional SKIPs.
`--require-msolve` requires every invoked solver leg to complete with a recognized verdict;
timeout, nonzero exit, missing output, and malformed/unrecognized output fail the run. The
HEAVY-only legs are invoked, and hence required, only under `HEAVY=1`.

Current recorded no-`msolve` run: **145 passed, 0 failed, 3 skipped, wall time 87.0 s**
(`HEAVY=0`, `msolve=NOT RERUN/SKIPPED`); an independent `python -O` rerun recorded the
same status in **90.1 s**. The skipped groups were the HEAVY cap-`d=4` leg, live msolve parser
validation plus every load-bearing msolve-dependent proof leg, and the formal constant-top
Task C msolve probe. No solver-backed conclusion was rerun. `§0` engine; `§1` wall +
general rungs; `§2` normalisations; `§3` degree-free node table; `§4` diff-1 radical certificate; `§5` diff-1 closure and the
symbolic `Q₂` congruence; `§6` diff-2 `(Q₄,Q₃)`-locus witness; `§7` bounded controls; `§8`
Task B and strict parser; `§9` bounded Task C. Runtime is environment-dependent. The final
banner excludes skipped checks from unconditional success statements and lists every skip.
