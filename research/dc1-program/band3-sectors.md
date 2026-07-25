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
> **The diff-1 sector is EMPTY at arbitrary degree for `κ ≠ 0`, and also for
> `κ = 0` with `b₁ ≠ 0`.**
>
> **diff-2 — RESTORATION IS REFUTED.** Explicit exact solution of `Q₄ = Q₃ = 0`
> at **symbolic `(r, κ)`**:
> `a₂ = 3(E-r)(E-r-1)`, `b₁ = 2κ(E-r-1)²`, `a₁ = -3`, `b₀ = 0`, with
> `a₂(r-1) = a₂(r+2) = 6 ≠ 0`. The cascade forces **exactly** the degraded
> `(E-r)(E-r-1) | a₂` and nothing more; `Q₃` contributes only the branch
> disjunction `a₂(r-1) · b₁(r+1) = 0`. **The clean-divisibility route to diff-2 is
> dead — a different mechanism is required.** diff-2 remains open at arbitrary
> degree, with bounded emptiness at cap `d ≤ 3` (committed) and `d = 4` (HEAVY).

Exact certificate: [`verify_band3_sectors.py`](verify_band3_sectors.py) —
**117 checks, all passing, ~43 s default run** (3 `HEAVY` skips, all listed); heavier legs behind `HEAVY=1`.
Every load-bearing upstream fact (the crossed-product ladder engine
`Q_m = [D,X]_m`, `Q₀ = (T−1)G`, `G(0)=0` under membership, the `Q₅` wall, the
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
> **The diff-1 broken class is closed** — *for `κ ≠ 0`, and (§4.1) for `κ = 0` with
> `b₁ ≠ 0`. The sub-branch `κ = 0` **and** `b₁ = 0` remains open, and is inherited and
> `h`-independent: it affects the cube-separated classes equally (§4.1).*

### 4.1 The `κ = 0` branch, and one inherited gap

`κ = 0` (`b₂ = 0`): `Q₄ = h^{[1]}h^{[2]}·W` with `W = h b₁^{[3]} − h^{[3]}b₁`, and
`W(r) = −6b₁(r)`, `W(r+1) = −12b₁(r+1)`, so `h | b₁`; writing `b₁ = hψ` gives
`ψ^{[3]} = ψ`, hence `ψ = c` constant. For `c ≠ 0`, `Q₃ = 0` forces `a₃ | F` with
`F = h^{[2]}a₂ − h a₂^{[1]}`, and `F(r) = 2a₂(r)`, `F(r+1) = 6a₂(r+1)`,
`F(r-2) = −6a₂(r-1)`, `F'(r) = 3a₂(r) + a₂(r+1) + 2a₂'(r)` — the same clean jet.
`Q₂ = 0` then forces `h h^{[1]} | V = h^{[1]}a₁ − h a₁^{[1]}`, and
`V(r+1) = 2a₁(r+1)`, `V(r-1) = −2a₁(r)`, so `h | a₁`. The closure of §4 runs.
**All degree-free (verifier `§5`).**

> **Named residual, and it is INHERITED, not new.** The sub-branch `κ = 0` **and**
> `b₁ = 0` is not closed here. In it `Q₃ = a₃(T³−1)b₀` only forces `b₀` constant
> and puts no condition on `a₂`. This gap is **`h`-independent**: the corpus's
> cube-separated `κ = 0` chain ([`shifted-power-residuals.md`](shifted-power-residuals.md)
> §1.2, "`Q₄ ⇒ b₁ = c h`, `Q₃ ⇒ h h^{[1]} | a₂`") silently assumes `c ≠ 0` in the
> second step. **This memo flags that as a corpus gap; it is not introduced by
> diff-1 and it is not closed by this memo.**

---

## 5. diff-2 — restoration is REFUTED (explicit exact witness)

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

> **Corollary (refutation).** For diff-2 the cascade `Q₄ ∧ Q₃` forces **exactly**
> the degraded `(E-r)(E-r-1) | a₂` of `broken-separation.md` §2.2 and nothing
> more. The clean `h h^{[1]} | a₂` is **not** restored, and no proper factor
> beyond `(E-r)(E-r-1)` is forced either. **The clean-divisibility route to
> diff-2 is closed off; a different mechanism is required.**

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

**Controls, committed:**
- **`κ₂ = 0` slice nonempty, explicit point.** `U = x + 2∂`, `X = U³ − ∂`,
  `D = U`. The verifier *computes* `[D,X]` in the crossed product and checks it
  equals `1` with `a₃ = 1`, `b₃ = 0`, `b₂ = 0` and genuine membership.
- **msolve parser validated in-file**, before any load-bearing call: a known
  **unit** ideal `(x−1, x−2)` must report `[-1]` (empty) and a known **feasible**
  ideal `(x²−2, y−x)` must report nonempty. Both trap guards are **unit-tested**
  (a body carrying `2*x^2/3` and a body carrying `**` are each rejected; a clean
  integer body passes) — no check in this file is true by construction.

**Engine finding (recorded for the corpus, and it corrects the previous
diagnosis).** SymPy's Gröbner engine does **not** finish the A\*-band3 cap-`d = 1`
system within 20 minutes, although the *rigid* diff-1/diff-2 sectors at `d = 3`
(44 variables) finish in ~4 s; the constant top `a₃ = 1` removes all rigidity.
**msolve does the same cap-`d = 1` system in 0.1 s.** So the binding constraint on
the previous attempt was **engine selection**, not (as previously supposed)
rational-coefficient misparsing — trap #2 is real and is guarded here, but it was
not what made `d = 3` unreachable. Every A\*-band3 emptiness leg in this verifier
is msolve-only.

**Results (run of record):**

| leg (cap `d`, `κ₂ = 1`, integer system) | engine / tier | result |
|---|---|---|
| `d = 1` emptiness over `ℚ` | msolve char-0, **committed default** (0.1 s) | `[-1]` = **EMPTY over `C̄`** |
| `d = 1` `κ₂ = 0` slice | msolve char-0, committed (0.3 s) | `[1,22,-1,[]]` = **NONEMPTY** (positive-dimensional) |
| `d = 1` unit ideal mod `p = 2³⁰+3` | msolve `-g`, committed (0.1 s) | `[1]` — **corroboration only** |
| `d = 1` positive cascade alone | explicit point, instant | `a₂ = 5`, `b₁ = 7` const ⇒ **NONEMPTY** with `κ₂ = 1` |
| `d = 1` drop `Q₀ = 1` | msolve char-0, `HEAVY=1`, 1200 s cap | slow (positive-dimensional); redundant with the two controls above |
| `d = 2` (branch `a₋₃ = 0` and full) | msolve char-0, `HEAVY=1`, 3000 s cap | **not completed in this session** (>150 s probe; the corpus records ≈35 min) |
| **`d = 3`** (branch `a₋₃ = 0` and full) | msolve char-0, `HEAVY=1`, 3000 s cap | **NOT ACHIEVED** — the leg is implemented and gated, and prints SKIP if it times out |

> **Honest statement of Task B.** The `d = 3` target is **not** met. What this
> memo adds over [`astar-band3.md`](astar-band3.md) is: (i) the `d = 1` char-0
> `ℚ` emptiness is now a **committed default-run** certificate rather than an
> "additional/manual" computation; (ii) the `κ₂ = 0` nonemptiness at the *same*
> cap and the positive-cascade nonemptiness are committed **controls**, so the
> `d = 1` exclusion is a genuine separation and not a broken encoding; (iii) the
> msolve parser is **validated in-file** (known unit ⇒ `[-1]`, known feasible ⇒
> nonempty) and both traps are **unit-tested**; (iv) the branch split
> (`a₋₃ = 0` / full) is implemented for `d = 2, 3`; (v) the real obstacle is
> named — engine, not encoding.

The `F_p` row is **corroboration only**. Unit mod `p` does *not* imply unit over
`ℚ` — `(p x − 1)` is a counterexample — and the verifier says so at the check.

## 7. TASK C — the slope-forcing probe: the route is DEAD at `d = 1`

The exact W2 analogue for the constant top: is the Lemma-P moment slope `G(1)`
forced to `0` by the **tail alone** (all `Q_m = 0` for `m ≠ 0`, membership,
`κ₂ = 1`, **no** `Q₀` condition)? If it were, `Q₀ = 1` — which needs `G(1) = 1` —
would be impossible, and that would be a degree-free target.

Radical-correct test (verifier `§9`, committed at `d = 1`, ~34 s): Rabinowitsch on
`{tail} ∪ {1 − t·G(1)}`. msolve char-0 returns `[1,23,-1,[]]` — a
**positive-dimensional, hence nonempty** variety.

> **TASK C VERDICT (cap `d = 1`, `κ₂ = 1`): `G(1)` is NOT forced to `0` by the
> tail alone.** There is a genuine point of the tail with `G(1) ≠ 0`. This is an
> *existence* statement, so it is radical-correct as a negative: **the slope route
> is dead at this cap**, and the constant-top obstruction does **not** reduce to
> the moment slope. `d = 2` is implemented behind `HEAVY=1`.

This is a useful negative: it tells the program not to look for the `(κ₂-closure)`
kill in the moment-slope covector, which is where the W2 analogy pointed.

## 8. What this changes in the band-3 ledger

| class | status before | **status now** |
|---|---|---|
| cube-separated / 2-separated `h` (incl. diff-3, `(E-r)²`) | closed, arbitrary degree | unchanged |
| **diff-1** `(E-r)(E-r-1)` | degraded `Q₄` forcing; restoration **unresolved**; bounded emptiness only | **CLOSED at arbitrary degree** (`κ ≠ 0`; and `κ = 0, b₁ ≠ 0`). Restoration **forced**, radical certificate exponent 2 |
| **diff-2** `(E-r)(E-r-2)` | degraded `Q₄` forcing; restoration **unresolved** | restoration **REFUTED** (explicit exact witness at symbolic `r, κ`); forced divisor is **exactly** `(E-r)(E-r-1)`; **open** at arbitrary degree, bounded-empty `d ≤ 4` |
| `κ = 0 ∧ b₁ = 0` sub-branch | (not previously isolated) | **named open gap, inherited and `h`-independent** — also affects the corpus's cube-separated `κ = 0` chain |

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
5. **diff-1 closure at arbitrary degree** for `κ ≠ 0`, and for `κ = 0` with
   `b₁ ≠ 0`: `h | a₁`, `G = h^{[-1]}M`, `Q₀ = 1 ⇒ h^{[-1]} | E` impossible (`§5`).
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
- A\*-band3 `κ₂ ≠ 0`: cap `d = 1` **EMPTY over `ℚ`** (msolve char-0 `[-1]`,
  committed default, 0.1 s), with committed controls (`κ₂ = 0` slice nonempty at
  the same cap; positive cascade permits `κ₂ ≠ 0` by an explicit point; explicit
  tame pair `[D,X] = 1`). `F_p` unit is corroboration only. **`d = 2` and `d = 3`
  are implemented and `HEAVY`-gated but were NOT completed here.**
- Task C: `G(1)` is **NOT** forced to `0` by the tail alone at cap `d = 1`
  (Rabinowitsch system nonempty, committed, 27 s) — the slope route is dead at
  that cap. `d = 2` `HEAVY`-gated, not run.

**Open / NOT claimed:**
1. **diff-2 at arbitrary degree.** The clean-divisibility route is now *refuted*,
   so the surviving branch `b₁(r+1) = 0` needs a genuinely different argument —
   the natural next targets are the `Q₂` derivative-node system on that branch and
   a `G`-factorisation through a divisor other than `h^{[-1]}`.
2. **`κ = 0 ∧ b₁ = 0`** in *every* shifted-cube class, cube-separated included.
3. A\*-band3 `(κ₂-closure)` at arbitrary degree; whether `κ₂` survives arbitrary
   composite tame words; whether every `κ₂ = 0` pair is tame.
4. Imbalanced coprime walls; general-`k` negative tail; **W2**.

No Weyl pair, no counterexample; DC1/JC2 untouched.

---

## 10. Verification

```sh
uv run --with sympy python research/dc1-program/verify_band3_sectors.py
HEAVY=1 uv run --with sympy python research/dc1-program/verify_band3_sectors.py
```

Default run of record: **117 checks executed, 117 passed, 0 failed, 3 skipped,
wall time 43.0 s** (`HEAVY=0`, `msolve=yes`). `§0` engine; `§1` wall + general rungs; `§2` normalisation; `§3`
degree-free node table; `§4` diff-1 radical certificate; `§5` diff-1
arbitrary-degree closure (+ the `κ = 0` branches); `§6` diff-2 witness; `§7`
bounded emptiness + non-vacuity controls; `§8` Task B (msolve parser validation,
`κ₂ = 0` control, `HEAVY` emptiness legs); `§9` Task C probe. Runtime is
environment-dependent. The final banner distinguishes *all checks passed* from
*all executed checks passed; optional checks skipped*, and lists every skip with
its reason.
