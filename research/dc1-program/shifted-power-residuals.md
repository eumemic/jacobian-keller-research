# Shifted-power residuals at band 3: h | a_1 derived, 2-separation, and the kappa_2 invariant

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED — PARTIAL RESULT WITH NAMED GAPS**

This memo closes and sharpens the three named residuals of the band-3
shifted-cube wall — the band-3 instance of band-reduction **Gap 2** — left open
in [`shifted-power-descent.md`](shifted-power-descent.md) §6/§7. It supplies one
of the two inputs the conditional band-reduction framework
([`band-reduction.md`](band-reduction.md) §9, post-`727ce8a` form) is waiting on.

The one-line summary of each residual:

> **RESIDUAL 1 (the sharp one) — DISCHARGED.** The divisibility `h | a₁`,
> previously *imported* from [`../band3/quantum-shifted-cube.md`](../band3/quantum-shifted-cube.md)
> §2.2 (the audit-flagged conditionality of `shifted-power-descent.md` §3), is now
> **derived in-file** from the ladder rung `Q₃`: with the `Q₄`-forced shapes
> `a₂ = h h^{[1]} g`, `b₁ = h β` substituted, `Q₃ = 0` with a **polynomial `b₀`**
> forces `h | a₁` — **exactly, with no slack** — for `κ ≠ 0` (and via `Q₂` for
> `κ = 0`). The closed-form factorization `G = h^{[-1]} M` therefore holds with
> **every** divisibility derived; the conditionality is discharged.
>
> **RESIDUAL 2 (non-cube-separated) — one class CLOSED arbitrary-degree.** The
> h-forcing needs only **2-separation** `gcd(h,h^{[1]}) = gcd(h,h^{[2]}) = 1`,
> **not** cube-separation. The `b₁` step needs *no coprimality at all* (a rational
> `ψ` with `(T³−1)ψ` polynomial is polynomial); `a₂` and `a₁` need only `j = 1,2`.
> Hence the entire **diff-3 class** (roots differing by exactly `3`, so `2`-separated
> but not cube-separated) is **forced `h` constant at arbitrary degree** — one of
> the three residual classes is eliminated. The genuinely broken classes **diff-1**
> and **diff-2** (the `a₂` divisibility fails) have committed bounded emptiness at
> cap `D = 2`; cap `D = 3,4` is an optional `HEAVY` `msolve` payload (`κ` free).
>
> **RESIDUAL 3 (constant-h, `κ₂ ≠ 0`) — one-step obstruction only.** No listed
> gauge-preserving generator removes `κ₂`: `D ↦ D − λX` preserves `b₃ = 0` only for
> `λ = 0`; higher transvections overshoot to band `≥ 6`; pair-exchange relabels the
> band-3 top; Fourier breaks membership. This does **not** establish invariance under
> an arbitrary composite tame word, which may leave and later return to the displayed
> gauge/band sector. The positive cascade permits `κ₂ ≠ 0`; both composite tame escape
> and the **A\*-band3 negative-tail closure** remain open.

Exact certificate:
[`verify_shifted_power_residuals.py`](verify_shifted_power_residuals.py) — default
run exact SymPy over QQ (`ALL SHIFTED POWER RESIDUAL CHECKS PASSED`); heavier
`msolve` corroboration behind `HEAVY=1`. Every load-bearing upstream fact (the
crossed-product ladder engine `Q_m = [D,X]_m`, `Q₀ = (T−1)G`, the closed-form
potential, the `Q₄` shapes) is **re-derived in file**, not merely cited.

Conventions frozen from the corpus: `A₁[x^{-1}] = ⊕_k x^k C[E]`, `E = x∂`,
`(x^a f)(x^b g) = x^{a+b} f(E+b) g(E)`, `f^{[n]}(E) = f(E+n)`, `T f = f^{[1]}`,
`S_n = 1+T+⋯+T^{n-1}`, `Q_m = ∑_{k+l=m}[b_l^{[k]}a_k − a_k^{[l]}b_l]`,
`[D,X]=1 ⇔ Q_m=δ_{m0}`, membership `(E)_j = E(E−1)⋯(E−j+1) | a_{-j},b_{-j}`, gauge
`b_k = 0`, `G = ∑_{k≥1}∑_{j=0}^{k-1}(a_k^{[j-k]}b_{-k}^{[j]} − b_k^{[j-k]}a_{-k}^{[j]})`,
`Q₀ = (T−1)G`. Sector: gauge `b₃ = 0`, `a₃ = h(E)h(E+1)h(E+2)`,
`b₂ = κ·h(E)h(E+1)` (the `Q₅` wall).

---

## 1. RESIDUAL 1 — `h | a₁` derived in-file from `Q₃`

### 1.1 The rung

In gauge `b₃ = 0`, the four `(k,l)` with `k+l = 3` and `k,l ∈ [−3,3]` (with
`b₃ = 0` killing `(0,3)`) give the exact three-term identity

```
Q₃ = a₃·[ (T³−1)b₀ + (β^{[2]}g − g^{[1]}β) ]  +  κ h^{[1]}( h^{[2]} a₁ − h a₁^{[2]} ),
```

after substituting the `Q₄`-forced shapes `a₂ = h h^{[1]} g`, `b₁ = h β`,
`b₂ = κ h h^{[1]}` (`a₃ = h h^{[1]} h^{[2]}`). The verifier checks this identity
against the engine at symbolic `h` (`§1`, arbitrary degree).

### 1.2 The forcing

The first bracket carries the full factor `a₃ = h h^{[1]} h^{[2]}`. A coefficient
of `D` is a genuine polynomial, so `b₀ ∈ C[E]`; hence `a₃ (T³−1)b₀` is divisible
by `a₃`. Rearranging `Q₃ = 0`,

```
a₃·[ (T³−1)b₀ + (β^{[2]}g − g^{[1]}β) ]  =  −κ h^{[1]}( h^{[2]} a₁ − h a₁^{[2]} ),
```

the left side is divisible by `a₃`, so the right side is too:
`a₃ | κ h^{[1]}( h^{[2]} a₁ − h a₁^{[2]} )`. Cancel `h^{[1]}` (nonzero):

```
h h^{[2]}  |  κ( h^{[2]} a₁ − h a₁^{[2]} ).
```

Reduce mod `h`: `h^{[2]} a₁ − h a₁^{[2]} ≡ h^{[2]} a₁`, so `h | κ h^{[2]} a₁`; with
`gcd(h, h^{[2]}) = 1` this is `h | κ a₁`, and for `κ ≠ 0` (char 0, `κ` a unit)

```
   h | a₁      (write a₁ = h·p,  p ∈ C[E] free).
```

**No slack.** The single reduction mod `h` already gives the full divisibility; the
`h^{[2]}`-reduction is redundant. The derivation uses only `gcd(h, h^{[2]}) = 1`
(the `j = 2` separation), nothing about `j = 1` or `j = 3`.

The verifier machine-checks the forcing as **necessity** at two genuinely
cube-separated tops — `h = E(2E−1)` (roots `{0, ½}`) and `h = (3E−1)(4E−1)` (roots
`{⅓, ¼}`): solving `Q₃ = 0` for a *fully generic* `a₁` and *fully generic* `b₀`
returns a locus on which `h | a₁` holds and the quotient `p` is free and nonzero
(`§1`). The `κ = 0` branch is handled separately: `b₂ = 0`, `Q₄ ⇒ b₁ = c h`,
`Q₃ ⇒ h h^{[1]} | a₂`, and then `Q₂ = 0` forces `h | a₁` (re-derived at symbolic
`h`, `§1`).

### 1.3 The conditionality is discharged

With `h | a₁` **derived**, the three sub-leading divisibilities
`h h^{[1]} | a₂`, `h | b₁`, `h | a₁` are all forced by the ladder itself
(`Q₄, Q₃`), so the closed-form central factorization

```
   G = h^{[-1]} · M,     M(0) = 0   (membership),
```

holds with **no imported hypothesis** (`§1`, re-verified at linear `h`). The
affine kill (`Q₀ = 1 ⇒ G = E ⇒ h^{[-1]} | E ⇒ deg h ≤ 1`, and `deg h = 1` killed by
`M(0) = 0` vs `M = 1/α`) then runs exactly as in
[`shifted-power-descent.md`](shifted-power-descent.md) §3. **The §3 conditionality
note is retired: "`Q₀ = 1 ⇒ h` constant (cube-separated)" is now unconditional on
the ladder.**

---

## 2. RESIDUAL 2 — 2-separation suffices; the diff-3 class is closed

### 2.1 The `b₁` step needs no coprimality

The `shifted-power-descent`/`quantum-shifted-cube` derivation of `h | b₁` was
stated using `gcd(h, h^{[3]}) = 1`. It does **not** need it. After `a₂ = h h^{[1]} g`
is established, `Q₄ = 0` reduces to

```
   b₁^{[3]} h − h^{[3]} b₁ = κ h h^{[3]}(g − g^{[2]}).
```

Divide by `h h^{[3]}` (both nonzero): with `ψ := b₁/h`,

```
   (T³ − 1) ψ = −κ (T² − 1) g       (the S₃ quantum midpoint for ψ).
```

The right side is a polynomial, and:

> **Lemma (denominator killing).** If `ψ ∈ C(E)` and `(Tⁿ − 1)ψ ∈ C[E]`, then
> `ψ ∈ C[E]`.

*Proof.* Write `ψ = N/D` in lowest terms. `(Tⁿ−1)ψ = (N^{[n]}D − N D^{[n]})/(D D^{[n]})`
polynomial forces `D | N^{[n]}D − N D^{[n]}`, hence `D | N D^{[n]}`, hence (as
`gcd(N,D)=1`) `D | D^{[n]}`. Equal degrees give `D^{[n]} = c·D`, i.e. the root
multiset of `D` is invariant under translation by `n` — impossible for a nonempty
finite multiset unless `D` is constant. ∎

So `ψ = b₁/h` is a polynomial: **`h | b₁`, using no separation hypothesis.** The
verifier checks the lemma's engine (no nonconstant divisor `D | D^{[3]}` for the
diff-3 tops; and the mechanism `(T³−1)(1/E) = 1/(E+3) − 1/E ∉ C[E]`) in `§2`.

### 2.2 Which separations each divisibility needs

| divisibility | rung | coprimality used |
|---|---|---|
| `h h^{[1]} | a₂` | `Q₄` | `gcd(h,h^{[1]})=1` **and** `gcd(h,h^{[2]})=1` (`j=1,2`) |
| `h | b₁` | `Q₄` | **none** (the lemma above) |
| `h | a₁` | `Q₃` | `gcd(h,h^{[2]})=1` (`j=2`) |

So the whole cascade — hence `G = h^{[-1]}M`, `M(0)=0`, and the affine kill — needs
only **2-separation** `gcd(h,h^{[1]}) = gcd(h,h^{[2]}) = 1`. Cube-separation
(`j=1,2,3`) was never necessary.

### 2.3 The diff-3 class is closed (arbitrary degree)

A **diff-3** top has two roots differing by exactly `3` (necessarily `deg h ≥ 2`)
but no roots differing by `1` or `2` — it is `2`-separated but **not**
cube-separated. By §2.2 the full cascade runs: `Q₄` forces `h h^{[1]} | a₂` and
`h | b₁`, `Q₃` forces `h | a₁`, `G = h^{[-1]}M` holds term-by-term with `M(0)=0`,
and `Q₀ = 1` forces `h^{[-1]} | E`, impossible for `deg h ≥ 2`. **Hence no
nonconstant diff-3 top carries a genuine pair.**

The verifier machine-checks the whole closure at **three** diff-3 tops —
`E(E−3)`, `(E−1)(E−4)`, `(2E−1)(2E−7)` — as necessity solves plus the
`G = h^{[-1]}M` identity plus the affine kill (`§2`). The arbitrary-degree content
is the coprimality/denominator prose of §2.1–§2.2; the machine checks are the
per-top necessity, exactly the tier of the original cube-separated `h`-forcing.

> **The closed shifted-cube class enlarges from `cube-separated` to `2-separated`
> `h`**, adding the entire diff-3 family. This eliminates one of the three residual
> classes of `shifted-power-descent.md` §7.1.

### 2.4 diff-1, diff-2 genuinely break — stronger bounded evidence

For **diff-1** (`gcd(h,h^{[1]}) ≠ 1`, e.g. `E(E−1)`) and **diff-2**
(`gcd(h,h^{[2]}) ≠ 1`, e.g. `E(E−2)`, `(E−1)(E−3)`) the `a₂` step genuinely fails:
solving `Q₄ = 0` with generic `a₂, b₁` gives a locus on which `h h^{[1]} ∤ a₂`
(the slack that removes the closed-form factorization). The verifier exhibits this
failure (`§2`). Only **bounded** evidence excludes these tops:

- **committed (exact SymPy, cap `D = 2`, `κ` symbolic → all `κ`):** the positive
  cascade `Q₄=Q₃=Q₂=Q₁=0` with `Q₀=1` and membership is the **unit ideal** for
  `E(E−1)` and `E(E−2)` (`§2`).
- **HEAVY (`msolve`, `κ` a free variable → all `κ`):** the same sector is the unit
  ideal at **cap `D = 3` and cap `D = 4`** for `E(E−1)`, `E(E−2)`, `(E−1)(E−3)`
  (behind `HEAVY=1`; SKIPs cleanly if `msolve` absent).

This is corroboration, **not** an arbitrary-degree exclusion; diff-1/diff-2 remain
a genuine sub-wall, open at arbitrary degree, now with cap-`D=4` emptiness.

---

## 3. RESIDUAL 3 — the constant-h gauged wall constant `κ₂`

With `h` constant, normalize `a₃ = 1`; the `Q₅` wall gives `b₂ = κ₂` (a constant).
The tame family
`U = x+c₀+c₁∂`, `X = U³ − ∂/κ − A`, `D = λX + κU + β` lives in the gauged
`κ₂ = 0` branch (`D ↦ D − λX` is the **B0-band3 collapse** killing `b₃, b₂, b₋₃`
at once).

### 3.1 One-step gauge obstruction (arbitrary degree; composite words open)

The only transvection `D ↦ D − p'(X)` reaching band `2` with band `≤ 3` is
`p'(X) = λX` (band `3`, since `a₃ = 1`). It sends

```
   b₃ ↦ −λ a₃ = −λ,        b₂ ↦ κ₂ − λ a₂.
```

Preserving the gauge `b₃ = 0` forces `λ = 0`, leaving `b₂ = κ₂` **unchanged**. Any
higher transvection `p'(X) = X², …` has band `≥ 6` (overshoots band `2`).
Pair-exchange `(X,D) ↦ (D,−X)` keeps `band = max(2,3) = 3` — it merely relabels
which generator carries the band-3 top. Fourier `E ↦ −E−1` sends the constant top
`a₃ = 1` to an `x^{-3}`-coefficient `1`, **not** divisible by `E(E−1)(E−2)`,
breaking `A₁`-membership (as in [`../band3/astar-band3.md`](../band3/astar-band3.md)
§6). All four facts are machine-checked (`§3`).

> **Scope.** These checks obstruct each displayed move when applied directly while
> retaining the stated low-band/gauge description. They do not exclude a composite
> tame word that temporarily exits this sector and later returns with `κ₂=0`. Thus
> `κ₂` is not proved to be a tame-orbit invariant here.

### 3.2 The positive cascade permits `κ₂ ≠ 0`

The constant-top positive rung is `Q₄ = κ₂(a₂ − a₂^{[2]}) + (b₁^{[3]} − b₁)`
(exact identity, `§3`), which has `κ₂ ≠ 0` solutions (e.g. `a₂ = E`,
`b₁ = ⅔κ₂E`). So `κ₂ ≠ 0` is **not excluded at the positive level**; its
disposition is a genuine **negative-tail** question.

A bounded emptiness probe of the `κ₂ ≠ 0` sector is inconclusive: with `κ₂`
fixed nonzero the positive cascade `+ Q₀=1` does *not* collapse to the unit ideal
(the Gröbner basis stays large; the exact full-tail sweep is too costly at the
caps attempted, matching the identical remark for the mirror A\* sector in
[`../band3/astar-band3.md`](../band3/astar-band3.md) §8). **No emptiness claim is
made here.** The honest disposition is §4 point 4: the direct gauge-reduction move
stalls, while composite tame escape and the A\*-band3 negative-tail closure are open.

---

## 4. Bookkeeping — what Gap 2's band-3 instance still needs

After this memo, the band-3 shifted-cube branch of band-reduction **Gap 2** stands
as follows (in the conditional-framework terms of `band-reduction.md` §9).

**Discharged / newly closed:**

1. **RESIDUAL 1 (h | a₁):** derived in-file from `Q₃` (arbitrary degree, no slack).
   The propagation `Q₀ = 1 ⇒ h` constant is now **unconditional on the ladder** for
   `2`-separated `h` — the audit-flagged conditionality of `shifted-power-descent.md`
   §3 is retired.
2. **RESIDUAL 2, diff-3 class:** closed at arbitrary degree; the closed
   shifted-cube class is **`2`-separated `h`**, not merely cube-separated.

**Still needed for full band-3 shifted-power closure (open, precisely delimited):**

3. **RESIDUAL 2, diff-1 & diff-2 classes** (`gcd(h,h^{[1]}) ≠ 1` or
   `gcd(h,h^{[2]}) ≠ 1`): the `a₂` divisibility genuinely fails and `G = h^{[-1]}M`
   is unavailable. Bounded emptiness only: cap `D = 2` committed, with cap `D = 3,4`
   optional under `HEAVY` when `msolve` is available. A new sub-wall, open at arbitrary degree.
4. **RESIDUAL 3, `κ₂ ≠ 0`:** the direct gauge-reduction move stalls (§3.1), but
   composite tame-word escape is not excluded. Independently, the exact tail statement is:

   > *(κ₂-closure)* There is **no** genuine Weyl pair `[D,X] = 1` with `a₃ = 1`
   > (constant top), gauge `b₃ = 0`, `b₂ = κ₂ ≠ 0`, and membership-valid negative
   > tail solving `Q_{-1} = ⋯ = Q_{-6} = 0`, `Q₀ = 1`.

   This is the **top-mirror of the open quantum A\*-I / negative-tail closure**:
   `b₂ = κ₂ ≠ 0` (a top-wall constant) is the mirror, under the falling-factorial
   reflection `E ↦ −E−1`, of A\*'s bottom-wall residue `μ̃ ≠ 0`. Reflection breaks
   `A₁`-membership for the constant top `a₃ = 1` (§3.1, exactly the obstruction of
   [`../band3/astar-band3.md`](../band3/astar-band3.md) §6 (i)), so the classical
   reflection route does not transcribe and the sector needs its own negative-tail
   analysis — the same inhomogeneous `μ`-source tail that keeps quantum A\*-I open.
   Supplying *(κ₂-closure)* is exactly what remains of RESIDUAL 3.

Beyond band 3, unchanged: the **imbalanced coprime walls** (`q ∤ k`, `q ≥ 2`) reduce
only modulo the composite-move escape (**Gap 1** = the DC1 core); the general-`k`
negative tail is untouched; **W2** remains open independently. **No Weyl pair and no
counterexample is constructed; DC1/JC2 untouched.**

---

## 5. Honest ledger — proved / bounded / open

**Proved (exact algebra, machine-checked identities; arbitrary degree where
stated):**
- Engine `Q_m = [D,X]_m`, `Q₀ = (T−1)G` (`§0`).
- **RESIDUAL 1:** the exact three-term `Q₃` identity; the coprimality forcing
  `h | a₁` (arbitrary degree, `gcd(h,h^{[2]})=1`, `κ ≠ 0`; `κ = 0` via `Q₂`),
  machine-verified as necessity at two cube-separated tops; the chain closes
  `G = h^{[-1]}M`, `M(0)=0` with `h|a₁` derived. **No slack.**
- **RESIDUAL 2:** the denominator-killing lemma (`h | b₁` needs no coprimality);
  the separation table (`a₂`: `j=1,2`; `a₁`: `j=2`; `b₁`: none); the **diff-3
  closure** at arbitrary degree (`2`-separated `h`), machine-verified at three
  diff-3 tops.
- **RESIDUAL 3:** the one-step gauge obstruction (none of the displayed direct
  transvection / pair-exchange / Fourier moves removes `κ₂` while retaining the sector);
  the constant-top `Q₄` identity and a `κ₂ ≠ 0` positive solution.

**Bounded / finite evidence (exact scope):**
- diff-1, diff-2 sector emptiness: cap `D = 2` committed (exact SymPy, all `κ`),
  cap `D = 3,4` HEAVY (`msolve`, all `κ`).
- `κ₂ ≠ 0` constant-h sector: **no emptiness claim** — the bounded probe is
  inconclusive (Gröbner too costly, mirroring `astar-band3.md` §8). Only the direct
  one-step gauge obstruction and positive-cascade consistency are established.

**Refuted (machine-checked) — a corpus correction:**
- The diff-3 tops (e.g. `E(E−3)`), listed as an **open** arbitrary-degree residual
  in [`../band3/quantum-shifted-cube.md`](../band3/quantum-shifted-cube.md) §7/§10
  and [`shifted-power-descent.md`](shifted-power-descent.md) §7.1, are in fact
  **closed**: the `h`-forcing never needed `gcd(h,h^{[3]}) = 1`, only
  `2`-separation. The cube-separation hypothesis was not tight.

**Open / NOT claimed:**
1. diff-1, diff-2 non-separated tops at arbitrary degree (bounded evidence only).
2. `κ₂ ≠ 0` constant-h: composite tame-word escape and the A\*-band3 negative-tail closure.
3. Imbalanced coprime walls (`q ∤ k`, `q ≥ 2`); general-`k` negative tail; **W2**.
No Weyl pair, no counterexample; DC1/JC2 untouched.

---

## 6. Verification

```sh
uv run --with sympy python research/dc1-program/verify_shifted_power_residuals.py
HEAVY=1 uv run --with sympy python research/dc1-program/verify_shifted_power_residuals.py
```

Exact SymPy over QQ: `§0` engine; `§1` RESIDUAL 1 (`Q₃` identity, `h|a₁` forcing +
necessity, `κ=0` via `Q₂`, chain closes); `§2` RESIDUAL 2 (denominator lemma,
diff-3 closure ×3 tops, diff-1/diff-2 slack + cap-`D=2` emptiness); `§3`
RESIDUAL 3 (gauge obstruction, positive-cascade permits `κ₂≠0`). Behind `HEAVY=1`,
`msolve` corroboration at cap `D = 3,4` for diff-1/diff-2/diff-2b. A successful run
ends `ALL SHIFTED POWER RESIDUAL CHECKS PASSED`.
