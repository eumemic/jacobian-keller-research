# Band-k weapons: lifting the band-2 arsenal to general band index

**INDEPENDENTLY DERIVED — MACHINE-CHECKED IDENTITIES (k = 2, 3, 4) — NOT PEER REVIEWED**

This memo lifts the band-2 arsenal (W1–W9) to a symbolic band index `k`, proving
what generalizes verbatim, and — this is the point — pinning down precisely
where the lift **breaks**. Two mechanisms generalize as clean arbitrary-`k`
theorems (W1 top Wronskian; W4 the m=0 moment). One generalizes classically but
acquires a genuinely new solution class quantumly at `k ≥ 3` (W2, the
gatekeeper — the sharpest finding here). One (W6, the bottom staggered rigidity)
is shown to be the **mirror of W2**, present in *both* faces, correcting an
over-strong "quantum-only" framing in the band-2 corpus. W5 (the lattice
conjecture) is given a precise formulation with the band-2 instance proved as a
theorem of the framework.

All conventions are frozen exactly as in the band-2 corpus:

- **Classical.** `C[x,ξ]`, `τ = xξ`, `{G,F} = G_ξ F_x − G_x F_ξ`, `{ξ,x}=1`.
  `F = Σ_{i=−k}^{k} x^i a_i(τ)`, `G = Σ_l x^l b_l(τ)`; `{x^l b, x^i a} =
  x^{i+l}(i a b′ − l a′ b)`; primes `= d/dτ`; **membership** `τ^j | a_{−j},
  b_{−j}`; ladder equations `C_m = Σ_{i+l=m}(i a_i b_l′ − l a_i′ b_l) = δ_{m0}`.
- **Quantum.** `A_1[x^{−1}] = ⊕_k x^k C[E]`, `E = x∂`, `(x^a f(E))(x^b g(E)) =
  x^{a+b} f(E+b) g(E)`, `f^{[r]}(E) = f(E+r)`; **membership** `E(E−1)⋯(E−j+1) |
  a_{−j}, b_{−j}`; ladder equations `Q_m = Σ_{i+l=m}[b_l^{[i]} a_i − a_i^{[l]}
  b_l] = δ_{m0}`.

Everything below is checked exactly by
[`verify_bandk_weapons.py`](verify_bandk_weapons.py) at `k = 2, 3, 4` (94 `PASS`
lines, ends `ALL BAND-K WEAPON CHECKS PASSED`). Citations to the band-2 corpus
are pinned to commits: band-2 assembly `84978b9`, full classical `f637b1a`,
partial cascades `91a053a`, classical A\* closure and quantum mirror in the
`84978b9` tree (`classical-Astar.md`, `quantum-mirror.md`, `quantum-M4.md`).

---

## Headline

Two of the band-2 levers are **band-index-agnostic theorems** (W1, W4). The
gatekeeper W2 splits: **classically it lifts verbatim** — the top coefficient of
a genuine band-`k` pair with `a_k ≠ 0` is forced to be a scalar times a perfect
`k`-th power, `a_k = c·h^k` — but **quantumly the naive lift is false for `k ≥
3`**: the ladder-`(2k−1)` gatekeeper equation admits solutions in which `a_k` is
*not* a shifted `k`-th power `c·h(E)h(E+1)⋯(E+k−1)`. The obstruction is exact and
structural: over the root-necklace ring `Z[σ]` the operators `S_k = 1+σ+⋯+σ^{k−1}`
and `S_{k−1}` are coprime, so for `k = 2` (`S_{k−1} = S_1 = 1`) the cofactor is
forced effective and `a_2` must be a shifted square, but for `k ≥ 3` the nontrivial
`S_{k−1}` opens a gap. The witness `u = E(E+3)`, `a_3 = E(E+2)(E+4)` solves the
band-3 quantum gatekeeper with `a_3` a product of three **non-consecutive** roots.
This is the precise point at which a naive quantum band-`k` induction fails, and
it is a band-3 phenomenon (absent at band 2).

---

## W1 — top Wronskian ⇒ proportionality (THEOREM, any k, both faces)

**Theorem W1.** For any `k`, the extreme ladder coefficient is a pure Wronskian:

- Classical: `C_{2k} = k(a_k b_k′ − a_k′ b_k)`. Hence if `a_k ≠ 0`, then
  `(b_k/a_k)′ = 0`, so `b_k = λ a_k`, `λ ∈ C`.
- Quantum: `Q_{2k} = b_k^{[k]} a_k − a_k^{[k]} b_k`. Hence, if `a_k ≠ 0`,
  `b_k/a_k` is `k`-periodic as a rational function, so (char 0) constant:
  `b_k = λ a_k`.

*Proof.* Only the term `(i,l) = (k,k)` has `i+l = 2k` with both indices in
`[−k,k]`; the displayed one-term brackets give the formulas directly. Vanishing
of a Wronskian (classical) or of `b_k^{[k]}a_k − a_k^{[k]}b_k` (quantum) is
exactly the statement that `b_k/a_k` has zero derivative / is `k`-periodic; a
rational function with an additive period has no poles it cannot repeat
infinitely, hence is constant. ∎

This is the mirror-symmetric bottom fact too: if `a_{−k}≠0`, then
`C_{−2k} = −k(a_{−k}b_{−k}′ − a_{−k}′b_{−k})` and `Q_{−2k}` likewise give
`b_{−k} = μ a_{−k}` (used in W6). If `a_{−k}=0`, the extreme bottom equation is
vacuous and belongs to a separate vanishing branch.
The band-2 case is W1 as used in `M4_proof_memo.md` (`f637b1a`) `C_4`, and
`quantum-M4.md` (`84978b9`) `Q_4`. **Verified k = 2, 3, 4** (`verify §1`).

---

## W2 — the gatekeeper (the load-bearing lift; CLASSICAL verbatim, QUANTUM breaks at k≥3)

Gauge the top by `W1` (`G ↦ G − λF`, so `b_k = 0`). The next ladder down is the
homogeneous gatekeeper controlling every level's collapse.

### W2 classical (THEOREM, any k)

**Theorem W2c.** In the gauge `b_k = 0`, with `u := b_{k−1}`,
`C_{2k−1} = k a_k u′ − (k−1) a_k′ u`. For `a_k ≠ 0`, the equation `C_{2k−1} = 0`
has a **nonzero polynomial** solution `u` **iff** `a_k = c·h^k` for some `h ∈
C[τ]`, `c ∈ C^*`; and then `u = e·h^{k−1}`, with first integral `u^k = c′
a_k^{k−1}`.

*Proof.* `k a_k u′ = (k−1) a_k′ u` gives `(u^k / a_k^{k−1})′ = 0`, so `u^k = c′
a_k^{k−1}`. In the UFD `C[τ]`, comparing the exponent `e` of any irreducible in
`a_k` and `f` in `u` gives `k f = (k−1) e`; since `gcd(k, k−1) = 1`, `k | e`, so
`a_k = c·h^k` with `h` the product of the irreducibles to their `e/k`. Conversely
`a_k = c h^k` makes `u = e h^{k−1}` a solution. ∎

This is exactly the M4/J2 square lemma at `k = 2` (`a_2 = c h^2`, `2 a_2 u′ − a_2′
u = 0`; `M4_proof_memo.md` `f637b1a` §4). The generalization is the **`k`-th
power classification**; the negative direction (a non-`k`-th-power `a_k` kills
`u`, hence forces `b_{k−1} = 0`) is what makes it a gatekeeper. This one
homogeneous equation does **not** by itself force the lower coefficients to
vanish; any further descent must be proved rung by rung. **Verified k = 2, 3, 4**, including a machine witness that `a_k = τ`
(not a `k`-th power) admits no nonzero `u` (`verify §2`).

### W2 quantum (the break: THEOREM at k=2, FALSE naive lift at k≥3)

**Theorem W2q (structure).** In the gauge `b_k = 0`, with `u := b_{k−1}`,
`Q_{2k−1} = u^{[k]} a_k − a_k^{[k−1]} u`. Every nonzero solution satisfies the
exact **telescoping first integral**
```
  Π_{j=0}^{k−1} u^{[j]}  =  C · Π_{j=0}^{k−2} a_k^{[j]}       (C ∈ C^*),
```
obtained because `R := (Π u^{[j]})/(Π a_k^{[j]})` has `R^{[1]}/R =
(u^{[k]}/u)(a_k/a_k^{[k−1]}) = 1` on the equation locus, hence `R` is constant.

**Sufficiency (any k).** If `a_k = c·Π_{j=0}^{k−1} h^{[j]}` is a **shifted `k`-th
power**, then `u = e·Π_{j=0}^{k−2} h^{[j]}` solves `Q_{2k−1} = 0`.

**The break (k ≥ 3).** For `k ≥ 3` the shifted-`k`-th-power class is a *proper
subset* of the solution set: there are nonzero `u` with `a_k` **not** a shifted
`k`-th power. Explicit band-3 witness:
```
  u = E(E+3),   a_3 = E(E+2)(E+4):   u^{[3]} a_3 = a_3^{[2]} u   holds,
```
yet `a_3 = E(E+2)(E+4)` has three **non-consecutive** roots `{0, −2, −4}`, so it
is **not** `c·h(E)h(E+1)h(E+2)` for any `h, c` (that would need consecutive roots
`{r, r−1, r−2}`). The invariant `Π_{j=0}^{2} u^{[j]} = Π_{j=0}^{1} a_3^{[j]}`
still holds (`C = 1`).

**Why (the exact mechanism).** Pass to the root-necklace divisor. Per residue
class mod 1, write `δ(p) ∈ Z_{≥0}[σ]` for the multiset of roots (`σ = T^{−1}`,
`δ(p^{[r]}) = σ^r δ(p)`). The equation is equivalent to
```
  S_k(σ) δ(u) = S_{k−1}(σ) δ(a_k),      S_r(σ) := 1 + σ + ⋯ + σ^{r−1}.
```
Since `gcd(S_k, S_{k−1}) = 1` in `Q[σ]` (their roots are distinct roots of unity;
`gcd(σ^k−1, σ^{k−1}−1) = σ−1`, cancelled by the common `S`), a solution forces
`δ(u) = S_{k−1} g` and `δ(a_k) = S_k g` for a cofactor `g ∈ Z[σ]`. A nonzero `u`
exists **iff** both `S_{k−1} g` and `S_k g` are effective (non-negative). For
`k = 2`, `S_{k−1} = S_1 = 1`, so `g = δ(u)` is forced effective, forcing `a_2 =
S_2 g`-shape, i.e. `a_2 = c h(E)h(E+1)` a shifted square (this **is** Lemma J2q,
`quantum-M4.md` `84978b9`). For `k ≥ 3`, `S_{k−1}` is nontrivial and a
**non-effective** `g` can have both `S_{k−1} g` and `S_k g` effective. In fact the
same explicit cofactor works for every `k ≥ 3`:
```
  g = 1 − σ + σ²,
  S_{k−1}g = 1 + σ² + σ³ + ⋯ + σ^{k−2} + σ^k,
  S_k g     = 1 + σ² + σ³ + ⋯ + σ^{k−1} + σ^{k+1}.
```
(An indicated sum is empty when its lower limit exceeds its upper limit.) Thus
both products are effective while `g` is not. Taking `δ(u)=S_{k−1}g` and
`δ(a_k)=S_k g` in one root coset, and defining the monic polynomials by
`δ(∏_j(E+j)^{c_j})=∑_j c_jσ^j`, gives an exotic gatekeeper solution for every
`k ≥ 3`. The displayed formulas follow by multiplying
`(1+⋯+σ^{r−1})(1−σ+σ²)` and collecting coefficients; the necklace identity is
then `S_kδ(u)=S_{k−1}δ(a_k)` exactly. At `k=3` this recovers
`δ(u)=1+σ³`, `δ(a_3)=1+σ²+σ⁴` and the witness above.

**Consequences for the program.** The clean quantum band-2 gatekeeper "either
`b_1 = 0` or `a_2` shifted-square" (`quantum-M4.md`, `quantum-mirror.md`) does
**not** transcribe to band 3: at `k = 3` the gatekeeper equation `Q_5 = 0` can be
satisfied with `b_2 = u ≠ 0` and `a_3` outside the shifted-cube class. Since `u =
b_{k−1}` sits at a **positive** ladder level (no membership constraint), this
branch is fully admissible at the gatekeeper stage — it is not disposed of by
membership there, as the band-2 argument's kill was. **This is the sharpest
concrete obstruction a naive quantum band-3 induction must confront**, and it is
genuinely a band-3 (not band-2) phenomenon. Whether the exotic branch survives
the *full* band-3 cascade + membership is the open question it exposes (see
"Open branches"). **Verified k = 2, 3, 4** for the equation form, sufficiency,
the invariant, the exotic witness, the non-shifted-cube certificate, the
necklace identities, and the `S_k/S_{k−1}` coprimality (`verify §3`).

### W2 unifying picture (the confluent limit)

Classical W2 is the **confluent (`σ → 1`) degeneration** of quantum W2: as all
roots merge, `S_r(1) = r`, the necklace condition collapses to the exponent
condition `k f = (k−1) e`, and the effective-cofactor gap closes (there is only
one point, so effectivity is automatic). Differentiation is shift-invariant,
which is exactly why the classical face sees only multiplicities and never the
additive spread that opens the quantum `k ≥ 3` exotica.

---

## W4 — the m=0 moment (THEOREM, any k, both faces)

**Theorem W4c (classical).** For any `k`,
```
  C_0 = Ψ′,      Ψ := Σ_{i=1}^{k} i (a_i b_{−i} − a_{−i} b_i),
```
an **exact total `τ`-derivative**, with **no proportionalities assumed**. Hence
`{G,F} = 1` forces `Ψ = τ + const`, and membership (`τ^i | a_{−i}, b_{−i}`)
makes every product vanish at `τ = 0`, so `Ψ(0) = 0` and `Ψ = τ` exactly.

*Proof.* At `m = 0` only `l = −i` survives; the summand is `i a_i b_{−i}′ − (−i)
a_i′ b_{−i} = i(a_i b_{−i})′`, and pairing `i` with `−i` gives `Σ_{i≥1}
i(a_i b_{−i} − a_{−i} b_i)`, whose derivative is `C_0`. ∎

This is the general mechanism the band-2 memos observed as "telescoping": the
moment equation is a **pure divergence** at the level of the raw ladder system;
the cascades are needed only to *solve* it, not to *close* it. At `k = 2` this is
`C_0 = 2(a_2 b_{−2})′ + (a_1 b_{−1})′ − (a_{−1} b_1)′ − 2(a_{−2} b_2)′`
(`M4_proof_memo.md` `f637b1a`).

**Theorem W4q (quantum).** For any `k`, `Q_0 = (T−1)G` with the explicit staggered
potential
```
  G = Σ_{i=1}^{k} Σ_{r=1}^{i} ( a_i^{[−r]} b_{−i}^{[i−r]} − a_{−i}^{[i−r]} b_i^{[−r]} ).
```
Hence `[D,X] = 1` forces `(T−1)G = 1`, so `G = E + c_0`; membership gives `G(0) =
0`, so `G = E`.

*Proof.* For each `i`, `(T−1) Σ_{r=1}^{i} a_i^{[−r]} b_{−i}^{[i−r]}` telescopes
(`f(r) := a_i^{[−r]} b_{−i}^{[i−r]}`, term `r` contributes `f(r−1) − f(r)`) to
`f(0) − f(i) = a_i b_{−i}^{[i]} − a_i^{[−i]} b_{−i}`, the index-`i` term of `Q_0`;
the second sum gives the index-`(−i)` term. ∎

At `k = 2`, in the gauge of `quantum-mirror.md` (`84978b9`, `a_2 = 1`, `b_1 = κ`),
this `G` collapses exactly to that memo's central telescoping potential `G =
w^{[1]} + w + p^{[−1]} v − κ u` (§2 there) — machine-confirmed. **Verified k = 2,
3, 4** (`verify §4`).

---

## W5 — the lattice conjecture (FORMULATION + band-2 instance as a framework theorem)

The strategic heart. In band 2 every resistant branch died of **integer
arithmetic on degree vectors** (`classical-Astar.md`, `84978b9`). We formulate the
general mechanism and prove the band-2 instance inside the framework.

### Formulation

After the W1/W2 top cascade and the W4 moment normalization, the residual
system is `C_{−1}, …, C_{−(2k)}` (classical) / `Q_{−1}, …, Q_{−2k}` (quantum) in
the surviving coefficient polynomials, whose **degrees** form a vector `d =
(d_1, …, d_n) ∈ Z_{≥0}^n`. Reading each residual equation at its generic top
degree yields, **away from leading-term cancellation loci**, a `Z`-affine
relation
```
  Σ_i α_{ji} d_i = β_j.
```
The coefficients `α` come from the shift/derivative degree bookkeeping (each
factor contributes its degree; a derivative or a `[·]`-shift lowers/raises by a
fixed integer), and the inhomogeneous `β_j` is nonzero for exactly the relations
fed by the **W4 moment's unit** (the lone `+1`). Collect these into `L d = β`,
`L ∈ Z^{r×n}`, `β ∈ Z^r`.

> **Lattice conjecture (band k).** Outside the tame strata, the affine system
> `L d = β` has **no solution** `d ∈ Z_{≥0}^n`. The obstruction is a
> congruence: some integer row combination `γ^⊤ L d = γ^⊤ β` reduces to `g·(lin)
> = m` with `g := gcd(γ^⊤ L) ∤ m` (equivalently, `β ∉ L(Z^n)` in the relevant
> quotient). The tame strata are exactly the cancellation loci excluded above,
> where a leading term drops and the corresponding coefficient is forced constant
> — recovering the tame normal form.

The leading-term cancellation loci (where an `α`-relation fails because the top
coefficient vanishes) are handled by descent: each is lower-dimensional and
either forces a coefficient constant (tame) or reduces to a smaller residual
system. This is the "split cancellation loci explicitly" discipline made
structural.

### Band-2 instance (THEOREM of the framework; = `classical-Astar.md` A\*)

For the classical square sector, degrees `d = (P, V, W) = (deg p, deg v, deg w)`.
The two first integrals `Φ`, `I₂` (`classical-Astar.md` Lemma 2.1) give, at
leading order,
```
  (Φ = 0):   2V = P + W,        (I₂ = 0):   V + W = 2P + 1,
```
i.e. `L d = β` with `L = [[−1,2,−1],[−2,1,1]]`, `β = (0,1)^⊤`. Row₁ + Row₂ =
`(−3, 3, 0)`, so `3(V − P) = 1`: `g = gcd(3,3,0) = 3 ∤ 1`. **Infeasible over `Z`
(a fortiori over `Z_{≥0}`)** — branch A\* is empty. This is the exact **mod-3
kill** of `classical-Astar.md`, now realized as the generic congruence
obstruction of the lattice framework. **Verified** (`verify §5`).

### Toward band 3

The `+1` in `V + W = 2P + 1` is the W4 moment unit (the `∫₀ᵗ` from the second
integral); the modulus `3` is the gcd of the reduced row. At band 2, `3 = k+1 =
2k−1`, so the band-2 instance does not distinguish which quantity controls the
general modulus. Resolving that requires the full band-3 residual first-integral
set (the sibling M5-band-3 derivation) — this memo supplies the **proven balance
inputs** for it (W1/W2/W4/W6 at band 3) and the framework into which they feed,
but does **not** claim the band-3 modulus. Status: **formulation + band-2
theorem proved; general-`k` congruence a conjecture; band-3 modulus open.**

---

## W3 — conditional trailing-coefficient identity (any k)

**Conditional identity W3.** The trailing coefficient `a_{−k}` enters the
equation `C_{−(k−1)}` **only** through
```
  −k a_{−k} b_1′ − a_{−k}′ b_1  =  −(k−1) a_{−k} b_1′  −  (a_{−k} b_1)′.
```
This decomposition is an arbitrary-`k` identity. If an **independent descent
argument** has already proved `b_1′=0`, its residue vanishes and the
`a_{−k}`-dependent part is the exact derivative `−(a_{−k}b_1)′`.

Nothing in W2 alone supplies that hypothesis: W2 controls the homogeneous
freedom `b_{k−1}` at the top wall, not the intervening coefficients
`b_{k−2},…,b_1`. Nor does the displayed identity show that all the other terms
of `C_{−(k−1)}` form an exact derivative. Therefore no arbitrary-`k` W3 first
integral, no free-`a_{−k}` theorem, and no whole-cascade conclusion is claimed.
At `k=2`, the separately proved gauged parametrization in
`classical-Astar.md` yields its genuine W3 integral; extending that extra input
to larger `k` remains open. **The conditional decomposition is verified at
`k=2,3,4`** (`verify §6`).

---

## W6 — Lemma R is the mirror of W2 (audit correction; balance present in BOTH faces)

**Theorem W6.** On the nonzero bottom branch `a_{−k}≠0`, the bottom-adjacent
equation `C_{−(2k−1)}` / `Q_{−(2k−1)}` is the mirror image of the W2 gatekeeper.
With the bottom proportionality `b_{−k} = μ a_{−k}` (W1 mirror) and
`φ := μ a_{−(k−1)} − b_{−(k−1)}`, `s := a_{−k}`:

- **Quantum:** `Q_{−(2k−1)} = s^{[−(k−1)]} φ − s φ^{[−k]}`. The `E^{S+Φ−1}` leading
  coefficient is `(kΦ − (k−1)S)·lc(s)·lc(φ)`, so vanishing forces the **staggered
  rigidity** `(k−1)·deg(a_{−k}) = k·deg(φ)`; since `gcd(k, k−1) = 1`, `k | deg
  a_{−k}`. Beyond leading order, this is the necklace equation of W2q — with the
  **same `k ≥ 3` exotica**.
- **Classical:** `C_{−(2k−1)} = k s φ′ − (k−1) s′ φ`. Its `E^{S+Φ−1}` leading
  coefficient is **also** `(kΦ − (k−1)S)`, and the full equation integrates to the
  **mirror `k`-th power** first integral `φ^k = c·s^{k−1}` — exactly W2c reflected.

**Audit correction.** The band-2 memo `quantum-mirror.md` (`84978b9`) §5 frames
the `(2Φ − S)` rigidity as **quantum-only** ("differentiation is shift-invariant,
so classical `C_{−3}` does not force `deg a_{−2} = 2 deg φ`"). That is too strong.
The reduced classical `C_{−3} = 2 s φ′ − s′ φ` has leading coefficient `(2Φ − S)`
and integrates to `φ² = c s`, so it **does** force `S = 2Φ`. The band-2 A\*
closure simply *routed around* `C_{−3}` (using the `Φ, I₂` integrals), which is a
choice of proof, not an absence of the balance. `quantum-mirror.md`'s **theorem**
(branch emptiness) stands; only its "quantum-only" *characterization of the
mechanism* is corrected here. The genuine classical/quantum distinction at the
bottom is identical to W2: classical = clean `k`-th power `φ^k = c s^{k−1}`;
quantum = staggered necklace with `k ≥ 3` exotica. **Verified k = 2, 3, 4**
(`verify §7`), including the classical integral `φ^k = c s^{k−1}`.

Whether a Lemma-R-style *staggered* rigidity gives the quantum face an extra lever
over the classical one is therefore a **red herring at the level of the degree
balance** — the balance is shared. The real asymmetry is the solution *class*
(power vs necklace), and it cuts the *same way* at top (W2) and bottom (W6).

---

## Rigor ledger

**THEOREM (arbitrary k, structural, machine-checked at k = 2, 3, 4):**
- W1 top/bottom (shifted-)Wronskian ⇒ proportionality (both faces).
- W2c classical gatekeeper: nonzero `u` ⇔ `a_k = c h^k` (`k`-th power), UFD
  exponent argument.
- W2q quantum gatekeeper: the equation form, the telescoping first integral, and
  shifted-`k`-th-power **sufficiency**; the necklace criterion `S_k δ(u) = S_{k−1}
  δ(a_k)` with the `S_{k}/S_{k−1}` coprimality (all `k`).
- W4 moment: `C_0 = Ψ′` and `Q_0 = (T−1)G` with explicit potentials (both faces).
- W3 algebraic decomposition of the `a_{−k}`-part of `C_{−(k−1)}`; its
  exact-derivative consequence is conditional on an independently proved
  `b_1′=0` and is not a general first-integral theorem.
- W6 bottom = W2 mirror: staggered balance `(k−1)S = kΦ` in both faces; classical
  integral `φ^k = c s^{k−1}`.

**VERIFIED at k = 2, 3, 4 as distinguishing computations:**
- The quantum-W2 **break** at band 3: `u = E(E+3)`, `a_3 = E(E+2)(E+4)` and its
  non-shifted-cube certificate. The written product formulas for
  `g=1−σ+σ²` prove the exotic family for every `k≥3`; finite checks merely
  corroborate the first cases.

**CONJECTURE:**
- W5 lattice conjecture at general `k` (infeasibility of `L d = β` over `Z_{≥0}`
  outside tame strata, via a congruence `g ∤ m`). The band-2 A\* instance is
  proved inside the framework; the general modulus and the band-3 modulus are open.

**AUDIT CORRECTION (proved here):**
- The bottom staggered degree balance `(k−1)S = kΦ` is present in the **classical**
  face too (via `C_{−(2k−1)} = k s φ′ − (k−1) s′ φ`, integral `φ^k = c s^{k−1}`),
  not quantum-only as `quantum-mirror.md` §5 states. That memo's emptiness theorem
  is unaffected.

**NOT claimed:** any band-`k` classification, any full band-3 theorem, JC2, or
DC1. The membership disposition of the quantum-W2 exotic branch, the band-3
first-integral set, and the general lattice modulus are all open.

---

## What genuinely differs at band 3 (vs band 2)

1. **Quantum W2 acquires exotic solutions** (`k ≥ 3`): the gatekeeper no longer
   forces `a_k` into the shifted-`k`-th-power class. Band 2 is the degenerate
   case `S_{k−1} = 1` where no gap exists. *(New, sharp, machine-witnessed.)*
2. **W6 is unmasked as the W2 mirror**, and the "quantum-only rigidity" reading of
   the band-2 corpus is corrected: the staggered degree balance is shared by both
   faces. *(Audit correction.)*
3. **W1 and W4 are untouched** — band-index-agnostic theorems with clean
   arbitrary-`k` proofs and explicit potentials.
4. **The lattice mechanism (W5) survives as a formulation**, with the band-2
   congruence kill re-derived as its generic instance; the band-3 modulus is the
   next quantitative target and depends on the full band-3 residual system.

---

## Verification

```sh
uv run --with sympy python research/band3/verify_bandk_weapons.py
```

Exact SymPy; 94 `PASS` lines at `k = 2, 3, 4`; ends `ALL BAND-K WEAPON CHECKS
PASSED`. The checks establish the displayed identities and the distinguishing
witnesses; the arbitrary-`k` completeness of the THEOREM items is the written
structural argument (leading-coefficient bookkeeping, UFD/necklace divisibility,
telescoping), not the finite-`k` computation.
