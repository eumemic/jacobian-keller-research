# The classical band-3 cascade

**INDEPENDENTLY DERIVED — MACHINE-CHECKED IDENTITIES — NOT PEER REVIEWED — BAND-SCOPED**

This memo derives the complete classical band-3 coefficient system and pushes the
band-2 arsenal down the cascade as far as it forces, mapping exactly where the
band-2 weapons break. Every displayed algebraic identity is machine-verified by
`verify_classical_band3.py` (same directory); a successful run ends
`ALL CLASSICAL BAND3 CHECKS PASSED` (66 checks). The degree-balance and
membership arguments are written proofs; the bounded-degree kernel sweeps are
regression corroboration only.

Conventions are frozen and identical to the band-2 program
(`research/band2-classical-proved/M4_proof_memo.md`, in the full-classical
package at commit `f637b1a`; square-sector closure `research/band2-square-sector/classical-Astar.md`
at `91a053a`; quantum assembly `84978b9`). Over `C`, put `τ = xξ`,
`{G,F} = G_ξ F_x − G_x F_ξ`, so `{ξ,x} = 1`.

## 0. Setup, support, membership

Write, in `C[x,ξ]` (not merely the Laurent ring),
```
F = Σ_{k=-3}^{3} x^k a_k(τ),      G = Σ_{l=-3}^{3} x^l b_l(τ),
```
with all coefficients in `C[τ]`. Genuine polynomial membership is
```
τ | a_{-1}, b_{-1};   τ² | a_{-2}, b_{-2};   τ³ | a_{-3}, b_{-3},
```
because `x^k τ^j = x^{k+j} ξ^j`. "Band 3" means support ⊆ `{-3,…,3}`; an extreme
coefficient may vanish. Primes are `d/dτ`.

For one-term pieces `{x^ℓ b(τ), x^k a(τ)} = x^{k+ℓ}(k a b' − ℓ a' b)`, so the
coefficient of `x^m` in `{G,F}` is
```
C_m = Σ_{k+ℓ=m} (k a_k b_ℓ' − ℓ a_k' b_ℓ),      C_m = δ_{m0},   m ∈ [-6,6].
```

Allowed normalizations (as in band 2): (i) `G → G − λF`; (ii) additive
constants; (iii) diagonal symplectic scaling `(x,ξ) → (ρx, ρ^{-1}ξ)`, which
sends `a_k → ρ^k a_k`, so it multiplies `a_3` by `ρ³` and `a_2` by `ρ²`.
Orientation operations (pair exchange, reflection) are used exactly as in
`classical-Astar.md`; the reflection `R f(x,ξ)=f(ξ,x)` acts by
`(Rf)_k = τ^{-k} a_{-k}` and preserves band 3.

## 1. The thirteen equations

Machine-checked identically against a direct two-variable Poisson computation
with fully generic coefficient functions (verifier §1):

```
C_6  = 3 a3 b3' − 3 a3' b3.

C_5  = 3 a3 b2' − 2 a3' b2 + 2 a2 b3' − 3 a2' b3.

C_4  = 3 a3 b1' −   a3' b1 + a1 b3' − 3 a1' b3 + 2 a2 b2' − 2 a2' b2.

C_3  = 3 a3 b0'          − 3 a0' b3 + 2 a2 b1' −   a2' b1 + a1 b2' − 2 a1' b2.

C_2  = 3 a3 b_{-1}' +   a3' b_{-1} − a_{-1} b3' − 3 a_{-1}' b3
       + 2 a2 b0' − 2 a0' b2 + a1 b1' − a1' b1.

C_1  = 3 a3 b_{-2}' + 2 a3' b_{-2} − 2 a_{-2} b3' − 3 a_{-2}' b3
       + 2 a2 b_{-1}' + a2' b_{-1} − a_{-1} b2' − 2 a_{-1}' b2
       + a1 b0' − a0' b1.

C_0  = 3 a3 b_{-3}' + 3 a3' b_{-3} − 3 a_{-3} b3' − 3 a_{-3}' b3
       + 2 a2 b_{-2}' + 2 a2' b_{-2} − 2 a_{-2} b2' − 2 a_{-2}' b2
       + a1 b_{-1}' + a1' b_{-1} − a_{-1} b1' − a_{-1}' b1.

C_{-1} = −a_{-1} b0' + a0' b_{-1} + a1 b_{-2}' + 2 a1' b_{-2}
         − 2 a_{-2} b1' − a_{-2}' b1 + 2 a2 b_{-3}' + 3 a2' b_{-3}
         − 3 a_{-3} b2' − 2 a_{-3}' b2.

C_{-2} = −a_{-1} b_{-1}' + a_{-1}' b_{-1} + a1 b_{-3}' − 2 a_{-2} b0'
         + 2 a0' b_{-2} + 3 a1' b_{-3} − b1 a_{-3}' − 3 a_{-3} b1'.

C_{-3} = −a_{-1} b_{-2}' + 2 a_{-1}' b_{-2} − 2 a_{-2} b_{-1}' + a_{-2}' b_{-1}
         − 3 a_{-3} b0' + 3 a0' b_{-3}.

C_{-4} = −a_{-1} b_{-3}' + a_{-3}' b_{-1} − 2 a_{-2} b_{-2}' + 2 a_{-2}' b_{-2}
         − 3 a_{-3} b_{-1}' + 3 a_{-1}' b_{-3}.

C_{-5} = 2 a_{-3}' b_{-2} − 3 a_{-3} b_{-2}' − 2 a_{-2} b_{-3}' + 3 a_{-2}' b_{-3}.

C_{-6} = 3 a_{-3}' b_{-3} − 3 a_{-3} b_{-3}'.
```

**Reflection antisymmetry.** Substituting `k → −k, ℓ → −ℓ` gives
`C_{-m}(a,b) = −C_m(ã, b̃)` with `ã_k = a_{-k}`, `b̃_ℓ = b_{-ℓ}`. Thus the
negative half of the cascade is the reflected image of the positive half; every
top-side result has a bottom-side mirror.

## 2. Top and bottom Wronskians (W1)

Assume `a_3 ≠ 0`. Then `C_6 = 3(a_3 b_3' − a_3' b_3) = 0` gives `(b_3/a_3)' = 0`,
hence
```
b_3 = λ_3 a_3,   λ_3 ∈ C.                                        (W1)
```
By reflection, `a_{-3} ≠ 0` and `C_{-6} = 0` give `b_{-3} = μ_3 a_{-3}`.
The single gauge move `G → G − λ_3 F` sets `b_3 = 0`. There is **one** gauge
parameter, and it kills the level-3 multiple only; the level-2 relation below is
then a *consequence*, not a second gauge (see §4).

## 3. THE WALL: C_5 and the cube classification (W2 → 2/3-power)

This is the first equation coupling two different weight levels. Substitute
`b_3 = λ_3 a_3`. A direct reduction (verifier §3) shows that with
`u_2 := b_2 − λ_3 a_2`,
```
C_5 = L[u_2],       L[u] := 3 a_3 u' − 2 a_3' u.                  (H5)
```
The gauge parameter drops out; `C_5 = 0 ⇔ L[u_2] = 0`. This is the exact
band-3 analogue of the band-2 wall `2 a_2 u' − a_2' u = 0` (M4, "H3"), with the
operator advancing from the 1/2-power to the **2/3-power** form.

**Lemma 3.1 (cube integrating factor).** For `a_3 ≠ 0`,
```
(u³ / a_3²)' = u² · L[u] / a_3³        (machine-checked identity),
```
so `L[u] = 0 ⇔ (u³/a_3²)' = 0 ⇔ u³ = c a_3²` for some `c ∈ C^*` (or `u=0`).

**Proposition 3.2 (exact homogeneous-freedom criterion).** Over `C`, a *nonzero*
polynomial solution of `L[u] = 0` exists **iff `a_3` is a scalar cube**,
`a_3 = c h³` with `h ∈ C[τ]`, `c ∈ C^*`. In that case the full polynomial
solution space is the line `{u = e h² : e ∈ C}`.

*Proof.* Nonzero `u` gives `u³ = c a_3²`, `c ≠ 0`. Factor over `C`:
`a_3 = γ ∏_i (τ − r_i)^{m_i}`. Then `u³` and `a_3²` have the same roots with
`3·mult_{r_i}(u) = 2 m_i`. Since `gcd(2,3) = 1`, `3 ∣ m_i` for every `i`;
writing `m_i = 3 p_i`, `a_3 = γ (∏_i(τ−r_i)^{p_i})³ = c h³` with `h = ∏_i(τ−r_i)^{p_i}`,
`c = γ`. Conversely `a_3 = c h³` gives `L[e h²] = 3 c h³·2 e h h' − 2·3 c h² h'·e h² = 0`
(machine-checked), and `u³ = c' a_3² = c'' h⁶` forces `u ∈ C·h²`. ∎

The general pattern (verifier §4): the level-`ℓ` descent operator is
```
L_ℓ[u] := 3 a_3 u' − ℓ a_3' u,      (u³/a_3^ℓ)' = u² L_ℓ[u] / a_3^{ℓ+1},
```
so `L_ℓ[u] = 0 ⇔ u³ = c a_3^ℓ`. The polynomial-kernel criterion is
sign-asymmetric:
- `ℓ > 0`: nonzero solution iff `a_3^ℓ` is a scalar cube — iff `a_3` is a scalar
  cube (when `gcd(ℓ,3)=1`, i.e. `ℓ ∈ {1,2}`), giving `u = e h^ℓ`; or for any
  `a_3` (when `3 ∣ ℓ`), giving `u = e·a_3^{ℓ/3}`.
- `ℓ = 0`: `u = ` const always.
- `ℓ < 0`: `u³ a_3^{|ℓ|} = ` const, so `deg a_3 ≥ 1` forces `u = 0`; only
  constant `a_3` admits a nonzero (constant) solution.

Over `C` a nonzero constant is a cube, so **"a_3 not a cube" forces
`deg a_3 ≥ 1`.** The wall `L = L_2`; the `C_4` operator is `M := L_1`; the
negative rungs of §6 use `L_{-1}, L_{-2}` under `deg a_3 ≥ 1`.

## 4. C_4: the three-level coupling and the mixed (a_3, a_2) sector

Substitute `b_3=λ_3 a_3`, `b_2=λ_3 a_2+u_2`, `b_1=λ_3 a_1+u_1`. Then
(verifier §5), gauge-invariantly,
```
C_4 = M[u_1] + 2 W(a_2, u_2),   M[u] = 3 a_3 u' − a_3' u,
      W(f,g) = f g' − f' g.                                      (C4*)
```
So `C_4 = 0` is an inhomogeneous first-order ODE for `u_1` driven by the
Wronskian of `a_2` with the wall datum `u_2`. This is the first place the
level-2 coefficient `a_2` — a *free* function of `F`, not a top coefficient —
enters. Band 2 had no analogue: there the top `a_2 = h²` was rigid, whereas
here the rigidity sits at level 3 (`a_3 = c h³`) and `a_2` is only coupled in.

**The mixed coupling.** Work in the cubic sector `a_3 = c h³`, gauge `b_3 = 0`
(so `λ_3 = 0`, `u_2 = b_2`, `u_1 = b_1`). By §3, `b_2 = e h²` (`e ∈ C` free).
Using `2 a_2 h' − a_2' h = −h³ (a_2/h²)'` and `M[hφ] = 3 a_3 h φ'`, equation
`(C4*)` integrates in closed form (verifier §9):
```
b_1 = κ_1 h + (2e / 3c) · a_2 / h,      κ_1 ∈ C.                  (4.1)
```
**Consequence (new phenomenon).** For `e ≠ 0`, polynomiality of `b_1` forces
```
h | a_2.
```
The verifier exhibits `h = τ`, `a_2 = 1`, `e ≠ 0`: the forced `b_1` has a genuine
pole, so no polynomial pair exists there. When `e = 0`, (4.1) collapses to the
clean band-2-style `b_1 = κ_1 h`, and the mixed coupling is silent.

## 5. The m = 0 moment (W4 at arbitrary band)

Because index pairs `(k, −k)` make the summand a total derivative,
```
C_0 = d/dτ [ M(τ) ],   M(τ) := Σ_{k≥1} k (a_k b_{-k} − a_{-k} b_k)
    = 3(a_3 b_{-3} − a_{-3} b_3) + 2(a_2 b_{-2} − a_{-2} b_2) + (a_1 b_{-1} − a_{-1} b_1)
```
(verifier §8). Hence `C_0 = 1 ⇔ M' = 1 ⇔ M(τ) = τ + const`. Every product in `M`
vanishes at `τ = 0` by membership (`τ^j ∣ a_{-j}, b_{-j}`), so `const = 0` and
```
M(τ) = τ.                                                        (W4)
```
This is the exact band-3 central identity, generalizing band-2's
`(2h²w + hpv − κh a_{-1})' = 1`; it holds at every band `K` with
`M = Σ_{k=1}^{K} k(a_k b_{-k} − a_{-k} b_k)`.

## 6. Theorem A: the non-cube top sector is empty

> **Theorem A.** Let `F, G ∈ C[x,ξ]` be a band-3 Keller pair with `a_3 ≠ 0`.
> If `a_3` is **not** a scalar cube (equivalently, over `C`, `deg a_3 ≥ 1` and
> `a_3 ≠ c h³`), then no such pair exists: the sector is empty.

*Proof.* Orient so `a_3 ≠ 0` and gauge `b_3 = 0` (§2). All reductions below are
machine-checked identities (verifier §6); the conclusions use §3–§4 and
membership.

1. `C_5 = L[b_2] = 0`. Non-cube ⇒ (Prop. 3.2) `b_2 = 0`.
2. `C_4 = M[b_1] = 0` (the `2W(a_2,b_2)` term vanishes since `b_2=0`).
   Non-cube ⇒ `b_1 = 0` (`M = L_1`, `gcd(1,3)=1`).
3. `C_3 = 3 a_3 b_0' = 0` ⇒ `b_0 = ` const (`a_3 ≠ 0`).
4. `C_2 = L_{-1}[b_{-1}] = 3 a_3 b_{-1}' + a_3' b_{-1} = 0`
   ⇒ `(a_3 b_{-1}³)' = 0` ⇒ `a_3 b_{-1}³ = ` const. Since `deg a_3 ≥ 1`,
   `b_{-1} = 0`.
5. `C_1 = L_{-2}[b_{-2}] = 3 a_3 b_{-2}' + 2 a_3' b_{-2} = 0`
   ⇒ `(a_3² b_{-2}³)' = 0` ⇒ `b_{-2} = 0` (same reason).
6. `C_0 = 3(a_3 b_{-3})' = 1` ⇒ `a_3 b_{-3} = τ/3 + ` const; membership
   `τ³ ∣ b_{-3}` kills the constant (`a_3 b_{-3}(0)=0`), giving `a_3 b_{-3} = τ/3`.
   But `τ³ ∣ b_{-3}` makes `a_3 b_{-3}` have `τ`-order `≥ 3`, contradicting the
   `τ`-order-`1` right side. (If `b_{-3}=0` the left side is `0 ≠ τ/3`.) ∎

This is the band-3 twin of the band-2 nonsquare-empty theorem (M4 §5). Note the
symmetric statement holds by reflection: `a_{-3} ≠ 0` non-cube is also empty.
The endgame at step 6 is a *strict* `τ`-order inequality — no leading-term
cancellation is possible — so the argument is cancellation-safe.

## 7. The cubic sector: reduction and residual system

By §6, a band-3 pair with `a_3 ≠ 0` has `a_3 = c h³` (`c ≠ 0`). Gauge `b_3 = 0`.
The wall gives `b_2 = e h²`; (4.1) gives `b_1 = κ_1 h + (2e/3c)a_2/h` (needing
`h ∣ a_2` when `e ≠ 0`). Two structural sub-splits organize the sector.

### 7.1 Constant vs nonconstant h

Diagonal scaling multiplies `a_3` by `ρ³`, hence `c` by `ρ³` and `h` by `ρ`
(up to redistributing the constant); it can normalize `c = 1` but cannot make a
nonconstant `h` constant. In band 2, the central identity forced `h ∣ τ` and
excluded nonconstant `h` in one line (`classical-Astar.md` §6). At band 3 that
one-liner **does not transfer**: the moment `M = τ` (W4) carries the `e`-term
`−(2e/3c) a_2 a_{-1}/h`, which is polynomial precisely under the already-forced
`h ∣ a_2`; no "`h ∣ τ`" divisibility is produced. Excluding nonconstant `h` is
therefore an **open branch** requiring a genuine negative-cascade argument
(candidate weapon: W5 leading-coefficient integer-lattice obstruction on the
moment `M = τ` together with `C_{-5}`, the bottom wall — not yet carried out).

### 7.2 e = 0 vs e ≠ 0

- **e = 0** (`b_2 = 0`): the top data is `b_2=0`, `b_1 = κ_1 h` — structurally
  identical to the band-2 square-sector opening (`b_2=0`, `b_1=κh`), but now
  sitting *above* an extra free level (`a_2`). The descent `C_3, C_2, C_1`
  integrates cleanly (below).
- **e ≠ 0** (`b_2 ≠ 0`): genuinely new. Forces `h ∣ a_2`, and the source term in
  `(C4*)` propagates the extra datum `e` through the entire descent. No band-3
  Keller pair with `e ≠ 0` is known; whether the sector is empty or rigid is
  **open** (the generic bounded-degree solver did not terminate within budget,
  so not even regression evidence is claimed here).

### 7.3 Explicit constant-h residual system (h = 1, a_3 = c)

Here `h | a_2` is automatic, `b_2 = e`, and the positive cascade integrates to
(all machine-checked, verifier §9):
```
b_1  = κ_1 + (2e/3c) a_2,
b_0  = (2e a_1 + κ_1 a_2)/(3c) − e a_2²/(9c²) + β,
b_{-1} = (1/3c)[ 2e a_0 + κ_1 a_1 − (2e/3c) a_1 a_2
                 − (κ_1/3c) a_2² + (4e/27c²) a_2³ ] + γ.
```
The central identity is the moment `M = τ`:
```
3c·b_{-3} − 2e·a_{-2} − κ_1 a_{-1} + a_1 b_{-1} + 2 a_2 b_{-2}
   − (2e/3c) a_2 a_{-1} = τ.
```
The residual system is then `C_{-1} = C_{-2} = C_{-3} = 0` together with
`C_{-6}` (bottom proportionality `b_{-3} = μ_3 a_{-3}`), `C_{-5}` (bottom wall),
`C_{-4}`, and the memberships `τ ∣ b_{-1}`, `τ² ∣ a_{-2}, b_{-2}`,
`τ³ ∣ a_{-3}, b_{-3}`. Unlike band 2, this system has **two** free
lower F-functions (`a_2` and `a_1`) surviving into the negative cascade; closing
it is exactly the additional work band 3 demands over band 2. The `e = 0` slice
of this system is the one that most closely mirrors the band-2 square sector.

## 8. First integrals of the negative cascade (W3 status)

The clean, provable first integrals are:

1. **Wronskian ratios** (§2): `b_3/a_3 = λ_3`, `b_{-3}/a_{-3} = μ_3` — first
   integrals of `C_6`, `C_{-6}`.
2. **The moment** `M = τ` (§5) — a two-sided first integral (W4), with the
   very-bottom coefficients `a_{-3}, b_{-3}` entering linearly.

For a band-2-style **W3 free-trailing-coefficient** integral (`Φ' = κ C_{-1}`
with `a_{-2}` free, `classical-Astar.md` Lemma 2.1), the honest band-3 finding is
an **obstruction**: writing
```
C_{-1} = d/dτ [ Σ_k k a_k b_{-1-k} ] + Σ_k a_k' b_{-1-k},
```
the residual `Σ_k a_k' b_{-1-k}` is not a total derivative of the raw
coefficients, so no gauge-free W3 integral of `C_{-1}` exists at band 3. The
band-2 W3 integral survives only after the full gauged parametrization, which at
band 3 is entangled with the mixed `(a_3, a_2)` coupling of §4. In the
constant-h `e = 0` slice (§7.3) the descent integrations *are* the first
integrals (`b_0, b_{-1}` as explicit antiderivatives, central identity `M = τ`);
producing the trailing-free analogue for `e ≠ 0` is part of the open §7.2
branch. Thus the W3 weapon is present in reduced form but does not generalize
verbatim.

## 9. Solver-facing case tree

```
band-3 pair {G,F}=1, supports ⊆ [-3,3], genuine membership
│
├─ a_3 = a_{-3} = 0  and lower extremes 0 ......... routes DOWN to band ≤ 2
│      (a_2=0 or a_{-2}=0 reductions handed to the band-2 theorem:
│       full classical band-2 = commit f637b1a; square sector = 91a053a)
│
├─ orient (pair exchange / reflection, classical-Astar §1) so a_3 ≠ 0
│  │
│  ├─ a_3 NOT a scalar cube (deg a_3 ≥ 1, a_3 ≠ c h³)
│  │        ⇒ EMPTY.  [Theorem A, §6 — PROVED at arbitrary degree]
│  │
│  └─ a_3 = c h³ (cube sector).  Gauge b_3=0; b_2=e h²; b_1 per (4.1).
│     │
│     ├─ e = 0  (b_2 = 0)
│     │   ├─ h constant  →  constant-h e=0 residual system (§7.3, e=0);
│     │   │                  mirrors band-2 square sector above a free a_2.
│     │   │                  [descent integrations PROVED; closure OPEN]
│     │   │                  Contains the tame shear family
│     │   │                  F = ξ + α x³ + β x² + γ x + δ, G = −x  (witness, §10).
│     │   └─ h nonconstant →  OPEN branch (§7.1): band-2 one-line
│     │                        exclusion does NOT transfer.
│     │
│     └─ e ≠ 0  (b_2 ≠ 0)
│         ⇒ forces h | a_2  (§4, new).  Genuinely new mixed sector.
│         [no known pair; empty-or-rigid is OPEN, §7.2]
```

Sectors that reduce to band 2 are cited above. The genuinely new sectors are the
**cubic sector `a_3 = c h³`** (all of it) and, within it, the **mixed `e ≠ 0`
coupling** and the **nonconstant-`h`** branch.

## 10. Witness

`F = ξ + α x³ + β x² + γ x + δ`, `G = −x` is a tame (triangular) polynomial
automorphism with `{G,F} = 1` (verifier §10). Its ladder is
`a_3=α, a_2=β, a_1=γ, a_0=δ, a_{-1}=τ`, `b_1=−1`; it satisfies all thirteen
`C_m = δ_{m0}` and lies in the constant-h, `e = 0` cubic sub-sector with
`κ_1 h = −1`, `h = α^{1/3}`. It anchors the case tree and exercises
`a_3, a_2, a_1, a_0 ≠ 0` simultaneously.

## 11. Claim disposition

**Proved at arbitrary degree (with machine-checked identities):**
- the complete thirteen-equation system and its match to the direct Poisson
  bracket (§1);
- top/bottom Wronskian proportionality (§2);
- the wall `C_5 = L[u_2]` and the exact cube criterion for homogeneous freedom,
  with the full descent-operator law `L_ℓ` (§3);
- the gauge-invariant three-level reduction `C_4 = M[u_1] + 2W(a_2,u_2)` and the
  mixed-coupling formula (4.1) with its `h ∣ a_2` consequence (§4);
- the m=0 moment identity `M(τ) = τ` at arbitrary band (§5);
- **Theorem A**: the non-cube top sector is empty (§6);
- the explicit constant-h cubic positive-cascade integrations (§7.3).

**Computed / regression only:** bounded-degree kernel triviality for non-cube
`a_3` and one-dimensionality for cube `a_3` (verifier §7); the non-polynomial
`b_1` witness for `h ∤ a_2` (verifier §9).

**Open branches (band-2 weapons do not close them here):**
- the cubic `e ≠ 0` mixed sector (empty or rigid — unknown);
- nonconstant `h` in the cubic sector (the band-2 one-line exclusion fails);
- full closure of the constant-h `e = 0` residual system (§7.3), i.e. the band-3
  analogue of the band-2 square-sector classification;
- the `a_3 = 0` sub-bands, routed to band ≤ 2 but not assembled here.

**Not claimed:** a full band-3 theorem, JC2, DC1, the quantum band-3 mirror, or
completeness inferred from computation.

## 12. Verification

From this directory:
```
uv run --with sympy python verify_classical_band3.py
```
The exact SymPy script checks: the thirteen equations against direct
differentiation with generic coefficients; every displayed `C_m`; the wall
reduction, cube integrating factor, and descent-operator laws; the `C_4`
three-level reduction; each rung of the Theorem-A cascade; the moment identity
and its membership vanishing; the mixed-coupling formula and `h ∣ a_2`
necessity; the constant-h cascade integrations; and the tame witness. A
successful run prints `ALL CLASSICAL BAND3 CHECKS PASSED`. These checks certify
the algebra only; the degree/membership arguments of §3, §6, §7.1 are the
proofs.
