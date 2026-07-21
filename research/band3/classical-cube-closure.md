# The classical band-3 cube sector: constant-h, e = 0 closure

**INDEPENDENTLY DERIVED — MACHINE-CHECKED IDENTITIES — NOT PEER REVIEWED — BAND-SCOPED**

This memo carries out the band-3 analogue of the band-2 square-sector program
(`research/band2-square-sector/classical-Astar.md` at `91a053a`;
`research/band2-m5-partial/classical-square-sector-partial.md` at `44c66d5`) for
the **constant-`h`, `e = 0`** slice of the band-3 cubic sector, the slice the
Wave-A cascade memo flagged as "the one that most closely mirrors the band-2
square sector" and left open (`research/band3/classical-band3-cascade.md` §7.3,
§11, commit `99fe6ee`). It builds directly on that memo and on the band-`k`
weapons memo (`research/band3/band-k-weapons.md`, `99fe6ee`), cited as [CASC] and
[WEAP]; the tame catalog (`research/band3/band3-tame-catalog.md`, `99fe6ee`) is
[CAT]. Every displayed algebraic identity is machine-verified by
`verify_classical_cube_closure.py` (same directory); a successful run prints
`ALL CLASSICAL CUBE CLOSURE CHECKS PASSED` (35 checks). Degree-lattice and
membership statements are written arguments; the bounded checks are regression
corroboration only.

Conventions are frozen and identical to [CASC]: over `C`, `τ = xξ`,
`{G,F} = G_ξ F_x − G_x F_ξ`, `{ξ,x} = 1`,
`F = Σ_{k=-3}^{3} x^k a_k(τ)`, `G = Σ_l x^l b_l(τ)`, membership
`τ^j | a_{-j}, b_{-j}` (`j = 1,2,3`), ladder equations
`C_m = Σ_{k+l=m}(k a_k b_l' − l a_k' b_l) = δ_{m0}`, primes `= d/dτ`. The thirteen
`C_m` are re-checked against the direct two-variable Poisson bracket (verifier §1).

## 0. The sector and the normalization

By [CASC] §2–§6 (Theorem A), an oriented band-3 pair with `a_3 ≠ 0` has `a_3`
a scalar cube `a_3 = c h³`; gauge `b_3 = 0`; the wall `C_5` gives `b_2 = e h²`;
`C_4` gives `b_1 = κ_1 h + (2e/3c) a_2/h`. This memo takes the **constant-`h`**
branch (`h = 1`, so `a_3 = c` constant) and the **`e = 0`** branch
(`b_2 = 0`), whence `b_1 = κ` a constant (writing `κ := κ_1`).

**Normalization `c = 1` is a genuine scaling.** The diagonal symplectic
`(x,ξ) → (ρx, ρ^{-1}ξ)` sends `a_k → ρ^k a_k` (verifier §0), so `a_3 = c → ρ³ c`;
choosing `ρ³ = 1/c` normalizes `c = 1` while keeping `h = 1`. We carry `c`
symbolically as a weight-3 grading bookkeeper — **every identity below holds for
all `c ≠ 0`** and specializes at `c = 1`. (Additive normalizations `β` on `b_0`
etc. are the analogue of band 2's; they are retained as free scalars.)

The free data is the F-ladder
```
q := a_2,  p := a_1,  r := a_0,  al := a_{-1},  s := a_{-2},  sig := a_{-3},
```
with membership `τ | al`, `τ² | s`, `τ³ | sig`, together with the scalars
`κ, β, γ, δ` (`γ, δ` are the integration constants of `b_{-1}`, `b_{-2}` below).

## 1. The reduction: the entire G-ladder is explicit in F

With `a_3 = c`, `b_3 = b_2 = 0`, `b_1 = κ`, the positive cascade and the moment
integrate **in closed form** (verifier §2 — each is an identity):
```
b_0    = (κ/3c) a_2 + β,                                   [C_3 = 0]
b_{-1} = (κ/3c) a_1 − (κ/9c²) a_2² + γ,                    [C_2 = 0]
b_{-2} = (κ/3c) a_0 − (2κ/9c²) a_1 a_2 + (5κ/81c³) a_2³
             − (γ/3c) a_2 + δ,                             [C_1 = 0]
b_{-3} = (1/3c)[ τ − 2 a_2 b_{-2} − a_1 b_{-1} + κ a_{-1} ].  [C_0 = 1, moment W4]
```
The `b_{-2}` line is genuinely closed-form: the `C_1` right-hand side is an
**exact `τ`-derivative** of a polynomial in `a_0, a_1, a_2` (verifier §2), so no
integral survives. Thus **`b_0, b_{-1}, b_{-2}` depend only on `a_0, a_1, a_2`**,
and `b_{-3}` additionally on `a_{-1}`, via the W4 moment identity `M = τ`.

The residual system is `C_{-1} = ⋯ = C_{-6} = 0` together with the memberships,
now equations purely in the F-ladder.

**This is the band-3 analogue of band-2's `(2.1)–(2.6)`** (`classical-square-sector-partial.md`),
with the crucial structural difference [CASC] flagged: band 2 had one free upper
level (`a_1`) above the normalized top `a_2 = 1`; band 3 has **two** free upper
levels (`a_2` and `a_0`, besides `a_1`) above the normalized top `a_3 = c`. That
extra freedom is the source of every complication below — and, as §7 shows, of
the *disappearance* of the band-2 obstruction.

## 2. C_{-1} determines a_{-2}; C_{-2} determines a_{-3}

The two trailing coefficients enter the cascade one rung at a time, each as a
pure derivative (verifier §3):
```
∂C_{-1}/∂a_{-2}' = −κ,   ∂C_{-1}/∂a_{-2} = 0,   C_{-1} does not contain a_{-3};
∂C_{-2}/∂a_{-3}' = −κ,   ∂C_{-2}/∂a_{-3} = 0.
```
Hence, for `κ ≠ 0`:
- `C_{-1} = 0` reads `κ a_{-2}' = (an explicit polynomial in a_2,a_1,a_0,a_{-1})`
  and **determines `a_{-2}`** as its antiderivative;
- `C_{-2} = 0` then **determines `a_{-3}`** likewise.

(That `κ ≠ 0` is forced is the band-3 mirror of band-2's `κ = 0` exclusion,
`classical-square-sector-partial.md` §1; the constant-`h` `κ = 0` sub-branch is
routed but not re-proved here.) After these two determinations the residual is
`C_{-3} = C_{-4} = C_{-5} = C_{-6} = 0` — the **closing constraints** — plus the
memberships on `a_{-2}, a_{-3}` and the induced memberships on `b_{-1}, b_{-2}, b_{-3}`.

## 3. Φ₁ : the W3-replacement first integral (the honest band-3 W3)

[CASC] §8 records that the band-2 W3 free-trailing-coefficient integral is
**obstructed** at band 3 in the raw ladder: `C_{-1} = d/dτ[Σ_k k a_k b_{-1-k}] +
Σ_k a_k' b_{-1-k}` and the residual is not a total derivative. The resolution is
that after the **full constant-`h`, `e = 0` parametrization** of §1 the trailing
integral reappears — with a single nonlocal generator. Define
`Q₁ := ∫ a_2` (`Q₁' = a_2`). Then (verifier §4)

> **Lemma 3.1 (band-3 first integral of `C_{-1}`).**
> ```
> Φ₁ := 2δ a_1 + γ a_0 − κ a_{-2}
>       − (5δ/3c) a_2² − (4γ/3c) a_1 a_2 + (2κ/3c) a_{-1} a_2 + (2κ/3c) a_1 a_0
>       + (τ a_2)/c − Q₁/(3c)
>       + (14γ/27c²) a_2³ − (5κ/9c²) a_1² a_2 − (5κ/9c²) a_2² a_0
>       + (40κ/81c³) a_1 a_2³ − (22κ/243c⁴) a_2⁵
> ```
> satisfies `Φ₁' = C_{-1}` identically.

The obstruction that forces the nonlocal term is a *constant* variational
density: the Euler operator gives `E_{a_2}(C_{-1}) = −1/(3c)` (verifier §4), and
the combination `τ a_2 − Q₁` (whose derivative is `τ a_2'`) is exactly what
absorbs it. This is the precise band-3 analogue of band-2's Lemma 2.1, in which
`Φ` carried the nonlocal `P₁ = ∫ a_1`; **at band 3 the nonlocal generator is
`Q₁ = ∫ a_2`** — one level up, matching the cube's extra free level. On any
solution `Φ₁ = Φ₁(0)` is constant; because `a_{-2}, a_{-1}, Q₁` vanish at `τ = 0`
(memberships and `Q₁(0)=0`), `Φ₁(0)` is a fixed polynomial in the values
`a_2(0), a_1(0), a_0(0)`. Lemma 3.1 **is** the equation determining `a_{-2}` of §2,
in integrated form.

## 4. The branch tree

The bottom two rungs give the reflection-mirror of the top cascade
([WEAP] W1-mirror, W6), now as constraints because `b_{-2}, b_{-3}` are already
determined by §1.

**Bottom Wronskian (`C_{-6}`).** `C_{-6} = 3(a_{-3}' b_{-3} − a_{-3} b_{-3}')`;
`C_{-6} = 0` with `a_{-3} ≠ 0` gives `b_{-3} = μ_3 a_{-3}`, `μ_3 ∈ C`
(verifier §5). This is the cross-coupling datum of [CAT] §1.

**Bottom wall (`C_{-5}`, the cube-class mirror).** Substituting `b_{-3} = μ_3 a_{-3}`,
with `φ := μ_3 a_{-2} − b_{-2}`,
```
C_{-5} = 3 a_{-3} φ' − 2 a_{-3}' φ  =  L_bot[φ],
```
the exact reflection of the top wall `C_5 = L[u_2] = 3 a_3 u_2' − 2 a_3' u_2`
([CASC] §3). Its integrating factor is the **cube** one (verifier §5):
`(φ³/a_{-3}²)' = φ² L_bot[φ] / a_{-3}³`. Hence `C_{-5} = 0` with `a_{-3} ≠ 0`
forces `φ³ = c' a_{-3}²`, so either
- `a_{-3}` is a scalar cube and `φ = e_bot·(cube root)²` (bottom cube class), or
- `a_{-3}` is not a scalar cube and `φ = 0`, i.e. **`b_{-2} = μ_3 a_{-2}`** — the
  mirror of the top `e = 0` (`b_2 = 0`).

The solver-facing tree (mirroring [CAT] §4, [CASC] §9):
```
constant-h, e=0 cube sector  (a_3 = c const, b_3=b_2=0, b_1=κ)
│
├─ a_{-3} = 0  ................ ONESIDED-TOP.  C_{-6}=0 automatically; b_{-3}
│     determined by the moment with τ³|b_{-3}.  Contains the tame witness
│     F = ξ + α x³ + β x² + γ x + δ, G = −x.  Routes toward band ≤ 2 on the
│     negative side.  [families exhibited; full onesided closure routed]
│
└─ a_{-3} ≠ 0  ⇒ b_{-3} = μ_3 a_{-3}  (C_{-6})
   │
   ├─ μ_3 = 0  (b_{-3} = 0) ..... B0-band3.  Contains B3-2, B3-3 (§6).
   │      Bottom wall: φ = −b_{-2}; C_{-5} forces b_{-2} = 0 unless a_{-3} a
   │      scalar cube.  [tame stratum; positive controls live here]
   │
   └─ μ_3 ≠ 0  (b_{-3} ≠ 0) ..... A*-band3  (resistant).  Bottom wall forces
          b_{-2} = μ_3 a_{-2} (a_{-3} not a cube) or the bottom cube class.
          [§7: closed by the tropical degree gap except on one leading-coeff
           locus; bounded-empty by [CAT] §5]
```
`A0-band3` (`μ_3 = 0` but the mid-level `b_{-2} = μ_3 a_{-2}` proportionality
nontrivial) is the intermediate rung, the band-3 image of band-2's `A0`.

## 5. The `a_2 = const` reduction (the stratum of the positive controls)

Every positive control and the witness has `a_2` **constant** (§6). Setting
`a_2 = Q_0` (const) makes §2 fully explicit (verifier §6):
- `a_{-2}` becomes an explicit polynomial in `a_1, a_0, a_{-1}` (the antiderivative
  of `C_{-1}` closes locally, since the `Q₁ = ∫a_2` term degenerates to
  `Q_0 τ`);
- `a_{-3}` becomes explicit in `a_1, a_0, a_{-1}` **with one nonlocal generator
  `P₁ = ∫ a_1`** — exactly the band-2 nonlocal `∫ a_1`. So the band-2 nonlocal
  reappears one rung lower, at the level that determines `a_{-3}`.

The closing constraints reduce to `C_{-3} = C_{-4} = C_{-5} = C_{-6} = 0` in
`(a_1, a_0, a_{-1})` and the scalars, with memberships. This is the tightest
band-3 system that still contains all known solutions.

## 6. Positive controls and the witness

Gauged into this sector (subtract `λ_3 F`), the catalog pairs [CAT] §4 land as
follows (verifier §7, all `{G,F} = 1` re-checked exactly):

- **B3-2** `→` `c = 1, κ = −1, e = 0`, ladder
  `a_2 = 1, a_1 = 3τ+1, a_0 = 2τ, a_{-1} = 3τ²+2τ, a_{-2} = τ², a_{-3} = τ³`;
  gauged `b_0 = 0, b_{-1} = −τ, b_{-2} = 0, b_{-3} = 0`. Branch **B0-band3**
  (`μ_3 = 0`), degrees `(deg a_1, deg a_0, deg a_{-1}) = (1,1,2)`.
- **B3-3** `→` `c = 1, κ = −1, e = 0, a_2 = 0`, ladder
  `a_1 = 3τ+2, a_{-1} = 3τ²+3τ, a_{-3} = τ³` (odd-only); gauged
  `b_{-1} = −τ, b_{-2} = b_{-3} = 0`. Branch **B0-band3**, degrees `(1,·,2)`.
- **Tame witness** `F = ξ + α x³ + β x² + γ x + δ, G = −x`: `a_{-2} = a_{-3} = 0`,
  `a_{-1} = τ`, `κ = −1`. Branch **onesided-top**.

The verifier checks that the §1 reduction formulas reproduce the B3-2 and B3-3
G-ladders exactly from their `(κ,β,γ,δ)`, and that all three satisfy every
`C_m = δ_{m0}`. **Both catalog B0-band3 controls therefore appear inside the
classification, with `b_{-2} = b_{-3} = 0`.**

## 7. The W5 verdict: a tropical gap, not a congruence — and the band-3 modulus

The strategic question [WEAP] W5 poses is whether the residual **degree** system
obstructs the resistant branch, as it did for band-2's `A*` (there `Φ` and `I₂`
gave `2V = P + W` and `V + W = 2P + 1`, whose sum `3(V − P) = 1` is infeasible
over `Z`: the **mod-3 kill**). [WEAP] flags the band-3 modulus as the key unknown:
`k+1 = 4` vs `2k−1 = 5` vs other.

The resistant branch is `A*-band3`: `a_{-3} ≠ 0` and `μ_3 ≠ 0`, so by the bottom
Wronskian (§4) `b_{-3} = μ_3 a_{-3}`, forcing the exact degree equality
```
deg b_{-3}  =  deg a_{-3}.                                      (RES)
```
Now compare how the two sides scale with `P := deg a_1` in the `a_2 = const`
reduction (§5). Write `R := deg a_0`, `L := deg a_{-1}`.

The explicit `a_{-3}` (from `C_{-2}`) and `b_{-3}` (from the moment) carry three
top terms whose `τ`-degrees dominate **variable by variable** (verifier §8):
```
              a_1-term       a_0-term      a_-1-term
   a_{-3} :  a_1³  (3P)     a_0²  (2R)    a_{-1}a_1 (L+P)   [coeffs −5/27c², 1/3c, 2/3c]
   b_{-3} :  a_1²  (2P)     a_0   (R)     a_{-1}   (L)
```
so `deg a_{-3} ≥ max(3P, 2R, L+P)` and `deg b_{-3} ≤ max(2P, R, L, 1)`. For
`P ≥ 1` each entry of the lower bound strictly exceeds the matching entry of the
upper bound (`3P > 2P`, `2R > R` for `R ≥ 1`, `L+P > L`), so
`max(3P,2R,L+P) > max(2P,R,L,1)` **unconditionally**, giving
```
deg b_{-3} < deg a_{-3},
```
which **contradicts (RES)** — *unless* the top of `a_{-3}` cancels. Its top sits at
degree `max(3P, 2R, L+P)`, and cancels only where that maximum is a **tie**
between two of the three terms, with the leading coefficients summing to zero:
```
3P = 2R        with   −(5/27c²) lc(a_1)³ + (1/3c) lc(a_0)² = 0,   or
L  = 2P        with   −(5/27c²) lc(a_1)³ + (2/3c) lc(a_{-1}) lc(a_1) = 0,   or
2R = L + P     with   (1/3c) lc(a_0)² + (2/3c) lc(a_{-1}) lc(a_1) = 0.
```
This is a **finite union of leading-coefficient cancellation loci** (codimension
one each); off all of them the degree gap closes `A*-band3` at arbitrary degree.

> **Finding (band-3 modulus = "other": tropical, not congruential).** For the
> constant-`h`, `e = 0` sector (in the `a_2 = const` reduction of §5) the
> resistant branch `A*-band3` is obstructed by a **max-degree (tropical)
> inequality** `deg b_{-3} < deg a_{-3}`, with the gap
> equal to the cube's signature `3·deg a_1 − 2·deg a_1 = deg a_1` — the "`3` vs
> `2`" of `a_3 = c h³` (moment `a_1 b_{-1} ∼ a_1²` against `a_{-3} ∼ a_1³`). This
> is **not** a linear congruence `g·(linear) = m`, so the band-2 mod-3 kill has
> **no band-3 analogue as a modulus**: the answer to [WEAP]'s "`k+1 = 4` vs
> `2k−1 = 5` vs other" is **other** — the obstruction is the degree *gap* `3:2`
> inherited from the cube, not a residue class. (The independent **wall
> divisibility** `3 ∣ deg a_{±3}`, from the cube integrating factor `u³ = c a_3²`,
> `gcd(3,2)=1`, verifier §5, is consistent with the tame stratum and is not itself
> an emptiness obstruction.)

The obstruction closes `A*-band3` **off** the cancellation loci at arbitrary
degree by pure degree-counting. What survives is exactly the finite union of
leading-coefficient loci above — the boundaries where the top term of `a_{-3}` is
killed by a coefficient tie. This is the residual **leading-coefficient** problem;
it is *not* a lattice congruence, so no modulus disposes of it. Bounded evidence
that it too is empty is strong: [CAT] §5 enumerated 1216 genuine band-3 tame pairs
with **0** resistant, and every positive control (§6) is `B0-band3` with
`b_{-3} = 0`.

**Why band 2's congruence evaporates.** Band 2 normalized its top (`a_2 = h²`) and
had no free upper F-function; the residual lattice was rigid enough for the
`Φ, I₂` rows to collide mod 3. Band 3 keeps **two** free upper levels, `a_2` and
`a_0`, above the normalized `a_3 = c`. The extra freedom converts the band-2
*congruence* into a band-3 *degree gap*: the same cube that made the sector hard
to integrate (the `a_1³` growth of `a_{-3}`) is what out-scales `b_{-3}` and kills
the resistant branch by inequality rather than by residue. This sharpens the
"open" verdict of [CASC] §11: the resistant branch is now closed **except** on one
explicit leading-coefficient locus, and the mechanism is identified.

## 8. Claim disposition

**Proved at arbitrary degree (machine-checked identities):**
- the thirteen `C_m` vs the direct Poisson bracket (verifier §1);
- the scaling law justifying `c = 1` (verifier §0);
- the complete reduction: `b_0, b_{-1}, b_{-2}` explicit in `a_0,a_1,a_2` and
  `b_{-3}` via the moment, with `C_3=C_2=C_1=0`, `C_0=1` identically (§1);
- `C_{-1}` determines `a_{-2}`, `C_{-2}` determines `a_{-3}`, each entering
  linearly as `−κ·(deriv)` (§2);
- **Φ₁**, the first integral `Φ₁' = C_{-1}` with nonlocal `Q₁ = ∫a_2` — the
  honest band-3 replacement for the obstructed W3 (§3);
- the branch tree: bottom Wronskian `b_{-3} = μ_3 a_{-3}` (`C_{-6}`), the bottom
  cube-wall `C_{-5} = L_bot[φ]` with its integrating factor and cube-class
  criterion, and the wall divisibility `3 ∣ deg a_{±3}` (§4, §7);
- the `a_2 = const` reduction with the nonlocal `P₁ = ∫a_1` in `a_{-3}` (§5);
- the positive controls B3-2, B3-3 reproduced inside the classification, both
  B0-band3 with `b_{-2}=b_{-3}=0`, and the tame witness in the onesided branch (§6).

**Determined (this memo's headline quantitative results):**
- the **band-3 modulus verdict**: the resistant branch `A*-band3` is obstructed by
  a **tropical degree gap** `deg b_{-3} < deg a_{-3}` (the cube's `3:2` scaling in
  `deg a_1`), **not** a linear congruence; the answer to [WEAP]'s "`k+1 = 4` vs
  `2k−1 = 5` vs other" is **other** — no residue-class modulus exists at band 3.
  The independent wall divisibility is `3 ∣ deg a_{±3}` (`= k`), consistent with
  the tame stratum (§7). *(The degree-scaling facts — `deg a_{-3} ≥ 3 deg a_1`,
  `deg b_{-3} ≤ 2 deg a_1` in `a_1` — are machine-checked, verifier §8; the
  termwise-domination conclusion is the written argument.)*
- **A*-band3 closed off a finite set of cancellation loci**: for `deg a_1 ≥ 1`,
  no resistant pair exists unless `(deg a_1, deg a_0, deg a_{-1})` lies on one of
  the ties `3P=2R`, `L=2P`, `2R=L+P` with the matching leading-coefficient
  relation (§7) — the only places `a_{-3}`'s top term is killed.

**Computed / regression only:**
- the Newton degree-forms of the closing constraints and the balanceable stratum
  `deg a_{-1} = 2 deg a_1` that carries the controls (this memo's lattice sweep);
- emptiness of the resistant branch at bounded degree — [CAT] §5's 1216-pair
  enumeration, **0** resistant (the residual leading-coefficient loci are *not*
  independently swept here).

**Open (a genuinely new argument, not a band-2 weapon, is required):**
- the residual **leading-coefficient loci** of `A*-band3` (the finite union of
  ties `3P=2R`, `L=2P`, `2R=L+P` with the matching `lc` relations) — a coefficient
  problem, not a congruence, so no modulus disposes of it (§7);
- full closure of the `a_2` **nonconstant** sub-branch (all controls have `a_2`
  constant; whether `deg a_2 ≥ 1` is empty or admits a new family is not settled —
  the `Φ₁` top term `a_2⁵` forces cancellation loci but the casework is not
  completed);
- the `κ = 0` and full onesided (`a_{-3} = 0`) sub-branches, routed to band ≤ 2;
- the `e ≠ 0` and nonconstant-`h` sectors ([CASC] §7.1, §7.2 — out of scope here).

**Not claimed:** a full band-3 theorem, closure of the resistant branch, JC2,
DC1, or the quantum mirror.

## 9. Verification

```
uv run --with sympy python research/band3/verify_classical_cube_closure.py
```
Exact SymPy; 35 `PASS` lines; ends `ALL CLASSICAL CUBE CLOSURE CHECKS PASSED`.
The script certifies the algebra (reduction, `Φ₁`, the branch-tree operators and
integrating factors, the `a_2 = const` reduction, the controls, the degree-scaling
coefficients `a_1³, a_0², a_{-1}a_1` of `a_{-3}`, and the controls' degree pattern
and cube divisibility); the tropical-gap verdict of §7 and the membership
arguments are the written proofs.
