# Closing the A\*-band3 resistant branch (μ ≠ λ), both faces

**INDEPENDENTLY DERIVED — MACHINE-CHECKED IDENTITIES — NOT PEER REVIEWED — BAND-SCOPED**

This memo attacks the `μ ≠ λ` resistant branch **A\*-band3** — the band-3 analogue
of the band-2 square-sector A\* branch — on both the classical (Jacobian, JC2)
and quantum (Weyl/Dixmier, DC1) faces. The catalog established A\*-band3 as
**counterexample-or-nothing** (`band3-tame-catalog.md`, commit `99fe6ee`: the
single-cubic-shear normal form forces `μ = λ`, and decoupled `±3` sources need
band ≥ 9, so no tame word reaches it). The task is to close it from the
coefficient equations, upgrading the band-2 A\* playbook
(`research/band2-square-sector/classical-Astar.md`, `84978b9`).

Everything displayed is machine-verified by `verify_astar_band3.py` (same
directory); a successful run ends `ALL ASTAR BAND3 CHECKS PASSED` (48 checks,
~6 s). Conventions are frozen exactly as in the Wave-A band-3 memos (commit
`99fe6ee`): classical `{G,F}=G_ξF_x−G_xF_ξ`, `τ=xξ`,
`F=Σ_{k=-3}^{3}x^k a_k(τ)`, `G=Σ x^k b_k(τ)`,
`C_m=Σ_{k+l=m}(k a_k b_l' − l a_k' b_l)=δ_{m0}`, membership `τ^j∣a_{-j},b_{-j}`;
quantum `A_1[x^{-1}]=⊕_k x^k C[E]`, `(x^a f)(x^b g)=x^{a+b}f(E+b)g(E)`,
`Q_m=Σ_{k+l=m}(b_l^{[k]}a_k − a_k^{[l]}b_l)=δ_{m0}`, membership
`E(E−1)⋯(E−r+1)∣a_{-r},b_{-r}`.

## 0. Headline

> **The band-2 A\* playbook does not transfer verbatim — and the reason is
> instructive.** Band-2 A\* died of a `mod 3` congruence (`3V=3P+1`) read off two
> clean first integrals `Φ, I₂` of the negative cascade. At band 3 **those first
> integrals do not exist**: no constant-coefficient combination of `C_{-1..-5}`
> is an exact `τ`-derivative (the W3-obstruction of `classical-band3-cascade.md`
> §8, here made quantitative), so there is no band-3 `L·d=β` lattice of the
> band-2 shape and no "band-3 modulus" in that sense.
>
> **A different, stronger weapon closes half of A\*: reflection.** The geometric
> reflection `R:(x,ξ)↦(ξ,x)` sends a genuine A\* pair to a genuine *gauged*
> band-3 Keller pair whose top coefficient is `τ^{-3}a_{-3}`. **Theorem A**
> (the non-cube top sector is empty, `classical-band3-cascade.md` §6, PROVED)
> applied to the reflected pair yields:
>
> **Theorem (classical A\*-I empty).** In branch A\*-band3, if `a_{-3}` is **not**
> a scalar cube then no pair exists. Consequently every classical A\*-band3 pair
> has **both** extremes scalar cubes and is reflection-symmetric.
>
> This is a genuine theorem (machine-checked identities + Theorem A). The
> residual **A\*-II** (both extremes cubes) is reduced to a coupled two-sided
> cube cascade; its `nonconstant-ĥ` part coincides with the still-open
> nonconstant-`h` sector (§7.1 there), and its minimal (doubly-constant) sub-case
> is empty on exact bounded sweeps. **No counterexample survives** (the one
> pole-free candidate the solver emits fails `C_0` for `μ≠0`). On the **quantum**
> face the reflection needs a falling-factorial twist and — because the quantum
> "Theorem A" is itself open (the `Q_5` wall does not force a shifted cube,
> `E(E−2)(E−4)`) — the classical route does **not** transcribe; quantum A\*-I is
> *reduced to the open quantum non-cube question*, with the same bounded
> emptiness evidence and inhomogeneous `μ₃`-tail structure.

## 1. The branch, exactly

Orient (pair exchange / reflection) so `a_3 ≠ 0`; spend the single gauge
`G↦G−λ₃F` to set `b_3 = 0` (`C_6`, W1). By **Theorem A**
(`classical-band3-cascade.md` §6) `a_3` is then a scalar cube `a_3 = c h³`,
`c≠0`. The bottom Wronskian `C_{-6}=0` gives `b_{-3}=μ̃ a_{-3}` with `μ̃∈C` the
gauged residue (`μ̃ = μ₃ − λ₃`); the single gauge is already spent, so `μ̃` is
**not** removable. The branch split (catalog §4) is

```
onesided-top (a_{-3}=0) · B0-band3 (a_{-3}≠0, μ̃=0) · A*-band3 (a_{-3}≠0, μ̃≠0).
```

> **A\*-band3.** A genuine band-3 Keller pair, gauge `b_3=0`, with
> `a_3=c h³ ≠ 0`, `a_{-3}≠0`, and `b_{-3}=μ̃ a_{-3}`, `μ̃≠0`.

`verify §2` checks `b_{-3}=μ̃a_{-3}` solves `C_{-6}` identically. The **bottom
wall** `C_{-5}` reduces exactly (reflection of the H5 wall reduction) to
```
C_{-5} = −L̃[u₂],   u₂ := b_{-2} − μ̃ a_{-2},   L̃[u] := 3 a_{-3} u' − 2 a_{-3}' u
```
(`verify §2`), with the cube integrating factor `(u³/a_{-3}²)' = u²L̃[u]/a_{-3}³`.
Since `τ³∣a_{-3}` forces `deg a_{-3}≥3≥1`, the cube criterion (Prop. 3.2 reflected)
splits A\* into

- **A\*-I** (`a_{-3}` **not** a cube): forced `u₂=0`, i.e. `b_{-2}=μ̃ a_{-2}`.
- **A\*-II** (`a_{-3}=c̃ ĥ³`, `τ∣ĥ`, a cube): `u₂ = ẽ ĥ²`, `ẽ∈C`.

## 2. The reflection lemma (the new weapon)

The geometric reflection `R:(x,ξ)↦(ξ,x)`, `RF(x,ξ)=F(ξ,x)`, is a `C`-algebra
automorphism of `C[x,ξ]` with Jacobian `−1`. Two exact facts (`verify §1`,
generic coefficients, and on the genuine pair B3-2):
```
(RF)_j = τ^{-j} a_{-j},        {RG,RF} = −{G,F} ∘ R.
```
Both are proved identities. `RF` is genuine because `(RF)_{-j}=τ^{j}a_j` is
divisible by `τ^j`; `(RF)_j=τ^{-j}a_{-j}` is a genuine polynomial because
`τ^j∣a_{-j}`. So **reflection carries genuine band-3 pairs to genuine band-3
pairs, flipping the bracket sign and reversing the band.** In particular
`(RF)_3 = τ^{-3}a_{-3}` is a nonzero polynomial (of degree `deg a_{-3}−3`), and
`{RG,RF}=−1`.

## 3. Theorem: branch A\*-I is empty

> **Theorem A\*-I.** No classical band-3 Keller pair lies in A\* with `a_{-3}`
> not a scalar cube. Equivalently, in A\*-band3 both extremes `a_3, a_{-3}` are
> scalar cubes.

*Proof.* Suppose `(F,G)` is such a pair. Reflect and form the **gauge**
```
H := μ̃ RF − RG.
```
Then (all `verify §3`):
1. `H_j = τ^{-j}(μ̃ a_{-j} − b_{-j})`, so `H` is a genuine band-3 element and
   `H_{-j}=τ^{j}(μ̃ a_j − b_j)` is `τ^j`-divisible (membership holds);
2. `H_3 = τ^{-3}(μ̃ a_{-3} − b_{-3}) = 0`, using `b_{-3}=μ̃ a_{-3}` — the
   reflected top is **gauged flat**;
3. `{H,RF} = {μ̃ RF − RG, RF} = −{RG,RF} = +1`.

Thus `(RF, H)` is a genuine band-3 Keller pair with `{H,RF}=1`, top gauged
`H_3=0`, and top coefficient `(RF)_3=τ^{-3}a_{-3}≠0`. Because `τ³` is a cube,
`a_{-3}=τ³·(RF)_3` is a scalar cube **iff** `(RF)_3` is (`verify §3`). If `a_{-3}`
is **not** a cube then `(RF)_3` is not a cube, and **Theorem A** (whose six rung
reductions are re-verified in `verify §3` and whose closing degree/membership
argument is `classical-band3-cascade.md` §6) makes the pair `(RF,H)` **empty** —
contradicting its construction. Hence `a_{-3}` is a cube. The mirror statement for
`a_3` is the original Theorem A. ∎

The reflected top is itself in A\* (`H_{-3}=μ̃ (RF)_{-3}`, `verify` structure), so
reflection maps A\* to A\* — it does not iterate to a contradiction, but it
**closes the non-cube-bottom half outright** and reduces A\* to the cube/cube
sector. Note Theorem A\*-I is unconditional in `h` (it does not use the
constant-`h` slice).

## 4. Branch A\*-II: reduction and bounded evidence

With both extremes cubes, run the top cube cascade (`classical-band3-cascade.md`
§7.3) *and* its reflected image. In A\*-I (`u₂=0`) the next rung is the reflected
Theorem-A M-rung
```
C_{-4} = M̃[u₁],   u₁ := μ̃ a_{-1} − b_{-1},   M̃[u]=3 a_{-3} u' − a_{-3}' u
```
(`verify §4`) — the exact bottom mirror of `C_4=M[u_1]`. In A\*-II the reflected
cube cascade gives the mirror mixed-coupling relation, coupling the top data
`(a_2,a_1,a_0)` to the bottom data `(a_{-1},a_{-2})` through `M=τ` (W4).

**Constant-`h` slice (`a_3=c`, `e=0`).** The §7.3 parametrization solves
`C_5,C_4,C_3,C_2` identically (`verify §4`), and the moment `M=τ` holds with
`M(0)=0` by membership (`verify §4`). Within this slice:

- **Doubly-minimal A\*** (`a_{-3}=c̃ τ³`, so reflected-`ĥ` is also constant): the
  bottom wall `C_{-5}=0` is built in via `b_{-2}=μ̃a_{-2}+ẽτ²`, `C_{-6}=0`, and
  `M=τ` holds by construction (`verify §5`). Exact bounded-degree sweeps over
  degree profiles `(deg a_2,a_1,a_0,a_{-2})` and several `(μ̃,κ₁)` find **no
  genuine A\* pair** (`verify §5`; the standalone sweep covered 84 configurations,
  0 hits). This is the closest analogue of the band-2 constant-`h` A\* — but a
  **computation, not a proof**.
- **Nonconstant-`ĥ` A\*-II** (`deg ĥ≥2`): the reflected pair is a nonconstant-`h`
  cube-sector pair, whose emptiness is the **OPEN** branch of
  `classical-band3-cascade.md` §7.1. So this part of A\*-II is **routed** to that
  open sector, not closed.

**No counterexample.** The only pole-free profile the linear solver emits with
`μ̃` symbolic (`a_2=a_{-2}=b_{-2}=0`, `a_1` constant) is spurious: its moment
gives `C_0 = 1 + 9 μ̃ τ²`, which fails `C_0=1` for every `μ̃≠0` (`verify §5`).
Verified to destruction: it is a `μ̃=0` (B0) datum, not an A\* pair.

## 5. Why the band-2 W5 lattice does not transfer (and the redundancy that does)

Band-2 A\* was killed by two exact first integrals — `Φ` of `C_{-1}` and `I₂` of
`C_{-2}` — whose leading degrees gave `2V=P+W`, `V+W=2P+1`, hence `3V=3P+1`,
infeasible `mod 3` (`classical-Astar.md` §3; W5 framework
`band-k-weapons.md`). At band 3 this mechanism **fails at the source**:

- **No first integral exists.** A machine search over a rich basis (algebraic
  monomials to degree 3 in the surviving coefficients, plus all relevant
  antiderivatives `∫a_1,∫a_2,∫a_2²,∫a_1a_2,∫a_0,∫a_2³,…` and `τ`-weightings)
  finds **no** exact `τ`-derivative among `C_{-1}`, `C_{-2}`, nor any
  constant-coefficient combination of `C_{-1..-5}` — every candidate multiplier
  is forced to `0`. (The same search **does** reproduce band-2's `Φ`.) This is
  the W3-obstruction of `classical-band3-cascade.md` §8 made quantitative: the
  residual `Σ_k a_k' b_{-1-k}` is not a total derivative of the raw coefficients.
  Consequently there is **no band-3 `L·d=β` of the band-2 shape**, and the
  "determine the band-3 modulus" sub-goal resolves negatively: the clean
  congruence lattice is obstructed, and the closure runs through **reflection +
  Theorem A** instead. (Machine evidence for the search lives in the development
  scripts; the *positive* structural facts it certifies — the reflection lemma,
  the reduced rungs — are in `verify`.)
- **The redundancy that *does* transfer.** Band-2's A\* proof "never used `C_{-3}`."
  The band-3 doubly-minimal A\* system stays empty on the sweep **after dropping
  `C_{-3}`** — and even after dropping `C_1` (`verify §5`). So the load-bearing
  set is `{moment, C_{-1}, C_{-2}, C_{-4}, memberships}`, exactly one level-shifted
  from band 2. The bottom wall `C_{-5}` and Wronskian `C_{-6}` are consumed by the
  A\*-I reflection step; the middle rungs `C_{-3}, C_1` are redundant.

## 6. The quantum face (DC1)

Setup: gauge `b_3=0` (`Q_6`, 3-periodicity); `Q_{-6}=0` gives `b_{-3}=μ₃ a_{-3}`;
A\* is `μ₃≠0` (`verify §6`). The positive cascade for the trivial-cube top
`a_3=1`, `b_2=0` has `Q_4 = b_1^{[3]}−b_1`, so `b_1` is a constant (the 3-fold
periodicity replacing band-2's 2-fold), and `Q_3,Q_2` determine `b_0,b_{-1}`
(`verify §6`) — the exact quantum analogue of §4.

**The genuinely new feature: the inhomogeneous tail.** With `b_{-3}=μ₃a_{-3}`,
```
Q_{-5} = [ b_{-2}^{[-3]}a_{-3} − a_{-3}^{[-2]}b_{-2} ] + μ₃·[ a_{-3}^{[-2]}a_{-2} − a_{-2}^{[-3]}a_{-3} ],
Q_{-4} = [ b_{-1}^{[-3]}a_{-3} − a_{-3}^{[-1]}b_{-1} ] + μ₃·[ a_{-3}^{[-1]}a_{-1} − a_{-1}^{[-3]}a_{-3} ] + [ b_{-2}^{[-2]}a_{-2} − a_{-2}^{[-2]}b_{-2} ]
```
(`verify §6`, exactly the `quantum-band3-cascade.md` §5 decompositions): the
bottom wall and the level-`−4` rung acquire a **`μ₃`-proportional source** — the
top/bottom cross-coupling with **no band-2 shadow**. This is the quantum obstacle
the archived milestone flagged, here reproduced exactly.

**Why the classical route does not transcribe.** The classical kill used
reflection to import Theorem A. Quantumly (i) the reflection `E↦−E−1` needs a
falling-factorial twist — the naive substitution sends the trivial-cube top
`a_3=1` to a reflected `x^{-3}`-coefficient `1`, which is **not** divisible by
`E(E−1)(E−2)`, breaking `A_1`-membership (`verify §6`); and, decisively,
(ii) the quantum analogue of **Theorem A is open**: the `Q_5` wall does **not**
force a shifted cube — `a_3=E(E−2)(E−4)` solves it (`quantum-band3-cascade.md`
§3.2). So even a correctly-twisted quantum reflection would only reduce quantum
A\*-I to *"is the reflected non-shifted-cube-necklace top killed downstream?"* —
precisely the open question of `quantum-band3-cascade.md` §3.5. **Quantum A\*-I is
therefore reduced to the open quantum non-cube question, not closed.**

**Bounded evidence.** Exact quantum sweeps (trivial-cube top, minimal
`a_{-3}=c̃ E(E−1)(E−2)`, `b_{-3}=μ₃a_{-3}`, numeric `μ₃≠0`, degree profiles
`(dP,dN)`) find **no genuine A\* pair** (`verify §6`; standalone sweep 24
configurations, 0 hits). The machinery is sound: with `μ₃=0` the *same* solver
returns genuine B0 pairs (`verify §6`).

## 7. Solver-facing case tree

```
A*-band3  (gauge b3=0; a3=c h^3≠0; a_-3≠0; b_-3=μ̃ a_-3, μ̃≠0)
│
├─ A*-I  (a_-3 NOT a cube)
│     ⇒ EMPTY.  [Theorem A*-I, §3 — PROVED via reflection R + Theorem A,
│                all identities machine-checked]
│
└─ A*-II (a_-3 = c̃ ĥ^3 a cube;  by §3 also a3 = c h^3, both cubes)
      │  bottom wall: b_-2 = μ̃ a_-2 + ẽ ĥ^2
      │
      ├─ ĥ constant (doubly-minimal, a_-3 = c̃ τ^3)
      │     constant-h reflected sector; residual {moment,C_-1,C_-2,C_-4};
      │     EMPTY on exact bounded sweeps (84 configs, 0 hits) — COMPUTED, not proved
      │
      └─ ĥ nonconstant (deg ĥ ≥ 2)
            reflected pair is a NONCONSTANT-h cube-sector pair
            ⇒ ROUTED to the open nonconstant-h branch (classical-band3 §7.1)

QUANTUM A*-band3  (gauge b3=0; b_-3=μ₃ a_-3, μ₃≠0; INHOMOGENEOUS μ₃-tail)
   reflection needs falling-factorial twist; quantum Theorem A OPEN
   ⇒ A*-I ROUTED to the open quantum non-cube question (quantum-band3 §3.5)
   ⇒ EMPTY on exact bounded sweeps (24 configs, 0 hits) — COMPUTED
```

## 8. Claim disposition

**Proved (machine-checked identities + cited PROVED inputs):**
- the reflection lemma `(RF)_j=τ^{-j}a_{-j}`, `{RG,RF}=−{G,F}∘R` (§2);
- the A\* setup (`b_{-3}=μ̃a_{-3}` solves `C_{-6}`), the bottom-wall reduction
  `C_{-5}=−L̃[u₂]` and its cube split, and `C_{-4}=M̃[u₁]` in A\*-I (§1,§4);
- **Theorem A\*-I: branch A\*-I is empty**, hence A\* forces both extremes to be
  scalar cubes (§3), via reflection `R` + gauge `H=μ̃RF−RG` + Theorem A (whose
  rung reductions are re-verified here; degree/membership argument in
  `classical-band3-cascade.md` §6, commit `99fe6ee`);
- the quantum A\* setup, the positive-cascade `Q_4=b_1^{[3]}−b_1`, and the
  **inhomogeneous** `Q_{-5},Q_{-4}` `μ₃`-source decompositions (§6);
- the quantum reflection membership-break for constant `a_3` (§6).

**Computed / corroboration only (exact SymPy, bounded):**
- no genuine A\* pair on classical doubly-minimal sweeps (84 configs) or quantum
  sweeps (24 configs); the redundancy of `C_{-3}` and `C_1`; the destruction of
  the single spurious pole-free candidate (`C_0=1+9μ̃τ²`);
- the non-existence of a band-2-style first integral / `L·d=β` lattice at band 3
  (basis search; the W3-obstruction made quantitative).

**Open (precisely delimited — A\* is NOT fully closed):**
- **A\*-II nonconstant-`ĥ`** (classical): coincides with the open nonconstant-`h`
  cube sector (`classical-band3-cascade.md` §7.1);
- **A\*-II doubly-minimal** (classical): believed empty (sweeps), but lacks a
  degree/first-integral proof — the band-3 first-integral obstruction blocks the
  band-2 argument;
- **Quantum A\*-I**: reduced to the open quantum non-shifted-cube question
  (`quantum-band3-cascade.md` §3.5); the whole quantum A\* is counterexample-or-
  nothing with bounded emptiness evidence only.

**Not claimed:** a full A\*-band3 emptiness theorem on either face; the mixed
`e≠0` sector; JC2; DC1; any unbounded-degree completeness.

## 9. Verification

```sh
uv run --with sympy python research/band3/verify_astar_band3.py
```
Exact SymPy. §0 machinery (`C_m`, `Q_m` vs direct brackets); §1 reflection lemma
(generic + B3-2); §2 A\* setup + bottom-wall reduction + cube factor; §3 the
**A\*-I emptiness theorem** (H-construction, `H_3=0`, `{H,RF}=1`, cube-equivalence,
Theorem-A rungs); §4 reflected `C_{-4}`, positive cascade, moment; §5
doubly-minimal residual system, bounded sweeps, `C_{-3}/C_1` redundancy, spurious-
candidate destruction; §6 quantum setup, positive cascade, inhomogeneous tail,
reflection break, bounded sweep + `μ₃=0` validation. A successful run ends
`ALL ASTAR BAND3 CHECKS PASSED` (48 checks).
