# W2 theory: reflection does not close it, the common-root-at-0 census, and the completeness of its point annihilators

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER REVIEWED — BAND-SCOPED**

This memo attacks the exceptional member **W2** — the `r=-4` top of the band-3
step-2 exotic arithmetic-progression (AP) family — from three flanks: the
**symmetry transport** (does the quantum Fourier reflection inherit W2's death
from its reflection partner `r'=-1`?), the **common-root-at-0 class** at general
band `k` (the counterexample habitat's census), and the **non-point
annihilators** (is W2's point-functional basis complete, or is there a hidden
functional that closes it?). The identities used locally are re-checked in this
memo's verifier; the full arbitrary-degree lambda-wave inputs remain sibling
certificates and are rerun separately.

W2 datum, gauge `b_3=0`, quantum band-3 conventions
(`Q_m = sum_{k+l=m}[b_l^{[k]} a_k - a_k^{[l]} b_l]`, `f^{[n]}(E)=f(E+n)`,
membership `(E)_j = E(E-1)...(E-j+1) | a_{-j}, b_{-j}`):

```text
a_3 = E(E+2)(E+4)   (roots {0,-2,-4}),     b_2 = E(E+3)   (roots {0,-3}),     b_3 = 0.
```

Exact certificate: [`verify_w2_theory.py`](verify_w2_theory.py) (ends
`ALL W2 THEORY CHECKS PASSED`, 48 checks). Base commit `d8189fc`; re-verified the
AP-lambda and lambda-general-k certificates, both pass unchanged.

## 0. Headline

> **Three independent flanks confirm W2's specialness and reduce the W2 filler
> quotient obstruction to the slope scalar `R(1)=1`, which is achievable.** The
> complete dual of that quotient is silent at slope one; this does not exclude
> every possible linear invariant of a larger system. The later bounded verdict
> [`w2-verdict.md`](w2-verdict.md) proves the encoded combined slope+tail system
> empty over `QQ` at raw cap `d=3` and reports/optionally reproduces the same at
> `d=4`; arbitrary degree remains open.
>
> 0. **The slope gate is open.** `E-R in Im Phi(W2) <=> R(1)=1`.
>    The positive cascade forces `R(1)=0` at raw degree `d<=2`. At `d=3`, the
>    Rabinowitsch check proves only that `R(1)` is **not forced to zero**
>    (`R(1) not in sqrt(I)`). Value-one achievability comes instead from the exact
>    one-parameter family in W2-DECISIVE:
>    `R(1)=-(8/9)am1_3`, so `am1_3=-9/8` gives `R(1)=1`. The sibling verifier
>    reconstructs and checks that family; this verifier independently repeats the
>    radical control.
> 1. **The tested reflection does not close W2.** The quantum Fourier automorphism
>    `phi` (`x->-partial`, `E->-E-1`) maps W2's **top-wall** completion problem to a
>    **bottom-wall** problem, NOT to
>    the (lambda-killed) `r'=-1` top-wall problem. The membership anchor stays
>    rigidly at `0`; the band-fixed `E->-E-1` that *would* transport W2 to `r'=-1`
>    is not an algebra automorphism and moves the anchor to `-1`. **W2's silence
>    survives reflection** — consistent with its independently-found silence.
> 2. **Census.** A canonical common-root hatch exists for every band `k>=3`: the
>    step-`(k-1)` / step-`k` AP pair, with a unique common root forced onto the
>    anchor by `gcd(k-1,k)=1`; band 3 gives W2. Uniqueness of one hatch per band is
>    proved only in the stated AP/bounded census scopes, not in general.
> 3. **Completeness.** `Im Phi(W2)` is the **principal squarefree ideal**
>    `(E(E-1)(E+1))`, so its cokernel is semisimple and its entire dual is point
>    functionals: `{ev_{-1}, ev_0, ev_1}` is **complete** (`dim = codim = 3`).
>    There is no additional annihilator in the complete dual of this filler
>    quotient. This closes the hidden-functional route **within that quotient**;
>    it is not a theorem about every possible linear invariant of the full system.

## 1. The pivot fact and finite regressions

`Im Phi(W2)` is the admissible two-filler image
`Phi(C,V) = K_3[(E)_3 C] - H_2[(E)_2 V]`,
`K_3[c]=sum_{j=0}^2 a_3^{[j-3]} c^{[j]}`, `H_2[v]=sum_{j=0}^1 b_2^{[j-2]} v^{[j]}`
(`quantum-exotic-cokernel.md`, `quantum-ap-lambda.md`). Set `D=E(E-1)(E+1)`.
The degree-free proof and verifier regressions establish:

- **`Im Phi(W2) = D·F[E]` exactly.** The arbitrary-degree proof is the triangular
  quotient argument in [`quantum-ap-filler-image.md`](quantum-ap-filler-image.md)
  §4: every basis filler is divisible by `D`; the quotient operator sends `E^n`
  to degree `n+1` with leading coefficient `-2`; and an explicit filler maps to
  `D`, supplying the remaining constant. This verifier rechecks divisibility and
  finite ranks through degree 10 as regressions of that written proof. Since `(D)`
  has codimension `deg D = 3`, `codim Im Phi(W2)=3`.
- **Annihilators of `(D)` are exactly `{ev_{-1}, ev_0, ev_1}`** (`p in (D)` iff
  `p(-1)=p(0)=p(1)=0`).
- **Slope reduction.** At `r=-4` both filler blocks vanish at `E=1`, so
  `R(1)=G(1)=Q_0(0)`; with the proved cascade constraint `R(-1)=-R(1)`, every
  annihilator `alpha ev_0 + beta ev_1 + gamma ev_{-1}` sends `E-R` to
  `(beta-gamma)(1-R(1))`. Hence

```text
E - R  in  Im Phi(W2)   <=>   R(1) = 1        (R(0)=0 automatic; R(-1)=-1 forced consistent).
```

This is the complete **filler-quotient gate**: central completion requires and,
within that quotient, is equivalent to slope `R(1)=1`. Section 4 (re-verifying the
sibling W2-DECISIVE flank) shows that the positive cascade and membership can
produce that slope at `d=3`. A DC1 candidate would additionally have to satisfy
the full negative tail and every remaining Weyl equation simultaneously.

## 2. Task 1 — reflection transport: the anchor theorem

**Premise.** The map `E -> -E-1` sends the AP top `{r,r+2,r+4}` to
`{-r-5,-r-3,-r-1}`, i.e. pairs AP members `r <-> -r-5` (fixed centre `-5/2`,
never an integer AP). W2 (`r=-4`) pairs with `r'=-1`. Since `r'=-1 != -4`, the
symbolic-`r` theorem (`quantum-ap-lambda.md`, `d8189fc`, re-run) **kills the
`r'=-1` top-wall system** at arbitrary degree (`lambda_{-1}(E)=r'+4=3 != 0`). If a
symmetry transported W2's problem onto `r'=-1`'s, W2 would inherit that death.

**It does not.** The transport breaks, and the break is exactly the membership
anchor. Three exact facts (`verify §1`):

**(a) The band-fixed `E->-E-1` is not an algebra automorphism, and it moves the
anchor to `-1`.** `E=x·partial` forces any automorphism with `phi(x)=x` and
`phi(E)=-E-1` to have `phi(partial)=-partial-x^{-1} not in A_1`. So `E->-E-1`
alone (band preserved) is illegal. Applied formally to a coefficient it sends the
anchor-0 falling factorial `(E)_3` (roots `{0,1,2}`) to `+-(E+1)(E+2)(E+3)` (roots
`{-1,-2,-3}`) — the anchor slides to `-1`. This is the "or do they?" the task
flags: **under the naive band-fixed map the anchor would move to `-1`, and that is
precisely why it cannot be the genuine symmetry.** The same slide is visible in
the sub: band-fixed `b_2(W2)=E(E+3)` reflects to `(E-2)(E+1)` (roots `{-1,2}`),
whereas the genuine `r'=-1` sub is `E(E-3)` (roots `{0,3}`) — off by exactly a
shift of `1` (`reflected b_2(W2)(E) = b_2(r'=-1)(E+1)`). The naive transport is
inconsistent at the sub before the top is even considered.

**(b) The genuine `phi` keeps the anchor at `0` but reverses the band.** `phi` is
the unique band-reversing symmetry (`phi^4=id`, `phi^2=`parity is band-preserving;
band-2 memo §2.2, commit `84978b9`). On a positive band,
`(phi F)_{-k} = (-1)^k (E)_k a_k(-E-1)` — the reflected coefficient is *forced* to
carry the anchor-0 factor `(E)_k`. For W2's top,

```text
(phi F)_{-3} = (E)_3 · a_3(r'=-1) = E(E-1)(E-2)·(E+1)(E-1)(E-3)      (degree 6).
```

The `r'=-1` top-roots `{-1,1,3}` **do** appear (correct `E->-E-1` root transport),
but they are **trapped at the BOTTOM (band `-3`), inflated by the anchor-0 twist
`(E)_3` into a degree-6 coefficient** — never a clean degree-3 top. `(phi F)_{-3}`
is divisible by `(E)_3` (anchor `0`) and **not** by `(E+1)(E+2)(E+3)` (anchor
`-1`): `phi` preserves anchor `0`.

**(c) Therefore the two natural reflection candidates do not map W2's top-problem
to `r'=-1`'s top-problem.** The genuine `phi` reverses the band, so
`{W2 top-completions} ≅_phi {phi F bottom-completions}` — W2's *own mirror*
problem (`phi^2` returns to W2 up to a scalar), which the sibling W2-TAIL
independently confirms is the reflected **bottom wall** with the same
`Phi_3(S)=1+S+S^2` necklace criterion (`w2-negative-tail.md`). The band-fixed
coefficient reflection is not an automorphism. **These reflection routes give no
leverage on W2.**

> **Proposition (the tested reflections do not close W2).** The band-fixed
> `E->-E-1` substitution that matches the top root-sets is not an automorphism and
> shifts the membership anchor to `-1` (mismatching `r'=-1`). The genuine Fourier
> `phi` fixes the anchor at `0` but reverses the band, sending W2's degree-3 top to
> a degree-6 bottom `(E)_3·(r'=-1 top)`. This proves failure of these two natural
> transports; it is not a classification of all automorphisms of `A_1`.

**Contrast with the cases where reflection *did* work.** Classical A\*-I
(`astar-band3.md` §3) and quantum band-2 O2/O3 (`quantum-a2-zero.md`) used `phi`
successfully — because there the *mirror* problem was **independently closed**
(by Theorem A, resp. M4q+completion). For W2 the mirror bottom-problem is exactly
as open as W2 itself (W2-TAIL: the bottom is *feasible*), so reflection only
restates the question. `astar-band3.md` §6 already flagged that the quantum
reflection "needs a falling-factorial twist"; this memo pins that twist as the
`(E)_3` anchor factor and shows it is fatal to the W2 transport.

## 3. Task 2 — the common-root-at-0 census (the counterexample habitat)

The band-`k` top wall `Q_{2k-1}=0`, `u(E+k)a(E)=a(E+k-1)u(E)` (`u:=b_{k-1}`), is
equivalent on the root-necklace (`sigma=T^{-1}`, root at `E=-j <-> sigma^j`) to

```text
S_k(sigma) δ(u) = S_{k-1}(sigma) δ(a),      S_r = 1+sigma+...+sigma^{r-1},
```

solved exactly by cofactors `g` with `δ(a)=S_k g`, `δ(u)=S_{k-1} g` **both
effective** (`band-k-weapons.md`, re-checked). The **common-root-at-0** condition
`rho*=0` (top and sub sharing a root at the anchor) is `[δ(a)]_0>0` and
`[δ(u)]_0>0`.

**The habitat criterion is a *unique* common root, not merely a root at 0.** If a
top has common roots at `0` *and* at some `p != 0`, the `p`-functional has
`lambda_p(E)=p != 0` and still obstructs `Q_0=1` (modulo the residual identity
`lambda_p(R)=0`). So the E-pairing degenerates completely — the genuine escape
locus — **iff the top has a unique common root**, translated onto the anchor.

**The canonical general-`k` escape hatch (PROVED, arbitrary `k`).**

```text
a_k = prod_{i=0}^{k-1} (E + i(k-1))   (step-(k-1) AP, k roots),
u   = prod_{j=0}^{k-2} (E + j k)      (step-k AP, k-1 roots).
```

`verify §2` proves for all `k`:

- **Wall-admissible.** `S_k δ(u)` and `S_{k-1} δ(a)` are *both* the multiset
  `{0,1,...,k(k-1)-1}` (base-`k` resp. base-`(k-1)` digit ranges collapse), so both
  equal `S_{k(k-1)}` — the wall holds for **every** `k` (checked `k=3..7`, proof is
  the base-representation identity).
- **Unique common root at the anchor.** `i(k-1)=jk` with `gcd(k-1,k)=1` forces
  `k | i`, so `i=0`: the ONLY common root is at `E=0`.
- **Exotic** (`a_k` non-consecutive) for `k>=3`.

Band 3 is `a_3=E(E+2)(E+4)`, `u=E(E+3)` = **W2 exactly**. Band 4 is
`{0,-3,-6,-9}` / `{0,-4,-8}`; band 5 is `{0,-4,...,-16}` / `{0,-5,...,-15}`; etc.

**Finite higher-band evidence.** `verify §2` samples the band-4 canonical hatch.
Through the stated truncation, the sampled images have gcd `E(E-1)` and the
sampled quotients span through degree 14; that gcd is squarefree and contains the
anchor root `0`. This is consistent with `Im Phi=(E(E-1))` and the slope gate
`R(1)=1`, but it is **not** an arbitrary-degree equality or gate theorem. W2's own
sampled gcd/rank behavior is reproduced by the same routine, while the exact W2
image equality rests on the separate written triangular proof. Thus the Band-4
calculation is finite evidence for a higher-`k` analogy, not a proof for Band 4 or
for every hatch.

**Band-3 uniqueness.** Within the AP family, W2 (`r=-4`) is the unique
common-root-at-0 member (proved; matches `lambda-general-k.md` §6.2). Within a
bounded non-AP census (support `<=6`, coeffs in `{-1,0,1}`) W2 is the **only**
unique-common-root band-3 hatch; the deg-6 (`g(1)=2`) tops all carry multiple
common roots. *(bounded evidence for non-AP tops.)*

## 4. Task 3 — the point annihilators are complete at W2

`lambda-general-k.md` §2 warned that for a **single** block point functionals can
undercount (W1: `codim Im L_K = 6` but only `4` finitely-supported annihilators).
At W2 this does **not** happen, for a structural reason (`verify §3`):

> **`Im Phi(W2)` is the principal ideal `(D)`, `D=E(E-1)(E+1)` squarefree.** Hence
> `F[E]/(D)` is semisimple (`≅ F^3` by CRT), and the *entire* dual of `F[E]/(D)`
> is spanned by the coordinate = point evaluations at the roots `{-1,0,1}`. The
> map `lambda -> (lambda(1),lambda(E),lambda(E^2))` is an isomorphism onto `F^3`,
> under which `ev_{-1},ev_0,ev_1` are the (invertible-Vandermonde, `det=2`) rows.
> Therefore `Ann(Im Phi(W2)) = span{ev_{-1},ev_0,ev_1}` **exactly**:
> `dim Ann = codim = 3`, and there is **no non-point (infinite-support)
> annihilator.**

The single block `Im L_K(W2)` is *not* an ideal (`E·g` escapes its span,
`verify §3`) — which is precisely why *it* can have non-point annihilators; only
the **sum** `Im Phi = L_K + L_H` collapses to the principal ideal `(D)`. **The
non-point route cannot close W2:** every annihilator is a combination of
`ev_{-1},ev_0,ev_1`, and their pairing with `E-R` is `(beta-gamma)(1-R(1))` — zero
iff `R(1)=1`. No hidden functional exists.

## 5. Synthesis — the escape is real; the habitat says where else to look

Combining the three flanks with the siblings and the re-verified slope:

- **W2-DECISIVE / this memo §0,§4:** the slope `R(1)` is forced `0` at `d<=2`.
  At `d=3`, this memo's independent Rabinowitsch check proves only that `R(1)` is
  not forced to zero. The sibling W2-DECISIVE certificate proves the stronger
  free-modulus and achievability statements by its explicit family
  `R(1)=-(8/9)am1_3`, including the specialization `am1_3=-9/8`.
- **W2-TAIL (`w2-negative-tail.md`):** the negative-tail subsystem `Q_-1..Q_-6`
  is a **proper (feasible) ideal** at the tested raw caps `d=1,2` when the moment
  equation is omitted. This separate feasibility does not establish compatibility
  with slope one. The reflected **bottom wall** carries the *tame consecutive
  cube* `(E)_3` (not an exotic AP) — the mirror image of this memo's Task 1.

Reflection does not transport the desired closure, the **complete dual of the W2
filler quotient** is silent at slope one, and the moment slope is achievable.
These statements do not classify every possible linear invariant of the full
cascade-plus-tail system. The later bounded computation
[`w2-verdict.md`](w2-verdict.md) proves the encoded conjunction empty at raw cap
`d=3` and reports/optionally reproduces an exact `d=4` result; arbitrary-degree
W2 remains open.

> **Corrected escape-locus conjecture (band-`k`, evidence bands 3–5).** For every
> band-`k` wall-admissible top, `coker Phi` carries a functional obstructing
> `Q_0=1` — **unless** the top is (a translate onto the anchor of) a
> **unique-common-root** exotic hatch. At W2, the E-pairing degenerates
> completely, `Im Phi` is a principal squarefree ideal whose roots include the
> anchor, and the entire **filler-quotient obstruction** reduces to the single
> achievable slope scalar `R(1)=1`. It is conjectural that every higher-band
> hatch has the same principal-image and one-slope behavior. This does not dispose
> of the negative tail. The hatches are sparse in the bounded censuses; one
> canonical family exists per band (the step-`(k-1)`/step-`k` AP; band 3 = W2),
> but general uniqueness is unproved.

This **sharpens** `lambda-general-k.md` §6.2 in two ways. (i) The escape is not
"common root at 0" but "**unique** common root at 0": multiple common roots leave
a live obstruction (a nonzero-`rho*` functional), so the habitat is far thinner
than a bare root-at-anchor condition would suggest. (ii) The escape is genuine, not
an artifact — the slope gate `R(1)=1` it reduces to is open **and now known
achievable** at band 3.

**Where else to look (bounded census guidance).** The step-`(k-1)`/step-`k`
construction supplies a canonical hatch for every `k>=4`. At Band 4, finite
sampling finds gcd `E(E-1)` and quotient spanning through degree 14, consistent
with the analogous gate `R(1)=1`; no arbitrary-degree image equality or gate is
proved. It is **conjectural** that there is only one hatch per band, that every
hatch reduces to one common gate, or that this family supplies a missing uniform
induction ingredient. Higher-degree and non-AP alternatives remain open.

**The live frontier (updated by the bounded verdict).** A DC1 candidate requires
the combined cascade, `Q_0=1`, negative tail, and membership simultaneously. The
later [`w2-verdict.md`](w2-verdict.md) proves that encoded conjunction is the unit
ideal at raw cap `d=3` on both branches, so no candidate materializes there. It
reports/optionally reproduces the same exact verdict at `d=4`, but that result is
not self-contained without a completed optional run. The explicit slope witness
supplies branch-B separate feasibility at `d=3`, while committed tail feasibility
controls are at `d=1,2`; do not read these as a d=3-both-branches separate-
feasibility theorem. Exact `QQ` raw cap `d=4` remains unverified here unless that
optional run is completed; arbitrary-degree W2 and other Band-3 tops remain open.

## 6. Claim disposition

**Proved (exact algebra, arbitrary degree / arbitrary `k` as stated):**
- Task 1 reflection proposition: the band-fixed `E->-E-1` substitution is not an
  automorphism and moves the anchor to `-1`; genuine Fourier `phi` fixes anchor
  `0`, reverses band, and sends the top to the degree-6 bottom
  `(E)_3·(r'=-1 top)`. These exact identities exclude those two natural
  transports, not every possible automorphism.
- Task 2 general-`k` escape hatch: step-`(k-1)`/step-`k` AP is wall-admissible
  (base-representation collapse, all `k`), has a unique common root at the anchor
  (`gcd(k-1,k)=1`, all `k`), and is exotic (`k>=3`). Band 3 = W2; W2 unique in the
  band-3 AP family.
- Task 3 completeness: `Im Phi(W2)=(E(E-1)(E+1))` principal squarefree; point
  basis `{ev_{-1},ev_0,ev_1}` complete (`dim=codim=3`); no non-point annihilator.
- Pivot fact: `Im Phi(W2)=D·F[E]`, `E-R in Im Phi <=> R(1)=1`.

**Re-verified sibling result (exact, this memo §0,§4):** the slope `R(1)` is
forced `0` at `d<=2`. At `d=3`, Rabinowitsch proves only `R(1) not in sqrt(I)`,
i.e. the slope is not identically zero on the cascade variety. The explicit
W2-DECISIVE family proves the stronger claim: `R(1)=-(8/9)am1_3`, so `R(1)=1`
is achievable. Thus W2 is not closed by the moment slope alone.

**Bounded / exact-finite evidence:** band-3 uniqueness of the escape hatch among
non-AP tops (census window `<=6`, coeffs `{-1,0,1}`); band-4,5 census (canonical
hatch found; higher-degree hatches not excluded); the hatch `Im Phi`
degenerations computed to finite generator degree; the `d=3` slope freedom
(finite raw degree, exact).

**Exceptional loci (explicit):** the reflection pairing `r<->-r-5` has centre
`-5/2` (no integer AP) and places W2 (`r=-4`) in the orbit `{-4,-1}`; only
`r=-4` escapes the AP lambda obstruction, while its reflection partner `r=-1` is
killed. The general-`k` hatch exists at every `k>=3` (needs only
`gcd(k-1,k)=1`).

**Open / not claimed:** combined feasibility above the self-contained exact raw
cap `d=3` in [`w2-verdict.md`](w2-verdict.md); exact `QQ` raw cap `d=4` unless the
optional msolve check is run (a reported mod-65003 `d=5` result is not a `QQ`
proof); arbitrary-degree W2; completeness of the
non-AP Band-3 census; uniqueness/classification of higher-band hatches or a
uniform gate; other Band-3 tops, DC1, and JC2. No Weyl pair or counterexample is
constructed.

## 7. Verification

```sh
uv run --with sympy python research/band3/verify_w2_theory.py
```

Exact SymPy, ends `ALL W2 THEORY CHECKS PASSED` (48 checks): §0 pivot fact
(`Im Phi=(D)`, slope reduction); §1 reflection transport (anchor theorem); §2
common-root-at-0 census (necklace criterion, general-`k` hatch with base-rep and
`gcd` proofs, band-3 = W2, finite Band-4 gcd/rank regression); §3 point-annihilator
completeness; §4 the `d=2` forced-zero control and the `d=3` Rabinowitsch result
that `R(1)` is not forced to zero. Value-one achievability is supplied by the
explicit family in `verify_w2_decisive.py`, not inferred from radical
nonmembership. Re-run `verify_quantum_ap_lambda.py` and
`verify_lambda_general_k.py` to confirm the lambda-wave inputs.
