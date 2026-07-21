# W2 theory: reflection does not close it, the common-root-at-0 census, and the completeness of its point annihilators

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER REVIEWED — BAND-SCOPED**

This memo attacks the exceptional member **W2** — the `r=-4` top of the band-3
step-2 exotic arithmetic-progression (AP) family — from three flanks: the
**symmetry transport** (does the quantum Fourier reflection inherit W2's death
from its reflection partner `r'=-1`?), the **common-root-at-0 class** at general
band `k` (the counterexample habitat's census), and the **non-point
annihilators** (is W2's point-functional basis complete, or is there a hidden
functional that closes it?). Every lambda-wave fact relied on is **re-verified in
file**, not merely cited.

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

> **Three independent flanks all confirm W2's specialness is genuine and reduce
> it to the single slope scalar `R(1)=1` — and that scalar is now shown
> ACHIEVABLE.** Reflection, the non-point functionals, and the moment slope each
> FAIL to close W2; what remains is the combined slope+tail feasibility.
>
> 0. **The slope gate is open (re-verified).** `E-R in Im Phi(W2) <=> R(1)=1`.
>    The positive cascade forces `R(1)=0` at raw degree `d<=2` but at `d=3`
>    `R(1)` becomes a **genuine free modulus** (Rabinowitsch radical-membership
>    test: `R(1) not in sqrt(I)`), so **`R(1)=1` is achievable**. This
>    independently corroborates the sibling W2-DECISIVE (whose verifier crashes on
>    a downstream slice-display bug); the decisive test and the `d=2` control both
>    reproduce here.
> 1. **Reflection does not close W2.** The quantum Fourier automorphism `phi`
>    (`x->-partial`, `E->-E-1`) is the *only* band-reversing symmetry of `A_1`; it
>    maps W2's **top-wall** completion problem to a **bottom-wall** problem, NOT to
>    the (lambda-killed) `r'=-1` top-wall problem. The membership anchor stays
>    rigidly at `0`; the band-fixed `E->-E-1` that *would* transport W2 to `r'=-1`
>    is not an algebra automorphism and moves the anchor to `-1`. **W2's silence
>    survives reflection** — consistent with its independently-found silence.
> 2. **Census.** The counterexample habitat (`rho*=0`) is sparse and structured:
>    per band `k>=3` there is exactly **one canonical minimal escape hatch**, the
>    **step-`(k-1)` / step-`k` AP** pair, with a **unique** common root forced onto
>    the anchor by `gcd(k-1,k)=1`. Band 3 is **exactly W2**. Generic exotics have
>    *several* common roots — a nonzero one keeps the obstruction alive — so they
>    are **not** hatches.
> 3. **Completeness.** `Im Phi(W2)` is the **principal squarefree ideal**
>    `(E(E-1)(E+1))`, so its cokernel is semisimple and its entire dual is point
>    functionals: `{ev_{-1}, ev_0, ev_1}` is **complete** (`dim = codim = 3`).
>    There is **no** non-point annihilator; the route "find a hidden functional
>    that pairs `E-R` nontrivially" is closed. W2 lives or dies by `R(1)=1`.

## 1. The pivot fact, re-verified (not cited)

`Im Phi(W2)` is the admissible two-filler image
`Phi(C,V) = K_3[(E)_3 C] - H_2[(E)_2 V]`,
`K_3[c]=sum_{j=0}^2 a_3^{[j-3]} c^{[j]}`, `H_2[v]=sum_{j=0}^1 b_2^{[j-2]} v^{[j]}`
(`quantum-exotic-cokernel.md`, `quantum-ap-lambda.md`). Set `D=E(E-1)(E+1)`.
`verify §0` re-derives, exactly:

- **`Im Phi(W2) = D·F[E]` exactly.** Every basis filler is divisible by `D`
  (inclusion; the arbitrary-degree inclusion is the committed `lambda_r`
  annihilation identity, re-run), and the images `Phi/D` span `F[E]` (surjectivity
  onto `(D)`, checked through degree 10 with a triangular leading-term pattern).
  Since `(D)` has codimension `deg D = 3`, `codim Im Phi(W2)=3`.
- **Annihilators of `(D)` are exactly `{ev_{-1}, ev_0, ev_1}`** (`p in (D)` iff
  `p(-1)=p(0)=p(1)=0`).
- **Slope reduction.** At `r=-4` both filler blocks vanish at `E=1`, so
  `R(1)=G(1)=Q_0(0)`; with the proved cascade constraint `R(-1)=-R(1)`, every
  annihilator `alpha ev_0 + beta ev_1 + gamma ev_{-1}` sends `E-R` to
  `(beta-gamma)(1-R(1))`. Hence

```text
E - R  in  Im Phi(W2)   <=>   R(1) = 1        (R(0)=0 automatic; R(-1)=-1 forced consistent).
```

This is the whole game: W2 becomes a DC1 counterexample candidate **iff the
positive cascade + membership can produce slope `R(1)=1`** — which §4 (re-verifying
the sibling W2-DECISIVE flank) shows **is achievable at `d=3`**.

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

**(c) Therefore no symmetry maps W2's top-problem to `r'=-1`'s top-problem.**
`phi` reverses the band, so `{W2 top-completions} ≅_phi {phi F bottom-completions}`
— W2's *own mirror* problem (`phi^2` returns to W2 up to a scalar), which the
sibling W2-TAIL independently confirms is the reflected **bottom wall** with the
*same* `Phi_3(S)=1+S+S^2` necklace criterion (`w2-negative-tail.md`). The only
band-preserving symmetries (scaling, `phi^2`, pair-exchange) cannot turn a
bottom-problem into a top-problem, and pair-exchange kills the band-3 structure
(`b_3=0` becomes the top). **Reflection gives no leverage on W2.**

> **Theorem (reflection does not close W2).** There is no automorphism of `A_1`
> carrying W2's top-wall completion problem (top `{0,-2,-4}`, sub `{0,-3}`,
> anchor-0 membership) to the `r'=-1` top-wall completion problem (top `{-1,1,3}`,
> sub `{0,3}`, anchor-0 membership). The band-fixed `E->-E-1` that matches the top
> root-sets is not an automorphism and shifts the membership anchor to `-1`
> (mismatching `r'=-1`); the genuine Fourier `phi` fixes the anchor at `0` but
> reverses the band, sending W2's degree-3 top to a degree-6 bottom `(E)_3·(r'=-1
> top)`. W2's exceptionality is **stable under reflection.**

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

**These hatches degenerate `Im Phi` the same way as W2** (`verify §2`, exact
finite computation): the band-4 hatch has `Im Phi = (E(E-1))` (principal,
squarefree, codim 2, anchor `0` a root), and `E-R in Im Phi <=> R(1)=1` — the
*same* slope gate. W2's own `Im Phi=(E(E-1)(E+1))` (codim 3) is reproduced by the
identical generic routine. So the general-`k` hatch is a bona-fide W2-analogue.

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

- **W2-DECISIVE / this memo §0,§4:** the slope `R(1)` is forced `0` at `d<=2` but
  is a **free modulus at `d=3`**, so `R(1)=1` is **achievable**. (My independent
  Rabinowitsch radical-membership re-verification reproduces the committed `d=2`
  control *and* the `d=3` freedom; the only defect in the sibling certificate is a
  cosmetic slice-display formula that makes its script exit early.)
- **W2-TAIL (`w2-negative-tail.md`):** the negative tail `Q_-1..Q_-6` is a
  **proper (feasible) ideal**; it imposes no obstruction beyond the slope, and the
  reflected **bottom wall** carries the *tame consecutive cube* `(E)_3` (not an
  exotic AP) — the mirror image of this memo's Task 1.

So **every** functional-method obstruction on W2 has now failed: reflection
(Task 1), the completed point-annihilator dual (Task 3), and the moment slope
(§0/§4). W2 is **not closed** by any of them.

> **Corrected escape-locus conjecture (band-`k`, evidence bands 3–5).** For every
> band-`k` wall-admissible top, `coker Phi` carries a functional obstructing
> `Q_0=1` — **unless** the top is (a translate onto the anchor of) a
> **unique-common-root** exotic hatch. At such a hatch the E-pairing degenerates
> completely, `Im Phi` collapses to a principal squarefree ideal whose roots
> include the anchor, and the entire obstruction reduces to a **single slope
> scalar `R(1)=1` — which is achievable** (proved at band 3). The hatches are
> sparse: one canonical minimal family per band (the step-`(k-1)`/step-`k` AP;
> band 3 = W2), plus at most finitely many higher-degree unique-common-root
> exotics per bounded window.

This **sharpens** `lambda-general-k.md` §6.2 in two ways. (i) The escape is not
"common root at 0" but "**unique** common root at 0": multiple common roots leave
a live obstruction (a nonzero-`rho*` functional), so the habitat is far thinner
than a bare root-at-anchor condition would suggest. (ii) The escape is genuine, not
an artifact — the slope gate `R(1)=1` it reduces to is open **and now known
achievable** at band 3.

**Where else to look (the census answer).** The band-`k` step-`(k-1)`/step-`k`
hatches (`k>=4`) are the only other minimal habitats, and each reduces to the
*identical* slope gate `R(1)=1` (band 4 verified: `Im Phi = (E(E-1))`,
`E-R in Im Phi <=> R(1)=1`). Because the reduction is uniform across the family,
the achievability shown at band 3 makes these the natural next counterexample
candidates — the DC1 search should climb the hatch tower.

**The live frontier (no counterexample yet).** A DC1 counterexample candidate
requires the **combined** system — positive cascade **and** `Q_0=1` (slope `=1`)
**and** the negative tail `Q_-1..Q_-6` **and** membership — satisfied
*simultaneously* at a common degree, then `[D,X]=1` verified. The slope is
achievable at `d=3` and the tail is feasible separately, but their **intersection**
(slope=1 `∧` tail-feasible) has not yet been assembled or `[D,X]`-verified by any
flank (the escalation task is in progress team-wide). **No full candidate pair
materialized in this run, so the escalation orders were not triggered and no
counterexample is claimed** — but the obstruction has been narrowed to exactly one
finite feasibility question.

## 6. Claim disposition

**Proved (exact algebra, arbitrary degree / arbitrary `k` as stated):**
- Task 1 reflection theorem: no automorphism transports W2's top-problem to the
  `r'=-1` top-problem; band-fixed `E->-E-1` not an automorphism and anchor->-1;
  genuine `phi` fixes anchor `0`, reverses band, top->degree-6 bottom `(E)_3·(r'=-1
  top)`. (All identities exact.)
- Task 2 general-`k` escape hatch: step-`(k-1)`/step-`k` AP is wall-admissible
  (base-representation collapse, all `k`), has a unique common root at the anchor
  (`gcd(k-1,k)=1`, all `k`), and is exotic (`k>=3`). Band 3 = W2; W2 unique in the
  band-3 AP family.
- Task 3 completeness: `Im Phi(W2)=(E(E-1)(E+1))` principal squarefree; point
  basis `{ev_{-1},ev_0,ev_1}` complete (`dim=codim=3`); no non-point annihilator.
- Pivot fact: `Im Phi(W2)=D·F[E]`, `E-R in Im Phi <=> R(1)=1`.

**Re-verified sibling result (exact, this memo §0,§4):** the slope `R(1)` is
forced `0` at `d<=2` and is a free modulus at `d=3` (`R(1) not in sqrt(I)` by
Rabinowitsch), so **`R(1)=1` is achievable** — W2 is not closed by the moment
slope. (Corrects the sibling's crashing certificate; the decisive test and its
`d=2` control both reproduce here.)

**Bounded / exact-finite evidence:** band-3 uniqueness of the escape hatch among
non-AP tops (census window `<=6`, coeffs `{-1,0,1}`); band-4,5 census (canonical
hatch found; higher-degree hatches not excluded); the hatch `Im Phi`
degenerations computed to finite generator degree; the `d=3` slope freedom
(finite raw degree, exact).

**Exceptional loci (explicit):** the reflection pairing `r<->-r-5` has centre
`-5/2` (no integer AP); the lambda-escape AP members are exactly `r in {-4,-1}`
(`W2` and its reflection partner, which *is* killed). The general-`k` hatch exists
at every `k>=3` (needs only `gcd(k-1,k)=1`).

**Open / not claimed:** the **combined** feasibility (positive cascade `∧`
slope `=1` `∧` negative tail `Q_-1..Q_-6` `∧` membership) at a common degree, and
the `[D,X]=1` verification of any resulting pair — the live escalation frontier;
whether `R(1)=1` persists as achievable together with a feasible tail; the
arbitrary-degree slope law; completeness of the non-AP band-3 census and the
band-`k` hatch enumeration beyond the stated windows; DC1; JC2. No Weyl pair or
counterexample is constructed.

## 7. Verification

```sh
uv run --with sympy python research/band3/verify_w2_theory.py
```

Exact SymPy, ends `ALL W2 THEORY CHECKS PASSED` (48 checks): §0 pivot fact
(`Im Phi=(D)`, slope reduction); §1 reflection transport (anchor theorem); §2
common-root-at-0 census (necklace criterion, general-`k` hatch with base-rep and
`gcd` proofs, band-3 = W2, band-4 degeneration); §3 point-annihilator
completeness; §4 the re-verified slope gate (`d=2` forced control, `d=3` free
modulus via Rabinowitsch radical membership). Re-run `verify_quantum_ap_lambda.py` and `verify_lambda_general_k.py`
to confirm the lambda-wave inputs (both pass at `d8189fc`).
