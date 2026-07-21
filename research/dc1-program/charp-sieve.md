# The char-p sieve: mod-p screening of DC1 counterexample candidates

**INDEPENDENTLY DERIVED — EXACT SYMPY MACHINE-CHECKED — NOT PEER REVIEWED.**
This promotes the archived provisional Result **P2**
(`archive-import/provisional/dixmier-band-program/dc1-opening.md`) to a stated,
audited, and computationally backed screening theorem. It does **not** claim DC1,
JC2, or a counterexample. Backing script:
[`verify_charp_sieve.py`](verify_charp_sieve.py) (exact arithmetic; 140 `PASS`
lines; ends `ALL CHARP SIEVE CHECKS PASSED`). Sieve implementation:
[`sieve_dc1_candidate.py`](sieve_dc1_candidate.py).

> **Provenance discipline (hostile-audit standard).** Every load-bearing
> assertion below is tagged **[proved]** (self-contained proof here),
> **[computed]** (exact machine check on the stated instances — corroboration,
> never completeness), or **[cited]** (established in the literature; used as a
> black box, with the reference). The bridge's central-reduction step is
> **[cited]** for its general form and **[computed]** on the concrete pairs the
> sieve runs on. The sieve's exclusion power rests on the **[proved]**
> degree-2/tame-Kummer corollary.

---

## 0. What this is, in one paragraph

DC1 asks whether every endomorphism of the first Weyl algebra
`A_1 = C<x, d | [d,x]=1>` is an automorphism. An endomorphism is a pair `(X, D)`
with `[D, X] = 1`; it is automatically injective (`A_1` is simple), so only
surjectivity is at stake. The **Tsuchimoto/BKK bridge** sends a pair `(X, D)`
with algebraic coefficients, for almost every prime `p`, to a **Jacobian-1
polynomial self-map `psi_p` of the affine plane** `A^2_{F_p} = Spec Z`, where
`Z = F_p[x^p, d^p]` is the center of `A_1 (x) F_p`. The **geometric degree**
`d = deg psi_p` is bounded independently of `p`, and for `p > d` every
subextension of the induced function-field extension is **tame**, so the classical
cyclic (Kummer) obstruction applies to `psi_p`. The sharp consequence: **a Keller
center map of geometric degree 2 is impossible for `p > 2`**, and the monodromy
pair `(G, H)` must satisfy `H ->> G^ab` (no nontrivial abelian subcover). The
sieve computes `psi_p`, checks `Jac = 1`, computes `d`, and emits this verdict. It
**screens** candidates (can exclude a hypothetical counterexample) but proves
nothing about DC1 on its own.

Conventions match the band corpus: `A_1[x^{-1}] = (+)_k x^k C[E]`, `E = x d`,
`(x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E)`, membership
`x^{-j} c(E) in A_1 <=> E(E-1)...(E-j+1) | c(E)`,
`Q_m = sum_{i+l=m}[b_l(E+i)a_i(E) - a_i(E+l)b_l(E)] = delta_{m0}`. For the char-p
work the more convenient model is the **normal-ordered** presentation
`A_1 (x) F_p = (+) F_p . x^i d^j` with `d x = x d + 1`; the two are the same
algebra.

---

## 1. The center of `A_1 (x) F_p` [proved]

**Lemma 1 (center).** Over any field `k` of characteristic `p > 0`,
`Z := Z(A_1 (x) k) = k[x^p, d^p]`, and a normal-ordered monomial `x^i d^j` is
central **iff** `p | i` and `p | j`.

*Proof.* Write `w = sum c_{ij} x^i d^j`. From `d x = x d + 1` one gets the two
inner derivations
```
  [d, w] = sum c_{ij} * i * x^{i-1} d^j ,      [x, w] = - sum c_{ij} * j * x^i d^{j-1} ,
```
because `[d, x^i] = i x^{i-1}` and `[x, d^j] = -j d^{j-1}` (and `[d,d^j]=[x,x^i]=0`).
The monomials `x^{i-1} d^j` (resp. `x^i d^{j-1}`) are `k`-linearly independent, so
`[d,w]=0 <=> (p | i whenever c_{ij} != 0)` and `[x,w]=0 <=> (p | j whenever
c_{ij} != 0)`. An element is central iff it commutes with the generators `x, d`,
hence iff every occurring monomial has `p | i` and `p | j`. Conversely `x^p, d^p`
are central: `[d, x^p] = p x^{p-1} = 0` and `[x, d^p] = -p d^{p-1} = 0` in
characteristic `p`. ∎

The mechanism is a **falling-factorial collapse**: `(ad_x)^p (d^j) = (j)_p d^{j-p}`
with `(j)_p = j(j-1)...(j-p+1) ≡ 0 (mod p)` for every `j` (one of the `p`
consecutive factors is divisible by `p`), so `(ad_x)^p = 0`, i.e. `x^p` is
central. Same for `d^p`. **[proved]**, **[computed]** (`verify §1,§2`, incl.
`d^2 x^2 = x^2 d^2 + 4 x d + 2`).

**Structure [cited].** `A_1 (x) k` is a free `Z`-module of rank `p^2` with basis
`{x^i d^j : 0 <= i,j < p}`, and it is an **Azumaya algebra** of PI-degree `p` over
`Z ≅ k[u, v]` (`u = x^p`, `v = d^p`), i.e. `Spec Z = A^2_k` (Revoy 1973; Tsuchimoto
2005 §1; Belov-Kanel–Kontsevich 2007 §2). This identifies the center with the
**affine plane**.

---

## 2. The central-reduction bridge (Tsuchimoto/BKK), made explicit

### 2.1 The recipe

Let `(X, D)` be an endomorphism pair, `[D, X] = 1`, with coefficients in a number
ring `R` (`Z` in the examples). For a prime `p` not dividing any coefficient
denominator, reduce to `X_p, D_p in A_1 (x) F_p`. Then:

> **Bridge recipe.** `psi_p := (P, Q)` with `P = X_p^{\,p}` and `Q = D_p^{\,p}`,
> each read as an element of `Z = F_p[u, v]` via `u = x^p`, `v = d^p`. This is a
> polynomial self-map `psi_p : A^2_{F_p} -> A^2_{F_p}`.

For the recipe to be **defined**, `X_p^{\,p}` and `D_p^{\,p}` must be central. This
is the content of the bridge:

**Proposition 2 (central p-th powers).** If `[D, X] = 1` in `A_1 (x) F_p`, then
`X^p in Z` and `D^p in Z`.

*Partial proof + status.* Since `[D, X] = 1`, `ad_D` is a derivation with
`ad_D(X) = 1`, so `ad_D(X^p) = sum_{r=0}^{p-1} X^r (ad_D X) X^{p-1-r} = p X^{p-1}
= 0`; likewise `ad_X(D^p) = 0`. Also `ad_X(X^p) = 0` and `ad_D(D^p) = 0` trivially.
Hence `X^p` commutes with the subalgebra `<X, D>`. Moreover, because `x^p, d^p` are
central (Lemma 1), `ad_X` and `ad_D` are **`Z`-linear** derivations of the rank-`p^2`
Azumaya algebra `A_1 (x) F_p`; the operator `(ad_X)^p = ad_{X^p}` is its
**p-curvature**, and the assertion `X^p in Z` is exactly the vanishing of that
p-curvature. For `X = x` this vanishing is the falling-factorial collapse of §1;
for a pair arising as the reduction of a characteristic-0 endomorphism it is
**[cited]** — Tsuchimoto (2005), "Endomorphisms of Weyl algebra and p-curvatures,"
Osaka J. Math. 42, Prop. 2.3ff, and Belov-Kanel–Kontsevich (2007), Moscow Math. J.
7. We do **not** claim an original elementary proof for the general pair; we
**[computed]** the vanishing exactly for every tame pair the sieve runs
(identity, quadratic/cubic shears, rotation, the band-2 pairs, the composite) at
`p = 3, 5, 7, 11` (`verify §2`). Note it is genuinely a *theorem about
`[D,X]=1`*: for a non-Weyl pair the p-th power can be **non-central** — e.g.
`X = x`, `D = d + xd` (so `[D,X] = 1 + x`) has `D^p` non-central mod `5, 7`
(`verify §4`), which is exactly the failure the sieve traps. ∎(partial)

### 2.2 Why `psi_p` is Jacobian-1

**Proposition 3 (Keller).** `Jac(psi_p) = P_u Q_v - P_v Q_u = 1` in `F_p[u, v]`.

*Status.* The center `Z = F_p[u, v]` carries a canonical Poisson bracket with
`{u, v} = 1`, obtained by reducing `(1/p)[·,·]` from the flat `Z`-lift `A_1` over
`Z_{(p)}`; the p-th-power map is a **Poisson map** (Belov-Kanel–Kontsevich 2007),
so from `[X, D] = -1` one gets `{P, Q} = {X^p, D^p} = 1`, and
`Jac(psi_p) = {P, Q} = 1`. We use the *result* `Jac = 1` and **[computed]** it
exactly on every tame pair at `p = 3, 5, 7, 11` (`verify §3`); the sieve
**re-checks it exactly for each candidate and each prime** — a violation would
flag a computational inconsistency (or a non-genuine pair). **[cited]** for the
structural reason, **[computed]** per candidate.

### 2.3 The geometric degree is bounded independently of `p` [proved]

**Lemma 4 (degree bound).** If `X, D` have `(x, d)`-bidegree `<= delta`
(max of `i + j` over occurring monomials), then `deg_{u,v} psi_p <= delta` for
every good `p`, and the geometric degree `d = [F_p(u,v) : F_p(P,Q)]` satisfies
`d <= (deg P)(deg Q) <= delta^2`, independent of `p`.

*Proof.* `X^p` has `(x,d)`-degree `<= p*delta`; being central, every occurring
exponent is a multiple of `p`, so in `(u, v) = (x^p, d^p)` its degree is `<=
delta`. The geometric degree of a dominant plane map is `<=` the product of the
`(u,v)`-degrees of its components (Bézout). Both bounds are `p`-free. ∎
**[computed]** (`verify §7`, `p = 3, 5, 7, 11, 13`): every tame candidate has
`deg_{u,v} psi_p <= bideg(X, D)` at all sampled primes, with `d` constant.

### 2.4 Worked examples (the recipe in action) [computed]

All exact (`verify §3`). Balanced integer coefficients; `u = x^p, v = d^p`.

| candidate | `(X, D)` | `psi_p` (`p >= 5`) | `Jac` | geom. deg |
|---|---|---|---|---|
| identity | `(x, d)` | `u |-> u,  v |-> v` | `1` | `1` |
| shear | `(x, d + x^2)` | `u |-> u,  v |-> u^2 + v` | `1` | `1` |
| rotation | `(d, -x)` | `u |-> v,  v |-> -u` | `1` | `1` |
| band-2 | `(x + d^2, d)` | `u |-> u + v^2,  v |-> v` | `1` | `1` |
| composite | `(d, d^2 - x)` | `u |-> v,  v |-> v^2 - u` | `1` | `1` |

**Small-prime corrections are real and harmless to the degree.** The shear gives
`v |-> u^2 + v - 1` at `p = 3` (an Artin–Schreier correction term appears once
`2 = p - 1`), versus `u^2 + v` for `p >= 5`; the map stays a geometric-degree-1
automorphism. This is exactly the "screen at good primes, watch small primes"
discipline. `verify §3` regression-pins both forms.

---

## 3. The tame cyclic obstruction (P2's core), stated and audited

### 3.1 The squarefree-pullback lemma in characteristic `p`

**Lemma 5 (squarefree pullback; char-free over perfect fields) [cited, note's
Lemma 5.1].** Let `k` be perfect and `psi : A^2_k -> A^2_k` a dominant map with
`Jac(psi) in k^*` (Keller), hence étale on the locus where it is quasi-finite; in
characteristic `p` this locus is tame exactly where the local extension degrees
are prime to `p`. Then for a reduced (squarefree) principal divisor `div(f)`
downstairs, the pullback `div(f ∘ psi)` is reduced.

*Why.* `A^2` is factorial with unit group `k^*`, so Weil divisors are principal;
an étale (unramified) map does not raise multiplicities, so a squarefree divisor
pulls back squarefree. The only failure is **wild** ramification, which requires a
local degree divisible by `p`. **[cited]**; the char-free statement over perfect
fields is the note's Lemma 5.1, used verbatim.

### 3.2 The tame Kummer step [proved, modulo Lemma 5]

**Proposition 6 (no tame cyclic subcover).** Over `k = F_p-bar`, let `psi` be a
Keller plane map of geometric degree `d`, and let `M / k(psi)` be a subextension
of `k(u,v) / k(psi)` that is **cyclic of degree `n` with `p ∤ n`** (tame). Then
`M = k(psi)`; i.e. there is **no nontrivial tame cyclic subcover**.

*Proof.* `k` contains the `n`-th roots of unity (`p ∤ n`), so by Kummer theory
`M = k(psi)(f^{1/n})` for some `f in k(psi)^*`. Pull `f` up to `k(u,v) =
k[u,v]_{(0)}`. The cover `M` is unramified over `A^2` away from `div(f)`; by
Lemma 5 the ramification divisor of `f^{1/n}` upstairs is reduced, and by
Zariski–Nagata purity ramification occurs only in codimension 1. But
`k[u,v]` is a UFD with units `k^*`, and a tame degree-`n` cover of `A^2_{F_p-bar}`
is classified by `Hom(pi_1^{t}(A^2), Z/n)`; the **tame abelianized fundamental
group of the affine plane is trivial** (`pi_1^{t,ab}(A^2_{F_p-bar}) = 0` — the
tame analogue of the simple-connectivity of `A^2_C`; cf. SGA1/Grothendieck–Murre
tame-fundamental-group theory). Hence `f = c * g^n` with `c in k^*, g in k(psi)^*`,
so `f^{1/n} in k(psi)` and `M = k(psi)`. ∎

The single wild loophole is `p | n`: the **Artin–Schreier** covers
`v |-> v^p - v` (or `v - v^p`) are degree-`p`, unramified over `A^1`, tame theory
does not see them, and they are precisely the classical char-p pathology
(`t |-> t - t^p`, derivative `1`, non-injective). The sieve exhibits this: the map
`(u, v - v^p)` has `Jac = 1` and geometric degree `p` — a Keller self-map of the
plane that is **not** an automorphism (`verify §5`). So the theorem and the known
counterexamples **partition cleanly along tame/wild**, as P2 asserts. **[computed]**
consistency check.

### 3.3 The obstruction and its consequences

Let `psi_p` be the bridge image of a candidate, geometric degree `d`, monodromy
pair `(G, H)`: `G = Gal(N / k(psi_p))` for `N` the Galois closure of
`k(u,v)/k(psi_p)`, `H = Gal(N / k(u,v))` the point stabilizer, `[G : H] = d`. For
`p > d`, all subextension degrees are `< p`, hence prime to `p`, hence tame, so
Proposition 6 applies to every cyclic subextension. Translating "no nontrivial
tame cyclic subcover" into group theory (intermediate fields `<->` subgroups
between `H` and `G`; abelian-over-base subcovers `<->` normal subgroups
`N' ⊇ [G,G]` with `H ⊆ N' ⊊ G`):

**Theorem 7 (tame monodromy obstruction) [proved from Prop. 6].** For `p > d`:

- **(C1) `d != 2`.** A degree-2 extension is automatically Galois and cyclic
  (`C_2`), tame for `p > 2`; Proposition 6 forbids it. So `deg psi_p != 2` for
  `p > 2`. *(This mirrors the classical plane-Keller fact that quadratic
  extensions are Galois and excluded — cf. `research/jc2-frontier-audit.md`
  item 3.)*
- **(C2) `H . [G, G] = G`,** equivalently the composite `H ↪ G ->> G^ab` is
  **surjective**: `H ->> G^ab`. There is **no nontrivial abelian subcover** of
  `k(u,v)/k(psi_p)`.
- **(C3) No prime-degree Galois part.** For any prime `ell < p`, `psi_p` has no
  degree-`ell` **Galois (cyclic)** subcover: a `C_ell`-cover would give
  `H = 1 ⊊ G = C_ell` with `H ↛ G^ab`, contradicting (C2). Non-Galois prime
  degrees are *not* excluded by monodromy alone — e.g. `d = 3` with `G = S_3`
  passes, since `S_3^ab = C_2` and the index-3 subgroup `H = <transposition>`
  surjects onto `C_2`. (Degree 3 is killed by *other* means classically —
  Orevkov — not by this obstruction.)
- **(C4) Compatible family.** `d` and the monodromy type are constant for almost
  all `p` (Lemma 4 + spreading-out), so a genuine counterexample would give a
  **`p`-uniform** family of tame Keller maps of one fixed degree `d >= 3` with
  perfect-modulo-`H` monodromy — a strong, checkable rigidity.
- **(C5) `G` has no cyclic subquotient exposed by `H`.** In particular `psi_p`
  cannot be a composite of an automorphism with any tame cyclic cover; the
  "abelian layers" a solvable monodromy would provide are all forbidden by (C2).

**Consequence for DC1 (the screen).** If a candidate `(X, D)` has `[D,X] = 1` and
its bridge image `psi_p` has geometric degree `d`, then:
- `d = 1` at good primes  ⟺  `psi_p` is an automorphism, **consistent with**
  `(X,D)` being an automorphism of `A_1` (Tsuchimoto/BKK: `phi in Aut(A_1)` iff
  `psi_p in Aut(A^2)` for almost all `p`). The screen does **not** fire; this is
  **not** a proof of automorphy.
- `d = 2` at some `p > 2`: **EXCLUDED** — no such char-0 endomorphism exists
  (Theorem 7 C1). The candidate is spurious.
- `d >= 3` at good primes: not excluded by the degree screen; must additionally
  meet (C2)/(C3). Full monodromy is *not* mechanized here (see §5).

---

## 4. The sieve: what it computes

[`sieve_dc1_candidate.py`](sieve_dc1_candidate.py) takes a pair `(X, D)` over `Q`
(normal-ordered dicts `{(i,j): coeff}`) and a list of primes and returns, per
prime:

1. **precondition** `[D, X] = 1` over `Q` (exact; rejects loudly otherwise);
2. **bad-prime skip** if `p` divides a coefficient denominator;
3. **bridge** `psi_p = (P, Q)` via `X^p, D^p` reduced into `Z = F_p[u,v]`; raises
   `BridgeFailure` if a p-th power is **non-central** (traps non-Weyl inputs);
4. **`Jac(psi_p)`** reduced mod `p`, with the exact `Jac == 1` check;
5. **geometric degree** `d` via elimination resultants `Res_v(P-a, Q-b)` and
   `Res_u(P-a, Q-b)`, coefficients reduced mod `p`, degrees cross-checked;
6. **verdict** from Theorem 7: `automorphism` (`d=1`), `tame` + degree-2
   **EXCLUSION** (`d=2, p>2`), `possibly-wild` (`p <= d`, tame obstruction does not
   apply), or `tame, not excluded by degree screen` (`d >= 3, p > d`).

The aggregate screen flags **EXCLUDED** if any tame prime excludes; otherwise it
reports the sampled geometric degrees. Example run
(`uv run --with sympy python sieve_dc1_candidate.py`) screens the tame library
(all `d = 1`, `Jac = 1`, not excluded) and the broken pairs (rejected).

### Validation summary (`verify_charp_sieve.py`, 140 exact checks)

- **Weyl soundness**: `d x = x d + 1`, associativity, `d^2 x^2 = x^2 d^2 + 4 x d + 2`.
- **Center theorem**: central `<=>` exponents divisible by `p`; `X^p, D^p` central
  for every tame pair at `p = 3,5,7,11` (bridge defined).
- **Bridge**: `Jac = 1` and `d = 1` for identity, shear, rotation, band-2,
  composite at `p = 3,5,7,11`; regression on the `p = 3` Artin–Schreier correction.
- **Broken pairs fail loudly**: `[D,X]=2` and `[D,X]=1+x` rejected at the gate;
  a non-Weyl pair makes `D^p` non-central, raising `BridgeFailure`.
- **Degree engine**: `identity->1`, `(u^2,v)->2`, `(u^3,v)->3`, Artin–Schreier
  `(u, v - v^p)->p` (with `Jac = 1`).
- **Verdict truth table**: `d=1` automorphism; `d=2, p>2` EXCLUDED; `d=2, p=2` and
  `d=p, p` possibly-wild; `d=3, p=7` tame-not-excluded; `d=None` inconclusive.
- **Degree bound**: `deg_{u,v} psi_p <= bideg(X,D)` at `p = 3,5,7,11,13`.
- **End-to-end**: tame library not excluded; `half-shear` skips `p=2`; the
  exclusion path fires on a synthetic geometric-degree-2 map.

---

## 5. What the sieve can and cannot exclude

**CAN exclude [proved/mechanized]:**
- Any candidate whose bridge image has **geometric degree 2** at a prime `p > 2`
  (Theorem 7 C1) — impossible, candidate is not a genuine endomorphism.
- Any **non-endomorphism** (`[D,X] != 1`) — rejected at the precondition.
- Any input that is **not a Weyl pair mod `p`** — trapped as a non-central p-th
  power (`BridgeFailure`).
- More generally, any candidate whose center map, at a tame prime, exhibits a
  **cyclic (Galois) subcover of degree prime to `p`** — but see the caveat below on
  monodromy computation.

**CANNOT exclude / does NOT prove [explicitly disclaimed]:**
- **DC1 itself.** A candidate that reduces to an automorphism at every sampled
  prime (`d = 1`) is *consistent with* being an automorphism but is **not proven**
  to be one. The correspondence "`phi in Aut(A_1)` iff `psi_p in Aut(A^2)` for
  almost all `p`" is the deep BKK statement; a *finite* prime sweep can never
  certify "almost all `p`."
- **Wild-degree candidates at their own prime.** When `p <= d` the covering degree
  is not prime to `p`; an Artin–Schreier subcover may exist and the tame
  obstruction is silent. The sieve reports `possibly-wild` and defers to a larger
  prime. Since `d` is bounded independent of `p` (Lemma 4), any fixed candidate is
  eventually screened at every `p > d` — but a *single* prime `p = d` proves
  nothing.
- **Full monodromy `(G, H)`.** The sieve mechanizes the **degree** obstruction
  (C1) sharply and reports geometric degree exactly. It does **not** compute the
  Galois group `G` of the center map's function-field extension; conditions (C2)/
  (C3) beyond degree 2 are stated but not decided by the current code. A genuine
  degree-`>= 3` candidate would need a monodromy computation (a Galois-group /
  fundamental-group calculation) that is out of scope here. This is the honest
  boundary of the screen.
- **A single anomalous prime.** An exclusion should be confirmed **stable across
  several tame primes**. A lone degree-2 reading at one prime could reflect bad
  reduction; genuine exclusions recur at all `p > d`. The sieve prints every
  prime so stability is visible.

**Relation to the campaign.** The moment-unit principle
(`research/band3/band-k-weapons.md` W4q: `Q_0 = (T-1)G`, membership `G(0)=0`, so
`Q_0 = 1 <=> G = E`) is the *quantum-side* obstruction — it lives in
`A_1[x^{-1}]` and constrains the ladder structure of a putative counterexample
directly. The char-p sieve is the *arithmetic-side* obstruction — it constrains
the counterexample's **shadow on the center** `A^2_{F_p}`. They are independent
screens on the same object: a DC1 counterexample must survive **both** the
band-`k` rigidity cascade (bands 1–2 proved, `84978b9`/`b9f9cf3`; band 3 nearly
complete, `99fe6ee..9fa9f74`) **and** the tame monodromy obstruction on every
`psi_p`. The band cascade kills structure from inside `A_1`; the sieve kills it
from the arithmetic of its mod-`p` reductions. Neither alone is DC1.

---

## 6. Rigor ledger

**[proved] here (self-contained):**
- Lemma 1 (center `= F_p[x^p,d^p]`; central `<=>` exponents `≡ 0 mod p`).
- Lemma 4 (geometric degree of `psi_p` bounded by `bideg(X,D)`, `p`-free).
- Proposition 6 & Theorem 7 C1–C3 (tame Kummer step and its degree/monodromy
  consequences), **modulo** the cited Lemma 5 and `pi_1^{t,ab}(A^2)=0`.

**[cited] (literature black boxes):**
- Azumaya structure / PI-degree `p` of `A_1 (x) F_p` (Revoy 1973; Tsuchimoto
  2005; BKK 2007).
- Proposition 2 general form (central p-th powers / p-curvature vanishing) —
  Tsuchimoto 2005; BKK 2007.
- Proposition 3 general form (`psi_p` Poisson, hence `Jac = 1`) — BKK 2007.
- Lemma 5 (squarefree pullback, char-free over perfect fields) — the note's
  Lemma 5.1; `pi_1^{t,ab}(A^2_{F_p-bar}) = 0` — SGA1 tame π_1 theory.
- `phi in Aut(A_1) <=> psi_p in Aut(A^2)` for almost all `p` — BKK 2007.

**[computed] (exact SymPy, corroboration only — never completeness):**
- Proposition 2 & 3 on the concrete tame pairs at `p = 3,5,7,11`.
- The full worked-example table; the degree-bound sweep `p <= 13`.
- The verdict truth table and the Artin–Schreier tame/wild consistency check.
- Every identity in `verify_charp_sieve.py` (140 `PASS`).

**NOT claimed:** DC1; JC2; a counterexample; that any specific candidate is an
automorphism; a monodromy computation beyond geometric degree; completeness of any
prime sweep. The sieve is a **necessary-condition screen** on hypothetical
counterexamples.

---

## 7. Reproduce

```sh
uv run --with sympy python research/dc1-program/verify_charp_sieve.py
# -> 140 PASS lines, ends: ALL CHARP SIEVE CHECKS PASSED

uv run --with sympy python research/dc1-program/sieve_dc1_candidate.py
# -> screens the tame library + broken pairs (self-demo)
```

Screen your own candidate:

```python
from sieve_dc1_candidate import sieve_candidate, elt
# X = x + d^2,  D = d   (normal-ordered: coeff, x-exponent, d-exponent)
X = elt((1, 1, 0), (1, 0, 2))
D = elt((1, 0, 1))
sieve_candidate(X, D, primes=(3, 5, 7, 11, 13, 17), name="my candidate")
```

**References.**
Y. Tsuchimoto, *Endomorphisms of Weyl algebra and p-curvatures*,
Osaka J. Math. **42**(2) (2005) 435–452 (MR2147727); precursor:
*Preliminaries on Dixmier conjecture*, Mem. Fac. Sci. Kochi Univ. Ser. A Math.
**24** (2003) 43–59.
A. Belov-Kanel (Kanel-Belov), M. Kontsevich, *The Jacobian conjecture is stably
equivalent to the Dixmier conjecture*, Moscow Math. J. **7**(2) (2007) 209–218,
doi:10.17323/1609-4514-2007-7-2-209-218 (arXiv:math/0512171).
P. Revoy, *Algèbres de Weyl en caractéristique p*, C. R. Acad. Sci. Paris Sér.
A-B **276** (1973) A225–A228 (center/Azumaya structure; cf. also Tsuchimoto 2005
§1). SGA1 (Grothendieck–Murre) for the tame fundamental group. Archived source:
`archive-import/provisional/dixmier-band-program/dc1-opening.md` Result P2 and the
note's Lemma 5.1 / Remark 5.7. Campaign cross-refs: `research/band3/band-k-weapons.md`
(W4q moment unit), band-2 assembly `84978b9` / `b9f9cf3`, band-3 `99fe6ee..9fa9f74`.
