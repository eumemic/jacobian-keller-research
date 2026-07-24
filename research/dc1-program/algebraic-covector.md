# The algebraic-node covector calculus: Thm A′ over F(datum) via trace forms, and the two-block necklace of the W2 negative tail

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED — BAND-SCOPED — MIXED PROVED/BOUNDED**

This memo builds the **degree-free tool** the W2 endgame gap needs
([`residual-identity.md`](residual-identity.md) §6,
[`slope-forcing-degree-free.md`](slope-forcing-degree-free.md) §6): the covector on the
**algebraic `(a_2,b_1)` necklace** — the nodes are the roots of the cascade-solved data,
not integers, so the integer-node adjoint criterion
([`lambda-general-k.md`](lambda-general-k.md) Thm A′) does not directly apply. It does
five things:

1. **Extends Thm A′ to algebraic nodes (PROVED, node-free).** The moving-sum adjoint
   identity `lambda(S_n g) = (S_n^* lambda)(g)` holds verbatim for a **symbolic** node
   `rho` — nothing in the adjoint algebra ever used integrality. And a Galois-symmetric
   node-functional `lam = sum_{p(rho)=0} g(rho) ev_rho` over the roots of a datum
   polynomial `p` is a **trace form** `Tr_{F[E]/(p)}(g·)`, a rational function of the
   coefficients of `p` computed with **no root named** (polynomial after monic
   normalization; in general denominators are powers of the leading coefficient),
   closed under `S_n^*` — so the algebraic-necklace support conditions are
   degree-free-computable over the corresponding coefficient field via companion
   traces / resultants (lead 1).
2. **Derives the two-block structure of the negative tail (PROVED, degree-free).** The
   filler-linear parts of `Q_-1,Q_-2,Q_-3` are explicit **two-term** operators whose
   tops are the cascade-solved data `(b_1,b_0,b_-1)` and `(a_2,a_1,a_0)` (lead 2).
3. **Sharpens the necklace picture (CORRECTION).** In the **pure negative tail**
   `Q_-1..Q_-3` (branch B) the fixed tops `a_3,b_2` **do not appear at all**: the tail
   necklace is **entirely algebraic**. The `(a_3,b_2)` fixed necklace belongs to the
   `Q_0`/`Phi` joint operator ([`joint-covector.md`](joint-covector.md)), **not** the
   pure tail. Integer nodes enter the tail only through the **membership windows** `{0,1}`
   and `{0,1,2}`.
4. **Reads the `d=3` certificate through the calculus (BOUNDED, exact).** The depth-3
   cokernel is 16-dimensional; the W-kill is reproduced; `a_2(0)` is not forced; `W`
   enters the consistency conditions; `Q_-1` alone does not force `W` (lead 3).
5. **Tests the fixed trace-form recipe and localizes the obstruction (OPEN, lead 4/5).**
   On a rung where the two terms use independently variable filler evaluations, a single
   algebraic node would need a shared root of `top(E)` and `top(E-3)`, generically none;
   the sought covector must therefore **couple** terms and rungs across the **varying**
   tops `(a_2,a_1,a_0)/(b_1,b_0,b_-1)`. That coupling is the exact residual gap: a fixed
   finite trace-form recipe producing a unit multiple of `W` at **every** `d` is **not**
   obtained. At `d=4`, full filler-column rank and random mod-`p` samples give supporting
   evidence only; they do **not** construct a W-forcing covector or prove one exists.

W2 datum, gauge `b_3=0`, quantum band-3 conventions
(`Q_m=sum_(k+l=m)[b_l^[k]a_k - a_k^[l]b_l]`, `f^[n](E)=f(E+n)`,
membership `(E)_j=E(E-1)...(E-j+1) | a_-j,b_-j`):

```text
a_3 = E(E+2)(E+4)   (roots {0,-2,-4}),   b_2 = E(E+3),   b_3 = 0.
branch B:  a_-3=0,  b_-3=(E)_3 C;   branch-B fillers  a_-2=(E)_2 V,  b_-3=(E)_3 C.
```

Exact certificate: [`verify_algebraic_covector.py`](verify_algebraic_covector.py). Runtime is
environment-dependent; `HEAVY=1` adds the `d=4` linear-route mod-`p` leg and the optional
depth-3 `msolve` kill. The final summary distinguishes a run in which all checks passed from
one in which all executed checks passed but optional checks were skipped. Every load-bearing
upstream fact — the crossed-product engine, `Q_0=(T-1)G`, the slope gate, the both-ends Lemma
P, the factorization `R(1)=a_2(0)·W` — is **re-derived in file**.

## 0. Headline

> **THE TOOL IS BUILT; THE RECIPE REMAINS OPEN, WITH THE OBSTRUCTION LOCALIZED.**
> The adjoint criterion Thm A′ extends **node-free** to algebraic nodes, and
> Galois-symmetric annihilators descend to `F(datum)` as **trace forms** computable
> without naming a root (PROVED, machine-checked on symbolic blocks). The depth-3 negative
> tail is two explicit **two-term** block operators on a **purely algebraic** necklace
> (PROVED, degree-free). But on the tested generic b-block rung a single algebraic
> node does not suffice; the sought annihilator must **couple** terms and rungs across
> the varying tops `(a_2,a_1,a_0)/(b_1,b_0,b_-1)`. A **fixed finite trace-form recipe**
> yielding a unit multiple of `W` at
> every degree is **not** obtained. The residual identity
> `W ∈ sqrt(cascade+Q_-1..Q_-3)` at arbitrary `d` is therefore still **open**; what is new
> is the degree-free calculus in which it lives and the precise coupling that blocks it.

## 1. Algebraic-node Thm A′ — the tool (PROVED, node-free)

Thm A′ ([`lambda-general-k.md`](lambda-general-k.md) §2) reads a filler-image annihilator
off the adjoint of the moving sum `S_n=1+T+...+T^{n-1}`. Its proof
`lambda(S_n g)=(S_n^* lambda)(g)` is an identity in the test polynomial `g` and the node
`rho` — **integrality is never used**. So it holds verbatim for a symbolic (algebraic)
node, machine-checked as an identity in generic degree-5 `g` at `n=2,3,4` (`verify §S2(i)`).

**Trace-form descent (the "compute without naming roots").** For a squarefree datum
polynomial `p` (symbolic coefficients) and any weight/test `h`,

```text
   sum_{p(rho)=0} h(rho)  =  Tr_{F[E]/(p)}(h)  =  trace of mult-by-h on the companion basis,
```

a **rational function in the coefficients of `p,h`**, with denominators coming from
monic normalization (and hence polynomial in the coefficients of `h` and a monic `p`),
computed with no root named (`verify §S2(ii)`: the companion trace has free symbols
`⊆ coeffs(p)∪coeffs(h)`, and equals the actual root sum on random squarefree
specializations). **Tier note (audit-demoted): the equality clause
"companion trace = root-sum functional" is machine-checked in-file only NUMERICALLY at
fixed degree (`deg p=3, deg h=2`, 4 random squarefree specializations) — the fact is the
standard trace theorem, but the in-file evidence tier is bounded, not a symbolic proof
object. SEPARABILITY is essential and must be carried as a hypothesis: for NON-squarefree
`p` the companion trace returns the WITH-MULTIPLICITY sum, which differs from the
distinct-node covector (audit witness: `p=(E-2)^2(E+3)`, `h=E^2+1` — trace `20` ≠
distinct-node sum). Any use on a degenerate (repeated-root) necklace must either verify
squarefreeness or switch to the jet/multiplicity-extended criterion of
[`broken-separation.md`](broken-separation.md).** Hence, for squarefree `p`, a
Galois-symmetric node-functional
`lam=sum_{p(rho)=0} g(rho) ev_rho` **descends to `F(datum)`** and is computable by a
companion trace / resultant. Trace forms are **closed under `S_n^*`**:
`Tr(S_n g)=sum_{i=0}^{n-1} Tr(g(·+i))` as a coefficient identity (`verify §S2(iii)`), so
the shift-window support conditions stay trace-computable degree-freely. This formalizes
lead 1: the adjoint machinery works over `F(datum)`-bar and the annihilator conditions
descend as trace forms.

## 2. The two blocks of the negative tail (PROVED, degree-free)

In `Q_m=sum_(k+l=m)[b_l^[k]a_k-a_k^[l]b_l]` the branch-B fillers are `a_-2` (`k=-2`) and
`b_-3` (`l=-3`). Extracting their (linear) contributions to each tail rung gives the
explicit **two-term operators** (`verify §S3`, generic data):

```text
   a-block(Q_m) = b_{m+2}(E-2) a_-2(E) - b_{m+2}(E) a_-2(E+m+2),     top  b_{m+2},
   b-block(Q_m) = b_-3(E+m+3) a_{m+3}(E) - a_{m+3}(E-3) b_-3(E),     top  a_{m+3}.
```

| rung | a-block top | b-block top |
|---|---|---|
| `Q_-1` | `b_1` | `a_2` |
| `Q_-2` | `b_0` | `a_1` |
| `Q_-3` | `b_-1` | `a_0` |

All six tops are **cascade-solved data** (algebraic). This is the level-incidence
structural key of [`residual-identity.md`](residual-identity.md) §1, now written out as the
operators the covector must annihilate.

**The pure negative tail is on a purely algebraic necklace (CORRECTION).** With symbolic
`a_3,b_2` the machine confirms **none of `Q_-1,Q_-2,Q_-3` contains an `a_3` or `b_2`
symbol** (`verify §S3`). The `(a_3,b_2)` fixed necklace `{-4,..,3}` of
[`joint-covector.md`](joint-covector.md) §2 lives on the `Q_0`/`Phi` block
(`Phi=K_3[b_-3]-H_2[a_-2]`, whose tops **are** `a_3,b_2`); it does **not** occur in the
pure tail `Q_-1..Q_-3`. This **sharpens/corrects**
[`residual-identity.md`](residual-identity.md) §3, which described the pure-tail
consistency covectors as having "fixed part on the `(a_3,b_2)` root necklace": that fixed
part belongs to the **joint** `[Phi;M_-1]` operator, not to `Q_-1..Q_-3`. The only integer
nodes in the pure tail are the **membership windows**.

**Membership-window covectors annihilate a block (the `ev_0` analogue).** For `rho` a root
of the filler membership — `{0,1}` for `a_-2=(E)_2 V`, `{0,1,2}` for `b_-3=(E)_3 C` — the
single-row covector `ev_rho` on `lambda_-3` kills that block's contribution
(`verify §S3`, `d=2,3,4`), because the falling-factorial membership factor vanishes. These
are the tail's silent point functionals.

## 3. The block adjoint criterion (PROVED, node-free)

Pairing a point covector with a block and using membership makes the annihilation an
explicit **coefficient-of-filler-value** condition. For a symbolic node `rho`
(`verify §S4`):

```text
   ev_rho( b-block(Q_m) ) = (rho+m+3)_3 C(rho+m+3) a_{m+3}(rho) - (rho)_3 C(rho) a_{m+3}(rho-3),
   ev_rho( a-block(Q_m) ) = (rho)_2 V(rho) b_{m+2}(rho-2) - (rho+m+2)_2 V(rho+m+2) b_{m+2}(rho).
```

Collecting the coefficient of each surviving filler value `C(sigma)` / `V(sigma)` gives the
annihilation equations; a covector `lambda=(lambda_-1,lambda_-2,lambda_-3)` annihilates a
block iff every such coefficient (off the membership window) vanishes. When the nodes range
over a root necklace `rho∈roots(p)`, each coefficient is a **trace form** in the datum
(§1) — the algebraic-necklace support conditions are degree-free-computable. This is the
algebraic-node version of the moving-sum criterion, realized on the actual two-term tail
blocks.

## 4. The `d=3` certificate through the calculus (BOUNDED, exact)

On the branch-B parametrized `d=3` cascade (9 free coordinates,
`R(1)=-(4/9)a_2(0)·am1_3`, so `W=-(4/9)am1_3`), the depth-3 tail is **24 scalar equations,
all linear in the 8 fillers**; the filler map `M` has **full column rank 8**, cokernel
dimension **16** (`verify §S5`). Read through §1–§3:

- **The W-kill is reproduced.** `am1_3 ∈ sqrt(cascade+Q_-1..Q_-5)` over `QQ` and mod
  `65003` (sympy exact radical certificate); `msolve '^'` gives the depth-3 refinement
  `am1_3 ∈ sqrt(cascade+Q_-1..Q_-3)` under `HEAVY`. So the tail forces `W=0`.
- **`a_2(0)` is not forced.** Explicit `cascade+tail` witness with `a_2(0)≠0`, `R(1)=0`
  (`verify §S5`) — the kill is the **factor** `W`, `a_2(0)` free (as in
  [`residual-identity.md`](residual-identity.md) §2).
- **The specialized covectors are honest at the tested datum point.** After the free
  cascade coordinates are specialized, all 16 numerical cokernel covectors satisfy
  `mu·M=0` exactly at that specialization (annihilate every filler column). This is not
  a symbolic identity in the unspecialized free cascade coordinates. The filler map `M` is **independent of `am1_3`**: `W`
  lives purely in the residual `N`. And `W` **enters** the consistency conditions (some
  cokernel covector pairs `N` with a nonzero `am1_3`-coefficient), so tail-solvability
  genuinely constrains `W`.
- **Depth-3 is needed.** `Q_-1` alone does **not** force `W`: an explicit `cascade+Q_-1`
  solvable point has `R(1)=16/9≠0` (`verify §S5`). The kill needs the coupled depth-3 tail.

This is the explicit structural form of the `d=3` certificate (lead 3): the consistency
covectors of the two-block operator, on the algebraic necklace, with the membership windows
silent.

## 5. The fixed recipe and the coupling obstruction (OPEN, localized — lead 4/5)

**Why a single algebraic node is insufficient for the tested generic b-block.** On a
rung where the two filler evaluations are independent, `ev_rho` on `lambda_m` kills the
first term iff `a_{m+3}(rho)=0` and the second iff `a_{m+3}(rho-3)=0`; it kills that
whole rung only at a **shared** root of `a_{m+3}(E)` and `a_{m+3}(E-3)`. The verifier
checks the generic `a_2` instance, where `gcd(a_2(E),a_2(E-3))=1` (`verify §S3`). This
argument is not a blanket theorem for every block and rung: coincident filler evaluations
can cancel, and lower tops require separate analysis. For the tested generic b-block, the
annihilator must **couple** its two terms and ultimately the rungs with different tops
`a_2,a_1,a_0` (and analogously on the a-block side). In the integer-node band-3 AP family the
analogous coupling telescoped because the tops were shifts of one polynomial
([`lambda-general-k.md`](lambda-general-k.md) Thm C); here the tops are **independent
datum**, and no fixed telescoping is available.

**The recipe test.** A fixed finite trace-form recipe — a fixed node-selection rule at the
block root-sets through fixed shift windows, with trace-form weights, producing a covector
whose `N`-pairing is a unit multiple of `W` at every `d` — is **not obtained**. What is
verified degree by degree:

- `d=1,2,3`: `R(1)=a_2(0)·W` exact on the parametrized cascade; forcing vacuous for `d≤2`
  (cascade alone kills the slope) (`verify §S6`).
- `d=3`: the covector kill (§4).
- `d=4` (`HEAVY`): at sampled cascade points the depth-3 tail is 30 rows linear in 10
  fillers and its filler map has **full column rank 10**; among the sampled points with
  `W≠0`, the tail was never solvable. This is finite mod-`p` evidence only. It proves
  neither `W∈sqrt(I)` over the parameter space nor the existence of a symbolic W-forcing
  cokernel covector.

At `d=3` the exact covector kill exists. At `d=4` the nonempty cokernel and samples only
motivate the same search; whether its algebraic-necklace part pairs the residual to a unit
multiple of `W` remains unproved. The coupled two-term structure across varying tops does
not, in these tests, reduce to one fixed trace-form rule. **That coupling is the named residual gap** — the exact obstruction between
the built calculus and the arbitrary-`d` residual identity.

## 6. Evidence ledger — proved / bounded / refuted / open

**Proved (node-free / degree-free, char 0, machine-checked identities):**
- Engine `Q_m=[D,X]_m` (`m∈[-6,6]`), `Q_0=(T-1)G` (generic degree-2 coefficients); slope
  gate; both-ends Lemma P `R(1)=a_1(0)b_-1(1)+a_2(0)b_-2(2)-a_-1(1)b_1(0)`,
  filler-independent; the factorization `R(1)=a_2(0)·W` (`§S0,§S1`).
- **Algebraic-node Thm A′:** the moving-sum adjoint identity `lambda(S_n g)=(S_n^* lambda)(g)`
  for a **symbolic** node `rho` (`n=2,3,4`, generic degree-5 `g`); the **trace-form
  descent** `sum_{p(rho)=0} h(rho)=Tr_{F[E]/(p)}(h)` — a rational function in
  `coeffs(p,h)` (polynomial after monic normalization), equal to the root sum with
  multiplicity and **no root named** (in-file symbolic scope `deg p=3, deg h=2`, plus
  random squarefree root-sum confirmation); trace forms **closed under `S_n^*`** (`§S2`).
- **Two-block structure:** the explicit two-term operators `a-block`, `b-block` with tops
  `(b_1,b_0,b_-1)`, `(a_2,a_1,a_0)` (level incidence, generic data); the negative-tail
  necklace is **entirely algebraic** — `a_3,b_2` occur in **no** `Q_-1,Q_-2,Q_-3`;
  membership-window covectors annihilate a block (`d=2,3,4`); on the tested generic
  `a_2` b-block rung with independent filler evaluations, a single algebraic node cannot
  (`gcd(a_2(E),a_2(E-3))=1`) (`§S3`).
- **Block adjoint criterion:** the symbolic-node coefficient-of-`C(sigma)`/`V(sigma)`
  identities, whence the algebraic-necklace support conditions are trace forms (`§S4`).

**Bounded-finite (exact scope stated):**
- `d=3`: at the tested exact datum specialization, the depth-3 tail has 24 filler-linear
  equations, filler map full column rank 8, and cokernel dimension **16**; all 16 resulting
  covectors annihilate the specialized filler columns exactly. `M` is independent of
  `am1_3`; `W` enters the specialized consistency conditions; `Q_-1` alone does not force
  `W` (witness `R(1)=16/9`) (`§S5`).
- `d=3` W-kill (control, reproduced): `am1_3 ∈ sqrt(cascade+Q_-1..Q_-5)` over `QQ`+`GF(65003)`
  (sympy exact); depth-3 `am1_3 ∈ sqrt(cascade+Q_-1..Q_-3)` (`msolve '^'`, `HEAVY`);
  `a_2(0)` not forced (explicit witness) (`§S5,§S6`).
- `R(1)=a_2(0)·W` exact on the parametrized cascade at `d=1,2,3`; forcing vacuous `d≤2`
  (`§S6`).
- `d=4` (`HEAVY`): at sampled mod-`p` cascade points, the depth-3 filler map has full
  column rank **10**, and no sampled `W≠0` point was tail-solvable. This is supporting
  evidence only, not a symbolic W-forcing covector or radical certificate (`§S6`).

**Refuted / corrected (sharpening):**
- The description ([`residual-identity.md`](residual-identity.md) §3) that the **pure-tail**
  `Q_-1..Q_-3` consistency covectors have "fixed part on the `(a_3,b_2)` root necklace" is
  corrected: `a_3,b_2` appear in **no** `Q_-1,Q_-2,Q_-3` (`§S3`). The `(a_3,b_2)` fixed
  necklace is the `Q_0`/`Phi` block's ([`joint-covector.md`](joint-covector.md) §2); the
  pure tail is on a **purely algebraic** necklace with only membership-window integer nodes.
  (Not a change to any forcing conclusion — the kill and its tiers are unchanged.)

**Open / not claimed:**
- **A fixed finite trace-form recipe** producing a covector whose `N`-pairing is a unit
  multiple of `W` at **every** `d` — blocked by the two-term coupling across the varying
  tops `(a_2,a_1,a_0)/(b_1,b_0,b_-1)` (§5). The residual identity
  `W ∈ sqrt(cascade+Q_-1..Q_-3)` at arbitrary `d` therefore remains **open**; the calculus
  and the localized obstruction are the contribution.
- Everything the parent memos leave open: the arbitrary-degree slope forcing
  ([`slope-forcing-degree-free.md`](slope-forcing-degree-free.md) §6,
  [`residual-identity.md`](residual-identity.md) §6); no Weyl pair; all of Band 3, DC1, JC2.

## 7. Verification

```sh
uv run --with sympy python research/dc1-program/verify_algebraic_covector.py
# HEAVY (d=4 linear route mod p + depth-3 msolve kill):
HEAVY=1 uv run --with sympy python research/dc1-program/verify_algebraic_covector.py
```

`S0` engine; `S1` slope gate + both-ends Lemma P + factorization; `S2` algebraic-node
Thm A′ (symbolic-node adjoint + trace-form descent + `S_n^*` closure); `S3` the two blocks
(two-term operators, purely-algebraic necklace, membership-window covectors, tested generic
coupling obstruction); `S4` the block adjoint criterion; `S5` the bounded `d=3` certificate
(cokernel 16 at a specialization, W-kill, `a_2(0)` free, depth-3 needed); `S6` bounded
checks (`d=1,2,3` exact identities at their stated caps; `d=4` sampling plus optional
`msolve`, `HEAVY`). Runtime and executed-check counts depend on the environment and enabled
optional tools. The verifier reports either that all checks passed with no skips, or that all
executed checks passed with explicit skips; it never folds skipped optional legs into an
unqualified all-passed banner. Sampling does not certify the broader mathematical conjectures.
