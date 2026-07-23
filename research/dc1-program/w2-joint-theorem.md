# The W2 joint obstruction: the slope+tail incompatibility, localized and reduced

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED — BAND-SCOPED**

This memo attacks architecture step 1: turn the bounded `d<=4` Gröbner kill of the
W2 combined system ([`../band3/w2-verdict.md`](../band3/w2-verdict.md)) into a
degree-free theorem. It does **not** fully close W2 at arbitrary degree, but it
does three things that materially change the problem:

1. **Corrects the target.** The stated goal — "on the variety {positive cascade,
   negative tail, membership} the slope `R(1)=G(1)` vanishes identically" — is
   **FALSE**. `R(1)` is a *nonvanishing free modulus* on the cascade+tail variety.
   The joint kill is **not** a slope-forcing statement; it is a **filler**
   obstruction — the fillers `Q_0=1` demands (to make `G=E`) clash with the tail.
2. **Localizes the kill (generically) and linearizes it.** For the **generic**
   slope-`1` datum, the obstruction is already present in **`{Q_0=1, Q_-1}`** — the
   moment-unit condition and the **single** first negative-tail equation — and on
   the free-filler branch it is a **linear (Fredholm) inconsistency** in the two
   fillers. (A special slope-`1` sub-locus survives `Q_-1`; the *uniform* kill is
   the full tail `Q_-1..Q_-5`, both branches.)
3. **Reduces arbitrary-degree W2 to one named lemma** — the *Joint Filler Covector
   Lemma* — and supplies a new degree-free structural tool (the *both-ends Lemma
   P*) plus strong generic + degree-free-in-filler evidence for it.

W2 datum, gauge `b_3=0`, quantum band-3 conventions
(`Q_m = sum_(k+l=m)[b_l^[k] a_k - a_k^[l] b_l]`, `f^[n](E)=f(E+n)`,
membership `(E)_j=E(E-1)...(E-j+1) | a_-j, b_-j`):

```text
a_3 = E(E+2)(E+4)   (roots {0,-2,-4}),     b_2 = E(E+3),     b_3 = 0.
```

Exact certificate: [`verify_w2_joint.py`](verify_w2_joint.py) (ends
`ALL W2 JOINT CHECKS PASSED`). Base commit `30d8c59`; every load-bearing upstream
fact (the crossed-product engine, `Im Phi = D·F[E]`, the slope gate) is
**re-verified in file**, not merely cited.

## 0. Setup

After the positive cascade `Q_4=Q_3=Q_2=Q_1=0` solves `b_1,b_0,b_-1,b_-2`, and the
bottom Wronskian `Q_-6=0` splits into the two exhaustive branches
(`../band3/w2-negative-tail.md`)

```text
branch A:  a_-3 = (E)_3 · am3  (generic),  b_-3 = mu_3 · a_-3;
branch B:  a_-3 = 0,                       b_-3 = (E)_3 · C   (free filler).
```

the two admissible **fillers** are `a_-2 = (E)_2 V` and `b_-3` (`=mu_3 a_-3` on A,
`=(E)_3 C` on B). Writing `G` for the moment potential (`Q_0=(T-1)G`), membership
gives `G(0)=0`, and `Q_0=1 <=> G=E <=> E-R = Phi(C,V) in Im Phi`, where
`R` is the non-filler residual and `Phi(C,V)=K_3[(E)_3 C]-H_2[(E)_2 V]`. The pivot
`Im Phi(W2)=D·F[E]`, `D=E(E-1)(E+1)`, gives the slope gate

```text
Q_0 = 1   <=>   R(1)=1  and  R(-1)=-1     (R(0)=0 automatic).
```

## 1. The both-ends Lemma P (new, degree-free)

Evaluating the band-`k` moment potential `G` (`../dc1-program/moment-unit-general-k.md`,
Lemma P) at `E=1` gives the top-anchored boundary pairing; at `E=-1` its **mirror**:

```text
G(1)  = sum_(i=1)^k [ a_i(0) b_-i(i)   - a_-i(i) b_i(0) ],
G(-1) = sum_(i=1)^k [ a_i(-1-i) b_-i(-1) - a_-i(-1) b_i(-1-i) ]   (E -> -E-1 reflection of G(1)).
```

**Both are proved here as arbitrary-degree identities** (`verify §S2`). At W2 the
**level-3 terms drop on both ends**, for a structural reason:

```text
a_3(0)=a_3(-4)=0,   b_2(0)=b_2(-3)=0     =>
  R(1)  = G(1)  = a_1(0)b_-1(1)  + a_2(0)b_-2(2)  - a_-1(1)b_1(0),
  R(-1) = G(-1) = a_1(-2)b_-1(-1) - a_-1(-1)b_1(-2) + a_2(-3)b_-2(-1).
```

So the W2 slope pair is a **3-term boundary pairing at each end**, and — decisively
— `R(1)` and `R(-1)` are **independent of the fillers** `a_-2, b_-3` (the level-3
`b_-3` and `a_-2`-carrying terms are exactly the ones that vanish). This is why the
slope alone cannot see the tail: the tail lives in the filler data, the slope does
not.

## 2. The correct mechanism (not slope-forcing)

> **Refutation (machine-checked).** `R(1)=G(1)` does **NOT** vanish on the
> cascade+tail variety. At `d=3` it is a **positive-dimensional free modulus**
> there: `R(1) not in sqrt(cascade+tail)` by a Rabinowitsch test — msolve returns a
> positive-dimensional parametrization over `QQ` and over `GF(p)` (`verify §S8`
> reproduces the mod-`p` witness; `verify §S7` exhibits slope-`1` positive data on
> the cascade). Hence the proposed target — "cascade ∧ tail ⇒ `R(1)=0`" — is
> **false**, and the joint kill is *not* a slope statement.

The reason: `R(1)` is filler-independent (§1), while the tail `Q_-1..Q_-5` and the
moment-unit condition `Q_0=1` both constrain the **fillers**. The incompatibility is
between the fillers that `Q_0=1` forces (to make `G=E`) and the fillers the tail
allows. Concretely (`verify §S3`), **`a_-2` does not enter the positive cascade** —
it is a pure filler — and on **branch B**

```text
Q_0-1, Q_-1, Q_-2, Q_-3, Q_-4   are LINEAR in the two fillers (a_-2, b_-3);
Q_-5                            is the sole BILINEAR one.
```

## 3. The generic localization and the Fredholm gap

> **Generic localization (bounded, both branches).** At the explicit slope-`1`
> positive datum, the ideal generated by `{Q_0=1, Q_-1}` **in the fillers** is
> already the **unit ideal** — on branch A *and* branch B — while `{Q_0=1}` alone
> is feasible. The remaining tail equations `Q_-2,...,Q_-5` are **not used**
> *for this datum* (`verify §S5`).
>
> **This localization is generic, not uniform.** The symbolic
> `{cascade + Q_0=1 + Q_-1}` at `d=3` is **not** the unit ideal on either branch
> (`verify §S8`, msolve): a special slope-`1` sub-locus **survives** `Q_-1`, and
> the **uniform** kill is the full tail — `{cascade + Q_0=1 + Q_-1..Q_-5}` is the
> unit ideal, both branches (`verify §S8`; `w2-verdict.md` at `d=3,4`).

So the picture is *stratified*: `Q_-1` kills the **generic** slope-`1` datum, and
the deeper tail equations pick off the residual sub-loci. The kill is a **filler**
statement — *no filler completion realizes `Q_0=1` and the tail simultaneously* —
not a slope statement.

On branch B the fillers enter `{Q_0=1, Q_-1..Q_-4}` **linearly**, so the generic
kill is a Fredholm (rank) statement. Writing the combined operator
`L : (a_-2, b_-3) -> (coeffs of E-R-Phi ,  coeffs of Q_-1)` and `rhs` the
cascade/const terms:

> **Fredholm gap (branch B, generic).** `rank[L | rhs] - rank[L] = 1`: the `rhs`
> lies **outside the column span of `L`** — no filler completion exists. The gap
> is **exactly 1 and stable across filler degrees `dV=4,5,6,7,8`** at the explicit
> datum (`verify §S6`), so **that datum admits no tail-clearing filler of any
> degree**. Across **6 distinct slope-`1` data** at `d=3`, `L` has full column
> rank and the gap is `1` at every point (`verify §S7`) — generic in the positive
> data (rank only drops under specialization, so gap `1` at these points forces
> gap `>= 1` at the generic point of the slope-`1` fiber). The `Q_-1`-surviving
> sub-locus is exactly where this generic full-rank fails; the full tail supplies
> the extra rows that close it.

## 4. The named residual: the Joint Filler Covector Lemma

The arbitrary-degree W2 joint theorem is now reduced to exactly one statement.

> **Joint Filler Covector Lemma (branch B).** For every W2 slope-`1` positive
> datum, at every degree, the combined **linear** filler operator
> `L = [Phi ; M_{-1} ; M_{-2} ; M_{-3} ; M_{-4}]` acting on `(a_-2, b_-3)` has its
> `rhs` (the cascade/const parts of `E-R` and `Q_-1..Q_-4`) **outside its column
> span**: there is a covector `lam` with `lam·L = 0` and `lam·rhs != 0`. (`M_{-j}`
> is the filler-linear part of `Q_-j`; `Q_-5` is the one bilinear tail equation and
> is adjoined only on the residual sub-locus where the linear rows leave a kernel.)
> The **generic** case already uses only `[Phi ; M_{-1}]` (`§3`). **Branch A** is
> the same statement per fixed gauge `mu_3` (then `b_-3=mu_3 a_-3` is linear in
> `a_-3`, and the tail is linear in `(a_-2, a_-3)`).

- **Status.** Proved at bounded/generic degree: the explicit datum degree-free in
  the filler (`dV=4..8`, `{Q_0=1,Q_-1}`); 6 distinct slope-`1` data at `d=3`
  (`{Q_0=1,Q_-1}` gap `=1`); branch-A and branch-B **full-tail** unit ideal at
  `d=3` (both branches) and `d=4` (committed verdict). **Open**: arbitrary
  *positive-data* degree; the residual `Q_-1`-surviving sub-locus (needs the deeper
  tail rows / the bilinear `Q_-5`); and a degree-free construction of `lam` (the
  analogue of the `lambda_r` cofactor functional, now on the **joint**
  `Phi ⊕ M_{-1} ⊕ …` cokernel rather than `Phi` alone).

Why the single functional `lambda_{-4}` failed but the *joint* one need not: the
top cofactor functional `lambda_{-4}(f)=f(-1)+f(1)-2f(0)` annihilates `Im Phi` and
pairs `E-R` to `-(R(1)+R(-1))=0` (the cascade relation) — it is blind to the slope
(`w2-theory.md`: the point-annihilator dual `{ev_-1,ev_0,ev_1}` is complete and
`Q_0=1` alone is feasible). Adjoining the tail row `M_{-1}` **enlarges the cokernel**:
the combined operator `[Phi ; M_{-1}]` has a covector that pairs the `rhs`
nontrivially. The Fredholm gap `=1` is exactly the assertion that this enlarged
cokernel is one dimension bigger than the tail row can fill.

## 5. Evidence ledger — proved vs bounded

**Proved (arbitrary coefficient degree, char 0, machine-checked identities):**
- `Q_m=[D,X]_m` for `m in [-6,6]`; `Q_0=(T-1)G` (`§S0`).
- `Im Phi(W2) ⊂ D·F[E]` (every filler block divisible by `D=E(E-1)(E+1)`); slope
  gate `D | (E-R) <=> R(1)=1, R(-1)=-1` (`§S1`).
- **Both-ends Lemma P** and its W2 specialization with the level-3 drop; `R(±1)`
  filler-independent (`§S2`).
- **Branch-B linearity**: `a_-2` is a pure filler (absent from the cascade);
  `Q_0-1, Q_-1..Q_-4` linear and `Q_-5` bilinear in the two fillers (`§S3`).

**Bounded / generic evidence:**
- No constant-cofactor **linear Nullstellensatz certificate** for the (infeasible)
  `d=2` full system — the obstruction is genuinely **nonlinear** (`§S4`).
- **Generic localization** `{Q_0=1, Q_-1}=` unit ideal, both branches, at the
  explicit slope-`1` datum; `{Q_0=1}` alone feasible (`§S5`).
- **Fredholm gap `=1`**, stable across filler degrees `dV=4..8` at the explicit
  datum — that datum has no tail-clearing filler of any degree (`§S6`).
- **Gap `=1` at 6 distinct slope-`1` data**, `L` full column rank — generic in the
  positive data at `d=3` (`§S7`).
- **UNIFORM kill = the full tail.** `{cascade + Q_0=1 + Q_-1}` is **not** the unit
  ideal (either branch, `d=3`) — a slope-`1` sub-locus survives `Q_-1`; the FULL
  system (cascade `∧` `Q_0=1` `∧` `Q_-1..Q_-5`) **is** the unit ideal at `d=3`,
  both branches (`§S8`), reproducing `w2-verdict.md` (`d=3,4`).

**Refuted:** the stated target "`R(1)=0` on cascade+tail" (§2, `§S8` msolve); the
existence of a degree-free **linear** Nullstellensatz certificate (`§S4`); the
over-strong claim that `Q_-1` alone uniformly kills (`§S8`: it kills only the
generic datum).

**Open / not claimed:** the Joint Filler Covector Lemma at arbitrary positive-data
degree; the degree-free covector `lam`; branch-A at arbitrary degree; the
arbitrary-degree W2 joint theorem itself; all of Band 3, DC1, JC2. No Weyl pair and
no counterexample is constructed. The infinite-dimensional
`Im L_K ∩ Im L_H` ([`two-filler-cross-cancellation.md`](two-filler-cross-cancellation.md))
is untouched.

## 6. What it needs to generalize (feeding architecture step 2)

The reduction hands step 2 a **single, uniform, computable object**: the joint
cokernel of `[Phi ; M_{-1}]`. Three concrete openings:

1. **A degree-free `lam`.** The `r != -4` closure used the cofactor functional
   `lambda_r` with `Im Phi ⊂ ker lambda_r` (adjoint / moving-sum criterion,
   `lambda-general-k.md` Thm A'). The joint kill needs the analogous covector for
   the *stacked* operator `Phi ⊕ M_{-1}`. Because `M_{-1}` is itself a two-block
   filler map (top `a_2` and sub `b_1` in place of `a_3, b_2`), the same adjoint
   machinery applies — the target is a moving-sum covector on the **union** of the
   two root necklaces. This is the sharp next computation.
2. **Uniformity across the hatch tower.** `w2-theory.md` §3 makes W2 the band-3
   instance of the step-`(k-1)`/step-`k` hatch, each reducing to the same slope
   gate. The both-ends Lemma P is band-`k` verbatim (it is Lemma P at `E=±1`), and
   the level-`k` drop is the general `a_k(0)=a_k(-2(k-1))=0`, `b_{k-1}(0)=…=0`
   pattern; so the generic `{Q_0=1, Q_{-1}}` kill and the full-tail uniform kill
   are the natural band-`k` conjectures to test next.
3. **The singular-hatch theorem** (architecture step 2) wants a *uniform* reason
   every hatch is silent. The joint covector is that reason if it exists
   degree-free: it is the functional obstructing `Q_0=1 ∧ (Q_-1,…,Q_-5)=0` on the
   filler space (generically already `Q_0=1 ∧ Q_-1=0`).

## 7. Verification

```sh
uv run --with sympy python research/dc1-program/verify_w2_joint.py
```

Runs `§S0` (crossed-product engine, `Q_0=(T-1)G`), `§S1` (`Im Phi ⊂ D·F[E]`, slope
gate), `§S2` (both-ends Lemma P + W2 level-3 drop + filler-independence), `§S3`
(branch-B linearity, `a_-2` a pure filler), `§S4` (no linear Nullstellensatz
certificate), `§S5` (the generic `{Q_0=1, Q_-1}` kill at the explicit datum, both
branches), `§S6` (the branch-B Fredholm gap, stable across filler degree), `§S7`
(6 distinct slope-`1` data, generic gap `=1`), `§S8` (msolve, `d=3`: `Q_-1` alone
is **not** a uniform kill — `cascade+Q_0=1+Q_-1` is not the unit ideal on either
branch; the **full** `cascade+Q_0=1+Q_-1..Q_-5` **is** the unit ideal, both
branches; and the machine-checked refutation `R(1) not in sqrt(cascade+tail)`).
A successful run ends `ALL W2 JOINT CHECKS PASSED`. The exact-`QQ` `d=4`
full-system kill is the committed `w2-verdict.md` msolve certificate.
