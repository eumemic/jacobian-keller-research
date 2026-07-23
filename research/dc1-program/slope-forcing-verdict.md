# Slope-forcing verdict: the W2 negative tail forces R(1)=0 (d=3, both branches)

**INDEPENDENTLY DERIVED — EXACT ALGEBRA / MACHINE-CHECKED IDENTITIES — NOT PEER
REVIEWED**

This memo adjudicates the **FLAGGED single-engine `msolve` verdict** of
[`joint-covector.md`](joint-covector.md) §5 (the second `**`→`^` flip) and its
sibling [`w2-joint-theorem.md`](w2-joint-theorem.md) §2. The question, verbatim:

> On the variety `{`positive cascade `Q_4..Q_1=0`, negative tail `Q_-1..Q_-5=0`
> (**no** `Q_0` condition), membership, both exhaustive `Q_-6` branches — A
> (`a_-3=(E)_3 am3`, `b_-3=mu_3 a_-3`) and B (`a_-3=0`, `b_-3=(E)_3 C`)`}` at raw
> coefficient cap `d=3`, is `R(1)=G(1)` **identically zero**? Equivalently: is
> `R(1) ∈ sqrt(cascade+tail+membership)`?

W2 datum, gauge `b_3=0`, quantum band-3 conventions
(`Q_m = sum_(k+l=m)[b_l^[k] a_k - a_k^[l] b_l]`, `f^[n](E)=f(E+n)`,
membership `(E)_j | a_-j, b_-j`): `a_3=E(E+2)(E+4)`, `b_2=E(E+3)`.

Exact certificate: [`verify_slope_forcing.py`](verify_slope_forcing.py) (ends
`ALL SLOPE FORCING CHECKS PASSED`; default sympy legs, `HEAVY=1` adds the msolve
full-system cross-check). Every load-bearing upstream fact — the crossed-product
engine, `Q_0=(T-1)G`, the slope gate, `R(1)=G(1)` filler-independence, the positive
cascade — is **re-derived in file**.

## Verdict

> **YES — `R(1)=0` identically on `cascade+tail` at `d=3`, on BOTH branches.** That is,
> `R(1) ∈ sqrt(cascade+tail)`: the negative tail is **SLOPE-FORCING**. The flagged
> `^` verdict of `joint-covector.md` §5 is **CONFIRMED**; the original `**`-parser
> claim of `w2-joint-theorem.md` §2 — *"`R(1)` is a nonvanishing free modulus on
> `cascade+tail`; the joint kill is **not** slope-forcing"* — is a `msolve` parser
> artifact and is **REFUTED**.
>
> The membership is **NON-VACUOUS**: `cascade+tail` is non-empty (the origin lies on
> it, with `R(1)=0`) and positive-dimensional. `R(1)` is a genuine *non-constant free
> modulus on the cascade alone* (it attains `1`), so `R(1)≡0` on `cascade+tail` is a
> real forcing, not a triviality.

Two **independent Gröbner engines** agree, both branches, over two large primes and
over `QQ`:

| engine | branch | mod 65003 | mod 32003 | over QQ |
|---|---|---|---|---|
| **sympy** (parametrized reduction, this file) | B | UNIT | UNIT | **UNIT** |
| **sympy** (parametrized reduction, this file) | A | UNIT | UNIT | **UNIT** |
| **msolve 0.10.1 `^`** (full un-reduced d=3 system) | B | UNIT (25s) | UNIT (27s) | **UNIT (31s)** |
| **msolve 0.10.1 `^`** (full un-reduced d=3 system) | A | UNIT (134s) | UNIT (128s) | **UNIT (144s)** |

("UNIT" = the Rabinowitsch ideal `cascade + tail + {1 - t·R(1)}` is the unit ideal,
i.e. `R(1)` vanishes on the whole variety — a machine-checked **QQ Nullstellensatz
radical certificate**, §4.)

## 1. The reduction, and why the slope collapses to one monomial

The positive cascade `Q_4..Q_1=0` is triangular and forward-solves `b_1,b_0,b_-1,b_-2`
(the single `b_0` constant kernel retained). After substitution the free data is
`{a_2,a_1,a_0,a_-1}` (raw cap `d`) plus the fillers `a_-2=(E)_2 V`, and `b_-3=(E)_3 C`
(branch B) or `a_-3=(E)_3 am3, mu_3` (branch A). The cascade conditions cut a
**9-dimensional** variety (10 conditions, 8 independent), and admit an explicit
rational **parametrization** verified in file to satisfy every cascade condition
identically. On it the moment slope **collapses to a single monomial**:

```text
   R(1) = G(1) = Q_0(0)  ≡  -108 · a2_0 · am1_3      (EXACT identity on the cascade),
```

where `a2_0 = a_2(0)` and `am1_3` is the top coefficient of the membership quotient of
`a_-1`. (`R(1)` is *filler-independent* — the both-ends Lemma P, re-derived in file —
so it is a function of the positive data alone.) So the entire question becomes:

> **does the tail force `a2_0 · am1_3 = 0`?**  — and the answer is **yes**.

Because `a2_0` and `am1_3` are two of the nine **free** cascade coordinates, `R(1)` is
manifestly a non-constant free modulus on the cascade (it attains `1` at
`a2_0=1, am1_3=-1/108`). The tail — a *filler*-coupled system — nevertheless forces
their product to zero.

## 2. Feasibility (the side question): `cascade+tail` is NON-EMPTY

If `cascade+tail` were empty the radical membership would be vacuous. It is not:

- **The origin** (all free data and fillers `=0`) lies on `cascade+tail`, with
  `R(1)=0`, on **both** branches (`verify §S2`). So `cascade+tail` is non-empty and
  the forcing is non-vacuous.
- It is **positive-dimensional**: `cascade` alone is a proper (feasible) ideal
  (control, `msolve`/sympy), and the `d=1` two-parameter `a_-2` tail family of
  [`../band3/w2-negative-tail.md`](../band3/w2-negative-tail.md) §5 embeds. On all of
  it `R(1)=0`.

The flagged memo asserted feasibility only from *`msolve` non-termination* on
`cascade+tail` — which is **not** evidence of feasibility (it is a slow
positive-dimensional parametrization). The origin settles it cleanly.

## 3. The witness hunt (refutation direction) fails

A single `cascade+tail` point with `R(1)≠0` would refute membership and restore the
"free modulus" reading. There is none:

- `cascade + tail + {R(1)=1}` is the **unit ideal** (empty) — no slope-`1` tail point
  (`msolve` 1.1s; sympy). This also explains, cleanly, why the explicit slope-`1`
  datum of [`w2-decisive.md`](../band3/w2-decisive.md) **fails the tail on both
  branches** (recorded there): a slope-`1` point has `R(1)=1≠0`, so it *cannot* be on
  `cascade+tail`.
- The Rabinowitsch test (§Verdict) is exactly "no `cascade+tail` point has `R(1)≠0`",
  unit on both engines.

## 4. The exact certificate

Because `R(1)=-108·a2_0·am1_3` on the cascade, radical membership is the machine
statement **`(a2_0·am1_3) ∈ sqrt(tail-ideal)`**. The **positive certificate** produced:

- **QQ radical certificate (exact, char 0).** The reduced Gröbner basis of
  `parametrized-tail + {1 - t·a2_0·am1_3}` over `QQ` is `[1]` (sympy), both branches —
  a machine-checked Hilbert-Nullstellensatz certificate that `(a2_0·am1_3)^m` (hence
  `R(1)^m`) lies in the cascade+tail ideal for some finite `m` (`verify §S4/§S5`).
  The same is `[1]` over `GF(65003)`, `GF(32003)`, and — on the full un-reduced
  `d=3` system — via `msolve` over `QQ` and both primes (`§S7`, `HEAVY`).

*(Honest scope: the **minimal** power `m` is not separately extracted. The only
tractable route is the saturated elimination Gröbner basis
`(elim⟨Q_-1..Q_-4⟩ : det^∞)`, whose generators are high-degree after the cascade
parametrization and whose GB does not terminate in a sane budget; `m` is left bounded,
finite and `≥1`. The radical membership itself is fully certified above.)*

## 5. Graded depth: `kmin = 3` (and the reconciliation with the Fredholm mechanism)

The minimal tail depth that forces `R(1)=0` is **`k=3`**:

- `R(1) ∈ sqrt(cascade + Q_-1,Q_-2,Q_-3)` — **UNIT** (msolve, 100s). So `Q_-4,Q_-5`
  are **not needed** for slope-forcing.
- `R(1) ∉ sqrt(cascade + Q_-1)`: an **explicit witness** on `cascade+Q_-1` with
  `R(1)=-108≠0` (also `-432`, `verify §S6`) — so `Q_-1` alone does **not** force the
  slope (rigorous). `cascade+Q_-1,Q_-2` is likewise not unit within the `msolve`
  budget (evidence, not proof of `kmin>2`). So rigorously **`1 < kmin ≤ 3`**, and the
  msolve stratification pins it at `3`.

**Reconciliation.** This does **not** contradict the Fredholm-gap picture of
[`joint-covector.md`](joint-covector.md) §4 / `w2-joint-theorem.md` §3, which is about
`cascade + Q_0=1 + Q_-1` (the *moment-unit* system). Those are two genuinely different
mechanisms, and both are real:

- **Filler obstruction (depth 1, with `Q_0=1`).** On the slope-`1` locus (`R(1)=1`
  forced by `Q_0=1` via the slope gate), the *linear filler* operator `[Phi ; M_-1]`
  has its `rhs` outside its column span — no filler realizes `Q_0=1` and `Q_-1=0`
  simultaneously. This kills the **generic** slope-`1` datum and needs only `Q_-1`.
  It is **not** slope-forcing (`Q_-1` alone does not force `R(1)=0`, per the witness
  above).
- **Slope-forcing (depth 3, no `Q_0`).** The full negative tail through `Q_-3` forces
  `R(1)=0` outright, independently of any `Q_0` condition. This kills **every** datum
  of slope `≠0`, uniformly.

So `w2-joint-theorem.md` §2's *headline* ("not slope-forcing; `R(1)` a free modulus")
is wrong, but its §3–§4 *Fredholm localization at `Q_0=1+Q_-1`* stands as a separate,
compatible fact. The W2 full-system kill (`cascade ∧ Q_0=1 ∧ tail`, unit at `d≤4`,
[`../band3/w2-verdict.md`](../band3/w2-verdict.md)) is now **explained the simple way**:
the depth-3 tail forces `R(1)=0` while `Q_0=1` demands `R(1)=1` — an immediate
contradiction.

## 6. Consequence for the degree-free W2 target

The degree-free W2 target **simplifies**. The object to prove at arbitrary
positive-data degree is now a **slope statement**:

```text
   cascade ∧ (Q_-1=Q_-2=Q_-3=0)   ⇒   R(1)=0        (W2 slope-forcing lemma),
```

equivalently `a_2(0)·am1_3 ∈ sqrt(cascade+Q_-1..Q_-3)` — rather than the joint-filler
covector / Fredholm cokernel apparatus for the full tail. Because `R(1)=G(1)` is the
3-term boundary pairing of the both-ends Lemma P (degree-free) and `Q_-1,Q_-2,Q_-3`
are the first three reflected-tail walls, this is a sharply posed, low-height target.

## 7. Evidence ledger — proved / exact / bounded / refuted

**Proved (arbitrary coefficient degree, char 0, machine identities):**
- `Q_m=[D,X]_m` for `m∈[-6,6]`; `Q_0=(T-1)G` (`§S0`).
- slope gate `D|(E-R) ⇔ R(1)=1,R(-1)=-1`; `R(1)=Q_0(0)=G(1)`
  **filler-independent** (both-ends Lemma P) (`§S1`).

**Exact (d=3, char 0 + machine-checked):**
- Cascade parametrization valid (all conditions vanish identically), `dim=9`, and
  the identity `R(1) ≡ -108·a2_0·am1_3` on the cascade (`§S3`).
- `R(1) ∈ sqrt(cascade+tail)` **over `QQ`**, both branches, by the parametrized
  Rabinowitsch unit ideal (sympy) (`§S4`), with a two-prime cross-check — a
  machine-checked Nullstellensatz radical certificate (`§S5`).
- `cascade+tail` non-empty (origin), witness hunt for `R(1)≠0` fails (`§S2,§S6`).

**Bounded (d=3, `msolve` `^`, `HEAVY=1`):**
- The **full un-reduced** `d=3` system `cascade+tail+{1-t·R(1)}` is the unit ideal over
  `GF(65003)`, `GF(32003)` **and `QQ`**, both branches (`§S7`) — the independent
  second engine on the un-parametrized variety.

**Refuted / corrected:**
- `w2-joint-theorem.md` §2: *"`R(1)` not in `sqrt(cascade+tail)`; free modulus; not
  slope-forcing"* — a `**`-parser artifact. Corrected: `R(1) ∈ sqrt(cascade+tail)`,
  the tail is slope-forcing (depth 3).

**Open / not claimed:**
- The slope-forcing lemma at **arbitrary positive-data degree** (`d>4`); this file is a
  bounded `d=3` certificate (exact `QQ` on the parametrized reduction; `msolve` on the
  full system). No Weyl pair, no DC1 counterexample; all of Band 3, DC1, JC2 remain
  open. The infinite-dimensional `Im L_K ∩ Im L_H` is untouched.

## 8. Verification

```sh
uv run --with sympy python research/dc1-program/verify_slope_forcing.py
# full-system msolve cross-check (needs msolve 0.10.1 on PATH):
HEAVY=1 uv run --with sympy python research/dc1-program/verify_slope_forcing.py
```

`§S0` engine, `§S1` slope gate + `R(1)` filler-independence, `§S2` feasibility
(origin), `§S3` parametrization + `R(1)=-108 a2_0 am1_3`, `§S4` the independent sympy
Rabinowitsch (`GF(p1),GF(p2),QQ`, both branches), `§S5` the QQ Nullstellensatz
radical certificate, `§S6` witness hunt + graded `kmin`, `§S7` (`HEAVY`) the `msolve`
full-system cross-check.
Ends `ALL SLOPE FORCING CHECKS PASSED`.
