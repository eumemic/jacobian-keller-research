# W2 negative tail: the bottom wall, the mirror cube-class gate, and the slope/tail split

**INDEPENDENTLY DERIVED AND MACHINE-VERIFIED — NOT PEER REVIEWED — BAND-SCOPED**

QUANTUM MIRROR assault, band 3: the full `Q_m` system below the fillers for
**W2**, the `r=-4` member of the normalized degree-3 exotic AP family. Every other
member of that AP family dies upstream by the symbolic-`r` `lambda_r` theorem
(`quantum-ap-lambda.md`, commit `d8189fc`); higher-degree and non-AP Band-3 tops
remain open. This memo studies the negative tail for W2.

W2 datum, gauge `b_3 = 0`, quantum band-3 conventions
(`Q_m = sum_(k+l=m)[b_l^[k] a_k - a_k^[l] b_l]`, membership
`(E)_j = E(E-1)...(E-j+1) | a_-j, b_-j`):

```text
a_3 = E(E+2)(E+4)   (roots {0,-2,-4}),     b_2 = E(E+3),     b_3 = 0.
```

Everything below is checked exactly by
[`verify_w2_tail.py`](verify_w2_tail.py) (ends `ALL W2 TAIL CHECKS PASSED`,
~60 checks, ~27 s). All lambda-wave facts this memo builds on are **re-verified
in-file**, not merely cited.

## 0. Headline

> After the positive cascade solves `b_1, b_0, b_-1, b_-2` (from
> `Q_4=Q_3=Q_2=Q_1=0`) and the bottom Wronskian `Q_-6` sets `b_-3 = mu_3 a_-3`,
> **every** coefficient is expressed in the free data. Then `Q_0-1` and the whole
> negative tail `Q_-1,...,Q_-5` become **pure polynomial constraints**. A Gröbner
> decomposition at `d=1` splits the W2 infeasibility cleanly:
>
> ```text
>   (a) positive cascade + Q_0=1  (no tail)   = UNIT IDEAL   -- the slope kills;
>   (b) positive cascade + negative tail       = PROPER IDEAL -- the tail FEASIBLE.
> ```
>
> **At the tested raw caps `d=1,2`, the negative-tail subsystem is feasible without
> the moment equation.** This bounded separate-feasibility statement does not say
> the tail is compatible with slope one. The later bounded verdict
> [`w2-verdict.md`](w2-verdict.md) proves the encoded slope+tail conjunction empty
> at raw cap `d=3` and reports an optional exact check at `d=4`; arbitrary degree
> remains open.
>
> The bottom mirrors the top structurally — `Q_-5` is the reflected wall obeying
> the same `Phi_3(S)=1+S+S^2` necklace criterion — but with one decisive
> asymmetry: the membership-forced roots `{0,1,2}` of `a_-3` are the
> **consecutive (tame) cube**, `= Phi_3` itself, not an exotic step-2 AP. **The
> bottom wall does not independently force W2 exotic.** Separate tail feasibility
> at `d=1,2` does not imply that the tail creates no obstruction to slope-one
> compatibility.

## 1. The complete negative-tail system (arbitrary degree, machine-checked)

`verify` §0 checks `Q_m` equals the direct crossed-product commutator
`[D,X]_m` for **every** `m in [-6,6]` with fully generic degree-2 coefficients
(so the negative tail is on the same exact footing as the positive cascade), and
`Q_0 = (T-1)G`.

### 1.1 `Q_-6` — bottom proportionality (mirror of the top gauge)

`Q_-6 = b_-3^[-3] a_-3 - a_-3^[-3] b_-3`. On the branch `a_-3 != 0` the ratio
`b_-3/a_-3` is `(-3)`-periodic, hence constant:

```text
   b_-3 = mu_3 a_-3,      mu_3 in C.                              (verify §1)
```

`mu_3` carries membership `(E)_3` from `a_-3` automatically. The single available
gauge `D -> D - lambda_3 X` was spent flattening the **top** (`b_3=0`), so `mu_3`
survives **independent of** `lambda_3` — the `lambda_3–mu_3` cross-coupling the
cascade memo flagged (`quantum-band3-cascade.md` §5.2, commit `050a4c0`). If
`a_-3 = 0`, `Q_-6` is vacuous and `mu_3` is unavailable: a separate vanishing
branch.

### 1.2 `Q_-5` — the bottom wall in phi-gauge normal form

With `b_-3 = mu_3 a_-3`, direct substitution gives (verify §2)

```text
   Q_-5 = phi^[-3] a_-3 - a_-3^[-2] phi,      phi := b_-2 - mu_3 a_-2.  (BOTTOM-WALL)
```

The `mu_3 a_-2` correction is exactly what turns the raw-`b_-2` wall plus its
`mu_3`-source into the clean homogeneous wall — the sign in `phi` is fixed by
`verify §2` (`raw_wall + mu_source`). This is the **hardened cascade normal
form** for the bottom wall.

### 1.3 `Q_-4` decomposition and the bottom Lemma-R

```text
   Q_-4 = [ b_-1^[-3] a_-3 - a_-3^[-1] b_-1 ]
        + mu_3 [ a_-3^[-1] a_-1 - a_-1^[-3] a_-3 ]
        + [ b_-2^[-2] a_-2 - a_-2^[-2] b_-2 ]      (verify §3, cascade memo §5.3).
```

The staggered leading-coefficient Lemma-R for the bottom wall
`(a,b)/(a',b') = (0,-3)/(-2,0)` gives `coeff(E^{p+q-1}) = (2p-3q) lc(a_-3)
lc(phi)`, forcing the **mirror degree law `2 deg a_-3 = 3 deg phi`** (verify §3).

## 2. The reflection and the mirror Wall Lemma (arbitrary degree)

The quantum Fourier reflection `E -> -E-1` is a genuine automorphism flipping all
shift signs. `verify §4` checks it **exactly** carries the top wall to the bottom
wall:

```text
   refl( b_2^[3] a_3 - a_3^[2] b_2 )  =  phi^[-3] a_-3 - a_-3^[-2] phi.
```

Hence the entire top-wall Wall Lemma (`quantum-band3-cascade.md` §3.1) transports
to the bottom. In root-multiset data (positions grouped by `mod-Z` coset,
`S` = shift `+1`), the top relation `S(1+S) A = Phi_3(S) B` reflects to

```text
   Phi_3(S) B_phi = (S+1) A_-3,                                   (BOTTOM-M)
```

so a nonzero `phi` solving `(BOTTOM-WALL)` exists **iff**, in every `mod-Z` coset,
(i) the root multiset `A_-3` of `a_-3` is divisible by `Phi_3(S)` and (ii) the
forced quotient `B_phi = (S+1) A_-3 / Phi_3` is a genuine (nonnegative)
multiset. Because `gcd(S+1, Phi_3)=1`, (i) is equivalently `Phi_3 | A_-3`
directly. Summing exponents in `(BOTTOM-M)` reproduces `2 deg a_-3 = 3 deg phi`.
This is the exact reflection of the top-wall story.

## 3. The bottom cube-class gate — and why the base is TAME

The bottom wall is a **genuine gate** on `a_-3` (verify §5):

```text
   a_-3 = (E)_3         roots {0,1,2}=Phi_3:  bottom wall SOLVABLE, phi ~ E(E-1)  (shifted cube);
   a_-3, roots {0,1,2,3,4,5}=Phi_3(1+S^3):    SOLVABLE (two consecutive cubes);
   a_-3, roots {0,1,2,4}  (not Phi_3-divisible): bottom wall GATE FIRES, only phi=0.
```

So not every membership-respecting `a_-3` extends — the roots `{0,1,2,4}` are
killed exactly as the top-wall realizability killed `{1,5,6}`. This is the
reflection of the top-wall refutation.

**The decisive asymmetry.** Membership forces `(E)_3 = E(E-1)(E-2) | a_-3`, i.e.
`{0,1,2} subset roots(a_-3)`. But `{0,1,2}` is the **consecutive** cube
`h h^[1] h^[2]` with `h=E-2` — it **equals** `Phi_3` and is a genuine shifted
cube (the tame necklace), **not** an exotic step-2 AP like the W2 top's
`{0,-2,-4}`. Consequently the membership-minimal `a_-3 = (E)_3` **always** passes
the bottom gate; the genuine positive-control pair
`X=U^3-d/kappa, D=kappa U` realizes exactly this with `a_-3 = c_1^3 (E)_3`
(verify §9, all `Q_m=delta_{m0}` including the whole tail). **The bottom wall
therefore does not independently force W2 to be exotic** — the bottom half can be
as tame as the positive control. All the exotic pressure sits on the top, in the
slope.

## 4. Pivot facts, re-verified (not cited)

`verify §6` re-derives the load-bearing lambda-wave facts from scratch:

- `Im Phi(W2) = D * F[E]`, `D = E(E-1)(E+1)`: forward inclusion (every basis
  filler divisible by `D`) and reverse inclusion (the quotient operator
  `B(V) = (2-E)V - (E+2)V^[1]` is triangular, degree `n -> n+1`, leading coeff
  `-2`, plus the explicit preimage of `1`).
- `E - R in Im Phi <=> D | E-R <=> R(0)=0` (automatic from membership),
  `R(1)=1`, `R(-1)=-1`. (`(E-R)(1)=1-R(1)`, `(E-R)(-1)=-1-R(-1)`.)
- The proved cascade constraint `R(1)+R(-1)=0` is **consistent** with
  `(R(1),R(-1))=(1,-1)`, so W2 hinges on the one scalar `R(1)=1`.

**Correction / refinement on uniqueness.** `L_K: C |-> K_3[(E)_3 C]` and
`L_H: V |-> H_2[(E)_2 V]` are each **individually injective**, but the difference
`Phi = L_K - L_H` has a **nonzero kernel** — the two-filler cross-cancellation
lattice (`../dc1-program/two-filler-cross-cancellation.md`, commit `0889f8a`;
exhibited at low degree in `verify §6`). Hence given `E-R in Im Phi`, the filler
data `(C,V)` — equivalently `(b_-3, a_-2)` — is determined **only modulo
`ker Phi`**, not absolutely. The negative tail is solved over that whole affine
family, which is *extra freedom*, not a determinacy. The phrase "unique given the
divisibility" holds only for each block separately; the difference map does not
inherit it.

## 5. The combined system and the slope/tail split (bounded evidence)

After the positive cascade and `b_-3=mu_3 a_-3`, the full W2 DC1-face system is a
finite set of scalar polynomial equations in the free data
`{a_2,a_1,a_0,a_-1; a_-2, a_-3, mu_3; b_0`-kernel`}`:

```text
   positive-cascade conditions   (constrain a_2,a_1,a_0,a_-1),
   Q_0 - 1 = 0                    (the slope; <=> R(1)=1 and R(-1)=-1),
   Q_-1 = Q_-2 = Q_-3 = Q_-4 = Q_-5 = 0   (the negative tail).
```

Gröbner over `Q` (verify §7–§8):

| `d` | system | ideal |
|---|---|---|
| 1 | positive cascade only | proper (no false kill) |
| 1 | positive cascade + tail (no `Q_0=1`) | **proper — tail feasible** |
| 1 | positive cascade + `Q_0=1` (no tail) | **unit — slope kills** |
| 1 | full (cascade + `Q_0=1` + tail + membership) | unit (infeasible) |
| 2 | positive cascade + tail | **proper — tail feasible** |

At `d=1` the positive cascade forces `R=0` (single branch), so the slope equation
alone already makes the combined system infeasible; after omitting that equation,
the tail is satisfiable with a nondegenerate positive solution (an explicit two-
parameter `a_-2` family survives). At `d=2`, where `codim Im Phi` jumps `2 -> 3`
(`r=-4` is the unique rank-drop locus), the tail subsystem without the moment
equation remains feasible. These are **bounded separate-feasibility computations**
at `d=1,2`; they do not prove that the tail is compatible with slope one at those
caps or at arbitrary degree.

## 6. Escalation readiness

`verify §10` ships `run_full_system(top, d, extra_constraints)`, which builds the
full combined system at cap `d` and reports feasibility. The later decisive memo
supplies an exact slope-`1` point, and [`w2-verdict.md`](w2-verdict.md) subsequently
checks the full fiber: unit at raw cap `d=3` on both branches in committed SymPy,
and reported/optionally reproduced unit at `d=4` with msolve. Thus no candidate is
available at those tested caps. The harness here still self-tests only its stated
low-cap cases; it did not itself prove the later verdict.

## 7. Honest ledger — proved vs bounded

**Proved, arbitrary coefficient degree (char 0), machine-checked identities:**
- `Q_m = [D,X]_m` for all `m in [-6,6]` (§1); `Q_0=(T-1)G`.
- `Q_-6 => b_-3=mu_3 a_-3` (bottom proportionality, §1.1).
- `Q_-5 = phi^[-3] a_-3 - a_-3^[-2] phi`, `phi=b_-2-mu_3 a_-2` (bottom wall
  phi-gauge normal form, §1.2); the `Q_-4` decomposition and bottom Lemma-R
  degree law `2 deg a_-3 = 3 deg phi` (§1.3).
- The reflection `E->-E-1` carrying the top wall to the bottom wall exactly, and
  the mirror necklace criterion `Phi_3 B_phi = (S+1) A_-3` (§2).
- `Im Phi(W2)=D*F[E]` (both inclusions); `E-R in Im Phi <=> R(1)=1, R(-1)=-1`;
  the cross-cancellation kernel of `Phi` (§4).

**Bounded evidence:**
- The bottom cube-class gate by concrete necklace examples: `(E)_3` and
  `{0,1,2,3,4,5}` solvable, `{0,1,2,4}` gated (§3, verify §5).
- `R=0` and `R(1)+R(-1)=0` at `d=1` from the solved cascade (verify §7).
- The full-system unit-ideal + slope/tail Gröbner decomposition at `d=1`, and
  separate tail-subsystem feasibility without `Q_0=1` at `d=2` (verify §8).

**Explicit exceptional loci.** The `mu_3` gauge is unavailable on `a_-3=0`
(separate vanishing branch). If the filler forces `C=0` (`b_-3=0`), `E-R` must lie
in `Im L_H` alone (codimension 4), a strictly stronger condition; likewise
`mu_3=0` collapses `b_-3=0` with `a_-3` possibly nonzero — each a separate branch
of the tail, not covered by the `mu_3 != 0`, `a_-3 != 0` main line above.

**Open / not claimed.** Arbitrary-degree tail feasibility (only `d=1,2` here);
combined feasibility beyond the self-contained exact raw cap `d=3` in the bounded
verdict; exact `QQ` raw cap `d=4` unless the optional run is completed;
arbitrary-degree W2; higher-degree/non-AP Band-3 tops, DC1, and JC2. The moment
slope itself is achievable by the later `w2-decisive.md`, but is incompatible with
the encoded tail at raw cap `d=3`. No Weyl pair or counterexample is constructed. In unrestricted degree
`Im L_K intersect Im L_H` remains infinite-dimensional
(`../dc1-program/two-filler-cross-cancellation.md`); nothing here weakens that.

## 8. Verification

```sh
uv run --with sympy python research/band3/verify_w2_tail.py
```

runs §0 (crossed-product engine, all `m in [-6,6]`), §1 (bottom proportionality,
membership), §2 (bottom wall phi-gauge normal form), §3 (`Q_-4` decomposition,
Lemma-R), §4 (reflection, mirror necklace criterion), §5 (the bottom cube-class
gate by example), §6 (pivot facts re-derived, `Phi`-kernel note), §7 (`R=0`,
`R(1)+R(-1)=0` at `d=1`), §8 (combined-system feasibility + slope/tail
decomposition, `d=1,2`), §9 (positive-control genuine pair), §10 (escalation
harness). A successful run ends `ALL W2 TAIL CHECKS PASSED`.
