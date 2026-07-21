# Quantum band-2, the $a_2=0$ sector: reduction to the square sector and band 1

**INDEPENDENTLY DERIVED AND MACHINE-VERIFIED — NOT PEER REVIEWED — BAND-SCOPED**

QUANTUM MIRROR assault, $a_2=0$ edge. This memo mirrors, in $A_1$, the classical
$a_2=0$ analysis assembled upstream in
[`../band2-classical-full/full_classical_band2_theorem.md`](../band2-classical-full/full_classical_band2_theorem.md)
(commit `f637b1a`, Sections 1–3: the orientation operations and the independent
band-1 lemma). It supplies the **$a_2=0$ reduction infrastructure** of the
quantum band-2 program: a bracket- and membership-preserving problem-reduction
carrying **every** genuine $A_1$ pair $[D,X]=1$ with both ladder supports in
$[-2,2]$ and $a_2=0$ to one of two smaller problems —

- **(i)** the quantum $a_2\ne0$ **square sector** (worked elsewhere:
  [`quantum-mirror.md`](quantum-mirror.md) at commit `ad43ab5` closes the
  resistant constant-$h$ branch; `quantum-completion.md` supplies the repaired
  complete shifted-square classification), or
- **(ii)** the quantum **band-1** classification (now AUDITED:
  [`quantum-band1.md`](quantum-band1.md) Theorem 1 + Corollary 2, verdict
  CONFIRMED WITH REPAIRS — see §3 and §6).

### What this memo closes, and what it does not

- **Closed here, unconditionally (arbitrary degree):** the reduction
  infrastructure. The verifier checks exactly that the two routing operations
  preserve $[D,X]=1$ and $A_1$-membership (§2); the written Boolean case split
  (§4) gives an *exhaustive explicit partition* of the $a_2=0$ locus with no
  configuration left unrouted.
- **Closed here via the audited band-1 theorem:** **Sector Z** (§3), the
  band-1 collapse (input: `quantum-band1.md`, audited CONFIRMED WITH REPAIRS).
- **Closed in the full assembly, not by this memo alone:** Sectors O1–O3 are
  routed to $a_2\ne0$; `quantum-M4.md` and `quantum-completion.md` now close that
  target. This memo itself proves only the routing equivalence.

Everything asserted here is proved at arbitrary degree. The exact algebraic
identities, membership bookkeeping, and routing witnesses are checked by
[`verify_quantum_a2zero.py`](verify_quantum_a2zero.py) (ends
`ALL QUANTUM A2-ZERO CHECKS PASSED`); Boolean exhaustiveness is the direct
written case split in §4.

## Verdict

> **Reduction Theorem ($a_2=0$ edge).** Let $X=\sum_{k=-2}^{2}x^ka_k(E)$,
> $D=\sum_{l=-2}^{2}x^lb_l(E)$ be genuine elements of $A_1$ (with the
> falling-factorial memberships $E\mid a_{-1},b_{-1}$ and $E(E-1)\mid
> a_{-2},b_{-2}$) satisfying $[D,X]=1$, and suppose $a_2=0$. Then the four
> Boolean conditions of §5 **partition** the locus (mutually exclusive,
> exhaustive), and:
>
> - **Sector Z** ($a_{-2}=b_2=b_{-2}=0$): both $X,D$ are band-1; the pair is
>   an affine symplectic pair by the audited band-1 theorem
>   (`quantum-band1.md` Thm 1). *(Closed.)*
> - **Sectors O1–O3** (respectively $b_2\ne0$; $\,b_2=0,a_{-2}\ne0$;
>   $\,b_2=a_{-2}=0,b_{-2}\ne0$): quantum pair-exchange and/or the quantum Fourier
>   reflection carry the pair to an equivalent genuine $A_1$ pair
>   $[\,\widetilde D,\widetilde X\,]=1$ with the same support, same generated
>   subalgebra, and $\widetilde a_2\ne0$. *(Routed here; closed only after citing M4q and the repaired shifted-square
>   theorem in the full assembly.)*

This is the quantum analogue of the classical **orientation lemma** (Section 3
of the full classical memo, `f637b1a`) plus the classical **band-1 lemma**
(Section 2). As in the quantum-mirror closure (`ad43ab5`), the deformed argument
is in one respect *cleaner* than the classical one: the degree-reversing symmetry
is realised by a genuine algebra **automorphism** (the Fourier transform), so the
bracket transfers with **no sign flip**, whereas the classical swap
$R(x,\xi)=(\xi,x)$ has determinant $-1$ and must be compensated (classical (1.3):
$\{RG,RF\}=-R\{G,F\}$).

## 1. Setup and conventions

The crossed-product calculus and the $Q_m$ system below follow the audited
partial cascade
[`../band2-m5-partial/quantum-shifted-square-sector-partial.md`](../band2-m5-partial/quantum-shifted-square-sector-partial.md)
(commit `91a053a`) and [`quantum-mirror.md`](quantum-mirror.md) (`ad43ab5`); they
are re-derived and machine-checked independently here (`verify` §0).

Work in $A_1[x^{-1}]=\bigoplus_{k\in\mathbb Z}x^k\mathbb C[E]$, $E=x\partial$, with
$$(x^af(E))(x^bg(E))=x^{a+b}f(E+b)g(E),\qquad
\partial=x^{-1}E,\qquad f^{[r]}(E)=f(E+r).$$
For band-2 $X=\sum_{-2}^{2}x^ka_k$, $D=\sum_{-2}^{2}x^lb_l$ the ladder-$m$
coefficient of $[D,X]$ is
$$Q_m=\sum_{k+l=m}\bigl(b_l^{[k]}a_k-a_k^{[l]}b_l\bigr),\qquad Q_m=\delta_{m0}.$$

**Membership (quantum polynomiality).** $x^{-r}c(E)\in A_1$ iff the falling
factorial
$$E^{\underline r}:=E(E-1)\cdots(E-r+1)\ \text{ divides }\ c(E).$$
Thus $E\mid a_{-1},b_{-1}$ and $E(E-1)\mid a_{-2},b_{-2}$. This is the quantum
deformation of the classical $\tau^r\mid a_{-r}$ (`verify` §0, §2). Note
$\partial^{\,j}=x^{-j}E^{\underline j}$, so a genuine negative-level piece
$x^{-j}f(E)=\partial^{\,j}g(E)$ with $f=E^{\underline j}g$.

**Scope.** This is a theorem about the reduction of the fixed support-contained
$a_2=0$ edge to the square and band-1 sectors. It does **not** claim DC1: a
general Weyl pair need not admit the band-2 bound.

## 2. Two exact operations (`verify` §1, §2, §3)

### 2.1 Pair-exchange

If $[D,X]=1$, put $(X^\ast,D^\ast)=(D,-X)$. Then
$$[D^\ast,X^\ast]=[-X,D]=-[X,D]=[D,X]=1,$$
an **operator identity** checked for fully generic band-2 $X,D$ (`verify` §1).
The new top coefficient is $a_2^\ast=b_2$. Pair-exchange preserves $A_1$
membership (both entries are unchanged as operators), band-2 support, and the
generated subalgebra $\langle X,D\rangle=\langle X^\ast,D^\ast\rangle$ — so it
preserves the $A_1$-generation property that DC1 asks about.

### 2.2 Quantum Fourier reflection

Let $\varphi$ be the Fourier automorphism of $A_1$,
$$\boxed{\varphi:\quad x\mapsto-\partial,\quad \partial\mapsto x,\quad
E=x\partial\mapsto -\partial x=-(E+1).}$$
$\varphi$ is a genuine **algebra automorphism**: $[\varphi\partial,\varphi
x]=[x,-\partial]=1$, and $\varphi(AB)=\varphi(A)\varphi(B)$ is verified on all
ladder-generator pairs (`verify` §2). It satisfies $\varphi^4=\mathrm{id}$ and
$\varphi^2=$ parity ($x^kf\mapsto(-1)^kx^kf$), both machine-checked, so it is
invertible. Because $\varphi(x)=-\partial=-x^{-1}E$ has ladder degree $-1$,
$\varphi$ **reverses ladder degree**, $A_1^{(k)}\to A_1^{(-k)}$.

**Closed-form reflection dictionary** (`verify` §2). For $X=\sum x^ka_k$,
$$\begin{aligned}
(\varphi X)_{-2}&=E(E-1)\,a_2(-E-1), &
(\varphi X)_{-1}&=-E\,a_1(-E-1), &
(\varphi X)_{0}&=a_0(-E-1),\\[2pt]
(\varphi X)_{1}&=\frac{a_{-1}(-E-1)}{-E-1}, &
(\varphi X)_{2}&=\frac{a_{-2}(-E-1)}{(E+1)(E+2)}. & &
\end{aligned}$$
For $k\ge0$ this is $\varphi(x^ka_k)=(-1)^k\partial^{\,k}a_k(-E-1)
=(-1)^kx^{-k}E^{\underline k}a_k(-E-1)$; for $k=-j<0$ it is
$\varphi(\partial^{\,j}g)=x^jg(-E-1)$ with $a_{-j}=E^{\underline j}g$. The
divisions in $(\varphi X)_1,(\varphi X)_2$ are exact **precisely by
membership** $E\mid a_{-1}$, $E(E-1)\mid a_{-2}$ — the falling-factorial
bookkeeping the classical reflection did not need. $\varphi$ maps genuine $A_1$
elements to genuine $A_1$ elements (positive images carry no constraint; the
negative images carry the required falling factors automatically — `verify` §2).

**Bracket transfer (no sign).** The commutator of two genuine band-2 elements is
again a genuine $A_1$ element (membership at every negative level, `verify` §3),
so $\varphi([D,X])$ is defined and, $\varphi$ being an automorphism,
$$\boxed{[\varphi D,\varphi X]=\varphi([D,X])}$$
**exactly** (checked for generic band-2 $X,D$, `verify` §3). Since $\varphi$
fixes scalars ($\varphi(1)=1$), a Keller pair maps to a Keller pair with **no
compensating sign** — the point where the deformation improves on the classical
determinant-$(-1)$ swap.

**Orientation fact (nonvanishing, no cancellation).** By membership,
$a_{-2}=E(E-1)\,g$ with $g\in\mathbb C[E]$, and $g\ne0\iff a_{-2}\ne0$ (division
by a fixed nonzero polynomial). Then $(\varphi X)_2=g(-E-1)$. The substitution
$g\mapsto g(-E-1)$ is composition with the *invertible* affine map $E\mapsto
-E-1$, hence a ring automorphism of $\mathbb C[E]$; it sends nonzero to nonzero
with $\deg g(-E-1)=\deg g$ and no leading-term cancellation. Therefore
$$(\varphi X)_2\ne0\quad\Longleftrightarrow\quad a_{-2}\ne0\qquad(\text{`verify` §2}).$$
Thus $\varphi$ moves a nonzero level-$(-2)$ coefficient into the top slot.

**Shifted-square transfer.** If $a_{-2}$ is in the shifted-square class *modulo
the membership factor*, $a_{-2}=E(E-1)\,c\,h(E)h(E-1)$, then
$$(\varphi X)_2=c\,H(E)H(E-1),\qquad H(E)=h(-E-2)$$
(an exact identity, `verify` §2), again shifted-square — matching the J3q
reflection lemma
([`../band2-j3-provisional/theorem_J3_band2_nonsquare.md`](../band2-j3-provisional/theorem_J3_band2_nonsquare.md),
commit `5a76673`: "$c\,h(E)h(E-1)$"). So when the source is shifted-square the
reflection lands *inside* that class, not merely near it. (For a general
nonzero $a_{-2}$ the image top is a general nonzero polynomial — see §7.)

## 3. Quantum band-1 rigidity (Sector Z) — repaired input

If $a_2=a_{-2}=b_2=b_{-2}=0$ then both supports lie in $\{-1,0,1\}$: the pair is
band-1. Two facts (`verify` §0):

1. **Collapse.** With all four extremes zero, the band-2 system $Q_m$ is
   *identically equal* to the band-1 system for $-2\le m\le2$, and
   $Q_{\pm3}=Q_{\pm4}=0$ identically. So Sector Z is exactly the band-1 problem,
   with no residual band-2 coupling.
2. **Classification.** By the **audited quantum band-1 theorem**
   ([`quantum-band1.md`](quantum-band1.md), Theorem 1 and Corollary 2 — the
   repaired, citable form of the provisional Theorem P3; audit at
   [`../../archive-import/provisional/dixmier-band-program/AUDIT-band1.md`](../../archive-import/provisional/dixmier-band-program/AUDIT-band1.md),
   verdict CONFIRMED WITH REPAIRS): every band-1 pair $[D,X]=1$ **in $A_1$**
   is an affine symplectic pair, and such a pair generates $A_1$. That
   $A_1$-level statement is exactly what Sector Z consumes; nothing here uses
   P3's retracted $B\setminus A_1$ "affine + polar" classification (the audit
   exhibits a band-1, non-affine, non-polar localized pair
   $X=x^2\partial+x+x^{-1},\ D=x^{-1}$, which membership excludes from $A_1$),
   nor P3's mis-identified use of the Fourier map inside the localization
   (the audit's Lemma D closes the degenerate vanishing patterns instead).

**Sector Z is therefore closed by the audited band-1 theorem.** (Novelty
framing per the audit: the band-1 $A_1$ statement is a mild instance of
Bavula–Levandovskyy 2020 / Han–Tan 2024; it is consumed here as
infrastructure, not claimed as new.)

## 4. Quantum orientation lemma (Sector O)

### Lemma

Let $[D,X]=1$ be a support-contained band-2 pair with $a_2=0$ that is **not** in
band 1 (i.e. the extreme triple $(b_2,a_{-2},b_{-2})\ne(0,0,0)$). Then §2's
operations produce an equivalent genuine $A_1$ pair, still
$[\,\cdot\,,\cdot\,]=1$ with band-2 support and the same generated subalgebra,
whose top coefficient is nonzero.

### Proof (explicit exhaustive partition, `verify` §4)

The three Boolean conditions below **partition** $\{a_2=0,\ \text{not band-1}\}$:
they are mutually exclusive by construction (each assumes the negation of the
previous predicate) and exhaustive because $(b_2,a_{-2},b_{-2})\ne(0,0,0)$. No
degree comparison is used; only the zero/nonzero status of each extreme
polynomial, treated as an explicit branch.

- **O1 — $b_2\ne0$** (no condition on $a_{-2},b_{-2}$). Pair-exchange
  $(X,D)\mapsto(D,-X)$ gives new top $a_2^\ast=b_2\ne0$ (§2.1, exact relabelling).
  *(Witness: $X=x,\ D=\partial+x^2$, $[D,X]=1$, $a_2=0$, $b_2=1$; exchange gives
  $(\partial+x^2,-x)$ with $a_2=1$, `verify` §4.)*
- **O2 — $b_2=0,\ a_{-2}\ne0$.** Apply $\varphi$: $(\varphi X)_2\ne0$ by the
  orientation fact (nonvanishing via invertible substitution, §2.2). *(Witness:
  $X=x+\partial^2,\ D=\partial$, $[D,X]=1$, $a_2=0$, $a_{-2}=E(E-1)$; $\varphi$
  gives $(x^2-\partial,\ x)$ with $a_2=1$, `verify` §4.)*
- **O3 — $b_2=0,\ a_{-2}=0,\ b_{-2}\ne0$.** Pair-exchange first (new
  $a_{-2}^\ast=b_{-2}\ne0$, new $a_2^\ast=b_2=0$), then $\varphi$ as in O2; the
  composite fills the top slot (`verify` §4).

Each step preserves $[D,X]=1$ (§2.1, §2.2), band-2 support, membership, and the
generated subalgebra. $\qquad\blacksquare$

The oriented pair has $a_2\ne0$. The later theorem `quantum-M4.md` proves at
arbitrary degree that its top coefficient is a shifted square, and the repaired
`quantum-completion.md` plus `quantum-mirror.md` classify that sector. Thus the
routing proved here is now a completed input to `quantum-band2-theorem.md`; it
was open when this memo was first drafted.

## 5. The complete $a_2=0$ sector tree (current ledger)

| Sector | Defining condition | Routing | Current status |
|---|---|---|---|
| Z | $a_{-2}=b_2=b_{-2}=0$ | band-1 collapse | Closed by repaired `quantum-band1.md` |
| O1 | $b_2\ne0$ | pair exchange | Routed here; target closed by M4q + completion |
| O2 | $b_2=0,a_{-2}\ne0$ | Fourier $\varphi$ | Routed here; target closed by M4q + completion |
| O3 | $b_2=a_{-2}=0,b_{-2}\ne0$ | exchange then Fourier | Routed here; target closed by M4q + completion |

The partition and routing are the arbitrary-degree contribution of this memo.
The classification of the targets is cited, not re-proved or machine-certified
here. No localized involution $x\mapsto x^{-1}$ is used: orientation is by the
genuine Fourier automorphism of $A_1$.

## 6. Band-1 input consumed

Sector Z uses exactly Theorem 1 and Corollary 2 of `quantum-band1.md`: genuine
$A_1$ band-1 pairs are affine symplectic and generate $A_1$. The former
localized “affine + polar” classification remains retracted and is irrelevant.

## 7. Retrospective status

The open items listed in the original draft have since been supplied: the
band-1 repair, M4q shifted-square theorem, and complete shifted-square sector.
Their assembly is `quantum-band2-theorem.md`. This update is retrospective; the
present memo's verifier checks exact orientation identities and witnesses, while
the written Boolean split supplies exhaustiveness.

## 8. Contrast with the classical template

| step | classical (full memo §1–3) | quantum (here) |
|---|---|---|
| degree-reversal | swap $R(x,\xi)=(\xi,x)$, $\det=-1$ | Fourier $\varphi$, $x\mapsto-\partial$, a genuine automorphism |
| bracket under it | $\{RG,RF\}=-R\{G,F\}$ (sign) | $[\varphi D,\varphi X]=\varphi([D,X])$ (**no sign**) |
| reflected coeff | $(Rf)_k=t^{-k}a_{-k}$ | $(\varphi X)_2=a_{-2}(-E-1)/((E+1)(E+2))$, etc. |
| membership divisor | $\tau^r\mid a_{-r}$ | $E^{\underline r}\mid a_{-r}$ (falling factorial) |
| band-1 input | independent lemma (memo §2) | Theorem **P3** (conditional) |
| orientation target | M4 ($a_2$ square) + M5 | shifted-square sector + quantum-$M4$ |

As in [`quantum-mirror.md`](quantum-mirror.md) §5 (Lemma R), the difference
setting supplies rigidity the differential setting lacks; here it takes the mild
but real form that the degree-reversing symmetry is a *genuine automorphism*, so
the orientation half of band 2 transfers to $A_1$ without the classical sign
compensation, at the cost only of the falling-factorial membership bookkeeping.

## 9. Machine verification (`verify_quantum_a2zero.py`)

- **§0** — crossed-product engine: $Q_m=$ direct commutator; band-1 collapse
  ($a_2=a_{-2}=b_2=b_{-2}=0\Rightarrow$ band-1 system; $Q_{\pm3}=Q_{\pm4}=0$).
- **§1** — pair-exchange operator identity $[-X,D]=[D,X]$; top relabelling.
- **§2** — $\varphi$ on generators; $\varphi(AB)=\varphi(A)\varphi(B)$;
  $\varphi^4=\mathrm{id}$, $\varphi^2=$ parity; the closed-form reflection
  dictionary; orientation fact $a_{-2}\ne0\Rightarrow(\varphi X)_2\ne0$;
  shifted-square transfer; membership preservation.
- **§3** — genuineness of $[D,X]$; the exact automorphism identity
  $[\varphi D,\varphi X]=\varphi([D,X])$; $\varphi(1)=1$.
- **§4** — the three orientation leaves O1–O3 on concrete genuine $[D,X]=1$
  witnesses, each oriented to $a_2\ne0$ with the bracket preserved.

Run from the repository root:
```sh
uv run --with sympy python research/band2-square-sector/verify_quantum_a2zero.py
```
Ends `ALL QUANTUM A2-ZERO CHECKS PASSED`. These are exact identities and
arbitrary-degree structural facts; the completeness of the $a_2=0$ reduction is
the written case split of §4 plus the collapse of §3, not a bounded-degree
search.
