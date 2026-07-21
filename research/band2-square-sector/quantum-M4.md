# Quantum $M4$: no genuine band-2 Weyl pair has non-shifted-square top

**INDEPENDENTLY DERIVED AND MACHINE-VERIFIED — NOT PEER REVIEWED — BAND-SCOPED**

QUANTUM MIRROR assault, the non-shifted-square sector. This memo **proves** the
quantum analogue of the classical band-2 non-square theorem
([`../band2-classical-proved/M4_proof_memo.md`](../band2-classical-proved/M4_proof_memo.md)
§5, commit `f637b1a`). It is the missing brick identified by the $a_2=0$
reduction ([`quantum-a2-zero.md`](quantum-a2-zero.md) §7): the one-sided
non-shifted-square branch that the two-sided obstruction **J3q**
([`../band2-j3-provisional/theorem_J3_band2_nonsquare.md`](../band2-j3-provisional/theorem_J3_band2_nonsquare.md),
commit `5a76673`, classical audit folded in `ad43ab5`) does not reach.

## Verdict

> **Theorem $M4q$.** Work in the Weyl algebra $A_1$ over $\mathbb C$ with the
> genuine (falling-factorial) membership below. There is **no** pair
> $[D,X]=1$ with
> $$X=\sum_{k=-2}^{2}x^ka_k(E),\qquad D=\sum_{l=-2}^{2}x^lb_l(E),\qquad a_k,b_l\in\mathbb C[E],$$
> both ladder supports contained in $[-2,2]$, all genuine $A_1$ memberships
> holding, whose top coefficient $a_2(E)$ is nonzero and **outside the
> shifted-square class** — i.e. $a_2\ne c\,h(E)h(E+1)$ for every $c\in\mathbb C^*$
> and $h\in\mathbb C[E]$.
>
> Equivalently: **every** genuine band-2 Weyl pair $[D,X]=1$ with $a_2\ne0$ has
> $a_2$ a shifted square.

The argument places **no** hypothesis on $a_{-2}$ (it may vanish or not; it need
not be a shifted square). Thus $M4q$ covers the genuinely one-sided case
$a_{-2}=0$ **and**, as a by-product, re-derives the $a_2$-side of the two-sided
J3q obstruction; where J3q closes its case by the reflected symmetry between the
two extremes, $M4q$ closes uniformly by **genuine Weyl membership at the negative
extreme** — exactly as polynomiality did in the classical proof.

Everything below is checked exactly by
[`verify_quantum_M4.py`](verify_quantum_M4.py) (ends
`ALL QUANTUM M4 CHECKS PASSED`).

## 1. Setup and conventions

Work in the Ore localization $A_1[x^{-1}]=\bigoplus_{k\in\mathbb Z}x^k\mathbb C[E]$,
$E=x\partial$, with
$$(x^af(E))(x^bg(E))=x^{a+b}f(E+b)g(E),\qquad \partial=x^{-1}E,\qquad
f^{[r]}(E):=f(E+r).$$
For band-2 $X=\sum_{-2}^2x^ka_k$, $D=\sum_{-2}^2x^lb_l$ the ladder-$m$ coefficient
of $[D,X]=DX-XD$ is
$$Q_m=\sum_{k+l=m}\bigl(b_l^{[k]}a_k-a_k^{[l]}b_l\bigr),\qquad [D,X]=1\iff Q_m=\delta_{m0}.$$
(`verify` §0 checks $Q_m$ against the direct crossed-product commutator for every
$m\in[-4,4]$.)

**Genuine membership (quantum polynomiality).** $x^{-r}c(E)\in A_1$ iff the
falling factorial $E^{\underline r}=E(E-1)\cdots(E-r+1)$ divides $c(E)$. In band 2
this is
$$E\mid a_{-1},\,b_{-1},\qquad E(E-1)\mid a_{-2},\,b_{-2}.$$
This is the deformation of the classical $\tau^{\,r}\mid a_{-r}$; it is the one
place the theorem uses that we are in $A_1$ and not merely in the localization.

**Shifted-square class.** A nonzero $q\in\mathbb C[E]$ is *shifted-square* if
$q=c\,h(E)h(E+1)$ with $c\in\mathbb C^*$, $h\in\mathbb C[E]$. Over $\mathbb C$ every
nonzero **constant** is shifted-square ($c=c\cdot1\cdot1$), so a
non-shifted-square $a_2$ has $\deg a_2\ge1$. We use $\deg a_2\ge1$ repeatedly;
it is a consequence of the hypothesis, not an extra assumption.

**Allowed normalization.** Only $D\mapsto D-\lambda X$ ($\lambda\in\mathbb C$) is
used. It preserves $[D,X]=1$ (since $[X,X]=0$), band-2 support, the generated
subalgebra, and — because both $b_{-1},a_{-1}$ carry $E$ and both $b_{-2},a_{-2}$
carry $E(E-1)$ — genuine membership of the gauged $b_{-1},b_{-2}$ (`verify` §1).
We henceforth work in the **gauge $b_2=0$** and rename the gauged coefficients
$b_l$; the memberships $E\mid b_{-1}$ and $E(E-1)\mid b_{-2}$ still hold.

## 2. The descent (mirrors classical M4 §5)

Assume, for contradiction, a genuine pair as in the Theorem, gauged to $b_2=0$.

### $Q_4$ — top proportionality (the gauge)

$Q_4=b_2^{[2]}a_2-a_2^{[2]}b_2=0$ reads $b_2^{[2]}/a_2^{[2]}=b_2/a_2$, so the
rational $b_2/a_2$ is $2$-periodic, hence constant (rational periodicity lemma —
a pole would have infinitely many translates): $b_2=\lambda a_2$. Subtracting
$\lambda X$ realizes the gauge $b_2=0$. *(This is the source of the $\lambda$ the
classical proof carries; in the gauge it is absorbed.)*

### $Q_3$ — homogeneous ladder 3 kills $b_1$ (uses non-shifted-square)

In the gauge, $Q_3=b_1^{[2]}a_2-a_2^{[1]}b_1$ (`verify` §2), so $Q_3=0$ is exactly
the **J2q homogeneous equation** with $h=b_1$:
$$h(E+2)\,a_2(E)=a_2(E+1)\,h(E).$$

> **Lemma J2q** (quantum shifted-square lemma; audited J3 memo, `5a76673`,
> re-verified `ad43ab5`). For $a_2\ne0$ the polynomial solutions are $h=0$, or
> $a_2=c\,h(E)h(E+1)$ (after rescaling $h$), $c\in\mathbb C^*$.

*Re-verification of J2q (`verify` §2).* For $h\ne0$ set
$r=a_2/(h\,h^{[1]})\in\mathbb C(E)$. The equation is **exactly** $r^{[1]}=r$: as an
identity of rational functions,
$r^{[1]}-r=-\bigl(h^{[2]}a_2-a_2^{[1]}h\bigr)/(h\,h^{[1]}h^{[2]})$. Rational
periodicity forces $r$ to be a nonzero constant $c$, i.e. $a_2=c\,h\,h^{[1]}$. The
converse ($a_2$ shifted-square $\Rightarrow$ a nonzero $h$ solves) is also
checked. $\square$

Since $a_2$ is **non-shifted-square**, the case $a_2=c\,h\,h^{[1]}$ is excluded, so
the only solution is $h=b_1=0$.

### $Q_2$ — ladder 2 makes $b_0$ constant

With $b_2=b_1=0$, $Q_2=a_2\,(b_0^{[2]}-b_0)$ (`verify` §3). As $\mathbb C[E]$ is an
integral domain and $a_2\ne0$, $b_0^{[2]}=b_0$: a $2$-periodic polynomial. If
$\deg b_0=d\ge1$ the leading coefficient of $b_0^{[2]}-b_0$ is $2d\cdot\mathrm{lc}(b_0)\ne0$
(degree exactly $d-1$), so $d-1$ would be the degree of $0$ — impossible. Hence
$b_0=\gamma\in\mathbb C$ constant.

### $Q_1$ — ladder 1 kills $v:=b_{-1}$ (telescoping invariant)

With $b_2=b_1=0$ and $b_0=\gamma$ constant, all $a_1,a_0,a_{-1},a_{-2}$ terms in
$Q_1$ cancel, leaving (`verify` §4)
$$Q_1=v^{[2]}a_2-a_2^{[-1]}v,\qquad v=b_{-1},$$
i.e. the staggered-shift homogeneous equation $a_2(E)\,v(E+2)=a_2(E-1)\,v(E)$.
This is the quantum analogue of the classical $2a_2v'+a_2'v=0$. Introduce the
**triple-product invariant**
$$\boxed{\,P(E):=a_2(E-1)\,v(E)\,v(E+1)\,}$$
(the J3q-style $G$-invariant). One has the exact identity (`verify` §4)
$$P(E+1)-P(E)=v(E+1)\bigl(a_2(E)v(E+2)-a_2(E-1)v(E)\bigr)=v^{[1]}\cdot Q_1.$$
Thus $Q_1=0\Rightarrow P$ is $1$-periodic $\Rightarrow$ $P$ constant (leading
coefficient of $P^{[1]}-P$ is $\deg P\cdot\mathrm{lc}(P)$, nonzero unless
$\deg P=0$).

Enumerate the two branches explicitly:
- **$v=0$:** done.
- **$v\ne0$:** in the integral domain $\mathbb C[E]$ degrees add, so
  $\deg P=\deg a_2+2\deg v\ge\deg a_2\ge1$ and $P\ne0$; a nonzero constant $P$
  forces each factor to be a unit, in particular $a_2(E-1)$ constant, i.e. $a_2$
  constant — contradicting $\deg a_2\ge1$.

Hence $v=b_{-1}=0$. *(No membership is used here, matching the classical proof,
where $v=0$ preceded the polynomiality step.)*

### $Q_0=1$ — the central equation forces $a_2$ **linear**, $w$ a nonzero constant

With $v=0$ and $b_1=b_2=0$, the central equation is (`verify` §5)
$$Q_0=w^{[2]}a_2-a_2^{[-2]}w=1,\qquad w=b_{-2}.$$
Set $\Pi(E):=a_2(E-2)\,w(E)$. Then $\Pi(E+2)=a_2(E)\,w(E+2)$, so **exactly**
$$\Pi(E+2)-\Pi(E)=1.$$
This is the quantum analogue of the classical integral $a_2w=\tau/2+\gamma$: a
falling-factorial (here: $a_2$-shifted-product) divisibility. If $\deg\Pi=d$, the
leading coefficient of $\Pi^{[2]}-\Pi$ is $2d\cdot\mathrm{lc}(\Pi)$ with degree
exactly $d-1$ (for $d\ge1$); matching the constant $1$ forces $d=1$ and
$\mathrm{lc}(\Pi)=\tfrac12$. The branches $\Pi$ constant or $\Pi=0$ give
$0\ne1$ and are excluded. So $\Pi$ is exactly linear:
$$a_2(E-2)\,w(E)=\tfrac12E+\beta.$$
Degrees add in $\mathbb C[E]$: $\deg a_2+\deg w=1$ with $\deg a_2\ge1$, so
$$\boxed{\deg a_2=1\ (a_2\ \text{LINEAR}),\qquad \deg w=0,\ w\ne0\ (\text{nonzero constant}).}$$

### $Q_{-1}\dots Q_{-4}$ — the negative tail (consistency)

With $v=0$ and $w$ a nonzero constant (so $w^{[r]}=w$), the remaining equations
reduce to (`verify` §6)
$$Q_{-1}=w(a_1-a_1^{[-2]}),\quad Q_{-2}=w(a_0-a_0^{[-2]}),\quad
Q_{-3}=w(a_{-1}-a_{-1}^{[-2]}),\quad Q_{-4}=w(a_{-2}-a_{-2}^{[-2]}).$$
As $w\ne0$, each is a $2$-periodicity statement, forcing $a_1,a_0,a_{-1},a_{-2}$
all constant. This confirms the reduced system is consistent but is **not needed**
for the contradiction, which arrives at the next line.

### Membership — the kill

The gauged $b_{-2}=w$ is a **nonzero constant**, yet genuine membership demands
$E(E-1)\mid b_{-2}$. A nonzero constant is not divisible by $E(E-1)$
(degree $0<2$; `verify` §7). Contradiction. Therefore no such pair exists.
$\qquad\blacksquare$

*(The classical proof reached the same contradiction, $\tau^2\mid b_{-2}=w$; there
too the negative tail forcing $a_1,\dots,a_{-2}$ constant, and $E\mid a_{-1}$,
$E(E-1)\mid a_{-2}\Rightarrow a_{-1}=a_{-2}=0$, are corollaries, not links in the
contradiction chain.)*

## 3. Rigor ledger

- **Completeness.** The descent is the exact sequential use of
  $Q_4,Q_3,Q_2,Q_1,Q_0$ (plus, redundantly, $Q_{-1}\dots Q_{-4}$); no
  bounded-degree ansatz or solver-exhaustiveness enters the proof. §4 (below) is
  corroboration only.
- **No unwarranted degree equalities.** Every "$n$-periodic polynomial is
  constant" step is discharged by an explicit leading-coefficient computation
  ($n$-step difference of a degree-$d$ polynomial has degree exactly $d-1$ and
  leading coefficient $n\,d\cdot\mathrm{lc}$, nonzero for $d\ge1$). Every degree
  claim about a **product** ($P$ in $Q_1$, $\Pi$ in $Q_0$) uses degree additivity
  in the integral domain $\mathbb C[E]$, where no cancellation is possible.
- **Every zero branch enumerated.** $Q_1$: the branches $v=0$ and $v\ne0$ are
  split explicitly; the $v\ne0$ branch is closed by $P\ne0$. $Q_0$: the branches
  $\Pi=0$, $\deg\Pi=0$, $\deg\Pi\ge1$ are split; only $\deg\Pi=1$ survives, and it
  yields the contradiction.
- **Where each hypothesis is used.** Non-shifted-square $a_2$: at $Q_3$ (J2q gives
  $b_1=0$) and via $\deg a_2\ge1$ (at $Q_1$ and $Q_0$). Genuine membership: only
  at the final kill $E(E-1)\mid b_{-2}=w$.
- **Scope.** A theorem about the fixed support-contained band-2 sector with
  $a_2\ne0$ non-shifted-square. It does **not** claim DC1/JC2, does not classify
  the shifted-square sector, and imposes nothing on $a_{-2}$.

## 4. Machine verification (`verify_quantum_M4.py`)

- **§0** crossed-product engine; $Q_m$ equals the direct commutator for all
  $m\in[-4,4]$; explicit closed forms of $Q_4,Q_0$.
- **§1** gauge legitimacy $[D-\lambda X,X]_m=[D,X]_m-\lambda[X,X]_m$, $[X,X]_m=0$;
  membership survives the gauge.
- **§2** $Q_3=b_1^{[2]}a_2-a_2^{[1]}b_1$; J2q engine as an exact rational identity
  ($r^{[1]}-r$ vs the equation residual) plus the shifted-square converse.
- **§3** $Q_2=a_2(b_0^{[2]}-b_0)$; $2$-step-difference leading coefficient
  $=2d\,\mathrm{lc}$ (degrees $1$–$4$).
- **§4** $Q_1=v^{[2]}a_2-a_2^{[-1]}v$; the invariant identity
  $P^{[1]}-P=v^{[1]}Q_1$; $1$-step-difference leading coefficient $=d\,\mathrm{lc}$;
  product-degree additivity $\deg P=\deg a_2+2\deg v$.
- **§5** $Q_0=\Pi^{[2]}-\Pi$, $\Pi=a_2^{[-2]}w$; $2$-step leading coefficient;
  product-degree additivity $\deg\Pi=\deg a_2+\deg w$.
- **§6** negative tail $Q_{-1..-4}=w(a_j-a_j^{[-2]})$.
- **§7** membership contradiction (no nonzero constant is divisible by $E(E-1)$).
- **§8** bounded-degree **exhaustive emptiness** (Gröbner basis $=\langle1\rangle$)
  for concrete non-shifted-square tops $a_2\in\{E,\,2E-3,\,E^2{+}1,\,E^2{-}E{+}1,\,
  E^3,\,E^4{+}1\}$ with all other coefficients generic of degree $\le2$ and $\le3$
  and full membership imposed — all EMPTY. Corroboration only; completeness is §2.
- **§9** positive control: the explicit **genuine** pair
  $X=(x+\partial)^2-\tfrac12\partial$, $D=2x+2\partial$ satisfies $[D,X]=1$ with
  every membership holding and $a_2=1$ **shifted-square** — so the $Q_m$ machinery
  does accept real pairs, and the emptiness for non-shifted-square $a_2$ is a
  genuine distinction driven by membership, not an artifact.

Run from the repository root:
```sh
uv run --with sympy python research/band2-square-sector/verify_quantum_M4.py
```

## 5. Contrast with the classical template, and place in the program

| step | classical M4 (`f637b1a`, §5) | quantum $M4q$ (here) |
|---|---|---|
| $Q_4$ | $(b_2/a_2)'=0\Rightarrow b_2=\lambda a_2$ | $b_2/a_2$ $2$-periodic $\Rightarrow$ constant |
| $Q_3$ | $2a_2u'-a_2'u=0$, nonsquare $\Rightarrow u=0$ | J2q $h^{[2]}a_2=a_2^{[1]}h$, non-shifted-sq $\Rightarrow b_1=0$ |
| $Q_1$ | $2a_2v'+a_2'v=0\Rightarrow(a_2v^2)'=0$ | $P=a_2^{[-1]}v\,v^{[1]}$ is $1$-periodic $\Rightarrow$ const |
| $Q_0$ | $2(a_2w)'=1\Rightarrow a_2w=\tau/2+\gamma$ | $\Pi^{[2]}-\Pi=1$, $\Pi=a_2^{[-2]}w\Rightarrow a_2$ linear |
| kill | $\tau^2\mid b_{-2}=w$ nonzero const | $E(E-1)\mid b_{-2}=w$ nonzero const |

The mechanism is a faithful mirror: the classical first integral $(a_2v^2)$ becomes
the telescoping product $a_2^{[-1]}v\,v^{[1]}$, and the classical antiderivative
$a_2w=\tau/2+\gamma$ becomes the staggered difference $\Pi^{[2]}-\Pi=1$ with
$\Pi=a_2(E-2)w$. As in [`quantum-mirror.md`](quantum-mirror.md) §5, the difference
setting is if anything more rigid: the leading-coefficient bookkeeping of the
$n$-step differences replaces differentiation, and the same nonzero-constant
argument closes both product invariants.

**Place in the band-2 assembly.** Combined with:
- the $a_2\ne0$ **shifted-square** sector (constant-$h$ resistant branch closed in
  [`quantum-mirror.md`](quantum-mirror.md); the rest tracked in
  `quantum-completion.md`), and
- the $a_2=0$ **reduction** ([`quantum-a2-zero.md`](quantum-a2-zero.md), routing
  every $a_2=0$ pair to the square sector or to band-1 rigidity P3),

$M4q$ supplies the non-shifted-square leaf the reduction feeds (its §7 remainder):
the one-sided $a_2\ne0$ pair whose top need not be shifted-square, hence not
covered by J3q. With this leaf proved, the only openness inherited by the full
support-contained quantum band-2 theorem is P3 (band-1 rigidity, provisional) and
the non-constant-$h$ shifted-square sector — not the non-square sector, which is
now closed.

## 6. Status of claims

- **Proved, arbitrary degree (structural, not degree-bounded):** the descent of
  §2 and the Theorem $M4q$. The non-shifted-square band-2 sector with $a_2\ne0$
  contains no genuine $A_1$ pair, for **any** $a_{-2}$.
- **Consumed (audited) inputs:** the rational periodicity lemma and Lemma J2q
  (J3 memo `5a76673`/`ad43ab5`); the crossed-product $Q_m$ system (re-derived and
  machine-checked here, `verify` §0).
- **Corroboration only:** the bounded Gröbner sweep (`verify` §8) and the positive
  control witness (`verify` §9).
- **Out of scope:** the shifted-square sector; band-1 (P3); any global DC1/JC2
  claim; non-closed coefficient fields (over $\mathbb C$ non-shifted-square forces
  $\deg a_2\ge1$, which the proof uses).
- **No counterexample exists** in the non-shifted-square sector (it is empty).
