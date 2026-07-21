# Quantum band-2 shifted-square sector: complete classification

**INDEPENDENTLY DERIVED AND MACHINE-VERIFIED — NOT PEER REVIEWED**

QUANTUM MIRROR assault, M5 campaign. This memo closes **all** remaining
sub-branches of the quantum band-2 shifted-square sector, complementing the
resistant branch closed in `quantum-mirror.md` (commit `ad43ab5`, endorsed by
upstream audit). Every displayed identity and every branch closer is checked
exactly by `verify_quantum_completion.py` (ends
`ALL QUANTUM COMPLETION CHECKS PASSED`).

## 0. Precise scope — what this memo proves, and what it does not

Work in $A_1[x^{-1}]=\bigoplus_{k\in\mathbb Z}x^k\mathbb C[E]$,
$(x^af(E))(x^bg(E))=x^{a+b}f(E+b)g(E)$, $f^{[r]}(E)=f(E+r)$, $E=x\partial$,
$\partial=x^{-1}E$. For band-2 $X=\sum_{-2}^2x^ka_k$, $D=\sum_{-2}^2x^lb_l$, the
ladder-$m$ coefficient of $[D,X]$ is
$Q_m=\sum_{k+l=m}(b_l^{[k]}a_k-a_k^{[l]}b_l)$, and $[D,X]=1\iff Q_m=\delta_{m0}$.
Genuine Weyl membership means $E\mid a_{-1},b_{-1}$ and $E(E-1)\mid a_{-2},b_{-2}$.

**Sector hypothesis (input).** $a_2=h\,h^{[1]}\neq0$ — a *shifted square*. This is
the defining hypothesis of the sector, taken from the audited partial cascade
`../band2-m5-partial/quantum-shifted-square-sector-partial.md` (commit
`91a053a`), whose ladder-4/3/2 reductions ($b_2\to0$ by the gauge $D\mapsto
D-\lambda X$; $b_1=\kappa h$; for $\kappa\neq0$, $h\mid a_1$, $a_1=hp$, and the
quantum midpoint $B^{[1]}+B=\kappa p+\gamma$) are used verbatim.

**What is proved here (unconditional, arbitrary degree).** For a pair
$[D,X]=1$ in this sector with genuine membership, **exactly one** of:
- $\kappa=0$ — **impossible for every $h$** (§4);
- $\kappa\ne0$ and $h$ is nonconstant — **impossible** (§5);
- $\kappa\neq0$, $h$ constant (normalize $h=1$), and $b_{-2}\neq0$ — **impossible**
  (§6 for $a_{-2}=0$; the case $a_{-2}\neq0$ is `quantum-mirror.md`);
- $\kappa\neq0$, $h=1$, $b_{-2}=0$ — the pair is a member of the **tame family**
  of §8 (§7).

Equivalently: **the quantum band-2 shifted-square sector consists exactly of the
tame automorphism family of §8.** Combined with `quantum-mirror.md` this is a
complete, unconditional classification *of this sector*.

**What is not proved in this memo.** The following inputs or complementary
sectors are proved in the cited companion memos; none is silently assumed in the
sector argument here:
1. the shifted-square reduction for a nonzero extreme coefficient,
   $a_2\neq0\Rightarrow a_2=h h^{[1]}$, is M4q (`quantum-M4.md`);
2. the genuine-$A_1$ band-1 classification is Theorem 1 of
   `quantum-band1.md`;
3. the sectors with $a_2=0$ on the $X$ side are routed by
   `quantum-a2-zero.md` through Fourier reflection, pair exchange, or band 1;
4. the case $\kappa=0$ is disposed of here directly, not assumed.

Thus this memo proves the shifted-square sector unconditionally. The full
fixed-band-2 theorem follows only after combining it with those repaired inputs,
as done in `quantum-band2-theorem.md`; it is not an additional theorem proved by
this memo alone.

## 1. Reduced data and the branch partition

After the audited positive cascade, in the gauge $b_2=0$, $b_1=\kappa h$: write
$p=a_1/h$ (when $\kappa\neq0$), $u=a_{-1}$, $s=a_{-2}$, $B=b_0$, $v=b_{-1}$,
$w=b_{-2}$. The three exact integrals established in `quantum-mirror.md` (commit `ad43ab5`)
hold at $h=1$: the quantum midpoint $B^{[1]}+B=\kappa p+\gamma$; the ladder-1
integral $\kappa a_0=(v^{[1]}+v)+\tfrac1\kappa(B^2-\gamma B)+A$; and the central
telescoping
$$w^{[1]}+w+p^{[-1]}v-\kappa u=E,\tag{$\star$}$$
whose constant is $0$ by the point conditions $s(0)=s(1)=v(0)=0$. The remaining
equations are $Q_{-1}=Q_{-2}=Q_{-3}=0$ ($Q_{-4}=0$ gives $w^{[-2]}s=s^{[-2]}w$).

Partition first by $\kappa$. Branch **a** is $\kappa=0$ for arbitrary $h$ and is
disposed of in §4. Only on the complementary branch $\kappa\neq0$ may one use
the audited deduction $a_1=hp$ and the general-$h$ integral (2.1). There branch
**d** is nonconstant $h$ (§5); after it is excluded, normalize constant $h$ to
$1$ and split by the negative tails: **b** ($w=0,s\neq0$), **c**
($s=0,w\neq0$), resistant **$*$** ($s,w\neq0$, hence $w=\mu s$), and the base
($s=w=0$).

## 2. The general-$h$ central integral (only after $\kappa\ne0$)

**Lemma (general-$h$ central telescoping, $\kappa\ne0$).** With
$a_2=hh^{[1]}$, $a_1=hp$ (the latter deduction is available here because
$\kappa\ne0$), $b_1=\kappa h$, $b_2=0$,
$$Q_0=(T-1)\,G,\qquad G=h^{[-1]}\bigl(h\,w^{[1]}+h^{[-2]}w+p^{[-1]}v-\kappa u\bigr).$$
[Exact identity, machine-checked.] Hence $Q_0=1$ gives $G=E+c$, i.e.
$$h^{[-1]}\bigl(h\,w^{[1]}+h^{[-2]}w+p^{[-1]}v-\kappa u\bigr)=E+c.\tag{2.1}$$
At $h=1$ this is $(\star)$. This is the quantum transfer of the classical
integrated identity $(2h^2w+hpv-\kappa h\,a_{-1})'=1$ [oriented classical square
sector, commit `44c66d5`]; unlike the classical case it has no one-line closed
form for general $h$, but (2.1) is all §§4–5 need.

## 3. Two shifted-product invariants

From $Q_{-3}=0$, specialized:
- $w=0$: $Q_{-3}=s\,v^{[-2]}-s^{[-1]}v$. If $v\neq0$, the ratio $s/(v\,v^{[-1]})$
  is $1$-periodic, hence constant: $\boxed{s=c\,v\,v^{[-1]}}$ (§6).
- $s=0$: $Q_{-3}=w^{[-1]}u-u^{[-2]}w$. If $u\neq0$, likewise
  $\boxed{w=c\,u\,u^{[-1]}}$ (§7).

[Both solved forms machine-checked.] These are the quantum shifted-product
analogues of the classical $s=cv^2$ (branch B0) and $a_{-1}^2=cw$ (branch A0) of
`classical-Astar.md` (commit `ad43ab5`).

## 4. Branch a: $\kappa=0$ is impossible for all $h$

Here $b_1=\kappa h=0$ and $b_2=0$. $Q_2=a_2(B^{[2]}-B)=0$ with $a_2\neq0$ forces
$B$ $2$-periodic, hence $B=\beta$ constant. $Q_1=h\,(h^{[1]}v^{[2]}-h^{[-1]}v)=0$
with $h\neq0$ forces $h^{[-1]}v$ to be $2$-periodic, hence a constant $\nu$.
Enumerating: if $v\neq0$ and $h$ is nonconstant, $h(E-1)v(E)=\nu$ forces both
factors constant, impossible; if $h=1$, $v=\nu$ and $E\mid v$ give $\nu=0$. In
all cases $v=0$.

With $v=0$, $b_1=b_2=0$, $B$ constant, $Q_0=a_2w^{[2]}-a_2^{[-2]}w=(T-1)\bigl[
h^{[-1]}(h\,w^{[1]}+h^{[-2]}w)\bigr]=1$, so
$h^{[-1]}(h\,w^{[1]}+h^{[-2]}w)=E+c$. If $w=0$ then $Q_0=0\neq1$; so $w\neq0$ and
$E(E-1)\mid w$ gives $\deg w\geq2$.
- $h$ nonconstant: by the same degree/divisibility argument as §5, $h$ affine and $h\,w^{[1]}+h^{[-2]}w=(E+c)/h^{[-1]}$
  is a *constant*; but its degree is $1+\deg w\geq3$ with leading coefficient
  $2\,\mathrm{lc}(h)\,\mathrm{lc}(w)\neq0$ (no cancellation — same sign):
  contradiction.
- $h=1$: $w^{[1]}+w=E$ has the unique solution $w=\tfrac E2-\tfrac14$, and
  $w(0)=-\tfrac14\neq0$ violates $E(E-1)\mid w$.

**Branch a is empty (all $h$).** (Quantum mirror of `classical-Astar.md` §7 and
[commit `91a053a`] §7.)

## 5. Branch d ($\kappa\ne0$): nonconstant $h$ is impossible

From (2.1), $h^{[-1]}=h(E-1)$ divides the degree-$1$ polynomial $E+c$. A
polynomial of degree $\geq1$ dividing a degree-$1$ polynomial has degree exactly
$1$; so **$h$ is affine**, $h=\alpha(E-\rho)$ with $\alpha\neq0$. Then
$h^{[-1]}\mid(E+c)$ with both degree $1$ forces the quotient to be the constant
$1/\alpha$, so from (2.1)
$$\kappa u=h\,w^{[1]}+h^{[-2]}w+p^{[-1]}v-\tfrac1\alpha .$$
Evaluate at $E=0$. Membership gives $w(0)=w(1)=0$ (from $E(E-1)\mid w$) and
$v(0)=0$; hence $h(0)w(1)=h(-2)w(0)=p(-1)v(0)=0$ and
$$\kappa\,u(0)=-\tfrac1\alpha\neq0.$$
But $E\mid u$ forces $u(0)=0$: contradiction. **Branch d is empty.** (This is
the quantum mirror of the classical one-line kill, `classical-Astar.md` §6.)

## 6. Branch b ($w=0$, $s\neq0$): exactly the tame family

Since $s=c\,v\,v^{[-1]}$ (§3) and $s\neq0$, we have $v\neq0$ and $c\neq0$.
With $w=0$, $u=\tfrac1\kappa(p^{[-1]}v-E)$ from $(\star)$; the remaining equations
are $Q_{-2}(b)=(v^{[-1]}u-u^{[-1]}v)+s(B^{[-2]}-B)=0$ and
$Q_{-1}(b)=\kappa(s-s^{[1]})+v(a_0-a_0^{[-1]})+u(B^{[-1]}-B)=0$. Set $P=\deg p$.

- **$P\geq1$.** The leading coefficient of $Q_{-2}(b)$ (at degree $2\deg
  v+P-1$) is
  $$\mathrm{lc}(Q_{-2}(b))=-\tfrac{P}{\kappa}\,\mathrm{lc}(v)^2\,\mathrm{lc}(p)\,(c\kappa^2-1)$$
  [machine-checked identity]. As $P,\kappa,\mathrm{lc}(v),\mathrm{lc}(p)\neq0$,
  vanishing forces $c=1/\kappa^2$. But then the leading coefficient of
  $Q_{-1}(b)$ (at degree $P$) is $\tfrac{P}{2}\,\mathrm{lc}(p)\neq0$ [machine-checked],
  so $Q_{-1}(b)\neq0$: contradiction. Hence **$P\geq1$ is empty.**
- **$P=0$** ($p$ constant, $B$ constant). Then $s(B^{[-2]}-B)=0$ and
  $Q_{-2}(b)=\bigl[(E-1)v-E\,v^{[-1]}\bigr]/\kappa=0$ [machine-checked], so
  $v/E$ is $1$-periodic, i.e. **$v=v_1E$ is linear** (using $E\mid v$, $v\neq0$).
  The residual finite system $Q_{-1}(b)=Q_{-2}(b)=0$ then forces
  $c=1/\kappa^2$ (the alternative $v_1=0$ is excluded by $s\neq0$)
  [machine-checked].

So branch b is exactly $\{p\ \text{const},\ v=v_1E,\ c=1/\kappa^2\}$, a subfamily
of §8 (with $c_1\neq0$).

### The base branch $s=w=0$

Here $Q_{-2}=v^{[-1]}u-u^{[-1]}v=0$ and the central integral is
$p^{[-1]}v-\kappa u=E$.

If $v\ne0$, the first identity says $u/v$ is a 1-periodic rational function,
hence $u=\lambda v$. The central identity becomes
$v(p^{[-1]}-\kappa\lambda)=E$. Therefore $v\mid E$; membership $E\mid v$ gives
$v=v_1E$ with $v_1\ne0$, and then $p^{[-1]}-\kappa\lambda=1/v_1$, so $p$ is
constant. The midpoint $B^{[1]}+B=\kappa p+\gamma$ makes $B$ constant, and the
exact ladder equation reduces to
$$Q_{-1}=\frac{2Ev_1^2}{\kappa},$$
which cannot vanish in characteristic zero. Thus $v=0$.

The central identity now gives $u=-E/\kappa$. In the ladder-$(-1)$ equation,
the terms involving $v,s,w$ vanish and
$$Q_{-1}=u(B^{[-1]}-B)=-\frac E\kappa(B^{[-1]}-B).$$
Thus $B$ is constant. The midpoint identity then forces $p$ constant, and the
$a_0$ integral makes $a_0$ constant. All nine ladder equations
then hold, and these coefficients are precisely §8 with $c_1=0$:
$U=x+c_0$, $X=U^2-\partial/\kappa-A$, and
$D=\lambda X+\kappa U+\beta$ after restoring the gauge and constants. The
verifier checks every displayed reduction and all ladder coefficients.

## 7. Branch c ($s=0$, $w\neq0$): impossible

$w=c\,u\,u^{[-1]}$ (§3) and $w\neq0$ give $c\neq0$ and $u\neq0$.
- If $v\neq0$: $Q_{-2}=v^{[-1]}u-u^{[-1]}v=0$ makes $u/v$ $1$-periodic, so
  $u=\lambda v$. Substituting into $(\star)$ (with $s=0$, $w=cuu^{[-1]}$),
  the left side factors as $v\cdot\bigl[c\lambda^2(v^{[1]}+v^{[-1]})+p^{[-1]}
  -\kappa\lambda\bigr]$, so the equation reads $v\cdot(\text{poly})=E$: hence
  $v\mid E$, and with $E\mid v$, **$v=v_1E$**, $v_1\ne0$. Dividing the central
  equation by $v_1E$ gives the exact formula
  $$p^{[-1]}=\frac1{v_1}+\kappa\lambda-c\lambda^2v_1(2E),$$
  because $v^{[1]}+v^{[-1]}=2v_1E$. Thus $p$ has degree at most one; there are
  no unexamined higher-degree profiles. The saturated $(\deg p,\deg v)=(1,1)$
  Gröbner calculation for $Q_{-1}=Q_{-2}=0$ is therefore exhaustive; its
  saturation includes $\kappa c\lambda v_1\ne0$, exactly the assumptions
  $\kappa\ne0$, $w=c\lambda^2v v^{[-1]}\ne0$, and $v\ne0$. It yields the unit
  ideal. The additional degree profiles in the verifier are regression checks
  only.
- If $v=0$ (so $u\neq0$): $(\star)$ becomes $w^{[1]}+w-\kappa u=E$, i.e.
  $u\cdot\bigl[c(u^{[1]}+u^{[-1]})-\kappa\bigr]=E$, forcing $u\mid E$ and (with
  $E\mid u$) $u=u_1E$; the $E^2$-coefficient is then $2c\,u_1^2=0$, impossible
  ($c,u_1\neq0$) [machine-checked].

**Branch c is empty.** (Quantum mirror of `classical-Astar.md` §4.)

## 8. The tame family (normal form)

$$\boxed{\;U=x+c_0+c_1\partial,\qquad X=U^2-\tfrac1\kappa\partial-A,\qquad
D=\lambda X+\kappa U+\beta,\;}$$
with $c_0,c_1,A,\kappa,\lambda,\beta\in\mathbb C$, $\kappa\neq0$. Since
$[U,X]=[U,-\tfrac1\kappa\partial]=\tfrac1\kappa$, one has
$[D,X]=\kappa[U,X]=1$ [machine-checked as $Q_m=\delta_{m0}$]. Its reduced
($b_2=0$) coefficients are
$$a_2=1,\quad a_{-2}=c_1^2E(E-1),\quad b_1=\kappa,\quad b_{-1}=\kappa c_1E,\quad
b_{-2}=0,$$
so $w=0$, $s=c_1^2E(E-1)=\tfrac1{\kappa^2}\,v\,v^{[-1]}$ (i.e. $c=1/\kappa^2$),
and all memberships hold [machine-checked]. When $c_1\neq0$ this is branch b;
when $c_1=0$ it is the base ($s=w=0$). There are universal recovery identities
$$U=\frac{D-\lambda X-\beta}{\kappa},\qquad
\partial=\kappa(U^2-X-A),\qquad x=U-c_0-c_1\partial.$$
Thus generation is explicit. More strongly, the family is given by an explicit
tame word. Starting from $(x,\partial)$, apply the affine symplectic map
$(x,\partial)\mapsto(U,V)=(x+c_0+c_1\partial,\partial)$, then the triangular
map $(U,V)\mapsto(U,V-\kappa U^2+\kappa A)$, then the affine symplectic map
$$(U,V')\mapsto\bigl(-V'/\kappa,
\lambda(-V'/\kappa)+\kappa U+\beta\bigr)=(X,D).$$
Each linear part has determinant one, so this is a tame automorphism word, not
merely an inference from generation. The verifier checks the word and the three
recovery identities symbolically for general parameters.

## 9. Assembly

**Theorem (quantum shifted-square sector, complete — unconditional).** Let
$[D,X]=1$ in $A_1[x^{-1}]$ with support $[-2,2]$, $a_2=h\,h^{[1]}\neq0$, and
genuine Weyl membership. Then $h$ is constant, $\kappa\neq0$, $b_{-2}=0$, and the
pair is — after the reversible normalizations (gauge $D\mapsto D-\lambda X$,
diagonal scaling to $h=1$, additive constants) — a member of the tame family
§8. In particular $(X,D)$ generate $A_1$.

*Proof.* Section 4 gives $\kappa\neq0$ for every $h$; only then §5 gives $h$
constant. Section 7 and `quantum-mirror.md` give $b_{-2}=0$ in the cases
$a_{-2}=0$ and $a_{-2}\neq0$, respectively. Section 6 and the base branch then
identify the pair with the tame family. $\qquad\blacksquare$

**Orientation note.** A previous draft suggested the localized involution
$x\mapsto x^{-1},E\mapsto-E$ as an orientation mechanism. That map does not
preserve $A_1$ and is retracted for the present theorem. The genuine orientation
used by the assembly is the Fourier automorphism
$x\mapsto-\partial,\partial\mapsto x$, with its exact falling-factorial
coefficient dictionary proved in `quantum-a2-zero.md`. Together with pair
exchange, M4q, and the audited band-1 theorem, that genuine $A_1$ operation is
what extends this sector result to full fixed band 2.

## 10. Verification (`verify_quantum_completion.py`)

Exact SymPy checks: the general-$h$ central integral (§2); both shifted-product
invariants (§3); the branch-a cascade for all $h$ (§4); the scoped
$\kappa\ne0$ branch-d membership kill (§5); the branch-c
$u=\lambda v$, $v\mid E$ reduction, the $v=0$ sub-case $2cu_1^2=0$, and Gröbner
emptiness (§7); the exact base-branch split and contradiction; the branch-b leading-coefficient identities
$\mathrm{lc}(Q_{-2}(b))$ and $\mathrm{lc}(Q_{-1}(b))|_{c=1/\kappa^2}$, the $P=0$
reduction, and the finite-system forcing of $c=1/\kappa^2$ (§6); a regression
confirmation of the resistant branch; and the tame family $[D,X]=1$ with its
coefficients and memberships (§8). Runs to `ALL QUANTUM COMPLETION CHECKS
PASSED`.

## 11. Status of claims

- **Proved, arbitrary degree (unconditional):** §§2–8 and the sector Theorem of
  §9. Every argument uses exact identities or explicit finite reductions;
  degree claims are stated as the vanishing of a named leading coefficient, and
  zero-polynomial cases ($v=0$, $u=0$, $w=0$) are enumerated explicitly.
- **Exact symbolic checks:** displayed telescoping, coefficient, recovery, and
  tame-word identities.
- **Arbitrary-degree prose:** divisibility/periodicity and leading-degree
  arguments, plus the exhaustive branch partition.
- **Bounded computation:** resistant profiles and the extra branch-c degree
  profiles are regression only. The branch-c $(1,1)$ saturated Gröbner
  calculation is exhaustive only because the preceding exact central equation
  proves $\deg p\le1$ and $\deg v=1$.
- **Companion inputs (not proved here):** M4q's shifted-square reduction and the
  band-1 theorem are proved in `quantum-M4.md` and `quantum-band1.md`;
  `quantum-a2-zero.md` supplies the complementary routing. Their assembly is
  stated and proved in `quantum-band2-theorem.md`.
- **Out of scope for this memo:** direct analysis of the $a_2=0$ sectors.
- No counterexample exists in this sector: every pair is the tame family §8.
