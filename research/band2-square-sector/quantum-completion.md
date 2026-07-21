# Quantum band-2 shifted-square sector: complete constant-$h$ classification

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
- $h$ is nonconstant — **impossible** (§4);
- $\kappa=0$ — **impossible** (§5);
- $\kappa\neq0$, $h$ constant (normalize $h=1$), and $b_{-2}\neq0$ — **impossible**
  (§6 for $a_{-2}=0$; the case $a_{-2}\neq0$ is `quantum-mirror.md`);
- $\kappa\neq0$, $h=1$, $b_{-2}=0$ — the pair is a member of the **tame family**
  of §8 (§7).

Equivalently: **the quantum band-2 shifted-square sector consists exactly of the
tame automorphism family of §8.** Combined with `quantum-mirror.md` this is a
complete, unconditional classification *of this sector*.

**What is NOT proved here (the remainder).** The following are inputs or lie
outside this sector; none is silently used:
1. that a nonzero extreme coefficient $a_2$ *is* a shifted square (the quantum
   analogue of the classical M4 theorem / the J3 obstruction) — an **input**,
   under audit;
2. the quantum **band-1** classification, needed to assemble an *orientation*
   argument (pair-exchange / reflection) up to full band 2 — an **input**, under
   audit (P3); §9 states the assembly conditionally on 1–2;
3. the sectors with $a_2=0$ on the $X$ side (routed through lower bands);
4. the case $\kappa=0$ is disposed of here directly, not assumed.

So the headline "every band-2 quantum pair is a tame automorphism" is **only
conditional** (on inputs 1–2). Unconditional is the sector classification above.

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

Partition by $(\deg h,\ \kappa,\ s,\ w)$ into: **d** ($h$ nonconstant), **a**
($\kappa=0$), and — for $\kappa\neq0,h=1$ — the negative-tail cases split by
whether $s,w$ vanish: **b** ($w=0,s\neq0$), **c** ($s=0,w\neq0$), the resistant
**$*$** ($s\neq0,w\neq0$, i.e. $w=\mu s,\mu\neq0$), and the base ($s=w=0$).

## 2. The general-$h$ central integral (engine for §§4–5)

**Lemma (general-$h$ central telescoping).** With $a_2=hh^{[1]}$, $a_1=hp$,
$b_1=\kappa h$, $b_2=0$,
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

## 4. Branch d: nonconstant $h$ is impossible

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

## 5. Branch a: $\kappa=0$ is impossible

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
- $h$ nonconstant: as in §4, $h$ affine and $h\,w^{[1]}+h^{[-2]}w=(E+c)/h^{[-1]}$
  is a *constant*; but its degree is $1+\deg w\geq3$ with leading coefficient
  $2\,\mathrm{lc}(h)\,\mathrm{lc}(w)\neq0$ (no cancellation — same sign):
  contradiction.
- $h=1$: $w^{[1]}+w=E$ has the unique solution $w=\tfrac E2-\tfrac14$, and
  $w(0)=-\tfrac14\neq0$ violates $E(E-1)\mid w$.

**Branch a is empty (all $h$).** (Quantum mirror of `classical-Astar.md` §7 and
[commit `91a053a`] §7.)

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

## 7. Branch c ($s=0$, $w\neq0$): impossible

$w=c\,u\,u^{[-1]}$ (§3) and $w\neq0$ give $c\neq0$ and $u\neq0$.
- If $v\neq0$: $Q_{-2}=v^{[-1]}u-u^{[-1]}v=0$ makes $u/v$ $1$-periodic, so
  $u=\lambda v$. Substituting into $(\star)$ (with $s=0$, $w=cuu^{[-1]}$),
  the left side factors as $v\cdot\bigl[c\lambda^2(v^{[1]}+v^{[-1]})+p^{[-1]}
  -\kappa\lambda\bigr]$, so the equation reads $v\cdot(\text{poly})=E$: hence
  $v\mid E$, and with $E\mid v$, **$v=v_1E$ is linear.** Then $u,w,p,B,a_0$ are
  all of bounded degree; the residual finite system $Q_0{=}1,Q_{-1}{=}Q_{-2}{=}0$
  has no solution with $w\neq0$ [machine-checked, Gröbner + saturation].
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
when $c_1=0$ it is the base ($s=w=0$). As $U=x+c_0+c_1\partial$ is an affine
symplectic generator of $A_1$, $(X,D)$ generate $A_1$: the pair is a **tame
automorphism image**.

## 9. Assembly

**Theorem (quantum shifted-square sector, complete — unconditional).** Let
$[D,X]=1$ in $A_1[x^{-1}]$ with support $[-2,2]$, $a_2=h\,h^{[1]}\neq0$, and
genuine Weyl membership. Then $h$ is constant, $\kappa\neq0$, $b_{-2}=0$, and the
pair is — after the reversible normalizations (gauge $D\mapsto D-\lambda X$,
diagonal scaling to $h=1$, additive constants) — a member of the tame family
§8. In particular $(X,D)$ generate $A_1$.

*Proof.* §4 (⇒ $h$ constant), §5 (⇒ $\kappa\neq0$), §7 and `quantum-mirror.md`
(⇒ $b_{-2}=0$: the cases $a_{-2}=0$ and $a_{-2}\neq0$), §6 and the base case
(⇒ tame). $\qquad\blacksquare$

**Conditional extension to full band 2.** Following the classical assembly
[full classical band-2 theorem, commit `f637b1a`], an unoriented band-2 pair is
reduced to $a_2\neq0$ by pair-exchange $(F,G)\mapsto(G,-F)$ and the quantum
reflection $\psi:x\mapsto x^{-1},E\mapsto-E$ (an automorphism of $A_1[x^{-1}]$,
so $[\psi D,\psi X]=1$; $\psi$ sends the index-$(-2)$ coefficient to index $2$).
*Conditional on* (i) the quantum band-1 classification and (ii) the quantum
shifted-square input "$a_2\neq0\Rightarrow a_2=hh^{[1]}$" — both under audit (P3)
— every band-2 quantum Keller pair with genuine membership is a tame automorphism
image. This memo does **not** prove (i) or (ii); it proves the sector Theorem
above, which is the quantum counterpart of the classical square-sector closure.

## 10. Verification (`verify_quantum_completion.py`)

Exact SymPy checks: the general-$h$ central integral (§2); both shifted-product
invariants (§3); the branch-d membership kill (§4); the branch-a cascade,
telescoping, and $h=1$/nonconstant-$h$ contradictions (§5); the branch-c
$u=\lambda v$, $v\mid E$ reduction, the $v=0$ sub-case $2cu_1^2=0$, and Gröbner
emptiness (§7); the branch-b leading-coefficient identities
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
- **Regression only:** the bounded Gröbner sweeps (branch c; resistant).
- **Input / under audit (NOT proved here):** the shifted-square reduction
  "$a_2\neq0\Rightarrow a_2=hh^{[1]}$" and the quantum band-1 classification;
  the full-band-2 assembly of §9 is conditional on these.
- **Out of scope:** the $a_2=0$ sectors.
- No counterexample exists in this sector: every pair is the tame family §8.
