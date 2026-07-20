# Quantum mirror: the resistant constant-$h$ branch of band 2 is empty

**INDEPENDENTLY DERIVED AND MACHINE-VERIFIED — NOT PEER REVIEWED**

QUANTUM MIRROR assault, M5 campaign. This memo **closes** the quantum resistant
branch left open by
`research/band2-m5-partial/quantum-shifted-square-sector-partial.md`
(the "unresolved candidate branch $h=1,\ \kappa\mu s\ne0$", memo §6, §8).

## Verdict

> **Theorem.** In $A_1[x^{-1}]$, in the shifted-square sector of band 2 with
> constant top ($a_2=1$, i.e. $h=1$), after the audited gauge and positive
> cascade, the resistant branch
> $$\kappa\ne0,\quad \mu\ne0,\quad s=a_{-2}\ne0$$
> (with genuine Weyl memberships $E\mid u,v$ and $E(E-1)\mid s,w$) is **EMPTY**.
> No such pair $[D,X]=1$ exists, at any degree.

Consequently this branch yields **no DC1 counterexample candidate**: the sole
resistant locus of the constant-$h$ quantum square sector contains no pair at
all, so there is nothing to test for $A_1$-generation. The two-sided
constant-$h$ square sector reduces to the tame $s=0$ families already recorded
upstream (the swapped-shear and $Z$-generator families), which are genuine
$[D,X]=1$ pairs that **do** generate $A_1$.

This is a genuine advance over the audited partial cascade
(`research/band2-m5-partial/quantum-shifted-square-sector-partial.md`), which
left the entire branch open. The corresponding **classical** resistant branch
**A\*** was likewise left open by the audited
`research/band2-m5-partial/classical-square-sector-partial.md`; it is closed in
parallel, by a different (integrated first-integral) route, in the sibling memo
`classical-Astar.md` (same directory). The quantum closure here rests on an
extra rigidity that integer shifts provide and differentiation does not
(Lemma R, $S=2\deg\phi$); see §5.

Everything below is checked exactly by `verify_quantum_mirror.py`
(ends `ALL QUANTUM MIRROR CHECKS PASSED`).

## 1. Setup and conventions

Work in $A_1[x^{-1}]=\bigoplus_{k\in\mathbb Z}x^k\mathbb C[E]$ with
$(x^af(E))(x^bg(E))=x^{a+b}f(E+b)g(E)$ and $f^{[r]}(E)=f(E+r)$. For band-2
$X=\sum_{-2}^2x^ka_k$, $D=\sum_{-2}^2x^lb_l$ the ladder-$m$ coefficient of
$[D,X]$ is
$$Q_m=\sum_{k+l=m}\bigl(b_l^{[k]}a_k-a_k^{[l]}b_l\bigr),\qquad Q_m=\delta_{m0}.$$
`verify` §0 reproduces the audited memo tail equations (4.1),(4.3)–(4.6),(4.2)
from this $Q_m$ (note (4.1) is $Q_1/h$; all others are $Q_m$ verbatim), and
$Q_4=Q_3=0$ identically in the gauge $b_2=0,\ b_1=\kappa h$.

Branch data ($h=1$, so $a_2=1$; audited positive cascade):
$$a_1=p,\ a_0,\ a_{-1}=u,\ a_{-2}=s,\qquad
b_2=0,\ b_1=\kappa,\ b_0=B,\ b_{-1}=v,\ b_{-2}=w=\mu s,$$
memberships $E\mid u,v$ and $E(E-1)\mid s,w$; $T$ denotes $f\mapsto f^{[1]}$.
Write $P=\deg p$, $V=\deg v$, $S=\deg s\ (\ge2)$, $U=\deg u$, $A_0=\deg a_0$.

## 2. Three exact arbitrary-degree integrals (`verify` §1)

These make the audited "positive cascade" fully explicit as closed identities.

- **Midpoint** (ladder 2): $\;B^{[1]}+B=\kappa p+\gamma.$ Since $T+1$ is
  invertible on $\mathbb C[E]$, $B$ is determined; $\deg B=P$,
  $\operatorname{lc}(B)=\kappa\operatorname{lc}(p)/2$.
- **Ladder-1 integral** (from $Q_1=0$, using
  $p(B^{[1]}-B)=\tfrac1\kappa(T-1)(B^2-\gamma B)$):
  $$\boxed{\;\kappa a_0-(v^{[1]}+v)-\tfrac1\kappa\bigl(B^2-\gamma B\bigr)=A\in\mathbb C.\;}$$
  Hence $A_0=\deg a_0=\max(V,2P)$ with $\operatorname{lc}(a_0)$ the leading
  coefficient of the RHS (see below). This is the quantum replacement for the
  classical $a_0=\tfrac{p^2}4+\tfrac{2hv}\kappa-A$.
- **Central telescoping** (from $Q_0=1$): $Q_0=(T-1)G$ with
  $$G=w^{[1]}+w+p^{[-1]}v-\kappa u,\qquad (T-1)G=1\ \Rightarrow\ G=E+c_0.$$
  Point conditions $s(0)=s(1)=v(0)=0$ (so $w(0)=w(1)=0$, and $u(0)=0$ is then
  automatic) give $G(0)=0$, hence $c_0=0$:
  $$\boxed{\;w^{[1]}+w+p^{[-1]}v-\kappa u=E.\;}\tag{$\star$}$$
  Equivalently $\kappa u=\mu(s^{[1]}+s)+p^{[-1]}v-E$; this **determines** $u$,
  and $E\mid u$ is automatic.

So the free data are $p,v,s$ and the scalars $\kappa,\gamma,\mu,A$; $B,a_0,u,w$
are determined; and the only remaining equations are $Q_{-1}=Q_{-2}=Q_{-3}=0$
($Q_{-4}=0$ already gave $w=\mu s$).

## 3. Structural leading-coefficient lemmas (`verify` §2)

Set $\phi:=\mu u-v$ and $\Phi:=\deg\phi$.

- **Lemma R (quantum rigidity).** $Q_{-3}=s^{[-1]}\phi-s\phi^{[-2]}$. If
  $\phi\ne0$ the coefficient of $E^{S+\Phi-1}$ equals $s_S\,\phi_\Phi\,(2\Phi-S)$;
  vanishing forces
  $$\boxed{\,S=2\Phi\,}\qquad(\phi\ne0).$$
  (Leading coefficients: $s_S=\operatorname{lc}s\ne0$, $\phi_\Phi=\operatorname{lc}\phi\ne0$.)
- **Lemma B.** $Q_{-2}=s\,\Psi+(v^{[-1]}u-u^{[-1]}v)$ with
  $\Psi=\mu(a_0-a_0^{[-2]})+(B^{[-2]}-B)$; the bracket has degree $\le U+V-1$
  and $E^{U+V-1}$-coefficient $(U-V)\,u_U v_V$ (so it drops when $U=V$).
- **Lemma C.** $Q_{-1}=\mu\bigl(s^{[1]}p-p^{[-2]}s\bigr)+\kappa(s-s^{[1]})
  +v(a_0-a_0^{[-1]})+u(B^{[-1]}-B)$; the first bracket has degree $S+P-1$ with
  $E^{S+P-1}$-coefficient $(S+2P)\,s_S p_P$ (nonzero when $P\ge1$).
- **Lemma D.** $\deg(a_0-a_0^{[-2]})=A_0-1$ with leading coefficient
  $2A_0\operatorname{lc}(a_0)$ (for $A_0\ge1$).

## 4. Emptiness proof

Assume a solution with $\kappa\mu\ne0$, $s\ne0$. We derive a contradiction in
every case. All leaves are machine-checked in `verify` §3.

### 4.1 Case $\phi=0$ (i.e. $v=\mu u$)

Then $v^{[-1]}u-u^{[-1]}v=0$ (Lemma B), so $Q_{-2}=s\Psi=0$; as $s\ne0$,
$\Psi=0$, i.e. $\mu(a_0-a_0^{[-2]})=B-B^{[-2]}$.

- If $P=0$: $B$ constant, so $\Psi=\mu(a_0-a_0^{[-2]})=0$, i.e. $a_0$ is
  $2$-periodic, hence constant; by the ladder-1 integral (with $B$ constant)
  $v^{[1]}+v$ is constant, so $v$ is constant, so $v=0$ (as $E\mid v$), whence
  $u=v/\mu=0$. Then $(\star)$ reads $\mu(s^{[1]}+s)=E$, impossible since
  $\deg(s^{[1]}+s)=S\ge2\ne1$.
- If $P\ge1$: $B-B^{[-2]}$ has degree $P-1$, so $\Psi=0$ forces
  $\deg(a_0-a_0^{[-2]})\le P-1$, hence $A_0\le P$ (Lemma D), hence
  $\deg\Psi\le P-1$. Polynomiality of $u=(\mu(s^{[1]}+s)-E)/(\kappa-\mu p^{[-1]})$
  forces $\deg u=S-P$, so $V=\deg v=\deg u=S-P$ and $S>V$. In Lemma C the terms
  $\kappa(s-s^{[1]})$, $v(a_0-a_0^{[-1]})$, $u(B^{[-1]}-B)$ then have degrees
  $\le S-1$, $\le V+(A_0-1)\le S-1$, $\le(S-P)+(P-1)=S-1$ — all $<S+P-1$ — while
  $\mu(s^{[1]}p-p^{[-2]}s)$ has degree $S+P-1$ and coefficient
  $(S+2P)\mu s_S p_P\ne0$. Hence $Q_{-1}\ne0$, contradiction.

### 4.2 Case $\phi\ne0$: $S=2\Phi$ (Lemma R)

Central term-degrees in $\kappa u=\mu(s^{[1]}+s)+p^{[-1]}v-E$ are $S,\ P+V,\ 1$.

**Sub-case $P=0$** ($p=p_0$ const, $B$ const, $A_0=V$):
- $V<S$: top of $\kappa u$ is the $s$-term, $U=S$, $\deg(\mu u)=S>V$, so
  $\Phi=S$ and $S=2S$ — contradiction.
- $V=S$: $U=S$; in $Q_{-2}$ the bracket coefficient at $E^{2S-1}$ vanishes
  ($U=V$, Lemma B) while $s\Psi$ contributes $4\mu S\,s_S v_S/\kappa\ne0$
  (Lemmas D, ladder-1). So $Q_{-2}\ne0$ — contradiction.
- $V>S$: $U=V$; $\phi_V=v_V(\mu p_0-\kappa)/\kappa$. If $\mu p_0\ne\kappa$ then
  $\Phi=V$, $S=2V>2S$ — contradiction. If $\mu p_0=\kappa$ then
  $\phi=\tfrac\mu\kappa(\mu(s^{[1]}+s)-E)$ has degree $S$, so $\Phi=S$ and
  $S=2S$ — contradiction.

**Sub-case $P\ge1$** (split on $P+V$ vs $S$, using $\Phi\le\max(U,V)$):
- $P+V>S$: $U=P+V$ (unique top), $\Phi=P+V>V$, so $S=2(P+V)>S$ — contradiction.
- $P+V<S$: $U=S$, $\deg(\mu u)=S>V$, so $\Phi=S$ and $S=2S$ — contradiction.
- $P+V=S$:
  - if $2\mu s_S+p_Pv_V\ne0$: $U=S$, $\Phi=S$ ($>V$ as $P\ge1$), $S=2S$ —
    contradiction.
  - if $2\mu s_S+p_Pv_V=0$ (so $U\le S-1$): split on $2P$ vs $V\,(=S-P)$:
    - $2P\ne V$: then $\operatorname{lc}(a_0)\ne0$ ($=p_P^2/4$ if $2P>V$, or
      $2v_V/\kappa$ if $2P<V$), and the term $T_1=s\mu(a_0-a_0^{[-2]})$ of
      $Q_{-2}$ has degree $S+A_0-1$ strictly above $\deg T_2=S+P-1$ and
      $\deg T_3\le U+V-1$, with coefficient $2A_0\mu s_S\operatorname{lc}(a_0)\ne0$.
      So $Q_{-2}\ne0$ — contradiction.
    - $2P=V$ (hence $S=3P$): the **thin locus $L$**. If
      $\operatorname{lc}(a_0)=p_P^2/4+2v_V/\kappa\ne0$, $T_1$ still dominates
      $Q_{-2}$ (degree $5P-1$ vs $\le5P-2$) — contradiction. If
      $\operatorname{lc}(a_0)=0$, then $A_0\le2P-1$ and $U\le S-1=3P-1$, so in
      Lemma C the three tail terms have degree $\le4P-2$ while
      $\mu(s^{[1]}p-p^{[-2]}s)$ has degree $4P-1$ and coefficient
      $(S+2P)\mu s_S p_P=5P\mu s_S p_P\ne0$. So $Q_{-1}\ne0$ — contradiction.
      (`verify` §3 confirms this at $P=2$: the $E^{7}$-coefficient of $Q_{-1}$
      on $L$ equals $5\kappa p_P^4/8\ne0$.)

Every case is contradictory, so the branch is empty. $\qquad\blacksquare$

## 5. The mechanism, contrasted with the classical route

The two closers are: (i) **Lemma R** ($S=2\Phi$), a rigidity produced by the
staggered integer shift $\phi^{[-2]}$ in $Q_{-3}$ — it collapses all the
"generic" cases immediately; and (ii) the **leading term of $Q_{-1}$**,
$\mu(s^{[1]}p-p^{[-2]}s)$ with coefficient $(S+2P)s_Sp_P$, which kills the two
degenerate loci ($\phi=0,\,P\ge1$ and the thin $L$) where $a_0$'s top cancels.

This is structurally different from the classical A\* closure. The classical
$C_{-3}=2a_{-1}'b_{-2}-a_{-1}b_{-2}'+\dots$ has "leading" content a Wronskian
$a_{-1}b_{-2}'-a_{-1}'b_{-2}$ that carries **no** analogue of $2\Phi-S$:
differentiation is shift-invariant, so classical $C_{-3}$ does **not** force
$\deg a_{-2}=2\deg\phi$. The classical branch is instead closed through explicit
first integrals of the residual system (see `classical-Astar.md`). The quantum
closure is comparatively short precisely because Lemma R does the bulk of the
work — a rigidity available only in the difference (quantum) setting.

## 6. Machine verification (`verify_quantum_mirror.py`)

- §0 — $Q_m$ reproduce audited tail equations (general $h$); gauge $Q_4=Q_3=0$.
- §1 — midpoint rewrite, ladder-1 integral, central telescoping, $c_0=0$,
  automatic $E\mid u$.
- §2 — Lemmas R, B, C, D as exact leading-coefficient identities (several
  degrees each).
- §3 — every case leaf of §4: the $\phi=0$ collapse, $P=0$ cases (a)/(b),
  $Q_{-2}$ $T_1$-dominance, and the locus-$L$ kill (clean value $5\kappa p_P^4/8$).
- §4 — bounded-degree **exhaustive emptiness** by Gröbner + Rabinowitsch
  saturation of $\kappa\mu\ne0$ and $s\ne0$, over a grid
  ($\deg p\le2,\ \deg v\le3,\ \deg s\le4$ and more); all EMPTY. A separate
  extended sweep (not in the fast verifier) confirms EMPTY up to
  $\deg p\le3,\ \deg v\le4,\ \deg s\le6$, including the minimal instance of
  locus $L$ ($S=6,V=4,P=2$).

## 7. Aside: the general-$h$ central telescoping (recorded separately)

Attack line 3 asked whether the classical integrated identity
$\bigl(2h^2w+hpv-\kappa h\,a_{-1}\bigr)'=1$ transfers to the quantum sector. It
does, qualitatively: for **general** $h$ (still gauge $b_2=0,\ b_1=\kappa h$,
$a_2=hh^{[1]}$, $a_1=hp$), $Q_0$ is again a perfect unit difference of a
**local** potential,
$$Q_0=(T-1)G,\qquad G\ \text{a $\mathbb C$-combination of shifted products of }h,w,p,v,u,$$
found and checked by an exact linear solve (residual $0$; script kept in scratch,
not a deliverable). At $h=1$ it collapses to the identity of §2,
$G=w^{[1]}+w+p^{[-1]}v-\kappa u$. Unlike the classical $2h^2w+hpv-\kappa hu$, the
general-$h$ potential does **not** have a one-line closed form — the staggered
shifts $h^{[i]}h^{[j]}$ mix, exactly as the audited memo warned. Extracting a
canonical $G$ and reading off the nonconstant-$h$ divisibility (the analogue of
"$h\mid t+e$", i.e. "nonconstant $h$ is affine") is the natural next step; it is
**not** needed for, and lies beyond, the constant-$h$ closure proved here.

## 8. Status of claims

- **Proved, arbitrary degree:** the three integrals of §2; Lemmas R, B, C, D;
  the emptiness theorem of §4 (structural, not degree-bounded).
- **Corroboration only:** the bounded Gröbner sweeps of §6.4.
- **Out of scope:** $\kappa=0$; the shifted-square sector with *nonconstant* $h$
  (the audited "nonconstant $h$ affine" lemma remains open — Lemma R and the
  $Q_{-1}$ mechanism here are stated for $h=1$); non-square sectors; any global
  DC1/JC2 statement. This memo closes exactly one branch: the resistant
  constant-$h$ one.
- **No counterexample exists in this branch** (it is empty), so no
  $A_1$-generation test arises here.
