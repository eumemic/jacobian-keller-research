# Independent attack on the classical square sector of two-sided band 2

**INDEPENDENTLY DERIVED AND AUDITED — NOT PEER REVIEWED — PARTIAL CASCADE ONLY**

This memo does not establish a full classical band-2 classification, JC2, or DC1. The resistant arbitrary-degree branch identified below remains unresolved, and no counterexample was found.

## Executive verdict

Work over \(\mathbb C\), put \(t=\tau=x\xi\), and use
\[
\{G,F\}=G_\xi F_x-G_xF_\xi=1,
\qquad
F=\sum_{j=-2}^2x^ja_j(t),\quad G=\sum_{j=-2}^2x^jb_j(t).
\]
Assume \(a_2=h^2\ne0\). After subtracting the forced top multiple of \(F\) from \(G\), the complete top cascade begins
\[
b_2=0,\qquad b_1=\kappa h.
\]
This memo concentrates on the requested genuine homogeneous branch \(\kappa\ne0\), but records the lost branch \(\kappa=0\) separately.

The following arbitrary-degree statements are proved below:

1. For \(\kappa\ne0\), ladder 2 forces
   \[
   h\mid a_1,
   \quad a_1=hp,
   \quad b_0=\frac\kappa2p+\beta.
   \]
2. Ladder 1 then integrates exactly. With \(v=b_{-1}\),
   \[
   a_0=\frac{p^2}{4}+\frac{2hv}{\kappa}-A
   \]
   for a scalar \(A\).
3. The central equation integrates exactly. With \(w=b_{-2}\),
   \[
   a_{-1}=\frac{pv}{\kappa}+\frac{2hw}{\kappa}-\frac{t+e}{\kappa h}
   \]
   for a scalar \(e\). Hence polynomiality of the *coefficient* \(a_{-1}\) already implies
   \[
   \boxed{h\mid t+e}.
   \]
   Thus a nonconstant \(h\) is necessarily affine, with a single simple root. This claim does not require a degree bound or the negative equations.
4. The four remaining equations are displayed exactly below. They have a clean branch split at \(w=b_{-2}=0\) versus \(w\ne0\). In the latter branch, ladder \(-4\) gives
   \[
   a_{-2}=\widetilde\mu\,b_{-2}.
   \]
   The branch \(w\ne0,\widetilde\mu\ne0\) is a plausible match for the reported resistant “case A, \(\widetilde\mu\ne0\)” branch, subject to checking the author's parameter definitions (interchanging the ratio would replace \(\widetilde\mu\) by \(1/\widetilde\mu\)). I do **not** have an arbitrary-degree proof solving this residual branch.

There is a natural canonical family which solves the whole localized system:
\[
Y=xh(t),\qquad Z=x^{-1}\frac{t+e}{h(t)},\qquad \{Z,Y\}=1,
\]
\[
U=Y+c_0+c_1Z,
\qquad
F=U^2-\frac1\kappa Z-A,
\qquad
G=\lambda F+\kappa U+\beta.
\]
For constant \(h\), the genuine polynomial members necessarily have \(e=0\), and these are tame; for nonconstant \(h\), divisibility forces \(h=\alpha(t+e)\), so \(Z=1/(\alpha x)\), and every member fails polynomial membership. This proves the advertised conclusions **for this family**, not exhaustion of the resistant residual branch.

Accordingly, the defensible classification is:

- nonconstant \(h\Rightarrow h\mid t+e\): **proved**;
- the displayed nonconstant localized families are polar and fail membership: **proved**;
- every nonconstant solution belongs to those families: **not proved** because of the unsolved \(w\ne0,\widetilde\mu\ne0\) residual branch;
- constant \(h\) gives the displayed tame polynomial families: **proved as a construction**;
- constant \(h\) gives *exactly* those tame families: **not proved** for the same residual branch;
- therefore no complete classical square-sector theorem follows from this audit.

## 1. The nine equations

For one-term pieces,
\[
\{x^\ell b(t),x^ka(t)\}=x^{k+\ell}(ka b'-\ell a'b).
\]
Thus \(C_m=\sum_{k+\ell=m}(ka_kb_\ell'-\ell a_k'b_\ell)\), and \(C_m=\delta_{m0}\). Explicitly:
\[
\begin{aligned}
C_4={}&2a_2b_2'-2a_2'b_2,\\
C_3={}&2a_2b_1'-a_2'b_1+a_1b_2'-2a_1'b_2,\\
C_2={}&2a_2b_0'+a_1b_1'-a_1'b_1-2a_0'b_2,\\
C_1={}&2a_2b_{-1}'+a_2'b_{-1}+a_1b_0'-a_0'b_1-a_{-1}b_2'-2a_{-1}'b_2,\\
C_0={}&2a_2b_{-2}'+2a_2'b_{-2}+a_1b_{-1}'+a_1'b_{-1}
-a_{-1}b_1'-a_{-1}'b_1-2a_{-2}b_2'-2a_{-2}'b_2,\\
C_{-1}={}&a_1b_{-2}'+2a_1'b_{-2}+a_0'b_{-1}-a_{-1}b_0'
-2a_{-2}b_1'-a_{-2}'b_1,\\
C_{-2}={}&2a_0'b_{-2}+a_{-1}'b_{-1}-a_{-1}b_{-1}'-2a_{-2}b_0',\\
C_{-3}={}&2a_{-1}'b_{-2}-a_{-1}b_{-2}'+a_{-2}'b_{-1}-2a_{-2}b_{-1}',\\
C_{-4}={}&2a_{-2}'b_{-2}-2a_{-2}b_{-2}'.
\end{aligned}
\]
The extra \(a_2'b_{-1}\), etc., are essential; omitting them changes the integrated identities.

## 2. Top reduction, including zero branches

Since \(a_2=h^2\ne0\), \(C_4=0\) gives \(b_2=\lambda h^2\). Replace \(G\) by \(G-\lambda F\). Thereafter
\[
b_2=0.
\]
Equation \(C_3=0\) becomes
\[
2h^2b_1'-2hh'b_1=0,
\]
so \((b_1/h)'=0\) in \(\mathbb C(t)\), hence
\[
b_1=\kappa h.
\]
There are three distinct situations:

- \(h=0\): outside the assumed sector \(a_2\ne0\);
- \(\kappa=0\): a genuine branch, in which \(C_2\) gives only \(b_0=\beta\) and does not force \(h\mid a_1\);
- \(\kappa\ne0\): the retained homogeneous branch analyzed next.

The original top multiple \(\lambda\) can be restored at the end by replacing the reduced \(G\) by \(G+\lambda F\).

## 3. Positive cascade for \(\kappa\ne0\)

### Ladder 2

Equation \(C_2=0\) is
\[
2h^2b_0'+\kappa(a_1h'-a_1'h)=0.
\]
Set \(r=a_1/h\in\mathbb C(t)\). Then
\[
2b_0'=\kappa r'.
\]
Therefore \(r=2b_0/\kappa+c\) is polynomial. Absorb the scalar \(c\) into a polynomial \(p\), obtaining
\[
\boxed{a_1=hp,\qquad b_0=\frac\kappa2p+\beta.}\tag{3.1}
\]
This is an arbitrary-degree divisibility result, not a bounded ansatz.

### Ladder 1

Write \(v=b_{-1}\). Substitution of (3.1) into \(C_1=0\) gives
\[
2(hv)'+\frac\kappa2(pp'-2a_0')=0.
\]
After one integration,
\[
\boxed{a_0=\frac{p^2}{4}+\frac{2hv}{\kappa}-A,}\tag{3.2}
\]
where \(A\in\mathbb C\). The sign/name of \(A\) is conventional.

### Central equation

Write \(w=b_{-2}\). Substituting (3.1)--(3.2) in \(C_0=1\) makes all \(p'v,pv'\) terms telescope:
\[
\left(2h^2w+hpv-\kappa h a_{-1}\right)'=1.
\]
Hence
\[
2h^2w+hpv-\kappa h a_{-1}=t+e,
\]
or
\[
\boxed{a_{-1}=\frac{pv}{\kappa}+\frac{2hw}{\kappa}-\frac{t+e}{\kappa h}.}\tag{3.3}
\]
Every term except the last is polynomial. Because \(a_{-1}\in\mathbb C[t]\),
\[
\boxed{h\mid t+e.}\tag{3.4}
\]
No membership condition and no negative-tail equation has yet been used. If \(h\) is nonconstant, then
\[
h=\alpha(t+e),\qquad \alpha\in\mathbb C^*.
\]
If \(h\) is constant, (3.4) is automatic.

Polynomial membership is stronger than coefficient polynomiality:
\[
t\mid a_{-1},v,
\qquad
t^2\mid a_{-2},w.
\]
These constraints remain to be imposed on the residual system.

## 4. Exact residual negative cascade

Let
\[
s=a_{-2},
\]
and define \(a_0,a_{-1}\) by (3.2)--(3.3). Then the unsolved equations are exactly
\[
\begin{aligned}
R_{-1}:={}&a_1w'+2a_1'w+a_0'v-a_{-1}b_0'-2s(\kappa h)'-s'\kappa h=0,\\
R_{-2}:={}&2a_0'w+a_{-1}'v-a_{-1}v'-2s b_0'=0,\\
R_{-3}:={}&2a_{-1}'w-a_{-1}w'+s'v-2sv'=0,\\
R_{-4}:={}&2s'w-2sw'=0.
\end{aligned}\tag{4.1}
\]
This compact form is preferable to a long expanded expression and is exact after the proved integrations.

The branch tree is:

### B0: \(w=0\)

Then \(R_{-4}\) is automatic. The remaining equations are
\[
a_0'v-a_{-1}b_0'-2s(\kappa h)'-s'\kappa h=0,
\]
\[
a_{-1}'v-a_{-1}v'-2s b_0'=0,
\]
\[
s'v-2sv'=0.
\]
This includes \(s=0\), \(v=0\), and unit branches. One must not infer a ratio relation from \(R_{-4}\) here.

### A: \(w\ne0\)

Since \((s/w)'=0\),
\[
\boxed{s=\widetilde\mu w,\qquad \widetilde\mu\in\mathbb C.}\tag{4.2}
\]
This splits again:

- A0: \(\widetilde\mu=0\), so \(a_{-2}=0\) although \(b_{-2}\ne0\);
- A*: \(\widetilde\mu\ne0\), so both negative extremes are nonzero.

If another derivation defines \(b_{-2}=\mu a_{-2}\), then its \(\mu\) equals \(1/\widetilde\mu\) on A*. Therefore “case A with \(\widetilde\mu\ne0\)” almost certainly means the branch A* here: nonzero \(b_{-2}\), nonzero proportional \(a_{-2}\), after the top gauge and positive integrations. This is the branch that my hand derivation does not close.

The membership conditions in A* imply \(t^2\mid w\), hence automatically \(t^2\mid s\). They do not by themselves force \(w=0\).

## 5. Canonical localized family and exact symbolic check

Let
\[
q(t)=\frac{t+e}{h(t)},\qquad
Y=xh(t),\qquad Z=x^{-1}q(t).
\]
The one-term bracket rule gives
\[
\{Z,Y\}=(hq)'=1.
\]
For arbitrary \(c_0,c_1,A,\beta,\lambda\in\mathbb C\), set
\[
U=Y+c_0+c_1Z,
\]
\[
\boxed{F=U^2-\frac1\kappa Z-A,\qquad
G=\lambda F+\kappa U+\beta.}\tag{5.1}
\]
Then
\[
\{G,F\}=\kappa\{U,F\}
=-\{U,Z\}=1.
\]
Indeed, \(\{Z,Y\}=1\) gives \(\{U,Z\}=\{Y,Z\}=-1\), fixing the displayed sign.

Its coefficients in the reduced gauge \(\lambda=0\) are
\[
\begin{aligned}
a_2&=h^2,& a_1&=2c_0h,& a_0&=c_0^2+2c_1hq-A,\\
a_{-1}&=2c_0c_1q-q/\kappa,& a_{-2}&=c_1^2q^2,\\
b_1&=\kappa h,& b_0&=\kappa c_0+\beta,& b_{-1}&=\kappa c_1q,& b_{-2}&=0.
\end{aligned}
\]
Thus this lies in branch B0, and it realizes the first integrals above.

### Nonconstant \(h\)

Then (3.4) gives \(h=\alpha(t+e)\), so \(q=1/\alpha\) is a nonzero constant and
\[
Z=\frac1{\alpha x}.
\]
Membership fails in every subcase:

- if \(c_1\ne0\), then \(a_{-2}=c_1^2/\alpha^2\) is a nonzero constant, not divisible by \(t^2\);
- if \(c_1=0\), then \(a_{-1}=-1/(\alpha\kappa)\) is a nonzero constant, not divisible by \(t\).

No choice of \(c_0\) repairs the second assertion when \(c_1=0\), and when \(c_1\ne0\), failure at degree \(-2\) is already decisive. Hence every nonconstant member of (5.1) is a Laurent polar pair and not a genuine polynomial pair.

This is sharp as a localized construction, but is not an exhaustion proof for branch A*.

### Constant \(h=d\in\mathbb C^*\)

Now
\[
Y=dx,
\qquad
Z=x^{-1}\frac{t+e}{d}=\frac\xi d+\frac e{dx}.
\]
Genuine polynomial membership in (5.1) forces \(e=0\): if \(c_1\ne0\), the \(e^2x^{-2}\) term fails degree \(-2\) membership, while if \(c_1=0\), the term \(e/(d\kappa x)\) fails degree \(-1\) membership. Thus
\[
Z=\xi/d.
\]
Then \((Y,Z)=(dx,\xi/d)\) is a linear symplectic coordinate pair with \(\{Z,Y\}=1\), and
\[
U=dx+c_0+(c_1/d)\xi
\]
is affine linear. Formula (5.1) is the standard tame triangular construction: an affine symplectic change, followed by adding a quadratic polynomial of one coordinate and affine recombination. Direct symbolic expansion gives \(\{G,F\}=1\).

So the genuine polynomial constant-\(h\) members of (5.1), which necessarily have \(e=0\), are tame. What remains unproved is that every constant-\(h\) solution of (4.1), especially A*, reduces to (5.1) or another tame normal form.

## 6. Why the resistant branch is genuinely not dispatched here

When \(h=d\in\mathbb C^*\), the positive integrations become
\[
a_1=dp,
\quad b_0=\frac\kappa2p+\beta,
\quad a_0=\frac{p^2}{4}+\frac{2dv}{\kappa}-A,
\]
\[
a_{-1}=\frac{pv}{\kappa}+\frac{2dw}{\kappa}-\frac{t+e}{\kappa d}.
\]
In A*, put \(s=\widetilde\mu w\) with \(\widetilde\mu\ne0\). Equations \(R_{-1},R_{-2},R_{-3}\) remain a coupled nonlinear differential-polynomial system in \(p,v,w\), subject to
\[
t\mid v,a_{-1},
\qquad t^2\mid w.
\]
The top and central first integrals do not make this system triangular. In particular, neither \(w=0\) nor \(v=0\) follows from a single product invariant. A bounded-degree solver may find only the B0 tame family, but that does not exclude arbitrary-degree A* solutions.

Within this memo, the parameter correspondence is:

- \(w=b_{-2}\) is the branch-defining nonzero negative coefficient;
- \(\widetilde\mu=a_{-2}/b_{-2}\);
- if the reported derivation instead first used \(b_{-2}=\mu a_{-2}\), then \(\mu=1/\widetilde\mu\);
- “\(\widetilde\mu\ne0\)” means both negative extremes survive;
- the branch occurs only after \(\kappa\ne0\) has retained the positive homogeneous mode.

This makes A* a plausible match for the reported branch, but that identification requires checking the author's formulas and parameter definitions.

## 7. The omitted \(\kappa=0\) branch

For completeness, after the same top subtraction, \(\kappa=0\) gives
\[
b_2=b_1=0,
\qquad b_0=\beta.
\]
Then
\[
C_1=2h(hv)'=0,
\]
so \(hv=\nu\in\mathbb C\). If \(h\) is nonconstant and \(v\) polynomial, then \(\nu=0\), hence \(v=0\). The central equation reduces to
\[
2(h^2w)'=1,
\qquad h^2w=\frac t2+e_0.
\]
Thus \(h^2\mid t+2e_0\). A nonconstant \(h\) is impossible because a nonconstant square cannot divide a nonzero affine polynomial. Therefore \(h\) is constant in this branch. This is a useful complete exclusion of nonconstant \(h\) when \(\kappa=0\).

For constant \(h=d\), central normalization makes \(w\) affine with nonzero slope, which cannot satisfy the genuine membership condition \(t^2\mid w\). Hence the entire \(\kappa=0\) branch contains no genuine polynomial pair with \(a_2\ne0\). It may contain localized pairs. Thus every genuine square-sector pair must indeed lie in \(\kappa\ne0\), but that conclusion comes from this separate argument, not from dividing by \(\kappa\) at the outset.

## 8. Proof versus computation

### Proved algebraically, without degree bounds

- all nine coefficient equations;
- top proportionality and the full parameter \(b_1=\kappa h\);
- complete treatment of \(\kappa=0\) for genuine polynomial membership;
- for \(\kappa\ne0\), \(h\mid a_1\), formulas (3.1)--(3.3), and \(h\mid t+e\);
- exact residual system (4.1) and its branch split;
- exact canonical family (5.1);
- membership exclusion of its nonconstant polar members;
- tameness of its constant polynomial members.

### Symbolically checked identities

Using exact SymPy differentiation and simplification, I checked:

- the simplifications leading to the residual equations;
- \(\{G,F\}=1\) identically for (5.1), both with nonconstant affine \(h\) and with constant \(h\);
- the explicit Laurent coefficients and their membership failures.

These checks support identities only and are not completeness evidence.

### Not proved

- exhaustion of branches B0, A0, or A* by (5.1);
- nonexistence of A* under the membership divisibilities;
- exact tame classification for constant \(h\);
- consequently, full band-2 classical rigidity.

No counterexample to polynomial rigidity was found. The localized nonconstant family (5.1) is a counterexample to any claim that the coefficient ODEs alone kill nonconstant \(h\), but not to the genuine polynomial statement. The principal mathematical gap is A*, equivalently the likely reported case A with nonzero \(\widetilde\mu\).

## Final assessment

The reported split is half theorem and half conjectural classification. The strongest clean theorem recovered independently is:

> In the square sector with \(a_2=h^2\ne0\), every genuine polynomial pair has \(\kappa\ne0\), and then \(h\mid t+e\). Thus nonconstant \(h\) is affine. A canonical complete-looking localized family is polar and excluded by membership when \(h\) is nonconstant, while its constant-\(h\) specialization is tame.

What is missing is an arbitrary-degree solution or exclusion of the residual branch \(b_{-2}\ne0\), \(a_{-2}=\widetilde\mu b_{-2}\), \(\widetilde\mu\ne0\). Until that branch is closed, “nonconstant gives only polar families” and “constant gives exactly tame families” should not be stated as proved classifications.
