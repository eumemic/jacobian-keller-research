# Complete classification of the oriented classical a₂-square sector under band-2 support containment

**INDEPENDENTLY DERIVED AND AUDITED — NOT PEER REVIEWED — COMPLETE ONLY IN THE STATED SECTOR**

## Executive verdict and exact scope

Work over \(\mathbb C\), put \(t=x\xi\), and use the Poisson convention
\[
\{G,F\}=G_\xi F_x-G_xF_\xi=1.
\]
Suppose
\[
F=\sum_{j=-2}^2x^ja_j(t),\qquad G=\sum_{j=-2}^2x^jb_j(t),
\]
with the memo's genuine polynomial membership conditions
\[
t\mid a_{-1},b_{-1},\qquad t^2\mid a_{-2},b_{-2},
\]
and assume that the ladder-\(2\) coefficient of \(F\) is a nonzero square,
\[
a_2=h^2\ne0.
\]
After subtracting a scalar multiple of \(F\) from \(G\) and applying the diagonal symplectic scaling that normalizes the resulting constant \(h\) to \(1\), every such pair has the form
\[
U=x+c_0+c_1\xi,
\]
\[
\boxed{F=U^2-\frac1\kappa\xi-A,\qquad
G=\lambda F+\kappa U+\beta,}
\tag{0.1}
\]
where \(\kappa\in\mathbb C^*\) and \(c_0,c_1,A,\beta,\lambda\in\mathbb C\). Conversely, every parameter choice in (0.1) gives a normalized pair whose ladder supports are contained in \([-2,2]\), which satisfies genuine polynomial membership and obeys \(\{G,F\}=1\); undoing the scaling gives every original-coordinate pair.

This is a complete arbitrary-degree classification of the **oriented classical \(a_2\)-square sector \(a_2=h^2\ne0\)** under the support-containment and genuine membership hypotheses above. It is not a full band-2 theorem: it does not establish orientation or case coverage beyond this sector. It does not resolve the quantum shifted-square sector, which remains partial, and it proves neither JC2 nor DC1.

## 1. Coefficient equations and top cascade

For one-term pieces,
\[
\{x^\ell b(t),x^ka(t)\}=x^{k+\ell}(ka b'-\ell a'b).
\]
Thus the coefficient of \(x^m\) in \(\{G,F\}\) is
\[
C_m=\sum_{k+\ell=m}(ka_kb_\ell'-\ell a_k'b_\ell),
\qquad C_m=\delta_{m0}.
\]
Explicitly,
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
\tag{1.1}
\]

Since \(a_2=h^2\), equation \(C_4=0\) gives \(b_2=\lambda h^2\). Replace \(G\) by \(G-\lambda F\). Then \(b_2=0\), and \(C_3=0\) gives
\[
(b_1/h)'=0,
\qquad
\boxed{b_2=0,\quad b_1=\kappa h.}
\tag{1.2}
\]
The parameter \(\lambda\) is restored at the end.

### The branch \(\kappa=0\) has no genuine members

If \(\kappa=0\), then \(b_1=0\), and \(C_2=0\) gives \(b_0=\beta\). Writing \(v=b_{-1}\) and \(w=b_{-2}\), equation \(C_1=0\) says \((hv)'=0\). If \(h\) is nonconstant, polynomiality forces \(v=0\); if \(h\) is constant, membership \(t\mid v\) forces the same conclusion. The central equation is then
\[
2(h^2w)'=1,
\qquad h^2w=t/2+e_0.
\]
A nonconstant square cannot divide a nonzero affine polynomial, while for constant \(h\) the resulting affine \(w\) cannot satisfy \(t^2\mid w\). Thus genuine membership excludes \(\kappa=0\) entirely. Henceforth \(\kappa\ne0\).

## 2. Positive and central cascade for \(\kappa\ne0\)

Equation \(C_2=0\) gives, without a degree bound,
\[
\boxed{a_1=hp,\qquad b_0=\frac\kappa2p+\beta}
\tag{2.1}
\]
for a polynomial \(p\). With \(v=b_{-1}\), equation \(C_1=0\) integrates to
\[
\boxed{a_0=\frac{p^2}{4}+\frac{2hv}{\kappa}-A.}
\tag{2.2}
\]
With \(w=b_{-2}\), substitution in \(C_0=1\) gives the central identity
\[
\boxed{(2h^2w+hpv-\kappa h a_{-1})'=1,}
\tag{2.3}
\]
so
\[
2h^2w+hpv-\kappa h a_{-1}=t+e.
\tag{2.4}
\]
Membership makes every term on the left vanish at \(t=0\), and therefore \(e=0\). Coefficient polynomiality in (2.4) also gives \(h\mid t\). If \(h\) is nonconstant, then \(h=\alpha t\). But membership now puts every term on the left of (2.4) in \((t^2)\): use \(t^2\mid w\), \(t\mid v\), and \(t\mid a_{-1}\). This contradicts equality to \(t\). Therefore \(h\) is constant.

Write \(h=d\in\mathbb C^*\). Under the pullback substitution defined by the linear symplectic change \(x\mapsto x/d\), \(\xi\mapsto d\xi\), the quantity \(t=x\xi\), the bracket convention, support containment, and membership are preserved, while \(h\) is normalized to \(1\). We use that normalization from now on. Equations (2.1)--(2.4) become
\[
a_1=p,\quad b_0=\frac\kappa2p+\beta,\quad
a_0=\frac{p^2}{4}+\frac{2v}{\kappa}-A,
\tag{2.5}
\]
\[
a_{-1}=\frac{pv+2w-t}{\kappa}.
\tag{2.6}
\]
Put \(s=a_{-2}\).

## 3. Exclusion of the branch \(w=b_{-2}\ne0\)

Equation \(C_{-4}=0\) gives
\[
(s/w)'=0,
\qquad s=\mu w
\tag{3.1}
\]
for an arbitrary \(\mu\in\mathbb C\), including \(\mu=0\). Set
\[
P=p-\kappa\mu,
\qquad
Q=t-Pv-2w.
\tag{3.2}
\]
Direct substitution of (2.5)--(2.6) and (3.1) into \(C_{-1}=C_{-2}=C_{-3}=0\), with no division by \(P,v,Q\), gives
\[
P w' + (t/2+w)P' + (2/\kappa)vv'=0,
\tag{3.3}
\]
\[
(\kappa P w+v^2)P' + (t+2w)v' +2v w' -v=0,
\tag{3.4}
\]
\[
Qw'-2wQ'=0.
\tag{3.5}
\]
These equations retain all zero subbranches and are the independently certified arbitrary-degree reduction replacing the supplied checker's bounded ansatz.

### 3.1 The subcase \(Q=0\)

Here \(t=Pv+2w\). Substitution in (3.4) factors it as
\[
w(\kappa PP'+4v')=0.
\]
Because \(w\ne0\),
\[
v'=-\frac\kappa4PP',
\qquad
v=\gamma-\frac\kappa8P^2
\tag{3.6}
\]
for a scalar \(\gamma\). Equation (3.3) reduces to
\[
\kappa P^3P'-12PP'v+4P+8tP'=0.
\tag{3.7}
\]
If \(P\) is nonconstant of degree \(m\), then (3.6) has degree \(2m\) with leading coefficient \(-\kappa A^2/8\), where \(A\) is the leading coefficient of \(P\). The two degree-\(4m-1\) terms in (3.7) have combined leading coefficient
\[
\kappa mA^4+\frac{12\kappa mA^4}{8}
=\frac{5\kappa mA^4}{2}\ne0,
\]
a contradiction. Hence \(P\) is constant. Equation (3.7) then gives \(P=0\), and \(Q=0\) gives \(w=t/2\), contradicting \(t^2\mid w\). Thus \(Q=0\) is impossible.

### 3.2 The subcase \(Q\ne0\)

Since \(w\ne0\), equation (3.5) implies
\[
w=cQ^2,
\qquad c\in\mathbb C^*.
\tag{3.8}
\]
Indeed \((w/Q^2)'=0\) in \(\mathbb C(t)\). Membership \(t^2\mid w\) and unique factorization imply \(t\mid Q\).

Let
\[
q=\deg Q\ge1,
\quad m=\deg P,
\quad n=\deg v,
\]
and let \(C,A,B\) be the respective nonzero leading coefficients whenever the polynomial is nonzero. The defining identity for \(Q\) is
\[
t=Q+Pv+2cQ^2.
\tag{3.9}
\]
If \(m+n>2q\), the leading term of \(Pv\) in (3.9) cannot cancel; if \(m+n<2q\), the leading term of \(2cQ^2\) cannot cancel. Consequently \(P,v\ne0\),
\[
m+n=2q,
\tag{3.10}
\]
and cancellation at that degree gives
\[
AB=-2cC^2.
\tag{3.11}
\]
Membership excludes \(n=0\). If \(m=0\), then \(n=2q\), and the term \((2/\kappa)vv'\) in (3.3) has degree \(4q-1\), strictly above the other terms, which is impossible. Thus \(m,n>0\).

In (3.3), the combined leading coefficient of \(Pw'+wP'\) is \(cAC^2(2q+m)\ne0\), while \((2/\kappa)vv'\) has nonzero leading coefficient \((2n/\kappa)B^2\). Cancellation therefore requires
\[
m+2q=2n.
\tag{3.12}
\]
Together with (3.10), this gives
\[
m=\frac{2q}{3},\qquad n=\frac{4q}{3}.
\tag{3.13}
\]
Taking the leading coefficients in (3.3) and (3.4), respectively, yields
\[
cAC^2(2q+m)+\frac{2n}{\kappa}B^2=0,
\tag{3.14}
\]
\[
mA(\kappa cAC^2+B^2)+2cBC^2(n+2q)=0.
\tag{3.15}
\]
Use (3.11) and (3.13). Equation (3.14) becomes
\[
\kappa A^3+4cC^2=0,
\tag{3.16}
\]
whereas (3.15) becomes
\[
-\kappa A^3+16cC^2=0.
\tag{3.17}
\]
Adding (3.16) and (3.17) gives \(20cC^2=0\), contrary to \(cC\ne0\). Hence \(Q\ne0\) is also impossible.

The contradiction in both subcases proves
\[
\boxed{b_{-2}=w=0.}
\tag{3.18}
\]
This proof covers every \(\mu\), including the formerly resistant \(\mu\ne0\) branch and the zero branch \(\mu=0\), without a degree bound or computational exhaustion.

## 4. Complete classification when \(w=0\)

With \(w=0\), equations \(C_{-1}=C_{-2}=C_{-3}=0\) reduce exactly to
\[
4vv'+\kappa t p'-2\kappa^2s'=0,
\tag{4.1}
\]
\[
(v^2-\kappa^2s)p'+tv'-v=0,
\tag{4.2}
\]
\[
vs'-2sv'=0.
\tag{4.3}
\]
We classify all zero and nonzero branches.

### 4.1 The branch \(v=0\)

Equations (4.1)--(4.2) give
\[
\kappa tp'-2\kappa^2s'=0,
\qquad
-\kappa^2sp'=0.
\]
If \(s\ne0\), then \(p'=0\), hence \(s'=0\); membership \(t^2\mid s\) then forces \(s=0\), a contradiction. Thus \(s=0\), and (4.1) gives \(p'=0\). This is the final family with \(c_1=0\).

### 4.2 The branch \(s=0\) and \(v\ne0\)

Write \(v=tu\) by membership. Equations (4.1)--(4.2) become
\[
4vv'+\kappa tp'=0,
\qquad
u'+u^2p'=0.
\tag{4.4}
\]
If \(u\) is nonconstant, the second equation is impossible by degrees: its left summands have degrees \(\deg u-1\) and \(2\deg u+\deg p-1\), and \(p'=0\) would force \(u'=0\). If \(u\) is a nonzero constant, the second equation gives \(p'=0\), while the first gives \(4u^2=0\). Thus this branch has no solutions.

### 4.3 The branch \(v,s\ne0\)

Equation (4.3) gives
\[
s=Cv^2,
\qquad C\in\mathbb C^*.
\tag{4.5}
\]
Put
\[
D=1-\kappa^2C.
\]
Then (4.1)--(4.2) are
\[
4Dvv'+\kappa tp'=0,
\qquad
Dv^2p'+tv'-v=0.
\tag{4.6}
\]
Write \(v=tu\). The second equation gives
\[
u'=-Du^2p'.
\tag{4.7}
\]
If \(D\ne0\) and \(u\) is nonconstant, the two sides of (4.7) have degrees \(\deg u-1\) and \(2\deg u+\deg p-1\), an impossibility; if \(p'=0\), (4.7) already says \(u'=0\). If \(u\) is a nonzero constant, (4.7) gives \(p'=0\), and the first equation in (4.6), divided by \(t\), gives \(4Du^2=0\), again impossible. Therefore \(D\ne0\) has no solutions.

If \(D=0\), then \(C=1/\kappa^2\). Equations (4.6) give \(p'=0\) and \(tv'-v=0\). Therefore, for scalars \(c_0,c_1\),
\[
\boxed{p=2c_0,\qquad v=\kappa c_1t,\qquad s=c_1^2t^2.}
\tag{4.8}
\]
The case \(c_1=0\) agrees with the separately treated zero branch. Thus (4.8) is the complete solution of (4.1)--(4.3).

## 5. Reconstruction and bracket verification

Substituting (4.8) into (2.5)--(2.6) gives
\[
\begin{aligned}
a_2&=1,&a_1&=2c_0,&a_0&=c_0^2+2c_1t-A,\\
a_{-1}&=2c_0c_1t-t/\kappa,&a_{-2}&=c_1^2t^2,\\
b_1&=\kappa,&b_0&=\kappa c_0+\beta,&b_{-1}&=\kappa c_1t,&b_{-2}&=0.
\end{aligned}
\]
These are exactly the coefficients of (0.1) with
\[
U=x+c_0+c_1\xi.
\]
For the stated convention, \(\{U,\xi\}=-1\). Hence
\[
\{G,F\}
=\kappa\{U,U^2-\xi/\kappa-A\}
=-\{U,\xi\}=1.
\]
This also verifies the sign in (0.1). The displayed coefficients show genuine membership and common \([-2,2]\) support directly.

## 6. Claim disposition

### Proved algebraically at arbitrary degree

- the top cascade \(b_2=0\), \(b_1=\kappa h\) after subtracting \(\lambda F\);
- complete exclusion of \(\kappa=0\) under genuine membership;
- the positive integrations and central identity (2.3);
- \(e=0\), exclusion of nonconstant \(h\), and normalization to \(h=1\);
- the exact \(w\ne0\) system (3.3)--(3.5), including arbitrary \(\mu\);
- complete exclusion of both \(Q=0\) and \(Q\ne0\) by algebra and leading terms;
- complete treatment of every zero/nonzero branch of (4.1)--(4.3);
- reconstruction of (0.1) and direct verification of bracket, support, and membership;
- consequently, complete classification of the oriented classical \(a_2\)-square sector under this memo's exact hypotheses.

### Not claimed

- a full classical band-2 theorem or coverage of other orientations/sectors;
- a quantum classification;
- a proof of JC2 or DC1;
- peer review or completeness inferred from computation.

The earlier supplied checker remains unsuitable as a classification certificate. The proof here is independent of its bounded-degree searches and closes the classical resistant branch by an arbitrary-degree argument. The quantum caveats in this package are unchanged.
