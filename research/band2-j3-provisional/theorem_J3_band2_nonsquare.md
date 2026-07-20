# Audited J3 statement: the two-sided non-square sector of band 2

> **INDEPENDENTLY AUDITED — NOT PEER REVIEWED.** The independent audit found the core algebraic argument correct after the statement and proof repairs incorporated below. This is a narrowly scoped obstruction theorem, not a classification of band-2 pairs and not a proof of JC2 or DC1.

## Conventions and statement

Classically, work over \(\mathbb C\) in the Laurent Poisson algebra
\[
L=\mathbb C[x^{\pm1},\xi]=\bigoplus_{k\in\mathbb Z}x^k\mathbb C[\tau],
\qquad \tau=x\xi,
\]
with \(\{f,g\}=f_\xi g_x-f_xg_\xi\), so \(\{\xi,x\}=1\). Quantumly, work in the Ore localization
\[
A_1[x^{-1}]=\bigoplus_{k\in\mathbb Z}x^k\mathbb C[E],
\qquad [\partial,x]=1,\quad E=x\partial,
\]
using \([D,X]=DX-XD\). Thus \(f(E)x^j=x^jf(E+j)\) and
\[
(x^if(E))(x^jg(E))=x^{i+j}f(E+j)g(E).
\]
All displayed direct sums have finite ladder support. The **ladder support** of \(\sum x^ka_k\) is \(\{k:a_k\ne0\}\). A **common band-2 presentation** means that the supports of both \(X\) and \(D\) are contained in \([-2,2]\), precisely with no other ladder exponents. “Two-sided” below means \(a_2a_{-2}\ne0\), not merely that some positive and negative terms occur.

A nonzero classical polynomial is in the **square class** if it is \(c h(\tau)^2\); a nonzero quantum polynomial is in the **shifted-square class** if it is \(c h(E)h(E+1)\), where \(c\in\mathbb C^*\) and \(h\) is a polynomial. Nonzero constants belong to both classes.

**Theorem J3 (classical Laurent form).** Let
\[
X=\sum_{k=-2}^{2}x^ka_k(\tau),\qquad D=\sum_{l=-2}^{2}x^lb_l(\tau)
\]
with \(a_k,b_l\in\mathbb C[\tau]\), \(a_2a_{-2}\ne0\), and \(\{D,X\}=c_0\in\mathbb C^*\). Then \(a_2\) and \(a_{-2}\) cannot both lie outside the square class.

**Theorem J3q (quantum localized form).** Let
\[
X=\sum_{k=-2}^{2}x^ka_k(E),\qquad D=\sum_{l=-2}^{2}x^lb_l(E)
\]
with \(a_k,b_l\in\mathbb C[E]\), \(a_2a_{-2}\ne0\), and \([D,X]=1\). Then \(a_2\) and \(a_{-2}\) cannot both lie outside the shifted-square class.

The assertions concern only this two-sided common band presentation. They neither address one-sided sectors nor say that each nonzero extreme coefficient is a (shifted) square.

## Three elementary lemmas

**Rational periodicity lemma.** If \(r\in\mathbb C(t)\) and \(r(t+s)=r(t)\) for some \(s\in\mathbb C^*\), then \(r\) is constant.

**Proof.** A pole of \(r\) would give poles at all of its distinct translates by \(ns\), contrary to rationality. Hence \(r\) is a polynomial. If it has positive degree \(d\), then \(r(t+s)-r(t)\) has degree \(d-1\) with nonzero leading coefficient, a contradiction. ∎

**Classical homogeneous lemma.** For nonzero \(a\in\mathbb C[t]\), polynomial solutions of
\[
2ah'-a'h=0
\]
are either \(h=0\), or \(a=cq^2\) for some \(c\in\mathbb C^*\) and polynomial \(q\) (indeed one may take \(q=h\), after changing the scalar).

**Proof.** The zero solution is retained. If \(h\ne0\), direct differentiation gives \((h^2/a)'=0\) in \(\mathbb C(t)\). Thus \(h^2/a\in\mathbb C^*\), proving the assertion. ∎

**Quantum shifted-square lemma.** For nonzero \(a\in\mathbb C[E]\), polynomial solutions of
\[
h(E+2)a(E)=a(E+1)h(E)
\]
are either \(h=0\), or \(a(E)=c h(E)h(E+1)\) after changing \(h\) by a nonzero scalar, with \(c\in\mathbb C^*\).

**Proof.** Retain the zero solution. For \(h\ne0\), set
\[
r(E)=\frac{a(E)}{h(E)h(E+1)}\in\mathbb C(E).
\]
The displayed equation says \(r(E+1)=r(E)\). Rational periodicity makes \(r\) a nonzero constant. ∎

The reflected quantum equation gives \(c h(E)h(E-1)\), which is the same shifted-square class after replacing \(h(E)\) by \(h(E-1)\).

## Classical proof

Assume for contradiction that both extremes are non-squares. For the coefficient \(C_m\) of ladder degree \(m\) in \(\{D,X\}\), direct use of the chosen orientation gives
\[
C_m=\sum_{k+l=m}\bigl(k a_kb_l'-l a_k'b_l\bigr),
\qquad C_m=\delta_{m0}c_0.
\]

At \(m=4\), \(2(a_2b_2'-a_2'b_2)=0\), so \((b_2/a_2)'=0\) in \(\mathbb C(\tau)\) and \(b_2=\lambda_2a_2\). The reflected equation gives \(b_{-2}=\mu_2a_{-2}\).

At \(m=3\), put \(h=b_1-\lambda_2a_1\). The equation is
\[
2a_2h'-a_2'h=0.
\]
The classical homogeneous lemma and the non-square hypothesis force \(h=0\), hence \(b_1=\lambda_2a_1\). Reflection gives \(b_{-1}=\mu_2a_{-1}\).

At \(m=2\) and \(m=-2\), respectively,
\[
b_0'=\lambda_2a_0',\qquad b_0'=\mu_2a_0'.
\]
The \(m=1\) and \(m=-1\) equations then reduce to
\[
(\mu_2-\lambda_2)U=0,\quad U=a_2'a_{-1}+2a_2a_{-1}',
\]
\[
(\mu_2-\lambda_2)L=0,\quad L=a_{-2}'a_1+2a_{-2}a_1'.
\]
The central equation is exactly
\[
(\mu_2-\lambda_2)\bigl(2a_2a_{-2}+a_1a_{-1}\bigr)'=c_0.
\]
Thus \(\mu_2\ne\lambda_2\), so \(a_0'=b_0'=0\), \(U=L=0\), and
\(M=2a_2a_{-2}+a_1a_{-1}\) is linear nonconstant.

Now
\[
a_{-1}U=(a_2a_{-1}^2)',\qquad a_1L=(a_{-2}a_1^2)'.
\]
If \(a_{-1}\ne0\), then \(a_2a_{-1}^2\) is constant. It cannot be zero in the domain \(\mathbb C[\tau]\), and a nonzero constant product forces every factor to be a unit. Hence \(a_2\) is constant, contrary to its being a non-square. Therefore \(a_{-1}=0\). The same reasoning gives \(a_1=0\). Consequently \(2a_2a_{-2}=M\) is linear nonconstant. Degree additivity makes one extreme constant, again a square-class element. This contradiction proves J3. ∎

## Quantum proof

Again assume that both extremes are outside the shifted-square class. If \(Q_m\) is the ladder-degree-\(m\) coefficient of \([D,X]\), the multiplication convention gives the exact general formula
\[
Q_m(E)=\sum_{k+l=m}\left(b_l(E+k)a_k(E)-a_k(E+l)b_l(E)\right),
\qquad Q_m=\delta_{m0}.
\]

At \(m=4\),
\[
b_2(E+2)a_2(E)=a_2(E+2)b_2(E).
\]
Thus \(b_2/a_2\) is 2-periodic in \(\mathbb C(E)\), and rational periodicity gives \(b_2=\lambda_2a_2\). Reflection gives \(b_{-2}=\mu_2a_{-2}\).

At \(m=3\), with \(h=b_1-\lambda_2a_1\), the homogeneous equation is
\[
h(E+2)a_2(E)=a_2(E+1)h(E).
\]
The quantum lemma and the non-shifted-square hypothesis force \(h=0\). Reflection yields \(b_{-1}=\mu_2a_{-1}\).

At \(m=2\), \(b_0-\lambda_2a_0\) is a 2-periodic polynomial and hence a constant \(\gamma\). Independently, the \(m=-2\) equation gives
\[
b_0-\mu_2a_0=\gamma'
\]
for a possibly different constant \(\gamma'\). The \(m=1\) and \(m=-1\) equations are exactly
\[
(\mu_2-\lambda_2)U_q=0,
\quad U_q(E)=a_{-1}(E+2)a_2(E)-a_2(E-1)a_{-1}(E),
\]
\[
(\mu_2-\lambda_2)L_q=0,
\quad L_q(E)=a_{-2}(E+1)a_1(E)-a_1(E-2)a_{-2}(E).
\]
Define
\[
W_2(E)=a_2(E-2)a_{-2}(E),\qquad
W_1(E)=a_1(E-1)a_{-1}(E).
\]
The exact central identity is
\[
(\mu_2-\lambda_2)\bigl(W_2(E+2)-W_2(E)+W_1(E+1)-W_1(E)\bigr)=1.
\]
It first implies \(\mu_2\ne\lambda_2\). Only now do the two independent constant equations imply that \(a_0\), and then \(b_0\), is constant; and the cross equations imply \(U_q=L_q=0\).

Set
\[
G(E)=a_2(E)a_{-1}(E+1)a_{-1}(E+2),
\]
\[
H(E)=a_{-2}(E)a_1(E-1)a_1(E-2).
\]
Direct multiplication gives the exact identities
\[
G(E)-G(E-1)=a_{-1}(E+1)U_q(E),
\]
\[
H(E+1)-H(E)=a_1(E-1)L_q(E).
\]
Hence \(G\) and \(H\) are 1-periodic polynomials and therefore constants. If \(G\ne0\), all three polynomial factors are units, so \(a_2\) is constant, contradicting the hypothesis. If \(G=0\), the domain property and \(a_2\ne0\) force \(a_{-1}=0\). Thus in every allowed branch \(a_{-1}=0\); similarly \(a_1=0\).

Therefore \(W_1=0\), and the central identity becomes the exact formula
\[
W_2(E+2)-W_2(E)=\frac{1}{\mu_2-\lambda_2}\ne0.
\]
A polynomial with nonzero constant 2-step difference is linear nonconstant. Since \(W_2=a_2(E-2)a_{-2}(E)\), degree additivity makes one extreme constant, hence a shifted-square-class element. This contradiction proves J3q. ∎

## Polynomial-pair corollaries and limits

A genuine polynomial Keller pair in \(\mathbb C[x,\xi]\) that, in the above Laurent coordinates, has **both** entries supported in \([-2,2]\) satisfies J3. Likewise, a genuine Weyl pair \(D,X\in A_1\) with \([D,X]=1\) and the same common localized band presentation satisfies J3q.

These are corollaries of the Laurent/localized theorems, not identifications of all Laurent or localized pairs with Keller or Weyl pairs. In the quantum case, membership in \(A_1\) additionally requires
\[
E(E-1)\cdots(E-r+1)\mid a_{-r}(E)
\]
for each negative ladder \(-r\), and analogously for \(D\); no membership classification is asserted here.

The result excludes only the sector where both extreme coefficients are present and both are non-(shifted)-squares. It does not cover one-sided support, classify the remaining square sectors or all band-2 pairs, establish invertibility, or prove JC2/DC1.

## Computational support

`verify_J3.py` runs 13 exact checks of selected reduced classical and quantum identities and finite polynomial examples. In particular, it checks reduced ladder components, central and product identities, two bounded sample equations, and periodicity for polynomials of degree at most four. It does not verify the full unreduced coefficient systems, arbitrary-degree lemmas, branch completeness, localization or membership facts, or the completeness of either proof. The proofs above are self-contained independently of those finite checks.
