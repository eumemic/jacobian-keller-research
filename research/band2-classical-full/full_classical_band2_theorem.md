# Full classical band-2 theorem

**INDEPENDENTLY ASSEMBLED AND INTERNALLY AUDITED — NOT PEER REVIEWED — BAND-SCOPED**

## Theorem and exact scope

Work over \(\mathbb C\), put \(t=x\xi\), and use
\[
\{G,F\}=G_\xi F_x-G_xF_\xi.
\]
Let
\[
F=\sum_{k=-2}^{2}x^ka_k(t),\qquad
G=\sum_{k=-2}^{2}x^kb_k(t)
\]
be elements of \(\mathbb C[x,\xi]\), not merely of
\(\mathbb C[x,x^{-1},\xi]\), and suppose \(\{G,F\}=1\). Equivalently, in
this presentation the genuine polynomial membership conditions hold:
\[
t\mid a_{-1},b_{-1},\qquad t^2\mid a_{-2},b_{-2}.
\]
Thus the ladder support of **each** entry is contained in \([-2,2]\); an
extreme coefficient is allowed to vanish.

Then \((F,G)\) is a polynomial coordinate pair. More explicitly, unless the
pair is already in band 1, the orientation operations in Section 3 and a
diagonal symplectic scaling produce an equivalent transformed pair
\((\widetilde F,\widetilde G)\) for which there are constants
\(\kappa\ne0,c_0,c_1,A,\beta,\lambda\in\mathbb C\) such that
\[
\widetilde U=x+c_0+c_1\xi,
\qquad
\widetilde F=\widetilde U^2-\frac1\kappa\xi-A,
\qquad
\widetilde G=\lambda\widetilde F+\kappa\widetilde U+\beta.
\tag{0.1}
\]
The original unoriented pair need not literally have this displayed form;
undoing the reversible operations returns it. In all cases the map
\((x,\xi)\mapsto(F,G)\) is a polynomial automorphism of \(\mathbb A^2_\mathbb C\).

This is a theorem about the fixed support-contained band-2 class. It does
**not** claim the two-dimensional Jacobian conjecture (JC2), because a general
Keller pair need not admit this band bound in the fixed \(t=x\xi\) ladder
decomposition.

## 1. Frozen identities and operations

For one-term pieces,
\[
\{x^\ell b(t),x^ka(t)\}
=x^{k+\ell}(ka b'-\ell a'b).
\tag{1.1}
\]
The complete nine coefficient equations and their signs are recorded and
proved in the M4 memo cited below.

We use only the following elementary transformations.

1. **Pair exchange.** If \(\{G,F\}=1\), put
   \[
   (F^\ast,G^\ast)=(G,-F).
   \]
   Then \(\{G^\ast,F^\ast\}=\{-F,G\}=1\). This preserves polynomial
   membership and support containment. It also preserves the property of being
   a polynomial coordinate pair, since it is an invertible linear change of
   the target coordinates.

2. **Reflection.** Define
   \[
   (Rf)(x,\xi)=f(\xi,x).
   \]
   If \(f=\sum_kx^ka_k(t)\), then \(t=x\xi\) is fixed and
   \[
   Rf=\sum_k x^k\,t^{-k}a_{-k}(t),
   \qquad
   \boxed{(Rf)_k=t^{-k}a_{-k}}.
   \tag{1.2}
   \]
   Genuine polynomial membership makes (1.2) polynomial for every
   \(k\in[-2,2]\), and reflection preserves that support interval. The linear
   substitution has determinant \(-1\), so with the present Poisson convention
   it gives
   \[
   \{RG,RF\}=-R\{G,F\}.
   \tag{1.3}
   \]
   Consequently the oriented reflected pair is \((RF,-RG)\), because
   \(\{-RG,RF\}=1\). Reflection is an automorphism of the source polynomial
   ring, so it also preserves coordinate-pair status.

3. Target shears and translations, together with source diagonal symplectic
   scalings \((x,\xi)\mapsto(\alpha x,\alpha^{-1}\xi)\), are polynomial
   automorphisms and preserve the relevant bracket, support, and membership.

## 2. Independent band-1 theorem

### Lemma (classical band-1 classification)

Suppose both supports are contained in \([-1,1]\) and \(\{G,F\}=1\). Then the
pair is affine symplectic up to triangular target changes, hence is a polynomial
automorphism.

### Proof

Write
\[
F=xa_1(t)+a_0(t)+x^{-1}a_{-1}(t),\qquad
G=xb_1(t)+b_0(t)+x^{-1}b_{-1}(t),
\]
with \(t\mid a_{-1},b_{-1}\). The five equations obtained from (1.1) are
\[
\begin{aligned}
&a_1b_1'-a_1'b_1=0,\\
&a_1b_0'-a_0'b_1=0,\\
&a_1b_{-1}'+a_1'b_{-1}-a_{-1}b_1'-a_{-1}'b_1=1,\\
&a_0'b_{-1}-a_{-1}b_0'=0,\\
&a_{-1}'b_{-1}-a_{-1}b_{-1}'=0.
\end{aligned}
\tag{2.1}
\]

If \(a_1\ne0\), the top equation gives \(b_1=\lambda a_1\). Replace \(G\) by
\(G-\lambda F\), so \(b_1=0\). The second equation gives \(a_1b_0'=0\), hence
\(b_0\) is constant. The fourth equation then gives \(a_0'b_{-1}=0\), and the
central equation becomes
\[
(a_1b_{-1})'=1.
\tag{2.2}
\]
In particular \(b_{-1}\ne0\), so \(a_0\) is constant. Since
\(a_1b_{-1}=t+c\) and \(t\mid b_{-1}\), evaluation at \(t=0\) gives \(c=0\).
A product of nonzero polynomials equal to \(t\) has degree split \((0,1)\) or
\((1,0)\); divisibility forces
\[
a_1=\alpha\in\mathbb C^*,\qquad b_{-1}=t/\alpha.
\]
The bottom equation gives \(a_{-1}/b_{-1}=\mu\in\mathbb C\). Therefore
\[
F=\alpha x+c_0+(\mu/\alpha)\xi,
\qquad
G=\lambda F+c_1+(1/\alpha)\xi,
\]
an affine polynomial automorphism with bracket one.

If \(a_1=0\) but \(b_1\ne0\), pair exchange reduces to the preceding case. If
\(a_1=b_1=0\), apply the oriented reflection from Section 1 when at least one
negative coefficient is nonzero. If all four \(\pm1\) coefficients vanish,
then both entries depend only on \(t\), so their bracket is zero, contrary to
\(\{G,F\}=1\). These cases exhaust band 1. \(\square\)

## 3. Pair-exchange/reflection orientation lemma

### Lemma

If a support-contained band-2 Keller pair is not in band 1, the operations in
Section 1 produce an oriented pair, still denoted \((F,G)\), with
\(a_2\ne0\).

### Proof

Some coefficient among \(a_2,b_2,a_{-2},b_{-2}\) is nonzero. If it belongs to
\(G\), use pair exchange. If it has index \(-2\), use the oriented reflection;
formula (1.2) sends it to a nonzero index-\(2\) polynomial coefficient. The
Poisson signs in (1.3) and pair exchange retain \(\{G,F\}=1\). \(\square\)

## 4. Assembly of M4 and M5

Take a pair not covered by the band-1 lemma and orient it by Section 3, so
\(a_2\ne0\).

The proved M4 theorem
[`../band2-classical-proved/M4_proof_memo.md`](../band2-classical-proved/M4_proof_memo.md)
applies under exactly the present bracket, support-containment, and genuine
membership hypotheses. Over \(\mathbb C\), it proves that \(a_2\) is a square
in \(\mathbb C[t]\). Hence
\[
a_2=h^2\ne0.
\]

The proved M5 classical theorem
[`../band2-m5-partial/classical-square-sector-partial.md`](../band2-m5-partial/classical-square-sector-partial.md)
then applies under the same hypotheses. Its arbitrary-degree argument first
proves that the square root \(h\) is constant and then, after its stated
reversible normalizations, classifies the pair by the normal form (0.1).

That normalized form has an explicit polynomial inverse. From coordinates
\((P,Q)=(\widetilde F,\widetilde G)\),
\[
\widetilde U=\frac{Q-\lambda P-\beta}{\kappa},
\qquad
\xi=\kappa(\widetilde U^2-P-A),
\qquad
x=\widetilde U-c_0-c_1\xi.
\tag{4.1}
\]
Undoing pair exchange, reflection, target shears/translations, and diagonal
scaling composes this inverse with polynomial automorphisms. Thus every pair in
the theorem is a polynomial automorphism.

Together with the independent band-1 lemma, this exhausts every polynomial
Keller pair over \(\mathbb C\) satisfying the support-containment and genuine
membership hypotheses stated above. \(\square\)

## 5. Exact regression verifier

Run from the repository root:

```sh
python3 research/band2-classical-full/verify_full_classical_band2.py
```

The exact SymPy script checks the band-1 component identities and target
shear, representative examples of the reflection coefficient mapping and
Poisson sign, and the M5 normal-form bracket and inverse. These checks are
useful regression support only. They do
not establish theorem completeness; completeness is the written orientation
split followed by the proved M4 and M5 theorems and the independent band-1
argument above.
