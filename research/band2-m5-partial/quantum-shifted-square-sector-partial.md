# Independent audit: quantum shifted-square sector in band 2

**INDEPENDENTLY DERIVED AND AUDITED — NOT PEER REVIEWED — PARTIAL CASCADE ONLY**

This memo does not establish a full quantum band-2 classification, JC2, or DC1. The resistant arbitrary-degree branch identified below remains unresolved, and no counterexample was found.

## Verdict

The quantum shifted-square sector is **not presently classified** by the supplied M5 material or by the derivation below. Two arbitrary-degree positive-ladder reductions can be proved cleanly:

1. after subtracting a scalar multiple of `X` from `D`, the full ladder-3 homogeneous freedom is
   \[
   b_1=\kappa h;
   \]
2. when \(\kappa\ne0\), ladder 2 forces the genuine polynomial divisibility
   \[
   h\mid a_1,
   \]
   and gives an exact finite-difference parametrization of \(b_0\).

These are the quantum analogues of the first classical square-sector reductions. They do **not** by themselves imply that nonconstant \(h\) divides a linear polynomial, do not classify the two-sided negative tail, and do not prove that constant \(h\) yields only tame families. The remaining equations form a nonlinear finite-difference system with at least the branches \(\kappa=0\), \(\kappa\ne0\), and the independent negative homogeneous parameter \(\mu\). No arbitrary-degree argument supplied to me eliminates the candidate \(h\in\mathbb C^*,\ \kappa\mu\ne0,\ s\ne0\) branch, nor establishes whether that locus is nonempty.

Accordingly, the defensible M5 conclusion is a **partial cascade with a major unresolved branch**, not a quantum band-2 theorem. Degree-bounded SymPy solving and isolated quantum examples cannot close that gap.

## 1. Setup and exact coefficient system

Work in
\[
A_1[x^{-1}]=\bigoplus_{k\in\mathbb Z}x^k\mathbb C[E],
\qquad
(x^af(E))(x^bg(E))=x^{a+b}f(E+b)g(E).
\]
Write \(f^{[r]}(E)=f(E+r)\). For
\[
X=\sum_{k=-2}^2x^ka_k(E),\qquad
D=\sum_{l=-2}^2x^lb_l(E),
\]
the ladder-\(m\) coefficient of \([D,X]\) is
\[
Q_m=\sum_{k+l=m}\bigl(b_l^{[k]}a_k-a_k^{[l]}b_l\bigr),
\qquad Q_m=\delta_{m0}.
\]
This fixes all shift and sign conventions used below.

Assume
\[
a_2=h(E)h(E+1),\qquad h\ne0,
\]
since \(a_2\ne0\). Thus the formal zero branch \(h=0\) is incompatible with the stated sector and must not be silently treated as a solution.

## 2. Ladder 4 and the retained ladder-3 parameter

From
\[
Q_4=b_2^{[2]}a_2-a_2^{[2]}b_2=0
\]
the rational quotient \(b_2/a_2\) is 2-periodic. A rational function in characteristic zero with a nonzero additive period is constant, so
\[
b_2=\lambda a_2.
\]
Replace \(D\) by \(D-\lambda X\). This preserves \([D,X]=1\) and gives the useful gauge
\[
b_2=0.
\]
At ladder 3, if \(c=b_1-\lambda a_1\) before gauging, then
\[
c(E+2)a_2(E)=a_2(E+1)c(E).
\]
Substitution of \(a_2=h h^{[1]}\) and cancellation in the polynomial domain gives
\[
c^{[2]}h=h^{[2]}c.
\]
Hence \(c/h\) is 2-periodic as a rational function and therefore constant. The complete homogeneous solution is
\[
\boxed{b_1=\kappa h}
\]
in the gauge \(b_2=0\), with \(\kappa\in\mathbb C\) arbitrary. This includes the essential zero branch \(\kappa=0\); dividing by \(\kappa\) before separating it loses solutions.

## 3. Ladder 2: an arbitrary-degree quantum divisibility lemma

Put \(B=b_0\). Ladder 2 is
\[
a_2(B^{[2]}-B)+\kappa\bigl(h^{[1]}a_1-h a_1^{[1]}\bigr)=0.
\]
With the rational function \(r=a_1/h\), division by \(h h^{[1]}\) gives
\[
B^{[2]}-B=\kappa(r^{[1]}-r). \tag{3.1}
\]

### Branch 3A: \(\kappa=0\)

Equation (3.1) says \(B^{[2]}=B\). Polynomial periodicity gives
\[
\boxed{b_0=\beta\in\mathbb C.}
\]
There is no conclusion that \(h\mid a_1\) in this branch.

### Branch 3B: \(\kappa\ne0\)

The left side of (3.1) is polynomial, so \(r^{[1]}-r\) is polynomial. A rational function whose unit finite difference is polynomial must itself be polynomial: otherwise a pole of \(r\) propagates along an infinite additive orbit when the poles of \(r(E+1)-r(E)\) are required to cancel, impossible for a rational function with finitely many poles. Therefore
\[
\boxed{a_1=h p\quad\text{for some }p\in\mathbb C[E].}
\]
Equation (3.1) then becomes
\[
(T^2-1)B=\kappa(T-1)p,
\]
where \(T f=f^{[1]}\). Thus
\[
(T-1)\bigl((T+1)B-\kappa p\bigr)=0,
\]
and polynomial periodicity yields
\[
\boxed{B(E+1)+B(E)=\kappa p(E)+\gamma,\qquad \gamma\in\mathbb C.} \tag{3.2}
\]
The operator \(T+1\) is invertible on \(\mathbb C[E]\), so (3.2) determines one polynomial \(B\) for each \((p,\gamma)\). This is the precise quantum replacement for the classical midpoint formula; replacing it by \(b_0=\kappa p/2+\beta\) is generally false when \(p\) is nonconstant.

This establishes a real finite-difference divisibility theorem, but it is \(h\mid a_1\), not yet “\(h\) divides a linear polynomial.”

## 4. Exact remaining cascade

The following gives an arbitrary-degree reduction without pretending to solve the nonlinear tail. Retain the gauge
\[
a_2=h h^{[1]},\quad b_2=0,\quad b_1=\kappa h,
\]
and write
\[
u=a_{-1},\quad v=b_{-1},\quad s=a_{-2},\quad w=b_{-2}.
\]
Then ladder 1 is exactly
\[
h^{[1]}v^{[2]}-h^{[-1]}v
+p(B^{[1]}-B)+\kappa(a_0-a_0^{[1]})=0 \tag{4.1}
\]
in the branch \(\kappa\ne0\), where \(a_1=hp\). Without that branch assumption, replace \(p\) by \(a_1/h\) only as a rational expression, or retain the unreduced polynomial equation
\[
a_2v^{[2]}-a_2^{[-1]}v+B^{[1]}a_1-a_1B
+\kappa h a_0-\kappa a_0^{[1]}h=0.
\]

At the negative extreme, the exact equation is
\[
w^{[-2]}s-s^{[-2]}w=0. \tag{4.2}
\]
The remaining four equations, before any division or proportionality substitution, are
\[
w^{[-1]}u+v^{[-2]}s-u^{[-2]}w-s^{[-1]}v=0, \tag{4.3}
\]
\[
w a_0+v^{[-1]}u+B^{[-2]}s
-a_0^{[-2]}w-u^{[-1]}v-sB=0, \tag{4.4}
\]
\[
w^{[1]}a_1+v a_0+B^{[-1]}u+\kappa h^{[-2]}s
-a_1^{[-2]}w-a_0^{[-1]}v-uB-\kappa s^{[1]}h=0, \tag{4.5}
\]
and the central normalization
\[
w^{[2]}a_2+v^{[1]}a_1+\kappa h^{[-1]}u
-a_2^{[-2]}w-a_1^{[-1]}v-\kappa u^{[1]}h=1. \tag{4.6}
\]
Equations (4.1)--(4.6), together with (3.2) in the \(\kappa\ne0\) branch, are the exact reduced quantum system. They preserve constants, units, and the zero branches.

If \(s\ne0\), equation (4.2) and rational periodicity imply
\[
\boxed{w=\mu s},\qquad \mu\in\mathbb C. \tag{4.7}
\]
Here \(\mu=0\) is a genuine branch, and substitution into (4.3)--(4.6) gives the proportional negative-tail system. If \(s=0\), equation (4.2) imposes no condition on \(w\); the branch \(s=0,w\ne0\) must therefore be retained separately rather than absorbed into (4.7).

## 5. Membership constraints and what they do not prove

Genuine Weyl membership requires
\[
E\mid u,v,
\qquad E(E-1)\mid s,w.
\]
In the branch \(s\ne0\), equation (4.7) makes the membership condition for \(w\) automatic once it holds for \(s\), including \(\mu=0\). In the separate branch \(s=0,w\ne0\), membership of \(w\) must be imposed independently. These divisibilities belong to the full remaining system; they cannot be inferred from isolated localized examples.

A simple quantum boundary family illustrates the polar obstruction but does not classify the two-sided sector. Take \(p=0\), \(B\) and \(a_0\) constant, \(v=s=w=0\), \(u=A\) constant, and linear \(h\). Then (4.1) is satisfied, and the central equation becomes
\[
\kappa A\bigl(h(E-1)-h(E)\bigr)=1.
\]
If \(h(E)=cE+d\) with \(c\ne0\), choosing \(A=-1/(\kappa c)\) gives a localized solution, and it fails \(E\mid u\). But it has \(a_{-2}=0\), so it is a boundary example, not evidence that every nonconstant-\(h\), two-sided solution is polar or excluded by membership.

No derivation above forces a nonconstant \(h\) to divide an affine polynomial. Such a result, if true, has to emerge from (4.1), (4.3)--(4.6) after all zero and unit branches are separated. The classical integrated product identity does not transfer verbatim: the quantum equations involve staggered shifts and the quantum midpoint operator \((T+1)^{-1}\).

## 6. Constant \(h\)

If \(h\in\mathbb C^*\), scalar normalization can make \(h=1\). Then \(a_2=1\), and in the \(\kappa\ne0\) branch the positive cascade reduces to
\[
a_1=p,
\qquad B^{[1]}+B=\kappa p+\gamma,
\]
followed by the nonlinear tail (4.1), (4.3)--(4.6). This is a polynomial system with arbitrary-degree input unless further lemmas bound degrees; its nonempty components and their dimensions have not been determined.

Two familiar families can be checked exactly inside this branch:

- a swapped affine/shear family;
- the family obtained from a linear Weyl generator \(Z=x+r\partial\), with \(X\) polynomial in \(Z\) plus \(\partial\) and \(D=-Z+\beta\).

Exact verification proves that these displayed families satisfy the commutator relation. It does **not** prove that they exhaust all solutions of the reduced tail. In particular, no arbitrary-degree proof supplied here eliminates the branch
\[
 h=1,\qquad \kappa\ne0,\qquad \mu\ne0,\qquad s\ne0.
\]
This is the quantum counterpart of the resistant classical branch and must be listed as unresolved.

## 7. Assessment of machine evidence

The available M5 checker is predominantly classical. Its quantum section checks only:

1. one swapped-shear family;
2. one \(Z\)-family;
3. one bounded instance of the ladder-2 divisibility phenomenon, with a particular linear \(h\) and degree-two ansatz.

Those checks are useful regression tests. They do not verify the general ladder-3 lemma, the arbitrary-degree rational-pole argument in Section 3, the reduced tail (4.3)--(4.6), Weyl membership across all branches, nonconstant-\(h\) exclusion, or constant-\(h\) exhaustion.

Likewise, a solver at degree budgets one or two cannot establish an arbitrary-degree classification. Conditional/free SymPy solution dictionaries require saturation by the nonvanishing conditions; absence of a returned dictionary is not a proof of an empty algebraic branch. A script that catches a degree-two exception and still prints overall success does not even certify its own advertised bounded degree-two test.

## 8. Final classification of claims

- **Top proportionality and full ladder-3 homogeneous parameter:** proved for arbitrary degree.
- **Quantum ladder-2 implication \(h\mid a_1\) when \(\kappa\ne0\):** proved for arbitrary degree.
- **Quantum midpoint equation for \(b_0\):** proved exactly.
- **\(\kappa=0\) branch:** retained; it does not yield \(h\mid a_1\).
- **Negative proportionality \(b_{-2}=\mu a_{-2}\):** proved for arbitrary degree in the branch \(a_{-2}\ne0\), including \(\mu=0\); the separate branch \(a_{-2}=0,b_{-2}\ne0\) remains explicit.
- **Complete reduced finite-difference tail:** displayed exactly before branch specialization, not solved.
- **Nonconstant \(h\) forces affine \(h\) and only polar localized families:** not proved.
- **Membership kills every nonconstant-\(h\), two-sided family:** not proved.
- **Constant \(h\) gives exactly the tame families:** examples verified, exhaustion not proved.
- **Unresolved candidate branch:** \(h=1,\ \kappa\mu s\ne0\) subject to (4.1), (4.3)--(4.6); nonemptiness and dimension are not established.
- **Degree-bounded computation as exhaustive evidence:** invalid.

No counterexample was found.