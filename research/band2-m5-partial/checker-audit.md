# Independent adversarial audit of the M5 checker

**INDEPENDENTLY AUDITED — NOT PEER REVIEWED — NOT A CLASSIFICATION CERTIFICATE**

The supplied checker is not included or endorsed by this package. This audit preserves the useful identities it tests while rejecting any inference of arbitrary-degree exhaustion or full band-2 classification. It proves neither JC2 nor DC1, and no counterexample was found.

## Executive verdict

The checker does **not** establish the advertised classical band-2 square-sector classification. It contains useful symbolic regression identities, but it conflates those identities with completeness. The resistant case-A branch is tested only through inconsistent low-degree ansatzes, with no theorem bounding degrees and no complete algebraic emptiness certificate. The degree-2 test may be skipped without recording a failure. The displayed tame and polar Poisson checks also use the wrong sign relative to the checker’s own definition of the Poisson bracket.

The strongest defensible conclusions are:

1. the component formulas in S1 and S4 are correct for the stated reduced parametrizations and sign convention;
2. the two displayed tame families are genuine constant-bracket examples, but with classical constants opposite to those asserted;
3. the polar example does violate polynomial membership, but its bracket likewise has the opposite sign;
4. S5 proves neither arbitrary-degree exhaustion nor even a sound, certified algebraic exhaustion of all branches at its advertised degree budgets;
5. S6 consists of two exact quantum examples and one very small generic-parameter divisibility instance, not a quantum classification or meaningful spot check of the resistant tail.

At the time of this audit, the checker-based classification was therefore blocked at least by the arbitrary-degree branch
\[
h=1,\qquad \widetilde\mu\ne0,\qquad a_{-2}\ne0,\qquad c\ne0,
\]
and by the missing logical bridge from J3's statement that, for the designated member \(X\) in the two-sided case \(a_2a_{-2}\ne0\), at least one of \(a_2,a_{-2}\) lies in the square class (equivalently over \(\mathbb C\), is a literal square after absorbing a nonzero scalar) to the exhaustive oriented case split assumed by this checker.

**Current status note.** The companion classical memo now closes this resistant branch, and all of its zero subbranches, by a later independent arbitrary-degree proof. That result does not rehabilitate the supplied checker: the bounded ansatz, sign defects, silent skip, and absent orientation/case bridge documented here remain exactly as audited. The new theorem is scoped only to the oriented classical \(a_2\)-square sector under the companion memo's exact support and membership hypotheses; the quantum resistant branch and broader band-2 case coverage remain unresolved.

## 1. Independent Poisson derivation and orientation

For
\[
F=x^i f(\tau),\qquad G=x^j g(\tau),\qquad \tau=x\xi,
\]
one has
\[
\partial_xF=x^{i-1}(if+\tau f'),\qquad
\partial_\xi F=x^{i+1}f'.
\]
With the checker’s definition
\[
\operatorname{PB}(D,X)=D_xX_\xi-D_\xi X_x,
\]
the ladder-\(m\) component is therefore
\[
[\operatorname{PB}(D,X)]_m
 =\sum_{i+j=m}\bigl(jb_j a_i'-i b_j'a_i\bigr).
\]
This independently fixes every sign used below.

For example, if \(D=\kappa x+\beta\) and \(X=x^2+q_1x+q_0+A_1\xi\), then
\[
\operatorname{PB}(D,X)=\kappa A_1,
\]
not \(-\kappa A_1\). If \(Z=x+r\xi\), \(D=-Z+\beta\), and
\(X=Z^2+q_1Z+q_0+\xi\), then
\[
\operatorname{PB}(D,X)=-1,
\]
not \(+1\).

The quantum commutator has the opposite semiclassical orientation here: with \([x,\partial]=-1\), the checker correctly obtains \([\kappa x,A_1\partial]=-\kappa A_1\) and \([-Z,\partial]=+1\). The classical S2/S3 expectations appear to have copied this quantum orientation despite defining the classical bracket as \(D_xX_\xi-D_\xi X_x\).

## 2. S1: case-B reduced parametrization

### What is correctly checked

Under
\[
a_2=h^2,\quad a_1=hp,\quad a_0=p^2/4+hq+\delta,
\]
\[
b_2=0,\quad b_1=\kappa h,\quad b_0=\kappa p/2+\beta,
\quad b_{-1}=\kappa q/2,
\]
the script correctly verifies the stated vanishing positive components and the residual formulas. In particular, with
\[
I=h\left(2h b_{-2}+\frac\kappa2pq-\kappa a_{-1}\right),
\]
the independent component formula gives
\[
[\operatorname{PB}(D,X)]_0=-I'.
\]
The signs in S1 are internally correct.

The checks are performed with generic quadratic polynomials. Since these identities involve only polynomial arithmetic and first derivatives, this is good evidence for, and readily generalizes to, formal differential identities of arbitrary polynomials. It is not, however, evidence that every pair has the reduced parametrization.

If the constant bracket is nonzero, integration gives \(I=c\tau+e\), up to the sign chosen for the bracket constant. Because \(h\mid I\), a nonconstant polynomial \(h\) must divide a nonzero affine polynomial and hence must be affine. This particular degree argument is valid once all its hypotheses have actually been established.

### What is not checked

S1 does not prove:

- that every relevant pair falls into case B;
- that the displayed parametrization follows in every zero/nonzero parameter branch;
- that J3’s one-square conclusion allows the orientation assumed here;
- that the residual equations have only the later displayed families;
- completeness in arbitrary coefficient degree.

Calling a generic degree-2 identity check “fully generic” is harmless only for the identity itself, not for classification.

## 3. S2: polar family and membership obstruction

Substituting the displayed polar family into the independently derived central identity gives
\[
I=\sigma h,
\qquad
\operatorname{PB}(D,X)=-\sigma h_1.
\]
Thus S2’s assertion \(+\sigma h_1\) is false under the script’s own `PB` definition. Running the relevant sections produces a failure here.

The membership observation is mathematically correct under the declared nonvanishing assumptions: \(a_{-1}=-\sigma/\kappa\) is a nonzero constant, whereas polynomial membership at ladder \(-1\) requires \(\tau\mid a_{-1}\). In ordinary \(x,\xi\) variables, the offending term is a genuine localized \(x^{-1}\) term.

The implementation does not truly prove nonvanishing: `simplify(-sig/kap) != 0` is a syntactic symbolic comparison. The mathematical statement is valid only after separately assuming \(\sigma\kappa\ne0\).

Most importantly, one explicit nonmember polar family is not a proof that all nonconstant-\(h\) solutions are polar, nor that membership eliminates every nonconstant-\(h\), two-sided branch. That reduction is absent from the checker.

## 4. S3: tame families

Both displayed families are genuine polynomial constant-bracket pairs and satisfy the required negative-ladder membership:

- the \(A_1\xi\) term has localized coefficient \(a_{-1}=A_1\tau\), divisible by \(\tau\);
- in the \(Z\)-family, the \(\xi^2\) term gives a ladder \(-2\) coefficient proportional to \(\tau^2\), and ladder \(-1\) coefficients are divisible by \(\tau\).

But both asserted classical signs are wrong:
\[
\operatorname{PB}(\kappa x+\beta,
 x^2+q_1x+q_0+A_1\xi)=+\kappa A_1,
\]
\[
\operatorname{PB}(-Z+\beta,
 Z^2+q_1Z+q_0+\xi)=-1.
\]
Consequently the relation between \(\kappa A_1\) and the chosen bracket constant must also be sign-adjusted.

Even after that repair, S3 verifies only sufficiency: these are examples. It contains no uniqueness or exhaustion argument. The label “the two tame families” is therefore a classification claim unsupported by this section.

## 5. S4: case-A reduced parametrization

With \(b_{-2}=\widetilde\mu a_{-2}\), the ladder \(-4\) component vanishes identically, and the script’s formulas for ladders \(-3,-2,-1\) agree with the independent general component formula. Its central antiderivative
\[
E_{0,A}=2\widetilde\mu h^2a_{-2}
 +\frac\kappa2 hpq-\kappa a_{-1}h
\]
satisfies
\[
[\operatorname{PB}(D,X)]_0=-E_{0,A}'.
\]
These are valid symbolic identities.

As in S1, using generic quadratics checks formulas but not that every pair reaches this reduction, that all divisions used upstream are valid in every branch, or that the remaining equations have been solved. S4 should be described as a reduced-system regression test, not a classification step.

## 6. S5: resistant branch

### 6.1 No arbitrary-degree implication

S5 supplies no lemma bounding
\(\deg a_1\), \(\deg q\), \(\deg a_{-1}\), or \(\deg a_{-2}\). A solve at one or two selected ansatz sizes cannot exhaust arbitrary-degree polynomial solutions. The equations are nonlinear differential-polynomial equations, and no leading-term argument, recurrence, resultant, Gröbner basis, or projective-at-infinity analysis is supplied to show that high-degree branches reduce to these budgets.

After eliminating \(a_{-1}\) using the integrated central equation, one still obtains a nonlinear arbitrary-degree system. Writing \(p=a_1\), \(q=q\), \(v=a_{-2}\), and
\[
a_{-1}=\frac{2\widetilde\mu v}{\kappa}
 +\frac{pq}{2}-\frac{c\tau+e}{\kappa},
\]
the remaining three equations contain terms such as
\(v v'\), \(p q v'\), \(p v q'\), and \(qv p'\). Nothing in S5 bounds their polynomial degrees. The absence of a returned low-degree solution is therefore not evidence that an arbitrary-degree branch is empty.

In this historical checker audit, I did not obtain an explicit high-degree counterexample with
\(\widetilde\mu a_{-2}c\ne0\). That was not needed to invalidate the checker-based proof: the checker bore the burden of proving exhaustion, and it had no degree theorem. At that stage the resistant branch was unresolved rather than disproved; the companion classical memo now excludes it by a separate arbitrary-degree argument.

### 6.2 The advertised “degree budgets” are inconsistent

For `deg == 1`, the ansatz is
\[
\deg a_1,\deg q,\deg a_{-1}\le1,
\qquad a_{-2}=S_0\tau^2.
\]
For `deg == 2`, it is
\[
\deg a_1,\deg q,\deg a_{-1}\le2,
\qquad a_{-2}=S_0\tau^2+S_1\tau^3.
\]
Thus:

- if “budget” means the degree of the full coefficient, \(a_{-2}\) is allowed degrees 2 and 3 while the other coefficients are allowed degrees 1 and 2;
- if it means the degree after factoring the mandated \(\tau^2\), the quotient budgets are 0 and 1, not 1 and 2;
- a genuine quotient-degree-2 budget would require
  \(a_{-2}=\tau^2(S_0+S_1\tau+S_2\tau^2)\), including a \(\tau^4\) term.

Hence the truncation does not match any uniform degree-budget interpretation.

### 6.3 `sympy.solve` is not an exhaustion certificate

The script calls general-purpose `sympy.solve(..., dict=True)` on a nonlinear parametric system and then treats its returned dictionaries as a complete decomposition. SymPy does not promise that this call returns every component of an arbitrary nonlinear parametric variety in a form suitable for proof.

Specific missing safeguards include:

- no Gröbner-basis or elimination certificate;
- no saturation by \(\kappa\widetilde\mu c\) and by the ideal expressing \(a_{-2}\ne0\);
- no proof that singular branches lost through generic divisions inside `solve` are restored;
- no recording or checking of denominator/nonvanishing side conditions;
- no substitution of each dictionary back into every coefficient equation;
- no separate satisfiability check for free or conditional branches;
- no independent completeness count or algebraic-dimension analysis.

The degree-1 run returns 14 dictionaries, many duplicated and many leaving substantial free parameters. The post-filter merely asks whether the dictionary makes \(\widetilde\mu\), \(c\), or both displayed coefficients of \(a_{-2}\) identically zero. That filter is conservative for a straightforward free symbol, but it cannot repair missing solution components or implicit conditions from the solver. It proves only a property of the dictionaries that happened to be returned.

### 6.4 Degree 2 can be skipped while the script still passes

The degree-2 solve is wrapped in

```python
try:
    ...
    check(...)
except Exception:
    print("[INFO] degree-2 exhaustion skipped ...")
```

The exception branch appends no failed check to `RES`. Therefore an exception in the advertised degree-2 exhaustion permits the final line to print `ALL M5 CHECKS PASSED`, provided the unrelated checks pass. This is a direct control-flow defect.

External timeout or process termination is not caught, but ordinary SymPy exceptions are. In this environment the full solve is sufficiently expensive to run for minutes, itself illustrating that this is not a stable certificate. Separately, the current checker already fails S2 and S3 because of the sign errors, so its advertised all-PASS output is not reproduced here.

## 7. Completeness relative to J3

The supplied context says J3 establishes only that, for the designated member \(X\) in the two-sided case \(a_2a_{-2}\ne0\), at least one of \(a_2,a_{-2}\) lies in the square class (equivalently over \(\mathbb C\), is a literal square after absorbing a nonzero scalar). That is weaker than an exhaustive oriented classification. The checker assumes specific positive-square parametrizations and then splits according to a negative proportionality branch, but it does not prove:

1. that the square extreme can always be placed at the positive side without changing membership, bracket orientation, or the case equations;
2. that every possible zero branch in the reductions has been retained;
3. that the two reduced cases cover all consequences of the one-square statement;
4. that solutions of the residual equations are exactly polar or one of the two tame families.

Accordingly, J3 cannot be cited as filling the checker’s completeness gap.

## 8. Claimed tame classification and polar obstruction

The checker establishes the following limited facts:

- two displayed tame families are polynomial Keller examples, after correcting classical signs;
- one displayed nonconstant-\(h\) polar localized family violates membership.

It does **not** establish:

- every constant-\(h\) solution is tame;
- every nonconstant-\(h\) solution is polar;
- every nonconstant-\(h\), two-sided solution violates membership;
- no high-degree resistant case-A family exists.

The polar obstruction is family-specific. To become a classification obstruction, one needs a theorem deriving the polar form from the full residual equations, including all zero/nonzero cases. No such theorem is encoded or tested.

## 9. S6: exactly how little the quantum checks establish

The quantum multiplication routine implements
\[
(x^af(E))(x^bg(E))=x^{a+b}f(E+b)g(E),
\]
so its two family checks are exact for the displayed finite dictionaries. They show only:

1. the quantum swapped-shear example has commutator \(-\kappa A_1\);
2. the displayed quantum \(Z\)-family has commutator \(+1\).

They do not check a generic band-2 ansatz, negative-tail membership, case-A/case-B completeness, or exhaustion.

The divisibility test fixes \(h(E)=E+1\), takes \(a_1\) and \(b_0\) of degree at most two, and solves one ladder-2 equation. For generic symbolic \(\kappa\), SymPy returns
\[
d_0=d_1-d_2,
\]
so
\[
a_1(E)=(E+1)(d_2E+d_1-d_2).
\]
That is a valid small regression instance. It does not prove arbitrary-degree divisibility, does not test a general \(h\), and does not cover \(\kappa=0\): because \(\kappa\) is treated as a generic symbolic coefficient rather than as an unknown with an explicit branch split, the solve effectively uses the generic \(\kappa\ne0\) branch. At \(\kappa=0\), the displayed equation constrains \(b_0\) but does not force \(h\mid a_1\).

Even a general ladder-2 divisibility theorem would only be an early cascade lemma. It would not solve the shifted nonlinear negative tail or prove a quantum classification.

## 10. Claim-by-claim disposition

- **S1 identities:** verified, with the checker’s classical sign convention. They are reduced-system identities, not completeness.
- **S1 affine-\(h\) divisibility consequence:** valid under nonzero constant bracket and the established integrated identity.
- **S2 bracket sign:** false as written; the result is \(-\sigma h_1\).
- **S2 membership failure:** correct under \(\sigma\kappa\ne0\), but only for the displayed family.
- **S3 swapped shear:** constant-bracket example is correct; asserted classical sign is wrong.
- **S3 Z-family:** constant-bracket example is correct; asserted classical sign is wrong.
- **S3 “two tame families” as an exhaustive classification:** unproved.
- **S4 identities:** verified reduced-system identities, not exhaustion.
- **S5 bounded degree 1:** at most an uncertified report about returned SymPy dictionaries.
- **S5 bounded degree 2:** expensive/unstable, inconsistently truncated, and allowed to skip without failure.
- **S5 arbitrary-degree resistant branch:** not addressed by the supplied checker; subsequently closed in the companion classical memo by an independent arbitrary-degree proof.
- **S6 two quantum examples:** exact example checks only.
- **S6 deformed divisibility:** one linear-\(h\), degree-2, generic-\(\kappa\) instance only.
- **Overall classical classification by the supplied checker:** not established.
- **Overall quantum classification by the supplied checker:** not remotely established by these spot checks.

## 11. Historical blockers to a checker-based classification claim

The following were required before the supplied checker could support such a claim. The later companion proof addresses the classical oriented \(a_2\)-square mathematics independently; it does not repair these checker defects.

1. Supply a rigorous exhaustive case split from J3’s square-class disjunction, preserving orientation and membership.
2. Prove all reductions without losing zero branches or state the necessary branch hypotheses explicitly.
3. Prove an arbitrary-degree bound for the resistant \(h=1\), \(\widetilde\mu a_{-2}c\ne0\) branch, or solve it at arbitrary degree.
4. If bounded algebra is retained, use a reproducible algebraic certificate: appropriate saturation, elimination/Gröbner data, explicit side conditions, and branch verification.
5. Define degree budgets consistently with the forced factors \(\tau\) and \(\tau^2\).
6. Make any skipped advertised test a hard failure.
7. Correct the classical Poisson signs and the associated normalization constants.
8. Prove that the residual solutions are exhausted by the polar and tame families; checking those families is not enough.
9. Describe the quantum section as examples/regression tests unless the full arbitrary-degree shifted system is actually analyzed.

## Final assessment

This is a useful symbolic identity checker embedded in an invalid classification certificate. Its core S1/S4 component algebra is valuable, and its quantum examples are legitimate. But the checker's central completeness claim fails for both mathematical and computational reasons. At the time supplied, its resistant arbitrary-degree branch was a decisive proof blocker; the low-degree solver is not a complete algebraic decision procedure, its degree budgets are misaligned with membership, and degree 2 can be silently omitted. The classical sign assertions in S2/S3 are also directly false under the checker’s own bracket definition. The companion classical memo's later independent proof supersedes only the branch's mathematical open status within its stated oriented \(a_2\)-square sector; it does not alter any historical finding about this checker.
