# A concrete counterexample-oriented attack on JC2

> **PROVISIONAL RESEARCH MEMO — NOT PEER REVIEWED.** This document records possible tests and research directions. It does not claim a theorem, counterexample, proof of JC2, or independently validated result.

**Scope.** Work over a characteristic-zero field, reduced to \(\mathbb C\). “Keller” means \(J(P,Q)=P_xQ_y-P_yQ_x\in\mathbb C^*\). This memo proposes tests, not a claimed breakthrough.

## 1. Reductions that genuinely constrain a JC2 counterexample

### Verified facts

1. **Geometric/field criteria.** A Keller map is étale. If it is proper, it is an automorphism; if \([\mathbb C(x,y):\mathbb C(P,Q)]=1\) (the birational case), it is an automorphism. Thus a counterexample must be nonproper and have geometric degree \(>1\). The birational result goes back to Keller; see Bass–Connell–Wright (BCW), *Bull. AMS* 7 (1982), 287–330, DOI [10.1090/S0273-0979-1982-15032-7](https://doi.org/10.1090/S0273-0979-1982-15032-7).

2. **Degree reduction is stable, not plane-preserving.** BCW reduce the general JC to maps \(X+H\) with \(H\) cubic homogeneous, while increasing dimension; Drużkowski further reduces stably to cubic-linear maps. These are powerful for all-dimensional JC but do **not** say that JC2 reduces to cubic maps of \(\mathbb A^2\). Moreover degree-2 maps are known invertible (S. Wang, *J. Pure Appl. Algebra* 10 (1980), 103–105). Hence neither reduction supplies a low-degree plane search space.

3. **Plane degree exclusions.** Moh proved JC2 for \(\max(\deg P,\deg Q)\le100\): *J. Reine Angew. Math.* 340 (1983), 140–212, “On the Jacobian conjecture and the configurations of roots.” Magnus’ theorem gives invertibility when \(\gcd(\deg P,\deg Q)=1\); later refinements cover broader degree pairs. For a conservative search, impose both degrees \(>100\) and \(\gcd(\deg P,\deg Q)>1\). (These are exclusion filters, not normal forms.)

4. **One-line criterion.** If a Keller map \(F:\mathbb C^2\to\mathbb C^2\) restricts injectively to some affine line, then it is invertible (Gwoździewicz, *Bull. Soc. Sci. Lett. Łódź* 39 (1989), 1–10, “Injectivity on one line”). A candidate counterexample must therefore fail injectivity on every affine line.

5. **JC/DC direction.** Belov-Kanel–Kontsevich, arXiv:[math/0512171](https://arxiv.org/abs/math/0512171), prove stable equivalence and, at fixed indices, \(JC_{2n}\Rightarrow DC_n\). Therefore \(JC_2\Rightarrow DC_1\), and a counterexample to \(DC_1\) disproves JC2. Proving DC1—or a finite-band fragment—does not prove JC2.

### Deduction

The useful counterexample target is a nonproper, generically finite plane Keller map of degree \(>1\), beyond the known degree filters, with no injective line. Stable cubic reductions are better viewed as certification tools after stabilization, not as a direct plane ansatz.

## 2. Why the verified \(\mathbb A^3\) telescope does not descend automatically

### Verified facts from the supplied context

The three-dimensional family is noninjective Keller, has every generic degree \(\ge3\), and exhibits explicit \(S_4/S_5\) monodromy.

### Deduction

A coordinate slice \(z=c\) yields \((F_1,F_2)|_{z=c}\), whose Jacobian is the **minor** \(\partial(F_1,F_2)/\partial(x,y)\), not \(\det DF\); constant \(\det DF\) does not make that minor constant. Projection also identifies points but generally does not preserve étaleness. Eliminating \(z\) usually gives an algebraic correspondence or a hypersurface relation, not two polynomials in \(x,y\). Finally, a polynomial invariant quotient by a fiber symmetry can ramify; nontrivial \(S_4/S_5\) monodromy is not itself a free polynomial quotient of \(\mathbb A^3\). Thus the third variable may be exactly the state variable that makes telescoping and determinant cancellation coexist.

### Three precise descent/degeneration tests

1. **Invariant graph test.** For \(F=(F_1,F_2,F_3)\), solve for polynomials \(h(x,y),A(u,v),B(u,v)\) such that
   \[
   F_3(x,y,h)=h(F_1(x,y,h),F_2(x,y,h)),\quad
   (A,B)=(F_1,F_2)|_{z=h},\quad J(A,B)=1.
   \]
   Noninjectivity descends only if a known colliding pair lies on this graph. Gröbner reduction modulo degree bounds makes this finite and falsifiable.

2. **Affine foliation/quotient test.** Seek a polynomial submersion \(\pi:\mathbb A^3\to\mathbb A^2\) and \(G:\mathbb A^2\to\mathbb A^2\) satisfying \(\pi\circ F=G\circ\pi\), with \(J(G)=1\). Start with \(\pi=(\ell_1,\ell_2)\) affine-linear; coefficient comparison decides existence exactly. Require two known colliding points to have distinct or equal projected behavior appropriate to proving noninjectivity of \(G\), rather than merely of \(F\).

3. **Flat specialization test.** Introduce \(t\) by weighted scaling/conjugation and seek \(F_t\) Keller for \(t\ne0\), with a regular limit \(F_0=(G_1(x,y),G_2(x,y),z)\). Verify coefficientwise \(\det DF_t=1\), flatness of the collision/elimination ideal over \(\mathbb C[t]\), and persistence at \(t=0\) of two distinct points in one fiber. Without the last two checks, generic degree can drop and collisions can escape to infinity.

## 3. The finite-band \(A_1\) program, only as a counterexample route

Let \(A_1=\mathbb C\langle x,\partial\rangle/([\partial,x]-1)\), \(E=x\partial\). Since \([E,x^i\partial^j]=(i-j)x^i\partial^j\), a finite \(\mathrm{ad}(E)\)-band is a finite weight support.

### Deduction

Band rigidity is useful here only if its failure produces an explicit non-surjective endomorphism \(x\mapsto P,\partial\mapsto Q\) with \([Q,P]=1\). Positive rigidity merely eliminates a DC1 counterexample class and says nothing affirmative about JC2.

The smallest genuinely band-2 Laurent-normalized system is
\[
P=x+a(E)+b(E)\partial,
\qquad Q=\partial+c(E)x^{-1}+d(E)x^{-2},
\]
computed in the localization \(A_1[x^{-1}]\), using
\(x f(E)=f(E-1)x\), \(\partial f(E)=f(E+1)\partial\).
Clear powers of \(x\), impose \([Q,P]=1\), and equate each \(x^k\)-coefficient; additionally impose divisibility conditions making \(P,Q\in A_1\): \(c(E)\) divisible by \(E\), and \(d(E)\) by \(E(E-1)\). Bound \(\deg a,b,c,d\le2\) first. This is the smallest system containing both a second negative weight and nonlinear \(E\)-dependence; lower choices collapse toward band 1 or affine/triangular automorphisms.

### Speculation

A non-automorphic solution may first appear only at higher coefficient degree. Any solution must be tested for surjectivity; satisfying \([Q,P]=1\) constructs an endomorphism, not automatically a DC1 counterexample.

## 4. First milestone (days, exact symbolic computation)

### Selected milestone: affine-quotient obstruction for the explicit 3D family

**Inputs:** the exact formulas for each \(F\) in the verified \(S_4/S_5\) examples and exact colliding pairs.

**Unknowns:** two independent affine forms \(\ell_i=r_i x+s_i y+t_i z+u_i\), and plane polynomials \(G_1,G_2\) of degree \(\le\deg F\).

**Equations:** expand
\[
\ell_i(F(x,y,z))=G_i(\ell_1(x,y,z),\ell_2(x,y,z))\quad(i=1,2),
\]
coefficientwise; saturate by a nonzero \(2\times2\) minor of the \(2\times3\) linear-part matrix; then impose \(J(G)=1\). Solve over \(\mathbb Q\) by Gröbner basis, separately on minor charts.

**Falsifiable success criterion:** either (A) exhibit exact \((\pi,G)\) with \(J(G)=1\) and an exact collision for \(G\), which would be a JC2 counterexample requiring independent verification; or (B) certify the saturated ideals are unit ideals on every chart, proving no affine two-dimensional semiconjugate quotient exists for those examples at the stated degree bound. Outcome B is the realistic, publishable obstruction and directs effort to nonlinear graphs or flat degenerations.

No outcome here is presently established.