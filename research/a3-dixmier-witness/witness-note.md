> **Status: certified computational note with self-contained proofs.**
> This note certifies and promotes the provisional archive note
> `archive-import/provisional/weyl-a3/endomorphism-formulas.md`. Every statement
> labeled **[verified]** is checked by exact rational arithmetic in the
> accompanying `verify_a3_witness.py` (pure Python, `fractions.Fraction`, no
> dependencies); during certification all polynomial-level identities (Jacobian
> data, archived entries, bracket coefficients, fiber, family membership) were
> additionally reproduced with a second, independent implementation (SymPy exact
> symbolic arithmetic), and the two agree. The Weyl-operator-level identities
> are computed directly in the pure-Python Weyl implementation and are, in
> addition, equivalent to polynomial-level identities by the displayed
> computation in Section 3. Every statement labeled **[proved]** has a complete proof in
> this note. Statements labeled **[cited]** rest on the named reference and are
> not re-proved here. No novelty or priority is claimed anywhere in this note.

# An explicit injective, non-surjective endomorphism of the Weyl algebra A_3

## Status summary

- **Verified (exact arithmetic):** all formulas of the archived note reproduce
  from `F` alone with **zero discrepancies** — `det JF = -2`, the nine
  polynomial entries of `JF^{-1} = -adj(JF)/2`, the three `phi(x_i)`, the full
  operator relations `[phi(d_i), phi(x_j)] = delta_ij`,
  `[phi(d_i), phi(d_j)] = 0`, `[phi(x_i), phi(x_j)] = 0` (in an independently
  implemented normal-ordered Weyl algebra), the equivalent vector-field bracket
  identities, the three-point fiber of `F` over `(-1/4, 0, 0)`, and membership
  of `F` in the in-repo paper's telescoping family.
- **Proved (complete proofs below):** `A_3` is simple, so `phi` is injective
  (Section 5); if `phi_F` is an automorphism of `A_n` then `F` is an invertible
  polynomial map (Section 6, the classical "DC_n implies JC_n" transfer, with a
  full self-contained proof via exact order preservation); hence `phi` is not
  surjective (Section 7), and the coordinate generators `x_1, x_2, x_3` are
  explicit elements outside its image.
- **Cited (not re-proved here):** geometric degree of `F` equals 3 (in-repo
  `paper/main.tex`, Theorem "Telescoping family", applied with the verified
  member data); the literature provenance of the transfer statement (van den
  Essen, Belov-Kanel--Kontsevich, Tsuchimoto, Adjamagbo--van den Essen). The
  conclusion of this note does **not** depend on the cited geometric-degree
  theorem: non-injectivity of `F` is verified directly from the explicit fiber.
- **Not claimed:** novelty, priority, literature positioning, or peer-reviewed
  status. The interpretation — a counterexample to the Dixmier conjecture for
  `A_3` — is exactly the standard transfer of the non-injective Keller map `F`
  and stands or falls with the machine-checked computations and the proofs
  below; both are laid out so they can be independently reviewed.

## 1. Frozen conventions

Work over `C`; every displayed constant is rational, and everything below is
valid over any field `K` of characteristic zero. The Weyl algebra `A_n` is the
associative unital `C`-algebra with generators `x_1..x_n, d_1..d_n` and
relations

```
[d_i, x_j] = delta_ij,   [x_i, x_j] = 0,   [d_i, d_j] = 0.
```

Orientation is fixed as `[d_i, x_j] = d_i x_j - x_j d_i = delta_ij` (the same
convention as the archived note; it agrees with Belov-Kanel--Kontsevich, whose
`omega_{i,n+j} = -delta_ij` says `[x_i, d_j] = -delta_ij`). We write `n = 3`
and `x = x_1, y = x_2, z = x_3` when convenient.

**Normal form.** Every element of `A_n` is uniquely `sum_alpha p_alpha(x)
d^alpha` with `d^alpha = d_1^{a_1}...d_n^{a_n}`. Spanning follows from the
rewriting rule below; uniqueness follows from faithfulness of the action on
`C[x]` (with `x_i` acting by multiplication and `d_i` by `partial_i`): if
`P = sum p_alpha partial^alpha` has nonempty support, pick `alpha_0` minimal in
the support for the componentwise partial order; then `P(x^{alpha_0}) =
alpha_0! p_{alpha_0} != 0` in characteristic zero.

**Products.** The normal-ordering (Leibniz) rule, which is also exactly what
`verify_a3_witness.py` implements, is

```
(p d^alpha)(q d^beta) = sum_{gamma <= alpha} binom(alpha, gamma) p (partial^gamma q) d^{alpha + beta - gamma}.
```

**Order filtration.** Let `F_r` be the span of the `p(x) d^alpha` with
`|alpha| <= r` ("order at most `r`"). The rule above gives `F_r F_s <=
F_{r+s}`, and modulo `F_{r+s-1}` only the `gamma = 0` term survives:
`(p d^alpha)(q d^beta) == p q d^{alpha+beta}`. For `a` in `F_r \ F_{r-1}`
define the symbol `sigma_r(a) = sum_{|alpha| = r} p_alpha xi^alpha`, a
`xi`-homogeneous element of the commutative polynomial ring `C[x][xi_1..xi_n]`.
Then `sigma` is multiplicative, `C[x][xi]` is a domain, hence

```
ord(ab) = ord(a) + ord(b),   sigma(ab) = sigma(a) sigma(b),
```

and `A_n` is a domain. Note `F_0 = C[x]`: the order-zero elements are exactly
the polynomials in `x`.

## 2. The Keller map F and its verified properties

The polynomial map `F = (F_1, F_2, F_3) : C^3 -> C^3` is

```
F_1 = (1+xy)^3 z + y^2 (1+xy)(4+3xy)
F_2 = y + 3x (1+xy)^2 z + 3x y^2 (4+3xy)
F_3 = 2x - 3x^2 y - x^3 z
```

(total degree 7; this is the map the archived note calls "the announced
degree-7 Keller counterexample F").

- **[verified]** `det JF = -2` identically, where `JF_{ij} = partial F_i /
  partial x_j`. Hence `F` is a Keller map and `M := JF^{-1} = -adj(JF)/2` is a
  matrix of polynomials (with rational coefficients); **[verified]**
  `JF M = M JF = I` and `det M = -1/2`.
- **[verified]** The fiber of `F` over `(-1/4, 0, 0)` contains the three
  distinct rational points

  ```
  (0, 0, -1/4),   (1, -3/2, 13/2),   (-1, 3/2, 13/2),
  ```

  each mapped exactly to `(-1/4, 0, 0)`. Consequently **`F` is not injective**,
  so `F` admits no inverse — not a polynomial inverse, and not even a
  set-theoretic one. This is the only property of `F` beyond `det JF = -2`
  that the main conclusion uses.
- **[verified]** The three fiber points have pairwise distinct first
  coordinates and pairwise distinct second coordinates, and two distinct third
  coordinates. (Used in Section 7 to exhibit explicit non-image elements.)
- **[verified]** `F` is the member `(E, h) = (2 - 3 nu, 2 - 3 t)` of the
  telescoping family of the in-repo paper (`paper/main.tex`, Theorem
  "Telescoping family"): with `s = x^2 z`, `t = xy`, `u = 1+t`,
  `g = h(t) - s = 2 - 3t - s`, `A(nu) = nu^2 - nu^3`, `C(nu) = 4 nu - 3 nu^2`,
  the exact polynomial identities

  ```
  A(ug) + u g^2 = x^2 g^2 F_1,   C(ug) + 2g = x g F_2,   F_3 = x g
  ```

  hold, and the three scalar polynomiality conditions hold at `a = h(0) = 2`:
  `A(2) = -4 = -a^2`, `C(2) = -4 = -2a`, `1 + E(2) + h'(0)(E(2)+2)/a = 1 - 4 +
  (-3)(-2)/2 = 0`.
- **[cited]** By that theorem (degree of `E` is `n = 1`), the geometric degree
  of `F` is `n + 2 = 3`: the generic fiber has exactly three points, and the
  verified triple fiber above is an explicit rational instance. This citation
  is context only; nothing below depends on it.

## 3. The endomorphism phi

Define, following the archived note,

```
phi(x_i) = F_i,          phi(d_j) = Y_j := sum_k M_{kj} d_k,     M = JF^{-1},
```

so the coefficient of `d_k` in `phi(d_j)` is the `(k, j)` entry of `JF^{-1}`.

- **[verified]** All nine archived coefficient polynomials (the displayed
  `phi(d_1), phi(d_2), phi(d_3)` of the archive note) equal the recomputed
  `(JF^{-1})_{kj}`, entry by entry, and the three archived `phi(x_i)` equal
  `F_i`. **No discrepancy of any kind was found in the archived formulas.**
- **[verified]** In the independently implemented Weyl algebra of Section 1,
  the fifteen defining relations map to relations:

  ```
  [phi(d_i), phi(x_j)] = delta_ij   (9 identities)
  [phi(d_i), phi(d_j)] = 0          (3 nontrivial pairs, computed directly in A_3)
  [phi(x_i), phi(x_j)] = 0          (3 pairs)
  ```

**Well-definedness [proved].** `A_3` is the quotient of the free algebra on
the six generators by the two-sided ideal generated by the fifteen relators.
Any assignment of generator images extends uniquely to a homomorphism from the
free algebra, and it factors through `A_3` precisely when every relator maps to
zero — which is the verified list above. Hence `phi` is a unital `C`-algebra
endomorphism of `A_3` with `phi(1) = 1`, and `phi(g(x_1,x_2,x_3)) =
g(F_1,F_2,F_3)` for every polynomial `g`.

**Equivalence with vector-field commutativity [proved and verified both
ways].** For first-order operators `a = sum_k A_k d_k`, `b = sum_k B_k d_k`
with polynomial coefficients, the normal-ordering rule gives

```
ab = sum_{k,l} A_k B_l d_k d_l + sum_l ( sum_k A_k partial_k B_l ) d_l,
```

and the second-order double sums of `ab` and `ba` coincide (swap the summation
indices and use `d_k d_l = d_l d_k`). Hence

```
[a, b] = sum_l ( sum_k ( A_k partial_k B_l - B_k partial_k A_l ) ) d_l,
```

which is exactly the Lie bracket of the vector fields `sum_k A_k partial_k`
and `sum_k B_k partial_k`. So `[phi(d_i), phi(d_j)] = 0` in `A_3` is
equivalent to commutativity of the derivations `Y_j = sum_k M_{kj} partial_k`
of `C[x]`. **[verified]** both ways: the twenty-seven bracket coefficient
identities `sum_k ( M_{ki} partial_k M_{lj} - M_{kj} partial_k M_{li} ) = 0`
directly, and the three operator commutators inside the Weyl implementation.

**Remark 3.1 (general Keller maps; not needed for the sequel).** For every
Keller map `F` the fields `Y_j` commute, so `phi_F` is always an endomorphism.
Proof: by Lemma 6.1 below, `C(x)` is an algebraic — hence, in characteristic
zero, separable — extension of the subfield `C(F) = C(F_1,...,F_n)`. First,
`Y_j(F_i) = sum_k (JF)_{ik} M_{kj} = delta_ij`, so by the chain rule `Y_j`
maps `C(F)` into `C(F)` and restricts on it to the `j`-th coordinate partial
derivative. Therefore `[Y_i, Y_j]` is a derivation of `C(x)` vanishing on
`C(F)` (equality of mixed partials). A derivation `D` of `C(x)` vanishing on
`C(F)` vanishes: for `v` in `C(x)` with minimal polynomial `m` over `C(F)`,
`0 = D(m(v)) = m'(v) D(v)` and `m'(v) != 0` by separability. Hence
`[Y_i, Y_j] = 0`. (For the specific `F` of this note the commutativity is
machine-verified, so the sequel does not rely on this remark.)

## 4. Statement being certified

> **Theorem.** `phi` is an injective, non-surjective `C`-algebra endomorphism
> of `A_3`. Explicitly, none of `x_1, x_2, x_3` lies in the image of `phi`.

Injectivity is Section 5; non-surjectivity is Sections 6–7. Interpretation is
discussed in Section 8.

## 5. Injectivity: A_n is simple [proved]

**Proposition 5.1.** Over a field of characteristic zero, `A_n` is simple.

*Proof.* Let `I != 0` be a two-sided ideal and `0 != a = sum p_alpha d^alpha`
in `I`. From the normal-ordering rule, `[p d^alpha, x_i] = alpha_i p
d^{alpha - e_i}`. Choose `beta` with `p_beta != 0` and `|beta|` maximal, and
apply `(ad x_1)^{beta_1} ... (ad x_n)^{beta_n}` to `a` (each `ad x_i = [.,
x_i]` preserves `I`). A monomial `p_alpha d^alpha` survives only if `alpha >=
beta` componentwise; together with `|alpha| <= |beta|` this forces `alpha =
beta`. The result is `beta! p_beta`, a nonzero element of `I ∩ C[x]`
(characteristic zero: `beta! != 0`). Now for `q` in `C[x]`, `[d_i, q] =
partial_i q`; choosing `gamma` of maximal total degree in the support of
`p_beta` and applying `(ad d_1)^{gamma_1} ... (ad d_n)^{gamma_n}` — with
`ad d_i = [d_i, .]` — yields `gamma!` times the corresponding nonzero
coefficient, a nonzero scalar in `I`. Hence `1` is in `I` and `I = A_n`. ∎

**Corollary 5.2.** `phi` is injective: its kernel is a two-sided ideal of the
simple algebra `A_3` not containing `1` (as `phi(1) = 1 != 0`), hence zero. ∎

## 6. If phi_F is an automorphism, then F is invertible [proved]

This section makes precise and re-proves, self-contained, the classical easy
implication behind "DC_n implies JC_n". Provenance **[cited]**: Belov-Kanel
and Kontsevich, arXiv:math/0512171 (Mosc. Math. J. 7 (2007) 209–218), state on
p. 2: *"It is well-known that DC_n implies JC_n ... The endomorphism
phi*_diff of the Weyl algebra preserves the degree of differential operators.
Restricting phi*_diff to zero order differential operators, we obtain the
usual pullback phi* of functions"* — citing van den Essen, *Polynomial
Automorphisms and the Jacobian Conjecture*, Progr. Math. 190, Birkhäuser 2000
(their Remark 1 points at Theorem 10.4.2 there) and Bass–Connell–Wright. The
proof below is the complete version of exactly that argument, in this note's
conventions, so the note stands alone.

Throughout this section `F` is any Keller map of `C^n` (`det JF` a nonzero
constant), `M = JF^{-1} = adj(JF)/det JF` its polynomial inverse-Jacobian
matrix, and `phi = phi_F` is defined on generators as in Section 3 and is
assumed to be an endomorphism of `A_n` (machine-verified for the `F` of this
note; true in general by Remark 3.1).

**Lemma 6.1 (Jacobian criterion / dominance).** If `det JF != 0` then
`F_1, ..., F_n` are algebraically independent over `C`; equivalently, the
substitution homomorphism `tau : C[w_1..w_n] -> C[x]`, `g -> g(F)`, is
injective.

*Proof.* Suppose `P(F) = 0` with `P != 0` of minimal total degree. Applying
`partial/partial x_k` and the chain rule: `sum_i (partial_i P)(F) (JF)_{ik} =
0` for every `k`. Since `det JF != 0`, the row vector `((partial_i P)(F))_i`
vanishes: `(partial_i P)(F) = 0` for all `i`. Each `partial_i P` has smaller
degree, so minimality gives `partial_i P = 0` for all `i`; in characteristic
zero `P` is then constant, and `P(F) = 0` forces `P = 0` — contradiction. ∎

**Lemma 6.2 (exact order preservation).** For every `a != 0` in `A_n`,
`ord(phi(a)) = ord(a)`. More precisely, if `ord(a) = r` and `sigma_r(a) =
sum_{|alpha| = r} p_alpha xi^alpha`, then

```
sigma_r(phi(a)) = sum_{|alpha| = r} p_alpha(F) eta^alpha,    eta_j := sum_k M_{kj} xi_k,
```

which is nonzero.

*Proof.* Write `a = sum_{|alpha| <= r} p_alpha d^alpha`. Since `phi` is a
homomorphism, `phi(a) = sum_alpha p_alpha(F) Y^alpha` with `Y^alpha =
Y_1^{alpha_1} ... Y_n^{alpha_n}`. Each `Y_j` lies in `F_1` with `sigma_1(Y_j)
= eta_j != 0` (the `j`-th column of the invertible matrix `M` is nonzero), so
by multiplicativity of symbols in the domain `C[x][xi]`, `ord(Y^alpha) =
|alpha|` and `sigma_{|alpha|}(Y^alpha) = eta^alpha`. Hence `phi(a)` lies in
`F_r`, and modulo `F_{r-1}` only the `|alpha| = r` terms contribute, giving
the displayed formula `Q(eta)` where `Q(xi) := sum_{|alpha| = r} p_alpha(F)
xi^alpha` in `C[x][xi]`.

It remains to show `Q(eta) != 0`. First `Q != 0`: the monomials `xi^alpha` are
distinct, some `p_alpha != 0`, and then `p_alpha(F) != 0` by Lemma 6.1.
Second, `xi_j -> eta_j = sum_k M_{kj} xi_k` is the linear substitution by the
matrix `M^T`, which is invertible over `C(x)`; it therefore extends to an
automorphism of the polynomial algebra `C(x)[xi_1..xi_n]`, and automorphisms
are injective, so `Q != 0` implies `Q(eta) != 0`. Thus `sigma_r(phi(a)) != 0`
and `ord(phi(a)) = r`. ∎

**Proposition 6.3 (the transfer).** If `phi_F` is surjective — equivalently,
an automorphism, since injectivity is automatic by Corollary 5.2 — then there
is a polynomial map `G : C^n -> C^n` with `G ∘ F = F ∘ G = id`. In
particular `F` is bijective.

*Proof.* Let `psi = phi^{-1}`. For every `a != 0`, Lemma 6.2 applied to the
element `psi(a)` gives `ord(a) = ord(phi(psi(a))) = ord(psi(a))`: **the
inverse also preserves order exactly.** Since `ord(x_i) = 0` and `F_0 = C[x]`
(Section 1), the elements

```
G_i := psi(x_i)
```

are honest polynomials in `x_1, ..., x_n`, defining a polynomial map `G =
(G_1, ..., G_n)`. Applying `phi` and Section 3's substitution formula:

```
x_i = phi(psi(x_i)) = phi(G_i) = G_i(F_1, ..., F_n),
```

i.e. `G ∘ F = id` as an identity of polynomial maps. For the other
composition put `H := F ∘ G`. Then `H_i(F(x)) = F_i(G(F(x))) = F_i(x)` using
`G(F(x)) = x`, so `tau(H_i(w) - w_i) = 0` with `tau` the injective
substitution homomorphism of Lemma 6.1; hence `H_i(w) = w_i`, i.e. `F ∘ G =
id`. ∎

**Corollary 6.4 (contrapositive, as used here).** If `F` is not invertible —
for instance, not injective as a map `C^n -> C^n` — then `phi_F` is not
surjective. ∎

**Remark 6.5 (conventions audit).** Two places where orientation and
direction matter, both fixed above: (i) with `[d_i, x_j] = +delta_ij` and
`phi(d_j) = sum_k (JF^{-1})_{kj} d_k` (columns of `JF^{-1}`), the relation
check is `[phi(d_i), phi(x_j)] = Y_i(F_j) = (JF · JF^{-1})_{ji} = delta_ij`,
which is the **[verified]** list of Section 3; (ii) the automorphism
hypothesis produces a **left** inverse `G ∘ F = id` first (because `psi` is
applied to the coordinates and `phi` substitutes `F`), and the two-sided
inverse then follows from Lemma 6.1, not from any density or topology.

## 7. Non-surjectivity of phi and explicit non-image elements

**Theorem 7.1.** `phi` is injective and not surjective; in particular it is a
`C`-algebra endomorphism of `A_3` that is not an automorphism.

*Proof.* Injectivity is Corollary 5.2. If `phi` were surjective, Proposition
6.3 would make `F` invertible, in particular injective; but `F` identifies the
three distinct verified points of Section 2 over `(-1/4, 0, 0)`. ∎

**Proposition 7.2 (explicit witnesses).** None of `x_1, x_2, x_3` lies in the
image of `phi`.

*Proof.* Suppose `x_i = phi(P)`. By Lemma 6.2, `ord(P) = ord(x_i) = 0`, so
`P` is a polynomial and `x_i = P(F_1, F_2, F_3)`. Evaluating at the three
verified fiber points `p_1, p_2, p_3` of Section 2, whose `F`-images coincide,
the function `x_i` would take a single value on `{p_1, p_2, p_3}`. But the
`i`-th coordinates there are `(0, 1, -1)` for `i = 1`, `(0, -3/2, 3/2)` for
`i = 2`, and `(-1/4, 13/2, 13/2)` for `i = 3` — in every case not all equal.
∎

**Remark 7.3 (higher n).** Padding `F` by identity coordinates,
`F' = (F_1, F_2, F_3, x_4, ..., x_n)`, is again a non-injective Keller map for
every `n >= 3`, and the same construction and proofs produce an injective
non-surjective endomorphism of `A_n` for every `n >= 3`. This routine
extension is recorded for completeness; it is not separately machine-checked.
Nothing in this note bears on `A_1` or `A_2`.

## 8. Interpretation, stated conservatively

The Dixmier conjecture `DC_n` (Dixmier 1968, Problème 1 for `n = 1`; see the
references) asserts that every endomorphism of `A_n` over a characteristic-zero
field is an automorphism. Subject to (i) the exact computations reproduced by
`verify_a3_witness.py` (all of which also reproduced under an independent
SymPy implementation during certification) and (ii) the complete proofs in
Sections 5–6 above, Theorem 7.1 exhibits an explicit counterexample to `DC_3`.

Three cautions, in the repository's usual register:

1. **This is a transfer, not an independent phenomenon.** The witness is the
   standard Weyl-algebra image of the non-injectivity of the Keller map `F` —
   the map announced publicly as a Jacobian-conjecture counterexample (see the
   in-repo `paper/main.tex` and its reference [Gallagher26]; that announcement
   is recent and, to our knowledge, not yet peer-reviewed). What this note
   adds is certification: the two facts about `F` that the argument consumes —
   `det JF = -2` and the explicit triple fiber — are themselves machine-verified
   here by exact arithmetic from the displayed formulas, so the `A_3` statement
   does not inherit any unverified claim about `F`.
2. **What replaces the archived note's flagged gap.** The archive flagged the
   non-surjectivity appeal to "the standard easy implication DC_3 => JC_3" as
   needing independent review. Section 6 supplies the precise statement and a
   complete proof (exact order preservation of the filtration, order
   preservation of the inverse, restriction to order zero, and the
   left-to-two-sided inverse upgrade via algebraic independence), in the
   orientation and composition conventions of this repository. The archived
   note's mathematical content survives review unchanged; its formulas contain
   no errors.
3. **Scope.** No claim is made about `DC_1`, `DC_2`, novelty, or priority.
   The literature transfer statement itself is classical (van den Essen Thm
   10.4.2; Belov-Kanel--Kontsevich; Tsuchimoto; Adjamagbo--van den Essen) and
   is re-proved here only so that this note is self-contained.

## 9. Reproducibility

Run:

```
python3 verify_a3_witness.py
```

The checker is pure Python 3 (standard library only; exact
`fractions.Fraction` arithmetic; an independent normal-ordered Weyl-algebra
implementation with Leibniz-rule multiplication). It recomputes everything
from the three displayed polynomials `F_1, F_2, F_3` alone and compares
against the archived formulas, which are hardcoded verbatim from the archive
note. A successful run ends with the line `ALL A3 WITNESS CHECKS PASSED`; at
the time of writing, all 19 checks pass. No floating-point or numeric-only
step occurs anywhere.

## References

- J. Dixmier, *Sur les algèbres de Weyl*, Bull. Soc. Math. France 96 (1968),
  209–242.
- A. van den Essen, *Polynomial Automorphisms and the Jacobian Conjecture*,
  Progress in Mathematics 190, Birkhäuser, Basel, 2000. (Theorem 10.4.2 for
  the filtration-preserving invertibility statement.)
- A. Belov-Kanel, M. Kontsevich, *The Jacobian conjecture is stably equivalent
  to the Dixmier conjecture*, Mosc. Math. J. 7 (2007), no. 2, 209–218;
  arXiv:math/0512171. (P. 2: the "well-known" implication `DC_n => JC_n` and
  its degree-filtration sketch; Remark 1 cites van den Essen Thm 10.4.2.)
- Y. Tsuchimoto, *Endomorphisms of Weyl algebra and p-curvatures*, Osaka J.
  Math. 42 (2005), no. 2, 435–452.
- P. K. Adjamagbo, A. van den Essen, *A proof of the equivalence of the
  Dixmier, Jacobian and Poisson conjectures*, Acta Math. Vietnam. 32 (2007),
  no. 2–3; arXiv:math/0608009.
- H. Bass, E. H. Connell, D. Wright, *The Jacobian conjecture: reduction of
  degree and formal expansion of the inverse*, Bull. Amer. Math. Soc. (N.S.) 7
  (1982), 287–330.
- In-repo: `paper/main.tex` (Theorem "Telescoping family": polynomiality
  criterion, `det = -2`, geometric degree `n+2`; the `(E,h) = (2-3nu, 2-3t)`
  member is this note's `F`), and
  `archive-import/provisional/weyl-a3/endomorphism-formulas.md` (the certified
  source note).
