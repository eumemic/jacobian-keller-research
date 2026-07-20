> **Status: Provisional research archive.**
> Preserved for its derivations and research history. Computational checks pass in the archived environment, but the mathematical arguments and literature claims have not been independently reviewed. This is not a claim of a proof or counterexample to an open conjecture.

# Prong two, opened: the Dixmier route to the two-dimensional conjecture

**Objective.** Find an injective, non-surjective endomorphism of the first Weyl
algebra A1 = C⟨x, ∂⟩/([∂,x]=1), i.e. disprove DC1. Since JC2 ⇒ DC1
(Tsuchimoto 2005; Belov-Kanel–Kontsevich 2007), a DC1 counterexample disproves
the two-dimensional Jacobian Conjecture. Mind the direction: *proving* DC1 would
close this route but would **not** prove JC2.

An endomorphism of A1 is a pair (X, D) with [D, X] = 1; it is automatically
injective (A1 is simple), so non-surjectivity is the whole game.

## Result P1 (proved): graded rigidity of A1 — the quantum Euler obstruction

Let E = x∂ and consider endomorphisms homogeneous for the weight grading
(x, ∂) ↦ (λx, λ⁻¹∂): necessarily φ(x) = x·g(E), φ(∂) = h(E)·∂ with g, h ∈ C[E].
Using f(E)x = x f(E+1) and f(E)∂ = ∂ f(E−1), the relation [φ(∂), φ(x)] = 1
becomes a difference equation for P(E) := g(E)h(E):

    (E+1)·P(E) − E·P(E−1) = 1.

If deg P = k ≥ 1 with leading coefficient a, the left side has leading term
a(k+1)·E^k ≠ 0, so P must be constant, and the equation forces P ≡ 1. Hence the
only graded endomorphisms are the diagonal automorphisms x ↦ cx, ∂ ↦ c⁻¹∂.

This is the exact quantization of Theorem 4's Euler-positivity mechanism
(the classical coefficient 1 + a·deg h + b·deg g becomes k+1). Consequence: a
DC1 counterexample, like a JC2 counterexample, cannot be graded — any quantized
telescoping tower must be built from *filtered-but-not-graded* pieces, with the
interesting cancellation happening below the top order.

## Result P2 (proved modulo write-up): tame obstruction of mod-p shadows

The squarefree-pullback lemma (Lemma 5.1 of the note) is characteristic-free
over perfect fields, and the Kummer step needs only cyclic degree prime to p.
So over F̄p, the cyclic obstruction holds for all **tame** cyclic subextensions;
the loophole is exactly wild Artin–Schreier degree-p covers — which is where the
classical char-p pathologies (t ↦ t − t^p, derivative 1, non-injective) live.
Consistency check passed: the theorem and the known counterexamples partition
cleanly along tame/wild.

Consequence for DC1: a char-0 counterexample φ reduces, for almost all p, to an
endomorphism of A1 ⊗ Fp, which induces a Jacobian-1 endomorphism ψ_p of the
center ≅ A²_{Fp} (this is the Tsuchimoto/BKK bridge). The geometric degree d of
ψ_p is bounded independently of p, so for p > d every subextension is tame and
the full cyclic obstruction applies to ψ_p. In particular ψ_p cannot have
degree 2 for p > 2, and its monodromy pair (G, H) must satisfy H ↠ G^ab.
Any construction attempt can be screened against this for free.

## Attack plan

1. **Formalize P1, P2** (short; P2's statement already appears as Remark 5.7 in
   the note, flagged "pursued elsewhere").
2. **Localization telescoping.** Quantize the tower's middle layer: the classical
   monomial substitution becomes conjugation by quantum-torus elements in the
   localization A1[x⁻¹]. Compose: (shear in A1[x⁻¹]) ∘ (quantum-monomial
   conjugation) ∘ (shear), engineered so denominators telescope away, mirroring
   how g² cleared in the classical h1 = S/g². First milestone: exhibit *any*
   φ ∈ End(A1[x⁻¹]) whose associated-graded shadow is the classical degree-3
   base map, then measure exactly what obstructs pulling φ into A1. The
   obstruction, if structural, is itself a theorem.
3. **Calibration search.** Exhaustive solve of [D, X] = 1 in low Bernstein
   filtration degree (≤ 4, undetermined coefficients), modulo tame automorphisms:
   confirms the workspace and gives either evidence for DC1 in low degree or a
   miracle. (Expected: all solutions tame; the value is the machinery and the
   precise shape of the variety of solutions.)
4. **Decision point.** If a structural obstruction to step 2 emerges, invert the
   effort: try to prove DC1's graded/filtered rigidity in general — which would
   kill this route to ¬JC2 and redirect everything into the prong-one
   (leading-form degeneration + Theorem 4) program, whose known gap is the
   rank-drop of limits.

This archived memo ends here; step 2 was identified as the next concrete computation.
