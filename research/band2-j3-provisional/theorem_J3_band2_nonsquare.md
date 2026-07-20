# Provisional J3/M4 argument: the two-sided non-square sector of band 2

> **PROVISIONAL — INDEPENDENT AUDIT PENDING — NOT PEER REVIEWED.** The argument below is preserved as submitted for audit. It is not presented as an established theorem or claimed result. The accompanying computation checks selected symbolic identities only; it is computational support, not a proof, proof review, or independent validation.

The accompanying `verify_J3.py` script reports 13 exact symbolic checks.
Throughout, "shifted square" means: classically, c·h(τ)² with c ∈ C*, h ∈ C[τ]
(constants included, h = 1); quantumly, c·h(E)h(E+1). "Non-square" means not of
this form — in particular nonconstant.

## Statement

**Theorem J3 (classical).** Let X = Σ_{k=−2}^{2} x^k a_k(τ), D = Σ x^l b_l(τ)
(τ = xξ) satisfy {D, X} ∈ C*, with a₂ ≠ 0 and a₋₂ ≠ 0. Then at least one of
a₂, a₋₂ is a square (c·h²). Equivalently: there is **no** Keller pair — indeed
no Laurent pair in C[x^{±1}, ξ] — of band width ≤ 2 on both sides whose two
extreme coefficients are both non-squares.

**Theorem J3q (quantum).** Same statement in A₁[x⁻¹] = ⊕ x^k C[E] for pairs
with [D, X] = 1: at least one of a₂(E), a₋₂(E) is a shifted square c·h(E)h(E+1).

The theorem needs no membership/divisibility hypothesis: the emptiness already
holds in the localized algebras, so — unlike band 1 — the non-square sector has
no polar solutions either.

## Proof (classical)

Write the bracket components: for each m, Σ_{k+l=m}(k a_k b_l′ − l a_k′ b_l)
equals δ_{m,0}·(const). The cascade:

**m = ±4.** 2(a₂b₂′ − a₂′b₂) = 0 and the mirror: Wronskians vanish, so
b₂ = λ₂a₂, b₋₂ = μ₂a₋₂ with λ₂, μ₂ ∈ C (a₂, a₋₂ ≠ 0).

**m = ±3.** 2a₂b₁′ − a₂′b₁ = λ₂(2a₁′a₂ − a₁a₂′), for which b₁ = λ₂a₁ is a
particular solution; the homogeneous solutions satisfy h² = c·a₂ (Lemma J2), so
**non-square a₂ forces b₁ = λ₂a₁ exactly** — with the *same* constant λ₂.
Mirror: b₋₁ = μ₂a₋₁.

**m = ±2.** The (1,1) terms cancel identically once b₁ = λ₂a₁, leaving
2a₂(b₀′ − λ₂a₀′) = 0 and 2a₋₂(b₀′ − μ₂a₀′) = 0, so b₀′ = λ₂a₀′ = μ₂a₀′.

**m = ±1.** The four surviving terms collapse (using b₀′ = λ₂a₀′, resp. μ₂a₀′)
to (λ₂ − μ₂)·U = 0 and (λ₂ − μ₂)·L = 0, where
U := a₂′a₋₁ + 2a₂a₋₁′ and L := a₋₂′a₁ + 2a₋₂a₁′.

**m = 0.** The pairs (±2, ∓2) and (±1, ∓1) telescope to
(λ₂ − μ₂)·d/dτ(2a₂a₋₂ + a₁a₋₁) = const ≠ 0. Hence **λ₂ ≠ μ₂** (which
retroactively gives a₀, b₀ constant from m = ±2), and the **moment**
M(τ) := 2a₂a₋₂ + a₁a₋₁ is **linear nonconstant** in τ.

**Endgame.** From λ₂ ≠ μ₂: U = L = 0. The product invariants
a₋₁·U = (a₂a₋₁²)′ and a₁·L = (a₋₂a₁²)′ give a₂a₋₁² = const and a₋₂a₁² = const.
If a₋₁ ≠ 0, then a₂ = const/a₋₁² is a polynomial only if a₋₁ is constant,
making a₂ constant — a square, contradiction; so a₋₁ = 0, and likewise a₁ = 0.
Then M = 2a₂a₋₂ is linear, so deg a₂ + deg a₋₂ = 1, so one factor is constant —
a square. Contradiction. ∎

## Proof (quantum)

Identical skeleton with the difference calculus; every classical object deforms
by unit shifts:

m = ±4: b₂(E+2)a₂(E) = a₂(E+2)b₂(E) makes b₂/a₂ a 2-periodic rational function,
hence constant: b₂ = λ₂a₂; mirror b₋₂ = μ₂a₋₂.
m = ±3: particular solution b₁ = λ₂a₁; homogeneous solutions exist iff
a₂ = c·h(E)h(E+1) (Lemma J2q), so non-shifted-square a₂ forces b₁ = λ₂a₁;
mirror b₋₁ = μ₂a₋₁.
m = ±2: a₂(E)·[(b₀ − λ₂a₀)(E+2) − (b₀ − λ₂a₀)(E)] = 0, so b₀ − λ₂a₀ is a
2-periodic polynomial, hence constant; mirror with μ₂.
m = ±1: (μ₂ − λ₂)U_q = 0 and (μ₂ − λ₂)L_q = 0 with
U_q = a₋₁(E+2)a₂(E) − a₂(E−1)a₋₁(E), L_q = a₋₂(E+1)a₁(E) − a₁(E−2)a₋₂(E).
m = 0: (μ₂ − λ₂)·[W₂(E+2) − W₂(E) + W₁(E+1) − W₁(E)] = 1 with
W₂ = a₂(E−2)a₋₂(E), W₁ = a₁(E−1)a₋₁(E); hence λ₂ ≠ μ₂, a₀ constant, and the
telescoped moment is linear: W₂ + (lower shift bookkeeping) forces
**W₂ linear** once W₁ = 0.

Endgame: the deformed product invariants — machine-verified as exact
telescopings —
G(E) := a₂(E)·a₋₁(E+1)·a₋₁(E+2) satisfies G(E) − G(E−1) = a₋₁(E+1)·U_q(E),
H(E) := a₋₂(E)·a₁(E−1)·a₁(E−2) satisfies H(E+1) − H(E) = a₁(E−1)·L_q(E),
so U_q = 0 makes G a 1-periodic polynomial, i.e. constant. A nonzero constant
that is a product of three polynomials makes all three constant (a₂ constant =
shifted square, contradiction); G = 0 forces a₋₁ = 0 (a₂ ≠ 0). Likewise a₁ = 0.
Then W₂(E+2) − W₂(E) = 1/(μ₂−λ₂) ≠ 0 makes W₂ = a₂(E−2)a₋₂(E) linear
nonconstant, so one factor is constant — a shifted square. Contradiction. ∎

Note the dictionary held **verbatim at every step**: Wronskian ↔ periodic
ratio; h² = c·a₂ ↔ h(E)h(E+1) = c·a₂; (a₂a₋₁²)′ ↔ ΔG with the triple product
a₂(E)a₋₁(E+1)a₋₁(E+2); the moment (2a₂a₋₂ + a₁a₋₁)′ ↔ the telescoped
W₂, W₁ sum. Nothing quantum-exotic appeared.

## Consequences and the shape of the induction

1. **Band-2 confinement.** With J1/P3 (band 1), every band-2 Dixmier or Keller
pair has each nonzero extreme coefficient a (shifted) square. All remaining
band-2 life — including the shears (x, ξ + x²), whose extreme coefficient is
the square 1 — lives in the square sector. M5 = classify that sector.

2. **Provisional claimed consequence (audit pending).** The draft argues that JC₂ and DC₁ hold vacuously on the
two-sided non-square sector of width ≤ 4, at arbitrary polynomial degree —
orthogonal to every degree bound in the literature (Moh ≤ 100, GGV ≤ 108
excluding (72,108)).

3. **The emerging inductive mechanism (band-k conjecture).** The proof pattern
visibly generalizes: top equations should force a single λ on the whole
positive side and a single μ on the negative side; the cross equations should
force product invariants a_j·(a_{−j'} powers with shifts) constant; and the
m = 0 equation always telescopes the **moment** M = Σ_{j≥1} j·a_j a_{−j}
(suitably shifted, quantumly) to a linear function. In the fully non-square
regime this is self-contradictory exactly as above. What genuinely needs proof
for k ≥ 3 is the collapse of the cross equations — the m = ±3 coupling already
showed collapse is not automatic and is precisely governed by square-type
degeneracies. So the induction's shape is: *non-degenerate top data is
impossible at every band; all candidates descend into power/square strata;
strata descent = the difference analogue of Newton-polygon reduction.* This is
the precise point of contact with Abhyankar–Moh–GGV: their common-leading-form
normal forms are the global version of our square sectors.

## M5 scope (next)

Square sector of band 2: substitute a₂ = h² (classically; a₂ = h(E)h(E+1)
quantumly), allow the homogeneous parameter κ·h in b₁, and push the cascade.
Sub-cases: two-sided square, one-sided (a₋₂ = 0, includes shears), and the
h-degree hierarchy. Outcomes: (a) full band-2 classification = tame moves ⇒
band-2 theorem for JC₂/DC₁ and a template inductive step; (b) a non-tame
solution variety ⇒ the first candidate counterexample geometry, to be pushed
through the membership divisibilities (classically τ² | a₋₂-type; quantumly
E(E−1) | a₋₂-type) exactly as Theorem 1's constraints were satisfied in 3D.
