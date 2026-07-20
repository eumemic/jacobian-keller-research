> **Status: Provisional research archive.**
> Preserved for its derivations and research history. Computational checks pass in the archived environment, but the mathematical arguments and literature claims have not been independently reviewed. This is not a claim of a proof or counterexample to an open conjecture.

# Milestone 2 executed: the quantized tower meets A1

**Outcome in one line.** The tower's bottom layer quantizes *exactly*; the
obstruction to compressing the tower into A1 is now measured precisely (a
falling-factorial divisibility, rigid in the first ladder band); and the band
structure reorganizes DC1 into a width-induction whose base case is proven here
— and whose classical limit is literally JC2. The two prongs are one problem.

All computations below are machine-verified by `verify_prong2.py` (9 checks).

## 1. The crossed-product calculus

Localize: B := A1[x⁻¹] = ⊕_{k∈Z} x^k C[E], E = x∂, with

    (x^a f(E)) · (x^b g(E)) = x^{a+b} f(E+b) g(E),      ∂ = x⁻¹E.

**Membership criterion (the quantum polynomiality constraints).** For j ≥ 0,

    x^{-j} c(E) ∈ A1  ⟺  E(E−1)⋯(E−j+1) divides c(E).

(Witness: ∂^j = x^{-j}·E(E−1)⋯(E−j+1), verified for j = 3.) This is the exact
analogue of Theorem 1's constraints (i)–(iii): classically, s^i t^j / x^m is
polynomial iff 2i+j ≥ m; here τ^j → the falling factorial E^{(j)} under
quantization. Note the classical limit of the criterion is τ^j | c(τ), τ = xξ:
the divisor deforms from a power to a falling factorial.

**The quantum Keller system.** Writing X = Σ_k x^k a_k(E), D = Σ_l x^l b_l(E),
the relation [D, X] = 1 becomes, for every m ∈ Z,

    Σ_{k+l=m} [ b_l(E+k) a_k(E) − a_k(E+l) b_l(E) ] = δ_{m,0}.       (QK)

Its classical limit (symbols along the graded decomposition of C[x,ξ] by
τ = xξ) is, via the verified identity {x^l b, x^k a} = x^{k+l}(k a b′ − l a′ b):

    Σ_{k+l=m} ( k·a_k b_l′ − l·a_k′ b_l )(τ) = δ_{m,0}.              (CK)

(CK) is exactly the Jacobian-1 condition for a plane map written in Laurent
coordinates along τ-levels — i.e. **JC2 is the system (CK); DC1 is its
difference deformation (QK)**. Every element of A1 (or C[x,ξ]) has finite
x-support, so both conjectures are exhausted by the filtration

    Band_k := { pairs with supp(X), supp(D) ⊆ {−k, …, k} },

and "rigidity of Band_k for all k" is *equivalent* to the full conjecture on
each side. This turns "find a trick" into a structured induction on spectral
width (supp = spectrum of ad(E), which acts by x^k f ↦ k·x^k f).

## 2. Finding 1: the bottom layer quantizes exactly

Let N, W satisfy [W, N] = 1 and set (Weyl-ordered with W on the right)

    S̃ = A(N) + N²W,   T̃ = C(N) + 2NW,   A′ = νE(ν),  C′ = 2E(ν).

Then **[T̃, S̃] = 2N²W identically** — verified for a fully symbolic cubic E(ν),
with no lower-order quantum corrections. The classical bracket identity
{S, T} = 2ν²w that powers Theorem 1 survives quantization on the nose. So the
tower is not destroyed by ħ; what fails is the *assembly*, which classically
used a 2-dimensional base of invariants (s, t) and a fiber coordinate x. In B
the invariant subalgebra is C[E] — one dimension, not two — and the only fiber
is the ad(E)-ladder. The next finding measures exactly what that shortfall
costs in the first band.

## 3. Finding 2 (Theorem P3): band-1 rigidity, and the residue obstruction

**Theorem P3.** Let X, D ∈ A1[x⁻¹] with [D, X] = 1 and both supported in
x-degrees {−1, 0, 1}. Then, up to the degenerate subcases listed below, there
are constants with

    X = δx + α − ε∂ + c·x⁻¹/δ′,    D = λX + β′ + (…)x⁻¹,

and precisely one obstruction to A1-membership survives: the residue
**a_{−1}(0)** of the x⁻¹-coefficient. The A1 solutions are exactly the affine
symplectic pairs (images of (x, ∂) under X = ax+b∂+e, D = cx+d∂+f, ad−bc = 1);
the B∖A1 solutions form the one-parameter "polar" families such as (x, ∂+c/x).
In particular **every pair generating a proper subalgebra of A1 must have a
generator with ad(E)-support outside {−1, 0, 1}.**

*Proof.* Write X = xa₁ + a₀ + x⁻¹a₋₁, D = xb₁ + b₀ + x⁻¹b₋₁ (all coefficients
in C[E]) and expand (QK) by degree. Degree ±2: b₁(E+1)a₁(E) = a₁(E+1)b₁(E) and
the mirror equation; a rational function with r(E+1) = r(E) is constant, so
b₁ = λa₁ and b₋₁ = μa₋₁ (machine-checked: the solution space of the degree-2
equation for generic cubic a₁ is exactly the line C·a₁). Degree ±1:
a₁·Δb₀ = b₁·Δa₀ and a₋₁·Δ⁻b₀ = b₋₁·Δ⁻a₀ give Δb₀ = λΔa₀ = μΔa₀; the degree-0
equation will force λ ≠ μ, whence a₀, b₀ are constants α, β. Degree 0: the
surviving terms telescope —

    [D, X]₀ = (λ−μ)·( V(E) − V(E+1) ),     V(E) := a₁(E−1)·a₋₁(E)

(verified symbolically), so V(E) = −E/(λ−μ) + c is *linear*. A product of two
polynomials of degree sum 1 has one constant factor, giving the two branches
(a₁ constant, a₋₁ linear) and its Fourier mirror; back-substituting yields the
affine forms, with [D,X] = 1 equivalent to ad−bc = 1, and A1-membership
equivalent to a₋₁(0) = 0, i.e. c = 0, which kills exactly the polar parameter.
Degenerate subcases (a₁ = 0 or a₋₁ = 0) collapse the same way after applying
the Fourier automorphism x ↦ −∂, ∂ ↦ x, which negates ad(E)-degrees. ∎

**Reading of P3.** Classically, the tower had two free *function* parameters
(E, h) against three *scalar* constraints — room to spare, hence Theorem 1.
Quantumly in Band₁, the difference equations reduce everything to finitely many
scalars before the membership constraint is even imposed; the one surviving
parameter is precisely the polar residue, and membership kills it. The
dimension shortfall (base C[E] instead of C[s,t]) is now a counted fact, not a
metaphor. P3 covers unbounded Bernstein degree (a₁ may have any degree), so it
is a Newton-strip statement, not a low-degree statement; it is presumably
accessible to classical methods (Dixmier's ad-semisimple analysis) — literature
check pending — but its role here is the induction base.

## 4. Finding 3: the two prongs converge

For any pair X, D ∈ A1 with [D, X] = 1 and Bernstein degrees m + n > 2, the top
symbols Poisson-commute ({σD, σX} = σ_{m+n−2}([D,X]) = 0), hence are powers of
a common homogeneous polynomial (Dixmier's dichotomy). So **every DC1 candidate
lives entirely in the maximally degenerate leading-form stratum** — the exact
stratum that is the open gap of the prong-one program (limits of Keller maps
whose leading forms drop rank, where Theorem 4 cannot yet be applied). The two
prongs are the differential and difference faces of one problem:

    JC2  =  (CK) rigid in every band   |   DC1  =  (QK) rigid in every band,

with matching degenerate-leading-form cores, matching membership divisors
(τ^j vs falling factorial), and matching Euler-positivity base mechanisms
(Theorem 4 / P1 / P3). Dimension bookkeeping agrees: the minimal *faithful*
quantum home of the 3D tower is A3, via the cotangent lift of F — which is our
explicit witness that DC3 is false — and compressing 3 → 1 is exactly what the
band rigidity blocks so far.

## 5. Milestone 3 (next): Band₂

Supports ⊆ {−2, …, 2}. The degree ±4 equations again force b₂ = λ₂a₂,
b₋₂ = μ₂a₋₂ (2-periodic rational ⇒ constant). The first genuinely new
phenomenon is degree ±3, where cross-terms couple *different* proportionality
constants:

    b₁(E+2)a₂(E) − a₂(E+1)b₁(E) + b₂(E+1)a₁(E) − a₁(E+2)b₂(E) = 0,

and term-by-term collapse is no longer forced. Solve Band₂ completely: either
it is rigid (then formulate the width induction and hunt for the inductive
mechanism — a difference analogue of the Abhyankar–Moh/Newton-polygon step), or
it is not (then the solution variety contains the first candidate counterexample
geometry, to be pushed through the membership divisibilities). Parallel task:
check Dixmier (1968), Joseph (1975), Guccione–Guccione–Valqui for prior band or
strip results, and screen any Band₂ solution against Result P2's tame mod-p
obstruction for free.
