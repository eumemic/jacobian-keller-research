# The characteristic-p candidate sieve for DC1

**INDEPENDENTLY DERIVED — EXACT SYMPY CHECKS — NOT PEER REVIEWED — PROTOTYPE ONLY.**

This note documents a necessary-candidate screen, not an exclusion theorem. The
implementation is [`sieve_dc1_candidate.py`](sieve_dc1_candidate.py), with exact
regressions in [`verify_charp_sieve.py`](verify_charp_sieve.py). A finite prime
sweep is corroboration only.

## 1. Exact local algebra

In characteristic `p`, normal-ordered Weyl monomials satisfy

```text
(d^j x^i) = sum_k binom(j,k) binom(i,k) k! x^(i-k)d^(j-k).
```

The center of `A_1(F_p)` is exactly `F_p[x^p,d^p]`: commuting with `d` forces
every `x` exponent to be divisible by `p`, and commuting with `x` forces the same
for every `d` exponent. Thus a central element is read uniquely as a polynomial
in `u=x^p`, `v=d^p`.

For a rational Weyl pair `(X,D)` with `[D,X]=1`, at a prime not dividing a
coefficient denominator, the prototype computes `X^p,D^p` exactly modulo `p`.
It accepts the induced center map

```text
psi_p=(P,Q),  P=X^p, Q=D^p in F_p[u,v]
```

only after checking both powers are central. The general bridge and its Poisson
compatibility are literature inputs (Tsuchimoto; Belov-Kanel–Kontsevich), not
new theorems proved by these finite computations.

The verifier checks the center calculation, Weyl arithmetic, centrality on its
listed examples, and `Jac(P,Q)=1` exactly. It also retains the exact local
Artin–Schreier calculation: `(u,v-v^p)` has Jacobian one and geometric degree
`p`, but is not an automorphism. This is the essential wild-characteristic
warning.

## 2. Degree computation and production gates

The prototype forms both resultants

```text
Res_v(P-a,Q-b),    Res_u(P-a,Q-b)
```

and reduces their coefficients modulo `p`. It reports a geometric degree only
when both remaining-variable degrees are available and equal. If either is
missing or the two disagree, the degree is `None`; it never takes their maximum.

A production obstruction verdict is gated on all of the following:

1. `[D,X]=1` and denominator/good-prime checks pass;
2. `X^p,D^p` are central;
3. the induced center map has Jacobian exactly one;
4. both elimination degrees are available and agree;
5. good reduction and constancy of geometric degree are established, rather than
   inferred from a bounded sample; and
6. every hypothesis of an external finite tame-extension theorem is proved for
   the candidate.

The implementation verifies items 1–4 per sampled prime and passes those facts
explicitly to the verdict interface; every evidence parameter otherwise defaults
to `False` or `None`. A direct verdict call without that evidence is therefore
`INCONCLUSIVE`. Reports expose the per-prime `local_checks` separately and set
`global_good_reduction_established=False`. The implementation does not establish
items 5–6 and never emits a theorem-level exclusion. A degree-two reading at
`p>2` is at most `OBSTRUCTION_CANDIDATE (non-theorem)`, and only when every local
prerequisite is explicitly verified as the Boolean value `True`. A candidate at
one sampled good prime cannot override a bridge failure or failed/missing evidence
at another attempted good prime; mixed evidence is `INCONCLUSIVE`. Primes dividing
a coefficient denominator are outside the good-reduction domain and are explicitly
skipped rather than treated as evidence. Other failures or unsupported cases are
also `INCONCLUSIVE`.

## 3. Missing global geometry

Turning the local Kummer intuition into an exclusion requires a theorem that is
not established here. At minimum one must prove and track:

- generic finiteness and separability of the center map;
- finiteness of the normalization over `A^2`;
- absence of affine ramification;
- control of ramification and components at infinity;
- applicability of purity in the chosen finite model;
- good reduction and constant geometric degree across the relevant primes;
- control of the Galois closure and intermediate/subextension degrees; and
- tameness of the relevant quotient or cyclic subextension.

A constant Jacobian by itself does not supply this package. In particular, the
squarefree-pullback and tame Kummer local algebra may be valid once the required
finite étale/tame setting has been constructed, but those local statements do
not construct that setting. They therefore remain conditional ingredients, not
a proved tame-extension theorem for arbitrary Keller center maps.

## 4. What the checks establish

The exact verifier covers:

- normal-order Weyl arithmetic and the center criterion;
- bridge centrality and Jacobian one on the named tame examples;
- degree-one results for those examples;
- denominator skipping and loud rejection of broken Weyl pairs;
- exact degree readings for synthetic polynomial maps;
- the Artin–Schreier wild example;
- suppression of verdicts for a non-Keller center map;
- suppression of a degree when the two elimination computations disagree; and
- non-theorem candidate semantics for synthetic `(u^2,v)`.

The synthetic map `(u^2,v)` is only a pipeline test. It is not the center map of
a supplied Weyl candidate and proves no exclusion theorem.

A sampled degree-one center map is consistent with automorphy but does not prove
the original Weyl endomorphism is an automorphism. A bounded prime sweep cannot
certify an almost-all-primes statement.

## 5. Global status

Nothing here proves DC1, JC2, a Band-4 or Band-5 classification, or a
counterexample. If the characteristic-zero bridge is discussed globally, the
safe implication used by this project is only `JC2 => DC1`; no converse or
unrestricted equivalence is asserted here.

The characteristic-p pipeline is best viewed as an auditable source of necessary
computations and possible obstruction candidates while the finite tame-extension
geometry remains unresolved.

## 6. Reproduce

```sh
uv run --with sympy python research/dc1-program/verify_charp_sieve.py
uv run --with sympy python research/dc1-program/sieve_dc1_candidate.py
```

References: Y. Tsuchimoto, *Endomorphisms of Weyl algebra and p-curvatures*,
Osaka J. Math. 42 (2005); A. Belov-Kanel and M. Kontsevich, *The Jacobian
conjecture is stably equivalent to the Dixmier conjecture*, Moscow Math. J. 7
(2007); P. Revoy, *Algèbres de Weyl en caractéristique p* (1973).
