# M4 proof memo: classical band-2 non-square sector

## Status summary

- **Proved:** the complete nine classical band-2 coefficient equations below, with a fixed bracket orientation.
- **Proved:** for nonzero `a_2`, the degree-3 equation has a nonzero homogeneous deviation from proportionality exactly when `a_2` is a scalar times a polynomial square. Over an algebraically closed field this is exactly “`a_2` is a square,” up to the harmless choice of square root of the scalar.
- **Proved:** over an algebraically closed field of characteristic zero, a polynomial Keller pair in band 2 cannot have nonzero nonsquare `a_2`. Thus the proposed `a_2`-nonsquare sector is empty, so the rigidity statement there holds vacuously rather than by an affine classification.
- **Conditional:** over a non-algebraically-closed field, the same proof covers nonsquare `a_2` of positive degree once the degree-3 homogeneous deviation vanishes. A nonsquare constant has an additional residual branch and is not covered by that reduction.
- **Computational evidence only:** the exact verifier checks all signs, all nine equations, the degree-3 reduction, and representative square/nonsquare/localized/affine examples. It is not used as a completeness argument.
- **Status at the time of this memo:** the square sector was open here. A later companion result classifies the oriented \(a_2\)-square sector \(a_2=h^2\ne0\) over \(\mathbb C\) when both entries' ladder supports are contained in \([-2,2]\) and genuine polynomial membership holds. The still later [`../band2-classical-full/`](../band2-classical-full/) package supplies independent band-1 and orientation lemmas and assembles those two proved sectors into the full classical support-contained band-2 theorem over \(\mathbb C\). The nonsquare-constant residual system over nonclosed fields and the quantum analogue remain outside that theorem. The fixed support bound means JC2 is not claimed.

## 1. Scope, package inventory, and checker audit

The supplied ZIP had only two regular files and no paths, links, or special entries. It was inspected before extraction and copied with restrictive permissions into `earlier-package/`:

- `prong2_milestone2_band_rigidity.md`: definitions of the crossed product, quantum and classical coefficient systems, a claimed quantum band-1 theorem, and a band-2 prompt.
- `verify_prong2.py`: nine exact SymPy checks. Its classical scope is only four spot checks of the one-term bracket dictionary. It does not derive band-2 equations or test the J2/M4 claims. Its “band-1 rigidity” check uses one generic cubic for one proportionality step and symbolic checks for a telescoping identity and two examples; therefore it does not itself prove the claimed classification.

The old checker passes unchanged. The package text and checker consistently use `[D,X]` quantumly and the classical orientation fixed below. The warning that an earlier script had opposite orientation is addressed by independently differentiating the full two-variable band-2 expressions in the new verifier.

## 2. Frozen conventions

Let `K` be a field of characteristic zero. The main theorem below assumes `K` algebraically closed (in the intended application, `K=C`). Work in `K[x,xi]`, set `tau=x*xi`, and define

`{G,F} := G_xi F_x - G_x F_xi`.

Thus `{xi,x}=1`. A Keller pair in this memo means `{G,F}=1` in this orientation.

Write

`F = sum_{k=-2}^2 x^k a_k(tau)`,  `G = sum_{l=-2}^2 x^l b_l(tau)`,

with all coefficients in `K[tau]`. This is a Laurent presentation, but polynomiality is imposed:

- for `k>=0`, no additional condition is needed;
- `tau | a_{-1},b_{-1}`;
- `tau^2 | a_{-2},b_{-2}`.

These conditions are necessary and sufficient because `x^k tau^j=x^{k+j}xi^j`. “Band 2” means support is contained in `{-2,-1,0,1,2}`; it does not mean both extreme coefficients are nonzero.

Allowed normalizations in the proof are only:

1. replace `G` by `G-lambda F` (preserves the bracket and support);
2. add constants to either member;
3. diagonal symplectic scaling `(x,xi)->(rho*x,rho^{-1}*xi)`.

The third multiplies `a_2` by `rho^2`, so scalar-square class is invariant. General affine symplectic changes need not preserve this fixed band decomposition and are not silently used. Interchanging the pair reverses the bracket unless accompanied by a sign, and is also not used.

For one-term pieces, direct differentiation gives

`{x^l b(tau),x^k a(tau)} = x^{k+l}(k a b' - l a' b)`.

Consequently the coefficient of `x^m` is

`C_m = sum_{k+l=m}(k a_k b_l' - l a_k' b_l)`,

and the Keller equations are `C_m=delta_{m0}`.

## 3. Complete classical band-2 equations

Primes denote `d/dtau`. No cross-terms are omitted.

```
C_4  = 2 a2 b2' - 2 a2' b2 = 0.

C_3  = 2 a2 b1' - a2' b1
       + a1 b2' - 2 a1' b2 = 0.

C_2  = 2 a2 b0' + a1 b1' - a1' b1 - 2 a0' b2 = 0.

C_1  = 2 a2 b-1' + a2' b-1
       + a1 b0' - a0' b1
       - a-1 b2' - 2 a-1' b2 = 0.

C_0  = 2 a2 b-2' + 2 a2' b-2
       + a1 b-1' + a1' b-1
       - a-1 b1' - a-1' b1
       - 2 a-2 b2' - 2 a-2' b2 = 1.

C_-1 = a1 b-2' + 2 a1' b-2
       + a0' b-1 - a-1 b0'
       - 2 a-2 b1' - a-2' b1 = 0.

C_-2 = 2 a0' b-2 + a-1' b-1
       - a-1 b-1' - 2 a-2 b0' = 0.

C_-3 = 2 a-1' b-2 - a-1 b-2'
       + a-2' b-1 - 2 a-2 b-1' = 0.

C_-4 = 2 a-2' b-2 - 2 a-2 b-2' = 0.
```

The new verifier constructs `F,G`, differentiates them in `x,xi`, and compares the result to all nine displayed expressions identically. This also audits the orientation independently of the old script.

## 4. Exactly what the degree-3 equation proves

Assume `a_2 != 0`. From `C_4=0`, `(b_2/a_2)'=0` in `K(tau)`, hence

`b_2=lambda a_2`, with `lambda in K`.

Set `u=b_1-lambda a_1`. Substitution into `C_3=0` gives exactly

`2 a_2 u' - a_2' u = 0`.                                    (H3)

### Polynomial solutions of (H3)

If `u=0`, (H3) always holds. Thus “the homogeneous sector is nonempty” must mean that a **nonzero deviation** `u` exists; literally, the homogeneous solution space is never empty.

For nonzero `u`, (H3) implies

`(u^2/a_2)'=0`, hence `u^2=c a_2` for some `c in K^*`. Conversely, if `a_2=d h^2` with `d in K^*` and `h in K[tau]`, then every `u=e h` solves (H3), and these are all the polynomial solutions (with `e` arbitrary).

Therefore:

- a nonzero polynomial solution exists iff `a_2` is a **scalar-square**, `d h^2`;
- over algebraically closed `K`, this is equivalent to `a_2` being a literal square;
- over general `K`, “literal square” is too strong: a nonsquare unit times a square still has a nonzero solution;
- a nonzero constant always has nonzero constant solutions, whether or not it is a literal square in `K`;
- `a_2=0` does not permit division and `C_4,C_3` are different equations. Also, zero is usually called a square but has no square-class. It must be treated as a lower-positive-support edge case, not included in the theorem;
- after a field extension, scalar units may acquire square roots and irreducibles may factor, so literal square status is not generally extension-invariant. The intrinsic statement is the existence of `u` satisfying `u^2=c a_2` over the chosen coefficient field.

Diagonal symplectic scaling multiplies `a_2` by a square unit, preserving scalar-square status. The normalization `G->G-lambda F` replaces `b_2` by zero and `b_1` by `u`, making (H3) transparent.

## 5. Non-square attack and proof

### Theorem (proved band-sector statement)

Let `K` be algebraically closed of characteristic zero. Let `F,G in K[x,xi]` have support contained in `{-2,-1,0,1,2}` in the above `tau` decomposition and satisfy `{G,F}=1`. If `a_2 != 0`, then `a_2` is a square in `K[tau]`.

Equivalently, the sector with nonzero nonsquare top coefficient `a_2` contains no polynomial Keller pair.

### Proof

Suppose `a_2` is nonzero and nonsquare. Since `K` is algebraically closed, every nonzero constant is a square, so `deg(a_2)>0`.

From `C_4`, write `b_2=lambda a_2`. Section 4 and nonsquareness force `b_1=lambda a_1`. Then `C_2=0` gives

`2 a_2(b_0'-lambda a_0')=0`,

so `b_0=lambda a_0+c` for `c in K`.

Put `v=b_{-1}-lambda a_{-1}`. Equation `C_1=0` reduces to

`2a_2 v' + a_2'v=0`.

Multiplying by `v` gives `(a_2 v^2)'=0`. If `v!=0`, then `a_2 v^2` is a nonzero constant, impossible because `deg(a_2)>0` and `v` is a polynomial. Hence `v=0`, so `b_{-1}=lambda a_{-1}`.

Finally put `w=b_{-2}-lambda a_{-2}`. The remaining equations reduce, in descending order, to

```
2(a_2 w)' = 1,                         (R0)
a_1 w' + 2a_1' w = 0,                  (R-1)
2a_0' w = 0,                            (R-2)
2a_{-1}'w - a_{-1}w' = 0,              (R-3)
a_{-2}'w - a_{-2}w' = 0.               (R-4)
```

Equation (R0) gives

`a_2 w = tau/2 + gamma`.

Since `w` is a polynomial, `a_2` divides a nonconstant linear polynomial. Because `deg(a_2)>0`, `a_2` is linear and `w` is a nonzero constant. Then (R-1)--(R-4) imply respectively that `a_1,a_0,a_{-1},a_{-2}` are constants. Polynomiality requires `tau|a_{-1}` and `tau^2|a_{-2}`, hence `a_{-1}=a_{-2}=0`.

But `b_{-2}=lambda a_{-2}+w=w` is a nonzero constant, contradicting the necessary polynomiality condition `tau^2|b_{-2}`. This contradiction proves the theorem. QED.

The completeness step is the exact sequential use of `C_4,C_3,...,C_-4`; no bounded-degree ansatz or solver exhaustiveness is involved.

### Localized near-miss and role of polynomiality

The contradiction is sharp. In the Laurent ring, take

`F=x^2 tau`,  `G=(1/2)x^-2`.

Then `{G,F}=1`, `a_2=tau` is nonsquare, and it realizes the reduced solution `w=1/2`. It fails to lie in `K[x,xi]` precisely because `tau^2` does not divide the degree `-2` coefficient `1/2`. Thus polynomiality, not merely the coefficient ODEs, kills the residual branch.

## 6. Residual systems and falsifiable next calculations

### Nonclosed-field constant edge (open here)

If `a_2` is a nonsquare constant in a nonclosed field, (H3) allows arbitrary constant `u`, and the proof's proportionality cascade does not start. After replacing `G` by `G-lambda F`, the exact residual system begins with

`b_2=0`, `b_1=u in K`, `C_2=-a_1'u=0`,

followed by `C_1,...,C_-4` from Section 3 with these substitutions and the divisibilities at negative degrees. A falsifiable next calculation is to eliminate `b_0,b_-1,b_-2` from this residual differential-polynomial system, split `u=0` and `u!=0`, and test whether `tau^2|b_-2` forces a contradiction or an affine descent over `K`.

This edge disappears over `C`, the intended classical coefficient field.

### Square sector (historical next step; subsequently resolved in the companion M5 memo under its stated orientation)

At the time of this memo, the remaining calculation was to write `a_2=h^2` and, after subtracting `lambda F`, write `b_2=0`, `b_1=c h`, then substitute this parametrization into `C_2,...,C_-4` together with `tau|a_-1,b_-1` and `tau^2|a_-2,b_-2`. The later companion M5 memo carries out that arbitrary-degree calculation for the oriented `a_2=h^2!=0` sector over `C`, with both entries' ladder supports contained in `[-2,2]`, and classifies its genuine polynomial members. This historical reduction did not itself establish that later result or supply the remaining orientation/case bridge; the later full classical package supplies that bridge under its fixed support and membership hypotheses.

## 7. Verifier

From the repository root, run:

`python3 research/band2-classical-proved/verify_classical_band2.py`

It checks:

1. the full direct Poisson expansion against the coefficient formula;
2. every displayed `C_m`, including both extremes and all cross-terms;
3. the exact degree-3 reduction;
4. square and nonsquare representatives;
5. the localized nonsquare near-miss and its failed divisibility;
6. an affine symplectic baseline.

All checks pass. This verifier supports the algebraic identities only; the theorem's completeness is the written descent in Section 5.
