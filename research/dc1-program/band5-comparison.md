# Band-5 bounded comparison

**INDEPENDENTLY DERIVED — EXACT FINITE COMPUTATIONS — NOT PEER REVIEWED.**

This note compares selected Band-5 root-necklace and moment-unit computations
with Bands 3 and 4. It is not an exhaustive Band-5 classification and proves
neither DC1 nor JC2. The exact implementation is
[`verify_band5_comparison.py`](verify_band5_comparison.py).

## 1. Wall structure

After gauge `b_5=0`, the top wall is

```text
Q_9=b_4(E+5)a_5(E)-a_5(E+4)b_4(E)=0.
```

In one root coset the necklace equation is `S_5 delta(b_4)=S_4 delta(a_5)`.
For prime `5`, `S_5=Phi_5`; for composite `4`, `S_4=Phi_2 Phi_4`. This verifies
that primality affects factorization of the divisibility gate, while the
non-effective-cofactor phenomenon already exists for every tested `k>=3`.
The universal cofactor `1-S+S^2` supplies one exact Band-5 exotic example.

The step-two arithmetic progression is wall-admissible precisely for odd `k` by
the elementary divisibility of `(S^k+1)/(S+1)`. This local statement does not
classify all top roots.

## 2. Bounded count 13

The number 13 is the result of a bounded, normalized, single-coset integer-root
scan at minimal degree. Roots are distinct integers with minimum translated to
zero; the verifier compares finite maximum-root windows (including 20 and 30).
Both windows return 13 at `k=5`.

Agreement across these windows is a robustness check, not a proof that no top
outside the window exists. The count excludes higher-degree tops, multiple root
cosets, noninteger roots, and repeated-root configurations. The observed counts
`1,4,13` for `k=3,4,5` match a conjectural numerical formula only on those finite
runs. No `k=6` value is entered in the proved or verified ledger.

## 3. What the moment solver actually tests

The moment computation tests five fixed instances, not all 13 scan outputs:

```text
{0,2,3,4,6}
{0,2,4,6,8}
{0,1,4,7,8}
{0,1,3,4,7}
{1,3,4,5,7}
```

These are four normalized shapes plus one translate. The exact `Q_0` test
matrix is:

- for all five `d=1` systems, `Q_0=1` yields the unit ideal and `Q_0=0`
  yields a proper ideal;
- for the sole `d=2` system, with top `{0,2,3,4,6}`, only the `Q_0=1`
  unit-ideal result is tested. No `d=2`, `Q_0=0` feasibility test is made.

For each fixed top and wall solution, the solver imposes the positive cascade
`Q_8,...,Q_1`, solves the listed `b` coefficients with raw degree cap

```text
4d+5,
```

and then appends the coefficient equations for `Q_0` according to this matrix.

These results are exact only within the encoded fixed-top, fixed-degree ansätze
and solver cap. They are not arbitrary-degree statements and do not test all 13
bounded-scan tops. Propriety of a `Q_0=0` ideal means only that this partial
polynomial system has a common zero; it does not satisfy all negative ladder
equations and is not a complete homogeneous Weyl pair.

A genuine Band-5 tame pair is retained as a positive control. The removed
“direct commutator” loop duplicated the definition of `Q_m` and added no
independent validation.

## 4. Status

The computations support, but do not prove, a moment-unit pattern in selected
Band-5 polynomial ansätze. Outstanding issues include arbitrary free degree,
exceptional pivot and denominator loci, unrestricted top geometry, multiple
cosets, and realization of the complete ladder system. There is no completed
Band-5 rung or unrestricted counterexample exclusion.

Run:

```sh
uv run --with sympy python research/dc1-program/verify_band5_comparison.py
```
