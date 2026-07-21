# Band-4 moment-unit experiment

**INDEPENDENTLY DERIVED — EXACT FINITE COMPUTATIONS — NOT PEER REVIEWED.**

This is a bounded experiment in the fixed Band-4 polynomial ansatz implemented
by [`verify_band4_experiment.py`](verify_band4_experiment.py). It is not a
complete Band-4 classification, does not prove a span bound, and proves neither
DC1 nor JC2.

## 1. Fixed conventions and exact identities

Use `A_1[x^{-1}]=sum_i x^i C[E]`, `E=x∂`, with

```text
(x^a f(E))(x^b g(E)) = x^(a+b) f(E+b)g(E),
Q_m = sum_(i+l=m) (b_l(E+i)a_i(E)-a_i(E+l)b_l(E)).
```

Membership requires `E(E-1)...(E-r+1)` to divide each coefficient at level
`-r`. After the top gauge `b_4=0`, the wall is

```text
Q_7=b_3(E+4)a_4(E)-a_4(E+3)b_3(E)=0.
```

On one integer-root coset its necklace equation is `S_4 B=S S_3 A`, where
`S_4=1+S+S^2+S^3=Phi_2 Phi_4`. Since `gcd(S_4,S_3)=1`, the full `S_4` divisibility
is required; divisibility by only one cyclotomic factor is insufficient.

The central identity `Q_0=(T-1)G` and membership implication `G(0)=0` are checked
exactly. Thus `Q_0=1` would require the normalized potential `G=E`.

## 2. Precisely bounded root scan

The reported four exotic tops come only from the normalized, distinct,
single-coset integer-root scan

```text
roots=[0,a,b,c],   1 <= a < b < c <= 15.
```

Within this window the wall-admissible non-shifted-fourth-power tops are

```text
{0,2,3,5}, {0,1,3,6}, {0,3,5,6}, {0,3,6,9}.
```

This finite enumeration neither establishes that every Band-4 top has this form
nor proves that roots outside the window cannot occur. The tetromino/triomino
picture is explanatory only; it is not used as a proved global span bound.
Higher degree, repeated roots, noninteger roots, and multiple root cosets remain
outside the scan.

## 3. Scope of the `Q_0` computations

For each of the four fixed tops and its fixed normalized wall solution `b_3`, the
verifier introduces free lower coefficients of raw degree `d`, imposes membership
on negative levels, and forward-solves the positive cascade `Q_6,...,Q_1`. Each
solved `b` is provisioned with raw degree cap

```text
3d+6.
```

The exact Gröbner results are:

- all four fixed tops at `d=1` and `d=2`;
- selected integer translates of `{0,3,6,9}` at `d=1`; and
- fixed `b_3=0` systems at the tested degree.

In these encoded systems, positive cascade plus `Q_0=1` has unit Gröbner basis,
while positive cascade plus `Q_0=0` has a proper ideal. These are exact statements
about the displayed finite polynomial ansätze, not arbitrary-free-degree or
uniform-translation theorems.

Propriety of the `Q_0=0` ideal means only that this truncated system has a common
zero. It does not supply the remaining negative ladder equations and therefore
does not construct a complete homogeneous Weyl pair.

The verifier also runs a genuine tame Band-4 pair as an anti-false-kill control.
The removed duplicate “direct commutator” loop had merely repeated the definition
of `Q_m`; the retained tests check independent wall, potential, cascade, and
control identities.

## 4. Status

The experiment supports a moment-unit pattern inside the tested
positive-cascade-plus-moment ansätze. Residual work includes arbitrary free
degree, exceptional solver loci, unrestricted tops and root cosets, and assembly
with every negative ladder equation. There is no completed Band-4 rung and no
counterexample exclusion theorem.

Run:

```sh
uv run --with sympy python research/dc1-program/verify_band4_experiment.py
```
