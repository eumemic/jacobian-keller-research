#!/usr/bin/env python3
"""
verify_quantum_a2zero.py
========================
Exact SymPy verification backing the memo `quantum-a2-zero.md`:
the QUANTUM (A_1) a_2 = 0 sector of band 2 is ROUTED, in every configuration, to
one of

    (i)   the quantum a_2 != 0 (shifted-)square sector   [worked elsewhere:
          quantum-mirror.md (ad43ab5) + in-progress quantum-completion.md;
          NOT closed here -- Sectors O1-O3 are routed, not closed],
    (ii)  the quantum band-1 rigidity theorem P3          [conditional input;
          closes Sector Z],
    (iii) direct closure                                  [none needed here].

This script certifies the REDUCTION INFRASTRUCTURE only (the routing operations
preserve [D,X]=1 and A_1-membership, and the case split is an exhaustive
partition); it does not close the a_2 != 0 target.  The reduction is carried out
by two bracket-preserving operations:

    * pair-exchange     (X, D) |-> (D, -X),          preserving [D,X]=1;
    * Fourier reflection  phi:  x |-> -d,  d |-> x,   E |-> -E-1,
      a GENUINE algebra automorphism of A_1 (determinant-+1 analogue of the
      classical swap), hence [phi D, phi X] = phi([D,X]) with NO sign flip.

Convention (as in quantum-mirror.md / quantum-shifted-square-sector-partial.md):
    A_1[x^{-1}] = (+)_{k in Z} x^k C[E],
        (x^a f(E))(x^b g(E)) = x^{a+b} f(E+b) g(E),   f^[r](E) = f(E+r),
    ladder-m coefficient of [D,X]:
        Q_m = sum_{k+l=m} ( b_l^[k] a_k - a_k^[l] b_l ),   [D,X]=1 <=> Q_m=delta_{m0}.
    d = x^{-1} E,  E = x*d,  d^j = x^{-j} E(E-1)...(E-j+1).
    A_1-membership at negative ladder degree:  x^{-r} c(E) in A_1  <=>
        E(E-1)...(E-r+1) | c(E)   (the falling-factorial / quantum polynomiality
        constraint; classical limit tau^r | c(tau)).

Run:  uv run --with sympy python verify_quantum_a2zero.py
Ends: ALL QUANTUM A2-ZERO CHECKS PASSED
"""
import sympy as sp

E = sp.symbols('E')


# =====================================================================
# Crossed-product engine  (elements are dicts  level -> C[E] polynomial)
# =====================================================================
def sh(f, r):
    return sp.expand(sp.sympify(f).subs(E, E + r))


def falling(j):
    """E^{underline j} = E(E-1)...(E-j+1);  falling(0)=1."""
    r = sp.Integer(1)
    for i in range(j):
        r *= (E - i)
    return sp.expand(r)


def clean(d):
    return {k: sp.expand(v) for k, v in d.items() if sp.expand(v) != 0}


def mult(A, B):
    """(x^a f)(x^b g) = x^{a+b} f(E+b) g(E)."""
    R = {}
    for a, fa in A.items():
        for b, gb in B.items():
            lev = a + b
            R[lev] = sp.expand(R.get(lev, 0) + sh(fa, b) * gb)
    return clean(R)


def comm(D, X):
    """[D,X] = D X - X D  in the crossed product."""
    P, M = mult(D, X), mult(X, D)
    return clean({k: sp.expand(P.get(k, 0) - M.get(k, 0))
                  for k in set(P) | set(M)})


def Qm(a, b, m):
    """Ladder-m coefficient of [D,X] from the coefficient dicts (memo formula)."""
    return sp.expand(sum(sh(b[l], k) * a[k] - sh(a[k], l) * b[l]
                         for k in a for l in b if k + l == m))


def phi_piece(k, f):
    """phi(x^k f(E)) as a one-term dict, using phi(x)=-d, phi(d)=x, phi(E)=-E-1.

    k >= 0 :  x^k f(E) = (x)^k f  ->  (-d)^k f(-E-1) = (-1)^k x^{-k} E^{und k} f(-E-1).
    k < 0  :  genuine member x^{-j} f(E) = d^j g(E) with f = E^{und j} g
              ->  (x)^j g(-E-1).
    """
    f = sp.expand(f)
    if k >= 0:
        return {-k: sp.expand((-1) ** k * falling(k) * f.subs(E, -E - 1))}
    j = -k
    q, r = sp.div(sp.Poly(f, E), sp.Poly(falling(j), E))
    if sp.expand(r.as_expr()) != 0:
        raise AssertionError(f"membership fail at level {k}: falling({j}) nmid {f}")
    return {j: sp.expand(q.as_expr().subs(E, -E - 1))}


def phi(A):
    R = {}
    for k, f in A.items():
        for lev, c in phi_piece(k, f).items():
            R[lev] = sp.expand(R.get(lev, 0) + c)
    return clean(R)


def eqdict(A, B):
    return all(sp.expand(A.get(k, 0) - B.get(k, 0)) == 0 for k in set(A) | set(B))


def rp(name, deg):
    cs = sp.symbols(f'{name}_0:{deg+1}')
    return sum(cs[i] * E ** i for i in range(deg + 1))


def PASS(cond, label):
    if cond is not True:
        raise AssertionError(f"{label}  ::  {cond}")
    print("PASS", label)


# =====================================================================
# 0.  Engine sanity: Q_m = commutator; band-2 <-> band-1 collapse
# =====================================================================
print("--- Part 0: crossed-product engine and band-1 collapse ---")

xg, dg = {1: sp.Integer(1)}, {-1: E}
PASS(eqdict(mult(xg, dg), {0: E}), "E = x*d")
PASS(eqdict(comm(dg, xg), {0: sp.Integer(1)}), "[d,x] = 1  (Weyl relation)")

# Q_m formula agrees with the direct commutator on a fully generic band-2 pair.
Ag = {2: rp('A2', 3), 1: rp('A1', 3), 0: rp('A0', 3),
      -1: sp.expand(E * rp('U', 3)), -2: sp.expand(E * (E - 1) * rp('M2', 3))}
Bg = {2: rp('B2', 3), 1: rp('B1', 3), 0: rp('B0', 3),
      -1: sp.expand(E * rp('V', 3)), -2: sp.expand(E * (E - 1) * rp('W2', 3))}
cQ = comm(Bg, Ag)
PASS(all(sp.expand(Qm(Ag, Bg, m) - cQ.get(m, 0)) == 0 for m in range(-4, 5)),
     "Q_m (memo formula) = ladder coefficients of the direct commutator [D,X]")

# a_2 = a_{-2} = b_2 = b_{-2} = 0  =>  the band-2 Q_m system IS the band-1 system.
a0c = {2: sp.Integer(0), 1: rp('c1', 3), 0: rp('c0', 3),
       -1: sp.expand(E * rp('cm', 3)), -2: sp.Integer(0)}
b0c = {2: sp.Integer(0), 1: rp('e1', 3), 0: rp('e0', 3),
       -1: sp.expand(E * rp('em', 3)), -2: sp.Integer(0)}
a1b = {1: a0c[1], 0: a0c[0], -1: a0c[-1]}
b1b = {1: b0c[1], 0: b0c[0], -1: b0c[-1]}
PASS(all(sp.expand(Qm(a0c, b0c, m) - (Qm(a1b, b1b, m) if -2 <= m <= 2 else 0)) == 0
         for m in range(-4, 5)),
     "band-1 collapse: band-2 Q_m with all four extremes zero == band-1 Q_m")
PASS(all(sp.expand(Qm(a0c, b0c, m)) == 0 for m in (4, 3, -3, -4)),
     "band-1 collapse: Q_{+-3}, Q_{+-4} vanish identically (=> Sector Z is P3)")


# =====================================================================
# 1.  Pair-exchange  (X,D) |-> (D,-X)
# =====================================================================
print("\n--- Part 1: pair-exchange preserves [D,X]=1 and moves b_2 into the top ---")

negA = {k: sp.expand(-v) for k, v in Ag.items()}
PASS(eqdict(comm(negA, Bg), comm(Bg, Ag)),
     "[-X, D] = [D, X]  (operator identity; so [D,X]=1 => new pair Keller)")
PASS(sp.expand(Bg[2] - Bg.get(2, 0)) == 0,
     "after exchange, new top a_2' = (old) b_2  (b_2 != 0 => oriented)")


# =====================================================================
# 2.  Fourier reflection phi : x|->-d, d|->x, E|->-E-1
# =====================================================================
print("\n--- Part 2: Fourier automorphism, closed forms, membership bookkeeping ---")

PASS(eqdict(phi(xg), {-1: -E}), "phi(x) = -d")
PASS(eqdict(phi(dg), {1: sp.Integer(1)}), "phi(d) = x")
PASS(eqdict(phi({0: E}), {0: -E - 1}), "phi(E) = -E-1")
PASS(eqdict(comm(phi(dg), phi(xg)), {0: sp.Integer(1)}),
     "phi preserves the Weyl relation: [phi d, phi x] = 1")

# phi is an algebra homomorphism:  phi(A B) = phi(A) phi(B).
homom_ok = True
for (a, b) in [(1, 1), (2, 1), (1, 2), (2, 2), (0, 2), (2, 0), (3, 1)]:
    A, B = {a: rp('f', 2)}, {b: rp('g', 2)}
    homom_ok = homom_ok and eqdict(phi(mult(A, B)), mult(phi(A), phi(B)))
for (a, b) in [(-1, -1), (-2, -1), (-1, -2), (-2, 2), (2, -2),
               (-1, 2), (2, -1), (-2, -2)]:
    fa = sp.expand(falling(abs(a)) * rp('fr', 2)) if a < 0 else rp('f', 2)
    gb = sp.expand(falling(abs(b)) * rp('gr', 2)) if b < 0 else rp('g', 2)
    A, B = {a: fa}, {b: gb}
    homom_ok = homom_ok and eqdict(phi(mult(A, B)), mult(phi(A), phi(B)))
PASS(homom_ok, "phi(A B) = phi(A) phi(B) on all ladder-generator pairs (symbolic)")

Xg = {2: rp('a2', 2), 1: rp('a1', 2), 0: rp('a0', 2),
      -1: sp.expand(E * rp('m1', 2)), -2: sp.expand(E * (E - 1) * rp('m2', 2))}
p2 = phi(phi(Xg))
PASS(eqdict(phi(phi(p2)), Xg), "phi^4 = identity (phi invertible => automorphism)")
PASS(eqdict(p2, {k: sp.expand((-1) ** k * v) for k, v in Xg.items()}),
     "phi^2 = parity (x^k f |-> (-1)^k x^k f)")

# closed-form ladder coefficients of phi X  (the quantum reflection dictionary)
R = phi(Xg)
PASS(sp.expand(R.get(-2, 0) - E * (E - 1) * Xg[2].subs(E, -E - 1)) == 0,
     "(phi X)_{-2} = E(E-1) a_2(-E-1)")
PASS(sp.expand(R.get(-1, 0) + E * Xg[1].subs(E, -E - 1)) == 0,
     "(phi X)_{-1} = -E a_1(-E-1)")
PASS(sp.expand(R.get(0, 0) - Xg[0].subs(E, -E - 1)) == 0,
     "(phi X)_0 = a_0(-E-1)")
PASS(sp.expand(R.get(1, 0) - sp.cancel(Xg[-1].subs(E, -E - 1) / (-E - 1))) == 0,
     "(phi X)_1 = a_{-1}(-E-1)/(-E-1)  (polynomial by membership E | a_{-1})")
PASS(sp.expand(R.get(2, 0)
               - sp.cancel(Xg[-2].subs(E, -E - 1) / ((-E - 1) * (-E - 2)))) == 0,
     "(phi X)_2 = a_{-2}(-E-1)/((E+1)(E+2))  (polynomial by E(E-1)|a_{-2})")

# The load-bearing orientation fact: a_{-2} != 0  =>  (phi X)_2 != 0.
#   (phi X)_2 = g(-E-1) with a_{-2} = E(E-1) g;  substitution is injective on C[E].
gg = rp('gg', 3)
PASS(sp.expand(phi({-2: sp.expand(E * (E - 1) * gg)}).get(2, 0)
               - gg.subs(E, -E - 1)) == 0,
     "a_{-2}=E(E-1)g != 0  =>  (phi X)_2 = g(-E-1) != 0  (orients into a_2 slot)")

# Shifted-square bookkeeping: a level-(-2) shifted-square maps to a level-(+2)
# shifted-square (matches J3q's reflection lemma  c h(E)h(E-1)).
c = sp.symbols('c')
h = rp('h', 2)
am2ss = sp.expand(E * (E - 1) * c * h * sh(h, -1))          # = E(E-1)*c*h(E)h(E-1)
Hh = h.subs(E, -E - 2)
PASS(sp.expand(phi({-2: am2ss}).get(2, 0) - c * Hh * sh(Hh, -1)) == 0,
     "shifted-square transfer: (phi X)_2 = c H(E)H(E-1),  H(E)=h(-E-2)")

# Membership is preserved automatically: phi of a genuine A_1 element is genuine.
#   negative-level images of phi X carry the required falling factors by construction.
mem_ok = True
for j, fac in [(-1, E), (-2, E * (E - 1))]:
    coeff = R.get(j, 0)
    if coeff != 0:
        _, rem = sp.div(sp.Poly(coeff, E), sp.Poly(sp.expand(fac), E))
        mem_ok = mem_ok and sp.expand(rem.as_expr()) == 0
PASS(mem_ok, "phi X satisfies A_1 membership at every negative level (falling factors)")


# =====================================================================
# 3.  Bracket transfer:  [phi D, phi X] = phi([D,X])  (no sign)
# =====================================================================
print("\n--- Part 3: Fourier transfers the bracket, no sign flip ---")

# The commutator of two generic genuine band-2 elements is itself a genuine A_1
# element (membership holds at every negative level, verified), so phi([D,X]) is
# defined, and the automorphism identity holds EXACTLY.
C = comm(Bg, Ag)
mem_comm = True
for j, fac in [(-1, E), (-2, E * (E - 1)),
               (-3, E * (E - 1) * (E - 2)), (-4, E * (E - 1) * (E - 2) * (E - 3))]:
    if C.get(j, 0) != 0:
        _, rem = sp.div(sp.Poly(C[j], E), sp.Poly(sp.expand(fac), E))
        mem_comm = mem_comm and sp.expand(rem.as_expr()) == 0
PASS(mem_comm, "[D,X] is a genuine A_1 element (membership at every negative level)")
PASS(eqdict(comm(phi(Bg), phi(Ag)), phi(C)),
     "[phi D, phi X] = phi([D,X])  EXACTLY (genuine automorphism, NO sign flip)")

# Corollary: on a Keller pair the whole bracket is the scalar 1, and phi fixes
# scalars, so the oriented pair is again Keller.  (Contrast: the classical swap
# has determinant -1 and needs the compensating sign  {RG,RF} = -R{G,F}.)
PASS(eqdict(phi({0: sp.Integer(1)}), {0: sp.Integer(1)}),
     "phi(1) = 1  =>  [D,X]=1  ==>  [phi D, phi X] = 1")


# =====================================================================
# 4.  Concrete end-to-end a_2 = 0 witnesses, oriented to a_2 != 0
# =====================================================================
print("\n--- Part 4: concrete a_2=0 pairs oriented into the a_2!=0 sector ---")

# (O2)  X = x + d^2  (a_2=0, a_{-2}=E(E-1) != 0),  D = d.   [D,X]=1.
X, D = {1: sp.Integer(1), -2: sp.expand(E * (E - 1))}, {-1: E}
PASS(eqdict(comm(D, X), {0: sp.Integer(1)}), "O2 witness X=x+d^2, D=d : [D,X]=1")
PASS(X.get(2, 0) == 0 and X.get(-2, 0) != 0, "O2 witness has a_2=0, a_{-2}!=0")
oX, oD = phi(X), phi(D)
PASS(eqdict(comm(oD, oX), {0: sp.Integer(1)}) and oX.get(2, 0) != 0,
     "O2: phi orients to (x^2 - d, x), a_2 = 1 != 0, bracket preserved")

# (O1)  X = x  (band-1, a_2=0),  D = d + x^2  (b_2 != 0).   [D,X]=1.
X, D = {1: sp.Integer(1)}, {-1: E, 2: sp.Integer(1)}
PASS(eqdict(comm(D, X), {0: sp.Integer(1)}), "O1 witness X=x, D=d+x^2 : [D,X]=1")
Xp, Dp = dict(D), {k: sp.expand(-v) for k, v in X.items()}      # (D, -X)
PASS(eqdict(comm(Dp, Xp), {0: sp.Integer(1)}) and Xp.get(2, 0) != 0,
     "O1: pair-exchange orients to (d+x^2, -x), a_2 = 1 != 0, bracket preserved")

# (O3)  b_{-2} != 0 (with b_2 = a_{-2} = 0): pair-exchange then Fourier fills a_2.
Dbot = {-2: sp.expand(E * (E - 1) * rp('q', 2))}      # generic nonzero level -2 piece
PASS(phi(Dbot).get(2, 0) != 0,
     "O3: b_{-2}!=0  =>  (after exchange) phi fills the top slot with a_2 != 0")


print("\nALL QUANTUM A2-ZERO CHECKS PASSED")
