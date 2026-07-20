> **Status: Provisional research archive.**
> Preserved for its derivations and research history. Computational checks pass in the archived environment, but the mathematical arguments and literature claims have not been independently reviewed. This is not a claim of a proof or counterexample to an open conjecture.

# An explicit injective, non-surjective endomorphism of the Weyl algebra A_3

Weyl algebra A_3 = C< x1,x2,x3, d1,d2,d3 > with [di,xj] = delta_ij, [xi,xj]=[di,dj]=0.

Built from the Keller map F (det JF = -2) via phi(xi)=Fi, phi(dj) = sum_k (JF^-1)_{kj} dk.
Below, u = 1+x1*x2 (writing x=x1, y=x2, z=x3).

## Images of the coordinate generators

phi(x1) = (x*y + 1)*(x**2*y**2*z + 3*x*y**3 + 2*x*y*z + 4*y**2 + z)
phi(x2) = 3*x**3*y**2*z + 9*x**2*y**3 + 6*x**2*y*z + 12*x*y**2 + 3*x*z + y
phi(x3) = -x*(x**2*z + 3*x*y - 2)

## Images of the derivations   phi(dj) = A_j*d1 + B_j*d2 + C_j*d3

phi(d1):
  coeff of d1: x**3*(3*x**3*y*z + 9*x**2*y**2 + 3*x**2*z + 3*x*y - 4)
  coeff of d2: 3*x*(x**3*y*z + 3*x**2*y**2 + x**2*z + x*y - 1)
  coeff of d3: -9*x**5*y*z**2 - 45*x**4*y**2*z - 9*x**4*z**2 - 54*x**3*y**3 - 30*x**3*y*z - 27*x**2*y**2 + 9*x**2*z + 21*x*y + 1

phi(d2):
  coeff of d1: -x**2*(3*x**4*y**2*z + 9*x**3*y**3 + 6*x**3*y*z + 12*x**2*y**2 + 3*x**2*z - x*y - 3)/2
  coeff of d2: -(3*x**4*y**2*z + 9*x**3*y**3 + 6*x**3*y*z + 12*x**2*y**2 + 3*x**2*z - 2)/2
  coeff of d3: (9*x**5*y**2*z**2 + 45*x**4*y**3*z + 18*x**4*y*z**2 + 54*x**3*y**4 + 75*x**3*y**2*z + 9*x**3*z**2 + 81*x**2*y**3 + 21*x**2*y*z + 6*x*y**2 - 6*x*z - 16*y)/2

phi(d3):
  coeff of d1: -(x*y + 1)**2*(3*x**4*y**2*z + 9*x**3*y**3 + 6*x**3*y*z + 12*x**2*y**2 + 3*x**2*z - x*y - 1)/2
  coeff of d2: -3*(x*y + 1)**2*(x**2*y**2*z + 3*x*y**3 + 2*x*y*z + 4*y**2 + z)/2
  coeff of d3: (9*x**5*y**4*z**2 + 45*x**4*y**5*z + 36*x**4*y**3*z**2 + 54*x**3*y**6 + 165*x**3*y**4*z + 54*x**3*y**2*z**2 + 189*x**2*y**5 + 216*x**2*y**3*z + 36*x**2*y*z**2 + 222*x*y**4 + 117*x*y**2*z + 9*x*z**2 + 89*y**3 + 21*y*z)/2

## Compressed forms (t = xy, s = x^2 z, u = 1+t, P = 3s*u^2 + 3t^2*(3t+4))

phi(d1) = x^3*(3su + 9t^2 + 3t - 4) d1 + 3x*(su + 3t^2 + t - 1) d2 + [C_1] d3
phi(d2) = -(x^2/2)*(P - t - 3) d1 - (1/2)*(P - 2) d2 + [C_2] d3
phi(d3) = -(u^2/2)*(P - u) d1 - (3/2)*u*F1 d2 + [C_3] d3
  (C_1, C_2, C_3 are the d3-coefficients listed in full above)

## Verified relations (sympy, exact symbolic computation)

- [phi(dj), phi(xi)] = delta_ij for all i,j   (equivalent to JF * JF^-1 = I)
- [phi(di), phi(dj)] = 0 for all i,j          (vector-field commutators vanish identically)
- [phi(xi), phi(xj)] = 0                       (multiplication operators)
- det JF = -2, so JF^-1 = -adj(JF)/2 is a polynomial matrix and phi lands in A_3

## Why injective but not surjective

Injective: A_3 is a simple algebra and phi(1)=1, so ker(phi)=0.
Not surjective: the standard easy implication DC_3 => JC_3 shows that if phi_F is an
automorphism then the Keller map F is invertible. F is generically 3-to-1 (e.g. the fiber
over (-1/4,0,0) is {(0,0,-1/4), (1,-3/2,13/2), (-1,3/2,13/2)}), hence not invertible,
so phi is a proper endomorphism: a counterexample to the Dixmier conjecture for A_3.