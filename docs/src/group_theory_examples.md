# Group Theory Examples

Worked examples showing how group theory enables symbolic diagonalization of matrices with $n \geq 5$.

## Pentagon Graph (Dihedral $D_5$)

**Problem**: 5×5 adjacency matrix - normally requires degree-5 polynomial (unsolvable by radicals).

```julia
A = [0 1 0 0 1; 1 0 1 0 0; 0 1 0 1 0; 0 0 1 0 1; 1 0 0 1 0]
eigvals(A)  # [2, φ, φ, -1/φ, -1/φ] where φ = (1+√5)/2
```

**Why it works**: $D_5$ symmetry (5 rotations + 5 reflections) forces factorization into degree 1 + 2 + 2 polynomials. The golden ratio appears from pentagonal geometry.

**Application**: Cyclopentadienyl anion (C₅H₅⁻) molecular orbitals.

---

## Petersen Graph (Strongly Regular)

**Problem**: 10×10 adjacency matrix - degree-10 polynomial.

**Parameters**: $(n, k, \lambda, \mu) = (10, 3, 0, 1)$

**Eigenvalues** (from formula, no polynomial solving):
- $\lambda_1 = 3$ (multiplicity 1)
- $\lambda_2 = 1$ (multiplicity 5)
- $\lambda_3 = -2$ (multiplicity 4)

**Why it works**: Strong regularity constrains to exactly 3 distinct eigenvalues:
$$\lambda_{2,3} = \frac{(\lambda-\mu) \pm \sqrt{(\lambda-\mu)^2 + 4(k-\mu)}}{2}$$

---

## Hypercube $Q_n$

**Problem**: $2^n \times 2^n$ adjacency matrix (exponentially large).

**Eigenvalues** (closed-form):
$$\lambda_k = n - 2k, \quad \text{multiplicity } \binom{n}{k}$$

**Examples**:
- $Q_3$ (8×8): $\{3, 1^{(3)}, -1^{(3)}, -3\}$
- $Q_4$ (16×16): $\{4, 2^{(4)}, 0^{(6)}, -2^{(4)}, -4\}$
- $Q_5$ (32×32): $\{5, 3^{(5)}, 1^{(10)}, -1^{(10)}, -3^{(5)}, -5\}$

**Why it works**: Cayley graph of $(\mathbb{Z}_2)^n$, diagonalized by Walsh-Hadamard basis.

---

## Benzene (Circulant)

**Problem**: 6×6 Hückel Hamiltonian for π-electrons.

```julia
@variables α β
H = [α β 0 0 0 β; β α β 0 0 0; 0 β α β 0 0; 
     0 0 β α β 0; 0 0 0 β α β; β 0 0 0 β α]
```

**Eigenvalues** (DFT formula):
$$\lambda_k = \alpha + 2\beta\cos(2\pi k/6), \quad k = 0, \ldots, 5$$

Explicitly: $\{\alpha + 2\beta, (\alpha + \beta)^{(2)}, (\alpha - \beta)^{(2)}, \alpha - 2\beta\}$

**Application**: Explains aromatic stability - filled bonding orbitals.

---

## Quantum Spin-$j$ Systems

**Problem**: $(2j+1) \times (2j+1)$ angular momentum matrices.

**Eigenvalues** (from SU(2) representation theory):
$$\{-j, -j+1, \ldots, j-1, j\}$$

Always equally spaced by 1, regardless of which component ($J_x$, $J_y$, or $J_z$).

**Example**: Spin-5/2 (6×6 matrix) has eigenvalues $\{-5/2, -3/2, -1/2, 1/2, 3/2, 5/2\}$.

**Why it works**: SU(2) Lie algebra structure determines spectrum completely.

---

## Solvable Degree-5 Polynomials

While generic quintics are unsolvable, special cases work:

| Type | Example | Galois Group | Method |
|------|---------|--------------|--------|
| Pure power | $\lambda^5 - a = 0$ | $\mathbb{Z}_5$ | Roots of unity |
| Reducible | $(\lambda^2+1)(\lambda^3-2)=0$ | Product | Factor and solve |
| Dihedral | Special quintics | $D_5$ | Resolvent + quadratics |

**Pure quintic**: $\lambda^5 = a$ gives $\lambda_k = \sqrt[5]{a} \cdot e^{2\pi ik/5}$

---

## Key Insight

**Symmetry constrains eigenvalue structure.**

| Symmetry | Effect |
|----------|--------|
| Cyclic $\mathbb{Z}_n$ | DFT diagonalizes (circulant) |
| Dihedral $D_n$ | Factorizes into small pieces |
| Strong regularity | Only 3 eigenvalues |
| Tensor product | Eigenvalues multiply |
| SU(2) | Equally spaced spectrum |

These transform intractable problems into closed-form solutions.
