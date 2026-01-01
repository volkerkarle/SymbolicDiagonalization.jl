# Mathematical Background

Mathematical foundations for symbolic matrix diagonalization.

## The Eigenvalue Problem

Given matrix $\mathbf{A}$, find $\lambda$ and $\mathbf{v}$ such that:
$$\mathbf{A}\mathbf{v} = \lambda\mathbf{v}$$

This requires:
$$\det(\lambda\mathbf{I} - \mathbf{A}) = 0$$

The **characteristic polynomial** is degree $n$, so eigenvalues are polynomial roots.

## The Abel-Ruffini Theorem

**Theorem** (1824): No general algebraic formula exists for polynomial roots of degree $\geq 5$.

| Degree | Method | Status |
|--------|--------|--------|
| 1 | Linear | $\lambda = -c_0$ |
| 2 | Quadratic | $\lambda = \frac{-b \pm \sqrt{b^2-4c}}{2}$ |
| 3 | Cardano | Nested radicals |
| 4 | Ferrari | Reduces to cubic + quadratics |
| ≥5 | None | **Impossible in general** |

**Implication**: Generic $n \times n$ matrices with $n \geq 5$ cannot have symbolic eigenvalues.

**Exception**: Specific polynomials with special structure (e.g., $\lambda^5 - 1 = 0$) may still be solvable.

## Why Structure Enables Solutions

Matrices with structure have **special Galois groups**:

| Structure | Galois Group | Solvable? |
|-----------|--------------|-----------|
| Circulant | Cyclic $\mathbb{Z}_n$ | Always |
| Dihedral symmetry | Dihedral $D_n$ | Always |
| Generic matrix | Symmetric $S_n$ | No ($n \geq 5$) |

## Closed-Form Root Formulas

### Quadratic ($ax^2 + bx + c = 0$)

$$x = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}$$

### Cubic ($x^3 + px + q = 0$, depressed form)

$$x = \sqrt[3]{-\frac{q}{2} + \sqrt{\frac{q^2}{4} + \frac{p^3}{27}}} + \sqrt[3]{-\frac{q}{2} - \sqrt{\frac{q^2}{4} + \frac{p^3}{27}}}$$

### Quartic

Reduces to solving a cubic resolvent, then two quadratics. Expressions are very large (~13.5 MB for fully symbolic 4×4).

## Pattern-Specific Theory

### Circulant Matrices

**Why closed-form exists**: All circulant matrices are polynomials in the shift operator $\mathbf{S}$. The DFT matrix diagonalizes $\mathbf{S}$, hence all circulants.

$$\lambda_k = \sum_{j=0}^{n-1} c_j \omega^{jk}, \quad \omega = e^{2\pi i/n}$$

### Symmetric Toeplitz Tridiagonal

**Why closed-form exists**: Known eigenvector basis (discrete sine transform). Eigenvalue equation reduces to trigonometric identity.

$$\lambda_k = a + 2b\cos\left(\frac{k\pi}{n+1}\right)$$

### Kronecker Products

**Why closed-form exists**: Tensor product structure factorizes eigenspaces.

$$\lambda(\mathbf{A} \otimes \mathbf{B}) = \{\lambda_i(\mathbf{A}) \cdot \lambda_j(\mathbf{B})\}$$

### Permutation Matrices

**Why closed-form exists**: Finite order implies eigenvalues are roots of unity. Cycle structure determines which roots appear.

## Group Theory Foundations

### Representation Theory

When $\mathbf{A}$ commutes with symmetry group $G$:
1. Eigenvectors organize into irreducible representations
2. Eigenvalues within an irrep are equal (degeneracy)
3. Character tables predict eigenspace dimensions

### Key Groups

| Group | Structure | Application |
|-------|-----------|-------------|
| $\mathbb{Z}_n$ | Cyclic rotations | Circulant matrices |
| $D_n$ | Rotation + reflection | Polygon graphs |
| $S_n$ | All permutations | Generic case (bad) |

### Strongly Regular Graphs

Parameters $(n, k, \lambda, \mu)$ determine exactly 3 eigenvalues:
- $\lambda_1 = k$
- $\lambda_{2,3} = \frac{(\lambda-\mu) \pm \sqrt{(\lambda-\mu)^2 + 4(k-\mu)}}{2}$

Example: Petersen graph $(10,3,0,1)$ → $\{3, 1^{(5)}, -2^{(4)}\}$

### Hypercube Graphs

$n$-dimensional hypercube $Q_n$:
$$\lambda_k = n - 2k, \quad \text{multiplicity } \binom{n}{k}$$

## Pattern Discovery Principles

A pattern has closed-form eigenvalues when:

1. **Known diagonalizing transform** (e.g., DFT for circulants)
2. **Reducible structure** (e.g., block-diagonal splits problem)
3. **Tensor product structure** (eigenvalues are products)
4. **Symmetry constraints** (representation theory)
5. **Solvable Galois group** (cyclic, dihedral, etc.)

## References

- Horn & Johnson, *Matrix Analysis* (2012)
- Davis, *Circulant Matrices* (1979)
- Godsil & Royle, *Algebraic Graph Theory* (2001)
- Serre, *Linear Representations of Finite Groups* (1977)
