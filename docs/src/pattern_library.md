# Pattern Library

Catalog of matrix patterns with closed-form eigenvalue solutions.

## Quick Reference

| Pattern | Formula/Method | Complexity |
|---------|---------------|------------|
| Diagonal | Direct: $\lambda_i = d_i$ | $O(n)$ |
| Triangular | Diagonal entries | $O(n)$ |
| Block-diagonal | Union of block eigenvalues | $\sum O(n_i)$ |
| Circulant | DFT: $\lambda_k = \sum c_j \omega^{jk}$ | $O(n)$ |
| Sym. Toeplitz tridiag | $\lambda_k = a + 2b\cos(k\pi/(n+1))$ | $O(n)$ |
| Kronecker $A \otimes B$ | Products: $\lambda_i(A) \cdot \lambda_j(B)$ | $O(m+n)$ |
| Nested Kronecker $A_1 \otimes \cdots \otimes A_k$ | Product of all factors | $O(\sum n_i)$ |
| SO(2) Kronecker $R(\theta) \otimes R(\phi)$ | $e^{\pm i(\theta \pm \phi)}$ | $O(1)$ |
| SO(n) for n ≤ 4 | Trace-based symbolic formulas | $O(1)$ |
| Permutation | Roots of unity from cycle structure | $O(n)$ |
| Persymmetric | Half-size decomposition | $O(n/2)$ |
| Anti-diagonal | $\pm$ pairs | $O(n)$ |
| Jordan block | Single eigenvalue $\lambda$ with multiplicity $n$ | $O(1)$ |

---

## Basic Patterns

### Diagonal

**Eigenvalues**: Diagonal entries $\{d_1, d_2, \ldots, d_n\}$

```julia
@variables a b c
D = [a 0 0; 0 b 0; 0 0 c]
eigvals(D)  # [a, b, c]
```

### Triangular

**Eigenvalues**: Diagonal entries (upper or lower triangular)

```julia
@variables a b c
U = [a 1 2; 0 b 3; 0 0 c]
eigvals(U)  # [a, b, c]
```

### Jordan Blocks

**Eigenvalues**: $\lambda$ with algebraic multiplicity $n$, geometric multiplicity 1

```julia
@variables λ
J = [λ 1 0; 0 λ 1; 0 0 λ]
eigvals(J)  # [λ, λ, λ]
```

**Note**: Not diagonalizable for $n > 1$.

---

## Structure-Based Patterns

### Block-Diagonal

**Eigenvalues**: $\lambda(B) = \lambda(A_1) \cup \lambda(A_2) \cup \cdots$

```julia
@variables a b c d
B = [a b 0 0; b a 0 0; 0 0 c d; 0 0 d c]
eigvals(B)  # [a+b, a-b, c+d, c-d]
```

Reduces degree-$n$ problem to smaller independent problems.

### Persymmetric

**Definition**: $Q[i,j] = Q[n+1-j, n+1-i]$ (symmetric about anti-diagonal)

Splits into two half-sized eigenvalue problems via transformation with flip matrix.

```julia
@variables a b c d e f
P = [a b c d; b e f c; c f e b; d c b a]
eigvals(P)  # Solves two 2×2 problems
```

---

## Symmetry-Based Patterns

### Circulant Matrices

**Definition**: Each row is cyclic shift of previous row

**Eigenvalues** (any size $n$):
$$\lambda_k = \sum_{j=0}^{n-1} c_j \omega^{jk}, \quad \omega = e^{2\pi i/n}$$

```julia
@variables a b c
C = [a b c; c a b; b c a]
eigvals(C)  # DFT of first row
```

Works for arbitrarily large $n$ - bypasses Abel-Ruffini.

### Block Circulant

Block-level circulant structure. Eigenvalues from $n$ problems of size $k$:
$$\lambda(BC) = \bigcup_{m=0}^{n-1} \lambda\left(\sum_j A_j \omega^{mj}\right)$$

Requires block size $k \leq 4$ for closed form.

### Anti-Diagonal

**Definition**: Non-zero only on anti-diagonal with symmetry

**Eigenvalues**: Come in $\pm$ pairs (center unpaired for odd $n$)

```julia
@variables a b c
A = [0 0 0 0 a; 0 0 0 b 0; 0 0 c 0 0; 0 b 0 0 0; a 0 0 0 0]
eigvals(A)  # {c, ±a, ±b}
```

---

## Tensor Product Patterns

### Kronecker Products

**Eigenvalues**: All products of factor eigenvalues
$$\lambda(A \otimes B) = \{\lambda_i(A) \cdot \lambda_j(B)\}$$

```julia
@variables a b c d
A = [a 0; 0 b]
B = [c 0; 0 d]
K = kron(A, B)
eigvals(K)  # [ac, ad, bc, bd]
```

Reduces $(mn) \times (mn)$ to $m \times m$ and $n \times n$ problems.

### Nested Kronecker Products

**Eigenvalues**: For $A_1 \otimes A_2 \otimes \cdots \otimes A_k$, all products of factor eigenvalues:
$$\lambda(A_1 \otimes \cdots \otimes A_k) = \{\lambda_{i_1}(A_1) \cdot \lambda_{i_2}(A_2) \cdots \lambda_{i_k}(A_k)\}$$

```julia
using Symbolics

# Create 5 independent 2×2 symmetric matrices
matrices = Matrix{Num}[]
for i in 1:5
    a = Symbolics.variable(Symbol("a$i"))
    b = Symbolics.variable(Symbol("b$i"))
    c = Symbolics.variable(Symbol("c$i"))
    push!(matrices, [a b; b c])
end

# Build 32×32 matrix with 15 parameters
K = reduce(kron, matrices)
eigvals(K)  # 32 eigenvalues as products
```

**Performance**: Recursive detection and solving handles arbitrary depth efficiently:

| Depth | Matrix Size | Parameters | Eigenvalues | Time |
|-------|-------------|------------|-------------|------|
| 3 | 8×8 | 9 | 8 | ~12s |
| 5 | 32×32 | 15 | 32 | ~12s |
| 7 | 128×128 | 21 | 128 | ~23s |
| 10 | 1024×1024 | 30 | 1024 | ~33s |

The algorithm recursively factors $K = A \otimes B$, solves each factor independently,
then forms all pairwise products. This bypasses the exponential complexity of direct
symbolic computation on the full matrix.

### SO(2) Kronecker Products

**Definition**: Kronecker product of SO(2) rotation matrices
$$R(\theta) = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}$$

**Eigenvalues** for $R(\theta) \otimes R(\phi)$:
$$\lambda = e^{\pm i(\theta \pm \phi)} = \cos(\theta \pm \phi) \pm i \sin(\theta \pm \phi)$$

The four eigenvalues are:
- $\cos(\theta + \phi) + i\sin(\theta + \phi)$
- $\cos(\theta + \phi) - i\sin(\theta + \phi)$
- $\cos(\theta - \phi) + i\sin(\theta - \phi)$
- $\cos(\theta - \phi) - i\sin(\theta - \phi)$

**Convenience constructors**: Use `R2(θ)` to create rotation matrices and `so2_kron(angles)` for Kronecker products:

```julia
@variables α β γ

# 2-fold: 4×4 matrix
K2 = so2_kron([α, β])
eigvals(K2)
# [cos(-α - β) + im*sin(-α - β), cos(α - β) + im*sin(α - β),
#  cos(-α + β) + im*sin(-α + β), cos(α + β) + im*sin(α + β)]

# 3-fold: 8×8 matrix with 8 clean eigenvalues
K3 = so2_kron([α, β, γ])
eigvals(K3)  # cos(±α ± β ± γ) + i·sin(±α ± β ± γ) for all 8 sign combinations
```

**Direct eigenvalue computation**: For faster computation without building the matrix:

```julia
so2_kron_eigenvalues([α, β, γ])  # Returns 8 eigenvalues directly
```

**Same-angle case**: $R(\theta) \otimes R(\theta)$ uses double-angle formulas:

```julia
eigvals(kron(R2(θ), R2(θ)))
# [cos(2θ) + im*sin(2θ), 1, 1, cos(2θ) - im*sin(2θ)]
```

**Implementation**: Uses direct extraction of (cos, sin) pairs from matrix entries,
avoiding the `sqrt(cos²)` problem that produces messy symbolic expressions.

**Nesting**: Works recursively for $R(\theta_1) \otimes R(\theta_2) \otimes \cdots \otimes R(\theta_k)$
with $2^k$ clean eigenvalues.

---

## Lie Group Patterns

### SO(n) Rotation Matrices

**Definition**: Orthogonal matrices with determinant +1 (special orthogonal group)

**Eigenvalue Structure**:
- All eigenvalues lie on the unit circle: $|\lambda| = 1$
- They come in conjugate pairs: $e^{\pm i\theta_j}$
- Odd dimensions have eigenvalue 1 (rotation axis)

**Symbolic closed-form solutions** for $n \leq 4$:

| Dimension | Eigenvalues | Notes |
|-----------|-------------|-------|
| SO(2) | $\cos\theta \pm i\sin\theta$ | Direct extraction |
| SO(3) | $1, \cos\theta \pm i\sin\theta$ | Trace-based |
| SO(4) | $\cos\theta_1 \pm i\sin\theta_1, \cos\theta_2 \pm i\sin\theta_2$ | Block-diagonal preferred |

**Note**: SO(5)+ are not supported symbolically. For numeric SO(n) matrices, use `LinearAlgebra.eigvals()`.

**Block-diagonal detection**: For SO(4) matrices that are block-diagonal (two independent SO(2) rotations), the package prioritizes block-diagonal detection to produce the cleanest eigenvalues.

```julia
# Symbolic SO(3) rotation
@variables θ
R3 = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
eigvals(R3)  # [1, cos(θ) + im*sin(θ), cos(θ) - im*sin(θ)]

# Block-diagonal SO(4) - cleanest output
@variables θ φ
R4 = [cos(θ) -sin(θ) 0 0; sin(θ) cos(θ) 0 0; 0 0 cos(φ) -sin(φ); 0 0 sin(φ) cos(φ)]
eigvals(R4)  # [cos(θ) + im*sin(θ), cos(θ) - im*sin(θ), cos(φ) + im*sin(φ), cos(φ) - im*sin(φ)]
```

### Other Supported Lie Groups

| Group | Dimension | Description | Eigenvalue Property |
|-------|-----------|-------------|---------------------|
| SU(2) | 2×2 | Special unitary | $e^{\pm i\theta}$ on unit circle |
| SU(3) | 3×3 | Special unitary | Product = 1 |
| Sp(2) | 2×2 | Symplectic (≅ SL(2)) | Reciprocal pairs |
| Sp(4) | 4×4 | Symplectic | Reciprocal pairs |

---

## Tridiagonal Patterns

### Symmetric Toeplitz Tridiagonal

**Definition**: Constant diagonals $a$ (main), $b$ (off-diagonal)

**Eigenvalues** (any size $n$):
$$\lambda_k = a + 2b\cos\left(\frac{k\pi}{n+1}\right), \quad k = 1, \ldots, n$$

**Eigenvectors**: $v_k(j) = \sin(jk\pi/(n+1))$

```julia
@variables a b
T = [a b 0 0; b a b 0; 0 b a b; 0 0 b a]
eigvals(T)
```

### Special 5×5 Patterns

Two patterns with closed-form solutions discovered empirically:

**Pattern [b,d,b,b]**:
```
[a b 0 0 0; b a d 0 0; 0 d a b 0; 0 0 b a b; 0 0 0 b a]
```

**Pattern [b,b,d,b]**:
```
[a b 0 0 0; b a b 0 0; 0 b a d 0; 0 0 d a b; 0 0 0 b a]
```

**Both have identical eigenvalues**: $\{a \pm \sqrt{2b^2 + d^2}, a \pm b, a\}$

---

## Permutation Patterns

### Permutation Matrices

**Definition**: Exactly one 1 per row/column, rest zeros

**Eigenvalues**: Roots of unity from cycle decomposition
- $k$-cycle contributes: $\{1, \omega, \omega^2, \ldots, \omega^{k-1}\}$ where $\omega = e^{2\pi i/k}$

```julia
# 3-cycle (1→2→3→1)
P = [0 1 0; 0 0 1; 1 0 0]
eigvals(P)  # {1, ω, ω²} where ω = e^{2πi/3}
```

---

## Graph Patterns

### Dihedral Symmetry

Matrices commuting with rotation and reflection. Eigenspaces decompose by representation theory, often reducing to quadratic equations.

### Strongly Regular Graphs

Parameters $(n, k, \lambda, \mu)$ determine exactly 3 distinct eigenvalues:
- $\lambda_1 = k$
- $\lambda_{2,3} = \frac{(\lambda-\mu) \pm \sqrt{(\lambda-\mu)^2 + 4(k-\mu)}}{2}$

**Example**: Petersen graph $(10,3,0,1)$ → eigenvalues $\{3, 1^{(5)}, -2^{(4)}\}$

### Hypercube $Q_n$

**Eigenvalues**: $\lambda_k = n - 2k$ with multiplicity $\binom{n}{k}$

```julia
# Q_3 (8×8 cube adjacency)
eigvals(Q3)  # [3, 1, 1, 1, -1, -1, -1, -3]
```

### BCCB (Block-Circulant with Circulant Blocks)

Diagonalized by 2D DFT. Common in image processing and 2D convolution.

---

## Detection Notes

| Pattern | Detection Method |
|---------|-----------------|
| Diagonal | All off-diagonal zero |
| Triangular | Zeros above/below diagonal |
| Block-diagonal | Connected components of non-zero structure |
| Circulant | Row $i$ = cyclic shift of row $i-1$ |
| Toeplitz tridiag | Constant diagonals, symmetric |
| Kronecker | Block scaling pattern |
| SO(2) Kronecker | Block structure + orthogonality + trig pattern |
| SO(n) | $A^T A = I$ and $\det(A) = 1$ |
| Permutation | Exactly one 1 per row/column |

---

## Adding New Patterns

1. Verify numerically (multiple test cases)
2. Verify symbolically (characteristic polynomial)
3. Explain mathematical basis
4. Implement detector and solver
5. Add tests and documentation

See [Contributing](contributing.md) for details.
