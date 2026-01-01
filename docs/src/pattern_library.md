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
| Permutation | Exactly one 1 per row/column |

---

## Adding New Patterns

1. Verify numerically (multiple test cases)
2. Verify symbolically (characteristic polynomial)
3. Explain mathematical basis
4. Implement detector and solver
5. Add tests and documentation

See [Contributing](contributing.md) for details.
