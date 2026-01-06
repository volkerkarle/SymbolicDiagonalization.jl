# Pattern Library

SymbolicDiagonalization.jl automatically detects and exploits matrix structure to find closed-form eigenvalues. This page documents all supported patterns.

## Structural Patterns

These patterns are detected in any matrix regardless of specific form.

### Block-Diagonal

Matrices with independent blocks are solved recursively:

```@example patterns
using SymbolicDiagonalization, Symbolics, LinearAlgebra
using Main: LaTeX

@variables a b c d

B = [a  b  0  0;
     b  a  0  0;
     0  0  c  d;
     0  0  d  c]
nothing # hide
```

```@example patterns
LaTeX(B)
```

**Eigenvalues (each block solved independently):**

```@example patterns
LaTeX(eigvals(B))
```

### Persymmetric (Symmetric about Anti-Diagonal)

Persymmetric matrices decompose into two half-size problems:

```@example patterns
# A persymmetric matrix satisfies J*A*J = A where J is the exchange matrix
@variables a b c

P = [a  b  c;
     b  a  b;
     c  b  a]  # Symmetric circulant is also persymmetric
nothing # hide
```

```@example patterns
LaTeX(P)
```

```@example patterns
LaTeX(eigvals(P))
```

---

## Finite Group Patterns

Matrices with group-theoretic structure.

### Circulant Matrices (Cyclic Group Zₙ)

Any n×n circulant matrix has eigenvalues via the Discrete Fourier Transform:

```math
\lambda_k = \sum_{j=0}^{n-1} c_j \omega^{jk}, \quad \omega = e^{2\pi i/n}
```

```@example patterns
@variables a b c d e

# 5×5 circulant
C = [a b c d e;
     e a b c d;
     d e a b c;
     c d e a b;
     b c d e a]
nothing # hide
```

```@example patterns
LaTeX(C)
```

**Eigenvalues (via DFT of first row):**

```@example patterns
LaTeX(eigvals(C))
```

### Symmetric Circulant (Dihedral Group Dₙ)

When the first row is palindromic, eigenvalues simplify:

```@example patterns
@variables a b

# Symmetric circulant (first row = [a, b, b])
S = [a b b;
     b a b;
     b b a]
nothing # hide
```

```@example patterns
LaTeX(S)
```

```@example patterns
LaTeX(eigvals(S))
```

### Quaternion Group Q₈

#### Single Quaternion Matrix

```@example patterns
@variables q0 qi qj qk

# 2×2 quaternion representation
Q = [q0 + qi*im  qj + qk*im;
    -qj + qk*im  q0 - qi*im]
nothing # hide
```

```@example patterns
LaTeX(Q)
```

```@example patterns
LaTeX(eigvals(Q))
```

---

## Lie Groups

### SO(2): 2D Rotation Group

```@example patterns
@variables θ
R2 = SO2_rotation(θ)
nothing # hide
```

```@example patterns
LaTeX(R2)
```

**Eigenvalues e^{±iθ}:**

```@example patterns
LaTeX(eigvals(R2))
```

### SO(3): 3D Rotation Group

Axis-aligned rotations:

```@example patterns
Rx = SO3_Rx(θ)
nothing # hide
```

```@example patterns
LaTeX(Rx)
```

**Eigenvalues {1, e^{±iθ}}:**

```@example patterns
LaTeX(eigvals(Rx))
```

### SO(4): 4D Rotation Group

Double rotation structure:

```@example patterns
@variables θ φ

R4 = [cos(θ) -sin(θ) 0 0;
      sin(θ)  cos(θ) 0 0;
      0 0 cos(φ) -sin(φ);
      0 0 sin(φ)  cos(φ)]
nothing # hide
```

```@example patterns
LaTeX(R4)
```

**Eigenvalues {e^{±iθ}, e^{±iφ}}:**

```@example patterns
LaTeX(eigvals(R4))
```

### SU(2): Special Unitary Group

```@example patterns
Uz = SU2_Uz(θ)
nothing # hide
```

```@example patterns
LaTeX(Uz)
```

**Eigenvalues e^{±iθ/2}:**

```@example patterns
LaTeX(eigvals(Uz))
```

---

## Kronecker Products

### General Kronecker Products

Eigenvalues of A ⊗ B are products of individual eigenvalues. A 4×4 matrix (degree-4 polynomial) reduces to two quadratics:

```@example patterns
@variables a b c d

A = [a b; b a]
B = [c d; d c]
K = kron(A, B)
nothing # hide
```

```@example patterns
LaTeX(K)
```

**Eigenvalues (products of factor eigenvalues):**

```@example patterns
LaTeX(eigvals(K))
```

### SO(2) Kronecker Products

Trigonometric simplification gives clean angle-sum/difference forms:

```@example patterns
@variables θ φ

K = kron(SO2_rotation(θ), SO2_rotation(φ))
vals = eigvals(K)
nothing # hide
```

**Eigenvalues e^{i(±θ±φ)}:**

```@example patterns
LaTeX(vals)
```

### SU(2) Kronecker Products

Half-angle formulas for spin-1/2 representations:

```@example patterns
@variables α β
K = SU2_kron([α, β])
vals = eigvals(K)
nothing # hide
```

**Eigenvalues with half-angles:**

```@example patterns
LaTeX(vals)
```

### Nested Kronecker Products

Arbitrary depth Kronecker products:

```julia
# 5-fold Kronecker product: 32×32 with 15 parameters
@variables a1 b1 c1  a2 b2 c2  a3 b3 c3  a4 b4 c4  a5 b5 c5

matrices = [[a1 b1; b1 c1], [a2 b2; b2 c2], [a3 b3; b3 c3],
            [a4 b4; b4 c4], [a5 b5; b5 c5]]
K = reduce(kron, matrices)  # 32×32
eigvals(K)  # 32 symbolic eigenvalues
```

---

## Tridiagonal Patterns

### Symmetric Toeplitz Tridiagonal

Constant diagonals with closed-form eigenvalues:

```@example patterns
@variables a b

# Toeplitz tridiagonal [b, a, b]
T = [a b 0 0;
     b a b 0;
     0 b a b;
     0 0 b a]
nothing # hide
```

```@example patterns
LaTeX(T)
```

**Eigenvalues λₖ = a + 2b·cos(kπ/(n+1)):**

```@example patterns
LaTeX(eigvals(T))
```

---

## Pattern Detection Priority

When multiple patterns apply, the package uses this priority:

1. **Diagonal** - Trivial eigenvalues
2. **Block-diagonal** - Recursive decomposition
3. **Lie groups** - SO(2), SO(3), SO(4), SU(2), SU(3)
4. **Kronecker products** - Detect and factor
5. **Circulant** - DFT formula
6. **Tridiagonal** - Chebyshev formulas
7. **Persymmetric** - Half-size reduction
8. **General** - Cardano/Ferrari formulas (up to 4×4)

## Adding Custom Patterns

See [Mathematical Background](mathematical_background.md) for extending the pattern library.
