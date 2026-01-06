# SymbolicDiagonalization.jl

Symbolic matrix diagonalization for Julia.

## The Problem

The Abel-Ruffini theorem proves no closed-form solution exists for polynomials of degree 5+. General symbolic 5×5+ matrices cannot be diagonalized. However, matrices with exploitable structure can be solved regardless of size.

## Usage

```@example main
using SymbolicDiagonalization, Symbolics, LinearAlgebra
using Main: LaTeX

@variables a b c d
M = [a b 0 0; b a 0 0; 0 0 c d; 0 0 d c]
nothing # hide
```

**Block-diagonal matrix:**

```@example main
LaTeX(M)
```

**Eigenvalues (each block solved independently):**

```@example main
LaTeX(eigvals(M))
```

## Supported Patterns

### Structural

| Pattern | Method |
|---------|--------|
| Block-diagonal | Recursive decomposition |
| Persymmetric | Half-size reduction |

### Finite Groups

| Pattern | Method |
|---------|--------|
| Circulant (Zₙ) | DFT diagonalization |
| Symmetric Circulant (Dₙ) | Real DFT |
| Permutation (Sₙ) | Cycle decomposition |
| Quaternion Q₈ | Character theory |

### Kronecker Products

| Pattern | Method |
|---------|--------|
| A ⊗ B | λ(A)·λ(B) |
| Nested A₁⊗...⊗Aₙ | Recursive factorization |
| SO(2) Kronecker | Angle sum/difference |
| SU(2) Kronecker | Half-angle formulas |

### Transform Matrices

| Pattern | Eigenvalues |
|---------|-------------|
| Hadamard 2ⁿ | ±√(2ⁿ) |
| DFT n×n | Fourth roots of n |

### Lie Groups

| Group | Eigenvalues |
|-------|-------------|
| SO(2) | e^{±iθ} |
| SO(3) | 1, e^{±iθ} |
| SO(4) | e^{±iθ₁}, e^{±iθ₂} |
| SU(2) | e^{±iθ/2} |
| SU(3) | Cubic formula |

### Tridiagonal

| Pattern | Formula |
|---------|---------|
| Symmetric Toeplitz | λₖ = a + 2b·cos(kπ/(n+1)) |
| Path Laplacian | λₖ = 2 - 2cos(πk/n) |
| Cycle Laplacian | λₖ = 2 - 2cos(2πk/n) |
| Cartan Aₙ | λₖ = 4sin²(πk/2(n+1)) |

## Examples

### Circulant (any size)

```@example main
@variables a b c d
C = [a b c d; d a b c; c d a b; b c d a]
nothing # hide
```

```@example main
LaTeX(C)
```

**Eigenvalues via DFT:**

```@example main
LaTeX(eigvals(C))
```

### Kronecker Product

```@example main
@variables θ φ
K = kron(SO2_rotation(θ), SO2_rotation(φ))
nothing # hide
```

**Eigenvalues e^{i(±θ±φ)}:**

```@example main
LaTeX(eigvals(K))
```

### Nested Kronecker (1024×1024)

```julia
# 10-fold Kronecker: 1024×1024 matrix, 30 parameters
# Degree-1024 polynomial reduced to 10 quadratics
matrices = [[Symbolics.variable(Symbol("a$i")) Symbolics.variable(Symbol("b$i"));
             Symbolics.variable(Symbol("b$i")) Symbolics.variable(Symbol("c$i"))] for i in 1:10]
K = reduce(kron, matrices)
eigvals(K)  # ~33 seconds
```

### Path Laplacian

```@example main
L = path_laplacian(5)
nothing # hide
```

```@example main
LaTeX(L)
```

**Eigenvalues λₖ = 2 - 2cos(πk/n):**

```@example main
LaTeX(eigvals(L))
```

## API

### Main Functions

```julia
eigvals(M)                    # Eigenvalues (LinearAlgebra interface)
eigen(M)                      # Eigen decomposition
symbolic_eigenvalues(M)       # Returns (values, poly, λ)
symbolic_eigenpairs(M)        # Returns (pairs, poly, λ)
symbolic_diagonalize(M)       # Returns (P, D, pairs) where M = P*D*P⁻¹
```

### Constructors

```julia
# Rotations
SO2_rotation(θ)
SO3_Rx(θ), SO3_Ry(θ), SO3_Rz(θ)
SU2_Ux(θ), SU2_Uy(θ), SU2_Uz(θ)

# Kronecker
SO2_kron([θ, φ, ...])
SU2_kron([θ, φ, ...])

# Graphs
path_laplacian(n)
cycle_laplacian(n)

# Transforms
hadamard_matrix(n)  # 2ⁿ×2ⁿ
dft_matrix(n)

# Cartan
cartan_matrix_A(n)
cartan_matrix_G2()
```

## Limitations

- **5×5+ without structure**: No solution exists (Abel-Ruffini)
- **Large expressions**: 4×4 quartic formulas produce large expressions

## Source

[github.com/volkerkarle/SymbolicDiagonalization.jl](https://github.com/volkerkarle/SymbolicDiagonalization.jl)
