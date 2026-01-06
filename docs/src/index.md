# SymbolicDiagonalization.jl

Symbolic matrix diagonalization for Julia.

## The Problem

The eigenvalue problem requires solving the characteristic polynomial. For degree ≥ 5, no general formula exists (Abel-Ruffini theorem). This package finds closed-form symbolic eigenvalues by exploiting matrix structure.

## Key Feature: Truly Symbolic Eigenvalues

Unlike numerical solvers that return floating-point approximations, SymbolicDiagonalization returns **exact symbolic expressions**:

| Matrix Type | Numerical Solver | SymbolicDiagonalization |
|-------------|-----------------|-------------------------|
| Circulant 6×6 | `-3.0 - 5.196...im` | `-3 - 3√3·im` |
| Hadamard 4×4 | `2.0, -2.0` | `2, -2` (exact integers) |
| Hadamard 2×2 | `1.414...` | `√2` (symbolic) |
| DFT 4×4 | `2.0, -2.0, 2.0im` | `2, -2, 2im` (exact) |
| Path Laplacian | `0.381...` | `(3/2) - (1/2)√5` |

This preserves mathematical structure for downstream symbolic computation.

### Special Angle Simplification

Eigenvalues at special angles (multiples of π/6, π/5, π/4, π/3, π/2) are automatically simplified to their exact algebraic values:

| Angle | cos value | sin value |
|-------|-----------|-----------|
| π/6 | √3/2 | 1/2 |
| π/5 | (1+√5)/4 | √(10-2√5)/4 |
| π/4 | √2/2 | √2/2 |
| π/3 | 1/2 | √3/2 |
| 2π/5 | (√5-1)/4 | √(10+2√5)/4 |
| π/2 | 0 | 1 |

For example, path Laplacian eigenvalues `2 - 2cos(kπ/5)` simplify to golden ratio expressions:
- `2 - 2cos(π/5) → (3/2) - (1/2)√5`
- `2 - 2cos(2π/5) → (5/2) - (1/2)√5`

## Installation

```julia
using Pkg
Pkg.add("SymbolicDiagonalization")
```

## Usage

```@example main
using SymbolicDiagonalization, Symbolics, LinearAlgebra
using Main: LaTeX  # Import LaTeX wrapper from make.jl

@variables a b c
M = [a b; b c]
nothing # hide
```

**Input matrix:**

```@example main
LaTeX(M)
```

**Eigenvalues:**

```@example main
LaTeX(eigvals(M))
```

**Eigenvectors:**

```@example main
LaTeX(eigen(M).vectors)
```

## Supported Patterns

Matrices larger than 4×4 require exploitable structure. The package detects:

### Finite Groups

| Pattern | Group | Example | Eigenvalue Formula |
|---------|-------|---------|-------------------|
| Circulant | Zₙ | `[a b c; c a b; b c a]` | DFT of first row (exact roots of unity: `1, -1, ±i, ±√3/2`, etc.) |
| Symmetric Circulant | Dₙ | Palindromic first row | λₖ = c₀ + 2·Σcⱼcos(2πjk/n) |
| Permutation | Sₙ | `[0 1 0; 0 0 1; 1 0 0]` | Roots of unity from cycles |
| Quaternion (single) | Q₈ | 2×2 quaternion blocks | λ = a ± i√(b²+c²+d²) |
| Q₈ Regular Rep | Q₈ | 8×8 Q₈-invariant | Character theory: 5 distinct eigenvalues |
| Hypercube Qₙ | (Z₂)ⁿ | 2ⁿ×2ⁿ adjacency | λₖ = n - 2k (exact integers) |

### Transform Matrices

| Pattern | Size | Constructor | Eigenvalue Formula |
|---------|------|-------------|-------------------|
| Hadamard | 2ⁿ×2ⁿ | `hadamard_matrix(n)` | λ = ±2^(n/2) (exact: integers for even n, √2 multiples for odd n) |
| DFT (Fourier) | n×n | `dft_matrix(n)` | λ ∈ {±√n, ±i√n} (exact integers for perfect square n) |

### Coxeter/Weyl Groups

Only types with closed-form symbolic eigenvalue formulas are automatically diagonalized.
Other types have matrix constructors available for reference.

| Type | Lie Group | Constructor | Symbolic Eigenvalues |
|------|-----------|-------------|---------------------|
| Aₙ | SU(n+1) | `cartan_matrix_A(n)` | ✓ λₖ = 4sin²(πk/2(n+1)) |
| G₂ | Exceptional | `cartan_matrix_G2()` | ✓ λ = 2 ± √3 |
| Bₙ | SO(2n+1) | `cartan_matrix_B(n)` | Constructor only |
| Cₙ | Sp(2n) | `cartan_matrix_C(n)` | Constructor only |
| Dₙ | SO(2n) | `cartan_matrix_D(n)` | Constructor only |
| E₆,₇,₈ | Exceptional | `cartan_matrix_E(n)` | Constructor only |
| F₄ | Exceptional | `cartan_matrix_F4()` | Constructor only |

### Graph Laplacians

| Graph | Constructor | Eigenvalue Formula |
|-------|-------------|-------------------|
| Path Pₙ | `path_laplacian(n)` | λₖ = 2 - 2cos(πk/n) |
| Cycle Cₙ | `cycle_laplacian(n)` | λₖ = 2 - 2cos(2πk/n) |

### Lie Groups

| Pattern | Eigenvalues | Eigenvectors |
|---------|-------------|--------------|
| SO(2) | e^{±iθ} | [1, ±i] (fixed) |
| SO(3) | 1, e^{±iθ} | Axis + rotation plane |
| SO(4) | e^{±iθ₁}, e^{±iθ₂} | Projection operators |
| SU(2) | e^{±iθ/2} | From Pauli structure |
| SU(3) | Cubic formula | - |
| Sp(2), Sp(4) | Reciprocal pairs | - |

### Tensor Products

| Pattern | Formula |
|---------|---------|
| A ⊗ B | λᵢ(A) · λⱼ(B) |
| SO(2)⊗ᵏ | e^{i(±θ₁±θ₂±...)} |
| SU(2)⊗ᵏ | e^{i(±θ₁±θ₂±...)/2} |

### Structural

| Pattern | Method |
|---------|--------|
| Block-diagonal | Solve blocks independently |
| Sym. Toeplitz tridiagonal | λₖ = a + 2b·cos(kπ/(n+1)) |
| Persymmetric | Half-size decomposition |

## Examples

### Circulant (any size)

```@example main
@variables a b c d
C = [a b c d; d a b c; c d a b; b c d a]
nothing # hide
```

**Input matrix:**

```@example main
LaTeX(C)
```

**Eigenvalues (via DFT):**

```@example main
LaTeX(eigvals(C))
```

### SO(2) Rotation

```@example main
@variables θ
R = SO2_rotation(θ)
nothing # hide
```

**Rotation matrix:**

```@example main
LaTeX(R)
```

**Eigenvalues** $e^{\pm i\theta}$:

```@example main
LaTeX(eigvals(R))
```

### SO(2) Kronecker Product

```@example main
@variables θ φ
K = kron(SO2_rotation(θ), SO2_rotation(φ))
nothing # hide
```

**Kronecker product eigenvalues** $e^{i(\pm\theta \pm \phi)}$:

```@example main
LaTeX(eigvals(K))
```

### Quaternion Group Q₈ Regular Representation

```@example main
# Q₈-invariant 8×8 matrix with 8 coefficients (one per group element)
M = Q8_invariant_matrix(1.0, 0.5, 0.3, 0.2, 0.4, 0.1, 0.25, 0.15)
nothing # hide
```

**Eigenvalues (5 distinct via character theory):**

```@example main
LaTeX(eigvals(M))
```

### Nested Kronecker (1024×1024)

```julia
matrices = [[Symbolics.variable(Symbol("a\$i")) Symbolics.variable(Symbol("b\$i"));
             Symbolics.variable(Symbol("b\$i")) Symbolics.variable(Symbol("c\$i"))] for i in 1:10]
K = reduce(kron, matrices)  # 1024×1024, 30 parameters
eigvals(K)  # ~33 seconds
```

## Breaking the Abel-Ruffini Barrier

The Abel-Ruffini theorem states that polynomials of degree 5+ have no general closed-form solution. Yet these matrices---which would require solving degree-10, degree-16, or even degree-1024 polynomials---yield exact symbolic eigenvalues instantly:

### 10×10 Circulant Matrix

```julia
# Degree-10 characteristic polynomial - no closed form in general!
n = 10
first_row = collect(1:n)
C = [first_row[mod(j - i, n) + 1] for i in 1:n, j in 1:n]
eigvals(C)  # Instant via DFT formula
# Works for 100×100, 1000×1000, any size!
```

The Discrete Fourier Transform diagonalizes *all* circulant matrices, regardless of size.

### 16×16 Hypercube Graph Q₄

```julia
# 4-dimensional hypercube: 16 vertices, each connected to 4 neighbors
Q4 = zeros(Int, 16, 16)
for i in 0:15, j in 0:15
    count_ones(xor(i, j)) == 1 && (Q4[i+1, j+1] = 1)
end
eigvals(Q4)  # λₖ = 4-2k with multiplicity C(4,k)
# Works for Q₆ (64×64), Q₁₀ (1024×1024), any dimension!
```

The Walsh-Hadamard basis diagonalizes all hypercube graphs. Eigenvalues: $\lambda_k = n - 2k$ with multiplicity $\binom{n}{k}$.

### 20×20 Path Laplacian

```julia
L = path_laplacian(20)  # Tridiagonal [1,-1,0,...], [-1,2,-1,...], ...
eigvals(L)              # λₖ = 2 - 2cos(πk/20) for k=0,...,19
# Works for 50×50, 100×100, any size!
```

Chebyshev polynomial theory gives closed-form eigenvalues for all symmetric Toeplitz tridiagonal matrices.

### 16×16 Hadamard Matrix

```julia
H = hadamard_matrix(4)  # 2⁴ = 16×16 matrix of ±1 entries
eigvals(H)              # Only 2 distinct eigenvalues: ±4
# Works for H₃₂, H₆₄, H₁₀₂₄, any power of 2!
```

The self-similar structure of Sylvester-Hadamard matrices means a degree-16 polynomial factors completely!

### Double Rotation Kronecker (Symbolic)

```julia
@variables θ₁ θ₂
R = kron(SO2_rotation(θ₁), SO2_rotation(θ₂))
eigvals(R)  # 4 eigenvalues: e^{i(±θ₁±θ₂)}
```

Kronecker products reduce to smaller subproblems: $\lambda(A \otimes B) = \lambda(A) \cdot \lambda(B)$.

### Why Structure Beats Algebra

| Matrix | Size | Polynomial Degree | Why It Works |
|--------|------|------------------|--------------|
| Circulant | 10×10+ | 10+ | DFT diagonalizes all Zₙ-invariant matrices |
| Hypercube Q₄ | 16×16 | 16 | Walsh-Hadamard basis from (Z₂)⁴ symmetry |
| Hadamard H₁₆ | 16×16 | 16 | Only 2 distinct eigenvalues (±4) |
| Path Laplacian | 20×20 | 20 | Chebyshev polynomials → trigonometric formula |
| Cartan A₇ | 7×7 | 7 | Root system theory → sine formula |
| Kronecker 2^⊗10 | 1024×1024 | 1024 | Factors into 10 quadratic problems |

## API

### Main Functions

```julia
eigvals(M)                    # Eigenvalues (LinearAlgebra interface)
eigen(M)                      # Eigen decomposition
symbolic_eigenvalues(M)       # Returns (values, poly, λ)
symbolic_eigenpairs(M)        # Returns (pairs, poly, λ) where pairs = [(λ₁, [v₁]), ...]
symbolic_diagonalize(M)       # Returns (P, D, pairs) where M = P*D*P⁻¹
```

### Options

```julia
eigvals(M; 
    timeout=300,              # Seconds before timeout
    max_terms=10000,          # Expression complexity limit
    structure=:auto           # :auto, :symmetric, :hermitian, :unitary
)
```

### Lie Group Constructors

```julia
# SO(2)
SO2_rotation(θ)               # 2×2 rotation matrix
SO2_kron([θ, φ, ...])         # Kronecker product of rotations

# SO(3)
SO3_Rx(θ), SO3_Ry(θ), SO3_Rz(θ)

# SU(2)
SU2_Ux(θ), SU2_Uy(θ), SU2_Uz(θ)
SU2_kron([θ, φ, ...])
pauli_x(), pauli_y(), pauli_z()

# SU(3)
gellmann_matrices()           # 8 generators

# Coxeter/Weyl Groups
cartan_matrix_A(n)            # Type Aₙ (SU(n+1))
cartan_matrix_B(n)            # Type Bₙ (SO(2n+1))
cartan_matrix_C(n)            # Type Cₙ (Sp(2n))
cartan_matrix_D(n)            # Type Dₙ (SO(2n))
cartan_matrix_E(n)            # Type Eₙ (n ∈ {6,7,8})
cartan_matrix_F4()            # Type F₄
cartan_matrix_G2()            # Type G₂

coxeter_number(:A, n)         # Coxeter number h
coxeter_exponents(:A, n)      # Exponents [m₁, m₂, ...]
coxeter_element_eigenvalues(:A, n)  # e^{2πi·mⱼ/h}

path_laplacian(n)             # Path graph Pₙ Laplacian
cycle_laplacian(n)            # Cycle graph Cₙ Laplacian
householder_reflection(v)     # Householder matrix

# Transform Matrices
hadamard_matrix(n)            # 2ⁿ×2ⁿ Sylvester-Hadamard matrix
dft_matrix(n)                 # n×n DFT matrix (ω = e^{2πi/n})
dft_matrix(n, normalized=true)  # Unitary DFT (F/√n)

# Finite Group Matrices
Q8_invariant_matrix(c1, cm1, ci, cmi, cj, cmj, ck, cmk)  # Q₈ regular rep
```

### Simplification Functions

```julia
# Simplify eigenvalue expressions
simplify_eigenvalue(expr)     # Full simplification (trig + special angles + algebraic)
simplify_eigenvalues(vals)    # Apply to vector of eigenvalues

# Specific simplifications
simplify_special_angles(expr) # cos(π/3)→1/2, cos(2π/5)→(√5-1)/4, etc.
aggressive_simplify(expr)     # Combines all simplification rules
trig_simplify(expr)           # Trigonometric identities (sin²+cos²=1, etc.)

# Supported special angles: 0, π/6, π/5, π/4, π/3, 2π/5, π/2, 3π/5, 2π/3, 3π/4, 4π/5, 5π/6, π
# Both division form (π/n) and multiplication form ((k//n)*π) are handled
```

## Limitations

- **5×5+ without structure**: No solution exists (Abel-Ruffini)
- **Large expressions**: 4×4 quartic produces ~50MB expressions
- **No threading**: Symbolics.jl is not thread-safe

## Source

[github.com/volkerkarle/SymbolicDiagonalization.jl](https://github.com/volkerkarle/SymbolicDiagonalization.jl)
