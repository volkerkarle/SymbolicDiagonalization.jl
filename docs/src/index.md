# SymbolicDiagonalization.jl

Symbolic matrix diagonalization for Julia using [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl).

## Features

- Closed-form eigenvalue solvers for degrees 1-4 (linear, quadratic, Cardano, Ferrari)
- Automatic structure detection for larger matrices (16+ patterns)
- Lie group detection (SO(2)-SO(4), SU(2), SU(3), Sp(2), Sp(4)) with symbolic eigenvalues
- Aggressive trigonometric simplification (sqrt(1-cos²θ) → sin(θ))
- SO(2) Kronecker products with automatic trig simplification
- Diagonal shift optimization for full 6-parameter 3×3 symmetric matrices
- Nested Kronecker products (A₁ ⊗ A₂ ⊗ ... ⊗ Aₙ) - scales to 1024×1024 with 30 parameters
- Drop-in `LinearAlgebra` interface: `eigen()`, `eigvals()`

## Installation

```julia
using Pkg
Pkg.add("SymbolicDiagonalization")
```

## Quick Start

```@example main
using SymbolicDiagonalization
using Symbolics
using LinearAlgebra

@variables a b c
mat = [a 1 0; 0 b 1; 0 0 c]

E = eigen(mat)        # Eigen object with .values and .vectors
E.values
```

```@example main
eigvals(mat)          # eigenvalues only
```

## Highlights

### Full 6-Parameter 3×3 Symmetric Matrices

```julia
@variables a b c d e f
A = [a b c; b d e; c e f]  # 6 independent parameters
eigvals(A)  # Works via diagonal shift optimization (~107s)
```

### Nested Kronecker Products

```julia
# 10-fold Kronecker product: 1024×1024 matrix with 30 parameters!
matrices = [[Symbolics.variable(Symbol("a$i")) Symbolics.variable(Symbol("b$i"));
             Symbolics.variable(Symbol("b$i")) Symbolics.variable(Symbol("c$i"))] for i in 1:10]
K = reduce(kron, matrices)
eigvals(K)  # 1024 symbolic eigenvalues in ~33 seconds
```

### SO(2) Rotation Kronecker Products

When Kronecker products of SO(2) rotation matrices are detected, eigenvalues are 
automatically simplified to clean trigonometric form:

```julia
@variables θ φ
R(x) = [cos(x) -sin(x); sin(x) cos(x)]

eigvals(kron(R(θ), R(φ)))
# [cos(θ+φ) + im*sin(θ+φ), cos(θ-φ) + im*sin(θ-φ),
#  cos(θ-φ) - im*sin(θ-φ), cos(θ+φ) - im*sin(θ+φ)]

eigvals(kron(R(θ), R(θ)))  # Same-angle case
# [cos(2θ) + im*sin(2θ), 1, 1, cos(2θ) - im*sin(2θ)]
```

## Supported Patterns

For matrices larger than 4×4, the package detects and exploits special structure:

| Pattern | Description |
|---------|-------------|
| Diagonal | Direct extraction |
| Triangular | Diagonal elements |
| Block-diagonal | Recursive solving |
| Circulant | DFT-based (any size n) |
| Symmetric Toeplitz tridiagonal | Cosine formula (any size n) |
| Kronecker products A ⊗ B | Product of factor eigenvalues |
| Nested Kronecker A₁ ⊗ A₂ ⊗ ... ⊗ Aₙ | Recursive decomposition (any depth) |
| SO(2) Kronecker R(θ) ⊗ R(φ) | Clean trig form: cos(θ±φ) + i·sin(θ±φ) |
| Permutation | Roots of unity from cycle structure |
| Persymmetric | Half-size decomposition |

See [Pattern Library](pattern_library.md) for the complete list with examples.

## Performance

| Matrix | Size | Parameters | Time |
|--------|------|------------|------|
| 3×3 symmetric (full) | 3×3 | 6 | ~107s |
| 3×3 ⊗ 2×2 Kronecker | 6×6 | 9 | ~165s |
| 2×2^⊗5 nested Kronecker | 32×32 | 15 | ~12s |
| 2×2^⊗10 nested Kronecker | 1024×1024 | 30 | ~33s |

## Limitations

- **General 5×5+ matrices**: No closed-form solution exists (Abel-Ruffini theorem) - requires exploitable structure
- **Expression size**: Fully symbolic 4×4 quartics can produce large expressions
- **Simplification**: Results may not be in minimal form

## Documentation

**Getting Started**
- [User Guide](user_guide.md) - Installation, workflows, troubleshooting
- [API Reference](api_reference.md) - Complete function signatures

**Patterns & Theory**
- [Pattern Library](pattern_library.md) - All 16+ patterns with examples
- [Mathematical Background](mathematical_background.md) - Theory and algorithms
- [Group Theory Examples](group_theory_examples.md) - Symmetry-based approaches

**Development**
- [Implementation Details](implementation.md) - Internals and algorithms
- [Contributing](contributing.md) - Development setup and guidelines
