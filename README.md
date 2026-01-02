# SymbolicDiagonalization.jl

<p align="center">
  <img src="logo.png" alt="SymbolicDiagonalization.jl Logo" width="300"/>
</p>

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://volkerkarle.github.io/SymbolicDiagonalization.jl)
[![CI](https://github.com/volkerkarle/SymbolicDiagonalization.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/volkerkarle/SymbolicDiagonalization.jl/actions/workflows/ci.yml)
[![Documentation Build](https://github.com/volkerkarle/SymbolicDiagonalization.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/volkerkarle/SymbolicDiagonalization.jl/actions/workflows/documentation.yml)
[![codecov](https://codecov.io/gh/volkerkarle/SymbolicDiagonalization.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/volkerkarle/SymbolicDiagonalization.jl)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Julia Version](https://img.shields.io/badge/julia-v1.12-9558b2.svg)](https://julialang.org/)

**Symbolic eigenvalue computation for structured matrices**

A Julia package for symbolic matrix diagonalization using closed-form root solvers and automatic structure detection.

## The Problem

The Abel-Ruffini theorem proves that no general closed-form solution exists for polynomials of degree 5 or higher. This means general 5x5+ matrices can't be solved symbolically. However, many matrices have exploitable structure that allows closed-form solutions regardless of size.

## The Solution

**SymbolicDiagonalization.jl** automatically detects and exploits matrix structure:

- **Closed-form root solvers** for degrees 1-4 (linear, quadratic, Cardano, Ferrari)
- **Automatic structure detection** (block-diagonal, persymmetric, Hermitian)
- **16+ special pattern solvers** for arbitrary-sized matrices (circulant, Kronecker, Toeplitz tridiagonal, permutation, etc.)
- **Diagonal shift optimization** for full 6-parameter 3×3 symmetric matrices
- **Nested Kronecker products** (A₁ ⊗ A₂ ⊗ ... ⊗ Aₙ) - solve 1024×1024 matrices with 30 parameters!
- **Clean LinearAlgebra.jl interface** (`eigen()`, `eigvals()`)

## Quick Start

```julia
using Symbolics, SymbolicDiagonalization, LinearAlgebra

@variables a b c
mat = [a 1 0; 0 b 1; 0 0 c]

E = eigen(mat)        # Eigen object with .values and .vectors
λ = eigvals(mat)      # eigenvalues only (faster)
```

## Installation

```julia
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Examples

### Block-Diagonal Matrices

```julia
@variables a b c d
mat = [a  b  0  0;
       b  a  0  0;
       0  0  c  d;
       0  0  d  c]

eigvals(mat)  # [a+b, a-b, c+d, c-d]
```

### Circulant Matrices (any size)

```julia
@variables a b c
C = [a b c;
     c a b;
     b c a]

eigvals(C)  # Uses DFT: works for any n, even n=100
```

### Kronecker Products

```julia
@variables a b c d
A = [a 0; 0 b]
B = [c 0; 0 d]
M = kron(A, B)

eigvals(M)  # {ac, ad, bc, bd}
```

### Full 6-Parameter 3×3 Symmetric Matrix

```julia
@variables a b c d e f
A = [a b c; b d e; c e f]  # 6 independent parameters

eigvals(A)  # Works! Uses diagonal shift optimization
```

### Nested Kronecker Products (Scalable to 1000+ dimensions)

```julia
# 5-fold Kronecker product: 32×32 matrix with 15 parameters
@variables a1 b1 c1  a2 b2 c2  a3 b3 c3  a4 b4 c4  a5 b5 c5

matrices = [[a1 b1; b1 c1], [a2 b2; b2 c2], [a3 b3; b3 c3],
            [a4 b4; b4 c4], [a5 b5; b5 c5]]
K = reduce(kron, matrices)  # 32×32

eigvals(K)  # 32 symbolic eigenvalues in ~12 seconds!

# Works for 10-fold (1024×1024, 30 parameters) in ~33 seconds
```

## Documentation

**[Read the full documentation](https://volkerkarle.github.io/SymbolicDiagonalization.jl)** for:

- [User Guide](https://volkerkarle.github.io/SymbolicDiagonalization.jl/user_guide/) - Practical examples and workflows
- [API Reference](https://volkerkarle.github.io/SymbolicDiagonalization.jl/api_reference/) - Complete function reference
- [Pattern Library](https://volkerkarle.github.io/SymbolicDiagonalization.jl/pattern_library/) - All 16+ supported matrix patterns
- [Mathematical Background](https://volkerkarle.github.io/SymbolicDiagonalization.jl/mathematical_background/) - Theory and proofs
- [Contributing](https://volkerkarle.github.io/SymbolicDiagonalization.jl/contributing/) - Development guide

## Capabilities

| Feature | Details |
|---------|---------|
| All matrices up to 4×4 | Closed-form solutions via root formulas |
| Full 6-parameter 3×3 symmetric | Diagonal shift optimization |
| Block-diagonal decomposition | Automatic detection and recursion |
| 16+ pattern solvers | Circulant, Kronecker, tridiagonal, permutation, etc. |
| Nested Kronecker A₁⊗A₂⊗...⊗Aₙ | Scales to 1024×1024 with 30 parameters |
| 306 passing tests | Comprehensive test coverage |

## Performance Benchmarks

| Matrix | Size | Parameters | Time |
|--------|------|------------|------|
| 3×3 symmetric (full) | 3×3 | 6 | ~107s |
| 3×3 ⊗ 2×2 Kronecker | 6×6 | 9 | ~165s |
| 2×2^⊗5 nested Kronecker | 32×32 | 15 | ~12s |
| 2×2^⊗10 nested Kronecker | 1024×1024 | 30 | ~33s |

## Testing

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

## Contributing

Contributions welcome! See the [contributing guide](https://volkerkarle.github.io/SymbolicDiagonalization.jl/contributing/) for details.

## License

GPL v3 - see [LICENSE](LICENSE) for details.
