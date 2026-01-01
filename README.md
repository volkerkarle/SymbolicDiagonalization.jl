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

## Documentation

**[Read the full documentation](https://volkerkarle.github.io/SymbolicDiagonalization.jl)** for:

- [User Guide](https://volkerkarle.github.io/SymbolicDiagonalization.jl/user_guide/) - Practical examples and workflows
- [API Reference](https://volkerkarle.github.io/SymbolicDiagonalization.jl/api_reference/) - Complete function reference
- [Pattern Library](https://volkerkarle.github.io/SymbolicDiagonalization.jl/pattern_library/) - All 16+ supported matrix patterns
- [Mathematical Background](https://volkerkarle.github.io/SymbolicDiagonalization.jl/mathematical_background/) - Theory and proofs
- [Contributing](https://volkerkarle.github.io/SymbolicDiagonalization.jl/contributing/) - Development guide

## Status

**Functional prototype with 16+ special patterns**

| Works | Limitations |
|-------|-------------|
| All matrices up to 4x4 | General 5x5+ requires structure |
| Block-diagonal decomposition | Expression explosion for symbolic 4x4 |
| 16+ pattern solvers (circulant, Kronecker, etc.) | Basic structure detection |
| 242 passing tests | |

## Testing

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

## Contributing

Contributions welcome! See the [contributing guide](https://volkerkarle.github.io/SymbolicDiagonalization.jl/contributing/) for details.

## License

GPL v3 - see [LICENSE](LICENSE) for details.
