# SymbolicDiagonalization.jl

<p align="center">
  <img src="logo.png" alt="SymbolicDiagonalization.jl Logo" width="300"/>
</p>

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://volkerkarle.github.io/SymbolicDiagonalization.jl)
[![CI](https://github.com/volkerkarle/SymbolicDiagonalization.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/volkerkarle/SymbolicDiagonalization.jl/actions/workflows/ci.yml)
[![Documentation Build](https://github.com/volkerkarle/SymbolicDiagonalization.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/volkerkarle/SymbolicDiagonalization.jl/actions/workflows/documentation.yml)
[![codecov](https://codecov.io/gh/volkerkarle/SymbolicDiagonalization.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/volkerkarle/SymbolicDiagonalization.jl)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**Symbolic eigenvalue computation for structured matrices**

## The Problem

The Abel-Ruffini theorem proves no closed-form solution exists for polynomials of degree 5+. This means general symbolic 5×5+ matrices cannot be diagonalized. However, many matrices have exploitable structure.

## The Solution

This package automatically detects matrix structure and returns closed-form symbolic eigenvalues:

```julia
using Symbolics, SymbolicDiagonalization, LinearAlgebra

# 1024×1024 nested Kronecker product - degree 1024 polynomial!
# But structure reduces it to 10 quadratic problems
@variables a1 b1 c1  a2 b2 c2  a3 b3 c3  a4 b4 c4  a5 b5 c5  a6 b6 c6  a7 b7 c7  a8 b8 c8  a9 b9 c9  a10 b10 c10

matrices = [[a1 b1; b1 c1], [a2 b2; b2 c2], [a3 b3; b3 c3], [a4 b4; b4 c4], [a5 b5; b5 c5],
            [a6 b6; b6 c6], [a7 b7; b7 c7], [a8 b8; b8 c8], [a9 b9; b9 c9], [a10 b10; b10 c10]]
K = reduce(kron, matrices)  # 1024×1024 with 30 symbolic parameters

eigvals(K)  # 1024 symbolic eigenvalues in ~33 seconds
```

## Supported Patterns

| Pattern | Size | Method |
|---------|------|--------|
| Any matrix | ≤4×4 | Quadratic/Cardano/Ferrari formulas |
| Block-diagonal | Any | Recursive decomposition |
| Circulant | Any | DFT diagonalization |
| Kronecker A⊗B | Any | λ(A)·λ(B) factorization |
| Nested Kronecker | Any | Recursive factorization |
| Symmetric Toeplitz tridiagonal | Any | Chebyshev formula |
| Hadamard | 2ⁿ | ±√(2ⁿ) |
| DFT | Any | Fourth roots of n |
| Permutation | Any | Cycle decomposition |
| SO(2), SO(3), SO(4) | 2,3,4 | Rotation angle extraction |
| SU(2), SU(3) | 2,3 | Lie algebra structure |
| Cartan type Aₙ | Any | Root system formula |
| Path/Cycle Laplacian | Any | Trigonometric formula |

## Quick Start

```julia
using Symbolics, SymbolicDiagonalization, LinearAlgebra

@variables a b c d
```

**Block-diagonal (recursive decomposition):**
```julia
M = [a b 0 0; b a 0 0; 0 0 c d; 0 0 d c]
eigvals(M)  # [a+b, a-b, c+d, c-d]
```

**Circulant (DFT formula works for any n):**
```julia
C = [a b c; c a b; b c a]
eigvals(C)  # Closed-form for any size
```

**Kronecker product:**
```julia
A = [a 0; 0 b]
B = [c 0; 0 d]
eigvals(kron(A, B))  # [ac, ad, bc, bd]
```

**Rotation matrices:**
```julia
@variables θ φ
eigvals(kron(SO2_rotation(θ), SO2_rotation(φ)))  # e^{i(±θ±φ)}
```

## Installation

```julia
using Pkg
Pkg.add("SymbolicDiagonalization")
```

## Documentation

**[Full documentation](https://volkerkarle.github.io/SymbolicDiagonalization.jl)**

- [User Guide](https://volkerkarle.github.io/SymbolicDiagonalization.jl/user_guide/) - Examples and workflows
- [Pattern Library](https://volkerkarle.github.io/SymbolicDiagonalization.jl/pattern_library/) - All supported patterns
- [Mathematical Background](https://volkerkarle.github.io/SymbolicDiagonalization.jl/mathematical_background/) - Theory

## License

GPL v3 - see [LICENSE](LICENSE).
