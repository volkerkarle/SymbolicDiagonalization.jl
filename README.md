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
- **19+ special pattern solvers** for arbitrary-sized matrices (circulant, Kronecker, Hadamard, DFT, Toeplitz tridiagonal, permutation, Q₈ regular representation, etc.)
- **Lie group detection** (SO(2)-SO(4), SU(2), SU(3), Sp(2), Sp(4)) with symbolic eigenvalues
- **SO(2) Kronecker products** with automatic trig simplification (e.g., `cos(θ+φ) + i·sin(θ+φ)`)
- **SU(2) Kronecker products** with half-angle eigenvalues (e.g., `cos((α+β)/2) + i·sin((α+β)/2)`)
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

### SO(2) Rotation Kronecker Products

```julia
@variables θ φ
# Use SO2_rotation for clean 2×2 rotation matrices
R_θ = SO2_rotation(θ)
R_φ = SO2_rotation(φ)

# Kronecker product of rotations → clean trig eigenvalues
eigvals(kron(R_θ, R_φ))
# [cos(θ+φ) + im*sin(θ+φ), cos(θ-φ) + im*sin(θ-φ), 
#  cos(θ-φ) - im*sin(θ-φ), cos(θ+φ) - im*sin(θ+φ)]

# Same-angle case automatically simplifies
eigvals(kron(SO2_rotation(θ), SO2_rotation(θ)))
# [cos(2θ) + im*sin(2θ), 1, 1, cos(2θ) - im*sin(2θ)]
```

### SU(2) Kronecker Products

```julia
@variables α β
# SU(2) rotations use half-angles for spin-1/2 representation
K = SU2_kron([α, β])  # 4×4 matrix
eigvals(K)
# [cos((α+β)/2) + im*sin((α+β)/2), cos((α-β)/2) + im*sin((α-β)/2),
#  cos((α-β)/2) - im*sin((α-β)/2), cos((α+β)/2) - im*sin((α+β)/2)]
```

### Aggressive Symbolic Simplification

Clean eigenvalue expressions via automatic trigonometric simplification:

```julia
@variables θ
using SymbolicDiagonalization: aggressive_simplify

# sqrt(1 - cos²θ) automatically becomes sin(θ)
aggressive_simplify(sqrt(1 - cos(θ)^2))  # sin(θ)

# SO(3) rotation gives clean eigenvalues
Rz = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
eigvals(Rz)  # [1, cos(θ) + im*sin(θ), cos(θ) - im*sin(θ)]
```

### Closed-Form Eigenvectors for Lie Groups

For rotation and symplectic matrices, the package computes eigenvectors analytically without nullspace computation:

```julia
@variables θ

# SO(2) has fixed eigenvectors regardless of angle
R = SO2_rotation(θ)
pairs, _, _ = symbolic_eigenpairs(R)
# pairs[1] = (cos(θ) + im*sin(θ), [[1, -im]])  # e^{iθ} with eigenvector [1, -i]
# pairs[2] = (cos(θ) - im*sin(θ), [[1, im]])   # e^{-iθ} with eigenvector [1, i]

# SO(3) axis-aligned rotations
Rz = SO3_Rz(θ)
pairs, _, _ = symbolic_eigenpairs(Rz)
# pairs[1] = (1, [[0, 0, 1]])                  # eigenvalue 1 (rotation axis)
# pairs[2] = (cos(θ) + im*sin(θ), [[1, im, 0]])
# pairs[3] = (cos(θ) - im*sin(θ), [[1, -im, 0]])

# SO(4) block-diagonal: two independent rotations
@variables φ
R4 = [cos(θ) -sin(θ) 0 0; sin(θ) cos(θ) 0 0; 0 0 cos(φ) -sin(φ); 0 0 sin(φ) cos(φ)]
pairs, _, _ = symbolic_eigenpairs(R4)
# Eigenvectors are padded: [1, ±i, 0, 0] and [0, 0, 1, ±i]

# Full diagonalization: A = P * D * P⁻¹
P, D, pairs = symbolic_diagonalize(R)
```

Supported groups: **SO(2), SO(3), SO(4), SU(2), Sp(2), Sp(4)** (non-diagonal cases only)

### Full 6-Parameter 3×3 Symmetric Matrix

```julia
@variables a b c d e f
A = [a b c; b d e; c e f]  # 6 independent parameters

eigvals(A)  # Works! Uses diagonal shift optimization
```

### Hadamard Matrices (any power of 2)

```julia
H = hadamard_matrix(3)  # 8×8 Sylvester-Hadamard matrix
eigvals(H)              # [±2√2], each with multiplicity 4

# Works for any size 2ⁿ
H32 = hadamard_matrix(5)  # 32×32
eigvals(H32)              # [±√32], each with multiplicity 16
```

### DFT (Fourier) Matrices (any size)

```julia
F = dft_matrix(8)       # 8×8 Fourier matrix
eigvals(F)              # {±√8, ±i√8} with multiplicities (3,2,2,1)

# Normalized (unitary) version
F_norm = dft_matrix(8, normalized=true)
eigvals(F_norm)         # {1, -1, i, -i}
```

### Quaternion Group Q₈ Regular Representation

```julia
# Q₈ is the 8-element quaternion group {±1, ±i, ±j, ±k}
# The regular representation gives 8×8 matrices with closed-form eigenvalues
M = Q8_invariant_matrix(1.0, 0.5, 0.3, 0.2, 0.4, 0.1, 0.25, 0.15)
eigvals(M)
# Character theory formulas:
# - 4 one-dimensional irreps (multiplicity 1 each)
# - 1 two-dimensional irrep with eigenvalue c₁ - c₋₁ (multiplicity 4)
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
- [Pattern Library](https://volkerkarle.github.io/SymbolicDiagonalization.jl/pattern_library/) - All 19+ supported matrix patterns
- [Mathematical Background](https://volkerkarle.github.io/SymbolicDiagonalization.jl/mathematical_background/) - Theory and proofs
- [Contributing](https://volkerkarle.github.io/SymbolicDiagonalization.jl/contributing/) - Development guide

## Capabilities

| Feature | Details |
|---------|---------|
| All matrices up to 4×4 | Closed-form solutions via root formulas |
| Full 6-parameter 3×3 symmetric | Diagonal shift optimization |
| Block-diagonal decomposition | Automatic detection and recursion |
| 19+ pattern solvers | Circulant, Kronecker, Hadamard, DFT, tridiagonal, permutation, Q₈, etc. |
| Lie group detection | SO(2)-SO(4), SU(2), SU(3), Sp(2), Sp(4) with symbolic eigenvalues |
| **Lie group eigenvectors** | Closed-form eigenvectors for SO(2)-SO(4), SU(2), Sp(2), Sp(4) (non-diagonal) |
| SO(2) Kronecker products | `cos(θ±φ) + i·sin(θ±φ)` form via trig simplification |
| SU(2) Kronecker products | `cos((α±β)/2) + i·sin((α±β)/2)` with half-angle formulas |
| Nested Kronecker A₁⊗A₂⊗...⊗Aₙ | Scales to 1024×1024 with 30 parameters |
| 382 passing tests | Comprehensive test coverage |

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
