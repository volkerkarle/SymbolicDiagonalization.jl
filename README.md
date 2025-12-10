# SymbolicDiagonalization.jl

A Julia package for symbolic matrix diagonalization. Compute eigenvalues, eigenvectors, and full diagonalizations directly on `Symbolics.jl` expressions using closed-form root solvers and intelligent structure detection.

## Quick Start

```julia
using Symbolics, SymbolicDiagonalization, LinearAlgebra

@variables a b c
mat = [a 1 0; 0 b 1; 0 0 c]  # upper-triangular symbolic matrix

# Standard LinearAlgebra interface
E = eigen(mat)        # Eigen object with .values and .vectors
λ = eigvals(mat)      # eigenvalues only (faster)
```

## Features

### Closed-Form Solutions

- **Degrees 1-4**: Linear, quadratic, Cardano (cubic), Ferrari (quartic)
- **Characteristic polynomials**: Fraction-free Bareiss determinant
- **Symbolic nullspace**: Direct computation of eigenvectors

### Intelligent Structure Detection

- **Block-diagonal matrices**: Automatic decomposition into independent subproblems
- **Persymmetric matrices**: Specialized optimization for symmetric-about-antidiagonal structure
- **Special 5×5 tridiagonal patterns**: Closed-form solutions for patterns like `[b, d, b, b]`
  - Eigenvalues: `a ± √(2b² + d²)`, `a ± b`, `a`
  - See [docs/research/](docs/research/) for pattern discovery methodology

### Flexible API

- **LinearAlgebra interface**: `eigen()` and `eigvals()` for familiar usage
- **Original API**: Advanced control with `symbolic_eigenvalues()`, `symbolic_eigenpairs()`, `symbolic_diagonalize()`
- **Optional Oscar backend**: Use Oscar/Nemo for algebraic numbers without radical expansion (numeric/rational inputs)

## Installation

The package is currently local. From the package root:

```julia
julia --project -e 'using Pkg; Pkg.instantiate()'
```

To develop:

```julia
julia --project -e 'using Pkg; Pkg.develop(path=".")'
```

## Examples

### Basic Usage

```julia
using SymbolicDiagonalization, Symbolics, LinearAlgebra

# Simple 2×2 matrix
@variables a b
A = [a 0; 0 b]
E = eigen(A)
E.values   # [a, b]
E.vectors  # [1 0; 0 1]

# Eigenvalues only (faster)
λ = eigvals(A)  # [a, b]
```

### Structured Matrices

```julia
# Block-diagonal 4×4
@variables a b c d
mat = [a  b  0  0;
       b  a  0  0;
       0  0  c  d;
       0  0  d  c]

vals = eigvals(mat; structure=:hermitian)
# Result: [a+b, a-b, c+d, c-d] - clean and fast!

# 5×5 special tridiagonal pattern
@variables a b d
mat = [a  b  0  0  0;
       b  a  d  0  0;
       0  d  a  b  0;
       0  0  b  a  b;
       0  0  0  b  a]

vals = eigvals(mat)
# Closed-form: a ± √(2b² + d²), a ± b, a
```

### Advanced API

```julia
# Get characteristic polynomial too
vals, poly, λ = symbolic_eigenvalues(mat)

# Get eigenvalue-eigenvector pairs
pairs, _, _ = symbolic_eigenpairs(mat)
for (eigenval, eigenvecs) in pairs
    println("λ = $eigenval")
    println("v = $eigenvecs")
end

# Full diagonalization: A = P D P⁻¹
P, D, _ = symbolic_diagonalize(mat)
```

## API Reference

### LinearAlgebra Interface (Recommended)

- **`eigen(A; kwargs...)`** - Returns `Eigen` object with `.values` and `.vectors` fields
- **`eigvals(A; kwargs...)`** - Returns eigenvalues only (skips eigenvector computation)

### Original API (Advanced)

- **`symbolic_eigenvalues(A; kwargs...)`** - Returns eigenvalues, characteristic polynomial, and λ symbol
- **`symbolic_eigenpairs(A; kwargs...)`** - Returns eigenvalue-eigenvector pairs
- **`symbolic_diagonalize(A; kwargs...)`** - Returns diagonalization matrices `P`, `D`

### Common Options

- **`structure`** - Hint for matrix type: `:auto` (default), `:hermitian`, `:symmetric`, `:unitary`
- **`expand`** - Expand polynomial (default: `true`); set to `false` for large matrices
- **`complexity_threshold`** - Warn when symbolic variable count exceeds this (default: 5)
- **`backend`** - `:symbolics` (default) or `:oscar` (requires Oscar.jl for algebraic numbers)

## Performance Guide

### What Works Well

**✅ Diagonal/Triangular matrices**

- Instant eigenvalue extraction
- Any size

**✅ Block-diagonal matrices**

- Decomposes into independent subproblems
- 2×2 blocks use quadratic formula → clean results

**✅ Sparse symbolic matrices**

- Few parameters (≤5 variables recommended)
- Many zero entries

**✅ Numeric 4×4 matrices**

- ~10K characters per eigenvalue
- ~3-4 seconds computation

### Limitations

**⚠️ Fully symbolic matrices**

- 3×3: Usually practical
- 4×4: Works but produces very large expressions (~13.5 MB per eigenvalue)
- 5×5+: Only closed-form with special structure (see research docs)

**❌ Degree ≥5 without special structure**

- No closed-form solutions available
- Use numeric methods instead

### Optimization Tips

1. **Minimize symbolic parameters**: Fewer distinct variables = faster computation
2. **Use `expand=false`**: Skip polynomial expansion for eigenvalues-only queries
3. **Exploit structure**: Block-diagonal structure is automatically detected and optimized
4. **Numeric inputs**: Use `backend=:oscar` for algebraic numbers without radicals

## Testing

Run the test suite:

```julia
julia --project -e 'using Pkg; Pkg.test()'
```

Tests cover:

- Characteristic polynomials (Bareiss determinant)
- Structure detection (diagonal, triangular, block-diagonal, persymmetric)
- Closed-form root solvers (degrees 1-4)
- Special patterns (5×5 tridiagonal)
- Non-diagonalizable matrices (error handling)
- Numeric vs symbolic cases

## Documentation

Build local documentation:

```julia
julia --project=docs docs/make.jl
```

## Research

Pattern discovery methodology and results are documented in:

- [docs/research/RESEARCH_SUMMARY.md](docs/research/RESEARCH_SUMMARY.md) - Executive summary
- [docs/research/DISCOVERY_METHODOLOGY.md](docs/research/DISCOVERY_METHODOLOGY.md) - Systematic exploration guide
- [docs/research/PATTERN_DISCOVERIES.md](docs/research/PATTERN_DISCOVERIES.md) - Detailed pattern catalog

Explore new patterns using:

- [examples/explore_patterns.jl](examples/explore_patterns.jl) - Automated pattern search tool

## Contributing

This is a research project exploring symbolic eigenvalue computation. Contributions welcome, especially:

- New special pattern detection
- Performance optimizations
- Extended test coverage
- Documentation improvements

## License

See LICENSE file for details.
