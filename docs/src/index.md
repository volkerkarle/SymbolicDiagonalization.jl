# SymbolicDiagonalization.jl

SymbolicDiagonalization.jl provides symbolic matrix diagonalization for Julia using `Symbolics.jl`. Compute eigenvalues, eigenvectors, and full diagonalizations using closed-form root solvers and intelligent structure detection.

## Quick Start

```@example main
using SymbolicDiagonalization
using LinearAlgebra
using Symbolics

@variables a b c
mat = [a 1 0; 0 b 1; 0 0 c]  # upper-triangular matrix

# Standard LinearAlgebra interface
E = eigen(mat)        # Eigen object with .values and .vectors
λ = eigvals(mat)      # eigenvalues only (faster)
```

## API Overview

### LinearAlgebra Interface (Recommended)

- **`eigen(A; kwargs...)`** - Returns `Eigen` object with `.values` and `.vectors` fields
- **`eigvals(A; kwargs...)`** - Returns eigenvalues only (skips eigenvector computation)

### Original API (Advanced)

- **`characteristic_polynomial(A; var)`** → `(poly, coeffs, λ)`
- **`symbolic_eigenvalues(A; kwargs...)`** → `(vals, poly, λ)`
- **`symbolic_eigenpairs(A; kwargs...)`** → eigenvalue-eigenvector pairs
- **`symbolic_diagonalize(A; kwargs...)`** → `(P, D, pairs)` (throws if not diagonalizable)

### Common Options

- **`structure`** - Hint for matrix type: `:auto`, `:hermitian`, `:symmetric`, `:unitary`, `:none`
- **`expand`** - Expand polynomial (default: `true`); set to `false` for large matrices
- **`complexity_threshold`** - Warn when symbolic variable count exceeds this (default: 5)
- **`backend`** - `:symbolics` (default) or `:oscar` (requires Oscar.jl)

## Examples

### Basic Usage

```@example main
using SymbolicDiagonalization, LinearAlgebra, Symbolics

# Diagonal matrix
@variables a b
A = [a 0; 0 b]
E = eigen(A)
E.values   # [a, b]
E.vectors  # [1 0; 0 1]
```

### Block-Diagonal Matrix

```@example main
# 4×4 block-diagonal (automatically detected)
@variables a b c d λ
mat = [a  b  0  0;
       b  a  0  0;
       0  0  c  d;
       0  0  d  c]

vals = eigvals(mat; structure=:hermitian)
# Result: [a+b, a-b, c+d, c-d] - clean and fast!
```

### Special 5×5 Pattern

```@example main
# Special tridiagonal pattern with closed-form eigenvalues
@variables a b d
mat = [a  b  0  0  0;
       b  a  d  0  0;
       0  d  a  b  0;
       0  0  b  a  b;
       0  0  0  b  a]

vals = eigvals(mat)
# Closed-form: a ± √(2b² + d²), a ± b, a
```

## Implementation Details

### Characteristic Polynomial
- Computed via fraction-free Bareiss determinant on `λI - A`
- Coefficients extracted by differentiating at `λ = 0`
- Robust on symbolic rings

### Closed-Form Root Solvers
- Linear, quadratic, cubic (Cardano), quartic (Ferrari)
- Quartic solver produces large expressions for fully symbolic matrices
- Degrees ≥ 5 require special structure or throw error (Abel-Ruffini theorem)

### Structure Detection
The package automatically detects and exploits:

1. **Diagonal/Triangular** - Direct eigenvalue extraction from diagonal
2. **Block-diagonal** - Recursive decomposition into independent subproblems
3. **Persymmetric** - Symmetric matrices with `Q[i,j] = Q[n+1-j,n+1-i]` split into half-sized blocks
4. **Special 5×5 tridiagonal patterns** - Patterns `[b,d,b,b]` and `[b,b,d,b]` have closed-form solutions

### Eigenvector Computation
- Nullspace computed via `_rref`-based row reduction
- `symbolic_diagonalize` checks eigenvector independence
- Errors with explanation when matrix is not diagonalizable

### Optional Oscar Backend
- `backend=:oscar` uses Oscar/Nemo for algebraic numbers
- Works with numeric/rational coefficients
- Falls back to symbolic backend with warning for symbolic coefficients

## Performance Guide

### What Works Well

**Diagonal/Triangular matrices** - Any size, instant extraction

**Block-diagonal matrices** - Automatic decomposition into subproblems

**Sparse symbolic matrices** - Few parameters (≤5 variables recommended)

**Numeric 4×4 matrices** - ~10K chars per eigenvalue, ~3-4 seconds

### Limitations

**Fully symbolic 4×4** - Works but produces very large expressions (~13.5 MB per eigenvalue)

**Degree ≥5 without special structure** - No closed-form solutions (Abel-Ruffini theorem)

**Dense symbolic matrices** - Expression size grows rapidly

### Optimization Tips

1. **Minimize symbolic parameters** - Fewer variables = faster computation
2. **Use `expand=false`** - Skip polynomial expansion for eigenvalues-only
3. **Exploit structure** - Block-diagonal is automatically detected
4. **Numeric inputs** - Use `backend=:oscar` for algebraic numbers

## Research

This package includes systematic pattern discovery for special matrices with closed-form eigenvalues:

- **[research/RESEARCH_SUMMARY.md](../research/RESEARCH_SUMMARY.md)** - Executive summary
- **[research/DISCOVERY_METHODOLOGY.md](../research/DISCOVERY_METHODOLOGY.md)** - Pattern exploration methodology
- **[research/PATTERN_DISCOVERIES.md](../research/PATTERN_DISCOVERIES.md)** - Detailed pattern catalog

Explore new patterns using **[examples/explore_patterns.jl](../../examples/explore_patterns.jl)**.

## Building the Documentation

From the project root:

```bash
julia --project=docs docs/make.jl
```

This renders the documentation locally using Documenter.jl.
