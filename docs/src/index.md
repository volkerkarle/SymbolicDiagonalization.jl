# SymbolicDiagonalization.jl

**Status: Early Work in Progress**

SymbolicDiagonalization.jl provides symbolic matrix diagonalization for Julia using `Symbolics.jl` with closed-form root solvers and structure detection.

## Vision

The Abel-Ruffini theorem limits closed-form solutions to polynomials of degree ≤4, meaning general 5×5+ matrices cannot be solved symbolically. However, many real-world matrices have exploitable structure.

**Goal**: Build automatic structure detection to solve larger symbolic problems by recognizing and exploiting special patterns (block-diagonal, persymmetric, tridiagonal, circulant, Toeplitz, etc.).

**Current state**: Proof-of-concept. Basic block-diagonal detection and one special 5×5 pattern implemented. Most of the vision remains unimplemented.

## Quick Start

```@example main
using SymbolicDiagonalization
using LinearAlgebra
using Symbolics

@variables a b c
mat = [a 1 0; 0 b 1; 0 0 c]

E = eigen(mat)        # Eigen object with .values and .vectors
λ = eigvals(mat)      # eigenvalues only (faster)
```

## API Overview

### LinearAlgebra Interface (Recommended)

- **`eigen(A; kwargs...)`** - Returns `Eigen` object with `.values` and `.vectors` fields
- **`eigvals(A; kwargs...)`** - Returns eigenvalues only (skips eigenvector computation)

### Direct API (Advanced)

- **`characteristic_polynomial(A; var)`** → `(poly, coeffs, λ)`
- **`symbolic_eigenvalues(A; kwargs...)`** → `(vals, poly, λ)`
- **`symbolic_eigenpairs(A; kwargs...)`** → eigenvalue-eigenvector pairs
- **`symbolic_diagonalize(A; kwargs...)`** → `(P, D, pairs)` (throws if not diagonalizable)

### Common Options

- **`structure`** - Matrix type hint: `:auto`, `:hermitian`, `:symmetric`, `:unitary`
- **`expand`** - Expand polynomial (default: `true`)
- **`complexity_threshold`** - Warn when symbolic variable count exceeds this (default: 5)
- **`timeout`** - Maximum computation time in seconds (default: 300)
- **`max_terms`** - Expression complexity limit (default: 10000)

## Examples

### Block-Diagonal Matrix (Works Well)

```@example main
@variables a b c d λ
mat = [a  b  0  0;
       b  a  0  0;
       0  0  c  d;
       0  0  d  c]

vals = eigvals(mat)
# Result: [a+b, a-b, c+d, c-d]
```

### Special 5×5 Pattern (Experimental)

```@example main
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
- Works with symbolic rings

### Closed-Form Root Solvers (Degrees 1-4)
- Linear, quadratic, cubic (Cardano), quartic (Ferrari)
- Quartic produces very large expressions for fully symbolic matrices
- Degrees ≥5 require structure detection or throw error

### Structure Detection (Work in Progress)

**Currently implemented** (basic, not robust):
1. **Block-diagonal** - Simple zero-pattern detection, recursive solving
2. **Persymmetric** - Detects `Q[i,j] = Q[n+1-j,n+1-i]`, splits into half-sized blocks
3. **Special 5×5 tridiagonal** - ONE specific pattern `[b,d,b,b]` with known eigenvalues

**Missing** (high priority):
- Circulant matrices
- Toeplitz/Hankel matrices
- General tridiagonal families
- Arrow matrices
- Other structured patterns from literature

### Eigenvector Computation
- Nullspace via RREF-based row reduction
- `symbolic_diagonalize` verifies eigenvector independence
- Error handling for non-diagonalizable matrices

## Current Limitations

- **No general 5×5+ solver**: Requires detectable structure
- **Structure detection is rudimentary**: Misses most patterns
- **Expression explosion**: Fully symbolic 4×4 → ~13.5 MB per eigenvalue
- **Minimal simplification**: Results often not in simplest form
- **Edge cases**: Many not handled properly

## What Needs Work

### High Priority
- Robust structure detection algorithms
- More special patterns (circulant, Toeplitz, Hankel, tridiagonal families)
- Better expression simplification
- Comprehensive edge case testing

### Medium Priority
- Eigenvalue multiplicity handling
- Symbolic condition numbers
- Numerical fallback integration
- Complete pattern documentation

### Future
- Machine learning for pattern recognition
- User-defined pattern libraries
- Symbolic perturbation theory
- Parallel block processing

## Pattern Documentation

Pattern discovery notes and methodology:

- **[notes/RESEARCH_SUMMARY.md](../notes/RESEARCH_SUMMARY.md)** - Pattern discovery summary
- **[notes/DISCOVERY_METHODOLOGY.md](../notes/DISCOVERY_METHODOLOGY.md)** - Pattern exploration methodology
- **[notes/PATTERN_DISCOVERIES.md](../notes/PATTERN_DISCOVERIES.md)** - Detailed pattern catalog

Explore new patterns using **[examples/explore_patterns.jl](../../examples/explore_patterns.jl)**.

## Building the Documentation

From the project root:

```bash
julia --project=docs docs/make.jl
```

This renders the documentation locally using Documenter.jl.
