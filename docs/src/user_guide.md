# User Guide

Practical guide to using SymbolicDiagonalization.jl effectively.

## Installation

```julia
using Pkg
Pkg.add("SymbolicDiagonalization")
```

Or from the repository:

```bash
cd SymbolicDiagonalization.jl
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Basic Usage

```julia
using Symbolics
using SymbolicDiagonalization

@variables a b c d

# Create a symbolic matrix
M = [a b; b a]

# Eigenvalues only (faster)
eigvals(M)  # [a+b, a-b]

# Eigenvalues + eigenvectors
E = eigen(M)
E.values    # eigenvalues
E.vectors   # eigenvectors as columns
```

## Common Workflows

### Small Matrices (≤ 4×4)

Direct closed-form solution via root formulas:

```julia
@variables a b c d
M = [a b 0 0; b a 0 0; 0 0 c d; 0 0 d c]

eigvals(M)  # [a+b, a-b, c+d, c-d]
```

### Block-Diagonal Matrices

Automatically detected and solved recursively:

```julia
@variables a b c d e f
M = [a b 0 0 0 0;
     b a 0 0 0 0;
     0 0 c d 0 0;
     0 0 d c 0 0;
     0 0 0 0 e f;
     0 0 0 0 f e]

eigvals(M)  # Solves three 2×2 blocks
```

### Circulant Matrices

DFT-based solution for any size:

```julia
@variables a b c d
C = [a b c d; d a b c; c d a b; b c d a]

eigvals(C)  # Works for 100×100 too
```

### Tridiagonal Matrices

Closed-form cosine formula:

```julia
@variables a b
n = 5
T = diagm(0 => fill(a, n), 1 => fill(b, n-1), -1 => fill(b, n-1))

eigvals(T)  # λₖ = a + 2b·cos(kπ/(n+1))
```

### Kronecker Products

Products of factor eigenvalues:

```julia
@variables a b c d
A = [a 0; 0 b]
B = [c 0; 0 d]
K = kron(A, B)  # 4×4

eigvals(K)  # [ac, ad, bc, bd]
```

## Performance Tips

| Tip | Example |
|-----|---------|
| Use structure hints | `eigvals(M, structure=:hermitian)` |
| Eigenvalues only | `eigvals(M)` vs `eigen(M)` |
| Partial substitution | `substitute(M, Dict(c => 0))` |
| Increase timeout | `eigvals(M, timeout=600)` |
| Disable expansion | `eigvals(M, expand=false)` |

## Troubleshooting

### "Cannot solve degree n ≥ 5 without structure"

Matrix is 5×5+ with no detected structure. Solutions:
1. Add zeros to create block structure
2. Check for circulant, tridiagonal, or Kronecker patterns
3. Substitute numeric values for some variables
4. Use numerical methods via `LinearAlgebra.eigvals` on Float64 matrix

### Expression Complexity Error

Expressions exceed term limit. Solutions:
1. Increase limit: `eigvals(M, max_terms=50000)`
2. Add structure (zeros, blocks)
3. Disable expansion: `eigvals(M, expand=false)`

### Computation Timeout

Exceeded time limit. Solutions:
1. Increase timeout: `eigvals(M, timeout=900)`
2. Simplify matrix (fewer variables, more structure)
3. Use partial numeric substitution

### Memory Issues

Fully symbolic 4×4 quartic produces very large expressions. Solutions:
1. Add structure to avoid quartic formula
2. Use block decomposition
3. Partial numeric substitution

## Verification

```julia
# Verify diagonalization
P, D, pairs = symbolic_diagonalize(M)
using Symbolics: simplify
@assert all(iszero, simplify.(M - P * D * inv(P)))

# Verify numerically
M_test = substitute(M, Dict(a=>1, b=>2))
λ_sym = substitute.(eigvals(M), Ref(Dict(a=>1, b=>2)))
λ_num = LinearAlgebra.eigvals(Float64.(M_test))
@assert isapprox(sort(real.(λ_sym)), sort(λ_num))
```

## Physical Examples

### Vibrating String

```julia
@variables k m
n = 5
K = diagm(0 => fill(2k, n), 1 => fill(-k, n-1), -1 => fill(-k, n-1))
λ = eigvals(K)  # ωₖ² = λₖ/m
```

### Tight-Binding Model

```julia
@variables t ε
n = 4
H = diagm(0 => fill(ε, n), 1 => fill(-t, n-1), -1 => fill(-t, n-1))
E = eigvals(H)  # Eₖ = ε - 2t·cos(kπ/(n+1))
```

### Cycle Graph Laplacian

```julia
@variables d
L = [2d -d 0 0 0 -d; -d 2d -d 0 0 0; 0 -d 2d -d 0 0;
     0 0 -d 2d -d 0; 0 0 0 -d 2d -d; -d 0 0 0 -d 2d]
λ = eigvals(L)  # λₖ = 2d(1 - cos(2πk/n))
```

## See Also

- [API Reference](api_reference.md) - Complete function signatures
- [Pattern Library](pattern_library.md) - All 13 patterns with details
- [Mathematical Background](mathematical_background.md) - Theory
