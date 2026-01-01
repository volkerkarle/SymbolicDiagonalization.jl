# SymbolicDiagonalization.jl

Symbolic matrix diagonalization for Julia using [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl).

## Features

- Closed-form eigenvalue solvers for degrees 1-4 (linear, quadratic, Cardano, Ferrari)
- Automatic structure detection for larger matrices (13 patterns)
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

## Supported Patterns

For matrices larger than 4×4, the package detects and exploits special structure:

| Pattern | Description |
|---------|-------------|
| Diagonal | Direct extraction |
| Triangular | Diagonal elements |
| Block-diagonal | Recursive solving |
| Circulant | DFT-based (any size n) |
| Symmetric Toeplitz tridiagonal | Cosine formula (any size n) |
| Persymmetric | Half-size decomposition |
| Kronecker products | Product of factor eigenvalues |
| Permutation | Roots of unity |

See [Pattern Library](pattern_library.md) for the complete list with examples.

## Limitations

- **General 5×5+ matrices**: No closed-form solution exists (Abel-Ruffini theorem) - requires exploitable structure
- **Expression size**: Fully symbolic 4×4 quartics can produce large expressions
- **Simplification**: Results may not be in minimal form

## Documentation

**Getting Started**
- [User Guide](user_guide.md) - Installation, workflows, troubleshooting
- [API Reference](api_reference.md) - Complete function signatures

**Patterns & Theory**
- [Pattern Library](pattern_library.md) - All 13 patterns with examples
- [Mathematical Background](mathematical_background.md) - Theory and algorithms
- [Group Theory Examples](group_theory_examples.md) - Symmetry-based approaches

**Development**
- [Implementation Details](implementation.md) - Internals and algorithms
- [Contributing](contributing.md) - Development setup and guidelines
