# SymbolicDiagonalization.jl

**Status: Functional Prototype**

SymbolicDiagonalization.jl provides symbolic matrix diagonalization for Julia using `Symbolics.jl` with closed-form root solvers and structure detection.

## Vision

The Abel-Ruffini theorem limits closed-form solutions to polynomials of degree $\leq 4$, meaning general $5 \times 5$+ matrices cannot be solved symbolically. However, many real-world matrices have exploitable structure.

**Goal**: Build automatic structure detection to solve larger symbolic problems by recognizing and exploiting special patterns (block-diagonal, persymmetric, tridiagonal, circulant, Kronecker products, etc.).

**Current state**: Functional prototype with 9 special pattern solvers working for any size $n$. Comprehensive test coverage (172 tests, all passing). Ready for experimental use.

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
- **`expand`** - Expand polynomial coefficients (default: `true`)
- **`simplify`** - Simplify eigenvalues/eigenvectors (default: `false`)
- **`timeout`** - Maximum computation time in seconds (default: 300)
- **`max_terms`** - Expression complexity limit (default: 10000)
- **`complexity_threshold`** - Warn when symbolic variable count exceeds this (default: 5)
- **`check_diagonalizable`** - Verify eigenvectors are independent (default: `true`)
- **`var`** - Custom variable for characteristic polynomial (default: auto-generated)

## Examples

### Block-Diagonal Matrix

```@example main
@variables a b c d λ
mat = [a  b  0  0;
       b  a  0  0;
       0  0  c  d;
       0  0  d  c]

vals = eigvals(mat)
# Result: [a+b, a-b, c+d, c-d]
```

### Circulant Matrix (Any Size)

```@example main
@variables a b c
mat = [a  b  c;
       c  a  b;
       b  c  a]

vals = eigvals(mat)
# Uses DFT: eigenvalues are linear combinations with roots of unity
```

### Kronecker Product

```@example main
@variables a b c d
A = [a 0; 0 b]
B = [c 0; 0 d]
mat = kron(A, B)  # 4×4 matrix

vals = eigvals(mat)
# Result: [a*c, a*d, b*c, b*d] (products of eigenvalues)
```

### Symmetric Toeplitz Tridiagonal

```@example main
@variables a b
mat = [a  b  0;
       b  a  b;
       0  b  a]

vals = eigvals(mat)
# Closed-form using cosine formula: a + 2b*cos(kπ/(n+1))
```

### Permutation Matrix

```@example main
# 3-cycle permutation
mat = [0  1  0;
       0  0  1;
       1  0  0]

vals = eigvals(mat)
# Result: roots of unity [1, ω, ω²] where ω = exp(2πi/3)
```

## Implementation Details

### Characteristic Polynomial
- Computed via fraction-free Bareiss determinant on $\lambda I - A$
- Coefficients extracted by differentiating at $\lambda = 0$
- Works with symbolic rings

### Closed-Form Root Solvers (Degrees 1-4)
- Linear, quadratic, cubic (Cardano), quartic (Ferrari)
- Quartic produces very large expressions for fully symbolic matrices
- Degrees $\geq 5$ require structure detection or throw error

### Structure Detection (13 Patterns Implemented)

**Successfully implemented patterns:**
1. **Diagonal** - Trivial $O(n)$ eigenvalue extraction
2. **Triangular** - Diagonal elements are eigenvalues
3. **Block-Diagonal** - Recursive detection and solving (works with any block sizes)
4. **Persymmetric** - $Q[i,j] = Q[n+1-j,n+1-i]$, splits into half-sized eigenproblems
5. **Circulant** - Uses DFT, works for any size $n$
6. **Block Circulant** - Block-level DFT approach
7. **Symmetric Toeplitz Tridiagonal** - Closed-form cosine formula for any $n$
8. **Anti-Diagonal** - Eigenvalues come in $\pm$ pairs
9. **Permutation Matrices** - Eigenvalues are roots of unity
10. **Kronecker Products** - Eigenvalues are products of factor eigenvalues
11. **Special $5 \times 5$ pattern [b,d,b,b]** - Known closed-form solution
12. **Special $5 \times 5$ pattern [b,b,d,b]** - Known closed-form solution
13. **Jordan Blocks** - Repeated eigenvalue with known structure

**Detection robustness**: Currently basic pattern matching. Room for improvement with more sophisticated algorithms and tolerance handling for near-patterns.

### Eigenvector Computation
- Nullspace via RREF-based row reduction
- `symbolic_diagonalize` verifies eigenvector independence
- Error handling for non-diagonalizable matrices

## Current Capabilities and Limitations

### What Works Well

- **Small matrices ($\leq 4 \times 4$)**: Full symbolic diagonalization using closed-form root solvers
- **Structured matrices (any size)**: 13 special patterns with efficient $O(n)$ to $O(n^2)$ algorithms
- **Mixed symbolic/numeric**: Handles partially symbolic matrices effectively
- **Comprehensive testing**: 172 passing tests (37.4s) covering diverse scenarios

### Known Limitations

- **General $5 \times 5$+ matrices**: No closed-form solver (Abel-Ruffini theorem) - requires structure detection
- **Expression complexity**: Fully symbolic $4 \times 4$ quartic eigenvalues can be very large (~13.5 MB)
- **Simplification**: Basic simplification only; results may not be in minimal form
- **Pattern detection robustness**: Exact pattern matching; may miss near-patterns or numerically perturbed structures

## Development Priorities

### High Priority

- **Enhanced structure detection**: More robust algorithms with tolerance handling for near-patterns
- **Additional patterns**: Arrow matrices, general Hankel, more tridiagonal families
- **Expression simplification**: Deeper integration with `Symbolics.jl` simplification
- **Performance optimization**: Caching, memoization for repeated subproblems

### Medium Priority

- **Eigenvalue multiplicity**: Better handling of repeated eigenvalues and generalized eigenvectors
- **Numerical fallback**: Automatic hybrid symbolic-numeric mode for borderline cases
- **Documentation**: More examples, tutorials, pattern discovery guides
- **Symbolic conditioning**: Analyze stability and condition numbers symbolically

### Long-Term Vision

- **Machine learning**: Pattern recognition using trained models
- **User pattern libraries**: Custom pattern registration and sharing
- **Symbolic perturbation**: First-order eigenvalue sensitivity analysis
- **Parallel processing**: Multi-threaded block decomposition

## Documentation

### User Guides
- **[User Guide](user_guide.md)** - Installation, basic usage, workflow examples
- **[API Reference](api_reference.md)** - Complete function signatures and keyword arguments
- **[Pattern Library](pattern_library.md)** - All 13 implemented patterns with examples and complexity analysis

### Developer Resources
- **[Implementation Details](implementation.md)** - Algorithm descriptions, complexity analysis, performance considerations
- **[Mathematical Background](mathematical_background.md)** - Theory behind root solvers and pattern-specific algorithms
- **[Contributing Guide](contributing.md)** - Development setup, testing, adding new patterns

### Pattern Discovery

Additional resources for pattern research and exploration:
- **Research Summary** (`notes/RESEARCH_SUMMARY.md`) - Pattern discovery overview
- **Discovery Methodology** (`notes/DISCOVERY_METHODOLOGY.md`) - Pattern exploration techniques
- **Pattern Discoveries** (`notes/PATTERN_DISCOVERIES.md`) - Detailed catalog of investigated patterns
- **Explore Patterns Script** (`examples/explore_patterns.jl`) - Interactive pattern exploration tool

## Building the Documentation

From the project root:

```bash
julia --project=docs docs/make.jl
```

This renders the documentation locally using Documenter.jl.

### Viewing the Documentation

The documentation is built as static HTML in `docs/build/`. To view it properly with working navigation links, use a local web server:

```bash
# Option 1: Use the provided script
cd docs
./serve.sh

# Option 2: Use Python's built-in HTTP server
cd docs/build
python3 -m http.server 8000

# Option 3: Use Julia's HTTP server
cd docs/build
julia -e 'using HTTP; HTTP.serve(HTTP.Files("."))'
```

Then open your browser to `http://localhost:8000`

**Note**: Opening `index.html` directly with `file://` may cause navigation links to not work properly in some browsers due to security restrictions.
