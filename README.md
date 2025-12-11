# SymbolicDiagonalization.jl

**Status: Early Work in Progress**

A Julia package for symbolic matrix diagonalization using closed-form root solvers.

## Vision

The ultimate goal here is to build a practical symbolic eigenvalue solver that goes beyond simple 3x3 matrices. Closed-form solutions only exist for polynomials up to degree 4 (Abel-Ruffini theorem); this means we can't solve general 5×5+ matrices symbolically. However, many real-world matrices have exploitable structure. If we can automatically detect and exploit these structures, we can potentially solve much larger symbolic problems.

**What we're building**:

- Automatic structure detection (block-diagonal, persymmetric, special patterns)
- Library of special solvable patterns (tridiagonal, circulant, Toeplitz, etc.)
- Tools to discover new solvable patterns
- Clean integration with Symbolics.jl ecosystem

**Current state**: Very basic. We have rudimentary block-diagonal detection and one special 5×5 pattern. Most of the vision is unimplemented.

## Quick Start

```julia
using Symbolics, SymbolicDiagonalization, LinearAlgebra

@variables a b c
mat = [a 1 0; 0 b 1; 0 0 c]

E = eigen(mat)        # Eigen object with .values and .vectors
λ = eigvals(mat)      # eigenvalues only (faster)
```

## What Works

### Closed-Form Root Solvers (1-4 variables)

- Linear, quadratic, cubic (Cardano), quartic (Ferrari)
- Characteristic polynomial via Bareiss determinant
- Works but produces huge expressions for fully symbolic 4×4

### Basic Structure Detection (Work in Progress)

- **Block-diagonal**: Detects and splits into independent subproblems
- **Persymmetric**: Splits certain symmetric-about-antidiagonal matrices
- **Special 5×5 tridiagonal**: Specific patterns `[b,d,b,b]` with known eigenvalues

## Installation

From package root:

```julia
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Examples

### Block-Diagonal (Works Well)

```julia
@variables a b c d
mat = [a  b  0  0;
       b  a  0  0;
       0  0  c  d;
       0  0  d  c]

vals = eigvals(mat)  # [a+b, a-b, c+d, c-d]
```

### Special 5×5 Pattern (Experimental)

```julia
@variables a b d
mat = [a  b  0  0  0;
       b  a  d  0  0;
       0  d  a  b  0;
       0  0  b  a  b;
       0  0  0  b  a]

vals = eigvals(mat)  # a ± √(2b² + d²), a ± b, a
```

## API

### LinearAlgebra Interface

- `eigen(A; kwargs...)` - Returns `Eigen` object with `.values` and `.vectors`
- `eigvals(A; kwargs...)` - Eigenvalues only (faster)

### Direct API

- `symbolic_eigenvalues(A; kwargs...)` - Returns `(values, poly, λ)`
- `symbolic_eigenpairs(A; kwargs...)` - Returns `[(λ, [v₁, v₂, ...]), ...]`
- `symbolic_diagonalize(A; kwargs...)` - Returns `(P, D, pairs)`

### Options

- `structure` - Hint: `:auto`, `:hermitian`, `:symmetric`, `:unitary`
- `expand` - Expand polynomials (default: `true`)
- `complexity_threshold` - Warn if too many variables (default: 5)
- `timeout` - Max computation time in seconds (default: 300)
- `max_terms` - Expression complexity limit (default: 10000)

## Current Limitations

- **No general 5×5+ solver**: Only works if structure detected
- **Structure detection is basic**: Misses most exploitable patterns
- **Huge expressions**: Fully symbolic 4×4 produces ~13.5 MB symbols per eigenvalue
- **No simplification**: Results are often not in simplest form
- **Fragile**: Many edge cases not handled

## What Needs Work

### High Priority

- [ ] Robust structure detection algorithms
- [ ] More special patterns (circulant, Toeplitz, Hankel, tridiagonal families)
- [ ] Better expression simplification
- [ ] Comprehensive testing of edge cases

### Medium Priority

- [ ] Eigenvalue multiplicity handling
- [ ] Symbolic condition number estimation
- [ ] Integration with numerical fallback
- [ ] Documentation of all special patterns

### Future Ideas

- [ ] Machine learning for pattern recognition
- [ ] User-defined pattern libraries
- [ ] Symbolic perturbation theory
- [ ] Parallel processing for large block systems

## Development

### Testing

```julia
julia --project -e 'using Pkg; Pkg.test()'
```

### Exploring Patterns

See [examples/explore_patterns.jl](examples/explore_patterns.jl) for tools to discover new solvable patterns.

See [docs/notes/](docs/notes/) for documentation of discovered patterns.

### Documentation

```julia
julia --project=docs docs/make.jl
```

## Contributing

This is experimental software. Contributions, ideas, and feedback very welcome:

- Implementing known special patterns from linear algebra literature
- Improving structure detection algorithms
- Finding new solvable patterns
- Performance optimization
- Better testing

## License

See LICENSE file for details.
